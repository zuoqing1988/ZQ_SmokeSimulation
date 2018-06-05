#include "ZQ_OpticalFlow3DCuda.h"

using namespace ZQ;

float ZQ_OpticalFlow3DCuda::Coarse2Fine_HS_L2(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
											  const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid3DCuda GPyramid1;
	ZQ_GaussianPyramid3DCuda GPyramid2;
	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");
	float tmp_cuda_cost_time = 0;
	float ratio = 
		GPyramid1.ConstructPyramid(Im1,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GPyramid2.ConstructPyramid(Im2,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	if(opt.displayRunningInfo)
		printf("done!\n");

	ZQ_DImage3D<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();
		int depth = GPyramid1.Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,depth,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height,depth);
			v.allocate(width,height,depth);
			w.allocate(width,height,depth);
		}
		else
		{
			u.imresize(width,height,depth);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height,depth);
			v.Multiplywith(1.0/ratio);
			w.imresize(width,height,depth);
			w.Multiplywith(1.0/ratio);
		}

		cuda_cost_time += 
			ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,depth,nChannels,
			opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);		
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_HS_L1(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
											  const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid3DCuda GPyramid1;
	ZQ_GaussianPyramid3DCuda GPyramid2;
	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");
	float tmp_cuda_cost_time = 0;
	float ratio = 
		GPyramid1.ConstructPyramid(Im1,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GPyramid2.ConstructPyramid(Im2,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	if(opt.displayRunningInfo)
		printf("done!\n");

	ZQ_DImage3D<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();
		int depth = GPyramid1.Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,depth,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height,depth);
			v.allocate(width,height,depth);
			w.allocate(width,height,depth);
		}
		else
		{
			u.imresize(width,height,depth);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height,depth);
			v.Multiplywith(1.0/ratio);
			w.imresize(width,height,depth);
			w.Multiplywith(1.0/ratio);
		}

		cuda_cost_time += 
			ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L1(u.data(),v.data(),w.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,depth,nChannels,
			opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_HS_DL1(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
											   const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid3DCuda GPyramid1;
	ZQ_GaussianPyramid3DCuda GPyramid2;
	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");
	float tmp_cuda_cost_time = 0;
	float ratio = 
		GPyramid1.ConstructPyramid(Im1,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GPyramid2.ConstructPyramid(Im2,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	if(opt.displayRunningInfo)
		printf("done!\n");

	ZQ_DImage3D<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();
		int depth = GPyramid1.Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,depth,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height,depth);
			v.allocate(width,height,depth);
			w.allocate(width,height,depth);
		}
		else
		{
			u.imresize(width,height,depth);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height,depth);
			v.Multiplywith(1.0/ratio);
			w.imresize(width,height,depth);
			w.Multiplywith(1.0/ratio);
		}


		cuda_cost_time += 
			ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,depth,nChannels,
			opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_ADMM_L2(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
												const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid3DCuda GPyramid1;
	ZQ_GaussianPyramid3DCuda GPyramid2;
	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");
	float tmp_cuda_cost_time = 0;
	float ratio = 
		GPyramid1.ConstructPyramid(Im1,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GPyramid2.ConstructPyramid(Im2,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	if(opt.displayRunningInfo)
		printf("done!\n");

	ZQ_DImage3D<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();
		int depth = GPyramid1.Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,depth,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height,depth);
			v.allocate(width,height,depth);
			w.allocate(width,height,depth);

			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,depth,nChannels,
				opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);

		}
		else
		{
			u.imresize(width,height,depth);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height,depth);
			v.Multiplywith(1.0/ratio);
			w.imresize(width,height,depth);
			w.Multiplywith(1.0/ratio);
		}

		for(int alt = 0; alt < opt.nAlterations;alt++)
		{
			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u.data(),v.data(),w.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,depth,nChannels,
				opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_ADMM_DL1(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
												 const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid3DCuda GPyramid1;
	ZQ_GaussianPyramid3DCuda GPyramid2;
	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");
	float tmp_cuda_cost_time = 0;
	float ratio = 
		GPyramid1.ConstructPyramid(Im1,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GPyramid2.ConstructPyramid(Im2,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	if(opt.displayRunningInfo)
		printf("done!\n");

	ZQ_DImage3D<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();
		int depth = GPyramid1.Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,depth,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height,depth);
			v.allocate(width,height,depth);
			w.allocate(width,height,depth);

			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,depth,nChannels,
				opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);

		}
		else
		{
			u.imresize(width,height,depth);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height,depth);
			v.Multiplywith(1.0/ratio);
			w.imresize(width,height,depth);
			w.Multiplywith(1.0/ratio);
		}

		for(int alt = 0; alt < opt.nAlterations;alt++)
		{
			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,depth,nChannels,
				opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
				opt.nSORIterations,opt.nPoissonIterations);
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Inc_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													  std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}


		//OneResolution_OneDir_Inc_L2(u,v,w,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Inc_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													   std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}
	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}


		//OneResolution_OneDir_Inc_DL1(u,v,w,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Dec_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													  std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}


		//OneResolution_OneDir_Dec_L2(u,v,w,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0;i--)
			{
				if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Dec_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													   std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}


		//OneResolution_OneDir_Dec_DL1(u,v,w,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time +=
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0;i--)
			{
				if(i == vel_num-1)
				{
					cuda_cost_time +=
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time +=
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Inc_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													  std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}

		//OneResolution_TwoDir_Inc_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}

			if(!opt.isReflect)
				continue;

			for(int i = vel_num-1;i >= 0 ;i--)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}



float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Inc_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													   std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}


	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}

		//OneResolution_TwoDir_Inc_DL1(u,v,w,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}

			if(!opt.isReflect)
				continue;

			for(int i = vel_num-1;i >= 0 ;i--)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Dec_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													  std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}

		//OneResolution_TwoDir_Dec_L2(u,v,w,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0 ;i--)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}


			if(!opt.isReflect)
				continue;

			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}



float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Dec_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
													   std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(w.size() != image_num-1)
	{
		w.clear();
		ZQ_DImage3D<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			w.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage3D<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid3DCuda> GPyramids(image_num);

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage3D<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();
		int depth = GPyramids[0].Image(k).depth();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height,depth);
				v[i].allocate(width,height,depth);
				w[i].allocate(width,height,depth);
				warpIm[i].allocate(width,height,depth,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height,depth);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height,depth);
				v[i].Multiplywith(1.0/ratio);
				w[i].imresize(width,height,depth);
				w[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,depth,nChannels);
			}
		}

		//OneResolution_TwoDir_Dec_DL1(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
			}
			break;
		}

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0 ;i--)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}

			if(!opt.isReflect)
				continue;

			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u[i].data(),v[i].data(),w[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),
						u[i-1].data(),v[i-1].data(),w[i-1].data(),u[i+1].data(),v[i+1].data(),w[i+1].data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Inc_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;



	/*   Handle images from low resolution to high resolution     */

	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         L2 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}

		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_OneDir_Inc_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}


		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					buffer->GetUVW(u,v,w,k,i);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						width,height,depth,nChannels,opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}
		}
	}

	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Inc_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;



	/*   Handle images from low resolution to high resolution     */

	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         DL1 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}

		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_OneDir_Inc_DL1(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}


		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					buffer->GetUVW(u,v,w,k,i);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						width,height,depth,nChannels,opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				
				/************  output the final result     **************/
				if(k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}
		}
	}

	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Dec_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;



	/*   Handle images from low resolution to high resolution     */

	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         L2 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}

		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_OneDir_Dec_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}


		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0 ;i--)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == vel_num-1)
				{
					buffer->GetUVW(u,v,w,k,i);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						width,height,depth,nChannels,opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else 
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}
		}
	}

	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);
	

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_OneDir_Dec_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;



	/*   Handle images from low resolution to high resolution     */

	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         DL1 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}

		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_OneDir_Dec_DL1(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}


		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0 ;i--)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				
				if(i == vel_num-1)
				{
					buffer->GetUVW(u,v,w,k,i);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						width,height,depth,nChannels,opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}
		}
	}

	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);

	return cuda_cost_time;
}



float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Inc_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;



	/*   Handle images from low resolution to high resolution     */

	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         L2 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}

		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_TwoDir_Inc_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}


		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(!opt.isReflect && k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}

			if(!opt.isReflect)
				continue;

			for(int i = vel_num-1;i >= 0 ;i--)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}	
		}
	}


	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Inc_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;



	/*   Handle images from low resolution to high resolution     */

	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         DL1 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}

		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_TwoDir_Inc_DL1(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}


		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = 0;i < vel_num;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(!opt.isReflect && k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}

			
			
			if(!opt.isReflect)
				continue;

			for(int i = vel_num-1;i >= 0 ;i--)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				
				/************  output the final result     **************/
				if(k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}	
		}
	}

	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Dec_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;



	/*   Handle images from low resolution to high resolution     */

	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         L2 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}

		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_TwoDir_Dec_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_L2(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}


		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0 ;i--)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(!opt.isReflect && k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}

			

			if(!opt.isReflect)
				continue;

			for(int i = 0;i < vel_num;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if( k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}
		}
	}

	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);

	return cuda_cost_time;
}


float ZQ_OpticalFlow3DCuda::Coarse2Fine_TwoDir_Dec_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt)
{
	typedef ZQ_DImage3D<float> DImage3D;

	int nLevels = 0;
	float cuda_cost_time = 0;
	float tmp_cuda_cost_time = 0;
	if(buffer == 0 || !buffer->BuildPyramids(opt.ratioForPyramid,opt.minWidthForPyramid,nLevels,tmp_cuda_cost_time))
	{
		printf("Build Pyramids fail\n");
		return 0;
	}
	cuda_cost_time += tmp_cuda_cost_time;

	int image_num = buffer->GetImageNum();
	DImage3D im1,im2;
	DImage3D Im1,Im2,warpIm2;
	int width,height,depth,nChannels;
	DImage3D u,v,w;

	

	/*   Handle images from low resolution to high resolution     */
	
	for(int k = nLevels-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		/*************         DL1 for the start  resolution    ***************/

		if(k == nLevels-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				u.allocate(width,height,depth,1);
				v.allocate(width,height,depth,1);
				w.allocate(width,height,depth,1);

				cuda_cost_time += ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);

				buffer->WriteUVW(u,v,w,k,i);
			}
		}
		
		/**********    Upsample from the coarser resolution      *************/
		if(k < nLevels-1)
		{
			buffer->GetImage(im1,k,0);
			width = im1.width();
			height = im1.height();
			depth = im1.depth();

			for(int t = 0;t < image_num-1;t++)
			{
				buffer->GetUVW(u,v,w,k+1,t);
				u.imresize(width,height,depth);
				u.Multiplywith(1.0/opt.ratioForPyramid);
				v.imresize(width,height,depth);
				v.Multiplywith(1.0/opt.ratioForPyramid);
				w.imresize(width,height,depth);
				w.Multiplywith(1.0/opt.ratioForPyramid);
				buffer->WriteUVW(u,v,w,k,t);
			}
		}

		//OneResolution_TwoDir_Dec_DL1(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		for(int i = 0;i < vel_num;i++)
		{
			if(opt.initType == ZQ_OpticalFlowOptions::NONE_AS_INIT)
				break;

			buffer->GetImage(im1,k,i);
			buffer->GetImage(im2,k,i+1);

			im2feature(Im1,im1,isSmooth,opt);
			im2feature(Im2,im2,isSmooth,opt);
			warpIm2.allocate(Im2);

			width = Im1.width();
			height = Im1.height();
			depth = Im1.depth();
			nChannels = Im1.nchannels();

			DImage3D u,v,w;
			buffer->GetUVW(u,v,w,k,i);

			switch (opt.initType)
			{
			case ZQ_OpticalFlowOptions::L2_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
				break;
			case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),width,height,depth,nChannels,
					opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				break;
			}

			buffer->WriteUVW(u,v,w,k,i);
		}
		

		for(int alt_it = 0;alt_it < opt.nAlterations;alt_it++)
		{
			for(int i = vel_num-1;i >= 0 ;i--)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(!opt.isReflect && k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}


			if(!opt.isReflect)
				continue;

			for(int i = 0;i < vel_num;i++)
			{
				buffer->GetImage(im1,k,i);
				buffer->GetImage(im2,k,i+1);

				im2feature(Im1,im1,isSmooth,opt);
				im2feature(Im2,im2,isSmooth,opt);
				warpIm2.allocate(Im2);

				width = Im1.width();
				height = Im1.height();
				depth = Im1.depth();
				nChannels = Im1.nchannels();

				if(i == 0)
				{
					DImage3D u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_First(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else if(i == vel_num-1)
				{
					DImage3D u_pre,v_pre,w_pre;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Last(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,
						opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}
				else
				{
					DImage3D u_pre,v_pre,w_pre,u_nex,v_nex,w_nex;
					buffer->GetUVW(u,v,w,k,i);
					buffer->GetUVW(u_pre,v_pre,w_pre,k,i-1);
					buffer->GetUVW(u_nex,v_nex,w_nex,k,i+1);

					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow3D::OpticalFlow3D_ADMM_DL1_Middle(u.data(),v.data(),w.data(),warpIm2.data(),Im1.data(),Im2.data(),
						u_pre.data(),v_pre.data(),w_pre.data(),u_nex.data(),v_nex.data(),w_nex.data(),width,height,depth,nChannels,opt.alpha,opt.beta,opt.gamma,
						opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
					buffer->WriteUVW(u,v,w,k,i);
				}

				/************  output the final result     **************/
				if(k == 0 && alt_it == opt.nAlterations-1)
				{
					DImage3D flow,warp,other;
					flow.assemble(u,v,w);
					warpIm2.separate(1,warp,other);
					buffer->WriteFinalResult(flow,warp,i);
				}
			}
		}
	}

	buffer->ClearPyramids(nLevels);
	buffer->ClearIntermediateUVW(nLevels);

	return cuda_cost_time;
}

