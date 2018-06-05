#include "ZQ_OpticalFlow2DCuda.h"
//#include "ZQ_ImageIO.h"

using namespace ZQ;

float ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L2(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid2DCuda GPyramid1;
	ZQ_GaussianPyramid2DCuda GPyramid2;
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

	ZQ_DImage<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height);
			v.allocate(width,height);
		}
		else
		{
			u.imresize(width,height);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height);
			v.Multiplywith(1.0/ratio);
		}

		
		cuda_cost_time += 
			ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,nChannels,
			opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);		

		/*IplImage* flowimg =  ZQ_ImageIO::SaveFlowToColorImage(u, v, false, 0, 16, 0, false);
		char name[200];
		sprintf(name, "flow_%d.jpg", k);
		cvSaveImage(name, flowimg);
		cvReleaseImage(&flowimg);
		sprintf(name, "warp_%d.jpg", k);
		ZQ_ImageIO::saveImage(warpIm2, name);*/
	}

	return cuda_cost_time;
}


float ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L2_Occupy(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_DImage<float>& Occupy, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid2DCuda GPyramid1;
	ZQ_GaussianPyramid2DCuda GPyramid2;
	ZQ_GaussianPyramid2DCuda GOccupy;
	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");
	float tmp_cuda_cost_time = 0;
	float ratio = 
		GPyramid1.ConstructPyramid(Im1,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GPyramid2.ConstructPyramid(Im2,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GOccupy.ConstructPyramid(Occupy,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	if(opt.displayRunningInfo)
		printf("done!\n");

	ZQ_DImage<float> Image1,Image2,Occ;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);
		Occ = GOccupy.Image(k);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height);
			v.allocate(width,height);
		}
		else
		{
			u.imresize(width,height);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height);
			v.Multiplywith(1.0/ratio);
		}

		cuda_cost_time += 
			ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2_Occupy(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),Occ.data(),width,height,nChannels,
			opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);		
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L1(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid2DCuda GPyramid1;
	ZQ_GaussianPyramid2DCuda GPyramid2;
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

	ZQ_DImage<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height);
			v.allocate(width,height);
		}
		else
		{
			u.imresize(width,height);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height);
			v.Multiplywith(1.0/ratio);
			
		}

		cuda_cost_time += 
			ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L1(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,nChannels,
			opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_DL1(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid2DCuda GPyramid1;
	ZQ_GaussianPyramid2DCuda GPyramid2;
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

	ZQ_DImage<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height);
			v.allocate(width,height);
		}
		else
		{
			u.imresize(width,height);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height);
			v.Multiplywith(1.0/ratio);
		}


		cuda_cost_time += 
			ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,nChannels,
			opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_ADMM_L2(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid2DCuda GPyramid1;
	ZQ_GaussianPyramid2DCuda GPyramid2;
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

	ZQ_DImage<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height);
			v.allocate(width,height);

			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,nChannels,
				opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);

		}
		else
		{
			u.imresize(width,height);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height);
			v.Multiplywith(1.0/ratio);
		}

		for(int alt = 0; alt < opt.nAlterations;alt++)
		{
			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,nChannels,
				opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
		}
	}

	return cuda_cost_time;
}


float ZQ_OpticalFlow2DCuda::Coarse2Fine_ADMM_L2_Occupy(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_DImage<float>& Occupy, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid2DCuda GPyramid1;
	ZQ_GaussianPyramid2DCuda GPyramid2;
	ZQ_GaussianPyramid2DCuda GOccupy;
	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");
	float tmp_cuda_cost_time = 0;
	float ratio = 
		GPyramid1.ConstructPyramid(Im1,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GPyramid2.ConstructPyramid(Im2,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	GOccupy.ConstructPyramid(Occupy,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;
	if(opt.displayRunningInfo)
		printf("done!\n");

	ZQ_DImage<float> Image1,Image2,Occ;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);
		Occ = GOccupy.Image(k);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height);
			v.allocate(width,height);

			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2_Occupy(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),Occ.data(),width,height,nChannels,
				opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);

		}
		else
		{
			u.imresize(width,height);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height);
			v.Multiplywith(1.0/ratio);
		}

		for(int alt = 0; alt < opt.nAlterations;alt++)
		{
			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Occupy(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),Occ.data(),width,height,nChannels,
				opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_ADMM_DL1(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	ZQ_GaussianPyramid2DCuda GPyramid1;
	ZQ_GaussianPyramid2DCuda GPyramid2;
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

	ZQ_DImage<float> Image1,Image2;

	for(int k = GPyramid1.nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramid1.Image(k).width();
		int height = GPyramid1.Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		im2feature(Image1,GPyramid1.Image(k),isSmooth,opt);
		im2feature(Image2,GPyramid2.Image(k),isSmooth,opt);

		int nChannels = Image1.nchannels();
		warpIm2.allocate(width,height,nChannels);

		if(k == GPyramid1.nlevels()-1)
		{
			u.allocate(width,height);
			v.allocate(width,height);

			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,nChannels,
				opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);

		}
		else
		{
			u.imresize(width,height);
			u.Multiplywith(1.0/ratio);
			v.imresize(width,height);
			v.Multiplywith(1.0/ratio);
		}

		for(int alt = 0; alt < opt.nAlterations;alt++)
		{
			cuda_cost_time += 
				ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1(u.data(),v.data(),warpIm2.data(),Image1.data(),Image2.data(),width,height,nChannels,
				opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
				opt.nSORIterations,opt.nPoissonIterations);
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Inc_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
			}
		}


		//OneResolution_OneDir_Inc_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,
						opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Inc_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
			}
		}


		//OneResolution_OneDir_Inc_DL1(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Dec_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
			}
		}


		//OneResolution_OneDir_Dec_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nSORIterations,
						opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Dec_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;

	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
			}
		}


		//OneResolution_OneDir_Dec_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time +=
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
						opt.alpha,opt.beta,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time +=
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}

			}
		}
	}

	return cuda_cost_time;
}

float ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Inc_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
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
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}


float ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Inc_L2_Occupy(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_DImage<float>& Occupy, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);
	ZQ_GaussianPyramid2DCuda GOccupy;

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}

	GOccupy.ConstructPyramid(Occupy,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage<float>> Images(image_num);
	ZQ_DImage<float> Occ;

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);
		Occ = GOccupy.Image(k);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
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
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}



float ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Inc_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
			}
		}

		//OneResolution_TwoDir_Inc_DL1(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nInnerFixedPointIterations,opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}


float ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Dec_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
			}
		}

		//OneResolution_TwoDir_Dec_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}


			if(!opt.isReflect)
				continue;

			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}


float ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Dec_L2_Occupy(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_DImage<float>& Occupy, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);
	ZQ_GaussianPyramid2DCuda GOccupy;

	if(opt.displayRunningInfo)
		printf("Constructing pyramid...");

	float tmp_cuda_cost_time = 0;
	float ratio = 0;
	for(int i = 0;i < image_num;i++)
	{
		ratio = GPyramids[i].ConstructPyramid(Im[i],tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
		cuda_cost_time += tmp_cuda_cost_time;
	}
	GOccupy.ConstructPyramid(Occupy,tmp_cuda_cost_time,opt.ratioForPyramid,opt.minWidthForPyramid);
	cuda_cost_time += tmp_cuda_cost_time;

	if(opt.displayRunningInfo)
		printf("done!\n");

	/*   Handle images from low resolution to high resolution     */
	std::vector<ZQ_DImage<float>> Images(image_num);
	ZQ_DImage<float> Occ;

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

		bool isSmooth = true;
		if(k == 0)
			isSmooth = false;

		for(int i = 0;i < image_num;i++)
			im2feature(Images[i],GPyramids[i].Image(k),isSmooth,opt);
		Occ = GOccupy.Image(k);

		int nChannels = Images[0].nchannels();

		if(k == GPyramids[0].nlevels()-1)
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
			}
		}

		//OneResolution_TwoDir_Dec_L2(u,v,warpIm,Images,opt);

		int vel_num = image_num-1;

		switch (opt.initType)
		{
		case ZQ_OpticalFlowOptions::NONE_AS_INIT:
			break;
		case ZQ_OpticalFlowOptions::L2_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_L2_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}


			if(!opt.isReflect)
				continue;

			for(int i = 0;i < vel_num;i++)
			{
				if(i == 0)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_First_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Last_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_Middle_Occupy(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),Occ.data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}



float ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Dec_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt)
{
	float cuda_cost_time = 0;
	int image_num = Im.size();

	if(u.size() != image_num-1)
	{
		u.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			u.push_back(tmp);
	}

	if(v.size() != image_num-1)
	{
		v.clear();
		ZQ_DImage<float> tmp;
		for(int i = 0;i < image_num-1;i++)
			v.push_back(tmp);
	}

	if(warpIm.size() != image_num-1)
	{
		warpIm.clear();
		ZQ_DImage<float> tmp;
		for (int i = 0;i < image_num-1;i++)
		{
			warpIm.push_back(tmp);
		}
	}

	/*  Constructing Pyramids */
	std::vector<ZQ_GaussianPyramid2DCuda> GPyramids(image_num);

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
	std::vector<ZQ_DImage<float>> Images(image_num);

	for(int k = GPyramids[0].nlevels()-1;k >= 0;k--)
	{
		if(opt.displayRunningInfo)
			printf("Pyramid level %d \n",k);

		int width = GPyramids[0].Image(k).width();
		int height = GPyramids[0].Image(k).height();

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
				u[i].allocate(width,height);
				v[i].allocate(width,height);
				warpIm[i].allocate(width,height,nChannels);
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
		}
		else
		{
			for(int i = 0;i < image_num-1;i++)
			{
				u[i].imresize(width,height);
				u[i].Multiplywith(1.0/ratio);
				v[i].imresize(width,height);
				v[i].Multiplywith(1.0/ratio);
				warpIm[i].allocate(width,height,nChannels);
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
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
					opt.alpha,opt.beta,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,opt.nSORIterations);
			}
			break;
		case ZQ_OpticalFlowOptions::ADMM_AS_INIT:
			for(int i = 0;i < vel_num;i++)
			{
				cuda_cost_time += 
					ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),width,height,nChannels,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,
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
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_First(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i+1].data(),v[i+1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);

				}
				else if(i == vel_num-1)
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Last(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
				else
				{
					cuda_cost_time += 
						ZQ_CUDA_OpticalFlow2D::OpticalFlow2D_ADMM_DL1_Middle(u[i].data(),v[i].data(),warpIm[i].data(),Images[i].data(),Images[i+1].data(),u[i-1].data(),v[i-1].data(),
						u[i+1].data(),v[i+1].data(),width,height,nChannels,opt.alpha,opt.beta,opt.gamma,opt.lambda,opt.nADMMIterations,opt.nOuterFixedPointIterations,opt.nInnerFixedPointIterations,
						opt.nSORIterations,opt.nAdvectFixedPointIterations,opt.nPoissonIterations);
				}
			}
		}
	}

	return cuda_cost_time;
}

