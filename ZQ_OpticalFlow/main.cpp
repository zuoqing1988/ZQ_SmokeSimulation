#include "ZQ_DoubleImage.h"
#include "ZQ_ImageIO.h"
#include "ZQ_OpticalFlow2DCuda.h"
#include <stdio.h>
#include <time.h>

using namespace ZQ;

typedef ZQ_DImage<float> DImage;

void main_Pair(const ZQ_OpticalFlowOptions& opt, const char* in_fold, const char* prefix, const char* suffix, const int image_num, const int base_id,const char* out_fold);

void main_Seq(const ZQ_OpticalFlowOptions& opt, const char* in_fold, const char* prefix, const char* suffix, const int image_num, const int base_id,const char* out_fold);


int main(int argc, const char** argv)
{
	clock_t t1 = clock();

	if (argc < 4)
	{
		printf(" .exe image_name1 image_name2 flow_name [paras]\n");
		return 0;
	}

	const char* image_name1 = argv[1];
	const char* image_name2 = argv[2];
	const char* flow_name = argv[3];

	ZQ_OpticalFlowOptions opt;
	if (!opt.HandleParas(argc - 4, argv + 4))
	{
		printf("option args error\n");
		return 0;
	}

	DImage im1, im2, warpIm, u, v;
	
	
	if (!ZQ_ImageIO::loadImage(im1, image_name1, 0))
	{
		printf("failed to load %s\n", image_name1);
		return 0;
	}
	if (!ZQ_ImageIO::loadImage(im2, image_name2, 0))
	{
		printf("failed to load %s\n", image_name2);
		return 0;
	}

	double max_rad = opt.maxRad;
	bool user_input = max_rad > 0;


	ZQ_CUDA_OpticalFlow2D::InitDevice2D(opt.cudaDeviceID);
	float cuda_cost_time = 0;
	switch (opt.methodType)
	{
	case ZQ_OpticalFlowOptions::METHOD_HS_L2:
		cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L2(u, v, warpIm, im1, im2, opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_HS_DL1:
		cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_DL1(u, v, warpIm, im1, im2, opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_HS_L1:
		cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L1(u, v, warpIm, im1, im2, opt);
		break;
	default:
		printf("this method is no supported!\n");
		return 0;
		break;
	}


	DImage flow;
	flow.assemble(u, v);

	flow.saveImage(flow_name);

	char buf[200];
	sprintf_s(buf, "%s.png", flow_name);
	int wheelSize = 64;
	int color_type = 0;
	IplImage* flow_img = ZQ_ImageIO::SaveFlowToColorImage(u, v, user_input, max_rad, wheelSize, color_type);


	cvSaveImage(buf, flow_img);
	cvReleaseImage(&flow_img);

	DImage warp, other;
	warpIm.separate(1, warp, other);
	sprintf_s(buf, "warp_%s", image_name2);
	ZQ_ImageIO::saveImage(warp, buf);

	clock_t t2 = clock();

	printf("cuda_cost_time / total_cost = %f/%f \n", cuda_cost_time*0.001, (t2 - t1)*0.001);
	return 0;
}

void main1(int argc, const char** argv)
{
	if(argc < 7)
	{
		printf(" .exe inputfold prefix suffix imagenum baseid outputfold [paras]\n");
		return ;
	}

	const char* input_fold = argv[1];
	const char* prefix = argv[2];
	const char* suffix = argv[3];
	int image_num = atoi(argv[4]);
	int base_id = atoi(argv[5]);
	const char* output_fold = argv[6];

	ZQ_OpticalFlowOptions opt;
	if(!opt.HandleParas(argc-7,argv+7))
	{
		printf("option args error\n");
		return;
	}

	if(opt.methodType == ZQ_OpticalFlowOptions::METHOD_HS_L2 
		|| opt.methodType == ZQ_OpticalFlowOptions::METHOD_HS_L1
		|| opt.methodType == ZQ_OpticalFlowOptions::METHOD_HS_DL1
		|| opt.methodType == ZQ_OpticalFlowOptions::METHOD_ADMM_L2
		|| opt.methodType == ZQ_OpticalFlowOptions::METHOD_ADMM_DL1
		|| opt.methodType == ZQ_OpticalFlowOptions::METHOD_HS_L2_OCCUPY
		|| opt.methodType == ZQ_OpticalFlowOptions::METHOD_ADMM_L2_OCCUPY)
	{
		ZQ_CUDA_OpticalFlow2D::InitDevice2D(opt.cudaDeviceID);
		main_Pair(opt,input_fold,prefix,suffix,image_num,base_id,output_fold);
	}
	else
	{
		ZQ_CUDA_OpticalFlow2D::InitDevice2D(opt.cudaDeviceID);
		main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,output_fold);
	}

}

void main2()
{
	const char* opt_argv[] = {
		"methodtype", "HS_L2",
		"nAltIter", "3",
		"nAdmmIter", "20",
		"nOuterIter", "10",
		"nInnerIter", "3",
		"nPoissonIter", "30",
		"nAdvectIter", "3",
		"nSORIter", "30",
		"initType","ADMM_AS_INIT",
		"alpha","0.06",
		"beta","0",
		"gamma","0.1",
		"lambda","0.01",
		"display",
		"cubic"
	};


	ZQ_OpticalFlowOptions opt;
	if(!opt.HandleParas(sizeof(opt_argv)/sizeof(char*),opt_argv))
	{
		printf("option args error\n");
		return;
	}

	const char* input_fold = "INPUT";
	const char* prefix = "par_";
	const char* suffix = "png";
	const int image_num = 4;
	const int base_id = 0;
	
	const char* HS_L2_fold = "HS_L2";
	const char* HS_DL1_fold = "HS_DL1";
	const char* HS_L1_fold = "HS_L1";
	const char* ADMM_L2_fold = "ADMM_L2";
	const char* ADMM_DL1_fold = "ADMM_DL1";
	const char* ONEDIR_INC_L2_fold = "ONEDIR_INC_L2";
	const char* ONEDIR_INC_DL1_fold = "ONEDIR_INC_DL1";
	const char* ONEDIR_DEC_L2_fold = "ONEDIR_DEC_L2";
	const char* ONEDIR_DEC_DL1_fold = "ONEDIR_DEC_DL1";
	const char* TWODIR_INC_L2_fold = "TWODIR_INC_L2";
	const char* TWODIR_INC_DL1_fold = "TWODIR_INC_DL1";
	const char* TWODIR_DEC_L2_fold = "TWODIR_DEC_L2";
	const char* TWODIR_DEC_DL1_fold = "TWODIR_DEC_DL1";



	opt.methodType = ZQ_OpticalFlowOptions::METHOD_HS_L2;
	opt.alpha = 0.06;
	opt.beta = 0.001;
	opt.nOuterFixedPointIterations = 10;
	main_Pair(opt,input_fold,prefix,suffix,image_num,base_id,HS_L2_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_HS_DL1;
	opt.alpha = 0.2;
	opt.beta = 0.001;
	opt.nOuterFixedPointIterations = 10;
	opt.nInnerFixedPointIterations = 2;
	main_Pair(opt,input_fold,prefix,suffix,image_num,base_id,HS_DL1_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_HS_L1;
	opt.alpha = 0.012;
	opt.beta = 0.001;
	opt.nOuterFixedPointIterations = 10;
	opt.nInnerFixedPointIterations = 2;
	main_Pair(opt,input_fold,prefix,suffix,image_num,base_id,HS_L1_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_ADMM_L2;
	opt.alpha = 0.06;
	opt.beta = 0.001;
	opt.nADMMIterations = 20;
	opt.nOuterFixedPointIterations = 5;
	main_Pair(opt,input_fold,prefix,suffix,image_num,base_id,ADMM_L2_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_ADMM_DL1;
	opt.alpha = 0.2;
	opt.beta = 0.001;
	opt.nADMMIterations = 20;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 2;
	main_Pair(opt,input_fold,prefix,suffix,image_num,base_id,ADMM_DL1_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_ONEDIR_INC_L2;
	opt.alpha = 0.06;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,ONEDIR_INC_L2_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_ONEDIR_INC_DL1;
	opt.alpha = 0.2;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,ONEDIR_INC_DL1_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_ONEDIR_DEC_L2;
	opt.alpha = 0.06;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,ONEDIR_DEC_L2_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_ONEDIR_DEC_DL1;
	opt.alpha = 0.2;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,ONEDIR_DEC_DL1_fold);


	opt.methodType = ZQ_OpticalFlowOptions::METHOD_TWODIR_INC_L2;
	opt.alpha = 0.06;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,TWODIR_INC_L2_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_TWODIR_INC_DL1;
	opt.alpha = 0.2;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,TWODIR_INC_DL1_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_TWODIR_DEC_L2;
	opt.alpha = 0.06;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,TWODIR_DEC_L2_fold);

	opt.methodType = ZQ_OpticalFlowOptions::METHOD_TWODIR_DEC_DL1;
	opt.alpha = 0.2;
	opt.nAlterations = 3;
	opt.nOuterFixedPointIterations = 5;
	opt.nInnerFixedPointIterations = 1;
	main_Seq(opt,input_fold,prefix,suffix,image_num,base_id,TWODIR_DEC_DL1_fold);
}


void main_Pair(const ZQ_OpticalFlowOptions& opt, const char* in_fold, const char* prefix, const char* suffix, const int image_num, const int base_id,const char* out_fold)
{
	if(image_num < 2)
	{
		printf("error: image num < 2\n");
		return ;
	}

	char buf[2000];

	std::vector<DImage> Images;
	for(int i = 0;i < image_num;i++)
	{
		DImage Im;
		sprintf_s(buf,"%s\\%s%d.%s",in_fold,prefix,i+base_id,suffix);
		bool flag = ZQ_ImageIO::loadImage(Im,buf,0);
		if(!flag)
		{
			printf("fail to load image %s\n",buf);
			return ;
		}
		Images.push_back(Im);
	}

	double max_rad = opt.maxRad;
	bool user_input = max_rad > 0;

	std::vector<DImage> u(image_num-1),v(image_num-1),warpIm(image_num-1);

	DImage mask;

	float cuda_cost_time = 0;
	clock_t start = clock();
	for(int i = 0;i < image_num-1;i++)
	{
		switch(opt.methodType)
		{
		case ZQ_OpticalFlowOptions::METHOD_HS_L2:
			cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L2(u[i],v[i],warpIm[i],Images[i],Images[i+1],opt);
			break;
		case ZQ_OpticalFlowOptions::METHOD_HS_DL1:
			cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_DL1(u[i],v[i],warpIm[i],Images[i],Images[i+1],opt);
			break;
		case ZQ_OpticalFlowOptions::METHOD_HS_L1:
			cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L1(u[i],v[i],warpIm[i],Images[i],Images[i+1],opt);
			break;
		case ZQ_OpticalFlowOptions::METHOD_ADMM_L2:
			cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_ADMM_L2(u[i],v[i],warpIm[i],Images[i],Images[i+1],opt);
			break;
		case ZQ_OpticalFlowOptions::METHOD_ADMM_DL1:
			cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_ADMM_DL1(u[i],v[i],warpIm[i],Images[i],Images[i+1],opt);
			break;
		case ZQ_OpticalFlowOptions::METHOD_HS_L2_OCCUPY:
			{
				if(i == 0)
				{
					if(!ZQ_ImageIO::loadImage(mask,opt.maskFile,0))
					{
						printf("failed to load %s\n",opt.maskFile);
						return ;
					}
				}
				cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_HS_L2_Occupy(u[i],v[i],warpIm[i],Images[i],Images[i+1],mask,opt);
			}
			break;

		case ZQ_OpticalFlowOptions::METHOD_ADMM_L2_OCCUPY:
			{
				if(i == 0)
				{
					if(!ZQ_ImageIO::loadImage(mask,opt.maskFile,0))
					{
						printf("failed to load %s\n",opt.maskFile);
						return ;
					}
				}
				cuda_cost_time += ZQ_OpticalFlow2DCuda::Coarse2Fine_ADMM_L2_Occupy(u[i],v[i],warpIm[i],Images[i],Images[i+1],mask,opt);
			}
			break;
		}

	}	

	clock_t end = clock();

	for(int i = 0;i < image_num-1;i++)
	{

		DImage flow;
		flow.assemble(u[i],v[i]);
		sprintf_s(buf,"%s\\flow%d_%d.di2",out_fold,i+base_id,i+1+base_id);
		flow.saveImage(buf);


		int wheelSize = 64;
		int color_type = 1;
		IplImage* flow_img =  ZQ_ImageIO::SaveFlowToColorImage(u[i],v[i],user_input,max_rad,wheelSize,color_type);

		sprintf_s(buf,"%s\\flow%d_%d.%s",out_fold,i+base_id,i+1+base_id,suffix);
		cvSaveImage(buf,flow_img);
		cvReleaseImage(&flow_img);

		DImage warp,other;
		warpIm[i].separate(1,warp,other);
		sprintf_s(buf,"%s\\warp%d_%d.%s",out_fold,i+1+base_id,i+base_id,suffix);
		ZQ_ImageIO::saveImage(warp,buf);

	}

	printf("cuda_cost_time / total_cost = %f/%f \n",cuda_cost_time*0.001,(end-start)*0.001 );

	return ;
}


void main_Seq(const ZQ_OpticalFlowOptions& opt, const char* in_fold, const char* prefix, const char* suffix, const int image_num, const int base_id,const char* out_fold)
{
	if(image_num < 3)
	{
		printf("error: image num < 3\n");
		return ;
	}
	
	char buf[2000];

	std::vector<DImage> Images;
	for(int i = 0;i < image_num;i++)
	{
		DImage Im;
		sprintf_s(buf,"%s\\%s%d.%s",in_fold,prefix,i+base_id,suffix);
		bool flag = ZQ_ImageIO::loadImage(Im,buf,0);
		if(!flag)
		{
			printf("fail to load image %s\n",buf);
			return ;
		}
		Images.push_back(Im);

	}
	
	
	double max_rad = opt.maxRad;
	bool user_input = max_rad > 0;

	std::vector<DImage> u,v,warpIm;
	DImage mask;

	float cuda_cost_time = 0;
	clock_t start = clock();

	switch(opt.methodType)
	{

	case ZQ_OpticalFlowOptions::METHOD_ONEDIR_INC_L2:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Inc_L2(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_ONEDIR_INC_DL1:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Inc_DL1(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_ONEDIR_DEC_L2:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Dec_L2(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_ONEDIR_DEC_DL1:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_OneDir_Dec_DL1(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_TWODIR_INC_L2:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Inc_L2(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_TWODIR_INC_DL1:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Inc_DL1(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_TWODIR_DEC_L2:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Dec_L2(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_TWODIR_DEC_DL1:
		cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Dec_DL1(u,v,warpIm,Images,opt);
		break;
	case ZQ_OpticalFlowOptions::METHOD_TWODIR_INC_L2_OCCUPY:
		{
			if(!ZQ_ImageIO::loadImage(mask,opt.maskFile,0))
			{
				printf("failed to load %s\n",opt.maskFile);
				return;
			}
			cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Inc_L2_Occupy(u,v,warpIm,Images,mask,opt);
		}
		break;
	case ZQ_OpticalFlowOptions::METHOD_TWODIR_DEC_L2_OCCUPY:
		{
			if(!ZQ_ImageIO::loadImage(mask,opt.maskFile,0))
			{
				printf("failed to load %s\n",opt.maskFile);
				return;
			}
			cuda_cost_time = ZQ_OpticalFlow2DCuda::Coarse2Fine_TwoDir_Dec_L2_Occupy(u,v,warpIm,Images,mask,opt);
		}
		break;
	}
	
	clock_t end = clock();

	float total_time = 0.001*(end-start);
	
	for(int i = 0;i < image_num-1;i++)
	{
		
		
		DImage flow;
		flow.assemble(u[i],v[i]);
		sprintf_s(buf,"%s\\flow%d_%d.di2",out_fold,i+base_id,i+1+base_id);
		flow.saveImage(buf);
		
			
		int wheelSize = 64;
		int color_type = 1;
		IplImage* flow_img =  ZQ_ImageIO::SaveFlowToColorImage(u[i],v[i],user_input,max_rad,wheelSize,color_type);

		sprintf_s(buf,"%s\\flow%d_%d.%s",out_fold,i+base_id,i+1+base_id,suffix);
		cvSaveImage(buf,flow_img);
		cvReleaseImage(&flow_img);
		
		DImage warp,other;
		warpIm[i].separate(1,warp,other);
		sprintf_s(buf,"%s\\warp%d_%d.%s",out_fold,i+1+base_id,i+base_id,suffix);
		ZQ_ImageIO::saveImage(warp,buf);

	}

	printf("cuda_cost_time / total_time = %f/%f\n", cuda_cost_time*0.001,total_time);
	return ;
}
