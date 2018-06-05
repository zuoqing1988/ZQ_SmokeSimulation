#include "ZQ_GaussianPyramid2DCuda.h"
#include <math.h>

using namespace ZQ;

ZQ_GaussianPyramid2DCuda::ZQ_GaussianPyramid2DCuda(void)
{
}

ZQ_GaussianPyramid2DCuda::~ZQ_GaussianPyramid2DCuda(void)
{
	ImPyramid.clear();
}

float ZQ_GaussianPyramid2DCuda::ConstructPyramid(const ZQ_DImage<float> &image, float& cuda_cost_time, const float ratio, const int minWidth)
{
	cuda_cost_time = 0;

	if(ratio > 0.9 || ratio < 0.4)
		this->ratio = 0.75;
	else
		this->ratio = ratio;

	// first decide how many levels
	nLevels = log((float)minWidth/image.width())/log(this->ratio)+1;
	ImPyramid.clear();
	for(int i = 0;i < nLevels;i++)
	{
		ZQ_DImage<float> tmp;
		ImPyramid.push_back(tmp);
	}
	ImPyramid[0].copyData(image);
	float baseSigma = (1.0/(this->ratio)-1);
	int n = log(0.25f)/log(this->ratio);
	float nSigma=baseSigma*n;
	for(int i = 1;i < nLevels;i++)
	{
		ZQ_DImage<float> foo;
		if(i <= n)
		{
			int cur_width = image.width();
			int cur_height = image.height();
			int nChannels = image.nchannels();

			float sigma = baseSigma*i;
			foo.allocate(cur_width,cur_height,nChannels);
			cuda_cost_time += 
				ZQ_CUDA_ImageProcessing2D::GaussianSmoothing2_2D(foo.data(),image.data(),sigma,sigma*3,cur_width,cur_height,nChannels);

			double cur_ratio = pow(this->ratio,i);

			int dst_width = cur_width*cur_ratio;
			int dst_height = cur_height*cur_ratio;

			ImPyramid[i].allocate(dst_width,dst_height,nChannels);

			cuda_cost_time += 
				ZQ_CUDA_ImageProcessing2D::ResizeImage2D_Bicubic(ImPyramid[i].data(),foo.data(),cur_width,cur_height,dst_width,dst_height,nChannels);
		}
		else
		{
			int cur_width = ImPyramid[i-n].width();
			int cur_height = ImPyramid[i-n].height();
			int nChannels = ImPyramid[i-n].nchannels();

			foo.allocate(cur_width,cur_height,nChannels);

			cuda_cost_time +=
				ZQ_CUDA_ImageProcessing2D::GaussianSmoothing2_2D(foo.data(),ImPyramid[i-n].data(),nSigma,nSigma*3,cur_width,cur_height,nChannels);

			int dst_width = pow(this->ratio,i)*image.width();
			int dst_height = pow(this->ratio,i)*image.height();

			ImPyramid[i].allocate(dst_width,dst_height,nChannels);

			cuda_cost_time += 
				ZQ_CUDA_ImageProcessing2D::ResizeImage2D_Bicubic(ImPyramid[i].data(),foo.data(),cur_width,cur_height,dst_width,dst_height,nChannels);
		}
	}

	return this->ratio;
}
