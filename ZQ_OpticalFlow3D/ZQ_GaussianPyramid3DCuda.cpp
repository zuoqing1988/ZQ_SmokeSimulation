#include "ZQ_GaussianPyramid3DCuda.h"
#include <math.h>

using namespace ZQ;

ZQ_GaussianPyramid3DCuda::ZQ_GaussianPyramid3DCuda(void)
{
}

ZQ_GaussianPyramid3DCuda::~ZQ_GaussianPyramid3DCuda(void)
{
	ImPyramid.clear();
}

float ZQ_GaussianPyramid3DCuda::ConstructPyramid(const ZQ_DImage3D<float> &image, float& cuda_cost_time, const float ratio, const int minWidth)
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
		ZQ_DImage3D<float> tmp;
		ImPyramid.push_back(tmp);
	}
	ImPyramid[0].copyData(image);
	float baseSigma = (1.0/(this->ratio)-1);
	int n = log(0.25f)/log(this->ratio);
	float nSigma=baseSigma*n;
	for(int i = 1;i < nLevels;i++)
	{
		ZQ_DImage3D<float> foo;
		if(i <= n)
		{
			int cur_width = image.width();
			int cur_height = image.height();
			int cur_depth = image.depth();
			int nChannels = image.nchannels();

			float sigma = baseSigma*i;
			foo.allocate(cur_width,cur_height,cur_depth,nChannels);
			cuda_cost_time += 
				ZQ_CUDA_ImageProcessing3D::GaussianSmoothing2_3D(foo.data(),image.data(),sigma,sigma*3,cur_width,cur_height,cur_depth,nChannels);

			double cur_ratio = pow(this->ratio,i);

			int dst_width = cur_width*cur_ratio;
			int dst_height = cur_height*cur_ratio;
			int dst_depth = cur_depth*cur_ratio;

			ImPyramid[i].allocate(dst_width,dst_height,dst_depth,nChannels);

			cuda_cost_time += 
				ZQ_CUDA_ImageProcessing3D::ResizeImage3D_Tricubic(ImPyramid[i].data(),foo.data(),cur_width,cur_height,cur_depth,dst_width,dst_height,dst_depth,nChannels);
		}
		else
		{
			int cur_width = ImPyramid[i-n].width();
			int cur_height = ImPyramid[i-n].height();
			int cur_depth = ImPyramid[i-n].depth();
			int nChannels = ImPyramid[i-n].nchannels();

			foo.allocate(cur_width,cur_height,cur_depth,nChannels);

			cuda_cost_time +=
				ZQ_CUDA_ImageProcessing3D::GaussianSmoothing2_3D(foo.data(),ImPyramid[i-n].data(),nSigma,nSigma*3,cur_width,cur_height,cur_depth,nChannels);

			int dst_width = pow(this->ratio,i)*image.width();
			int dst_height = pow(this->ratio,i)*image.height();
			int dst_depth = pow(this->ratio,i)*image.depth();

			ImPyramid[i].allocate(dst_width,dst_height,dst_depth,nChannels);

			cuda_cost_time += 
				ZQ_CUDA_ImageProcessing3D::ResizeImage3D_Tricubic(ImPyramid[i].data(),foo.data(),cur_width,cur_height,cur_depth,dst_width,dst_height,dst_depth,nChannels);
		}
	}

	return this->ratio;
}
