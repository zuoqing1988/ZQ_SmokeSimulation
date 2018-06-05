#ifndef _ZQ_GAUSSIAN_PYRAMID_2D_CUDA_H_
#define _ZQ_GAUSSIAN_PYRAMID_2D_CUDA_H_

#pragma once


#include "ZQ_DoubleImage.h"
#include "ZQ_CUDA_ImageProcessing2D.h"
#include <vector>
#include <iostream>

namespace ZQ
{
	class ZQ_GaussianPyramid2DCuda
	{
	private:
		std::vector<ZQ::ZQ_DImage<float>> ImPyramid;
		int nLevels;
		float ratio;
	public:
		ZQ_GaussianPyramid2DCuda(void);
		~ZQ_GaussianPyramid2DCuda(void);
		float ConstructPyramid(const ZQ::ZQ_DImage<float>& image, float& cuda_cost_time, const float ratio = 0.5, const int minWidth = 16);
		int nlevels() const {return nLevels;};
		ZQ::ZQ_DImage<float>& Image(int index) { return ImPyramid[index]; };
	};
}

#endif