#ifndef _ZQ_GAUSSIAN_PYRAMID_3D_CUDA_H_
#define _ZQ_GAUSSIAN_PYRAMID_3D_CUDA_H_

#pragma once


#include "ZQ_DoubleImage3D.h"
#include "ZQ_CUDA_ImageProcessing3D.h"
#include <vector>
#include <iostream>

namespace ZQ
{
	class ZQ_GaussianPyramid3DCuda
	{
	private:
		std::vector<ZQ::ZQ_DImage3D<float>> ImPyramid;
		int nLevels;
		float ratio;
	public:
		ZQ_GaussianPyramid3DCuda(void);
		~ZQ_GaussianPyramid3DCuda(void);
		float ConstructPyramid(const ZQ::ZQ_DImage3D<float>& image, float& cuda_cost_time, const float ratio = 0.5, const int minWidth = 16);
		int nlevels() const {return nLevels;};
		ZQ::ZQ_DImage3D<float>& Image(int index) { return ImPyramid[index]; };
	};
}

#endif