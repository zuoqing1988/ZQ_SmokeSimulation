#ifndef _ZQ_OPTICAL_FLOW_2D_CUDA_H_
#define _ZQ_OPTICAL_FLOW_2D_CUDA_H_

#pragma once

#include "ZQ_GaussianPyramid2DCuda.h"
#include "ZQ_CUDA_OpticalFlow2D_Utils.h"
#include "ZQ_OpticalFlowOptions.h"
#include "ZQ_OpticalFlow.h"

namespace ZQ
{
	class ZQ_OpticalFlow2DCuda : ZQ_OpticalFlow
	{
	public:
		ZQ_OpticalFlow2DCuda();
		~ZQ_OpticalFlow2DCuda();

	public:

		static float Coarse2Fine_HS_L2(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_HS_L2_Occupy(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_DImage<float>& occupy, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_HS_L1(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_HS_DL1(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_ADMM_L2(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_ADMM_L2_Occupy(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_DImage<float>& occupy, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_ADMM_DL1(ZQ_DImage<float>& u, ZQ_DImage<float>& v, ZQ_DImage<float>& warpIm2, const ZQ_DImage<float>& Im1, const ZQ_DImage<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Inc_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Inc_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Dec_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Dec_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Inc_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Inc_L2_Occupy(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_DImage<float>& occupy, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Inc_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Dec_L2(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Dec_L2_Occupy(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_DImage<float>& occupy, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Dec_DL1(std::vector<ZQ_DImage<float>>& u, std::vector<ZQ_DImage<float>>& v, std::vector<ZQ_DImage<float>>& warpIm, const std::vector<ZQ_DImage<float>>& Im, const ZQ_OpticalFlowOptions& opt);
	};
}


#endif