#ifndef _ZQ_OPTICAL_FLOW_3D_CUDA_H_
#define _ZQ_OPTICAL_FLOW_3D_CUDA_H_

#pragma once

#include "ZQ_GaussianPyramid3DCuda.h"
#include "ZQ_CUDA_OpticalFlow3D_Utils.h"
#include "ZQ_OpticalFlowOptions.h"
#include "ZQ_OpticalFlow3D.h"
#include "ZQ_OpticalFlow3DBuffer.h"

namespace ZQ
{
	class ZQ_OpticalFlow3DCuda : ZQ_OpticalFlow3D
	{
	public:
		ZQ_OpticalFlow3DCuda();
		~ZQ_OpticalFlow3DCuda();

	public:

		static float Coarse2Fine_HS_L2(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
												const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_HS_L1(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
												const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_HS_DL1(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
												const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_ADMM_L2(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
												const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_ADMM_DL1(ZQ_DImage3D<float>& u, ZQ_DImage3D<float>& v, ZQ_DImage3D<float>& w, ZQ_DImage3D<float>& warpIm2, 
												const ZQ_DImage3D<float>& Im1, const ZQ_DImage3D<float>& Im2, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Inc_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Inc_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w,  
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Dec_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w, 
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Dec_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w,  
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Inc_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w,  
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Inc_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w,  
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Dec_L2(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w,  
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Dec_DL1(std::vector<ZQ_DImage3D<float>>& u, std::vector<ZQ_DImage3D<float>>& v, std::vector<ZQ_DImage3D<float>>& w,  
												std::vector<ZQ_DImage3D<float>>& warpIm, const std::vector<ZQ_DImage3D<float>>& Im, const ZQ_OpticalFlowOptions& opt);

		// Out Of Core
		static float Coarse2Fine_OneDir_Inc_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Inc_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Dec_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_OneDir_Dec_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Inc_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Inc_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Dec_L2(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);

		static float Coarse2Fine_TwoDir_Dec_DL1(ZQ_OpticalFlow3DBuffer* buffer, const ZQ_OpticalFlowOptions& opt);
	};
}


#endif