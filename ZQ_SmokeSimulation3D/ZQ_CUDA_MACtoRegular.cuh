#ifndef _ZQ_CUDA_MAC_TO_REGULAR_CUH_
#define _ZQ_CUDA_MAC_TO_REGULAR_CUH_
#pragma once

namespace ZQ_CUDA_MACtoRegular
{
	__global__
	void MAC_to_Regular4_Kernel(const float* mac_u, const float* mac_v, const float* mac_w, float4* vel4, int width, int height, int depth);

	__global__
	void Regular4_to_MAC_u_Kernel(const float4* vel4, float* mac_u, int width, int height, int depth);

	__global__
	void Regular4_to_MAC_v_Kernel(const float4* vel4, float* mac_v, int width, int height, int depth);

	__global__
	void Regular4_to_MAC_w_Kernel(const float4* vel4, float* mac_w, int width, int height, int depth);

	/************************/

	void cu_MAC_to_Regular4(const float* mac_u, const float* mac_v, const float* mac_w, float4* vel4, int width, int height, int depth);
	
	void cu_Regular4_to_MAC(const float4* vel4, float* mac_u, float* mac_v, float* mac_w, int width, int height, int depth);
}

#endif