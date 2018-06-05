#ifndef _ZQ_CUDA_ATTENUATION_3D_CUH_
#define _ZQ_CUDA_ATTENUATION_3D_CUH_

namespace ZQ_CUDA_Attenuation3D
{
	__global__
	void Atten_u_Kernel(float* mac_u, const bool* occupy, const float velAtten, const int width, const int height, const int depth);
	
	__global__
	void Atten_v_Kernel(float* mac_v, const bool* occupy, const float velAtten, const int width, const int height, const int depth);
	
	__global__
	void Atten_w_Kernel(float* mac_w, const bool* occupy, const float velAtten, const int width, const int height, const int depth);

	__global__
	void Atten_temperature_density_Kernel(float* temperature, float* density, const float tempAtten, const float densityAtten, const int width, const int height, const int depth);

	/********************************************/

	void cu_Attenuation3D(float* mac_u, float* mac_v, float* mac_w, float* temperature, float* density, const bool* occupy, 
		const float velAtten, const float tempAtten, const float densityAtten, const int width, const int height, const int depth);
}

#endif