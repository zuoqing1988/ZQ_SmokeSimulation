#ifndef _ZQ_CUDA_ATTENUATION_3D_H_
#define _ZQ_CUDA_ATTENUATION_3D_H_

namespace ZQ_CUDA_Attenuation3D
{
	extern "C"
	float Attenuation3D(float* mac_u, float* mac_v, float* mac_w, float* temperature, float* density, const bool* occupy,
		const float velAtten, const float tempAtten, const float densityAtten, const int width, const int height, const int depth);
}

#endif