#ifndef _ZQ_CUDA_ADD_FORCE_3D_H_
#define _ZQ_CUDA_ADD_FORCE_3D_H_

namespace ZQ_CUDA_AddForce3D
{
	enum AddForceType
	{
		ADD_FORCE_ENTIRE,
		ADD_FORCE_SLICE
	};

	extern "C"
	float AddForce3D(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* temperature, const float deltah, const float deltat,
		const float buoyCoeff, const float confineCoeff, const float Tamb, const int width, const int height, const int depth, enum AddForceType type);

}

#endif