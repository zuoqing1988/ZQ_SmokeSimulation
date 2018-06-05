#ifndef _ZQ_CUDA_ADVECTION_3D_H_
#define _ZQ_CUDA_ADVECTION_3D_H_
#pragma once

namespace ZQ_CUDA_Advection3D
{
	enum AdvectVelocityType{
		ADV_VEL_INREG_OUTREG = 0,
		ADV_VEL_INREG_OUTREG_BFECC,
		ADV_VEL_INMAC_OUTMAC,
		ADV_VEL_INMAC_OUTMAC_BFECC,
		ADV_VEL_INREG_OUTMAC,
		ADV_VEL_INREG_OUTMAC_BFECC
	};

	enum AdvectScalarType{
		ADV_SCA_MAC,
		ADV_SCA_MAC_BFECC,
		ADV_SCA_REG,
		ADV_SCA_REG_BFECC
	};

	extern "C" void ZQ_Cuda_Prepare_Advection(const unsigned int width, const unsigned int height, const unsigned int depth, const float voxelSize, const unsigned int steps, const float deltatt);


	/*	u		: [(width+1)*height*depth]
	*	v		: [width*(height+1)*depth]
	*	w		: [width*height*(depth+1)]
	*	occupy	: [width*height*depth]
	*/
	extern "C" float Advect_Velocity(float* u, float* v, float* w, const bool* occupy, bool is_open, enum AdvectVelocityType type);


	/*	u		:	[(width+1)*height*depth]	:	IN
	*	v		:	[width*(height+1)*depth]	:	IN
	*	w		:	[width*height*(depth+1)]	:	IN
	*	occupy	:	[width*height*depth]		:	IN
	*	input_temperature	:	[width*height*depth]		:	IN
	*	input_density	:	[width*height*depth]		:	IN
	*	output_temperature	:	[width*height*depth]	:	OUT
	*	output_density	:	[width*height*depth]	:	OUT
	*/
	extern "C" float Advect_Scalar(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density,
						float* output_temperature, float* output_density, enum AdvectScalarType type);

}

#endif