#ifndef _ZQ_CUDA_ADVECTION_3D_CUH_
#define _ZQ_CUDA_ADVECTION_3D_CUH_

#include "ZQlibCudaDefines.cuh"

namespace ZQ_CUDA_Advection3D
{
	/**************************** Kernel functions Begin **************************/
	__global__ void Velocity_Negative_Kernel(float* in_out_data, int len);
	__global__ void Velocity_Negative_4channels_Kernel(float4* in_out_data, int len);
	__global__ void Input_Increment_Kernel(float* input, const float* input_star, int len);
	__global__ void Input_Increment_4channels_Kernel(float4* input, const float4* input_star, int len);
	__global__ void Advect_Velocity_inRegular_outRegular_Kernel(float4 * d_output);
	__global__ void Advect_Velocity_inRegular_outMAC_u_Kernel(float* d_mac_u);
	__global__ void Advect_Velocity_inRegular_outMAC_v_Kernel(float* d_mac_v);
	__global__ void Advect_Velocity_inRegular_outMAC_w_Kernel(float* d_mac_w);
	__global__ void Advect_Velocity_inMAC_outMAC_u_Kernel(float* d_mac_u);
	__global__ void Advect_Velocity_inMAC_outMAC_v_Kernel(float* d_mac_v);
	__global__ void Advect_Velocity_inMAC_outMAC_w_Kernel(float* d_mac_w);
	__global__ void Advect_Scalar_Regular_Velocity_Kernel(float* d_temperature);
	__global__ void Advect_Scalar_MAC_Velocity_Kernel(float2* d_density);
	__global__ void Apply_Advect_Velocity_Result_Open_u_Kernel(const float* adv_u, const bool* occupy, float* u, int width, int height, int depth);
	__global__ void Apply_Advect_Velocity_Result_Open_v_Kernel(const float* adv_v, const bool* occupy, float* v, int width, int height, int depth);
	__global__ void Apply_Advect_Velocity_Result_Open_w_Kernel(const float* adv_w, const bool* occupy, float* w, int width, int height, int depth);
	__global__ void Apply_Advect_Velocity_Result_Closed_u_Kernel(const float* adv_u, const bool* occupy, float* u, int width, int height, int depth);
	__global__ void Apply_Advect_Velocity_Result_Closed_v_Kernel(const float* adv_v, const bool* occupy, float* v, int width, int height, int depth);
	__global__ void Apply_Advect_Velocity_Result_Closed_w_Kernel(const float* adv_w, const bool* occupy, float* w, int width, int height, int depth);
	/**************************** Kernel functions End **************************/


	/*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*/
	void cu_Velocity_Negative(float* u, float* v, float* w, int width, int height, int depth);

	/*	vel		: [width*height*depth]	:	IN & OUT
	*/
	void cu_Velocity_Negative_4channels(float4* vel, int width, int height, int depth);

	void cu_Input_Increment(float* input, const float* input_star, int len);

	void cu_Input_Increment_4channels(float4* input, const float4* input_star, int len);

	/**************************/
	/*	adv_u	:	[(width+1)*height*depth]	:	IN
	*	adv_v	:	[width*(height+1)*depth]	:	IN
	*	adv_w	:	[width*height*(depth+1)]	:	IN
	*	occupy	:	[width*height*depth]		:	IN
	*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*/
	void cu_Apply_Advect_Velocity_Result_Open(const float* adv_u, const float* adv_v, const float* adv_w, const bool* occupy, float* u, float* v, float* w);
	
	/**************************/
	/*	adv_u	:	[(width+1)*height*depth]	:	IN
	*	adv_v	:	[width*(height+1)*depth]	:	IN
	*	adv_w	:	[width*height*(depth+1)]	:	IN
	*	occupy	:	[width*height*depth]		:	IN
	*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*/
	void cu_Apply_Advect_Velocity_Result_Closed(const float* adv_u, const float* adv_v, const float* adv_w, const bool* occupy, float* u, float* v, float* w);

	
	/*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*	occupy	:	[width*height*depth]		:	IN
	*/
	void cu_Advect_Velocity_inRegular_outRegular(float* u, float* v, float* w, const bool* occupy, bool is_open);

	/*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*	occupy	:	[width*height*depth]		:	IN
	*/
	void cu_Advect_Velocity_inRegular_outRegular_BFECC(float* u, float* v, float* w, const bool* occupy, bool is_open);

	/*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*	occupy	:	[width*height*depth]		:	IN
	*/
	void cu_Advect_Velocity_inMAC_outMAC(float* u, float* v, float* w, const bool* occupy, bool is_open);

	/*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*	occupy	:	[width*height*depth]		:	IN
	*/
	void cu_Advect_Velocity_inMAC_outMAC_BFECC(float* u, float* v, float* w, const bool* occupy, bool is_open);

	/*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*	occupy	:	[width*height*depth]		:	IN
	*/
	void cu_Advect_Velocity_inRegular_outMAC(float* u, float* v, float* w, const bool* occupy, bool is_open);

	/*	u		:	[(width+1)*height*depth]	:	IN & OUT
	*	v		:	[width*(height+1)*depth]	:	IN & OUT
	*	w		:	[width*height*(depth+1)]	:	IN & OUT	
	*	occupy	:	[width*height*depth]		:	IN
	*/
	void cu_Advect_Velocity_inRegular_outMAC_BFECC(float* u, float* v, float* w, const bool* occupy, bool is_open);


	/*	u		:	[(width+1)*height*depth]	:	IN
	*	v		:	[width*(height+1)*depth]	:	IN
	*	w		:	[width*height*(depth+1)]	:	IN
	*	occupy	:	[width*height*depth]		:	IN
	*	input_temperature	:	[width*height*depth]		:	IN
	*	input_density	:	[width*height*depth]		:	IN
	*	output_temperature	:	[width*height*depth]	:	OUT
	*	output_density	:	[width*height*depth]	:	OUT
	*/
	void cu_Advect_Scalar_Regular_Velocity(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density, 
								float* output_temperature, float* output_density);

	/*	u		:	[(width+1)*height*depth]	:	IN
	*	v		:	[width*(height+1)*depth]	:	IN
	*	w		:	[width*height*(depth+1)]	:	IN
	*	occupy	:	[width*height*depth]		:	IN
	*	input_temperature	:	[width*height*depth]		:	IN
	*	input_density	:	[width*height*depth]		:	IN
	*	output_temperature	:	[width*height*depth]	:	OUT
	*	output_density	:	[width*height*depth]	:	OUT
	*/
	void cu_Advect_Scalar_Regular_Velocity_BFECC(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density, 
								float* output_temperature, float* output_density);

	/*	u		:	[(width+1)*height*depth]	:	IN
	*	v		:	[width*(height+1)*depth]	:	IN
	*	w		:	[width*height*(depth+1)]	:	IN
	*	occupy	:	[width*height*depth]		:	IN
	*	input_temperature	:	[width*height*depth]		:	IN
	*	input_density	:	[width*height*depth]		:	IN
	*	output_temperature	:	[width*height*depth]	:	OUT
	*	output_density	:	[width*height*depth]	:	OUT
	*/
	void cu_Advect_Scalar_MAC_Velocity(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density, 
								float* output_temperature, float* output_density);
								
	/*	u		:	[(width+1)*height*depth]	:	IN
	*	v		:	[width*(height+1)*depth]	:	IN
	*	w		:	[width*height*(depth+1)]	:	IN
	*	occupy	:	[width*height*depth]		:	IN
	*	input_temperature	:	[width*height*depth]		:	IN
	*	input_density	:	[width*height*depth]		:	IN
	*	output_temperature	:	[width*height*depth]	:	OUT
	*	output_density	:	[width*height*depth]	:	OUT
	*/
	void cu_Advect_Scalar_MAC_Velocity_BFECC(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density, 
								float* output_temperature, float* output_density);		
}


#endif