#ifndef _ZQ_CUDA_ADVECTION_3D_CU_
#define _ZQ_CUDA_ADVECTION_3D_CU_

#include "ZQ_CUDA_Advection3D.h"
#include "ZQ_CUDA_Advection3D.cuh"
#include "ZQ_CUDA_MACtoRegular.cuh"

namespace ZQ_CUDA_Advection3D
{
	texture<float4,3,cudaReadModeElementType> tex_velocity_regular;
	texture<float,3,cudaReadModeElementType> tex_velocity_MAC_u;
	texture<float,3,cudaReadModeElementType> tex_velocity_MAC_v;
	texture<float,3,cudaReadModeElementType> tex_velocity_MAC_w;
	texture<unsigned char,3,cudaReadModeElementType> tex_occupy;
	texture<float,3,cudaReadModeElementType> tex_inputVelocity_MAC_u;
	texture<float,3,cudaReadModeElementType> tex_inputVelocity_MAC_v;
	texture<float,3,cudaReadModeElementType> tex_inputVelocity_MAC_w;
	texture<float4,3,cudaReadModeElementType> tex_inputVelocity_regular;
	texture<float,3,cudaReadModeElementType> tex_temperature; 
	texture<float,3,cudaReadModeElementType> tex_density; 

	unsigned h_width;
	unsigned h_height;
	unsigned h_depth;
	unsigned int h_steps;
	float h_voxelSize;
	float h_deltatt;
	
	__constant__ unsigned int d_width;
	__constant__ unsigned int d_height;
	__constant__ unsigned int d_depth;
	__constant__ unsigned int d_steps;
	__constant__ float d_voxelSize;
	__constant__ float d_deltatt;
	
	/****************************************************************************************/
	
	__global__ 
	void Velocity_Negative_Kernel(float* in_out_data, int len)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		if(x >= len)
			return;
		in_out_data[x] = -in_out_data[x];
	}

	__global__ 
	void Velocity_Negative_4channels_Kernel(float4* in_out_data, int len)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		if(x >= len)
			return;
		in_out_data[x] = -in_out_data[x];
	}

	__global__ 
	void Input_Increment_Kernel(float* input, const float* input_star, int len)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		if(x >= len)
			return;

		input[x] = input[x]*1.5 - input_star[x]*0.5;
	}

	__global__ 
	void Input_Increment_4channels_Kernel(float4* input, const float4* input_star, int len)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		if(x >= len)
			return;

		input[x] = input[x]*1.5 - input_star[x]*0.5;
	}

	__global__ 
	void Advect_Velocity_inRegular_outRegular_Kernel(float4 * d_output)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x >= d_width || y >= d_height)
			return;
			
		for(int z = 0;z < d_depth;z++)
		{
			float4 pos = make_float4(x+0.5f,y+0.5f,z+0.5f,0.0f);
			float4 lastpos = pos;
			float3 velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
			float4 lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = velCoord;
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height&& pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);

				lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord = make_float3(lastpos.x/d_width,lastpos.y/d_height,lastpos.z/d_depth);
			float4 tempvel = tex3D(tex_inputVelocity_regular,out_coord.x,out_coord.y,out_coord.z);
			d_output[z*d_height*d_width+y*d_width+x] = tempvel;
		}
	}

	__global__ 
	void Advect_Velocity_inRegular_outMAC_u_Kernel(float * d_mac_u)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x > d_width || y >= d_height)
			return;

		for(int z = 0;z < d_depth;z++)
		{
			float4 pos = make_float4((float)x,y+0.5f,z+0.5f,0.0f);
			float4 lastpos = pos;
			float3 velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
			float4 lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = velCoord;
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);

				lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord = make_float3(lastpos.x/d_width,lastpos.y/d_height,lastpos.z/d_depth);
			float4 tempvel = tex3D(tex_inputVelocity_regular,out_coord.x,out_coord.y,out_coord.z);
			d_mac_u[z*d_height*(d_width+1)+y*(d_width+1)+x] = tempvel.x;
		}
	}

	__global__
	void Advect_Velocity_inRegular_outMAC_v_Kernel(float * d_mac_v)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x >= d_width || y > d_height)
			return;

		for(int z = 0;z < d_depth;z++)
		{
			float4 pos = make_float4(x+0.5f,(float)y,z+0.5f,0.0f);
			float4 lastpos = pos;
			float3 velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
			float4 lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = velCoord;
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);

				lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord = make_float3(lastpos.x/d_width,lastpos.y/d_height,lastpos.z/d_depth);
			float4 tempvel = tex3D(tex_inputVelocity_regular,out_coord.x,out_coord.y,out_coord.z);
			d_mac_v[z*(d_height+1)*d_width+y*d_width+x] = tempvel.y;
		}
	}
	
	__global__
	void Advect_Velocity_inRegular_outMAC_w_Kernel(float * d_mac_w)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x >= d_width || y >= d_height)
			return;

		for(int z = 0;z <= d_depth;z++)
		{
			float4 pos = make_float4(x+0.5f,y+0.5f,(float)z,0.0f);
			float4 lastpos = pos;
			float3 velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
			float4 lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = velCoord;
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);

				lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord = make_float3(lastpos.x/d_width,lastpos.y/d_height,lastpos.z/d_depth);
			float4 tempvel = tex3D(tex_inputVelocity_regular,out_coord.x,out_coord.y,out_coord.z);
			d_mac_w[z*d_height*d_width+y*d_width+x] = tempvel.z;
		}
	}

	__global__
	void Advect_Velocity_inMAC_outMAC_u_Kernel(float * d_mac_u)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x > d_width || y >= d_height)
			return;

		for(int z = 0;z < d_depth;z++)
		{
			float4 pos = make_float4((float)x,y+0.5f,z+0.5f,0.0f);
			float4 lastpos = pos;
			float3 velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
			float3 velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
			float3 velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));
			float4 lastvel = make_float4(
				tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
				tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
				tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
				velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
				velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));

				lastvel = make_float4(
					tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
					tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
					tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord_u = make_float3((lastpos.x+0.5f)/(d_width+1),lastpos.y/d_height,lastpos.z/d_depth);
			float tempvel = tex3D(tex_inputVelocity_MAC_u,out_coord_u.x,out_coord_u.y,out_coord_u.z);
			d_mac_u[z*d_height*(d_width+1)+y*(d_width+1)+x] = tempvel;
		}
	}

	__global__
	void Advect_Velocity_inMAC_outMAC_v_Kernel(float * d_mac_v)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x >= d_width || y > d_height)
			return;

		for(int z = 0;z < d_depth;z++)
		{
			float4 pos = make_float4(x+0.5f,(float)y,z+0.5f,0.0f);
			float4 lastpos = pos;
			float3 velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
			float3 velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
			float3 velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));
			float4 lastvel = make_float4(
				tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
				tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
				tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
				velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
				velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));

				lastvel = make_float4(
					tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
					tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
					tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord_v = make_float3(lastpos.x/d_width,(lastpos.y+0.5f)/(d_height+1),lastpos.z/d_depth);
			float tempvel = tex3D(tex_inputVelocity_MAC_v,out_coord_v.x,out_coord_v.y,out_coord_v.z);
			d_mac_v[z*(d_height+1)*d_width+y*d_width+x] = tempvel;
		}
	}
	
	__global__
	void Advect_Velocity_inMAC_outMAC_w_Kernel(float * d_mac_w)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x >= d_width || y >= d_height)
			return;

		for(int z = 0;z <= d_depth;z++)
		{
			float4 pos = make_float4(x+0.5f,y+0.5f,(float)z,0.0f);
			float4 lastpos = pos;
			float3 velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
			float3 velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
			float3 velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));
			float4 lastvel = make_float4(
				tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
				tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
				tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
				velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
				velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));

				lastvel = make_float4(
					tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
					tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
					tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord_w = make_float3(lastpos.x/d_width,lastpos.y/d_height,(lastpos.z+0.5f)/(d_depth+1));
			float tempvel = tex3D(tex_inputVelocity_MAC_w,out_coord_w.x,out_coord_w.y,out_coord_w.z);
			d_mac_w[z*d_height*d_width+y*d_width+x] = tempvel;
		}
	}
		
	__global__ 
	void Advect_Scalar_Regular_Velocity_Kernel(float* d_output_temperature, float* d_output_density)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x >= d_width || y >= d_height)
			return;

		for(int z = 0;z < d_depth;z++)
		{
			float4 pos = make_float4(x+0.5f,y+0.5f,z+0.5f,0.0f);
			float4 lastpos = pos;
			float3 velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
			float4 lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = velCoord;
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);

				lastvel = tex3D(tex_velocity_regular,velCoord.x,velCoord.y,velCoord.z);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord = make_float3(lastpos.x/d_width,lastpos.y/d_height,lastpos.z/d_depth);
			float temperature = tex3D(tex_temperature,out_coord.x,out_coord.y,out_coord.z);
			float density = tex3D(tex_density,out_coord.x,out_coord.y,out_coord.z);
			d_output_temperature[z*d_height*d_width+y*d_width+x] = temperature;
			d_output_density[z*d_height*d_width+y*d_width+x] = density;
		}
	}

	__global__ 
	void Advect_Scalar_MAC_Velocity_Kernel(float* d_output_temperature, float* d_output_density)
	{
		int x = threadIdx.x + blockIdx.x * blockDim.x;
		int y = threadIdx.y + blockIdx.y * blockDim.y;

		if(x >= d_width || y >= d_height)
			return;

		for(int z = 0;z < d_depth;z++)
		{
			float4 pos = make_float4(x+0.5f,y+0.5f,z+0.5f,0.0f);
			float4 lastpos = pos;
			float3 velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
			float3 velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
			float3 velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));
			float4 lastvel = make_float4(
				tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
				tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
				tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);

			unsigned int istep = 0;
			do 
			{
				float3 occupyCoord = make_float3(pos.x/d_width,pos.y/d_height,pos.z/d_depth);
				if(!(pos.x >= 0 && pos.x <= d_width && pos.y >= 0 && pos.y <= d_height && pos.z >= 0 && pos.z <= d_depth))
					break;
				if(tex3D(tex_occupy,occupyCoord.x,occupyCoord.y,occupyCoord.z) != 0)
					break;

				lastpos = pos;
				pos -= lastvel * d_deltatt / d_voxelSize;
				velCoord_u = make_float3((pos.x+0.5f)/(d_width+1),pos.y/d_height,pos.z/d_depth);
				velCoord_v = make_float3(pos.x/d_width,(pos.y+0.5f)/(d_height+1),pos.z/d_depth);
				velCoord_w = make_float3(pos.x/d_width,pos.y/d_height,(pos.z+0.5f)/(d_depth+1));

				lastvel = make_float4(
					tex3D(tex_velocity_MAC_u,velCoord_u.x,velCoord_u.y,velCoord_u.z),
					tex3D(tex_velocity_MAC_v,velCoord_v.x,velCoord_v.y,velCoord_v.z),
					tex3D(tex_velocity_MAC_w,velCoord_w.x,velCoord_w.y,velCoord_w.z),0.0f);
				istep ++;
			} while (istep < d_steps);

			float3 out_coord = make_float3(lastpos.x/d_width,lastpos.y/d_height,lastpos.z/d_depth);
			float temperature = tex3D(tex_temperature,out_coord.x,out_coord.y,out_coord.z);
			float density = tex3D(tex_density,out_coord.x,out_coord.y,out_coord.z);
			d_output_temperature[z*d_height*d_width+y*d_width+x] = temperature;
			d_output_density[z*d_height*d_width+y*d_width+x] = density;
		}
	}

	__global__ 
	void Apply_Advect_Velocity_Result_Open_u_Kernel(const float* adv_u, const bool* occupy, float* u, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x > width || y >= height)
			return;

		int offset = y*width+x;
		int u_offset = y*(width+1)+x;
		for (int z = 0; z < depth; z++)
		{
			if(x == 0)
			{
				if(!occupy[offset])
					u[u_offset] = adv_u[u_offset];
			}
			else if(x == width)
			{
				if(!occupy[offset-1])
					u[u_offset] = adv_u[u_offset];
			}
			else
			{
				if(!occupy[offset] && !occupy[offset-1])
					u[u_offset] = adv_u[u_offset];
			}
			offset += height*width;
			u_offset += height*(width+1);
		}
	}

	__global__ 
	void Apply_Advect_Velocity_Result_Open_v_Kernel(const float* adv_v, const bool* occupy, float* v, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x >= width || y > height)
			return;

		int offset = y*width+x;
		int v_offset = y*width+x;
		for (int z = 0; z < depth; z++)
		{
			if(y == 0)
			{
				if(!occupy[offset])
					v[v_offset] = adv_v[v_offset];
			}
			else if(y == height)
			{
				if(!occupy[offset-width])
					v[v_offset] = adv_v[v_offset];
			}
			else
			{
				if(!occupy[offset] && !occupy[offset-width])
					v[v_offset] = adv_v[v_offset];
			}
			offset += height*width;
			v_offset += (height+1)*width;
		}
	}

	__global__ 
	void Apply_Advect_Velocity_Result_Open_w_Kernel(const float* adv_w, const bool* occupy, float* w, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x >= width || y >= height)
			return;

		int offset = y*width+x;
		int w_offset = y*width+x;
		if(!occupy[offset])
			w[w_offset] = adv_w[w_offset];
		offset += height*width;
		w_offset += height*width;

		for (int z = 1; z < depth; z++)
		{
			if(!occupy[offset] && !occupy[offset-height*width])
				w[w_offset] = adv_w[w_offset];
			offset += height*width;
			w_offset += height*width;
		}

		if(!occupy[offset-height*width])
			w[w_offset] = adv_w[w_offset];
	}

	__global__ 
	void Apply_Advect_Velocity_Result_Closed_u_Kernel(const float* adv_u, const bool* occupy, float* u, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x == 0 || x >= width || y >= height)
			return;

		int offset = y*width+x;
		int u_offset = y*(width+1)+x;
		for (int z = 0; z < depth; z++)
		{
			if(!occupy[offset] && !occupy[offset-1])
				u[u_offset] = adv_u[u_offset];
			offset += height*width;
			u_offset += height*(width+1);
		}
	}

	__global__ 
	void Apply_Advect_Velocity_Result_Closed_v_Kernel(const float* adv_v, const bool* occupy, float* v, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x >= width || y == 0 || y >= height)
			return;

		int offset = y*width+x;
		int v_offset = y*width+x;
		for (int z = 0; z < depth; z++)
		{
			if(!occupy[offset] && !occupy[offset-width])
				v[v_offset] = adv_v[v_offset];
			offset += height*width;
			v_offset += (height+1)*width;
		}
	}

	__global__ 
	void Apply_Advect_Velocity_Result_Closed_w_Kernel(const float* adv_w, const bool* occupy, float* w, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x >= width || y >= height)
			return;

		int offset = y*width+x;
		int w_offset = y*width+x;
		offset += height*width;
		w_offset += height*width;

		for (int z = 1; z < depth; z++)
		{
			if(!occupy[offset] && !occupy[offset-height*width])
				w[w_offset] = adv_w[w_offset];
			offset += height*width;
			w_offset += height*width;
		}
	}

	/****************************************************************************************/

	void cu_Copy_to_tex_velocity_regular(const float4* vel, cudaArray** velocity_array)
	{
		tex_velocity_regular.normalized = true;
		tex_velocity_regular.filterMode = cudaFilterModeLinear;
		tex_velocity_regular.addressMode[0] = cudaAddressModeClamp;
		tex_velocity_regular.addressMode[1] = cudaAddressModeClamp;
		tex_velocity_regular.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf4 = cudaCreateChannelDesc<float4>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height,h_depth);

		checkCudaErrors( cudaMalloc3DArray(velocity_array, &channelDescf4, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)vel, texSize.width*sizeof(float4), texSize.width, texSize.height);
		copyParams.dstArray = *velocity_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_velocity_regular,*velocity_array,channelDescf4) );
	}

	void cu_Free_tex_velocity_regular(cudaArray** velocity_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_velocity_regular) );
		checkCudaErrors( cudaFreeArray(*velocity_array) );
		*velocity_array = 0;
	}

	void cu_Copy_to_tex_velocity_MAC_u(const float* u, cudaArray** u_array)
	{
		tex_velocity_MAC_u.normalized = true;
		tex_velocity_MAC_u.filterMode = cudaFilterModeLinear;
		tex_velocity_MAC_u.addressMode[0] = cudaAddressModeClamp;
		tex_velocity_MAC_u.addressMode[1] = cudaAddressModeClamp;
		tex_velocity_MAC_u.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width+1,h_height,h_depth);

		checkCudaErrors( cudaMalloc3DArray(u_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)u, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *u_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_velocity_MAC_u,*u_array,channelDescf) );
	}

	void cu_Free_tex_velocity_MAC_u(cudaArray** u_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_velocity_MAC_u) );
		checkCudaErrors( cudaFreeArray(*u_array) );
		*u_array = 0;
	}

	void cu_Copy_to_tex_velocity_MAC_v(const float* v, cudaArray** v_array)
	{
		tex_velocity_MAC_v.normalized = true;
		tex_velocity_MAC_v.filterMode = cudaFilterModeLinear;
		tex_velocity_MAC_v.addressMode[0] = cudaAddressModeClamp;
		tex_velocity_MAC_v.addressMode[1] = cudaAddressModeClamp;
		tex_velocity_MAC_v.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height+1,h_depth);

		checkCudaErrors( cudaMalloc3DArray(v_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)v, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *v_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_velocity_MAC_v,*v_array,channelDescf) );
	}

	void cu_Free_tex_velocity_MAC_v(cudaArray** v_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_velocity_MAC_v) );
		checkCudaErrors( cudaFreeArray(*v_array) );
		*v_array = 0;
	}

	void cu_Copy_to_tex_velocity_MAC_w(const float* w, cudaArray** w_array)
	{
		tex_velocity_MAC_w.normalized = true;
		tex_velocity_MAC_w.filterMode = cudaFilterModeLinear;
		tex_velocity_MAC_w.addressMode[0] = cudaAddressModeClamp;
		tex_velocity_MAC_w.addressMode[1] = cudaAddressModeClamp;
		tex_velocity_MAC_w.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height,h_depth+1);

		checkCudaErrors( cudaMalloc3DArray(w_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)w, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *w_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_velocity_MAC_w,*w_array,channelDescf) );
	}

	void cu_Free_tex_velocity_MAC_w(cudaArray** w_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_velocity_MAC_w) );
		checkCudaErrors( cudaFreeArray(*w_array) );
		*w_array = 0;
	}

	void cu_Copy_to_tex_occupy(const bool* occupy, cudaArray** occupy_array)
	{
		tex_occupy.normalized = true;                      
		tex_occupy.filterMode = cudaFilterModePoint;     
		tex_occupy.addressMode[0] = cudaAddressModeClamp; 
		tex_occupy.addressMode[1] = cudaAddressModeClamp;
		tex_occupy.addressMode[2] = cudaAddressModeClamp;	

		cudaExtent texSize = make_cudaExtent(h_width,h_height,h_depth);
		cudaChannelFormatDesc channelDescb = cudaCreateChannelDesc<uchar1>();

		checkCudaErrors( cudaMalloc3DArray(occupy_array, &channelDescb, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)occupy, texSize.width*sizeof(uchar1), texSize.width, texSize.height);
		copyParams.dstArray = *occupy_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_occupy,*occupy_array,channelDescb) );
	}

	void cu_Free_tex_occupy(cudaArray** occupy_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_occupy) );
		checkCudaErrors( cudaFreeArray(*occupy_array) );
		*occupy_array = 0;
	}

	void cu_Copy_to_tex_inputVelocity_regular(const float4* vel, cudaArray** inputVelocity_array)
	{
		tex_inputVelocity_regular.normalized = true;
		tex_inputVelocity_regular.filterMode = cudaFilterModeLinear;
		tex_inputVelocity_regular.addressMode[0] = cudaAddressModeClamp;
		tex_inputVelocity_regular.addressMode[1] = cudaAddressModeClamp;
		tex_inputVelocity_regular.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf4 = cudaCreateChannelDesc<float4>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height,h_depth);

		checkCudaErrors( cudaMalloc3DArray(inputVelocity_array, &channelDescf4, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)vel, texSize.width*sizeof(float4), texSize.width, texSize.height);
		copyParams.dstArray = *inputVelocity_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_inputVelocity_regular,*inputVelocity_array,channelDescf4) );
	}

	void cu_Free_tex_inputVelocity_regular(cudaArray** inputVelocity_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_inputVelocity_regular) );
		checkCudaErrors( cudaFreeArray(*inputVelocity_array) );
		*inputVelocity_array = 0;
	}

	void cu_Copy_to_tex_inputVelocity_MAC_u(const float* u, cudaArray** u_array)
	{
		tex_inputVelocity_MAC_u.normalized = true;
		tex_inputVelocity_MAC_u.filterMode = cudaFilterModeLinear;
		tex_inputVelocity_MAC_u.addressMode[0] = cudaAddressModeClamp;
		tex_inputVelocity_MAC_u.addressMode[1] = cudaAddressModeClamp;
		tex_inputVelocity_MAC_u.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width+1,h_height,h_depth);

		checkCudaErrors( cudaMalloc3DArray(u_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)u, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *u_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_inputVelocity_MAC_u,*u_array,channelDescf) );
	}

	void cu_Free_tex_inputVelocity_MAC_u(cudaArray** u_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_inputVelocity_MAC_u) );
		checkCudaErrors( cudaFreeArray(*u_array) );
		*u_array = 0;
	}

	void cu_Copy_to_tex_inputVelocity_MAC_v(const float* v, cudaArray** v_array)
	{
		tex_inputVelocity_MAC_v.normalized = true;
		tex_inputVelocity_MAC_v.filterMode = cudaFilterModeLinear;
		tex_inputVelocity_MAC_v.addressMode[0] = cudaAddressModeClamp;
		tex_inputVelocity_MAC_v.addressMode[1] = cudaAddressModeClamp;
		tex_inputVelocity_MAC_v.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height+1,h_depth);

		checkCudaErrors( cudaMalloc3DArray(v_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)v, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *v_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_inputVelocity_MAC_v,*v_array,channelDescf) );
	}

	void cu_Free_tex_inputVelocity_MAC_v(cudaArray** v_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_inputVelocity_MAC_v) );
		checkCudaErrors( cudaFreeArray(*v_array) );
		*v_array = 0;
	}

	void cu_Copy_to_tex_inputVelocity_MAC_w(const float* w, cudaArray** w_array)
	{
		tex_inputVelocity_MAC_w.normalized = true;
		tex_inputVelocity_MAC_w.filterMode = cudaFilterModeLinear;
		tex_inputVelocity_MAC_w.addressMode[0] = cudaAddressModeClamp;
		tex_inputVelocity_MAC_w.addressMode[1] = cudaAddressModeClamp;
		tex_inputVelocity_MAC_w.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height,h_depth+1);

		checkCudaErrors( cudaMalloc3DArray(w_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)w, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *w_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_inputVelocity_MAC_w,*w_array,channelDescf) );
	}

	void cu_Free_tex_inputVelocity_MAC_w(cudaArray** w_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_inputVelocity_MAC_w) );
		checkCudaErrors( cudaFreeArray(*w_array) );
		*w_array = 0;
	}

	void cu_Copy_to_tex_temperature(const float* temperature, cudaArray** temperature_array)
	{
		tex_temperature.normalized = true;
		tex_temperature.filterMode = cudaFilterModeLinear;
		tex_temperature.addressMode[0] = cudaAddressModeClamp;
		tex_temperature.addressMode[1] = cudaAddressModeClamp;
		tex_temperature.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height,h_depth);

		checkCudaErrors( cudaMalloc3DArray(temperature_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)temperature, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *temperature_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_temperature,*temperature_array,channelDescf) );
	}

	void cu_Free_tex_temperature(cudaArray** temperature_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_temperature) );
		checkCudaErrors( cudaFreeArray(*temperature_array) );
		*temperature_array = 0;
	}

	void cu_Copy_to_tex_density(const float* density, cudaArray** density_array)
	{
		tex_density.normalized = true;
		tex_density.filterMode = cudaFilterModeLinear;
		tex_density.addressMode[0] = cudaAddressModeClamp;
		tex_density.addressMode[1] = cudaAddressModeClamp;
		tex_density.addressMode[2] = cudaAddressModeClamp;

		cudaChannelFormatDesc channelDescf = cudaCreateChannelDesc<float>();
		cudaExtent texSize = make_cudaExtent(h_width,h_height,h_depth);

		checkCudaErrors( cudaMalloc3DArray(density_array, &channelDescf, texSize) );

		// copy data to 3D array
		cudaMemcpy3DParms copyParams = {0};
		copyParams.srcPtr   = make_cudaPitchedPtr((void*)density, texSize.width*sizeof(float), texSize.width, texSize.height);
		copyParams.dstArray = *density_array;
		copyParams.extent   = texSize;
		copyParams.kind     = cudaMemcpyDeviceToDevice;
		checkCudaErrors( cudaMemcpy3D(&copyParams) );

		checkCudaErrors( cudaBindTextureToArray(tex_density,*density_array,channelDescf) );
	}

	void cu_Free_tex_density(cudaArray** density_array)
	{
		checkCudaErrors( cudaUnbindTexture(tex_density) );
		checkCudaErrors( cudaFreeArray(*density_array) );
		*density_array = 0;
	}

	/**********************************************************************************/

	void cu_Velocity_Negative(float* u, float* v, float* w, int width, int height, int depth)
	{
		int u_len = (width+1)*height*depth;
		int v_len = width*(height+1)*depth;
		int w_len = width*height*(depth+1);
		dim3 blockSize(BLOCK_SIZE*BLOCK_SIZE,1);
		dim3 gridSize_u((u_len+blockSize.x-1)/blockSize.x,1);
		dim3 gridSize_v((v_len+blockSize.x-1)/blockSize.x,1);
		dim3 gridSize_w((w_len+blockSize.x-1)/blockSize.x,1);

		Velocity_Negative_Kernel<<<gridSize_u,blockSize>>>(u,u_len);
		Velocity_Negative_Kernel<<<gridSize_v,blockSize>>>(v,v_len);
		Velocity_Negative_Kernel<<<gridSize_w,blockSize>>>(w,w_len);
	}

	void cu_Velocity_Negative_4channels(float4* vel, int width, int height, int depth)
	{
		int vel_len = width*height*depth;
		dim3 blockSize(BLOCK_SIZE*BLOCK_SIZE,1);
		dim3 gridSize((vel_len+blockSize.x-1)/blockSize.x,1);

		Velocity_Negative_4channels_Kernel<<<gridSize,blockSize>>>(vel,vel_len);
	}

	void cu_Input_Increment(float* input, const float* input_star, int len)
	{
		dim3 blockSize(BLOCK_SIZE*BLOCK_SIZE,1);
		dim3 gridSize((len+blockSize.x-1)/blockSize.x,1);

		Input_Increment_Kernel<<<gridSize,blockSize>>>(input,input_star,len);
	}

	void cu_Input_Increment_4channels(float4* input, const float4* input_star, int len)
	{
		dim3 blockSize(BLOCK_SIZE*BLOCK_SIZE,1);
		dim3 gridSize((len+blockSize.x-1)/blockSize.x,1);

		Input_Increment_4channels_Kernel<<<gridSize,blockSize>>>(input,input_star,len);
	}

	void cu_Apply_Advect_Velocity_Result_Open(const float* adv_u, const float* adv_v, const float* adv_w, const bool* occupy, float* u, float* v, float* w)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize_u(((h_width+1)+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		dim3 gridSize_v((h_width+blockSize.x-1)/blockSize.x,((h_height+1)+blockSize.y-1)/blockSize.y);
		dim3 gridSize_w((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);

		 Apply_Advect_Velocity_Result_Open_u_Kernel<<<gridSize_u,blockSize>>>(adv_u,occupy,u,h_width,h_height,h_depth);
		 Apply_Advect_Velocity_Result_Open_v_Kernel<<<gridSize_v,blockSize>>>(adv_v,occupy,v,h_width,h_height,h_depth);
		 Apply_Advect_Velocity_Result_Open_w_Kernel<<<gridSize_w,blockSize>>>(adv_w,occupy,w,h_width,h_height,h_depth);
	}

	void cu_Apply_Advect_Velocity_Result_Closed(const float* adv_u, const float* adv_v, const float* adv_w, const bool* occupy, float* u, float* v, float* w)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize_u(((h_width+1)+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		dim3 gridSize_v((h_width+blockSize.x-1)/blockSize.x,((h_height+1)+blockSize.y-1)/blockSize.y);
		dim3 gridSize_w((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);

		 Apply_Advect_Velocity_Result_Closed_u_Kernel<<<gridSize_u,blockSize>>>(adv_u,occupy,u,h_width,h_height,h_depth);
		 Apply_Advect_Velocity_Result_Closed_v_Kernel<<<gridSize_v,blockSize>>>(adv_v,occupy,v,h_width,h_height,h_depth);
		 Apply_Advect_Velocity_Result_Closed_w_Kernel<<<gridSize_w,blockSize>>>(adv_w,occupy,w,h_width,h_height,h_depth);
	}

	void cu_Advect_Velocity_inRegular_outRegular(float* u, float* v, float* w, const bool* occupy, bool is_open)
	{
		float4* vel = 0;
		checkCudaErrors( cudaMalloc((void**)&vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(u,v,w,vel,h_width,h_height,h_depth);
		
		cudaArray* velocity_array = 0;
		cudaArray* inputVelocity_array = 0;
		cudaArray* occupy_array = 0;
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_inputVelocity_regular(vel,&inputVelocity_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Velocity_inRegular_outRegular_Kernel<<<gridSize,blockSize>>>(vel);

		float* out_u = 0;
		float* out_v = 0;
		float* out_w = 0;
		checkCudaErrors( cudaMalloc((void**)&out_u,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMemset(out_u,0,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_v,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMemset(out_v,0,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_w,sizeof(float)*h_width*h_height*(h_depth+1)) );
		checkCudaErrors( cudaMemset(out_w,0,sizeof(float)*h_width*h_height*(h_depth+1)) );
		
		ZQ_CUDA_MACtoRegular::cu_Regular4_to_MAC(vel,out_u,out_v,out_w,h_width,h_height,h_depth);

		if(is_open)
			cu_Apply_Advect_Velocity_Result_Open(out_u,out_v,out_w,occupy,u,v,w);
		else
			cu_Apply_Advect_Velocity_Result_Closed(out_u,out_v,out_w,occupy,u,v,w);

		checkCudaErrors( cudaFree(vel) );
		checkCudaErrors( cudaFree(out_u) );
		checkCudaErrors( cudaFree(out_v) );
		checkCudaErrors( cudaFree(out_w) );
		vel = 0;
		out_u = 0;
		out_v = 0;
		out_w = 0;

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);
		cu_Free_tex_occupy(&occupy_array);
	}

	void cu_Advect_Velocity_inRegular_outRegular_BFECC(float* u, float* v, float* w, const bool* occupy, bool is_open)
	{
		float4* vel = 0;
		checkCudaErrors( cudaMalloc((void**)&vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(u,v,w,vel,h_width,h_height,h_depth);

		float4* in_out_vel = 0;
		checkCudaErrors( cudaMalloc((void**)&in_out_vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(in_out_vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		
		cudaArray* velocity_array = 0;
		cudaArray* inputVelocity_array = 0;
		cudaArray* occupy_array = 0;
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_inputVelocity_regular(vel,&inputVelocity_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Velocity_inRegular_outRegular_Kernel<<<gridSize,blockSize>>>(in_out_vel);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);

		cu_Velocity_Negative_4channels(vel,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_inputVelocity_regular(in_out_vel,&inputVelocity_array);
		Advect_Velocity_inRegular_outRegular_Kernel<<<gridSize,blockSize>>>(in_out_vel);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);

		cu_Velocity_Negative_4channels(vel,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Input_Increment_4channels(vel,in_out_vel,h_width*h_height*h_depth);
		cu_Copy_to_tex_inputVelocity_regular(vel,&inputVelocity_array);

		Advect_Velocity_inRegular_outRegular_Kernel<<<gridSize,blockSize>>>(in_out_vel);


		float* out_u = 0;
		float* out_v = 0;
		float* out_w = 0;
		checkCudaErrors( cudaMalloc((void**)&out_u,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMemset(out_u,0,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_v,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMemset(out_v,0,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_w,sizeof(float)*h_width*h_height*(h_depth+1)) );
		checkCudaErrors( cudaMemset(out_w,0,sizeof(float)*h_width*h_height*(h_depth+1)) );
		
		ZQ_CUDA_MACtoRegular::cu_Regular4_to_MAC(in_out_vel,out_u,out_v,out_w,h_width,h_height,h_depth);

		if(is_open)
			cu_Apply_Advect_Velocity_Result_Open(out_u,out_v,out_w,occupy,u,v,w);
		else
			cu_Apply_Advect_Velocity_Result_Closed(out_u,out_v,out_w,occupy,u,v,w);

		checkCudaErrors( cudaFree(vel) );
		checkCudaErrors( cudaFree(in_out_vel));
		checkCudaErrors( cudaFree(out_u) );
		checkCudaErrors( cudaFree(out_v) );
		checkCudaErrors( cudaFree(out_w) );
		vel = 0;
		in_out_vel = 0;
		out_u = 0;
		out_v = 0;
		out_w = 0;

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);
		cu_Free_tex_occupy(&occupy_array);
	}

	void cu_Advect_Velocity_inMAC_outMAC(float* u, float* v, float* w, const bool* occupy, bool is_open)
	{
		cudaArray* velocity_MAC_u_array = 0;
		cudaArray* velocity_MAC_v_array = 0;
		cudaArray* velocity_MAC_w_array = 0;
		cudaArray* occupy_array = 0;
		cudaArray* inputVelocity_MAC_u_array = 0;
		cudaArray* inputVelocity_MAC_v_array = 0;
		cudaArray* inputVelocity_MAC_w_array = 0;
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);
		cu_Copy_to_tex_inputVelocity_MAC_u(u,&inputVelocity_MAC_u_array);
		cu_Copy_to_tex_inputVelocity_MAC_v(v,&inputVelocity_MAC_v_array);
		cu_Copy_to_tex_inputVelocity_MAC_w(w,&inputVelocity_MAC_w_array);

		float* out_u = 0;
		float* out_v = 0;
		float* out_w = 0;
		checkCudaErrors( cudaMalloc((void**)&out_u,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMemset(out_u,0,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_v,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMemset(out_v,0,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_w,sizeof(float)*h_width*h_height*(h_depth+1)) );
		checkCudaErrors( cudaMemset(out_w,0,sizeof(float)*h_width*h_height*(h_depth+1)) );

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize_u(((h_width+1)+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		dim3 gridSize_v((h_width+blockSize.x-1)/blockSize.x,((h_height+1)+blockSize.y-1)/blockSize.y);
		dim3 gridSize_w((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Velocity_inMAC_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inMAC_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inMAC_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);
		
		if(is_open)
			cu_Apply_Advect_Velocity_Result_Open(out_u,out_v,out_w,occupy,u,v,w);
		else
			cu_Apply_Advect_Velocity_Result_Closed(out_u,out_v,out_w,occupy,u,v,w);

		checkCudaErrors( cudaFree(out_u) );
		checkCudaErrors( cudaFree(out_v) );
		checkCudaErrors( cudaFree(out_w) );
		out_u = 0;
		out_v = 0;
		out_w = 0;

		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_occupy(&occupy_array);
		cu_Free_tex_inputVelocity_MAC_u(&inputVelocity_MAC_u_array);
		cu_Free_tex_inputVelocity_MAC_v(&inputVelocity_MAC_v_array);
		cu_Free_tex_inputVelocity_MAC_w(&inputVelocity_MAC_w_array);
	}

	void cu_Advect_Velocity_inMAC_outMAC_BFECC(float* u, float* v, float* w, const bool* occupy, bool is_open)
	{
		cudaArray* velocity_MAC_u_array = 0;
		cudaArray* velocity_MAC_v_array = 0;
		cudaArray* velocity_MAC_w_array = 0;
		cudaArray* occupy_array = 0;
		cudaArray* inputVelocity_MAC_u_array = 0;
		cudaArray* inputVelocity_MAC_v_array = 0;
		cudaArray* inputVelocity_MAC_w_array = 0;
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);
		cu_Copy_to_tex_inputVelocity_MAC_u(u,&inputVelocity_MAC_u_array);
		cu_Copy_to_tex_inputVelocity_MAC_v(v,&inputVelocity_MAC_v_array);
		cu_Copy_to_tex_inputVelocity_MAC_w(w,&inputVelocity_MAC_w_array);

		float* out_u = 0;
		float* out_v = 0;
		float* out_w = 0;
		checkCudaErrors( cudaMalloc((void**)&out_u,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMemset(out_u,0,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_v,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMemset(out_v,0,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_w,sizeof(float)*h_width*h_height*(h_depth+1)) );
		checkCudaErrors( cudaMemset(out_w,0,sizeof(float)*h_width*h_height*(h_depth+1)) );

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize_u(((h_width+1)+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		dim3 gridSize_v((h_width+blockSize.x-1)/blockSize.x,((h_height+1)+blockSize.y-1)/blockSize.y);
		dim3 gridSize_w((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Velocity_inMAC_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inMAC_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inMAC_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);

		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_inputVelocity_MAC_u(&inputVelocity_MAC_u_array);
		cu_Free_tex_inputVelocity_MAC_v(&inputVelocity_MAC_v_array);
		cu_Free_tex_inputVelocity_MAC_w(&inputVelocity_MAC_w_array);

		cu_Velocity_Negative(u,v,w,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Copy_to_tex_inputVelocity_MAC_u(out_u,&inputVelocity_MAC_u_array);
		cu_Copy_to_tex_inputVelocity_MAC_v(out_v,&inputVelocity_MAC_v_array);
		cu_Copy_to_tex_inputVelocity_MAC_w(out_w,&inputVelocity_MAC_w_array);

		Advect_Velocity_inMAC_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inMAC_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inMAC_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);

		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_inputVelocity_MAC_u(&inputVelocity_MAC_u_array);
		cu_Free_tex_inputVelocity_MAC_v(&inputVelocity_MAC_v_array);
		cu_Free_tex_inputVelocity_MAC_w(&inputVelocity_MAC_w_array);

		cu_Velocity_Negative(u,v,w,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Input_Increment(u,out_u,(h_width+1)*h_height*h_depth);
		cu_Input_Increment(v,out_v,h_width*(h_height+1)*h_depth);
		cu_Input_Increment(w,out_w,h_width*h_height*(h_depth+1));
		cu_Copy_to_tex_inputVelocity_MAC_u(u,&inputVelocity_MAC_u_array);
		cu_Copy_to_tex_inputVelocity_MAC_v(v,&inputVelocity_MAC_v_array);
		cu_Copy_to_tex_inputVelocity_MAC_w(w,&inputVelocity_MAC_w_array);

		Advect_Velocity_inMAC_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inMAC_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inMAC_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);
		
		if(is_open)
			cu_Apply_Advect_Velocity_Result_Open(out_u,out_v,out_w,occupy,u,v,w);
		else
			cu_Apply_Advect_Velocity_Result_Closed(out_u,out_v,out_w,occupy,u,v,w);

		checkCudaErrors( cudaFree(out_u) );
		checkCudaErrors( cudaFree(out_v) );
		checkCudaErrors( cudaFree(out_w) );
		out_u = 0;
		out_v = 0;
		out_w = 0;

		cu_Free_tex_occupy(&occupy_array);
		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_inputVelocity_MAC_u(&inputVelocity_MAC_u_array);
		cu_Free_tex_inputVelocity_MAC_v(&inputVelocity_MAC_v_array);
		cu_Free_tex_inputVelocity_MAC_w(&inputVelocity_MAC_w_array);

	}

	void cu_Advect_Velocity_inRegular_outMAC(float* u, float* v, float* w, const bool* occupy, bool is_open)
	{
		float4* vel = 0;
		checkCudaErrors( cudaMalloc((void**)&vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(u,v,w,vel,h_width,h_height,h_depth);

		cudaArray* velocity_array = 0;
		cudaArray* inputVelocity_array = 0;
		cudaArray* occupy_array = 0;
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_inputVelocity_regular(vel,&inputVelocity_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);

		float* out_u = 0;
		float* out_v = 0;
		float* out_w = 0;
		checkCudaErrors( cudaMalloc((void**)&out_u,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMemset(out_u,0,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_v,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMemset(out_v,0,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_w,sizeof(float)*h_width*h_height*(h_depth+1)) );
		checkCudaErrors( cudaMemset(out_w,0,sizeof(float)*h_width*h_height*(h_depth+1)) );

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize_u(((h_width+1)+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		dim3 gridSize_v((h_width+blockSize.x-1)/blockSize.x,((h_height+1)+blockSize.y-1)/blockSize.y);
		dim3 gridSize_w((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Velocity_inRegular_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inRegular_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inRegular_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);
		
		if(is_open)
			cu_Apply_Advect_Velocity_Result_Open(out_u,out_v,out_w,occupy,u,v,w);
		else
			cu_Apply_Advect_Velocity_Result_Closed(out_u,out_v,out_w,occupy,u,v,w);

		checkCudaErrors( cudaFree(vel) );
		checkCudaErrors( cudaFree(out_u) );
		checkCudaErrors( cudaFree(out_v) );
		checkCudaErrors( cudaFree(out_w) );
		vel = 0;
		out_u = 0;
		out_v = 0;
		out_w = 0;

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);
		cu_Free_tex_occupy(&occupy_array);
	}

	void cu_Advect_Velocity_inRegular_outMAC_BFECC(float* u, float* v, float* w, const bool* occupy, bool is_open)
	{
		float4* vel = 0;
		float4* in_out_vel = 0;
		checkCudaErrors( cudaMalloc((void**)&vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&in_out_vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(in_out_vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(u,v,w,vel,h_width,h_height,h_depth);

		cudaArray* velocity_array = 0;
		cudaArray* inputVelocity_array = 0;
		cudaArray* occupy_array = 0;
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_inputVelocity_regular(vel,&inputVelocity_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);

		float* out_u = 0;
		float* out_v = 0;
		float* out_w = 0;
		checkCudaErrors( cudaMalloc((void**)&out_u,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMemset(out_u,0,sizeof(float)*(h_width+1)*h_height*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_v,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMemset(out_v,0,sizeof(float)*h_width*(h_height+1)*h_depth) );
		checkCudaErrors( cudaMalloc((void**)&out_w,sizeof(float)*h_width*h_height*(h_depth+1)) );
		checkCudaErrors( cudaMemset(out_w,0,sizeof(float)*h_width*h_height*(h_depth+1)) );

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize_u(((h_width+1)+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		dim3 gridSize_v((h_width+blockSize.x-1)/blockSize.x,((h_height+1)+blockSize.y-1)/blockSize.y);
		dim3 gridSize_w((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Velocity_inRegular_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inRegular_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inRegular_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);
		
		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(out_u,out_v,out_w,in_out_vel,h_width,h_height,h_depth);

		cu_Velocity_Negative_4channels(vel,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_inputVelocity_regular(in_out_vel,&inputVelocity_array);

		Advect_Velocity_inRegular_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inRegular_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inRegular_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);

		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(out_u,out_v,out_w,in_out_vel,h_width,h_height,h_depth);

		cu_Velocity_Negative_4channels(vel,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		
		cu_Input_Increment_4channels(vel,in_out_vel,h_width*h_height*h_depth);

		cu_Copy_to_tex_inputVelocity_regular(vel,&inputVelocity_array);

		Advect_Velocity_inRegular_outMAC_u_Kernel<<<gridSize_u,blockSize>>>(out_u);
		Advect_Velocity_inRegular_outMAC_v_Kernel<<<gridSize_v,blockSize>>>(out_v);
		Advect_Velocity_inRegular_outMAC_w_Kernel<<<gridSize_w,blockSize>>>(out_w);


		
		if(is_open)
			cu_Apply_Advect_Velocity_Result_Open(out_u,out_v,out_w,occupy,u,v,w);
		else
			cu_Apply_Advect_Velocity_Result_Closed(out_u,out_v,out_w,occupy,u,v,w);

		checkCudaErrors( cudaFree(vel) );
		checkCudaErrors( cudaFree(in_out_vel) );
		checkCudaErrors( cudaFree(out_u) );
		checkCudaErrors( cudaFree(out_v) );
		checkCudaErrors( cudaFree(out_w) );
		vel = 0;
		in_out_vel = 0;
		out_u = 0;
		out_v = 0;
		out_w = 0;

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_inputVelocity_regular(&inputVelocity_array);
		cu_Free_tex_occupy(&occupy_array);
	}

	void cu_Advect_Scalar_Regular_Velocity(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density, 
								float* output_temperature, float* output_density)
	{
		
		float4* vel = 0;
		checkCudaErrors( cudaMalloc((void**)&vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(u,v,w,vel,h_width,h_height,h_depth);

		cudaArray* velocity_array = 0;
		cudaArray* occupy_array = 0;
		cudaArray* temperature_array = 0;
		cudaArray* density_array = 0;
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);
		cu_Copy_to_tex_temperature(input_temperature,&temperature_array);
		cu_Copy_to_tex_density(input_density,&density_array);

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Scalar_Regular_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_occupy(&occupy_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);

		checkCudaErrors( cudaFree(vel));
		vel = 0;
	}

	void cu_Advect_Scalar_Regular_Velocity_BFECC(const float* u, const float* v, const float* w, const bool* occupy, float* input_temperature, float* input_density, 
								float* output_temperature, float* output_density)
	{
		
		float4* vel = 0;
		checkCudaErrors( cudaMalloc((void**)&vel,sizeof(float4)*h_width*h_height*h_depth) );
		checkCudaErrors( cudaMemset(vel,0,sizeof(float4)*h_width*h_height*h_depth) );
		ZQ_CUDA_MACtoRegular::cu_MAC_to_Regular4(u,v,w,vel,h_width,h_height,h_depth);

		cudaArray* velocity_array = 0;
		cudaArray* occupy_array = 0;
		cudaArray* temperature_array = 0;
		cudaArray* density_array = 0;
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);
		cu_Copy_to_tex_temperature(input_temperature,&temperature_array);
		cu_Copy_to_tex_density(input_density,&density_array);

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Scalar_Regular_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);

		cu_Velocity_Negative_4channels(vel,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Copy_to_tex_temperature(output_temperature,&temperature_array);
		cu_Copy_to_tex_density(output_density,&density_array);

		Advect_Scalar_Regular_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);

		cu_Velocity_Negative_4channels(vel,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_regular(vel,&velocity_array);
		cu_Input_Increment(input_temperature,output_temperature,h_width*h_height*h_depth);
		cu_Input_Increment(input_density,output_density,h_width*h_height*h_depth);
		cu_Copy_to_tex_temperature(input_temperature,&temperature_array);
		cu_Copy_to_tex_density(input_density,&density_array);

		Advect_Scalar_Regular_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_regular(&velocity_array);
		cu_Free_tex_occupy(&occupy_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);

		checkCudaErrors( cudaFree(vel));
		vel = 0;
	}

	void cu_Advect_Scalar_MAC_Velocity(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density, 
								float* output_temperature, float* output_density)
	{
		cudaArray* velocity_MAC_u_array = 0;
		cudaArray* velocity_MAC_v_array = 0;
		cudaArray* velocity_MAC_w_array = 0;
		cudaArray* occupy_array = 0;
		cudaArray* temperature_array = 0;
		cudaArray* density_array = 0;
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);
		cu_Copy_to_tex_temperature(input_temperature,&temperature_array);
		cu_Copy_to_tex_density(input_density,&density_array);

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Scalar_MAC_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_occupy(&occupy_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);
	}

	void cu_Advect_Scalar_MAC_Velocity_BFECC(float* u, float* v, float* w, const bool* occupy, float* input_temperature, float* input_density, 
								float* output_temperature, float* output_density)
	{
		cudaArray* velocity_MAC_u_array = 0;
		cudaArray* velocity_MAC_v_array = 0;
		cudaArray* velocity_MAC_w_array = 0;
		cudaArray* occupy_array = 0;
		cudaArray* temperature_array = 0;
		cudaArray* density_array = 0;
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Copy_to_tex_occupy(occupy,&occupy_array);
		cu_Copy_to_tex_temperature(input_temperature,&temperature_array);
		cu_Copy_to_tex_density(input_density,&density_array);

		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((h_width+blockSize.x-1)/blockSize.x,(h_height+blockSize.y-1)/blockSize.y);
		Advect_Scalar_MAC_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);

		cu_Velocity_Negative(u,v,w,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Copy_to_tex_temperature(output_temperature,&temperature_array);
		cu_Copy_to_tex_density(output_density,&density_array);

		Advect_Scalar_MAC_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);

		cu_Velocity_Negative(u,v,w,h_width,h_height,h_depth);
		cu_Copy_to_tex_velocity_MAC_u(u,&velocity_MAC_u_array);
		cu_Copy_to_tex_velocity_MAC_v(v,&velocity_MAC_v_array);
		cu_Copy_to_tex_velocity_MAC_w(w,&velocity_MAC_w_array);
		cu_Input_Increment(input_temperature,output_temperature,h_width*h_height*h_depth);
		cu_Input_Increment(input_density,output_density,h_width*h_height*h_depth);
		cu_Copy_to_tex_temperature(input_temperature,&temperature_array);
		cu_Copy_to_tex_density(input_density,&density_array);

		Advect_Scalar_MAC_Velocity_Kernel<<<gridSize,blockSize>>>(output_temperature,output_density);

		cu_Free_tex_velocity_MAC_u(&velocity_MAC_u_array);
		cu_Free_tex_velocity_MAC_v(&velocity_MAC_v_array);
		cu_Free_tex_velocity_MAC_w(&velocity_MAC_w_array);
		cu_Free_tex_occupy(&occupy_array);
		cu_Free_tex_temperature(&temperature_array);
		cu_Free_tex_density(&density_array);
	}

	/****************************************************************************************/

	extern "C"
	void ZQ_Cuda_Prepare_Advection(const unsigned int width, const unsigned int height, const unsigned int depth, const float voxelSize, const unsigned int steps, const float deltatt)
	{
		h_width = width;
		h_height = height;
		h_depth = depth;
		h_steps = steps;
		h_voxelSize = voxelSize;
		h_deltatt = deltatt;

		checkCudaErrors( cudaMemcpyToSymbol(ZQ_CUDA_Advection3D::d_width,&width,sizeof(unsigned int)));
		checkCudaErrors( cudaMemcpyToSymbol(ZQ_CUDA_Advection3D::d_height,&height,sizeof(unsigned int)));
		checkCudaErrors( cudaMemcpyToSymbol(ZQ_CUDA_Advection3D::d_depth,&depth,sizeof(unsigned int)));
		checkCudaErrors( cudaMemcpyToSymbol(ZQ_CUDA_Advection3D::d_steps,&steps,sizeof(unsigned int)));
		checkCudaErrors( cudaMemcpyToSymbol(ZQ_CUDA_Advection3D::d_voxelSize,&voxelSize,sizeof(float)));
		checkCudaErrors( cudaMemcpyToSymbol(ZQ_CUDA_Advection3D::d_deltatt,&deltatt,sizeof(float)));	

		//int tmp = 0;
		//scanf("%d",&tmp);
		
		//checkCudaErrors( cudaMemcpyFromSymbol(&h_width,d_width,sizeof(unsigned int)));
		//checkCudaErrors( cudaMemcpyFromSymbol(&h_height,d_height,sizeof(unsigned int)));
		//checkCudaErrors( cudaMemcpyFromSymbol(&h_steps,d_steps,sizeof(unsigned int)));
		//checkCudaErrors( cudaMemcpyFromSymbol(&h_voxelSize,d_voxelSize,sizeof(float)));
		//checkCudaErrors( cudaMemcpyFromSymbol(&h_deltatt,d_deltatt,sizeof(float)));
		//printf("width = %d\n",h_width);
		//printf("height = %d\n",h_height);
		//printf("steps = %d\n",h_steps);
		//printf("voxelSize = %f\n",h_voxelSize);
		//printf("deltatt = %f\n",h_deltatt);
		
		//scanf("%d",&tmp);
	}

	extern "C"
	float Advect_Velocity(float* u, float* v, float* w, const bool* occupy, bool is_open, enum AdvectVelocityType type)
	{
		float time = 0;
		cudaEvent_t start,stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start,0);

		/////
		float* d_u = 0;
		float* d_v = 0;
		float* d_w = 0;
		bool* d_occupy = 0;
		checkCudaErrors( cudaMalloc((void**)&d_u,sizeof(float)*(h_width+1)*h_height*h_depth));
		checkCudaErrors( cudaMalloc((void**)&d_v,sizeof(float)*h_width*(h_height+1)*h_depth));
		checkCudaErrors( cudaMalloc((void**)&d_w,sizeof(float)*h_width*h_height*(h_depth+1)));
		checkCudaErrors( cudaMalloc((void**)&d_occupy,sizeof(bool)*h_width*h_height*h_depth));
		checkCudaErrors( cudaMemcpy(d_u,u,sizeof(float)*(h_width+1)*h_height*h_depth,cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(d_v,v,sizeof(float)*h_width*(h_height+1)*h_depth,cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(d_w,w,sizeof(float)*h_width*h_height*(h_depth+1),cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMemcpy(d_occupy,occupy,sizeof(bool)*h_width*h_height*h_depth,cudaMemcpyHostToDevice));

		switch(type)
		{
		case ADV_VEL_INREG_OUTREG:
			cu_Advect_Velocity_inRegular_outRegular(d_u,d_v,d_w,d_occupy,is_open);
			break;
		case ADV_VEL_INREG_OUTREG_BFECC:
			cu_Advect_Velocity_inRegular_outRegular_BFECC(d_u,d_v,d_w,d_occupy,is_open);
			break;
		case ADV_VEL_INMAC_OUTMAC:
			cu_Advect_Velocity_inMAC_outMAC(d_u,d_v,d_w,d_occupy,is_open);
			break;
		case ADV_VEL_INMAC_OUTMAC_BFECC:
			cu_Advect_Velocity_inMAC_outMAC_BFECC(d_u,d_v,d_w,d_occupy,is_open);
			break;
		case ADV_VEL_INREG_OUTMAC:
			cu_Advect_Velocity_inRegular_outMAC(d_u,d_v,d_w,d_occupy,is_open);
			break;
		case ADV_VEL_INREG_OUTMAC_BFECC:
			cu_Advect_Velocity_inRegular_outMAC_BFECC(d_u,d_v,d_w,d_occupy,is_open);
			break;
		}
		
		checkCudaErrors( cudaMemcpy(u,d_u,sizeof(float)*(h_width+1)*h_height*h_depth,cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(v,d_v,sizeof(float)*h_width*(h_height+1)*h_depth,cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(w,d_w,sizeof(float)*h_width*h_height*(h_depth+1),cudaMemcpyDeviceToHost) );

		checkCudaErrors( cudaFree(d_u) );
		checkCudaErrors( cudaFree(d_v) );
		checkCudaErrors( cudaFree(d_w) );
		checkCudaErrors( cudaFree(d_occupy));
		d_u = 0;
		d_v = 0;
		d_w = 0;
		d_occupy = 0;
		
		cudaEventRecord(stop,0);
		cudaEventSynchronize(start);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time,start,stop);
		return time;
	}

	extern "C"
	float Advect_Scalar(const float* u, const float* v, const float* w, const bool* occupy, const float* input_temperature, const float* input_density, 
					float* output_temperature, float* output_density, enum AdvectScalarType type)
	{
		float time = 0;
		cudaEvent_t start,stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start,0);
		
		float* d_u = 0;
		float* d_v = 0;
		float* d_w = 0;
		bool* d_occupy = 0;
		float* d_input_temperature = 0;
		float* d_input_density = 0;
		float* d_output_temperature = 0;
		float* d_output_density = 0;
		checkCudaErrors( cudaMalloc((void**)&d_u,sizeof(float)*(h_width+1)*h_height*h_depth));
		checkCudaErrors( cudaMemcpy(d_u,u,sizeof(float)*(h_width+1)*h_height*h_depth,cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMalloc((void**)&d_v,sizeof(float)*h_width*(h_height+1)*h_depth));
		checkCudaErrors( cudaMemcpy(d_v,v,sizeof(float)*h_width*(h_height+1)*h_depth,cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMalloc((void**)&d_w,sizeof(float)*h_width*h_height*(h_depth+1)));
		checkCudaErrors( cudaMemcpy(d_w,w,sizeof(float)*h_width*h_height*(h_depth+1),cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMalloc((void**)&d_occupy,sizeof(bool)*h_width*h_height*h_depth));
		checkCudaErrors( cudaMemcpy(d_occupy,occupy,sizeof(bool)*h_width*h_height*h_depth,cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMalloc((void**)&d_input_temperature,sizeof(float)*h_width*h_height*h_depth));
		checkCudaErrors( cudaMemcpy(d_input_temperature,input_temperature,sizeof(float)*h_width*h_height*h_depth,cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMalloc((void**)&d_input_density,sizeof(float)*h_width*h_height*h_depth));
		checkCudaErrors( cudaMemcpy(d_input_density,input_density,sizeof(float)*h_width*h_height*h_depth,cudaMemcpyHostToDevice));
		checkCudaErrors( cudaMalloc((void**)&d_output_temperature,sizeof(float)*h_width*h_height*h_depth));
		checkCudaErrors( cudaMemset(d_output_temperature,0,sizeof(float)*h_width*h_height*h_depth));
		checkCudaErrors( cudaMalloc((void**)&d_output_density,sizeof(float)*h_width*h_height*h_depth));
		checkCudaErrors( cudaMemset(d_output_density,0,sizeof(float)*h_width*h_height*h_depth));

		switch(type)
		{
		case ADV_SCA_MAC:
			cu_Advect_Scalar_MAC_Velocity(d_u,d_v,d_w,d_occupy,d_input_temperature,d_input_density,d_output_temperature,d_output_density);
			break;
		case ADV_SCA_REG:
			cu_Advect_Scalar_Regular_Velocity(d_u,d_v,d_w,d_occupy,d_input_temperature,d_input_density,d_output_temperature,d_output_density);
			break;
		}
		
		checkCudaErrors( cudaMemcpy(output_temperature,d_output_temperature,sizeof(float)*h_width*h_height*h_depth,cudaMemcpyDeviceToHost));
		checkCudaErrors( cudaMemcpy(output_density,d_output_density,sizeof(float)*h_width*h_height*h_depth,cudaMemcpyDeviceToHost));

		checkCudaErrors( cudaFree(d_u));
		checkCudaErrors( cudaFree(d_v));
		checkCudaErrors( cudaFree(d_w));
		checkCudaErrors( cudaFree(d_occupy));
		checkCudaErrors( cudaFree(d_input_temperature));
		checkCudaErrors( cudaFree(d_input_density));
		checkCudaErrors( cudaFree(d_output_temperature));
		checkCudaErrors( cudaFree(d_output_density));
		d_u = 0;
		d_v = 0;
		d_w = 0;
		d_occupy = 0;
		d_input_temperature = 0;
		d_input_density = 0;
		d_output_temperature = 0;
		d_output_density = 0;
		
		cudaEventRecord(stop,0);
		cudaEventSynchronize(start);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time,start,stop);
		return time;
	}
}


#endif