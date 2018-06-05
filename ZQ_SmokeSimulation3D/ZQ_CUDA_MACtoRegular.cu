#ifndef _ZQ_CUDA_MAC_TO_REGULAR_CU_
#define _ZQ_CUDA_MAC_TO_REGULAR_CU_

#include "ZQlibCudaDefines.cuh"
#include "ZQ_CUDA_MACtoRegular.cuh"

namespace ZQ_CUDA_MACtoRegular
{
	__global__
	void MAC_to_Regular4_Kernel(const float* mac_u, const float* mac_v, const float* mac_w, float4* vel4, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x >= width || y >= height)
			return;

		for (int z = 0; z < depth; z++)
		{
			int offset = z*height*width + y*width + x;
			float u = 0.5f*(mac_u[z*height*(width + 1) + y*(width + 1) + x] + mac_u[z*height*(width + 1) + y*(width + 1) + x + 1]);
			float v = 0.5f*(mac_v[z*(height + 1)*width + y*width + x] + mac_v[z*(height + 1)*width + (y + 1)*width + x]);
			float w = 0.5f*(mac_w[z*height*width + y*width + x] + mac_w[(z + 1)*height*width + y*width + x]);
			vel4[offset] = make_float4(u, v, w, 0);
		}
	}

	__global__
		void Regular4_to_MAC_u_Kernel(const float4* vel4, float* mac_u, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x > width || y >= height)
			return;

		for (int z = 0; z < depth; z++)
		{
			if (x == 0)
				mac_u[z*height*(width + 1) + y*(width + 1) + x] = vel4[z*height*width + y*width + x].x;
			else if (x == width)
				mac_u[z*height*(width + 1) + y*(width + 1) + x] = vel4[z*height*width + y*width + x - 1].x;
			else
				mac_u[z*height*(width + 1) + y*(width + 1) + x] = 0.5*(vel4[z*height*width + y*width + x].x + vel4[z*height*width + y*width + x - 1].x);
		}
	}

	__global__
		void Regular4_to_MAC_v_Kernel(const float4* vel4, float* mac_v, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x >= width || y > height)
			return;

		for (int z = 0; z < depth; z++)
		{
			if (y == 0)
				mac_v[z*(height + 1)*width + y*width + x] = vel4[z*height*width + y*width + x].y;
			else if (y == height)
				mac_v[z*(height + 1)*width + y*width + x] = vel4[z*height*width + (y - 1)*width + x].y;
			else
				mac_v[z*(height + 1)*width + y*width + x] = 0.5*(vel4[z*height*width + y*width + x].y + vel4[z*height*width + (y - 1)*width + x].y);
		}
	}

	__global__
		void Regular4_to_MAC_w_Kernel(const float4* vel4, float* mac_w, int width, int height, int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x + tx;
		int y = by*blockDim.y + ty;
		if (x >= width || y >= height)
			return;

		mac_w[y*width + x] = vel4[y*width + x].z;
		mac_w[depth*height*width + y*width + x] = vel4[(depth - 1)*height*width + y*width + x].z;
		for (int z = 1; z < depth; z++)
		{
			mac_w[z*height*width + y*width + x] = 0.5*(vel4[z*height*width + y*width + x].z + vel4[(z - 1)*height*width + y*width + x].z);
		}
	}

	/******************/
	
	void cu_MAC_to_Regular4(const float* mac_u, const float* mac_v, const float* mac_w, float4* vel4, int width, int height, int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);

		MAC_to_Regular4_Kernel<<<gridSize,blockSize>>>(mac_u,mac_v,mac_w,vel4,width,height,depth);
	}

	void cu_Regular4_to_MAC(const float4* vel4, float* mac_u, float* mac_v, float* mac_w, int width, int height, int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 u_gridSize((width+1+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		dim3 v_gridSize((width+blockSize.x-1)/blockSize.x,(height+1+blockSize.y-1)/blockSize.y);
		dim3 w_gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		Regular4_to_MAC_u_Kernel<<<u_gridSize,blockSize>>>(vel4,mac_u,width,height,depth);
		Regular4_to_MAC_v_Kernel<<<v_gridSize,blockSize>>>(vel4,mac_v,width,height,depth);
		Regular4_to_MAC_w_Kernel<<<w_gridSize,blockSize>>>(vel4,mac_w,width,height,depth);
	}
}

#endif