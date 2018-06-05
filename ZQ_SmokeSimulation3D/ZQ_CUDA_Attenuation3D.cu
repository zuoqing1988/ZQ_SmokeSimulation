#ifndef _ZQ_CUDA_ATTENUATION_3D_CU_
#define _ZQ_CUDA_ATTENUATION_3D_CU_

#include "ZQlibCudaDefines.cuh"

namespace ZQ_CUDA_Attenuation3D
{
	__global__
	void Atten_u_Kernel(float* mac_u, const bool* occupy, const float velAtten, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x > width || y >= height)
			return ;
		
		for(int z = 0;z < depth;z++)
		{
			if(x == 0)
			{
				if(!occupy[z*height*width+y*width+0])
					mac_u[z*height*(width+1)+y*(width+1)+0] *= velAtten;
			}
			else if(x == width)
			{
				if(!occupy[z*height*width+y*width+width-1])
					mac_u[z*height*(width+1)+y*(width+1)+width] *= velAtten;
			}
			else
			{
				if(!occupy[z*height*width+y*width+x-1] && !occupy[z*height*width+y*width+x])
					mac_u[z*height*(width+1)+y*(width+1)+x] *= velAtten;
			}
		}
	}
	
	__global__
	void Atten_v_Kernel(float* mac_v, const bool* occupy, const float velAtten, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width || y > height)
			return ;
		
		for(int z = 0;z < depth;z++)
		{
			if(y == 0)
			{
				if(!occupy[z*height*width+x])
					mac_v[z*(height+1)*width+x] *= velAtten;
			}
			else if(y == height)
			{
				if(!occupy[z*height*width+(height-1)*width+x])
					mac_v[z*(height+1)*width+height*width+x] *= velAtten;
			}
			else
			{
				if(!occupy[z*height*width+(y-1)*width+x] && !occupy[z*height*width+y*width+x])
					mac_v[z*(height+1)*width+y*width+x] *= velAtten;
			}
		}
	}
	
	__global__
	void Atten_w_Kernel(float* mac_w, const bool* occupy, const float velAtten, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width || y >= height)
			return ;
		
		if(!occupy[y*width+x])
			mac_w[y*width+x] *= velAtten;
		if(!occupy[(depth-1)*height*width+y*width+x])
			mac_w[depth*height*width+y*width+x] *= velAtten;
		for(int z = 1;z < depth;z++)
		{
			if(!occupy[(z-1)*height*width+y*width+x] && !occupy[z*height*width+y*width+x])
				mac_w[z*height*width+y*width+x] *= velAtten;
		}
	}
	
	__global__
	void Atten_temperature_density_Kernel(float* temperature, float* density, const float tempAtten, const float densityAtten, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width || y >= height)
			return ;
			
		for(int z = 0;z < depth;z++)
		{
			temperature[z*height*width+y*width+x] *= tempAtten;
			density[z*height*width+y*width+x] *= densityAtten;
		}
	}

	/*************************************************/

	void cu_Attenuation3D(float* mac_u, float* mac_v, float* mac_w, float* temperature, float* density, const bool* occupy, 
						const float velAtten, const float tempAtten, const float densityAtten, const int width, const int height, const int depth)
	{	
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		dim3 u_gridSize((width+1+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		dim3 v_gridSize((width+blockSize.x-1)/blockSize.x,(height+1+blockSize.y-1)/blockSize.y);
		dim3 w_gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		
		Atten_u_Kernel<<<u_gridSize,blockSize>>>(mac_u,occupy,velAtten,width,height,depth);
		Atten_v_Kernel<<<v_gridSize,blockSize>>>(mac_v,occupy,velAtten,width,height,depth);
		Atten_w_Kernel<<<w_gridSize,blockSize>>>(mac_w,occupy,velAtten,width,height,depth);
		Atten_temperature_density_Kernel<<<gridSize,blockSize>>>(temperature,density,tempAtten,densityAtten,width,height,depth);
		
	}

	/***********************************************/

	extern "C"
	float Attenuation3D(float* mac_u, float* mac_v, float* mac_w, float* temperature, float* density, const bool* occupy, 
						const float velAtten, const float tempAtten, const float densityAtten, const int width, const int height, const int depth)
	{
		float time = 0;
		cudaEvent_t start,stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start,0);
		
		float* mac_u_d = 0;
		float* mac_v_d = 0;
		float* mac_w_d = 0;
		float* temperature_d = 0;
		float* density_d = 0;
		bool* occupy_d = 0;
		
		checkCudaErrors( cudaMalloc((void**)&mac_u_d,sizeof(float)*(width+1)*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&mac_v_d,sizeof(float)*width*(height+1)*depth) );
		checkCudaErrors( cudaMalloc((void**)&mac_w_d,sizeof(float)*width*height*(depth+1)) );
		checkCudaErrors( cudaMalloc((void**)&temperature_d,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&density_d,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&occupy_d,sizeof(bool)*width*height*depth) );
		checkCudaErrors( cudaMemcpy(mac_u_d,mac_u,sizeof(float)*(width+1)*height*depth,cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(mac_v_d,mac_v,sizeof(float)*width*(height+1)*depth,cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(mac_w_d,mac_w,sizeof(float)*width*height*(depth+1),cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(temperature_d,temperature,sizeof(float)*width*height*depth,cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(density_d,density,sizeof(float)*width*height*depth,cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(occupy_d,occupy,sizeof(bool)*width*height*depth,cudaMemcpyHostToDevice) );
		
		cu_Attenuation3D(mac_u_d,mac_v_d,mac_w_d,temperature_d,density_d,occupy_d,velAtten,tempAtten,densityAtten,width,height,depth);
		
		checkCudaErrors( cudaMemcpy(mac_u,mac_u_d,sizeof(float)*(width+1)*height*depth,cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(mac_v,mac_v_d,sizeof(float)*width*(height+1)*depth,cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(mac_w,mac_w_d,sizeof(float)*width*height*(depth+1),cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(density,density_d,sizeof(float)*width*height*depth,cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(temperature,temperature_d,sizeof(float)*width*height*depth,cudaMemcpyDeviceToHost) );
		
		checkCudaErrors( cudaFree(mac_u_d) );
		checkCudaErrors( cudaFree(mac_v_d) );
		checkCudaErrors( cudaFree(mac_w_d) );
		checkCudaErrors( cudaFree(temperature_d) );
		checkCudaErrors( cudaFree(density_d) );
		checkCudaErrors( cudaFree(occupy_d) );
		mac_u_d = 0;
		mac_v_d = 0;
		mac_w_d = 0;
		temperature_d = 0;
		density_d = 0;
		occupy_d = 0;
		
		cudaEventRecord(stop,0);
		cudaEventSynchronize(start);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time,start,stop);
		return time;
	}
}

#endif