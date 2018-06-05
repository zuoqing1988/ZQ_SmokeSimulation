#ifndef _ZQ_CUDA_POISSON_EDITING_3D_CU_
#define _ZQ_CUDA_POISSON_EDITING_3D_CU_


#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 16
#endif

__global__
void poisson_editing3d_RedBlack_Kernel(const bool* mask, const float* laplace, float* output, const int width, const int height, const int depth, const int nChannels, bool redkernel)
{
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	int x = bx*blockDim.x+tx;
	int y = by*blockDim.y+ty;
	if(x >= width || y >= height)
		return ;

	int start_z = redkernel ? 0 : 1;

	int SLICE = width*height;
	for(int z = start_z;z < depth;z += 2)
	{
		int offset_single = z*height*width+y*width+x;
		for(int c = 0;c < nChannels;c++)
		{
			int offset = offset_single*nChannels+c;
			if(mask[offset])
			{
				float coeff = 6;
				float sigma = output[offset+SLICE] + output[offset-SLICE]
				+ output[offset+width] + output[offset-width]
				+ output[offset+1] + output[offset-1]
				- laplace[offset];
				output[offset] = sigma/coeff;
			}
		}
	}
}

////////////////////////
void cu_PoissonEditing3D(const bool* mask, const float* laplace, float* output, const int width, const int height, const int depth, const int nChannels, const int nIteration)
{
	dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
	dim3 gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y)/blockSize.y);

	for(int i = 0;i < nIteration;i++)
	{
		poisson_editing3d_RedBlack_Kernel<<<gridSize,blockSize>>>(mask,laplace,output,width,height,depth,nChannels,true);
		poisson_editing3d_RedBlack_Kernel<<<gridSize,blockSize>>>(mask,laplace,output,width,height,depth,nChannels,false);
	}
}

extern "C"
float ZQ_CUDA_PoissonEditing3D(const bool* mask, const float* laplace, float* output, const int width, const int height, const int depth, const int nChannels, const int nIteration)
{
	float time = 0;
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);

	bool* mask_d = 0;
	float* laplace_d = 0;
	float* output_d = 0;
	checkCudaErrors( cudaMalloc((void**)&mask_d,sizeof(bool)*width*height*depth) );
	checkCudaErrors( cudaMalloc((void**)&laplace_d,sizeof(float)*width*height*depth*nChannels) );
	checkCudaErrors( cudaMalloc((void**)&output_d,sizeof(float)*width*height*depth*nChannels) );

	checkCudaErrors( cudaMemcpy(mask_d,mask,sizeof(bool)*width*height*depth,cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(laplace_d,laplace,sizeof(float)*width*height*depth*nChannels,cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy(output_d,output,sizeof(float)*width*height*depth*nChannels,cudaMemcpyHostToDevice) );

	cu_PoissonEditing3D(mask_d,laplace_d,output_d,width,height,depth,nChannels,nIteration);

	checkCudaErrors( cudaMemcpy(output,output_d,sizeof(float)*width*height*depth*nChannels,cudaMemcpyDeviceToHost) );

	checkCudaErrors( cudaFree(mask_d) );
	checkCudaErrors( cudaFree(laplace_d) );
	checkCudaErrors( cudaFree(output_d) );
	mask_d = 0;
	laplace_d = 0;
	output_d = 0;

	cudaEventRecord(stop,0);
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);
	return time;
}

#endif