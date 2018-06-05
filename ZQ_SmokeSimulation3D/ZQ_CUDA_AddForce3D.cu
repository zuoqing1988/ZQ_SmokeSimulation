#ifndef _ZQ_CUDA_ADD_FORCE_3D_CU_
#define _ZQ_CUDA_ADD_FORCE_3D_CU_

#include "ZQlibCudaDefines.cuh"
#include "ZQ_CUDA_AddForce3D.cuh"
#include "ZQ_CUDA_AddForce3D.h"
#include "ZQ_CUDA_PoissonSolver3D.cuh"

namespace ZQ_CUDA_AddForce3D
{
	__global__
	void Compute_Vorticity_Vector_Scale_Kernel(float* vortVector, float* vortScale, const float* u, const float* v, const float* w, const float deltah, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width-1 || x <= 0 || y >= height-1 || y <= 0)
			return ;
		
		for(int z = 1; z < depth-1;z++)
		{
				int offset = z*height*width+y*width+x;
				vortVector[offset*3+0] = (w[offset+width] - w[offset-width] - v[offset+height*width] + v[offset-height*width])/(2*deltah);
				vortVector[offset*3+1] = (u[offset+height*width] - u[offset-height*width] - w[offset+1] + w[offset-1])/(2*deltah);
				vortVector[offset*3+2] = (v[offset+1] - v[offset-1] - u[offset+width] + u[offset-width])/(2*deltah);
				vortScale[offset] = sqrt(vortVector[offset*3+0]*vortVector[offset*3+0]
										+vortVector[offset*3+1]*vortVector[offset*3+1]
										+vortVector[offset*3+2]*vortVector[offset*3+2]);
		}	
	}
	
	
	__global__
	void Compute_Vorticity_Vector_Scale_Kernel(const int z, float* vortVector, float* vortScale, const float* u, const float* v, const float* w, const float deltah, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		x = x + 1;
		y = y + 1;
		if(x >= width-1|| y >= height-1) // x in [1,width-2],y in [1,height-2]
			return ;
		
		//z in [1,deoth-2] 
		
		int offset = z*height*width+y*width+x;
		vortVector[offset*3+0] = (w[offset+width] - w[offset-width] - v[offset+height*width] + v[offset-height*width])/(2*deltah);
		vortVector[offset*3+1] = (u[offset+height*width] - u[offset-height*width] - w[offset+1] + w[offset-1])/(2*deltah);
		vortVector[offset*3+2] = (v[offset+1] - v[offset-1] - u[offset+width] + u[offset-width])/(2*deltah);
		vortScale[offset] = sqrt(vortVector[offset*3+0]*vortVector[offset*3+0]
								+vortVector[offset*3+1]*vortVector[offset*3+1]
								+vortVector[offset*3+2]*vortVector[offset*3+2]);
		
	}
	
	__global__
	void Compute_Gradient_of_VorticityScale_Kernel(float* gradVort, float* vortScale, const int width, const int height, const int depth)
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
			int offset = z*height*width+y*width+x;
			if(x == 0)
				gradVort[offset*3+0] = vortScale[offset+1] - vortScale[offset];
			else if(x == width-1)
				gradVort[offset*3+0] = vortScale[offset] - vortScale[offset-1];
			else
				gradVort[offset*3+0] = 0.5f*(vortScale[offset+1] - vortScale[offset-1]);

			if(y == 0)
				gradVort[offset*3+1] = vortScale[offset+width] - vortScale[offset];
			else if(y == height-1)
				gradVort[offset*3+1] = vortScale[offset] - vortScale[offset-width];
			else
				gradVort[offset*3+1] = 0.5f*(vortScale[offset+width] - vortScale[offset-width]);

			if(z == 0)
				gradVort[offset*3+2] = vortScale[offset+height*width] - vortScale[offset];
			else if(z == depth-1)
				gradVort[offset*3+2] = vortScale[offset] - vortScale[offset-height*width];
			else
				gradVort[offset*3+2] = 0.5f*(vortScale[offset+height*width] - vortScale[offset-height*width]);

			float len = sqrt(gradVort[offset*3+0]*gradVort[offset*3+0]
							+ gradVort[offset*3+1]*gradVort[offset*3+1]
							+ gradVort[offset*3+2]*gradVort[offset*3+2]);
			if(len != 0)
			{
				gradVort[offset*3+0] /= len;
				gradVort[offset*3+1] /= len;
				gradVort[offset*3+2] /= len;
			}
		}
	}
	
	__global__
	void Compute_Gradient_of_VorticityScale_Kernel(const int z, float* gradVort, float* vortScale, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width || y >= height)
			return ;
		
		// z in [0,depth-1]
		
		int offset = z*height*width+y*width+x;
		if(x == 0)
			gradVort[offset*3+0] = vortScale[offset+1] - vortScale[offset];
		else if(x == width-1)
			gradVort[offset*3+0] = vortScale[offset] - vortScale[offset-1];
		else
			gradVort[offset*3+0] = 0.5f*(vortScale[offset+1] - vortScale[offset-1]);

		if(y == 0)
			gradVort[offset*3+1] = vortScale[offset+width] - vortScale[offset];
		else if(y == height-1)
			gradVort[offset*3+1] = vortScale[offset] - vortScale[offset-width];
		else
			gradVort[offset*3+1] = 0.5f*(vortScale[offset+width] - vortScale[offset-width]);

		if(z == 0)
			gradVort[offset*3+2] = vortScale[offset+height*width] - vortScale[offset];
		else if(z == depth-1)
			gradVort[offset*3+2] = vortScale[offset] - vortScale[offset-height*width];
		else
			gradVort[offset*3+2] = 0.5f*(vortScale[offset+height*width] - vortScale[offset-height*width]);

		float len = sqrt(gradVort[offset*3+0]*gradVort[offset*3+0]
						+ gradVort[offset*3+1]*gradVort[offset*3+1]
						+ gradVort[offset*3+2]*gradVort[offset*3+2]);
		if(len != 0)
		{
			gradVort[offset*3+0] /= len;
			gradVort[offset*3+1] /= len;
			gradVort[offset*3+2] /= len;
		}
	}
	
	__global__
	void Compute_Force_Kernel(float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, 
					const float buoyCoeff, const float deltah, const float Tamb, const int width, const int height, const int depth)
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
			int offset = z*height*width+y*width+x;	
				
			force[offset*3+0] = confineCoeff*deltah*(gradVort[offset*3+1]*vortVector[offset*3+2]-gradVort[offset*3+2]*vortVector[offset*3+1]);
			force[offset*3+1] = confineCoeff*deltah*(gradVort[offset*3+2]*vortVector[offset*3+0]-gradVort[offset*3+0]*vortVector[offset*3+2]);
			force[offset*3+2] = confineCoeff*deltah*(gradVort[offset*3+0]*vortVector[offset*3+1]-gradVort[offset*3+1]*vortVector[offset*3+0]);
			force[offset*3+1] += buoyCoeff*(temperature[offset]-Tamb);
		}
	}	
	
	__global__
	void Compute_Force_Kernel(const int z, float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, 
					const float buoyCoeff, const float deltah, const float Tamb, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width || y >= height)
			return ;
		
		//z in [0,depth-1]	
		
		int offset = z*height*width+y*width+x;	
				
		force[offset*3+0] = confineCoeff*deltah*(gradVort[offset*3+1]*vortVector[offset*3+2]-gradVort[offset*3+2]*vortVector[offset*3+1]);
		force[offset*3+1] = confineCoeff*deltah*(gradVort[offset*3+2]*vortVector[offset*3+0]-gradVort[offset*3+0]*vortVector[offset*3+2]);
		force[offset*3+2] = confineCoeff*deltah*(gradVort[offset*3+0]*vortVector[offset*3+1]-gradVort[offset*3+1]*vortVector[offset*3+0]);
		force[offset*3+1] += buoyCoeff*(temperature[offset]-Tamb);
		
	}
	
	__global__
	void AddForce_u_Kernel(float* mac_u, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
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
				if(!occupy[z*height*width+y*width+x])
				{
					mac_u[z*height*(width+1)+y*(width+1)+x] += deltat*force[(z*height*width+y*width+x)*3];
				}
			}
			else if(x == width)
			{
				if(!occupy[z*height*width+y*width+x-1])
				{
					mac_u[z*height*(width+1)+y*(width+1)+x] += deltat*force[(z*height*width+y*width+x-1)*3];
				}
			}
			else
			{
				if(!occupy[z*height*width+y*width+x-1] && !occupy[z*height*width+y*width+x])
					mac_u[z*height*(width+1)+y*(width+1)+x] += deltat*0.5f*(force[(z*height*width+y*width+x)*3]+force[(z*height*width+y*width+x-1)*3]);
			}
		}
	}
	
	__global__
	void AddForce_u_Kernel(const int z, float* mac_u, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x > width || y >= height)
			return ;
		
		//z in [0,depth-1]
		
		if(x == 0)
		{
			if(!occupy[z*height*width+y*width+x])
			{
				mac_u[z*height*(width+1)+y*(width+1)+x] += deltat*force[(z*height*width+y*width+x)*3];
			}
		}
		else if(x == width)
		{
			if(!occupy[z*height*width+y*width+x-1])
			{
				mac_u[z*height*(width+1)+y*(width+1)+x] += deltat*force[(z*height*width+y*width+x-1)*3];
			}
		}
		else
		{
			if(!occupy[z*height*width+y*width+x-1] && !occupy[z*height*width+y*width+x])
				mac_u[z*height*(width+1)+y*(width+1)+x] += deltat*0.5f*(force[(z*height*width+y*width+x)*3]+force[(z*height*width+y*width+x-1)*3]);
		}
	}
	
	__global__
	void AddForce_v_Kernel(float* mac_v, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
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
				if(!occupy[z*height*width+y*width+x])
				{
					mac_v[z*(height+1)*width+y*width+x] += deltat*force[(z*height*width+y*width+x)*3+1];
				}
			}
			else if(y == height)
			{
				if(!occupy[z*height*width+(y-1)*width+x])
				{
					mac_v[z*(height+1)*width+y*width+x] += deltat*force[(z*height*width+(y-1)*width+x)*3+1];
				}
			}
			else
			{
				if(!occupy[z*height*width+(y-1)*width+x] && !occupy[z*height*width+y*width+x])
					mac_v[z*(height+1)*width+y*width+x] += deltat*0.5f*(force[(z*height*width+y*width+x)*3+1]+force[(z*height*width+(y-1)*width+x)*3+1]);
			}
		}
	}
	
	__global__
	void AddForce_v_Kernel(const int z, float* mac_v, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width || y > height)
			return ;
		
		//z in [0,depth-1]	
		if(y == 0)
		{
			if(!occupy[z*height*width+y*width+x])
			{
				mac_v[z*(height+1)*width+y*width+x] += deltat*force[(z*height*width+y*width+x)*3+1];
			}
		}
		else if(y == height)
		{
			if(!occupy[z*height*width+(y-1)*width+x])
			{
				mac_v[z*(height+1)*width+y*width+x] += deltat*force[(z*height*width+(y-1)*width+x)*3+1];
			}
		}
		else
		{
			if(!occupy[z*height*width+(y-1)*width+x] && !occupy[z*height*width+y*width+x])
				mac_v[z*(height+1)*width+y*width+x] += deltat*0.5f*(force[(z*height*width+y*width+x)*3+1]+force[(z*height*width+(y-1)*width+x)*3+1]);
		}
	}
	
	__global__
	void AddForce_w_Kernel(float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
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
		{
			mac_w[y*width+x] += deltat*force[(y*width+x)*3+2];
		}
		
		if(!occupy[(depth-1)*height*width+y*width+x])
		{
			mac_w[depth*height*width+y*width+x] += deltat*force[((depth-1)*height*width+y*width+x)*3+2];
		}
			
		for(int z = 0;z < depth;z++)
		{
			if(!occupy[(z-1)*height*width+y*width+x] && !occupy[z*height*width+y*width+x])
				mac_w[z*height*width+y*width+x] += deltat*0.5f*(force[(z*height*width+y*width+x)*3+2]+force[((z-1)*height*width+y*width+x)*3+2]);
		}
	}
	
	__global__
	void AddForce_w_Kernel(const int z, float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
	{
		int bx = blockIdx.x;
		int by = blockIdx.y;
		int tx = threadIdx.x;
		int ty = threadIdx.y;

		int x = bx*blockDim.x+tx;
		int y = by*blockDim.y+ty;
		if(x >= width || y >= height)
			return ;
			
		if(z == 0)
		{
			if(!occupy[y*width+x])
			{
				mac_w[y*width+x] += deltat*force[(y*width+x)*3+2];

			}
		}
		else if(z == depth)
		{
			if(!occupy[(depth-1)*height*width+y*width+x])
			{
				mac_w[depth*height*width+y*width+x] += deltat*force[((depth-1)*height*width+y*width+x)*3+2];
			}
		}
		else
		{	
			if(!occupy[(z-1)*height*width+y*width+x] && !occupy[z*height*width+y*width+x])
				mac_w[z*height*width+y*width+x] += deltat*0.5f*(force[(z*height*width+y*width+x)*3+2]+force[((z-1)*height*width+y*width+x)*3+2]);
		}
	}
	
	
	
	/*****************************************************/

	void cu_Compute_Vorticity_Vector_Scale(float* vortVector, float* vortScale, float* u, float* v, float* w, const float deltah, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width+blockSize.x-1)/blockSize.x, (height+blockSize.y-1)/blockSize.y);
		
		Compute_Vorticity_Vector_Scale_Kernel<<<gridSize,blockSize>>>(vortVector,vortScale,u,v,w,deltah,width,height,depth);
	}
	
	void cu_Compute_Vorticity_Vector_Scale2(float* vortVector, float* vortScale, float* u, float* v, float* w, const float deltah, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width-2+blockSize.x-1)/blockSize.x, (height-2+blockSize.y-1)/blockSize.y);
		
		for(int z = 1;z < depth-1;z++)
			Compute_Vorticity_Vector_Scale_Kernel<<<gridSize,blockSize>>>(z,vortVector,vortScale,u,v,w,deltah,width,height,depth);
	}
	
	void cu_Compute_Gradient_of_VorticityScale(float* gradVort, float* vortScale, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		
		Compute_Gradient_of_VorticityScale_Kernel<<<gridSize,blockSize>>>(gradVort, vortScale, width, height, depth);
	}
	
	void cu_Compute_Gradient_of_VorticityScale2(float* gradVort, float* vortScale, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		
		for(int z = 0;z < depth;z++)
			Compute_Gradient_of_VorticityScale_Kernel<<<gridSize,blockSize>>>(z, gradVort, vortScale, width, height, depth);
	}
	
	void cu_Compute_Force(float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, const float buoyCoeff, 
			const float deltah, const float Tamb, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		
		Compute_Force_Kernel<<<gridSize,blockSize>>>(force,gradVort,vortVector,temperature,confineCoeff,buoyCoeff,deltah,Tamb,width,height,depth);
	}
	
	void cu_Compute_Force2(float* force, const float* gradVort, const float* vortVector, const float* temperature, const float confineCoeff, const float buoyCoeff, 
			const float deltah, const float Tamb, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		
		for(int z = 0;z < depth;z++)
			Compute_Force_Kernel<<<gridSize,blockSize>>>(z,force,gradVort,vortVector,temperature,confineCoeff,buoyCoeff,deltah,Tamb,width,height,depth);
	}
	
	void cu_AddForce_u_v_w(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 u_gridSize((width+1+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		dim3 v_gridSize((width+blockSize.x-1)/blockSize.x,(height+1+blockSize.y-1)/blockSize.y);
		dim3 w_gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		
		AddForce_u_Kernel<<<u_gridSize,blockSize>>>(mac_u,occupy,force,deltat,width,height,depth);
		AddForce_v_Kernel<<<v_gridSize,blockSize>>>(mac_v,occupy,force,deltat,width,height,depth);
		AddForce_w_Kernel<<<w_gridSize,blockSize>>>(mac_w,occupy,force,deltat,width,height,depth);
	}
	
	void cu_AddForce_u_v_w2(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* force, const float deltat, const int width, const int height, const int depth)
	{
		dim3 blockSize(BLOCK_SIZE,BLOCK_SIZE);
		dim3 u_gridSize((width+1+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		dim3 v_gridSize((width+blockSize.x-1)/blockSize.x,(height+1+blockSize.y-1)/blockSize.y);
		dim3 w_gridSize((width+blockSize.x-1)/blockSize.x,(height+blockSize.y-1)/blockSize.y);
		
		for(int z = 0;z < depth;z++)
			AddForce_u_Kernel<<<u_gridSize,blockSize>>>(z,mac_u,occupy,force,deltat,width,height,depth);
		for(int z = 0;z < depth;z++)
			AddForce_v_Kernel<<<v_gridSize,blockSize>>>(z,mac_v,occupy,force,deltat,width,height,depth);
		
		for(int z = 0;z <= depth;z++)
			AddForce_w_Kernel<<<w_gridSize,blockSize>>>(z,mac_w,occupy,force,deltat,width,height,depth);
	}

	/****************************************************************/

	 
	void cu_AddForce3D(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* temperature, const float deltah, const float deltat, 
						const float buoyCoeff, const float confineCoeff, const float Tamb, const int width, const int height, const int depth)
	{	
		float* u = 0;
		float* v = 0;
		float* w = 0;
		checkCudaErrors( cudaMalloc((void**)&u,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&v,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&w,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(u,0,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(v,0,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(w,0,sizeof(float)*width*height*depth) );
		
		ZQ_CUDA_PoissonSolver3D::cu_MAC_to_Regular_vel(u,v,w,mac_u,mac_v,mac_w,width,height,depth);
		
		float* vortVector = 0;
		float* vortScale = 0;
		checkCudaErrors( cudaMalloc((void**)&vortVector,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMalloc((void**)&vortScale,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(vortVector,0,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMemset(vortScale,0,sizeof(float)*width*height*depth) );
		
		cu_Compute_Vorticity_Vector_Scale(vortVector,vortScale,u,v,w,deltah, width,height,depth);
		
		float* gradVort = 0;
		checkCudaErrors( cudaMalloc((void**)&gradVort,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMemset(gradVort,0,sizeof(float)*width*height*depth*3) );
		
		cu_Compute_Gradient_of_VorticityScale(gradVort,vortScale,width,height,depth);	
		
		float* force = 0;
		checkCudaErrors( cudaMalloc((void**)&force,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMemset(force,0,sizeof(float)*width*height*depth*3) );
		
		cu_Compute_Force(force,gradVort,vortVector, temperature, confineCoeff, buoyCoeff, deltah, Tamb, width, height, depth);
		
		cu_AddForce_u_v_w(mac_u,mac_v,mac_w,occupy,force,deltat,width,height,depth);
	
		checkCudaErrors( cudaFree(u) );
		checkCudaErrors( cudaFree(v) );
		checkCudaErrors( cudaFree(w) );
		checkCudaErrors( cudaFree(vortVector) );
		checkCudaErrors( cudaFree(vortScale) );
		checkCudaErrors( cudaFree(gradVort) );
		checkCudaErrors( cudaFree(force) );
		
	}

	void cu_AddForce3D2(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* temperature, const float deltah, const float deltat, 
						const float buoyCoeff, const float confineCoeff, const float Tamb, const int width, const int height, const int depth)
	{	
		float* u = 0;
		float* v = 0;
		float* w = 0;
		checkCudaErrors( cudaMalloc((void**)&u,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&v,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&w,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(u,0,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(v,0,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(w,0,sizeof(float)*width*height*depth) );
		
		ZQ_CUDA_PoissonSolver3D::cu_MAC_to_Regular_vel(u,v,w,mac_u,mac_v,mac_w,width,height,depth);
		
		float* vortVector = 0;
		float* vortScale = 0;
		checkCudaErrors( cudaMalloc((void**)&vortVector,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMalloc((void**)&vortScale,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemset(vortVector,0,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMemset(vortScale,0,sizeof(float)*width*height*depth) );
		
		cu_Compute_Vorticity_Vector_Scale2(vortVector,vortScale,u,v,w,deltah, width,height,depth);
		
		float* gradVort = 0;
		checkCudaErrors( cudaMalloc((void**)&gradVort,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMemset(gradVort,0,sizeof(float)*width*height*depth*3) );
		
		cu_Compute_Gradient_of_VorticityScale2(gradVort,vortScale,width,height,depth);	
		
		float* force = 0;
		checkCudaErrors( cudaMalloc((void**)&force,sizeof(float)*width*height*depth*3) );
		checkCudaErrors( cudaMemset(force,0,sizeof(float)*width*height*depth*3) );
		
		cu_Compute_Force2(force,gradVort,vortVector, temperature, confineCoeff, buoyCoeff, deltah, Tamb, width, height, depth);
		
		cu_AddForce_u_v_w2(mac_u,mac_v,mac_w,occupy,force,deltat,width,height,depth);
	
		checkCudaErrors( cudaFree(u) );
		checkCudaErrors( cudaFree(v) );
		checkCudaErrors( cudaFree(w) );
		checkCudaErrors( cudaFree(vortVector) );
		checkCudaErrors( cudaFree(vortScale) );
		checkCudaErrors( cudaFree(gradVort) );
		checkCudaErrors( cudaFree(force) );
	}

	/**********************************************************/
	extern "C" 
	float AddForce3D(float* mac_u, float* mac_v, float* mac_w, const bool* occupy, const float* temperature, const float deltah, const float deltat, 
						const float buoyCoeff, const float confineCoeff, const float Tamb, const int width, const int height, const int depth, enum AddForceType type)
	{
		float time = 0;
		cudaEvent_t start,stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start,0);
		
		float* mac_u_d = 0;
		float* mac_v_d = 0;
		float* mac_w_d = 0;
		bool* occupy_d = 0;
		float* temperature_d = 0;
		
		checkCudaErrors( cudaMalloc((void**)&mac_u_d,sizeof(float)*(width+1)*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&mac_v_d,sizeof(float)*width*(height+1)*depth) );
		checkCudaErrors( cudaMalloc((void**)&mac_w_d,sizeof(float)*width*height*(depth+1)) );
		checkCudaErrors( cudaMalloc((void**)&occupy_d,sizeof(bool)*width*height*depth) );
		checkCudaErrors( cudaMalloc((void**)&temperature_d,sizeof(float)*width*height*depth) );
		checkCudaErrors( cudaMemcpy(mac_u_d,mac_u,sizeof(float)*(width+1)*height*depth,cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(mac_v_d,mac_v,sizeof(float)*width*(height+1)*depth,cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(mac_w_d,mac_w,sizeof(float)*width*height*(depth+1),cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(occupy_d,occupy,sizeof(bool)*width*height*depth,cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy(temperature_d,temperature,sizeof(float)*width*height*depth,cudaMemcpyHostToDevice) );
		
		switch(type)
		{
		case AddForceType::ADD_FORCE_ENTIRE:
			cu_AddForce3D(mac_u_d, mac_v_d, mac_w_d, occupy_d, temperature_d, deltah, deltat, buoyCoeff, confineCoeff, Tamb, width, height, depth);
			break;
		case AddForceType::ADD_FORCE_SLICE:
			cu_AddForce3D(mac_u_d, mac_v_d, mac_w_d, occupy_d, temperature_d, deltah, deltat, buoyCoeff, confineCoeff, Tamb, width, height, depth);
			break;
		}
	
		checkCudaErrors( cudaMemcpy(mac_u,mac_u_d,sizeof(float)*(width+1)*height*depth,cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(mac_v,mac_v_d,sizeof(float)*width*(height+1)*depth,cudaMemcpyDeviceToHost) );
		checkCudaErrors( cudaMemcpy(mac_w,mac_w_d,sizeof(float)*width*height*(depth+1),cudaMemcpyDeviceToHost) );
		
		checkCudaErrors( cudaFree(mac_u_d) );
		checkCudaErrors( cudaFree(mac_v_d) );
		checkCudaErrors( cudaFree(mac_w_d) );
		checkCudaErrors( cudaFree(occupy_d) );
		checkCudaErrors( cudaFree(temperature_d) );
		
		cudaEventRecord(stop,0);
		cudaEventSynchronize(start);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time,start,stop);
		return time;
		
	}
}

#endif