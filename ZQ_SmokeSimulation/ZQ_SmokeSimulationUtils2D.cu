#ifndef _ZQ_SMOKE_SIMULATION_UTILS_2D_CU_
#define _ZQ_SMOKE_SIMULATION_UTILS_2D_CU_

#include "ZQlibCudaDefines.cuh"

namespace ZQ_CUDA_MatrixMul
{
	extern "C" void MatrixMul(const float* A, const float* B, const int hA, const int wA, const int wB, float* C);
}

namespace ZQ_SmokeUtils2D
{
	extern "C" void SubProjectionCuda(const int subRow, const int subCol, const int fast_count, const float* subMatrix, const float* input, float* output)
	{
		ZQ_CUDA_MatrixMul::MatrixMul(subMatrix, input, subRow, subCol, fast_count, output);
	}
}

#endif