#include "ZQ_MathBase.h"
#include "ZQ_CUDA_MatrixMul.h"
#include <stdio.h>

const int h1 = 2, w1 = 3;
const int h2 = w1, w2 = 4;

int main()
{
	float cpu_A[h1*w1] = { 1, 2, 3, 4, 5, 6 };
	float cpu_B[h2*w2] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
	float cpu_C[h1*w2];
	float gpu_C[h1*w2];

	

	ZQ::ZQ_MathBase::MatrixMul(cpu_A, cpu_B, h1, w1, w2, cpu_C);
	ZQ_CUDA_MatrixMul::MatrixMul(cpu_A, cpu_B, h1, w1, w2, gpu_C);

	for (int i = 0; i < h1; i++)
	{
		for (int j = 0; j < w2; j++)
		{
			printf("%12.1f", cpu_C[i*w2 + j]);
		}
		printf("\n");
	}

	for (int i = 0; i < h1; i++)
	{
		for (int j = 0; j < w2; j++)
		{
			printf("%12.1f", gpu_C[i*w2 + j]);
		}
		printf("\n");
	}

	return 0;
}