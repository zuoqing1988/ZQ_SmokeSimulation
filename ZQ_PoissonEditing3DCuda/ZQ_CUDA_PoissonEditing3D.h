#ifndef _ZQ_CUDA_POISSON_EDITING_3D_H_
#define _ZQ_CUDA_POISSON_EDITING_3D_H_
#pragma once

extern "C"
float ZQ_CUDA_PoissonEditing3D(const bool* mask, const float* laplace, float* output, const int width, const int height, const int depth, const int nChannels, const int nIteration);

#endif