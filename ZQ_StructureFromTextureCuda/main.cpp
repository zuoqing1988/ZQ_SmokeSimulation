#include "ZQ_ImageIO.h"
#include "ZQ_CUDA_StructureFromTexture.h"
#include <time.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <vector_types.h>
#include <helper_math.h>

using namespace ZQ;
int exp1(int argc, const char** argv);
int exp2(int argc, const char** argv);
int exp3(int argc, const char** argv);
int exp4(int argc, const char** argv);
int exp5(int argc, const char** argv);
int exp6(int argc, const char** argv);
int exp7(int argc, const char** argv);
int exp8(int argc, const char** argv);
int exp9(int argc, const char** argv);
int for_nvprof(int argc, const char** argv);

int main(int argc, const char** argv)
{
	if (_strcmpi(argv[1], "exp1") == 0)
		return exp1(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp2") == 0)
		return exp2(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp3") == 0)
		return exp3(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp4") == 0)
		return exp4(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp5") == 0)
		return exp5(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp6") == 0)
		return exp6(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp7") == 0)
		return exp7(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp8") == 0)
		return exp8(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "exp9") == 0)
		return exp9(argc - 1, argv + 1);
	else if (_strcmpi(argv[1], "for_nvprof") == 0)
		return for_nvprof(argc - 1, argv + 1);
	return EXIT_FAILURE;
}

int exp1(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();
	output = input;

	float lambda = 0.02;
	//int nOuterIter = 4;
	int nSORIter = 200;
	float eps = 1e-3;
	//int fsize = 3;
	//float sigma = 5;
	
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);

	sprintf(cur_fold, "%s\\out-0.5-1-2-1", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int nOutIter = 0; nOutIter <= 20; nOutIter++)
	{
		ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda, nOutIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2, 1, eps);

		sprintf(name, "%s\\iteration-%02d.png", cur_fold, nOutIter);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}

	}
		

	sprintf(cur_fold, "%s\\out-0.5-1-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int nOutIter = 0; nOutIter <= 20; nOutIter++)
	{
		ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda, nOutIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2, 0.5, eps);

		sprintf(name, "%s\\iteration-%02d.png", cur_fold, nOutIter);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}

	}

	sprintf(cur_fold, "%s\\out-1.0-1-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int nOutIter = 0; nOutIter <= 20; nOutIter++)
	{
		ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda, nOutIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1.0, 1, 2, 0.5, eps);

		sprintf(name, "%s\\iteration-%02d.png", cur_fold, nOutIter);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}

	}

	sprintf(cur_fold, "%s\\out-1.0-1-1.5-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int nOutIter = 0; nOutIter <= 20; nOutIter++)
	{
		ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*3, nOutIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1.0, 1, 1.5, 0.5, eps);

		sprintf(name, "%s\\iteration-%02d.png", cur_fold, nOutIter);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}

	}
	
	return EXIT_SUCCESS;

}

int exp2(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 2)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	float lambda = 0.001;
	int nOuterIter = 3;
	int nSORIter = 20;
	float eps = 1e-3;
	int fsize = 3;
	float sigma = 5;

	output = input;
	if (!ZQ_ImageIO::saveImage(output, "orginal.png"))
	{
		return EXIT_FAILURE;
	}
	ZQ_DImage<float> tmp = input;
	float cost_time = 0;
	cost_time += ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width,height,nChannels, 
		lambda, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1, eps);
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl1-type1-improve.png"))
	{
		return EXIT_FAILURE;
	}

	cost_time = 0;
	cost_time+= ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		lambda * 3, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2, 1, eps);
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl0.5-type1-improve.png"))
	{
		return EXIT_FAILURE;
	}

	cost_time = 0;
	cost_time+= ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		lambda * 5, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.8, 1, 2, 1, eps);
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl0.8-type1-improve.png"))
	{
		return EXIT_FAILURE;
	}

	tmp = input;
	cost_time = 0;
	for (int it = 0; it < nOuterIter; it++)
	{
		cost_time += ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), tmp.data(), width, height, nChannels,
			lambda, 1, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1, eps);
		tmp = output;
	}
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl1-type1-improve10.png"))
	{
		return EXIT_FAILURE;
	}

	cost_time = 0;
	cost_time+=ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		lambda * 20, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 1, eps);
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl1-type2-improve.png"))
	{
		return EXIT_FAILURE;
	}

	tmp = input;
	cost_time = 0;
	for (int it = 0; it < nOuterIter; it++)
	{
		cost_time += ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), tmp.data(),  width, height, nChannels,
			lambda, 1, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 1, eps);
		tmp = output;
	}
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl1-type2-improve10.png"))
	{
		return EXIT_FAILURE;
	}

	cost_time = 0;
	cost_time += ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		lambda, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 0, eps);
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl1-type3-improve.png"))
	{
		return EXIT_FAILURE;
	}

	tmp = input;
	cost_time = 0;
	for (int it = 0; it < nOuterIter; it++)
	{
		cost_time += ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), tmp.data(), width, height, nChannels, 
			lambda*0.1, 1, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 0, eps);
		tmp = output;
	}
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-rtvl1-type3-improve10.png"))
	{
		return EXIT_FAILURE;
	}

	cost_time = 0;
	cost_time += ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels, 
		lambda * 30, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 1, eps);
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-wls-type2-improve.png"))
	{
		return EXIT_FAILURE;
	}

	tmp = input;
	cost_time = 0;
	for (int it = 0; it < nOuterIter; it++)
	{
		cost_time += ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), tmp.data(), width, height, nChannels, 
			lambda * 3, 1, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 1, eps);
		tmp = output;
	}
	printf("%f\n", cost_time);
	if (!ZQ_ImageIO::saveImage(output, "out-wls-type2-improve10.png"))
	{
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}


int exp3(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	float lambda = 0.001;
	int nOuterIter = 3;
	int nSORIter = 20;
	float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;
	if (!ZQ_ImageIO::saveImage(output, "orginal.png"))
	{
		return EXIT_FAILURE;
	}
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);
	
	sprintf(cur_fold, "%s\\out-0.0-1-1-1", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 1, 1, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0.0-1-1.5-1", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 1.5, 1, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}
	
	sprintf(cur_fold, "%s\\out-0.0-1-2-1-RTV", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i+=5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1, eps);
		
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}
	
	sprintf(cur_fold, "%s\\out-0.5-1-2-1", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2, 1, eps);
		
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}
	
	sprintf(cur_fold, "%s\\out-1-1-2-1", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 1, eps);
		
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-1-1-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 0.5, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}
	

	sprintf(cur_fold, "%s\\out-1-1-2-0", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 0, eps);
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-0-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 0.5, eps);
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-0-2-1-WLS", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 1, eps);
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-0-2-1.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 1.5, eps);
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-0-2-2", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 2, eps);
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}
	

	return EXIT_SUCCESS;
}


int exp4(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	float lambda = 0.001;
	int nOuterIter = 3;
	int nSORIter = 20;
	float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;
	if (!ZQ_ImageIO::saveImage(output, "orginal.png"))
	{
		return EXIT_FAILURE;
	}
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);

	sprintf(cur_fold, "%s\\out-0.0-1-2-1.0-RTV", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0.2-1-2-0.8", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.2, 1, 2, 0.8, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0.5-1-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2, 0.5, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0.8-1-2-0.2", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.8, 1, 2, 0.2, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-1.0-1-2-0.0", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 1, 1, 2, 0, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}


	return EXIT_SUCCESS;
}


int exp5(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	float lambda = 0.001;
	int nOuterIter = 3;
	int nSORIter = 20;
	float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;
	if (!ZQ_ImageIO::saveImage(output, "orginal.png"))
	{
		return EXIT_FAILURE;
	}
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);


	sprintf(cur_fold, "%s\\out-0-1-2-2.0", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 2.0, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-1.8", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1.8, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-1.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1.5, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-1.2", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1.2, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-1.0-RTV", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-0.8", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 0.8, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 0.5, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-0.2", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 0.2, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0-1-2-0.0", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 100; i += 5)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 0, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}


	return EXIT_SUCCESS;
}


int exp6(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	float lambda = 0.001;
	int nOuterIter = 3;
	int nSORIter = 20;
	float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;
	if (!ZQ_ImageIO::saveImage(output, "orginal.png"))
	{
		return EXIT_FAILURE;
	}
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);

	float p = 2, s = 1;
	for (float q = 0; q <= 2.0; q += 0.5)
	{
		for (float r = 0; r <= 2.0; r += 0.5)
		{
			sprintf(cur_fold, "%s\\out-%.1f-%.0f-%.0f-%.1f", out_fold, r, s, p, q);
			sprintf(cmdline, "mkdir %s", cur_fold);
			system(cmdline);
			for (int i = 0; i < 100; i += 5)
			{
				float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
					lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, r, s, p, q, eps);

				sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
				if (!ZQ_ImageIO::saveImage(output, name))
				{
					return EXIT_FAILURE;
				}
			}
		}
	}

	return EXIT_SUCCESS;
}


int exp7(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	float lambda = 0.01;
	int nOuterIter = 3;
	int nSORIter = 20;
	float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;
	if (!ZQ_ImageIO::saveImage(output, "orginal.png"))
	{
		return EXIT_FAILURE;
	}
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);

	sprintf(cur_fold, "%s\\out-0.0-1-2-1-RTV", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 10; i ++)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0.5-1-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (int i = 0; i < 10; i ++)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2, 0.5, eps);

		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	

	sprintf(cur_fold, "%s\\out-0-0-2-1-WLS", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int i = 0; i < 10; i ++)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 1, eps);
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}

	sprintf(cur_fold, "%s\\out-0.5-0-2-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (int i = 0; i < 10; i ++)
	{
		float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
			lambda*i, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 0, 2, 0.5, eps);
		sprintf(name, "%s\\%.3f.png", cur_fold, lambda*i);
		if (!ZQ_ImageIO::saveImage(output, name))
		{
			return EXIT_FAILURE;
		}
	}


	return EXIT_SUCCESS;
}




int exp8(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	float lambda = 0.02;
	int nOuterIter = 3;
	int nSORIter = 20;
	float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;
	
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);

	sprintf(cur_fold, "%s\\out-0.45-1.1-1.8-0.55", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		lambda, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.45, 1.1, 1.8, 0.55, eps);

	printf("time = %f\n", cost_time);
	sprintf(name, "%s\\%.3f.png", cur_fold, lambda);
	if (!ZQ_ImageIO::saveImage(output, name))
	{
		return EXIT_FAILURE;
	}
	sprintf(cur_fold, "%s\\out-0.0-1.0-2.0-1.0-RTV", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		lambda, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 1, 2, 1, eps);

	printf("time = %f\n", cost_time);
	sprintf(name, "%s\\%.3f.png", cur_fold, lambda);
	if (!ZQ_ImageIO::saveImage(output, name))
	{
		return EXIT_FAILURE;
	}
	sprintf(cur_fold, "%s\\out-0.0-0.0-2.0-1.0-WLS", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		lambda*30, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0, 0, 2, 0.8, eps*0.1);

	printf("time = %f\n", cost_time);
	sprintf(name, "%s\\%.3f.png", cur_fold, lambda);
	if (!ZQ_ImageIO::saveImage(output, name))
	{
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int exp9(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	//float lambda = 0.02;
	int nOuterIter = 3;
	int nSORIter = 20;
	//float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;

	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);

	sprintf(cur_fold, "%s\\out-0.5-1.0-2.0-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);

	for (float lambda = 0.01; lambda < 1; lambda += 0.01)
	{
		for (float eps = 1e-6; eps <= 1e-3; eps *= 10)
		{
			float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
				lambda, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2.0, 0.5, eps);

			printf("time = %f\n", cost_time);
			sprintf(name, "%s\\%.3f-%.6f.png", cur_fold, lambda, eps);
			if (!ZQ_ImageIO::saveImage(output, name))
			{
				return EXIT_FAILURE;
			}
		}
	}
	
	sprintf(cur_fold, "%s\\out-0.5-0.1-2.0-0.5", out_fold);
	sprintf(cmdline, "mkdir %s", cur_fold);
	system(cmdline);
	for (float lambda = 0.01; lambda < 1; lambda += 0.01)
	{
		for (float eps = 1e-6; eps <= 1e-3; eps *= 10)
		{
			float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
				lambda, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 0.1, 2.0, 0.5, eps);

			printf("time = %f\n", cost_time);
			sprintf(name, "%s\\%.3f-%.6f.png", cur_fold, lambda, eps);
			if (!ZQ_ImageIO::saveImage(output, name))
			{
				return EXIT_FAILURE;
			}
		}
	}
	
	
	return EXIT_SUCCESS;
}

int for_nvprof(int argc, const char** argv)
{
	ZQ_DImage<float> input, output;

	if (argc < 5)
		return EXIT_FAILURE;
	const char* inputfile = argv[1];
	const char* out_fold = argv[2];
	int fsize = atoi(argv[3]);
	float sigma = atof(argv[4]);
	if (!ZQ_ImageIO::loadImage(input, inputfile, 1))
	{
		printf("failed to load %s\n", inputfile);
		return EXIT_FAILURE;
	}

	int width = input.width();
	int height = input.height();
	int nChannels = input.nchannels();

	//float lambda = 0.02;
	int nOuterIter = 3;
	int nSORIter = 20;
	//float eps = 1e-3;
	//int fsize = 1;
	//float sigma = 2;

	output = input;
	printf("!\n");
	char name[200];
	char cur_fold[200];
	char cmdline[200];
	sprintf(cmdline, "mkdir %s", out_fold);
	system(cmdline);

	checkCudaErrors(cudaProfilerStart());
	printf("begin...\n");
	float cost_time = ZQ_CUDA_StructureFromTexture::StructureFromTextureImprovedWLS(output.data(), input.data(), width, height, nChannels,
		0.02, nOuterIter, nSORIter, fsize, sigma, fsize, sigma, fsize, sigma, 0.5, 1, 2.0, 0.5, 1e-3);
	
	printf("time = %f\n", cost_time);
	checkCudaErrors(cudaProfilerStop());
	sprintf(name, "%s\\out.png", out_fold);
	if (!ZQ_ImageIO::saveImage(output, name))
	{
		return EXIT_FAILURE;
	}

	

	return EXIT_SUCCESS;
}