#include "ZQ_CUDA_BlendTwoImages.h"
#include "ZQ_ImageIO.h"
#include <time.h>

using namespace ZQ;
typedef float BaseType;
typedef ZQ_DImage<BaseType> DImage;

int main(const int argc, const char** argv)
{
	clock_t t1 = clock();

	if (argc != 6)
	{
		printf(".exe imagefile1 imagefile2 flowfile weight1 outfile\n");
		return 0;
	}
	const char* imagefile1 = argv[1];
	const char* imagefile2 = argv[2];
	const char* flowfile = argv[3];
	float weight1 = atof(argv[4]);
	const char* outfile = argv[5];
	
	DImage img1, img2, flow, u, v, img3;
	if (!ZQ_ImageIO::loadImage(img1, imagefile1, 1))
	{
		printf("failed to load image %s\n", imagefile1);
		return 0;
	}
	if (!ZQ_ImageIO::loadImage(img2, imagefile2, 1))
	{
		printf("failed to load image %s\n", imagefile2);
		return 0;
	}
	if (!flow.loadImage(flowfile))
	{
		printf("failed to load flow %s\n", flowfile);
		return 0;
	}

	int width = img1.width();
	int height = img1.height();
	int nChannels = img1.nchannels();
	if (!img2.matchDimension(width, height, nChannels))
	{
		printf("dimension does not match\n");
		return 0;
	}

	if (!flow.matchDimension(width, height, 2))
	{
		printf("dimension does not match\n");
		return 0;
	}

	flow.separate(1, u, v);

	img3.allocate(width, height, nChannels);
	BaseType*& img1_data = img1.data();
	BaseType*& img2_data = img2.data();
	BaseType*& img3_data = img3.data();
	BaseType*& u_data = u.data();
	BaseType*& v_data = v.data();
	float cuda_cost_time = ZQ_CUDA_BlendTwoImages::Cutil_BlendTwoImagesByMedFilt(width, height, nChannels, img1_data, img2_data, u_data, v_data, weight1, img3_data, 0, 0);

	if (!ZQ_ImageIO::saveImage(img3, outfile))
	{
		printf("failed to save image %s\n", outfile);
		return 0;
	}

	clock_t t2 = clock();

	printf("cuda_cost_time / total_cost = %f/%f \n", cuda_cost_time*0.001, (t2 - t1)*0.001);
	return 0;
}