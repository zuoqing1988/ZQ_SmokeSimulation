#include "ZQ_CUDA_Wavelet.h"
#include "ZQ_DoubleImage.h"
#include "ZQ_Wavelet.h"
#include "ZQ_ImageIO.h"

using namespace ZQ;
typedef ZQ_DImage<double> DImage;
void main()
{
	const char* filename = "1.jpg";
	IplImage* img = cvLoadImage(filename,0);
	if(img == 0)
	{
		printf("failed to load image %s\n",filename);
		return;
	}

	int width = img->width;
	int height = img->height;

	float* input = new float[width*height];
	for(int h = 0;h < height;h++)
	{
		for(int w = 0;w < width;w++)
			input[h*width+w] = cvGetReal2D(img,h,w);
	}

	float* output = 0;
	int out_width,out_height;
	ZQ_CUDA_Wavelet::DWT2(input,width,height,5,output,out_width,out_height);

	DImage image,cpu_out_image;
	ZQ_ImageIO::loadImage(image,filename,0);

	ZQ_Wavelet<double> m_wave;
	m_wave.DiscreteWaveletImageNLevels(image,"db1",5,ZQ_Wavelet<double>::PADDING_ZPD);
	m_wave.GetWaveletImage(cpu_out_image);



	IplImage* out_img = cvCreateImage(cvSize(out_width,out_height),IPL_DEPTH_8U,1);
	IplImage* cpu_out_img = cvCreateImage(cvSize(out_width,out_height),IPL_DEPTH_8U,1);
	for(int h = 0; h < out_height;h++)
	{
		for(int w = 0;w < out_width;w++)
		{
			cvSetReal2D(out_img,h,w,fabs(output[h*out_width+w]));
			cvSetReal2D(cpu_out_img,h,w,fabs(cpu_out_image.data()[h*out_width+w]*255));
		}
	}

	cvNamedWindow("ori");
	cvNamedWindow("gpu_dwt");
	cvNamedWindow("cpu_dwt");
	cvShowImage("ori",img);
	cvShowImage("gpu_dwt",out_img);
	cvShowImage("cpu_dwt",cpu_out_img);
	cvWaitKey(0);
	cvReleaseImage(&img);
	cvReleaseImage(&out_img);
	cvReleaseImage(&cpu_out_img);
	delete []input;
	delete []output;

}