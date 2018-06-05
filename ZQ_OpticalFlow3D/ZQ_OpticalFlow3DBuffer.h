#ifndef _ZQ_OPTICAL_FLOW3D_BUFFER_H_
#define _ZQ_OPTICAL_FLOW3D_BUFFER_H_

#include "ZQ_DoubleImage3D.h"
#include "ZQ_GaussianPyramid3DCuda.h"

namespace ZQ
{

#ifndef STRING_BUFFER_LENGTH
#define STRING_BUFFER_LENGTH 2000
#endif

#ifndef IMAGE3D_BUFFER_SIZE
#define IMAGE3D_BUFFER_SIZE 3
#endif

	class ZQ_OpticalFlow3DBuffer
	{
		typedef ZQ_DImage3D<float> DImage3D;
	public:
		ZQ_OpticalFlow3DBuffer();
		~ZQ_OpticalFlow3DBuffer();

	private:
		char input_fold[STRING_BUFFER_LENGTH];
		char output_fold[STRING_BUFFER_LENGTH];
		char image_prefix[STRING_BUFFER_LENGTH];
		int image_num;
		int base_index;

	private://for image
		int buf_image_index[IMAGE3D_BUFFER_SIZE];
		int buf_image_level[IMAGE3D_BUFFER_SIZE];
		int buf_image_live_time[IMAGE3D_BUFFER_SIZE];
		DImage3D buf_image[IMAGE3D_BUFFER_SIZE];

	private://for u,v,w
		int buf_uvw_index[IMAGE3D_BUFFER_SIZE];
		int buf_uvw_level[IMAGE3D_BUFFER_SIZE];
		int buf_uvw_live_time[IMAGE3D_BUFFER_SIZE];
		DImage3D buf_u[IMAGE3D_BUFFER_SIZE];
		DImage3D buf_v[IMAGE3D_BUFFER_SIZE];
		DImage3D buf_w[IMAGE3D_BUFFER_SIZE];

	public:

		int GetImageNum(){return image_num;}

		bool SetFilePara(const char* in_fold, const char* out_fold, const char* img_prefix, int img_num, int base_idx);

		bool BuildPyramids(double ratio,int minWidth, int& levels, float& cuda_cost_time);

		bool ClearPyramids(int levels);

		bool ClearIntermediateUVW(int levels);

		bool GetImage(DImage3D& out_img, int level,int i);

		bool GetUVW(DImage3D& u, DImage3D& v, DImage3D& w, int level, int i);

		bool WriteUVW(DImage3D& u, DImage3D& v, DImage3D& w, int level, int i);

		bool WriteFinalResult(DImage3D& flow, DImage3D& warp,int i);
	};
}

#endif