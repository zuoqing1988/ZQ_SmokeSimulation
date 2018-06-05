#ifndef _ZQ_MERGE_IMAGE_H_
#define _ZQ_MERGE_IMAGE_H_

#include "ZQ_MergeImageOptions.h"
#include "ZQ_DoubleImage.h"
#include "ZQ_ImageIO.h"
#include "ZQ_MergeImage.h"
#include "ZQ_MathBase.h"
#include "ZQ_CPURenderer2DWorkspace.h""
#include "ZQ_TextureSampler.h"
#include "ZQ_TaucsBase.h"

namespace ZQ
{
	class ZQ_MergeImage
	{
		typedef ZQ_DImage<float> DImage;
		typedef ZQ_Vec2D Vec2;
	public:
		static bool Go(const ZQ_MergeImageOptions& opt);
	private:
		static bool Go_MergeDirectly(const ZQ_MergeImageOptions& opt);
		static bool Go_MergeDensity(const ZQ_MergeImageOptions& opt);
		static bool Go_MergeSourcePatch(const ZQ_MergeImageOptions& opt);
		static bool Go_ImageBlur(const ZQ_MergeImageOptions& opt);
		static bool Go_MergeLowHigh(const ZQ_MergeImageOptions& opt);

	public:
		static void MakeMatrix(const Vec2& scale, const float rot_rad, const Vec2& trans, float output_mat[9]);

	public:
		static bool MergeDirectly(DImage& output, const std::vector<DImage>& inputs,bool display);

		static bool MergeDensity(DImage& output, const std::vector<DImage>& inputs, const std::vector<Vec2>& target_size, const std::vector<float>& angles, const std::vector<Vec2>& trans, const bool blend_mode, const bool yAixsUp, const bool display);
	private:
		static bool _renderToBackground(ZQ::ZQ_CPURenderer2DWorkspace& renderer, const DImage& source, const Vec2& target_size, const float rot_rad, const Vec2& trans);

	public:
		static bool MergeSourcePatch(DImage& output, const DImage& source, const DImage& patch, const DImage& mask);

		static void DecomposeByBlur(DImage& low, DImage& high, const DImage& input, const float sigma, const int fsize);


	private: /* IO */
		static bool _load(DImage& img, const char* file, const int isColor);

		static bool _save(const DImage&img, const char* file);

	};
}

#endif
