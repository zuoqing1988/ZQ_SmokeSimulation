#include "ZQ_SmokeGrid3D.h"
#include "ZQ_Wavelet.h"
#include "ZQ_CompressedImage.h"

ZQ_SmokeGrid3D::ZQ_SmokeGrid3D(const unsigned int width, const unsigned int height, const unsigned int depth)
{
	this->width = width;
	this->height = height;
	this->depth = depth;
	occupy.allocate(width,height,depth);
	pressure.allocate(width,height,depth);
	temperature.allocate(width,height,depth);
	density.allocate(width,height,depth);

	u.allocate(width+1,height,depth);
	v.allocate(width,height+1,depth);
	w.allocate(width,height,depth+1);
}

ZQ_SmokeGrid3D::~ZQ_SmokeGrid3D()
{

}

void ZQ_SmokeGrid3D::BuildFaceUnOccupyRatio()
{
	unOccupyRatioU.allocate(width+1,height,depth);
	unOccupyRatioV.allocate(width,height+1,depth);
	unOccupyRatioW.allocate(width,height,depth+1);
}

void ZQ_SmokeGrid3D::BuildVolumeUnOccupyRatio()
{
	unOccupyVolume.allocate(width,height,depth);
}

ZQ_SmokeGrid3D* ZQ_SmokeGrid3D::Clone() const
{
	ZQ_SmokeGrid3D* tmp = new ZQ_SmokeGrid3D(width,height,depth);
	tmp->occupy = occupy;
	tmp->pressure = pressure;
	tmp->temperature = temperature;
	tmp->density = density;
	tmp->u = u;
	tmp->v = v;
	tmp->w = w;
	tmp->unOccupyRatioU = unOccupyRatioU;
	tmp->unOccupyRatioV = unOccupyRatioV;
	tmp->unOccupyRatioW = unOccupyRatioW;
	tmp->unOccupyVolume = unOccupyVolume;

	return tmp;
}

bool ZQ_SmokeGrid3D::ExportDensity(const char* file, float quality) const
{
	if(file == 0)
		return false;
	int len = strlen(file);

	char suffix_di3[] = ".di3";
	char suffix_zqwv[] = ".zqwv";
	char suffix_zqci[] = ".zqci";
	int len_di3 = strlen(suffix_di3);
	int len_zqwv = strlen(suffix_zqwv);
	int len_zqci = strlen(suffix_zqci);

	if(len >= len_di3 && _strcmpi(file+len-len_di3,suffix_di3) == 0)
	{
		return density.saveImage(file);
	}
	else if(len >= len_zqwv && _strcmpi(file+len-len_zqwv,suffix_zqwv) == 0)
	{
		return ZQ::ZQ_Wavelet<float>::SaveWaveletFile(file,density,1000,quality);
	}
	else if(len >= len_zqci && _strcmpi(file+len-len_zqci,suffix_zqci) == 0)
	{
		return ZQ::ZQ_CompressedImage::SaveCompressedImage(file,density,1000,quality);
	}
	return false;
}

bool ZQ_SmokeGrid3D::ExportGuide(const char* file) const
{
	if(file == 0)
		return false;
	FILE* out = fopen(file,"wb");
	if(out == 0)
		return false;
	if(fwrite(u.data(),sizeof(float),(width+1)*height*depth,out) != (width+1)*height*depth)
	{
		fclose(out);
		return false;
	}

	if(fwrite(v.data(),sizeof(float),width*(height+1)*depth,out) != width*(height+1)*depth)
	{
		fclose(out);
		return false;
	}

	if(fwrite(w.data(),sizeof(float),width*height*(depth+1),out) != width*height*(depth+1))
	{
		fclose(out);
		return false;
	}

	if(fwrite(density.data(),sizeof(float),width*height*depth,out) != width*height*depth)
	{
		fclose(out);
		return false;
	}
	if(fwrite(temperature.data(),sizeof(float),width*height*depth,out) != width*height*depth)
	{
		fclose(out);
		return false;
	}
	fclose(out);

	return true;
}

bool ZQ_SmokeGrid3D::ExportOccupy(const char* file) const
{
	return true;
}

bool ZQ_SmokeGrid3D::ExportVelocity(const char* file) const
{
	return true;
}
