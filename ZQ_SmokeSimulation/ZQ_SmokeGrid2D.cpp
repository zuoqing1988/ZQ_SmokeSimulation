#include "ZQ_SmokeGrid2D.h"

ZQ_SmokeGrid2D::ZQ_SmokeGrid2D(const unsigned int width, const unsigned int height)
{
	this->width = width;
	this->height = height;
	occupy.allocate(width,height);
	pressure.allocate(width,height);
	temperature.allocate(width,height);
	density.allocate(width,height);

	u.allocate(width+1,height);
	v.allocate(width,height+1);
}

ZQ_SmokeGrid2D::~ZQ_SmokeGrid2D()
{

}

void ZQ_SmokeGrid2D::BuildFaceUnOccupyRatio()
{
	unOccupyRatioU.allocate(width+1,height);
	unOccupyRatioV.allocate(width,height+1);
}

void ZQ_SmokeGrid2D::BuildVolumeUnOccupyRatio()
{
	unOccupyVolume.allocate(width,height);
}

ZQ_SmokeGrid2D* ZQ_SmokeGrid2D::Clone() const
{
	ZQ_SmokeGrid2D* tmp = new ZQ_SmokeGrid2D(width,height);
	tmp->occupy = occupy;
	tmp->pressure = pressure;
	tmp->temperature = temperature;
	tmp->density = density;
	tmp->u = u;
	tmp->v = v;
	tmp->unOccupyRatioU = unOccupyRatioU;
	tmp->unOccupyRatioV = unOccupyRatioV;
	tmp->unOccupyVolume = unOccupyVolume;

	return tmp;
}

bool ZQ_SmokeGrid2D::ExportDensity(const char* file) const
{
	return true;
}

bool ZQ_SmokeGrid2D::ExportGuide(const char* file) const
{
	if(file == 0)
		return false;
	FILE* out = fopen(file,"wb");
	if(out == 0)
		return false;
	if(fwrite(u.data(),sizeof(float),height*(width+1),out) != height*(width+1))
	{
		fclose(out);
		return false;
	}

	if(fwrite(v.data(),sizeof(float),width*(height+1),out) != width*(height+1))
	{
		fclose(out);
		return false;
	}

	if(fwrite(density.data(),sizeof(float),width*height,out) != width*height)
	{
		fclose(out);
		return false;
	}
	if(fwrite(temperature.data(),sizeof(float),width*height,out) != width*height)
	{
		fclose(out);
		return false;
	}
	fclose(out);

	return true;
}

bool ZQ_SmokeGrid2D::ExportOccupy(const char* file) const
{
	return true;
}

bool ZQ_SmokeGrid2D::ExportVelocity(const char* file) const
{
	return true;
}
