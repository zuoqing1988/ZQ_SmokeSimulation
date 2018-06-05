#include "ZQ_SmokeGuideGrid2D.h"

ZQ_SmokeGuideGrid2D::ZQ_SmokeGuideGrid2D(const unsigned int width, const unsigned int height)
{
	this->width = width;
	this->height = height;
	u.allocate(width+1,height);
	v.allocate(width,height+1);
	temperature.allocate(width,height);
	density.allocate(width,height);
}

ZQ_SmokeGuideGrid2D::~ZQ_SmokeGuideGrid2D()
{

}

bool ZQ_SmokeGuideGrid2D::LoadFromFile(const char *file)
{
	if(file == 0)
		return false;
	FILE* in = fopen(file,"rb");
	if(in == 0)
		return false;
	if(fread(u.data(),sizeof(float),height*(width+1),in) != height*(width+1))
	{
		fclose(in);
		return false;
	}

	if(fread(v.data(),sizeof(float),width*(height+1),in) != width*(height+1))
	{
		fclose(in);
		return false;
	}

	if(fread(density.data(),sizeof(float),width*height,in) != width*height)
	{
		fclose(in);
		return false;
	}
	if(fread(temperature.data(),sizeof(float),width*height,in) != width*height)
	{
		fclose(in);
		return false;
	}
	fclose(in);
	return true;
}