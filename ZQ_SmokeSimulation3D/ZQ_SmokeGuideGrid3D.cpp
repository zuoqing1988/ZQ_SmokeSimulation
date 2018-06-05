#include "ZQ_SmokeGuideGrid3D.h"

ZQ_SmokeGuideGrid3D::ZQ_SmokeGuideGrid3D(const unsigned int width, const unsigned int height, const unsigned int depth)
{
	this->width = width;
	this->height = height;
	this->depth = depth;
	u.allocate(width+1,height,depth);
	v.allocate(width,height+1,depth);
	w.allocate(width,height,depth+1);
	temperature.allocate(width,height,depth);
	density.allocate(width,height,depth);
}

ZQ_SmokeGuideGrid3D::~ZQ_SmokeGuideGrid3D()
{

}

bool ZQ_SmokeGuideGrid3D::LoadFromFile(const char *file)
{
	if(file == 0)
		return false;
	FILE* in = fopen(file,"rb");
	if(in == 0)
		return false;
	if(fread(u.data(),sizeof(float),(width+1)*height*depth,in) != (width+1)*height*depth)
	{
		fclose(in);
		return false;
	}

	if(fread(v.data(),sizeof(float),width*(height+1)*depth,in) != width*(height+1)*depth)
	{
		fclose(in);
		return false;
	}

	if(fread(w.data(),sizeof(float),width*height*(depth+1),in) != width*height*(depth+1))
	{
		fclose(in);
		return false;
	}

	if(fread(density.data(),sizeof(float),width*height*depth,in) != width*height*depth)
	{
		fclose(in);
		return false;
	}
	if(fread(temperature.data(),sizeof(float),width*height*depth,in) != width*height*depth)
	{
		fclose(in);
		return false;
	}
	fclose(in);
	return true;
}