#ifndef _ZQ_SMOKE_GUIDE_GRID_2D_H_
#define _ZQ_SMOKE_GUIDE_GRID_2D_H_
#pragma once

#include "ZQ_DoubleImage.h"

class ZQ_SmokeGuideGrid2D
{
public:
	ZQ_SmokeGuideGrid2D(const unsigned int width, const unsigned int height);
	~ZQ_SmokeGuideGrid2D();
private:
	unsigned int width, height;
	ZQ::ZQ_DImage<float> u;
	ZQ::ZQ_DImage<float> v;
	ZQ::ZQ_DImage<float> density;
	ZQ::ZQ_DImage<float> temperature;

public:
	unsigned int GetWidth() const {return width;}
	unsigned int GetHeight() const {return height;}

	bool LoadFromFile(const char* file);
	const float*& GetUptr() const{return u.data();}
	const float*& GetVptr() const{return v.data();}
	const float*& GetDensityPtr() const {return density.data();}
	const float*& GetTemperaturePtr() const {return temperature.data();}
};

#endif