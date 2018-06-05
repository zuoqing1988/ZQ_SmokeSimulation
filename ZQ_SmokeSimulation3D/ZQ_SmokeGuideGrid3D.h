#ifndef _ZQ_SMOKE_GUIDE_GRID_3D_H_
#define _ZQ_SMOKE_GUIDE_GRID_3D_H_
#pragma once

#include "ZQ_DoubleImage3D.h"

class ZQ_SmokeGuideGrid3D
{
public:
	ZQ_SmokeGuideGrid3D(const unsigned int width, const unsigned int height, const unsigned int depth);
	~ZQ_SmokeGuideGrid3D();
private:
	unsigned int width, height,depth;
	ZQ::ZQ_DImage3D<float> u;
	ZQ::ZQ_DImage3D<float> v;
	ZQ::ZQ_DImage3D<float> w;
	ZQ::ZQ_DImage3D<float> density;
	ZQ::ZQ_DImage3D<float> temperature;

public:
	unsigned int GetWidth() const {return width;}
	unsigned int GetHeight() const {return height;}
	unsigned int GetDepth() const {return depth;}

	bool LoadFromFile(const char* file);
	const float*& GetUptr() const{return u.data();}
	const float*& GetVptr() const{return v.data();}
	const float*& GetWptr() const{return w.data();}
	const float*& GetDensityPtr() const {return density.data();}
	const float*& GetTemperaturePtr() const {return temperature.data();}
};

#endif