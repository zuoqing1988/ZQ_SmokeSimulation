#ifndef _ZQ_SMOKE_GRID_2D_H_
#define _ZQ_SMOKE_GRID_2D_H_
#pragma once

#include "ZQ_DoubleImage.h"

class ZQ_SmokeGrid2D
{
public:
	ZQ_SmokeGrid2D(const unsigned int width, const unsigned int height);
	~ZQ_SmokeGrid2D();

private:
	unsigned int width, height;

	ZQ::ZQ_DImage<bool> occupy;
	ZQ::ZQ_DImage<float> pressure;
	ZQ::ZQ_DImage<float> temperature;
	ZQ::ZQ_DImage<float> density;

	ZQ::ZQ_DImage<float> u;
	ZQ::ZQ_DImage<float> v;
	
	ZQ::ZQ_DImage<float> unOccupyRatioU;
	ZQ::ZQ_DImage<float> unOccupyRatioV;	
	ZQ::ZQ_DImage<float> unOccupyVolume;

public:
	int GetWidth() const {return width;}
	int GetHeight() const {return height;}

	bool*& GetOccupyPtr() {return occupy.data();}
	const bool*& GetOccupyPtr() const {return occupy.data();}
	float*& GetPressurePtr() {return pressure.data();}
	const float*& GetPressurePtr() const {return pressure.data();}
	float*& GetTemperaturePtr() {return temperature.data();}
	const float*& GetTemperaturePtr() const {return temperature.data();}
	float*& GetDensityPtr() {return density.data();}
	const float*& GetDensityPtr() const {return density.data();}
	float*& GetVelocityUptr() {return u.data();}
	const float*& GetVelocityUptr() const {return u.data();}
	float*& GetVelocityVptr() {return v.data();}
	const float*& GetVelocityVptr() const {return v.data();}
	float*& GetFaceUnOccupyRatioUPtr() {return unOccupyRatioU.data();}
	const float*& GetFaceUnOccupyRatioUPtr() const {return unOccupyRatioU.data();}
	float*& GetFaceUnOccupyRatioVPtr() {return unOccupyRatioV.data();}
	const float*& GetFaceUnOccupyRatioVPtr() const {return unOccupyRatioV.data();}
	float*& GetVolumeUnOccupyRatioPtr() {return unOccupyVolume.data();}
	const float*& GetVolumeUnOccupyRatioPtr() const {return unOccupyVolume.data();}

	void BuildFaceUnOccupyRatio();
	void ClearFaceUnOccupyRatio() {unOccupyRatioU.clear(); unOccupyRatioV.clear();}
	void BuildVolumeUnOccupyRatio();
	void ClearVolumeUnOccupyRatio() {unOccupyVolume.clear();}

	ZQ_SmokeGrid2D* Clone() const;

	bool ExportDensity(const char* file) const;
	bool ExportOccupy(const char* file) const;
	bool ExportVelocity(const char* file) const;
	bool ExportGuide(const char* file) const;

};

#endif