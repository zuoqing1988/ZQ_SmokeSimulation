#ifndef _ZQ_SMOKE_GRID_3D_H_
#define _ZQ_SMOKE_GRID_3D_H_
#pragma once

#include "ZQ_DoubleImage3D.h"

class ZQ_SmokeGrid3D
{
public:
	ZQ_SmokeGrid3D(const unsigned int width, const unsigned int height, const unsigned int depth);
	~ZQ_SmokeGrid3D();

private:
	unsigned int width, height,depth;

	ZQ::ZQ_DImage3D<bool> occupy;
	ZQ::ZQ_DImage3D<float> pressure;
	ZQ::ZQ_DImage3D<float> temperature;
	ZQ::ZQ_DImage3D<float> density;

	ZQ::ZQ_DImage3D<float> u;
	ZQ::ZQ_DImage3D<float> v;
	ZQ::ZQ_DImage3D<float> w;

	ZQ::ZQ_DImage3D<float> unOccupyRatioU;
	ZQ::ZQ_DImage3D<float> unOccupyRatioV;	
	ZQ::ZQ_DImage3D<float> unOccupyRatioW;
	ZQ::ZQ_DImage3D<float> unOccupyVolume;

public:
	int GetWidth() const {return width;}
	int GetHeight() const {return height;}
	int GetDepth() const {return depth;}

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
	float*& GetVelocityWptr() {return w.data();}
	const float*& GetVelocityWptr() const {return w.data();}
	float*& GetFaceUnOccupyRatioUPtr() {return unOccupyRatioU.data();}
	const float*& GetFaceUnOccupyRatioUPtr() const {return unOccupyRatioU.data();}
	float*& GetFaceUnOccupyRatioVPtr() {return unOccupyRatioV.data();}
	const float*& GetFaceUnOccupyRatioVPtr() const {return unOccupyRatioV.data();}
	float*& GetFaceUnOccupyRatioWPtr() {return unOccupyRatioW.data();}
	const float*& GetFaceUnOccupyRatioWPtr() const {return unOccupyRatioW.data();}
	float*& GetVolumeUnOccupyRatioPtr() {return unOccupyVolume.data();}
	const float*& GetVolumeUnOccupyRatioPtr() const {return unOccupyVolume.data();}

	void BuildFaceUnOccupyRatio();
	void ClearFaceUnOccupyRatio() {unOccupyRatioU.clear(); unOccupyRatioV.clear(); unOccupyRatioW.clear();}
	void BuildVolumeUnOccupyRatio();
	void ClearVolumeUnOccupyRatio() {unOccupyVolume.clear();}

	ZQ_SmokeGrid3D* Clone() const;

	bool ExportDensity(const char* file, float quality = 0.9999) const;
	bool ExportOccupy(const char* file) const;
	bool ExportVelocity(const char* file) const;
	bool ExportGuide(const char* file) const;

};

#endif