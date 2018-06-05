#ifndef _ZQ_SMOKE_SIMULATION_3D_H_
#define _ZQ_SMOKE_SIMULATION_3D_H_
#pragma once

#include "ZQ_SmokeSimulationUtils3D.h"
#include "ZQ_SmokeOctreeGrid3D.h"
#include "windows.h"

#ifndef _m_pi_
#define _m_pi_ 3.1415926535f
#endif

class ZQ_SmokeSimulation3D
{
public:
	ZQ_SmokeSimulation3D();
	~ZQ_SmokeSimulation3D();

private:
	ZQ_SmokeSimulationPara3D para;
	ZQ_SmokeOctreeGrid3D* octreeGrid;
	ZQ_SmokeGrid3D* fineGrid;
	ZQ_SmokeGrid3D* coarseGrid;
	int* interestingRegionDist;
	taucs_ccs_matrix* coarseA;
	float* submatrix;
	int subrow;
	int subcol;

	std::vector<ZQ_SmokeMovingObject3D*> mvobjs;

public:
	bool Init(const char* file);
	bool Run(const char* fold);
	bool UpdateOneFrame(const ZQ_SmokeGuideGrid3D* guideGrid);

private:
	bool _initOctree();
	bool _initRegular();

	bool _rebuild();
	bool _projection(const ZQ_SmokeGuideGrid3D* guideGrid);
	bool _totalProjection();
	bool _globalProjection();


};

#endif