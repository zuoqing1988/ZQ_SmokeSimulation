#ifndef _ZQ_SMOKE_SIMULATION_2D_H_
#define _ZQ_SMOKE_SIMULATION_2D_H_
#pragma once

#include "ZQ_SmokeSimulationUtils2D.h"
#include "ZQ_SmokeOctreeGrid2D.h"
#include "ZQ_PCGSolver.h"
#include "windows.h"

#ifndef _m_pi_
#define _m_pi_ 3.1415926535f
#endif

class ZQ_SmokeSimulation2D
{
public:
	ZQ_SmokeSimulation2D();
	~ZQ_SmokeSimulation2D();

private:
	ZQ_SmokeSimulationPara2D para;
	ZQ_SmokeOctreeGrid2D* octreeGrid;
	ZQ_SmokeGrid2D* fineGrid;
	ZQ_SmokeGrid2D* coarseGrid;
	int* interestingRegionDist;
	taucs_ccs_matrix* coarseA;
	float* submatrix;
	int subrow;
	int subcol;

	std::vector<ZQ_SmokeMovingObject2D*> mvobjs;

public:
	bool Init(const char* file);
	bool Run(const char* fold);
	bool UpdateOneFrame(const ZQ_SmokeGuideGrid2D* guideGrid);

private:
	bool _initOctree();
	bool _initRegular();

	bool _rebuild();
	bool _projection(const ZQ_SmokeGuideGrid2D* guideGrid);
	bool _totalProjection();
	bool _globalProjection();


};

#endif