#ifndef _ZQ_SMOKE_SIMULATION_UTILS_2D_H_
#define _ZQ_SMOKE_SIMULATION_UTILS_2D_H_
#pragma once

#include "ZQ_SmokeGrid2D.h"
#include "ZQ_SmokeMovingObject2D.h"
#include "ZQ_SmokeGuideGrid2D.h"

#include "ZQ_CUDA_PoissonSolver2D.h"
#include "ZQ_CUDA_Advection2D.h"
#include "ZQ_taucs.h"

enum GlobalSolverType
{
	CLOSED_POISSON,
	CLOSED_FLUX,
	OPEN_POISSON,
	OPEN_FLUX,
	COARSE_OPEN_POISSON_SOR,
	COARSE_CLOSED_POISSON_SOR,
	COARSE_OPEN_FLUX_SOR,
	COARSE_CLOSED_FLUX_SOR,
	TOTAL_OPEN_POISSON_SOR,
	TOTAL_CLOSED_POISSON_SOR,
	TOTAL_OPEN_FLUX_SOR,
	TOTAL_CLOSED_FLUX_SOR,
	CLOSED_OCTREE_POISSON,
	CLOSED_OCTREE_FLUX,
	OPEN_OCTREE_POISSON,
	OPEN_OCTREE_FLUX,
	OPEN_OCTREE_POISSON_SOR,
	CLOSED_OCTREE_POISSON_SOR
};


class ZQ_SmokeSource2D
{
public:
	float cx,cy;
	float half_xlen,half_ylen;
	float reinjectDensity;
	float reinjectTemperature;
};

class ZQ_SmokeSimulationPara2D
{
public:
	enum InternalObject{INT_OBJ_NULL,INT_OBJ_SPHERE,INT_OBJ_BOX,INT_OBJ_OVAL,INT_OBJ_SPECIAL1};

public:
	ZQ_SmokeSimulationPara2D();
	~ZQ_SmokeSimulationPara2D();

public:
	int width, height;
	int coarseWidth,coarseHeight;
	int subSize;
	float gravityCoeff;						//gravity coefficient, unused
	float buoyancyCoeff;					//buoyancy coefficient
	float confineCoeff;						//confinement force coefficient
	std::vector<ZQ_SmokeSource2D> sources;	//smoke sources
	float velAttenCoeff;					//velocity attenuation coefficient
	float Tamb;								//ambient temperature
	float voxelSize;						//voxelSize
	float stepTime;							//step time
	float tempAttenCoeff;					//temperature attenuation coefficient
	float densityAttenCoeff;				//density attenuation coefficient
	int steps;								//steps for advection
	int frameCount;							//frames
	int maxIter;							//max iterations for solvers
	int innerIter;							//inner iteration for solvers
	GlobalSolverType globalGridSolver;		//global grid solver
	int subGridSolver;						//sub grid solver
	int octree_thresh1;						//octree thresh level 1
	int octree_thresh2;						//octree thresh level 2
	int octree_thresh3;						//octree thresh level 3
	bool useBFECC;							//using BFECC or not
	bool useMovingObject;					//using moving object or not
	char movingObjectFile[200];				//moving object file name
	bool useGuidingControl;					//using guiding control or not
	float guidingCoeff;						//guiding coefficient
	bool guidingNaiveBlend;					//naive blend for guiding
	char guidingFold[200];					//guiding information fold
	bool exportGuide;						//export current simulation as guide grid
	bool exportOccupy;						//export occupy information
	bool randDensity;						//use random density for source 
	float maxDensity;						//max density
	float maxTemperature;					//max temperature
	InternalObject intobj;					//internal object, disabled when using moving object 

public:
	bool LoadFromFile(const char* file, bool display = true);
};

//CUDA UTILS
extern "C" void SubProjectionCuda(const int subRow, const int subCol, const int fast_count, const float* subMatrix, const float* input, float* output);


int ZQ_InitMatOpenPoisson2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A);
int ZQ_InitMatClosedPoisson2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A);
int ZQ_InitMatOpenFlux2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A);
int ZQ_InitMatClosedFlux2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A);
bool ZQ_InitSubMatPoisson2D(const int subSize, float** subMatrix,int* row, int* col);
bool ZQ_SolveSubGrid2D(const ZQ_SmokeSimulationPara2D& para, const bool* occupy, const float* input, float* output);
bool ZQ_InitFineGrid2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_InitCoarseGrid2D(const ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid);
bool ZQ_InitSimulation2D(ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, taucs_ccs_matrix** coarseA,
						 float** subMatrix, int* subRow, int* subCol);
bool ZQ_ApplyMovingObject2D(ZQ_SmokeGrid2D* fineGrid, const std::vector<ZQ_SmokeMovingObject2D*> mvobjs, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_AddForce2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_Attenuation2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_VelocityAdvect2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_VelocityAdvectBFECC2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_MapFine2Coarse2D(const ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_SolveOpenPressure2D(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA);
bool ZQ_SolveOpenPressure2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_SolveClosedPressure2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_SolveOpenFlux2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_SolveClosedFlux2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_SolveClosedPressure2D(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA);
bool ZQ_SolveCoarseOpenPoisson2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid);
bool ZQ_SolveCoarseClosedPoisson2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid);
bool ZQ_SolveCoarseOpenFlux2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid);
bool ZQ_SolveCoarseClosedFlux2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid);
bool ZQ_AdjustOpenVelocity2D(const ZQ_SmokeGrid2D* coarseGrid, ZQ_SmokeGrid2D* deltaCoarseGrid);
bool ZQ_AdjustClosedVelocity2D(const ZQ_SmokeGrid2D* coarseGrid, ZQ_SmokeGrid2D* deltaCoarseGrid);
bool ZQ_SolveFlux2D(const ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA,  ZQ_SmokeGrid2D* deltaCoarseGrid);

/*deltaCoarseGrid is in and out*/
bool ZQ_ApplyGuidingVelocity2D(const ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeGuideGrid2D* guideGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaCoarseGrid);

bool ZQ_MapCoarse2Fine2D(const ZQ_SmokeGrid2D* deltaCoarseGrid, ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_SubProjection2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para, const float* subMatrix, int subRow, int subCol, bool useGPUflag = false);
bool ZQ_TotalProjection2D(ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para,  const taucs_ccs_matrix* coarseA,
						  const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid2D* guideGrid);
bool ZQ_DensityTemperatureAdvect2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_DensityTemperatureAdvectBFECC2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_ReinjectTemperatureDensity2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para);
bool ZQ_UpdateOneFrame2D(ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA,
						 const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid2D* guideGrid);

#endif