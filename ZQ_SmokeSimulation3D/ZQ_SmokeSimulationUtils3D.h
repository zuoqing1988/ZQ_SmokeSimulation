#ifndef _ZQ_SMOKE_SIMULATION_UTILS_3D_H_
#define _ZQ_SMOKE_SIMULATION_UTILS_3D_H_
#pragma once

#include "ZQ_SmokeGrid3D.h"
#include "ZQ_SmokeMovingObject3D.h"
#include "ZQ_SmokeGuideGrid3D.h"
#include "ZQ_CUDA_PoissonSolver3D.h"
#include "ZQ_CUDA_Advection3D.h"
#include "ZQ_CUDA_AddForce3D.h"
#include "ZQ_CUDA_Attenuation3D.h"
#include "ZQ_taucs.h"


extern "C" float SubProjection3DCuda(const int subRow, const int subCol, const int fast_count, const float* subMatrix, const float* input, float* output);


enum GlobalSolverType3D
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


class ZQ_SmokeSource3D
{
public:
	float cx,cy,cz;
	float half_xlen,half_ylen,half_zlen;
	float reinjectDensity;
	float reinjectTemperature;
};

class ZQ_SmokeSimulationPara3D
{
public:
	enum InternalObject{ INT_OBJ_NULL, INT_OBJ_SPHERE, INT_OBJ_BOX, INT_OBJ_OVAL, INT_OBJ_SPECIAL1 };
	enum AdvectVelocityType{ 
		ADV_VEL_INMAC_OUTMAC, 
		ADV_VEL_INMAC_OUTMAC_BFECC, 
		ADV_VEL_INREG_OUTREG, 
		ADV_VEL_INREG_OUTREG_BFECC, 
		ADV_VEL_INREG_OUTMAC, 
		ADV_VEL_INREG_OUTMAC_BFECC 
	};
	enum AdvectScalarType{ 
		ADV_SCA_MAC, 
		ADV_SCA_MAC_BFECC,
		ADV_SCA_REG,
		ADV_SCA_REG_BFECC
	};
	enum ExportType{
		EXPORT_DI3,
		EXPORT_ZQWV,
		EXPORT_ZQCI
	};
	enum AddForceType{
		ADD_FORCE_ENTIRE,
		ADD_FORCE_SLICE
	};

public:
	ZQ_SmokeSimulationPara3D();
	~ZQ_SmokeSimulationPara3D();

public:
	int width, height, depth;
	int coarseWidth,coarseHeight,coarseDepth;
	int subSize;
	float gravityCoeff;						//gravity coefficient, unused
	float buoyancyCoeff;					//buoyancy coefficient
	float confineCoeff;						//confinement force coefficient
	std::vector<ZQ_SmokeSource3D> sources;	//smoke sources
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
	GlobalSolverType3D globalGridSolver;		//global grid solver
	int subGridSolver;						//sub grid solver
	int octree_thresh1;						//octree thresh level 1
	int octree_thresh2;						//octree thresh level 2
	int octree_thresh3;						//octree thresh level 3
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
	AdvectVelocityType advVelType;			//advect velocity type
	AdvectScalarType advScaType;			//advect scalar type
	ExportType exportType;					// export file type
	float export_quality;					// export file quality
	AddForceType addforce_type;				// addforce type

public:
	bool LoadFromFile(const char* file, bool display = true);
};

void MACtoRegular4(int width, int height, int depth, const float* mac_u, const float* mac_v, const float* mac_w, float* vel4);
void Regular4toMAC(int width, int height, int depth, const float* vel4, float* mac_u, float* mac_v, float* mac_w);
void ApplyAdvectVelocityResult(int width, int height, int depth, const bool* occupy, const float* adv_u, const float* adv_v, const float* adv_w, float* u, float* v, float* w);

int ZQ_InitMatOpenPoisson3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A);
int ZQ_InitMatClosedPoisson3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A);
int ZQ_InitMatOpenFlux3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A);
int ZQ_InitMatClosedFlux3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A);
bool ZQ_InitSubMatPoisson3D(const int subSize, float** subMatrix,int* row, int* col);
bool ZQ_SolveSubGrid3D(const ZQ_SmokeSimulationPara3D& para, const bool* occupy, const float* input, float* output);
bool ZQ_InitFineGrid3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_InitCoarseGrid3D(const ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid);
bool ZQ_InitSimulation3D(ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, taucs_ccs_matrix** coarseA,
						 float** subMatrix, int* subRow, int* subCol);
bool ZQ_ApplyMovingObject3D(ZQ_SmokeGrid3D* fineGrid, const std::vector<ZQ_SmokeMovingObject3D*> mvobjs, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_AddForce3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_Attenuation3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_AdvectVelocity3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_MapFine2Coarse3D(const ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_SolveOpenPressure3D(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA);
bool ZQ_SolveOpenPressure3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_SolveClosedPressure3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_SolveOpenFlux3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_SolveClosedFlux3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_SolveClosedPressure3D(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA);
bool ZQ_SolveCoarseOpenPoisson3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid);
bool ZQ_SolveCoarseClosedPoisson3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid);
bool ZQ_SolveCoarseOpenFlux3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid);
bool ZQ_SolveCoarseClosedFlux3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid);
bool ZQ_AdjustOpenVelocity3D(const ZQ_SmokeGrid3D* coarseGrid, ZQ_SmokeGrid3D* deltaCoarseGrid);
bool ZQ_AdjustClosedVelocity3D(const ZQ_SmokeGrid3D* coarseGrid, ZQ_SmokeGrid3D* deltaCoarseGrid);
bool ZQ_SolveFlux3D(const ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA,  ZQ_SmokeGrid3D* deltaCoarseGrid);

/*deltaCoarseGrid is in and out*/
bool ZQ_ApplyGuidingVelocity3D(const ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeGuideGrid3D* guideGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaCoarseGrid);

bool ZQ_MapCoarse2Fine3D(const ZQ_SmokeGrid3D* deltaCoarseGrid, ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_SubProjection3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para, const float* subMatrix, int subRow, int subCol, bool useGPUflag = false);
bool ZQ_TotalProjection3D(ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para,  const taucs_ccs_matrix* coarseA,
						  const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid3D* guideGrid);

bool ZQ_AdvectScalar3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_ReinjectTemperatureDensity3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para);
bool ZQ_UpdateOneFrame3D(ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA,
						 const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid3D* guideGrid);

#endif