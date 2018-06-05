#ifndef _ZQ_SMOKE_OCTREE_GRID_2D_H_
#define _ZQ_SMOKE_OCTREE_GRID_2D_H_

#include "ZQ_SmokeSimulationUtils2D.h"

class ZQ_SmokeOctreeGrid2D
{
	class Info
	{
	public:
		int index;
		int level1;
		int i1,j1;
		int level2;
		int i2,j2;
		Info(){memset(this,0,sizeof(Info));}
		Info(int _lv1, int _i1, int _j1, int _lv2, int _i2, int _j2) : level1(_lv1), i1(_i1), j1(_j1),
			level2(_lv2), i2(_i2), j2(_j2)
		{index = -1;}
	};
public:
	ZQ_SmokeOctreeGrid2D();
	~ZQ_SmokeOctreeGrid2D();

private:
	enum CUDA_OpenOctreePoissonSOR_Type
	{
		CUDA_OPEN_OCTREE_POISSON_SOR_TYPE1,
		CUDA_OPEN_OCTREE_POISSON_SOR_TYPE2,
		CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3
	};

	static const CUDA_OpenOctreePoissonSOR_Type cuda_openoctree_poisson_type = CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3;

	static const bool guiding_specialRegion_useCUDA = true;

	//int uCount;
	//int vCount;
	//int lamdaCount;
	int k_index;
	int dim;
	double* deltax;
	float* deltaU;
	float* deltaV;
	int* Uidx2;
	int* Vidx2;
	int* Uidx4;
	int* Vidx4;
	int* Uidx8;
	int* Vidx8;
	std::vector<Info> Uinfos;
	std::vector<Info> Vinfos;
	std::vector<Info> lamdainfos;

	/*** for CUDA open octree poisson**********/
	std::vector<int> level0_index_red;
	std::vector<int> level0_neighborinfo_red;
	std::vector<int> level0_index_black;
	std::vector<int> level0_neighborinfo_black;
	std::vector<int> level1_index_red;
	std::vector<int> level1_neighborinfo_red;
	std::vector<int> level1_index_black;
	std::vector<int> level1_neighborinfo_black;
	std::vector<int> level2_index_red;
	std::vector<int> level2_neighborinfo_red;
	std::vector<int> level2_index_black;
	std::vector<int> level2_neighborinfo_black;
	std::vector<int> level3_index_red;
	std::vector<int> level3_neighborinfo_red;
	std::vector<int> level3_index_black;
	std::vector<int> level3_neighborinfo_black;
	/******************************************/

	bool* leafLevel0;
	bool* leafLevel1;
	bool* leafLevel2;
	bool* leafLevel3;

	int* lamda_index[4];

	taucs_ccs_matrix* A;

	float* subMatrix2;
	float* subMatrix4;
	float* subMatrix8;

	static const int min_default_thresh1 = 0;
	static const int min_default_thresh2 = 0;
	static const int min_default_thresh3 = 0;
	int level_dis_thresh1;
	int level_dis_thresh2;
	int level_dis_thresh3;

	bool* specialRegion_mask;
	int* specialRegion_label;
	int specialRegion_num; //when using CUDA, this var must have value
	std::vector<taucs_ccs_matrix*> specialRegion_A;


public:
	bool InitFromGlobalGrid(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist = 0);
	bool ReBuildFromGlobalGrid(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist = 0);
	bool GlobalProjection(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);
	bool ApplyGuidingVelocity(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const ZQ_SmokeGuideGrid2D* guideGrid);
	bool SubProjection(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);

private:
	
	
	bool _buildGlobalMatrix(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);
	bool _solveGlobal(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);

	//build global matrix
	bool _buildGlobalMatrixOpenPoissonSOR(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);
	bool _buildGlobalMatrixClosedPoissonSOR(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);
	bool _buildGlobalMatrixOpenPoisson(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);
	bool _buildGlobalMatrixClosedPoisson(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);
	bool _buildGlobalMatrixOpenFlux(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);
	bool _buildGlobalMatrixClosedFlux(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);



	//solve global
	bool _solveGlobalOpenPoissonSOR(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);
	bool _solveGlobalClosedPoissonSOR(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);
	bool _solveGlobalOpenPoisson(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);
	bool _solveGlobalClosedPoisson(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);
	bool _solveGlobalOpenFlux(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);
	bool _solveGlobalClosedFlux(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para);

	// build sub
	bool _buildSubMatrix();

	void _seedFilling(const int width, const int height, int* distanceField);
	void _find_special_label(const int width, const int height, const bool* special_mask, int* special_label, int& region_num);
	bool _buildOctreeInfo(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist);

	void _selectInfos_lambda(const int width, const int height, const bool* leafLevel0, const bool* leafLevel1, const bool* leafLevel2, const bool* leafLevel3,
		std::vector<Info>& lamdainfos);

	void _selectInfosForCPUSolvers(const int width, const int height, const bool* leafLevel0, const bool* leafLevel1, const bool* leafLevel2, const bool* leafLevel3,
		std::vector<Info>& Uinfos, std::vector<Info>& Vinfos, std::vector<Info>& lamdainfos);
	
	void _selectInfos_for_CUDAOpenPoissonRedBlack2(const int width, const int height, const bool* leafLevel0, const bool* leafLevel1, 
		const bool* leafLevel2, const bool* leafLevel3, std::vector<Info>& lamdainfos,
		std::vector<int>& level0_index_red,	std::vector<int>& level0_index_black,
		std::vector<int>& level1_index_red,	std::vector<int>& level1_index_black,
		std::vector<int>& level2_index_red,	std::vector<int>& level2_index_black,
		std::vector<int>& level3_index_red,	std::vector<int>& level3_index_black);

	void _selectInfos_for_CUDAOpenPoissonRedBlack3(const int width, const int height, const bool* leafLevel0, const bool* leafLevel1,
		const bool* leafLevel2, const bool* leafLevel3, std::vector<Info>& lamdainfos, 
		std::vector<int>& level0_index_red,		std::vector<int>& level0_neighborinfo_red,
		std::vector<int>& level0_index_black,	std::vector<int>& level0_neighborinfo_black,
		std::vector<int>& level1_index_red,		std::vector<int>& level1_neighborinfo_red,
		std::vector<int>& level1_index_black,	std::vector<int>& level1_neighborinfo_black,
		std::vector<int>& level2_index_red,		std::vector<int>& level2_neighborinfo_red,
		std::vector<int>& level2_index_black,	std::vector<int>& level2_neighborinfo_black,
		std::vector<int>& level3_index_red,		std::vector<int>& level3_neighborinfo_red,
		std::vector<int>& level3_index_black,	std::vector<int>& level3_neighborinfo_black
		);
	
	

	
};

#endif