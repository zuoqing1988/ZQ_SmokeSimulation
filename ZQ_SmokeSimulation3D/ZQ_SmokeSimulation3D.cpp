#include "ZQ_SmokeSimulation3D.h"
#include "ZQ_PCGSolver.h"
#include "ZQ_TaucsBase.h"

using namespace ZQ;

ZQ_SmokeSimulation3D::ZQ_SmokeSimulation3D()
{
	octreeGrid = 0;
	fineGrid = 0;
	coarseGrid = 0;
	coarseA = 0;
	submatrix = 0;
	subrow = 0;
	subcol = 0;

	interestingRegionDist = 0;
}

ZQ_SmokeSimulation3D::~ZQ_SmokeSimulation3D()
{
	if(octreeGrid)
	{
		delete octreeGrid;
		octreeGrid = 0;
	}
	if(fineGrid)
	{
		delete fineGrid;
		fineGrid = 0;
	}
	if(coarseGrid)
	{
		delete coarseGrid;
		coarseGrid = 0;
	}
	if(interestingRegionDist)
	{
		delete []interestingRegionDist;
		interestingRegionDist = 0;
	}
	if(coarseA)
	{
		ZQ_TaucsBase::ZQ_taucs_ccs_free(coarseA);
		coarseA = 0;
	}
	if(submatrix)
	{
		delete []submatrix;
		submatrix = 0;
	}
	for(int i = 0;i < mvobjs.size();i++)
	{
		if(mvobjs[i])
			delete mvobjs[i];
	}
	mvobjs.clear();
}

bool ZQ_SmokeSimulation3D::Init(const char *file)
{
	LARGE_INTEGER tf;
	QueryPerformanceFrequency(&tf);


	if(!para.LoadFromFile(file,true))
		return false;
	if(para.useMovingObject)
	{
		if(!LoadMovingObjects3D(para.movingObjectFile,mvobjs))
			return false;
	}

	if(para.globalGridSolver == OPEN_OCTREE_FLUX || para.globalGridSolver == OPEN_OCTREE_POISSON
		|| para.globalGridSolver == CLOSED_OCTREE_FLUX || para.globalGridSolver == CLOSED_OCTREE_POISSON 
		|| para.globalGridSolver == OPEN_OCTREE_POISSON_SOR || para.globalGridSolver == CLOSED_OCTREE_POISSON_SOR)
	{
		LARGE_INTEGER init_octree_st, init_octree_end;
		QueryPerformanceCounter(&init_octree_st);
	
		if(!_initOctree())
		{
			return false;
		}
		QueryPerformanceCounter(&init_octree_end);
		printf("init octree cost:%f\n",1.0*(init_octree_end.QuadPart-init_octree_st.QuadPart)/tf.QuadPart);
	}
	else
	{
		LARGE_INTEGER init_regular_st, init_regular_end;
		QueryPerformanceCounter(&init_regular_st);
		if(!_initRegular())
		{
			return false;
		}
		QueryPerformanceCounter(&init_regular_end);
		printf("init regular cost:%f\n",1.0*(init_regular_end.QuadPart-init_regular_st.QuadPart)/tf.QuadPart);
	}
	return true;
}

bool ZQ_SmokeSimulation3D::Run(const char* fold)
{
	int width = para.width;
	int height = para.height;
	int depth = para.depth;

	char name[500] = "";

	ZQ_SmokeGuideGrid3D* guideGrid = new ZQ_SmokeGuideGrid3D(para.coarseWidth,para.coarseHeight,para.coarseDepth);
	for(int i = 0;i < para.frameCount;i++)
	{
		printf("frame [%d]\n",i);
		if(para.useGuidingControl)
		{
			char guidefile[200];
			sprintf(guidefile,"%s\\%d.guidedat",para.guidingFold,i);
			if(!guideGrid->LoadFromFile(guidefile))
			{
				printf("load guide file %s fail\n",guidefile);
				delete guideGrid;
				return false;
			}
		}
		if(!UpdateOneFrame(guideGrid))
			return false;
		switch (para.exportType)
		{
		case ZQ_SmokeSimulationPara3D::EXPORT_DI3:
			sprintf(name, "%s\\frame%d.di3", fold, i);
			break;
		case ZQ_SmokeSimulationPara3D::EXPORT_ZQCI:
			sprintf(name, "%s\\frame%d.zqci", fold, i);
			break;
		case ZQ_SmokeSimulationPara3D::EXPORT_ZQWV:
			sprintf(name, "%s\\frame%d.zqwv", fold, i);
			break;
		}
		
		fineGrid->ExportDensity(name,para.export_quality);
		if(para.exportGuide)
		{
			sprintf(name,"%s\\%d.guidedat",fold,i);
			fineGrid->ExportGuide(name);
		}
		if(para.exportOccupy)
		{
			sprintf(name,"%s\\%d.occdat",fold,i);
			fineGrid->ExportOccupy(name);
		}	
	}
	delete guideGrid;
	return true;
}

bool ZQ_SmokeSimulation3D::UpdateOneFrame(const ZQ_SmokeGuideGrid3D* guideGrid)
{
	LARGE_INTEGER tf;
	QueryPerformanceFrequency(&tf);

	static float moving_frame = 0;
	moving_frame +=1;

	/********************** apply boundary begin ***************************************/
	LARGE_INTEGER app_boundary_st, app_boundary_end;
	QueryPerformanceCounter(&app_boundary_st);
	if(para.useMovingObject)
	{
		if(moving_frame > 1) 
		{
			if(!ZQ_ApplyMovingObject3D(fineGrid,mvobjs,para))
				return false;
		}
	}
	
	QueryPerformanceCounter(&app_boundary_end);
	printf("app boundary cost:%f\n",1.0*(app_boundary_end.QuadPart-app_boundary_st.QuadPart)/tf.QuadPart);
	/***********************  apply boundary end *********************************/

	/********************** add force begin **********************************/
	LARGE_INTEGER add_force_st, add_force_end;
	QueryPerformanceCounter(&add_force_st);
	if(!ZQ_AddForce3D(fineGrid,para))
		return false;
	QueryPerformanceCounter(&add_force_end);
	printf("add force cost:%f\n",1.0*(add_force_end.QuadPart-add_force_st.QuadPart)/tf.QuadPart);
	
	LARGE_INTEGER attenuation_st, attenuation_end;
	QueryPerformanceCounter(&attenuation_st);
	if(!ZQ_Attenuation3D(fineGrid,para))
		return false;
	QueryPerformanceCounter(&attenuation_end);
	printf("attenuation cost:%f\n",1.0*(attenuation_end.QuadPart-attenuation_st.QuadPart)/tf.QuadPart);
	/****************************   add force end   ********************************/


	/************************   velocity advection begin    ***********************/
	LARGE_INTEGER advect_velocity_st, advect_velocity_end;
	QueryPerformanceCounter(&advect_velocity_st);
	
	if(!ZQ_AdvectVelocity3D(fineGrid,para))
		return false;

	QueryPerformanceCounter(&advect_velocity_end);
	printf("advect velocity cost:%f\n",1.0*(advect_velocity_end.QuadPart-advect_velocity_st.QuadPart)/tf.QuadPart);
	/************************   velocity advection end    ***************************/

	/*************************  projection begin ***************************************/
	
	if(!_projection(guideGrid))
	{
		return false;
	}
	/*************************  projection end ***************************************/


	/********************  scalar fields advection begin ***************************/
	LARGE_INTEGER advect_tempden_st, advect_tempden_end;
	QueryPerformanceCounter(&advect_tempden_st);
	

	if(!ZQ_AdvectScalar3D(fineGrid,para))
		return false;
	
	QueryPerformanceCounter(&advect_tempden_end);
	printf("advect tempden cost:%f\n",1.0*(advect_tempden_end.QuadPart-advect_tempden_st.QuadPart)/tf.QuadPart);

	/********************  scalar fields advection end ***************************/

	if(!ZQ_ReinjectTemperatureDensity3D(fineGrid,para))
		return false;
	if(para.useMovingObject)
	{
		for(int i = 0;i < mvobjs.size();i++)
		{
			if(mvobjs[i])
			{
				mvobjs[i]->UpdateOneStep();
			}
		}
	}
	return true;
}



bool ZQ_SmokeSimulation3D::_initOctree()
{
	fineGrid = new ZQ_SmokeGrid3D(para.width,para.height,para.depth);
	if(!ZQ_InitFineGrid3D(fineGrid,para))
	{
		delete fineGrid;
		fineGrid = 0;
		return false;
	}
	if(para.useMovingObject)
	{
		if(!ZQ_ApplyMovingObject3D(fineGrid,mvobjs,para))
			return false;
	}

	interestingRegionDist = new int[para.width*para.height*para.depth];
	memset(interestingRegionDist,0,sizeof(int)*para.width*para.height*para.depth);
	
	octreeGrid = new ZQ_SmokeOctreeGrid3D();


	static float moving_frame = 0;
	moving_frame +=1;

	LARGE_INTEGER tf;
	QueryPerformanceFrequency(&tf);
	LARGE_INTEGER init_octree_st, init_octree_end;
	QueryPerformanceCounter(&init_octree_st);

	if(!octreeGrid->InitFromGlobalGrid(fineGrid,para,0/*interestingRegionDist*/))
	{
		delete octreeGrid;
		octreeGrid = 0;
		delete fineGrid;
		fineGrid = 0;
		delete []interestingRegionDist;
		interestingRegionDist = 0;
		return false;
	}

	QueryPerformanceCounter(&init_octree_end);
	printf("init octree cost:%f\n",1.0*(init_octree_end.QuadPart-init_octree_st.QuadPart)/tf.QuadPart);

	return true;
}

bool ZQ_SmokeSimulation3D::_initRegular()
{
	fineGrid = new ZQ_SmokeGrid3D(para.width,para.height,para.depth);
	coarseGrid = new ZQ_SmokeGrid3D(para.coarseWidth,para.coarseHeight,para.coarseDepth);
	if(para.useMovingObject)
	{
		if(!ZQ_InitFineGrid3D(fineGrid,para))
		{
			delete fineGrid;
			fineGrid = 0;
			delete coarseGrid;
			coarseGrid = 0;
			return false;
		}
		if(!ZQ_ApplyMovingObject3D(fineGrid,mvobjs,para))
		{
			delete fineGrid;
			fineGrid = 0;
			delete coarseGrid;
			coarseGrid = 0;
			return false;
		}
		if(!ZQ_InitCoarseGrid3D(fineGrid,coarseGrid))
		{
			delete fineGrid;
			fineGrid = 0;
			delete coarseGrid;
			coarseGrid = 0;
			return false;
		}

		int matrixsz = 0;
		if(para.globalGridSolver == OPEN_POISSON)
		{
			if((matrixsz =  ZQ_InitMatOpenPoisson3D(coarseGrid,&coarseA)) == 0)
			{
				delete fineGrid;
				fineGrid = 0;
				delete coarseGrid;
				coarseGrid = 0;
				return false;
			}
		}
		else if(para.globalGridSolver == CLOSED_POISSON)
		{
			if((matrixsz =  ZQ_InitMatClosedPoisson3D(coarseGrid,&coarseA)) == 0)
			{
				delete fineGrid;
				fineGrid = 0;
				delete coarseGrid;
				coarseGrid = 0;
				return false;
			}
		}
		else if(para.globalGridSolver == OPEN_FLUX)
		{
			if((matrixsz = ZQ_InitMatOpenFlux3D(coarseGrid,&coarseA)) == 0)
			{
				delete fineGrid;
				fineGrid = 0;
				delete coarseGrid;
				coarseGrid = 0;
				return false;
			}
		}
		else if(para.globalGridSolver == CLOSED_FLUX)
		{
			if((matrixsz = ZQ_InitMatClosedFlux3D(coarseGrid,&coarseA)) == 0)
			{
				delete fineGrid;
				fineGrid = 0;
				delete coarseGrid;
				coarseGrid = 0;
				return false;
			}
		}
		if(para.subSize >= 2)
		{	
			if(!ZQ_InitSubMatPoisson3D(para.subSize, &submatrix,&subrow,&subcol))
			{
				delete fineGrid;
				fineGrid = 0;
				delete coarseGrid;
				coarseGrid = 0;
				return false;
			}
		}
	}
	else
	{
		if(!ZQ_InitSimulation3D(fineGrid,coarseGrid,para,&coarseA,&submatrix,&subrow,&subcol))
		{
			delete fineGrid;
			fineGrid = 0;
			delete coarseGrid;
			coarseGrid = 0;
			return false;
		}
	}

	return true;
}

bool ZQ_SmokeSimulation3D::_projection(const ZQ_SmokeGuideGrid3D* guideGrid)
{
	LARGE_INTEGER tf;
	QueryPerformanceFrequency(&tf);


	/**************************** TOTAL  Projection  **********************************************/
	/*appling guiding is not allowed */
	if(para.globalGridSolver == TOTAL_OPEN_POISSON_SOR || para.globalGridSolver == TOTAL_CLOSED_POISSON_SOR
		|| para.globalGridSolver == TOTAL_OPEN_FLUX_SOR || para.globalGridSolver == TOTAL_CLOSED_FLUX_SOR)
	{
		if(!_totalProjection())
		{
			return false;
		}

		return true;
	}
	
	/*CPU projection*/
	if(para.globalGridSolver == OPEN_FLUX || para.globalGridSolver == OPEN_POISSON
		|| para.globalGridSolver == CLOSED_FLUX || para.globalGridSolver == CLOSED_POISSON)
	{
		if(para.useMovingObject)
		{
			LARGE_INTEGER rebuild_st, rebuild_end;
			QueryPerformanceCounter(&rebuild_st);

			int matrixsz = 0;
			if(coarseA)
			{
				ZQ_TaucsBase::ZQ_taucs_ccs_free(coarseA);
				coarseA = 0;
			}

			if(!ZQ_InitCoarseGrid3D(fineGrid,coarseGrid))
				return false;
			if(para.globalGridSolver == OPEN_POISSON)
			{
				if((matrixsz =  ZQ_InitMatOpenPoisson3D(coarseGrid,&coarseA)) == 0)
					return false;
			}
			else if(para.globalGridSolver == CLOSED_POISSON)
			{
				if((matrixsz =  ZQ_InitMatClosedPoisson3D(coarseGrid,&coarseA)) == 0)
					return false;
			}
			else if(para.globalGridSolver == OPEN_FLUX)
			{
				if((matrixsz = ZQ_InitMatOpenFlux3D(coarseGrid,&coarseA)) == 0)
					return false;
			}
			else if(para.globalGridSolver == CLOSED_FLUX)
			{
				if((matrixsz = ZQ_InitMatClosedFlux3D(coarseGrid,&coarseA)) == 0)
					return false;
			}
			QueryPerformanceCounter(&rebuild_end);
			printf("rebuild cost:%f\n",1.0*(rebuild_end.QuadPart-rebuild_st.QuadPart)/tf.QuadPart);
		}
		if (!ZQ_TotalProjection3D(fineGrid,coarseGrid,para,coarseA,submatrix,subrow,subcol,guideGrid))
			return false;

		return true;
	}
	
	/*GPU projection*/
	if(para.globalGridSolver == COARSE_OPEN_POISSON_SOR || para.globalGridSolver == COARSE_CLOSED_POISSON_SOR 
		|| para.globalGridSolver == COARSE_OPEN_FLUX_SOR || para.globalGridSolver == COARSE_CLOSED_FLUX_SOR)
	{
		if(para.useMovingObject)
		{
			if(!ZQ_InitCoarseGrid3D(fineGrid,coarseGrid))
				return false;
		}
		if (!ZQ_TotalProjection3D(fineGrid, coarseGrid, para, coarseA, submatrix, subrow, subcol, guideGrid))
			return false;

		return true;
	}
	

	if(para.globalGridSolver == OPEN_OCTREE_POISSON || para.globalGridSolver == OPEN_OCTREE_FLUX
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON || para.globalGridSolver == CLOSED_OCTREE_FLUX
		|| para.globalGridSolver == OPEN_OCTREE_POISSON_SOR || para.globalGridSolver == CLOSED_OCTREE_POISSON_SOR)
	{
		if(para.useMovingObject)
		{
			LARGE_INTEGER rebuild_st, rebuild_end;
			QueryPerformanceCounter(&rebuild_st);

			if(!octreeGrid->ReBuildFromGlobalGrid(fineGrid,para,0/*interestingRegion*/))
				return false;
			QueryPerformanceCounter(&rebuild_end);
			printf("rebuild cost:%f\n",1.0*(rebuild_end.QuadPart-rebuild_st.QuadPart)/tf.QuadPart);
		}

		LARGE_INTEGER solve_global_st, solve_global_end;
		QueryPerformanceCounter(&solve_global_st);
		if(!octreeGrid->GlobalProjection(fineGrid,para))
			return false;
		QueryPerformanceCounter(&solve_global_end);
		printf("solve global cost:%f\n", 1.0*(solve_global_end.QuadPart-solve_global_st.QuadPart)/tf.QuadPart);

		if(para.useGuidingControl)
		{
			LARGE_INTEGER app_guide_st, app_guide_end;
			QueryPerformanceCounter(&app_guide_st);
			if(!octreeGrid->ApplyGuidingVelocity(fineGrid,para,guideGrid))
				return false;
			QueryPerformanceCounter(&app_guide_end);
			printf("app guide cost:%f\n",1.0*(app_guide_end.QuadPart-app_guide_st.QuadPart)/tf.QuadPart);
		}

		LARGE_INTEGER solve_sub_st, solve_sub_end;
		QueryPerformanceCounter(&solve_sub_st);
		if(!octreeGrid->SubProjection(fineGrid,para))
			return false;
		QueryPerformanceCounter(&solve_sub_end);
		printf("solve sub cost:%f\n",1.0*(solve_sub_end.QuadPart-solve_sub_st.QuadPart)/tf.QuadPart);

		return true;
	}

	return false;
}

bool ZQ_SmokeSimulation3D::_totalProjection()
{
	LARGE_INTEGER tf,solve_global_st, solve_global_end;
	QueryPerformanceFrequency(&tf);
	QueryPerformanceCounter(&solve_global_st);

	switch (para.globalGridSolver)
	{
	case TOTAL_OPEN_FLUX_SOR:
		ZQ_SolveOpenFlux3DSOR(fineGrid,para);
		break;
	case TOTAL_CLOSED_FLUX_SOR:
		ZQ_SolveClosedFlux3DSOR(fineGrid,para);
		break;
	case TOTAL_OPEN_POISSON_SOR:
		ZQ_SolveOpenPressure3DSOR(fineGrid,para);
		break;
	case TOTAL_CLOSED_POISSON_SOR:
		ZQ_SolveClosedPressure3DSOR(fineGrid,para);
		break;
	default:
		return false;
	}
	QueryPerformanceCounter(&solve_global_end);
	printf("solve global cost:%f\n",1.0*(solve_global_end.QuadPart-solve_global_st.QuadPart)/tf.QuadPart);
	return true;
}