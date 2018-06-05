#include "ZQ_SmokeOctreeGrid2D.h"
#include <math.h>
#include "ZQ_SparseMatrix.h"
#include "ZQ_PCGSolver.h"
#include "ZQ_PoissonSolver.h"
#include "ZQ_TaucsBase.h"

ZQ_SmokeOctreeGrid2D::ZQ_SmokeOctreeGrid2D()
{
	//uCount = 0;
	//vCount = 0;
	//lamdaCount = 0;
	
	dim = 0;
	deltax = 0;
	deltaU = 0;
	deltaV = 0;
	Uidx2 = 0;
	Vidx2 = 0;
	Uidx4 = 0;
	Vidx4 = 0;
	Uidx8 = 0;
	Vidx8 = 0;

	leafLevel0 = 0;
	leafLevel1 = 0;
	leafLevel2 = 0;
	leafLevel3 = 0;

	lamda_index[0] = 0;
	lamda_index[1] = 0;
	lamda_index[2] = 0;
	lamda_index[3] = 0;
	
	A = 0;

	subMatrix2 = 0;
	subMatrix4 = 0;
	subMatrix8 = 0;

	level_dis_thresh1 = 0;
	level_dis_thresh2 = 0;
	level_dis_thresh3 = 0;

	specialRegion_mask = 0;
	specialRegion_label = 0;
}

ZQ_SmokeOctreeGrid2D::~ZQ_SmokeOctreeGrid2D()
{
	if(deltax)
	{
		delete []deltax;
		deltax = 0;
	}
	if(deltaU)
	{
		delete []deltaU;
		deltaU = 0;
	}
	if(deltaV)
	{
		delete []deltaV;
		deltaV = 0;
	}
	if(Uidx2)
	{
		delete []Uidx2;
		Uidx2 = 0;
	}
	if(Uidx4)
	{
		delete []Uidx4;
		Uidx4 = 0;
	}
	if(Uidx8)
	{
		delete []Uidx8;
		Uidx8 = 0;
	}
	if(Vidx2)
	{
		delete []Vidx2;
		Vidx2 = 0;
	}
	if(Vidx4)
	{
		delete []Vidx4;
		Vidx4 = 0;
	}
	if(Vidx8)
	{
		delete []Vidx8;
		Vidx8 = 0;
	}
	if(leafLevel0)
	{
		delete []leafLevel0;
		leafLevel0 = 0;
	}
	if(leafLevel1)
	{
		delete []leafLevel1;
		leafLevel1 = 0;
	}
	if(leafLevel2)
	{
		delete []leafLevel2;
		leafLevel2 = 0;
	}
	if(leafLevel3)
	{
		delete []leafLevel3;
		leafLevel3 = 0;
	}
	if(lamda_index[0])
	{
		delete []lamda_index[0];
		lamda_index[0] = 0;
	}
	if(lamda_index[1])
	{
		delete []lamda_index[1];
		lamda_index[1] = 0;
	}
	if(lamda_index[2])
	{
		delete []lamda_index[2];
		lamda_index[2] = 0;
	}
	if(lamda_index[3])
	{
		delete []lamda_index[3];
		lamda_index[3] = 0;
	}
	if(subMatrix2)
	{
		delete []subMatrix2;
		subMatrix2 = 0;
	}
	if(subMatrix4)
	{
		delete []subMatrix4;
		subMatrix4 = 0;
	}
	if(subMatrix8)
	{
		delete []subMatrix8;
		subMatrix8 = 0;
	}

	if(specialRegion_mask)
	{
		delete []specialRegion_mask;
		specialRegion_mask = 0;
	}

	if(specialRegion_label)
	{
		delete []specialRegion_label;
		specialRegion_label = 0;
	}

	if(A)
	{
		ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(A);
		A = 0;
	}
	for(int i = 0;i < specialRegion_A.size();i++)
	{
		if(specialRegion_A[i])
		{
			ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(specialRegion_A[i]);
			specialRegion_A[i] = 0;
		}
	}

}

bool ZQ_SmokeOctreeGrid2D::InitFromGlobalGrid(const ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para, const int *interestingRegionDist)
{
	if(!_buildGlobalMatrix(srcGrid,para,interestingRegionDist))
		return false;
	if(!_buildSubMatrix())
		return false;
	return true;
}

bool ZQ_SmokeOctreeGrid2D::ReBuildFromGlobalGrid(const ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para, const int *interestingRegionDist)
{
	return _buildGlobalMatrix(srcGrid,para,interestingRegionDist);
}

bool ZQ_SmokeOctreeGrid2D::GlobalProjection(ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	return _solveGlobal(srcGrid,para);
}

bool ZQ_SmokeOctreeGrid2D::ApplyGuidingVelocity(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const ZQ_SmokeGuideGrid2D* guideGrid)
{

	//IMPORTANT FLAG
	bool blend_flag = para.guidingNaiveBlend;

	if(guideGrid == 0 || srcGrid == 0)
		return false;
	int guide_width = guideGrid->GetWidth();
	int guide_height = guideGrid->GetHeight();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	if(width % guide_width != 0 || height % guide_height != 0 || width/guide_width != height/guide_height)
		return false;

	int scale = width / guide_width;
	int levelWidth1 = width/2;
	int levelHeight1 = height/2;
	int levelWidth2 = width/4;
	int levelHeight2 = height/4;
	int levelWidth3 = width/8;
	int levelHeight3 = height/8;

	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	const float*& guideUptr = guideGrid->GetUptr();
	const float*& guideVptr = guideGrid->GetVptr();
	
	float* deltaU = new float[(width+1)*height];
	memset(deltaU,0,sizeof(float)*(width+1)*height);
	float* deltaV = new float[(height+1)*width];
	memset(deltaV,0,sizeof(float)*(height+1)*width);

	int size;
	int i_shift,j_shift;
	float guidingCoeff = para.guidingCoeff;
	for(int i = 0;i < Uinfos.size();i++)
	{
		int lvl = __min(Uinfos[i].level1, Uinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Uinfos[i].level1)
		{

			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		if(blend_flag)
		{
			float val;
			for(int sj = 0;sj < size;sj++)
			{
				float coord_x = (float)i_shift/width*guide_width;
				float coord_y = (float)(j_shift+sj+0.5f)/height*guide_height-0.5f;

				ZQ::ZQ_ImageProcessing::BilinearInterpolate(guideUptr,guide_width+1,guide_height,1,coord_x,coord_y,&val,false);
				deltaU[(j_shift+sj)*(width+1)+i_shift] = guidingCoeff * (val - Uptr[(j_shift+sj)*(width+1)+i_shift]);
			}
		}
		else
		{
			float val;
			float total_coarse = 0;
			float total_fine = 0;
			for(int sj = 0;sj < size;sj++)
			{
				float coord_x = (float)i_shift/width*guide_width;
				float coord_y = (float)(j_shift+sj+0.5f)/height*guide_height-0.5f;

				ZQ::ZQ_ImageProcessing::BilinearInterpolate(guideUptr,guide_width+1,guide_height,1,coord_x,coord_y,&val,false);
				total_coarse += val;
				total_fine += Uptr[(j_shift+sj)*(width+1)+i_shift];
			}
			float delta_tmp = guidingCoeff * (total_coarse - total_fine) / size;
			for(int sj = 0;sj < size;sj++)
			{
				deltaU[(j_shift+sj)*(width+1)+i_shift] = delta_tmp;
			}
		}
	}

	for(int i = 0;i < Vinfos.size();i++)
	{
		int lvl = __min(Vinfos[i].level1, Vinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		if(blend_flag)
		{
			float val;
			for(int si = 0;si < size;si++)
			{
				float coord_x = (i_shift+si+0.5f)/width*guide_width-0.5f;
				float coord_y = (float)j_shift/height*guide_height;
				ZQ::ZQ_ImageProcessing::BilinearInterpolate(guideVptr,guide_width,guide_height+1,1,coord_x,coord_y,&val,false);
				deltaV[j_shift*width+i_shift+si] = guidingCoeff * (val-Vptr[j_shift*width+i_shift+si]);
			}
		}
		else
		{
			float val;
			float total_coarse = 0;
			float total_fine = 0;
			for(int si = 0;si < size;si++)
			{
				float coord_x = (i_shift+si+0.5f)/width*guide_width-0.5f;
				float coord_y = (float)j_shift/height*guide_height;
				ZQ::ZQ_ImageProcessing::BilinearInterpolate(guideVptr,guide_width,guide_height+1,1,coord_x,coord_y,&val,false);
				total_coarse += val;
				total_fine += Vptr[j_shift*width+i_shift+si];
			}

			float delta_tmp = guidingCoeff * (total_coarse - total_fine) / size;
			for(int si = 0;si < size;si++)
			{
				deltaV[j_shift*width+i_shift+si] = delta_tmp;
			}
		}
	}

	bool*& occupyPtr = srcGrid->GetOccupyPtr();

	for(int j = 0;j < height;j++)
	{
		for(int i = 1;i < width;i++)
		{
			if(specialRegion_mask[j*width+i-1] && specialRegion_mask[j*width+i])
				deltaU[j*(width+1)+i] = 0;
		}
	}
	for(int i = 0;i < width;i++)
	{
		for(int j = 1;j < height;j++)
		{
			if(specialRegion_mask[(j-1)*width+i] && specialRegion_mask[j*width+i])
				deltaV[j*width+i] = 0;
		}
	}
	

	if(guiding_specialRegion_useCUDA) //means using CUDA
	{
		for(int iregion = 0;iregion < specialRegion_num;iregion++)
		{
			bool* tmp_occupy = new bool[width*height];
			int first_x = -1, first_y = -1;
			float div_per_volume = 0.0f;
			int div_count = 0;
			int min_x = width;
			int min_y = height;
			int max_x = -1;
			int max_y = -1;
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					tmp_occupy[j*width+i] = (specialRegion_label[j*width+i] != iregion+1);
					if(!tmp_occupy[j*width+i])
					{
						if(first_x < 0 && first_y < 0)
						{
							first_x = i;
							first_y = j;
						}

						div_count ++;
						div_per_volume += Uptr[j*(width+1)+i+1] - Uptr[j*(width+1)+i] + Vptr[(j+1)*width+i] - Vptr[j*width+i];

						min_x = __min(min_x,i);
						min_y = __min(min_y,j);
						max_x = __max(max_x,i);
						max_y = __max(max_y,j);
					}
				}
			}
			div_per_volume /= div_count;

			
			int cur_maxIter = __min(para.maxIter,div_count);
			if(false)
			{
				ZQ_CUDA_PoissonSolver2D::SolveClosedPoissonRedBlackwithOccupy2D_MAC(Uptr,Vptr,tmp_occupy,/*first_x,first_y,*/div_per_volume,width,height,cur_maxIter);
			}
			else
			{
				int offset_x = min_x;
				int offset_y = min_y;
				int size_x = max_x - min_x + 1;
				int size_y = max_y - min_y + 1;

				bool* cur_occupy = new bool[size_x*size_y];
				float* cur_mac_u = new float[(size_x+1)*size_y];
				float* cur_mac_v = new float[size_x*(size_y+1)];
				for(int i = 0;i < size_y;i++)
				{
					memcpy(cur_occupy+i*size_x,tmp_occupy+(offset_y+i)*width+offset_x,sizeof(bool)*size_x);
					memcpy(cur_mac_u+i*(size_x+1),Uptr+(offset_y+i)*(width+1)+offset_x,sizeof(float)*(size_x+1));
					memcpy(cur_mac_v+i*size_x,Vptr+(offset_y+i)*width+offset_x,sizeof(float)*size_x);
				}
				memcpy(cur_mac_v+size_y*size_x,Vptr+(offset_y+size_y)*width+offset_x,sizeof(float)*size_x);

				first_x -= offset_x;
				first_y -= offset_y;
				ZQ_CUDA_PoissonSolver2D::SolveClosedPoissonRedBlackwithOccupy2D_MAC(cur_mac_u,cur_mac_v,cur_occupy,/*first_x,first_y,*/div_per_volume,size_x,size_y,cur_maxIter);

				for(int i = 0;i < size_y;i++)
				{
					memcpy(Uptr+(offset_y+i)*(width+1)+offset_x,cur_mac_u+i*(size_x+1),sizeof(float)*(size_x+1));
					memcpy(Vptr+(offset_y+i)*width+offset_x,cur_mac_v+i*size_x,sizeof(float)*size_x);
				}
				memcpy(Vptr+(offset_y+size_y)*width+offset_x,cur_mac_v+size_y*size_x,sizeof(float)*size_x);
				delete []cur_mac_u;
				delete []cur_mac_v;
				delete []cur_occupy;
			}
			
			

			
			delete []tmp_occupy;
		}
	}
	else
	{
		for(int iregion = 0;iregion < specialRegion_A.size();iregion++)
		{
			if(specialRegion_A[iregion])
			{
				double* divergence = new double[specialRegion_A[iregion]->m];
				double* x = new double[specialRegion_A[iregion]->n];
				double* x0 = new double[specialRegion_A[iregion]->n];
				memset(x0,0,sizeof(double)*specialRegion_A[iregion]->n);
				int idx = 0;
				for(int j = 0;j < height;j++)
				{
					for(int i = 0;i < width;i++)
					{
						if(specialRegion_label[j*width+i] == iregion+1 && !occupyPtr[j*width+i])
						{
							divergence[idx++] = deltaU[j*(width+1)+i+1] - deltaU[j*(width+1)+i]
							+ deltaV[(j+1)*width+i] - deltaV[j*width+i];
						}
					}
				}
				ZQ::ZQ_PCGSolver solver;
				double tol = 1e-9;
				int max_iter = para.maxIter;
				int it = 0;
				solver.PCG_sparse_unsquare(specialRegion_A[iregion],divergence,x0,max_iter,tol,x,it,false);

				float* pressure = new float[width*height];
				memset(pressure,0,sizeof(float)*width*height);
				idx = 0;
				for(int j = 0;j < height;j++)
				{
					for(int i = 0;i < width;i++)
					{
						if(specialRegion_label[j*width+i] == iregion+1 && !occupyPtr[j*width+i])
						{
							if(idx == 0)
								idx ++;
							else
							{
								pressure[j*width+i] = x[idx-1];
								idx++;
							}
						}
					}
				}
				for(int j = 0;j < height;j++)
				{
					for(int i = 1;i < width;i++)
					{
						if(specialRegion_label[j*width+(i-1)] == iregion+1 && specialRegion_label[j*width+i] == iregion+1)
							deltaU[j*(width+1)+i] = -(pressure[j*width+i] - pressure[j*width+i-1]);
					}
				}
				for(int j = 1;j < height;j++)
				{
					for(int i = 0;i < width;i++)
					{
						if(specialRegion_label[(j-1)*width+i] == iregion+1 && specialRegion_label[j*width+i] == iregion+1)
							deltaV[j*width+i] = -(pressure[j*width+i] - pressure[(j-1)*width+i]);
					}
				}

				delete []pressure;
				delete []x0;
				delete []x;
				delete []divergence;
			}
		}

		for(int i = 0;i < height*(width+1);i++)
			Uptr[i] += deltaU[i];
		for(int i = 0;i < width*(height+1);i++)
			Vptr[i] += deltaV[i];
		delete []deltaU;
		delete []deltaV;
		deltaU = 0;
		deltaV = 0;

	}
	

	return true;
}

bool ZQ_SmokeOctreeGrid2D::SubProjection(ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(srcGrid == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	
	int size,i_shift,j_shift;
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	
	//************** sub projection ************//
	int count1 = 0;
	int count2 = 0;
	int count3 = 0;
	for(int i = 0;i < lamdainfos.size();i++)
	{
		if(lamdainfos[i].level1 == 1)
			count1 ++;
		else if(lamdainfos[i].level1 == 2)
			count2 ++;
		else if(lamdainfos[i].level1 == 3)
			count3 ++;
	}

	int inputRow1 = 2*2*(2+1);
	int inputCol1 = count1;
	int outputRow1 = 2*2*(2-1);
	int outputCol1 = count1;
	int inputRow2 = 2*4*(4+1);
	int inputCol2 = count2;
	int outputRow2 = 2*4*(4-1);
	int outputCol2 = count2;
	int inputRow3 = 2*8*(8+1);
	int inputCol3 = count3;
	int outputRow3 = 2*8*(8-1);
	int outputCol3 = count3;
	float* inputVelocity1 = 0;
	float* inputVelocity2 = 0;
	float* inputVelocity3 = 0;
	float* outputVelocity1 = 0;
	float* outputVelocity2 = 0;
	float* outputVelocity3 = 0;
	if(count1 > 0)
	{
		inputVelocity1 = new float[inputRow1*inputCol1];
		outputVelocity1 = new float[outputRow1*outputCol1];
		memset(inputVelocity1,0,sizeof(float)*inputRow1*inputCol1);
		memset(outputVelocity1,0,sizeof(float)*outputRow1*outputCol1);

	}
	if(count2 > 0)
	{
		inputVelocity2 = new float[inputRow2*inputCol2];
		outputVelocity2 = new float[outputRow2*outputCol2];
		memset(inputVelocity2,0,sizeof(float)*inputRow2*inputCol2);
		memset(outputVelocity2,0,sizeof(float)*outputRow2*outputCol2);

	}
	
	if(count3 > 0)
	{
		inputVelocity3 = new float[inputRow3*inputCol3];
		outputVelocity3 = new float[outputRow3*outputCol3];
		memset(inputVelocity3,0,sizeof(float)*inputRow3*inputCol3);
		memset(outputVelocity3,0,sizeof(float)*outputRow3*outputCol3);

	}
	
	int col_idx1 = 0;
	int col_idx2 = 0;
	int col_idx3 = 0;


	for(int i = 0;i < lamdainfos.size();i++)
	{
		if(lamdainfos[i].level1 == 1)
		{
			size = 2;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			int row_idx = 0;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 0;si <= size;si++)
				{
					inputVelocity1[row_idx*inputCol1+col_idx1] = Uptr[(j_shift+sj)*(width+1)+(i_shift+si)];
					row_idx ++;
				}
			}
			for(int sj = 0;sj <= size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					inputVelocity1[row_idx*inputCol1+col_idx1] = Vptr[(j_shift+sj)*width+(i_shift+si)];
					row_idx ++;
				}
			}
			col_idx1 ++;
		}
		else if(lamdainfos[i].level1 == 2)
		{
			size = 4;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			int row_idx = 0;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 0;si <= size;si++)
				{
					inputVelocity2[row_idx*inputCol2+col_idx2] = Uptr[(j_shift+sj)*(width+1)+(i_shift+si)];
					row_idx ++;
				}
			}
			for(int sj = 0;sj <= size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					inputVelocity2[row_idx*inputCol2+col_idx2] = Vptr[(j_shift+sj)*width+(i_shift+si)];
					row_idx ++;
				}
			}
			col_idx2 ++;
		}
		else if(lamdainfos[i].level1 == 3)
		{
			size = 8;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			int row_idx = 0;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 0;si <= size;si++)
				{
					inputVelocity3[row_idx*inputCol3+col_idx3] = Uptr[(j_shift+sj)*(width+1)+(i_shift+si)];
					row_idx ++;
				}
			}
			for(int sj = 0;sj <= size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					inputVelocity3[row_idx*inputCol3+col_idx3] = Vptr[(j_shift+sj)*width+(i_shift+si)];
					row_idx ++;
				}
			}
			col_idx3 ++;
		}
	}

	bool returnflag = true;
	if(count1 > 0)
	{
		SubProjectionCuda(outputRow1,inputRow1,inputCol1,subMatrix2, inputVelocity1, outputVelocity1);
	}
	if(count2 > 0)
	{
		SubProjectionCuda(outputRow2,inputRow2,inputCol2,subMatrix4, inputVelocity2, outputVelocity2);
	}

	if(count3 > 0)
	{
		SubProjectionCuda(outputRow3,inputRow3,inputCol3,subMatrix8, inputVelocity3, outputVelocity3);
	}

	col_idx1 = 0;
	col_idx2 = 0;
	col_idx3 = 0;

	for(int i = 0;i < lamdainfos.size();i++)
	{
		if(lamdainfos[i].level1 == 1)
		{
			size = 2;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			int row_idx = 0;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 1;si < size;si++)
				{
					Uptr[(j_shift+sj)*(width+1)+(i_shift+si)] = outputVelocity1[row_idx*outputCol1+col_idx1];
					row_idx ++;
				}
			}
			for(int sj = 1;sj < size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					Vptr[(j_shift+sj)*width+(i_shift+si)] = outputVelocity1[row_idx*outputCol1+col_idx1];
					row_idx ++;
				}
			}
			col_idx1 ++;
		}
		else if(lamdainfos[i].level1 == 2)
		{
			size = 4;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			int row_idx = 0;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 1;si < size;si++)
				{
					Uptr[(j_shift+sj)*(width+1)+(i_shift+si)] = outputVelocity2[row_idx*outputCol2+col_idx2];
					row_idx ++;
				}
			}
			for(int sj = 1;sj < size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					Vptr[(j_shift+sj)*width+(i_shift+si)] = outputVelocity2[row_idx*outputCol2+col_idx2];
					row_idx ++;
				}
			}
			col_idx2 ++;
		}
		else if(lamdainfos[i].level1 == 3)
		{
			size = 8;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			int row_idx = 0;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 1;si < size;si++)
				{
					Uptr[(j_shift+sj)*(width+1)+(i_shift+si)] = outputVelocity3[row_idx*outputCol3+col_idx3];
					row_idx ++;
				}
			}
			for(int sj = 1;sj < size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					Vptr[(j_shift+sj)*width+(i_shift+si)] = outputVelocity3[row_idx*outputCol3+col_idx3];
					row_idx ++;
				}
			}
			col_idx3 ++;
		}
	}

	delete []inputVelocity1;
	delete []inputVelocity2;
	delete []inputVelocity3;
	delete []outputVelocity1;
	delete []outputVelocity2;
	delete []outputVelocity3;
	return returnflag;
}

bool ZQ_SmokeOctreeGrid2D::_buildGlobalMatrix(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist)
{
	switch(para.globalGridSolver)
	{
	case OPEN_OCTREE_FLUX:
		return _buildGlobalMatrixOpenFlux(srcGrid,para,interestingRegionDist);
		break;
	case OPEN_OCTREE_POISSON:
		return _buildGlobalMatrixOpenPoisson(srcGrid,para,interestingRegionDist);
		break;
	case CLOSED_OCTREE_FLUX:
		return _buildGlobalMatrixClosedFlux(srcGrid,para,interestingRegionDist);
		break;
	case CLOSED_OCTREE_POISSON:
		return _buildGlobalMatrixClosedPoisson(srcGrid,para,interestingRegionDist);
		break;
	case OPEN_OCTREE_POISSON_SOR:
		return _buildGlobalMatrixOpenPoissonSOR(srcGrid,para,interestingRegionDist);
		break;
	case CLOSED_OCTREE_POISSON_SOR:
		return _buildGlobalMatrixClosedPoissonSOR(srcGrid,para,interestingRegionDist);
		break;
	}
	return false;
}

bool ZQ_SmokeOctreeGrid2D::_solveGlobal(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	switch(para.globalGridSolver)
	{
	case OPEN_OCTREE_FLUX:
		return _solveGlobalOpenFlux(srcGrid,para);
		break;
	case OPEN_OCTREE_POISSON:
		return _solveGlobalOpenPoisson(srcGrid,para);
		break;
	case CLOSED_OCTREE_FLUX:
		return _solveGlobalClosedFlux(srcGrid,para);
		break;
	case CLOSED_OCTREE_POISSON:
		return _solveGlobalClosedPoisson(srcGrid,para);
		break;
	case OPEN_OCTREE_POISSON_SOR:
		return _solveGlobalOpenPoissonSOR(srcGrid,para);
		break;
	case CLOSED_OCTREE_POISSON_SOR:
		return _solveGlobalClosedPoissonSOR(srcGrid,para);
		break;
	}

	return false;
}

bool ZQ_SmokeOctreeGrid2D::_buildGlobalMatrixOpenPoissonSOR(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	if(para.useGuidingControl)
		_selectInfosForCPUSolvers(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,lamdainfos);

	switch(cuda_openoctree_poisson_type)
	{
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE1:
		_selectInfos_lambda(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos);
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE2:
		_selectInfos_for_CUDAOpenPoissonRedBlack2(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos,
			level0_index_red,level0_index_black,
			level1_index_red,level1_index_black,
			level2_index_red,level2_index_black,
			level3_index_red,level3_index_black);
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		_selectInfos_for_CUDAOpenPoissonRedBlack3(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos,
			level0_index_red,level0_neighborinfo_red,level0_index_black,level0_neighborinfo_black,
			level1_index_red,level1_neighborinfo_red,level1_index_black,level1_neighborinfo_black,
			level2_index_red,level2_neighborinfo_red,level2_index_black,level2_neighborinfo_black,
			level3_index_red,level3_neighborinfo_red,level3_index_black,level3_neighborinfo_black);
		break;
	}
	
	return true;
}

bool ZQ_SmokeOctreeGrid2D::_buildGlobalMatrixClosedPoissonSOR(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	if(para.useGuidingControl)
		_selectInfosForCPUSolvers(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,lamdainfos);

	switch(cuda_openoctree_poisson_type)
	{
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE1:
		_selectInfos_lambda(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos);
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE2:
		_selectInfos_for_CUDAOpenPoissonRedBlack2(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos,
			level0_index_red,level0_index_black,
			level1_index_red,level1_index_black,
			level2_index_red,level2_index_black,
			level3_index_red,level3_index_black);
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		_selectInfos_for_CUDAOpenPoissonRedBlack3(para.width,para.height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos,
			level0_index_red,level0_neighborinfo_red,level0_index_black,level0_neighborinfo_black,
			level1_index_red,level1_neighborinfo_red,level1_index_black,level1_neighborinfo_black,
			level2_index_red,level2_neighborinfo_red,level2_index_black,level2_neighborinfo_black,
			level3_index_red,level3_neighborinfo_red,level3_index_black,level3_neighborinfo_black);
		break;
	}

	return true;
}

bool ZQ_SmokeOctreeGrid2D::_buildGlobalMatrixOpenPoisson(const ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para, const int *interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();

	_selectInfosForCPUSolvers(width,height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,lamdainfos);

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lamdaCount = lamdainfos.size();

	int cur_idx = 0;
	for(int i = 0;i < lamdaCount;i++)
	{
		lamdainfos[i].index = cur_idx++;
	}

	for(int i = 0;i < 4;i++)
	{
		if(lamda_index[i])
			delete []lamda_index[i];
		lamda_index[i] = 0;
	}
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	lamda_index[0] = new int[width*height];
	lamda_index[1] = new int[levelWidth1*levelHeight1];
	lamda_index[2] = new int[levelWidth2*levelHeight2];
	lamda_index[3] = new int[levelWidth3*levelHeight3];
	for(int i = 0;i < width*height;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3;i++)
	{
		lamda_index[3][i] = -1;
	}

	for(int i = 0;i < lamdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;

		switch(lvl)
		{
		case 0:
			lamda_index[lvl][jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}
	dim = lamdaCount;

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	double area[4] = {1,2,4,8};
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
				len1 = 4;
			}
			else
			{
				col1 = -1;
				len1 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Uinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
				len2 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(col1 >= 0 && col2 >= 0)
		{
			double val = area[lvl]/(len1+len2);
			mat.AddTo(col1,col2,val);
			mat.AddTo(col2,col2,-val);
			mat.AddTo(col2,col1,val);
			mat.AddTo(col1,col1,-val);
		}
		else if(col1 < 0 && col2 >= 0)
		{
			double val = area[lvl]/(2*len2);
			mat.AddTo(col2,col2,-val);
		}
		else if(col1 >= 0 && col2 < 0)
		{
			double val = area[lvl]/(2*len1);
			mat.AddTo(col1,col1,-val);
		}
		else
		{
			printf("error in build octree OPEN_POISSON\n");
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1,Vinfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Vinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
				len1 = 4;
			}
			else
			{
				col1 = -1;
				len1 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Vinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
				len2 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		if(col1 >= 0 && col2 >= 0)
		{
			double val = area[lvl]/(len1+len2);
			mat.AddTo(col1,col2,val);
			mat.AddTo(col2,col2,-val);
			mat.AddTo(col2,col1,val);
			mat.AddTo(col1,col1,-val);
		}
		else if(col1 < 0 && col2 >= 0)
		{
			double val = area[lvl]/(2*len2);
			mat.AddTo(col2,col2,-val);
		}
		else if(col1 >= 0 && col2 < 0)
		{
			double val = area[lvl]/(2*len1);
			mat.AddTo(col1,col1,-val);
		}
		else
		{
			printf("error in build octree OPEN_POISSON\n");
		}
	}


	if(A)
	{
		ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(A);
		A = 0;
	}
	A = mat.ExportCCS(TAUCS_DOUBLE);

	if(A == 0)
	{
		printf("init global matrix fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}

	return true;
}

bool ZQ_SmokeOctreeGrid2D::_buildGlobalMatrixClosedPoisson(const ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para, const int *interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();

	_selectInfosForCPUSolvers(width,height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,lamdainfos);

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lamdaCount = lamdainfos.size();

	int cur_idx = 0;
	for(int i = 0;i < lamdaCount;i++)
	{
		lamdainfos[i].index = cur_idx++;
	}

	for(int i = 0;i < 4;i++)
	{
		if(lamda_index[i])
			delete []lamda_index[i];
		lamda_index[i] = 0;
	}
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;

	lamda_index[0] = new int[width*height];
	lamda_index[1] = new int[levelWidth1*levelHeight1];
	lamda_index[2] = new int[levelWidth2*levelHeight2];
	lamda_index[3] = new int[levelWidth3*levelHeight3];
	for(int i = 0;i < width*height;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3;i++)
	{
		lamda_index[3][i] = -1;
	}

	//the last preesure is set 0 
	for(int i = 0;i < lamdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;

		switch(lvl)
		{
		case 0:
			lamda_index[lvl][jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}
	dim = lamdaCount;
	int row_num = dim+1;

	ZQ::ZQ_SparseMatrix<float> mat(row_num,dim);
	mat.SetValue(row_num-1,dim-1,1);

	double area[4] = {1,2,4,8};

	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
				len1 = 4;
			}
			else
			{
				col1 = -1;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Uinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(col1 >= 0 && col2 >= 0)
		{
			double val = area[lvl]/(len1+len2);///volume[lvl];
			//if(col2 != lamdaCount-1)
			{
				mat.AddTo(col1,col2,val);
				mat.AddTo(col2,col2,-val);
			}
			//if(col1 != lamdaCount-1)
			{
				mat.AddTo(col2,col1,val);
				mat.AddTo(col1,col1,-val);
			}
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1,Vinfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Vinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
				len1 = 4;
			}
			else
				col1 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Vinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		if(col1 >= 0 && col2 >= 0)
		{
			double val = area[lvl]/(len1+len2);///volume[lvl];
			//if(col2 != lamdaCount-1)
			{
				mat.AddTo(col2,col1,val);
				mat.AddTo(col2,col2,-val);
			}
			//if(col1 != lamdaCount-1)
			{
				mat.AddTo(col1,col2,val);
				mat.AddTo(col1,col1,-val);
			}
		}
	}

	if(A)
	{
		ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(A);
		A = 0;
	}

	A = mat.ExportCCS(TAUCS_DOUBLE);

	if(A == 0)
	{
		printf("init global matrix fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}

	return true;
}

bool ZQ_SmokeOctreeGrid2D::_buildGlobalMatrixOpenFlux(const ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para, const int *interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();

	_selectInfosForCPUSolvers(width,height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,lamdainfos);
	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lamdaCount = lamdainfos.size();

	int cur_idx = 0;
	for(int i = 0;i < uCount;i++)
	{
		Uinfos[i].index = cur_idx++;
	}

	for(int i = 0;i < vCount;i++)
	{
		Vinfos[i].index = cur_idx++;
	}
	for(int i = 0;i < lamdaCount;i++)
	{
		lamdainfos[i].index = cur_idx++;
	}
	dim = cur_idx;

	for(int i = 0;i < 4;i++)
	{
		if(lamda_index[i])
			delete []lamda_index[i];
		lamda_index[i] = 0;
	}
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;


	for(int i = 0;i < 4;i++)
	{
		if(lamda_index[i])
		{
			delete []lamda_index[i];
			lamda_index[i] = 0;
		}
	}
	lamda_index[0] = new int[width*height];
	lamda_index[1] = new int[levelWidth1*levelHeight1];
	lamda_index[2] = new int[levelWidth2*levelHeight2];
	lamda_index[3] = new int[levelWidth3*levelHeight3];
	for(int i = 0;i < width*height;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3;i++)
	{
		lamda_index[3][i] = -1;
	}

	for(int i = 0;i < lamdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;
		switch(lvl)
		{
		case 0:
			lamda_index[lvl][jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	double area[4] = {1,2,4,8};
	double weight[4] = {1,2,4,8};
	//double volume[4] = {1,4,16,64};
	//double weight[4] = {1,1,1,1};
	double volume[4] = {1,1,1,1};
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		int row = Uinfos[i].index;
		int col ;
		float wt1 = 1, wt2 = 1;

		double val;
		switch(Uinfos[i].level1)
		{
		case 0:
			col = lamda_index[0][Uinfos[i].j1*width+Uinfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col = lamda_index[3][Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
				val = -area[lvl]/volume[3];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			else 
				wt1 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		
		switch(Uinfos[i].level2)
		{
		case 0:
			col = lamda_index[0][Uinfos[i].j2*width+Uinfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col = lamda_index[3][Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
				val = area[lvl]/volume[3];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			else
				wt2 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		col = Uinfos[i].index;
		mat.SetValue(row,col,2*weight[lvl]);
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1,Vinfos[i].level2);
		int row = Vinfos[i].index;
		int col = Vinfos[i].index;
		float wt1 = 1, wt2 = 1;
		
		double val;
		switch(Vinfos[i].level1)
		{
		case 0:
			col = lamda_index[0][Vinfos[i].j1*width+Vinfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col = lamda_index[3][Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
				val = -area[lvl]/volume[3];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			else
				wt1 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Vinfos[i].level2)
		{
		case 0:
			col = lamda_index[0][Vinfos[i].j2*width+Vinfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j2 < levelWidth3)
			{
				col = lamda_index[3][Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
				val = area[lvl]/volume[3];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			else
				wt2 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		col = Vinfos[i].index;
		mat.SetValue(row,col,2*weight[lvl]);
	}

	if(A)
	{
		ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(A);
		A = 0;
	}
	
	A = mat.ExportCCS(TAUCS_DOUBLE);
	
	if(A == 0)
	{
		printf("init global matrix fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}

	return true;
}

bool ZQ_SmokeOctreeGrid2D::_buildGlobalMatrixClosedFlux(const ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para, const int* interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();

	_selectInfosForCPUSolvers(width,height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,lamdainfos);

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lambdaCount = lamdainfos.size();

	int cur_idx = 0;
	for(int i = 0;i < uCount;i++)
	{
		if((Uinfos[i].level1 == 3 && Uinfos[i].i1 < 0) || (Uinfos[i].level2 == 3 && Uinfos[i].i2 >= width/8))
			Uinfos[i].index = -1;
		else
			Uinfos[i].index = cur_idx++;
	}

	for(int i = 0;i < vCount;i++)
	{
		if((Vinfos[i].level1 == 3 && Vinfos[i].j1 < 0) || (Vinfos[i].level2 == 3 && Vinfos[i].j2 >= height/8))
			Vinfos[i].index = -1;
		else
			Vinfos[i].index = cur_idx++;
	}
	for(int i = 0;i < lambdaCount;i++)
	{
		lamdainfos[i].index = cur_idx++;
	}
	k_index = cur_idx++;
	dim = cur_idx;


	for(int i = 0;i < 4;i++)
	{
		if(lamda_index[i])
			delete []lamda_index[i];
		lamda_index[i] = 0;
	}
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;

	lamda_index[0] = new int[width*height];
	lamda_index[1] = new int[levelWidth1*levelHeight1];
	lamda_index[2] = new int[levelWidth2*levelHeight2];
	lamda_index[3] = new int[levelWidth3*levelHeight3];
	for(int i = 0;i < width*height;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3;i++)
	{
		lamda_index[3][i] = -1;
	}

	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;
		switch(lvl)
		{
		case 0:
			lamda_index[lvl][jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	double area[4] = {1,2,4,8};
	double weight[4] = {1,2,4,8};
	double volume[4] = {1,4,16,64};
	
	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int row = lamdainfos[i].index;
		int col = k_index;
		mat.SetValue(row,col,volume[lvl]);
		mat.SetValue(col,row,volume[lvl]);
	}
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		int row = Uinfos[i].index;
		int col = Uinfos[i].index;
		if(row < 0)
			continue;
		mat.SetValue(row,col,2*weight[lvl]);

		double val;
		switch(Uinfos[i].level1)
		{
		case 0:
			col = lamda_index[0][Uinfos[i].j1*width+Uinfos[i].i1];
			val = -area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			val = -area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			val = -area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col = lamda_index[3][Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
				val = -area[lvl];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Uinfos[i].level2)
		{
		case 0:
			col = lamda_index[0][Uinfos[i].j2*width+Uinfos[i].i2];
			val = area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			val = area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			val = area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col = lamda_index[3][Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
				val = area[lvl];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1,Vinfos[i].level2);
		int row = Vinfos[i].index;
		int col = Vinfos[i].index;
		if(row < 0)
			continue;
		mat.SetValue(row,col,2*weight[lvl]);
		double val;
		switch(Vinfos[i].level1)
		{
		case 0:
			col = lamda_index[0][Vinfos[i].j1*width+Vinfos[i].i1];
			val = -area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			val = -area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			val = -area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col = lamda_index[3][Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
				val = -area[lvl];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Vinfos[i].level2)
		{
		case 0:
			col = lamda_index[0][Vinfos[i].j2*width+Vinfos[i].i2];
			val = area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			val = area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			val = area[lvl];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col = lamda_index[3][Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
				val = area[lvl];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}

	if(A)
	{
		ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(A);
		A = 0;
	}

	A = mat.ExportCCS(TAUCS_DOUBLE);

	if(A == 0)
	{
		printf("init global matrix fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}
	return true;
}




bool ZQ_SmokeOctreeGrid2D::_solveGlobalOpenPoissonSOR(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(srcGrid == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();

	switch(cuda_openoctree_poisson_type)
	{
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE1:
		ZQ_CUDA_PoissonSolver2D::SolveOpenOctreePoissonRedBlack2D_MAC(Uptr,Vptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,width,height,para.maxIter);
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE2:
		{
			int level0_num_red = level0_index_red.size()/2;
			int level1_num_red = level1_index_red.size()/2;
			int level2_num_red = level2_index_red.size()/2;
			int level3_num_red = level3_index_red.size()/2;
			int level0_num_black = level0_index_black.size()/2;
			int level1_num_black = level1_index_black.size()/2;
			int level2_num_black = level2_index_black.size()/2;
			int level3_num_black = level3_index_black.size()/2;
			int* level0_index_red_ptr = level0_num_red > 0 ? &level0_index_red[0] : 0;
			int* level1_index_red_ptr = level1_num_red > 0 ? &level1_index_red[0] : 0;
			int* level2_index_red_ptr = level2_num_red > 0 ? &level2_index_red[0] : 0;
			int* level3_index_red_ptr = level3_num_red > 0 ? &level3_index_red[0] : 0;
			int* level0_index_black_ptr = level0_num_black > 0 ? &level0_index_black[0] : 0;
			int* level1_index_black_ptr = level1_num_black > 0 ? &level1_index_black[0] : 0;
			int* level2_index_black_ptr = level2_num_black > 0 ? &level2_index_black[0] : 0;
			int* level3_index_black_ptr = level3_num_black > 0 ? &level3_index_black[0] : 0;
			ZQ_CUDA_PoissonSolver2D::SolveOpenOctreePoissonRedBlack2_2D_MAC(Uptr,Vptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,width,height,para.maxIter,
				level0_num_red,level0_index_red_ptr,level0_num_black,level0_index_black_ptr,
				level1_num_red,level1_index_red_ptr,level1_num_black,level1_index_black_ptr,
				level2_num_red,level2_index_red_ptr,level2_num_black,level2_index_black_ptr,
				level3_num_red,level3_index_red_ptr,level3_num_black,level3_index_black_ptr);
		}
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		{
			const int index_channels = 4;
			const int neighborinfo_channels = 3;
			int level0_num_red = level0_index_red.size()/index_channels;
			int level1_num_red = level1_index_red.size()/index_channels;
			int level2_num_red = level2_index_red.size()/index_channels;
			int level3_num_red = level3_index_red.size()/index_channels;
			int level0_num_black = level0_index_black.size()/index_channels;
			int level1_num_black = level1_index_black.size()/index_channels;
			int level2_num_black = level2_index_black.size()/index_channels;
			int level3_num_black = level3_index_black.size()/index_channels;
			int level0_info_len_red = level0_neighborinfo_red.size()/neighborinfo_channels;
			int level1_info_len_red = level1_neighborinfo_red.size()/neighborinfo_channels;
			int level2_info_len_red = level2_neighborinfo_red.size()/neighborinfo_channels;
			int level3_info_len_red = level3_neighborinfo_red.size()/neighborinfo_channels;
			int level0_info_len_black = level0_neighborinfo_black.size()/neighborinfo_channels;
			int level1_info_len_black = level1_neighborinfo_black.size()/neighborinfo_channels;
			int level2_info_len_black = level2_neighborinfo_black.size()/neighborinfo_channels;
			int level3_info_len_black = level3_neighborinfo_black.size()/neighborinfo_channels;

			int* level0_index_red_ptr = level0_num_red > 0 ? &level0_index_red[0] : 0;
			int* level1_index_red_ptr = level1_num_red > 0 ? &level1_index_red[0] : 0;
			int* level2_index_red_ptr = level2_num_red > 0 ? &level2_index_red[0] : 0;
			int* level3_index_red_ptr = level3_num_red > 0 ? &level3_index_red[0] : 0;
			int* level0_index_black_ptr = level0_num_black > 0 ? &level0_index_black[0] : 0;
			int* level1_index_black_ptr = level1_num_black > 0 ? &level1_index_black[0] : 0;
			int* level2_index_black_ptr = level2_num_black > 0 ? &level2_index_black[0] : 0;
			int* level3_index_black_ptr = level3_num_black > 0 ? &level3_index_black[0] : 0;
			int* level0_info_red_ptr = level0_info_len_red > 0 ? &level0_neighborinfo_red[0] : 0;
			int* level1_info_red_ptr = level1_info_len_red > 0 ? &level1_neighborinfo_red[0] : 0;
			int* level2_info_red_ptr = level2_info_len_red > 0 ? &level2_neighborinfo_red[0] : 0;
			int* level3_info_red_ptr = level3_info_len_red > 0 ? &level3_neighborinfo_red[0] : 0;
			int* level0_info_black_ptr = level0_info_len_black > 0 ? &level0_neighborinfo_black[0] : 0;
			int* level1_info_black_ptr = level1_info_len_black > 0 ? &level1_neighborinfo_black[0] : 0;
			int* level2_info_black_ptr = level2_info_len_black > 0 ? &level2_neighborinfo_black[0] : 0;
			int* level3_info_black_ptr = level3_info_len_black > 0 ? &level3_neighborinfo_black[0] : 0;
			ZQ_CUDA_PoissonSolver2D::SolveOpenOctreePoissonRedBlack3_2D_MAC(Uptr,Vptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,width,height,para.maxIter,
				level0_num_red,		level0_index_red_ptr,	level0_info_len_red,	level0_info_red_ptr,
				level0_num_black,	level0_index_black_ptr,	level0_info_len_black,	level0_info_black_ptr,
				level1_num_red,		level1_index_red_ptr,	level1_info_len_red,	level1_info_red_ptr,
				level1_num_black,	level1_index_black_ptr,	level1_info_len_black,	level1_info_black_ptr,
				level2_num_red,		level2_index_red_ptr,	level2_info_len_red,	level2_info_red_ptr,
				level2_num_black,	level2_index_black_ptr,	level2_info_len_black,	level2_info_black_ptr,
				level3_num_red,		level3_index_red_ptr,	level3_info_len_red,	level3_info_red_ptr,
				level3_num_black,	level3_index_black_ptr,	level3_info_len_black,	level3_info_black_ptr);
		}
		break;
	}

	return true;
}


bool ZQ_SmokeOctreeGrid2D::_solveGlobalClosedPoissonSOR(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(srcGrid == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	bool*& occupy = srcGrid->GetOccupyPtr();

	float div_per_volume = 0.0f;
	int count = 0;

	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			if(!occupy[j*width+i])
			{
				count += 1;
				div_per_volume += Uptr[j*(width+1)+i+1] - Uptr[j*(width+1)+i] + Vptr[(j+1)*width+i] - Vptr[j*width+i];
			}
		}
	}
	if(count == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}

	div_per_volume /= count;

	switch(cuda_openoctree_poisson_type)
	{
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE1:
		ZQ_CUDA_PoissonSolver2D::SolveClosedOctreePoissonRedBlack2D_MAC(Uptr,Vptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,div_per_volume,width,height,para.maxIter);
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE2:
		{
			int level0_num_red = level0_index_red.size()/2;
			int level1_num_red = level1_index_red.size()/2;
			int level2_num_red = level2_index_red.size()/2;
			int level3_num_red = level3_index_red.size()/2;
			int level0_num_black = level0_index_black.size()/2;
			int level1_num_black = level1_index_black.size()/2;
			int level2_num_black = level2_index_black.size()/2;
			int level3_num_black = level3_index_black.size()/2;
			int* level0_index_red_ptr = level0_num_red > 0 ? &level0_index_red[0] : 0;
			int* level1_index_red_ptr = level1_num_red > 0 ? &level1_index_red[0] : 0;
			int* level2_index_red_ptr = level2_num_red > 0 ? &level2_index_red[0] : 0;
			int* level3_index_red_ptr = level3_num_red > 0 ? &level3_index_red[0] : 0;
			int* level0_index_black_ptr = level0_num_black > 0 ? &level0_index_black[0] : 0;
			int* level1_index_black_ptr = level1_num_black > 0 ? &level1_index_black[0] : 0;
			int* level2_index_black_ptr = level2_num_black > 0 ? &level2_index_black[0] : 0;
			int* level3_index_black_ptr = level3_num_black > 0 ? &level3_index_black[0] : 0;
			ZQ_CUDA_PoissonSolver2D::SolveClosedOctreePoissonRedBlack2_2D_MAC(Uptr,Vptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,div_per_volume,width,height,para.maxIter,
				level0_num_red,level0_index_red_ptr,level0_num_black,level0_index_black_ptr,
				level1_num_red,level1_index_red_ptr,level1_num_black,level1_index_black_ptr,
				level2_num_red,level2_index_red_ptr,level2_num_black,level2_index_black_ptr,
				level3_num_red,level3_index_red_ptr,level3_num_black,level3_index_black_ptr);
		}
		break;
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		{
			const int index_channels = 4;
			const int neighborinfo_channels = 3;
			int level0_num_red = level0_index_red.size()/index_channels;
			int level1_num_red = level1_index_red.size()/index_channels;
			int level2_num_red = level2_index_red.size()/index_channels;
			int level3_num_red = level3_index_red.size()/index_channels;
			int level0_num_black = level0_index_black.size()/index_channels;
			int level1_num_black = level1_index_black.size()/index_channels;
			int level2_num_black = level2_index_black.size()/index_channels;
			int level3_num_black = level3_index_black.size()/index_channels;
			int level0_info_len_red = level0_neighborinfo_red.size()/neighborinfo_channels;
			int level1_info_len_red = level1_neighborinfo_red.size()/neighborinfo_channels;
			int level2_info_len_red = level2_neighborinfo_red.size()/neighborinfo_channels;
			int level3_info_len_red = level3_neighborinfo_red.size()/neighborinfo_channels;
			int level0_info_len_black = level0_neighborinfo_black.size()/neighborinfo_channels;
			int level1_info_len_black = level1_neighborinfo_black.size()/neighborinfo_channels;
			int level2_info_len_black = level2_neighborinfo_black.size()/neighborinfo_channels;
			int level3_info_len_black = level3_neighborinfo_black.size()/neighborinfo_channels;

			int* level0_index_red_ptr = level0_num_red > 0 ? &level0_index_red[0] : 0;
			int* level1_index_red_ptr = level1_num_red > 0 ? &level1_index_red[0] : 0;
			int* level2_index_red_ptr = level2_num_red > 0 ? &level2_index_red[0] : 0;
			int* level3_index_red_ptr = level3_num_red > 0 ? &level3_index_red[0] : 0;
			int* level0_index_black_ptr = level0_num_black > 0 ? &level0_index_black[0] : 0;
			int* level1_index_black_ptr = level1_num_black > 0 ? &level1_index_black[0] : 0;
			int* level2_index_black_ptr = level2_num_black > 0 ? &level2_index_black[0] : 0;
			int* level3_index_black_ptr = level3_num_black > 0 ? &level3_index_black[0] : 0;
			int* level0_info_red_ptr = level0_info_len_red > 0 ? &level0_neighborinfo_red[0] : 0;
			int* level1_info_red_ptr = level1_info_len_red > 0 ? &level1_neighborinfo_red[0] : 0;
			int* level2_info_red_ptr = level2_info_len_red > 0 ? &level2_neighborinfo_red[0] : 0;
			int* level3_info_red_ptr = level3_info_len_red > 0 ? &level3_neighborinfo_red[0] : 0;
			int* level0_info_black_ptr = level0_info_len_black > 0 ? &level0_neighborinfo_black[0] : 0;
			int* level1_info_black_ptr = level1_info_len_black > 0 ? &level1_neighborinfo_black[0] : 0;
			int* level2_info_black_ptr = level2_info_len_black > 0 ? &level2_neighborinfo_black[0] : 0;
			int* level3_info_black_ptr = level3_info_len_black > 0 ? &level3_neighborinfo_black[0] : 0;
			ZQ_CUDA_PoissonSolver2D::SolveClosedOctreePoissonRedBlack3_2D_MAC(Uptr,Vptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,div_per_volume,width,height,para.maxIter,
				level0_num_red,		level0_index_red_ptr,	level0_info_len_red,	level0_info_red_ptr,
				level0_num_black,	level0_index_black_ptr,	level0_info_len_black,	level0_info_black_ptr,
				level1_num_red,		level1_index_red_ptr,	level1_info_len_red,	level1_info_red_ptr,
				level1_num_black,	level1_index_black_ptr,	level1_info_len_black,	level1_info_black_ptr,
				level2_num_red,		level2_index_red_ptr,	level2_info_len_red,	level2_info_red_ptr,
				level2_num_black,	level2_index_black_ptr,	level2_info_len_black,	level2_info_black_ptr,
				level3_num_red,		level3_index_red_ptr,	level3_info_len_red,	level3_info_red_ptr,
				level3_num_black,	level3_index_black_ptr,	level3_info_len_black,	level3_info_black_ptr);
		}
		break;
	}

	return true;
}

bool ZQ_SmokeOctreeGrid2D::_solveGlobalOpenPoisson(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(srcGrid == 0 || A == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	
	int dim = A->n;
	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[A->m];
	memset(b,0,sizeof(double)*A->m);
	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	int i_shift,j_shift;
	int Length[4] = {1,2,4,8};
	int size;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lambdaCount = lamdainfos.size();

	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			break;
		}
		
		i_shift = lamdainfos[i].i1*size;
		j_shift = lamdainfos[i].j1*size;
		double val = 0;
		int size = Length[lvl];
		for(int sj = 0;sj < size;sj++)
		{
			val -= Uptr[(j_shift+sj)*(width+1)+i_shift];
			val += Uptr[(j_shift+sj)*(width+1)+i_shift+size];
		}
		for(int si = 0;si < size;si++)
		{
			val -= Vptr[j_shift*width+i_shift+si];
			val += Vptr[(j_shift+size)*width+i_shift+si];
		}
		b[lamdainfos[i].index] = val;
	}

	double tol = 1e-9;
	int it = 0;
	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(A,b,x0,para.maxIter,tol,x,it,false);

	double* pressure = new double[dim];
	memcpy(pressure,x,sizeof(double)*dim);

	delete []x;
	delete []x0;
	delete []b;

	if(deltax)
		delete []deltax;
	deltax = new double[uCount+vCount];
	memset(deltax,0,sizeof(double)*(uCount+vCount));
		
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight3 = height/8;
	
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
				len1 = 4;
			}
			else
			{
				col1 = -1;
				len1 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		
		switch(Uinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
				len2 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		
		size = Length[lvl];
	
		if(lvl == Uinfos[i].level1)
		{
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		if( i_shift != 0 && i_shift != width)
		{
			double val = 1.0/(len1+len2);
			deltax[i] = -val*(pressure[col2]-pressure[col1]);
		}
		else if(i_shift == 0)
		{
			double val = 1.0/(2*len2);
			deltax[i] = -val*(pressure[col2]-0);
		}
		else if(i_shift == width)
		{
			double val = 1.0/(2*len1);
			deltax[i] = -val*(0-pressure[col1]);
		}
		
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1,Vinfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Vinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
				len1 = 4;
			}
			else
			{
				col1 = -1;
				len1 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Vinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
				len2 = 4;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		
		size = Length[lvl];
		
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		if(j_shift != 0 && j_shift != height)
		{
			double val = 1.0/(len1+len2);
			deltax[i+uCount] = -val*(pressure[col2]-pressure[col1]);
		}
		else if(j_shift == 0)
		{
			double val = 1.0/(2*len2);
			deltax[i+uCount] = -val*(pressure[col2]-0);
		}
		else if(j_shift == height)
		{
			double val = 1.0/(2*len1);
			deltax[i+uCount] = -val*(0-pressure[col1]);
		}
	}

	delete []pressure;

	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;
	
	deltaU = new float[(width+1)*height];
	memset(deltaU,0,sizeof(float)*(width+1)*height);
	deltaV = new float[(height+1)*width];
	memset(deltaV,0,sizeof(float)*(height+1)*width);
	
	
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1, Uinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Uinfos[i].level1)
		{

			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}

		for(int sj = 0;sj < size;sj++)
		{
			deltaU[(j_shift+sj)*(width+1)+i_shift] += deltax[i];
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1, Vinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		for(int si = 0;si < size;si++)
		{
			deltaV[j_shift*width+i_shift+si] += deltax[i+uCount];
		}
	}

	delete []deltax;
	deltax = 0;

	for(int i = 0;i < height*(width+1);i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < width*(height+1);i++)
		Vptr[i] += deltaV[i];

	delete []deltaU;
	delete []deltaV;
	
	deltaU = 0;
	deltaV = 0;
	
	return true;
}


bool ZQ_SmokeOctreeGrid2D::_solveGlobalClosedPoisson(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(srcGrid == 0 || A == 0)
		return false;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lambdaCount = lamdainfos.size();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();

	int dim = A->n;
	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[A->m];
	memset(b,0,sizeof(double)*A->m);
	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	int i_shift,j_shift;
	int Length[4] = {1,2,4,8};
	//double volume[4] = {1,4,16,64};
	int size;
	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			break;
		}

		i_shift = lamdainfos[i].i1*size;
		j_shift = lamdainfos[i].j1*size;
		double val = 0;
		int size = Length[lvl];
		for(int sj = 0;sj < size;sj++)
		{
			val -= Uptr[(j_shift+sj)*(width+1)+i_shift];
			val += Uptr[(j_shift+sj)*(width+1)+(i_shift+size)];
		}
		for(int si = 0;si < size;si++)
		{
			val -= Vptr[j_shift*width+i_shift+si];
			val += Vptr[(j_shift+size)*width+i_shift+si];
		}
		b[lamdainfos[i].index] = val;///volume[lamdainfos[i].level1];
	}

	double tol = 1e-9;
	int it = 0;
	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(A, b, x0, para.maxIter, tol, x, it, false);




	double* pressure = new double[dim+1];
	memcpy(pressure,x,sizeof(double)*dim);
	pressure[dim] = 0;

	delete []x;
	delete []x0;
	delete []b;

	if(deltax)
		delete []deltax;
	deltax = new double[uCount+vCount];
	memset(deltax,0,sizeof(double)*(uCount+vCount));

	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight3 = height/8;

	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
				len1 = 4;
			}
			else
			{
				col1 = -1;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Uinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		double val = 1.0/(len1+len2);
		size = Length[lvl];

		if(lvl == Uinfos[i].level1)
		{
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		if( i_shift != 0 && i_shift != width)
		{
			deltax[i] = -val*(pressure[col2]-pressure[col1]);
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1,Vinfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Vinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
				len1 = 4;
			}
			else
				col1 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Vinfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
				len2 = 4;
			}
			else
			{
				col2 = -1;
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		double val = 1.0/(len1+len2);
		size = Length[lvl];

		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		if(j_shift != 0 && j_shift != height)
		{
			deltax[i+uCount] = -val*(pressure[col2]-pressure[col1]);
		}
	}

	delete []pressure;

	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;

	deltaU = new float[(width+1)*height];
	memset(deltaU,0,sizeof(float)*(width+1)*height);
	deltaV = new float[(height+1)*width];
	memset(deltaV,0,sizeof(float)*(height+1)*width);


	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1, Uinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Uinfos[i].level1)
		{

			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}


		for(int sj = 0;sj < size;sj++)
		{
			deltaU[(j_shift+sj)*(width+1)+i_shift] += deltax[i];
		}

	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1, Vinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		for(int si = 0;si < size;si++)
		{
			deltaV[j_shift*width+i_shift+si] += deltax[i+uCount];
		}	
	}

	delete []deltax;
	deltax = 0;

	for(int i = 0;i < height*(width+1);i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < width*(height+1);i++)
		Vptr[i] += deltaV[i];

	delete []deltaU;
	delete []deltaV;

	deltaU = 0;
	deltaV = 0;

	return true;
}

bool ZQ_SmokeOctreeGrid2D::_solveGlobalOpenFlux(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(srcGrid == 0)
		return false;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lambdaCount = lamdainfos.size();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight3 = height/8;

	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[dim];
	double* oldVelocity = new double[uCount+vCount];

	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	memset(b,0,sizeof(double)*dim);

	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();

	int size;
	int i_shift,j_shift;
	double area[4] = {1,2,4,8};
	double weight[4] = {1,2,4,8};
	//double volume[4] = {1,4,16,64};
	//double weight[4] = {1,1,1,1};
	double volume[4] = {1,1,1,1};
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1, Uinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Uinfos[i].level1)
		{
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		double val = 0;
		for(int sj = 0;sj < size;sj++)
		{
			val += Uptr[(j_shift+sj)*(width+1)+i_shift];
		}
		
		b[Uinfos[i].index] = 2*(val/size)*weight[lvl];
		oldVelocity[Uinfos[i].index] = val / (size);
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1, Vinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		double val = 0;
		for(int si = 0;si < size;si++)
		{
			val += Vptr[j_shift*width+i_shift+si];
		}
	
		b[Vinfos[i].index] = 2*(val/size)*weight[lvl];
		oldVelocity[Vinfos[i].index] = val / (size);
	}

	bool*& occupy = srcGrid->GetOccupyPtr();
	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			break;
		}
		
		i_shift = lamdainfos[i].i1*size;
		j_shift = lamdainfos[i].j1*size;
		double val = 0;
		for(int sj = 0;sj < size;sj++)
		{
			if(occupy[(j_shift+sj)*width+i_shift] || (i_shift > 0 && occupy[(j_shift+sj)*width+i_shift-1]))
				val -= Uptr[(j_shift+sj)*(width+1)+i_shift];
			if(occupy[(j_shift+sj)*width+i_shift+size-1] || (i_shift+size < width && occupy[(j_shift+sj)*width+i_shift+size]))
				val += Uptr[(j_shift+sj)*(width+1)+i_shift+size];
		}
		for(int si = 0;si < size;si++)
		{
			if(occupy[j_shift*width+i_shift+si] || (j_shift > 0 && occupy[(j_shift-1)*width+i_shift+si]))
				val -= Vptr[j_shift*width+i_shift+si];
			if(occupy[(j_shift+size-1)*width+i_shift+si] || (j_shift+size < height && occupy[(j_shift+size)*width+i_shift+si]))
				val += Vptr[(j_shift+size)*width+i_shift+si];
		}
		b[lamdainfos[i].index] = val / volume[lvl];
	}

	double tol = 1e-9;
	int it = 0;

	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(A, b, x0, para.maxIter, tol, x, it, false);

	if(deltax)
	{
		delete []deltax;
	}
	deltax = new double[uCount+vCount];
	for(int i = 0;i < uCount+vCount;i++)
		deltax[i] = x[i] - oldVelocity[i];

	delete []x0;
	delete []x;
	delete []b;
	delete []oldVelocity;

	
	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;
	
	deltaU = new float[(width+1)*height];
	memset(deltaU,0,sizeof(float)*(width+1)*height);
	deltaV = new float[(height+1)*width];
	memset(deltaV,0,sizeof(float)*(height+1)*width);
	
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1, Uinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Uinfos[i].level1)
		{
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}

	
		for(int sj = 0;sj < size;sj++)
		{
			deltaU[(j_shift+sj)*(width+1)+i_shift] += deltax[Uinfos[i].index];
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1, Vinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}
	
		for(int si = 0;si < size;si++)
		{
			deltaV[j_shift*width+i_shift+si] += deltax[Vinfos[i].index];
		}
	}

	delete []deltax;
	deltax = 0;

	

	for(int i = 0;i < (width+1)*height;i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < (height+1)*width;i++)
		Vptr[i] += deltaV[i];
	
	delete []deltaU;
	delete []deltaV;
	deltaU = 0;
	deltaV = 0;

	return true;
}



bool ZQ_SmokeOctreeGrid2D::_solveGlobalClosedFlux(ZQ_SmokeGrid2D* srcGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(srcGrid == 0)
		return false;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int lambdaCount = lamdainfos.size();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight3 = height/8;

	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[dim];
	double* oldVelocity = new double[uCount+vCount];

	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	memset(b,0,sizeof(double)*dim);

	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();

	int size;
	int i_shift,j_shift;
	double area[4] = {1,2,4,8};
	double weight[4] = {1,2,4,8};
	double volume[4] = {1,4,16,64};
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1, Uinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Uinfos[i].level1)
		{
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		double val = 0;
		for(int sj = 0;sj < size;sj++)
		{
			val += Uptr[(j_shift+sj)*(width+1)+i_shift];
		}
		if(Uinfos[i].index >= 0)
		{
			b[Uinfos[i].index] = 2*(val/size)*area[lvl];
			oldVelocity[Uinfos[i].index] = val / (size);
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1, Vinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		double val = 0;
		for(int si = 0;si < size;si++)
		{
			val += Vptr[j_shift*width+i_shift+si];
		}
		if(Vinfos[i].index >= 0)
		{
			b[Vinfos[i].index] = 2*(val/size)*area[lvl];
			oldVelocity[Vinfos[i].index] = val / (size);
		}
	}

	bool*& occupy = srcGrid->GetOccupyPtr();
	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			break;
		}
		
		i_shift = lamdainfos[i].i1*size;
		j_shift = lamdainfos[i].j1*size;
		double val = 0;
		for(int sj = 0;sj < size;sj++)
		{
			if((i_shift > 0 && occupy[(j_shift+sj)*width+i_shift-1]) || i_shift == 0)
				val -= Uptr[(j_shift+sj)*(width+1)+i_shift];
			if((i_shift+size < width && occupy[(j_shift+sj)*width+i_shift+size]) || i_shift == width)
				val += Uptr[(j_shift+sj)*(width+1)+i_shift+size];
		}
		for(int si = 0;si < size;si++)
		{
			if((j_shift > 0 && occupy[(j_shift-1)*width+i_shift+si]) || j_shift == 0)
				val -= Vptr[j_shift*width+i_shift+si];
			if((j_shift+size < height && occupy[(j_shift+size)*width+i_shift+si]) || j_shift == height)
				val += Vptr[(j_shift+size)*width+i_shift+si];
		}
		b[lamdainfos[i].index] = val;
	}

	double tol = 1e-9;
	int it = 0;

	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(A,b,x0,para.maxIter,tol,x,it,false);

	if(deltax)
	{
		delete []deltax;
	}
	deltax = new double[uCount+vCount];
	for(int i = 0;i < uCount+vCount;i++)
		deltax[i] = x[i] - oldVelocity[i];

	delete []x0;
	delete []x;
	delete []b;
	delete []oldVelocity;

	
	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;
	
	deltaU = new float[(width+1)*height];
	memset(deltaU,0,sizeof(float)*(width+1)*height);
	deltaV = new float[(height+1)*width];
	memset(deltaV,0,sizeof(float)*(height+1)*width);
	
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1, Uinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Uinfos[i].level1)
		{
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		if(Uinfos[i].index >= 0)
		{
			for(int sj = 0;sj < size;sj++)
			{
				deltaU[(j_shift+sj)*(width+1)+i_shift] += deltax[Uinfos[i].index];
			}
		}
	}

	for(int i = 0;i < vCount;i++)
	{
		int lvl = __min(Vinfos[i].level1, Vinfos[i].level2);
		switch(lvl)
		{
		case 0:
			size = 1;
			break;
		case 1:
			size = 2;
			break;
		case 2:
			size = 4;
			break;
		case 3:
			size = 8;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
		if(lvl == Vinfos[i].level1)
		{
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}
	
		if(Vinfos[i].index >= 0)
		{
			for(int si = 0;si < size;si++)
			{
				deltaV[j_shift*width+i_shift+si] += deltax[Vinfos[i].index];
			}
		}
	}

	delete []deltax;
	deltax = 0;

	for(int i = 0;i < (width+1)*height;i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < (height+1)*width;i++)
		Vptr[i] += deltaV[i];


	delete []deltaU;
	delete []deltaV;
	deltaU = 0;
	deltaV = 0;

	return true;
}

bool ZQ_SmokeOctreeGrid2D::_buildSubMatrix()
{
	int row2,col2;
	if(subMatrix2)
	{
		delete []subMatrix2;
		subMatrix2 = 0;
	}
	if(subMatrix4)
	{
		delete []subMatrix4;
		subMatrix4 = 0;
	}
	if(subMatrix8)
	{
		delete []subMatrix8;
		subMatrix8 = 0;
	}
	if(!ZQ_InitSubMatPoisson2D(2,&subMatrix2,&row2,&col2))
	{
		printf("init subMatrix2 fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}

	if(!ZQ_InitSubMatPoisson2D(4,&subMatrix4,&row2,&col2))
	{
		printf("init subMatrix4 fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}
	if(!ZQ_InitSubMatPoisson2D(8,&subMatrix8,&row2,&col2))
	{
		printf("init subMatrix2 fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}
	
	return true;
}


void ZQ_SmokeOctreeGrid2D::_seedFilling(const int width, const int height, int *distanceField)
{
	int head = -1;
	int tail = -1;
	int* queue_i = new int[width*height];
	int* queue_j = new int[width*height];
	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			if(distanceField[j*width+i] == 0)
			{
				tail++;
				queue_i[tail] = i;
				queue_j[tail] = j;
			}
		}
	}
	while(head < tail)
	{
		head ++;
		int cur_i = queue_i[head];
		int cur_j = queue_j[head];
		int cur_dis = distanceField[cur_j*width+cur_i];
		for(int sj = __max(0,cur_j-1);sj <= __min(height-1,cur_j+1);sj++)
		{
			for(int si = __max(0,cur_i-1);si <= __min(width-1,cur_i+1);si++)
			{
				if(si == cur_i && sj == cur_j)
					continue;
				if(distanceField[sj*width+si] < 0 || distanceField[sj*width+si] > cur_dis+1)
				{
					distanceField[sj*width+si] = cur_dis+1;
					tail ++;
					queue_i[tail] = si;
					queue_j[tail] = sj;
				}
			}
		}
	}
	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			if(distanceField[j*width+i] < 0)
			{
				printf("distance[%d,%d] = %d\n",i,j,distanceField[j*width+i]);
			}
		}
	}
	delete []queue_i;
	delete []queue_j;
}




void ZQ_SmokeOctreeGrid2D::_find_special_label(const int width, const int height, const bool* special_mask, int *special_label, int &region_num)
{

	for(int i = 0;i < width*height;i++)
	{
		if(special_mask[i])
			special_label[i] = -1;
		else
			special_label[i] = 0;
	}


	/*IplImage* img = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,3);
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(special_mask[i*width+j])
			{
				cvSet2D(img,i,j,cvScalar(0,255,0));
			}
			else
			{
				cvSet2D(img,i,j,cvScalar(0,0,0));
			}
		}
	}

	cvNamedWindow("mask");
	cvShowImage("mask",img);
	cvWaitKey(0);
	cvReleaseImage(&img);
	cvDestroyWindow("mask");*/

	bool no_new_region = false;
	int num_region = 0;

	bool no_new_cell_for_region = false;

	while(!no_new_region)
	{
		no_new_region = true;

		int cur_i,cur_j;
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(special_label[j*width+i] == -1)
				{
					special_label[j*width+i] = num_region+1;
					no_new_region = false;
					cur_i = i;
					cur_j = j;
					break;
				}
			}
			if(!no_new_region)
				break;
		}

		if(!no_new_region)
		{
			int* queue_i = new int[width*height];
			int* queue_j = new int[width*height];
			int tail = -1;
			int head = 0;
			queue_i[head] = cur_i;
			queue_j[head] = cur_j;

			while(tail < head)
			{
				tail++;
				int iii = queue_i[tail];
				int jjj = queue_j[tail];
				if(jjj+1 < height)
				{
					if(special_label[(jjj+1)*width+iii] == -1)
					{
						head ++;
						queue_i[head] = iii;
						queue_j[head] = jjj+1;
						special_label[(jjj+1)*width+iii] = num_region+1;
					}
				}
				if(jjj-1 >= 0)
				{
					if(special_label[(jjj-1)*width+iii] == -1)
					{
						head ++;
						queue_i[head] = iii;
						queue_j[head] = jjj-1;
						special_label[(jjj-1)*width+iii] = num_region+1;
					}
				}
				if(iii+1 < width)
				{
					if(special_label[jjj*width+(iii+1)] == -1)
					{
						head ++;
						queue_i[head] = iii+1;
						queue_j[head] = jjj;
						special_label[jjj*width+(iii+1)] = num_region+1;
					}
				}
				if(iii-1 >= 0)
				{
					if(special_label[jjj*width+(iii-1)] == -1)
					{
						head ++;
						queue_i[head] = iii-1;
						queue_j[head] = jjj;
						special_label[jjj*width+(iii-1)] = num_region+1;
					}
				}
			}
			delete []queue_i;
			delete []queue_j;

			num_region ++;

		}
	}
	region_num = num_region;
}

bool ZQ_SmokeOctreeGrid2D::_buildOctreeInfo(const ZQ_SmokeGrid2D *srcGrid, const ZQ_SmokeSimulationPara2D& para, const int *interestingRegionDist)
{
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();

	const bool*& occupyPtr = srcGrid->GetOccupyPtr();

	int* distanceField = new int[width*height];
	for(int i = 0;i < width*height;i++)
		distanceField[i] = -1;

	bool hasSeed = false;

	for(int i = 0;i < height*width;i++)
	{
		if(occupyPtr[i])
		{
			hasSeed = true;
			distanceField[i] = 0;
		}
	}
	if(interestingRegionDist)
	{
		for(int i = 0;i < height*width;i++)
		{
			if(interestingRegionDist[i] > 0)
			{
				hasSeed = true;
				distanceField[i] = interestingRegionDist[i];
			}
		}
	}

	if(hasSeed)
	{
		_seedFilling(width,height,distanceField);
	}
	else
	{
		for(int i = 0;i < height*width;i++)
		{
			distanceField[i] = width+height;
		}
	}

	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;

	if(leafLevel0)
		delete []leafLevel0;
	if(leafLevel1)
		delete []leafLevel1;
	if(leafLevel2)
		delete []leafLevel2;
	if(leafLevel3)
		delete []leafLevel3;
	leafLevel0 = new bool[width*height];
	leafLevel1 = new bool[levelWidth1*levelHeight1];
	leafLevel2 = new bool[levelWidth2*levelHeight2];
	leafLevel3 = new bool[levelWidth3*levelHeight3];

	memset(leafLevel0,false,sizeof(bool)*width*height);
	memset(leafLevel1,false,sizeof(bool)*levelWidth1*levelHeight1);
	memset(leafLevel2,false,sizeof(bool)*levelWidth2*levelHeight2);
	memset(leafLevel3,false,sizeof(bool)*levelWidth3*levelHeight3);

	for(int j = 0;j < levelHeight3;j++)
	{
		for(int i = 0;i < levelWidth3;i++)
		{
			bool flag = true;
			for(int sj = 0;sj < 8;sj++)
			{
				if(!flag)
					break;
				for(int si = 0;si < 8;si ++)
				{
					if(!flag)
						break;
					if(distanceField[(j*8+sj)*width+(i*8+si)] <= level_dis_thresh3)
					{
						flag = false;
						break;
					}
				}
			}
			leafLevel3[j*levelWidth3+i] = flag;
		}
	}
	for(int j = 0;j < levelHeight2;j++)
	{
		for(int i = 0;i < levelWidth2;i++)
		{
			if(leafLevel3[(j/2)*levelWidth3+(i/2)])
				leafLevel2[j*levelWidth2+i] = false;
			else
			{
				bool flag = true;
				for(int sj = 0;sj < 4;sj++)
				{
					if(!flag)
						break;
					for(int si = 0;si < 4;si++)
					{
						if(!flag)
							break;

						if(distanceField[(j*4+sj)*width+(i*4+si)] <= level_dis_thresh2)
						{
							flag = false;
							break;
						}

					}
				}
				leafLevel2[j*levelWidth2+i] = flag;
			}
		}
	}

	for(int j = 0;j < levelHeight1;j++)
	{
		for(int i = 0;i < levelWidth1;i++)
		{
			if(leafLevel3[(j/4)*levelWidth3+(i/4)] || leafLevel2[(j/2)*levelWidth2+(i/2)])
				leafLevel1[j*levelWidth1+i] = false;
			else
			{
				bool flag = true;
				for(int sj = 0;sj < 2;sj++)
				{
					if(!flag)
						break;
					for(int si = 0;si < 2;si++)
					{
						if(!flag)
							break;
						if(distanceField[(j*2+sj)*width+(i*2+si)] <= level_dis_thresh1)
						{
							flag = false;
							break;
						}
					}
				}
				leafLevel1[j*levelWidth1+i] = flag;
			}
		}
	}
	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			if(!leafLevel1[(j/2)*levelWidth1+(i/2)] && !leafLevel2[(j/4)*levelWidth2+(i/4)]
			&& !leafLevel3[(j/8)*levelWidth3+(i/8)] && !occupyPtr[j*width+i])
			{
				leafLevel0[j*width+i] = true;
			}
		}
	}
	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			bool has_use = true;
			has_use = occupyPtr[j*width+i] || leafLevel0[j*width+i] || leafLevel1[(j/2)*levelWidth1+(i/2)]
			|| leafLevel2[(j/4)*levelWidth2+(i/4)] || leafLevel3[(j/8)*levelWidth3+(i/8)];
			if(!has_use)
			{
				printf("%3d,%3d is not use\n",i,j);
			}
		}
	}
	delete []distanceField;
	distanceField = 0;



	if(para.useGuidingControl)
	{
		for(int i = 0;i < specialRegion_A.size();i++)
		{
			if(specialRegion_A[i])
				ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(specialRegion_A[i]);
		}
		specialRegion_A.clear();
		if(specialRegion_mask)
		{
			delete []specialRegion_mask;
		}
		specialRegion_mask = new bool[width*height];
		memset(specialRegion_mask,0,sizeof(bool)*width*height);
		if(specialRegion_label)
		{
			delete []specialRegion_label;
		}
		specialRegion_label = new int[width*height];
		memset(specialRegion_label,0,sizeof(int)*width*height);

		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(leafLevel0[j*width+i] || leafLevel1[(j/2)*(width/2)+(i/2)])
					specialRegion_mask[j*width+i] = true;
				else
					specialRegion_mask[j*width+i] = false;
			}
		}
		
		_find_special_label(width,height,specialRegion_mask,specialRegion_label,specialRegion_num);

		bool* new_Occupy = new bool[width*height];

		if(!guiding_specialRegion_useCUDA)
		{
			for(int iregion = 0;iregion < specialRegion_num;iregion++)
			{
				for(int i = 0;i < width*height;i++)
					new_Occupy[i] = occupyPtr[i] || (specialRegion_label[i] != iregion+1);
				ZQ::ZQ_PoissonSolver::BuildClosedPoisson<double>(width,height,new_Occupy,&A,false);
				specialRegion_A.push_back(A);
			}
			delete []new_Occupy;
		}	
	}

	return true;
}

void ZQ_SmokeOctreeGrid2D::_selectInfos_lambda(const int width, const int height, const bool* leafLevel0, const bool* leafLevel1, const bool* leafLevel2, const bool* leafLevel3, 
											   std::vector<Info>& lamdainfos)
{
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	lamdainfos.clear();

	/*level3 lambda*/
	for(int j = 0;j < levelHeight3;j++)
	{
		for(int i = 0;i < levelWidth3;i++)
		{
			if(leafLevel3[j*levelWidth3+i])
			{
				Info tmpInfo(3,i,j,3,i,j);
				lamdainfos.push_back(tmpInfo);
			}
		}
	}

	/*level2 lambda*/
	for(int j = 0;j < levelHeight2;j++)
	{
		for(int i = 0;i < levelWidth2;i++)
		{
			if(leafLevel2[j*levelWidth2+i])
			{
				Info tmpInfo(2,i,j,2,i,j);
				lamdainfos.push_back(tmpInfo);
			}
		}
	}

	/*level1 lambda*/
	for(int j = 0;j < levelHeight1;j++)
	{
		for(int i = 0;i < levelWidth1;i++)
		{
			if(leafLevel1[j*levelWidth1+i])
			{
				Info tmpInfo(1,i,j,1,i,j);
				lamdainfos.push_back(tmpInfo);
			}
		}
	}

	/*level0 lambda*/
	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			if(leafLevel0[j*width+i])
			{
				Info tmpInfo(0,i,j,0,i,j);
				lamdainfos.push_back(tmpInfo);
			}
		}
	}
}

void ZQ_SmokeOctreeGrid2D::_selectInfosForCPUSolvers(const int width, const int height, const bool *leafLevel0, const bool *leafLevel1, const bool *leafLevel2, const bool *leafLevel3, 
													 std::vector<Info> &Uinfos, std::vector<Info> &Vinfos, std::vector<Info> &lamdainfos)
{
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;

	Uinfos.clear();
	Vinfos.clear();
	lamdainfos.clear();

	/* level3 U */
	for(int j = 0;j < levelHeight3;j ++)
	{
		for(int i = 0;i <= levelWidth3;i++)
		{
			if(i == 0)
			{
				if(leafLevel3[j*levelWidth3+i])
				{
					Info tmpInfo(3,i-1,j,3,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
			else if(i == levelWidth3)
			{
				if(leafLevel3[j*levelWidth3+i-1])
				{
					Info tmpInfo(3,i-1,j,3,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
			else
			{
				if(leafLevel3[j*levelWidth3+i] && leafLevel3[j*levelWidth3+i-1])
				{
					Info tmpInfo(3,i-1,j,3,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
		}
	}

	/*level2 U*/
	for(int j = 0;j < levelHeight2;j++)
	{
		for(int i = 0;i <= levelWidth2;i++)
		{
			if(i == 0)
			{
				if(leafLevel2[j*levelWidth2+i])
				{
					Info tmpInfo(3,i/2-1,j/2,2,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
			else if(i == levelWidth2)
			{
				if(leafLevel2[j*levelWidth2+i-1])
				{
					Info tmpInfo(2,i-1,j,3,i/2,j/2);
					Uinfos.push_back(tmpInfo);
				}
			}
			else
			{
				if(leafLevel2[j*levelWidth2+i-1] && leafLevel2[j*levelWidth2+i])
				{
					Info tmpInfo(2,i-1,j,2,i,j);
					Uinfos.push_back(tmpInfo);
				}
				else if( i%2 == 0 && leafLevel2[j*levelWidth2+i-1] && leafLevel3[(j/2)*levelWidth3+(i/2)])
				{
					Info tmpInfo(2,i-1,j,3,i/2,j/2);
					Uinfos.push_back(tmpInfo);
				}
				else if( i%2 == 0 && leafLevel3[j/2*levelWidth3+(i/2-1)] && leafLevel2[j*levelWidth2+i])
				{
					Info tmpInfo(3,i/2-1,j/2,2,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}			
		}
	}

	/*level1 U*/
	for(int j = 0; j < levelHeight1;j++)
	{
		for(int i = 0;i <= levelWidth1;i++)
		{
			if(i == 0)
			{
				if(leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(3,i/4-1,j/4,1,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
			else if(i == levelWidth1)
			{
				if(leafLevel1[j*levelWidth1+i-1])
				{
					Info tmpInfo(1,i-1,j,3,i/4,j/4);
					Uinfos.push_back(tmpInfo);
				}
			}
			else
			{
				if(leafLevel1[j*levelWidth1+i-1] && leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(1,i-1,j,1,i,j);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%2 == 0 && leafLevel1[j*levelWidth1+i-1] && leafLevel2[(j/2)*levelWidth2+(i/2)])
				{
					Info tmpInfo(1,i-1,j,2,i/2,j/2);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%4 == 0 && leafLevel1[j*levelWidth1+i-1] && leafLevel3[(j/4)*levelWidth3+(i/4)])
				{
					Info tmpInfo(1,i-1,j,3,i/4,j/4);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%2 == 0 && leafLevel2[j/2*levelWidth2+(i/2-1)] && leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(2,i/2-1,j/2,1,i,j);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%4 == 0 && leafLevel3[j/4*levelWidth3+(i/4-1)] && leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(3,i/4-1,j/4,1,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
		}
	}

	/*level0 U*/
	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i <= width;i++)
		{
			if(i == 0)
			{
				if(leafLevel0[j*width+i])
				{
					Info tmpInfo(3,i/8-1,j/8,0,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
			else if(i == width)
			{
				if(leafLevel0[j*width+i-1])
				{
					Info tmpInfo(0,i-1,j,3,i/8,j/8);
					Uinfos.push_back(tmpInfo);
				}
			}
			else 
			{
				if(leafLevel0[j*width+i-1] && leafLevel0[j*width+i])
				{
					Info tmpInfo(0,i-1,j,0,i,j);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%2 == 0 && leafLevel0[j*width+i-1] && leafLevel1[(j/2)*levelWidth1+(i/2)])
				{
					Info tmpInfo(0,i-1,j,1,i/2,j/2);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%4 == 0 && leafLevel0[j*width+i-1] && leafLevel2[(j/4)*levelWidth2+(i/4)])
				{
					Info tmpInfo(0,i-1,j,2,i/4,j/4);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%8 == 0 && leafLevel0[j*width+i-1] && leafLevel3[(j/8)*levelWidth3+(i/8)])
				{
					Info tmpInfo(0,i-1,j,3,i/8,j/8);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%2 == 0 && leafLevel1[j/2*levelWidth1+(i/2-1)] && leafLevel0[j*width+i])
				{
					Info tmpInfo(1,i/2-1,j/2,0,i,j);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%4 == 0 && leafLevel2[j/4*levelWidth2+(i/4-1)] && leafLevel0[j*width+i])
				{
					Info tmpInfo(2,i/4-1,j/4,0,i,j);
					Uinfos.push_back(tmpInfo);
				}
				else if(i%8 == 0 && leafLevel3[j/8*levelWidth3+(i/8-1)] && leafLevel0[j*width+i])
				{
					Info tmpInfo(3,i/8-1,j/8,0,i,j);
					Uinfos.push_back(tmpInfo);
				}
			}
		}
	}

	/*level3 V*/
	for(int j = 0;j <= levelHeight3;j++)
	{
		for(int i = 0;i < levelWidth3;i++)
		{
			if(j == 0)
			{
				if(leafLevel3[j*levelWidth3+i])
				{
					Info tmpInfo(3,i,j-1,3,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
			else if(j == levelHeight3)
			{
				if(leafLevel3[(j-1)*levelWidth3+i])
				{
					Info tmpInfo(3,i,j-1,3,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
			else
			{
				if(leafLevel3[(j-1)*levelWidth3+i] && leafLevel3[j*levelWidth3+i])
				{
					Info tmpInfo(3,i,j-1,3,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
		}
	}

	/*level2 V*/
	for(int j = 0;j <= levelHeight2;j++)
	{
		for(int i = 0;i < levelWidth2;i++)
		{
			if(j == 0)
			{
				if(leafLevel2[j*levelWidth2+i])
				{
					Info tmpInfo(3,i/2,j/2-1,2,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
			else if(j == levelHeight2)
			{
				if(leafLevel2[(j-1)*levelWidth2+i])
				{
					Info tmpInfo(2,i,j-1,3,i/2,j/2);
					Vinfos.push_back(tmpInfo);
				}
			}
			else
			{
				if(leafLevel2[(j-1)*levelWidth2+i] && leafLevel2[j*levelWidth2+i])
				{
					Info tmpInfo(2,i,j-1,2,i,j);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%2 == 0 && leafLevel2[(j-1)*levelWidth2+i] && leafLevel3[(j/2)*levelWidth3+(i/2)])
				{
					Info tmpInfo(2,i,j-1,3,i/2,j/2);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%2 == 0 && leafLevel3[(j/2-1)*levelWidth3+i/2] && leafLevel2[j*levelWidth2+i])
				{
					Info tmpInfo(3,i/2,j/2-1,2,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
		}
	}

	/*level1 V*/
	for(int j = 0;j <= levelHeight1;j++)
	{
		for(int i = 0;i < levelWidth1;i++)
		{
			if(j == 0)
			{
				if(leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(3,i/4,j/4-1,1,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
			else if(j == levelHeight1)
			{
				if(leafLevel1[(j-1)*levelWidth1+i])
				{
					Info tmpInfo(1,i,j-1,3,i/4,j/4);
					Vinfos.push_back(tmpInfo);
				}
			}
			else
			{
				if(leafLevel1[(j-1)*levelWidth1+i] && leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(1,i,j-1,1,i,j);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%2 == 0 && leafLevel1[(j-1)*levelWidth1+i] && leafLevel2[(j/2)*levelWidth2+(i/2)])
				{
					Info tmpInfo(1,i,j-1,2,i/2,j/2);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%4 == 0 && leafLevel1[(j-1)*levelWidth1+i] && leafLevel3[(j/4)*levelWidth3+(i/4)])
				{
					Info tmpInfo(1,i,j-1,3,i/4,j/4);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%2 == 0 && leafLevel2[(j/2-1)*levelWidth2+i/2] && leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(2,i/2,j/2-1,1,i,j);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%4 == 0 && leafLevel3[(j/4-1)*levelWidth3+i/4] && leafLevel1[j*levelWidth1+i])
				{
					Info tmpInfo(3,i/4,j/4-1,1,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
		}
	}

	/*level0 V*/
	for(int j = 0;j <= height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			if(j == 0)
			{
				if(leafLevel0[j*width+i])
				{
					Info tmpInfo(3,i/8,j/8-1,0,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
			else if(j == height)
			{
				if(leafLevel0[(j-1)*width+i])
				{
					Info tmpInfo(0,i,j-1,3,i/8,j/8);
					Vinfos.push_back(tmpInfo);
				}
			}
			else
			{
				if(leafLevel0[(j-1)*width+i] && leafLevel0[j*width+i])
				{
					Info tmpInfo(0,i,j-1,0,i,j);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%2 == 0 && leafLevel0[(j-1)*width+i] && leafLevel1[(j/2)*levelWidth1+(i/2)])
				{
					Info tmpInfo(0,i,j-1,1,i/2,j/2);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%4 == 0 && leafLevel0[(j-1)*width+i] && leafLevel2[(j/4)*levelWidth2+(i/4)])
				{
					Info tmpInfo(0,i,j-1,2,i/4,j/4);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%8 == 0 && leafLevel0[(j-1)*width+i] && leafLevel3[(j/8)*levelWidth3+(i/8)])
				{
					Info tmpInfo(0,i,j-1,3,i/8,j/8);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%2 == 0 && leafLevel1[(j/2-1)*levelWidth1+i/2] && leafLevel0[j*width+i])
				{
					Info tmpInfo(1,i/2,j/2-1,0,i,j);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%4 == 0 && leafLevel2[(j/4-1)*levelWidth2+i/4] && leafLevel0[j*width+i])
				{
					Info tmpInfo(2,i/4,j/4-1,0,i,j);
					Vinfos.push_back(tmpInfo);
				}
				else if(j%8 == 0 && leafLevel3[(j/8-1)*levelWidth3+i/8] && leafLevel0[j*width+i])
				{
					Info tmpInfo(3,i/8,j/8-1,0,i,j);
					Vinfos.push_back(tmpInfo);
				}
			}
		}
	}


	_selectInfos_lambda(width,height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos);
	
}

void ZQ_SmokeOctreeGrid2D::_selectInfos_for_CUDAOpenPoissonRedBlack2(
	const int width, const int height, const bool* leafLevel0, const bool* leafLevel1, const bool* leafLevel2, const bool* leafLevel3, std::vector<Info>& lambdainfos,
	std::vector<int>& level0_index_red, std::vector<int>& level0_index_black, 
	std::vector<int>& level1_index_red, std::vector<int>& level1_index_black, 
	std::vector<int>& level2_index_red, std::vector<int>& level2_index_black, 
	std::vector<int>& level3_index_red, std::vector<int>& level3_index_black)
{
	level0_index_red.clear();
	level1_index_red.clear();
	level2_index_red.clear();
	level3_index_red.clear();
	level0_index_black.clear();
	level1_index_black.clear();
	level2_index_black.clear();
	level3_index_black.clear();

	int levelWidth1 = width/2;
	int levelHeight1 = height/2;
	int levelWidth2 = width/4;
	int levelHeight2 = height/4;
	int levelWidth3 = width/8;
	int levelHeight3 = height/8;

	_selectInfos_lambda(width,height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos);

	/*level3 lambda*/
	for(int y = 0;y < levelHeight3;y++)
	{
		for(int x = 0;x < levelWidth3;x++)
		{
			if(leafLevel3[y*levelWidth3+x])
			{
				if((x+y)%2 == 0)
				{
					level3_index_red.push_back(x);
					level3_index_red.push_back(y);
				}
				else
				{
					level3_index_black.push_back(x);
					level3_index_black.push_back(y);
				}
			}
		}
	}

	/*level2 lambda*/
	for(int y = 0;y < levelHeight2;y++)
	{
		for(int x = 0;x < levelWidth2;x++)
		{
			if(leafLevel2[y*levelWidth2+x])
			{
				if((x+y)%2 == 0)
				{
					level2_index_red.push_back(x);
					level2_index_red.push_back(y);
				}
				else
				{
					level2_index_black.push_back(x);
					level2_index_black.push_back(y);
				}
			}
		}
	}

	/*level1 lambda*/
	for(int y = 0;y < levelHeight1;y++)
	{
		for(int x = 0;x < levelWidth1;x++)
		{
			if(leafLevel1[y*levelWidth1+x])
			{
				if((x+y)%2 == 0)
				{
					level1_index_red.push_back(x);
					level1_index_red.push_back(y);
				}
				else
				{
					level1_index_black.push_back(x);
					level1_index_black.push_back(y);
				}
			}
		}
	}

	/*level0 lambda*/
	for(int y = 0;y < height;y++)
	{
		for(int x = 0;x < width;x++)
		{
			if(leafLevel0[y*width+x])
			{
				if((x+y)%2 == 0)
				{
					level0_index_red.push_back(x);
					level0_index_red.push_back(y);
				}
				else
				{
					level0_index_black.push_back(x);
					level0_index_black.push_back(y);
				}
			}
		}
	}
}

void ZQ_SmokeOctreeGrid2D::_selectInfos_for_CUDAOpenPoissonRedBlack3(
	const int width, const int height, const bool* leafLevel0, const bool* leafLevel1, const bool* leafLevel2, const bool* leafLevel3, std::vector<Info>& lambdainfos,
	std::vector<int>& level0_index_red, std::vector<int>& level0_neighborinfo_red, std::vector<int>& level0_index_black, std::vector<int>& level0_neighborinfo_black, 
	std::vector<int>& level1_index_red, std::vector<int>& level1_neighborinfo_red, std::vector<int>& level1_index_black, std::vector<int>& level1_neighborinfo_black, 
	std::vector<int>& level2_index_red, std::vector<int>& level2_neighborinfo_red, std::vector<int>& level2_index_black, std::vector<int>& level2_neighborinfo_black, 
	std::vector<int>& level3_index_red, std::vector<int>& level3_neighborinfo_red, std::vector<int>& level3_index_black, std::vector<int>& level3_neighborinfo_black )
{
	level0_index_red.clear();
	level0_index_black.clear();
	level1_index_red.clear();
	level1_index_black.clear();
	level2_index_red.clear();
	level2_index_black.clear();
	level3_index_red.clear();
	level3_index_black.clear();
	level0_neighborinfo_red.clear();
	level0_neighborinfo_black.clear();
	level1_neighborinfo_red.clear();
	level1_neighborinfo_black.clear();
	level2_neighborinfo_red.clear();
	level2_neighborinfo_black.clear();
	level3_neighborinfo_red.clear();
	level3_neighborinfo_black.clear();
	
	
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;

	_selectInfos_lambda(width,height,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos);

	/****************  for level0  ******************************/
	int level0_red_offset = 0;
	int level0_black_offset = 0;
	for(int y = 0;y < height;y++)
	{
		for(int x = 0;x < width;x++)
		{
			if(!leafLevel0[y*width+x])
				continue;
			if((x+y)%2 == 0)
			{
				level0_index_red.push_back(x);
				level0_index_red.push_back(y);
				int cur_offset = level0_red_offset;
				
				//x-
				if(x == 0)
				{
					level0_neighborinfo_red.push_back(3);
					level0_neighborinfo_red.push_back(x/8-1);
					level0_neighborinfo_red.push_back(y/8);
					level0_red_offset ++;
				}
				else
				{
					if(leafLevel0[y*width+x-1])
					{
						level0_neighborinfo_red.push_back(0);
						level0_neighborinfo_red.push_back(x-1);
						level0_neighborinfo_red.push_back(y);
						level0_red_offset ++;
					}
					else if(x%2 == 0 && leafLevel1[(y/2)*levelWidth1+(x/2)-1])
					{
						level0_neighborinfo_red.push_back(1);
						level0_neighborinfo_red.push_back(x/2-1);
						level0_neighborinfo_red.push_back(y/2);
						level0_red_offset ++;
					}
					else if(x%4 == 0 && leafLevel2[(y/4)*levelWidth2+(x/4)-1])
					{
						level0_neighborinfo_red.push_back(2);
						level0_neighborinfo_red.push_back(x/4-1);
						level0_neighborinfo_red.push_back(y/4);
						level0_red_offset ++;
					}
					else if(x%8 == 0 && leafLevel3[(y/8)*levelWidth3+(x/8)-1])
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8-1);
						level0_neighborinfo_red.push_back(y/8);
						level0_red_offset ++;
					}
				}
				
				//x+
				if(x == width-1)
				{
					level0_neighborinfo_red.push_back(3);
					level0_neighborinfo_red.push_back((x+1)/8);
					level0_neighborinfo_red.push_back(y/8);
					level0_red_offset ++;
				}
				else
				{
					if(leafLevel0[y*width+(x+1)])
					{
						level0_neighborinfo_red.push_back(0);
						level0_neighborinfo_red.push_back(x+1);
						level0_neighborinfo_red.push_back(y);
						level0_red_offset ++;
					}
					else if((x+1)%2 == 0 && leafLevel1[(y/2)*levelWidth1+(x+1)/2])
					{
						level0_neighborinfo_red.push_back(1);
						level0_neighborinfo_red.push_back((x+1)/2);
						level0_neighborinfo_red.push_back(y/2);
						level0_red_offset ++;
					}
					else if((x+1)%4 == 0 && leafLevel2[(y/4)*levelWidth2+(x+1)/4])
					{
						level0_neighborinfo_red.push_back(2);
						level0_neighborinfo_red.push_back((x+1)/4);
						level0_neighborinfo_red.push_back(y/4);
						level0_red_offset ++;
					}
					else if((x+1)%8 == 0 && leafLevel3[(y/8)*levelWidth3+(x+1)/8])
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back((x+1)/8);
						level0_neighborinfo_red.push_back(y/8);
						level0_red_offset ++;
					}
				}

				//y-
				if(y == 0)
				{
					level0_neighborinfo_red.push_back(3);
					level0_neighborinfo_red.push_back(x/8);
					level0_neighborinfo_red.push_back(y/8-1);
					level0_red_offset ++;
				}
				else
				{
					if(leafLevel0[(y-1)*width+x])
					{
						level0_neighborinfo_red.push_back(0);
						level0_neighborinfo_red.push_back(x);
						level0_neighborinfo_red.push_back(y-1);
						level0_red_offset ++;
					}
					else if(y%2 == 0 && leafLevel1[(y/2-1)*levelWidth1+x/2])
					{
						level0_neighborinfo_red.push_back(1);
						level0_neighborinfo_red.push_back(x/2);
						level0_neighborinfo_red.push_back(y/2-1);
						level0_red_offset ++;
					}
					else if(y%4 == 0 && leafLevel2[(y/4-1)*levelWidth2+x/4])
					{
						level0_neighborinfo_red.push_back(2);
						level0_neighborinfo_red.push_back(x/4);
						level0_neighborinfo_red.push_back(y/4-1);
						level0_red_offset ++;
					}
					else if(y%8 == 0 && leafLevel3[(y/8-1)*levelWidth3+x/8])
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8);
						level0_neighborinfo_red.push_back(y/8-1);
						level0_red_offset ++;
					}
				}

				//y+
				if(y == height-1)
				{
					level0_neighborinfo_red.push_back(3);
					level0_neighborinfo_red.push_back(x/8);
					level0_neighborinfo_red.push_back((y+1)/8);
					level0_red_offset ++;
				}
				else
				{
					if(leafLevel0[(y+1)*width+x])
					{
						level0_neighborinfo_red.push_back(0);
						level0_neighborinfo_red.push_back(x);
						level0_neighborinfo_red.push_back(y+1);
						level0_red_offset ++;
					}
					else if((y+1)%2 == 0 && leafLevel1[(y+1)/2*levelWidth1+x/2])
					{
						level0_neighborinfo_red.push_back(1);
						level0_neighborinfo_red.push_back(x/2);
						level0_neighborinfo_red.push_back((y+1)/2);
						level0_red_offset ++;
					}
					else if((y+1)%4 == 0 && leafLevel2[(y+1)/4*levelWidth2+x/4])
					{
						level0_neighborinfo_red.push_back(2);
						level0_neighborinfo_red.push_back(x/4);
						level0_neighborinfo_red.push_back((y+1)/4);
						level0_red_offset ++;
					}
					else if((y+1)%8 == 0 && leafLevel3[(y+1)/8*levelWidth3+x/8])
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8);
						level0_neighborinfo_red.push_back((y+1)/8);
						level0_red_offset ++;
					}
				}

				level0_index_red.push_back(level0_red_offset-cur_offset);
				level0_index_red.push_back(cur_offset);
			}
			else
			{
				level0_index_black.push_back(x);
				level0_index_black.push_back(y);
				int cur_offset = level0_black_offset;

				//x-
				if(x == 0)
				{
					level0_neighborinfo_black.push_back(3);
					level0_neighborinfo_black.push_back(x/8-1);
					level0_neighborinfo_black.push_back(y/8);
					level0_black_offset ++;
				}
				else
				{
					if(leafLevel0[y*width+x-1])
					{
						level0_neighborinfo_black.push_back(0);
						level0_neighborinfo_black.push_back(x-1);
						level0_neighborinfo_black.push_back(y);
						level0_black_offset ++;
					}
					else if(x%2 == 0 && leafLevel1[(y/2)*levelWidth1+(x/2)-1])
					{
						level0_neighborinfo_black.push_back(1);
						level0_neighborinfo_black.push_back(x/2-1);
						level0_neighborinfo_black.push_back(y/2);
						level0_black_offset ++;
					}
					else if(x%4 == 0 && leafLevel2[(y/4)*levelWidth2+(x/4)-1])
					{
						level0_neighborinfo_black.push_back(2);
						level0_neighborinfo_black.push_back(x/4-1);
						level0_neighborinfo_black.push_back(y/4);
						level0_black_offset ++;
					}
					else if(x%8 == 0 && leafLevel3[(y/8)*levelWidth3+(x/8)-1])
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8-1);
						level0_neighborinfo_black.push_back(y/8);
						level0_black_offset ++;
					}
				}

				//x+
				if(x == width-1)
				{
					level0_neighborinfo_black.push_back(3);
					level0_neighborinfo_black.push_back((x+1)/8);
					level0_neighborinfo_black.push_back(y/8);
					level0_black_offset ++;
				}
				else
				{
					if(leafLevel0[y*width+(x+1)])
					{
						level0_neighborinfo_black.push_back(0);
						level0_neighborinfo_black.push_back(x+1);
						level0_neighborinfo_black.push_back(y);
						level0_black_offset ++;
					}
					else if((x+1)%2 == 0 && leafLevel1[(y/2)*levelWidth1+(x+1)/2])
					{
						level0_neighborinfo_black.push_back(1);
						level0_neighborinfo_black.push_back((x+1)/2);
						level0_neighborinfo_black.push_back(y/2);
						level0_black_offset ++;
					}
					else if((x+1)%4 == 0 && leafLevel2[(y/4)*levelWidth2+(x+1)/4])
					{
						level0_neighborinfo_black.push_back(2);
						level0_neighborinfo_black.push_back((x+1)/4);
						level0_neighborinfo_black.push_back(y/4);
						level0_black_offset ++;
					}
					else if((x+1)%8 == 0 && leafLevel3[(y/8)*levelWidth3+(x+1)/8])
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back((x+1)/8);
						level0_neighborinfo_black.push_back(y/8);
						level0_black_offset ++;
					}
				}

				//y-
				if(y == 0)
				{
					level0_neighborinfo_black.push_back(3);
					level0_neighborinfo_black.push_back(x/8);
					level0_neighborinfo_black.push_back(y/8-1);
					level0_black_offset ++;
				}
				else
				{
					if(leafLevel0[(y-1)*width+x])
					{
						level0_neighborinfo_black.push_back(0);
						level0_neighborinfo_black.push_back(x);
						level0_neighborinfo_black.push_back(y-1);
						level0_black_offset ++;
					}
					else if(y%2 == 0 && leafLevel1[(y/2-1)*levelWidth1+x/2])
					{
						level0_neighborinfo_black.push_back(1);
						level0_neighborinfo_black.push_back(x/2);
						level0_neighborinfo_black.push_back(y/2-1);
						level0_black_offset ++;
					}
					else if(y%4 == 0 && leafLevel2[(y/4-1)*levelWidth2+x/4])
					{
						level0_neighborinfo_black.push_back(2);
						level0_neighborinfo_black.push_back(x/4);
						level0_neighborinfo_black.push_back(y/4-1);
						level0_black_offset ++;
					}
					else if(y%8 == 0 && leafLevel3[(y/8-1)*levelWidth3+x/8])
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8);
						level0_neighborinfo_black.push_back(y/8-1);
						level0_black_offset ++;
					}
				}

				//y+
				if(y == height-1)
				{
					level0_neighborinfo_black.push_back(3);
					level0_neighborinfo_black.push_back(x/8);
					level0_neighborinfo_black.push_back((y+1)/8);
					level0_black_offset ++;
				}
				else
				{
					if(leafLevel0[(y+1)*width+x])
					{
						level0_neighborinfo_black.push_back(0);
						level0_neighborinfo_black.push_back(x);
						level0_neighborinfo_black.push_back(y+1);
						level0_black_offset ++;
					}
					else if((y+1)%2 == 0 && leafLevel1[(y+1)/2*levelWidth1+x/2])
					{
						level0_neighborinfo_black.push_back(1);
						level0_neighborinfo_black.push_back(x/2);
						level0_neighborinfo_black.push_back((y+1)/2);
						level0_black_offset ++;
					}
					else if((y+1)%4 == 0 && leafLevel2[(y+1)/4*levelWidth2+x/4])
					{
						level0_neighborinfo_black.push_back(2);
						level0_neighborinfo_black.push_back(x/4);
						level0_neighborinfo_black.push_back((y+1)/4);
						level0_black_offset ++;
					}
					else if((y+1)%8 == 0 && leafLevel3[(y+1)/8*levelWidth3+x/8])
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8);
						level0_neighborinfo_black.push_back((y+1)/8);
						level0_black_offset ++;
					}
				}

				level0_index_black.push_back(level0_black_offset-cur_offset);
				level0_index_black.push_back(cur_offset);
			}
		}
	}

	/************************** for level1 **************************************/
	int level1_red_offset = 0;
	int level1_black_offset = 0;
	for(int y = 0;y < levelHeight1;y++)
	{
		for(int x = 0;x < levelWidth1;x++)
		{
			if(!leafLevel1[y*levelWidth1+x])
				continue;

			if((x+y)%2 == 0)
			{
				level1_index_red.push_back(x);
				level1_index_red.push_back(y);
				int cur_offset = level1_red_offset;

				//x-
				if(x == 0)
				{
					level1_neighborinfo_red.push_back(3);
					level1_neighborinfo_red.push_back(x/4-1);
					level1_neighborinfo_red.push_back(y/4);
					level1_red_offset ++;
				}
				else
				{
					if(leafLevel1[y*levelWidth1+x-1])
					{
						level1_neighborinfo_red.push_back(1);
						level1_neighborinfo_red.push_back(x-1);
						level1_neighborinfo_red.push_back(y);
						level1_red_offset ++;
					}
					else if(x%2 == 0 && leafLevel2[(y/2)*levelWidth2+x/2-1])
					{
						level1_neighborinfo_red.push_back(2);
						level1_neighborinfo_red.push_back(x/2-1);
						level1_neighborinfo_red.push_back(y/2);
						level1_red_offset ++;
					}
					else if(x%4 == 0 && leafLevel3[(y/4)*levelWidth3+x/4-1])
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4-1);
						level1_neighborinfo_red.push_back(y/4);
						level1_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y*2)*width+x*2-1])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back(x*2-1);
							level1_neighborinfo_red.push_back(y*2);
							level1_red_offset ++;
						}
						if(leafLevel0[(y*2+1)*width+x*2-1])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back(x*2-1);
							level1_neighborinfo_red.push_back(y*2+1);
							level1_red_offset ++;
						}
					}
				}

				//x+
				if(x == levelWidth1-1)
				{
					level1_neighborinfo_red.push_back(3);
					level1_neighborinfo_red.push_back((x+1)/4);
					level1_neighborinfo_red.push_back(y/4);
					level1_red_offset ++;
				}
				else
				{
					if(leafLevel1[y*levelWidth1+x+1])
					{
						level1_neighborinfo_red.push_back(1);
						level1_neighborinfo_red.push_back(x+1);
						level1_neighborinfo_red.push_back(y);
						level1_red_offset ++;
					}
					else if((x+1)%2 == 0 && leafLevel2[(y/2)*levelWidth2+(x+1)/2])
					{
						level1_neighborinfo_red.push_back(2);
						level1_neighborinfo_red.push_back((x+1)/2);
						level1_neighborinfo_red.push_back(y/2);
						level1_red_offset ++;
					}
					else if((x+1)%4 == 0 && leafLevel3[(y/4)*levelWidth3+(x+1)/4])
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back((x+1)/4);
						level1_neighborinfo_red.push_back(y/4);
						level1_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y*2)*width+(x+1)*2])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back((x+1)*2);
							level1_neighborinfo_red.push_back(y*2);
							level1_red_offset ++;
						}
						if(leafLevel0[(y*2+1)*width+(x+1)*2])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back((x+1)*2);
							level1_neighborinfo_red.push_back(y*2+1);
							level1_red_offset ++;
						}
					}
				}
				
				//y-
				if(y == 0)
				{
					level1_neighborinfo_red.push_back(3);
					level1_neighborinfo_red.push_back(x/4);
					level1_neighborinfo_red.push_back(y/4-1);
					level1_red_offset ++;
				}
				else 
				{
					if(leafLevel1[(y-1)*levelWidth1+x])
					{
						level1_neighborinfo_red.push_back(1);
						level1_neighborinfo_red.push_back(x);
						level1_neighborinfo_red.push_back(y-1);
						level1_red_offset ++;
					}
					else if(y%2 == 0 && leafLevel2[(y/2-1)*levelWidth2+x/2])
					{
						level1_neighborinfo_red.push_back(2);
						level1_neighborinfo_red.push_back(x/2);
						level1_neighborinfo_red.push_back(y/2-1);
						level1_red_offset ++;
					}
					else if(y%4 == 0 && leafLevel3[(y/4-1)*levelWidth3+x/4])
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4);
						level1_neighborinfo_red.push_back(y/4-1);
						level1_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y*2-1)*width+x*2])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back(x*2);
							level1_neighborinfo_red.push_back(y*2-1);
							level1_red_offset ++;
						}
						if(leafLevel0[(y*2-1)*width+x*2+1])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back(x*2+1);
							level1_neighborinfo_red.push_back(y*2-1);
							level1_red_offset ++;
						}
					}
				}

				//y+
				if(y == levelHeight1-1)
				{
					level1_neighborinfo_red.push_back(3);
					level1_neighborinfo_red.push_back(x/4);
					level1_neighborinfo_red.push_back((y+1)/4);
					level1_red_offset ++;
				}
				else
				{
					if(leafLevel1[(y+1)*levelWidth1+x])
					{
						level1_neighborinfo_red.push_back(1);
						level1_neighborinfo_red.push_back(x);
						level1_neighborinfo_red.push_back(y+1);
						level1_red_offset ++;
					}
					else if((y+1)%2 == 0 && leafLevel2[(y+1)/2*levelWidth2+x/2])
					{
						level1_neighborinfo_red.push_back(2);
						level1_neighborinfo_red.push_back(x/2);
						level1_neighborinfo_red.push_back((y+1)/2);
						level1_red_offset ++;
					}
					else if((y+1)%4 == 0 && leafLevel3[(y+1)/4*levelWidth3+x/4])
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4);
						level1_neighborinfo_red.push_back((y+1)/4);
						level1_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y+1)*2*width+x*2])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back(x*2);
							level1_neighborinfo_red.push_back((y+1)*2);
							level1_red_offset ++;
						}
						if(leafLevel0[(y+1)*2*width+x*2+1])
						{
							level1_neighborinfo_red.push_back(0);
							level1_neighborinfo_red.push_back(x*2+1);
							level1_neighborinfo_red.push_back((y+1)*2);
							level1_red_offset ++;
						}
					}
				}
				level1_index_red.push_back(level1_red_offset-cur_offset);
				level1_index_red.push_back(cur_offset);
			}
			else
			{
				level1_index_black.push_back(x);
				level1_index_black.push_back(y);
				int cur_offset = level1_black_offset;

				//x-
				if(x == 0)
				{
					level1_neighborinfo_black.push_back(3);
					level1_neighborinfo_black.push_back(x/4-1);
					level1_neighborinfo_black.push_back(y/4);
					level1_black_offset ++;
				}
				else
				{
					if(leafLevel1[y*levelWidth1+x-1])
					{
						level1_neighborinfo_black.push_back(1);
						level1_neighborinfo_black.push_back(x-1);
						level1_neighborinfo_black.push_back(y);
						level1_black_offset ++;
					}
					else if(x%2 == 0 && leafLevel2[(y/2)*levelWidth2+x/2-1])
					{
						level1_neighborinfo_black.push_back(2);
						level1_neighborinfo_black.push_back(x/2-1);
						level1_neighborinfo_black.push_back(y/2);
						level1_black_offset ++;
					}
					else if(x%4 == 0 && leafLevel3[(y/4)*levelWidth3+x/4-1])
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4-1);
						level1_neighborinfo_black.push_back(y/4);
						level1_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y*2)*width+x*2-1])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back(x*2-1);
							level1_neighborinfo_black.push_back(y*2);
							level1_black_offset ++;
						}
						if(leafLevel0[(y*2+1)*width+x*2-1])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back(x*2-1);
							level1_neighborinfo_black.push_back(y*2+1);
							level1_black_offset ++;
						}
					}
				}

				//x+
				if(x == levelWidth1-1)
				{
					level1_neighborinfo_black.push_back(3);
					level1_neighborinfo_black.push_back((x+1)/4);
					level1_neighborinfo_black.push_back(y/4);
					level1_black_offset ++;
				}
				else
				{
					if(leafLevel1[y*levelWidth1+x+1])
					{
						level1_neighborinfo_black.push_back(1);
						level1_neighborinfo_black.push_back(x+1);
						level1_neighborinfo_black.push_back(y);
						level1_black_offset ++;
					}
					else if((x+1)%2 == 0 && leafLevel2[(y/2)*levelWidth2+(x+1)/2])
					{
						level1_neighborinfo_black.push_back(2);
						level1_neighborinfo_black.push_back((x+1)/2);
						level1_neighborinfo_black.push_back(y/2);
						level1_black_offset ++;
					}
					else if((x+1)%4 == 0 && leafLevel3[(y/4)*levelWidth3+(x+1)/4])
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back((x+1)/4);
						level1_neighborinfo_black.push_back(y/4);
						level1_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y*2)*width+(x+1)*2])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back((x+1)*2);
							level1_neighborinfo_black.push_back(y*2);
							level1_black_offset ++;
						}
						if(leafLevel0[(y*2+1)*width+(x+1)*2])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back((x+1)*2);
							level1_neighborinfo_black.push_back(y*2+1);
							level1_black_offset ++;
						}
					}
				}

				//y-
				if(y == 0)
				{
					level1_neighborinfo_black.push_back(3);
					level1_neighborinfo_black.push_back(x/4);
					level1_neighborinfo_black.push_back(y/4-1);
					level1_black_offset ++;
				}
				else 
				{
					if(leafLevel1[(y-1)*levelWidth1+x])
					{
						level1_neighborinfo_black.push_back(1);
						level1_neighborinfo_black.push_back(x);
						level1_neighborinfo_black.push_back(y-1);
						level1_black_offset ++;
					}
					else if(y%2 == 0 && leafLevel2[(y/2-1)*levelWidth2+x/2])
					{
						level1_neighborinfo_black.push_back(2);
						level1_neighborinfo_black.push_back(x/2);
						level1_neighborinfo_black.push_back(y/2-1);
						level1_black_offset ++;
					}
					else if(y%4 == 0 && leafLevel3[(y/4-1)*levelWidth3+x/4])
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4);
						level1_neighborinfo_black.push_back(y/4-1);
						level1_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y*2-1)*width+x*2])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back(x*2);
							level1_neighborinfo_black.push_back(y*2-1);
							level1_black_offset ++;
						}
						if(leafLevel0[(y*2-1)*width+x*2+1])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back(x*2+1);
							level1_neighborinfo_black.push_back(y*2-1);
							level1_black_offset ++;
						}
					}
				}

				//y+
				if(y == levelHeight1-1)
				{
					level1_neighborinfo_black.push_back(3);
					level1_neighborinfo_black.push_back(x/4);
					level1_neighborinfo_black.push_back((y+1)/4);
					level1_black_offset ++;
				}
				else
				{
					if(leafLevel1[(y+1)*levelWidth1+x])
					{
						level1_neighborinfo_black.push_back(1);
						level1_neighborinfo_black.push_back(x);
						level1_neighborinfo_black.push_back(y+1);
						level1_black_offset ++;
					}
					else if((y+1)%2 == 0 && leafLevel2[(y+1)/2*levelWidth2+x/2])
					{
						level1_neighborinfo_black.push_back(2);
						level1_neighborinfo_black.push_back(x/2);
						level1_neighborinfo_black.push_back((y+1)/2);
						level1_black_offset ++;
					}
					else if((y+1)%4 == 0 && leafLevel3[(y+1)/4*levelWidth3+x/4])
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4);
						level1_neighborinfo_black.push_back((y+1)/4);
						level1_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel0[(y+1)*2*width+x*2])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back(x*2);
							level1_neighborinfo_black.push_back((y+1)*2);
							level1_black_offset ++;
						}
						if(leafLevel0[(y+1)*2*width+x*2+1])
						{
							level1_neighborinfo_black.push_back(0);
							level1_neighborinfo_black.push_back(x*2+1);
							level1_neighborinfo_black.push_back((y+1)*2);
							level1_black_offset ++;
						}
					}
				}
				level1_index_black.push_back(level1_black_offset-cur_offset);
				level1_index_black.push_back(cur_offset);
			}
		}
	}

	/************************** for level2 ***************************/
	int level2_red_offset = 0;
	int level2_black_offset = 0;
	for(int y = 0;y < levelHeight2;y++)
	{
		for(int x = 0;x < levelWidth2;x++)
		{
			if(!leafLevel2[y*levelWidth2+x])
				continue;

			if((x+y)%2 == 0)
			{
				level2_index_red.push_back(x);
				level2_index_red.push_back(y);
				int cur_offset = level2_red_offset;

				//x-
				if(x == 0)
				{
					level2_neighborinfo_red.push_back(3);
					level2_neighborinfo_red.push_back(x/2-1);
					level2_neighborinfo_red.push_back(y/2);
					level2_red_offset ++;
				}
				else
				{
					if(leafLevel2[y*levelWidth2+x-1])
					{
						level2_neighborinfo_red.push_back(2);
						level2_neighborinfo_red.push_back(x-1);
						level2_neighborinfo_red.push_back(y);
						level2_red_offset ++;
					}
					else if(x%2 == 0 && leafLevel3[(y/2)*levelWidth3+x/2-1])
					{
						level2_neighborinfo_red.push_back(3);
						level2_neighborinfo_red.push_back(x/2-1);
						level2_neighborinfo_red.push_back(y/2);
						level2_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y*2)*levelWidth1+x*2-1])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back(x*2-1);
							level2_neighborinfo_red.push_back(y*2);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4)*width+x*4-1])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4-1);
								level2_neighborinfo_red.push_back(y*4);
								level2_red_offset ++;
							}
							if(leafLevel0[(y*4+1)*width+x*4-1])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4-1);
								level2_neighborinfo_red.push_back(y*4+1);
								level2_red_offset ++;
							}
						}

						if(leafLevel1[(y*2+1)*levelWidth1+x*2-1])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back(x*2-1);
							level2_neighborinfo_red.push_back(y*2+1);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4+2)*width+x*4-1])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4-1);
								level2_neighborinfo_red.push_back(y*4+2);
								level2_red_offset ++;
							}
							if(leafLevel0[(y*4+3)*width+x*4-1])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4-1);
								level2_neighborinfo_red.push_back(y*4+3);
								level2_red_offset ++;
							}
						}
					}
				}

				//x+
				if(x == levelWidth2-1)
				{
					level2_neighborinfo_red.push_back(3);
					level2_neighborinfo_red.push_back((x+1)/2);
					level2_neighborinfo_red.push_back(y/2);
					level2_red_offset ++;
				}
				else
				{
					if(leafLevel2[y*levelWidth2+x+1])
					{
						level2_neighborinfo_red.push_back(2);
						level2_neighborinfo_red.push_back(x+1);
						level2_neighborinfo_red.push_back(y);
						level2_red_offset ++;
					}
					else if((x+1)%2 == 0 && leafLevel3[(y/2)*levelWidth3+(x+1)/2])
					{
						level2_neighborinfo_red.push_back(3);
						level2_neighborinfo_red.push_back((x+1)/2);
						level2_neighborinfo_red.push_back(y/2);
						level2_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y*2)*levelWidth1+(x+1)*2])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back((x+1)*2);
							level2_neighborinfo_red.push_back(y*2);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4)*width+(x+1)*4])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back((x+1)*4);
								level2_neighborinfo_red.push_back(y*4);
								level2_red_offset ++;
							}
							if(leafLevel0[(y*4+1)*width+(x+1)*4])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back((x+1)*4);
								level2_neighborinfo_red.push_back(y*4+1);
								level2_red_offset ++;
							}
						}

						if(leafLevel1[(y*2+1)*levelWidth1+(x+1)*2])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back((x+1)*2);
							level2_neighborinfo_red.push_back(y*2+1);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4+2)*width+(x+1)*4])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back((x+1)*4);
								level2_neighborinfo_red.push_back(y*4+2);
								level2_red_offset ++;
							}
							if(leafLevel0[(y*4+3)*width+(x+1)*4])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back((x+1)*4);
								level2_neighborinfo_red.push_back(y*4+3);
								level2_red_offset ++;
							}
						}
					}
				}

				//y-
				if(y == 0)
				{
					level2_neighborinfo_red.push_back(3);
					level2_neighborinfo_red.push_back(x/2);
					level2_neighborinfo_red.push_back(y/2-1);
					level2_red_offset ++;
				}
				else
				{
					if(leafLevel2[(y-1)*levelWidth2+x])
					{
						level2_neighborinfo_red.push_back(2);
						level2_neighborinfo_red.push_back(x);
						level2_neighborinfo_red.push_back(y-1);
						level2_red_offset ++;
					}
					else if(y%2 == 0 && leafLevel3[(y/2-1)*levelWidth3+x/2])
					{
						level2_neighborinfo_red.push_back(3);
						level2_neighborinfo_red.push_back(x/2);
						level2_neighborinfo_red.push_back(y/2-1);
						level2_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y*2-1)*levelWidth1+x*2])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back(x*2);
							level2_neighborinfo_red.push_back(y*2-1);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4-1)*width+x*4])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4);
								level2_neighborinfo_red.push_back(y*4-1);
								level2_red_offset ++;
							}
							if(leafLevel0[(y*4-1)*width+x*4+1])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4+1);
								level2_neighborinfo_red.push_back(y*4-1);
								level2_red_offset ++;
							}
						}

						if(leafLevel1[(y*2-1)*levelWidth1+x*2+1])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back(x*2+1);
							level2_neighborinfo_red.push_back(y*2-1);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4-1)*width+x*4+2])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4+2);
								level2_neighborinfo_red.push_back(y*4-1);
								level2_red_offset ++;
							}
							if(leafLevel0[(y*4-1)*width+x*4+3])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4+3);
								level2_neighborinfo_red.push_back(y*4-1);
								level2_red_offset ++;
							}
						}
					}
				}

				//y+
				if(y == levelHeight2-1)
				{
					level2_neighborinfo_red.push_back(3);
					level2_neighborinfo_red.push_back(x/2);
					level2_neighborinfo_red.push_back((y+1)/2);
					level2_red_offset ++;
				}
				else
				{
					if(leafLevel2[(y+1)*levelWidth2+x])
					{
						level2_neighborinfo_red.push_back(2);
						level2_neighborinfo_red.push_back(x);
						level2_neighborinfo_red.push_back(y+1);
						level2_red_offset ++;
					}
					else if((y+1)%2 == 0 && leafLevel3[(y+1)/2*levelWidth3+x/2])
					{
						level2_neighborinfo_red.push_back(3);
						level2_neighborinfo_red.push_back(x/2);
						level2_neighborinfo_red.push_back((y+1)/2);
						level2_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y+1)*2*levelWidth1+x*2])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back(x*2);
							level2_neighborinfo_red.push_back((y+1)*2);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y+1)*4*width+x*4])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4);
								level2_neighborinfo_red.push_back((y+1)*4);
								level2_red_offset ++;
							}
							if(leafLevel0[(y+1)*4*width+x*4+1])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4+1);
								level2_neighborinfo_red.push_back((y+1)*4);
								level2_red_offset ++;
							}
						}

						if(leafLevel1[(y+1)*2*levelWidth1+x*2+1])
						{
							level2_neighborinfo_red.push_back(1);
							level2_neighborinfo_red.push_back(x*2+1);
							level2_neighborinfo_red.push_back((y+1)*2);
							level2_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y+1)*4*width+x*4+2])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4+2);
								level2_neighborinfo_red.push_back((y+1)*4);
								level2_red_offset ++;
							}
							if(leafLevel0[(y+1)*4*width+x*4+3])
							{
								level2_neighborinfo_red.push_back(0);
								level2_neighborinfo_red.push_back(x*4+3);
								level2_neighborinfo_red.push_back((y+1)*4);
								level2_red_offset ++;
							}
						}
					}
				}
				level2_index_red.push_back(level2_red_offset-cur_offset);
				level2_index_red.push_back(cur_offset);
			}
			else
			{
				level2_index_black.push_back(x);
				level2_index_black.push_back(y);
				int cur_offset = level2_black_offset;

				//x-
				if(x == 0)
				{
					level2_neighborinfo_black.push_back(3);
					level2_neighborinfo_black.push_back(x/2-1);
					level2_neighborinfo_black.push_back(y/2);
					level2_black_offset ++;
				}
				else
				{
					if(leafLevel2[y*levelWidth2+x-1])
					{
						level2_neighborinfo_black.push_back(2);
						level2_neighborinfo_black.push_back(x-1);
						level2_neighborinfo_black.push_back(y);
						level2_black_offset ++;
					}
					else if(x%2 == 0 && leafLevel3[(y/2)*levelWidth3+x/2-1])
					{
						level2_neighborinfo_black.push_back(3);
						level2_neighborinfo_black.push_back(x/2-1);
						level2_neighborinfo_black.push_back(y/2);
						level2_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y*2)*levelWidth1+x*2-1])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back(x*2-1);
							level2_neighborinfo_black.push_back(y*2);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4)*width+x*4-1])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4-1);
								level2_neighborinfo_black.push_back(y*4);
								level2_black_offset ++;
							}
							if(leafLevel0[(y*4+1)*width+x*4-1])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4-1);
								level2_neighborinfo_black.push_back(y*4+1);
								level2_black_offset ++;
							}
						}

						if(leafLevel1[(y*2+1)*levelWidth1+x*2-1])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back(x*2-1);
							level2_neighborinfo_black.push_back(y*2+1);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4+2)*width+x*4-1])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4-1);
								level2_neighborinfo_black.push_back(y*4+2);
								level2_black_offset ++;
							}
							if(leafLevel0[(y*4+3)*width+x*4-1])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4-1);
								level2_neighborinfo_black.push_back(y*4+3);
								level2_black_offset ++;
							}
						}
					}
				}

				//x+
				if(x == levelWidth2-1)
				{
					level2_neighborinfo_black.push_back(3);
					level2_neighborinfo_black.push_back((x+1)/2);
					level2_neighborinfo_black.push_back(y/2);
					level2_black_offset ++;
				}
				else
				{
					if(leafLevel2[y*levelWidth2+x+1])
					{
						level2_neighborinfo_black.push_back(2);
						level2_neighborinfo_black.push_back(x+1);
						level2_neighborinfo_black.push_back(y);
						level2_black_offset ++;
					}
					else if((x+1)%2 == 0 && leafLevel3[(y/2)*levelWidth3+(x+1)/2])
					{
						level2_neighborinfo_black.push_back(3);
						level2_neighborinfo_black.push_back((x+1)/2);
						level2_neighborinfo_black.push_back(y/2);
						level2_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y*2)*levelWidth1+(x+1)*2])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back((x+1)*2);
							level2_neighborinfo_black.push_back(y*2);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4)*width+(x+1)*4])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back((x+1)*4);
								level2_neighborinfo_black.push_back(y*4);
								level2_black_offset ++;
							}
							if(leafLevel0[(y*4+1)*width+(x+1)*4])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back((x+1)*4);
								level2_neighborinfo_black.push_back(y*4+1);
								level2_black_offset ++;
							}
						}

						if(leafLevel1[(y*2+1)*levelWidth1+(x+1)*2])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back((x+1)*2);
							level2_neighborinfo_black.push_back(y*2+1);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4+2)*width+(x+1)*4])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back((x+1)*4);
								level2_neighborinfo_black.push_back(y*4+2);
								level2_black_offset ++;
							}
							if(leafLevel0[(y*4+3)*width+(x+1)*4])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back((x+1)*4);
								level2_neighborinfo_black.push_back(y*4+3);
								level2_black_offset ++;
							}
						}
					}
				}

				//y-
				if(y == 0)
				{
					level2_neighborinfo_black.push_back(3);
					level2_neighborinfo_black.push_back(x/2);
					level2_neighborinfo_black.push_back(y/2-1);
					level2_black_offset ++;
				}
				else
				{
					if(leafLevel2[(y-1)*levelWidth2+x])
					{
						level2_neighborinfo_black.push_back(2);
						level2_neighborinfo_black.push_back(x);
						level2_neighborinfo_black.push_back(y-1);
						level2_black_offset ++;
					}
					else if(y%2 == 0 && leafLevel3[(y/2-1)*levelWidth3+x/2])
					{
						level2_neighborinfo_black.push_back(3);
						level2_neighborinfo_black.push_back(x/2);
						level2_neighborinfo_black.push_back(y/2-1);
						level2_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y*2-1)*levelWidth1+x*2])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back(x*2);
							level2_neighborinfo_black.push_back(y*2-1);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4-1)*width+x*4])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4);
								level2_neighborinfo_black.push_back(y*4-1);
								level2_black_offset ++;
							}
							if(leafLevel0[(y*4-1)*width+x*4+1])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4+1);
								level2_neighborinfo_black.push_back(y*4-1);
								level2_black_offset ++;
							}
						}

						if(leafLevel1[(y*2-1)*levelWidth1+x*2+1])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back(x*2+1);
							level2_neighborinfo_black.push_back(y*2-1);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y*4-1)*width+x*4+2])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4+2);
								level2_neighborinfo_black.push_back(y*4-1);
								level2_black_offset ++;
							}
							if(leafLevel0[(y*4-1)*width+x*4+3])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4+3);
								level2_neighborinfo_black.push_back(y*4-1);
								level2_black_offset ++;
							}
						}
					}
				}

				//y+
				if(y == levelHeight2-1)
				{
					level2_neighborinfo_black.push_back(3);
					level2_neighborinfo_black.push_back(x/2);
					level2_neighborinfo_black.push_back((y+1)/2);
					level2_black_offset ++;
				}
				else
				{
					if(leafLevel2[(y+1)*levelWidth2+x])
					{
						level2_neighborinfo_black.push_back(2);
						level2_neighborinfo_black.push_back(x);
						level2_neighborinfo_black.push_back(y+1);
						level2_black_offset ++;
					}
					else if((y+1)%2 == 0 && leafLevel3[(y+1)/2*levelWidth3+x/2])
					{
						level2_neighborinfo_black.push_back(3);
						level2_neighborinfo_black.push_back(x/2);
						level2_neighborinfo_black.push_back((y+1)/2);
						level2_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel1[(y+1)*2*levelWidth1+x*2])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back(x*2);
							level2_neighborinfo_black.push_back((y+1)*2);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y+1)*4*width+x*4])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4);
								level2_neighborinfo_black.push_back((y+1)*4);
								level2_black_offset ++;
							}
							if(leafLevel0[(y+1)*4*width+x*4+1])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4+1);
								level2_neighborinfo_black.push_back((y+1)*4);
								level2_black_offset ++;
							}
						}

						if(leafLevel1[(y+1)*2*levelWidth1+x*2+1])
						{
							level2_neighborinfo_black.push_back(1);
							level2_neighborinfo_black.push_back(x*2+1);
							level2_neighborinfo_black.push_back((y+1)*2);
							level2_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel0[(y+1)*4*width+x*4+2])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4+2);
								level2_neighborinfo_black.push_back((y+1)*4);
								level2_black_offset ++;
							}
							if(leafLevel0[(y+1)*4*width+x*4+3])
							{
								level2_neighborinfo_black.push_back(0);
								level2_neighborinfo_black.push_back(x*4+3);
								level2_neighborinfo_black.push_back((y+1)*4);
								level2_black_offset ++;
							}
						}
					}
				}
				level2_index_black.push_back(level2_black_offset-cur_offset);
				level2_index_black.push_back(cur_offset);
			}
		}
	}

	/********************** for level3   ****************************/
	int level3_red_offset = 0;
	int level3_black_offset = 0;
	for(int y = 0;y < levelHeight3;y++)
	{
		for(int x = 0;x < levelWidth3;x++)
		{
			if(!leafLevel3[y*levelWidth3+x])
				continue;
			
			if((x+y)%2 == 0)
			{
				level3_index_red.push_back(x);
				level3_index_red.push_back(y);
				int cur_offset = level3_red_offset;

				//x-
				if(x == 0)
				{
					level3_neighborinfo_red.push_back(3);
					level3_neighborinfo_red.push_back(x-1);
					level3_neighborinfo_red.push_back(y);
					level3_red_offset ++;
				}
				else
				{
					if(leafLevel3[y*levelWidth3+x-1])
					{
						level3_neighborinfo_red.push_back(3);
						level3_neighborinfo_red.push_back(x-1);
						level3_neighborinfo_red.push_back(y);
						level3_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel2[(y*2)*levelWidth2+x*2-1])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back(x*2-1);
							level3_neighborinfo_red.push_back(y*2);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4)*levelWidth1+x*4-1])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4-1);
								level3_neighborinfo_red.push_back(y*4);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+1)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8+1);
									level3_red_offset ++;
								}
							}

							if(leafLevel1[(y*4+1)*levelWidth1+x*4-1])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4-1);
								level3_neighborinfo_red.push_back(y*4+1);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+2)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8+2);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+3)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8+3);
									level3_red_offset ++;
								}
							}
						}

						if(leafLevel2[(y*2+1)*levelWidth2+x*2-1])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back(x*2-1);
							level3_neighborinfo_red.push_back(y*2+1);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[y*4+2]*levelWidth1+x*4-1)
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4-1);
								level3_neighborinfo_red.push_back(y*4+2);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+4)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8+4);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+5)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8+5);
									level3_red_offset ++;
								}
							}

							if(leafLevel1[(y*4+3)*levelWidth1+x*4-1])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4-1);
								level3_neighborinfo_red.push_back(y*4+3);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+6)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8+6);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+7)*width+x*8-1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8-1);
									level3_neighborinfo_red.push_back(y*8+7);
									level3_red_offset ++;
								}
							}
						}
					}
				}

				//x+
				if(x == levelWidth3-1)
				{
					level3_neighborinfo_red.push_back(3);
					level3_neighborinfo_red.push_back(x+1);
					level3_neighborinfo_red.push_back(y);
					level3_red_offset ++;
				}
				else
				{
					if(leafLevel3[y*levelWidth3+x+1])
					{
						level3_neighborinfo_red.push_back(3);
						level3_neighborinfo_red.push_back(x+1);
						level3_neighborinfo_red.push_back(y);
						level3_red_offset ++;
					}
					else //smaller 
					{
						if(leafLevel2[(y*2)*levelWidth2+(x+1)*2])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back((x+1)*2);
							level3_neighborinfo_red.push_back(y*2);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back((x+1)*4);
								level3_neighborinfo_red.push_back(y*4);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+1)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8+1);
									level3_red_offset ++;
								}
							}
							if(leafLevel1[(y*4+1)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back((x+1)*4);
								level3_neighborinfo_red.push_back(y*4+1);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+2)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8+2);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+3)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8+3);
									level3_red_offset ++;
								}
							}
						}
						if(leafLevel2[(y*2+1)*levelWidth2+(x+1)*2])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back((x+1)*2);
							level3_neighborinfo_red.push_back(y*2+1);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4+2)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back((x+1)*4);
								level3_neighborinfo_red.push_back(y*4+2);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+4)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8+4);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+5)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8+5);
									level3_red_offset ++;
								}
							}
							if(leafLevel1[(y*4+3)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back((x+1)*4);
								level3_neighborinfo_red.push_back(y*4+3);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+6)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8+6);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8+7)*width+(x+1)*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back((x+1)*8);
									level3_neighborinfo_red.push_back(y*8+7);
									level3_red_offset ++;
								}
							}
						}
					}
				}

				//y-
				if(y == 0)
				{
					level3_neighborinfo_red.push_back(3);
					level3_neighborinfo_red.push_back(x);
					level3_neighborinfo_red.push_back(y-1);
					level3_red_offset ++;
				}
				else
				{
					if(leafLevel3[(y-1)*levelWidth3+x])
					{
						level3_neighborinfo_red.push_back(3);
						level3_neighborinfo_red.push_back(x);
						level3_neighborinfo_red.push_back(y-1);
						level3_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel2[(y*2-1)*levelWidth2+x*2])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back(x*2);
							level3_neighborinfo_red.push_back(y*2-1);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4-1)*levelWidth1+x*4])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4);
								level3_neighborinfo_red.push_back(y*4-1);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+1);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
							}
							if(leafLevel1[(y*4-1)*levelWidth1+x*4+1])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4+1);
								level3_neighborinfo_red.push_back(y*4-1);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8+2])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+2);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+3])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+3);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
							}
						}
						if(leafLevel2[(y*2-1)*levelWidth2+x*2+1])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back(x*2+1);
							level3_neighborinfo_red.push_back(y*2-1);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4-1)*levelWidth1+x*4+2])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4+2);
								level3_neighborinfo_red.push_back(y*4-1);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8+4])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+4);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+5])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+5);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
							}
							if(leafLevel1[(y*4-1)*levelWidth1+x*4+3])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4+3);
								level3_neighborinfo_red.push_back(y*4-1);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8+6])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+6);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+7])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+7);
									level3_neighborinfo_red.push_back(y*8-1);
									level3_red_offset ++;
								}
							}
						}
					}
				}

				//y+
				if(y == levelHeight3-1)
				{
					level3_neighborinfo_red.push_back(3);
					level3_neighborinfo_red.push_back(x);
					level3_neighborinfo_red.push_back(y+1);
					level3_red_offset ++;
				}
				else
				{
					if(leafLevel3[(y+1)*levelWidth3+x])
					{
						level3_neighborinfo_red.push_back(3);
						level3_neighborinfo_red.push_back(x);
						level3_neighborinfo_red.push_back(y+1);
						level3_red_offset ++;
					}
					else //smaller
					{
						if(leafLevel2[(y+1)*2*levelWidth2+x*2])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back(x*2);
							level3_neighborinfo_red.push_back((y+1)*2);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y+1)*4*levelWidth1+x*4])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4);
								level3_neighborinfo_red.push_back((y+1)*4);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;
								}
								if(leafLevel0[(y+1)*8*width+x*8+1])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+1);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;
								}
							}
							if(leafLevel1[(y+1)*4*levelWidth1+x*4+1])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4+1);
								level3_neighborinfo_red.push_back((y+1)*4);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8+2])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+2);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;
								}
								if(leafLevel0[(y+1)*8*width+x*8+3])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+3);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;
								}
							}
						}
						if(leafLevel2[(y+1)*2*levelWidth2+x*2+1])
						{
							level3_neighborinfo_red.push_back(2);
							level3_neighborinfo_red.push_back(x*2+1);
							level3_neighborinfo_red.push_back((y+1)*2);
							level3_red_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y+1)*4*levelWidth1+x*4+2])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4+2);
								level3_neighborinfo_red.push_back((y+1)*4);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8+4])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+4);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;
								}
								if(leafLevel0[(y+1)*8*width+x*8+5])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+5);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;
								}
							}
							if(leafLevel1[(y+1)*4*levelWidth1+x*4+3])
							{
								level3_neighborinfo_red.push_back(1);
								level3_neighborinfo_red.push_back(x*4+3);
								level3_neighborinfo_red.push_back((y+1)*4);
								level3_red_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8+6])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+6);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;

								}
								if(leafLevel0[(y+1)*8*width+x*8+7])
								{
									level3_neighborinfo_red.push_back(0);
									level3_neighborinfo_red.push_back(x*8+7);
									level3_neighborinfo_red.push_back((y+1)*8);
									level3_red_offset ++;
								}
							}
						}
					}
				}
				level3_index_red.push_back(level3_red_offset-cur_offset);
				level3_index_red.push_back(cur_offset);
			}
			else
			{
				level3_index_black.push_back(x);
				level3_index_black.push_back(y);
				int cur_offset = level3_black_offset;

				//x-
				if(x == 0)
				{
					level3_neighborinfo_black.push_back(3);
					level3_neighborinfo_black.push_back(x-1);
					level3_neighborinfo_black.push_back(y);
					level3_black_offset ++;
				}
				else
				{
					if(leafLevel3[y*levelWidth3+x-1])
					{
						level3_neighborinfo_black.push_back(3);
						level3_neighborinfo_black.push_back(x-1);
						level3_neighborinfo_black.push_back(y);
						level3_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel2[(y*2)*levelWidth2+x*2-1])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back(x*2-1);
							level3_neighborinfo_black.push_back(y*2);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4)*levelWidth1+x*4-1])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4-1);
								level3_neighborinfo_black.push_back(y*4);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+1)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8+1);
									level3_black_offset ++;
								}
							}

							if(leafLevel1[(y*4+1)*levelWidth1+x*4-1])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4-1);
								level3_neighborinfo_black.push_back(y*4+1);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+2)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8+2);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+3)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8+3);
									level3_black_offset ++;
								}
							}
						}

						if(leafLevel2[(y*2+1)*levelWidth2+x*2-1])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back(x*2-1);
							level3_neighborinfo_black.push_back(y*2+1);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[y*4+2]*levelWidth1+x*4-1)
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4-1);
								level3_neighborinfo_black.push_back(y*4+2);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+4)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8+4);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+5)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8+5);
									level3_black_offset ++;
								}
							}

							if(leafLevel1[(y*4+3)*levelWidth1+x*4-1])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4-1);
								level3_neighborinfo_black.push_back(y*4+3);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+6)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8+6);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+7)*width+x*8-1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8-1);
									level3_neighborinfo_black.push_back(y*8+7);
									level3_black_offset ++;
								}
							}
						}
					}
				}

				//x+
				if(x == levelWidth3-1)
				{
					level3_neighborinfo_black.push_back(3);
					level3_neighborinfo_black.push_back(x+1);
					level3_neighborinfo_black.push_back(y);
					level3_black_offset ++;
				}
				else
				{
					if(leafLevel3[y*levelWidth3+x+1])
					{
						level3_neighborinfo_black.push_back(3);
						level3_neighborinfo_black.push_back(x+1);
						level3_neighborinfo_black.push_back(y);
						level3_black_offset ++;
					}
					else //smaller 
					{
						if(leafLevel2[(y*2)*levelWidth2+(x+1)*2])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back((x+1)*2);
							level3_neighborinfo_black.push_back(y*2);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back((x+1)*4);
								level3_neighborinfo_black.push_back(y*4);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+1)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8+1);
									level3_black_offset ++;
								}
							}
							if(leafLevel1[(y*4+1)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back((x+1)*4);
								level3_neighborinfo_black.push_back(y*4+1);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+2)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8+2);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+3)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8+3);
									level3_black_offset ++;
								}
							}
						}
						if(leafLevel2[(y*2+1)*levelWidth2+(x+1)*2])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back((x+1)*2);
							level3_neighborinfo_black.push_back(y*2+1);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4+2)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back((x+1)*4);
								level3_neighborinfo_black.push_back(y*4+2);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+4)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8+4);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+5)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8+5);
									level3_black_offset ++;
								}
							}
							if(leafLevel1[(y*4+3)*levelWidth1+(x+1)*4])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back((x+1)*4);
								level3_neighborinfo_black.push_back(y*4+3);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8+6)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8+6);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8+7)*width+(x+1)*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back((x+1)*8);
									level3_neighborinfo_black.push_back(y*8+7);
									level3_black_offset ++;
								}
							}
						}
					}
				}

				//y-
				if(y == 0)
				{
					level3_neighborinfo_black.push_back(3);
					level3_neighborinfo_black.push_back(x);
					level3_neighborinfo_black.push_back(y-1);
					level3_black_offset ++;
				}
				else
				{
					if(leafLevel3[(y-1)*levelWidth3+x])
					{
						level3_neighborinfo_black.push_back(3);
						level3_neighborinfo_black.push_back(x);
						level3_neighborinfo_black.push_back(y-1);
						level3_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel2[(y*2-1)*levelWidth2+x*2])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back(x*2);
							level3_neighborinfo_black.push_back(y*2-1);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4-1)*levelWidth1+x*4])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4);
								level3_neighborinfo_black.push_back(y*4-1);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+1);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
							}
							if(leafLevel1[(y*4-1)*levelWidth1+x*4+1])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4+1);
								level3_neighborinfo_black.push_back(y*4-1);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8+2])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+2);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+3])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+3);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
							}
						}
						if(leafLevel2[(y*2-1)*levelWidth2+x*2+1])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back(x*2+1);
							level3_neighborinfo_black.push_back(y*2-1);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y*4-1)*levelWidth1+x*4+2])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4+2);
								level3_neighborinfo_black.push_back(y*4-1);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8+4])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+4);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+5])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+5);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
							}
							if(leafLevel1[(y*4-1)*levelWidth1+x*4+3])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4+3);
								level3_neighborinfo_black.push_back(y*4-1);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y*8-1)*width+x*8+6])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+6);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
								if(leafLevel0[(y*8-1)*width+x*8+7])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+7);
									level3_neighborinfo_black.push_back(y*8-1);
									level3_black_offset ++;
								}
							}
						}
					}
				}

				//y+
				if(y == levelHeight3-1)
				{
					level3_neighborinfo_black.push_back(3);
					level3_neighborinfo_black.push_back(x);
					level3_neighborinfo_black.push_back(y+1);
					level3_black_offset ++;
				}
				else
				{
					if(leafLevel3[(y+1)*levelWidth3+x])
					{
						level3_neighborinfo_black.push_back(3);
						level3_neighborinfo_black.push_back(x);
						level3_neighborinfo_black.push_back(y+1);
						level3_black_offset ++;
					}
					else //smaller
					{
						if(leafLevel2[(y+1)*2*levelWidth2+x*2])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back(x*2);
							level3_neighborinfo_black.push_back((y+1)*2);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y+1)*4*levelWidth1+x*4])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4);
								level3_neighborinfo_black.push_back((y+1)*4);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;
								}
								if(leafLevel0[(y+1)*8*width+x*8+1])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+1);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;
								}
							}
							if(leafLevel1[(y+1)*4*levelWidth1+x*4+1])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4+1);
								level3_neighborinfo_black.push_back((y+1)*4);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8+2])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+2);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;
								}
								if(leafLevel0[(y+1)*8*width+x*8+3])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+3);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;
								}
							}
						}
						if(leafLevel2[(y+1)*2*levelWidth2+x*2+1])
						{
							level3_neighborinfo_black.push_back(2);
							level3_neighborinfo_black.push_back(x*2+1);
							level3_neighborinfo_black.push_back((y+1)*2);
							level3_black_offset ++;
						}
						else //smaller
						{
							if(leafLevel1[(y+1)*4*levelWidth1+x*4+2])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4+2);
								level3_neighborinfo_black.push_back((y+1)*4);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8+4])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+4);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;
								}
								if(leafLevel0[(y+1)*8*width+x*8+5])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+5);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;
								}
							}
							if(leafLevel1[(y+1)*4*levelWidth1+x*4+3])
							{
								level3_neighborinfo_black.push_back(1);
								level3_neighborinfo_black.push_back(x*4+3);
								level3_neighborinfo_black.push_back((y+1)*4);
								level3_black_offset ++;
							}
							else //smaller
							{
								if(leafLevel0[(y+1)*8*width+x*8+6])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+6);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;

								}
								if(leafLevel0[(y+1)*8*width+x*8+7])
								{
									level3_neighborinfo_black.push_back(0);
									level3_neighborinfo_black.push_back(x*8+7);
									level3_neighborinfo_black.push_back((y+1)*8);
									level3_black_offset ++;
								}
							}
						}
					}
				}
				level3_index_black.push_back(level3_black_offset-cur_offset);
				level3_index_black.push_back(cur_offset);
			}
		}
	}
}

