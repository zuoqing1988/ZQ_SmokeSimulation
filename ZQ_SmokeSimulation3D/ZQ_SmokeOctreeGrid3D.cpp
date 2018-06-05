#include "ZQ_SmokeOctreeGrid3D.h"
#include "ZQ_PCGSolver.h"
#include "ZQ_SparseMatrix.h"
#include "ZQ_PoissonSolver3D.h"
#include <math.h>
#include <time.h>
#include "windows.h"

ZQ_SmokeOctreeGrid3D::ZQ_SmokeOctreeGrid3D()
{	
	dim = 0;
	deltax = 0;
	deltaU = 0;
	deltaV = 0;
	deltaW = 0;
	Uidx2 = 0;
	Vidx2 = 0;
	Widx2 = 0;
	Uidx4 = 0;
	Vidx4 = 0;
	Widx4 = 0;
	Uidx8 = 0;
	Vidx8 = 0;
	Widx8 = 0;

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

ZQ_SmokeOctreeGrid3D::~ZQ_SmokeOctreeGrid3D()
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
	if(deltaW)
	{
		delete []deltaW;
		deltaW = 0;
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
	if(Widx2)
	{
		delete []Widx2;
		Widx2 = 0;
	}
	if(Widx4)
	{
		delete []Widx4;
		Widx4 = 0;
	}
	if(Widx8)
	{
		delete []Widx8;
		Widx8 = 0;
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

bool ZQ_SmokeOctreeGrid3D::InitFromGlobalGrid(const ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para, const int *interestingRegionDist)
{
	if(!_buildGlobalMatrix(srcGrid,para,interestingRegionDist))
		return false;
	if(!_buildSubMatrix())
		return false;
	return true;
}

bool ZQ_SmokeOctreeGrid3D::ReBuildFromGlobalGrid(const ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para, const int *interestingRegionDist)
{
	return _buildGlobalMatrix(srcGrid,para,interestingRegionDist);
}

bool ZQ_SmokeOctreeGrid3D::GlobalProjection(ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	return _solveGlobal(srcGrid,para);
}


bool ZQ_SmokeOctreeGrid3D::ApplyGuidingVelocity(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para, const ZQ_SmokeGuideGrid3D* guideGrid)
{

	//IMPORTANT FLAG
	bool blend_flag = para.guidingNaiveBlend;

	if(guideGrid == 0 || srcGrid == 0)
		return false;
	int guide_width = guideGrid->GetWidth();
	int guide_height = guideGrid->GetHeight();
	int guide_depth = guideGrid->GetDepth();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	if(width % guide_width != 0 || height % guide_height != 0  || depth % guide_depth != 0
		|| width/guide_width != height/guide_height || width/guide_width != depth/guide_depth)
		return false;

	int scale = width / guide_width;
	int levelWidth1 = width/2;
	int levelHeight1 = height/2;
	int levelDepth1 = depth/2;
	int levelWidth2 = width/4;
	int levelHeight2 = height/4;
	int levelDepth2 = depth/4;
	int levelWidth3 = width/8;
	int levelHeight3 = height/8;
	int levelDepth3 = depth/8;

	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();
	const float*& guideUptr = guideGrid->GetUptr();
	const float*& guideVptr = guideGrid->GetVptr();
	const float*& guideWptr = guideGrid->GetWptr();
	
	float* deltaU = new float[(width+1)*height*depth];
	memset(deltaU,0,sizeof(float)*(width+1)*height*depth);
	float* deltaV = new float[width*(height+1)*depth];
	memset(deltaV,0,sizeof(float)*width*(height+1)*depth);
	float* deltaW = new float[width*height*(depth+1)];
	memset(deltaW,0,sizeof(float)*width*height*(depth+1));

	int size;
	int i_shift,j_shift,k_shift;
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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		if(blend_flag)
		{
			float val;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					float coord_x = (float)i_shift/width*guide_width;
					float coord_y = (float)(j_shift+sj+0.5f)/height*guide_height-0.5f;
					float coord_z = (float)(k_shift+sk+0.5f)/depth*guide_depth-0.5f;

					ZQ::ZQ_ImageProcessing3D::TrilinearInterpolate(guideUptr,guide_width+1,guide_height,guide_depth,1,coord_x,coord_y,coord_z,&val,false);
					deltaU[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift] = 
						guidingCoeff * (val - Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift]);
				}
			}
			
		}
		else
		{
			float val;
			float total_coarse = 0;
			float total_fine = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					float coord_x = (float)i_shift/width*guide_width;
					float coord_y = (float)(j_shift+sj+0.5f)/height*guide_height-0.5f;
					float coord_z = (float)(k_shift+sk+0.5f)/depth*guide_depth-0.5f;

					ZQ::ZQ_ImageProcessing3D::TrilinearInterpolate(guideUptr,guide_width+1,guide_height,guide_depth,1,coord_x,coord_y,coord_z,&val,false);
					total_coarse += val;
					total_fine += Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
				}
			}
			
			float delta_tmp = guidingCoeff * (total_coarse - total_fine) / (size*size);
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					deltaU[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift] = delta_tmp;
				}
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		if(blend_flag)
		{
			float val;
			for(int sk = 0;sk < size;sk++)
			{
				for(int si = 0;si < size;si++)
				{
					float coord_x = (i_shift+si+0.5f)/width*guide_width-0.5f;
					float coord_y = (float)j_shift/height*guide_height;
					float coord_z = (k_shift+sk+0.5f)/depth*guide_depth-0.5f;
					ZQ::ZQ_ImageProcessing3D::TrilinearInterpolate(guideVptr,guide_width,guide_height+1,guide_depth,1,coord_x,coord_y,coord_z,&val,false);
					deltaV[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si] = 
						guidingCoeff * (val-Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si]);
				}
			}
		}
		else
		{
			float val;
			float total_coarse = 0;
			float total_fine = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int si = 0;si < size;si++)
				{
					float coord_x = (i_shift+si+0.5f)/width*guide_width-0.5f;
					float coord_y = (float)j_shift/height*guide_height;
					float coord_z = (k_shift+sk+0.5f)/depth*guide_depth-0.5f;
					ZQ::ZQ_ImageProcessing3D::TrilinearInterpolate(guideVptr,guide_width,guide_height+1,guide_depth,1,coord_x,coord_y,coord_z,&val,false);
					total_coarse += val;
					total_fine += Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
				}
			}

			float delta_tmp = guidingCoeff * (total_coarse - total_fine) / (size*size);
			for(int sk = 0;sk < size;sk++)
			{
				for(int si = 0;si < size;si++)
				{
					deltaV[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si] = delta_tmp;
				}
			}
		}
	}

	for(int i = 0;i < Winfos.size();i++)
	{
		int lvl = __min(Winfos[i].level1, Winfos[i].level2);
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
		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		if(blend_flag)
		{
			float val;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					float coord_x = (i_shift+si+0.5f)/width*guide_width-0.5f;
					float coord_y = (j_shift+sj+0.5f)/height*guide_height-0.5f;
					float coord_z = (float)k_shift/depth*guide_depth;
					ZQ::ZQ_ImageProcessing3D::TrilinearInterpolate(guideWptr,guide_width,guide_height,guide_depth+1,1,coord_x,coord_y,coord_z,&val,false);
					deltaW[k_shift*height*width+(j_shift+sj)*width+i_shift+si] = 
						guidingCoeff * (val-Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si]);
				}
			}
		}
		else
		{
			float val;
			float total_coarse = 0;
			float total_fine = 0;
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					float coord_x = (i_shift+si+0.5f)/width*guide_width-0.5f;
					float coord_y = (j_shift+sj+0.5f)/height*guide_height-0.5f;
					float coord_z = (float)k_shift/depth*guide_depth;
					ZQ::ZQ_ImageProcessing3D::TrilinearInterpolate(guideWptr,guide_width,guide_height,guide_depth+1,1,coord_x,coord_y,coord_z,&val,false);
					total_coarse += val;
					total_fine += Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
				}
			}

			float delta_tmp = guidingCoeff * (total_coarse - total_fine) / (size*size);
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					deltaW[k_shift*height*width+(j_shift+sj)*width+i_shift+si] = delta_tmp;
				}
			}
		}
	}

	bool*& occupyPtr = srcGrid->GetOccupyPtr();

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 1;i < width;i++)
			{
				if(specialRegion_mask[k*height*width+j*width+i-1] && specialRegion_mask[k*height*width+j*width+i])
					deltaU[k*height*(width+1)+j*(width+1)+i] = 0;
			}
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 1;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(specialRegion_mask[k*height*width+(j-1)*width+i] && specialRegion_mask[k*height*width+j*width+i])
					deltaV[k*(height+1)*width+j*width+i] = 0;
			}
		}
	}

	for(int k = 1;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(specialRegion_mask[(k-1)*height*width+j*width+i] && specialRegion_mask[k*height*width+j*width+i])
					deltaW[k*height*width+j*width+i] = 0;
			}
		}
	}
	
	

	if(guiding_specialRegion_useCUDA) //means using CUDA
	{
		for(int iregion = 0;iregion < specialRegion_num;iregion++)
		{
			bool* tmp_occupy = new bool[width*height*depth];
			int first_x = -1, first_y = -1, first_z = -1;
			float div_per_volume = 0.0f;
			int div_count = 0;
			int min_x = width;
			int min_y = height;
			int min_z = depth;
			int max_x = -1;
			int max_y = -1;
			int max_z = -1;
			for(int k = 0;k < depth;k++)
			{
				for(int j = 0;j < height;j++)
				{
					for(int i = 0;i < width;i++)
					{
						tmp_occupy[k*height*width+j*width+i] = (specialRegion_label[k*height*width+j*width+i] != iregion+1);
						if(!tmp_occupy[k*height*width+j*width+i])
						{
							if(first_x < 0 && first_y < 0 && first_z < 0)
							{
								first_x = i;
								first_y = j;
								first_z = k;
							}

							div_count ++;
							div_per_volume += Uptr[k*height*(width+1)+j*(width+1)+i+1] - Uptr[k*height*(width+1)+j*(width+1)+i] 
											+ Vptr[k*(height+1)*width+(j+1)*width+i] - Vptr[k*(height+1)*width+j*width+i]
											+ Wptr[(k+1)*height*width+j*width+i] - Wptr[k*height*width+j*width+i];

							min_x = __min(min_x,i);
							min_y = __min(min_y,j);
							min_z = __min(min_z,k);
							max_x = __max(max_x,i);
							max_y = __max(max_y,j);
							max_z = __max(max_z,k);
						}
					}
				}
			}
			
			div_per_volume /= div_count;

			
			int cur_maxIter = __min(para.maxIter,div_count);
			if(false)
			{
				ZQ_CUDA_PoissonSolver3D::SolveClosedPoissonRedBlackwithOccupy3D_MAC(Uptr,Vptr,Wptr,tmp_occupy,div_per_volume,width,height,depth,cur_maxIter);
			}
			else
			{
				int offset_x = min_x;
				int offset_y = min_y;
				int offset_z = min_z;
				int size_x = max_x - min_x + 1;
				int size_y = max_y - min_y + 1;
				int size_z = max_z - min_z + 1;

				bool* cur_occupy = new bool[size_x*size_y*size_z];
				float* cur_mac_u = new float[(size_x+1)*size_y*size_z];
				float* cur_mac_v = new float[size_x*(size_y+1)*size_z];
				float* cur_mac_w = new float[size_x*size_y*(size_z+1)];
				for(int k = 0;k < size_z;k++)
				{
					for(int j = 0;j < size_y;j++)
					{
						for(int i = 0;i < size_x;i++)
						{
							int rx = offset_x+i;
							int ry = offset_y+j;
							int rz = offset_z+k;
							cur_occupy[k*size_y*size_x+j*size_x+i] = tmp_occupy[rz*height*width+ry*width+rx];
						}
					}
				}
				for(int k = 0;k < size_z;k++)
				{
					for(int j = 0;j < size_y;j++)
					{
						for(int i = 0;i <= size_x;i++)
						{
							int rx = offset_x+i;
							int ry = offset_y+j;
							int rz = offset_z+k;
							cur_mac_u[k*size_y*(size_x+1)+j*(size_x+1)+i] = Uptr[rz*height*(width+1)+ry*(width+1)+rx];
						}
					}
				}

				for(int k = 0;k < size_z;k++)
				{
					for(int j = 0;j <= size_y;j++)
					{
						for(int i = 0;i < size_x;i++)
						{
							int rx = offset_x+i;
							int ry = offset_y+j;
							int rz = offset_z+k;
							cur_mac_v[k*(size_y+1)*size_x+j*size_x+i] = Vptr[rz*(height+1)*width+ry*width+rx];
						}
					}
				}

				for(int k = 0;k <= size_z;k++)
				{
					for(int j = 0;j < size_y;j++)
					{
						for(int i = 0;i < size_x;i++)
						{
							int rx = offset_x+i;
							int ry = offset_y+j;
							int rz = offset_z+k;
							cur_mac_w[k*size_y*size_x+j*size_x+i] = Wptr[rz*height*width+ry*width+rx];
						}
					}
				}

				first_x -= offset_x;
				first_y -= offset_y;
				first_z -= offset_z;
				ZQ_CUDA_PoissonSolver3D::SolveClosedPoissonRedBlackwithOccupy3D_MAC(cur_mac_u,cur_mac_v,cur_mac_w,cur_occupy,div_per_volume,size_x,size_y,size_z,cur_maxIter);

				for(int k = 0;k < size_z;k++)
				{
					for(int j = 0;j < size_y;j++)
					{
						for(int i = 0;i < size_x;i++)
						{
							int rx = offset_x+i;
							int ry = offset_y+j;
							int rz = offset_z+k;
							Uptr[rz*height*(width+1)+ry*(width+1)+rx] = cur_mac_u[k*size_y*(size_x+1)+j*(size_x+1)+i];
						}
					}
				}
				for(int k = 0;k < size_z;k++)
				{
					for(int j = 0;j <= size_y;j++)
					{
						for(int i = 0;i < size_x;i++)
						{
							int rx = offset_x+i;
							int ry = offset_y+j;
							int rz = offset_z+k;
							Vptr[rz*(height+1)*width+ry*width+rx] = cur_mac_v[k*(size_y+1)*size_x+j*size_x+i];
						}
					}
				}
				for(int k = 0;k <= size_z;k++)
				{
					for(int j = 0;j < size_y;j++)
					{
						for(int i = 0;i < size_x;i++)
						{
							int rx = offset_x+i;
							int ry = offset_y+j;
							int rz = offset_z+k;
							Wptr[rz*height*width+ry*width+rx] = cur_mac_w[k*size_y*size_x+j*size_x+i];
						}
					}
				}
				
				delete []cur_mac_u;
				delete []cur_mac_v;
				delete []cur_mac_w;
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
				for(int k = 0;k < depth;k++)
				{
					for(int j = 0;j < height;j++)
					{
						for(int i = 0;i < width;i++)
						{
							if(specialRegion_label[k*height*width+j*width+i] == iregion+1 && !occupyPtr[k*height*width+j*width+i])
							{
								divergence[idx++] = deltaU[k*height*(width+1)+j*(width+1)+i+1] - deltaU[k*height*(width+1)+j*(width+1)+i]
													+ deltaV[k*(height+1)*width+(j+1)*width+i] - deltaV[k*(height+1)*width+j*width+i]
													+ deltaW[(k+1)*height*width+j*width+i] - deltaW[k*height*width+j*width+i];
							}
						}
					}
				}
				
				double tol = 1e-9;
				int max_iter = para.maxIter;
				int it = 0;
				ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(specialRegion_A[iregion], divergence, x0, max_iter, tol, x, it, false);

				float* pressure = new float[width*height*depth];
				memset(pressure,0,sizeof(float)*width*height*depth);
				idx = 0;
				for(int k = 0;k < depth;k++)
				{
					for(int j = 0;j < height;j++)
					{
						for(int i = 0;i < width;i++)
						{
							if(specialRegion_label[k*height*width+j*width+i] == iregion+1 && !occupyPtr[k*height*width+j*width+i])
							{
								if(idx == 0)
									idx ++;
								else
								{
									pressure[k*height*width+j*width+i] = x[idx-1];
									idx++;
								}
							}
						}
					}
				}
				
				for(int k = 0;k < depth;k++)
				{
					for(int j = 0;j < height;j++)
					{
						for(int i = 1;i < width;i++)
						{
							if(specialRegion_label[k*height*width+j*width+(i-1)] == iregion+1 && specialRegion_label[k*height*width+j*width+i] == iregion+1)
								deltaU[k*height*(width+1)+j*(width+1)+i] = -(pressure[k*height*width+j*width+i] - pressure[k*height*width+j*width+i-1]);
						}
					}
				}
				
				for(int k = 0;k < depth;k++)
				{
					for(int j = 1;j < height;j++)
					{
						for(int i = 0;i < width;i++)
						{
							if(specialRegion_label[k*height*width+(j-1)*width+i] == iregion+1 && specialRegion_label[k*height*width+j*width+i] == iregion+1)
								deltaV[k*(height+1)*width+j*width+i] = -(pressure[k*height*width+j*width+i] - pressure[k*height*width+(j-1)*width+i]);
						}
					}
				}

				for(int k = 1;k < depth;k++)
				{
					for(int j = 0;j < height;j++)
					{
						for(int i = 0;i < width;i++)
						{
							if(specialRegion_label[(k-1)*height*width+j*width+i] == iregion+1 && specialRegion_label[k*height*width+j*width+i] == iregion+1)
								deltaW[k*height*width+j*width+i] = -(pressure[k*height*width+j*width+i] - pressure[(k-1)*height*width+j*width+i]);
						}
					}
				}
				

				delete []pressure;
				delete []x0;
				delete []x;
				delete []divergence;
			}
		}

		for(int i = 0;i < (width+1)*height*depth;i++)
			Uptr[i] += deltaU[i];
		for(int i = 0;i < width*(height+1)*depth;i++)
			Vptr[i] += deltaV[i];
		for(int i = 0;i < width*height*(depth+1);i++)
			Wptr[i] += deltaW[i];
		delete []deltaU;
		delete []deltaV;
		delete []deltaW;
		deltaU = 0;
		deltaV = 0;
		deltaW = 0;
	}
	

	return true;
}

bool ZQ_SmokeOctreeGrid3D::SubProjection(ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	LARGE_INTEGER tf;
	QueryPerformanceFrequency(&tf);
	LARGE_INTEGER pre_st, pre_end, aft_st, aft_end;
	QueryPerformanceCounter(&pre_st);
	
	if(srcGrid == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	
	int size,i_shift,j_shift,k_shift;
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();
	
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

	int inputRow1 = 3*2*2*(2+1);
	int inputCol1 = count1;
	int outputRow1 = 3*2*2*(2-1);
	int outputCol1 = count1;
	int inputRow2 = 3*4*4*(4+1);
	int inputCol2 = count2;
	int outputRow2 = 3*4*4*(4-1);
	int outputCol2 = count2;
	int inputRow3 = 3*8*8*(8+1);
	int inputCol3 = count3;
	int outputRow3 = 3*8*8*(8-1);
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
			k_shift = lamdainfos[i].k1 * size;
			int row_idx = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si <= size;si++)
					{
						inputVelocity1[row_idx*inputCol1+col_idx1] = Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)];
						row_idx ++;
					}
				}
			}
			
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj <= size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						inputVelocity1[row_idx*inputCol1+col_idx1] = Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)];
						row_idx ++;
					}
				}
			}

			for(int sk = 0;sk <= size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						inputVelocity1[row_idx*inputCol1+col_idx1] = Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)];
						row_idx ++;
					}
				}
			}
			
			col_idx1 ++;
		}
		else if(lamdainfos[i].level1 == 2)
		{
			size = 4;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			k_shift = lamdainfos[i].k1 * size;
			int row_idx = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si <= size;si++)
					{
						inputVelocity2[row_idx*inputCol2+col_idx2] = Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)];
						row_idx ++;
					}
				}
			}
			
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj <= size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						inputVelocity2[row_idx*inputCol2+col_idx2] = Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)];
						row_idx ++;
					}
				}
			}

			for(int sk = 0;sk <= size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						inputVelocity2[row_idx*inputCol2+col_idx2] = Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)];
						row_idx ++;
					}
				}
			}
			
			col_idx2 ++;
		}
		else if(lamdainfos[i].level1 == 3)
		{
			size = 8;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			k_shift = lamdainfos[i].k1 * size;
			int row_idx = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si <= size;si++)
					{
						inputVelocity3[row_idx*inputCol3+col_idx3] = Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)];
						row_idx ++;
					}
				}
			}
			
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj <= size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						inputVelocity3[row_idx*inputCol3+col_idx3] = Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)];
						row_idx ++;
					}
				}
			}

			for(int sk = 0;sk <= size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						inputVelocity3[row_idx*inputCol3+col_idx3] = Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)];
						row_idx ++;
					}
				}
			}
			
			col_idx3 ++;
		}
	}

	QueryPerformanceCounter(&pre_end);
	

	bool returnflag = true;
	float cuda_cost_time = 0.0f;
	if(count1 > 0)
	{
		float time1 = 0;
		SubProjection3DCuda(outputRow1, inputRow1, outputCol1, subMatrix2, inputVelocity1, outputVelocity1);
		cuda_cost_time += time1;
	}
	if(count2 > 0)
	{
		float time2 = 0;
		SubProjection3DCuda(outputRow2, inputRow2, outputCol2, subMatrix4, inputVelocity2, outputVelocity2);
		cuda_cost_time += time2;
	}

	if(count3 > 0)
	{
		float time3 = 0;
		SubProjection3DCuda(outputRow3, inputRow3, outputCol3, subMatrix8, inputVelocity3, outputVelocity3);
		cuda_cost_time += time3;
	}
	

	QueryPerformanceCounter(&aft_st);

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
			k_shift = lamdainfos[i].k1 * size;
			int row_idx = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 1;si < size;si++)
					{
						Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)] = outputVelocity1[row_idx*outputCol1+col_idx1];
						row_idx ++;
					}
				}
			}
			
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 1;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)] = outputVelocity1[row_idx*outputCol1+col_idx1];
						row_idx ++;
					}
				}
			}

			for(int sk = 1;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)] = outputVelocity1[row_idx*outputCol1+col_idx1];
					}
				}
			}
			
			col_idx1 ++;
		}
		else if(lamdainfos[i].level1 == 2)
		{
			size = 4;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			k_shift = lamdainfos[i].k1 * size;
			int row_idx = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 1;si < size;si++)
					{
						Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)] = outputVelocity2[row_idx*outputCol2+col_idx2];
						row_idx ++;
					}
				}
			}
			
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 1;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)] = outputVelocity2[row_idx*outputCol2+col_idx2];
						row_idx ++;
					}
				}
			}

			for(int sk = 1;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)] = outputVelocity2[row_idx*outputCol2+col_idx2];
						row_idx ++;
					}
				}
			}
			
			col_idx2 ++;
		}
		else if(lamdainfos[i].level1 == 3)
		{
			size = 8;
			i_shift = lamdainfos[i].i1 * size;
			j_shift = lamdainfos[i].j1 * size;
			k_shift = lamdainfos[i].k1 * size;
			int row_idx = 0;
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 1;si < size;si++)
					{
						Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)] = outputVelocity3[row_idx*outputCol3+col_idx3];
						row_idx ++;
					}
				}
			}
			
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 1;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)] = outputVelocity3[row_idx*outputCol3+col_idx3];
						row_idx ++;
					}
				}
			}

			for(int sk = 1;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					for(int si = 0;si < size;si++)
					{
						Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)] = outputVelocity3[row_idx*outputCol3+col_idx3];
						row_idx ++;
					}
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

	QueryPerformanceCounter(&aft_end);

	cuda_cost_time *= 0.001;
	double pre_time = 1.0*(pre_end.QuadPart-pre_st.QuadPart)/tf.QuadPart;
	double aft_time = 1.0*(aft_end.QuadPart-aft_st.QuadPart)/tf.QuadPart;
	double total_time = cuda_cost_time + pre_time+aft_time;
	printf("sub projection: cuda / total = %f / %f\n",cuda_cost_time, total_time);
	return returnflag;
}

bool ZQ_SmokeOctreeGrid3D::_buildGlobalMatrix(const ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para, const int* interestingRegionDist)
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

bool ZQ_SmokeOctreeGrid3D::_solveGlobal(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para)
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

bool ZQ_SmokeOctreeGrid3D::_buildGlobalMatrixOpenPoissonSOR(const ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para, const int* interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	//clock_t t1 = clock();
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	//clock_t t2 = clock();
	if(para.useGuidingControl)
		_selectInfosForCPUSolvers(para.width,para.height,para.depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,Winfos,lamdainfos);

	//clock_t t3 = clock();
	switch(cuda_openoctree_poisson_type)
	{
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		_selectInfos_for_CUDAOpenPoissonRedBlack3(para.width,para.height,para.depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos,
			level0_index_red,level0_neighborinfo_red,level0_index_black,level0_neighborinfo_black,
			level1_index_red,level1_neighborinfo_red,level1_index_black,level1_neighborinfo_black,
			level2_index_red,level2_neighborinfo_red,level2_index_black,level2_neighborinfo_black,
			level3_index_red,level3_neighborinfo_red,level3_index_black,level3_neighborinfo_black);
		break;
	}
	//clock_t t4 = clock();
	//printf("octree time part1: %.3f  part2: %.3f\n",0.001*(t2-t1),0.001*(t4-t3));
	
	return true;
}

bool ZQ_SmokeOctreeGrid3D::_buildGlobalMatrixClosedPoissonSOR(const ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para, const int* interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	if(para.useGuidingControl)
		_selectInfosForCPUSolvers(para.width,para.height,para.depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,Winfos,lamdainfos);

	switch(cuda_openoctree_poisson_type)
	{
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		_selectInfos_for_CUDAOpenPoissonRedBlack3(para.width,para.height,para.depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos,
			level0_index_red,level0_neighborinfo_red,level0_index_black,level0_neighborinfo_black,
			level1_index_red,level1_neighborinfo_red,level1_index_black,level1_neighborinfo_black,
			level2_index_red,level2_neighborinfo_red,level2_index_black,level2_neighborinfo_black,
			level3_index_red,level3_neighborinfo_red,level3_index_black,level3_neighborinfo_black);
		break;
	}

	return true;
}

bool ZQ_SmokeOctreeGrid3D::_buildGlobalMatrixOpenPoisson(const ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para, const int *interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();

	_selectInfosForCPUSolvers(width,height,depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,Winfos,lamdainfos);

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();
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
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;
	lamda_index[0] = new int[width*height*depth];
	lamda_index[1] = new int[levelWidth1*levelHeight1*levelDepth1];
	lamda_index[2] = new int[levelWidth2*levelHeight2*levelDepth2];
	lamda_index[3] = new int[levelWidth3*levelHeight3*levelDepth3];
	for(int i = 0;i < width*height*depth;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1*levelDepth1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2*levelDepth2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3*levelDepth3;i++)
	{
		lamda_index[3][i] = -1;
	}

	for(int i = 0;i < lamdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;
		int kk = lamdainfos[i].k1;

		switch(lvl)
		{
		case 0:
			lamda_index[lvl][kk*height*width+jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][kk*levelHeight1*levelWidth1+jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][kk*levelHeight2*levelWidth2+jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][kk*levelHeight3*levelWidth3+jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}
	dim = lamdaCount;

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	double area[4] = {1,4,16,64};
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].k1*height*width+Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].k1*levelHeight1*levelWidth1+Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].k1*levelHeight2*levelWidth2+Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].k1*levelHeight3*levelWidth3+Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
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
			col2 = lamda_index[0][Uinfos[i].k2*height*width+Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].k2*levelHeight1*levelWidth1+Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].k2*levelHeight2*levelWidth2+Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].k2*levelHeight3*levelWidth3+Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
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
			col1 = lamda_index[0][Vinfos[i].k1*height*width+Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].k1*levelHeight1*levelWidth1+Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].k1*levelHeight2*levelWidth2+Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].k1*levelHeight3*levelWidth3+Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
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
			col2 = lamda_index[0][Vinfos[i].k2*height*width+Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].k2*levelHeight1*levelWidth1+Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].k2*levelHeight2*levelWidth2+Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].k2*levelHeight3*levelWidth3+Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
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

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1,Winfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Winfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Winfos[i].k1*height*width+Winfos[i].j1*width+Winfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Winfos[i].k1*levelHeight1*levelWidth1+Winfos[i].j1*levelWidth1+Winfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Winfos[i].k1*levelHeight2*levelWidth2+Winfos[i].j1*levelWidth2+Winfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Winfos[i].k1 >= 0)
			{
				col1 = lamda_index[3][Winfos[i].k1*levelHeight3*levelWidth3+Winfos[i].j1*levelWidth3+Winfos[i].i1];
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

		switch(Winfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Winfos[i].k2*height*width+Winfos[i].j2*width+Winfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Winfos[i].k2*levelHeight1*levelWidth1+Winfos[i].j2*levelWidth1+Winfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Winfos[i].k2*levelHeight2*levelWidth2+Winfos[i].j2*levelWidth2+Winfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Winfos[i].k2 < levelDepth3)
			{
				col2 = lamda_index[3][Winfos[i].k2*levelHeight3*levelWidth3+Winfos[i].j2*levelWidth3+Winfos[i].i2];
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

bool ZQ_SmokeOctreeGrid3D::_buildGlobalMatrixClosedPoisson(const ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para, const int *interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();

	_selectInfosForCPUSolvers(width,height,depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,Winfos,lamdainfos);

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();
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
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;

	lamda_index[0] = new int[width*height*depth];
	lamda_index[1] = new int[levelWidth1*levelHeight1*levelDepth1];
	lamda_index[2] = new int[levelWidth2*levelHeight2*levelDepth2];
	lamda_index[3] = new int[levelWidth3*levelHeight3*levelDepth3];
	for(int i = 0;i < width*height*depth;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1*levelDepth1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2*levelDepth2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3*levelDepth3;i++)
	{
		lamda_index[3][i] = -1;
	}

	//the last pressure is set 0 
	for(int i = 0;i < lamdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;
		int kk = lamdainfos[i].k1;

		switch(lvl)
		{
		case 0:
			lamda_index[lvl][kk*height*width+jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][kk*levelHeight1*levelWidth1+jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][kk*levelHeight2*levelWidth2+jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][kk*levelHeight3*levelWidth3+jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}
	dim = lamdaCount-1;

	ZQ::ZQ_SparseMatrix<float> mat(lamdaCount,dim);

	double area[4] = {1,4,16,64};
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].k1*height*width+Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].k1*levelHeight1*levelWidth1+Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].k1*levelHeight2*levelWidth2+Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].k1*levelHeight3*levelWidth3+Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
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
			col2 = lamda_index[0][Uinfos[i].k2*height*width+Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].k2*levelHeight1*levelWidth1+Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].k2*levelHeight2*levelWidth2+Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].k2*levelHeight3*levelWidth3+Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
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
			double val = area[lvl]/(len1+len2);
			if(col2 != lamdaCount-1)
			{
				mat.AddTo(col1,col2,val);
				mat.AddTo(col2,col2,-val);
			}
			if(col1 != lamdaCount-1)
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
			col1 = lamda_index[0][Vinfos[i].k1*height*width+Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].k1*levelHeight1*levelWidth1+Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].k1*levelHeight2*levelWidth2+Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].k1*levelHeight3*levelWidth3+Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
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
			col2 = lamda_index[0][Vinfos[i].k2*height*width+Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].k2*levelHeight1*levelWidth1+Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].k2*levelHeight2*levelWidth2+Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].k2*levelHeight3*levelWidth3+Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
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
			double val = area[lvl]/(len1+len2);
			if(col2 != lamdaCount-1)
			{
				mat.AddTo(col2,col1,val);
				mat.AddTo(col2,col2,-val);
			}
			if(col1 != lamdaCount-1)
			{
				mat.AddTo(col1,col2,val);
				mat.AddTo(col1,col1,-val);
			}
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1,Winfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Winfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Winfos[i].k1*height*width+Winfos[i].j1*width+Winfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Winfos[i].k1*levelHeight1*levelWidth1+Winfos[i].j1*levelWidth1+Winfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Winfos[i].k1*levelHeight2*levelWidth2+Winfos[i].j1*levelWidth2+Winfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Winfos[i].k1 >= 0)
			{
				col1 = lamda_index[3][Winfos[i].k1*levelHeight3*levelWidth3+Winfos[i].j1*levelWidth3+Winfos[i].i1];
				len1 = 4;
			}
			else
				col1 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Winfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Winfos[i].k2*height*width+Winfos[i].j2*width+Winfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Winfos[i].k2*levelHeight1*levelWidth1+Winfos[i].j2*levelWidth1+Winfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Winfos[i].k2*levelHeight2*levelWidth2+Winfos[i].j2*levelWidth2+Winfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Winfos[i].k2 < levelDepth3)
			{
				col2 = lamda_index[3][Winfos[i].k2*levelHeight3*levelWidth3+Winfos[i].j2*levelWidth3+Winfos[i].i2];
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
			double val = area[lvl]/(len1+len2);
			if(col2 != lamdaCount-1)
			{
				mat.AddTo(col2,col1,val);
				mat.AddTo(col2,col2,-val);
			}
			if(col1 != lamdaCount-1)
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

bool ZQ_SmokeOctreeGrid3D::_buildGlobalMatrixOpenFlux(const ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para, const int *interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();

	_selectInfosForCPUSolvers(width,height,depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,Winfos,lamdainfos);
	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();
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

	for(int i = 0;i < wCount;i++)
	{
		Winfos[i].index = cur_idx++;
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
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;


	for(int i = 0;i < 4;i++)
	{
		if(lamda_index[i])
		{
			delete []lamda_index[i];
			lamda_index[i] = 0;
		}
	}
	lamda_index[0] = new int[width*height*depth];
	lamda_index[1] = new int[levelWidth1*levelHeight1*levelDepth1];
	lamda_index[2] = new int[levelWidth2*levelHeight2*levelDepth2];
	lamda_index[3] = new int[levelWidth3*levelHeight3*levelDepth3];
	for(int i = 0;i < width*height*depth;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1*levelDepth1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2*levelDepth2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3*levelDepth3;i++)
	{
		lamda_index[3][i] = -1;
	}

	for(int i = 0;i < lamdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;
		int kk = lamdainfos[i].k1;
		switch(lvl)
		{
		case 0:
			lamda_index[lvl][kk*height*width+jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][kk*levelHeight1*levelWidth1+jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][kk*levelHeight2*levelWidth2+jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][kk*levelHeight3*levelWidth3+jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	double area[4] = {1,4,16,64};
	double weight[4] = {1,4,16,64};
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
			col = lamda_index[0][Uinfos[i].k1*height*width+Uinfos[i].j1*width+Uinfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].k1*levelHeight1*levelWidth1+Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].k1*levelHeight2*levelWidth2+Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col = lamda_index[3][Uinfos[i].k1*levelHeight3*levelWidth3+Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
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
			col = lamda_index[0][Uinfos[i].k2*height*width+Uinfos[i].j2*width+Uinfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].k2*levelHeight1*levelWidth1+Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].k2*levelHeight2*levelWidth2+Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col = lamda_index[3][Uinfos[i].k2*levelHeight3*levelWidth3+Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
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
			col = lamda_index[0][Vinfos[i].k1*height*width+Vinfos[i].j1*width+Vinfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].k1*levelHeight1*levelWidth1+Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].k1*levelHeight2*levelWidth2+Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col = lamda_index[3][Vinfos[i].k1*levelHeight3*levelWidth3+Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
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
			col = lamda_index[0][Vinfos[i].k2*height*width+Vinfos[i].j2*width+Vinfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].k2*levelHeight1*levelWidth1+Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].k2*levelHeight2*levelWidth2+Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col = lamda_index[3][Vinfos[i].k2*levelHeight3*levelWidth3+Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
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

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1,Winfos[i].level2);
		int row = Winfos[i].index;
		int col = Winfos[i].index;
		float wt1 = 1, wt2 = 1;

		double val;
		switch(Winfos[i].level1)
		{
		case 0:
			col = lamda_index[0][Winfos[i].k1*height*width+Winfos[i].j1*width+Winfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Winfos[i].k1*levelHeight1*levelWidth1+Winfos[i].j1*levelWidth1+Winfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Winfos[i].k1*levelHeight2*levelWidth2+Winfos[i].j1*levelWidth2+Winfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Winfos[i].k1 >= 0)
			{
				col = lamda_index[3][Winfos[i].k1*levelHeight3*levelWidth3+Winfos[i].j1*levelWidth3+Winfos[i].i1];
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

		switch(Winfos[i].level2)
		{
		case 0:
			col = lamda_index[0][Winfos[i].k2*height*width+Winfos[i].j2*width+Winfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Winfos[i].k2*levelHeight1*levelWidth1+Winfos[i].j2*levelWidth1+Winfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Winfos[i].k2*levelHeight2*levelWidth2+Winfos[i].j2*levelWidth2+Winfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Winfos[i].k2 < levelDepth3)
			{
				col = lamda_index[3][Winfos[i].k2*levelHeight3*levelWidth3+Winfos[i].j2*levelWidth3+Winfos[i].i2];
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
		col = Winfos[i].index;
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

bool ZQ_SmokeOctreeGrid3D::_buildGlobalMatrixClosedFlux(const ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para, const int* interestingRegionDist)
{
	level_dis_thresh1 = para.octree_thresh1;
	level_dis_thresh2 = para.octree_thresh2;
	level_dis_thresh3 = para.octree_thresh3;
	if(!_buildOctreeInfo(srcGrid,para,interestingRegionDist))
		return false;

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();

	_selectInfosForCPUSolvers(width,height,depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,Uinfos,Vinfos,Winfos,lamdainfos);

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();
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

	for(int i = 0;i < wCount;i++)
	{
		if((Winfos[i].level1 == 3 && Winfos[i].k1 < 0)  || (Winfos[i].level2 == 3 && Winfos[i].k2 >= depth/8))
			Winfos[i].index = -1;
		else
			Winfos[i].index = cur_idx++;
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
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;

	lamda_index[0] = new int[width*height*depth];
	lamda_index[1] = new int[levelWidth1*levelHeight1*levelDepth1];
	lamda_index[2] = new int[levelWidth2*levelHeight2*levelDepth2];
	lamda_index[3] = new int[levelWidth3*levelHeight3*levelDepth3];
	for(int i = 0;i < width*height*depth;i++)
	{
		lamda_index[0][i] = -1;
	}
	for(int i = 0;i < levelWidth1*levelHeight1*levelDepth1;i++)
	{
		lamda_index[1][i] = -1;
	}
	for(int i = 0;i < levelWidth2*levelHeight2*levelDepth2;i++)
	{
		lamda_index[2][i] = -1;
	}
	for(int i = 0;i < levelWidth3*levelHeight3*levelDepth3;i++)
	{
		lamda_index[3][i] = -1;
	}

	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int ii = lamdainfos[i].i1;
		int jj = lamdainfos[i].j1;
		int kk = lamdainfos[i].k1;
		switch(lvl)
		{
		case 0:
			lamda_index[lvl][kk*height*width+jj*width+ii] = lamdainfos[i].index;
			break;
		case 1:
			lamda_index[lvl][kk*levelHeight1*levelWidth1+jj*levelWidth1+ii] = lamdainfos[i].index;
			break;
		case 2:
			lamda_index[lvl][kk*levelHeight2*levelWidth2+jj*levelWidth2+ii] = lamdainfos[i].index;
			break;
		case 3:
			lamda_index[lvl][kk*levelHeight3*levelWidth3+jj*levelWidth3+ii] = lamdainfos[i].index;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	double area[4] = {1,4,16,64};
	double weight[4] = {1,4,16,64};
	double volume[4] = {1,1,1,1};
	for(int i = 0;i < lambdaCount;i++)
	{
		int lvl = lamdainfos[i].level1;
		int row = lamdainfos[i].index;
		int col = k_index;
		mat.SetValue(row,col,1);
		mat.SetValue(col,row,1);
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
			col = lamda_index[0][Uinfos[i].k1*height*width+Uinfos[i].j1*width+Uinfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].k1*levelHeight1*levelWidth1+Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].k1*levelHeight2*levelWidth2+Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col = lamda_index[3][Uinfos[i].k1*levelHeight3*levelWidth3+Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
				val = -area[lvl]/volume[3];
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
			col = lamda_index[0][Uinfos[i].k2*height*width+Uinfos[i].j2*width+Uinfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Uinfos[i].k2*levelHeight1*levelWidth1+Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Uinfos[i].k2*levelHeight2*levelWidth2+Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col = lamda_index[3][Uinfos[i].k2*levelHeight3*levelWidth3+Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
				val = area[lvl]/volume[3];
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
			col = lamda_index[0][Vinfos[i].k1*height*width+Vinfos[i].j1*width+Vinfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].k1*levelHeight1*levelWidth1+Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].k1*levelHeight2*levelWidth2+Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col = lamda_index[3][Vinfos[i].k1*levelHeight3*levelWidth3+Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
				val = -area[lvl]/volume[3];
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
			col = lamda_index[0][Vinfos[i].k2*height*width+Vinfos[i].j2*width+Vinfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Vinfos[i].k2*levelHeight1*levelWidth1+Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Vinfos[i].k2*levelHeight2*levelWidth2+Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col = lamda_index[3][Vinfos[i].k2*levelHeight3*levelWidth3+Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
				val = area[lvl]/volume[3];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1,Winfos[i].level2);
		int row = Winfos[i].index;
		int col = Winfos[i].index;
		if(row < 0)
			continue;
		mat.SetValue(row,col,2*weight[lvl]);
		double val;
		switch(Winfos[i].level1)
		{
		case 0:
			col = lamda_index[0][Winfos[i].k1*height*width+Winfos[i].j1*width+Winfos[i].i1];
			val = -area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Winfos[i].k1*levelHeight1*levelWidth1+Winfos[i].j1*levelWidth1+Winfos[i].i1];
			val = -area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Winfos[i].k1*levelHeight2*levelWidth2+Winfos[i].j1*levelWidth2+Winfos[i].i1];
			val = -area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Winfos[i].k1 >= 0)
			{
				col = lamda_index[3][Winfos[i].k1*levelHeight3*levelWidth3+Winfos[i].j1*levelWidth3+Winfos[i].i1];
				val = -area[lvl]/volume[3];
				mat.SetValue(row,col,val);
				mat.SetValue(col,row,val);
			}
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Winfos[i].level2)
		{
		case 0:
			col = lamda_index[0][Winfos[i].k2*height*width+Winfos[i].j2*width+Winfos[i].i2];
			val = area[lvl]/volume[0];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 1:
			col = lamda_index[1][Winfos[i].k2*levelHeight1*levelWidth1+Winfos[i].j2*levelWidth1+Winfos[i].i2];
			val = area[lvl]/volume[1];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 2:
			col = lamda_index[2][Winfos[i].k2*levelHeight2*levelWidth2+Winfos[i].j2*levelWidth2+Winfos[i].i2];
			val = area[lvl]/volume[2];
			mat.SetValue(row,col,val);
			mat.SetValue(col,row,val);
			break;
		case 3:
			if(Winfos[i].k2 < levelDepth3)
			{
				col = lamda_index[3][Winfos[i].k2*levelHeight3*levelWidth3+Winfos[i].j2*levelWidth3+Winfos[i].i2];
				val = area[lvl]/volume[3];
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




bool ZQ_SmokeOctreeGrid3D::_solveGlobalOpenPoissonSOR(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(srcGrid == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();

	switch(cuda_openoctree_poisson_type)
	{
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		{
			const int index_channels = 5;
			const int neighborinfo_channels = 4;
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
			float cuda_cost_time = 
			ZQ_CUDA_PoissonSolver3D::SolveOpenOctreePoissonRedBlack3_3D_MAC(
			//CPU_SolveOpenOctreePoissonRedBlack3_MAC(
				Uptr,Vptr,Wptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,width,height,depth,para.maxIter,
				level0_num_red,		level0_index_red_ptr,	level0_info_len_red,	level0_info_red_ptr,
				level0_num_black,	level0_index_black_ptr,	level0_info_len_black,	level0_info_black_ptr,
				level1_num_red,		level1_index_red_ptr,	level1_info_len_red,	level1_info_red_ptr,
				level1_num_black,	level1_index_black_ptr,	level1_info_len_black,	level1_info_black_ptr,
				level2_num_red,		level2_index_red_ptr,	level2_info_len_red,	level2_info_red_ptr,
				level2_num_black,	level2_index_black_ptr,	level2_info_len_black,	level2_info_black_ptr,
				level3_num_red,		level3_index_red_ptr,	level3_info_len_red,	level3_info_red_ptr,
				level3_num_black,	level3_index_black_ptr,	level3_info_len_black,	level3_info_black_ptr);

			printf("open octree global projection cuda cost = %f\n",0.001*cuda_cost_time);
		}
		break;
	}

	return true;
}


bool ZQ_SmokeOctreeGrid3D::_solveGlobalClosedPoissonSOR(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(srcGrid == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();
	bool*& occupy = srcGrid->GetOccupyPtr();

	float div_per_volume = 0.0f;
	int count = 0;

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(!occupy[k*height*width+j*width+i])
				{
					count += 1;
					div_per_volume += Uptr[k*height*(width+1)+j*(width+1)+i+1] - Uptr[k*height*(width+1)+j*(width+1)+i] 
									+ Vptr[k*(height+1)*width+(j+1)*width+i] - Vptr[k*(height+1)*width+j*width+i]
									+ Wptr[(k+1)*height*width+j*width+i] - Wptr[k*height*width+j*width+i];
				}
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
	case CUDA_OPEN_OCTREE_POISSON_SOR_TYPE3:
		{
			const int index_channels = 5;
			const int neighborinfo_channels = 4;
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
			float cuda_cost_time = 
			ZQ_CUDA_PoissonSolver3D::SolveClosedOctreePoissonRedBlack3_3D_MAC(Uptr,Vptr,Wptr,leafLevel0,leafLevel1,leafLevel2,leafLevel3,div_per_volume,width,height,depth,para.maxIter,
				level0_num_red,		level0_index_red_ptr,	level0_info_len_red,	level0_info_red_ptr,
				level0_num_black,	level0_index_black_ptr,	level0_info_len_black,	level0_info_black_ptr,
				level1_num_red,		level1_index_red_ptr,	level1_info_len_red,	level1_info_red_ptr,
				level1_num_black,	level1_index_black_ptr,	level1_info_len_black,	level1_info_black_ptr,
				level2_num_red,		level2_index_red_ptr,	level2_info_len_red,	level2_info_red_ptr,
				level2_num_black,	level2_index_black_ptr,	level2_info_len_black,	level2_info_black_ptr,
				level3_num_red,		level3_index_red_ptr,	level3_info_len_red,	level3_info_red_ptr,
				level3_num_black,	level3_index_black_ptr,	level3_info_len_black,	level3_info_black_ptr);

			printf("closed octree global projection cuda cost = %f\n",0.001*cuda_cost_time);
		}
		break;
	}

	return true;
}

bool ZQ_SmokeOctreeGrid3D::_solveGlobalOpenPoisson(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(srcGrid == 0 || A == 0)
		return false;
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();
	
	int dim = A->n;
	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[A->m];
	memset(b,0,sizeof(double)*A->m);
	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	int i_shift,j_shift,k_shift;
	int Length[4] = {1,2,4,8};
	int size;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();
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
		k_shift = lamdainfos[i].k1*size;
		double val = 0;
		int size = Length[lvl];
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				val -= Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
				val += Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift+size];
			}
		}
		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				val -= Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
				val += Vptr[(k_shift+sk)*(height+1)*width+(j_shift+size)*width+i_shift+si];
			}
		}

		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				val -= Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
				val += Wptr[(k_shift+size)*height*width+(j_shift+sj)*width+i_shift+si];
			}
		}
		
		b[lamdainfos[i].index] = val;
	}

	double tol = 1e-9;
	int it = 0;
	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(A, b, x0, para.maxIter, tol, x, it, false);

	double* pressure = new double[dim];
	memcpy(pressure,x,sizeof(double)*dim);

	delete []x;
	delete []x0;
	delete []b;

	if(deltax)
		delete []deltax;
	deltax = new double[uCount+vCount+wCount];
	memset(deltax,0,sizeof(double)*(uCount+vCount+wCount));
		
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	int levelDepth3 = depth/8;
	
	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].k1*height*width+Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].k1*levelHeight1*levelWidth1+Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].k1*levelHeight2*levelWidth2+Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].k1*levelHeight3*levelWidth3+Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
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
			col2 = lamda_index[0][Uinfos[i].k2*height*width+Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].k2*levelHeight1*levelWidth1+Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].k2*levelHeight2*levelWidth2+Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].k2*levelHeight3*levelWidth3+Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
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
			col1 = lamda_index[0][Vinfos[i].k1*height*width+Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].k1*levelHeight1*levelWidth1+Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].k1*levelHeight2*levelWidth2+Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].k1*levelHeight3*levelWidth3+Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
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
			col2 = lamda_index[0][Vinfos[i].k2*height*width+Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].k2*levelHeight1*levelWidth1+Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].k2*levelHeight2*levelWidth2+Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].k2*levelHeight3*levelWidth3+Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
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

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1,Winfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Winfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Winfos[i].k1*height*width+Winfos[i].j1*width+Winfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Winfos[i].k1*levelHeight1*levelWidth1+Winfos[i].j1*levelWidth1+Winfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Winfos[i].k1*levelHeight2*levelWidth2+Winfos[i].j1*levelWidth2+Winfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Winfos[i].k1 >= 0)
			{
				col1 = lamda_index[3][Winfos[i].k1*levelHeight3*levelWidth3+Winfos[i].j1*levelWidth3+Winfos[i].i1];
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

		switch(Winfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Winfos[i].k2*height*width+Winfos[i].j2*width+Winfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Winfos[i].k2*levelHeight1*levelWidth1+Winfos[i].j2*levelWidth1+Winfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Winfos[i].k2*levelHeight2*levelWidth2+Winfos[i].j2*levelWidth2+Winfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Winfos[i].k2 < levelDepth3)
			{
				col2 = lamda_index[3][Winfos[i].k2*levelHeight3*levelWidth3+Winfos[i].j2*levelWidth3+Winfos[i].i2];
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

		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		if(k_shift != 0 && k_shift != depth)
		{
			double val = 1.0/(len1+len2);
			deltax[i+uCount+vCount] = -val*(pressure[col2]-pressure[col1]);
		}
		else if(k_shift == 0)
		{
			double val = 1.0/(2*len2);
			deltax[i+uCount+vCount] = -val*(pressure[col2]-0);
		}
		else if(k_shift == depth)
		{
			double val = 1.0/(2*len1);
			deltax[i+uCount+vCount] = -val*(0-pressure[col1]);
		}
	}

	delete []pressure;

	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;
	if(deltaW)
		delete []deltaW;
	
	deltaU = new float[(width+1)*height*depth];
	memset(deltaU,0,sizeof(float)*(width+1)*height*depth);
	deltaV = new float[width*(height+1)*depth];
	memset(deltaV,0,sizeof(float)*width*(height+1)*depth);
	deltaW = new float[width*height*(depth+1)];
	memset(deltaW,0,sizeof(float)*width*height*(depth+1));
	
	
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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				deltaU[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift] += deltax[i];
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				deltaV[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si] += deltax[i+uCount];
			}
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1, Winfos[i].level2);
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
		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				deltaW[k_shift*height*width+(j_shift+sj)*width+i_shift+si] += deltax[i+uCount+vCount];
			}
		}
	}

	delete []deltax;
	deltax = 0;

	for(int i = 0;i < (width+1)*height*depth;i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		Vptr[i] += deltaV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		Wptr[i] += deltaW[i];

	delete []deltaU;
	delete []deltaV;
	delete []deltaW;
	
	deltaU = 0;
	deltaV = 0;
	deltaW = 0;
	
	return true;
}


bool ZQ_SmokeOctreeGrid3D::_solveGlobalClosedPoisson(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(srcGrid == 0 || A == 0)
		return false;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();

	int lambdaCount = lamdainfos.size();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();

	int dim = A->n;
	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[A->m];
	memset(b,0,sizeof(double)*A->m);
	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	int i_shift,j_shift,k_shift;
	int Length[4] = {1,2,4,8};
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
		k_shift = lamdainfos[i].k1*size;
		double val = 0;
		int size = Length[lvl];
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				val -= Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
				val += Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+size)];
			}
		}
		
		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				val -= Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
				val += Vptr[(k_shift+sk)*(height+1)*width+(j_shift+size)*width+i_shift+si];
			}
		}

		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				val -= Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
				val += Wptr[(k_shift+size)*height*width+(j_shift+sj)*width+i_shift+si];
			}
		}
		
		b[lamdainfos[i].index] = val;
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
	deltax = new double[uCount+vCount+wCount];
	memset(deltax,0,sizeof(double)*(uCount+vCount+wCount));

	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	int levelDepth3 = depth/8;

	for(int i = 0;i < uCount;i++)
	{
		int lvl = __min(Uinfos[i].level1,Uinfos[i].level2);
		double len1 = 0;
		double len2 = 0;
		int col1,col2;
		switch(Uinfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Uinfos[i].k1*height*width+Uinfos[i].j1*width+Uinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Uinfos[i].k1*levelHeight1*levelWidth1+Uinfos[i].j1*levelWidth1+Uinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Uinfos[i].k1*levelHeight2*levelWidth2+Uinfos[i].j1*levelWidth2+Uinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Uinfos[i].i1 >= 0)
			{
				col1 = lamda_index[3][Uinfos[i].k1*levelHeight3*levelWidth3+Uinfos[i].j1*levelWidth3+Uinfos[i].i1];
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
			col2 = lamda_index[0][Uinfos[i].k2*height*width+Uinfos[i].j2*width+Uinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Uinfos[i].k2*levelHeight1*levelWidth1+Uinfos[i].j2*levelWidth1+Uinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Uinfos[i].k2*levelHeight2*levelWidth2+Uinfos[i].j2*levelWidth2+Uinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Uinfos[i].i2 < levelWidth3)
			{
				col2 = lamda_index[3][Uinfos[i].k2*levelHeight3*levelWidth3+Uinfos[i].j2*levelWidth3+Uinfos[i].i2];
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
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
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
			col1 = lamda_index[0][Vinfos[i].k1*height*width+Vinfos[i].j1*width+Vinfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Vinfos[i].k1*levelHeight1*levelWidth1+Vinfos[i].j1*levelWidth1+Vinfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Vinfos[i].k1*levelHeight2*levelWidth2+Vinfos[i].j1*levelWidth2+Vinfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Vinfos[i].j1 >= 0)
			{
				col1 = lamda_index[3][Vinfos[i].k1*levelHeight3*levelWidth3+Vinfos[i].j1*levelWidth3+Vinfos[i].i1];
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
			col2 = lamda_index[0][Vinfos[i].k2*height*width+Vinfos[i].j2*width+Vinfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Vinfos[i].k2*levelHeight1*levelWidth1+Vinfos[i].j2*levelWidth1+Vinfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Vinfos[i].k2*levelHeight2*levelWidth2+Vinfos[i].j2*levelWidth2+Vinfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Vinfos[i].j2 < levelHeight3)
			{
				col2 = lamda_index[3][Vinfos[i].k2*levelHeight3*levelWidth3+Vinfos[i].j2*levelWidth3+Vinfos[i].i2];
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		if(j_shift != 0 && j_shift != height)
		{
			deltax[i+uCount] = -val*(pressure[col2]-pressure[col1]);
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1,Winfos[i].level2);
		double len1 = 0,len2 = 0;
		int col1,col2;
		switch(Winfos[i].level1)
		{
		case 0:
			col1 = lamda_index[0][Winfos[i].k1*height*width+Winfos[i].j1*width+Winfos[i].i1];
			len1 = 0.5;
			break;
		case 1:
			col1 = lamda_index[1][Winfos[i].k1*levelHeight1*levelWidth1+Winfos[i].j1*levelWidth1+Winfos[i].i1];
			len1 = 1;
			break;
		case 2:
			col1 = lamda_index[2][Winfos[i].k1*levelHeight2*levelWidth2+Winfos[i].j1*levelWidth2+Winfos[i].i1];
			len1 = 2;
			break;
		case 3:
			if(Winfos[i].k1 >= 0)
			{
				col1 = lamda_index[3][Winfos[i].k1*levelHeight3*levelWidth3+Winfos[i].j1*levelWidth3+Winfos[i].i1];
				len1 = 4;
			}
			else
				col1 = -1;
			break;
		default:
			printf("something going wrong!(%s):(%d)\n",__FILE__,__LINE__);
			break;
		}

		switch(Winfos[i].level2)
		{
		case 0:
			col2 = lamda_index[0][Winfos[i].k2*height*width+Winfos[i].j2*width+Winfos[i].i2];
			len2 = 0.5;
			break;
		case 1:
			col2 = lamda_index[1][Winfos[i].k2*levelHeight1*levelWidth1+Winfos[i].j2*levelWidth1+Winfos[i].i2];
			len2 = 1;
			break;
		case 2:
			col2 = lamda_index[2][Winfos[i].k2*levelHeight2*levelWidth2+Winfos[i].j2*levelWidth2+Winfos[i].i2];
			len2 = 2;
			break;
		case 3:
			if(Winfos[i].k2 < levelDepth3)
			{
				col2 = lamda_index[3][Winfos[i].k2*levelHeight3*levelWidth3+Winfos[i].j2*levelWidth3+Winfos[i].i2];
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
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		if(k_shift != 0 && k_shift != depth)
		{
			deltax[i+uCount+vCount] = -val*(pressure[col2]-pressure[col1]);
		}
	}

	delete []pressure;

	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;

	if(deltaW)
		delete []deltaW;

	deltaU = new float[(width+1)*height*depth];
	memset(deltaU,0,sizeof(float)*(width+1)*height*depth);
	deltaV = new float[width*(height+1)*depth];
	memset(deltaV,0,sizeof(float)*width*(height+1)*depth);
	deltaW = new float[width*height*(depth+1)];
	memset(deltaW,0,sizeof(float)*width*height*(depth+1));


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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}


		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				deltaU[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift] += deltax[i];
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				deltaV[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si] += deltax[i+uCount];
			}
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1, Winfos[i].level2);
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
		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				deltaW[k_shift*height*width+(j_shift+sj)*width+i_shift+si] += deltax[i+uCount+vCount];
			}
		}
	}

	delete []deltax;
	deltax = 0;

	for(int i = 0;i < (width+1)*height*depth;i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		Vptr[i] += deltaV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		Wptr[i] += deltaW[i];

	delete []deltaU;
	delete []deltaV;
	delete []deltaW;

	deltaU = 0;
	deltaV = 0;
	deltaW = 0;
	return true;
}

bool ZQ_SmokeOctreeGrid3D::_solveGlobalOpenFlux(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(srcGrid == 0)
		return false;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();
	int lambdaCount = lamdainfos.size();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	int levelDepth3 = depth/8;

	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[dim];
	double* oldVelocity = new double[uCount+vCount+wCount];

	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	memset(b,0,sizeof(double)*dim);

	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();

	int size;
	int i_shift,j_shift,k_shift;
	double area[4] = {1,4,16,64};
	double weight[4] = {1,4,16,64};
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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		double val = 0;
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				val += Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
			}
		}
		
		b[Uinfos[i].index] = 2*(val/(size*size))*weight[lvl];
		oldVelocity[Uinfos[i].index] = val / (size*size);
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		double val = 0;
		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				val += Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
			}
		}
		
		b[Vinfos[i].index] = 2*(val/(size*size))*weight[lvl];
		oldVelocity[Vinfos[i].index] = val / (size*size);
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1, Winfos[i].level2);
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
		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		double val = 0;
		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				val += Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
			}
		}

		b[Winfos[i].index] = 2*(val/(size*size))*weight[lvl];
		oldVelocity[Winfos[i].index] = val / (size*size);
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
		
		k_shift = lamdainfos[i].k1*size;
		i_shift = lamdainfos[i].i1*size;
		j_shift = lamdainfos[i].j1*size;
		double val = 0;
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				if(occupy[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift] || (i_shift > 0 && occupy[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift-1]))
					val -= Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
				if(occupy[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift+size-1] || (i_shift+size < width && occupy[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift+size]))
					val += Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift+size];
			}
		}
		
		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				if(occupy[(k_shift+sk)*height*width+j_shift*width+i_shift+si] || (j_shift > 0 && occupy[(k_shift+sk)*height*width+(j_shift-1)*width+i_shift+si]))
					val -= Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
				if(occupy[(k_shift+sk)*height*width+(j_shift+size-1)*width+i_shift+si] || (j_shift+size < height && occupy[(k_shift+sk)*height*width+(j_shift+size)*width+i_shift+si]))
					val += Vptr[(k_shift+sk)*(height+1)*width+(j_shift+size)*width+i_shift+si];
			}
		}

		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				if(occupy[k_shift*height*width+(j_shift+sj)*width+i_shift+si] || (k_shift > 0 && occupy[(k_shift-1)*height*width+(j_shift+sj)*width+i_shift+si]))
					val -= Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
				if(occupy[(k_shift+size-1)*height*width+(j_shift+sj)*width+i_shift+si] || (k_shift+size < width && occupy[(k_shift+size)*height*width+(j_shift+sj)*width+i_shift+si]))
					val += Wptr[(k_shift+size)*height*width+(j_shift+sj)*width+i_shift+si];
			}
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
	deltax = new double[uCount+vCount+wCount];
	for(int i = 0;i < uCount+vCount+wCount;i++)
		deltax[i] = x[i] - oldVelocity[i];

	delete []x0;
	delete []x;
	delete []b;
	delete []oldVelocity;

	
	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;
	if(deltaW)
		delete []deltaW;
	
	deltaU = new float[(width+1)*height*depth];
	memset(deltaU,0,sizeof(float)*(width+1)*height*depth);
	deltaV = new float[width*(height+1)*depth];
	memset(deltaV,0,sizeof(float)*width*(height+1)*depth);
	deltaW = new float[width*height*(depth+1)];
	memset(deltaW,0,sizeof(float)*width*height*(depth+1));
	
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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}

	
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				deltaU[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift] += deltax[Uinfos[i].index];
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}
	
		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				deltaV[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si] += deltax[Vinfos[i].index];
			}
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1, Winfos[i].level2);
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
		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				deltaW[k_shift*height*width+(j_shift+sj)*width+i_shift+si] += deltax[Winfos[i].index];
			}
		}
	}

	delete []deltax;
	deltax = 0;

	

	for(int i = 0;i < (width+1)*height*depth;i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		Vptr[i] += deltaV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		Wptr[i] += deltaW[i];
	
	delete []deltaU;
	delete []deltaV;
	delete []deltaW;
	deltaU = 0;
	deltaV = 0;
	deltaW = 0;

	return true;
}



bool ZQ_SmokeOctreeGrid3D::_solveGlobalClosedFlux(ZQ_SmokeGrid3D* srcGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(srcGrid == 0)
		return false;

	int uCount = Uinfos.size();
	int vCount = Vinfos.size();
	int wCount = Winfos.size();
	int lambdaCount = lamdainfos.size();

	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;

	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[dim];
	double* oldVelocity = new double[uCount+vCount+wCount];

	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	memset(b,0,sizeof(double)*dim);

	float*& Uptr = srcGrid->GetVelocityUptr();
	float*& Vptr = srcGrid->GetVelocityVptr();
	float*& Wptr = srcGrid->GetVelocityWptr();

	int size;
	int i_shift,j_shift,k_shift;
	double area[4] = {1,4,16,64};
	double weight[4] = {1,4,16,64};
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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		double val = 0;
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				val += Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
			}
		}
		
		if(Uinfos[i].index >= 0)
		{
			b[Uinfos[i].index] = 2*(val/(size*size))*area[lvl];
			oldVelocity[Uinfos[i].index] = val / (size*size);
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}

		double val = 0;
		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				val += Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
			}
		}
		
		if(Vinfos[i].index >= 0)
		{
			b[Vinfos[i].index] = 2*(val/(size*size))*area[lvl];
			oldVelocity[Vinfos[i].index] = val / (size*size);
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1, Winfos[i].level2);
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
		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		double val = 0;
		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				val += Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
			}
		}

		if(Winfos[i].index >= 0)
		{
			b[Winfos[i].index] = 2*(val/(size*size))*area[lvl];
			oldVelocity[Winfos[i].index] = val / (size*size);
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
		k_shift = lamdainfos[i].k1*size;
		double val = 0;
		for(int sk = 0;sk < size;sk++)
		{
			for(int sj = 0;sj < size;sj++)
			{
				if((i_shift > 0 && occupy[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift-1]) || i_shift == 0)
					val -= Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
				if((i_shift+size < width && occupy[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift+size]) || i_shift == width)
					val += Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift+size];
			}
		}
		
		for(int sk = 0;sk < size;sk++)
		{
			for(int si = 0;si < size;si++)
			{
				if((j_shift > 0 && occupy[(k_shift+sk)*height*width+(j_shift-1)*width+i_shift+si]) || j_shift == 0)
					val -= Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
				if((j_shift+size < height && occupy[(k_shift+sk)*height*width+(j_shift+size)*width+i_shift+si]) || j_shift == height)
					val += Vptr[(k_shift+sk)*(height+1)*width+(j_shift+size)*width+i_shift+si];
			}
		}

		for(int sj = 0;sj < size;sj++)
		{
			for(int si = 0;si < size;si++)
			{
				if((k_shift > 0 && occupy[(k_shift-1)*height*width+(j_shift+sj)*width+i_shift+si]) || k_shift == 0)
					val -= Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
				if((k_shift+size < depth && occupy[(k_shift+size)*height*width+(j_shift+sj)*width+i_shift+si]) || k_shift == depth)
					val += Wptr[(k_shift+size)*height*width+(j_shift+sj)*width+i_shift+si];
			}
		}

		b[lamdainfos[i].index] = val/volume[lvl];
	}

	double tol = 1e-9;
	int it = 0;

	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(A, b, x0, para.maxIter, tol, x, it, false);

	if(deltax)
	{
		delete []deltax;
	}
	deltax = new double[uCount+vCount+wCount];
	for(int i = 0;i < uCount+vCount+wCount;i++)
		deltax[i] = x[i] - oldVelocity[i];

	delete []x0;
	delete []x;
	delete []b;
	delete []oldVelocity;

	
	if(deltaU)
		delete []deltaU;
	if(deltaV)
		delete []deltaV;
	if(deltaW)
		delete []deltaW;
	
	deltaU = new float[(width+1)*height*depth];
	memset(deltaU,0,sizeof(float)*(width+1)*height*depth);
	deltaV = new float[width*(height+1)*depth];
	memset(deltaV,0,sizeof(float)*width*(height+1)*depth);
	deltaW = new float[width*height*(depth+1)];
	memset(deltaW,0,sizeof(float)*width*height*(depth+1));
	
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
			k_shift = Uinfos[i].k1*size;
			j_shift = Uinfos[i].j1*size;
			i_shift = (Uinfos[i].i1+1)*size;
		}
		else
		{
			k_shift = Uinfos[i].k2*size;
			j_shift = Uinfos[i].j2*size;
			i_shift = Uinfos[i].i2*size;
		}
		if(Uinfos[i].index >= 0)
		{
			for(int sk = 0;sk < size;sk++)
			{
				for(int sj = 0;sj < size;sj++)
				{
					deltaU[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift] += deltax[Uinfos[i].index];
				}
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
			k_shift = Vinfos[i].k1*size;
			j_shift = (Vinfos[i].j1+1)*size;
			i_shift = Vinfos[i].i1*size;
		}
		else
		{
			k_shift = Vinfos[i].k2*size;
			j_shift = Vinfos[i].j2*size;
			i_shift = Vinfos[i].i2*size;
		}
	
		if(Vinfos[i].index >= 0)
		{
			for(int sk = 0;sk < size;sk++)
			{
				for(int si = 0;si < size;si++)
				{
					deltaV[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si] += deltax[Vinfos[i].index];
				}
			}	
		}
	}

	for(int i = 0;i < wCount;i++)
	{
		int lvl = __min(Winfos[i].level1, Winfos[i].level2);
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
		if(lvl == Winfos[i].level1)
		{
			k_shift = (Winfos[i].k1+1)*size;
			j_shift = Winfos[i].j1*size;
			i_shift = Winfos[i].i1*size;
		}
		else
		{
			k_shift = Winfos[i].k2*size;
			j_shift = Winfos[i].j2*size;
			i_shift = Winfos[i].i2*size;
		}

		if(Winfos[i].index >= 0)
		{
			for(int sj = 0;sj < size;sj++)
			{
				for(int si = 0;si < size;si++)
				{
					deltaW[k_shift*height*width+(j_shift+sj)*width+i_shift+si] += deltax[Winfos[i].index];
				}
			}	
		}
	}
	delete []deltax;
	deltax = 0;

	for(int i = 0;i < (width+1)*height*depth;i++)
		Uptr[i] += deltaU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		Vptr[i] += deltaV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		Wptr[i] += deltaW[i];

	delete []deltaU;
	delete []deltaV;
	delete []deltaW;
	deltaU = 0;
	deltaV = 0;
	deltaW = 0;

	return true;
}

bool ZQ_SmokeOctreeGrid3D::_buildSubMatrix()
{
	int row2,col2;
	int row4,col4;
	int row8,col8;
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
	if(!ZQ_InitSubMatPoisson3D(2,&subMatrix2,&row2,&col2))
	{
		printf("init subMatrix2 fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}

	if(!ZQ_InitSubMatPoisson3D(4,&subMatrix4,&row4,&col4))
	{
		printf("init subMatrix4 fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}
	if(!ZQ_InitSubMatPoisson3D(8,&subMatrix8,&row8,&col8))
	{
		printf("init subMatrix2 fail!(%s):(%d)\n",__FILE__,__LINE__);
		return false;
	}
	return true;
}


void ZQ_SmokeOctreeGrid3D::_seedFilling(const int width, const int height, const int depth, int *distanceField)
{
	int head = -1;
	int tail = -1;
	int* queue_i = new int[width*height*depth];
	int* queue_j = new int[width*height*depth];
	int* queue_k = new int[width*height*depth];
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(distanceField[k*height*width+j*width+i] == 0)
				{
					tail++;
					queue_i[tail] = i;
					queue_j[tail] = j;
					queue_k[tail] = k;
				}
			}
		}
	}
	
	while(head < tail)
	{
		head ++;
		int cur_i = queue_i[head];
		int cur_j = queue_j[head];
		int cur_k = queue_k[head];
		int cur_dis = distanceField[cur_k*height*width+cur_j*width+cur_i];

		//x-
		if(cur_i-1 >= 0 && distanceField[cur_k*height*width+cur_j*width+cur_i-1] < 0)
		{
			distanceField[cur_k*height*width+cur_j*width+cur_i-1] = cur_dis+1;
			tail ++;
			queue_i[tail] = cur_i-1;
			queue_j[tail] = cur_j;
			queue_k[tail] = cur_k;
		}

		//x+
		if(cur_i+1 < width && distanceField[cur_k*height*width+cur_j*width+cur_i+1] < 0)
		{
			distanceField[cur_k*height*width+cur_j*width+cur_i+1] = cur_dis+1;
			tail ++;
			queue_i[tail] = cur_i+1;
			queue_j[tail] = cur_j;
			queue_k[tail] = cur_k;
		}

		//y-
		if(cur_j-1 >= 0 && distanceField[cur_k*height*width+(cur_j-1)*width+cur_i] < 0)
		{
			distanceField[cur_k*height*width+(cur_j-1)*width+cur_i] = cur_dis+1;
			tail ++;
			queue_i[tail] = cur_i;
			queue_j[tail] = cur_j-1;
			queue_k[tail] = cur_k;
		}

		//y+
		if(cur_j+1 < height && distanceField[cur_k*height*width+(cur_j+1)*width+cur_i] < 0)
		{
			distanceField[cur_k*height*width+(cur_j+1)*width+cur_i] = cur_dis+1;
			tail ++;
			queue_i[tail] = cur_i;
			queue_j[tail] = cur_j+1;
			queue_k[tail] = cur_k;
		}

		//z-
		if(cur_k-1 >= 0 && distanceField[(cur_k-1)*height*width+cur_j*width+cur_i] < 0)
		{
			distanceField[(cur_k-1)*height*width+cur_j*width+cur_i] = cur_dis+1;
			tail ++;
			queue_i[tail] = cur_i;
			queue_j[tail] = cur_j;
			queue_k[tail] = cur_k-1;
		}

		//z+
		if(cur_k+1 < depth && distanceField[(cur_k+1)*height*width+cur_j*width+cur_i] < 0)
		{
			distanceField[(cur_k+1)*height*width+cur_j*width+cur_i] = cur_dis+1;
			tail ++;
			queue_i[tail] = cur_i;
			queue_j[tail] = cur_j;
			queue_k[tail] = cur_k+1;
		}

		/*for(int sk = __max(0,cur_k-1);sk <= __min(depth-1,cur_k+1);sk++)
		{
			for(int sj = __max(0,cur_j-1);sj <= __min(height-1,cur_j+1);sj++)
			{
				for(int si = __max(0,cur_i-1);si <= __min(width-1,cur_i+1);si++)
				{
					if(si == cur_i && sj == cur_j && sk == cur_k)
						continue;
					if(distanceField[sk*height*width+sj*width+si] < 0 )
					{
						distanceField[sk*height*width+sj*width+si] = cur_dis+1;
						tail ++;
						queue_i[tail] = si;
						queue_j[tail] = sj;
						queue_k[tail] = sk;
					}
				}
			}
		}*/

	}
	/*for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(distanceField[k*height*width+j*width+i] < 0)
				{
					printf("distance[%d,%d,%d] = %d\n",i,j,k,distanceField[k*height*width+j*width+i]);
				}
			}
		}
	}*/
	
	delete []queue_i;
	delete []queue_j;
	delete []queue_k;
}




void ZQ_SmokeOctreeGrid3D::_find_special_label(const int width, const int height, const int depth, const bool* special_mask, int *special_label, int &region_num)
{

	for(int i = 0;i < width*height*depth;i++)
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

		int cur_i,cur_j,cur_k;
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(special_label[k*height*width+j*width+i] == -1)
					{
						special_label[k*height*width+j*width+i] = num_region+1;
						no_new_region = false;
						cur_i = i;
						cur_j = j;
						cur_k = k;
						break;
					}
				}
				if(!no_new_region)
					break;
			}
		}
		

		if(!no_new_region)
		{
			int* queue_i = new int[width*height*depth];
			int* queue_j = new int[width*height*depth];
			int* queue_k = new int[width*height*depth];
			int tail = -1;
			int head = 0;
			queue_i[head] = cur_i;
			queue_j[head] = cur_j;
			queue_k[head] = cur_k;

			while(tail < head)
			{
				tail++;
				int iii = queue_i[tail];
				int jjj = queue_j[tail];
				int kkk = queue_k[tail];
				if(kkk+1 < depth)
				{
					if(special_label[(kkk+1)*height*width+jjj*width+iii] == -1)
					{
						head ++;
						queue_i[head] = iii;
						queue_j[head] = jjj;
						queue_k[head] = kkk+1;
					}
				}
				if(kkk-1 >= 0)
				{
					if(special_label[(kkk-1)*height*width+jjj*width+iii] == -1)
					{
						head ++;
						queue_i[head] = iii;
						queue_j[head] = jjj;
						queue_k[head] = kkk-1;
					}
				}
				if(jjj+1 < height)
				{
					if(special_label[kkk*height*width+(jjj+1)*width+iii] == -1)
					{
						head ++;
						queue_i[head] = iii;
						queue_j[head] = jjj+1;
						queue_k[head] = kkk;
						special_label[kkk*height*width+(jjj+1)*width+iii] = num_region+1;
					}
				}
				if(jjj-1 >= 0)
				{
					if(special_label[kkk*height*width+(jjj-1)*width+iii] == -1)
					{
						head ++;
						queue_i[head] = iii;
						queue_j[head] = jjj-1;
						queue_k[head] = kkk;
						special_label[kkk*height*width+(jjj-1)*width+iii] = num_region+1;
					}
				}
				if(iii+1 < width)
				{
					if(special_label[kkk*height*width+jjj*width+(iii+1)] == -1)
					{
						head ++;
						queue_i[head] = iii+1;
						queue_j[head] = jjj;
						queue_k[head] = kkk;
						special_label[kkk*height*width+jjj*width+(iii+1)] = num_region+1;
					}
				}
				if(iii-1 >= 0)
				{
					if(special_label[kkk*height*width+jjj*width+(iii-1)] == -1)
					{
						head ++;
						queue_i[head] = iii-1;
						queue_j[head] = jjj;
						queue_k[head] = kkk;
						special_label[kkk*height*width+jjj*width+(iii-1)] = num_region+1;
					}
				}
			}
			delete []queue_i;
			delete []queue_j;
			delete []queue_k;

			num_region ++;

		}
	}
	region_num = num_region;
}

bool ZQ_SmokeOctreeGrid3D::_buildOctreeInfo(const ZQ_SmokeGrid3D *srcGrid, const ZQ_SmokeSimulationPara3D& para, const int *interestingRegionDist)
{
	int width = srcGrid->GetWidth();
	int height = srcGrid->GetHeight();
	int depth = srcGrid->GetDepth();

	const bool*& occupyPtr = srcGrid->GetOccupyPtr();

	int* distanceField = new int[width*height*depth];
	for(int i = 0;i < width*height*depth;i++)
		distanceField[i] = -1;

	bool hasSeed = false;

	for(int i = 0;i < width*height*depth;i++)
	{
		if(occupyPtr[i])
		{
			hasSeed = true;
			distanceField[i] = 0;
		}
	}
	if(interestingRegionDist)
	{
		for(int i = 0;i < width*height*depth;i++)
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
		//clock_t t1 = clock();
		_seedFilling(width,height,depth,distanceField);
		//clock_t t2 = clock();
		//printf("seed filling cost: %.3f\n",0.001*(t2-t1));
	}
	else
	{
		for(int i = 0;i < width*height*depth;i++)
		{
			distanceField[i] = width+height+depth;
		}
	}

	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;

	if(leafLevel0)
		delete []leafLevel0;
	if(leafLevel1)
		delete []leafLevel1;
	if(leafLevel2)
		delete []leafLevel2;
	if(leafLevel3)
		delete []leafLevel3;
	leafLevel0 = new bool[width*height*depth];
	leafLevel1 = new bool[levelWidth1*levelHeight1*levelDepth1];
	leafLevel2 = new bool[levelWidth2*levelHeight2*levelDepth2];
	leafLevel3 = new bool[levelWidth3*levelHeight3*levelDepth3];

	memset(leafLevel0,false,sizeof(bool)*width*height*depth);
	memset(leafLevel1,false,sizeof(bool)*levelWidth1*levelHeight1*levelDepth1);
	memset(leafLevel2,false,sizeof(bool)*levelWidth2*levelHeight2*levelDepth2);
	memset(leafLevel3,false,sizeof(bool)*levelWidth3*levelHeight3*levelDepth3);

	for(int k = 0;k < levelDepth3;k++)
	{
		for(int j = 0;j < levelHeight3;j++)
		{
			for(int i = 0;i < levelWidth3;i++)
			{
				bool flag = true;
				for(int sk = 0;sk < 8;sk++)
				{
					if(!flag)
						break;
					for(int sj = 0;sj < 8;sj++)
					{
						if(!flag)
							break;
						for(int si = 0;si < 8;si ++)
						{
							if(!flag)
								break;
							if(distanceField[(k*8+sk)*height*width+(j*8+sj)*width+(i*8+si)] <= level_dis_thresh3)
							{
								flag = false;
								break;
							}
						}
					}
				}
				
				leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i] = flag;
			}
		}
	}
	
	for(int k = 0;k < levelDepth2;k++)
	{
		for(int j = 0;j < levelHeight2;j++)
		{
			for(int i = 0;i < levelWidth2;i++)
			{
				if(leafLevel3[(k/2)*levelHeight3*levelWidth3+(j/2)*levelWidth3+(i/2)])
					leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i] = false;
				else
				{
					bool flag = true;
					for(int sk = 0;sk < 4;sk++)
					{
						if(!flag)
							break;
						for(int sj = 0;sj < 4;sj++)
						{
							if(!flag)
								break;
							for(int si = 0;si < 4;si++)
							{
								if(!flag)
									break;

								if(distanceField[(k*4+sk)*height*width+(j*4+sj)*width+(i*4+si)] <= level_dis_thresh2)
								{
									flag = false;
									break;
								}

							}
						}
					}
					
					leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i] = flag;
				}
			}
		}
	}
	
	for(int k = 0;k < levelDepth1;k++)
	{
		for(int j = 0;j < levelHeight1;j++)
		{
			for(int i = 0;i < levelWidth1;i++)
			{
				if(leafLevel3[(k/4)*levelHeight3*levelWidth3+(j/4)*levelWidth3+(i/4)] || leafLevel2[(k/2)*levelHeight2*levelWidth2+(j/2)*levelWidth2+(i/2)])
					leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i] = false;
				else
				{
					bool flag = true;
					for(int sk = 0;sk < 2;sk++)
					{
						if(!flag)
							break;
						for(int sj = 0;sj < 2;sj++)
						{
							if(!flag)
								break;
							for(int si = 0;si < 2;si++)
							{
								if(!flag)
									break;
								if(distanceField[(k*2+sk)*height*width+(j*2+sj)*width+(i*2+si)] <= level_dis_thresh1)
								{
									flag = false;
									break;
								}
							}
						}
					}
					leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i] = flag;
				}
			}
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(!leafLevel1[(k/2)*levelHeight1*levelWidth1+(j/2)*levelWidth1+(i/2)] && !leafLevel2[(k/4)*levelHeight2*levelWidth2+(j/4)*levelWidth2+(i/4)]
				&& !leafLevel3[(k/8)*levelHeight3*levelWidth3+(j/8)*levelWidth3+(i/8)] && !occupyPtr[k*height*width+j*width+i])
				{
					leafLevel0[k*height*width+j*width+i] = true;
				}
			}
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				bool has_use = true;
				has_use = occupyPtr[k*height*width+j*width+i] 
				|| leafLevel0[k*height*width+j*width+i] 
				|| leafLevel1[(k/2)*levelHeight1*levelWidth1+(j/2)*levelWidth1+(i/2)]
				|| leafLevel2[(k/4)*levelHeight2*levelWidth2+(j/4)*levelWidth2+(i/4)] 
				|| leafLevel3[(k/8)*levelHeight3*levelWidth3+(j/8)*levelWidth3+(i/8)];
				if(!has_use)
				{
					printf("%3d,%3d,%3d is not use\n",i,j,k);
				}
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
		specialRegion_mask = new bool[width*height*depth];
		memset(specialRegion_mask,0,sizeof(bool)*width*height*depth);
		if(specialRegion_label)
		{
			delete []specialRegion_label;
		}
		specialRegion_label = new int[width*height*depth];
		memset(specialRegion_label,0,sizeof(int)*width*height*depth);

		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(leafLevel0[k*height*width+j*width+i] || leafLevel1[(k/2)*(height/2)*(width/2)+(j/2)*(width/2)+(i/2)])
						specialRegion_mask[k*height*width+j*width+i] = true;
					else
						specialRegion_mask[k*height*width+j*width+i] = false;
				}
			}
		}
		
		
		_find_special_label(width,height,depth,specialRegion_mask,specialRegion_label,specialRegion_num);

		bool* new_Occupy = new bool[width*height*depth];

		if(!guiding_specialRegion_useCUDA)
		{
			for(int iregion = 0;iregion < specialRegion_num;iregion++)
			{
				for(int i = 0;i < width*height*depth;i++)
					new_Occupy[i] = occupyPtr[i] || (specialRegion_label[i] != iregion+1);
				ZQ::ZQ_PoissonSolver3D::BuildClosedPoisson<double>(width,height,depth,new_Occupy,&A,false);
				specialRegion_A.push_back(A);
			}
			delete []new_Occupy;
		}	
	}

	return true;
}

void ZQ_SmokeOctreeGrid3D::_selectInfos_lambda(const int width, const int height, const int depth, const bool* leafLevel0, const bool* leafLevel1, const bool* leafLevel2, const bool* leafLevel3, 
											   std::vector<Info>& lamdainfos)
{
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;
	lamdainfos.clear();

	/*level3 lambda*/
	for(int k = 0;k < levelDepth3;k++)
	{
		for(int j = 0;j < levelHeight3;j++)
		{
			for(int i = 0;i < levelWidth3;i++)
			{
				if(leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i])
				{
					Info tmpInfo(3,i,j,k,3,i,j,k);
					lamdainfos.push_back(tmpInfo);
				}
			}
		}
	}
	

	/*level2 lambda*/
	for(int k = 0;k < levelDepth2;k++)
	{
		for(int j = 0;j < levelHeight2;j++)
		{
			for(int i = 0;i < levelWidth2;i++)
			{
				if(leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
				{
					Info tmpInfo(2,i,j,k,2,i,j,k);
					lamdainfos.push_back(tmpInfo);
				}
			}
		}
	}
	

	/*level1 lambda*/
	for(int k = 0;k < levelDepth1;k++)
	{
		for(int j = 0;j < levelHeight1;j++)
		{
			for(int i = 0;i < levelWidth1;i++)
			{
				if(leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
				{
					Info tmpInfo(1,i,j,k,1,i,j,k);
					lamdainfos.push_back(tmpInfo);
				}
			}
		}
	}
	

	/*level0 lambda*/
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(leafLevel0[k*height*width+j*width+i])
				{
					Info tmpInfo(0,i,j,k,0,i,j,k);
					lamdainfos.push_back(tmpInfo);
				}
			}
		}
	}
	
}

void ZQ_SmokeOctreeGrid3D::_selectInfosForCPUSolvers(const int width, const int height, const int depth, const bool *leafLevel0, 
													 const bool *leafLevel1, const bool *leafLevel2, const bool *leafLevel3, 
													 std::vector<Info> &Uinfos, std::vector<Info> &Vinfos, std::vector<Info> &Winfos, std::vector<Info> &lamdainfos)
{
	int levelWidth1 = width/2;
	int levelWidth2 = width/4;
	int levelWidth3 = width/8;
	int levelHeight1 = height/2;
	int levelHeight2 = height/4;
	int levelHeight3 = height/8;
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;

	Uinfos.clear();
	Vinfos.clear();
	Winfos.clear();
	lamdainfos.clear();

	/* level3 U */
	for(int k = 0;k < levelDepth3;k++)
	{
		for(int j = 0;j < levelHeight3;j ++)
		{
			for(int i = 0;i <= levelWidth3;i++)
			{
				if(i == 0)
				{
					if(leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i])
					{
						Info tmpInfo(3,i-1,j,k,3,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
				else if(i == levelWidth3)
				{
					if(leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i-1])
					{
						Info tmpInfo(3,i-1,j,k,3,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i] && leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i-1])
					{
						Info tmpInfo(3,i-1,j,k,3,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
			}
		}
	}
	

	/*level2 U*/
	for(int k = 0;k < levelDepth2;k++)
	{
		for(int j = 0;j < levelHeight2;j++)
		{
			for(int i = 0;i <= levelWidth2;i++)
			{
				if(i == 0)
				{
					if(leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(3,i/2-1,j/2,k/2,2,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
				else if(i == levelWidth2)
				{
					if(leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i-1])
					{
						Info tmpInfo(2,i-1,j,k,3,i/2,j/2,k/2);
						Uinfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i-1] && leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(2,i-1,j,k,2,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
					else if( i%2 == 0 && leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i-1] && leafLevel3[(k/2)*levelHeight3*levelWidth3+(j/2)*levelWidth3+(i/2)])
					{
						Info tmpInfo(2,i-1,j,k,3,i/2,j/2,k/2);
						Uinfos.push_back(tmpInfo);
					}
					else if( i%2 == 0 && leafLevel3[(k/2)*levelHeight3*levelWidth3+(j/2)*levelWidth3+(i/2-1)] && leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(3,i/2-1,j/2,k/2,2,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}			
			}
		}
	}
	

	/*level1 U*/
	for(int k = 0;k < levelDepth1;k++)
	{
		for(int j = 0; j < levelHeight1;j++)
		{
			for(int i = 0;i <= levelWidth1;i++)
			{
				if(i == 0)
				{
					if(leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(3,i/4-1,j/4,k/4,1,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
				else if(i == levelWidth1)
				{
					if(leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i-1])
					{
						Info tmpInfo(1,i-1,j,k,3,i/4,j/4,k/4);
						Uinfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i-1] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(1,i-1,j,k,1,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%2 == 0 && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i-1] && leafLevel2[(k/2)*levelHeight2*levelWidth2+(j/2)*levelWidth2+(i/2)])
					{
						Info tmpInfo(1,i-1,j,k,2,i/2,j/2,k/2);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%4 == 0 && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i-1] && leafLevel3[(k/4)*levelHeight3*levelWidth3+(j/4)*levelWidth3+(i/4)])
					{
						Info tmpInfo(1,i-1,j,k,3,i/4,j/4,k/4);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%2 == 0 && leafLevel2[(k/2)*levelHeight2*levelWidth2+(j/2)*levelWidth2+(i/2-1)] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(2,i/2-1,j/2,k/2,1,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%4 == 0 && leafLevel3[(k/4)*levelHeight3*levelWidth3+(j/4)*levelWidth3+(i/4-1)] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(3,i/4-1,j/4,k/4,1,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
			}
		}
	}
	

	/*level0 U*/
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i <= width;i++)
			{
				if(i == 0)
				{
					if(leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(3,i/8-1,j/8,k/8,0,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
				else if(i == width)
				{
					if(leafLevel0[k*height*width+j*width+i-1])
					{
						Info tmpInfo(0,i-1,j,k,3,i/8,j/8,k/8);
						Uinfos.push_back(tmpInfo);
					}
				}
				else 
				{
					if(leafLevel0[k*height*width+j*width+i-1] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(0,i-1,j,k,0,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%2 == 0 && leafLevel0[k*height*width+j*width+i-1] && leafLevel1[(k/2)*levelHeight1*levelWidth1+(j/2)*levelWidth1+(i/2)])
					{
						Info tmpInfo(0,i-1,j,k,1,i/2,j/2,k/2);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%4 == 0 && leafLevel0[k*height*width+j*width+i-1] && leafLevel2[(k/4)*levelHeight2*levelWidth2+(j/4)*levelWidth2+(i/4)])
					{
						Info tmpInfo(0,i-1,j,k,2,i/4,j/4,k/4);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%8 == 0 && leafLevel0[k*height*width+j*width+i-1] && leafLevel3[(k/8)*levelHeight3*levelWidth3+(j/8)*levelWidth3+(i/8)])
					{
						Info tmpInfo(0,i-1,j,k,3,i/8,j/8,k/8);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%2 == 0 && leafLevel1[(k/2)*levelHeight1*levelWidth1+(j/2)*levelWidth1+(i/2-1)] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(1,i/2-1,j/2,k/2,0,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%4 == 0 && leafLevel2[(k/4)*levelHeight2*levelWidth2+(j/4)*levelWidth2+(i/4-1)] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(2,i/4-1,j/4,k/4,0,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
					else if(i%8 == 0 && leafLevel3[(k/8)*levelHeight3*levelWidth3+(j/8)*levelWidth3+(i/8-1)] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(3,i/8-1,j/8,k/8,0,i,j,k);
						Uinfos.push_back(tmpInfo);
					}
				}
			}
		}
	}
	

	/*level3 V*/
	for(int k = 0;k < levelDepth3;k++)
	{
		for(int j = 0;j <= levelHeight3;j++)
		{
			for(int i = 0;i < levelWidth3;i++)
			{
				if(j == 0)
				{
					if(leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i])
					{
						Info tmpInfo(3,i,j-1,k,3,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
				else if(j == levelHeight3)
				{
					if(leafLevel3[k*levelHeight3*levelWidth3+(j-1)*levelWidth3+i])
					{
						Info tmpInfo(3,i,j-1,k,3,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel3[k*levelHeight3*levelWidth3+(j-1)*levelWidth3+i] && leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i])
					{
						Info tmpInfo(3,i,j-1,k,3,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
			}
		}
	}
	

	/*level2 V*/
	for(int k = 0;k < levelDepth2;k++)
	{
		for(int j = 0;j <= levelHeight2;j++)
		{
			for(int i = 0;i < levelWidth2;i++)
			{
				if(j == 0)
				{
					if(leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(3,i/2,j/2-1,k/2,2,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
				else if(j == levelHeight2)
				{
					if(leafLevel2[k*levelHeight2*levelWidth2+(j-1)*levelWidth2+i])
					{
						Info tmpInfo(2,i,j-1,k,3,i/2,j/2,k/2);
						Vinfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel2[k*levelHeight2*levelWidth2+(j-1)*levelWidth2+i] && leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(2,i,j-1,k,2,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%2 == 0 && leafLevel2[k*levelHeight2*levelWidth2+(j-1)*levelWidth2+i] && leafLevel3[(k/2)*levelHeight3*levelWidth3+(j/2)*levelWidth3+(i/2)])
					{
						Info tmpInfo(2,i,j-1,k,3,i/2,j/2,k/2);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%2 == 0 && leafLevel3[(k/2)*levelHeight3*levelWidth3+(j/2-1)*levelWidth3+i/2] && leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(3,i/2,j/2-1,k/2,2,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
			}
		}
	}
	

	/*level1 V*/
	for(int k = 0;k < levelDepth1;k++)
	{
		for(int j = 0;j <= levelHeight1;j++)
		{
			for(int i = 0;i < levelWidth1;i++)
			{
				if(j == 0)
				{
					if(leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(3,i/4,j/4-1,k/4,1,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
				else if(j == levelHeight1)
				{
					if(leafLevel1[k*levelHeight1*levelWidth1+(j-1)*levelWidth1+i])
					{
						Info tmpInfo(1,i,j-1,k,3,i/4,j/4,k/4);
						Vinfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel1[k*levelHeight1*levelWidth1+(j-1)*levelWidth1+i] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(1,i,j-1,k,1,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%2 == 0 && leafLevel1[k*levelHeight1*levelWidth1+(j-1)*levelWidth1+i] && leafLevel2[(k/2)*levelHeight2*levelWidth2+(j/2)*levelWidth2+(i/2)])
					{
						Info tmpInfo(1,i,j-1,k,2,i/2,j/2,k/2);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%4 == 0 && leafLevel1[k*levelHeight1*levelWidth1+(j-1)*levelWidth1+i] && leafLevel3[(k/4)*levelHeight3*levelWidth3+(j/4)*levelWidth3+(i/4)])
					{
						Info tmpInfo(1,i,j-1,k,3,i/4,j/4,k/4);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%2 == 0 && leafLevel2[(k/2)*levelHeight2*levelWidth2+(j/2-1)*levelWidth2+i/2] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(2,i/2,j/2-1,k/2,1,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%4 == 0 && leafLevel3[(k/4)*levelHeight3*levelWidth3+(j/4-1)*levelWidth3+i/4] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(3,i/4,j/4-1,k/4,1,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
			}
		}
	}
	

	/*level0 V*/
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j <= height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(j == 0)
				{
					if(leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(3,i/8,j/8-1,k/8,0,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
				else if(j == height)
				{
					if(leafLevel0[k*height*width+(j-1)*width+i])
					{
						Info tmpInfo(0,i,j-1,k,3,i/8,j/8,k/8);
						Vinfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel0[k*height*width+(j-1)*width+i] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(0,i,j-1,k,0,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%2 == 0 && leafLevel0[k*height*width+(j-1)*width+i] && leafLevel1[(k/2)*levelHeight1*levelWidth1+(j/2)*levelWidth1+(i/2)])
					{
						Info tmpInfo(0,i,j-1,k,1,i/2,j/2,k/2);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%4 == 0 && leafLevel0[k*height*width+(j-1)*width+i] && leafLevel2[(k/4)*levelHeight2*levelWidth2+(j/4)*levelWidth2+(i/4)])
					{
						Info tmpInfo(0,i,j-1,k,2,i/4,j/4,k/4);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%8 == 0 && leafLevel0[k*height*width+(j-1)*width+i] && leafLevel3[(k/8)*levelHeight3*levelWidth3+(j/8)*levelWidth3+(i/8)])
					{
						Info tmpInfo(0,i,j-1,k,3,i/8,j/8,k/8);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%2 == 0 && leafLevel1[(k/2)*levelHeight1*levelWidth1+(j/2-1)*levelWidth1+i/2] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(1,i/2,j/2-1,k/2,0,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%4 == 0 && leafLevel2[(k/4)*levelHeight2*levelWidth2+(j/4-1)*levelWidth2+i/4] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(2,i/4,j/4-1,k/4,0,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
					else if(j%8 == 0 && leafLevel3[(k/8)*levelHeight3*levelWidth3+(j/8-1)*levelWidth3+i/8] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(3,i/8,j/8-1,k/8,0,i,j,k);
						Vinfos.push_back(tmpInfo);
					}
				}
			}
		}
	}
	

	/*level3 W*/
	for(int k = 0;k <= levelDepth3;k++)
	{
		for(int j = 0;j < levelHeight3;j++)
		{
			for(int i = 0;i < levelWidth3;i++)
			{
				if(k == 0)
				{
					if(leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i])
					{
						Info tmpInfo(3,i,j,k-1,3,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
				else if(k == levelDepth3)
				{
					if(leafLevel3[(k-1)*levelHeight3*levelWidth3+j*levelWidth3+i])
					{
						Info tmpInfo(3,i,j,k-1,3,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel3[(k-1)*levelHeight3*levelWidth3+j*levelWidth3+i] && leafLevel3[k*levelHeight3*levelWidth3+j*levelWidth3+i])
					{
						Info tmpInfo(3,i,j,k-1,3,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
			}
		}
	}


	/*level2 W*/
	for(int k = 0;k <= levelDepth2;k++)
	{
		for(int j = 0;j < levelHeight2;j++)
		{
			for(int i = 0;i < levelWidth2;i++)
			{
				if(k == 0)
				{
					if(leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(3,i/2,j/2,k/2-1,2,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
				else if(k == levelDepth2)
				{
					if(leafLevel2[(k-1)*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(2,i,j,k-1,3,i/2,j/2,k/2);
						Winfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel2[(k-1)*levelHeight2*levelWidth2+j*levelWidth2+i] && leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(2,i,j,k-1,2,i,j,k);
						Winfos.push_back(tmpInfo);
					}
					else if(k%2 == 0 && leafLevel2[(k-1)*levelHeight2*levelWidth2+j*levelWidth2+i] && leafLevel3[(k/2)*levelHeight3*levelWidth3+(j/2)*levelWidth3+(i/2)])
					{
						Info tmpInfo(2,i,j,k-1,3,i/2,j/2,k/2);
						Winfos.push_back(tmpInfo);
					}
					else if(k%2 == 0 && leafLevel3[(k/2-1)*levelHeight3*levelWidth3+j/2*levelWidth3+i/2] && leafLevel2[k*levelHeight2*levelWidth2+j*levelWidth2+i])
					{
						Info tmpInfo(3,i/2,j/2,k/2-1,2,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
			}
		}
	}


	/*level1 W*/
	for(int k = 0;k <= levelDepth1;k++)
	{
		for(int j = 0;j < levelHeight1;j++)
		{
			for(int i = 0;i < levelWidth1;i++)
			{
				if(k == 0)
				{
					if(leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(3,i/4,j/4,k/4-1,1,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
				else if(k == levelDepth1)
				{
					if(leafLevel1[(k-1)*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(1,i,j,k-1,3,i/4,j/4,k/4);
						Winfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel1[(k-1)*levelHeight1*levelWidth1+j*levelWidth1+i] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(1,i,j,k-1,1,i,j,k);
						Winfos.push_back(tmpInfo);
					}
					else if(k%2 == 0 && leafLevel1[(k-1)*levelHeight1*levelWidth1+j*levelWidth1+i] && leafLevel2[(k/2)*levelHeight2*levelWidth2+(j/2)*levelWidth2+(i/2)])
					{
						Info tmpInfo(1,i,j,k-1,2,i/2,j/2,k/2);
						Winfos.push_back(tmpInfo);
					}
					else if(k%4 == 0 && leafLevel1[k*levelHeight1*levelWidth1+(j-1)*levelWidth1+i] && leafLevel3[(k/4)*levelHeight3*levelWidth3+(j/4)*levelWidth3+(i/4)])
					{
						Info tmpInfo(1,i,j,k-1,3,i/4,j/4,k/4);
						Winfos.push_back(tmpInfo);
					}
					else if(k%2 == 0 && leafLevel2[(k/2-1)*levelHeight2*levelWidth2+(j/2)*levelWidth2+i/2] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(2,i/2,j/2,k/2-1,1,i,j,k);
						Winfos.push_back(tmpInfo);
					}
					else if(k%4 == 0 && leafLevel3[(k/4-1)*levelHeight3*levelWidth3+(j/4)*levelWidth3+i/4] && leafLevel1[k*levelHeight1*levelWidth1+j*levelWidth1+i])
					{
						Info tmpInfo(3,i/4,j/4-1,k/4,1,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
			}
		}
	}


	/*level0 W*/
	for(int k = 0;k <= depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(k == 0)
				{
					if(leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(3,i/8,j/8,k/8-1,0,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
				else if(k == depth)
				{
					if(leafLevel0[(k-1)*height*width+j*width+i])
					{
						Info tmpInfo(0,i,j,k-1,3,i/8,j/8,k/8);
						Winfos.push_back(tmpInfo);
					}
				}
				else
				{
					if(leafLevel0[(k-1)*height*width+j*width+i] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(0,i,j,k-1,0,i,j,k);
						Winfos.push_back(tmpInfo);
					}
					else if(k%2 == 0 && leafLevel0[(k-1)*height*width+j*width+i] && leafLevel1[(k/2)*levelHeight1*levelWidth1+(j/2)*levelWidth1+(i/2)])
					{
						Info tmpInfo(0,i,j,k-1,1,i/2,j/2,k/2);
						Winfos.push_back(tmpInfo);
					}
					else if(k%4 == 0 && leafLevel0[(k-1)*height*width+j*width+i] && leafLevel2[(k/4)*levelHeight2*levelWidth2+(j/4)*levelWidth2+(i/4)])
					{
						Info tmpInfo(0,i,j,k-1,2,i/4,j/4,k/4);
						Winfos.push_back(tmpInfo);
					}
					else if(k%8 == 0 && leafLevel0[(k-1)*height*width+j*width+i] && leafLevel3[(k/8)*levelHeight3*levelWidth3+(j/8)*levelWidth3+(i/8)])
					{
						Info tmpInfo(0,i,j,k-1,3,i/8,j/8,k/8);
						Winfos.push_back(tmpInfo);
					}
					else if(k%2 == 0 && leafLevel1[(k/2-1)*levelHeight1*levelWidth1+(j/2)*levelWidth1+i/2] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(1,i/2,j/2,k/2-1,0,i,j,k);
						Winfos.push_back(tmpInfo);
					}
					else if(k%4 == 0 && leafLevel2[(k/4-1)*levelHeight2*levelWidth2+(j/4)*levelWidth2+i/4] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(2,i/4,j/4,k/4-1,0,i,j,k);
						Winfos.push_back(tmpInfo);
					}
					else if(k%8 == 0 && leafLevel3[(k/8-1)*levelHeight3*levelWidth3+(j/8)*levelWidth3+i/8] && leafLevel0[k*height*width+j*width+i])
					{
						Info tmpInfo(3,i/8,j/8,k/8-1,0,i,j,k);
						Winfos.push_back(tmpInfo);
					}
				}
			}
		}
	}


	_selectInfos_lambda(width,height,depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos);
	
}


void ZQ_SmokeOctreeGrid3D::_selectInfos_for_CUDAOpenPoissonRedBlack3(
	const int width, const int height, const int depth, const bool* leafLevel0, const bool* leafLevel1, const bool* leafLevel2, const bool* leafLevel3, std::vector<Info>& lambdainfos,
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
	int levelDepth1 = depth/2;
	int levelDepth2 = depth/4;
	int levelDepth3 = depth/8;

	_selectInfos_lambda(width,height,depth,leafLevel0,leafLevel1,leafLevel2,leafLevel3,lamdainfos);

	/****************  for level0  ******************************/
	int level0_red_offset = 0;
	int level0_black_offset = 0;
	for(int z = 0;z < depth;z++)
	{
		for(int y = 0;y < height;y++)
		{
			for(int x = 0;x < width;x++)
			{
				if(!leafLevel0[z*height*width+y*width+x])
					continue;
				if((x+y+z)%2 == 0)
				{
					level0_index_red.push_back(x);
					level0_index_red.push_back(y);
					level0_index_red.push_back(z);
					int cur_offset = level0_red_offset;

					//x-
					if(x == 0)
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8-1);
						level0_neighborinfo_red.push_back(y/8);
						level0_neighborinfo_red.push_back(z/8);
						level0_red_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+y*width+x-1])
						{
							level0_neighborinfo_red.push_back(0);
							level0_neighborinfo_red.push_back(x-1);
							level0_neighborinfo_red.push_back(y);
							level0_neighborinfo_red.push_back(z);
							level0_red_offset ++;
						}
						else if(x%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y/2)*levelWidth1+(x/2)-1])
						{
							level0_neighborinfo_red.push_back(1);
							level0_neighborinfo_red.push_back(x/2-1);
							level0_neighborinfo_red.push_back(y/2);
							level0_neighborinfo_red.push_back(z/2);
							level0_red_offset ++;
						}
						else if(x%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y/4)*levelWidth2+(x/4)-1])
						{
							level0_neighborinfo_red.push_back(2);
							level0_neighborinfo_red.push_back(x/4-1);
							level0_neighborinfo_red.push_back(y/4);
							level0_neighborinfo_red.push_back(z/4);
							level0_red_offset ++;
						}
						else if(x%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y/8)*levelWidth3+(x/8)-1])
						{
							level0_neighborinfo_red.push_back(3);
							level0_neighborinfo_red.push_back(x/8-1);
							level0_neighborinfo_red.push_back(y/8);
							level0_neighborinfo_red.push_back(z/8);
							level0_red_offset ++;
						}
					}

					//x+
					if(x == width-1)
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back((x+1)/8);
						level0_neighborinfo_red.push_back(y/8);
						level0_neighborinfo_red.push_back(z/8);
						level0_red_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+y*width+(x+1)])
						{
							level0_neighborinfo_red.push_back(0);
							level0_neighborinfo_red.push_back(x+1);
							level0_neighborinfo_red.push_back(y);
							level0_neighborinfo_red.push_back(z);
							level0_red_offset ++;
						}
						else if((x+1)%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y/2)*levelWidth1+(x+1)/2])
						{
							level0_neighborinfo_red.push_back(1);
							level0_neighborinfo_red.push_back((x+1)/2);
							level0_neighborinfo_red.push_back(y/2);
							level0_neighborinfo_red.push_back(z/2);
							level0_red_offset ++;
						}
						else if((x+1)%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y/4)*levelWidth2+(x+1)/4])
						{
							level0_neighborinfo_red.push_back(2);
							level0_neighborinfo_red.push_back((x+1)/4);
							level0_neighborinfo_red.push_back(y/4);
							level0_neighborinfo_red.push_back(z/4);
							level0_red_offset ++;
						}
						else if((x+1)%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y/8)*levelWidth3+(x+1)/8])
						{
							level0_neighborinfo_red.push_back(3);
							level0_neighborinfo_red.push_back((x+1)/8);
							level0_neighborinfo_red.push_back(y/8);
							level0_neighborinfo_red.push_back(z/8);
							level0_red_offset ++;
						}
					}

					//y-
					if(y == 0)
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8);
						level0_neighborinfo_red.push_back(y/8-1);
						level0_neighborinfo_red.push_back(z/8);
						level0_red_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+(y-1)*width+x])
						{
							level0_neighborinfo_red.push_back(0);
							level0_neighborinfo_red.push_back(x);
							level0_neighborinfo_red.push_back(y-1);
							level0_neighborinfo_red.push_back(z);
							level0_red_offset ++;
						}
						else if(y%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y/2-1)*levelWidth1+x/2])
						{
							level0_neighborinfo_red.push_back(1);
							level0_neighborinfo_red.push_back(x/2);
							level0_neighborinfo_red.push_back(y/2-1);
							level0_neighborinfo_red.push_back(z/2);
							level0_red_offset ++;
						}
						else if(y%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y/4-1)*levelWidth2+x/4])
						{
							level0_neighborinfo_red.push_back(2);
							level0_neighborinfo_red.push_back(x/4);
							level0_neighborinfo_red.push_back(y/4-1);
							level0_neighborinfo_red.push_back(z/4);
							level0_red_offset ++;
						}
						else if(y%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y/8-1)*levelWidth3+x/8])
						{
							level0_neighborinfo_red.push_back(3);
							level0_neighborinfo_red.push_back(x/8);
							level0_neighborinfo_red.push_back(y/8-1);
							level0_neighborinfo_red.push_back(z/8);
							level0_red_offset ++;
						}
					}

					//y+
					if(y == height-1)
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8);
						level0_neighborinfo_red.push_back((y+1)/8);
						level0_neighborinfo_red.push_back(z/8);
						level0_red_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+(y+1)*width+x])
						{
							level0_neighborinfo_red.push_back(0);
							level0_neighborinfo_red.push_back(x);
							level0_neighborinfo_red.push_back(y+1);
							level0_neighborinfo_red.push_back(z);
							level0_red_offset ++;
						}
						else if((y+1)%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y+1)/2*levelWidth1+x/2])
						{
							level0_neighborinfo_red.push_back(1);
							level0_neighborinfo_red.push_back(x/2);
							level0_neighborinfo_red.push_back((y+1)/2);
							level0_neighborinfo_red.push_back(z/2);
							level0_red_offset ++;
						}
						else if((y+1)%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y+1)/4*levelWidth2+x/4])
						{
							level0_neighborinfo_red.push_back(2);
							level0_neighborinfo_red.push_back(x/4);
							level0_neighborinfo_red.push_back((y+1)/4);
							level0_neighborinfo_red.push_back(z/4);
							level0_red_offset ++;
						}
						else if((y+1)%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y+1)/8*levelWidth3+x/8])
						{
							level0_neighborinfo_red.push_back(3);
							level0_neighborinfo_red.push_back(x/8);
							level0_neighborinfo_red.push_back((y+1)/8);
							level0_neighborinfo_red.push_back(z/8);
							level0_red_offset ++;
						}
					}

					//z-
					if(z == 0)
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8);
						level0_neighborinfo_red.push_back(y/8);
						level0_neighborinfo_red.push_back(z/8-1);
						level0_red_offset ++;
					}
					else
					{
						if(leafLevel0[(z-1)*height*width+y*width+x])
						{
							level0_neighborinfo_red.push_back(0);
							level0_neighborinfo_red.push_back(x);
							level0_neighborinfo_red.push_back(y);
							level0_neighborinfo_red.push_back(z-1);
							level0_red_offset ++;
						}
						else if(z%2 == 0 && leafLevel1[(z/2-1)*levelHeight1*levelWidth1+(y/2)*levelWidth1+x/2])
						{
							level0_neighborinfo_red.push_back(1);
							level0_neighborinfo_red.push_back(x/2);
							level0_neighborinfo_red.push_back(y/2);
							level0_neighborinfo_red.push_back(z/2-1);
							level0_red_offset ++;
						}
						else if(z%4 == 0 && leafLevel2[(z/4-1)*levelHeight2*levelWidth2+(y/4)*levelWidth2+x/4])
						{
							level0_neighborinfo_red.push_back(2);
							level0_neighborinfo_red.push_back(x/4);
							level0_neighborinfo_red.push_back(y/4);
							level0_neighborinfo_red.push_back(z/4-1);
							level0_red_offset ++;
						}
						else if(z%8 == 0 && leafLevel3[(z/8-1)*levelHeight3*levelWidth3+(y/8)*levelWidth3+x/8])
						{
							level0_neighborinfo_red.push_back(3);
							level0_neighborinfo_red.push_back(x/8);
							level0_neighborinfo_red.push_back(y/8);
							level0_neighborinfo_red.push_back(z/8-1);
							level0_red_offset ++;
						}
					}

					//z+
					if(z == depth-1)
					{
						level0_neighborinfo_red.push_back(3);
						level0_neighborinfo_red.push_back(x/8);
						level0_neighborinfo_red.push_back(y/8);
						level0_neighborinfo_red.push_back((z+1)/8);
						level0_red_offset ++;
					}
					else
					{
						if(leafLevel0[(z+1)*height*width+y*width+x])
						{
							level0_neighborinfo_red.push_back(0);
							level0_neighborinfo_red.push_back(x);
							level0_neighborinfo_red.push_back(y);
							level0_neighborinfo_red.push_back(z+1);
							level0_red_offset ++;
						}
						else if((z+1)%2 == 0 && leafLevel1[(z+1)/2*levelHeight1*levelWidth1+y/2*levelWidth1+x/2])
						{
							level0_neighborinfo_red.push_back(1);
							level0_neighborinfo_red.push_back(x/2);
							level0_neighborinfo_red.push_back(y/2);
							level0_neighborinfo_red.push_back((z+1)/2);
							level0_red_offset ++;
						}
						else if((z+1)%4 == 0 && leafLevel2[(z+1)/4*levelHeight2*levelWidth2+y/4*levelWidth2+x/4])
						{
							level0_neighborinfo_red.push_back(2);
							level0_neighborinfo_red.push_back(x/4);
							level0_neighborinfo_red.push_back(y/4);
							level0_neighborinfo_red.push_back((z+1)/4);
							level0_red_offset ++;
						}
						else if((z+1)%8 == 0 && leafLevel3[(z+1)/8*levelHeight3*levelWidth3+y/8*levelWidth3+x/8])
						{
							level0_neighborinfo_red.push_back(3);
							level0_neighborinfo_red.push_back(x/8);
							level0_neighborinfo_red.push_back(y/8);
							level0_neighborinfo_red.push_back((z+1)/8);
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
					level0_index_black.push_back(z);
					int cur_offset = level0_black_offset;

					//x-
					if(x == 0)
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8-1);
						level0_neighborinfo_black.push_back(y/8);
						level0_neighborinfo_black.push_back(z/8);
						level0_black_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+y*width+x-1])
						{
							level0_neighborinfo_black.push_back(0);
							level0_neighborinfo_black.push_back(x-1);
							level0_neighborinfo_black.push_back(y);
							level0_neighborinfo_black.push_back(z);
							level0_black_offset ++;
						}
						else if(x%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y/2)*levelWidth1+(x/2)-1])
						{
							level0_neighborinfo_black.push_back(1);
							level0_neighborinfo_black.push_back(x/2-1);
							level0_neighborinfo_black.push_back(y/2);
							level0_neighborinfo_black.push_back(z/2);
							level0_black_offset ++;
						}
						else if(x%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y/4)*levelWidth2+(x/4)-1])
						{
							level0_neighborinfo_black.push_back(2);
							level0_neighborinfo_black.push_back(x/4-1);
							level0_neighborinfo_black.push_back(y/4);
							level0_neighborinfo_black.push_back(z/4);
							level0_black_offset ++;
						}
						else if(x%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y/8)*levelWidth3+(x/8)-1])
						{
							level0_neighborinfo_black.push_back(3);
							level0_neighborinfo_black.push_back(x/8-1);
							level0_neighborinfo_black.push_back(y/8);
							level0_neighborinfo_black.push_back(z/8);
							level0_black_offset ++;
						}
					}

					//x+
					if(x == width-1)
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back((x+1)/8);
						level0_neighborinfo_black.push_back(y/8);
						level0_neighborinfo_black.push_back(z/8);
						level0_black_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+y*width+(x+1)])
						{
							level0_neighborinfo_black.push_back(0);
							level0_neighborinfo_black.push_back(x+1);
							level0_neighborinfo_black.push_back(y);
							level0_neighborinfo_black.push_back(z);
							level0_black_offset ++;
						}
						else if((x+1)%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y/2)*levelWidth1+(x+1)/2])
						{
							level0_neighborinfo_black.push_back(1);
							level0_neighborinfo_black.push_back((x+1)/2);
							level0_neighborinfo_black.push_back(y/2);
							level0_neighborinfo_black.push_back(z/2);
							level0_black_offset ++;
						}
						else if((x+1)%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y/4)*levelWidth2+(x+1)/4])
						{
							level0_neighborinfo_black.push_back(2);
							level0_neighborinfo_black.push_back((x+1)/4);
							level0_neighborinfo_black.push_back(y/4);
							level0_neighborinfo_black.push_back(z/4);
							level0_black_offset ++;
						}
						else if((x+1)%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y/8)*levelWidth3+(x+1)/8])
						{
							level0_neighborinfo_black.push_back(3);
							level0_neighborinfo_black.push_back((x+1)/8);
							level0_neighborinfo_black.push_back(y/8);
							level0_neighborinfo_black.push_back(z/8);
							level0_black_offset ++;
						}
					}

					//y-
					if(y == 0)
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8);
						level0_neighborinfo_black.push_back(y/8-1);
						level0_neighborinfo_black.push_back(z/8);
						level0_black_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+(y-1)*width+x])
						{
							level0_neighborinfo_black.push_back(0);
							level0_neighborinfo_black.push_back(x);
							level0_neighborinfo_black.push_back(y-1);
							level0_neighborinfo_black.push_back(z);
							level0_black_offset ++;
						}
						else if(y%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y/2-1)*levelWidth1+x/2])
						{
							level0_neighborinfo_black.push_back(1);
							level0_neighborinfo_black.push_back(x/2);
							level0_neighborinfo_black.push_back(y/2-1);
							level0_neighborinfo_black.push_back(z/2);
							level0_black_offset ++;
						}
						else if(y%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y/4-1)*levelWidth2+x/4])
						{
							level0_neighborinfo_black.push_back(2);
							level0_neighborinfo_black.push_back(x/4);
							level0_neighborinfo_black.push_back(y/4-1);
							level0_neighborinfo_black.push_back(z/4);
							level0_black_offset ++;
						}
						else if(y%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y/8-1)*levelWidth3+x/8])
						{
							level0_neighborinfo_black.push_back(3);
							level0_neighborinfo_black.push_back(x/8);
							level0_neighborinfo_black.push_back(y/8-1);
							level0_neighborinfo_black.push_back(z/8);
							level0_black_offset ++;
						}
					}

					//y+
					if(y == height-1)
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8);
						level0_neighborinfo_black.push_back((y+1)/8);
						level0_neighborinfo_black.push_back(z/8);
						level0_black_offset ++;
					}
					else
					{
						if(leafLevel0[z*height*width+(y+1)*width+x])
						{
							level0_neighborinfo_black.push_back(0);
							level0_neighborinfo_black.push_back(x);
							level0_neighborinfo_black.push_back(y+1);
							level0_neighborinfo_black.push_back(z);
							level0_black_offset ++;
						}
						else if((y+1)%2 == 0 && leafLevel1[(z/2)*levelHeight1*levelWidth1+(y+1)/2*levelWidth1+x/2])
						{
							level0_neighborinfo_black.push_back(1);
							level0_neighborinfo_black.push_back(x/2);
							level0_neighborinfo_black.push_back((y+1)/2);
							level0_neighborinfo_black.push_back(z/2);
							level0_black_offset ++;
						}
						else if((y+1)%4 == 0 && leafLevel2[(z/4)*levelHeight2*levelWidth2+(y+1)/4*levelWidth2+x/4])
						{
							level0_neighborinfo_black.push_back(2);
							level0_neighborinfo_black.push_back(x/4);
							level0_neighborinfo_black.push_back((y+1)/4);
							level0_neighborinfo_black.push_back(z/4);
							level0_black_offset ++;
						}
						else if((y+1)%8 == 0 && leafLevel3[(z/8)*levelHeight3*levelWidth3+(y+1)/8*levelWidth3+x/8])
						{
							level0_neighborinfo_black.push_back(3);
							level0_neighborinfo_black.push_back(x/8);
							level0_neighborinfo_black.push_back((y+1)/8);
							level0_neighborinfo_black.push_back(z/8);
							level0_black_offset ++;
						}
					}

					//z-
					if(z == 0)
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8);
						level0_neighborinfo_black.push_back(y/8);
						level0_neighborinfo_black.push_back(z/8-1);
						level0_black_offset ++;
					}
					else
					{
						if(leafLevel0[(z-1)*height*width+y*width+x])
						{
							level0_neighborinfo_black.push_back(0);
							level0_neighborinfo_black.push_back(x);
							level0_neighborinfo_black.push_back(y);
							level0_neighborinfo_black.push_back(z-1);
							level0_black_offset ++;
						}
						else if(z%2 == 0 && leafLevel1[(z/2-1)*levelHeight1*levelWidth1+(y/2)*levelWidth1+x/2])
						{
							level0_neighborinfo_black.push_back(1);
							level0_neighborinfo_black.push_back(x/2);
							level0_neighborinfo_black.push_back(y/2);
							level0_neighborinfo_black.push_back(z/2-1);
							level0_black_offset ++;
						}
						else if(z%4 == 0 && leafLevel2[(z/4-1)*levelHeight2*levelWidth2+(y/4)*levelWidth2+x/4])
						{
							level0_neighborinfo_black.push_back(2);
							level0_neighborinfo_black.push_back(x/4);
							level0_neighborinfo_black.push_back(y/4);
							level0_neighborinfo_black.push_back(z/4-1);
							level0_black_offset ++;
						}
						else if(z%8 == 0 && leafLevel3[(z/8-1)*levelHeight3*levelWidth3+(y/8)*levelWidth3+x/8])
						{
							level0_neighborinfo_black.push_back(3);
							level0_neighborinfo_black.push_back(x/8);
							level0_neighborinfo_black.push_back(y/8);
							level0_neighborinfo_black.push_back(z/8-1);
							level0_black_offset ++;
						}
					}

					//z+
					if(z == depth-1)
					{
						level0_neighborinfo_black.push_back(3);
						level0_neighborinfo_black.push_back(x/8);
						level0_neighborinfo_black.push_back(y/8);
						level0_neighborinfo_black.push_back((z+1)/8);
						level0_black_offset ++;
					}
					else
					{
						if(leafLevel0[(z+1)*height*width+y*width+x])
						{
							level0_neighborinfo_black.push_back(0);
							level0_neighborinfo_black.push_back(x);
							level0_neighborinfo_black.push_back(y);
							level0_neighborinfo_black.push_back(z+1);
							level0_black_offset ++;
						}
						else if((z+1)%2 == 0 && leafLevel1[(z+1)/2*levelHeight1*levelWidth1+y/2*levelWidth1+x/2])
						{
							level0_neighborinfo_black.push_back(1);
							level0_neighborinfo_black.push_back(x/2);
							level0_neighborinfo_black.push_back(y/2);
							level0_neighborinfo_black.push_back((z+1)/2);
							level0_black_offset ++;
						}
						else if((z+1)%4 == 0 && leafLevel2[(z+1)/4*levelHeight2*levelWidth2+y/4*levelWidth2+x/4])
						{
							level0_neighborinfo_black.push_back(2);
							level0_neighborinfo_black.push_back(x/4);
							level0_neighborinfo_black.push_back(y/4);
							level0_neighborinfo_black.push_back((z+1)/4);
							level0_black_offset ++;
						}
						else if((z+1)%8 == 0 && leafLevel3[(z+1)/8*levelHeight3*levelWidth3+y/8*levelWidth3+x/8])
						{
							level0_neighborinfo_black.push_back(3);
							level0_neighborinfo_black.push_back(x/8);
							level0_neighborinfo_black.push_back(y/8);
							level0_neighborinfo_black.push_back((z+1)/8);
							level0_black_offset ++;
						}
					}

					level0_index_black.push_back(level0_black_offset-cur_offset);
					level0_index_black.push_back(cur_offset);
				}
			}
		}
	}
	

	/************************** for level1 **************************************/
	int level1_red_offset = 0;
	int level1_black_offset = 0;
	for(int z = 0;z < levelDepth1;z++)
	{
		for(int y = 0;y < levelHeight1;y++)
		{
			for(int x = 0;x < levelWidth1;x++)
			{
				if(!leafLevel1[z*levelHeight1*levelWidth1+y*levelWidth1+x])
					continue;

				if((x+y+z)%2 == 0)
				{
					level1_index_red.push_back(x);
					level1_index_red.push_back(y);
					level1_index_red.push_back(z);
					int cur_offset = level1_red_offset;

					//x-
					if(x == 0)
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4-1);
						level1_neighborinfo_red.push_back(y/4);
						level1_neighborinfo_red.push_back(z/4);
						level1_red_offset ++;
					}
					else
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+y*levelWidth1+x-1])
						{
							level1_neighborinfo_red.push_back(1);
							level1_neighborinfo_red.push_back(x-1);
							level1_neighborinfo_red.push_back(y);
							level1_neighborinfo_red.push_back(z);
							level1_red_offset ++;
						}
						else if(x%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y/2)*levelWidth2+x/2-1])
						{
							level1_neighborinfo_red.push_back(2);
							level1_neighborinfo_red.push_back(x/2-1);
							level1_neighborinfo_red.push_back(y/2);
							level1_neighborinfo_red.push_back(z/2);
							level1_red_offset ++;
						}
						else if(x%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y/4)*levelWidth3+x/4-1])
						{
							level1_neighborinfo_red.push_back(3);
							level1_neighborinfo_red.push_back(x/4-1);
							level1_neighborinfo_red.push_back(y/4);
							level1_neighborinfo_red.push_back(z/4);
							level1_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2;cur_z < z*2+2;cur_z++)
							{
								for(int cur_y = y*2;cur_y < y*2+2;cur_y++)
								{
									if(leafLevel0[cur_z*height*width+cur_y*width+x*2-1])
									{
										level1_neighborinfo_red.push_back(0);
										level1_neighborinfo_red.push_back(x*2-1);
										level1_neighborinfo_red.push_back(cur_y);
										level1_neighborinfo_red.push_back(cur_z);
										level1_red_offset ++;
									}
								}
							}
						}
					}

					//x+
					if(x == levelWidth1-1)
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back((x+1)/4);
						level1_neighborinfo_red.push_back(y/4);
						level1_neighborinfo_red.push_back(z/4);
						level1_red_offset ++;
					}
					else
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+y*levelWidth1+x+1])
						{
							level1_neighborinfo_red.push_back(1);
							level1_neighborinfo_red.push_back(x+1);
							level1_neighborinfo_red.push_back(y);
							level1_neighborinfo_red.push_back(z);
							level1_red_offset ++;
						}
						else if((x+1)%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y/2)*levelWidth2+(x+1)/2])
						{
							level1_neighborinfo_red.push_back(2);
							level1_neighborinfo_red.push_back((x+1)/2);
							level1_neighborinfo_red.push_back(y/2);
							level1_neighborinfo_red.push_back(z/2);
							level1_red_offset ++;
						}
						else if((x+1)%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y/4)*levelWidth3+(x+1)/4])
						{
							level1_neighborinfo_red.push_back(3);
							level1_neighborinfo_red.push_back((x+1)/4);
							level1_neighborinfo_red.push_back(y/4);
							level1_neighborinfo_red.push_back(z/4);
							level1_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2; cur_z < z*2+2;cur_z++)
							{
								for(int cur_y = y*2; cur_y < y*2+2;cur_y++)
								{
									if(leafLevel0[cur_z*height*width+cur_y*width+(x+1)*2])
									{
										level1_neighborinfo_red.push_back(0);
										level1_neighborinfo_red.push_back((x+1)*2);
										level1_neighborinfo_red.push_back(cur_y);
										level1_neighborinfo_red.push_back(cur_z);
										level1_red_offset ++;
									}
								}
							}
						}
					}

					//y-
					if(y == 0)
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4);
						level1_neighborinfo_red.push_back(y/4-1);
						level1_neighborinfo_red.push_back(z/4);
						level1_red_offset ++;
					}
					else 
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+(y-1)*levelWidth1+x])
						{
							level1_neighborinfo_red.push_back(1);
							level1_neighborinfo_red.push_back(x);
							level1_neighborinfo_red.push_back(y-1);
							level1_neighborinfo_red.push_back(z);
							level1_red_offset ++;
						}
						else if(y%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y/2-1)*levelWidth2+x/2])
						{
							level1_neighborinfo_red.push_back(2);
							level1_neighborinfo_red.push_back(x/2);
							level1_neighborinfo_red.push_back(y/2-1);
							level1_neighborinfo_red.push_back(z/2);
							level1_red_offset ++;
						}
						else if(y%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y/4-1)*levelWidth3+x/4])
						{
							level1_neighborinfo_red.push_back(3);
							level1_neighborinfo_red.push_back(x/4);
							level1_neighborinfo_red.push_back(y/4-1);
							level1_neighborinfo_red.push_back(z/4);
							level1_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2; cur_z < z*2+2;cur_z++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[cur_z*height*width+(y*2-1)*width+cur_x])
									{
										level1_neighborinfo_red.push_back(0);
										level1_neighborinfo_red.push_back(cur_x);
										level1_neighborinfo_red.push_back(y*2-1);
										level1_neighborinfo_red.push_back(cur_z);
										level1_red_offset ++;
									}
								}
							}
						}
					}

					//y+
					if(y == levelHeight1-1)
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4);
						level1_neighborinfo_red.push_back((y+1)/4);
						level1_neighborinfo_red.push_back(z/4);
						level1_red_offset ++;
					}
					else
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+(y+1)*levelWidth1+x])
						{
							level1_neighborinfo_red.push_back(1);
							level1_neighborinfo_red.push_back(x);
							level1_neighborinfo_red.push_back(y+1);
							level1_neighborinfo_red.push_back(z);
							level1_red_offset ++;
						}
						else if((y+1)%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y+1)/2*levelWidth2+x/2])
						{
							level1_neighborinfo_red.push_back(2);
							level1_neighborinfo_red.push_back(x/2);
							level1_neighborinfo_red.push_back((y+1)/2);
							level1_neighborinfo_red.push_back(z/2);
							level1_red_offset ++;
						}
						else if((y+1)%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y+1)/4*levelWidth3+x/4])
						{
							level1_neighborinfo_red.push_back(3);
							level1_neighborinfo_red.push_back(x/4);
							level1_neighborinfo_red.push_back((y+1)/4);
							level1_neighborinfo_red.push_back(z/4);
							level1_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2;cur_z < z*2+2;cur_z++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[cur_z*height*width+(y+1)*2*width+cur_x])
									{
										level1_neighborinfo_red.push_back(0);
										level1_neighborinfo_red.push_back(cur_x);
										level1_neighborinfo_red.push_back((y+1)*2);
										level1_neighborinfo_red.push_back(cur_z);
										level1_red_offset ++;
									}
								}
							}
						}
					}

					//z-
					if(z == 0)
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4);
						level1_neighborinfo_red.push_back(y/4);
						level1_neighborinfo_red.push_back(z/4-1);
						level1_red_offset ++;
					}
					else 
					{
						if(leafLevel1[(z-1)*levelHeight1*levelWidth1+y*levelWidth1+x])
						{
							level1_neighborinfo_red.push_back(1);
							level1_neighborinfo_red.push_back(x);
							level1_neighborinfo_red.push_back(y);
							level1_neighborinfo_red.push_back(z-1);
							level1_red_offset ++;
						}
						else if(z%2 == 0 && leafLevel2[(z/2-1)*levelHeight2*levelWidth2+(y/2)*levelWidth2+x/2])
						{
							level1_neighborinfo_red.push_back(2);
							level1_neighborinfo_red.push_back(x/2);
							level1_neighborinfo_red.push_back(y/2);
							level1_neighborinfo_red.push_back(z/2-1);
							level1_red_offset ++;
						}
						else if(z%4 == 0 && leafLevel3[(z/4-1)*levelHeight3*levelWidth3+(y/4)*levelWidth3+x/4])
						{
							level1_neighborinfo_red.push_back(3);
							level1_neighborinfo_red.push_back(x/4);
							level1_neighborinfo_red.push_back(y/4);
							level1_neighborinfo_red.push_back(z/4-1);
							level1_red_offset ++;
						}
						else //smaller
						{
							for(int cur_y = y*2; cur_y < y*2+2;cur_y++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[(z*2-1)*height*width+cur_y*width+cur_x])
									{
										level1_neighborinfo_red.push_back(0);
										level1_neighborinfo_red.push_back(cur_x);
										level1_neighborinfo_red.push_back(cur_y);
										level1_neighborinfo_red.push_back(z*2-1);
										level1_red_offset ++;
									}
								}
							}
						}
					}

					//z+
					if(y == levelDepth1-1)
					{
						level1_neighborinfo_red.push_back(3);
						level1_neighborinfo_red.push_back(x/4);
						level1_neighborinfo_red.push_back(y/4);
						level1_neighborinfo_red.push_back((z+1)/4);
						level1_red_offset ++;
					}
					else
					{
						if(leafLevel1[(z+1)*levelHeight1*levelWidth1+y*levelWidth1+x])
						{
							level1_neighborinfo_red.push_back(1);
							level1_neighborinfo_red.push_back(x);
							level1_neighborinfo_red.push_back(y);
							level1_neighborinfo_red.push_back(z+1);
							level1_red_offset ++;
						}
						else if((z+1)%2 == 0 && leafLevel2[(z+1)/2*levelHeight2*levelWidth2+(y/2)*levelWidth2+x/2])
						{
							level1_neighborinfo_red.push_back(2);
							level1_neighborinfo_red.push_back(x/2);
							level1_neighborinfo_red.push_back(y/2);
							level1_neighborinfo_red.push_back((z+1)/2);
							level1_red_offset ++;
						}
						else if((z+1)%4 == 0 && leafLevel3[(z+1)/4*levelHeight3*levelWidth3+(y)/4*levelWidth3+x/4])
						{
							level1_neighborinfo_red.push_back(3);
							level1_neighborinfo_red.push_back(x/4);
							level1_neighborinfo_red.push_back(y/4);
							level1_neighborinfo_red.push_back((z+1)/4);
							level1_red_offset ++;
						}
						else //smaller
						{
							for(int cur_y = y*2;cur_y < y*2+2;cur_y++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[(z+1)*2*height*width+cur_y*width+cur_x])
									{
										level1_neighborinfo_red.push_back(0);
										level1_neighborinfo_red.push_back(cur_x);
										level1_neighborinfo_red.push_back(cur_y);
										level1_neighborinfo_red.push_back((z+1)*2);
										level1_red_offset ++;
									}
								}
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
					level1_index_black.push_back(z);
					int cur_offset = level1_black_offset;

					//x-
					if(x == 0)
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4-1);
						level1_neighborinfo_black.push_back(y/4);
						level1_neighborinfo_black.push_back(z/4);
						level1_black_offset ++;
					}
					else
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+y*levelWidth1+x-1])
						{
							level1_neighborinfo_black.push_back(1);
							level1_neighborinfo_black.push_back(x-1);
							level1_neighborinfo_black.push_back(y);
							level1_neighborinfo_black.push_back(z);
							level1_black_offset ++;
						}
						else if(x%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y/2)*levelWidth2+x/2-1])
						{
							level1_neighborinfo_black.push_back(2);
							level1_neighborinfo_black.push_back(x/2-1);
							level1_neighborinfo_black.push_back(y/2);
							level1_neighborinfo_black.push_back(z/2);
							level1_black_offset ++;
						}
						else if(x%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y/4)*levelWidth3+x/4-1])
						{
							level1_neighborinfo_black.push_back(3);
							level1_neighborinfo_black.push_back(x/4-1);
							level1_neighborinfo_black.push_back(y/4);
							level1_neighborinfo_black.push_back(z/4);
							level1_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2;cur_z < z*2+2;cur_z++)
							{
								for(int cur_y = y*2;cur_y < y*2+2;cur_y++)
								{
									if(leafLevel0[cur_z*height*width+cur_y*width+x*2-1])
									{
										level1_neighborinfo_black.push_back(0);
										level1_neighborinfo_black.push_back(x*2-1);
										level1_neighborinfo_black.push_back(cur_y);
										level1_neighborinfo_black.push_back(cur_z);
										level1_black_offset ++;
									}
								}
							}
						}
					}

					//x+
					if(x == levelWidth1-1)
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back((x+1)/4);
						level1_neighborinfo_black.push_back(y/4);
						level1_neighborinfo_black.push_back(z/4);
						level1_black_offset ++;
					}
					else
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+y*levelWidth1+x+1])
						{
							level1_neighborinfo_black.push_back(1);
							level1_neighborinfo_black.push_back(x+1);
							level1_neighborinfo_black.push_back(y);
							level1_neighborinfo_black.push_back(z);
							level1_black_offset ++;
						}
						else if((x+1)%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y/2)*levelWidth2+(x+1)/2])
						{
							level1_neighborinfo_black.push_back(2);
							level1_neighborinfo_black.push_back((x+1)/2);
							level1_neighborinfo_black.push_back(y/2);
							level1_neighborinfo_black.push_back(z/2);
							level1_black_offset ++;
						}
						else if((x+1)%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y/4)*levelWidth3+(x+1)/4])
						{
							level1_neighborinfo_black.push_back(3);
							level1_neighborinfo_black.push_back((x+1)/4);
							level1_neighborinfo_black.push_back(y/4);
							level1_neighborinfo_black.push_back(z/4);
							level1_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2; cur_z < z*2+2;cur_z++)
							{
								for(int cur_y = y*2; cur_y < y*2+2;cur_y++)
								{
									if(leafLevel0[cur_z*height*width+cur_y*width+(x+1)*2])
									{
										level1_neighborinfo_black.push_back(0);
										level1_neighborinfo_black.push_back((x+1)*2);
										level1_neighborinfo_black.push_back(cur_y);
										level1_neighborinfo_black.push_back(cur_z);
										level1_black_offset ++;
									}
								}
							}
						}
					}

					//y-
					if(y == 0)
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4);
						level1_neighborinfo_black.push_back(y/4-1);
						level1_neighborinfo_black.push_back(z/4);
						level1_black_offset ++;
					}
					else 
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+(y-1)*levelWidth1+x])
						{
							level1_neighborinfo_black.push_back(1);
							level1_neighborinfo_black.push_back(x);
							level1_neighborinfo_black.push_back(y-1);
							level1_neighborinfo_black.push_back(z);
							level1_black_offset ++;
						}
						else if(y%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y/2-1)*levelWidth2+x/2])
						{
							level1_neighborinfo_black.push_back(2);
							level1_neighborinfo_black.push_back(x/2);
							level1_neighborinfo_black.push_back(y/2-1);
							level1_neighborinfo_black.push_back(z/2);
							level1_black_offset ++;
						}
						else if(y%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y/4-1)*levelWidth3+x/4])
						{
							level1_neighborinfo_black.push_back(3);
							level1_neighborinfo_black.push_back(x/4);
							level1_neighborinfo_black.push_back(y/4-1);
							level1_neighborinfo_black.push_back(z/4);
							level1_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2; cur_z < z*2+2;cur_z++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[cur_z*height*width+(y*2-1)*width+cur_x])
									{
										level1_neighborinfo_black.push_back(0);
										level1_neighborinfo_black.push_back(cur_x);
										level1_neighborinfo_black.push_back(y*2-1);
										level1_neighborinfo_black.push_back(cur_z);
										level1_black_offset ++;
									}
								}
							}
						}
					}

					//y+
					if(y == levelHeight1-1)
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4);
						level1_neighborinfo_black.push_back((y+1)/4);
						level1_neighborinfo_black.push_back(z/4);
						level1_black_offset ++;
					}
					else
					{
						if(leafLevel1[z*levelHeight1*levelWidth1+(y+1)*levelWidth1+x])
						{
							level1_neighborinfo_black.push_back(1);
							level1_neighborinfo_black.push_back(x);
							level1_neighborinfo_black.push_back(y+1);
							level1_neighborinfo_black.push_back(z);
							level1_black_offset ++;
						}
						else if((y+1)%2 == 0 && leafLevel2[(z/2)*levelHeight2*levelWidth2+(y+1)/2*levelWidth2+x/2])
						{
							level1_neighborinfo_black.push_back(2);
							level1_neighborinfo_black.push_back(x/2);
							level1_neighborinfo_black.push_back((y+1)/2);
							level1_neighborinfo_black.push_back(z/2);
							level1_black_offset ++;
						}
						else if((y+1)%4 == 0 && leafLevel3[(z/4)*levelHeight3*levelWidth3+(y+1)/4*levelWidth3+x/4])
						{
							level1_neighborinfo_black.push_back(3);
							level1_neighborinfo_black.push_back(x/4);
							level1_neighborinfo_black.push_back((y+1)/4);
							level1_neighborinfo_black.push_back(z/4);
							level1_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z = z*2;cur_z < z*2+2;cur_z++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[cur_z*height*width+(y+1)*2*width+cur_x])
									{
										level1_neighborinfo_black.push_back(0);
										level1_neighborinfo_black.push_back(cur_x);
										level1_neighborinfo_black.push_back((y+1)*2);
										level1_neighborinfo_black.push_back(cur_z);
										level1_black_offset ++;
									}
								}
							}
						}
					}

					//z-
					if(z == 0)
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4);
						level1_neighborinfo_black.push_back(y/4);
						level1_neighborinfo_black.push_back(z/4-1);
						level1_black_offset ++;
					}
					else 
					{
						if(leafLevel1[(z-1)*levelHeight1*levelWidth1+y*levelWidth1+x])
						{
							level1_neighborinfo_black.push_back(1);
							level1_neighborinfo_black.push_back(x);
							level1_neighborinfo_black.push_back(y);
							level1_neighborinfo_black.push_back(z-1);
							level1_black_offset ++;
						}
						else if(z%2 == 0 && leafLevel2[(z/2-1)*levelHeight2*levelWidth2+(y/2)*levelWidth2+x/2])
						{
							level1_neighborinfo_black.push_back(2);
							level1_neighborinfo_black.push_back(x/2);
							level1_neighborinfo_black.push_back(y/2);
							level1_neighborinfo_black.push_back(z/2-1);
							level1_black_offset ++;
						}
						else if(z%4 == 0 && leafLevel3[(z/4-1)*levelHeight3*levelWidth3+(y/4)*levelWidth3+x/4])
						{
							level1_neighborinfo_black.push_back(3);
							level1_neighborinfo_black.push_back(x/4);
							level1_neighborinfo_black.push_back(y/4);
							level1_neighborinfo_black.push_back(z/4-1);
							level1_black_offset ++;
						}
						else //smaller
						{
							for(int cur_y = y*2; cur_y < y*2+2;cur_y++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[(z*2-1)*height*width+cur_y*width+cur_x])
									{
										level1_neighborinfo_black.push_back(0);
										level1_neighborinfo_black.push_back(cur_x);
										level1_neighborinfo_black.push_back(cur_y);
										level1_neighborinfo_black.push_back(z*2-1);
										level1_black_offset ++;
									}
								}
							}
						}
					}

					//z+
					if(y == levelDepth1-1)
					{
						level1_neighborinfo_black.push_back(3);
						level1_neighborinfo_black.push_back(x/4);
						level1_neighborinfo_black.push_back(y/4);
						level1_neighborinfo_black.push_back((z+1)/4);
						level1_black_offset ++;
					}
					else
					{
						if(leafLevel1[(z+1)*levelHeight1*levelWidth1+y*levelWidth1+x])
						{
							level1_neighborinfo_black.push_back(1);
							level1_neighborinfo_black.push_back(x);
							level1_neighborinfo_black.push_back(y);
							level1_neighborinfo_black.push_back(z+1);
							level1_black_offset ++;
						}
						else if((z+1)%2 == 0 && leafLevel2[(z+1)/2*levelHeight2*levelWidth2+(y/2)*levelWidth2+x/2])
						{
							level1_neighborinfo_black.push_back(2);
							level1_neighborinfo_black.push_back(x/2);
							level1_neighborinfo_black.push_back(y/2);
							level1_neighborinfo_black.push_back((z+1)/2);
							level1_black_offset ++;
						}
						else if((z+1)%4 == 0 && leafLevel3[(z+1)/4*levelHeight3*levelWidth3+(y)/4*levelWidth3+x/4])
						{
							level1_neighborinfo_black.push_back(3);
							level1_neighborinfo_black.push_back(x/4);
							level1_neighborinfo_black.push_back(y/4);
							level1_neighborinfo_black.push_back((z+1)/4);
							level1_black_offset ++;
						}
						else //smaller
						{
							for(int cur_y = y*2;cur_y < y*2+2;cur_y++)
							{
								for(int cur_x = x*2;cur_x < x*2+2;cur_x++)
								{
									if(leafLevel0[(z+1)*2*height*width+cur_y*width+cur_x])
									{
										level1_neighborinfo_black.push_back(0);
										level1_neighborinfo_black.push_back(cur_x);
										level1_neighborinfo_black.push_back(cur_y);
										level1_neighborinfo_black.push_back((z+1)*2);
										level1_black_offset ++;
									}
								}
							}
						}
					}

					level1_index_black.push_back(level1_black_offset-cur_offset);
					level1_index_black.push_back(cur_offset);
				}
			}
		}
	}
	

	/************************** for level2 ***************************/
	int level2_red_offset = 0;
	int level2_black_offset = 0;
	for(int z = 0;z < levelDepth2;z++)
	{
		for(int y = 0;y < levelHeight2;y++)
		{
			for(int x = 0;x < levelWidth2;x++)
			{
				if(!leafLevel2[z*levelHeight2*levelWidth2+y*levelWidth2+x])
					continue;

				if((x+y+z)%2 == 0)
				{
					level2_index_red.push_back(x);
					level2_index_red.push_back(y);
					level2_index_red.push_back(z);
					int cur_offset = level2_red_offset;

					//x-
					if(x == 0)
					{
						level2_neighborinfo_red.push_back(3);
						level2_neighborinfo_red.push_back(x/2-1);
						level2_neighborinfo_red.push_back(y/2);
						level2_neighborinfo_red.push_back(z/2);
						level2_red_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+y*levelWidth2+x-1])
						{
							level2_neighborinfo_red.push_back(2);
							level2_neighborinfo_red.push_back(x-1);
							level2_neighborinfo_red.push_back(y);
							level2_neighborinfo_red.push_back(z);
							level2_red_offset ++;
						}
						else if(x%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y/2)*levelWidth3+x/2-1])
						{
							level2_neighborinfo_red.push_back(3);
							level2_neighborinfo_red.push_back(x/2-1);
							level2_neighborinfo_red.push_back(y/2);
							level2_neighborinfo_red.push_back(z/2);
							level2_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+cur_y1*levelWidth1+x*2-1])
									{
										level2_neighborinfo_red.push_back(1);
										level2_neighborinfo_red.push_back(x*2-1);
										level2_neighborinfo_red.push_back(cur_y1);
										level2_neighborinfo_red.push_back(cur_z1);
										level2_red_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel0[cur_z2*height*width+cur_y2*width+x*4-1])
												{
													level2_neighborinfo_red.push_back(0);
													level2_neighborinfo_red.push_back(x*4-1);
													level2_neighborinfo_red.push_back(cur_y2);
													level2_neighborinfo_red.push_back(cur_z2);
													level2_red_offset ++;
												}
											}
										}
									}
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
						level2_neighborinfo_red.push_back(z/2);
						level2_red_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+y*levelWidth2+x+1])
						{
							level2_neighborinfo_red.push_back(2);
							level2_neighborinfo_red.push_back(x+1);
							level2_neighborinfo_red.push_back(y);
							level2_neighborinfo_red.push_back(z);
							level2_red_offset ++;
						}
						else if((x+1)%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y/2)*levelWidth3+(x+1)/2])
						{
							level2_neighborinfo_red.push_back(3);
							level2_neighborinfo_red.push_back((x+1)/2);
							level2_neighborinfo_red.push_back(y/2);
							level2_neighborinfo_red.push_back(z/2);
							level2_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+cur_y1*levelWidth1+(x+1)*2])
									{
										level2_neighborinfo_red.push_back(1);
										level2_neighborinfo_red.push_back((x+1)*2);
										level2_neighborinfo_red.push_back(cur_y1);
										level2_neighborinfo_red.push_back(cur_z1);
										level2_red_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel0[cur_z2*height*width+cur_y2*width+(x+1)*4])
												{
													level2_neighborinfo_red.push_back(0);
													level2_neighborinfo_red.push_back((x+1)*4);
													level2_neighborinfo_red.push_back(cur_y2);
													level2_neighborinfo_red.push_back(cur_z2);
													level2_red_offset ++;
												}
											}
										}
									}
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
						level2_neighborinfo_red.push_back(z/2);
						level2_red_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+(y-1)*levelWidth2+x])
						{
							level2_neighborinfo_red.push_back(2);
							level2_neighborinfo_red.push_back(x);
							level2_neighborinfo_red.push_back(y-1);
							level2_neighborinfo_red.push_back(z);
							level2_red_offset ++;
						}
						else if(y%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y/2-1)*levelWidth3+x/2])
						{
							level2_neighborinfo_red.push_back(3);
							level2_neighborinfo_red.push_back(x/2);
							level2_neighborinfo_red.push_back(y/2-1);
							level2_neighborinfo_red.push_back(z/2);
							level2_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+(y*2-1)*levelWidth1+cur_x1])
									{
										level2_neighborinfo_red.push_back(1);
										level2_neighborinfo_red.push_back(cur_x1);
										level2_neighborinfo_red.push_back(y*2-1);
										level2_neighborinfo_red.push_back(cur_z1);
										level2_red_offset ++;
									}
									else 
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												level2_neighborinfo_red.push_back(0);
												level2_neighborinfo_red.push_back(cur_x2);
												level2_neighborinfo_red.push_back(y*4-1);
												level2_neighborinfo_red.push_back(cur_z2);
												level2_red_offset ++;
											}
										}
									}
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
						level2_neighborinfo_red.push_back(z/2);
						level2_red_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+(y+1)*levelWidth2+x])
						{
							level2_neighborinfo_red.push_back(2);
							level2_neighborinfo_red.push_back(x);
							level2_neighborinfo_red.push_back(y+1);
							level2_neighborinfo_red.push_back(z);
							level2_red_offset ++;
						}
						else if((y+1)%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y+1)/2*levelWidth3+x/2])
						{
							level2_neighborinfo_red.push_back(3);
							level2_neighborinfo_red.push_back(x/2);
							level2_neighborinfo_red.push_back((y+1)/2);
							level2_neighborinfo_red.push_back(z/2);
							level2_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+(y+1)*2*levelWidth1+cur_x1])
									{
										level2_neighborinfo_red.push_back(1);
										level2_neighborinfo_red.push_back(cur_x1);
										level2_neighborinfo_red.push_back((y+1)*2);
										level2_neighborinfo_red.push_back(cur_z1);
										level2_red_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel0[cur_z2*height*width+(y+1)*4*width+cur_x2])
												{
													level2_neighborinfo_red.push_back(0);
													level2_neighborinfo_red.push_back(cur_x2);
													level2_neighborinfo_red.push_back((y+1)*4);
													level2_neighborinfo_red.push_back(cur_z2);
													level2_red_offset ++;
												}
											}
										}
									}
								}
							}
						}
					}

					//z-
					if(z == 0)
					{
						level2_neighborinfo_red.push_back(3);
						level2_neighborinfo_red.push_back(x/2);
						level2_neighborinfo_red.push_back(y/2);
						level2_neighborinfo_red.push_back(z/2-1);
						level2_red_offset ++;
					}
					else
					{
						if(leafLevel2[(z-1)*levelHeight2*levelWidth2+y*levelWidth2+x])
						{
							level2_neighborinfo_red.push_back(2);
							level2_neighborinfo_red.push_back(x);
							level2_neighborinfo_red.push_back(y);
							level2_neighborinfo_red.push_back(z-1);
							level2_red_offset ++;
						}
						else if(z%2 == 0 && leafLevel3[(z/2-1)*levelHeight3*levelWidth3+(y/2)*levelWidth3+x/2])
						{
							level2_neighborinfo_red.push_back(3);
							level2_neighborinfo_red.push_back(x/2);
							level2_neighborinfo_red.push_back(y/2);
							level2_neighborinfo_red.push_back(z/2-1);
							level2_red_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[(z*2-1)*levelHeight1*levelWidth1+cur_y1*levelWidth1+cur_x1])
									{
										level2_neighborinfo_red.push_back(1);
										level2_neighborinfo_red.push_back(cur_x1);
										level2_neighborinfo_red.push_back(cur_y1);
										level2_neighborinfo_red.push_back(z*2-1);
										level2_red_offset ++;
									}
									else 
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												level2_neighborinfo_red.push_back(0);
												level2_neighborinfo_red.push_back(cur_x2);
												level2_neighborinfo_red.push_back(cur_y2);
												level2_neighborinfo_red.push_back(z*4-1);
												level2_red_offset ++;
											}
										}
									}
								}
							}
						}
					}

					//z+
					if(z == levelDepth2-1)
					{
						level2_neighborinfo_red.push_back(3);
						level2_neighborinfo_red.push_back(x/2);
						level2_neighborinfo_red.push_back(y/2);
						level2_neighborinfo_red.push_back((z+1)/2);
						level2_red_offset ++;
					}
					else
					{
						if(leafLevel2[(z+1)*levelHeight2*levelWidth2+y*levelWidth2+x])
						{
							level2_neighborinfo_red.push_back(2);
							level2_neighborinfo_red.push_back(x);
							level2_neighborinfo_red.push_back(y);
							level2_neighborinfo_red.push_back(z+1);
							level2_red_offset ++;
						}
						else if((z+1)%2 == 0 && leafLevel3[(z+1)/2*levelHeight3*levelWidth3+(y/2)*levelWidth3+x/2])
						{
							level2_neighborinfo_red.push_back(3);
							level2_neighborinfo_red.push_back(x/2);
							level2_neighborinfo_red.push_back(y/2);
							level2_neighborinfo_red.push_back((z+1)/2);
							level2_red_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[(z+1)*2*levelHeight1*levelWidth1+cur_y1*levelWidth1+cur_x1])
									{
										level2_neighborinfo_red.push_back(1);
										level2_neighborinfo_red.push_back(cur_x1);
										level2_neighborinfo_red.push_back(cur_y1);
										level2_neighborinfo_red.push_back((z+1)*2);
										level2_red_offset ++;
									}
									else
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel0[(z+1)*4*height*width+cur_y2*width+cur_x2])
												{
													level2_neighborinfo_red.push_back(0);
													level2_neighborinfo_red.push_back(cur_x2);
													level2_neighborinfo_red.push_back(cur_y2);
													level2_neighborinfo_red.push_back((z+1)*4);
													level2_red_offset ++;
												}
											}
										}
									}
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
					level2_index_black.push_back(z);
					int cur_offset = level2_black_offset;

					//x-
					if(x == 0)
					{
						level2_neighborinfo_black.push_back(3);
						level2_neighborinfo_black.push_back(x/2-1);
						level2_neighborinfo_black.push_back(y/2);
						level2_neighborinfo_black.push_back(z/2);
						level2_black_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+y*levelWidth2+x-1])
						{
							level2_neighborinfo_black.push_back(2);
							level2_neighborinfo_black.push_back(x-1);
							level2_neighborinfo_black.push_back(y);
							level2_neighborinfo_black.push_back(z);
							level2_black_offset ++;
						}
						else if(x%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y/2)*levelWidth3+x/2-1])
						{
							level2_neighborinfo_black.push_back(3);
							level2_neighborinfo_black.push_back(x/2-1);
							level2_neighborinfo_black.push_back(y/2);
							level2_neighborinfo_black.push_back(z/2);
							level2_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+cur_y1*levelWidth1+x*2-1])
									{
										level2_neighborinfo_black.push_back(1);
										level2_neighborinfo_black.push_back(x*2-1);
										level2_neighborinfo_black.push_back(cur_y1);
										level2_neighborinfo_black.push_back(cur_z1);
										level2_black_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel0[cur_z2*height*width+cur_y2*width+x*4-1])
												{
													level2_neighborinfo_black.push_back(0);
													level2_neighborinfo_black.push_back(x*4-1);
													level2_neighborinfo_black.push_back(cur_y2);
													level2_neighborinfo_black.push_back(cur_z2);
													level2_black_offset ++;
												}
											}
										}
									}
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
						level2_neighborinfo_black.push_back(z/2);
						level2_black_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+y*levelWidth2+x+1])
						{
							level2_neighborinfo_black.push_back(2);
							level2_neighborinfo_black.push_back(x+1);
							level2_neighborinfo_black.push_back(y);
							level2_neighborinfo_black.push_back(z);
							level2_black_offset ++;
						}
						else if((x+1)%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y/2)*levelWidth3+(x+1)/2])
						{
							level2_neighborinfo_black.push_back(3);
							level2_neighborinfo_black.push_back((x+1)/2);
							level2_neighborinfo_black.push_back(y/2);
							level2_neighborinfo_black.push_back(z/2);
							level2_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+cur_y1*levelWidth1+(x+1)*2])
									{
										level2_neighborinfo_black.push_back(1);
										level2_neighborinfo_black.push_back((x+1)*2);
										level2_neighborinfo_black.push_back(cur_y1);
										level2_neighborinfo_black.push_back(cur_z1);
										level2_black_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel0[cur_z2*height*width+cur_y2*width+(x+1)*4])
												{
													level2_neighborinfo_black.push_back(0);
													level2_neighborinfo_black.push_back((x+1)*4);
													level2_neighborinfo_black.push_back(cur_y2);
													level2_neighborinfo_black.push_back(cur_z2);
													level2_black_offset ++;
												}
											}
										}
									}
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
						level2_neighborinfo_black.push_back(z/2);
						level2_black_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+(y-1)*levelWidth2+x])
						{
							level2_neighborinfo_black.push_back(2);
							level2_neighborinfo_black.push_back(x);
							level2_neighborinfo_black.push_back(y-1);
							level2_neighborinfo_black.push_back(z);
							level2_black_offset ++;
						}
						else if(y%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y/2-1)*levelWidth3+x/2])
						{
							level2_neighborinfo_black.push_back(3);
							level2_neighborinfo_black.push_back(x/2);
							level2_neighborinfo_black.push_back(y/2-1);
							level2_neighborinfo_black.push_back(z/2);
							level2_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+(y*2-1)*levelWidth1+cur_x1])
									{
										level2_neighborinfo_black.push_back(1);
										level2_neighborinfo_black.push_back(cur_x1);
										level2_neighborinfo_black.push_back(y*2-1);
										level2_neighborinfo_black.push_back(cur_z1);
										level2_black_offset ++;
									}
									else 
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												level2_neighborinfo_black.push_back(0);
												level2_neighborinfo_black.push_back(cur_x2);
												level2_neighborinfo_black.push_back(y*4-1);
												level2_neighborinfo_black.push_back(cur_z2);
												level2_black_offset ++;
											}
										}
									}
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
						level2_neighborinfo_black.push_back(z/2);
						level2_black_offset ++;
					}
					else
					{
						if(leafLevel2[z*levelHeight2*levelWidth2+(y+1)*levelWidth2+x])
						{
							level2_neighborinfo_black.push_back(2);
							level2_neighborinfo_black.push_back(x);
							level2_neighborinfo_black.push_back(y+1);
							level2_neighborinfo_black.push_back(z);
							level2_black_offset ++;
						}
						else if((y+1)%2 == 0 && leafLevel3[(z/2)*levelHeight3*levelWidth3+(y+1)/2*levelWidth3+x/2])
						{
							level2_neighborinfo_black.push_back(3);
							level2_neighborinfo_black.push_back(x/2);
							level2_neighborinfo_black.push_back((y+1)/2);
							level2_neighborinfo_black.push_back(z/2);
							level2_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[cur_z1*levelHeight1*levelWidth1+(y+1)*2*levelWidth1+cur_x1])
									{
										level2_neighborinfo_black.push_back(1);
										level2_neighborinfo_black.push_back(cur_x1);
										level2_neighborinfo_black.push_back((y+1)*2);
										level2_neighborinfo_black.push_back(cur_z1);
										level2_black_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel0[cur_z2*height*width+(y+1)*4*width+cur_x2])
												{
													level2_neighborinfo_black.push_back(0);
													level2_neighborinfo_black.push_back(cur_x2);
													level2_neighborinfo_black.push_back((y+1)*4);
													level2_neighborinfo_black.push_back(cur_z2);
													level2_black_offset ++;
												}
											}
										}
									}
								}
							}
						}
					}

					//z-
					if(z == 0)
					{
						level2_neighborinfo_black.push_back(3);
						level2_neighborinfo_black.push_back(x/2);
						level2_neighborinfo_black.push_back(y/2);
						level2_neighborinfo_black.push_back(z/2-1);
						level2_black_offset ++;
					}
					else
					{
						if(leafLevel2[(z-1)*levelHeight2*levelWidth2+y*levelWidth2+x])
						{
							level2_neighborinfo_black.push_back(2);
							level2_neighborinfo_black.push_back(x);
							level2_neighborinfo_black.push_back(y);
							level2_neighborinfo_black.push_back(z-1);
							level2_black_offset ++;
						}
						else if(z%2 == 0 && leafLevel3[(z/2-1)*levelHeight3*levelWidth3+(y/2)*levelWidth3+x/2])
						{
							level2_neighborinfo_black.push_back(3);
							level2_neighborinfo_black.push_back(x/2);
							level2_neighborinfo_black.push_back(y/2);
							level2_neighborinfo_black.push_back(z/2-1);
							level2_black_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[(z*2-1)*levelHeight1*levelWidth1+cur_y1*levelWidth1+cur_x1])
									{
										level2_neighborinfo_black.push_back(1);
										level2_neighborinfo_black.push_back(cur_x1);
										level2_neighborinfo_black.push_back(cur_y1);
										level2_neighborinfo_black.push_back(z*2-1);
										level2_black_offset ++;
									}
									else 
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												level2_neighborinfo_black.push_back(0);
												level2_neighborinfo_black.push_back(cur_x2);
												level2_neighborinfo_black.push_back(cur_y2);
												level2_neighborinfo_black.push_back(z*4-1);
												level2_black_offset ++;
											}
										}
									}
								}
							}
						}
					}

					//z+
					if(z == levelDepth2-1)
					{
						level2_neighborinfo_black.push_back(3);
						level2_neighborinfo_black.push_back(x/2);
						level2_neighborinfo_black.push_back(y/2);
						level2_neighborinfo_black.push_back((z+1)/2);
						level2_black_offset ++;
					}
					else
					{
						if(leafLevel2[(z+1)*levelHeight2*levelWidth2+y*levelWidth2+x])
						{
							level2_neighborinfo_black.push_back(2);
							level2_neighborinfo_black.push_back(x);
							level2_neighborinfo_black.push_back(y);
							level2_neighborinfo_black.push_back(z+1);
							level2_black_offset ++;
						}
						else if((z+1)%2 == 0 && leafLevel3[(z+1)/2*levelHeight3*levelWidth3+(y/2)*levelWidth3+x/2])
						{
							level2_neighborinfo_black.push_back(3);
							level2_neighborinfo_black.push_back(x/2);
							level2_neighborinfo_black.push_back(y/2);
							level2_neighborinfo_black.push_back((z+1)/2);
							level2_black_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel1[(z+1)*2*levelHeight1*levelWidth1+cur_y1*levelWidth1+cur_x1])
									{
										level2_neighborinfo_black.push_back(1);
										level2_neighborinfo_black.push_back(cur_x1);
										level2_neighborinfo_black.push_back(cur_y1);
										level2_neighborinfo_black.push_back((z+1)*2);
										level2_black_offset ++;
									}
									else
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel0[(z+1)*4*height*width+cur_y2*width+cur_x2])
												{
													level2_neighborinfo_black.push_back(0);
													level2_neighborinfo_black.push_back(cur_x2);
													level2_neighborinfo_black.push_back(cur_y2);
													level2_neighborinfo_black.push_back((z+1)*4);
													level2_black_offset ++;
												}
											}
										}
									}
								}
							}
						}
					}
					level2_index_black.push_back(level2_black_offset-cur_offset);
					level2_index_black.push_back(cur_offset);
				}
			}
		}
	}
	

	/********************** for level3   ****************************/
	int level3_red_offset = 0;
	int level3_black_offset = 0;
	for(int z = 0;z < levelDepth3;z++)
	{
		for(int y = 0;y < levelHeight3;y++)
		{
			for(int x = 0;x < levelWidth3;x++)
			{
				if(!leafLevel3[y*levelWidth3+x])
					continue;

				if((x+y+z)%2 == 0)
				{
					level3_index_red.push_back(x);
					level3_index_red.push_back(y);
					level3_index_red.push_back(z);
					int cur_offset = level3_red_offset;

					//x-
					if(x == 0)
					{
						level3_neighborinfo_red.push_back(3);
						level3_neighborinfo_red.push_back(x-1);
						level3_neighborinfo_red.push_back(y);
						level3_neighborinfo_red.push_back(z);
						level3_red_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+y*levelWidth3+x-1])
						{
							level3_neighborinfo_red.push_back(3);
							level3_neighborinfo_red.push_back(x-1);
							level3_neighborinfo_red.push_back(y);
							level3_neighborinfo_red.push_back(z);
							level3_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+cur_y1*levelWidth2+x*2-1])
									{
										level3_neighborinfo_red.push_back(2);
										level3_neighborinfo_red.push_back(x*2-1);
										level3_neighborinfo_red.push_back(cur_y1);
										level3_neighborinfo_red.push_back(cur_z1);
										level3_red_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+cur_y2*levelWidth1+x*4-1])
												{
													level3_neighborinfo_red.push_back(1);
													level3_neighborinfo_red.push_back(x*4-1);
													level3_neighborinfo_red.push_back(cur_y2);
													level3_neighborinfo_red.push_back(cur_z2);
													level3_red_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
														{
															if(leafLevel0[cur_z3*height*width+cur_y3*width+x*8-1])
															{
																level3_neighborinfo_red.push_back(0);
																level3_neighborinfo_red.push_back(x*8-1);
																level3_neighborinfo_red.push_back(cur_y3);
																level3_neighborinfo_red.push_back(cur_z3);
																level3_red_offset ++;
															}
														}
													}
												}
											}
										}
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
						level3_neighborinfo_red.push_back(z);
						level3_red_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+y*levelWidth3+x+1])
						{
							level3_neighborinfo_red.push_back(3);
							level3_neighborinfo_red.push_back(x+1);
							level3_neighborinfo_red.push_back(y);
							level3_neighborinfo_red.push_back(z);
							level3_red_offset ++;
						}
						else //smaller 
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+cur_y1*levelWidth2+(x+1)*2])
									{
										level3_neighborinfo_red.push_back(2);
										level3_neighborinfo_red.push_back((x+1)*2);
										level3_neighborinfo_red.push_back(cur_y1);
										level3_neighborinfo_red.push_back(cur_z1);
										level3_red_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+cur_y2*levelWidth1+(x+1)*4])
												{
													level3_neighborinfo_red.push_back(1);
													level3_neighborinfo_red.push_back((x+1)*4);
													level3_neighborinfo_red.push_back(cur_y2);
													level3_neighborinfo_red.push_back(cur_z2);
													level3_red_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
														{
															if(leafLevel0[cur_z3*height*width+cur_y3*width+(x+1)*8])
															{
																level3_neighborinfo_red.push_back(0);
																level3_neighborinfo_red.push_back((x+1)*8);
																level3_neighborinfo_red.push_back(cur_y3);
																level3_neighborinfo_red.push_back(cur_z3);
																level3_red_offset ++;
															}
														}
													}
												}
											}
										}
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
						level3_neighborinfo_red.push_back(z);
						level3_red_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+(y-1)*levelWidth3+x])
						{
							level3_neighborinfo_red.push_back(3);
							level3_neighborinfo_red.push_back(x);
							level3_neighborinfo_red.push_back(y-1);
							level3_neighborinfo_red.push_back(z);
							level3_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+(y*2-1)*levelWidth2+cur_x1])
									{
										level3_neighborinfo_red.push_back(2);
										level3_neighborinfo_red.push_back(cur_x1);
										level3_neighborinfo_red.push_back(y*2-1);
										level3_neighborinfo_red.push_back(cur_z1);
										level3_red_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+(y*4-1)*levelWidth1+cur_x2])
												{
													level3_neighborinfo_red.push_back(1);
													level3_neighborinfo_red.push_back(cur_x2);
													level3_neighborinfo_red.push_back(y*4-1);
													level3_neighborinfo_red.push_back(cur_z2);
													level3_red_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[cur_z3*height*width+(y*8-1)*width+cur_x3])
															{
																level3_neighborinfo_red.push_back(0);
																level3_neighborinfo_red.push_back(cur_x3);
																level3_neighborinfo_red.push_back(y*8-1);
																level3_neighborinfo_red.push_back(cur_z3);
																level3_red_offset ++;
															}
														}
													}
												}
											}
										}
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
						level3_neighborinfo_red.push_back(z);
						level3_red_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+(y+1)*levelWidth3+x])
						{
							level3_neighborinfo_red.push_back(3);
							level3_neighborinfo_red.push_back(x);
							level3_neighborinfo_red.push_back(y+1);
							level3_neighborinfo_red.push_back(z);
							level3_red_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+(y+1)*2*levelWidth2+cur_x1])
									{
										level3_neighborinfo_red.push_back(2);
										level3_neighborinfo_red.push_back(cur_x1);
										level3_neighborinfo_red.push_back((y+1)*2);
										level3_neighborinfo_red.push_back(cur_z1);
										level3_red_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+(y+1)*4*levelWidth1+cur_z2])
												{
													level3_neighborinfo_red.push_back(1);
													level3_neighborinfo_red.push_back(cur_x2);
													level3_neighborinfo_red.push_back((y+1)*4);
													level3_neighborinfo_red.push_back(cur_z2);
													level3_red_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[cur_z3*height*width+(y+1)*8*width+cur_x3])
															{
																level3_neighborinfo_red.push_back(0);
																level3_neighborinfo_red.push_back(cur_x3);
																level3_neighborinfo_red.push_back((y+1)*8);
																level3_neighborinfo_red.push_back(cur_z3);
																level3_red_offset ++;
															}
														}
													}
												}
											}
										}
									}
								}
							}
							
						}
					}

					//z-
					if(z == 0)
					{
						level3_neighborinfo_red.push_back(3);
						level3_neighborinfo_red.push_back(x);
						level3_neighborinfo_red.push_back(y);
						level3_neighborinfo_red.push_back(z-1);
						level3_red_offset ++;
					}
					else
					{
						if(leafLevel3[(z-1)*levelHeight3*levelWidth3+y*levelWidth3+x])
						{
							level3_neighborinfo_red.push_back(3);
							level3_neighborinfo_red.push_back(x);
							level3_neighborinfo_red.push_back(y);
							level3_neighborinfo_red.push_back(z-1);
							level3_red_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[(z*2-1)*levelHeight2*levelWidth2+cur_y1*levelWidth2+cur_x1])
									{
										level3_neighborinfo_red.push_back(2);
										level3_neighborinfo_red.push_back(cur_x1);
										level3_neighborinfo_red.push_back(cur_y1);
										level3_neighborinfo_red.push_back(z*2-1);
										level3_red_offset ++;
									}
									else
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[(z*4-1)*levelHeight1*levelWidth1+cur_y2*levelWidth1+cur_x2])
												{
													level3_neighborinfo_red.push_back(1);
													level3_neighborinfo_red.push_back(cur_x2);
													level3_neighborinfo_red.push_back(cur_y2);
													level3_neighborinfo_red.push_back(z*4-1);
													level3_red_offset ++;
												}
												else
												{
													for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[(z*8-1)*height*width+cur_y3*width+cur_x3])
															{
																level3_neighborinfo_red.push_back(0);
																level3_neighborinfo_red.push_back(cur_x3);
																level3_neighborinfo_red.push_back(cur_y3);
																level3_neighborinfo_red.push_back(z*8-1);
																level3_red_offset ++;
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}

					//z+
					if(z == levelDepth3-1)
					{
						level3_neighborinfo_red.push_back(3);
						level3_neighborinfo_red.push_back(x);
						level3_neighborinfo_red.push_back(y);
						level3_neighborinfo_red.push_back(z+1);
						level3_red_offset ++;
					}
					else
					{
						if(leafLevel3[(z+1)*levelHeight3*levelWidth3+y*levelWidth3+x])
						{
							level3_neighborinfo_red.push_back(3);
							level3_neighborinfo_red.push_back(x);
							level3_neighborinfo_red.push_back(y);
							level3_neighborinfo_red.push_back(z+1);
							level3_red_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[(z+1)*2*levelHeight2*levelWidth2+cur_y1*levelWidth2+cur_x1])
									{
										level3_neighborinfo_red.push_back(2);
										level3_neighborinfo_red.push_back(cur_x1);
										level3_neighborinfo_red.push_back(cur_y1);
										level3_neighborinfo_red.push_back((z+1)*2);
										level3_red_offset ++;
									}
									else
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[(z+1)*4*levelHeight1*levelWidth1+cur_y2*levelWidth1+cur_x2])
												{
													level3_neighborinfo_red.push_back(1);
													level3_neighborinfo_red.push_back(cur_x2);
													level3_neighborinfo_red.push_back(cur_y2);
													level3_neighborinfo_red.push_back((z+1)*4);
													level3_red_offset ++;
												}
												else
												{
													for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[(z+1)*8*height*width+cur_y3*width+cur_x3])
															{
																level3_neighborinfo_red.push_back(0);
																level3_neighborinfo_red.push_back(cur_x3);
																level3_neighborinfo_red.push_back(cur_y3);
																level3_neighborinfo_red.push_back((z+1)*8);
																level3_red_offset ++;
															}
														}
													}
												}
											}
										}
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
					level3_index_black.push_back(z);
					int cur_offset = level3_black_offset;

					//x-
					if(x == 0)
					{
						level3_neighborinfo_black.push_back(3);
						level3_neighborinfo_black.push_back(x-1);
						level3_neighborinfo_black.push_back(y);
						level3_neighborinfo_black.push_back(z);
						level3_black_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+y*levelWidth3+x-1])
						{
							level3_neighborinfo_black.push_back(3);
							level3_neighborinfo_black.push_back(x-1);
							level3_neighborinfo_black.push_back(y);
							level3_neighborinfo_black.push_back(z);
							level3_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+cur_y1*levelWidth2+x*2-1])
									{
										level3_neighborinfo_black.push_back(2);
										level3_neighborinfo_black.push_back(x*2-1);
										level3_neighborinfo_black.push_back(cur_y1);
										level3_neighborinfo_black.push_back(cur_z1);
										level3_black_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+cur_y2*levelWidth1+x*4-1])
												{
													level3_neighborinfo_black.push_back(1);
													level3_neighborinfo_black.push_back(x*4-1);
													level3_neighborinfo_black.push_back(cur_y2);
													level3_neighborinfo_black.push_back(cur_z2);
													level3_black_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
														{
															if(leafLevel0[cur_z3*height*width+cur_y3*width+x*8-1])
															{
																level3_neighborinfo_black.push_back(0);
																level3_neighborinfo_black.push_back(x*8-1);
																level3_neighborinfo_black.push_back(cur_y3);
																level3_neighborinfo_black.push_back(cur_z3);
																level3_black_offset ++;
															}
														}
													}
												}
											}
										}
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
						level3_neighborinfo_black.push_back(z);
						level3_black_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+y*levelWidth3+x+1])
						{
							level3_neighborinfo_black.push_back(3);
							level3_neighborinfo_black.push_back(x+1);
							level3_neighborinfo_black.push_back(y);
							level3_neighborinfo_black.push_back(z);
							level3_black_offset ++;
						}
						else //smaller 
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+cur_y1*levelWidth2+(x+1)*2])
									{
										level3_neighborinfo_black.push_back(2);
										level3_neighborinfo_black.push_back((x+1)*2);
										level3_neighborinfo_black.push_back(cur_y1);
										level3_neighborinfo_black.push_back(cur_z1);
										level3_black_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+cur_y2*levelWidth1+(x+1)*4])
												{
													level3_neighborinfo_black.push_back(1);
													level3_neighborinfo_black.push_back((x+1)*4);
													level3_neighborinfo_black.push_back(cur_y2);
													level3_neighborinfo_black.push_back(cur_z2);
													level3_black_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
														{
															if(leafLevel0[cur_z3*height*width+cur_y3*width+(x+1)*8])
															{
																level3_neighborinfo_black.push_back(0);
																level3_neighborinfo_black.push_back((x+1)*8);
																level3_neighborinfo_black.push_back(cur_y3);
																level3_neighborinfo_black.push_back(cur_z3);
																level3_black_offset ++;
															}
														}
													}
												}
											}
										}
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
						level3_neighborinfo_black.push_back(z);
						level3_black_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+(y-1)*levelWidth3+x])
						{
							level3_neighborinfo_black.push_back(3);
							level3_neighborinfo_black.push_back(x);
							level3_neighborinfo_black.push_back(y-1);
							level3_neighborinfo_black.push_back(z);
							level3_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+(y*2-1)*levelWidth2+cur_x1])
									{
										level3_neighborinfo_black.push_back(2);
										level3_neighborinfo_black.push_back(cur_x1);
										level3_neighborinfo_black.push_back(y*2-1);
										level3_neighborinfo_black.push_back(cur_z1);
										level3_black_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+(y*4-1)*levelWidth1+cur_x2])
												{
													level3_neighborinfo_black.push_back(1);
													level3_neighborinfo_black.push_back(cur_x2);
													level3_neighborinfo_black.push_back(y*4-1);
													level3_neighborinfo_black.push_back(cur_z2);
													level3_black_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[cur_z3*height*width+(y*8-1)*width+cur_x3])
															{
																level3_neighborinfo_black.push_back(0);
																level3_neighborinfo_black.push_back(cur_x3);
																level3_neighborinfo_black.push_back(y*8-1);
																level3_neighborinfo_black.push_back(cur_z3);
																level3_black_offset ++;
															}
														}
													}
												}
											}
										}
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
						level3_neighborinfo_black.push_back(z);
						level3_black_offset ++;
					}
					else
					{
						if(leafLevel3[z*levelHeight3*levelWidth3+(y+1)*levelWidth3+x])
						{
							level3_neighborinfo_black.push_back(3);
							level3_neighborinfo_black.push_back(x);
							level3_neighborinfo_black.push_back(y+1);
							level3_neighborinfo_black.push_back(z);
							level3_black_offset ++;
						}
						else //smaller
						{
							for(int cur_z1 = z*2; cur_z1 < z*2+2; cur_z1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[cur_z1*levelHeight2*levelWidth2+(y+1)*2*levelWidth2+cur_x1])
									{
										level3_neighborinfo_black.push_back(2);
										level3_neighborinfo_black.push_back(cur_x1);
										level3_neighborinfo_black.push_back((y+1)*2);
										level3_neighborinfo_black.push_back(cur_z1);
										level3_black_offset ++;
									}
									else
									{
										for(int cur_z2 = cur_z1*2; cur_z2 < cur_z1*2+2; cur_z2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[cur_z2*levelHeight1*levelWidth1+(y+1)*4*levelWidth1+cur_z2])
												{
													level3_neighborinfo_black.push_back(1);
													level3_neighborinfo_black.push_back(cur_x2);
													level3_neighborinfo_black.push_back((y+1)*4);
													level3_neighborinfo_black.push_back(cur_z2);
													level3_black_offset ++;
												}
												else
												{
													for(int cur_z3 = cur_z2*2; cur_z3 < cur_z2*2+2; cur_z3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[cur_z3*height*width+(y+1)*8*width+cur_x3])
															{
																level3_neighborinfo_black.push_back(0);
																level3_neighborinfo_black.push_back(cur_x3);
																level3_neighborinfo_black.push_back((y+1)*8);
																level3_neighborinfo_black.push_back(cur_z3);
																level3_black_offset ++;
															}
														}
													}
												}
											}
										}
									}
								}
							}

						}
					}

					//z-
					if(z == 0)
					{
						level3_neighborinfo_black.push_back(3);
						level3_neighborinfo_black.push_back(x);
						level3_neighborinfo_black.push_back(y);
						level3_neighborinfo_black.push_back(z-1);
						level3_black_offset ++;
					}
					else
					{
						if(leafLevel3[(z-1)*levelHeight3*levelWidth3+y*levelWidth3+x])
						{
							level3_neighborinfo_black.push_back(3);
							level3_neighborinfo_black.push_back(x);
							level3_neighborinfo_black.push_back(y);
							level3_neighborinfo_black.push_back(z-1);
							level3_black_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[(z*2-1)*levelHeight2*levelWidth2+cur_y1*levelWidth2+cur_x1])
									{
										level3_neighborinfo_black.push_back(2);
										level3_neighborinfo_black.push_back(cur_x1);
										level3_neighborinfo_black.push_back(cur_y1);
										level3_neighborinfo_black.push_back(z*2-1);
										level3_black_offset ++;
									}
									else
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[(z*4-1)*levelHeight1*levelWidth1+cur_y2*levelWidth1+cur_x2])
												{
													level3_neighborinfo_black.push_back(1);
													level3_neighborinfo_black.push_back(cur_x2);
													level3_neighborinfo_black.push_back(cur_y2);
													level3_neighborinfo_black.push_back(z*4-1);
													level3_black_offset ++;
												}
												else
												{
													for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[(z*8-1)*height*width+cur_y3*width+cur_x3])
															{
																level3_neighborinfo_black.push_back(0);
																level3_neighborinfo_black.push_back(cur_x3);
																level3_neighborinfo_black.push_back(cur_y3);
																level3_neighborinfo_black.push_back(z*8-1);
																level3_black_offset ++;
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}

					//z+
					if(z == levelDepth3-1)
					{
						level3_neighborinfo_black.push_back(3);
						level3_neighborinfo_black.push_back(x);
						level3_neighborinfo_black.push_back(y);
						level3_neighborinfo_black.push_back(z+1);
						level3_black_offset ++;
					}
					else
					{
						if(leafLevel3[(z+1)*levelHeight3*levelWidth3+y*levelWidth3+x])
						{
							level3_neighborinfo_black.push_back(3);
							level3_neighborinfo_black.push_back(x);
							level3_neighborinfo_black.push_back(y);
							level3_neighborinfo_black.push_back(z+1);
							level3_black_offset ++;
						}
						else //smaller
						{
							for(int cur_y1 = y*2; cur_y1 < y*2+2; cur_y1++)
							{
								for(int cur_x1 = x*2; cur_x1 < x*2+2; cur_x1++)
								{
									if(leafLevel2[(z+1)*2*levelHeight2*levelWidth2+cur_y1*levelWidth2+cur_x1])
									{
										level3_neighborinfo_black.push_back(2);
										level3_neighborinfo_black.push_back(cur_x1);
										level3_neighborinfo_black.push_back(cur_y1);
										level3_neighborinfo_black.push_back((z+1)*2);
										level3_black_offset ++;
									}
									else
									{
										for(int cur_y2 = cur_y1*2; cur_y2 < cur_y1*2+2; cur_y2++)
										{
											for(int cur_x2 = cur_x1*2; cur_x2 < cur_x1*2+2; cur_x2++)
											{
												if(leafLevel1[(z+1)*4*levelHeight1*levelWidth1+cur_y2*levelWidth1+cur_x2])
												{
													level3_neighborinfo_black.push_back(1);
													level3_neighborinfo_black.push_back(cur_x2);
													level3_neighborinfo_black.push_back(cur_y2);
													level3_neighborinfo_black.push_back((z+1)*4);
													level3_black_offset ++;
												}
												else
												{
													for(int cur_y3 = cur_y2*2; cur_y3 < cur_y2*2+2; cur_y3++)
													{
														for(int cur_x3 = cur_x2*2; cur_x3 < cur_x2*2+2; cur_x3++)
														{
															if(leafLevel0[(z+1)*8*height*width+cur_y3*width+cur_x3])
															{
																level3_neighborinfo_black.push_back(0);
																level3_neighborinfo_black.push_back(cur_x3);
																level3_neighborinfo_black.push_back(cur_y3);
																level3_neighborinfo_black.push_back((z+1)*8);
																level3_black_offset ++;
															}
														}
													}
												}
											}
										}
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

	///bucket sort

	//_sort_neighborinfo(level0_index_red,level0_neighborinfo_red);		// level0 no need sort
	//_sort_neighborinfo(level0_index_black,level0_neighborinfo_black);	// level0 no need sort

	// WARNING: sorting seems to slowing down the CUDA kernels
	//_sort_neighborinfo(level1_index_red,level1_neighborinfo_red); 
	//_sort_neighborinfo(level1_index_black,level1_neighborinfo_black);
	//_sort_neighborinfo(level2_index_red,level2_neighborinfo_red);
	//_sort_neighborinfo(level2_index_black,level2_neighborinfo_black);
	//_sort_neighborinfo(level3_index_red,level3_neighborinfo_red);
	//_sort_neighborinfo(level3_index_black,level3_neighborinfo_black);
	
}

void ZQ_SmokeOctreeGrid3D::_sort_neighborinfo(std::vector<int>& index, std::vector<int>& neighborinfo)
{
	const int index_channels = 5;
	const int neighborinfo_channels = 4;

	int index_num = index.size()/index_channels;
	int neighborinfo_num_with_channels = neighborinfo.size();

	if(index_num <= 32) //32 is the CUDA warp size
		return;

	int* index_dat = new int[index_num*index_channels];
	memcpy(index_dat,&index[0],sizeof(int)*index_channels*index_num);

	int* neighborinfo_dat = new int[neighborinfo_num_with_channels];
	memcpy(neighborinfo_dat,&neighborinfo[0],sizeof(int)*neighborinfo_num_with_channels);

	const int max_neighbor_num = 8*8*6; // six faces, each face have 8*8 neighbors

	std::vector<int> buckets[max_neighbor_num+1];

	for(int i = 0;i < index_num;i++)
	{
		int offset = i*index_channels;

		int cur_neighbor_num = index_dat[offset+3];
		buckets[cur_neighbor_num].push_back(i);
	}

	int cur_offset = 0;
	int cur_neighborinfo_offset = 0;
	for(int bb = max_neighbor_num;bb >= 0;bb--)
	{
		int cur_size = buckets[bb].size();
		for(int cc = 0;cc < cur_size;cc++)
		{
			int offset = buckets[bb][cc]*index_channels;
			int neighborinfo_offset = index_dat[offset+4];
			int neighborinfo_num = index_dat[offset+3];

			index[cur_offset  ] = index_dat[offset  ];
			index[cur_offset+1] = index_dat[offset+1];
			index[cur_offset+2] = index_dat[offset+2];
			index[cur_offset+3] = neighborinfo_num;
			index[cur_offset+4] = cur_neighborinfo_offset;
			cur_offset += index_channels;

			for(int dd = 0;dd < neighborinfo_num;dd++)
			{
				neighborinfo[cur_neighborinfo_offset*neighborinfo_channels  ] = neighborinfo_dat[neighborinfo_offset*neighborinfo_channels  ];
				neighborinfo[cur_neighborinfo_offset*neighborinfo_channels+1] = neighborinfo_dat[neighborinfo_offset*neighborinfo_channels+1];
				neighborinfo[cur_neighborinfo_offset*neighborinfo_channels+2] = neighborinfo_dat[neighborinfo_offset*neighborinfo_channels+2];
				neighborinfo[cur_neighborinfo_offset*neighborinfo_channels+3] = neighborinfo_dat[neighborinfo_offset*neighborinfo_channels+3];
				cur_neighborinfo_offset ++;
				neighborinfo_offset ++;
			}
		}
	}

	delete []index_dat;
	delete []neighborinfo_dat;
}

