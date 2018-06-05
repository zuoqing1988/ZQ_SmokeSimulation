#include "ZQ_SmokeSimulationUtils3D.h"
#include "windows.h"
#include "ZQ_taucs.h"
#include "ZQ_TaucsBase.h"
#include "ZQ_SparseMatrix.h"
#include "ZQ_PCGSolver.h"
#include "ZQ_Matrix.h"
#include "ZQ_SVD.h"
#include "ZQ_MathBase.h"
#include "ZQ_CUDA_MatrixMul.h"

extern "C"
float SubProjection3DCuda(const int subRow, const int subCol, const int fast_count, const float* subMatrix, const float* input, float* output)
{
	return ZQ_CUDA_MatrixMul::MatrixMul(subMatrix, input, subRow, subCol, fast_count, output);
}

ZQ_SmokeSimulationPara3D::ZQ_SmokeSimulationPara3D()
{
	width = height = depth = 128;
	subSize = 4;
	coarseWidth = width/subSize;
	coarseHeight = height/subSize;
	coarseDepth = depth/subSize;
	gravityCoeff = 0;
	buoyancyCoeff = 1;
	confineCoeff = 0;
	Tamb = 10;
	voxelSize = 1.0 / width;
	
	ZQ_SmokeSource3D source;
	source.cx = 0.5;
	source.cy = 0.1;
	source.cz = 0.5;
	source.half_xlen = 0.1;
	source.half_ylen = 0.05;
	source.half_zlen = 0.1;
	source.reinjectDensity = 1;
	source.reinjectTemperature = Tamb * 20;
	sources.push_back(source);

	velAttenCoeff = 0;
	densityAttenCoeff = 0;
	tempAttenCoeff = 0;
	stepTime = 0.01;
	steps = 20;
	frameCount = 100;
	maxIter = 100;
	innerIter = 30;
	globalGridSolver = COARSE_OPEN_POISSON_SOR;
	subGridSolver = CLOSED_POISSON;
	octree_thresh1 = 2;
	octree_thresh2 = 4;
	octree_thresh3 = 16;
	useMovingObject = false;
	strcpy(movingObjectFile,"");
	useGuidingControl = false;
	guidingCoeff = 0;
	guidingNaiveBlend = true;
	exportGuide = false;
	strcpy(guidingFold,"");
	exportOccupy = false;
	randDensity = false;
	maxDensity = 0.8f;
	maxTemperature = 20.0f;
	intobj = INT_OBJ_NULL;
	advVelType = ADV_VEL_INREG_OUTREG;
	advScaType = ADV_SCA_REG;
	exportType = EXPORT_DI3;
	export_quality = 0.9999;
	addforce_type = ADD_FORCE_ENTIRE;
}

ZQ_SmokeSimulationPara3D::~ZQ_SmokeSimulationPara3D()
{
}

bool ZQ_SmokeSimulationPara3D::LoadFromFile(const char *file, bool display)
{
	FILE* config = fopen(file,"r");
	if(config == 0)
		return false;
	sources.clear();
	char s[1000];
	do{
		strcpy(s,"");
		fgets(s,999,config);
		int len = strlen(s);
		if(s[len] == '\n')
			s[--len] = '\0';
		char arg1[100] = {0},arg2[100]={0};
		sscanf(s,"%s%s",arg1,arg2);
		if(_strcmpi(arg1,"resolution") == 0)
		{
			sscanf(s,"%s%d%d%d",arg1,&width,&height,&depth);
		}
		else if(_strcmpi(arg1,"subSize") == 0)
			subSize = atoi(arg2);
		else if(_strcmpi(arg1,"Tamb") == 0)
			Tamb = atof(arg2);
		else if(_strcmpi(arg1,"velAttenCoeff") == 0)
			velAttenCoeff = atof(arg2);
		else if(_strcmpi(arg1,"tempAttenCoeff") == 0)
			tempAttenCoeff = atof(arg2);
		else if(_strcmpi(arg1,"densityAttenCoeff") == 0)
			densityAttenCoeff = atof(arg2);
		else if(_strcmpi(arg1,"confineCoeff") == 0)
			confineCoeff = atof(arg2);
		else if(_strcmpi(arg1,"source") == 0)
		{
			ZQ_SmokeSource3D source;
			sscanf(s,"%s%f%f%f%f%f%f%f%f",arg1,&source.cx,&source.cy,&source.cz,&source.half_xlen,&source.half_ylen,&source.half_zlen,&source.reinjectDensity,&source.reinjectTemperature);
			sources.push_back(source);
		}
		else if(_strcmpi(arg1,"steps") == 0)
			steps = atoi(arg2);
		else if(_strcmpi(arg1,"stepTime") == 0)
			stepTime = atof(arg2);
		else if(_strcmpi(arg1,"frameCount") == 0)
			frameCount = atoi(arg2);
		else if(_strcmpi(arg1,"maxIter") == 0)
			maxIter = atoi(arg2);
		else if(_strcmpi(arg1,"innerIter") == 0)
			innerIter = atoi(arg2);
		else if(_strcmpi(arg1,"globalGridSolver") == 0)
		{
			if(_strcmpi(arg2,"CLOSED_POISSON") == 0)
				globalGridSolver = CLOSED_POISSON;
			else if(_strcmpi(arg2,"OPEN_POISSON") == 0)
				globalGridSolver = OPEN_POISSON;
			else if(_strcmpi(arg2,"CLOSED_FLUX") == 0)
				globalGridSolver = CLOSED_FLUX;
			else if(_strcmpi(arg2,"OPEN_FLUX") == 0)
				globalGridSolver = OPEN_FLUX;		
			else if(_strcmpi(arg2,"CLOSED_OCTREE_FLUX") == 0)
				globalGridSolver = CLOSED_OCTREE_FLUX;
			else if(_strcmpi(arg2,"OPEN_OCTREE_FLUX") == 0)
				globalGridSolver = OPEN_OCTREE_FLUX;
			else if(_strcmpi(arg2,"CLOSED_OCTREE_POISSON") == 0)
				globalGridSolver = CLOSED_OCTREE_POISSON;
			else if(_strcmpi(arg2,"OPEN_OCTREE_POISSON") == 0)
				globalGridSolver = OPEN_OCTREE_POISSON;
			else if(_strcmpi(arg2,"TOTAL_OPEN_POISSON_SOR") == 0)
				globalGridSolver = TOTAL_OPEN_POISSON_SOR;
			else if(_strcmpi(arg2,"TOTAL_CLOSED_POISSON_SOR") == 0)
				globalGridSolver = TOTAL_CLOSED_POISSON_SOR;
			else if(_strcmpi(arg2,"TOTAL_OPEN_FLUX_SOR") == 0)
				globalGridSolver = TOTAL_OPEN_FLUX_SOR;
			else if(_strcmpi(arg2,"TOTAL_CLOSED_FLUX_SOR") == 0)
				globalGridSolver = TOTAL_CLOSED_FLUX_SOR;
			else if(_strcmpi(arg2,"COARSE_OPEN_POISSON_SOR") == 0)
				globalGridSolver = COARSE_OPEN_POISSON_SOR;
			else if(_strcmpi(arg2,"COARSE_CLOSED_POISSON_SOR") == 0)
				globalGridSolver = COARSE_CLOSED_POISSON_SOR;
			else if(_strcmpi(arg2,"COARSE_OPEN_FLUX_SOR") == 0)
				globalGridSolver = COARSE_OPEN_FLUX_SOR;
			else if(_strcmpi(arg2,"COARSE_CLOSED_FLUX_SOR") == 0)
				globalGridSolver = COARSE_CLOSED_FLUX_SOR;
			else if(_strcmpi(arg2,"OPEN_OCTREE_POISSON_SOR") == 0)
				globalGridSolver = OPEN_OCTREE_POISSON_SOR;
			else if(_strcmpi(arg2,"CLOSED_OCTREE_POISSON_SOR") == 0)
				globalGridSolver = CLOSED_OCTREE_POISSON_SOR;
			else
				globalGridSolver = OPEN_FLUX;
		}
		else if(_strcmpi(arg1,"octree_thresh1") == 0)
			octree_thresh1 = atoi(arg2);
		else if(_strcmpi(arg1,"octree_thresh2") == 0)
			octree_thresh2 = atoi(arg2);
		else if(_strcmpi(arg1,"octree_thresh3") == 0)
			octree_thresh3 = atoi(arg2);
		else if(_strcmpi(arg1,"useMovingObject") == 0)
		{
			if(_strcmpi(arg2,"true") == 0)
				useMovingObject = true;
			else
				useMovingObject = false;
		}
		else if(_strcmpi(arg1,"movingObjectFile") == 0)
		{
			strcpy(movingObjectFile,arg2);
		}
		else if(_strcmpi(arg1,"useGuidingControl") == 0)
		{
			if(_strcmpi(arg2,"true") == 0)
				useGuidingControl = true;
			else
				useGuidingControl = false;
		}
		else if(_strcmpi(arg1,"guidingCoeff") == 0)
		{
			guidingCoeff = atof(arg2);
		}
		else if(_strcmpi(arg1,"guidingNaiveBlend") == 0)
		{
			if(_strcmpi(arg2,"true") == 0)
				guidingNaiveBlend = true;
			else
				guidingNaiveBlend = false;
		}
		else if(_strcmpi(arg1,"guidingFold") == 0)
		{
			strcpy(guidingFold,arg2);
		}
		else if(_strcmpi(arg1,"exportGuide") == 0)
		{
			if(strcmp(arg2,"true") == 0)
				exportGuide = true;
			else
				exportGuide = false;
		}
		else if(_strcmpi(arg1,"exportOccupy") == 0)
		{
			if(_strcmpi(arg2,"true") == 0)
				exportOccupy = true;
			else
				exportOccupy = false;
		}
		else if(_strcmpi(arg1,"randDensity") == 0)
		{
			if(_strcmpi(arg2,"true") == 0)
				randDensity = true;
			else
				randDensity = false;
		}
		else if(_strcmpi(arg1,"intobj") == 0)
		{
			if(_strcmpi(arg2,"sphere") == 0)
				intobj = INT_OBJ_SPHERE;
			else if(_strcmpi(arg2,"oval") == 0)
				intobj = INT_OBJ_OVAL;
			else if(_strcmpi(arg2,"box") == 0)
				intobj = INT_OBJ_BOX;
			else if(_strcmpi(arg2,"special1") == 0)
				intobj = INT_OBJ_SPECIAL1;
			else
				intobj = INT_OBJ_NULL;
		}
		else if (_strcmpi(arg1, "advVelType") == 0)
		{
			if (_strcmpi(arg2, "inMAC_outMAC") == 0)
				advVelType = ADV_VEL_INMAC_OUTMAC;
			else if (_strcmpi(arg2, "inMAC_outMAC_BFECC") == 0)
				advVelType = ADV_VEL_INMAC_OUTMAC_BFECC;
			else if (_strcmpi(arg2, "inReg_outReg") == 0)
				advVelType = ADV_VEL_INREG_OUTREG;
			else if (_strcmpi(arg2, "inReg_outReg_BFECC") == 0)
				advVelType = ADV_VEL_INREG_OUTREG_BFECC;
			else if (_strcmpi(arg2, "inReg_outMAC") == 0)
				advVelType = ADV_VEL_INREG_OUTMAC;
			else if (_strcmpi(arg2, "inReg_outMAC_BFECC") == 0)
				advVelType = ADV_VEL_INREG_OUTMAC_BFECC;
			else
			{
				printf("unknown para: %s: %s\n", arg1, arg2);
				return false;
			}
		}
		else if (_strcmpi(arg1, "advScaType") == 0)
		{
			if (_strcmpi(arg2, "MAC") == 0)
				advScaType = ADV_SCA_MAC;
			else if (_strcmpi(arg2, "MAC_BFECC") == 0)
				advScaType = ADV_SCA_MAC_BFECC;
			else if (_strcmpi(arg2, "Reg") == 0)
				advScaType = ADV_SCA_REG;
			else if (_strcmpi(arg2, "Reg_BFECC") == 0)
				advScaType = ADV_SCA_REG_BFECC;
			else
			{
				printf("unknown para: %s: %s\n", arg1, arg2);
				return false;
			}
		}
		else if (_strcmpi(arg1, "exportType") == 0)
		{
			if (_strcmpi(arg2, "DI3") == 0)
				exportType = EXPORT_DI3;
			else if (_strcmpi(arg2, "ZQWV") == 0)
				exportType = EXPORT_ZQWV;
			else if (_strcmpi(arg2, "ZQCI") == 0)
				exportType = EXPORT_ZQCI;
			else
			{
				printf("unknown para: %s: %s\n", arg1, arg2);
				return false;
			}
		}
		else if (_strcmpi(arg1, "exportQuality") == 0)
		{
			export_quality = atof(arg2);
		}
		else if (_strcmpi(arg1, "AddForceType") == 0)
		{
			if (_strcmpi(arg2, "ENTIRE") == 0)
				addforce_type = ADD_FORCE_ENTIRE;
			else if (_strcmpi(arg2, "SLICE") == 0)
				addforce_type = ADD_FORCE_SLICE;
			else
			{
				printf("unknown para: %s: %s\n", arg1, arg2);
				return false;
			}
		}
		else if(_strcmpi(arg1,"maxDensity") == 0)
		{
			maxDensity = atof(arg2);
		}
		else if(_strcmpi(arg1,"maxTemperature") == 0)
		{
			maxTemperature = atof(arg2);
		}
	}while(strcmp(s,"") != 0);
	fclose(config);
	if(width%subSize != 0 || height%subSize != 0 || depth%subSize != 0)
	{
		printf("Error: Grid[%d x %d x %d] cannot be divided by subSize[%d]\n",width,height,depth,subSize);
		return false;
	}
	coarseWidth = width/subSize;
	coarseHeight = height/subSize;
	coarseDepth = depth/subSize;
	voxelSize = 1.0/width;
	if(sources.size() == 0)
	{
		printf("ERROR: no smoke source exist\n");
		return false;
	}

	if(display)
	{
		printf("resolution:\t\t%d x %d x %d\n",width,height,depth);
		printf("subSize:\t\t%d\n",subSize);
		printf("coarse:\t\t%d x %d x %d\n",coarseWidth,coarseHeight,coarseDepth);
		printf("Tamb:\t\t%f\n",Tamb);
		printf("confineCoeff:\t\t%f\n",confineCoeff);
		printf("velAttenCoeff:\t\t%f\n",velAttenCoeff);
		printf("tempAttenCoeff:\t\t%f\n",tempAttenCoeff);
		printf("densityAttenCoeff:\t\t%f\n",densityAttenCoeff);
		for(int cc = 0;cc < sources.size();cc++)
		{
			printf("Source[%d]: center(%.3f,%.3f,%.3f), size(%.3f,%.3f,%.3f), density(%.3f), temperature(%.3f)\n",
				cc,sources[cc].cx,sources[cc].cy,sources[cc].cz,sources[cc].half_xlen,sources[cc].half_ylen,sources[cc].half_zlen,
							sources[cc].reinjectDensity,sources[cc].reinjectTemperature);
		}
		printf("voxelSize:\t\t%f\n",voxelSize);
		printf("steps:\t\t%d\n",steps);
		printf("stepTime:\t\t%f\n",stepTime);
		printf("frameCount:\t\t%d\n",frameCount);
	}
	return true;
}


void MACtoRegular4(int width, int height, int depth, const float* mac_u, const float* mac_v, const float* mac_w, float* vel4)
{
	for (int k = 0; k < depth; k++)
	{
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				int offset = k*height*width + j*width + i;
				int u_offset = k*height*(width + 1) + j*(width + 1) + i;
				int v_offset = k*(height + 1)*width + j*width + i;
				int w_offset = k*height*width + j*width + i;
				vel4[offset * 4 + 0] = 0.5*(mac_u[u_offset] + mac_u[u_offset + 1]);
				vel4[offset * 4 + 1] = 0.5*(mac_v[v_offset] + mac_v[v_offset + width]);
				vel4[offset * 4 + 2] = 0.5*(mac_w[w_offset] + mac_w[w_offset + height*width]);
				vel4[offset * 4 + 3] = 0;
			}
		}
	}
}

void Regular4toMAC(int width, int height, int depth, const float* vel4, float* mac_u, float* mac_v, float* mac_w)
{
	for (int k = 0; k < depth; k++)
	{
		for (int j = 0; j < height; j++)
		{
			int u_offset = k*height*(width + 1) + j*(width + 1);
			int offset = k*height*width + j*width;
			for (int i = 1; i < width; i++)
				mac_u[u_offset + i] = 0.5*(vel4[(offset + i) * 4 + 0] + vel4[(offset + i - 1) * 4 + 0]);
			mac_u[u_offset + 0] = vel4[offset * 4 + 0];
			mac_v[u_offset + width] = vel4[(offset + width - 1) * 4 + 0];
		}
	}

	for (int k = 0; k < depth; k++)
	{
		for (int i = 0; i < width; i++)
		{
			int v_offset = k*(height + 1)*width + i;
			int offset = k*height*width + i;
			for (int j = 1; j < height; j++)
				mac_v[v_offset + j*width] = 0.5*(vel4[(offset + j*width) * 4 + 1] + vel4[(offset + (j - 1)*width) * 4 + 1]);
			mac_v[v_offset + 0] = vel4[offset * 4 + 1];
			mac_v[v_offset + height*width] = vel4[(offset + (height - 1)*width) * 4 + 1];
		}
	}

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int w_offset = j*width + i;
			int offset = j*width + i;
			for (int k = 1; k < depth; k++)
				mac_w[w_offset + k*height*width] = 0.5*(vel4[(offset + k*height*width) * 4 + 2] + vel4[(offset + (k - 1)*height*width) * 4 + 2]);
			mac_w[w_offset + 0] = vel4[offset * 4 + 2];
			mac_w[w_offset + depth*height*width] = vel4[(offset + (depth - 1)*height*width) * 4 + 2];
		}
	}
}


void ApplyAdvectVelocityResult(int width, int height, int depth, const bool* occupy, const float* adv_u, const float* adv_v, const float* adv_w, float* u, float* v, float* w)
{
	for (int k = 0; k < depth; k++)
	{
		for (int j = 0; j < height; j++)
		{
			for (int i = 1; i < width; i++)
			{
				if (occupy[k*height*width + j*width + i] || occupy[k*height*width + j*width + i - 1])
					;
				else
					u[k*height*(width + 1) + j*(width + 1) + i] = adv_u[k*height*(width + 1) + j*(width + 1) + i];
			}
			if (!occupy[k*height*width + j*width + 0])
				u[k*height*(width + 1) + j*(width + 1) + 0] = adv_u[k*height*(width + 1) + j*(width + 1) + 0];
			if (!occupy[k*height*width + j*width + width - 1])
				u[k*height*(width + 1) + j*(width + 1) + width] = adv_u[k*height*(width + 1) + j*(width + 1) + width];
		}
	}

	for (int k = 0; k < depth; k++)
	{
		for (int i = 0; i < width; i++)
		{
			for (int j = 1; j < height; j++)
			{
				if (occupy[k*height*width + j*width + i] || occupy[k*height*width + (j - 1)*width + i])
					;
				else
					v[k*(height + 1)*width + j*width + i] = adv_v[k*(height + 1)*width + j*width + i];
			}
			if (!occupy[k*height*width + 0 * width + i])
				v[k*(height + 1)*width + 0 * width + i] = adv_v[k*(height + 1)*width + 0 * width + i];
			if (!occupy[k*height*width + (height - 1)*width + i])
				v[k*(height + 1)*width + height*width + i] = adv_v[k*(height + 1)*width + height*width + i];
		}
	}

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			for (int k = 1; k < depth; k++)
			{
				if (occupy[k*height*width + j*width + i] || occupy[(k - 1)*height*width + j*width + i])
					;
				else
					w[k*height*width + j*width + i] = adv_w[k*height*width + j*width + i];
			}
			if (!occupy[0 * height*width + j*width + i])
				w[j*width + i] = adv_w[j*width + i];
			if (!occupy[(depth - 1)*height*width + j*width + i])
				w[depth*height*width + j*width + i] = adv_w[depth*height*width + j*width + i];
		}
	}
}

int ZQ_InitMatOpenPoisson3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A)
{

	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();

	int* index = new int[width*height*depth];
	memset(index,0,sizeof(int)*width*height*depth);
	int idx = 1;
	
	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();

	for(int i = 0;i < depth*height*width;i++)
	{
		if(!occupyPtr[i])
			index[i] = idx++;
	}

	if(idx == 1)
	{
		delete []index;
		return 0;
	}


	int dim = idx-1;
	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	int KSLICE = height*width;
	int JSLICE = width;
	int ISLICE = 1;

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int offset = k*height*width+j*width+i;
				int row_id = index[offset]-1;
				if(row_id < 0)
					continue;

				std::vector<int> indices;
				std::vector<int> u_indices;
				std::vector<int> v_indices;
				std::vector<int> w_indices;
				if(k > 0)
				{
					indices.push_back(offset-KSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(-1);
					w_indices.push_back(k*height*width+j*width+i);
				}
				if(k < depth-1)
				{
					indices.push_back(offset+KSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(-1);
					w_indices.push_back((k+1)*height*width+j*width+i);
				}
				if(j > 0)
				{
					indices.push_back(offset-JSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(k*(height+1)*width+j*width+i);
					w_indices.push_back(-1);
				}
				if(j < height-1)
				{
					indices.push_back(offset+JSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(k*(height+1)*width+(j+1)*width+i);
					w_indices.push_back(-1);
				}
				if(i > 0)
				{
					indices.push_back(offset-ISLICE);
					u_indices.push_back(k*height*(width+1)+j*(width+1)+i);
					v_indices.push_back(-1);
					w_indices.push_back(-1);
				}
				if(i < width-1)
				{
					indices.push_back(offset+ISLICE);
					u_indices.push_back(k*height*(width+1)+j*(width+1)+i+1);
					v_indices.push_back(-1);
					w_indices.push_back(-1);
				}

				float count = 6-indices.size();
				for(int cc = 0;cc < indices.size();cc++)
				{
					int cur_index = index[indices[cc]];
					if(cur_index != 0)
					{
						float tmp_ratio = 0.0f;
						if(u_indices[cc] != -1)
							tmp_ratio = unoccupyUptr[u_indices[cc]];
						if(v_indices[cc] != -1)
							tmp_ratio = unoccupyVptr[v_indices[cc]];
						if(w_indices[cc] != -1)
							tmp_ratio = unoccupyWptr[w_indices[cc]];
						count += tmp_ratio;

						mat.SetValue(row_id,cur_index-1, tmp_ratio);
					}
				}

				mat.SetValue(row_id,row_id, -count);
			}
		}
	}
	


	delete []index;

	*A = mat.ExportCCS(TAUCS_DOUBLE);

	return dim;
}

int ZQ_InitMatClosedPoisson3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A)
{	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetWidth();
	int depth = coarseGrid->GetDepth();

	int* index = new int[width*height*depth];
	memset(index,0,sizeof(int)*width*height*depth);
	
	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();

	int idx = 1;
	bool haveFirst = false;
	int first = -1;
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(!occupyPtr[k*height*width+j*width+i])
				{
					if(!haveFirst)
					{
						haveFirst = true;
						first = k*height*width+j*width+i;
					}
					index[k*height*width+j*width+i] = idx++;
				}
			}
		}
	}
	
	if(!haveFirst || first == -1)
	{
		delete []index;
		return 0;
	}

	
	int row = idx-1;
	int dim = idx-2;

	ZQ::ZQ_SparseMatrix<float> mat(row,dim);

	int KSLICE = height*width;
	int JSLICE = width;
	int ISLICE = 1;

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int offset = k*height*width+j*width+i;
				int row_id = index[offset]-1;

				if(row_id < 0)
					continue;

				std::vector<int> indices;
				std::vector<int> u_indices;
				std::vector<int> v_indices;
				std::vector<int> w_indices;

				if(k > 0)
				{
					indices.push_back(offset-KSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(-1);
					w_indices.push_back(k*height*width+j*width+i);
				}
				if(k < depth-1)
				{
					indices.push_back(offset+KSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(-1);
					w_indices.push_back((k+1)*height*width+j*width+i);
				}
				if(j > 0)
				{
					indices.push_back(offset-JSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(k*(height+1)*width+j*width+i);
					w_indices.push_back(-1);
				}
				if(j < height-1)
				{
					indices.push_back(offset+JSLICE);
					u_indices.push_back(-1);
					v_indices.push_back(k*(height+1)*width+(j+1)*width+i);
					w_indices.push_back(-1);
				}
				if(i > 0)
				{
					indices.push_back(offset-ISLICE);
					u_indices.push_back(k*height*(width+1)+j*(width+1)+i);
					v_indices.push_back(-1);
					w_indices.push_back(-1);
				}
				if(i < width-1)
				{
					indices.push_back(offset+ISLICE);
					u_indices.push_back(k*height*(width+1)+j*(width+1)+i+1);
					v_indices.push_back(-1);
					w_indices.push_back(-1);
				}

				float count = 0.0f;
				for(int cc = 0;cc < indices.size();cc++)
				{
					int cur_index = index[indices[cc]];
					if(cur_index != 0)
					{
						float tmp_ratio = 0;
						if(u_indices[cc] != -1)
							tmp_ratio = unoccupyUptr[u_indices[cc]];
						if(v_indices[cc] != -1)
							tmp_ratio = unoccupyVptr[v_indices[cc]];
						if(w_indices[cc] != -1)
							tmp_ratio = unoccupyWptr[w_indices[cc]];
						count += tmp_ratio;
						if(cur_index != 1)
							mat.SetValue(row_id,cur_index-2, tmp_ratio);
					}
				}

				if(index[offset] != 1)
					mat.SetValue(row_id,index[offset]-2, -count);

			}
		}
	}
	

	delete []index;

	*A = mat.ExportCCS(TAUCS_DOUBLE);

	return dim;
}

int ZQ_InitMatOpenFlux3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A)
{
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();

	int* u_index = new int[(width+1)*height*depth];
	int* v_index = new int[width*(height+1)*depth];
	int* w_index = new int[width*height*(depth+1)];
	int* lamda_index = new int[width*height*depth];

	for(int i = 0;i < (width+1)*height*depth;i++)
		u_index[i] = -1;
	for(int i = 0;i < width*(height+1)*depth;i++)
		v_index[i] = -1;
	for(int i = 0;i < width*height*(depth+1);i++)
		w_index[i] = -1;
	for(int i = 0;i < width*height*depth;i++)
		lamda_index[i] = -1;

	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();

	int idx = 0;
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i <= width;i++)
			{
				if(unoccupyUptr[k*height*(width+1)+j*(width+1)+i] != 0)
					u_index[k*height*(width+1)+j*(width+1)+i] = idx++;
			}
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j <= height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(unoccupyVptr[k*(height+1)*width+j*width+i] != 0)
					v_index[k*(height+1)*width+j*width+i] = idx++;
			}
		}
	}

	for(int k = 0;k <= depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(unoccupyWptr[k*height*width+j*width+i] != 0)
					w_index[k*height*width+j*width+i] = idx++;
			}
		}
	}

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(!occupyPtr[k*height*width+j*width+i])
					lamda_index[k*height*width+j*width+i] = idx++;
			}
		}
	}
	

	int dim = idx;

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	idx = 0;
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i <= width;i++)
			{
				int cur_u_index = u_index[k*height*(width+1)+j*(width+1)+i];

				if(cur_u_index >= 0)
				{
					float ratio = unoccupyUptr[k*height*(width+1)+j*(width+1)+i];
					mat.SetValue(cur_u_index,cur_u_index, 2*ratio);
					if(i > 0)
						mat.SetValue(cur_u_index,lamda_index[k*height*width+j*width+i-1], -ratio);
					if(i < width)
						mat.SetValue(cur_u_index,lamda_index[k*height*width+j*width+i], ratio);
				}
			}
		}
	}
	

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j <= height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int cur_v_index = v_index[k*(height+1)*width+j*width+i];

				if(cur_v_index >= 0)
				{
					float ratio = unoccupyVptr[k*(height+1)*width+j*width+i];
					mat.SetValue(cur_v_index,cur_v_index, 2*ratio);
					if(j > 0)
						mat.SetValue(cur_v_index,lamda_index[k*height*width+(j-1)*width+i], -ratio);
					if(j < height)
						mat.SetValue(cur_v_index,lamda_index[k*height*width+j*width+i], ratio);
				}
			}
		}
	}
	
	for(int k = 0;k <= depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int cur_w_index = w_index[k*height*width+j*width+i];
				
				if(cur_w_index >= 0)
				{
					float ratio = unoccupyWptr[k*height*width+j*width+i];
					mat.SetValue(cur_w_index,cur_w_index, 2*ratio);
					if(k > 0)
						mat.SetValue(cur_w_index,lamda_index[(k-1)*height*width+j*width+i], -ratio);
					if(k < depth)
						mat.SetValue(cur_w_index,lamda_index[k*height*width+j*width+i], ratio);
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
				int cur_lambda_index = lamda_index[k*height*width+j*width+i];
				if(cur_lambda_index >= 0)
				{
					float ratio = 0;
					ratio = unoccupyUptr[k*height*(width+1)+j*(width+1)+i+1];
					if(ratio > 0)
						mat.SetValue(cur_lambda_index,u_index[k*height*(width+1)+j*(width+1)+i+1], -ratio);
					ratio = unoccupyUptr[k*height*(width+1)+j*(width+1)+i];
					if(ratio > 0)
						mat.SetValue(cur_lambda_index,u_index[k*height*(width+1)+j*(width+1)+i], ratio);
					ratio = unoccupyVptr[k*(height+1)*width+(j+1)*width+i];
					if(ratio > 0)
						mat.SetValue(cur_lambda_index,v_index[k*(height+1)*width+(j+1)*width+i], -ratio);
					ratio = unoccupyVptr[k*(height+1)*width+j*width+i];
					if(ratio > 0)
						mat.SetValue(cur_lambda_index,v_index[k*(height+1)*width+j*width+i], ratio);
					ratio = unoccupyWptr[(k+1)*height*width+j*width+i];
					if(ratio > 0)
						mat.SetValue(cur_lambda_index,w_index[(k+1)*height*width+j*width+i], -ratio);
					ratio = unoccupyWptr[k*height*width+j*width+i];
					if(ratio > 0)
						mat.SetValue(cur_lambda_index,w_index[k*height*width+j*width+i], ratio);
				}
			}
		}
	}
	

	*A = mat.ExportCCS(TAUCS_DOUBLE);

	delete []u_index;
	delete []v_index;
	delete []w_index;
	delete []lamda_index;

	return dim;
}

int ZQ_InitMatClosedFlux3D(const ZQ_SmokeGrid3D* coarseGrid, taucs_ccs_matrix** A)
{
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();
	
	int* u_index = new int[(width+1)*height*depth];
	int* v_index = new int[width*(height+1)*depth];
	int* w_index = new int[width*height*(depth+1)];
	int* lamda_index = new int[width*height*depth];
	int k_index = 0;

	for(int i = 0;i < (width+1)*height*depth;i++)
		u_index[i] = -1;

	for(int i = 0;i < width*(height+1)*depth;i++)
		v_index[i] = -1;

	for(int i = 0;i < width*height*(depth+1);i++)
		w_index[i] = -1;

	for(int i = 0;i < width*height*depth;i++)
		lamda_index[i] = -1;

	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();
	const float*& unoccupyVolumePtr = coarseGrid->GetVolumeUnOccupyRatioPtr();

	int idx = 0;
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 1;i < width;i++)
			{
				if(unoccupyUptr[k*height*(width+1)+j*(width+1)+i] != 0)
					u_index[k*height*(width+1)+j*(width+1)+i] = idx++;
			}
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 1;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(unoccupyVptr[k*(height+1)*width+j*width+i] != 0)
					v_index[k*(height+1)*width+j*width+i] = idx++;
			}
		}
	}
	
	
	for(int k = 1;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(unoccupyWptr[k*height*width+j*width+i] != 0)
					w_index[k*height*width+j*width+i] = idx++;
			}
		}
	}

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(!occupyPtr[k*height*width+j*width+i])
					lamda_index[k*height*width+j*width+i] = idx++;
			}
		}
	}
	
	k_index = idx++;
	int dim = idx;
	idx = 0;

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 1;i < width;i++)
			{
				int cur_u_index = u_index[k*height*(width+1)+j*(width+1)+i];
				if(cur_u_index >= 0)
				{
					float ratio = unoccupyUptr[k*height*(width+1)+j*(width+1)+i];
					mat.SetValue(cur_u_index,cur_u_index, 2*ratio);
					mat.SetValue(cur_u_index,lamda_index[k*height*width+j*width+i-1], -ratio);
					mat.SetValue(cur_u_index,lamda_index[k*height*width+j*width+i], ratio);
				}
			}
		}
	}
	

	for(int k = 0;k < depth;k++)
	{
		for(int j = 1;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int cur_v_index = v_index[k*(height+1)*width+j*width+i];
				if(cur_v_index >= 0)
				{
					float ratio = unoccupyVptr[k*(height+1)*width+j*width+i];
					mat.SetValue(cur_v_index,cur_v_index, 2*ratio);
					mat.SetValue(cur_v_index,lamda_index[k*height*width+(j-1)*width+i], -ratio);
					mat.SetValue(cur_v_index,lamda_index[k*height*width+j*width+i], ratio);
				}
			}
		}
	}

	for(int k = 1;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int cur_w_index = w_index[k*height*width+j*width+i];
				if(cur_w_index >= 0)
				{
					float ratio = unoccupyWptr[k*height*width+j*width+i];
					mat.SetValue(cur_w_index,cur_w_index,2*ratio);
					mat.SetValue(cur_w_index,lamda_index[(k-1)*height*width+j*width+i], -ratio);
					mat.SetValue(cur_w_index,lamda_index[k*height*width+j*width+i], ratio);
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
				int cur_lambda_index = lamda_index[k*height*width+j*width+i];
				if(cur_lambda_index >= 0)
				{
					float ratio = 0;
					ratio = unoccupyUptr[k*height*(width+1)+j*(width+1)+i+1];
					if(ratio > 0 && u_index[k*height*(width+1)+j*(width+1)+i+1] >= 0)
						mat.SetValue(cur_lambda_index,u_index[k*height*(width+1)+j*(width+1)+i+1], -ratio);

					ratio = unoccupyUptr[k*height*(width+1)+j*(width+1)+i];
					if(ratio > 0 && u_index[k*height*(width+1)+j*(width+1)+i] >= 0)
						mat.SetValue(cur_lambda_index,u_index[k*height*(width+1)+j*(width+1)+i], ratio);

					ratio = unoccupyVptr[k*(height+1)*width+(j+1)*width+i];
					if(ratio > 0 && v_index[k*(height+1)*width+(j+1)*width+i] >= 0)
						mat.SetValue(cur_lambda_index,v_index[k*(height+1)*width+(j+1)*width+i], -ratio);

					ratio = unoccupyVptr[k*(height+1)*width+j*width+i];
					if(ratio > 0 && v_index[k*(height+1)*width+j*width+i] >= 0)
						mat.SetValue(cur_lambda_index,v_index[k*(height+1)*width+j*width+i], ratio);

					ratio = unoccupyWptr[(k+1)*height*width+j*width+i];
					if(ratio > 0 && w_index[(k+1)*height*width+j*width+i] >= 0)
						mat.SetValue(cur_lambda_index,w_index[(k+1)*height*width+j*width+i], -ratio);

					ratio = unoccupyWptr[k*height*width+j*width+i];
					if(ratio > 0 && w_index[k*height*width+j*width+i] >= 0)
						mat.SetValue(cur_lambda_index,w_index[k*height*width+j*width+i], ratio);

					ratio = unoccupyVolumePtr[k*height*width+j*width+i];
					if(ratio > 0)
						mat.SetValue(cur_lambda_index,k_index, ratio);
				}
			}
		}
	}
	


	//k_index
	{
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(lamda_index[k*height*width+j*width+i] >= 0)
					{
						float ratio = unoccupyVolumePtr[k*height*width+j*width+i];
						mat.SetValue(k_index,lamda_index[k*height*width+j*width+i], ratio);
					}
				}
			}
		}
	}

	*A = mat.ExportCCS(TAUCS_DOUBLE);
	
	delete []u_index;
	delete []v_index;
	delete []w_index;
	delete []lamda_index;

	return dim;
}

bool ZQ_InitSubMatPoisson3D(const int subSize, float** subMatrix, int* row, int* col)
{
	typedef ZQ::ZQ_Matrix<double> ZQ_Mat;
	ZQ_Mat matP(subSize*subSize*subSize,subSize*subSize*subSize-1);
	ZQ_Mat matV(subSize*subSize*subSize,3*subSize*subSize*(subSize+1));
	ZQ_Mat matD(3*subSize*subSize*(subSize+1),subSize*subSize*subSize);


	int ushift = 0;
	int vshift = subSize*subSize*(subSize+1);
	int wshift = 2*subSize*subSize*(subSize+1);

	for(int k = 0;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				int vcol1 = ushift + k*subSize*(subSize+1)+j*(subSize+1)+i;
				int vcol2 = ushift + k*subSize*(subSize+1)+j*(subSize+1)+i+1;
				int vcol3 = vshift + k*(subSize+1)*subSize+j*subSize+i;
				int vcol4 = vshift + k*(subSize+1)*subSize+(j+1)*subSize+i;
				int vcol5 = wshift + k*subSize*subSize+j*subSize+i;
				int vcol6 = wshift + (k+1)*subSize*subSize+j*subSize+i;

				int vrow = k*subSize*subSize+j*subSize+i;
				matV.SetData(vrow,vcol1,-1);
				matV.SetData(vrow,vcol2,1);
				matV.SetData(vrow,vcol3,-1);
				matV.SetData(vrow,vcol4,1);
				matV.SetData(vrow,vcol5,-1);
				matV.SetData(vrow,vcol6,1);

				int prow = k*subSize*subSize+j*subSize+i;
				int pcol0 = prow - 1;
				int pcol1 = k*subSize*subSize+j*subSize+(i-1) - 1;
				int pcol2 = k*subSize*subSize+j*subSize+(i+1) - 1;
				int pcol3 = k*subSize*subSize+(j-1)*subSize+i - 1;
				int pcol4 = k*subSize*subSize+(j+1)*subSize+i - 1;
				int pcol5 = (k-1)*subSize*subSize+j*subSize+i - 1;
				int pcol6 = (k+1)*subSize*subSize+j*subSize+i - 1;

				if(prow == 0)
				{
					matP.SetData(0,pcol2,1);
					matP.SetData(0,pcol4,1);
					matP.SetData(0,pcol6,1);
					continue;
				}

				int count = 0;
				if(i == 0)
				{
					count ++;
					matP.SetData(prow,pcol2,1);
				}
				else if(i == subSize-1)
				{
					count ++;
					if(pcol1 != -1)
						matP.SetData(prow,pcol1,1);
				}
				else
				{
					count += 2;
					if(pcol1 != -1)
						matP.SetData(prow,pcol1,1);
					matP.SetData(prow,pcol2,1);
				}
				if(j == 0)
				{
					count ++;
					matP.SetData(prow,pcol4,1);
				}
				else if(j == subSize-1)
				{
					count ++;
					if(pcol3 != -1)
						matP.SetData(prow,pcol3,1);
				}
				else
				{
					count += 2;
					if(pcol3 != -1)
						matP.SetData(prow,pcol3,1);
					matP.SetData(prow,pcol4,1);
				}
				if(k == 0)
				{
					count ++;
					matP.SetData(prow,pcol6,1);
				}
				else if(k == subSize-1)
				{
					count ++;
					if(pcol5 != -1)
						matP.SetData(prow,pcol5,1);
				}
				else 
				{
					count += 2;
					if(pcol5 != -1)
						matP.SetData(prow,pcol5,1);
					matP.SetData(prow,pcol6,1);
				}

				matP.SetData(prow,pcol0,-count);
			}
		}
	}
	
	for(int k = 0;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 1;i < subSize;i++)
			{
				int drow = ushift + k*subSize*(subSize+1)+j*(subSize+1)+i;
				int dcol1 = k*subSize*subSize+j*subSize+i-1;
				int dcol2 = k*subSize*subSize+j*subSize+i;
				matD.SetData(drow,dcol1,-1);
				matD.SetData(drow,dcol2,1);
			}
		}
	}

	for(int k = 0;k < subSize;k++)
	{
		for(int j = 1;j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				int drow = vshift + k*(subSize+1)*subSize+j*subSize+i;
				int dcol1 = k*subSize*subSize+(j-1)*subSize+i;
				int dcol2 = k*subSize*subSize+j*subSize+i;
				matD.SetData(drow,dcol1,-1);
				matD.SetData(drow,dcol2,1);
			}
		}
	}

	for(int k = 1;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				int drow = wshift + k*subSize*subSize+j*subSize+i;
				int dcol1 = (k-1)*subSize*subSize+j*subSize+i;
				int dcol2 = k*subSize*subSize+j*subSize+i;
				matD.SetData(drow,dcol1,-1);
				matD.SetData(drow,dcol2,1);
			}
		}
	}
	
	
	ZQ_Mat invP(subSize*subSize*subSize-1,subSize*subSize*subSize);

	ZQ::ZQ_SVD::Invert(matP,invP);
	ZQ_Mat invPV = invP*matV;

	ZQ_Mat invPV1(subSize*subSize*subSize,3*subSize*subSize*(subSize+1));
	bool flag;
	for(int rr = 1; rr < subSize*subSize*subSize;rr++)
	{
		for(int cc = 0;cc < 3*subSize*subSize*(subSize+1);cc++)
			invPV1.SetData(rr,cc,invPV.GetData(rr-1,cc,flag));
	}

	ZQ_Mat DinvPV = matD*invPV1;

	int tmpheight = 3*subSize*subSize*(subSize+1);
	int tmpwidth = 3*subSize*subSize*(subSize+1);
	int height = 3*subSize*subSize*(subSize-1);
	int width = 3*subSize*subSize*(subSize+1);

	if(*subMatrix != 0)
		delete [](*subMatrix);
	*subMatrix = new float[width*height];
	float* tmpMat = new float[tmpheight*tmpwidth];

	for(int i = 0;i < tmpheight;i++)
	{
		for(int j = 0;j < tmpwidth;j++)
		{
			if(i == j)
				tmpMat[i*tmpwidth+j] = 1 - DinvPV.GetData(i,j,flag);
			else
				tmpMat[i*tmpwidth+j] = -DinvPV.GetData(i,j,flag);
		}
	}
	
	int useCount = 0;
	for(int k = 0;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 1;i < subSize;i++)
			{
				int wholeCount = k*subSize*(subSize+1)+j*(subSize+1)+i;
				memcpy((*subMatrix)+useCount*width,tmpMat+wholeCount*tmpwidth,sizeof(float)*width);
				useCount ++;
			}
		}
	}
	
	for(int k = 0;k < subSize;k++)
	{
		for(int j = 1; j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				int wholeCount = subSize*subSize*(subSize+1) + k*(subSize+1)*subSize+j*subSize+i;
				memcpy((*subMatrix)+useCount*width,tmpMat+wholeCount*tmpwidth,sizeof(float)*width);
				useCount ++;
			}
		}
	}

	for(int k = 1;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				int wholeCount = 2*subSize*subSize*(subSize+1) + k*subSize*subSize+j*subSize+i;
				memcpy((*subMatrix)+useCount*width,tmpMat+wholeCount*tmpwidth,sizeof(float)*width);
				useCount ++;
			}
		}
	}
	
	

	delete []tmpMat;
	*row = height;
	*col = width;

	return true;
}

bool ZQ_SolveSubGrid3D(const ZQ_SmokeSimulationPara3D& para, const bool* occupy, const float* input, float* output)
{
	if(occupy == 0 || input == 0 || output == 0)
		return false;
	int subSize = para.subSize;
	ZQ_SmokeGrid3D* tmpGrid = new ZQ_SmokeGrid3D(subSize,subSize,subSize);
	
	tmpGrid->BuildFaceUnOccupyRatio();
	float*& unoccupyUptr = tmpGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyVptr = tmpGrid->GetFaceUnOccupyRatioVPtr();
	float*& unoccupyWptr = tmpGrid->GetFaceUnOccupyRatioWPtr();
	bool*& occupyPtr = tmpGrid->GetOccupyPtr();
	
	for(int k = 0;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 0;i <= subSize;i++)
			{
				unoccupyUptr[k*subSize*(subSize+1)+j*(subSize+1)+i] = 1;
			}
		}
	}
	for(int k = 0;k < subSize;k++)
	{
		for(int j = 0;j <= subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				unoccupyVptr[k*(subSize+1)*subSize+j*subSize+i] = 1;
			}
		}
	}
	for(int k = 0;k <= subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				unoccupyWptr[k*subSize*subSize+j*subSize+i] = 1;
			}
		}
	}
	
	for(int k = 0;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				occupyPtr[k*subSize*subSize+j*subSize+i] = occupy[k*subSize*subSize+j*subSize+i];
				if(occupy[k*subSize*subSize+j*subSize+i])
				{
					unoccupyUptr[k*subSize*(subSize+1)+j*(subSize+1)+i] = 0;
					unoccupyUptr[k*subSize*(subSize+1)+j*(subSize+1)+i+1] = 0;
					unoccupyVptr[k*(subSize+1)*subSize+j*subSize+i] = 0;
					unoccupyVptr[k*(subSize+1)*subSize+(j+1)*subSize+i] = 0;
					unoccupyWptr[k*subSize*subSize+j*subSize+i] = 0;
					unoccupyWptr[(k+1)*subSize*subSize+j*subSize+i] = 0;
				}
			}
		}
	}
	

	tmpGrid->BuildVolumeUnOccupyRatio();
	float*& unoccupyVolumePtr = tmpGrid->GetVolumeUnOccupyRatioPtr();
	for(int k = 0;k < subSize;k++)
	{
		for(int j = 0;j < subSize;j++)
		{
			for(int i = 0;i < subSize;i++)
			{
				if(occupy[k*subSize*subSize+j*subSize+i])
					unoccupyVolumePtr[k*subSize*subSize+j*subSize+i] = 0;
				else
					unoccupyVolumePtr[k*subSize*subSize+j*subSize+i] = 1;
			}
		}
	}
	

	float*& U = tmpGrid->GetVelocityUptr();
	float*& V = tmpGrid->GetVelocityVptr();
	float*& W = tmpGrid->GetVelocityWptr();
	memcpy(U,input,sizeof(float)*subSize*subSize*(subSize+1));
	memcpy(V,input+subSize*subSize*(subSize+1),sizeof(float)*subSize*subSize*(subSize+1));
	memcpy(W,input+2*subSize*subSize*(subSize+1),sizeof(float)*subSize*subSize*(subSize+1));
	
	taucs_ccs_matrix* tmpA = 0;
	
	ZQ_SmokeSimulationPara3D sub_para;
	sub_para.maxIter = 2*(subSize*subSize*subSize+3*subSize*subSize*(subSize+1)+1);
	
	sub_para.globalGridSolver = CLOSED_POISSON;
	sub_para.width = subSize;
	sub_para.height = subSize;
	sub_para.depth = subSize;
	sub_para.subSize = 1;

	ZQ_SmokeGrid3D* deltaGrid = new ZQ_SmokeGrid3D(subSize,subSize,subSize);
	
	if(!ZQ_InitMatClosedPoisson3D(tmpGrid,&tmpA))
	{
		delete tmpGrid;
		delete deltaGrid;
		return false;
	}
	ZQ_SolveClosedPressure3D(tmpGrid,sub_para,tmpA);
	ZQ_AdjustClosedVelocity3D(tmpGrid,deltaGrid);
	

	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	float*& deltaW = deltaGrid->GetVelocityWptr();

	for(int i = 0;i < subSize*subSize*(subSize+1);i++)
	{
		output[i] = U[i]+deltaU[i];
		output[i+subSize*subSize*(subSize+1)] = V[i]+deltaV[i];
		output[i+2*subSize*subSize*(subSize+1)] = W[i]+deltaW[i];
	}
	if(tmpA)
		ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(tmpA);
	
	delete tmpGrid;
	delete deltaGrid;
	return true;
}

bool ZQ_InitFineGrid3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	int width = para.width;
	int height = para.height;
	int depth = para.depth;
	
	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	float*& pressurePtr = fineGrid->GetPressurePtr();

	memset(occupyPtr,0,sizeof(bool)*width*height*depth);
	memset(pressurePtr,0,sizeof(float)*width*height*depth);

	if(para.globalGridSolver == CLOSED_POISSON || para.globalGridSolver == CLOSED_FLUX 
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON || para.globalGridSolver == CLOSED_OCTREE_FLUX
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON_SOR
		|| para.globalGridSolver == TOTAL_CLOSED_POISSON_SOR || para.globalGridSolver == TOTAL_CLOSED_FLUX_SOR
		|| para.globalGridSolver == COARSE_CLOSED_POISSON_SOR || para.globalGridSolver == COARSE_CLOSED_FLUX_SOR)
	{
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				occupyPtr[k*height*width+j*width+0] = true;
			/*	occupyPtr[k*height*width+j*width+1] = true;
				occupyPtr[k*height*width+j*width+2] = true;
				occupyPtr[k*height*width+j*width+3] = true;
				occupyPtr[k*height*width+j*width+width-4] = true;
				occupyPtr[k*height*width+j*width+width-3] = true;
				occupyPtr[k*height*width+j*width+width-2] = true;*/
				occupyPtr[k*height*width+j*width+width-1] = true;
			}
		}
		for(int k = 0;k < depth;k++)
		{
			for(int i = 0;i < width;i++)
			{
				occupyPtr[k*height*width+0*width+i] = true;
				/*occupyPtr[k*height*width+1*width+i] = true;
				occupyPtr[k*height*width+2*width+i] = true;
				occupyPtr[k*height*width+3*width+i] = true;
				occupyPtr[k*height*width+(height-4)*width+i] = true;
				occupyPtr[k*height*width+(height-3)*width+i] = true;
				occupyPtr[k*height*width+(height-2)*width+i] = true;*/
				occupyPtr[k*height*width+(height-1)*width+i] = true;
			}
		}

		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				occupyPtr[0*height*width+j*width+i] = true;
				/*occupyPtr[1*height*width+j*width+i] = true;
				occupyPtr[2*height*width+j*width+i] = true;
				occupyPtr[3*height*width+j*width+i] = true;
				occupyPtr[(depth-4)*height*width+j*width+i] = true;
				occupyPtr[(depth-3)*height*width+j*width+i] = true;
				occupyPtr[(depth-2)*height*width+j*width+i] = true;*/
				occupyPtr[(depth-1)*height*width+j*width+i] = true;
			}
		}
	}
	
	switch(para.intobj)
	{
	case ZQ_SmokeSimulationPara3D::INT_OBJ_SPHERE:
		{
			for(int k = 0;k < depth;k++)
			{
				for(int j = 0;j < height;j++)
				{
					for(int i = 0;i < width;i++)
					{
						if((k+0.5-depth/2)*(k+0.5-depth/2) + (j+0.5-height/2)*(j+0.5-height/2) + (i+0.5-width/2)*(i+0.5-width/2) < width*width/400.0)
							occupyPtr[k*height*width+j*width+i] = true;
					}
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara3D::INT_OBJ_BOX:
		{
			for(int k = 0;k < depth;k++)
			{
				for(int j = 0;j < height;j++)
				{
					for(int i = 0;i < width;i++)
					{
						if(fabs(k+0.5-depth/2) < depth/20 && fabs(j+0.5-height/2) < height/20 && fabs(i+0.5-width/2) < width/10)
							occupyPtr[k*height*width+j*width+i] = true;
					}
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara3D::INT_OBJ_OVAL:
		{
			for(int k = 0;k < depth;k++)
			{
				for(int i = 0;i < height;i++)
				{
					for(int j = 0;j < width;j++)
					{
						if(sqrt((i+0.5-3.0*width/7)*(i+0.5-3.0*width/7) + (j+0.5-height/2)*(j+0.5-height/2) + (k+0.5-depth/2)*(k+0.5-depth/2))
							+ sqrt((i+0.5-4.0*width/7)*(i+0.5-4.0*width/7) + (j+0.5-height/2)*(j+0.5-height/2) + (k+0.5-depth/2)*(k+0.5-depth/2))
							< width / 5.0)
						{
							occupyPtr[k*height*width+j*width+i] = true;
						}
					}
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara3D::INT_OBJ_SPECIAL1:
		{
			for(int k = 0;k < depth;k++)
			{
				for(int j = 0;j < height;j++)
				{
					for(int i = 0;i < width;i++)
					{
						if((j == height/2+13 || j == height/2+14) && (i%4 == 0 || i%4 == 3) && i >=4 && i < width-4)
							occupyPtr[k*height*width+j*width+i] = true;
					}
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara3D::INT_OBJ_NULL: default:
		break;
	}
	
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	float*& Wptr = fineGrid->GetVelocityWptr();
	memset(Uptr,0,sizeof(float)*(width+1)*height*depth);
	memset(Vptr,0,sizeof(float)*width*(height+1)*depth);
	memset(Wptr,0,sizeof(float)*width*height*(depth+1));

	return true;
}

bool ZQ_InitCoarseGrid3D(const ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid)
{
	if(fineGrid == 0 || coarseGrid == 0)
		return false;

	int width = fineGrid->GetWidth();
	int height = fineGrid->GetHeight();
	int depth = fineGrid->GetDepth();

	int coarseWidth = coarseGrid->GetWidth();
	int coarseHeight = coarseGrid->GetHeight();
	int coarseDepth = coarseGrid->GetDepth();

	int subSize = width/coarseWidth;
	if(subSize*coarseWidth != width || subSize*coarseHeight != height || subSize*coarseDepth != depth)
		return false;

	coarseGrid->BuildFaceUnOccupyRatio();
	
	const bool*& fine_occupyPtr = fineGrid->GetOccupyPtr();
	float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();

	int subFaceSize = subSize*subSize;
	for(int k = 0;k < coarseDepth;k++)
	{
		for(int j = 0;j < coarseHeight;j++)
		{
			for(int i = 1; i < coarseWidth;i++)
			{
				int occupyNum = 0;
				int i_shift = i*subSize;
				int j_shift = j*subSize;
				int k_shift = k*subSize;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						if(fine_occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift] 
						|| fine_occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift-1])
							occupyNum ++;
					}					
				}
				unoccupyUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+i] = 1.0f - (float)occupyNum/subFaceSize;
			}
			
			{
				int k_shift = k*subSize;
				int j_shift = j*subSize;
				int occupyNum = 0;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						if(fine_occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+0])
							occupyNum ++;
					}
				}
				unoccupyUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+0] = 1.0f - (float)occupyNum/subFaceSize;
				
				occupyNum = 0;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						if(fine_occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+width-1])
							occupyNum ++;
					}
				}
				unoccupyUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+coarseWidth] = 1.0f - (float)occupyNum/subFaceSize;
			}
		}
	}
	

	for(int k = 0;k < coarseDepth;k++)
	{
		for(int i = 0;i < coarseWidth;i++)
		{
			for(int j = 1;j < coarseHeight;j++)
			{
				int occupyNum = 0;
				int k_shift = k*subSize;
				int i_shift = i*subSize;
				int j_shift = j*subSize;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(fine_occupyPtr[(k_shift+sk)*height*width+j_shift*width+i_shift+si]
						|| fine_occupyPtr[(k_shift+sk)*height*width+(j_shift-1)*width+i_shift+si])
							occupyNum ++;
					}
				}
				unoccupyVptr[k*(coarseHeight+1)*coarseWidth+j*coarseWidth+i] = 1.0f - (float)occupyNum/subFaceSize;
			}
			{
				int k_shift = k*subSize;
				int i_shift = i*subSize;
				int occupyNum = 0;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(fine_occupyPtr[(k_shift+sk)*height*width+0*width+i_shift+si])
							occupyNum ++;
					}
				}
				unoccupyVptr[k*(coarseHeight+1)*coarseWidth+0*coarseWidth+i] = 1.0f - (float)occupyNum/subFaceSize;
				occupyNum = 0;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(fine_occupyPtr[(k_shift+sk)*height*width+(height-1)*width+i_shift+si])
							occupyNum ++;
					}
				}
				unoccupyVptr[k*(coarseHeight+1)*coarseWidth+coarseHeight*coarseWidth+i] = 1.0f - (float)occupyNum/subFaceSize;
			}
		}
	}

	for(int j = 0;j < coarseHeight;j++)
	{
		for(int i = 0;i < coarseWidth;i++)
		{
			for(int k = 1;k < coarseDepth;k++)
			{
				int occupyNum = 0;
				int j_shift = j*subSize;
				int i_shift = i*subSize;
				int k_shift = k*subSize;
				for(int sj = 0;sj < subSize;sj++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(fine_occupyPtr[k_shift*height*width+(j_shift+sj)*width+i_shift+si]
						||fine_occupyPtr[(k_shift-1)*height*width+(j_shift+sj)*width+i_shift+si])
							occupyNum ++;
					}
				}
				unoccupyWptr[k*coarseHeight*coarseWidth+j*coarseWidth+i] = 1.0f - (float)occupyNum/subFaceSize;
			}
			{
				int j_shift = j*subSize;
				int i_shift = i*subSize;
				int occupyNum = 0;
				for(int sj = 0;sj < subSize;sj++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(fine_occupyPtr[(j_shift+sj)*width+i_shift+si])
							occupyNum ++;
					}
				}
				unoccupyWptr[j*coarseWidth+i] = 1.0f - (float)occupyNum/subFaceSize;
				occupyNum = 0;
				for(int sj = 0;sj < subSize;sj++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(fine_occupyPtr[(depth-1)*height*width+(j_shift+sj)*width+i_shift+si])
							occupyNum ++;
					}
				}
				unoccupyWptr[coarseDepth*coarseHeight*coarseWidth+j*coarseWidth+i] = 1.0f - (float)occupyNum/subFaceSize;
			}
		}
	}
	
	bool*& coarse_occupyPtr = coarseGrid->GetOccupyPtr();
	for(int k = 0;k < coarseDepth;k++)
	{
		for(int j = 0;j < coarseHeight;j++)
		{
			for(int i = 0;i < coarseWidth;i++)
			{
				if(unoccupyUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+i] == 0
					&& unoccupyUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+i+1] == 0
					&& unoccupyVptr[k*(coarseHeight+1)*coarseWidth+j*coarseWidth+i] == 0
					&& unoccupyVptr[k*(coarseHeight+1)*coarseWidth+(j+1)*coarseWidth+i] == 0
					&& unoccupyWptr[k*coarseHeight*coarseWidth+j*coarseWidth+i] == 0
					&& unoccupyWptr[(k+1)*coarseHeight*coarseWidth+j*coarseWidth+i] == 0)
				{
					coarse_occupyPtr[k*coarseHeight*coarseWidth+j*coarseWidth+i] = true;
				}
				else
					coarse_occupyPtr[k*coarseHeight*coarseWidth+j*coarseWidth+i] = false;
			}
		}
	}
	
	
	coarseGrid->BuildVolumeUnOccupyRatio();
	float*& unoccupyVolumePtr = coarseGrid->GetVolumeUnOccupyRatioPtr();
	for(int k = 0;k < coarseDepth;k++)
	{
		for(int j = 0;j < coarseHeight;j++)
		{
			for(int i = 0;i < coarseWidth;i++)
			{
				int occupyNum = 0;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						for(int si = 0;si < subSize;si++)
						{
							if(fine_occupyPtr[(k*subSize+sk)*height*width+(j*subSize+sj)*width+(i*subSize+si)])
								occupyNum ++;
						}
					}
				}
				
				unoccupyVolumePtr[k*coarseHeight*coarseWidth+j*coarseWidth+i] = 1.0f - (float)occupyNum/(subSize*subSize*subSize);
			}
		}
	}
	
	return true;
}

bool ZQ_InitSimulation3D(ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, taucs_ccs_matrix** coarseA,
						 float** subMatrix, int* subRow, int* subCol)
{
	if(!ZQ_InitFineGrid3D(fineGrid,para))
		return false;
	
	if(!ZQ_InitCoarseGrid3D(fineGrid,coarseGrid))
	{
		return false;
	}

	int subSize = para.subSize;
	
	int matrixsz = 0;
	if(para.globalGridSolver == OPEN_POISSON)
	{
		if((matrixsz = ZQ_InitMatOpenPoisson3D(coarseGrid,coarseA)) == 0)
			return false;
	}
	else if(para.globalGridSolver == CLOSED_POISSON)
	{
		if((matrixsz =  ZQ_InitMatClosedPoisson3D(coarseGrid,coarseA)) == 0)
			return false;
	}
	else if(para.globalGridSolver == OPEN_FLUX)
	{
		if((matrixsz = ZQ_InitMatOpenFlux3D(coarseGrid,coarseA)) == 0)
			return false;
	}
	else if(para.globalGridSolver == CLOSED_FLUX)
	{
		if((matrixsz = ZQ_InitMatClosedFlux3D(coarseGrid,coarseA)) == 0)
			return false;
	}

	if(para.subSize >= 2)
	{
		if(!ZQ_InitSubMatPoisson3D(subSize, subMatrix,subRow,subCol))
			return false;
	}
	return true;
}

bool ZQ_ApplyMovingObject3D(ZQ_SmokeGrid3D* fineGrid, const std::vector<ZQ_SmokeMovingObject3D*> mvobjs, const ZQ_SmokeSimulationPara3D& para)
{
	int count = mvobjs.size();
	int width = para.width;
	int height = para.height;
	int depth = para.depth;
	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	float deltah = para.voxelSize;
	float deltat = para.stepTime;
	float*& u = fineGrid->GetVelocityUptr();
	float*& v = fineGrid->GetVelocityVptr();
	float*& w = fineGrid->GetVelocityWptr();
	float*& density = fineGrid->GetDensityPtr();
	float*& temperature = fineGrid->GetTemperaturePtr();
	memset(occupyPtr,false,sizeof(bool)*width*height*depth);

	for(int i = 0;i < count;i++)
	{
		int cx,cy,cz;
		int posx,posy,posz;
		int sizex,sizey,sizez;
		int curu,curv,curw;
		bool* objoccupy = 0;
		if(mvobjs[i] && (objoccupy = (bool*)(mvobjs[i]->GetOccupyPtr())) != 0)
		{
			mvobjs[i]->GetCenter(cx,cy,cz);
			mvobjs[i]->GetPos(posx,posy,posz);
			mvobjs[i]->GetSize(sizex,sizey,sizez);
			mvobjs[i]->GetCurVel(curu,curv,curw);
			double scale = deltah / deltat;
			curu *= scale;
			curv *= scale;
			curw *= scale;

			for(int sk = 0;sk < sizez;sk++)
			{
				for(int sj = 0;sj < sizey;sj++)
				{
					for(int si = 0;si < sizex;si++)
					{
						if(objoccupy[sk*sizey*sizex+sj*sizex+si])
						{
							int ix = si+posx-cx;
							int iy = sj+posy-cy;
							int iz = sk+posz-cz;
							if(ix >= 0 && ix < width && iy >= 0 && iy < height && iz >= 0 && iz < depth)
							{
								occupyPtr[iz*height*width+iy*width+ix] = true;
								u[iz*height*(width+1)+iy*(width+1)+ix] = curu;
								u[iz*height*(width+1)+iy*(width+1)+(ix+1)] = curu;
								v[iz*(height+1)*width+iy*width+ix] = curv;
								v[iz*(height+1)*width+(iy+1)*width+ix] = curv;
								w[iz*height*width+iy*width+ix] = curw;
								w[(iz+1)*height*width+iy*width+ix] = curw;
							}					
						}
					}
				}
			}
		}
	}

	if(para.globalGridSolver == CLOSED_POISSON || para.globalGridSolver == CLOSED_FLUX 
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON || para.globalGridSolver == CLOSED_OCTREE_FLUX
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON_SOR
		|| para.globalGridSolver == TOTAL_CLOSED_POISSON_SOR || para.globalGridSolver == TOTAL_CLOSED_FLUX_SOR
		|| para.globalGridSolver == COARSE_CLOSED_POISSON_SOR || para.globalGridSolver == COARSE_CLOSED_FLUX_SOR)
	{
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				occupyPtr[k*height*width+j*width+0] = true;
			/*	occupyPtr[k*height*width+j*width+1] = true;
				occupyPtr[k*height*width+j*width+2] = true;
				occupyPtr[k*height*width+j*width+3] = true;
				occupyPtr[k*height*width+j*width+width-4] = true;
				occupyPtr[k*height*width+j*width+width-3] = true;
				occupyPtr[k*height*width+j*width+width-2] = true;*/
				occupyPtr[k*height*width+j*width+width-1] = true;
			}
		}
		for(int k = 0;k < depth;k++)
		{
			for(int i = 0;i < width;i++)
			{
				occupyPtr[k*height*width+0*width+i] = true;
				/*occupyPtr[k*height*width+1*width+i] = true;
				occupyPtr[k*height*width+2*width+i] = true;
				occupyPtr[k*height*width+3*width+i] = true;
				occupyPtr[k*height*width+(height-4)*width+i] = true;
				occupyPtr[k*height*width+(height-3)*width+i] = true;
				occupyPtr[k*height*width+(height-2)*width+i] = true;*/
				occupyPtr[k*height*width+(height-1)*width+i] = true;
			}
		}

		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				occupyPtr[0*height*width+j*width+i] = true;
				/*occupyPtr[1*height*width+j*width+i] = true;
				occupyPtr[2*height*width+j*width+i] = true;
				occupyPtr[3*height*width+j*width+i] = true;
				occupyPtr[(depth-4)*height*width+j*width+i] = true;
				occupyPtr[(depth-3)*height*width+j*width+i] = true;
				occupyPtr[(depth-2)*height*width+j*width+i] = true;*/
				occupyPtr[(depth-1)*height*width+j*width+i] = true;
			}
		}
	}
	return true;
}


bool ZQ_AddForce3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(fineGrid == 0)
		return false;

	ZQ_CUDA_AddForce3D::AddForceType type;
	switch (para.addforce_type)
	{
	case ZQ_SmokeSimulationPara3D::ADD_FORCE_ENTIRE:
		type = ZQ_CUDA_AddForce3D::ADD_FORCE_ENTIRE;
		break;
	case ZQ_SmokeSimulationPara3D::ADD_FORCE_SLICE:
		type = ZQ_CUDA_AddForce3D::ADD_FORCE_SLICE;
		break;
	default:
		return false;
		break;
	}

	int width = fineGrid->GetWidth();
	int height = fineGrid->GetHeight();
	int depth = fineGrid->GetDepth();

	float alpha = para.gravityCoeff;
	float beta  = para.buoyancyCoeff;
	float Tamb  = para.Tamb;
	float deltat = para.stepTime;
	float deltah = para.voxelSize;
	float eta = para.confineCoeff;


	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	float*& Wptr = fineGrid->GetVelocityWptr();
	bool*& occupy = fineGrid->GetOccupyPtr();
	float*& temperature = fineGrid->GetTemperaturePtr();


	float cuda_cost_time =
	ZQ_CUDA_AddForce3D::AddForce3D(Uptr,Vptr,Wptr,occupy,temperature,deltah,deltat,beta,eta,Tamb,width,height,depth,type);

	printf("add force cuda cost = %f\n",0.001*cuda_cost_time);
	return true;
}

bool ZQ_Attenuation3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(fineGrid == 0)
		return false;
	float attenVel = exp(-para.stepTime * para.velAttenCoeff);
	float attenDen = exp(-para.stepTime * para.densityAttenCoeff);
	float attenTemp = exp(-para.stepTime * para.tempAttenCoeff);
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	float*& Wptr = fineGrid->GetVelocityWptr();
	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	float*& densityPtr = fineGrid->GetDensityPtr();
	float*& temperaturePtr = fineGrid->GetTemperaturePtr();
	int width = para.width;
	int height = para.height;
	int depth = para.depth;

	float cuda_cost_time = 
	ZQ_CUDA_Attenuation3D::Attenuation3D(Uptr,Vptr,Wptr,temperaturePtr,densityPtr,occupyPtr,attenVel,attenTemp,attenDen,width,height,depth);

	printf("attenaution cuda cost = %f\n",0.001*cuda_cost_time);

	return true;
}

bool ZQ_AdvectVelocity3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if (fineGrid == 0)
		return false;
	
	ZQ_CUDA_Advection3D::AdvectVelocityType type;
	switch (para.advVelType)
	{
	case ZQ_SmokeSimulationPara3D::ADV_VEL_INMAC_OUTMAC:
		type = ZQ_CUDA_Advection3D::ADV_VEL_INMAC_OUTMAC;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_VEL_INMAC_OUTMAC_BFECC:
		type = ZQ_CUDA_Advection3D::ADV_VEL_INMAC_OUTMAC_BFECC;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_VEL_INREG_OUTREG:
		type = ZQ_CUDA_Advection3D::ADV_VEL_INREG_OUTREG;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_VEL_INREG_OUTREG_BFECC:
		type = ZQ_CUDA_Advection3D::ADV_VEL_INREG_OUTREG_BFECC;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_VEL_INREG_OUTMAC:
		type = ZQ_CUDA_Advection3D::ADV_VEL_INREG_OUTMAC;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_VEL_INREG_OUTMAC_BFECC:
		type = ZQ_CUDA_Advection3D::ADV_VEL_INREG_OUTMAC_BFECC;
		break;
	default:
		return false;
		break;
	}

	unsigned int width = para.width;
	unsigned int height = para.height;
	unsigned int depth = para.depth;
	unsigned int steps = para.steps;
	float voxelSize = para.voxelSize;
	float deltatt = para.stepTime / steps;

	ZQ_CUDA_Advection3D::ZQ_Cuda_Prepare_Advection(width, height, depth, voxelSize, steps, deltatt);

	float*& u = fineGrid->GetVelocityUptr();
	float*& v = fineGrid->GetVelocityVptr();
	float*& w = fineGrid->GetVelocityWptr();
	bool*& occupy = fineGrid->GetOccupyPtr();
	float cuda_cost_time = ZQ_CUDA_Advection3D::Advect_Velocity(u, v, w, occupy, true, type);
	
	printf("velocity advection cuda cost = %f\n", 0.001*cuda_cost_time);

	return true;
}

bool ZQ_MapFine2Coarse3D(const ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(fineGrid == 0 || coarseGrid == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int depth = para.depth;
	int coarseWidth = para.coarseWidth;
	int coarseHeight = para.coarseHeight;
	int coarseDepth = para.coarseDepth;
	int subSize = para.subSize;

	const float*& fine_Uptr = fineGrid->GetVelocityUptr();
	const float*& fine_Vptr = fineGrid->GetVelocityVptr();
	const float*& fine_Wptr = fineGrid->GetVelocityWptr();
	float*& coarse_Uptr = coarseGrid->GetVelocityUptr();
	float*& coarse_Vptr = coarseGrid->GetVelocityVptr();
	float*& coarse_Wptr = coarseGrid->GetVelocityWptr();

	int subFaceSize = subSize*subSize;
	for(int k = 0;k < coarseDepth;k++)
	{
		for(int j = 0;j < coarseHeight;j++)
		{
			for(int i = 0;i <= coarseWidth;i++)
			{
				int i_shift = i*subSize;
				int j_shift = j*subSize;
				int k_shift = k*subSize;
				float subU = 0;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						subU += fine_Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift];
					}
				}
				coarse_Uptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+i] = subU/subFaceSize;
			}
		}
	}
	
	for(int k = 0;k < coarseDepth;k++)
	{
		for(int j = 0;j <= coarseHeight;j++)
		{
			for(int i = 0;i < coarseWidth;i++)
			{
				int i_shift = i*subSize;
				int j_shift = j*subSize;
				int k_shift = k*subSize;
				float subV = 0;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int si = 0;si < subSize;si++)
					{
						subV += fine_Vptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si];
					}
				}
				coarse_Vptr[k*(coarseHeight+1)*coarseWidth+j*coarseWidth+i] = subV/subFaceSize;
			}
		}
	}
	for(int k = 0;k <= coarseDepth;k++)
	{
		for(int j = 0;j < coarseHeight;j++)
		{
			for(int i = 0;i < coarseWidth;i++)
			{
				int i_shift = i*subSize;
				int j_shift = j*subSize;
				int k_shift = k*subSize;
				float subW = 0;
				for(int sj = 0;sj < subSize;sj++)
				{
					for(int si = 0;si < subSize;si++)
					{
						subW += fine_Wptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si];
					}
				}
				coarse_Wptr[k*coarseHeight*coarseWidth+j*coarseWidth+i] = subW/subFaceSize;
			}
		}
	}
	
	
	return true;
}

bool ZQ_SolveOpenPressure3D(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA)
{
	if(coarseGrid == 0)
		return false;
	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetWidth();
	int depth = coarseGrid->GetDepth();
	bool*& occupy = coarseGrid->GetOccupyPtr();
	float*& Uptr = coarseGrid->GetVelocityUptr();
	float*& Vptr = coarseGrid->GetVelocityVptr();
	float*& Wptr = coarseGrid->GetVelocityWptr();
	float*& pressure = coarseGrid->GetPressurePtr();

	int idx = 1;
	int* index = new int[width*height*depth];
	memset(index,0,sizeof(int)*width*height*depth);
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(occupy[k*height*width+j*width+i])
					continue;
				index[k*height*width+j*width+i] = idx++;
			}
		}
	}
	

	int dim = idx-1;
	
	double* divVelocity = new double[dim];
	memset(divVelocity,0,sizeof(double)*dim);
	double* x = new double[dim];
	memset(x,0,sizeof(double)*dim);
	double* x0 = new double[dim];
	memset(x0,0,sizeof(double)*dim);

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int offset = k*height*width+j*width+i;
				int row = index[offset] - 1;
				if(row < 0)
					continue;
				divVelocity[row] = Uptr[k*height*(width+1)+j*(width+1)+i+1] - Uptr[k*height*(width+1)+j*(width+1)+i] 
								+ Vptr[k*(height+1)*width+(j+1)*width+i] - Vptr[k*(height+1)*width+j*width+i]
								+ Wptr[(k+1)*height*width+j*width+i] - Wptr[k*height*width+j*width+i];
			}
		}
	}
	

	
	int max_iter = para.maxIter;
	double tol = 1e-12;
	int it = 0;
	
	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(coarseA, divVelocity, x0, max_iter, tol, x, it, false);
	

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(index[k*height*width+j*width+i] != 0)
				{
					pressure[k*height*width+j*width+i] = x[index[k*height*width+j*width+i]-1];
				}
			}
		}
	}
	
	delete []x0;
	delete []x;
	delete []divVelocity;
	delete []index;

	return true;
}


bool ZQ_SolveOpenPressure3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(fineGrid == 0)
		return false;

	ZQ_CUDA_PoissonSolver3D::SolveOpenPoissonRedBlackwithOccupy3D_MAC(
		fineGrid->GetVelocityUptr(),
		fineGrid->GetVelocityVptr(),
		fineGrid->GetVelocityWptr(),
		fineGrid->GetOccupyPtr(),
		fineGrid->GetWidth(),
		fineGrid->GetHeight(),
		fineGrid->GetDepth(),
		para.maxIter
		);

	return true;
}

bool ZQ_SolveClosedPressure3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(fineGrid == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int depth = para.depth;
	int first_x = -1;
	int first_y = -1;
	int first_z = -1;
	float div_per_volume = 0.0f;
	int count = 0;

	bool*& occupy = fineGrid->GetOccupyPtr();
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	float*& Wptr = fineGrid->GetVelocityWptr();

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(!occupy[k*height*width+j*width+i])
				{
					if(first_x == -1 && first_y == -1 && first_z == -1)
					{
						first_x = i;
						first_y = j;
						first_z = k;
					}
					count ++;

					div_per_volume += Uptr[k*height*(width+1)+j*(width+1)+i+1] - Uptr[k*height*(width+1)+j*(width+1)+i] 
									+ Vptr[k*(height+1)*width+(j+1)*width+i] - Vptr[k*(height+1)*width+j*width+i]
									+ Wptr[(k+1)*height*width+j*width+i] - Wptr[k*height*width+j*width+i];

				}
			}
		}
	}
	

	if(first_x < 0 && first_y < 0 && first_z < 0)
		return false;

	if(count == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}
	div_per_volume /= count;


	ZQ_CUDA_PoissonSolver3D::SolveClosedPoissonRedBlackwithOccupy3D_MAC(
		Uptr,
		Vptr,
		Wptr,
		occupy,
		//first_x,
		//first_y,
		//first_z,
		div_per_volume,
		width,
		height,
		depth,
		para.maxIter
		);
	
	return true;
}

bool ZQ_SolveOpenFlux3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(fineGrid == 0)
		return false;


	ZQ_CUDA_PoissonSolver3D::SolveOpenFluxRedBlackwithOccupy3D_MAC(
		fineGrid->GetVelocityUptr(),
		fineGrid->GetVelocityVptr(),
		fineGrid->GetVelocityWptr(),
		fineGrid->GetOccupyPtr(),
		fineGrid->GetWidth(),
		fineGrid->GetHeight(),
		fineGrid->GetDepth(),
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);

	
	return true;
}


bool ZQ_SolveClosedFlux3DSOR(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(fineGrid == 0)
		return false;
	int width = para.width;
	int height = para.height;
	int depth = para.depth;
	float div_per_volume = 0.0f;

	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	float*& Wptr = fineGrid->GetVelocityWptr();
	bool*& occupy = fineGrid->GetOccupyPtr();

	int count = 0;
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(!occupy[k*height*width+j*width+i])
				{
					count ++;
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

	ZQ_CUDA_PoissonSolver3D::SolveClosedFluxRedBlackwithOccupy3D_MAC(
		Uptr,
		Vptr,
		Wptr,
		occupy,
		div_per_volume,
		width,
		height,
		depth,
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);


	return true;
}

bool ZQ_SolveClosedPressure3D(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA)
{
	if(coarseGrid == 0)
		return false;
	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();
	bool*& occupy = coarseGrid->GetOccupyPtr();
	float*& Uptr = coarseGrid->GetVelocityUptr();
	float*& Vptr = coarseGrid->GetVelocityVptr();
	float*& Wptr = coarseGrid->GetVelocityWptr();
	float*& pressure = coarseGrid->GetPressurePtr();

	int idx = 1;
	int* index = new int[width*height*depth];
	memset(index,0,sizeof(int)*width*height*depth);
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(occupy[k*height*width+j*width+i])
					continue;
				index[k*height*width+j*width+i] = idx++;
			}
		}
	}
	

	int dim = idx-2;
	int rownum = dim+1;

	double* divVelocity = new double[rownum];
	memset(divVelocity,0,sizeof(double)*rownum);
	double* x = new double[dim];
	memset(x,0,sizeof(double)*dim);
	double* x0 = new double[dim];
	memset(x0,0,sizeof(double)*dim);

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int row_id = index[k*height*width+j*width+i] - 1;
				if(row_id < 0)
					continue;
				divVelocity[row_id] = Uptr[k*height*(width+1)+j*(width+1)+i+1] - Uptr[k*height*(width+1)+j*(width+1)+i] 
									+ Vptr[k*(height+1)*width+(j+1)*width+i] - Vptr[k*(height+1)*width+j*width+i]
									+ Wptr[(k+1)*height*width+j*width+i] - Wptr[k*height*width+j*width+i];
			}
		}
	}
	

	int max_iter = para.maxIter;
	double tol = 1e-12;
	int it = 0;
	
	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(coarseA, divVelocity, x0, max_iter, tol, x, it, false);
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int offset = k*height*width+j*width+i;
				if(index[offset] > 0)
				{
					if(index[offset] == 1)
						pressure[offset] = 0;
					else
						pressure[offset] = x[index[offset]-2];
				}
			}
		}
	}
	
	delete []x0;
	delete []x;
	delete []divVelocity;
	delete []index;
	return true;
}

bool ZQ_SolveCoarseOpenPoisson3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	int depth = para.coarseDepth;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	float*& coarseW = coarseGrid->GetVelocityWptr();
	bool*& occupy = coarseGrid->GetOccupyPtr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& unoccupyW = coarseGrid->GetFaceUnOccupyRatioWPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	float*& deltaW = deltaGrid->GetVelocityWptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height*depth);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1)*depth);
	memcpy(deltaW,coarseW,sizeof(float)*width*height*(depth+1));

	ZQ_CUDA_PoissonSolver3D::SolveOpenPoissonRedBlackwithFaceRatio3D_MAC(
		deltaU,deltaV,deltaW,occupy,unoccupyU,unoccupyV,unoccupyW,width,height,depth,para.maxIter
		);

	for(int i = 0;i < (width+1)*height*depth;i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		deltaV[i] -= coarseV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		deltaW[i] -= coarseW[i];

	return true;
}

bool ZQ_SolveCoarseClosedPoisson3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	int depth = para.coarseDepth;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	float*& coarseW = coarseGrid->GetVelocityWptr();
	float*& unoccupyVolume = coarseGrid->GetVolumeUnOccupyRatioPtr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& unoccupyW = coarseGrid->GetFaceUnOccupyRatioWPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	float*& deltaW = deltaGrid->GetVelocityWptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height*depth);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1)*depth);
	memcpy(deltaW,coarseW,sizeof(float)*width*height*(depth+1));

	float div_per_volume = 0.0f;
	float volumeSize = 0.0f;
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				volumeSize += unoccupyVolume[k*height*width+j*width+i];
				div_per_volume += coarseU[k*height*(width+1)+j*(width+1)+i+1] - coarseU[k*height*(width+1)+j*(width+1)+i]
								+ coarseV[k*(height+1)*width+(j+1)*width+i] - coarseV[k*(height+1)*width+j*width+i]
								+ coarseW[(k+1)*height*width+j*width+i] - coarseW[k*height*width+j*width+i];
			}
		}
	}
	
	if(volumeSize == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}
	div_per_volume /= volumeSize;

	ZQ_CUDA_PoissonSolver3D::SolveClosedPoissonRedBlackwithFaceRatio3D_MAC(
		deltaU,deltaV,deltaW,unoccupyVolume,unoccupyU,unoccupyV,unoccupyW,div_per_volume,width,height,depth,para.maxIter
		);

	for(int i = 0;i < (width+1)*height*depth;i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		deltaV[i] -= coarseV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		deltaW[i] -= coarseW[i];

	return true;
}

bool ZQ_SolveCoarseOpenFlux3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	int depth = para.coarseDepth;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	float*& coarseW = coarseGrid->GetVelocityWptr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& unoccupyW = coarseGrid->GetFaceUnOccupyRatioWPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	float*& deltaW = deltaGrid->GetVelocityWptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height*depth);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1)*depth);
	memcpy(deltaW,coarseW,sizeof(float)*width*height*(depth+1));

	ZQ_CUDA_PoissonSolver3D::SolveOpenFluxRedBlackwithFaceRatio3D_MAC(
		deltaU,deltaV,deltaW,unoccupyU,unoccupyV,unoccupyW,width,height,depth,
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);

	for(int i = 0;i < (width+1)*height*depth;i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		deltaV[i] -= coarseV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		deltaW[i] -= coarseW[i];

	return true;
}

bool ZQ_SolveCoarseClosedFlux3DSOR(ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	int depth = para.coarseDepth;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	float*& coarseW = coarseGrid->GetVelocityWptr();
	float*& unoccupyVolume = coarseGrid->GetVolumeUnOccupyRatioPtr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& unoccupyW = coarseGrid->GetFaceUnOccupyRatioWPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	float*& deltaW = deltaGrid->GetVelocityWptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height*depth);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1)*depth);
	memcpy(deltaW,coarseW,sizeof(float)*width*height*(depth+1));

	float div_per_volume = 0.0f;
	float volume = 0.0f;

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				volume += unoccupyVolume[k*height*width+j*width+i];
				div_per_volume += coarseU[k*height*(width+1)+j*(width+1)+i+1] - coarseU[k*height*(width+1)+j*(width+1)+i]
								+ coarseV[k*(height+1)*width+(j+1)*width+i] - coarseV[k*(height+1)*width+j*width+i]
								+ coarseW[(k+1)*height*width+j*width+i] - coarseW[k*height*width+j*width+i];
			}
		}
	}
	

	if(volume == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}
	div_per_volume /= volume;

	ZQ_CUDA_PoissonSolver3D::SolveClosedFluxRedBlackwithFaceRatio3D_MAC(
		deltaU,deltaV,deltaW,unoccupyVolume,unoccupyU,unoccupyV,unoccupyW,div_per_volume,width,height,depth,
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);

	for(int i = 0;i < (width+1)*height*depth;i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1)*depth;i++)
		deltaV[i] -= coarseV[i];
	for(int i = 0;i < width*height*(depth+1);i++)
		deltaW[i] -= coarseW[i];

	return true;
}

bool ZQ_AdjustOpenVelocity3D(const ZQ_SmokeGrid3D* coarseGrid, ZQ_SmokeGrid3D* deltaCoarseGrid)
{
	if(coarseGrid == 0 || deltaCoarseGrid == 0)
		return false;

	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();

	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();
	const float*& pressure = coarseGrid->GetPressurePtr();
	float*& delta_Uptr = deltaCoarseGrid->GetVelocityUptr();
	float*& delta_Vptr = deltaCoarseGrid->GetVelocityVptr();
	float*& delta_Wptr = deltaCoarseGrid->GetVelocityWptr();

	memset(delta_Uptr,0,sizeof(float)*(width+1)*height*depth);
	memset(delta_Vptr,0,sizeof(float)*width*(height+1)*depth);
	memset(delta_Wptr,0,sizeof(float)*width*height*(depth+1));

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 1;i < width;i++)
			{
				if(unoccupyUptr[k*height*(width+1)+j*(width+1)+i] > 0)
					delta_Uptr[k*height*(width+1)+j*(width+1)+i] = -(pressure[k*height*width+j*width+i]-pressure[k*height*width+j*width+i-1]);
			}
			if(unoccupyUptr[k*height*(width+1)+j*(width+1)+0] > 0)
				delta_Uptr[k*height*(width+1)+j*(width+1)+0] = -(pressure[k*height*width+j*width+0]-0);
			if(unoccupyUptr[k*height*width+j*(width+1)+width] > 0)
				delta_Uptr[k*height*(width+1)+j*(width+1)+width] = -(0-pressure[k*height*width+j*width+width-1]);
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int i = 0;i < width;i++)
		{
			for(int j = 1;j < height;j++)
			{
				if(unoccupyVptr[k*(height+1)*width+j*width+i] > 0)
					delta_Vptr[k*(height+1)*width+j*width+i] = -(pressure[k*height*width+j*width+i]-pressure[k*height*width+(j-1)*width+i]);
			}
			if(unoccupyVptr[k*(height+1)*width+0*width+i] > 0)
				delta_Vptr[k*(height+1)*width+0*width+i] = -(pressure[k*height*width+0*width+i]-0);
			if(unoccupyVptr[k*(height+1)*width+height*width+i] > 0)
				delta_Vptr[k*(height+1)*width+height*width+i] = -(0-pressure[k*height*width+(height-1)*width+i]);
		}
	}

	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			for(int k = 1;k < depth;k++)
			{
				if(unoccupyWptr[k*height*width+j*width+i] > 0)
					delta_Wptr[k*height*width+j*width+i] = -(pressure[k*height*width+j*width+i]-pressure[(k-1)*height*width+j*width+i]);
			}
			if(unoccupyWptr[j*width+i] > 0)
				delta_Wptr[j*width+i] = -(pressure[j*width+i]-0);
			if(unoccupyWptr[depth*height*width+j*width+i] > 0)
				delta_Wptr[depth*height*width+j*width+i] = -(0-pressure[(depth-1)*height*width+j*width+i]);
		}
	}
	

	return true;
}

bool ZQ_AdjustClosedVelocity3D(const ZQ_SmokeGrid3D* coarseGrid, ZQ_SmokeGrid3D* deltaCoarseGrid)
{
	if(coarseGrid == 0 || deltaCoarseGrid == 0)
		return false;

	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();

	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();
	const float*& pressure = coarseGrid->GetPressurePtr();
	float*& delta_Uptr = deltaCoarseGrid->GetVelocityUptr();
	float*& delta_Vptr = deltaCoarseGrid->GetVelocityVptr();
	float*& delta_Wptr = deltaCoarseGrid->GetVelocityWptr();

	memset(delta_Uptr,0,sizeof(float)*(width+1)*height*depth);
	memset(delta_Vptr,0,sizeof(float)*width*(height+1)*depth);
	memset(delta_Wptr,0,sizeof(float)*width*height*(depth+1));

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 1;i < width;i++)
			{
				if(unoccupyUptr[k*height*(width+1)+j*(width+1)+i] > 0)
					delta_Uptr[k*height*(width+1)+j*(width+1)+i] = -(pressure[k*height*width+j*width+i]-pressure[k*height*width+j*width+i-1]);
			}
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int i = 0;i < width;i++)
		{
			for(int j = 1;j < height;j++)
			{
				if(unoccupyVptr[k*(height+1)*width+j*width+i] > 0)
					delta_Vptr[k*(height+1)*width+j*width+i] = -(pressure[k*height*width+j*width+i]-pressure[k*height*width+(j-1)*width+i]);
			}
		}
	}

	for(int j = 0;j < height;j++)
	{
		for(int i = 0;i < width;i++)
		{
			for(int k = 1;k < depth;k++)
			{
				if(unoccupyWptr[k*height*width+j*width+i] > 0)
					delta_Wptr[k*height*width+j*width+i] = -(pressure[k*height*width+j*width+i]-pressure[(k-1)*height*width+j*width+i]);
			}
		}
	}
	

	return true;
}

bool ZQ_SolveFlux3D(const ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA, ZQ_SmokeGrid3D* deltaCoarseGrid)
{
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();

	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyWptr = coarseGrid->GetFaceUnOccupyRatioWPtr();
	const bool*& occupy = coarseGrid->GetOccupyPtr();
	const float*& Uptr = coarseGrid->GetVelocityUptr();
	const float*& Vptr = coarseGrid->GetVelocityVptr();
	const float*& Wptr = coarseGrid->GetVelocityWptr();
	float*& delta_Uptr = deltaCoarseGrid->GetVelocityUptr();
	float*& delta_Vptr = deltaCoarseGrid->GetVelocityVptr();
	float*& delta_Wptr = deltaCoarseGrid->GetVelocityWptr();

	int* u_index = new int[(width+1)*height*depth];
	int* v_index = new int[width*(height+1)*depth];
	int* w_index = new int[width*height*(depth+1)];
	int* lamda_index = new int[width*height*depth];
	int k_index = 0;

	for(int i = 0;i < (width+1)*height*depth;i++)
		u_index[i] = -1;
	for(int i = 0;i < width*(height+1)*depth;i++)
		v_index[i] = -1;
	for(int i = 0;i < width*height*(depth+1);i++)
		w_index[i] = -1;
	for(int i = 0;i < width*height*depth;i++)
		lamda_index[i] = -1;


	int idx = 0;
	if(para.globalGridSolver == OPEN_FLUX)
	{
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i <= width;i++)
				{
					if(unoccupyUptr[k*height*(width+1)+j*(width+1)+i] != 0)
						u_index[k*height*(width+1)+j*(width+1)+i] = idx++;
				}
			}
		}
		
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j <= height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(unoccupyVptr[k*(height+1)*width+j*width+i] != 0)
						v_index[k*(height+1)*width+j*width+i] = idx++;
				}
			}
		}

		for(int k = 0;k <= depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(unoccupyWptr[k*height*width+j*width+i] != 0)
						w_index[k*height*width+j*width+i] = idx++;
				}
			}
		}
		
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(!occupy[k*height*width+j*width+i])
						lamda_index[k*height*width+j*width+i] = idx++;
				}
			}
		}
	}
	else
	{
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 1;i < width;i++)
				{
					if(unoccupyUptr[k*height*(width+1)+j*(width+1)+i] != 0)
						u_index[k*height*(width+1)+j*(width+1)+i] = idx++;
				}
			}
		}
		
		for(int k = 0;k < depth;k++)
		{
			for(int j = 1;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(unoccupyVptr[k*(height+1)*width+j*width+i] != 0)
						v_index[k*(height+1)*width+j*width+i] = idx++;
				}
			}
		}
		
		for(int k = 1;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(unoccupyWptr[k*height*width+j*width+i] != 0)
						w_index[k*height*width+j*width+i] = idx++;
				}
			}
		}

		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(!occupy[k*height*width+j*width+i])
						lamda_index[k*height*width+j*width+i] = idx++;
				}
			}
		}
		
		k_index = idx++;
	}
	
	int dim = idx;
	double* x0 = new double[dim];
	double* x = new double[dim];
	double* b = new double[dim];
	memset(x0,0,sizeof(double)*dim);
	memset(x,0,sizeof(double)*dim);
	memset(b,0,sizeof(double)*dim);
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				if(lamda_index[k*height*width+j*width+i] >= 0)
				{
					b[lamda_index[k*height*width+j*width+i]] = Uptr[k*height*(width+1)+j*(width+1)+i+1] - Uptr[k*height*(width+1)+j*(width+1)+i] 
															+ Vptr[k*(height+1)*width+(j+1)*width+i] - Vptr[k*(height+1)*width+j*width+i]
															+ Wptr[(k+1)*height*width+j*width+i] - Wptr[k*height*width+j*width+i];
				}
			}
		}
	}
	
	int it = 0;
	double tol = 1e-14;
	ZQ::ZQ_PCGSolver::PCG_sparse_unsquare(coarseA, b, x0, para.maxIter, tol, x, it, false);

	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i <= width;i++)
			{
				int offset = k*height*(width+1)+j*(width+1)+i;
				delta_Uptr[offset] = u_index[offset] >= 0 ? x[u_index[offset]] : 0;
			}
		}
	}
	
	for(int k = 0;k < depth;k++)
	{
		for(int j = 0;j <= height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int offset = k*(height+1)*width+i*width+j;
				delta_Vptr[offset] = v_index[offset] >= 0 ? x[v_index[offset]] : 0;
			}
		}
	}

	for(int k = 0;k <= depth;k++)
	{
		for(int j = 0;j < height;j++)
		{
			for(int i = 0;i < width;i++)
			{
				int offset = k*height*width+j*width+i;
				delta_Wptr[offset] = w_index[offset] >= 0 ? x[w_index[offset]] : 0;
			}
		}
	}
	
	
	delete []u_index;
	delete []v_index;
	delete []w_index;
	delete []lamda_index;
	delete []x0;
	delete []x;
	delete []b;
	return true;
}

bool ZQ_ApplyGuidingVelocity3D(const ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeGuideGrid3D* guideGrid, const ZQ_SmokeSimulationPara3D& para, ZQ_SmokeGrid3D* deltaCoarseGrid)
{
	if( coarseGrid == 0 || guideGrid == 0 || deltaCoarseGrid == 0)
		return false;
	if(!para.useGuidingControl)
		return false;

	float guideCoeff = para.guidingCoeff;
	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	int depth = coarseGrid->GetDepth();

	const float*& Uptr = coarseGrid->GetVelocityUptr();
	const float*& Vptr = coarseGrid->GetVelocityVptr();
	const float*& Wptr = coarseGrid->GetVelocityWptr();
	const float*& GUptr = guideGrid->GetUptr();
	const float*& GVptr = guideGrid->GetVptr();
	const float*& GWptr = guideGrid->GetWptr();
	float*& deltaU = deltaCoarseGrid->GetVelocityUptr();
	float*& deltaV = deltaCoarseGrid->GetVelocityVptr();
	float*& deltaW = deltaCoarseGrid->GetVelocityWptr();

	for(int i = 0;i < (width+1)*height*depth;i++)
		deltaU[i] += guideCoeff*(GUptr[i] - Uptr[i]);
	for(int i = 0;i < width*(height+1)*depth;i++)
		deltaV[i] += guideCoeff*(GVptr[i] - Vptr[i]);
	for(int i = 0;i < width*height*(depth+1);i++)
		deltaW[i] += guideCoeff*(GWptr[i] - Wptr[i]);

	return true;
}

bool ZQ_MapCoarse2Fine3D(const ZQ_SmokeGrid3D* deltaCoarseGrid, ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if(deltaCoarseGrid == 0 || fineGrid == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int depth = para.depth;
	int subSize = para.subSize;
	int coarseWidth = para.coarseWidth;
	int coarseHeight = para.coarseHeight;
	int coarseDepth = para.coarseDepth;
	
	
	/*adjust the fineGrid velocity by plus the deltaCoarseGrid*/
	float*& fineUptr = fineGrid->GetVelocityUptr();
	float*& fineVptr = fineGrid->GetVelocityVptr();
	float*& fineWptr = fineGrid->GetVelocityWptr();
	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	const float*& deltaUptr = deltaCoarseGrid->GetVelocityUptr();
	const float*& deltaVptr = deltaCoarseGrid->GetVelocityVptr();
	const float*& deltaWptr = deltaCoarseGrid->GetVelocityWptr();

	for(int k = 0;k < coarseDepth;k++)
	{
		for(int j = 0;j < coarseHeight;j++)
		{
			int k_shift = k*subSize;
			int j_shift = j*subSize;
			for(int i = 1;i < coarseWidth;i++)
			{
				int i_shift = i*subSize;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						if(!occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift] && !occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift-1])
							fineUptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift] += deltaUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+i];
					}
				}
			}
			int i_shift0 = 0;
			int i_shiftN = coarseWidth*subSize;
			for(int sk = 0;sk < subSize;sk++)
			{
				for(int sj = 0;sj < subSize;sj++)
				{
					if(!occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shift0])
						fineUptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shift0] += deltaUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+0];
					if(!occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+i_shiftN])
						fineUptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+i_shiftN] += deltaUptr[k*coarseHeight*(coarseWidth+1)+j*(coarseWidth+1)+coarseWidth];
				}
			}
		}
	}
	
	for(int k = 0;k < coarseDepth;k++)
	{
		for(int i = 0;i < coarseWidth;i++)
		{
			int k_shift = k*subSize;
			int i_shift = i*subSize;
			for(int j = 1;j < coarseHeight;j++)
			{
				int j_shift = j*subSize;
				for(int sk = 0;sk < subSize;sk++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(!occupyPtr[(k_shift+sk)*height*width+j_shift*width+i_shift+si] && !occupyPtr[(k_shift+sk)*height*width+(j_shift-1)*width+i_shift+si])
							fineVptr[(k_shift+sk)*(height+1)*width+j_shift*width+i_shift+si] += deltaVptr[k*(coarseHeight+1)*coarseWidth+j*coarseWidth+i];
					}
				}
			}
			int j_shift0 = 0;
			int j_shiftN = coarseHeight*subSize;
			for(int sk = 0;sk < subSize;sk++)
			{
				for(int si = 0;si < subSize;si++)
				{
					if(!occupyPtr[(k_shift+sk)*height*width+j_shift0*width+i_shift+si])
						fineVptr[(k_shift+sk)*(height+1)*width+j_shift0*width+i_shift+si] += deltaVptr[k*(coarseHeight+1)*coarseWidth+i];
					if(!occupyPtr[(k_shift+sk)*height*width+(j_shiftN-1)*width+i_shift+si])
						fineVptr[(k_shift+sk)*(height+1)*width+j_shiftN*width+i_shift+si] += deltaVptr[k*(coarseHeight+1)*coarseWidth+(coarseHeight-1)*coarseWidth+i];
				}
			}
		}
	}

	for(int j = 0;j < coarseHeight;j++)
	{
		for(int i = 0;i < coarseWidth;i++)
		{
			int j_shift = j*subSize;
			int i_shift = i*subSize;
			for(int k = 1;k < coarseDepth;k++)
			{
				int k_shift = k*subSize;
				for(int sj = 0;sj < subSize;sj++)
				{
					for(int si = 0;si < subSize;si++)
					{
						if(!occupyPtr[k_shift*height*width+(j_shift+sj)*width+i_shift+si] && !occupyPtr[(k_shift-1)*height*width+(j_shift+sj)*width+i_shift+si])
							fineWptr[k_shift*height*width+(j_shift+sj)*width+i_shift+si] += deltaWptr[k*coarseHeight*coarseWidth+j*coarseWidth+i];
					}
				}
			}
			int k_shift0 = 0;
			int k_shiftN = coarseDepth*subSize;
			for(int sj = 0;sj < subSize;sj++)
			{
				for(int si = 0;si < subSize;si++)
				{
					if(!occupyPtr[k_shift0*height*width+(j_shift+sj)*width+i_shift+si])
						fineWptr[k_shift0*height*width+(j_shift+sj)*width+i_shift+si] += deltaWptr[j*coarseWidth+i];
					if(!occupyPtr[(k_shiftN-1)*height*width+(j_shift+sj)*width+i_shift+si])
						fineWptr[k_shiftN*height*width+(j_shift+sj)*width+i_shift+si] += deltaWptr[coarseDepth*coarseHeight*coarseWidth+j*coarseWidth+i];
				}
			}
		}
	}
	

	return true;
}

bool ZQ_SubProjection3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para, const float* subMatrix, int subRow, int subCol, bool useGPUflag)
{
	if(fineGrid == 0 || subMatrix == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int depth = para.depth;
	int subSize = para.subSize;
	int coarseWidth = para.coarseWidth;
	int coarseHeight = para.coarseHeight;
	int coarseDepth = para.coarseDepth;

	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	float*& Wptr = fineGrid->GetVelocityWptr();
	if(useGPUflag)
	{
		bool* regular_flag = new bool[coarseWidth*coarseHeight*coarseDepth];
		memset(regular_flag,0,sizeof(bool)*coarseWidth*coarseHeight*coarseDepth);
		int fast_count = 0;
		
		for(int k = 0;k < coarseDepth;k++)
		{
			for(int j = 0;j < coarseHeight;j++)
			{
				for(int i = 0;i < coarseWidth;i++)
				{
					for(int sk = 0;sk < subSize;sk++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							for(int si = 0;si < subSize;si++)
							{
								if(occupyPtr[(k*subSize+sk)*height*width+(j*subSize+sj)*width+(i*subSize+si)])
									regular_flag[k*coarseHeight*coarseWidth+j*coarseWidth+i] = true;
							}
						}
					}
					
					if(!regular_flag[k*coarseHeight*coarseWidth+j*coarseWidth+i])
						fast_count++;
				}
			}
		}
		
		if(fast_count > 0)
		{
			float* input = new float[subCol*fast_count];
			memset(input,0,sizeof(float)*subCol*fast_count);
			float* output = new float[subRow*fast_count];
			memset(output,0,sizeof(float)*subRow*fast_count);
			int cur_count = 0;
			for(int k = 0;k < coarseDepth;k++)
			{
				for(int j = 0;j < coarseHeight;j++)
				{
					for(int i = 0;i < coarseWidth;i++)
					{
						if(!regular_flag[k*coarseHeight*coarseWidth+j*coarseWidth+i])
						{
							int cur_row = 0;
							for(int sk = 0;sk < subSize;sk++)
							{
								for(int sj = 0;sj < subSize;sj++)
								{
									for(int si = 0;si <= subSize;si++)
									{
										input[cur_row*fast_count+cur_count] = Uptr[(k*subSize+sk)*height*(width+1)+(j*subSize+sj)*(width+1)+(i*subSize+si)];
										cur_row ++;
									}
								}
							}
							
							for(int sk = 0;sk < subSize;sk++)
							{
								for(int sj = 0;sj <= subSize;sj++)
								{
									for(int si = 0;si < subSize;si++)
									{
										input[cur_row*fast_count+cur_count] = Vptr[(k*subSize+sk)*(height+1)*width+(j*subSize+sj)*width+(i*subSize+si)];
										cur_row ++;
									}
								}
							}

							for(int sk = 0;sk <= subSize;sk++)
							{
								for(int sj = 0;sj < subSize;sj++)
								{
									for(int si = 0;si < subSize;si++)
									{
										input[cur_row*fast_count+cur_count] = Wptr[(k*subSize+sk)*height*width+(j*subSize+sj)*width+(i*subSize+si)];
										cur_row ++;
									}
								}
							}
							
							cur_count ++;
						}
					}
				}
			}
			
			float cuda_cost =
			SubProjection3DCuda(subRow,subCol,fast_count,subMatrix,input,output);
			printf("sub projection cuda cost:%f\n", cuda_cost*0.001);
			
			cur_count = 0;
			for(int k = 0;k < coarseDepth;k++)
			{
				for(int j = 0;j < coarseHeight;j++)
				{
					for(int i = 0;i < coarseWidth;i++)
					{
						if(!regular_flag[k*coarseHeight*coarseWidth+j*coarseWidth+i])
						{
							int cur_row = 0;
							for(int sk = 0;sk < subSize;sk++)
							{
								for(int sj = 0; sj < subSize;sj++)
								{
									for(int si = 1; si < subSize; si++)
									{
										Uptr[(k*subSize+sk)*height*(width+1)+(j*subSize+sj)*(width+1)+(i*subSize+si)] = output[cur_row*fast_count+cur_count];
										cur_row++;
									}
								}
							}
							
							for(int sk = 0;sk < subSize;sk++)
							{
								for(int sj = 1;sj < subSize;sj++)
								{
									for(int si = 0; si < subSize;si++)
									{
										Vptr[(k*subSize+sk)*(height+1)*width+(j*subSize+sj)*width+(i*subSize+si)] = output[cur_row*fast_count+cur_count];
										cur_row++;
									}
								}
							}

							for(int sk = 1;sk < subSize;sk++)
							{
								for(int sj = 0;sj < subSize;sj++)
								{
									for(int si = 0;si < subSize;si++)
									{
										Wptr[(k*subSize+sk)*height*width+(j*subSize+sj)*width+(i*subSize+si)] = output[cur_row*fast_count+cur_count];
										cur_row++;
									}
								}
							}
							
							cur_count++;
						}
					}
				}
			
			}
			
			delete []input;
			delete []output;
		}
		
		for(int k = 0;k < coarseDepth;k++)
		{
			for(int j = 0;j < coarseHeight;j++)
			{
				for(int i = 0;i < coarseWidth;i++)
				{
					if(regular_flag[k*coarseHeight*coarseWidth+j*coarseWidth+i])
					{
						float* sub_input = new float[subSize*subSize*(subSize+1)*3];
						memset(sub_input,0,sizeof(float)*subSize*subSize*(subSize+1)*3);
						float* sub_output = new float[subSize*subSize*(subSize+1)*3];
						memset(sub_output,0,sizeof(float)*subSize*subSize*(subSize+1)*3);
						bool* sub_occupy = new bool[subSize*subSize*subSize];
						memset(sub_occupy,0,sizeof(bool)*subSize*subSize*subSize);

						int i_shift = i*subSize;
						int j_shift = j*subSize;
						int k_shift = k*subSize;
						int v_shift = subSize*subSize*(subSize+1);
						int w_shift = 2*subSize*subSize*(subSize+1);

						for(int sk = 0;sk < subSize;sk++)
						{
							for(int sj = 0;sj < subSize;sj++)
							{
								for(int si = 0;si <= subSize;si++)
								{
									sub_input[sk*subSize*(subSize+1)+sj*(subSize+1)+si] = Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)];
								}
							}
						}

						for(int sk = 0;sk < subSize;sk++)
						{
							for(int sj = 0;sj <= subSize;sj++)
							{
								for(int si = 0;si < subSize;si++)
								{
									sub_input[v_shift+sk*(subSize+1)*subSize+sj*subSize+si] = Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)];
								}
							}
						}
						
						for(int sk = 0;sk <= subSize;sk++)
						{
							for(int sj = 0;sj < subSize;sj++)
							{
								for(int si = 0;si < subSize;si++)
								{
									sub_input[w_shift+sk*subSize*subSize+sj*subSize+si] = Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)];
								}
							}
						}

						for(int sk = 0;sk < subSize;sk++)
						{
							for(int sj = 0;sj < subSize;sj++)
							{
								for(int si = 0;si < subSize;si++)
								{
									sub_occupy[sk*subSize*subSize+sj*subSize+si] = occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)];
								}
							}
						}
						
						memcpy(sub_output,sub_input,sizeof(float)*subSize*subSize*(subSize+1)*3);
						ZQ_SolveSubGrid3D(para,sub_occupy,sub_input,sub_output);
						for(int sk = 0;sk < subSize;sk++)
						{
							for(int sj = 0;sj < subSize;sj++)
							{
								for(int si = 1;si < subSize;si++)
								{
									Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)] = sub_output[sk*subSize*(subSize+1)+sj*(subSize+1)+si];
								}
							}
						}
						
						for(int sk = 0;sk < subSize;sk++)
						{
							for(int sj = 1;sj < subSize;sj++)
							{
								for(int si = 0;si < subSize;si++)
								{
									Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)] = sub_output[v_shift+sk*(subSize+1)*subSize+sj*subSize+si];
								}
							}
						}

						for(int sk = 1;sk < subSize;sk++)
						{
							for(int sj = 0;sj < subSize;sj++)
							{
								for(int si = 0;si < subSize;si++)
								{
									Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)] = sub_output[w_shift+sk*subSize*subSize+sj*subSize+si];
								}
							}
						}
						
						delete []sub_input;
						delete []sub_output;
						delete []sub_occupy;
					}
				}
			}
		}
		
		delete []regular_flag;
		regular_flag = 0;
	}
	else
	{
		for(int k = 0;k < coarseDepth;k++)
		{
			for(int j = 0;j < coarseHeight;j++)
			{
				for(int i = 0;i < coarseWidth;i++)
				{

					float* sub_input = new float[subSize*subSize*(subSize+1)*3];
					memset(sub_input,0,sizeof(float)*subSize*subSize*(subSize+1)*3);
					float* sub_output = new float[subSize*subSize*(subSize+1)*3];
					memset(sub_output,0,sizeof(float)*subSize*subSize*(subSize+1)*3);
					bool* sub_occupy = new bool[subSize*subSize*subSize];
					memset(sub_occupy,0,sizeof(bool)*subSize*subSize*subSize);

					int i_shift = i*subSize;
					int j_shift = j*subSize;
					int k_shift = k*subSize;
					int v_shift = subSize*subSize*(subSize+1);
					int w_shift = 2*subSize*subSize*(subSize+1);

					for(int sk = 0;sk < subSize;sk++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							for(int si = 0;si <= subSize;si++)
							{
								sub_input[sk*subSize*(subSize+1)+sj*(subSize+1)+si] = Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)];
							}
						}
					}

					for(int sk = 0;sk < subSize;sk++)
					{
						for(int sj = 0;sj <= subSize;sj++)
						{
							for(int si = 0;si < subSize;si++)
							{
								sub_input[v_shift+sk*(subSize+1)*subSize+sj*subSize+si] = Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)];
							}
						}
					}

					for(int sk = 0;sk <= subSize;sk++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							for(int si = 0;si < subSize;si++)
							{
								sub_input[w_shift+sk*subSize*subSize+sj*subSize+si] = Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)];
							}
						}
					}

					for(int sk = 0;sk < subSize;sk++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							for(int si = 0;si < subSize;si++)
							{
								sub_occupy[sk*subSize*subSize+sj*subSize+si] = occupyPtr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)];
							}
						}
					}

					memcpy(sub_output,sub_input,sizeof(float)*subSize*subSize*(subSize+1)*3);
					ZQ_SolveSubGrid3D(para,sub_occupy,sub_input,sub_output);
					for(int sk = 0;sk < subSize;sk++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							for(int si = 1;si < subSize;si++)
							{
								Uptr[(k_shift+sk)*height*(width+1)+(j_shift+sj)*(width+1)+(i_shift+si)] = sub_output[sk*subSize*(subSize+1)+sj*(subSize+1)+si];
							}
						}
					}

					for(int sk = 0;sk < subSize;sk++)
					{
						for(int sj = 1;sj < subSize;sj++)
						{
							for(int si = 0;si < subSize;si++)
							{
								Vptr[(k_shift+sk)*(height+1)*width+(j_shift+sj)*width+(i_shift+si)] = sub_output[v_shift+sk*(subSize+1)*subSize+sj*subSize+si];
							}
						}
					}

					for(int sk = 1;sk < subSize;sk++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							for(int si = 0;si < subSize;si++)
							{
								Wptr[(k_shift+sk)*height*width+(j_shift+sj)*width+(i_shift+si)] = sub_output[w_shift+sk*subSize*subSize+sj*subSize+si];
							}
						}
					}

					delete []sub_input;
					delete []sub_output;
					delete []sub_occupy;

				}
			}
		}
	}
	return true;
}


bool ZQ_TotalProjection3D(ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA,
								const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid3D* guideGrid)
{

	LARGE_INTEGER tf,solve_global_st, solve_global_end;
	QueryPerformanceFrequency(&tf);
	QueryPerformanceCounter(&solve_global_st);
	
	if(!ZQ_MapFine2Coarse3D(fineGrid,coarseGrid,para))
		return false;
	ZQ_SmokeGrid3D* deltaGrid = new ZQ_SmokeGrid3D(para.coarseWidth,para.coarseHeight,para.coarseDepth);

	if(para.globalGridSolver == CLOSED_POISSON)
	{
		if(!ZQ_SolveClosedPressure3D(coarseGrid,para,coarseA))
			return false;
		if(!ZQ_AdjustClosedVelocity3D(coarseGrid,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == OPEN_POISSON)
	{
		
		if(!ZQ_SolveOpenPressure3D(coarseGrid,para,coarseA))
			return false;
		
		if(!ZQ_AdjustOpenVelocity3D(coarseGrid,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
		
	}
	else if(para.globalGridSolver == OPEN_FLUX || para.globalGridSolver == CLOSED_FLUX)
	{
		if(!ZQ_SolveFlux3D(coarseGrid,para,coarseA,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_OPEN_POISSON_SOR)
	{
		if(!ZQ_SolveCoarseOpenPoisson3DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_CLOSED_POISSON_SOR)
	{
		if(!ZQ_SolveCoarseOpenPoisson3DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_OPEN_FLUX_SOR)
	{
		if(!ZQ_SolveCoarseOpenFlux3DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_CLOSED_FLUX_SOR)
	{
		if(!ZQ_SolveCoarseClosedFlux3DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}

	if(para.useGuidingControl)
	{
		if(!ZQ_ApplyGuidingVelocity3D(coarseGrid,guideGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	
	if(!ZQ_MapCoarse2Fine3D(deltaGrid,fineGrid,para))
	{
		delete deltaGrid;
		return false;
	}
	delete deltaGrid;
	QueryPerformanceCounter(&solve_global_end);
	printf("solve global cost:%f\n",1.0*(solve_global_end.QuadPart-solve_global_st.QuadPart)/tf.QuadPart);

	if(para.subSize > 1)
	{
		LARGE_INTEGER solve_sub_st, solve_sub_end;
		QueryPerformanceCounter(&solve_sub_st);
		if(!ZQ_SubProjection3D(fineGrid,para,subMatrix,subRow,subCol,true))
		{
			return false;
		}
		QueryPerformanceCounter(&solve_sub_end);
		printf("solve sub cost:%f\n",1.0*(solve_sub_end.QuadPart-solve_sub_st.QuadPart)/tf.QuadPart);
	}
	return true;
}

bool ZQ_AdvectScalar3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	if (fineGrid == 0)
		return false;

	ZQ_CUDA_Advection3D::AdvectScalarType type;
	switch (para.advScaType)
	{
	case ZQ_SmokeSimulationPara3D::ADV_SCA_MAC:
		type = ZQ_CUDA_Advection3D::ADV_SCA_MAC;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_SCA_MAC_BFECC:
		type = ZQ_CUDA_Advection3D::ADV_SCA_MAC_BFECC;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_SCA_REG:
		type = ZQ_CUDA_Advection3D::ADV_SCA_REG;
		break;
	case ZQ_SmokeSimulationPara3D::ADV_SCA_REG_BFECC:
		type = ZQ_CUDA_Advection3D::ADV_SCA_REG_BFECC;
		break;
	default:
		return false;
		break;
	}

	unsigned int width = para.width;
	unsigned int height = para.height;
	unsigned int depth = para.depth;
	unsigned int steps = para.steps;
	float voxelSize = para.voxelSize;
	float deltatt = para.stepTime / steps;

	bool*& occupy = fineGrid->GetOccupyPtr();
	float*& u = fineGrid->GetVelocityUptr();
	float*& v = fineGrid->GetVelocityVptr();
	float*& w = fineGrid->GetVelocityWptr();
	float*& temperaturePtr = fineGrid->GetTemperaturePtr();
	float*& densityPtr = fineGrid->GetDensityPtr();

	float cuda_cost_time = ZQ_CUDA_Advection3D::Advect_Scalar(u, v, w, occupy, temperaturePtr, densityPtr, temperaturePtr, densityPtr, type);

	printf("scalar advection cuda cost = %f\n", 0.001*cuda_cost_time);
	return true;
}


bool ZQ_ReinjectTemperatureDensity3D(ZQ_SmokeGrid3D* fineGrid, const ZQ_SmokeSimulationPara3D& para)
{
	int width = para.width;
	int height = para.height;
	int depth = para.depth;

	float*& density = fineGrid->GetDensityPtr();
	float*& temperature = fineGrid->GetTemperaturePtr();

	bool randDenFlag = para.randDensity;
	
	for(int cc = 0;cc < para.sources.size();cc++)
	{
		float cx = para.sources[cc].cx * width;
		float cy = para.sources[cc].cy * height;
		float cz = para.sources[cc].cz * depth;
		float half_xlen = para.sources[cc].half_xlen * width;
		float half_ylen = para.sources[cc].half_ylen * height;
		float half_zlen = para.sources[cc].half_zlen * depth;
		float reinject_den = para.sources[cc].reinjectDensity;
		float reinject_tem = para.sources[cc].reinjectTemperature;
		for(int k = 0;k < depth;k++)
		{
			for(int j = 0;j < height;j++)
			{
				for(int i = 0;i < width;i++)
				{
					if(fabs(i+0.5f-cx) <= half_xlen && fabs(j+0.5f-cy) <= half_ylen && fabs(k+0.5f-cz) <= half_zlen)
					{
						density[k*height*width+j*width+i] = reinject_den * (randDenFlag ? rand()%256/255.0f : 1.0f);
						temperature[k*height*width+j*width+i] = reinject_tem;
					}
				}
			}
		}
		
	}
	
	return true;

}

bool ZQ_UpdateOneFrame3D(ZQ_SmokeGrid3D* fineGrid, ZQ_SmokeGrid3D* coarseGrid, const ZQ_SmokeSimulationPara3D& para, const taucs_ccs_matrix* coarseA,
								const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid3D* guideGrid)
{
	if(!ZQ_AddForce3D(fineGrid,para))
		return false;
	if(!ZQ_Attenuation3D(fineGrid,para))
		return false;

	if(!ZQ_AdvectVelocity3D(fineGrid,para))
		return false;

	if(!ZQ_TotalProjection3D(fineGrid,coarseGrid,para,coarseA,subMatrix,subRow,subCol,guideGrid))
		return false;

	if(!ZQ_AdvectScalar3D(fineGrid,para))
		return false;

	if(!ZQ_ReinjectTemperatureDensity3D(fineGrid,para))
		return false;
	return true;
}