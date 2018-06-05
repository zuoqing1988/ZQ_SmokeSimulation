#include "ZQ_SmokeSimulationUtils2D.h"
#include "ZQ_SparseMatrix.h"
#include "ZQ_Matrix.h"
#include "ZQ_SVD.h"
#include "ZQ_TaucsBase.h"
#include "ZQ_PoissonSolver.h"
#include "windows.h"

ZQ_SmokeSimulationPara2D::ZQ_SmokeSimulationPara2D()
{
	width = height = 128;
	subSize = 4;
	coarseWidth = width/subSize;
	coarseHeight = height/subSize;
	gravityCoeff = 0;
	buoyancyCoeff = 1;
	confineCoeff = 0;
	Tamb = 10;
	voxelSize = 1.0 / width;
	
	ZQ_SmokeSource2D source;
	source.cx = 0.5;
	source.cy = 0.1;
	source.half_xlen = 0.1;
	source.half_ylen = 0.05;
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
	useBFECC = true;
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

}

ZQ_SmokeSimulationPara2D::~ZQ_SmokeSimulationPara2D()
{
}

bool ZQ_SmokeSimulationPara2D::LoadFromFile(const char *file, bool display)
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
			sscanf(s,"%s%d%d",arg1,&width,&height);
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
			ZQ_SmokeSource2D source;
			sscanf(s,"%s%f%f%f%f%f%f",arg1,&source.cx,&source.cy,&source.half_xlen,&source.half_ylen,&source.reinjectDensity,&source.reinjectTemperature);
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
		else if(_strcmpi(arg1,"BFECC") == 0)
		{
			if(strcmp(arg2,"true") == 0)
				useBFECC = true;
			else 
				useBFECC = false;
		}
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
	if(width%subSize != 0 || height%subSize != 0)
	{
		printf("Error: Grid[%d x %d] cannot be divided by subSize[%d]\n",width,height,subSize);
		return false;
	}
	coarseWidth = width/subSize;
	coarseHeight = height/subSize;
	voxelSize = 1.0/width;
	if(sources.size() == 0)
	{
		printf("ERROR: no smoke source exist\n");
		return false;
	}

	if(display)
	{
		printf("resolution:\t\t%d x %d\n",width,height);
		printf("subSize:\t\t%d\n",subSize);
		printf("coarse:\t\t%d x %d\n",coarseWidth,coarseHeight);
		printf("Tamb:\t\t%f\n",Tamb);
		printf("confineCoeff:\t\t%f\n",confineCoeff);
		printf("velAttenCoeff:\t\t%f\n",velAttenCoeff);
		printf("tempAttenCoeff:\t\t%f\n",tempAttenCoeff);
		printf("densityAttenCoeff:\t\t%f\n",densityAttenCoeff);
		for(int cc = 0;cc < sources.size();cc++)
		{
			printf("Source[%d]: center(%.3f,%.3f), size(%.3f,%.3f), density(%.3f), temperature(%.3f)\n",
				cc,sources[cc].cx,sources[cc].cy,sources[cc].half_xlen,sources[cc].half_ylen,sources[cc].reinjectDensity,sources[cc].reinjectTemperature);
		}
		printf("voxelSize:\t\t%f\n",voxelSize);
		printf("steps:\t\t%d\n",steps);
		printf("stepTime:\t\t%f\n",stepTime);
		printf("frameCount:\t\t%d\n",frameCount);
	}
	return true;
}



int ZQ_InitMatOpenPoisson2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A)
{

	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();

	int* index = new int[width*height];
	memset(index,0,sizeof(int)*width*height);
	int idx = 1;
	
	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();

	for(int i = 0;i < height*width;i++)
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

	int ISLICE = width;
	int JSLICE = 1;

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			int row_id = index[offset]-1;
			if(row_id < 0)
				continue;

			std::vector<int> indices;
			std::vector<int> u_indices;
			std::vector<int> v_indices;
			if(i > 0)
			{
				indices.push_back(offset-ISLICE);
				u_indices.push_back(-1);
				v_indices.push_back(i*width+j);
			}
			if(i < height-1)
			{
				indices.push_back(offset+ISLICE);
				u_indices.push_back(-1);
				v_indices.push_back((i+1)*width+j);
			}
			if(j > 0)
			{
				indices.push_back(offset-JSLICE);
				u_indices.push_back(i*(width+1)+j);
				v_indices.push_back(-1);
			}
			if(j < width-1)
			{
				indices.push_back(offset+JSLICE);
				u_indices.push_back(i*(width+1)+j+1);
				v_indices.push_back(-1);
			}

			float count = 4-indices.size();
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
					count += tmp_ratio;
					mat.SetValue(row_id,cur_index-1, tmp_ratio);
				}
			}

			mat.SetValue(row_id,row_id, -count);
		}
	}


	delete []index;

	*A = mat.ExportCCS(TAUCS_DOUBLE);

	return dim;
}

int ZQ_InitMatClosedPoisson2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A)
{	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetWidth();

	int* index = new int[width*height];
	memset(index,0,sizeof(int)*width*height);
	
	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();

	int idx = 1;
	bool haveFirst = false;
	int first = -1;
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(!occupyPtr[i*width+j])
			{
				if(!haveFirst)
				{
					haveFirst = true;
					first = i*width+j;
				}
				index[i*width+j] = idx++;
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

	int ISLICE = width;
	int JSLICE = 1;

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			int row_id = index[offset]-1;

			if(row_id < 0)
				continue;
			
			std::vector<int> indices;
			std::vector<int> u_indices;
			std::vector<int> v_indices;

			if(i > 0)
			{
				indices.push_back(offset-ISLICE);
				u_indices.push_back(-1);
				v_indices.push_back(i*width+j);
			}
			if(i < height-1)
			{
				indices.push_back(offset+ISLICE);
				u_indices.push_back(-1);
				v_indices.push_back((i+1)*width+j);
			}
			if(j > 0)
			{
				indices.push_back(offset-JSLICE);
				u_indices.push_back(i*(width+1)+j);
				v_indices.push_back(-1);
			}
			if(j < width-1)
			{
				indices.push_back(offset+JSLICE);
				u_indices.push_back(i*(width+1)+j+1);
				v_indices.push_back(-1);
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
					count += tmp_ratio;
					if(cur_index != 1)
						mat.SetValue(row_id,cur_index-2, tmp_ratio);
				}
			}

			if(index[offset] != 1)
				mat.SetValue(row_id,index[offset]-2, -count);

		}
	}

	delete []index;

	*A = mat.ExportCCS(TAUCS_DOUBLE);

	return dim;
}

int ZQ_InitMatOpenFlux2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A)
{
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();

	int* u_index = new int[height*(width+1)];
	int* v_index = new int[width*(height+1)];
	int* lamda_index = new int[width*height];

	for(int i = 0;i < height*(width+1);i++)
		u_index[i] = -1;
	for(int i = 0;i < width*(height+1);i++)
		v_index[i] = -1;
	for(int i = 0;i < width*height;i++)
		lamda_index[i] = -1;

	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();

	int idx = 0;
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j <= width;j++)
		{
			if(unoccupyUptr[i*(width+1)+j] != 0)
			{
				u_index[i*(width+1)+j] = idx++;
			}
		}
	}

	for(int i = 0;i <= height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(unoccupyVptr[i*width+j] != 0)
			{
				v_index[i*width+j] = idx++;
				
			}
		}
	}

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(!occupyPtr[i*width+j])
				lamda_index[i*width+j] = idx++;
		}
	}

	int dim = idx;

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);

	idx = 0;
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j <= width;j++)
		{
			int cur_u_index = u_index[i*(width+1)+j];

			if(cur_u_index >= 0)
			{
				float ratio = unoccupyUptr[i*(width+1)+j];
				mat.SetValue(cur_u_index,cur_u_index, 2*ratio);
				if(j > 0)
					mat.SetValue(cur_u_index,lamda_index[i*width+j-1], -ratio);
				if(j < width)
					mat.SetValue(cur_u_index,lamda_index[i*width+j], ratio);
			}
		}
	}

	for(int i = 0;i <= height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int cur_v_index = v_index[i*width+j];

			if(cur_v_index >= 0)
			{
				float ratio = unoccupyVptr[i*width+j];
				mat.SetValue(cur_v_index,cur_v_index, 2*ratio);
				if(i > 0)
					mat.SetValue(cur_v_index,lamda_index[(i-1)*width+j], -ratio);
				if(i < height)
					mat.SetValue(cur_v_index,lamda_index[i*width+j], ratio);
			}
		}
	}

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int cur_lambda_index = lamda_index[i*width+j];
			if(cur_lambda_index >= 0)
			{
				float ratio = 0;
				ratio = unoccupyUptr[i*(width+1)+j+1];
				if(ratio > 0)
					mat.SetValue(cur_lambda_index,u_index[i*(width+1)+j+1], -ratio);
				ratio = unoccupyUptr[i*(width+1)+j];
				if(ratio > 0)
					mat.SetValue(cur_lambda_index,u_index[i*(width+1)+j], ratio);
				ratio = unoccupyVptr[(i+1)*width+j];
				if(ratio > 0)
					mat.SetValue(cur_lambda_index,v_index[(i+1)*width+j], -ratio);
				ratio = unoccupyVptr[i*width+j];
				if(ratio > 0)
					mat.SetValue(cur_lambda_index,v_index[i*width+j], ratio);
			}
		}
	}

	*A = mat.ExportCCS(TAUCS_DOUBLE);

	delete []u_index;
	delete []v_index;
	delete []lamda_index;

	return dim;
}

int ZQ_InitMatClosedFlux2D(const ZQ_SmokeGrid2D* coarseGrid, taucs_ccs_matrix** A)
{
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	
	int* u_index = new int[height*(width+1)];
	int* v_index = new int[width*(height+1)];
	int* lamda_index = new int[width*height];
	int k_index = 0;

	for(int i = 0;i < height*(width+1);i++)
		u_index[i] = -1;

	for(int i = 0;i < width*(height+1);i++)
		v_index[i] = -1;

	for(int i = 0;i < width*height;i++)
		lamda_index[i] = -1;

	const bool*& occupyPtr = coarseGrid->GetOccupyPtr();
	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& unoccupyVolumePtr = coarseGrid->GetVolumeUnOccupyRatioPtr();

	int idx = 0;
	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{
			if(unoccupyUptr[i*(width+1)+j] != 0)
				u_index[i*(width+1)+j] = idx++;
		}
	}
	for(int i = 1;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(unoccupyVptr[i*width+j] != 0)
				v_index[i*width+j] = idx++;
		}
	}
	
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(!occupyPtr[i*width+j])
				lamda_index[i*width+j] = idx++;
		}
	}
	k_index = idx++;
	int dim = idx;
	idx = 0;

	ZQ::ZQ_SparseMatrix<float> mat(dim,dim);
	
	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{
			int cur_u_index = u_index[i*(width+1)+j];
			if(cur_u_index >= 0)
			{
				float ratio = unoccupyUptr[i*(width+1)+j];
				mat.SetValue(cur_u_index,cur_u_index, 2*ratio);
				mat.SetValue(cur_u_index,lamda_index[i*width+j-1], -ratio);
				mat.SetValue(cur_u_index,lamda_index[i*width+j], ratio);
			}
		}
	}

	for(int i = 1;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int cur_v_index = v_index[i*width+j];
			if(cur_v_index >= 0)
			{
				float ratio = unoccupyVptr[i*width+j];
				mat.SetValue(cur_v_index,cur_v_index, 2*ratio);
				mat.SetValue(cur_v_index,lamda_index[(i-1)*width+j], -ratio);
				mat.SetValue(cur_v_index,lamda_index[i*width+j], ratio);
			}
		}
	}

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int cur_lambda_index = lamda_index[i*width+j];
			if(cur_lambda_index >= 0)
			{
				float ratio = 0;
				ratio = unoccupyUptr[i*(width+1)+j+1];
				if(ratio > 0 && u_index[i*(width+1)+j+1] >= 0)
					mat.SetValue(cur_lambda_index,u_index[i*(width+1)+j+1], -ratio);

				ratio = unoccupyUptr[i*(width+1)+j];
				if(ratio > 0 && u_index[i*(width+1)+j] >= 0)
					mat.SetValue(cur_lambda_index,u_index[i*(width+1)+j], ratio);

				ratio = unoccupyVptr[(i+1)*width+j];
				if(ratio > 0 && v_index[(i+1)*width+j] >= 0)
					mat.SetValue(cur_lambda_index,v_index[(i+1)*width+j], -ratio);

				ratio = unoccupyVptr[i*width+j];
				if(ratio > 0 && v_index[i*width+j] >= 0)
					mat.SetValue(cur_lambda_index,v_index[i*width+j], ratio);

				ratio = unoccupyVolumePtr[i*width+j];
				if(ratio > 0)
					mat.SetValue(cur_lambda_index,k_index, ratio);
			}
		}
	}


	//k_index
	{
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(lamda_index[i*width+j] >= 0)
				{
					float ratio = unoccupyVolumePtr[i*width+j];
					mat.SetValue(k_index,lamda_index[i*width+j], ratio);
				}
			}
		}
	}

	*A = mat.ExportCCS(TAUCS_DOUBLE);
	
	delete []u_index;
	delete []v_index;
	delete []lamda_index;

	return dim;
}

bool ZQ_InitSubMatPoisson2D(const int subSize, float** subMatrix, int* row, int* col)
{
	typedef ZQ::ZQ_Matrix<double> ZQ_Mat;
	ZQ_Mat matP(subSize*subSize,subSize*subSize-1);
	ZQ_Mat matV(subSize*subSize,2*subSize*(subSize+1));
	ZQ_Mat matD(2*subSize*(subSize+1),subSize*subSize);


	int ushift = 0;
	int vshift = subSize*(subSize+1);

	for(int i = 0;i < subSize;i++)
	{
		for(int j = 0;j < subSize;j++)
		{
			int vcol1 = ushift + i*(subSize+1) + j;
			int vcol2 = ushift + i*(subSize+1) + (j+1);
			int vcol3 = vshift + i*subSize + j;
			int vcol4 = vshift + (i+1)*subSize + j;

			int vrow = i*subSize+j;
			matV.SetData(vrow,vcol1,-1);
			matV.SetData(vrow,vcol2,1);
			matV.SetData(vrow,vcol3,-1);
			matV.SetData(vrow,vcol4,1);

			int prow = i*subSize+j;
			int pcol0 = prow - 1;
			int pcol1 = i*subSize+(j-1) - 1;
			int pcol2 = i*subSize+(j+1) - 1;
			int pcol3 = (i-1)*subSize+j - 1;
			int pcol4 = (i+1)*subSize+j - 1;

			if(prow == 0)
			{
				matP.SetData(0,pcol2,1);
				matP.SetData(0,pcol4,1);
				continue;
			}
			
			int count = 0;
			if(j == 0)
			{
				count ++;
				matP.SetData(prow,pcol2,1);
			}
			else if(j == subSize-1)
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
			if(i == 0)
			{
				count ++;
				matP.SetData(prow,pcol4,1);
			}
			else if(i == subSize-1)
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

			matP.SetData(prow,pcol0,-count);
		}
	}

	for(int i = 0;i < subSize;i++)
	{
		for(int j = 1;j < subSize;j++)
		{
			int drow = ushift + i*(subSize+1)+j;
			int dcol1 = i*subSize + (j-1);
			int dcol2 = i*subSize + j;
			matD.SetData(drow,dcol1,-1);
			matD.SetData(drow,dcol2,1);
		}
	}

	for(int i = 1;i < subSize;i++)
	{
		for(int j = 0;j < subSize;j++)
		{
			int drow = vshift + i*subSize+j;
			int dcol1 = (i-1)*subSize+j;
			int dcol2 = i*subSize+j;
			matD.SetData(drow,dcol1,-1);
			matD.SetData(drow,dcol2,1);
		}
	}
	
	ZQ_Mat invP(subSize*subSize-1,subSize*subSize);

	ZQ::ZQ_SVD::Invert(matP,invP);
	ZQ_Mat invPV = invP*matV;

	ZQ_Mat invPV1(subSize*subSize,2*subSize*(subSize+1));
	bool flag;
	for(int rr = 1; rr < subSize*subSize;rr++)
	{
		for(int cc = 0;cc < 2*subSize*(subSize+1);cc++)
			invPV1.SetData(rr,cc,invPV.GetData(rr-1,cc,flag));
	}

	ZQ_Mat DinvPV = matD*invPV1;

	int tmpheight = 2*subSize*(subSize+1);
	int tmpwidth = 2*subSize*(subSize+1);
	int height = 2*subSize*(subSize-1);
	int width = 2*subSize*(subSize+1);

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
	for(int i = 0;i < subSize;i++)
	{
		for(int j = 1;j < subSize;j++)
		{
			int wholeCount = i*(subSize+1)+j;
			memcpy((*subMatrix)+useCount*width,tmpMat+wholeCount*tmpwidth,sizeof(float)*width);
			useCount ++;
		}
	}
	for(int i = 1; i < subSize;i++)
	{
		for(int j = 0;j < subSize;j++)
		{
			int wholeCount = subSize*(subSize+1) + i*subSize+j;
			memcpy((*subMatrix)+useCount*width,tmpMat+wholeCount*tmpwidth,sizeof(float)*width);
			useCount ++;
		}
	}
	

	delete []tmpMat;
	*row = height;
	*col = width;

	return true;
}

bool ZQ_SolveSubGrid2D(const ZQ_SmokeSimulationPara2D& para, const bool* occupy, const float* input, float* output)
{
	if(occupy == 0 || input == 0 || output == 0)
		return false;
	int subSize = para.subSize;
	ZQ_SmokeGrid2D* tmpGrid = new ZQ_SmokeGrid2D(subSize,subSize);
	
	tmpGrid->BuildFaceUnOccupyRatio();
	float*& unoccupyUptr = tmpGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyVptr = tmpGrid->GetFaceUnOccupyRatioVPtr();
	bool*& occupyPtr = tmpGrid->GetOccupyPtr();
	
	for(int i = 0;i <= subSize;i++)
	{
		for(int j = 0;j < subSize;j++)
		{
			unoccupyUptr[j*(subSize+1)+i] = 1;
			unoccupyVptr[i*subSize+j] = 1;
		}
	}
	for(int i = 0;i < subSize;i++)
	{
		for(int j = 0;j < subSize;j++)
		{
			occupyPtr[i*subSize+j] = occupy[i*subSize+j];

			if(occupy[i*subSize+j])
			{
				unoccupyUptr[i*(subSize+1)+j] = 0;
				unoccupyUptr[i*(subSize+1)+j+1] = 0;
				unoccupyVptr[i*subSize+j] = 0;
				unoccupyVptr[(i+1)*subSize+j] = 0;
			}
		}
	}

	tmpGrid->BuildVolumeUnOccupyRatio();
	float*& unoccupyVolumePtr = tmpGrid->GetVolumeUnOccupyRatioPtr();
	for(int i = 0;i < subSize;i++)
	{
		for(int j = 0;j < subSize;j++)
		{
			if(occupy[i*subSize+j])
				unoccupyVolumePtr[i*subSize+j] = 0;
			else
				unoccupyVolumePtr[i*subSize+j] = 1;
		}
	}

	float*& U = tmpGrid->GetVelocityUptr();
	float*& V = tmpGrid->GetVelocityVptr();
	memcpy(U,input,sizeof(float)*subSize*(subSize+1));
	memcpy(V,input+subSize*(subSize+1),sizeof(float)*subSize*(subSize+1));
	
	taucs_ccs_matrix* tmpA = 0;
	
	ZQ_SmokeSimulationPara2D sub_para;
	sub_para.maxIter = 2*(subSize*subSize+2*subSize*(subSize+1)+1);
	
	sub_para.globalGridSolver = CLOSED_POISSON;
	sub_para.width = subSize;
	sub_para.height = subSize;
	sub_para.subSize = 1;

	ZQ_SmokeGrid2D* deltaGrid = new ZQ_SmokeGrid2D(subSize,subSize);
	
	if(!ZQ_InitMatClosedPoisson2D(tmpGrid,&tmpA))
	{
		delete tmpGrid;
		delete deltaGrid;
		return false;
	}
	ZQ_SolveClosedPressure2D(tmpGrid,sub_para,tmpA);
	ZQ_AdjustClosedVelocity2D(tmpGrid,deltaGrid);
	

	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();

	for(int i = 0;i < subSize*(subSize+1);i++)
	{
		output[i] = U[i]+deltaU[i];
		output[i+subSize*(subSize+1)] = V[i]+deltaV[i];
	}
	if(tmpA)
		ZQ::ZQ_TaucsBase::ZQ_taucs_ccs_free(tmpA);
	
	delete tmpGrid;
	delete deltaGrid;
	return true;
}

bool ZQ_InitFineGrid2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	int width = para.width;
	int height = para.height;
	
	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	float*& pressurePtr = fineGrid->GetPressurePtr();

	memset(occupyPtr,0,sizeof(bool)*width*height);
	memset(pressurePtr,0,sizeof(float)*width*height);

	if(para.globalGridSolver == CLOSED_POISSON || para.globalGridSolver == CLOSED_FLUX 
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON || para.globalGridSolver == CLOSED_OCTREE_FLUX
		|| para.globalGridSolver == COARSE_CLOSED_FLUX_SOR || para.globalGridSolver == COARSE_CLOSED_POISSON_SOR
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON_SOR)
	{
		for(int j = 0;j < width;j++)
		{
			occupyPtr[0*width+j] = true;
			occupyPtr[1*width+j] = true;
			occupyPtr[2*width+j] = true;
			occupyPtr[3*width+j] = true;
			occupyPtr[(height-1)*width+j] = true;
			occupyPtr[(height-2)*width+j] = true;
			occupyPtr[(height-3)*width+j] = true;
			occupyPtr[(height-4)*width+j] = true;
		}
		for(int i = 0;i < height;i++)
		{
			occupyPtr[i*width+0] = true;
			occupyPtr[i*width+1] = true;
			occupyPtr[i*width+2] = true;
			occupyPtr[i*width+3] = true;
			occupyPtr[i*width+width-1] = true;
			occupyPtr[i*width+width-2] = true;
			occupyPtr[i*width+width-3] = true;
			occupyPtr[i*width+width-4] = true;
		}

	}
	
	switch(para.intobj)
	{
	case ZQ_SmokeSimulationPara2D::INT_OBJ_SPHERE:
		{
			for(int i = 0;i < height;i++)
			{
				for(int j = 0;j < width;j++)
				{
					if((i+0.5-height/2)*(i+0.5-height/2) + (j+0.5-width/2)*(j+0.5-width/2) < width*width/400.0)
						occupyPtr[i*width+j] = true;
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara2D::INT_OBJ_BOX:
		{
			for(int i = 0;i < height;i++)
			{
				for(int j = 0;j < width;j++)
				{
					if(fabs(i+0.5-height/2) < height/20 && fabs(j+0.5-width/2) < width/10)
						occupyPtr[i*width+j] = true;
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara2D::INT_OBJ_OVAL:
		{
			for(int i = 0;i < height;i++)
			{
				for(int j = 0;j < width;j++)
				{
					if(sqrt((j + 0.5 - 3.0*width/7)*(j + 0.5 - 3.0*width/7) + (i+0.5-height/2)*(i+0.5-height/2))
						+ sqrt((j + 0.5 - 4.0*width/7)*(j + 0.5 - 4.0*width/7) + (i+0.5-height/2)*(i+0.5-height/2))
						< width / 5.0)
					{
						occupyPtr[i*width+j] = true;
					}
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara2D::INT_OBJ_SPECIAL1:
		{
			for(int i = 0;i < height;i++)
			{
				for(int j = 0;j < width;j++)
				{
					if((i == height/2+13 || i == height/2+14) && (j%4 == 0 || j%4 == 3) && j >=4 && j < width-4)
						occupyPtr[i*width+j] = true;
				}
			}
		}
		break;
	case ZQ_SmokeSimulationPara2D::INT_OBJ_NULL: default:
		break;
	}
	
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	memset(Uptr,0,sizeof(float)*height*(width+1));
	memset(Vptr,0,sizeof(float)*width*(height+1));

	return true;
}

bool ZQ_InitCoarseGrid2D(const ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid)
{
	if(fineGrid == 0 || coarseGrid == 0)
		return false;

	int width = fineGrid->GetWidth();
	int height = fineGrid->GetWidth();

	int coarseWidth = coarseGrid->GetWidth();
	int coarseHeight = coarseGrid->GetHeight();

	int subSize = width/coarseWidth;
	if(subSize*coarseWidth != width || subSize*coarseHeight != height)
		return false;

	coarseGrid->BuildFaceUnOccupyRatio();
	
	const bool*& fine_occupyPtr = fineGrid->GetOccupyPtr();
	float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();

	int subFaceSize = subSize;
	for(int i = 0;i < coarseHeight;i++)
	{
		for(int j = 1; j < coarseWidth;j++)
		{
			int occupyNum = 0;
			int i_shift = i*subSize;
			int j_shift = j*subSize;
			for(int si = 0;si < subSize;si++)
			{
				if(fine_occupyPtr[(i_shift+si)*width+j_shift] || fine_occupyPtr[(i_shift+si)*width+j_shift-1])
					occupyNum ++;
			}

			unoccupyUptr[i*(coarseWidth+1)+j] = 1.0f - (float)occupyNum/subFaceSize;
		}
	}

	for(int i = 0;i < coarseHeight;i++)
	{
		int i_shift = i*subSize;
		int occupyNum = 0;
		for(int si = 0;si < subSize;si++)
		{
			if(fine_occupyPtr[(i_shift+si)*width+0])
				occupyNum ++;
		}
		unoccupyUptr[i*(coarseWidth+1)+0] = 1.0f - (float)occupyNum/subFaceSize;
		occupyNum = 0;
		for(int si = 0;si < subSize;si++)
		{
			if(fine_occupyPtr[(i_shift+si)*width+width-1])
				occupyNum ++;
		}
		unoccupyUptr[i*(coarseWidth+1)+coarseWidth] = 1.0f - (float)occupyNum/subFaceSize;
	}
	
	for(int i = 1;i < coarseHeight;i++)
	{
		for(int j = 0;j < coarseWidth;j++)
		{
			int i_shift = i*subSize;
			int j_shift = j*subSize;
			int occupyNum = 0;
			for(int sj = 0; sj < subSize;sj++)
			{					
				if(fine_occupyPtr[i_shift*width+(j_shift+sj)] || fine_occupyPtr[(i_shift-1)*width+(j_shift+sj)])
					occupyNum ++;
			}
			unoccupyVptr[i*coarseWidth+j] = 1.0f - (float)occupyNum/subFaceSize;
		}
	}
	for(int j = 0;j < coarseWidth;j++)
	{
		int j_shift = j*subSize;

		int occpuyNum = 0;
		for(int sj = 0; sj < subSize;sj++)
		{
			if(fine_occupyPtr[0*width+(j_shift+sj)])
				occpuyNum ++;
		}
		unoccupyVptr[0*coarseWidth+j] = 1.0f - (float)occpuyNum/subFaceSize;
		occpuyNum = 0;
		for(int sj = 0; sj < subSize;sj++)
		{
			if(fine_occupyPtr[(height-1)*width+(j_shift+sj)])
				occpuyNum ++;
		}
		unoccupyVptr[coarseHeight*coarseWidth+j] = 1.0f - (float)occpuyNum/subFaceSize;
	}
	
	bool*& coarse_occupyPtr = coarseGrid->GetOccupyPtr();
	for(int i = 0;i < coarseHeight;i++)
	{
		for(int j = 0;j < coarseWidth;j++)
		{
			if(unoccupyUptr[i*(coarseWidth+1)+j] == 0
				&& unoccupyUptr[i*(coarseWidth+1)+j+1] == 0
				&& unoccupyVptr[i*coarseWidth+j] == 0
				&& unoccupyVptr[(i+1)*coarseWidth+j] == 0)
			{
				coarse_occupyPtr[i*coarseWidth+j] = true;
			}
			else
				coarse_occupyPtr[i*coarseWidth+j] = false;
		}
	}
	
	coarseGrid->BuildVolumeUnOccupyRatio();
	float*& unoccupyVolumePtr = coarseGrid->GetVolumeUnOccupyRatioPtr();
	for(int i = 0;i < coarseHeight;i++)
	{
		for(int j = 0;j < coarseWidth;j++)
		{
			int occupyNum = 0;
			for(int si = 0;si < subSize;si++)
			{
				for(int sj = 0;sj < subSize;sj++)
				{
					if(fine_occupyPtr[(i*subSize+si)*width+(j*subSize+sj)])
						occupyNum ++;
				}
			}
			unoccupyVolumePtr[i*coarseWidth+j] = 1.0f - (float)occupyNum/(subSize*subSize);
		}
	}
	return true;
}

bool ZQ_InitSimulation2D(ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, taucs_ccs_matrix** coarseA,
						 float** subMatrix, int* subRow, int* subCol)
{
	if(!ZQ_InitFineGrid2D(fineGrid,para))
		return false;
	
	if(!ZQ_InitCoarseGrid2D(fineGrid,coarseGrid))
	{
		return false;
	}

	int subSize = para.subSize;
	
	int matrixsz = 0;
	if(para.globalGridSolver == OPEN_POISSON)
	{
		if((matrixsz = ZQ_InitMatOpenPoisson2D(coarseGrid,coarseA)) == 0)
			return false;
	}
	else if(para.globalGridSolver == CLOSED_POISSON)
	{
		if((matrixsz =  ZQ_InitMatClosedPoisson2D(coarseGrid,coarseA)) == 0)
			return false;
	}
	else if(para.globalGridSolver == OPEN_FLUX)
	{
		if((matrixsz = ZQ_InitMatOpenFlux2D(coarseGrid,coarseA)) == 0)
			return false;
	}
	else if(para.globalGridSolver == CLOSED_FLUX)
	{
		if((matrixsz = ZQ_InitMatClosedFlux2D(coarseGrid,coarseA)) == 0)
			return false;
	}

	if(para.subSize >= 2)
	{
		if(!ZQ_InitSubMatPoisson2D(subSize, subMatrix,subRow,subCol))
			return false;
	}
	return true;
}

bool ZQ_ApplyMovingObject2D(ZQ_SmokeGrid2D* fineGrid, const std::vector<ZQ_SmokeMovingObject2D*> mvobjs, const ZQ_SmokeSimulationPara2D& para)
{
	int count = mvobjs.size();
	int width = para.width;
	int height = para.height;
	bool*& occupy = fineGrid->GetOccupyPtr();
	float deltah = para.voxelSize;
	float deltat = para.stepTime;
	float*& u = fineGrid->GetVelocityUptr();
	float*& v = fineGrid->GetVelocityVptr();
	float*& density = fineGrid->GetDensityPtr();
	float*& temperature = fineGrid->GetTemperaturePtr();
	memset(occupy,false,sizeof(bool)*width*height);

	for(int i = 0;i < count;i++)
	{
		int cx,cy;
		int posx,posy;
		int sizex,sizey;
		int curu,curv;
		bool* objoccupy = 0;
		if(mvobjs[i] && (objoccupy = (bool*)(mvobjs[i]->GetOccupyPtr())) != 0)
		{
			mvobjs[i]->GetCenter(cx,cy);
			mvobjs[i]->GetPos(posx,posy);
			mvobjs[i]->GetSize(sizex,sizey);
			mvobjs[i]->GetCurVel(curu,curv);
			double scale = deltah / deltat;
			curu *= scale;
			curv *= scale;

			for(int sj = 0;sj < sizey;sj++)
			{
				for(int si = 0;si < sizex;si++)
				{
					if(objoccupy[sj*sizex+si])
					{
						int ix = si+posx-cx;
						int iy = sj+posy-cy;
						if(ix >= 0 && ix < width && iy >= 0 && iy < height)
						{
							occupy[iy*width+ix] = true;
							u[iy*(width+1)+ix] = curu;
							u[iy*(width+1)+(ix+1)] = curu;
							v[iy*width+ix] = curv;
							v[(iy+1)*width+ix] = curv;
						}					
					}
				}
			}
		}
	}
	return true;
}


bool ZQ_AddForce2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;
	int width = fineGrid->GetWidth();
	int height = fineGrid->GetHeight();

	float alpha = para.gravityCoeff;
	float beta  = para.buoyancyCoeff;
	float Tamb  = para.Tamb;
	float deltat = para.stepTime;
	float deltah = para.voxelSize;
	float eta = para.confineCoeff;
	float max_vel = 100.0 * deltah / para.stepTime;
	float* vortVector = new float[width*height];
	memset(vortVector,0,sizeof(float)*width*height);
	float* vortScale  = new float[width*height];
	memset(vortScale,0,sizeof(float)*width*height);
	float* gradVort = new float[width*height*2];
	memset(gradVort,0,sizeof(float)*width*height*2);
	float* u = new float[width*height];
	memset(u,0,sizeof(float)*width*height);
	float* v = new float[width*height];
	memset(v,0,sizeof(float)*width*height);

	float* fconf = new float[width*height*2];
	memset(fconf,0,sizeof(float)*width*height*2);
	float* fbuoy = new float[width*height];
	memset(fbuoy,0,sizeof(float)*width*height);

	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();


	ZQ::ZQ_PoissonSolver::MACtoRegularGrid(width,height,Uptr,Vptr,u,v);

	/*calculate vorticity vector and scale*/
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{

			if(i == 0 || i == height-1 || j == 0 || j == width-1)
			{
			}
			else
			{
				vortVector[i*width+j] = (v[i*width+(j+1)] - v[i*width+(j-1)] - u[(i+1)*width+j] + u[(i-1)*width+j]) / (2*deltah);
				vortScale[i*width+j] = fabs(vortVector[i*width+j]);
			}	
		}
	}


	/*calculate gradient of vorticity scale*/
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{

			if(j == 0)
				gradVort[(i*width+j)*2  ] = vortScale[i*width+(j+1)] - vortScale[i*width+j];
			else if(j == width-1)
				gradVort[(i*width+j)*2  ] = vortScale[i*width+j] - vortScale[i*width+(j-1)];
			else
				gradVort[(i*width+j)*2  ] = 0.5f*(vortScale[i*width+(j+1)] - vortScale[i*width+(j-1)]);

			if(i == 0)
				gradVort[(i*width+j)*2+1] = vortScale[(i+1)*width+j] - vortScale[i*width+j];
			else if(i == height-1)
				gradVort[(i*width+j)*2+1] = vortScale[i*width+j] - vortScale[(i-1)*width+j];
			else
				gradVort[(i*width+j)*2+1] = 0.5f*(vortScale[(i+1)*width+j] - vortScale[(i-1)*width+j]);

			float len = sqrt(gradVort[(i*width+j)*2  ]*gradVort[(i*width+j)*2  ]
			+ gradVort[(i*width+j)*2+1]*gradVort[(i*width+j)*2+1]);
			if(len == 0)
				memset(gradVort+(i*width+j)*2,0,sizeof(float)*2);
			else
			{
				gradVort[(i*width+j)*2  ] /= len;
				gradVort[(i*width+j)*2+1] /= len;
			}
		}
	}

	/*fconf = eta * deltat *(gradVort X vortVector)*/
	float*& temperature = fineGrid->GetTemperaturePtr();
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			fconf[(i*width+j)*2  ] = eta*deltah*(gradVort[(i*width+j)*2+1]*vortVector[i*width+j]);
			fconf[(i*width+j)*2+1] = eta*deltah*(-gradVort[(i*width+j)*2  ]*vortVector[i*width+j]);
			fbuoy[i*width+j] = beta*(temperature[i*width+j]);
		}
	}


	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			fconf[(i*width+j)*2+1] += fbuoy[i*width+j];
		}
	}

	bool*& occupy = fineGrid->GetOccupyPtr();
	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{


			if(occupy[i*width+j] || occupy[i*width+j-1])
				;
			else
				Uptr[i*(width+1)+j] += deltat*0.5f*(fconf[(i*width+j)*2]+fconf[(i*width+j-1)*2]);
		}
	}
	for(int i = 0;i < height;i++)
	{
		if(!occupy[i*width+0])
			Uptr[i*(width+1)+0] += deltat*fconf[(i*width+0)*2];
		if(!occupy[i*width+width-1])
			Uptr[i*(width+1)+width] += deltat*fconf[(i*width+width-1)*2];
	}

	for(int i = 1;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(occupy[i*width+j] || occupy[(i-1)*width+j])
				;
			else
				Vptr[i*width+j] += deltat*0.5f*(fconf[(i*width+j)*2+1]+fconf[((i-1)*width+j)*2+1]);
		}
	}
	for(int j = 0;j < width;j++)
	{
		if(!occupy[0*width+j])
			Vptr[0*width+j] += deltat*fconf[(0*width+j)*2+1];
		if(!occupy[(height-1)*width+j])
			Vptr[height*width+j] += deltat*fconf[((height-1)*width+j)*2+1];
	}

	if(para.globalGridSolver == CLOSED_POISSON || para.globalGridSolver == CLOSED_FLUX 
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON || para.globalGridSolver == CLOSED_OCTREE_FLUX
		|| para.globalGridSolver == COARSE_CLOSED_FLUX_SOR || para.globalGridSolver == COARSE_CLOSED_POISSON_SOR
		|| para.globalGridSolver == CLOSED_OCTREE_POISSON_SOR)
	{
		for(int i = 0;i < height;i++)
		{
			Uptr[i*(width+1)+0] = 0;
			Uptr[i*(width+1)+width] = 0;
		}
		for(int i = 0;i < width;i++)
		{
			Vptr[0*width+i] = 0;
			Vptr[height*width+i] = 0;
		}
	}
	
	delete []vortVector;
	delete []vortScale;
	delete []gradVort;
	delete []u;
	delete []v;
	delete []fconf;
	delete []fbuoy;
	return true;
}

bool ZQ_Attenuation2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;
	float attenVel = exp(-para.stepTime * para.velAttenCoeff);
	float attenDen = exp(-para.stepTime * para.densityAttenCoeff);
	float attenTemp = exp(-para.stepTime * para.tempAttenCoeff);
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	float*& densityPtr = fineGrid->GetDensityPtr();
	float*& temperaturePtr = fineGrid->GetTemperaturePtr();
	int width = para.width;
	int height = para.height;
	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{
			if(occupyPtr[i*width+j] || occupyPtr[i*width+j-1])
				;
			else
				Uptr[i*(width+1)+j] *= attenVel;
		}
		if(!occupyPtr[i*width+0])
			Uptr[i*(width+1)+0] *= attenVel;
		if(!occupyPtr[i*width+width-1])
			Uptr[i*(width+1)+width] *= attenVel;
	}
	
	for(int j = 0;j < width;j++)
	{
		for(int i = 1;i < height;i++)
		{
			if(occupyPtr[i*width+j] || occupyPtr[(i-1)*width+j])
				;
			else
				Vptr[i*width+j] *= attenVel;
		}
		if(!occupyPtr[0*width+j])
			Vptr[0*width+j] *= attenVel;
		if(!occupyPtr[(height-1)*width+j])
			Vptr[height*width+j] *= attenVel;
	}

	for(int i = 0;i < width*height;i++)
	{
		densityPtr[i] *= attenDen;
		temperaturePtr[i] *= attenTemp;
	}
	return true;
}

bool ZQ_VelocityAdvect2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;
	unsigned int width = para.width;
	unsigned int height = para.height;
	unsigned int steps = para.steps;
	float voxelSize = para.voxelSize;
	float deltatt = para.stepTime/steps;
	
	float*& vel_u = fineGrid->GetVelocityUptr();
	float*& vel_v = fineGrid->GetVelocityVptr();
	bool*& occupy = fineGrid->GetOccupyPtr();
	float* input_u = new float[(width+1)*height];
	float* input_v = new float[width*(height+1)];
	float* output_u = new float[(width+1)*height];
	float* output_v = new float[width*(height+1)];

	memcpy(input_u,vel_u,sizeof(float)*(width+1)*height);
	memcpy(input_v,vel_v,sizeof(float)*width*(height+1));
	memset(output_u,0,sizeof(float)*(width+1)*height);
	memset(output_v,0,sizeof(float)*width*(height+1));

	ZQ_CUDA_Advection2D::ZQ_Cuda_Prepare_Advection2D(width,height,voxelSize,steps,deltatt);
	ZQ_CUDA_Advection2D::ZQ_Cuda_Velocity_Advection2D_inMAC_outMAC(vel_u,vel_v,occupy,input_u,input_v,output_u,output_v);

	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{
			if(occupy[i*width+j] || occupy[i*width+j-1])
				;
			else
				vel_u[i*(width+1)+j] = output_u[i*(width+1)+j];
		}
		if(!occupy[i*width+0])
			vel_u[i*(width+1)+0] = output_u[i*(width+1)+0];
		if(!occupy[i*width+width-1])
			vel_u[i*(width+1)+width] = output_u[i*(width+1)+width];
	}
	for(int j = 0;j < width;j++)
	{
		for(int i = 1; i < height;i++)
		{
			if(occupy[i*width+j] || occupy[(i-1)*width+j])
				;
			else
				vel_v[i*width+j] = output_v[i*width+j];
		}
		if(!occupy[0*width+j])
			vel_v[0*width+j] = output_v[0*width+j];
		if(!occupy[(height-1)*width+j])
			vel_v[height*width+j] = output_v[height*width+j];
	}
	
	delete []input_u;
	delete []input_v;
	delete []output_u;
	delete []output_v;

	return true;
}

bool ZQ_VelocityAdvectBFECC2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;
	unsigned int width = para.width;
	unsigned int height = para.height;
	unsigned int steps = para.steps;
	float voxelSize = para.voxelSize;
	float deltatt = para.stepTime/steps;

	float*& vel_u = fineGrid->GetVelocityUptr();
	float*& vel_v = fineGrid->GetVelocityVptr();
	bool*& occupy = fineGrid->GetOccupyPtr();
	float* input_u = new float[(width+1)*height];
	float* input_v = new float[width*(height+1)];
	float* output_u = new float[(width+1)*height];
	float* output_v = new float[width*(height+1)];
	float* tmpoutput_u = new float[(width+1)*height];
	float* tmpoutput_v = new float[width*(height+1)];

	memcpy(input_u,vel_u,sizeof(float)*(width+1)*height);
	memcpy(input_v,vel_v,sizeof(float)*width*(height+1));

	memset(output_u,0,sizeof(float)*(width+1)*height);
	memset(output_v,0,sizeof(float)*width*(height+1));
	memset(tmpoutput_u,0,sizeof(float)*(width+1)*height);
	memset(tmpoutput_v,0,sizeof(float)*width*(height+1));

	
	ZQ_CUDA_Advection2D::ZQ_Cuda_Prepare_Advection2D(width,height,voxelSize,steps,deltatt);
	
	//first advection
	ZQ_CUDA_Advection2D::ZQ_Cuda_Velocity_Advection2D_inMAC_outMAC(vel_u,vel_v,occupy,input_u,input_v,tmpoutput_u,tmpoutput_v);

	//second advection
	for(int i = 0;i < (width+1)*height;i++)
		vel_u[i] = -vel_u[i];
	for(int i = 0;i < width*(height+1);i++)
		vel_v[i] = -vel_v[i];

	ZQ_CUDA_Advection2D::ZQ_Cuda_Velocity_Advection2D_inMAC_outMAC(vel_u,vel_v,occupy,tmpoutput_u,tmpoutput_v,output_u,output_v);
	
	//turn off BFECC close to boundary
	for(int i = 6;i < height-6;i++)
	{
		for(int j = 6;j <= width-6;j++)
		{
			int offset = i*(width+1)+j;
			input_u[offset] = input_u[offset] + 0.5f*(input_u[offset]-output_u[offset]);
		}
	}
	for(int i = 6;i <= height-6;i++)
	{
		for(int j = 6;j < width;j++)
		{
			int offset = i*width+j;
			input_v[offset] = input_v[offset] + 0.5f*(input_v[offset]-output_v[offset]);
		}
	}

	//third advection
	for(int i = 0;i < (width+1)*height;i++)
		vel_u[i] = -vel_u[i];
	for(int i = 0;i < width*(height+1);i++)
		vel_v[i] = -vel_v[i];

	ZQ_CUDA_Advection2D::ZQ_Cuda_Velocity_Advection2D_inMAC_outMAC(vel_u,vel_v,occupy,input_u,input_v,output_u,output_v);

	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{
			if(occupy[i*width+j] || occupy[i*width+j-1])
				;
			else
				vel_u[i*(width+1)+j] = output_u[i*(width+1)+j];
		}
		if(!occupy[i*width+0])
			vel_u[i*(width+1)+0] = output_u[i*(width+1)+0];
		if(!occupy[i*width+width-1])
			vel_u[i*(width+1)+width] = output_u[i*(width+1)+width];
	}
	for(int j = 0;j < width;j++)
	{
		for(int i = 1; i < height;i++)
		{
			if(occupy[i*width+j] || occupy[(i-1)*width+j])
				;
			else
				vel_v[i*width+j] = output_v[i*width+j];
		}
		if(!occupy[0*width+j])
			vel_v[0*width+j] = output_v[0*width+j];
		if(!occupy[(height-1)*width+j])
			vel_v[height*width+j] = output_v[height*width+j];
	}
	
	delete []input_u;
	delete []input_v;
	delete []output_u;
	delete []output_v;
	delete []tmpoutput_u;
	delete []tmpoutput_v;
	return true;
}

bool ZQ_MapFine2Coarse2D(const ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0 || coarseGrid == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int coarseWidth = para.coarseWidth;
	int coarseHeight = para.coarseHeight;
	int subSize = para.subSize;

	const float*& fine_Uptr = fineGrid->GetVelocityUptr();
	const float*& fine_Vptr = fineGrid->GetVelocityVptr();
	float*& coarse_Uptr = coarseGrid->GetVelocityUptr();
	float*& coarse_Vptr = coarseGrid->GetVelocityVptr();

	for(int i = 0;i < coarseHeight;i++)
	{
		for(int j = 0;j <= coarseWidth;j++)
		{
			int i_shift = i*subSize;
			int j_shift = j*subSize;
			float subU = 0;
			for(int si = 0;si < subSize;si++)
				subU += fine_Uptr[(i_shift+si)*(width+1)+j_shift];
			coarse_Uptr[i*(coarseWidth+1)+j] = subU/subSize;
		}
	}

	for(int i = 0;i <= coarseHeight;i++)
	{
		for(int j = 0;j < coarseWidth;j++)
		{
			int i_shift = i*subSize;
			int j_shift = j*subSize;
			float subV = 0;
			for(int sj = 0;sj < subSize;sj++)
				subV += fine_Vptr[i_shift*width+(j_shift+sj)];
			coarse_Vptr[i*coarseWidth+j] = subV/subSize;
		}
	}
	
	return true;
}

bool ZQ_SolveOpenPressure2D(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA)
{
	if(coarseGrid == 0)
		return false;
	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetWidth();
	bool*& occupy = coarseGrid->GetOccupyPtr();
	float*& Uptr = coarseGrid->GetVelocityUptr();
	float*& Vptr = coarseGrid->GetVelocityVptr();
	float*& pressure = coarseGrid->GetPressurePtr();

	int idx = 1;
	int* index = new int[width*height];
	memset(index,0,sizeof(int)*width*height);
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(occupy[i*width+j])
				continue;
			index[i*width+j] = idx++;
		}
	}

	int dim = idx-1;
	
	double* divVelocity = new double[dim];
	memset(divVelocity,0,sizeof(double)*dim);
	double* x = new double[dim];
	memset(x,0,sizeof(double)*dim);
	double* x0 = new double[dim];
	memset(x0,0,sizeof(double)*dim);

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			int row = index[offset] - 1;
			if(row < 0)
				continue;
			divVelocity[row] = Uptr[i*(width+1)+j+1] - Uptr[i*(width+1)+j] + Vptr[(i+1)*width+j] - Vptr[i*width+j];
		}
	}

	
	ZQ::ZQ_PCGSolver solver;
	int max_iter = para.maxIter;
	double tol = 1e-12;
	int it = 0;
	
	solver.PCG_sparse_unsquare(coarseA,divVelocity,x0,max_iter,tol,x,it,false);
	

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(index[i*width+j] != 0)
			{
				pressure[i*width+j] = x[index[i*width+j]-1];
			}
		}
	}
	delete []x0;
	delete []x;
	delete []divVelocity;
	delete []index;

	return true;
}


bool ZQ_SolveOpenPressure2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;

	ZQ_CUDA_PoissonSolver2D::SolveOpenPoissonRedBlackwithOccupy2D_MAC(
		fineGrid->GetVelocityUptr(),
		fineGrid->GetVelocityVptr(),
		fineGrid->GetOccupyPtr(),
		fineGrid->GetWidth(),
		fineGrid->GetHeight(),
		para.maxIter
		);

	/*ZQ::ZQ_PoissonSolver::SolveOpenPoissonSOR_MACGrid(
		fineGrid->GetVelocityUptr(),
		fineGrid->GetVelocityVptr(),
		fineGrid->GetWidth(),
		fineGrid->GetHeight(),
		fineGrid->GetOccupyPtr(),
		para.maxIter,ZQ_FLOAT,false);*/
	return true;
}

bool ZQ_SolveClosedPressure2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int first_x = -1;
	int first_y = -1;
	float div_per_volume = 0.0f;
	int count = 0;

	bool*& occupy = fineGrid->GetOccupyPtr();
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(!occupy[i*width+j])
			{
				if(first_x == -1 && first_y == -1)
				{
					first_x = j;
					first_y = i;
				}
				count ++;

				div_per_volume += Uptr[i*(width+1)+j+1] - Uptr[i*(width+1)+j] + Vptr[(i+1)*width+j] - Vptr[i*width+j];

			}
		}
	}

	if(first_x < 0 && first_y < 0)
		return false;

	if(count == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}
	div_per_volume /= count;


	ZQ_CUDA_PoissonSolver2D::SolveClosedPoissonRedBlackwithOccupy2D_MAC(
		Uptr,
		Vptr,
		occupy,
		//first_x,
		//first_y,
		div_per_volume,
		width,
		height,
		para.maxIter
		);
	
	return true;
}

bool ZQ_SolveOpenFlux2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;


	ZQ_CUDA_PoissonSolver2D::SolveOpenFluxRedBlackwithOccupy2D_MAC(
		fineGrid->GetVelocityUptr(),
		fineGrid->GetVelocityVptr(),
		fineGrid->GetOccupyPtr(),
		fineGrid->GetWidth(),
		fineGrid->GetHeight(),
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);

	
	return true;
}


bool ZQ_SolveClosedFlux2DSOR(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;
	int width = para.width;
	int height = para.height;
	float div_per_volume = 0.0f;

	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	bool*& occupy = fineGrid->GetOccupyPtr();

	int count = 0;
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(!occupy[i*width+j])
			{
				count ++;
				div_per_volume += Uptr[i*(width+1)+j+1] - Uptr[i*(width+1)+j] + Vptr[(i+1)*width+j] - Vptr[i*width+j];
			}
		}
	}

	if(count == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}
	div_per_volume /= count;

	ZQ_CUDA_PoissonSolver2D::SolveClosedFluxRedBlackwithOccupy2D_MAC(
		Uptr,
		Vptr,
		occupy,
		div_per_volume,
		width,
		height,
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);


	return true;
}

bool ZQ_SolveClosedPressure2D(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA)
{
	if(coarseGrid == 0)
		return false;
	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();
	bool*& occupy = coarseGrid->GetOccupyPtr();
	float*& Uptr = coarseGrid->GetVelocityUptr();
	float*& Vptr = coarseGrid->GetVelocityVptr();
	float*& pressure = coarseGrid->GetPressurePtr();

	int idx = 1;
	int* index = new int[width*height];
	memset(index,0,sizeof(int)*width*height);
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(occupy[i*width+j])
				continue;
			index[i*width+j] = idx++;
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

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int row_id = index[i*width+j] - 1;
			if(row_id < 0)
				continue;
			divVelocity[row_id] = Uptr[i*(width+1)+j+1] - Uptr[i*(width+1)+j] + Vptr[(i+1)*width+j] - Vptr[i*width+j];
		}
	}
	
	ZQ::ZQ_PCGSolver solver;
	int max_iter = para.maxIter;
	double tol = 1e-12;
	int it = 0;
	
	solver.PCG_sparse_unsquare(coarseA,divVelocity,x0,max_iter,tol,x,it,false);
	
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			if(index[offset] > 0)
			{
				if(index[offset] == 1)
					pressure[offset] = 0;
				else
					pressure[offset] = x[index[offset]-2];
			}
		}
	}
	delete []x0;
	delete []x;
	delete []divVelocity;
	delete []index;
	return true;
}

bool ZQ_SolveCoarseOpenPoisson2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	bool*& occupy = coarseGrid->GetOccupyPtr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1));

	ZQ_CUDA_PoissonSolver2D::SolveOpenPoissonRedBlackwithFaceRatio2D_MAC(
		deltaU,deltaV,occupy,unoccupyU,unoccupyV,width,height,para.maxIter
		);

	for(int i = 0;i < height*(width+1);i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1);i++)
		deltaV[i] -= coarseV[i];

	return true;
}

bool ZQ_SolveCoarseClosedPoisson2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	float*& unoccupyVolume = coarseGrid->GetVolumeUnOccupyRatioPtr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1));

	float div_per_volume = 0.0f;
	float volumeSize = 0.0f;
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			volumeSize += unoccupyVolume[i*width+j];
			div_per_volume += coarseU[i*(width+1)+j+1] - coarseU[i*(width+1)+j]
							+ coarseV[(i+1)*width+j] - coarseV[i*width+j];
		}
	}
	if(volumeSize == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}
	div_per_volume /= volumeSize;

	ZQ_CUDA_PoissonSolver2D::SolveClosedPoissonRedBlackwithFaceRatio2D_MAC(
		deltaU,deltaV,unoccupyVolume,unoccupyU,unoccupyV,div_per_volume,width,height,para.maxIter
		);

	for(int i = 0;i < height*(width+1);i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1);i++)
		deltaV[i] -= coarseV[i];

	return true;
}

bool ZQ_SolveCoarseOpenFlux2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1));

	ZQ_CUDA_PoissonSolver2D::SolveOpenFluxRedBlackwithFaceRatio2D_MAC(
		deltaU,deltaV,unoccupyU,unoccupyV,width,height,
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);

	for(int i = 0;i < height*(width+1);i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1);i++)
		deltaV[i] -= coarseV[i];

	return true;
}

bool ZQ_SolveCoarseClosedFlux2DSOR(ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaGrid)
{
	int width = para.coarseWidth;
	int height = para.coarseHeight;
	float*& coarseU = coarseGrid->GetVelocityUptr();
	float*& coarseV = coarseGrid->GetVelocityVptr();
	float*& unoccupyVolume = coarseGrid->GetVolumeUnOccupyRatioPtr();
	float*& unoccupyU = coarseGrid->GetFaceUnOccupyRatioUPtr();
	float*& unoccupyV = coarseGrid->GetFaceUnOccupyRatioVPtr();
	float*& deltaU = deltaGrid->GetVelocityUptr();
	float*& deltaV = deltaGrid->GetVelocityVptr();
	memcpy(deltaU,coarseU,sizeof(float)*(width+1)*height);
	memcpy(deltaV,coarseV,sizeof(float)*width*(height+1));

	float div_per_volume = 0.0f;
	float volume = 0.0f;

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			volume += unoccupyVolume[i*width+j];
			div_per_volume += coarseU[i*(width+1)+j+1] - coarseU[i*(width+1)+j]
							+ coarseV[(i+1)*width+j] - coarseV[i*width+j];
		}
	}

	if(volume == 0)
	{
		printf("something is wrong!%s,%d\n",__FILE__,__LINE__);
		return false;
	}
	div_per_volume /= volume;

	ZQ_CUDA_PoissonSolver2D::SolveClosedFluxRedBlackwithFaceRatio2D_MAC(
		deltaU,deltaV,unoccupyVolume,unoccupyU,unoccupyV,div_per_volume,width,height,
		(para.maxIter-1)/para.innerIter+1,
		para.innerIter
		);

	for(int i = 0;i < height*(width+1);i++)
		deltaU[i] -= coarseU[i];
	for(int i = 0;i < width*(height+1);i++)
		deltaV[i] -= coarseV[i];

	return true;
}

bool ZQ_AdjustOpenVelocity2D(const ZQ_SmokeGrid2D* coarseGrid, ZQ_SmokeGrid2D* deltaCoarseGrid)
{
	if(coarseGrid == 0 || deltaCoarseGrid == 0)
		return false;

	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();

	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& pressure = coarseGrid->GetPressurePtr();
	float*& delta_Uptr = deltaCoarseGrid->GetVelocityUptr();
	float*& delta_Vptr = deltaCoarseGrid->GetVelocityVptr();

	memset(delta_Uptr,0,sizeof(float)*(width+1)*height);
	memset(delta_Vptr,0,sizeof(float)*width*(height+1));

	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{
			if(unoccupyUptr[i*(width+1)+j] > 0)
				delta_Uptr[i*(width+1)+j] = -(pressure[i*width+j]-pressure[i*width+j-1]);
		}
		if(unoccupyUptr[i*(width+1)+0] > 0)
			delta_Uptr[i*(width+1)+0] = -(pressure[i*width+0]-0);
		if(unoccupyUptr[i*(width+1)+width] > 0)
			delta_Uptr[i*(width+1)+width] = -(0-pressure[i*width+width-1]);
	}

	for(int j = 0;j < width;j++)
	{
		for(int i = 1;i < height;i++)
		{
			if(unoccupyVptr[i*width+j] > 0)
				delta_Vptr[i*width+j] = -(pressure[i*width+j]-pressure[(i-1)*width+j]);
		}
		if(unoccupyVptr[0*width+j] > 0)
			delta_Vptr[0*width+j] = -(pressure[0*width+j]-0);
		if(unoccupyVptr[height*width+j] > 0)
			delta_Vptr[height*width+j] = -(0-pressure[(height-1)*width+j]);
	}

	return true;
}

bool ZQ_AdjustClosedVelocity2D(const ZQ_SmokeGrid2D* coarseGrid, ZQ_SmokeGrid2D* deltaCoarseGrid)
{
	if(coarseGrid == 0 || deltaCoarseGrid == 0)
		return false;

	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();

	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const float*& pressure = coarseGrid->GetPressurePtr();
	float*& delta_Uptr = deltaCoarseGrid->GetVelocityUptr();
	float*& delta_Vptr = deltaCoarseGrid->GetVelocityVptr();

	memset(delta_Uptr,0,sizeof(float)*(width+1)*height);
	memset(delta_Vptr,0,sizeof(float)*width*(height+1));

	for(int i = 0;i < height;i++)
	{
		for(int j = 1;j < width;j++)
		{
			if(unoccupyUptr[i*(width+1)+j] > 0)
				delta_Uptr[i*(width+1)+j] = -(pressure[i*width+j]-pressure[i*width+j-1]);
		}
	}

	for(int j = 0;j < width;j++)
	{
		for(int i = 1;i < height;i++)
		{
			if(unoccupyVptr[i*width+j] > 0)
				delta_Vptr[i*width+j] = -(pressure[i*width+j]-pressure[(i-1)*width+j]);
		}
	}

	return true;
}

bool ZQ_SolveFlux2D(const ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA, ZQ_SmokeGrid2D* deltaCoarseGrid)
{
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();

	const float*& unoccupyUptr = coarseGrid->GetFaceUnOccupyRatioUPtr();
	const float*& unoccupyVptr = coarseGrid->GetFaceUnOccupyRatioVPtr();
	const bool*& occupy = coarseGrid->GetOccupyPtr();
	const float*& Uptr = coarseGrid->GetVelocityUptr();
	const float*& Vptr = coarseGrid->GetVelocityVptr();
	float*& delta_Uptr = deltaCoarseGrid->GetVelocityUptr();
	float*& delta_Vptr = deltaCoarseGrid->GetVelocityVptr();

	int* u_index = new int[height*(width+1)];
	int* v_index = new int[width*(height+1)];
	int* lamda_index = new int[width*height];
	int k_index = 0;

	for(int i = 0;i < height*(width+1);i++)
		u_index[i] = -1;
	for(int i = 0;i < width*(height+1);i++)
		v_index[i] = -1;
	for(int i = 0;i < width*height;i++)
		lamda_index[i] = -1;


	int idx = 0;
	if(para.globalGridSolver == OPEN_FLUX)
	{
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j <= width;j++)
			{
				if(unoccupyUptr[i*(width+1)+j] != 0)
					u_index[i*(width+1)+j] = idx++;
			}
		}
		for(int i = 0;i <= height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(unoccupyVptr[i*width+j] != 0)
					v_index[i*width+j] = idx++;
			}
		}
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(!occupy[i*width+j])
					lamda_index[i*width+j] = idx++;
			}
		}
	}
	else
	{
		for(int i = 0;i < height;i++)
		{
			for(int j = 1;j < width;j++)
			{
				if(unoccupyUptr[i*(width+1)+j] != 0)
					u_index[i*(width+1)+j] = idx++;
			}
		}
		for(int i = 1;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(unoccupyVptr[i*width+j] != 0)
					v_index[i*width+j] = idx++;
			}
		}
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(!occupy[i*width+j])
					lamda_index[i*width+j] = idx++;
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
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			if(lamda_index[i*width+j] >= 0)
			{
				b[lamda_index[i*width+j]] = Uptr[i*(width+1)+j+1] - Uptr[i*(width+1)+j] + Vptr[(i+1)*width+j] - Vptr[i*width+j];
			}
		}
	}

	ZQ::ZQ_PCGSolver solver;
	int it = 0;
	double tol = 1e-14;
	solver.PCG_sparse_unsquare(coarseA,b,x0,para.maxIter,tol,x,it,false);

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j <= width;j++)
		{
			int offset = i*(width+1)+j;
			delta_Uptr[offset] = u_index[offset] >= 0 ? x[u_index[offset]] : 0;
		}
	}
	for(int i = 0;i <= height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			delta_Vptr[offset] = v_index[offset] >= 0 ? x[v_index[offset]] : 0;
		}
	}
	
	delete []u_index;
	delete []v_index;
	delete []lamda_index;
	delete []x0;
	delete []x;
	delete []b;
	return true;
}

bool ZQ_ApplyGuidingVelocity2D(const ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeGuideGrid2D* guideGrid, const ZQ_SmokeSimulationPara2D& para, ZQ_SmokeGrid2D* deltaCoarseGrid)
{
	if( coarseGrid == 0 || guideGrid == 0 || deltaCoarseGrid == 0)
		return false;
	if(!para.useGuidingControl)
		return false;

	float guideCoeff = para.guidingCoeff;
	
	int width = coarseGrid->GetWidth();
	int height = coarseGrid->GetHeight();

	const float*& Uptr = coarseGrid->GetVelocityUptr();
	const float*& Vptr = coarseGrid->GetVelocityVptr();
	const float*& GUptr = guideGrid->GetUptr();
	const float*& GVptr = guideGrid->GetVptr();
	float*& deltaU = deltaCoarseGrid->GetVelocityUptr();
	float*& deltaV = deltaCoarseGrid->GetVelocityVptr();

	for(int i = 0;i < height*(width+1);i++)
		deltaU[i] += guideCoeff*(GUptr[i] - Uptr[i]);
	for(int i = 0;i < width*(height+1);i++)
		deltaV[i] += guideCoeff*(GVptr[i] - Vptr[i]);

	return true;
}

bool ZQ_MapCoarse2Fine2D(const ZQ_SmokeGrid2D* deltaCoarseGrid, ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(deltaCoarseGrid == 0 || fineGrid == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int subSize = para.subSize;
	int coarseWidth = para.coarseWidth;
	int coarseHeight = para.coarseHeight;
	
	
	/*adjust the fineGrid velocity by plus the deltaCoarseGrid*/
	float*& fineUptr = fineGrid->GetVelocityUptr();
	float*& fineVptr = fineGrid->GetVelocityVptr();
	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	const float*& deltaUptr = deltaCoarseGrid->GetVelocityUptr();
	const float*& deltaVptr = deltaCoarseGrid->GetVelocityVptr();

	for(int i = 0;i < coarseHeight;i++)
	{
		int i_shift = i*subSize;
		for(int j = 1;j < coarseWidth;j++)
		{
			int j_shift = j*subSize;
			for(int si = 0;si < subSize;si++)
			{
				if(!occupyPtr[(i_shift+si)*width+j_shift] && !occupyPtr[(i_shift+si)*width+j_shift-1])
					fineUptr[(i_shift+si)*(width+1)+j_shift] += deltaUptr[i*(coarseWidth+1)+j];
			}
		}
		int j_shift0 = 0;
		int j_shiftN = coarseWidth*subSize;
		for(int si = 0;si < subSize;si++)
		{
			if(!occupyPtr[(i_shift+si)*width+j_shift0])
				fineUptr[(i_shift+si)*(width+1)+j_shift0] += deltaUptr[i*(coarseWidth+1)+0];
			if(!occupyPtr[(i_shift+si)*width+j_shiftN-1])
				fineUptr[(i_shift+si)*(width+1)+j_shiftN] += deltaUptr[i*(coarseWidth+1)+coarseWidth];
		}
	}
	
	for(int j = 0;j < coarseWidth;j++)
	{
		int j_shift = j*subSize;
		for(int i = 1;i < coarseHeight;i++)
		{
			int i_shift = i*subSize;
			for(int sj = 0;sj < subSize;sj++)
			{
				if(!occupyPtr[i_shift*width+(j_shift+sj)] && !occupyPtr[(i_shift-1)*width+(j_shift+sj)])
					fineVptr[i_shift*width+(j_shift+sj)] += deltaVptr[i*coarseWidth+j];
			}
		}
		int i_shift0 = 0;
		int i_shiftN = coarseHeight*subSize;
		for(int sj = 0;sj < subSize;sj++)
		{
			if(!occupyPtr[i_shift0*width+(j_shift+sj)])
				fineVptr[i_shift0*width+(j_shift+sj)] += deltaVptr[0*coarseWidth+j];
			if(!occupyPtr[(i_shiftN-1)*width+(j_shift+sj)])
				fineVptr[i_shiftN*width+(j_shift+sj)] += deltaVptr[coarseHeight*coarseWidth+j];
		}
	}

	return true;
}

bool ZQ_SubProjection2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para, const float* subMatrix, int subRow, int subCol, bool useGPUflag)
{
	if(fineGrid == 0 || subMatrix == 0)
		return false;

	int width = para.width;
	int height = para.height;
	int subSize = para.subSize;
	int coarseWidth = para.coarseWidth;
	int coarseHeight = para.coarseHeight;

	bool*& occupyPtr = fineGrid->GetOccupyPtr();
	float*& Uptr = fineGrid->GetVelocityUptr();
	float*& Vptr = fineGrid->GetVelocityVptr();
	if(useGPUflag)
	{
		bool* regular_flag = new bool[coarseWidth*coarseHeight];
		memset(regular_flag,0,sizeof(bool)*coarseWidth*coarseHeight);
		int fast_count = 0;
		
		for(int i = 0;i < coarseHeight;i++)
		{
			for(int j = 0;j < coarseWidth;j++)
			{
				for(int si = 0;si < subSize;si++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						if(occupyPtr[(i*subSize+si)*width+(j*subSize+sj)])
							regular_flag[i*coarseWidth+j] = true;
					}
				}
				if(!regular_flag[i*coarseWidth+j])
					fast_count++;
			}
		}
		if(fast_count > 0)
		{
			float* input = new float[subCol*fast_count];
			memset(input,0,sizeof(float)*subCol*fast_count);
			float* output = new float[subRow*fast_count];
			memset(output,0,sizeof(float)*subRow*fast_count);
			int cur_count = 0;
			for(int i = 0;i < coarseHeight;i++)
			{
				for(int j = 0;j < coarseWidth;j++)
				{
					if(!regular_flag[i*coarseWidth+j])
					{
						int cur_row = 0;
						for(int si = 0;si < subSize;si++)
						{
							for(int sj = 0;sj <= subSize;sj++)
							{
								input[cur_row*fast_count+cur_count] = Uptr[(i*subSize+si)*(width+1)+(j*subSize+sj)];
								cur_row ++;
							}
						}
						for(int si = 0;si <= subSize;si++)
						{
							for(int sj = 0;sj < subSize;sj++)
							{
								input[cur_row*fast_count+cur_count] = Vptr[(i*subSize+si)*width+(j*subSize+sj)];
								cur_row ++;
							}
						}
						cur_count ++;
					}
				}
			}
			
			SubProjectionCuda(subRow,subCol,fast_count,subMatrix,input,output);
			
			cur_count = 0;
			for(int i = 0;i < coarseHeight;i++)
			{
				for(int j = 0;j < coarseWidth;j++)
				{
					if(!regular_flag[i*coarseWidth+j])
					{
						int cur_row = 0;
						for(int si = 0; si < subSize;si++)
						{
							for(int sj = 1; sj < subSize; sj++)
							{
								Uptr[(i*subSize+si)*(width+1)+(j*subSize+sj)] = output[cur_row*fast_count+cur_count];
								cur_row++;
							}
						}
						for(int si = 1;si < subSize;si++)
						{
							for(int sj = 0; sj < subSize;sj++)
							{
								Vptr[(i*subSize+si)*width+(j*subSize+sj)] = output[cur_row*fast_count+cur_count];
								cur_row++;
							}
						}
						cur_count++;
					}
				}
			}
			delete []input;
			delete []output;
		}
		for(int i = 0;i < coarseHeight;i++)
		{
			for(int j = 0;j < coarseWidth;j++)
			{
				if(regular_flag[i*coarseWidth+j])
				{
					float* sub_input = new float[subSize*(subSize+1)*2];
					memset(sub_input,0,sizeof(float)*subSize*(subSize+1)*2);
					float* sub_output = new float[subSize*(subSize+1)*2];
					memset(sub_output,0,sizeof(float)*subSize*(subSize+1)*2);
					bool* sub_occupy = new bool[subSize*subSize];
					memset(sub_occupy,0,sizeof(bool)*subSize*subSize);

					int i_shift = i*subSize;
					int j_shift = j*subSize;
					int v_shift = subSize*(subSize+1);

					for(int si = 0;si < subSize;si++)
					{
						for(int sj = 0;sj <= subSize;sj++)
						{
							sub_input[si*(subSize+1)+sj] = Uptr[(i_shift+si)*(width+1)+(j_shift+sj)];
						}
					}
					
					for(int si = 0;si <= subSize;si++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							sub_input[v_shift+si*subSize+sj] = Vptr[(i_shift+si)*width+(j_shift+sj)];
						}
					}
					for(int si = 0;si < subSize;si++)
					{
						for(int sj = 0;sj < subSize;sj++)
						{
							sub_occupy[si*subSize+sj] = occupyPtr[(i_shift+si)*width+(j_shift+sj)];
						}
					}
					memcpy(sub_output,sub_input,sizeof(float)*subSize*(subSize+1)*2);
					ZQ_SolveSubGrid2D(para,sub_occupy,sub_input,sub_output);
					for(int si = 0;si < subSize;si++)
					{
						for(int sj = 1;sj < subSize;sj++)
						{
							Uptr[(i_shift+si)*(width+1)+(j_shift+sj)] = sub_output[si*(subSize+1)+sj];
						}
					}
					for(int sj = 0;sj < subSize;sj++)
					{
						for(int si = 1;si < subSize;si++)
						{
							Vptr[(i_shift+si)*width+(j_shift+sj)] = sub_output[v_shift+si*subSize+sj];
						}
					}
					delete []sub_input;
					delete []sub_output;
					delete []sub_occupy;
				}
			}
		}
		delete []regular_flag;
		regular_flag = 0;
	}
	else
	{
		for(int i = 0;i < coarseHeight;i++)
		{
			for(int j = 0;j < coarseWidth;j++)
			{
				float* sub_input = new float[subSize*(subSize+1)*2];
				memset(sub_input,0,sizeof(float)*subSize*(subSize+1)*2);
				float* sub_output = new float[subSize*(subSize+1)*2];
				memset(sub_output,0,sizeof(float)*subSize*(subSize+1)*2);
				bool* sub_occupy = new bool[subSize*subSize];
				memset(sub_occupy,0,sizeof(bool)*subSize*subSize);

				int i_shift = i*subSize;
				int j_shift = j*subSize;
				int v_shift = subSize*(subSize+1);

				for(int si = 0;si < subSize;si++)
				{
					for(int sj = 0;sj <= subSize;sj++)
					{
						sub_input[si*(subSize+1)+sj] = Uptr[(i_shift+si)*(width+1)+(j_shift+sj)];
					}
				}

				for(int si = 0;si <= subSize;si++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						sub_input[v_shift+si*subSize+sj] = Vptr[(i_shift+si)*width+(j_shift+sj)];
					}
				}
				for(int si = 0;si < subSize;si++)
				{
					for(int sj = 0;sj < subSize;sj++)
					{
						sub_occupy[si*subSize+sj] = occupyPtr[(i_shift+si)*width+(j_shift+sj)];
					}
				}
				memcpy(sub_output,sub_input,sizeof(float)*subSize*(subSize+1)*2);
				ZQ_SolveSubGrid2D(para,sub_occupy,sub_input,sub_output);
				for(int si = 0;si < subSize;si++)
				{
					for(int sj = 1;sj < subSize;sj++)
					{
						Uptr[(i_shift+si)*(width+1)+(j_shift+sj)] = sub_output[si*(subSize+1)+sj];
					}
				}
				for(int sj = 0;sj < subSize;sj++)
				{
					for(int si = 1;si < subSize;si++)
					{
						Vptr[(i_shift+si)*width+(j_shift+sj)] = sub_output[v_shift+si*subSize+sj];
					}
				}
				delete []sub_input;
				delete []sub_output;
				delete []sub_occupy;
			}
		}

	}
	return true;
}

bool ZQ_TotalProjection2D(ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA,
								const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid2D* guideGrid)
{

	LARGE_INTEGER tf,solve_global_st, solve_global_end;
	QueryPerformanceFrequency(&tf);
	QueryPerformanceCounter(&solve_global_st);
	
	if(!ZQ_MapFine2Coarse2D(fineGrid,coarseGrid,para))
		return false;
	ZQ_SmokeGrid2D* deltaGrid = new ZQ_SmokeGrid2D(para.coarseWidth,para.coarseHeight);

	if(para.globalGridSolver == CLOSED_POISSON)
	{
		if(!ZQ_SolveClosedPressure2D(coarseGrid,para,coarseA))
			return false;
		if(!ZQ_AdjustClosedVelocity2D(coarseGrid,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == OPEN_POISSON)
	{
		
		if(!ZQ_SolveOpenPressure2D(coarseGrid,para,coarseA))
			return false;
		
		if(!ZQ_AdjustOpenVelocity2D(coarseGrid,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
		
	}
	else if(para.globalGridSolver == OPEN_FLUX || para.globalGridSolver == CLOSED_FLUX)
	{
		if(!ZQ_SolveFlux2D(coarseGrid,para,coarseA,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_OPEN_POISSON_SOR)
	{
		if(!ZQ_SolveCoarseOpenPoisson2DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_CLOSED_POISSON_SOR)
	{
		if(!ZQ_SolveCoarseOpenPoisson2DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_OPEN_FLUX_SOR)
	{
		if(!ZQ_SolveCoarseOpenFlux2DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	else if(para.globalGridSolver == COARSE_CLOSED_FLUX_SOR)
	{
		if(!ZQ_SolveCoarseClosedFlux2DSOR(coarseGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}

	if(para.useGuidingControl)
	{
		if(!ZQ_ApplyGuidingVelocity2D(coarseGrid,guideGrid,para,deltaGrid))
		{
			delete deltaGrid;
			return false;
		}
	}
	
	if(!ZQ_MapCoarse2Fine2D(deltaGrid,fineGrid,para))
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
		if(!ZQ_SubProjection2D(fineGrid,para,subMatrix,subRow,subCol,true))
		{
			return false;
		}
		QueryPerformanceCounter(&solve_sub_end);
		printf("solve sub cost:%f\n",1.0*(solve_sub_end.QuadPart-solve_sub_st.QuadPart)/tf.QuadPart);
	}
	return true;
}

bool ZQ_DensityTemperatureAdvect2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;
	
	unsigned int width = para.width;
	unsigned int height = para.height;
	unsigned int steps = para.steps;
	float voxelSize = para.voxelSize;
	float deltatt = para.stepTime/steps;

	float* inputTempDen = new float[width*height*2];
	float* outputTempDen = new float[width*height*2];
	memset(inputTempDen,0,sizeof(float)*width*height*2);
	memset(outputTempDen,0,sizeof(float)*width*height*2);
	
	bool*& occupy = fineGrid->GetOccupyPtr();
	float*& u = fineGrid->GetVelocityUptr();
	float*& v = fineGrid->GetVelocityVptr();

	float*& temperaturePtr = fineGrid->GetTemperaturePtr();
	float*& densityPtr = fineGrid->GetDensityPtr();
	
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			inputTempDen[offset*2  ] = temperaturePtr[offset];
			inputTempDen[offset*2+1] = densityPtr[offset];
		}
	}
	
	ZQ_CUDA_Advection2D::ZQ_Cuda_Prepare_Advection2D(width,height,voxelSize,steps,deltatt);
	ZQ_CUDA_Advection2D::ZQ_Cuda_Scalar_Advection2D_MAC_Velocity(u,v,occupy,inputTempDen,outputTempDen);
	
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			
			temperaturePtr[offset] = outputTempDen[offset*2];
			densityPtr[offset] = outputTempDen[offset*2+1];
		}
	}
	
	delete []inputTempDen;
	delete []outputTempDen;
	return true;
}

bool ZQ_DensityTemperatureAdvectBFECC2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	if(fineGrid == 0)
		return false;

	unsigned int width = para.width;
	unsigned int height = para.height;
	unsigned int steps = para.steps;
	float voxelSize = para.voxelSize;
	float deltatt = para.stepTime/steps;

	float* inputTempDen = new float[width*height*2];
	float* outputTempDen = new float[width*height*2];
	float* tempTempDen = new float[width*height*2];
	memset(inputTempDen,0,sizeof(float)*width*height*2);
	memset(outputTempDen,0,sizeof(float)*width*height*2);
	memset(tempTempDen,0,sizeof(float)*width*height*2);

	bool*& occupy = fineGrid->GetOccupyPtr();
	float*& u = fineGrid->GetVelocityUptr();
	float*& v = fineGrid->GetVelocityVptr();

	float*& temperaturePtr = fineGrid->GetTemperaturePtr();
	float*& densityPtr = fineGrid->GetDensityPtr();

	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			inputTempDen[offset*2  ] = temperaturePtr[offset];
			inputTempDen[offset*2+1] = densityPtr[offset];
		}
	}

	
	ZQ_CUDA_Advection2D::ZQ_Cuda_Prepare_Advection2D(width,height,voxelSize,steps,deltatt);

	//first advection
	ZQ_CUDA_Advection2D::ZQ_Cuda_Scalar_Advection2D_MAC_Velocity(u,v,occupy,inputTempDen,tempTempDen);

	// second advection

	for(int i = 0;i < height*(width+1);i++)
		u[i] = -u[i];
	for(int i = 0;i < width*(height+1);i++)
		v[i] = -v[i];
 
	ZQ_CUDA_Advection2D::ZQ_Cuda_Scalar_Advection2D_MAC_Velocity(u,v,occupy,inputTempDen,outputTempDen);
	
	//turn off BFECC close to the boundary
	for(int i = 6;i < height-6;i++)
	{
		for(int j = 6;j < width-6;j++)
		{
			int offset = i*width+j;
			inputTempDen[offset*2] = inputTempDen[offset*2] + 0.5f*(inputTempDen[offset*2] - outputTempDen[offset*2]);
			inputTempDen[offset*2+1] = inputTempDen[offset*2+1] + 0.5f*(inputTempDen[offset*2+1] - outputTempDen[offset*2+1]);
			if(inputTempDen[offset*2] < 0)
				inputTempDen[offset*2] = 0;
			else if(inputTempDen[offset*2] > para.maxTemperature)
				inputTempDen[offset*2] = para.maxTemperature;
			if(inputTempDen[offset*2+1] < 0)
				inputTempDen[offset*2+1] = 0;
			else if(inputTempDen[offset*2+1] > para.maxDensity)
				inputTempDen[offset*2+1] = para.maxDensity;
		}
	}

	for(int i = 0;i < height*(width+1);i++)
		u[i] = -u[i];
	for(int i = 0;i < width*(height+1);i++)
		v[i] = -v[i];

	//third advection
	ZQ_CUDA_Advection2D::ZQ_Cuda_Scalar_Advection2D_MAC_Velocity(u,v,occupy,inputTempDen,outputTempDen);
	
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{
			int offset = i*width+j;
			temperaturePtr[offset] = outputTempDen[offset*2];
			densityPtr[offset] = outputTempDen[offset*2+1];
		}
	}
	
	delete []inputTempDen;
	delete []outputTempDen;
	delete []tempTempDen;
	return true;
}

bool ZQ_ReinjectTemperatureDensity2D(ZQ_SmokeGrid2D* fineGrid, const ZQ_SmokeSimulationPara2D& para)
{
	int width = para.width;
	int height = para.height;

	float*& density = fineGrid->GetDensityPtr();
	float*& temperature = fineGrid->GetTemperaturePtr();

	bool randDenFlag = para.randDensity;
	
	for(int cc = 0;cc < para.sources.size();cc++)
	{
		float cx = para.sources[cc].cx * width;
		float cy = para.sources[cc].cy * height;
		float half_xlen = para.sources[cc].half_xlen * width;
		float half_ylen = para.sources[cc].half_ylen * height;
		float reinject_den = para.sources[cc].reinjectDensity;
		float reinject_tem = para.sources[cc].reinjectTemperature;
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(fabs(j+0.5f-cx) <= half_xlen && fabs(i+0.5f-cy) <= half_ylen)
				{
					density[i*width+j] = reinject_den * (randDenFlag ? rand()%256/255.0f : 1.0f);
					temperature[i*width+j] = reinject_tem;
				}
			}
		}
	}
	
	return true;

}

bool ZQ_UpdateOneFrame2D(ZQ_SmokeGrid2D* fineGrid, ZQ_SmokeGrid2D* coarseGrid, const ZQ_SmokeSimulationPara2D& para, const taucs_ccs_matrix* coarseA,
								const float* subMatrix, int subRow, int subCol, const ZQ_SmokeGuideGrid2D* guideGrid)
{
	if(!ZQ_AddForce2D(fineGrid,para))
		return false;
	if(!ZQ_Attenuation2D(fineGrid,para))
		return false;

	if(!ZQ_VelocityAdvect2D(fineGrid,para))
		return false;

	if(!ZQ_TotalProjection2D(fineGrid,coarseGrid,para,coarseA,subMatrix,subRow,subCol,guideGrid))
		return false;

	if(!ZQ_DensityTemperatureAdvect2D(fineGrid,para))
		return false;

	if(!ZQ_ReinjectTemperatureDensity2D(fineGrid,para))
		return false;
	return true;
}