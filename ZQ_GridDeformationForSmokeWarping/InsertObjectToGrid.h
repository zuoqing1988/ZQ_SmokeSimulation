#ifndef _INSERT_OBJECT_TO_GRID_H_
#define _INSERT_OBJECT_TO_GRID_H_
#pragma once


#include <string.h>
#include <vector>
#include "ZQ_GridDeformation.h"

class InsertObjectToGridOptions
{
public:
	InsertObjectToGridOptions()
	{
		Reset();
	}
	~InsertObjectToGridOptions()
	{

	}
	enum CONST_VAL{FILE_LEN = 200};
	int max_iteration;
	int FPiteration;
	float line_weight;
	float angle_weight;
	float distance_weight;
	bool has_input_file;
	char input_file[FILE_LEN];
	bool has_map_file;
	char map_file[FILE_LEN];
	bool has_show_file;
	char show_file[FILE_LEN];

	void Reset()
	{
		max_iteration = 10000;
		FPiteration = 5;
		line_weight = 1;
		angle_weight = 0.1;
		distance_weight = 1;
		has_input_file = false;
		input_file[0] = '\0';
		has_map_file = false;
		map_file[0] = '\0';
		has_show_file = false;
		show_file[0] = '\0';
	}
	bool HandleArgs(int argc, const char** argv)
	{
		for(int k = 0;k < argc;k++)
		{
			if(_strcmpi(argv[k],"max_iteration") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of max_iteration ?\n");
					return false;
				}
				max_iteration = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"FPiteration") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of FPiteration ?\n");
					return false;
				}
				FPiteration = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"line_weight") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of line_weight ?\n");
					return false;
				}
				line_weight = atof(argv[k]);
			}
			else if(_strcmpi(argv[k],"angle_weight") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of angle_weight ?\n");
					return false;
				}
				angle_weight = atof(argv[k]);
			}
			else if(_strcmpi(argv[k],"distance_weight") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of distance_weight ?\n");
					return false;
				}
				distance_weight = atof(argv[k]);
			}
			else if(_strcmpi(argv[k],"input_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of input_file ?\n");
					return false;
				}
				strcpy(input_file,argv[k]);
				has_input_file = true;
			}
			else if(_strcmpi(argv[k],"map_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of map_file ?\n");
					return false;
				}
				strcpy(map_file,argv[k]);
				has_map_file = true;
			}
			else if(_strcmpi(argv[k],"show_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of show_file ?\n");
					return false;
				}
				strcpy(show_file,argv[k]);
				has_show_file = true;
			}
			else
			{
				printf("unknown parameter: %s\n",argv[k]);
				return false;
			}
		}
		return true;
	}
};

class InsertObjectToGrid
{
	class CtrlPt
	{
	public:
		CtrlPt(int _ix, int _iy, float _fx, float _fy):ix(_ix),iy(_iy),fx(_fx),fy(_fy){}
		int ix,iy;
		float fx,fy;
	};

public:
	InsertObjectToGrid(const int width, const int height)
	{
		this->width = width;
		this->height = height;
		objMask = new bool[width*height];
		contourMask = new bool[width*height];
		outCoord = new float[width*height*2];
		max_iteration = 1000;
		FPiteration = 5;
		line_weight = 1;
		angle_weight = 0.1;
		distance_weight = 0;
		memset(objMask,0,sizeof(bool)*width*height);
		memset(contourMask,0,sizeof(bool)*width*height);
		memset(outCoord,0,sizeof(float)*width*height*2);

		show_scale = 5;
		show_shift = 2;
		background = cvCreateImage(cvSize(width*show_scale,height*show_scale),IPL_DEPTH_8U,3);
		nSelected = 0;
		has_split_point = false;

		color_non_object = cvScalar(0,200,0);
		color_object = cvScalar(200,0,0);
		color_contour = cvScalar(0,0,200);
		color_selected = cvScalar(0,100,250);
		color_split = cvScalar(250,0,250);
		thick_line = 1;
		thick_point = 2;
		size_point = 1;
	}
	~InsertObjectToGrid()
	{
		if(objMask)
		{
			delete []objMask;
			objMask = 0;
		}
		if(contourMask)
		{
			delete []contourMask;
			contourMask = 0;
		}
		if(outCoord)
		{
			delete []outCoord;
			outCoord = 0;
		}
		cvReleaseImage(&background);
	}

private:
	int width,height;
	bool* objMask;
	bool* contourMask;
	float* outCoord;
	int max_iteration;
	int FPiteration;
	float line_weight;
	float angle_weight;
	float distance_weight;

	std::vector<int> final_contour_x;
	std::vector<int> final_contour_y;

	//for dialog
	CvScalar color_non_object/* = cvScalar(0,200,0)*/;
	CvScalar color_object/* = cvScalar(200,0,0)*/;
	CvScalar color_contour/* = cvScalar(0,0,200)*/;
	CvScalar color_selected/* = cvScalar(0,100,150)*/;
	CvScalar color_split/* = cvScalar(150,0,150)*/;
	int thick_line/* = 1*/;
	int thick_point/* = 1*/;
	int size_point/* = 1*/;

	//
	int show_scale;
	int show_shift;
	IplImage* background;
	int nSelected;
	enum CONST_VAL
	{
		NMAX_SELECTED = 2
	};
	int selectedX[NMAX_SELECTED];
	int selectedY[NMAX_SELECTED];
	bool has_split_point;
	int splitIdx[NMAX_SELECTED];
	float splitX[NMAX_SELECTED];
	float splitY[NMAX_SELECTED];



public:
	const float* GetOutCoord() const{return outCoord;}
	void SetMaxIteration(const int max_it){max_iteration = max_it;}
	void SetFPIteration(const int fp_it){FPiteration = fp_it;}
	void SetLineWeight(const float line_w){line_weight = line_w;}
	void SetAngleWeight(const float angle_w){angle_weight = angle_w;}
	void SetDistanceWeight(const float dist_w){distance_weight = dist_w;}
	void SetObjectMask(bool* mask)
	{
		memcpy(objMask,mask,sizeof(bool)*width*height);
	}

	bool ShowDialog()
	{
		_objMask2ContourMask(width,height,objMask,contourMask);
		std::vector<int> contour_x,contour_y;
		if(!_contourMask2Contour(width,height,contourMask,contour_x,contour_y))
		{
			printf("failed to extract contour\n");
			return false;
		}
		if(!_reduceContour(width,height,objMask,contour_x,contour_y,final_contour_x,final_contour_y))
		{
			printf("failed to reduce contour\n");
			return false;
		}

		/*******************/
	
		cvZero(background);
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(objMask[i*width+j])
				{
					cvCircle(background,cvPoint(j*show_scale+show_shift,i*show_scale+show_shift),size_point,color_object,thick_point);
				}
				else
				{
					cvCircle(background,cvPoint(j*show_scale+show_shift,i*show_scale+show_shift),size_point,color_non_object,thick_point);
				}
			}
		}

		for(int i = 0;i < contour_x.size()-1;i++)
		{
			float x1 = contour_x[i];
			float y1 = contour_y[i];
			float x2 = contour_x[i+1];
			float y2 = contour_y[i+1];
			cvLine(background,cvPoint(x1*show_scale+show_shift,y1*show_scale+show_shift),cvPoint(x2*show_scale+show_shift,y2*show_scale+show_shift),color_contour,thick_line);
		}

		const char* dialog_name = "Dialog";

		cvNamedWindow(dialog_name);
		cvSetMouseCallback(dialog_name,_mouseHandler,this);

		bool flag = false;

		do 
		{
			IplImage* show_img = cvCloneImage(background);
			_drawToImage(show_img);
			cvShowImage(dialog_name,show_img);
			cvReleaseImage(&show_img);
			int key = cvWaitKey(10);
			if(key == 'q' || key == 'Q')
			{
				flag = false;
				break;
			}
			else if(key == 'o' || key == 'O')
			{
				flag = true;
				
				break;
			}
		} while (true);
		
		cvDestroyWindow(dialog_name);
		if(!flag)
		{
			return false;
		}
		
		return _deform();
	}

private:

	bool _map_contour(const std::vector<int>& contour_x, const std::vector<int>& contour_y, const float pt1[2], const float pt2[2], std::vector<CtrlPt>& ctrl_pts)
	{
		int contour_size = contour_x.size();
		if(contour_size < 2 || contour_size != contour_y.size())
			return false;
		float* len_shift = new float[contour_size];
		float total_len = 0;
		for(int i = 0;i < contour_size-1;i++)
		{
			len_shift[i] = total_len;
			float cur_x = contour_x[i+1] - contour_x[i];
			float cur_y = contour_y[i+1] - contour_y[i];
			float cur_len = sqrt(cur_x*cur_x+cur_y*cur_y);
			total_len += cur_len;
		}
		len_shift[contour_size-1] = total_len;
		if(total_len == 0)
		{
			delete []len_shift;
			return false;
		}

		float dir[2] = {pt2[0]-pt1[0],pt2[1]-pt1[1]};
		ctrl_pts.clear();
		for(int i = 0;i < contour_size;i++)
		{
			float weight = len_shift[i]/total_len;
			ctrl_pts.push_back(CtrlPt(contour_x[i],contour_y[i],pt1[0]+weight*dir[0],pt1[1]+weight*dir[1]));
		}
		delete []len_shift;
		return true;
	}

	bool _deform()
	{
		if(!has_split_point)
			return false;

		int contour_size = final_contour_x.size();
		if(contour_size < 2)
			return false;

		int fixed_idx0 = splitIdx[0];
		int fixed_idx1 = splitIdx[1];
		float pts1[2] = {splitX[0],splitY[0]};
		float pts2[2] = {splitX[1],splitY[1]};
		if(fixed_idx0 > fixed_idx1)
		{
			int tmp = fixed_idx0;
			fixed_idx0 = fixed_idx1;
			fixed_idx1 = tmp;
			float tmpcoord[2];
			memcpy(tmpcoord,pts1,sizeof(float)*2);
			memcpy(pts1,pts2,sizeof(float)*2);
			memcpy(pts2,tmpcoord,sizeof(float)*2);
		}

		std::vector<int> left_contour_x;
		std::vector<int> left_contour_y;
		std::vector<int> right_contour_x;
		std::vector<int> right_contour_y;

		

		left_contour_x.insert(left_contour_x.begin(),final_contour_x.begin()+fixed_idx0,final_contour_x.begin()+fixed_idx1+1);
		left_contour_y.insert(left_contour_y.begin(),final_contour_y.begin()+fixed_idx0,final_contour_y.begin()+fixed_idx1+1);

		right_contour_x.insert(right_contour_x.begin(),final_contour_x.begin()+fixed_idx1,final_contour_x.begin()+contour_size-1);
		right_contour_x.insert(right_contour_x.end(),final_contour_x.begin(),final_contour_x.begin()+fixed_idx0+1);

		right_contour_y.insert(right_contour_y.begin(),final_contour_y.begin()+fixed_idx1,final_contour_y.begin()+contour_size-1);
		right_contour_y.insert(right_contour_y.end(),final_contour_y.begin(),final_contour_y.begin()+fixed_idx0+1);

		std::vector<CtrlPt> left_ctrl_pts;
		std::vector<CtrlPt> right_ctrl_pts;
		if(!_map_contour(left_contour_x,left_contour_y,pts1,pts2,left_ctrl_pts)
			|| !_map_contour(right_contour_x,right_contour_y,pts2,pts1,right_ctrl_pts))
		{
			return false;
		}

		bool* nouseful_flag = new bool[width*height];
		bool* fixed_flag = new bool[width*height];
		float* init_coord = new float[width*height*2];
		memset(nouseful_flag,0,sizeof(bool)*width*height);
		memset(fixed_flag,0,sizeof(bool)*width*height);
		memset(init_coord,0,sizeof(float)*width*height*2);

		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				init_coord[(i*width+j)*2+0] = j;
				init_coord[(i*width+j)*2+1] = i;
			}
		}

		for(int i = 0;i < width;i++)
		{
			fixed_flag[0*width+i] = true;
			init_coord[(0*width+i)*2+0] = i;
			init_coord[(0*width+i)*2+1] = 0;
			fixed_flag[(height-1)*width+i] = true;
			init_coord[((height-1)*width+i)*2+0] = i;
			init_coord[((height-1)*width+i)*2+1] = height-1;
		}
		for(int i = 0;i < height;i++)
		{
			fixed_flag[i*width+0] = true;
			init_coord[(i*width+0)*2+0] = 0;
			init_coord[(i*width+0)*2+1] = i;
			fixed_flag[i*width+width-1] = true;
			init_coord[(i*width+width-1)*2+0] = width-1;
			init_coord[(i*width+width-1)*2+1] = i;
		}

		for(int i = 0;i < left_ctrl_pts.size();i++)
		{
			int ix = left_ctrl_pts[i].ix;
			int iy = left_ctrl_pts[i].iy;
			float fx = left_ctrl_pts[i].fx;
			float fy = left_ctrl_pts[i].fy;
			fixed_flag[iy*width+ix] = true;
			init_coord[(iy*width+ix)*2+0] = fx;
			init_coord[(iy*width+ix)*2+1] = fy;
		}
		for(int i = 0;i < right_ctrl_pts.size();i++)
		{
			int ix = right_ctrl_pts[i].ix;
			int iy = right_ctrl_pts[i].iy;
			float fx = right_ctrl_pts[i].fx;
			float fy = right_ctrl_pts[i].fy;
			fixed_flag[iy*width+ix] = true;
			init_coord[(iy*width+ix)*2+0] = fx;
			init_coord[(iy*width+ix)*2+1] = fy;
		}
		for(int i = 0;i < width*height;i++)
		{
			nouseful_flag[i] = objMask[i];
		}


		ZQ::ZQ_GridDeformation<float> deform;
		ZQ::ZQ_GridDeformationOptions opt;
		if(distance_weight == 0)
		{
			opt.methodType == ZQ::ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_ENERGY;
		}
		else
		{
			opt.methodType = ZQ::ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_DISTANCE_ENERGY;
			opt.distance = 1;
		}
			
		opt.line_weight = line_weight;
		opt.angle_weight = angle_weight;
		opt.distance_weight = distance_weight;
		printf("build matrix begin...\n");
		if(!deform.BuildMatrix(width,height,nouseful_flag,fixed_flag,opt))
		{
			delete []nouseful_flag;
			delete []fixed_flag;
			delete []init_coord;
			printf("failed to build matrix\n");
			return false;
		}
		printf("build matrix done\n");

		printf("solve begin...\n");
		if(!deform.Deformation(init_coord,outCoord))
		{
			delete []nouseful_flag;
			delete []fixed_flag;
			delete []init_coord;
			printf("failed to deform\n");
			return false;
		}
		printf("solve done\n");
		delete []nouseful_flag;
		delete []fixed_flag;
		delete []init_coord;
		return true;
	}

	void _drawToImage(IplImage* img)
	{
		for(int i = 0;i < nSelected;i++)
		{
			cvCircle(img,cvPoint(selectedX[i],selectedY[i]),size_point,color_selected,thick_point);
		}
		if(has_split_point)
		{
			for(int i = 0;i < nSelected;i++)
			{
				//cvCircle(img,cvPoint(splitX[i]*show_scale+show_shift,splitY[i]*show_scale+show_shift),size_point,color_split,thick_point);
				float cur_x = final_contour_x[splitIdx[i]];
				float cur_y = final_contour_y[splitIdx[i]];
				cvCircle(img,cvPoint(cur_x*show_scale+show_shift,cur_y*show_scale+show_shift),size_point,color_split,thick_point);
			}
		}
	}

	void _calSplitPoints()
	{
		
		if(final_contour_x.size() == 0)
			has_split_point = false;
		else
		{
			for(int i = 0;i < nSelected;i++)
			{
				float cur_x = (float)(selectedX[i] - show_shift)/show_scale;
				float cur_y = (float)(selectedY[i] - show_shift)/show_scale;
				splitIdx[i] = 0;
				double min_dist2 = (cur_x-final_contour_x[0])*(cur_x-final_contour_x[0])+(cur_y-final_contour_y[0])*(cur_y-final_contour_y[0]);
				for(int j = 1;j < final_contour_x.size();j++)
				{
					double cur_dist2 = (cur_x-final_contour_x[j])*(cur_x-final_contour_x[j])+(cur_y-final_contour_y[j])*(cur_y-final_contour_y[j]);
					if(cur_dist2 < min_dist2)
					{
						splitIdx[i] = j;
						min_dist2 = cur_dist2;
					}
				}
				//splitX[i] = final_contour_x[splitIdx[i]];
				//splitY[i] = final_contour_y[splitIdx[i]];
				splitX[i] = cur_x;
				splitY[i] = cur_y;
			}
			has_split_point = true;
		}
	}

	static void _mouseHandler(int event, int x,int y ,int flags, void* para)
	{
		InsertObjectToGrid* mApp = (InsertObjectToGrid*)para;

		if(event == CV_EVENT_LBUTTONDOWN)
		{
			if(mApp->nSelected < NMAX_SELECTED)
			{
				mApp->selectedX[mApp->nSelected] = x;
				mApp->selectedY[mApp->nSelected] = y;
				mApp->nSelected++;
				
				mApp->_calSplitPoints();
			}		
		}
		else if(event == CV_EVENT_RBUTTONDOWN)
		{
			if(mApp->nSelected > 0)
			{
				mApp->nSelected --;
				mApp->_calSplitPoints();
			}
		}
	}

	static void _objMask2ContourMask(const int width, const int height, const bool* objMask, bool* contourMask)
	{
		memset(contourMask,0,sizeof(bool)*width*height);
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(objMask[i*width+j])
				{
					for(int hh = __max(0,i-1);hh <= __min(height-1,i+1);hh++)
					{
						for(int ww = __max(0,j-1);ww <= __min(width-1,j+1);ww++)
						{
							contourMask[hh*width+ww] = true;
						}
					}
				}
			}
		}
	}


	static bool _contourMask2Contour(const int width, const int height, const bool* contourMask, std::vector<int>& contour_x, std::vector<int>& contour_y)
	{
		bool find_the_first = false;
		int first_x = -1;
		int first_y = -1;
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(contourMask[i*width+j])
				{
					find_the_first = true;
					first_x = j;
					first_y = i;
					break;
				}
			}
			if(find_the_first)
				break;
		}
		if(!find_the_first)
			return false;

		contour_x.clear();
		contour_y.clear();

		contour_x.push_back(first_x);
		contour_y.push_back(first_y);

		if(first_y+1 >= height) 
			return false;
		contour_x.push_back(first_x);
		contour_y.push_back(first_y+1);
		const int DIR_LEFT = 0;
		const int DIR_RIGHT = 1;
		const int DIR_UP = 2;
		const int DIR_DOWN = 3;
		int cur_dir = DIR_UP;
		int cur_x = first_x;
		int cur_y = first_y+1;
		while(true)
		{
			switch(cur_dir)
			{
			case DIR_UP:
				{
					if(cur_x-1 >= 0 && contourMask[cur_y*width+cur_x-1])
					{
						contour_x.push_back(cur_x-1);
						contour_y.push_back(cur_y);
						cur_dir = DIR_LEFT;
						cur_x = cur_x-1;
					}
					else if(cur_y+1 < height && contourMask[(cur_y+1)*width+cur_x])
					{
						contour_x.push_back(cur_x);
						contour_y.push_back(cur_y+1);
						cur_dir = DIR_UP;
						cur_y = cur_y+1;
					}
					else if(cur_x+1 < width && contourMask[cur_y*width+cur_x+1])
					{
						contour_x.push_back(cur_x+1);
						contour_y.push_back(cur_y);
						cur_dir = DIR_RIGHT;
						cur_x = cur_x+1;
					}
					else
					{
						return false;
					}
				}
				break;
			case DIR_DOWN:
				{
					if(cur_x+1 < width && contourMask[cur_y*width+cur_x+1])
					{
						contour_x.push_back(cur_x+1);
						contour_y.push_back(cur_y);
						cur_dir = DIR_RIGHT;
						cur_x = cur_x+1;
					}
					else if(cur_y-1 >= 0 && contourMask[(cur_y-1)*width+cur_x])
					{
						contour_x.push_back(cur_x);
						contour_y.push_back(cur_y-1);
						cur_dir = DIR_DOWN;
						cur_y = cur_y-1;
					}
					else if(cur_x-1 >= 0 && contourMask[cur_y*width+cur_x-1])
					{
						contour_x.push_back(cur_x-1);
						contour_y.push_back(cur_y);
						cur_dir = DIR_LEFT;
						cur_x = cur_x-1;
					}
					else
					{
						return false;
					}
				}
				break;
			case DIR_LEFT:
				{
					if(cur_y-1 >= 0 && contourMask[(cur_y-1)*width+cur_x])
					{
						contour_x.push_back(cur_x);
						contour_y.push_back(cur_y-1);
						cur_dir = DIR_DOWN;
						cur_y = cur_y-1;
					}
					else if(cur_x-1 >= 0 && contourMask[cur_y*width+cur_x-1])
					{
						contour_x.push_back(cur_x-1);
						contour_y.push_back(cur_y);
						cur_dir = DIR_LEFT;
						cur_x = cur_x-1;
					}
					else if(cur_y+1 < height && contourMask[(cur_y+1)*width+cur_x])
					{
						contour_x.push_back(cur_x);
						contour_y.push_back(cur_y+1);
						cur_dir = DIR_UP;
						cur_y = cur_y+1;
					}
					else
					{
						return false;
					}
				}
				break;;
			case DIR_RIGHT:
				{
					if(cur_y+1 < height && contourMask[(cur_y+1)*width+cur_x])
					{
						contour_x.push_back(cur_x);
						contour_y.push_back(cur_y+1);
						cur_dir = DIR_UP;
						cur_y = cur_y+1;
					}
					else if(cur_x+1 < width && contourMask[cur_y*width+(cur_x+1)])
					{
						contour_x.push_back(cur_x+1);
						contour_y.push_back(cur_y);
						cur_dir = DIR_RIGHT;
						cur_x = cur_x+1;
					}
					else if(cur_y-1 >= 0 && contourMask[(cur_y-1)*width+cur_x])
					{
						contour_x.push_back(cur_x);
						contour_y.push_back(cur_y-1);
						cur_dir = DIR_DOWN;
						cur_y = cur_y-1;
					}
					else
					{
						return false;
					}
				}
				break;
			}
			if(cur_x == first_x && cur_y == first_y)
				break;
		}
		return true;
	}


	static bool _isRealContour(const int width, const int height, const bool* objMask, const int x, const int y)
	{
		if(x < 0 || x >= width || y < 0 || y >= height)
			return false;

		if(x-1 >= 0 && objMask[y*width+x-1])
			return true;
		if(x+1 < width && objMask[y*width+x+1])
			return true;
		if(y-1 >= 0 && objMask[(y-1)*width+x])
			return true;
		if(y+1 < height && objMask[(y+1)*width+x])
			return true;

		return false;
	}


	static bool _reduceContour(const int width, const int height, const bool* objMask, const std::vector<int>& contour_x, const std::vector<int>& contour_y, std::vector<int>& out_x, std::vector<int>& out_y)
	{
		int pts_num = contour_x.size();
		if(pts_num <= 2 || pts_num != contour_y.size())
			return false;

		std::vector<bool> isRealContour(pts_num);
		for(int i = 0;i < pts_num;i++)
		{
			isRealContour[i] = _isRealContour(width,height,objMask,contour_x[i],contour_y[i]);
		}

		int start_i = -1;
		for(int i = 0;i < pts_num;i++)
		{
			if(isRealContour[i])
			{
				start_i = i;
				break;
			}
		}
		if(start_i < 0)
			return false;

		int realPtsNum = pts_num-1;
		out_x.clear();
		out_y.clear();
		out_x.push_back(contour_x[start_i]);
		out_y.push_back(contour_y[start_i]);
		for(int i = start_i+1;i < realPtsNum;i++)
		{
			if(!isRealContour[i] && abs(contour_x[i+1]-contour_x[i-1]) == 1 && abs(contour_y[i+1]-contour_y[i-1]) == 1)
			{
			}
			else
			{
				out_x.push_back(contour_x[i]);
				out_y.push_back(contour_y[i]);
			}
		}
		for(int i = 0;i <= start_i-1;i++)
		{
			if(i == 0)
			{
				if(!isRealContour[0] && abs(contour_x[1]-contour_x[realPtsNum-1]) == 1 && abs(contour_y[1]-contour_y[realPtsNum-1]) == 1)
				{

				}
				else
				{
					out_x.push_back(contour_x[i]);
					out_y.push_back(contour_y[i]);
				}
			}
			else
			{
				if(!isRealContour[i] && abs(contour_x[i+1]-contour_x[i-1]) == 1 && abs(contour_y[i+1]-contour_y[i-1]) == 1)
				{
				}
				else
				{
					out_x.push_back(contour_x[i]);
					out_y.push_back(contour_y[i]);
				}
			}
		}

		out_x.push_back(contour_x[start_i]);
		out_y.push_back(contour_y[start_i]);
		return true;
	}
};



#endif