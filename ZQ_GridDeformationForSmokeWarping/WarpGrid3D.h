#ifndef _WARP_GRID_3D_H_
#define _WARP_GRID_3D_H_
#pragma once

#include "ZQ_DoubleImage3D.h"
#include "ZQ_GridDeformation3D.h"
#include "RenderGridWithTexture3D.h"
#include <time.h>

class WarpGrid3DOptions
{
	class NouseSegment3D
	{
	public:
		NouseSegment3D(){idx1 = idx2 = start = end = 0;}
		NouseSegment3D(int _idx1, int _idx2, int _st, int _end)
			: idx1(_idx1),idx2(_idx2),start(_st),end(_end){}
		int idx1, idx2, start, end;
	};
	class FixedSegment3D
	{
	public:
		FixedSegment3D(){idx1 = idx2 = start = end = 0; x1 = y1 = z1 = x2 = y2 = z2 = 0;}
		FixedSegment3D(int _idx1, int _idx2, int _st, int _end, float _x1, float _y1, float _z1, float _x2, float _y2, float _z2)
			:idx1(_idx1),idx2(_idx2),start(_st),end(_end),x1(_x1),y1(_y1),z1(_z1),x2(_x2),y2(_y2),z2(_z2){}

		int idx1,idx2;
		int start,end;
		float x1,y1,z1,x2,y2,z2;
	};
	class FixedLine3D
	{
	public:
		FixedLine3D(){idx1 = idx2 = 0; x1 = y1 = z1 = x2 = y2 = z2 = 0;}
		FixedLine3D(int _idx1, int _idx2, int _x1, int _y1, int _z1, int _x2, int _y2, int _z2)
			:idx1(_idx1),idx2(_idx2),x1(_x1),y1(_y1),z1(_z1),x2(_x2),y2(_y2),z2(_z2){}
		int idx1,idx2;
		float x1,y1,z1,x2,y2,z2;
	};
	class FixedPoint3D
	{
	public:
		FixedPoint3D(){coord_x = coord_y = coord_z = 0; pos_x = pos_y = pos_z = 0;}
		FixedPoint3D(int _cx, int _cy, int _cz, float _px, float _py, float _pz)
			:coord_x(_cx),coord_y(_cy),coord_z(_cz),pos_x(_px),pos_y(_py),pos_z(_pz){}

		int coord_x,coord_y,coord_z;
		float pos_x,pos_y,pos_z;
	};
public:
	WarpGrid3DOptions(){Reset();}
	~WarpGrid3DOptions(){}

	enum CONST_VAL{
		FILE_LEN = 200
	};
	bool has_show_file;
	char show_file[FILE_LEN];
	bool has_map_file;
	char map_file[FILE_LEN];
	bool has_map_mask_file;
	char map_mask_file[FILE_LEN];
	int GridWidth,GridHeight,GridDepth;
	int BoundingBoxmin[3],BoundingBoxmax[3];
	ZQ::ZQ_DImage3D<bool> nouseful_flag;
	ZQ::ZQ_DImage3D<bool> fixed_flag;
	ZQ::ZQ_DImage3D<float> init_coord;
	float line_weight;
	float angle_weight;
	float distance_weight;
	float distance;
	bool xloop;
	float max_iteration;
	float FPiteration;
	int show_scale;

	void Reset()
	{
		has_show_file = false;
		show_file[0] = '\0';
		has_map_file = false;
		map_file[0] = '\0';
		has_map_mask_file = false;
		map_mask_file[0] = '\0';
		GridWidth = GridHeight = 0;
		BoundingBoxmin[0] = BoundingBoxmin[1] = BoundingBoxmin[2] = 0;
		BoundingBoxmax[0] = BoundingBoxmax[1] = BoundingBoxmax[2] = 0;
		line_weight = 1;
		angle_weight = 0.1;
		distance_weight = 0;
		distance = 1;
		xloop = false;
		max_iteration = 1000;
		FPiteration = 5;
		show_scale = 10;

	}

	bool HandleArgs(const int argc, const char** argv)
	{
		std::vector<NouseSegment3D> nouse_seg_x;
		std::vector<NouseSegment3D> nouse_seg_y;
		std::vector<NouseSegment3D> nouse_seg_z;
		std::vector<FixedSegment3D> fixed_seg_x;
		std::vector<FixedSegment3D> fixed_seg_y;
		std::vector<FixedSegment3D> fixed_seg_z;
		std::vector<FixedLine3D> fixed_line_x;
		std::vector<FixedLine3D> fixed_line_y;
		std::vector<FixedLine3D> fixed_line_z;
		std::vector<FixedPoint3D> fixed_pts;
		for(int k = 0;k < argc;k++)
		{
			if(_strcmpi(argv[k],"show_file") == 0)
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
			else if(_strcmpi(argv[k],"map_mask_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of map_mask_file ?\n");
					return false;
				}
				strcpy(map_mask_file,argv[k]);
				has_map_mask_file = true;
			}
			else if(_strcmpi(argv[k],"show_scale") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of show_scale ?\n");
					return false;
				}
				show_scale = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"max_iteration") == 0)
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
			else if(_strcmpi(argv[k],"distance") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of distance ?\n");
					return false;
				}
				distance = atof(argv[k]);
			}
			else if(_strcmpi(argv[k],"xloop") == 0)
			{
				xloop = true;
			}
			else if(_strcmpi(argv[k],"resolution") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of resolution:x ?\n");
					return false;
				}
				GridWidth = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of resolution:y ?\n");
					return false;
				}
				GridHeight = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of resolution:z ?\n");
					return false;
				}
				GridDepth = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"boundingboxmin") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of boundingboxmin:x ?\n");
					return false;
				}
				BoundingBoxmin[0] = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of boundingboxmin:y ?\n");
					return false;
				}
				BoundingBoxmin[1] = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of boundingboxmin:z ?\n");
					return false;
				}
				BoundingBoxmin[2] = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"boundingboxmax") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of boundingboxmax:x ?\n");
					return false;
				}
				BoundingBoxmax[0] = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of boundingboxmax:y ?\n");
					return false;
				}
				BoundingBoxmax[1] = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of boundingboxmax:z ?\n");
					return false;
				}
				BoundingBoxmax[2] = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"nouseSegX") == 0)
			{
				int idx1,idx2,st,end;
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegX:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegX:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegX:start ?\n");
					return false;
				}
				st = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegX:end ?\n");
					return false;
				}
				end = atoi(argv[k]);
				nouse_seg_x.push_back(NouseSegment3D(idx1,idx2,st,end));
			}
			else if(_strcmpi(argv[k],"nouseSegY") == 0)
			{
				int idx1,idx2,st,end;
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegY:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegY:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegY:start ?\n");
					return false;
				}
				st = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegY:end ?\n");
					return false;
				}
				end = atoi(argv[k]);
				nouse_seg_y.push_back(NouseSegment3D(idx1,idx2,st,end));
			}
			else if(_strcmpi(argv[k],"nouseSegZ") == 0)
			{
				int idx1,idx2,st,end;
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegZ:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegZ:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegZ:start ?\n");
					return false;
				}
				st = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegZ:end ?\n");
					return false;
				}
				end = atoi(argv[k]);
				nouse_seg_z.push_back(NouseSegment3D(idx1,idx2,st,end));
			}
			else if(_strcmpi(argv[k],"fixedSegX") == 0)
			{
				int idx1,idx2,st,end;
				float x1,y1,z1,x2,y2,z2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:start ?\n");
					return false;
				}
				st = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:end ?\n");
					return false;
				}
				end = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:x1 ?\n");
					return false;
				}
				x1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:y1 ?\n");
					return false;
				}
				y1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:z1 ?\n");
					return false;
				}
				z1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:x2 ?\n");
					return false;
				}
				x2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:y2 ?\n");
					return false;
				}
				y2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:z2 ?\n");
					return false;
				}
				z2 = atof(argv[k]);
				fixed_seg_x.push_back(FixedSegment3D(idx1,idx2,st,end,x1,y1,z1,x2,y2,z2));
			}
			else if(_strcmpi(argv[k],"fixedSegY") == 0)
			{
				int idx1,idx2,st,end;
				float x1,y1,z1,x2,y2,z2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:start ?\n");
					return false;
				}
				st = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:end ?\n");
					return false;
				}
				end = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:x1 ?\n");
					return false;
				}
				x1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:y1 ?\n");
					return false;
				}
				y1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:z1 ?\n");
					return false;
				}
				z1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:x2 ?\n");
					return false;
				}
				x2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:y2 ?\n");
					return false;
				}
				y2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:z2 ?\n");
					return false;
				}
				z2 = atof(argv[k]);
				fixed_seg_y.push_back(FixedSegment3D(idx1,idx2,st,end,x1,y1,z1,x2,y2,z2));
			}
			else if(_strcmpi(argv[k],"fixedSegZ") == 0)
			{
				int idx1,idx2,st,end;
				float x1,y1,z1,x2,y2,z2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:start ?\n");
					return false;
				}
				st = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:end ?\n");
					return false;
				}
				end = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:x1 ?\n");
					return false;
				}
				x1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:y1 ?\n");
					return false;
				}
				y1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:z1 ?\n");
					return false;
				}
				z1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:x2 ?\n");
					return false;
				}
				x2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:y2 ?\n");
					return false;
				}
				y2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegZ:z2 ?\n");
					return false;
				}
				z2 = atof(argv[k]);
				fixed_seg_z.push_back(FixedSegment3D(idx1,idx2,st,end,x1,y1,z1,x2,y2,z2));
			}
			else if(_strcmpi(argv[k],"fixedLineX") == 0)
			{
				int idx1,idx2;
				float x1,y1,z1,x2,y2,z2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:x1 ?\n");
					return false;
				}
				x1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:y1 ?\n");
					return false;
				}
				y1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:z1 ?\n");
					return false;
				}
				z1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:x2 ?\n");
					return false;
				}
				x2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:y2 ?\n");
					return false;
				}
				y2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:z2 ?\n");
					return false;
				}
				z2 = atof(argv[k]);
				fixed_line_x.push_back(FixedLine3D(idx1,idx2,x1,y1,z1,x2,y2,z2));
			}
			else if(_strcmpi(argv[k],"fixedLineY") == 0)
			{
				int idx1,idx2;
				float x1,y1,z1,x2,y2,z2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:x1 ?\n");
					return false;
				}
				x1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:y1 ?\n");
					return false;
				}
				y1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:z1 ?\n");
					return false;
				}
				z1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:x2 ?\n");
					return false;
				}
				x2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:y2 ?\n");
					return false;
				}
				y2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:z2 ?\n");
					return false;
				}
				z2 = atof(argv[k]);
				fixed_line_y.push_back(FixedLine3D(idx1,idx2,x1,y1,z1,x2,y2,z2));
			}
			else if(_strcmpi(argv[k],"fixedLineZ") == 0)
			{
				int idx1,idx2;
				float x1,y1,z1,x2,y2,z2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:idx1 ?\n");
					return false;
				}
				idx1 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:idx2 ?\n");
					return false;
				}
				idx2 = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:x1 ?\n");
					return false;
				}
				x1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:y1 ?\n");
					return false;
				}
				y1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:z1 ?\n");
					return false;
				}
				z1 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:x2 ?\n");
					return false;
				}
				x2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:y2 ?\n");
					return false;
				}
				y2 = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineZ:z2 ?\n");
					return false;
				}
				z2 = atof(argv[k]);
				fixed_line_z.push_back(FixedLine3D(idx1,idx2,x1,y1,z1,x2,y2,z2));
			}
			else if(_strcmpi(argv[k],"fixedPoint") == 0)
			{
				int cx,cy,cz;
				float px,py,pz;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedPoint: coord_x ?\n");
					return false;
				}
				cx = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedPoint: coord_y ?\n");
					return false;
				}
				cy = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedPoint: coord_z ?\n");
					return false;
				}
				cz = atoi(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedPoint: pos_x ?\n");
					return false;
				}
				px = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedPoint: pos_y ?\n");
					return false;
				}
				py = atof(argv[k]);
				k++;
				if(k >= argc)
				{
					printf("the value of fixedPoint: pos_z ?\n");
					return false;
				}
				pz = atof(argv[k]);
				fixed_pts.push_back(FixedPoint3D(cx,cy,cz,px,py,pz));
			}
			else
			{
				printf("unknown parameter %s\n",argv[k]);
				return false;
			}
		}

		if(GridWidth < 2 || GridHeight < 2 || GridHeight < 2)
		{
			printf("invalid grid resolution\n");
			return false;
		}
		nouseful_flag.allocate(GridWidth,GridHeight,GridDepth);
		fixed_flag.allocate(GridWidth,GridHeight,GridDepth);
		init_coord.allocate(GridWidth,GridHeight,GridDepth,3);
		bool*& nouseful_flag_data = nouseful_flag.data();
		bool*& fixed_flag_data = fixed_flag.data();
		float*& init_coord_data = init_coord.data();
		for(int i = 0;i < nouse_seg_x.size();i++)
		{
			int idx_y = nouse_seg_x[i].idx1;
			int idx_z = nouse_seg_x[i].idx2;
			if(idx_y < 0 || idx_y >= GridHeight || idx_z < 0 || idx_z >= GridDepth)
			{
				printf("invalid nouseSegX:idx=(:,%d,%d)\n",idx_y,idx_z);
				return false;
			}
			int start = nouse_seg_x[i].start;
			int end = nouse_seg_x[i].end;
			if(start < 0 || start >= GridWidth || end < 0 || end >= GridWidth || start >= end)
			{
				printf("invalid nouseSegX:start=%d,end=%d\n",start,end);
				return false;
			}
			for(int w = start;w <= end;w++)
				nouseful_flag_data[idx_z*GridHeight*GridWidth+idx_y*GridWidth+w] = true;
		}
		for(int i = 0;i < nouse_seg_y.size();i++)
		{
			int idx_x = nouse_seg_y[i].idx1;
			int idx_z = nouse_seg_y[i].idx2;
			if(idx_x < 0 || idx_x >= GridWidth || idx_z < 0 || idx_z >= GridDepth)
			{
				printf("invalid nouseSegY:idx=(%d,:,%d)\n",idx_x,idx_z);
				return false;
			}
			int start = nouse_seg_y[i].start;
			int end = nouse_seg_y[i].end;
			if(start < 0 || start >= GridHeight || end < 0 || end >= GridHeight || start >= end)
			{
				printf("invalid nouseSegY:start=%d,end=%d\n",start,end);
				return false;
			}
			for(int h = start;h <= end;h++)
				nouseful_flag_data[idx_z*GridHeight*GridWidth+h*GridWidth+idx_x] = true;
		}
		for(int i = 0;i < nouse_seg_z.size();i++)
		{
			int idx_x = nouse_seg_z[i].idx1;
			int idx_y = nouse_seg_z[i].idx2;
			if(idx_x < 0 || idx_x >= GridWidth || idx_y < 0 || idx_y >= GridHeight)
			{
				printf("invalid nouseSegZ:idx=(%d,%d,:)\n",idx_x,idx_y);
				return false;
			}
			int start = nouse_seg_y[i].start;
			int end = nouse_seg_y[i].end;
			if(start < 0 || start >= GridDepth || end < 0 || end >= GridDepth || start >= end)
			{
				printf("invalid nouseSegZ:start=%d,end=%d\n",start,end);
				return false;
			}
			for(int d = start;d <= end;d++)
				nouseful_flag_data[d*GridHeight*GridWidth+idx_y*GridWidth+idx_x] = true;
		}
		for(int i = 0;i < fixed_seg_x.size();i++)
		{
			int idx_y = fixed_seg_x[i].idx1;
			int idx_z = fixed_seg_x[i].idx2;
			if(idx_y < 0 || idx_y >= GridHeight || idx_z < 0 || idx_z >= GridDepth)
			{
				printf("invalid fixedSegX:idx=(:,%d,%d)\n",idx_y,idx_z);
				return false;
			}
			int start = fixed_seg_x[i].start;
			int end = fixed_seg_x[i].end;
			if(start < 0 || start >= GridWidth || end < 0 || end >= GridWidth || start >= end)
			{
				printf("invalid fixedSegX:start=%d,end=%d\n",start,end);
				return false;
			}
			float x1 = fixed_seg_x[i].x1;
			float y1 = fixed_seg_x[i].y1;
			float z1 = fixed_seg_x[i].z1;
			float x2 = fixed_seg_x[i].x2;
			float y2 = fixed_seg_x[i].y2;
			float z2 = fixed_seg_x[i].z2;
			for(int w = start;w <= end;w++)
			{
				fixed_flag_data[idx_z*GridHeight*GridWidth+idx_y*GridWidth+w] = true;
				float weight = (float)(w-start)/(end-start);
				init_coord_data[(idx_z*GridHeight*GridWidth+idx_y*GridWidth+w)*3+0] = x1 + weight*(x2-x1);
				init_coord_data[(idx_z*GridHeight*GridWidth+idx_y*GridWidth+w)*3+1] = y1 + weight*(y2-y1);
				init_coord_data[(idx_z*GridHeight*GridWidth+idx_y*GridWidth+w)*3+2] = z1 + weight*(z2-z1);
			}
		}
		for(int i = 0;i < fixed_seg_y.size();i++)
		{
			int idx_x = fixed_seg_y[i].idx1;
			int idx_z = fixed_seg_y[i].idx2;
			if(idx_x < 0 || idx_x >= GridWidth || idx_z < 0 || idx_z >= GridDepth)
			{
				printf("invalid fixedSegY:idx=(%d,:,%d)\n",idx_x,idx_z);
				return false;
			}
			int start = fixed_seg_y[i].start;
			int end = fixed_seg_y[i].end;
			if(start < 0 || start >= GridHeight || end < 0 || end >= GridHeight || start >= end)
			{
				printf("invalid fixedSegY:start=%d,end=%d\n",start,end);
				return false;
			}
			float x1 = fixed_seg_y[i].x1;
			float y1 = fixed_seg_y[i].y1;
			float z1 = fixed_seg_y[i].z1;
			float x2 = fixed_seg_y[i].x2;
			float y2 = fixed_seg_y[i].y2;
			float z2 = fixed_seg_y[i].z2;
			for(int h = start;h <= end;h++)
			{
				fixed_flag_data[idx_z*GridHeight*GridWidth+h*GridWidth+idx_x] = true;
				float weight = (float)(h-start)/(end-start);
				init_coord_data[(idx_z*GridHeight*GridWidth+h*GridWidth+idx_x)*3+0] = x1 + weight*(x2-x1);
				init_coord_data[(idx_z*GridHeight*GridWidth+h*GridWidth+idx_x)*3+1] = y1 + weight*(y2-y1);
				init_coord_data[(idx_z*GridHeight*GridWidth+h*GridWidth+idx_x)*3+2] = z1 + weight*(z2-z1);
			}
		}
		for(int i = 0;i < fixed_seg_z.size();i++)
		{
			int idx_x = fixed_seg_z[i].idx1;
			int idx_y = fixed_seg_z[i].idx2;
			if(idx_x < 0 || idx_x >= GridWidth || idx_y < 0 || idx_y >= GridHeight)
			{
				printf("invalid fixedSegZ:idx=(%d,%d,:)\n",idx_x,idx_y);
				return false;
			}
			int start = fixed_seg_z[i].start;
			int end = fixed_seg_z[i].end;
			if(start < 0 || start >= GridDepth || end < 0 || end >= GridDepth || start >= end)
			{
				printf("invalid fixedSegZ:start=%d,end=%d\n",start,end);
				return false;
			}
			float x1 = fixed_seg_z[i].x1;
			float y1 = fixed_seg_z[i].y1;
			float z1 = fixed_seg_z[i].z1;
			float x2 = fixed_seg_z[i].x2;
			float y2 = fixed_seg_z[i].y2;
			float z2 = fixed_seg_z[i].z2;
			for(int d = start;d <= end;d++)
			{
				fixed_flag_data[d*GridHeight*GridWidth+idx_y*GridWidth+idx_x] = true;
				float weight = (float)(d-start)/(end-start);
				init_coord_data[(d*GridHeight*GridWidth+idx_y*GridWidth+idx_x)*3+0] = x1 + weight*(x2-x1);
				init_coord_data[(d*GridHeight*GridWidth+idx_y*GridWidth+idx_x)*3+1] = y1 + weight*(y2-y1);
				init_coord_data[(d*GridHeight*GridWidth+idx_y*GridWidth+idx_x)*3+2] = z1 + weight*(z2-z1);
			}
		}
		for(int i = 0;i < fixed_line_x.size();i++)
		{
			int idx_y = fixed_line_x[i].idx1;
			int idx_z = fixed_line_x[i].idx2;
			if(idx_y < 0 || idx_y >= GridHeight || idx_z < 0 || idx_z >= GridDepth)
			{
				printf("invalid fixedLineX:idx=(:,%d,%d)\n",idx_y,idx_z);
				return false;
			}
			float x1 = fixed_line_x[i].x1;
			float y1 = fixed_line_x[i].y1;
			float z1 = fixed_line_x[i].z1;
			float x2 = fixed_line_x[i].x2;
			float y2 = fixed_line_x[i].y2;
			float z2 = fixed_line_x[i].z2;
			for(int w = 0;w < GridWidth;w++)
			{
				fixed_flag_data[idx_z*GridHeight*GridWidth+idx_y*GridWidth+w] = true;
				float weight = (float)w/(GridWidth-1);
				init_coord_data[(idx_z*GridHeight*GridWidth+idx_y*GridWidth+w)*3+0] = x1 + weight*(x2-x1);
				init_coord_data[(idx_z*GridHeight*GridWidth+idx_y*GridWidth+w)*3+1] = x1 + weight*(x2-x1);
				init_coord_data[(idx_z*GridHeight*GridWidth+idx_y*GridWidth+w)*3+2] = x1 + weight*(x2-x1);
			}
		}
		for(int i = 0;i < fixed_line_y.size();i++)
		{
			int idx_x = fixed_line_y[i].idx1;
			int idx_z = fixed_line_y[i].idx2;
			if(idx_x < 0 || idx_x >= GridWidth || idx_z < 0 || idx_z >= GridDepth)
			{
				printf("invalid fixedLineY:idx=(%d,:,%d)\n",idx_x,idx_z);
				return false;
			}
			float x1 = fixed_line_x[i].x1;
			float y1 = fixed_line_x[i].y1;
			float z1 = fixed_line_x[i].z1;
			float x2 = fixed_line_x[i].x2;
			float y2 = fixed_line_x[i].y2;
			float z2 = fixed_line_x[i].z2;
			for(int h = 0;h < GridHeight;h++)
			{
				fixed_flag_data[idx_z*GridHeight*GridWidth+h*GridWidth+idx_x] = true;
				float weight = (float)h/(GridHeight-1);
				init_coord_data[(idx_z*GridHeight*GridWidth+h*GridWidth+idx_x)*3+0] = x1 + weight*(x2-x1);
				init_coord_data[(idx_z*GridHeight*GridWidth+h*GridWidth+idx_x)*3+1] = x1 + weight*(x2-x1);
				init_coord_data[(idx_z*GridHeight*GridWidth+h*GridWidth+idx_x)*3+2] = x1 + weight*(x2-x1);
			}
		}
		for(int i = 0;i < fixed_line_z.size();i++)
		{
			int idx_x = fixed_line_z[i].idx1;
			int idx_y = fixed_line_z[i].idx2;
			if(idx_x < 0 || idx_x >= GridWidth || idx_y < 0 || idx_y >= GridHeight)
			{
				printf("invalid fixedLineZ:idx=(%d,%d,:)\n",idx_x,idx_y);
				return false;
			}
			float x1 = fixed_line_x[i].x1;
			float y1 = fixed_line_x[i].y1;
			float z1 = fixed_line_x[i].z1;
			float x2 = fixed_line_x[i].x2;
			float y2 = fixed_line_x[i].y2;
			float z2 = fixed_line_x[i].z2;
			for(int d = 0;d < GridDepth;d++)
			{
				fixed_flag_data[d*GridHeight*GridWidth+idx_y*GridWidth+idx_x] = true;
				float weight = (float)d/(GridDepth-1);
				init_coord_data[(d*GridHeight*GridWidth+idx_y*GridWidth+idx_x)*3+0] = x1 + weight*(x2-x1);
				init_coord_data[(d*GridHeight*GridWidth+idx_y*GridWidth+idx_x)*3+1] = x1 + weight*(x2-x1);
				init_coord_data[(d*GridHeight*GridWidth+idx_y*GridWidth+idx_x)*3+2] = x1 + weight*(x2-x1);
			}
		}
		for(int i = 0;i < fixed_pts.size();i++)
		{
			int cx = fixed_pts[i].coord_x;
			int cy = fixed_pts[i].coord_y;
			int cz = fixed_pts[i].coord_z;
			if(cx < 0 || cx >= GridWidth || cy < 0 || cy >= GridHeight || cz < 0 || cz >= GridDepth)
			{
				printf("invalid fixedPoint: <x,y,z> = <%d,%d,%d>\n",cx,cy,cz);
				return false;
			}
			float px = fixed_pts[i].pos_x;
			float py = fixed_pts[i].pos_y;
			float pz = fixed_pts[i].pos_z;
			fixed_flag_data[cz*GridHeight*GridWidth+cy*GridWidth+cx] = true;
			init_coord_data[(cz*GridHeight*GridWidth+cy*GridWidth+cx)*3+0] = px;
			init_coord_data[(cz*GridHeight*GridWidth+cy*GridWidth+cx)*3+1] = py;
			init_coord_data[(cz*GridHeight*GridWidth+cy*GridWidth+cx)*3+2] = pz;
		}
		return true;
	}
};

class WarpGrid3D
{
public:
	static bool go(const WarpGrid3DOptions& opt)
	{
		int width = opt.GridWidth;
		int height = opt.GridHeight;
		int depth = opt.GridDepth;
		if(!opt.nouseful_flag.matchDimension(width,height,depth,1)
			|| !opt.fixed_flag.matchDimension(width,height,depth,1)
			|| !opt.init_coord.matchDimension(width,height,depth,3))
		{
			printf("dimensions don't match\n");
			return false;
		}

		if(opt.has_map_mask_file)
		{
			ZQ::ZQ_DImage3D<float> obj_mask(width,height,depth);
			for(int i = 0;i < width*height*depth;i++)
				obj_mask.data()[i] = opt.nouseful_flag.data()[i] ? 0 : 1;
			if(!obj_mask.saveImage(opt.map_mask_file))
			{
				printf("failed to save %s\n",opt.map_mask_file);
				return false;
			}
		}
		ZQ::ZQ_DImage3D<float> out_coord(width,height,depth,3);
		ZQ::ZQ_GridDeformation3D<float> deform;
		ZQ::ZQ_GridDeformation3DOptions w_opt;
		
		if(opt.xloop)
			w_opt.methodType = ZQ::ZQ_GridDeformation3DOptions::METHOD_ARAP_VERT_AS_CENTER_XLOOP;
		else
			w_opt.methodType = ZQ::ZQ_GridDeformation3DOptions::METHOD_ARAP_VERT_AS_CENTER;
		
		w_opt.line_weight = opt.line_weight;
		w_opt.angle_weight = opt.angle_weight;
		w_opt.distance_weight = opt.distance_weight;
		w_opt.distance = opt.distance;
		w_opt.iteration = opt.max_iteration;
		w_opt.FPIteration = opt.FPiteration;
		clock_t build_st = clock();
		if(!deform.BuildMatrix(width,height,depth,opt.nouseful_flag.data(),opt.fixed_flag.data(),w_opt))
		{
			printf("failed to build matrix for deformation\n");
			return false;
		}
		clock_t build_end = clock();
		printf("build cost: %f seconds\n",0.001*(build_end-build_st));

		clock_t deform_st = clock();
		if(!deform.Deformation(opt.init_coord.data(),out_coord.data(),false))
		{
			printf("failed to deform\n");
			return false;
		}
		clock_t deform_end = clock();
		printf("build cost: %f seconds\n",0.001*(deform_end-deform_st));

		if(opt.has_map_file)
		{
			if(!out_coord.saveImage(opt.map_file))
			{
				printf("failed to save %s\n",opt.map_file);
				return false;
			}
		}

		return true;
	}
};

#endif