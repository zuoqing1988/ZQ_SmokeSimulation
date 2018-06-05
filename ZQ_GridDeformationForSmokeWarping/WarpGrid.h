#ifndef _WARP_GRID_H_
#define _WARP_GRID_H_
#pragma once

#include "ZQ_DoubleImage.h"
#include "ZQ_GridDeformation.h"
#include "RenderGridWithTexture.h"
#include <time.h>

class WarpGridOptions
{
	class NouseSegment
	{
	public:
		NouseSegment(){idx = start = end = 0;}
		NouseSegment(int _idx, int _st, int _end)
			: idx(_idx),start(_st),end(_end){}
		int idx, start, end;
	};
	class FixedSegment
	{
	public:
		FixedSegment(){idx = start = end = 0; x1 = y1 = x2 = y2 = 0;}
		FixedSegment(int _idx, int _st, int _end, float _x1, float _y1, float _x2, float _y2)
			:idx(_idx),start(_st),end(_end),x1(_x1),y1(_y1),x2(_x2),y2(_y2){}

		int idx;
		int start,end;
		float x1,y1,x2,y2;
	};
	class FixedLine
	{
	public:
		FixedLine(){idx = 0; x1 = y1 = x2 = y2 = 0;}
		FixedLine(int _idx, int _x1, int _y1, int _x2, int _y2)
			:idx(_idx),x1(_x1),y1(_y1),x2(_x2),y2(_y2){}
		int idx;
		float x1,y1,x2,y2;
	};
	class FixedPoint
	{
	public:
		FixedPoint(){coord_x = coord_y = 0; pos_x = pos_y = 0;}
		FixedPoint(int _cx, int _cy, float _px, float _py)
			:coord_x(_cx),coord_y(_cy),pos_x(_px),pos_y(_py){}

		int coord_x,coord_y;
		float pos_x,pos_y;
	};
public:
	WarpGridOptions(){Reset();}
	~WarpGridOptions(){}

	enum CONST_VAL{
		FILE_LEN = 200
	};
	bool has_show_file;
	char show_file[FILE_LEN];
	bool has_map_file;
	char map_file[FILE_LEN];
	bool has_map_mask_file;
	char map_mask_file[FILE_LEN];
	int GridWidth,GridHeight;
	int BoundingBoxmin[2],BoundingBoxmax[2];
	ZQ::ZQ_DImage<bool> nouseful_flag;
	ZQ::ZQ_DImage<bool> fixed_flag;
	ZQ::ZQ_DImage<float> init_coord;
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
		BoundingBoxmin[0] = BoundingBoxmin[1] = 0;
		BoundingBoxmax[0] = BoundingBoxmax[1] = 0;
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
		std::vector<NouseSegment> nouse_seg_x;
		std::vector<NouseSegment> nouse_seg_y;
		std::vector<FixedSegment> fixed_seg_x;
		std::vector<FixedSegment> fixed_seg_y;
		std::vector<FixedLine> fixed_line_x;
		std::vector<FixedLine> fixed_line_y;
		std::vector<FixedPoint> fixed_pts;
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
			}
			else if(_strcmpi(argv[k],"nouseSegX") == 0)
			{
				int idx,st,end;
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegX:idx ?\n");
					return false;
				}
				idx = atoi(argv[k]);
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
				nouse_seg_x.push_back(NouseSegment(idx,st,end));
			}
			else if(_strcmpi(argv[k],"nouseSegY") == 0)
			{
				int idx,st,end;
				k++;
				if(k >= argc)
				{
					printf("the value of nouseSegY:idx ?\n");
					return false;
				}
				idx = atoi(argv[k]);
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
				nouse_seg_y.push_back(NouseSegment(idx,st,end));
			}
			else if(_strcmpi(argv[k],"fixedSegX") == 0)
			{
				int idx,st,end;
				float x1,y1,x2,y2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegX:idx ?\n");
					return false;
				}
				idx = atoi(argv[k]);
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
				fixed_seg_x.push_back(FixedSegment(idx,st,end,x1,y1,x2,y2));
			}
			else if(_strcmpi(argv[k],"fixedSegY") == 0)
			{
				int idx,st,end;
				float x1,y1,x2,y2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedSegY:idx ?\n");
					return false;
				}
				idx = atoi(argv[k]);
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
				fixed_seg_y.push_back(FixedSegment(idx,st,end,x1,y1,x2,y2));
			}
			else if(_strcmpi(argv[k],"fixedLineX") == 0)
			{
				int idx;
				float x1,y1,x2,y2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineX:idx ?\n");
					return false;
				}
				idx = atoi(argv[k]);
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
				fixed_line_x.push_back(FixedLine(idx,x1,y1,x2,y2));
			}
			else if(_strcmpi(argv[k],"fixedLineY") == 0)
			{
				int idx;
				float x1,y1,x2,y2;
				k++;
				if(k >= argc)
				{
					printf("the value of fixedLineY:idx ?\n");
					return false;
				}
				idx = atoi(argv[k]);
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
				fixed_line_y.push_back(FixedLine(idx,x1,y1,x2,y2));
			}
			else if(_strcmpi(argv[k],"fixedPoint") == 0)
			{
				int cx,cy;
				float px,py;
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
				fixed_pts.push_back(FixedPoint(cx,cy,px,py));
			}
			else
			{
				printf("unknown parameter %s\n",argv[k]);
				return false;
			}
		}

		if(GridWidth < 2 || GridHeight < 2)
		{
			printf("invalid grid resolution\n");
			return false;
		}
		nouseful_flag.allocate(GridWidth,GridHeight);
		fixed_flag.allocate(GridWidth,GridHeight);
		init_coord.allocate(GridWidth,GridHeight,2);
		bool*& nouseful_flag_data = nouseful_flag.data();
		bool*& fixed_flag_data = fixed_flag.data();
		float*& init_coord_data = init_coord.data();
		for(int i = 0;i < nouse_seg_x.size();i++)
		{
			int idx = nouse_seg_x[i].idx;
			if(idx < 0 || idx >= GridHeight)
			{
				printf("invalid nouseSegX:idx=(:,%d)\n",idx);
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
				nouseful_flag_data[idx*GridWidth+w] = true;
		}
		for(int i = 0;i < nouse_seg_y.size();i++)
		{
			int idx = nouse_seg_y[i].idx;
			if(idx < 0 || idx >= GridWidth)
			{
				printf("invalid nouseSegY:idx=(%d,:)\n",idx);
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
				nouseful_flag_data[h*GridWidth+idx] = true;
		}
		for(int i = 0;i < fixed_seg_x.size();i++)
		{
			int idx = fixed_seg_x[i].idx;
			if(idx < 0 || idx >= GridHeight)
			{
				printf("invalid fixedSegX:idx=(:,%d)\n",idx);
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
			float x2 = fixed_seg_x[i].x2;
			float y2 = fixed_seg_x[i].y2;
			for(int w = start;w <= end;w++)
			{
				fixed_flag_data[idx*GridWidth+w] = true;
				float weight = (float)(w-start)/(end-start);
				init_coord_data[(idx*GridWidth+w)*2+0] = x1 + weight*(x2-x1);
				init_coord_data[(idx*GridWidth+w)*2+1] = y1 + weight*(y2-y1);
			}
		}
		for(int i = 0;i < fixed_seg_y.size();i++)
		{
			int idx = fixed_seg_y[i].idx;
			if(idx < 0 || idx >= GridWidth)
			{
				printf("invalid fixedSegY:idx=(%d,:)\n",idx);
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
			float x2 = fixed_seg_y[i].x2;
			float y2 = fixed_seg_y[i].y2;
			for(int h = start;h <= end;h++)
			{
				fixed_flag_data[h*GridWidth+idx] = true;
				float weight = (float)(h-start)/(end-start);
				init_coord_data[(h*GridWidth+idx)*2+0] = x1 + weight*(x2-x1);
				init_coord_data[(h*GridWidth+idx)*2+1] = y1 + weight*(y2-y1);
			}
		}
		for(int i = 0;i < fixed_line_x.size();i++)
		{
			int idx = fixed_line_x[i].idx;
			if(idx < 0 || idx >= GridHeight)
			{
				printf("invalid fixedLineX:idx=%d\n",idx);
				return false;
			}
			float x1 = fixed_line_x[i].x1;
			float y1 = fixed_line_x[i].y1;
			float x2 = fixed_line_x[i].x2;
			float y2 = fixed_line_x[i].y2;
			for(int w = 0;w < GridWidth;w++)
			{
				fixed_flag_data[idx*GridWidth+w] = true;
				float weight = (float)w/(GridWidth-1);
				init_coord_data[(idx*GridWidth+w)*2+0] = x1 + weight*(x2-x1);
				init_coord_data[(idx*GridWidth+w)*2+1] = y1 + weight*(y2-y1);
			}
		}
		for(int i = 0;i < fixed_line_y.size();i++)
		{
			int idx = fixed_line_y[i].idx;
			if(idx < 0 || idx >= GridWidth)
			{
				printf("invalid fixedLineY:idx=%d\n",idx);
				return false;
			}
			float x1 = fixed_line_y[i].x1;
			float y1 = fixed_line_y[i].y1;
			float x2 = fixed_line_y[i].x2;
			float y2 = fixed_line_y[i].y2;
			for(int h = 0;h < GridHeight;h++)
			{
				fixed_flag_data[h*GridWidth+idx] = true;
				float weight = (float)h/(GridHeight-1);
				init_coord_data[(h*GridWidth+idx)*2+0] = x1 + weight*(x2-x1);
				init_coord_data[(h*GridWidth+idx)*2+1] = y1 + weight*(y2-y1);
			}
		}
		for(int i = 0;i < fixed_pts.size();i++)
		{
			int cx = fixed_pts[i].coord_x;
			int cy = fixed_pts[i].coord_y;
			if(cx < 0 || cx >= GridWidth || cy < 0 || cy >= GridHeight)
			{
				printf("invalid fixedPoint: <x,y> = <%d,%d>\n",cx,cy);
				return false;
			}
			float px = fixed_pts[i].pos_x;
			float py = fixed_pts[i].pos_y;
			fixed_flag_data[cy*GridWidth+cx] = true;
			init_coord_data[(cy*GridWidth+cx)*2+0] = px;
			init_coord_data[(cy*GridWidth+cx)*2+1] = py;
		}
		return true;
	}
};

class WarpGrid
{
public:
	static bool go(const WarpGridOptions& opt)
	{
		int width = opt.GridWidth;
		int height = opt.GridHeight;
		if(!opt.nouseful_flag.matchDimension(width,height,1)
			|| !opt.fixed_flag.matchDimension(width,height,1)
			|| !opt.init_coord.matchDimension(width,height,2))
		{
			printf("dimensions don't match\n");
			return false;
		}

		if(opt.has_map_mask_file)
		{
			ZQ::ZQ_DImage<float> obj_mask(width,height);
			for(int i = 0;i < width*height;i++)
				obj_mask.data()[i] = opt.nouseful_flag.data()[i] ? 0 : 1;
			if(!ZQ::ZQ_ImageIO::saveImage(obj_mask,opt.map_mask_file))
			{
				printf("failed to save %s\n",opt.map_mask_file);
				return false;
			}
		}
		ZQ::ZQ_DImage<float> out_coord(width,height,2);
		ZQ::ZQ_GridDeformation<float> deform;
		ZQ::ZQ_GridDeformationOptions w_opt;
		if(opt.distance_weight != 0)
		{
			if(opt.xloop)
				w_opt.methodType = ZQ::ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_DISTANCE_ENERGY_XLOOP;
			else
				w_opt.methodType = ZQ::ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_DISTANCE_ENERGY;
		}
		else
		{
			if(opt.xloop)
				w_opt.methodType = ZQ::ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_ENERGY_XLOOP;
			else
				w_opt.methodType = ZQ::ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_ENERGY;
		}
		w_opt.line_weight = opt.line_weight;
		w_opt.angle_weight = opt.angle_weight;
		w_opt.distance_weight = opt.distance_weight;
		w_opt.distance = opt.distance;
		w_opt.iteration = opt.max_iteration;
		w_opt.FPIteration = opt.FPiteration;
		clock_t build_st = clock();
		if(!deform.BuildMatrix(width,height,opt.nouseful_flag.data(),opt.fixed_flag.data(),w_opt))
		{
			printf("failed to build matrix for deformation\n");
			return false;
		}
		clock_t build_end = clock();
		printf("build cost: %f seconds\n",0.001*(build_end-build_st));

		clock_t deform_st = clock();
		if(!deform.Deformation(opt.init_coord.data(),out_coord.data()))
		{
			printf("failed to deform\n");
			return false;
		}
		clock_t deform_end = clock();
		printf("build cost: %f seconds\n",0.001*(deform_end-deform_st));

		if(opt.has_show_file)
		{
			
			int show_scale = opt.show_scale;
			int back_width = opt.BoundingBoxmax[0] - opt.BoundingBoxmin[0];
			int back_height = opt.BoundingBoxmax[1] - opt.BoundingBoxmin[1];
			int show_shift_x = (0.5-opt.BoundingBoxmin[0])*show_scale;
			int show_shift_y = (0.5-opt.BoundingBoxmin[1])*show_scale;

			IplImage* show_img = cvCreateImage(cvSize(back_width*show_scale,back_height*show_scale),IPL_DEPTH_8U,3);
			if(show_img == 0)
			{
				printf("failed to create background image\n");
				return false;
			}

			bool retflag = RenderGridWithTexture::render_wireframe(show_img,out_coord,opt.nouseful_flag,show_scale,show_shift_x,show_shift_y,opt.xloop);
			cvSaveImage(opt.show_file,show_img);
			cvReleaseImage(&show_img);
			if(!retflag)
			{
				printf("failed to render\n");
				return false;
			}
		}

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