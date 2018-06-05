#ifndef _RESAMPLE_MAP_RBF_3D_H_
#define _RESAMPLE_MAP_RBF_3D_H_
#pragma once

#include "ZQ_DoubleImage3D.h"
#include "ZQ_ScatteredInterpolationRBF.h"
#include "RenderGridWithTexture3D.h"

class ResampleMapRBF3DOptions
{
public:
	ResampleMapRBF3DOptions()
	{
		Reset();
	}
	~ResampleMapRBF3DOptions()
	{

	}

	enum CONST_VAL{FILE_LEN = 200};
	bool has_map_file;
	char map_file[FILE_LEN];
	bool has_show_file;
	char show_file[FILE_LEN];
	bool xloop;
	float xloop_thresh;
	int BoundingBoxmin[3];
	int BoundingBoxmax[3];
	int show_scale;
	int rbf_neighbor_num;
	float rbf_radius;
	bool display;
	int border_size;
	bool with_linear_mask;

	void Reset()
	{
		has_map_file = false;
		map_file[0] = '\0';
		has_show_file = false;
		show_file[0] = '\0';
		xloop = false;
		xloop_thresh = 0.05;
		BoundingBoxmin[0] = BoundingBoxmin[1] = BoundingBoxmin[2] = 0;
		BoundingBoxmax[0] = BoundingBoxmax[1] = BoundingBoxmax[2] = 0;
		show_scale = 10;
		rbf_neighbor_num = 4;
		rbf_radius = 3;
		display = false;
		border_size = 0;
		with_linear_mask = false;
	}

	bool HandleArgs(int argc, const char** argv)
	{
		for(int k = 0;k < argc;k++)
		{
			if(_strcmpi(argv[k],"map_file") == 0)
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
			else if(_strcmpi(argv[k],"xloop") == 0)
			{
				xloop = true;
			}
			else if(_strcmpi(argv[k],"xloop_thresh") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of xloop_thresh ?\n");
					return false;
				}
				xloop_thresh = atof(argv[k]);
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
			else if(_strcmpi(argv[k],"rbf_neighbor_num") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of rbf_neighbor_num ?\n");
					return false;
				}
				rbf_neighbor_num = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"rbf_radius") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of rbf_radius ?\n");
					return false;
				}
				rbf_radius = atof(argv[k]);
			}
			else if(_strcmpi(argv[k],"display") == 0)
			{
				display = true;
			}
			else if(_strcmpi(argv[k],"border_size") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of border_size ?\n");
					return false;
				}
				border_size = atoi(argv[k]);
			}
			else if(_strcmpi(argv[k],"with_linear_mask") == 0)
			{
				with_linear_mask = true;
			}
			else
			{
				printf("unknown parameters %s\n",argv[k]);
				return false;
			}
		}
		return true;
	}
};


class ResampleMapRBF3D
{
public:
	static bool go(const ResampleMapRBF3DOptions& opt)
	{
		if(!opt.has_map_file)
		{
			printf("error: need map_file!\n");
			return false;
		}
		if(!opt.has_show_file)
		{
			printf("error: need show_file!\n");
			return false;
		}
		if(opt.show_scale <= 0)
		{
			printf("error: invalid show_scale!\n");
			return false;
		}

		int show_scale = opt.show_scale;
		int back_width = opt.BoundingBoxmax[0] - opt.BoundingBoxmin[0];
		int back_height = opt.BoundingBoxmax[1] - opt.BoundingBoxmin[1];
		int back_depth = opt.BoundingBoxmax[2] - opt.BoundingBoxmin[2];
		float show_shift_x = (0-opt.BoundingBoxmin[0])*show_scale;
		float show_shift_y = (0-opt.BoundingBoxmin[1])*show_scale;
		float show_shift_z = (0-opt.BoundingBoxmin[2])*show_scale;
		if(back_width <= 0 || back_height <= 0 || back_depth <= 0)
		{
			printf("invalid boundingbox !\n");
			return false;
		}
		back_width *= show_scale;
		back_height *= show_scale;
		back_depth *= show_scale;

		ZQ::ZQ_DImage3D<float> map_img;
		if(!map_img.loadImage(opt.map_file))
		{
			printf("failed to load %s\n",opt.map_file);
			return false;
		}

		ZQ::ZQ_DImage3D<float> show_map(back_width,back_height,back_depth,3);
		ZQ::ZQ_DImage3D<float> show_mask(back_width,back_height,back_depth,1);

		if(opt.xloop)
		{
			if(!_resample_map_xloop(show_map,show_mask,map_img,show_scale,show_shift_x,show_shift_y,show_shift_z,
				opt.rbf_neighbor_num,opt.rbf_radius,opt.xloop_thresh,opt.border_size,opt.display))
				return false;
		}
		else
		{
			if(!_resample_map(show_map,show_mask,map_img,show_scale,show_shift_x,show_shift_y,show_shift_z,
				opt.rbf_neighbor_num,opt.rbf_radius,opt.border_size,opt.display))
				return false;
		}

		if(!show_map.saveImage(opt.show_file))
		{
			printf("failed to save %s\n",opt.show_file);
			return false;
		}

		if(opt.with_linear_mask)
		{

			int map_width = map_img.width();
			int map_height = map_img.height();
			int map_depth = map_img.depth();
			ZQ::ZQ_DImage3D<float> tex_linear(2,2,2,1);
			ZQ::ZQ_DImage3D<bool> obj_mask(map_width,map_height,map_depth,1);
			ZQ::ZQ_DImage3D<float> mask_linear(back_width,back_height,back_depth,1);
			float*& tex_linear_data = tex_linear.data();
			tex_linear_data[0] = 1; tex_linear_data[1] = 1; tex_linear_data[2] = 1; tex_linear_data[3] = 1;
			tex_linear_data[4] = 1; tex_linear_data[5] = 1; tex_linear_data[6] = 1; tex_linear_data[7] = 1;
			if(!RenderGridWithTexture3D::render_with_texture(mask_linear,map_img,tex_linear,obj_mask,show_scale,show_shift_x,show_shift_y,show_shift_z,opt.xloop,0))
			{
				return false;
			}
			float*& mask_linear_data = mask_linear.data();
			float*& show_mask_data = show_mask.data();
			for(int i = 0;i < back_width*back_height*back_depth;i++)
				show_mask_data[i] = (show_mask_data[i] > 0.5 && mask_linear_data[i] > 0.5) ? 1 : 0;
		}

		std::string mask_filename;
		mask_filename.append(opt.show_file);
		mask_filename.append("_mask.di3");
		if(!show_mask.saveImage(mask_filename.c_str()))
		{
			printf("failed to save %s\n",mask_filename.c_str());
			return false;
		}

		return true;
	}

	static bool _resample_map_xloop(ZQ::ZQ_DImage3D<float>& show_map, ZQ::ZQ_DImage3D<float>& mask, const ZQ::ZQ_DImage3D<float>& map_img, int show_scale, float show_shift_x, float show_shift_y, float show_shift_z,
		int rbf_neighbor_num, float rbf_radius, float xloop_thresh, int boarder_size, bool display)
	{
		if(3 != show_map.nchannels() || 3 != map_img.nchannels())
			return false;

		if(boarder_size < 0)
			return false;

		int show_width = show_map.width();
		int show_height = show_map.height();
		int show_depth = show_map.depth();
		int map_width = map_img.width();
		int map_height = map_img.height();
		int map_depth = map_img.depth();
		if(map_width <= 1 || map_height <= 1 + 2*boarder_size || map_depth <= 1 + 2*boarder_size)
			return false;

		if(!mask.matchDimension(show_width,show_height,show_depth,1))
			mask.allocate(show_width,show_height,show_depth);

		int npts = map_width*map_height*map_depth;
		const float*& map_img_data = map_img.data();
		ZQ::ZQ_DImage<float> pts_img(npts,3);
		ZQ::ZQ_DImage<float> vals_img(npts,4);

		const double m_pi = atan(1.0)*4;

		float*& pts = pts_img.data();
		float*& vals = vals_img.data();
		for(int k = 0;k < map_depth;k++)
		{
			for(int j = 0;j < map_height;j++)
			{
				for(int i = 0;i < map_width;i++)
				{
					int offset = k*map_height*map_width + j*map_width + i;
					pts[offset*3+0] = map_img_data[offset*3+0]*show_scale + show_shift_x;
					pts[offset*3+1] = map_img_data[offset*3+1]*show_scale + show_shift_y;
					pts[offset*3+2] = map_img_data[offset*3+2]*show_scale + show_shift_z;
					double angle = (double)i/map_width*2*m_pi;
					vals[offset*4+0] = cos(angle);
					vals[offset*4+1] = sin(angle);
					vals[offset*4+2] = (double)(j-boarder_size)/(map_height-1-2*boarder_size);
					vals[offset*4+3] = (double)(k-boarder_size)/(map_depth-1-2*boarder_size);
				}
			}
		}

		ZQ::ZQ_ScatteredInterpolationRBF<float> rbf;
		if(!rbf.SetLandmarks(npts,3,pts,vals,4))
			return false;
		if(!rbf.SolveCoefficient(rbf_neighbor_num,rbf_radius,npts*2,ZQ::ZQ_RBFKernel::COMPACT_CPC2,display))
			return false;

		ZQ::ZQ_DImage3D<float> tmp_coord(show_width,show_height,show_depth,4);
		float*& tmp_coord_data = tmp_coord.data();
		if(!rbf.GridData3D(show_width,show_height,show_depth,0,show_width-1,0,show_height-1,0,show_depth-1,tmp_coord_data))
		{
			return false;
		}

		float*& mask_data = mask.data();
		float*& show_map_data = show_map.data();

		for(int i = 0;i < show_width*show_height*show_depth;i++)
		{
			double cos_x = tmp_coord_data[i*4+0];
			double sin_x = tmp_coord_data[i*4+1];
			float coord_x = _get_the_angle(cos_x,sin_x)/(2*m_pi);
			float coord_y = tmp_coord_data[i*4+2];
			float coord_z = tmp_coord_data[i*4+3];
			float radius = sqrt(sin_x*sin_x+cos_x*cos_x);
			mask_data[i] = (coord_y <= 0 || coord_y >= 1 || coord_z <= 0 || coord_z >= 1 ||  fabs(radius-1) > xloop_thresh) ? 0 : 1;
		
			show_map_data[i*3+0] = coord_x;
			show_map_data[i*3+1] = coord_y;
			show_map_data[i*3+2] = coord_z;
		}
		return true;

	}

	static bool _resample_map(ZQ::ZQ_DImage3D<float>& show_map, ZQ::ZQ_DImage3D<float>& mask, const ZQ::ZQ_DImage3D<float>& map_img, int show_scale, float show_shift_x, float show_shift_y, float show_shift_z, 
		int rbf_neighbor_num, float rbf_radius, int boarder_size, bool display)
	{
		if(3 != show_map.nchannels() || 3 != map_img.nchannels())
			return false;

		if(boarder_size < 0)
			return false;

		int show_width = show_map.width();
		int show_height = show_map.height();
		int show_depth = show_map.depth();
		int map_width = map_img.width();
		int map_height = map_img.height();
		int map_depth = map_img.depth();
		if(map_width <= 1+2*boarder_size || map_height <= 1+2*boarder_size || map_depth <= 1+2*boarder_size)
			return false;

		if(!mask.matchDimension(show_width,show_height,show_depth,1))
			mask.allocate(show_width,show_height,show_depth);

		int npts = map_width*map_height*map_depth;
		const float*& map_img_data = map_img.data();
		ZQ::ZQ_DImage3D<float> pts_img(npts,1,1,3);
		ZQ::ZQ_DImage3D<float> vals_img(npts,1,1,3);
		float*& pts = pts_img.data();
		float*& vals = vals_img.data();
		for(int k = 0;k < map_depth;k++)
		{
			for(int j = 0;j < map_height;j++)
			{
				for(int i = 0;i < map_width;i++)
				{
					int offset = k*map_height*map_width + j*map_width + i;
					pts[offset*3+0] = map_img_data[offset*3+0]*show_scale + show_shift_x;
					pts[offset*3+1] = map_img_data[offset*3+1]*show_scale + show_shift_y;
					pts[offset*3+2] = map_img_data[offset*3+2]*show_scale + show_shift_z;
					vals[offset*3+0] = (double)(i-boarder_size)/(map_width-1-2*boarder_size);
					vals[offset*3+1] = (double)(j-boarder_size)/(map_height-1-2*boarder_size);
					vals[offset*3+2] = (double)(k-boarder_size)/(map_depth-1-2*boarder_size);
				}
			}
		}
		
		ZQ::ZQ_ScatteredInterpolationRBF<float> rbf;
		if(!rbf.SetLandmarks(npts,3,pts,vals,3))
			return false;
		if(!rbf.SolveCoefficient(rbf_neighbor_num,rbf_radius,npts*2,ZQ::ZQ_RBFKernel::COMPACT_CPC2,display))
			return false;

		float*& show_map_data = show_map.data();
		if(!rbf.GridData3D(show_width,show_height,show_depth,0,show_width-1,0,show_height-1,0,show_depth-1,show_map_data))
		{
			return false;
		}

		float*& mask_data = mask.data();

		for(int i = 0;i < show_width*show_height*show_depth;i++)
		{
			float coord_x = show_map_data[i*3+0];
			float coord_y = show_map_data[i*3+1];
			float coord_z = show_map_data[i*3+2];
			mask_data[i] = (coord_x <= 0 || coord_x >= 1 || coord_y <= 0 || coord_y >= 1 || coord_z <= 0 || coord_z >= 1) ? 0 : 1;
		}
		return true;
	}

private:
	static double _get_the_angle(double cos_x, double sin_x)
	{
		double len = sqrt(sin_x*sin_x+cos_x*cos_x);
		if(len == 0)
			return 0;
		cos_x /= len;
		sin_x /= len;

		const double m_pi = atan(1.0)*4;

		double angle = acos(cos_x);
		return sin_x >= 0 ? angle : (2*m_pi-angle);

	}

};


#endif