#ifndef _TESSELLATION_H_
#define _TESSELLATION_H_
#pragma once

#include "ZQ_DoubleImage.h"
#include "ZQ_GridDeformation.h"

class TessellationOptions
{
public:
	TessellationOptions()
	{
		Reset();
	}
	~TessellationOptions()
	{

	}

	enum CONST_VAL{FILE_LEN = 200};
	bool has_in_map_file;
	char in_map_file[FILE_LEN];
	bool has_in_mask_file;
	char in_mask_file[FILE_LEN];
	bool has_out_map_file;
	char out_map_file[FILE_LEN];
	bool has_out_mask_file;
	char out_mask_file[FILE_LEN];
	float line_weight;
	float angle_weight;
	float distance_weight;
	float distance;
	bool xloop;
	float max_iteration;
	float FPiteration;
	bool display;
	int tesse_scale;

	void Reset()
	{
		has_in_map_file = false;
		in_map_file[0] = '\0';
		has_in_mask_file = false;
		in_mask_file[0] = '\0';
		has_out_map_file = false;
		out_map_file[0] = '\0';
		has_out_mask_file = false;
		out_mask_file[0] = '\0';
		line_weight = 1;
		angle_weight = 0.1;
		distance_weight = 0;
		distance = 1;
		max_iteration = 1000;
		FPiteration = 5;
		xloop = false;
		display = false;
		tesse_scale = 4;
	}

	bool HandleArgs(int argc, const char** argv)
	{
		for(int k = 0;k < argc;k++)
		{
			if(_strcmpi(argv[k],"in_map_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of in_map_file ?\n");
					return false;
				}
				strcpy(in_map_file,argv[k]);
				has_in_map_file = true;
			}
			else if(_strcmpi(argv[k],"in_mask_fiel") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of in_mask_file ?\n");
					return false;
				}
				strcpy(in_mask_file,argv[k]);
				has_in_mask_file = true;
			}
			else if(_strcmpi(argv[k],"out_map_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of out_map_file ?\n");
					return false;
				}
				strcpy(out_map_file,argv[k]);
				has_out_map_file = true;
			}
			else if(_strcmpi(argv[k],"out_mask_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of out_mask_file ?\n");
					return false;
				}
				strcpy(out_mask_file,argv[k]);
				has_out_mask_file = true;
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

			else if(_strcmpi(argv[k],"display") == 0)
			{
				display = true;
			}
			else if(_strcmpi(argv[k],"tesse_scale") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of tesse_scale ?\n");
					return false;
				}
				tesse_scale = atoi(argv[k]);
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

class Tessellation
{
public:
	static bool Go(const TessellationOptions& opt)
	{
		if(!opt.has_in_map_file)
		{
			if(opt.display)
				printf("need in map file!\n");
			return false;
		}
		ZQ::ZQ_DImage<float> in_map,in_mask;
		ZQ::ZQ_DImage<bool> in_nouseful_flag;
		if(!in_map.loadImage(opt.in_map_file))
		{
			if(opt.display)
				printf("failed to open %s\n",opt.in_map_file);
			return false;
		}
		int width = in_map.width();
		int height = in_map.height();
		if(in_map.nchannels() != 2)
		{
			if(opt.display)
				printf("map should has 2 channels\n");
			return false;
		}

		bool*& in_nouseful_flag_data = in_nouseful_flag.data();
		if(!opt.has_in_mask_file)
		{
			in_nouseful_flag.allocate(width,height);
		}
		else
		{
			if(!ZQ::ZQ_ImageIO::loadImage(in_mask,opt.in_mask_file,0))
			{
				if(opt.display)
					printf("failed to open %s\n",opt.in_mask_file);
				return false;
			}
			if(!in_mask.matchDimension(width,height,1))
			{
				if(opt.display)
					printf("map and mask should match in dimension\n");
				return false;
			}
			in_nouseful_flag.allocate(width,height);
			float*& in_mask_data = in_mask.data();
			for(int i = 0;i < width*height;i++)
				in_nouseful_flag_data[i] = in_mask_data[i] < 0.5;
		}

		if(opt.tesse_scale <= 1)
		{
			if(opt.display)
				printf("tesse_scale should at least 2\n");
			return false;
		}

		if(!opt.has_out_map_file)
		{
			if(opt.display)
				printf("need out_map_file\n");
			return false;
		}

		/////////

		int out_width,out_height;
		ZQ::ZQ_DImage<float> init_coord;
		ZQ::ZQ_DImage<bool> nouseful_flag;
		ZQ::ZQ_DImage<bool> fixed_flag;
		bool*& nouseful_data = nouseful_flag.data();
		bool*& fixed_data = fixed_flag.data();
		float*& init_coord_data = init_coord.data();
		float*& in_map_data = in_map.data();
		int tesse_scale = opt.tesse_scale;
		if(opt.xloop)
		{
			out_width = width*opt.tesse_scale;
			out_height = (height-1)*opt.tesse_scale+1;
			init_coord.allocate(out_width,out_height,2);
			nouseful_flag.allocate(out_width,out_height,1);
			fixed_flag.allocate(out_width,out_height,1);
		}
		else
		{
			out_width = (width-1)*opt.tesse_scale+1;
			out_height = (height-1)*opt.tesse_scale+1;
			init_coord.allocate(out_width,out_height,2);
			nouseful_flag.allocate(out_width,out_height,1);
			fixed_flag.allocate(out_width,out_height,1);
		}
		for(int i = 0;i < out_height;i++)
		{
			for(int j = 0;j < out_width;j++)
			{
				int i0 = i/tesse_scale;
				int i1 = i0+1;
				int j0 = j/tesse_scale;
				int j1 = opt.xloop ? (j0+1)%width : j0+1;
				int off_i = i%tesse_scale;
				int off_j = j%tesse_scale;
				fixed_data[i*out_width+j] = (off_i == 0)&&(off_j==0);
				if(off_i == 0)
				{
					if(off_j == 0)
					{
						nouseful_data[i*out_width+j] = in_nouseful_flag_data[i0*width+j0];
						init_coord_data[(i*out_width+j)*2+0] = in_map_data[(i0*width+j0)*2+0];
						init_coord_data[(i*out_width+j)*2+1] = in_map_data[(i0*width+j0)*2+1];
					}
					else
					{
						nouseful_data[i*out_width+j] = in_nouseful_flag_data[i0*width+j0] && in_nouseful_flag_data[i0*width+j1];
						float weight2 = (float)off_j/(float)tesse_scale;
						float weight1 =1.0f-weight2;
						init_coord_data[(i*out_width+j)*2+0] = weight1*in_map_data[(i0*width+j0)*2+0] + weight2*in_map_data[(i0*width+j1)*2+0];
						init_coord_data[(i*out_width+j)*2+1] = weight1*in_map_data[(i0*width+j0)*2+1] + weight2*in_map_data[(i0*width+j1)*2+1];
					}
				}
				else
				{
					if(off_j == 0)
					{
						nouseful_data[i*out_width+j] = in_nouseful_flag_data[i0*width+j0] && in_nouseful_flag_data[i1*width+j0];
						float weight2 = (float)off_i/(float)tesse_scale;
						float weight1 =1.0f-weight2;
						init_coord_data[(i*out_width+j)*2+0] = weight1*in_map_data[(i0*width+j0)*2+0] + weight2*in_map_data[(i1*width+j0)*2+0];
						init_coord_data[(i*out_width+j)*2+1] = weight1*in_map_data[(i0*width+j0)*2+1] + weight2*in_map_data[(i1*width+j0)*2+1];
					}
					else
					{
						nouseful_data[i*out_width+j] = in_nouseful_flag_data[i0*width+j0] && in_nouseful_flag_data[i0*width+j1]
						&& in_nouseful_flag_data[i1*width+j0] && in_nouseful_flag_data[i1*width+j1];
						float s_i = (float)off_i/(float)tesse_scale;
						float s_j = (float)off_j/(float)tesse_scale;
						init_coord_data[(i*out_width+j)*2+0] = (1.0-s_i)*(1.0-s_j)*in_map_data[(i0*width+j0)*2+0] + (1.0-s_i)*s_j*in_map_data[(i0*width+j1)*2+0]
						+ s_i*(1.0-s_j)*in_map_data[(i1*width+j0)*2+0] + s_i*s_j*in_map_data[(i1*width+j1)*2+0];
						init_coord_data[(i*out_width+j)*2+1] = (1.0-s_i)*(1.0-s_j)*in_map_data[(i0*width+j0)*2+1] + (1.0-s_i)*s_j*in_map_data[(i0*width+j1)*2+1]
						+ s_i*(1.0-s_j)*in_map_data[(i1*width+j0)*2+1] + s_i*s_j*in_map_data[(i1*width+j1)*2+1];
					}
				}
			}
		}

		/////////

		if(opt.has_out_mask_file)
		{
			ZQ::ZQ_DImage<float> obj_mask(out_width,out_height);
			for(int i = 0;i < out_width*out_height;i++)
				obj_mask.data()[i] = nouseful_data[i] ? 0 : 1;
			if(!ZQ::ZQ_ImageIO::saveImage(obj_mask,opt.out_mask_file))
			{
				printf("failed to save %s\n",opt.out_mask_file);
				return false;
			}
		}
		ZQ::ZQ_DImage<float> out_coord(out_width,out_height,2);
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
		w_opt.distance = opt.distance/(float)tesse_scale;
		w_opt.iteration = opt.max_iteration;
		w_opt.FPIteration = opt.FPiteration;
		clock_t build_st = clock();
		if(!deform.BuildMatrix(out_width,out_height,nouseful_data,fixed_data,w_opt))
		{
			printf("failed to build matrix for deformation\n");
			return false;
		}
		clock_t build_end = clock();
		printf("build cost: %f seconds\n",0.001*(build_end-build_st));

		clock_t deform_st = clock();
		if(!deform.Deformation(init_coord_data,out_coord.data(),true))
		{
			printf("failed to deform\n");
			return false;
		}
		clock_t deform_end = clock();
		printf("deform cost: %f seconds\n",0.001*(deform_end-deform_st));


		if(!out_coord.saveImage(opt.out_map_file))
		{
			printf("failed to save %s\n",opt.out_map_file);
			return false;
		}
		return true;
	}
};

#endif