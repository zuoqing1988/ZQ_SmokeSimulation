#ifndef _RENDER_GRID_WITH_TEXTURE_H_
#define _RENDER_GRID_WITH_TEXTURE_H_
#pragma once

#include "ZQ_DoubleImage.h"
#include "ZQ_ScanLinePolygonFill.h"
#include "ZQ_TextureSampler.h"
#include "ZQ_ImageIO.h"
#include "ZQ_Vec2D.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

class RenderGridWithTextureOptions
{
public:
	RenderGridWithTextureOptions()
	{
		Reset();
	}
	~RenderGridWithTextureOptions()
	{

	}

	enum CONST_VAL{FILE_LEN = 200};
	bool has_map_file;
	char map_file[FILE_LEN];
	bool has_map_mask_file;
	char map_mask_file[FILE_LEN];
	bool has_tex_file;
	char tex_file[FILE_LEN];
	bool has_show_file;
	char show_file[FILE_LEN];
	bool xloop;
	int BoundingBoxmin[2];
	int BoundingBoxmax[2];
	int show_scale;
	int border_size;

	void Reset()
	{
		has_map_file = false;
		map_file[0] = '\0';
		has_map_mask_file = false;
		map_mask_file[0] = '\0';
		has_tex_file = false;
		tex_file[0] = '\0';
		has_show_file = false;
		show_file[0] = '\0';
		xloop = false;
		BoundingBoxmin[0] = BoundingBoxmin[1] = 0;
		BoundingBoxmax[0] = BoundingBoxmax[1] = 0;
		show_scale = 10;
		border_size = 0;
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
			else if(_strcmpi(argv[k],"tex_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of tex_file ?\n");
					return false;
				}
				strcpy(tex_file,argv[k]);
				has_tex_file = true;
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
			else
			{
				printf("unknown parameters %s\n",argv[k]);
				return false;
			}
		}
		return true;
	}
};

class RenderGridWithTexture
{
public:
	static bool go(const RenderGridWithTextureOptions& opt)
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
		float show_shift_x = (0.5-opt.BoundingBoxmin[0])*show_scale;
		float show_shift_y = (0.5-opt.BoundingBoxmin[1])*show_scale;
		if(back_width <= 0 || back_height <= 0)
		{
			printf("invalid boundingbox !\n");
			return false;
		}
		back_width *= show_scale;
		back_height *= show_scale;

		ZQ::ZQ_DImage<float> map_img;
		if(!map_img.loadImage(opt.map_file))
		{
			printf("failed to load %s\n",opt.map_file);
			return false;
		}
		ZQ::ZQ_DImage<bool> obj_mask(map_img.width(),map_img.height());
		if(opt.has_map_mask_file)
		{
			ZQ::ZQ_DImage<float> mask;
			if(!ZQ::ZQ_ImageIO::loadImage(mask,opt.map_mask_file,0))
			{
				printf("failed to load %s\n",opt.map_mask_file);
				return false;
			}
			if(mask.width() != map_img.width() || mask.height() != map_img.height())
			{
				printf("mask and map don't match in dimension\n");
				return false;
			}
			for(int i = 0;i < obj_mask.npixels();i++)
			{
				obj_mask.data()[i] = mask.data()[i] > 0.5;
			}
		}

		if(!opt.has_tex_file)
		{
			IplImage* render_img = cvCreateImage(cvSize(back_width,back_height),IPL_DEPTH_8U,3);
			if(render_img == 0)
			{
				printf("failed to create render image\n");
				return false;
			}
			cvZero(render_img);
			
			if(!render_wireframe(render_img,map_img,obj_mask,show_scale,show_shift_x,show_shift_y,opt.xloop))
				return false;
			cvSaveImage(opt.show_file,render_img);
			cvReleaseImage(&render_img);
		}
		else
		{
			ZQ::ZQ_DImage<float> tex_img;
			if(!ZQ::ZQ_ImageIO::loadImage(tex_img,opt.tex_file,1))
			{
				printf("failed to load %s\n",opt.tex_file);
				return false;
			}
			int nChannels = tex_img.nchannels();
			ZQ::ZQ_DImage<float> render_img(back_width,back_height,nChannels);
			if(!render_with_texture(render_img,map_img,tex_img,obj_mask,show_scale,show_shift_x,show_shift_y,opt.xloop,opt.border_size))
				return false;
			if(!ZQ::ZQ_ImageIO::saveImage(render_img,opt.show_file))
			{
				printf("failed to save %s\n",opt.show_file);
				return false;
			}
		}

		return true;
	}

	static bool render_wireframe(IplImage* render_img, const ZQ::ZQ_DImage<float>& map_img, const ZQ::ZQ_DImage<bool> obj_mask, 
		const int show_scale, float show_shift_x, float show_shift_y, bool xloop)
	{
		int width = map_img.width();
		int height = map_img.height();
		if(!map_img.matchDimension(width,height,2))
			return false;
		if(!obj_mask.matchDimension(width,height,1))
			return false;
		if(render_img == 0)
			return false;

		CvScalar color_line = cvScalar(0,200,0);
		CvScalar color_point = cvScalar(0,250,0);
		int thick_line = 1;
		int thick_point = 2;
		int size_point = 1;

		const bool*& mask_data = obj_mask.data();
		const float*& coord_data = map_img.data();
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(!mask_data[i*width+j])
				{
					float pos_x = coord_data[(i*width+j)*2+0];
					float pos_y = coord_data[(i*width+j)*2+1];
					cvCircle(render_img,cvPoint(pos_x*show_scale+show_shift_x,pos_y*show_scale+show_shift_y),size_point,color_point,thick_point);
				}
			}
		}
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width-1;j++)
			{
				if(!mask_data[i*width+j] && !mask_data[i*width+j+1])
				{
					float pos_x1 = coord_data[(i*width+j)*2+0];
					float pos_y1 = coord_data[(i*width+j)*2+1];
					float pos_x2 = coord_data[(i*width+j+1)*2+0];
					float pos_y2 = coord_data[(i*width+j+1)*2+1];
					cvLine(render_img,cvPoint(pos_x1*show_scale+show_shift_x,pos_y1*show_scale+show_shift_y),
						cvPoint(pos_x2*show_scale+show_shift_x,pos_y2*show_scale+show_shift_y),
						color_line,thick_point);
				}
			}
			if(xloop)
			{
				if(!mask_data[i*width+0] && !mask_data[i*width+width-1])
				{
					float pos_x1 = coord_data[(i*width+0)*2+0];
					float pos_y1 = coord_data[(i*width+0)*2+1];
					float pos_x2 = coord_data[(i*width+width-1)*2+0];
					float pos_y2 = coord_data[(i*width+width-1)*2+1];
					cvLine(render_img,cvPoint(pos_x1*show_scale+show_shift_x,pos_y1*show_scale+show_shift_y),
						cvPoint(pos_x2*show_scale+show_shift_x,pos_y2*show_scale+show_shift_y),
						color_line,thick_point);
				}
			}
		}
		for(int i = 0;i < height-1;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(!mask_data[i*width+j] && !mask_data[i*width+j+width])
				{
					float pos_x1 = coord_data[(i*width+j)*2+0];
					float pos_y1 = coord_data[(i*width+j)*2+1];
					float pos_x2 = coord_data[(i*width+j+width)*2+0];
					float pos_y2 = coord_data[(i*width+j+width)*2+1];
					cvLine(render_img,cvPoint(pos_x1*show_scale+show_shift_x,pos_y1*show_scale+show_shift_y),
						cvPoint(pos_x2*show_scale+show_shift_x,pos_y2*show_scale+show_shift_y),
						color_line,thick_point);
				}
			}
		}
		return true;
	}

	static bool render_with_texture(ZQ::ZQ_DImage<float>& render_img, const ZQ::ZQ_DImage<float>& map_img, const ZQ::ZQ_DImage<float>& tex_img, const ZQ::ZQ_DImage<bool>& obj_mask, 
		const int show_scale, float show_shift_x, float show_shift_y, bool xloop, int border_size)
	{
		int width = map_img.width();
		int height = map_img.height();
		if(width-1-2*border_size <= 0 || height-1-2*border_size <= 0)
			return false;

		if(!map_img.matchDimension(width,height,2))
			return false;
		int nChannels = tex_img.nchannels();
		
		if(!obj_mask.matchDimension(width,height,1))
			return false;
		if(render_img.nchannels() != nChannels)
			return false;

		const float*& coord_data = map_img.data();
		const float*& tex_data = tex_img.data();
		const bool*& mask_data = obj_mask.data();


		ZQ::ZQ_TextureSampler<float> tex_sampler;
		tex_sampler.BindImage(tex_img,xloop);

		
		if(xloop)
		{
			for(int i = border_size;i < height-1-border_size;i++)
			{
				for(int j = 0;j <= width-1;j++)
				{
					int offset_00 = i*width+j;
					int offset_01 = i*width+(j+1)%width;
					int offset_10 = (i+1)*width+j;
					int offset_11 = (i+1)*width+(j+1)%width;
					if(mask_data[offset_00] || mask_data[offset_01] || mask_data[offset_10] || mask_data[offset_11])
					{
						continue;
					}

					ZQ::ZQ_Vec2D pt00(coord_data[offset_00*2+0]*show_scale+show_shift_x,coord_data[offset_00*2+1]*show_scale+show_shift_y);
					ZQ::ZQ_Vec2D pt01(coord_data[offset_01*2+0]*show_scale+show_shift_x,coord_data[offset_01*2+1]*show_scale+show_shift_y);
					ZQ::ZQ_Vec2D pt10(coord_data[offset_10*2+0]*show_scale+show_shift_x,coord_data[offset_10*2+1]*show_scale+show_shift_y);
					ZQ::ZQ_Vec2D pt11(coord_data[offset_11*2+0]*show_scale+show_shift_x,coord_data[offset_11*2+1]*show_scale+show_shift_y);

					ZQ::ZQ_Vec2D texCoord00((float)j/width,(float)(i-border_size)/(height-1-2*border_size));
					ZQ::ZQ_Vec2D texCoord01((float)(j+1)/width,(float)(i-border_size)/(height-1-2*border_size));
					ZQ::ZQ_Vec2D texCoord10((float)j/width,(float)(i+1-border_size)/(height-1-2*border_size));
					ZQ::ZQ_Vec2D texCoord11((float)(j+1)/width,(float)(i+1-border_size)/(height-1-2*border_size));

					_render_triangle(render_img,tex_sampler,pt00,texCoord00,pt01,texCoord01,pt10,texCoord10);
					_render_triangle(render_img,tex_sampler,pt01,texCoord01,pt10,texCoord10,pt11,texCoord11);
				}
			}
		}
		else
		{
			for(int i = border_size;i < height-1-border_size;i++)
			{
				for(int j = border_size;j < width-1-border_size;j++)
				{
					int offset_00 = i*width+j;
					int offset_01 = i*width+(j+1)%width;
					int offset_10 = (i+1)*width+j;
					int offset_11 = (i+1)*width+(j+1)%width;
					if(mask_data[offset_00] || mask_data[offset_01] || mask_data[offset_10] || mask_data[offset_11])
					{
						continue;
					}

					ZQ::ZQ_Vec2D pt00(coord_data[offset_00*2+0]*show_scale+show_shift_x,coord_data[offset_00*2+1]*show_scale+show_shift_y);
					ZQ::ZQ_Vec2D pt01(coord_data[offset_01*2+0]*show_scale+show_shift_x,coord_data[offset_01*2+1]*show_scale+show_shift_y);
					ZQ::ZQ_Vec2D pt10(coord_data[offset_10*2+0]*show_scale+show_shift_x,coord_data[offset_10*2+1]*show_scale+show_shift_y);
					ZQ::ZQ_Vec2D pt11(coord_data[offset_11*2+0]*show_scale+show_shift_x,coord_data[offset_11*2+1]*show_scale+show_shift_y);

					ZQ::ZQ_Vec2D texCoord00((float)(j-border_size)/(width-1-2*border_size),(float)(i-border_size)/(height-1-2*border_size));
					ZQ::ZQ_Vec2D texCoord01((float)(j+1-border_size)/(width-1-2*border_size),(float)(i-border_size)/(height-1-2*border_size));
					ZQ::ZQ_Vec2D texCoord10((float)(j-border_size)/(width-1-2*border_size),(float)(i+1-border_size)/(height-1-2*border_size));
					ZQ::ZQ_Vec2D texCoord11((float)(j+1-border_size)/(width-1-2*border_size),(float)(i+1-border_size)/(height-1-2*border_size));

					_render_triangle(render_img,tex_sampler,pt00,texCoord00,pt01,texCoord01,pt10,texCoord10);
					_render_triangle(render_img,tex_sampler,pt01,texCoord01,pt10,texCoord10,pt11,texCoord11);
				}
			}
		}

		
		
		return true;
	}

private:

	static void _render_triangle(ZQ::ZQ_DImage<float>& render_img, const ZQ::ZQ_TextureSampler<float>& tex_sampler, const ZQ::ZQ_Vec2D& pt1, const ZQ::ZQ_Vec2D& texCoord1,
		const ZQ::ZQ_Vec2D& pt2, const ZQ::ZQ_Vec2D& texCoord2,const ZQ::ZQ_Vec2D& pt3, const ZQ::ZQ_Vec2D& texCoord3)
	{
		std::vector<ZQ::ZQ_Vec2D> poly,pixels;
		poly.push_back(pt1);
		poly.push_back(pt2);
		poly.push_back(pt3);
		ZQ::ZQ_ScanLinePolygonFill::ScanLinePolygonFill(poly,pixels);
		if(pixels.size() == 0)
			return;

		float*& render_data = render_img.data();
		int nChannels = render_img.nchannels();
		int width = render_img.width();
		int height = render_img.height();
		for(int i = 0;i < pixels.size();i++)
		{
			ZQ::ZQ_Vec2D cur_pt = pixels[i];
			float area3 = _area_of_triangle(pt1,pt2,cur_pt);
			float area1 = _area_of_triangle(pt2,pt3,cur_pt);
			float area2 = _area_of_triangle(pt3,pt1,cur_pt);
			int cur_x = cur_pt.x;
			int cur_y = cur_pt.y;
			if(cur_x < 0 || cur_x >= width || cur_y < 0 || cur_y >= height)
				continue;
			float sum_area = area1+area2+area3;
			if(sum_area == 0)
				continue;
			float w1 = area1/sum_area;
			float w2 = area2/sum_area;
			float w3 = area3/sum_area;
			ZQ::ZQ_Vec2D cur_texCoord = texCoord1*w1 + texCoord2*w2 + texCoord3*w3;
			float tex_coord_x = cur_texCoord.x;
			float tex_coord_y = cur_texCoord.y;
			tex_sampler.Sample_NormalizedCoord(tex_coord_x,tex_coord_y,render_data+(cur_y*width+cur_x)*nChannels,false);
		}
	}

	static float _area_of_triangle(const ZQ::ZQ_Vec2D& pt1,const ZQ::ZQ_Vec2D& pt2,const ZQ::ZQ_Vec2D& pt3)
	{
		float x1 = pt2.x - pt1.x;
		float y1 = pt2.y - pt1.y;
		float x2 = pt3.x - pt1.x;
		float y2 = pt3.y - pt1.y;
		return fabs(x1*y2-x2*y1);
	}

};

#endif