#ifndef _RENDER_GRID_WITH_TEXTURE_3D_H_
#define _RENDER_GRID_WITH_TEXTURE_3D_H_
#pragma once

#include "ZQ_DoubleImage3D.h"
#include "ZQ_ScanLinePolygonFill.h"
#include "ZQ_TextureSampler3D.h"
#include "ZQ_Vec3D.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

class RenderGridWithTexture3DOptions
{
public:
	RenderGridWithTexture3DOptions()
	{
		Reset();
	}
	~RenderGridWithTexture3DOptions()
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
	int BoundingBoxmin[3];
	int BoundingBoxmax[3];
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
		BoundingBoxmin[0] = BoundingBoxmin[1] = BoundingBoxmin[2] = 0;
		BoundingBoxmax[0] = BoundingBoxmax[1] = BoundingBoxmax[2] = 0;
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

class RenderGridWithTexture3D
{
public:
	static bool go(const RenderGridWithTexture3DOptions& opt)
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
		if(!opt.has_tex_file)
		{
			printf("error: need tex_file!\n");
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
		float show_shift_x = (0.5-opt.BoundingBoxmin[0])*show_scale;
		float show_shift_y = (0.5-opt.BoundingBoxmin[1])*show_scale;
		float show_shift_z = (0.5-opt.BoundingBoxmin[2])*show_scale;
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
		ZQ::ZQ_DImage3D<bool> obj_mask(map_img.width(),map_img.height(),map_img.depth());
		if(opt.has_map_mask_file)
		{
			ZQ::ZQ_DImage3D<float> mask;
			if(!mask.loadImage(opt.map_mask_file))
			{
				printf("failed to load %s\n",opt.map_mask_file);
				return false;
			}
			if(mask.width() != map_img.width() || mask.height() != map_img.height() || mask.depth() != map_img.depth())
			{
				printf("mask and map don't match in dimension\n");
				return false;
			}
			for(int i = 0;i < obj_mask.npixels();i++)
			{
				obj_mask.data()[i] = mask.data()[i] > 0.5;
			}
		}


		ZQ::ZQ_DImage3D<float> tex_img;
		if(!tex_img.loadImage(opt.tex_file))
		{
			printf("failed to load %s\n",opt.tex_file);
			return false;
		}
		int nChannels = tex_img.nchannels();
		ZQ::ZQ_DImage3D<float> render_img(back_width,back_height,back_depth,nChannels);
		if(!render_with_texture(render_img,map_img,tex_img,obj_mask,show_scale,show_shift_x,show_shift_y,show_shift_z,opt.xloop,opt.border_size))
			return false;
		if(!render_img.saveImage(opt.show_file))
		{
			printf("failed to save %s\n",opt.show_file);
			return false;
		}

		return true;
	}

	static bool render_with_texture(ZQ::ZQ_DImage3D<float>& render_img, const ZQ::ZQ_DImage3D<float>& map_img, const ZQ::ZQ_DImage3D<float>& tex_img, const ZQ::ZQ_DImage3D<bool>& obj_mask, 
		const int show_scale, float show_shift_x, float show_shift_y, float show_shift_z, bool xloop, int border_size)
	{
		int width = map_img.width();
		int height = map_img.height();
		int depth = map_img.depth();
		if(width-1-2*border_size <= 0 || height-1-2*border_size <= 0 || depth-1-2*border_size <= 0 || width == 0 || height == 0 || depth == 0)
			return false;

		if(!map_img.matchDimension(width,height,depth,3))
			return false;
		int nChannels = tex_img.nchannels();

		if(!obj_mask.matchDimension(width,height,depth,1))
			return false;
		if(render_img.nchannels() != nChannels)
			return false;

		const float*& coord_data = map_img.data();
		const float*& tex_data = tex_img.data();
		const bool*& mask_data = obj_mask.data();


		ZQ::ZQ_TextureSampler3D<float> tex_sampler;
		tex_sampler.BindImage(tex_img,xloop);


		if(xloop)
		{
			float xlen = width;
			float ylen = height-1-2*border_size;
			float zlen = depth-1-2*border_size;
			for(int k = border_size;k < depth-1-border_size;k++)
			{
				for(int j = border_size;j < height-1-border_size;j++)
				{
					for(int i = 0;i <= width-1;i++)
					{
						int offset_000 = k*height*width+j*width+i;
						int offset_001 = k*height*width+j*width+(i+1)%width;
						int offset_010 = k*height*width+(j+1)*width+i;
						int offset_011 = k*height*width+(j+1)*width+(i+1)%width;
						int offset_100 = (k+1)*height*width+j*width+i;
						int offset_101 = (k+1)*height*width+j*width+(i+1)%width;
						int offset_110 = (k+1)*height*width+(j+1)*width+i;
						int offset_111 = (k+1)*height*width+(j+1)*width+(i+1)%width;
						if(mask_data[offset_000] || mask_data[offset_001] || mask_data[offset_010] || mask_data[offset_011]
						|| mask_data[offset_100] || mask_data[offset_101] || mask_data[offset_110] || mask_data[offset_111])
							continue;

						ZQ::ZQ_Vec3D pt000(coord_data[offset_000*3+0]*show_scale+show_shift_x,coord_data[offset_000*3+1]*show_scale+show_shift_y,coord_data[offset_000*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt001(coord_data[offset_001*3+0]*show_scale+show_shift_x,coord_data[offset_001*3+1]*show_scale+show_shift_y,coord_data[offset_001*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt010(coord_data[offset_010*3+0]*show_scale+show_shift_x,coord_data[offset_010*3+1]*show_scale+show_shift_y,coord_data[offset_010*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt011(coord_data[offset_011*3+0]*show_scale+show_shift_x,coord_data[offset_011*3+1]*show_scale+show_shift_y,coord_data[offset_011*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt100(coord_data[offset_100*3+0]*show_scale+show_shift_x,coord_data[offset_100*3+1]*show_scale+show_shift_y,coord_data[offset_100*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt101(coord_data[offset_101*3+0]*show_scale+show_shift_x,coord_data[offset_101*3+1]*show_scale+show_shift_y,coord_data[offset_101*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt110(coord_data[offset_110*3+0]*show_scale+show_shift_x,coord_data[offset_110*3+1]*show_scale+show_shift_y,coord_data[offset_110*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt111(coord_data[offset_111*3+0]*show_scale+show_shift_x,coord_data[offset_111*3+1]*show_scale+show_shift_y,coord_data[offset_111*3+2]*show_scale+show_shift_z);

						ZQ::ZQ_Vec3D texCoord000((float)i/xlen,(float)(j-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord001((float)(i+1)/xlen,(float)(j-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord010((float)i/xlen,(float)(j+1-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord011((float)(i+1)/xlen,(float)(j+1-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord100((float)i/xlen,(float)(j-border_size)/ylen,(float)(k+1-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord101((float)(i+1)/xlen,(float)(j-border_size)/ylen,(float)(k+1-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord110((float)i/xlen,(float)(j+1-border_size)/ylen,(float)(k+1-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord111((float)(i+1)/xlen,(float)(j+1-border_size)/ylen,(float)(k+1-border_size)/zlen);

						_render_tetrahedron(render_img,tex_sampler,pt000,texCoord000,pt001,texCoord001,pt010,texCoord010,pt100,texCoord100);
						_render_tetrahedron(render_img,tex_sampler,pt100,texCoord100,pt101,texCoord101,pt110,texCoord110,pt010,texCoord010);
						_render_tetrahedron(render_img,tex_sampler,pt100,texCoord100,pt101,texCoord101,pt001,texCoord001,pt010,texCoord010);

						_render_tetrahedron(render_img,tex_sampler,pt011,texCoord011,pt001,texCoord001,pt010,texCoord010,pt111,texCoord111);
						_render_tetrahedron(render_img,tex_sampler,pt111,texCoord111,pt101,texCoord101,pt110,texCoord110,pt010,texCoord010);
						_render_tetrahedron(render_img,tex_sampler,pt111,texCoord111,pt101,texCoord101,pt001,texCoord001,pt010,texCoord010);

					}
				}
			}
		}
		else
		{
			float xlen = width-1-2*border_size;
			float ylen = height-1-2*border_size;
			float zlen = depth-1-2*border_size;
			for(int k = border_size;k < depth-1-border_size;k++)
			{
				for(int j = border_size;j < height-1-border_size;j++)
				{
					for(int i = border_size;i < width-1-border_size;i++)
					{
						int offset_000 = k*height*width+j*width+i;
						int offset_001 = k*height*width+j*width+i+1;
						int offset_010 = k*height*width+(j+1)*width+i;
						int offset_011 = k*height*width+(j+1)*width+i+1;
						int offset_100 = (k+1)*height*width+j*width+i;
						int offset_101 = (k+1)*height*width+j*width+i+1;
						int offset_110 = (k+1)*height*width+(j+1)*width+i;
						int offset_111 = (k+1)*height*width+(j+1)*width+i+1;

						if(mask_data[offset_000] || mask_data[offset_001] || mask_data[offset_010] || mask_data[offset_011]
						|| mask_data[offset_100] || mask_data[offset_101] || mask_data[offset_110] || mask_data[offset_111])
							continue;

						ZQ::ZQ_Vec3D pt000(coord_data[offset_000*3+0]*show_scale+show_shift_x,coord_data[offset_000*3+1]*show_scale+show_shift_y,coord_data[offset_000*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt001(coord_data[offset_001*3+0]*show_scale+show_shift_x,coord_data[offset_001*3+1]*show_scale+show_shift_y,coord_data[offset_001*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt010(coord_data[offset_010*3+0]*show_scale+show_shift_x,coord_data[offset_010*3+1]*show_scale+show_shift_y,coord_data[offset_010*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt011(coord_data[offset_011*3+0]*show_scale+show_shift_x,coord_data[offset_011*3+1]*show_scale+show_shift_y,coord_data[offset_011*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt100(coord_data[offset_100*3+0]*show_scale+show_shift_x,coord_data[offset_100*3+1]*show_scale+show_shift_y,coord_data[offset_100*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt101(coord_data[offset_101*3+0]*show_scale+show_shift_x,coord_data[offset_101*3+1]*show_scale+show_shift_y,coord_data[offset_101*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt110(coord_data[offset_110*3+0]*show_scale+show_shift_x,coord_data[offset_110*3+1]*show_scale+show_shift_y,coord_data[offset_110*3+2]*show_scale+show_shift_z);
						ZQ::ZQ_Vec3D pt111(coord_data[offset_111*3+0]*show_scale+show_shift_x,coord_data[offset_111*3+1]*show_scale+show_shift_y,coord_data[offset_111*3+2]*show_scale+show_shift_z);

						ZQ::ZQ_Vec3D texCoord000((float)(i-border_size)/xlen,(float)(j-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord001((float)(i+1-border_size)/xlen,(float)(j-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord010((float)(i-border_size)/xlen,(float)(j+1-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord011((float)(i+1-border_size)/xlen,(float)(j+1-border_size)/ylen,(float)(k-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord100((float)(i-border_size)/xlen,(float)(j-border_size)/ylen,(float)(k+1-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord101((float)(i+1-border_size)/xlen,(float)(j-border_size)/ylen,(float)(k+1-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord110((float)(i-border_size)/xlen,(float)(j+1-border_size)/ylen,(float)(k+1-border_size)/zlen);
						ZQ::ZQ_Vec3D texCoord111((float)(i+1-border_size)/xlen,(float)(j+1-border_size)/ylen,(float)(k+1-border_size)/zlen);


						_render_tetrahedron(render_img,tex_sampler,pt000,texCoord000,pt001,texCoord001,pt010,texCoord010,pt100,texCoord100);
						_render_tetrahedron(render_img,tex_sampler,pt100,texCoord100,pt101,texCoord101,pt110,texCoord110,pt010,texCoord010);
						_render_tetrahedron(render_img,tex_sampler,pt100,texCoord100,pt101,texCoord101,pt001,texCoord001,pt010,texCoord010);

						_render_tetrahedron(render_img,tex_sampler,pt011,texCoord011,pt001,texCoord001,pt010,texCoord010,pt111,texCoord111);
						_render_tetrahedron(render_img,tex_sampler,pt111,texCoord111,pt101,texCoord101,pt110,texCoord110,pt010,texCoord010);
						_render_tetrahedron(render_img,tex_sampler,pt111,texCoord111,pt101,texCoord101,pt001,texCoord001,pt010,texCoord010);

					}
				}
			}
		}


		return true;
	}

private:

	static void _render_tetrahedron(ZQ::ZQ_DImage3D<float>& render_img, const ZQ::ZQ_TextureSampler3D<float>& tex_sampler, const ZQ::ZQ_Vec3D& pt1, const ZQ::ZQ_Vec3D& texCoord1,
		const ZQ::ZQ_Vec3D& pt2, const ZQ::ZQ_Vec3D& texCoord2,const ZQ::ZQ_Vec3D& pt3, const ZQ::ZQ_Vec3D& texCoord3,const ZQ::ZQ_Vec3D& pt4, const ZQ::ZQ_Vec3D& texCoord4)
	{
		float*& render_img_data = render_img.data();
		int render_width = render_img.width();
		int render_height = render_img.height();
		int render_depth = render_img.depth();
		int nChannels = render_img.nchannels();
		ZQ::ZQ_Vec3D in_pts[4] = {pt1,pt2,pt3,pt4};
		ZQ::ZQ_Vec3D in_coords[4] = {texCoord1,texCoord2,texCoord3,texCoord4};
		int idx[4];
		float z_coord[4] = {pt1.z,pt2.z,pt3.z,pt4.z};
		_sort4_inc(z_coord,idx);
		ZQ::ZQ_Vec3D order_pts[4] = {in_pts[idx[0]],in_pts[idx[1]],in_pts[idx[2]],in_pts[idx[3]]};
		ZQ::ZQ_Vec3D order_coords[4] = {in_coords[idx[0]],in_coords[idx[1]],in_coords[idx[2]],in_coords[idx[3]]};

		ZQ::ZQ_Vec3D dir_0_1 = order_pts[1]-order_pts[0];
		ZQ::ZQ_Vec3D dir_0_2 = order_pts[2]-order_pts[0];
		ZQ::ZQ_Vec3D dir_0_3 = order_pts[3]-order_pts[0];
		ZQ::ZQ_Vec3D dir_1_2 = order_pts[2]-order_pts[1];
		ZQ::ZQ_Vec3D dir_1_3 = order_pts[3]-order_pts[1];
		ZQ::ZQ_Vec3D dir_2_3 = order_pts[3]-order_pts[2];
		
		//if dir_0_1.z == 0, it will not loop
		for(int z = __max(0,ceil(order_pts[0].z)); z <= __min(render_depth-1,ceil(order_pts[1].z-1)); z++)
		{
			float z_shift = z - order_pts[0].z;
			float w1 = z_shift/dir_0_1.z;
			float w2 = z_shift/dir_0_2.z;
			float w3 = z_shift/dir_0_3.z;
			
			ZQ::ZQ_Vec3D pt_z_0_1 = order_pts[0] + dir_0_1 * w1;
			ZQ::ZQ_Vec3D pt_z_0_2 = order_pts[0] + dir_0_2 * w2;
			ZQ::ZQ_Vec3D pt_z_0_3 = order_pts[0] + dir_0_3 * w3;

			ZQ::ZQ_Vec3D coord_z_0_1 = order_coords[0]*(1-w1) + order_coords[1]*w1;
			ZQ::ZQ_Vec3D coord_z_0_2 = order_coords[0]*(1-w2) + order_coords[2]*w2;
			ZQ::ZQ_Vec3D coord_z_0_3 = order_coords[0]*(1-w3) + order_coords[3]*w3;

			std::vector<ZQ::ZQ_Vec2D> pts2d,pixels;
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_0_1.x,pt_z_0_1.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_0_2.x,pt_z_0_2.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_0_3.x,pt_z_0_3.y));

			ZQ::ZQ_ScanLinePolygonFill::ScanLinePolygonFill(pts2d,pixels);

			float area_total = _area_of_triangle(pts2d[0],pts2d[1],pts2d[2]);
			for(int pp = 0;pp < pixels.size();pp++)
			{
				if(pixels[pp].x < 0 || pixels[pp].x >= render_width || pixels[pp].y < 0 || pixels[pp].y >= render_height)
					continue;

				float area0 = _area_of_triangle(pts2d[1],pts2d[2],pixels[pp]);
				float area1 = _area_of_triangle(pts2d[0],pts2d[2],pixels[pp]);
				float area2 = _area_of_triangle(pts2d[0],pts2d[1],pixels[pp]);
				ZQ::ZQ_Vec3D cur_tex_coord = coord_z_0_1*(area0/area_total) + coord_z_0_2*(area1/area_total) + coord_z_0_3*(area2/area_total);

				int xx = pixels[pp].x;
				int yy = pixels[pp].y;
				tex_sampler.Sample_NormalizedCoord(cur_tex_coord.x,cur_tex_coord.y,cur_tex_coord.z,render_img_data+(z*render_height*render_width+yy*render_width+xx)*nChannels,false);
			}
		}

		//if dir_1_2.z == 0, it will not loop
		for(int z = __max(0,ceil(order_pts[1].z)); z <= __min(render_depth-1,ceil(order_pts[2].z-1)); z++)
		{
			float z_shift_0 = z - order_pts[0].z;
			float z_shift_1 = z - order_pts[1].z;
			float w_0_2 = z_shift_0/dir_0_2.z;
			float w_0_3 = z_shift_0/dir_0_3.z;
			float w_1_2 = z_shift_1/dir_1_2.z;
			float w_1_3 = z_shift_1/dir_1_3.z;

			ZQ::ZQ_Vec3D pt_z_0_2 = order_pts[0] + dir_0_2 * w_0_2;
			ZQ::ZQ_Vec3D pt_z_0_3 = order_pts[0] + dir_0_3 * w_0_3;
			ZQ::ZQ_Vec3D pt_z_1_2 = order_pts[1] + dir_1_2 * w_1_2;
			ZQ::ZQ_Vec3D pt_z_1_3 = order_pts[1] + dir_1_3 * w_1_3;

			ZQ::ZQ_Vec3D coord_z_0_2 = order_coords[0]*(1-w_0_2) + order_coords[2]*w_0_2;
			ZQ::ZQ_Vec3D coord_z_0_3 = order_coords[0]*(1-w_0_3) + order_coords[3]*w_0_3;
			ZQ::ZQ_Vec3D coord_z_1_2 = order_coords[1]*(1-w_1_2) + order_coords[2]*w_1_2;
			ZQ::ZQ_Vec3D coord_z_1_3 = order_coords[1]*(1-w_1_3) + order_coords[3]*w_1_3;

			std::vector<ZQ::ZQ_Vec2D> pts2d,pixels;
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_0_2.x,pt_z_0_2.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_1_3.x,pt_z_1_3.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_0_3.x,pt_z_0_3.y));

			ZQ::ZQ_ScanLinePolygonFill::ScanLinePolygonFill(pts2d,pixels);

			float area_total = _area_of_triangle(pts2d[0],pts2d[1],pts2d[2]);
			for(int pp = 0;pp < pixels.size();pp++)
			{
				if(pixels[pp].x < 0 || pixels[pp].x >= render_width || pixels[pp].y < 0 || pixels[pp].y >= render_height)
					continue;

				float area0 = _area_of_triangle(pts2d[1],pts2d[2],pixels[pp]);
				float area1 = _area_of_triangle(pts2d[0],pts2d[2],pixels[pp]);
				float area2 = _area_of_triangle(pts2d[0],pts2d[1],pixels[pp]);
				ZQ::ZQ_Vec3D cur_tex_coord = coord_z_0_2*(area0/area_total) + coord_z_1_3*(area1/area_total) + coord_z_0_3*(area2/area_total);

				int xx = pixels[pp].x;
				int yy = pixels[pp].y;
				tex_sampler.Sample_NormalizedCoord(cur_tex_coord.x,cur_tex_coord.y,cur_tex_coord.z,render_img_data+(z*render_height*render_width+yy*render_width+xx)*nChannels,false);
			}

			pts2d.clear(); pixels.clear();
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_0_2.x,pt_z_0_2.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_1_3.x,pt_z_1_3.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_1_2.x,pt_z_1_2.y));

			ZQ::ZQ_ScanLinePolygonFill::ScanLinePolygonFill(pts2d,pixels);

			area_total = _area_of_triangle(pts2d[0],pts2d[1],pts2d[2]);
			for(int pp = 0;pp < pixels.size();pp++)
			{
				if(pixels[pp].x < 0 || pixels[pp].x >= render_width || pixels[pp].y < 0 || pixels[pp].y >= render_height)
					continue;

				float area0 = _area_of_triangle(pts2d[1],pts2d[2],pixels[pp]);
				float area1 = _area_of_triangle(pts2d[0],pts2d[2],pixels[pp]);
				float area2 = _area_of_triangle(pts2d[0],pts2d[1],pixels[pp]);
				ZQ::ZQ_Vec3D cur_tex_coord = coord_z_0_2*(area0/area_total) + coord_z_1_3*(area1/area_total) + coord_z_1_2*(area2/area_total);

				int xx = pixels[pp].x;
				int yy = pixels[pp].y;
				tex_sampler.Sample_NormalizedCoord(cur_tex_coord.x,cur_tex_coord.y,cur_tex_coord.z,render_img_data+(z*render_height*render_width+yy*render_width+xx)*nChannels,false);
			}
		}

		//if dir_2_3.z == 0, it will not loop
		for(int z = __max(0,ceil(order_pts[2].z)); z <= __min(render_depth-1,ceil(order_pts[3].z-1)); z++)
		{
			float z_shift = order_pts[3].z - z;
			float w0 = z_shift/dir_0_3.z;
			float w1 = z_shift/dir_1_3.z;
			float w2 = z_shift/dir_2_3.z;

			ZQ::ZQ_Vec3D pt_z_0_3 = order_pts[3] - dir_0_3 * w0;
			ZQ::ZQ_Vec3D pt_z_1_3 = order_pts[3] - dir_1_3 * w1;
			ZQ::ZQ_Vec3D pt_z_2_3 = order_pts[3] - dir_2_3 * w2;

			ZQ::ZQ_Vec3D coord_z_0_3 = order_coords[3]*(1-w0) + order_coords[0]*w0;
			ZQ::ZQ_Vec3D coord_z_1_3 = order_coords[3]*(1-w1) + order_coords[1]*w1;
			ZQ::ZQ_Vec3D coord_z_2_3 = order_coords[3]*(1-w2) + order_coords[2]*w2;

			std::vector<ZQ::ZQ_Vec2D> pts2d,pixels;
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_0_3.x,pt_z_0_3.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_1_3.x,pt_z_1_3.y));
			pts2d.push_back(ZQ::ZQ_Vec2D(pt_z_2_3.x,pt_z_2_3.y));

			ZQ::ZQ_ScanLinePolygonFill::ScanLinePolygonFill(pts2d,pixels);

			float area_total = _area_of_triangle(pts2d[0],pts2d[1],pts2d[2]);
			for(int pp = 0;pp < pixels.size();pp++)
			{
				if(pixels[pp].x < 0 || pixels[pp].x >= render_width || pixels[pp].y < 0 || pixels[pp].y >= render_height)
					continue;
				float area0 = _area_of_triangle(pts2d[1],pts2d[2],pixels[pp]);
				float area1 = _area_of_triangle(pts2d[0],pts2d[2],pixels[pp]);
				float area2 = _area_of_triangle(pts2d[0],pts2d[1],pixels[pp]);
				ZQ::ZQ_Vec3D cur_tex_coord = coord_z_0_3*(area0/area_total) + coord_z_1_3*(area1/area_total) + coord_z_2_3*(area2/area_total);

				int xx = pixels[pp].x;
				int yy = pixels[pp].y;
				tex_sampler.Sample_NormalizedCoord(cur_tex_coord.x,cur_tex_coord.y,cur_tex_coord.z,render_img_data+(z*render_height*render_width+yy*render_width+xx)*nChannels,false);
			}
		}

	}

	static void _sort4_inc(const float in_pts[4], int idx[4])
	{
		float val[4] = {in_pts[0],in_pts[1],in_pts[2],in_pts[3]};
		idx[0] = 0; idx[1] = 1; idx[2] = 2; idx[3] = 3;
		for(int pass = 0;pass < 3;pass++)
		{
			for(int i = 0;i < 3;i++)
			{
				if(val[i] > val[i+1])
				{
					float tmp = val[i];
					val[i] = val[i+1];
					val[i+1] = tmp;
					int ii = idx[i];
					idx[i] = idx[i+1];
					idx[i+1] = ii;
				}
			}
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