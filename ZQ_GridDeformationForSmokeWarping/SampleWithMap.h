#ifndef _SAMPLE_WITH_MAP_H_
#define _SAMPLE_WITH_MAP_H_
#pragma once

#include "ZQ_DoubleImage.h"
#include "ZQ_TextureSampler.h"
#include "ZQ_ImageIO.h"

class SampleWithMapOptions
{
public:
	SampleWithMapOptions()
	{
		Reset();
	}
	~SampleWithMapOptions()
	{

	}

	enum CONST_VAL{FILE_LEN = 200};
	bool has_map_file;
	char map_file[FILE_LEN];
	bool has_mask_file;
	char mask_file[FILE_LEN];
	bool has_tex_file;
	char tex_file[FILE_LEN];
	bool has_show_file;
	char show_file[FILE_LEN];
	bool xloop;
	bool normalized;
	bool cubic;
	
	void Reset()
	{
		has_map_file = false;
		map_file[0] = '\0';
		has_mask_file = false;
		mask_file[0] = '\0';
		has_tex_file = false;
		tex_file[0] = '\0';
		xloop = false;
		normalized = false;
		cubic = false;
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
			else if(_strcmpi(argv[k],"mask_file") == 0)
			{
				k++;
				if(k >= argc)
				{
					printf("the value of mask_file ?\n");
					return false;
				}
				strcpy(mask_file,argv[k]);
				has_mask_file = true;
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
			else if(_strcmpi(argv[k],"normalized") == 0)
			{
				normalized = true;
			}
			else if(_strcmpi(argv[k],"cubic") == 0)
			{
				cubic = true;
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

class SampleWithMap
{
public:
	static bool go(const SampleWithMapOptions& opt)
	{
		if(!opt.has_map_file)
		{
			printf("error: need map_file!\n");
			return false;
		}
		if(!opt.has_tex_file)
		{
			printf("error: need tex_file\n");
			return false;
		}
		if(!opt.has_show_file)
		{
			printf("error: need show_file!\n");
			return false;
		}


		ZQ::ZQ_DImage<float> map_img;
		if(!map_img.loadImage(opt.map_file))
		{
			printf("failed to load %s\n",opt.map_file);
			return false;
		}

		int map_width = map_img.width();
		int map_height = map_img.height();

		ZQ::ZQ_DImage<bool> obj_mask(map_width,map_height);
		if(opt.has_mask_file)
		{
			ZQ::ZQ_DImage<float> mask;
			if(!ZQ::ZQ_ImageIO::loadImage(mask,opt.mask_file,0))
			{
				printf("failed to load %s\n",opt.mask_file);
				return false;
			}
			if(mask.width() != map_img.width() || mask.height() != map_img.height())
			{
				printf("mask and map don't match in dimension\n");
				return false;
			}
			for(int i = 0;i < obj_mask.npixels();i++)
			{
				obj_mask.data()[i] = mask.data()[i] < 0.5;
			}
		}


		ZQ::ZQ_DImage<float> tex_img;
		if(!ZQ::ZQ_ImageIO::loadImage(tex_img,opt.tex_file,1))
		{
			printf("failed to load %s\n",opt.tex_file);
			return false;
		}

		int nChannels = tex_img.nchannels();
		ZQ::ZQ_DImage<float> show_img(map_width,map_height,nChannels);

		if(!_sample_with_map(show_img,map_img,tex_img,obj_mask,opt.normalized,opt.cubic,opt.xloop))
			return false;
		if(!ZQ::ZQ_ImageIO::saveImage(show_img,opt.show_file))
		{
			printf("failed to save %s\n",opt.show_file);
			return false;
		}
		return true;
	}

	static bool _sample_with_map(ZQ::ZQ_DImage<float>& show_img, const ZQ::ZQ_DImage<float>& map_img, const ZQ::ZQ_DImage<float>& tex_img, const ZQ::ZQ_DImage<bool>& obj_mask, bool normalized, bool cubic, bool xloop)
	{
		int map_width = map_img.width();
		int map_height = map_img.height();
		if(2 != map_img.nchannels())
			return false;
		if(!obj_mask.matchDimension(map_width,map_height,1))
			return false;

		int nChannels = tex_img.nchannels();
		if(!show_img.matchDimension(map_width,map_height,nChannels))
			show_img.allocate(map_width,map_height,nChannels);

		ZQ::ZQ_TextureSampler<float> sampler;
		sampler.BindImage(tex_img,xloop);

		const bool*& obj_mask_data = obj_mask.data();
		const float*& map_img_data = map_img.data();
		float*& show_img_data = show_img.data();

		for(int i = 0;i < map_height;i++)
		{
			for(int j = 0;j < map_width;j++)
			{
				int offset = i*map_width+j;
				if(obj_mask_data[offset])
					continue;

				if(normalized)
					sampler.Sample_NormalizedCoord(map_img_data[offset*2+0],map_img_data[offset*2+1],show_img_data+offset*nChannels,cubic);
				else
					sampler.Sample(map_img_data[offset*2+0],map_img_data[offset*2+1],show_img_data+offset*nChannels,cubic);
			}
		}
		
		return true;
	}
};

#endif