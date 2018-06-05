#include "ZQ_DoubleImage.h"
#include "ZQ_DoubleImage3D.h"
#include "ZQ_ImageIO.h"
#include "ZQ_PoissonEditing.h"
#include <vector>
#include <stdio.h>
#include <string.h>
#include "ZQ_CUDA_PoissonEditing3D.h"


using namespace ZQ;
typedef ZQ_DImage<float> DImage;
typedef ZQ_DImage3D<float> DImage3D;

class ROI
{
public:
	ROI(){off_x = off_y = width = height = 0;}
	~ROI(){}
	int off_x,off_y;
	int width,height;
	ZQ_DImage<bool> mask;
};

void _find_roi(const ZQ_DImage<bool>& mask, std::vector<ROI>& rois);
bool _solve_sub_problem(const std::vector<DImage>& input_images,const std::vector<DImage>& copyin_images,DImage3D& output3D, const ROI& roi);

bool PoissonEditing3DCuda(const ZQ_DImage3D<bool>& mask, const DImage3D& input, const DImage3D& copyin, DImage3D& output, const int nIteration);

void main(const int argc, const char** argv)
{
	if(argc != 8 && argc != 9)
	{
		printf(" args error\n");
		return;
	}
	int num = atoi(argv[1]);
	const char* input_fold = argv[2];
	const char* input_prefix = "";
	const char* input_suffix = argv[3];
	const char* copyin_fold = argv[4];
	const char* copyin_prefix = "";
	const char* copyin_suffix = argv[5];
	const char* mask_file = argv[6];
	const char* output_fold = argv[7];
	int base_id = 0;
	if(argc == 9)
		base_id = atoi(argv[8]);

	std::vector<DImage> input_images;
	std::vector<DImage> copyin_images;
	DImage mask;

	if(!ZQ_ImageIO::loadImage(mask,mask_file,0))
	{
		printf("failed to load %s\n",mask_file);
		return;
	}

	int width = mask.width();
	int height = mask.height();
	ZQ_DImage<bool> mask_flag(width,height);
	bool*& mask_flag_data = mask_flag.data();
	float*& mask_data = mask.data();
	for(int pp = 0;pp < width*height;pp++)
		mask_flag_data[pp] = mask_data[pp] >= 0.5;

	std::vector<ROI> rois;
	_find_roi(mask_flag,rois);

	char buf[500];
	for(int i = 0;i < num;i++)
	{
		DImage tmp;
		sprintf(buf,"%s\\%s%d.%s",input_fold,input_prefix,i+base_id,input_suffix);
		if(!ZQ_ImageIO::loadImage(tmp,buf,0))
		{
			printf("failed to load %s\n",buf);
			return;
		}
		if(!tmp.matchDimension(width,height,1))
		{
			printf("dimensions don't match!\n");
			return;
		}
		input_images.push_back(tmp);
	}

	for(int i = 0;i < num;i++)
	{
		DImage tmp;
		sprintf(buf,"%s\\%s%d.%s",copyin_fold,copyin_prefix,i+base_id,copyin_suffix);
		if(!ZQ_ImageIO::loadImage(tmp,buf,0))
		{
			printf("failed to load %s\n",buf);
			return;
		}
		if(!tmp.matchDimension(width,height,1))
		{
			printf("dimensions don't match!\n");
			return;
		}
		copyin_images.push_back(tmp);
	}

	ZQ_PoissonEditingOptions opt;
	opt.type = ZQ_PoissonEditingOptions::METHOD_NAIVE;
	opt.grad_scale = 1;
	opt.nSORIteration = 500;
	opt.display = true;

	DImage tmp_output;
	ZQ_PoissonEditing::PoissonEditing(mask,copyin_images[0],input_images[0],tmp_output,opt);
	input_images[0] = tmp_output;
	ZQ_PoissonEditing::PoissonEditing(mask,copyin_images[num-1],input_images[num-1],tmp_output,opt);
	input_images[num-1] = tmp_output;

	DImage3D output3D(width,height,num);
	float*& output_data = output3D.data();
	for(int i = 0;i < num;i++)
	{
		float*& inputData = input_images[i].data();
		memcpy(output_data+i*width*height,inputData,sizeof(float)*width*height);
	}

	for(int i = 0;i < rois.size();i++)
	{
		_solve_sub_problem(input_images,copyin_images,output3D,rois[i]);
	}

	
	for(int i = 0;i < num;i++)
	{
		DImage tmp(width,height);
		for(int pp = 0;pp < width*height;pp++)
			tmp.data()[pp] = output_data[i*width*height+pp];

		sprintf(buf,"%s\\%s%d.%s",output_fold,input_prefix,i+base_id,input_suffix);
		if(!ZQ_ImageIO::saveImage(tmp,buf))
		{
			printf("failed to save %s\n",buf);
		}
	}

	return;
}

void _find_roi(const ZQ_DImage<bool>& mask, std::vector<ROI>& rois)
{
	rois.clear();
	int width = mask.width();
	int height = mask.height();
	const bool*& mask_data = mask.data();

	if(width == 0 || height == 0)
		return;

	bool* visit_flag = new bool[width*height];
	memset(visit_flag,0,sizeof(bool)*width*height);
	bool done_flag = false;
	while(!done_flag)
	{
		done_flag = true;
		int find_h = 0;
		int find_w = 0;
		for(int h = 0;h < height;h++)
		{
			for(int w = 0;w < width;w++)
			{
				if(!visit_flag[h*width+w] && mask_data[h*width+w])	
				{
					done_flag = false;
					find_h = h;
					find_w = w;
					break;
				}
			}
			if(!done_flag)
				break;
		}

		if(done_flag)
			break;

		bool* cur_mask = new bool[width*height];
		bool* add_tail_flag = new bool[width*height];
		memset(add_tail_flag,0,sizeof(bool)*width*height);
		memset(cur_mask,0,sizeof(bool)*width*height);
		int* index_x = new int[width*height];
		int* index_y = new int[width*height];
		int head = -1;
		int tail = -1;
		tail++;
		index_x[tail] = find_w;
		index_y[tail] = find_h;
		add_tail_flag[find_h*width+find_w] = true;
		while(head < tail)
		{
			head++;
			int cur_x = index_x[head];
			int cur_y = index_y[head];
			visit_flag[cur_y*width+cur_x] = true;
			cur_mask[cur_y*width+cur_x] = true;
			if(cur_y+1 < height && !add_tail_flag[(cur_y+1)*width+cur_x] && mask_data[(cur_y+1)*width+cur_x])
			{
				tail++;
				index_x[tail] = cur_x;
				index_y[tail] = cur_y+1;
				add_tail_flag[(cur_y+1)*width+cur_x] = true;
			}
			if(cur_y-1 >= 0 && !add_tail_flag[(cur_y-1)*width+cur_x] && mask_data[(cur_y-1)*width+cur_x])
			{
				tail++;
				index_x[tail] = cur_x;
				index_y[tail] = cur_y-1;
				add_tail_flag[(cur_y-1)*width+cur_x] = true;
			}
			if(cur_x+1 < width && !add_tail_flag[cur_y*width+cur_x+1] && mask_data[cur_y*width+cur_x+1])
			{
				tail++;
				index_x[tail] = cur_x+1;
				index_y[tail] = cur_y;
				add_tail_flag[cur_y*width+cur_x+1] = true;
			}
			if(cur_x-1 >= 0 && !add_tail_flag[cur_y*width+cur_x-1] && mask_data[cur_y*width+cur_x-1])
			{
				tail++;
				index_x[tail] = cur_x-1;
				index_y[tail] = cur_y;
				add_tail_flag[cur_y*width+cur_x-1] = true;
			}
		}
		int min_x = -1;
		int min_y = -1;
		int max_x = width;
		int max_y = height;
		for(int h = 0;h < height;h++)
		{
			for(int w = 0;w < width;w++)
			{
				if(cur_mask[h*width+w])
				{
					if(w < min_x)
						min_x = w;
					if(w > max_x)
						max_x = w;
					if(h < min_y)
						min_y = h;
					if(h > max_y)
						max_y = h;
				}
			}
		}
		min_x = __max(min_x-1,0);
		min_y = __max(min_y-1,0);
		max_x = __min(max_x+1,width-1);
		max_y = __min(max_y+1,height-1);

		int roi_width = max_x-min_x+1;
		int roi_height = max_y-min_y+1;
		int roi_off_x = min_x;
		int roi_off_y = min_y;

		ROI cur_roi;
		rois.push_back(cur_roi);
		int roi_num = rois.size();
		rois[roi_num-1].width = roi_width;
		rois[roi_num-1].height = roi_height;
		rois[roi_num-1].off_x = roi_off_x;
		rois[roi_num-1].off_y = roi_off_y;
		rois[roi_num-1].mask.allocate(roi_width,roi_height);
		bool*& roi_mask_data = rois[roi_num-1].mask.data();

		for(int hh = 0;hh < roi_height;hh++)
		{
			for(int ww = 0;ww < roi_width;ww++)
			{
				roi_mask_data[hh*roi_width+ww] = cur_mask[(hh+roi_off_y)*width+(ww+roi_off_x)];
			}
		}
		delete []cur_mask;
		delete []index_x;
		delete []index_y;
		delete []add_tail_flag;
	}
	delete []visit_flag;
}

bool _solve_sub_problem(const std::vector<DImage>& input_images,const std::vector<DImage>& copyin_images,DImage3D& output3D, const ROI& roi)
{
	int num = input_images.size();
	if(num == 0 || copyin_images.size() != num)
		return false;

	int roi_width = roi.width;
	int roi_height = roi.height;
	int roi_off_x = roi.off_x;
	int roi_off_y = roi.off_y;
	
	int width = input_images[0].width();
	int height = input_images[0].height();
	int nChannels = input_images[0].nchannels();
	
	if(!output3D.matchDimension(width,height,num,nChannels))
		return false;

	for(int i = 0;i < num;i++)
	{
		if(!input_images[i].matchDimension(width,height,nChannels)
			|| !copyin_images[i].matchDimension(width,height,nChannels))
			return false;
	}

	DImage3D input3D(roi_width,roi_height,num,nChannels);
	DImage3D copyin3D(roi_width,roi_height,num,nChannels);
	ZQ_DImage3D<bool> mask3D(roi_width,roi_height,num,1);

	float*& input3D_data = input3D.data();
	float*& copyin3D_data = copyin3D.data();
	bool*& mask3D_data = mask3D.data();
	const bool*& maskData = roi.mask.data();
	for(int i = 0;i < num;i++)
	{
		const float*& inputData = input_images[i].data();
		const float*& copyinData = copyin_images[i].data();
		for(int hh = 0;hh < roi_height;hh++)
		{
			for(int ww = 0;ww < roi_width;ww++)
			{
				memcpy(input3D_data + (i*roi_width*roi_height + hh*roi_width + ww)*nChannels, inputData + ((hh+roi_off_y)*width+(ww+roi_off_x))*nChannels, sizeof(float)*nChannels);
				memcpy(copyin3D_data + (i*roi_width*roi_height + hh*roi_width + ww)*nChannels, copyinData + ((hh+roi_off_y)*width+(ww+roi_off_x))*nChannels, sizeof(float)*nChannels);
			}
		}
	}
	for(int i = 1;i < num-1;i++)
	{
		for(int pp = 0;pp < roi_width*roi_height;pp++)
		{
			mask3D_data[i*roi_width*roi_height+pp] = maskData[pp];
		}
	}

	DImage3D tmpout3D;
	int nIteration3D = 5000;
	if(!PoissonEditing3DCuda(mask3D,input3D,copyin3D,tmpout3D,nIteration3D))
	{
		printf("poisson edting fail\n");
		return false;
	}

	float*& output3D_data = output3D.data();
	float*& tmp3D_data = tmpout3D.data();
	for(int i = 1;i < num-1;i++)
	{
		for(int hh = 0;hh < roi_height;hh++)
		{
			for(int ww = 0;ww < roi_width;ww++)
			{
				if(mask3D_data[i*roi_width*roi_height + hh*roi_width + ww])
				{
					output3D_data[i*width*height + (roi_off_y+hh)*width + (roi_off_x+ww)] = tmp3D_data[i*roi_width*roi_height + hh*roi_width + ww];
				}
			}
		}
	}
	return true;
}

bool PoissonEditing3DCuda(const ZQ_DImage3D<bool>& mask, const DImage3D& input, const DImage3D& copyin, DImage3D& output, const int nIteration)
{
	int width = input.width();
	int height = input.height();
	int depth = input.depth();
	int nChannels = input.nchannels();

	if(!mask.matchDimension(width,height,depth,1) 
		|| !copyin.matchDimension(width,height,depth,nChannels))
		return false;

	output.copyData(input);

	DImage3D lap(width,height,depth,nChannels);

	const float*& copyin_data = copyin.data();
	const bool*& mask_data = mask.data();
	float*& output_data = output.data();
	float*& lap_data = lap.data();

	ZQ_ImageProcessing3D::Laplacian(copyin_data,lap_data,width,height,depth,nChannels);

	ZQ_CUDA_PoissonEditing3D(mask_data,lap_data,output_data,width,height,depth,nChannels,nIteration);

	return true;
}