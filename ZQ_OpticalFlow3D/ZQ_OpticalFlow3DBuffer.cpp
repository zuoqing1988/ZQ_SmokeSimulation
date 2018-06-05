#include "ZQ_OpticalFlow3DBuffer.h"

using namespace ZQ;

ZQ_OpticalFlow3DBuffer::ZQ_OpticalFlow3DBuffer()
{
	for(int i = 0;i < IMAGE3D_BUFFER_SIZE;i++)
	{
		buf_image_index[i] = -1;
		buf_image_level[i] = -1;
		buf_image_live_time[i] = -1;
	}

	for(int i = 0;i < IMAGE3D_BUFFER_SIZE;i++)
	{
		buf_uvw_index[i] = -1;
		buf_uvw_level[i] = -1;
		buf_uvw_live_time[i] = -1;
	}

}

ZQ_OpticalFlow3DBuffer::~ZQ_OpticalFlow3DBuffer()
{

}

bool ZQ_OpticalFlow3DBuffer::SetFilePara(const char* in_fold, const char* out_fold, const char* img_prefix, int img_num, int base_idx)
{
	strcpy(input_fold,in_fold);
	strcpy(output_fold,out_fold);
	strcpy(image_prefix,img_prefix);
	image_num = img_num;
	base_index = base_idx;

	return true;
}

bool ZQ_OpticalFlow3DBuffer::BuildPyramids(double ratio,int minWidth, int& levels, float& cuda_cost_time)
{
	char tmp_buf[STRING_BUFFER_LENGTH];
	ZQ_GaussianPyramid3DCuda GPyramid3D;

	cuda_cost_time = 0;

	DImage3D tmp_img;
	for(int i = 0;i < image_num;i++)
	{
		sprintf(tmp_buf,"%s\\%s%d.di3",input_fold,image_prefix,i+base_index);

		bool flag = tmp_img.loadImage(tmp_buf);
		if(!flag)
		{
			printf("load input image %s fail\n",tmp_buf);
			return false;
		}

		float tmp_cuda_cost = 0;
		GPyramid3D.ConstructPyramid(tmp_img,tmp_cuda_cost,ratio,minWidth);
		cuda_cost_time += tmp_cuda_cost;

		if(i == 0)
		{
			levels = GPyramid3D.nlevels();
		}

		for(int k = 1;k < levels;k++)
		{
			sprintf(tmp_buf,"%s\\%s%d_%d.di3",input_fold,image_prefix,i+base_index,k);
			tmp_img = GPyramid3D.Image(k);
			if(!tmp_img.saveImage(tmp_buf))
			{
				printf("save image %s fail, MAYBE DISK FULL\n",tmp_buf);
				return false;
			}
		}
	}
	return true;
}

bool ZQ_OpticalFlow3DBuffer::ClearPyramids(int levels)
{

	char tmp_buf[STRING_BUFFER_LENGTH];
	for(int i = 0;i < image_num;i++)
	{
		for(int k = 1;k < levels;k++)
		{
			sprintf(tmp_buf,"%s\\%s%d_%d.di3",input_fold,image_prefix,i+base_index,k);
			if(remove(tmp_buf) == -1)
			{
				printf("warning: remove file %s fail\n",tmp_buf);
			}
		}
	}
	return true;
}

bool ZQ_OpticalFlow3DBuffer::ClearIntermediateUVW(int levels)
{
	char tmp_buf[STRING_BUFFER_LENGTH];
	for(int i = 0;i < image_num-1;i++)
	{
		for(int k = 0;k < levels;k++)
		{
			sprintf(tmp_buf,"%s\\%s%d_%d.di3",output_fold,"u_",i+base_index,k);
			if(remove(tmp_buf) == -1)
			{
				printf("warning: remove file %s fail\n",tmp_buf);
			}
			sprintf(tmp_buf,"%s\\%s%d_%d.di3",output_fold,"v_",i+base_index,k);
			if(remove(tmp_buf) == -1)
			{
				printf("warning: remove file %s fail\n",tmp_buf);
			}
			sprintf(tmp_buf,"%s\\%s%d_%d.di3",output_fold,"w_",i+base_index,k);
			if(remove(tmp_buf) == -1)
			{
				printf("warning: remove file %s fail\n",tmp_buf);
			}
		}
	}
	return true;
}

bool ZQ_OpticalFlow3DBuffer::GetImage(DImage3D& out_img, int level,int i)
{
	if(i < 0 || i >= image_num)
		return false;

	bool is_in_buf = false;

	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		if(buf_image_index[t] == i && buf_image_level[t] == level)
		{
			buf_image_live_time[t] = 0;
			out_img.copyData(buf_image[t]);
			is_in_buf = true;
			break;
		}
	}


	if(!is_in_buf)
	{
		char name_buf[STRING_BUFFER_LENGTH];
		if(level == 0)
		{
			sprintf(name_buf,"%s\\%s%d.di3",input_fold,image_prefix,i+base_index);
		}
		else
		{
			sprintf(name_buf,"%s\\%s%d_%d.di3",input_fold,image_prefix,i+base_index,level);
		}


		if(!out_img.loadImage(name_buf))
		{
			return false;
		}

		//replace buffer

		bool has_empty = false;
		int empty_idx = -1;
		for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
		{
			if(buf_image_index[t] < 0)
			{
				has_empty = true;
				empty_idx = t;
				break;
			}
		}

		if(has_empty)
		{
			buf_image[empty_idx].copyData(out_img);
			buf_image_index[empty_idx] = i;
			buf_image_level[empty_idx] = level;
			buf_image_live_time[empty_idx] = 0;
		}
		else
		{
			int max_live_time = buf_image_live_time[0];
			int max_live_idx = 0;
			for(int t = 1;t < IMAGE3D_BUFFER_SIZE;t++)
			{
				if(max_live_time < buf_image_live_time[t])
				{
					max_live_time = buf_image_live_time[t];
					max_live_idx = t;
				}
			}

			buf_image[max_live_idx].copyData(out_img);
			buf_image_index[max_live_idx] = i;
			buf_image_level[max_live_idx] = level;
			buf_image_live_time[max_live_idx] = 0;
		}
	}

	//refresh live time

	int pos[IMAGE3D_BUFFER_SIZE];
	int live_time[IMAGE3D_BUFFER_SIZE];
	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		pos[t] = t;
		live_time[t] = buf_image_live_time[t];
	}

	for(int p = 0;p < IMAGE3D_BUFFER_SIZE-1;p++)
	{
		for(int t = 0;t < IMAGE3D_BUFFER_SIZE-1;t++)
		{
			if(live_time[t] > live_time[t+1])
			{
				int tmp_live = live_time[t];
				live_time[t] = live_time[t+1];
				live_time[t+1] = tmp_live;

				int tmp_pos = pos[t];
				pos[t] = pos[t+1];
				pos[t+1] = tmp_pos;
			}
		}
	}

	int last_live_time = -1;
	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		if(live_time[t] < 0)
			continue;

		if(last_live_time == -1)
		{
			live_time[t] = 1;
			last_live_time = 1;
		}
		else
		{
			live_time[t] = last_live_time+1;
			last_live_time = live_time[t];
		}
	}

	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		buf_image_live_time[pos[t]] = live_time[t];
	}

	return true;
}

bool ZQ_OpticalFlow3DBuffer::GetUVW(DImage3D& u, DImage3D& v, DImage3D& w, int level, int i)
{
	if(i < 0 || i >= image_num-1)
		return false;

	bool is_in_buf = false;

	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		if(buf_uvw_index[t] == i && buf_uvw_level[t] == level)
		{
			buf_uvw_live_time[t] = 0;
			u.copyData(buf_u[t]);
			v.copyData(buf_v[t]);
			w.copyData(buf_w[t]);
			is_in_buf = true;
			break;
		}
	}


	if(!is_in_buf)
	{
		char name_buf[STRING_BUFFER_LENGTH];

		sprintf(name_buf,"%s\\%s%d_%d.di3",output_fold,"u_",i+base_index,level);
		if(!u.loadImage(name_buf))
		{
			return false;
		}
		sprintf(name_buf,"%s\\%s%d_%d.di3",output_fold,"v_",i+base_index,level);
		if(!v.loadImage(name_buf))
		{
			return false;
		}
		sprintf(name_buf,"%s\\%s%d_%d.di3",output_fold,"w_",i+base_index,level);
		if(!w.loadImage(name_buf))
		{
			return false;
		}

		//replace buffer

		bool has_empty = false;
		int empty_idx = -1;
		for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
		{
			if(buf_uvw_index[t] < 0)
			{
				has_empty = true;
				empty_idx = t;
				break;
			}
		}

		if(has_empty)
		{
			buf_u[empty_idx].copyData(u);
			buf_v[empty_idx].copyData(v);
			buf_w[empty_idx].copyData(w);
			buf_uvw_index[empty_idx] = i;
			buf_uvw_level[empty_idx] = level;
			buf_uvw_live_time[empty_idx] = 0;
		}
		else
		{
			int max_live_time = buf_uvw_live_time[0];
			int max_live_idx = 0;
			for(int t = 1;t < IMAGE3D_BUFFER_SIZE;t++)
			{
				if(max_live_time < buf_uvw_live_time[t])
				{
					max_live_time = buf_uvw_live_time[t];
					max_live_idx = t;
				}
			}

			buf_u[max_live_idx].copyData(u);
			buf_v[max_live_idx].copyData(v);
			buf_w[max_live_idx].copyData(w);
			buf_uvw_index[max_live_idx] = i;
			buf_uvw_level[max_live_idx] = level;
			buf_uvw_live_time[max_live_idx] = 0;
		}
	}

	//refresh live time

	int pos[IMAGE3D_BUFFER_SIZE];
	int live_time[IMAGE3D_BUFFER_SIZE];
	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		pos[t] = t;
		live_time[t] = buf_uvw_live_time[t];
	}

	for(int p = 0;p < IMAGE3D_BUFFER_SIZE-1;p++)
	{
		for(int t = 0;t < IMAGE3D_BUFFER_SIZE-1;t++)
		{
			if(live_time[t] > live_time[t+1])
			{
				int tmp_live = live_time[t];
				live_time[t] = live_time[t+1];
				live_time[t+1] = tmp_live;

				int tmp_pos = pos[t];
				pos[t] = pos[t+1];
				pos[t+1] = tmp_pos;
			}
		}
	}

	int last_live_time = -1;
	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		if(live_time[t] < 0)
			continue;

		if(last_live_time == -1)
		{
			live_time[t] = 1;
			last_live_time = 1;
		}
		else
		{
			live_time[t] = last_live_time+1;
			last_live_time = live_time[t];
		}
	}

	for(int t = 0;t < IMAGE3D_BUFFER_SIZE;t++)
	{
		buf_uvw_live_time[pos[t]] = live_time[t];
	}

	return true;
}

bool ZQ_OpticalFlow3DBuffer::WriteUVW(DImage3D& u, DImage3D& v, DImage3D& w, int level, int i)
{
	char tmp_buf[STRING_BUFFER_LENGTH];
	sprintf(tmp_buf,"%s\\%s%d_%d.di3",output_fold,"u_",i+base_index,level);
	if(!u.saveImage(tmp_buf))
	{	
		printf("save image %s fail, MAYBE DISK FULL\n",tmp_buf);
		return false;
	}
	sprintf(tmp_buf,"%s\\%s%d_%d.di3",output_fold,"v_",i+base_index,level);
	if(!v.saveImage(tmp_buf))
	{	
		printf("save image %s fail, MAYBE DISK FULL\n",tmp_buf);
		return false;
	}
	sprintf(tmp_buf,"%s\\%s%d_%d.di3",output_fold,"w_",i+base_index,level);
	if(!w.saveImage(tmp_buf))
	{	
		printf("save image %s fail, MAYBE DISK FULL\n",tmp_buf);
		return false;
	}

	return true;
}

bool ZQ_OpticalFlow3DBuffer::WriteFinalResult(DImage3D& flow, DImage3D& warp,int i)
{
	char tmp_buf[STRING_BUFFER_LENGTH];
	sprintf(tmp_buf,"%s\\flow%d_%d.di3",output_fold,i+base_index,i+1+base_index);
	if(!flow.saveImage(tmp_buf))
	{	
		printf("save flow %s fail, MAYBE DISK FULL\n",tmp_buf);
		return false;
	}
	sprintf(tmp_buf,"%s\\warp%d_%d.di3",output_fold,i+1+base_index,i+base_index);
	if(!warp.saveImage(tmp_buf))
	{	
		printf("save image %s fail, MAYBE DISK FULL\n",tmp_buf);
		return false;
	}
	return true;
}