#ifndef _ZQ_SMOKE_EDITING_XLOOP_FIGGRID_3D_H_
#define _ZQ_SMOKE_EDITING_XLOOP_FIGGRID_3D_H_
#pragma once

#include "ZQ_SmokeEditingFigGrid3D.h"

namespace ZQ_SmokeEditing3D
{
	template<const int N, const int TYPE>
	class XloopDeformGrid3D: public DeformFigGrid3D
	{
	public:
		XloopDeformGrid3D(): type(TYPE)
		{
			arap = true;
			width = 32;
			height = 4;
			depth = 4;
			xloop = true;
			_apply();
		}
		~XloopDeformGrid3D(){}

	private:

		const int type;

		void _createTweakBar()
		{
			if(bar == 0)
			{
				const char* name = "XloopDeformGrid3D";
				TwBar* tmp_bar = TwGetBarByName(name);
				if(tmp_bar != 0)
					TwDeleteBar(tmp_bar);

				bar = TwNewBar(name);
				TwSetTopBar(bar);
				char buf[200];
				sprintf(buf,"%s refresh=0.1 movable=false resizable=false color='0 0 0' alpha='100' position='0 150'  size='200 600'",name);
				TwDefine(buf);
			}
		}

		void _init_ctrl_pts()
		{
			switch(type)
			{
			case 0:
			case 1:
			default:
				{
					ctrl_pts_num = N*4;
					ctrl_pts_x.clear();
					ctrl_pts_y.clear();
					ctrl_pts_z.clear();
					for(int i = 0;i < N;i++)
					{
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(0);
						ctrl_pts_z.push_back(0);
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(height-1);
						ctrl_pts_z.push_back(0);
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(height-1);
						ctrl_pts_z.push_back(depth-1);
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(0);
						ctrl_pts_z.push_back(depth-1);
					}
					ctrl_id = 0;
				}
				break;
			}

		}

		void _init_fixed_pts()
		{

			switch(type)
			{
			case 0:
				{
					for(int ii = 0;ii < width;ii++)
					{
						bool tmp_fix_flag = false;
						for(int nn = 0;nn < N;nn++)
						{
							if(ii == ctrl_pts_x[4*nn])
							{
								tmp_fix_flag = true;
								break;
							}	
						}
						if(tmp_fix_flag)
						{
							for(int kk = 0;kk <= depth-1;kk++)
							{
								fixed_flag.data()[kk*height*width+0*width+ii] = true;
								fixed_flag.data()[kk*height*width+(height-1)*width+ii] = true;
							}
							for(int jj = 0;jj <= height-1;jj++)
							{
								fixed_flag.data()[0*height*width+jj*width+ii] = true;
								fixed_flag.data()[(depth-1)*height*width+jj*width+ii] = true;
							}
						}
					}

				}
				break;
			case 1:
			default:
				{
					for(int pp = 0;pp < ctrl_pts_num;pp++)
					{
						int offset = ctrl_pts_z[pp]*height*width + ctrl_pts_y[pp]*width + ctrl_pts_x[pp];
						fixed_flag.data()[offset] = true;
					}
				}
				break;
			}

		}

		void _init_coords()
		{
			int w = width;
			int h = height;
			int d = depth;
			float*& coord_ptr = coords.data();

			const float m_pi = 3.1415926535;
			float radius = w/(2.0*m_pi);
			float shift_r = -0.5*(h-1);
			float shift_z = -0.5*(d-1);
			for(int ii = 0;ii < w;ii++)
			{
				float angle = (float)ii/w * 2 * m_pi;

				for(int jj = 0;jj < h;jj++)
				{
					float cur_r = jj+shift_r+radius;
					float x = cos(angle)*cur_r;
					float y = sin(angle)*cur_r;
					for(int kk = 0;kk < d;kk++)
					{
						coord_ptr[(kk*h*w+jj*w+ii)*3+0] = x;
						coord_ptr[(kk*h*w+jj*w+ii)*3+1] = y;
						coord_ptr[(kk*h*w+jj*w+ii)*3+2] = kk+shift_z;
					}
				}
			}
		}

		void _recompute_fixed()
		{
			switch(type)
			{
			case 0:
				{
					int W = coords.width();
					int H = coords.height();
					int D = coords.depth();
					float*& coord_ptr = coords.data();

					for(int nn = 0;nn < N;nn++)
					{
						int ii = ctrl_pts_x[nn*4];
						int offset_00 = ii;
						int offset_01 = (H-1)*W+ii;
						int offset_10 = (D-1)*H*W+ii;
						int offset_11 = (D-1)*H*W+(H-1)*W+ii;
						float pt00[3] = {coord_ptr[offset_00*3+0],coord_ptr[offset_00*3+1],coord_ptr[offset_00*3+2]};
						float pt01[3] = {coord_ptr[offset_01*3+0],coord_ptr[offset_01*3+1],coord_ptr[offset_01*3+2]};
						float pt10[3] = {coord_ptr[offset_10*3+0],coord_ptr[offset_10*3+1],coord_ptr[offset_10*3+2]};
						float pt11[3] = {coord_ptr[offset_11*3+0],coord_ptr[offset_11*3+1],coord_ptr[offset_11*3+2]};
						
						for(int jj = 1;jj < H-1;jj++)
						{
							float weight2 = (float)jj/(H-1);
							float weight1 = 1.0-weight2;
							for(int cc = 0;cc < 3;cc++)
							{
								coord_ptr[(0*H*W+jj*W+ii)*3+cc] = weight1*pt00[cc] + weight2*pt01[cc];
								coord_ptr[((D-1)*H*W+jj*W+ii)*3+cc] = weight1*pt10[cc] + weight2*pt11[cc];
							}
						}
						for(int kk = 1;kk < D-1;kk++)
						{
							float weight2 = (float)kk/(D-1);
							float weight1 = 1.0-weight2;
							for(int cc = 0;cc < 3;cc++)
							{
								coord_ptr[(kk*H*W+0*W+ii)*3+cc] = weight1*pt00[cc] + weight2*pt10[cc];
								coord_ptr[(kk*H*W+(H-1)*W+ii)*3+cc] = weight1*pt01[cc] + weight2*pt11[cc];
							}
						}
					}
				}
				break;
			}
		}
		void _addLocalTweakBar()
		{
			TwSetParam(bar,"ARAP","readonly",TW_PARAM_CSTRING,1,"true");
		}
	};
}

#endif