#ifndef _ZQ_SMOKE_EDITING_CORNERS_FIGGRID_3D_H_
#define _ZQ_SMOKE_EDITING_CORNERS_FIGGRID_3D_H_
#pragma once

#include "ZQ_SmokeEditingFigGrid3D.h"

namespace ZQ_SmokeEditing3D
{
	template<const int N, const int Type, const bool Hflag>
	class CornersFigGrid3D: public DeformFigGrid3D
	{
	public:
		CornersFigGrid3D():flag(Hflag),type(Type)
		{
			arap = true;
			if(flag)
			{
				width = 9;
				height = 4*N+1;
				depth = 9;
			}
			else
			{
				width = 4*N+1;
				height = 9;
				depth = 9;
			}
	
			_apply();
		}
		~CornersFigGrid3D(){}

	private:
		const bool flag;
		const int type;

		void _createTweakBar()
		{
			if(bar == 0)
			{
				char name[200];
				sprintf(name,"Corners8");
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
			ctrl_pts_x.clear();
			ctrl_pts_y.clear();
			ctrl_pts_z.clear();

			switch(type)
			{
			case 0:
			case 1:
			default:
				{
					ctrl_pts_num = N*4;
					for(int nn = 0;nn < N;nn++)
					{
						if(flag)
						{
							int iy = (height-1)*nn/(N-1);
							ctrl_pts_x.push_back(0);
							ctrl_pts_y.push_back(iy);
							ctrl_pts_z.push_back(0);
							ctrl_pts_x.push_back(0);
							ctrl_pts_y.push_back(iy);
							ctrl_pts_z.push_back(depth-1);
							ctrl_pts_x.push_back(width-1);
							ctrl_pts_y.push_back(iy);
							ctrl_pts_z.push_back(0);
							ctrl_pts_x.push_back(width-1);
							ctrl_pts_y.push_back(iy);
							ctrl_pts_z.push_back(depth-1);
						}
						else
						{
							int ix = (width-1)*nn/(N-1);
							ctrl_pts_x.push_back(ix);
							ctrl_pts_y.push_back(0);
							ctrl_pts_z.push_back(0);
							ctrl_pts_x.push_back(ix);
							ctrl_pts_y.push_back(0);
							ctrl_pts_z.push_back(depth-1);
							ctrl_pts_x.push_back(ix);
							ctrl_pts_y.push_back(height-1);
							ctrl_pts_z.push_back(0);
							ctrl_pts_x.push_back(ix);
							ctrl_pts_y.push_back(height-1);
							ctrl_pts_z.push_back(depth-1);
						}
					}
				}
				break;
			}
			
			ctrl_id = 0;
		}

		void _init_fixed_pts()
		{
			switch(type)
			{
			case 0:
				{
					for(int nn = 0;nn < N;nn++)
					{
						if(flag)
						{
							int iy = ctrl_pts_y[nn*4];
							for(int ww = 0;ww < width;ww++)
							{
								fixed_flag.data()[0*height*width+iy*width+ww] = true;
								fixed_flag.data()[(depth-1)*height*width+iy*width+ww] = true;
							}
							for(int dd = 0;dd < depth;dd++)
							{
								fixed_flag.data()[dd*height*width+iy*width+0] = true;
								fixed_flag.data()[dd*height*width+iy*width+width-1] = true;
							}
						}
						else
						{
							int ix = ctrl_pts_x[nn*4];
							for(int hh = 0;hh < height;hh++)
							{
								fixed_flag.data()[0*height*width+hh*width+ix] = true;
								fixed_flag.data()[(depth-1)*height*width+hh*width+ix] = true;
							}
							for(int dd = 0;dd < depth;dd++)
							{
								fixed_flag.data()[dd*height*width+0*width+ix] = true;
								fixed_flag.data()[dd*height*width+(height-1)*width+ix] = true;
							}
						}
					}
				}
				break;
			case 1:
			default:
				{
					for(int i = 0;i < ctrl_pts_x.size();i++)
					{
						int offset = ctrl_pts_z[i]*height*width+ctrl_pts_y[i]*width+ctrl_pts_x[i];
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

			coords.allocate(w,h,d,3);
			float*& coord_ptr = coords.data();
			float shift_x = -0.5*(w-1);
			float shift_y = -0.5*(h-1);
			float shift_z = -0.5*(d-1);
			for(int k = 0;k < d;k++)
			{
				for(int j = 0;j < h;j++)
				{
					for(int i = 0;i < w;i++)
					{
						coord_ptr[(k*h*w+j*w+i)*3+0] = i+shift_x;
						coord_ptr[(k*h*w+j*w+i)*3+1] = j+shift_y;
						coord_ptr[(k*h*w+j*w+i)*3+2] = k+shift_z;
					}
				}
			}
		}


		void _recompute_fixed()
		{
			int W = coords.width();
			int H = coords.height();
			int D = coords.depth();
			float*& coord_ptr = coords.data();

			switch(type)
			{
			case 0:
				{
					for(int nn = 0;nn < N;nn++)
					{
						if(flag)
						{
							int iy = ctrl_pts_y[nn*4];
							float pt0[3] = {coord_ptr[(0*H*W+iy*W+0)*3],coord_ptr[(0*H*W+iy*W+0)*3+1],coord_ptr[(0*H*W+iy*W+0)*3+2]};
							float pt1[3] = {coord_ptr[(0*H*W+iy*W+W-1)*3],coord_ptr[(0*H*W+iy*W+W-1)*3+1],coord_ptr[(0*H*W+iy*W+W-1)*3+2]};
							float pt2[3] = {coord_ptr[((D-1)*H*W+iy*W+0)*3],coord_ptr[((D-1)*H*W+iy*W+0)*3+1],coord_ptr[((D-1)*H*W+iy*W+0)*3+2]};
							float pt3[3] = {coord_ptr[((D-1)*H*W+iy*W+W-1)*3],coord_ptr[((D-1)*H*W+iy*W+W-1)*3+1],coord_ptr[((D-1)*H*W+iy*W+W-1)*3+2]};
							for(int ww = 1;ww < W-1;ww++)
							{
								float weight2 = (float)ww/(W-1);
								float weight1 = 1.0-weight2;
								coord_ptr[(0*H*W+iy*W+ww)*3+0] = weight1*pt0[0]+weight2*pt1[0];
								coord_ptr[(0*H*W+iy*W+ww)*3+1] = weight1*pt0[1]+weight2*pt1[1];
								coord_ptr[(0*H*W+iy*W+ww)*3+2] = weight1*pt0[2]+weight2*pt1[2];
								coord_ptr[((D-1)*H*W+iy*W+ww)*3+0] = weight1*pt2[0]+weight2*pt3[0];
								coord_ptr[((D-1)*H*W+iy*W+ww)*3+1] = weight1*pt2[1]+weight2*pt3[1];
								coord_ptr[((D-1)*H*W+iy*W+ww)*3+2] = weight1*pt2[2]+weight2*pt3[2];
							}
							for(int dd = 1;dd < D-1;dd++)
							{
								float weight2 = (float)dd/(D-1);
								float weight1 = 1.0-weight2;
								coord_ptr[(dd*H*W+iy*W+0)*3+0] = weight1*pt0[0]+weight2*pt2[0];
								coord_ptr[(dd*H*W+iy*W+0)*3+1] = weight1*pt0[1]+weight2*pt2[1];
								coord_ptr[(dd*H*W+iy*W+0)*3+2] = weight1*pt0[2]+weight2*pt2[2];
								coord_ptr[(dd*H*W+iy*W+W-1)*3+0] = weight1*pt1[0]+weight2*pt3[0];
								coord_ptr[(dd*H*W+iy*W+W-1)*3+1] = weight1*pt1[1]+weight2*pt3[1];
								coord_ptr[(dd*H*W+iy*W+W-1)*3+2] = weight1*pt1[2]+weight2*pt3[2];
							}
						}
						else
						{
							int ix = ctrl_pts_x[nn*4];
							float pt0[3] = {coord_ptr[(0*H*W+0*W+ix)*3],coord_ptr[(0*H*W+0*W+ix)*3+1],coord_ptr[(0*H*W+0*W+ix)*3+2]};
							float pt1[3] = {coord_ptr[(0*H*W+(H-1)*W+ix)*3],coord_ptr[(0*H*W+(H-1)*W+ix)*3+1],coord_ptr[(0*H*W+(H-1)*W+ix)*3+2]};
							float pt2[3] = {coord_ptr[((D-1)*H*W+0*W+ix)*3],coord_ptr[((D-1)*H*W+0*W+ix)*3+1],coord_ptr[((D-1)*H*W+0*W+ix)*3+2]};
							float pt3[3] = {coord_ptr[((D-1)*H*W+(H-1)*W+ix)*3],coord_ptr[((D-1)*H*W+(H-1)*W+ix)*3+1],coord_ptr[((D-1)*H*W+(H-1)*W+ix)*3+2]};
							for(int hh = 1;hh < H-1;hh++)
							{
								float weight2 = (float)hh/(H-1);
								float weight1 = 1.0-weight2;
								coord_ptr[(0*H*W+hh*W+ix)*3+0] = weight1*pt0[0]+weight2*pt1[0];
								coord_ptr[(0*H*W+hh*W+ix)*3+1] = weight1*pt0[1]+weight2*pt1[1];
								coord_ptr[(0*H*W+hh*W+ix)*3+2] = weight1*pt0[2]+weight2*pt1[2];
								coord_ptr[((D-1)*H*W+hh*W+ix)*3+0] = weight1*pt2[0]+weight2*pt3[0];
								coord_ptr[((D-1)*H*W+hh*W+ix)*3+1] = weight1*pt2[1]+weight2*pt3[1];
								coord_ptr[((D-1)*H*W+hh*W+ix)*3+2] = weight1*pt2[2]+weight2*pt3[2];
							}
							for(int dd = 1;dd < D-1;dd++)
							{
								float weight2 = (float)dd/(D-1);
								float weight1 = 1.0-weight2;
								coord_ptr[(dd*H*W+0*W+ix)*3+0] = weight1*pt0[0]+weight2*pt2[0];
								coord_ptr[(dd*H*W+0*W+ix)*3+1] = weight1*pt0[1]+weight2*pt2[1];
								coord_ptr[(dd*H*W+0*W+ix)*3+2] = weight1*pt0[2]+weight2*pt2[2];
								coord_ptr[(dd*H*W+(H-1)*W+ix)*3+0] = weight1*pt1[0]+weight2*pt3[0];
								coord_ptr[(dd*H*W+(H-1)*W+ix)*3+1] = weight1*pt1[1]+weight2*pt3[1];
								coord_ptr[(dd*H*W+(H-1)*W+ix)*3+2] = weight1*pt1[2]+weight2*pt3[2];
							}
						}
					}
				}
				break;
			case 1:
			default:
				{
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