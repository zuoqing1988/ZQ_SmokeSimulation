#ifndef _ZQ_SMOKE_EDITING_CORNERS_FIGGRID_H_
#define _ZQ_SMOKE_EDITING_CORNERS_FIGGRID_H_

#include "ZQ_SmokeEditingFigGrid.h"

namespace ZQ_SmokeEditing
{
	template<const unsigned int N, const bool Hflag>
	class CornersDeformFigGrid: public DeformFigGrid
	{
	public:
		CornersDeformFigGrid() : flag(Hflag)
		{
			if(flag)
			{
				width = 33;
				height = 16*N+1;
			}
			else
			{
				height = 33;
				width = 16*N+1;
			}
			
			_apply();
		}
		~CornersDeformFigGrid(){}

	private:
		const bool flag;
		
		void _createTweakBar()
		{
			if(bar == 0)
			{
				char name[200];
				sprintf(name,"Corners%d",N*2);
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
			ctrl_pts_num = N*2;
			ctrl_pts_x.clear();
			ctrl_pts_y.clear();
			if(flag)
			{
				for(int i = 0;i < N;i++)
				{
					ctrl_pts_x.push_back(0);
					ctrl_pts_y.push_back((height-1)*i/(N-1));

					ctrl_pts_x.push_back(width-1);
					ctrl_pts_y.push_back((height-1)*i/(N-1));
				}
			}
			else
			{
				for(int i = 0;i < N;i++)
				{
					ctrl_pts_x.push_back((width-1)*i/(N-1));
					ctrl_pts_y.push_back(0);

					ctrl_pts_x.push_back((width-1)*i/(N-1));
					ctrl_pts_y.push_back(height-1);
				}
			}
			ctrl_id = 0;
		}

		void _init_fixed_pts()
		{
			if(flag)
			{
				for(int i = 0;i < N;i++)
				{
					for(int j = 0;j < width;j++)
					{
						fixed_flag.data()[ctrl_pts_y[i*2]*width+j] = true;
					}
				}
			}
			else
			{
				for(int j = 0;j < N;j++)
				{
					for(int i = 0;i < height;i++)
					{
						fixed_flag.data()[i*width+ctrl_pts_x[j*2]] = true;
					}
				}
			}
		}

		void _init_coords()
		{
			int w = width;
			int h = height;
			
			coords.allocate(w,h,2);
			float*& coord_ptr = coords.data();
			float shift_x = -0.5*(w-1);
			float shift_y = -0.5*(h-1);
			for(int i = 0;i < h;i++)
			{
				for(int j = 0;j < w;j++)
				{
					coord_ptr[(i*w+j)*2+0] = j+shift_x;
					coord_ptr[(i*w+j)*2+1] = i+shift_y;
				}
			}

		}


		void _recompute_fixed()
		{
			int W = coords.width();
			int H = coords.height();
			float*& coord_ptr = coords.data();

			if(flag)
			{
				for(int i = 0;i < N;i++)
				{
					int y_id = ctrl_pts_y[i*2];
					float pt0[2] = {coord_ptr[(y_id*W)*2],coord_ptr[(y_id*W)*2+1]};
					float pt1[2] = {coord_ptr[(y_id*W+W-1)*2],coord_ptr[(y_id*W+W-1)*2+1]};
					for(int j = 1;j < W-1;j++)
					{
						float weight2 = (float)j/(W-1);
						float weight1 = 1.0-weight2;
						coord_ptr[(y_id*W+j)*2+0] = weight1*pt0[0]+weight2*pt1[0];
						coord_ptr[(y_id*W+j)*2+1] = weight1*pt0[1]+weight2*pt1[1];
					}
				}
			}
			else
			{
				for(int i = 0;i < N;i++)
				{
					int x_id = ctrl_pts_x[i*2];
					float pt0[2] = {coord_ptr[(0*W+x_id)*2],coord_ptr[(0*W+x_id)*2+1]};
					float pt1[2] = {coord_ptr[((H-1)*W+x_id)*2],coord_ptr[((H-1)*W+x_id)*2+1]};
					for(int j = 1;j < H-1;j++)
					{
						float weight2 = (float)j/(H-1);
						float weight1 = 1.0-weight2;
						coord_ptr[(j*W+x_id)*2+0] = weight1*pt0[0]+weight2*pt1[0];
						coord_ptr[(j*W+x_id)*2+1] = weight1*pt0[1]+weight2*pt1[1];
					}
				}
			}
		}
	};
}

#endif