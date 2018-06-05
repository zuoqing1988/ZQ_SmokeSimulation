#ifndef _ZQ_SMOKE_EDITING_XLOOP_FIGGRID_H_
#define _ZQ_SMOKE_EDITING_XLOOP_FIGGRID_H_

#include "ZQ_SmokeEditingFigGrid.h"

namespace ZQ_SmokeEditing
{
	template<const int N, const int TYPE>
	class XloopDeformGrid: public DeformFigGrid
	{
	public:
		XloopDeformGrid(): type(TYPE)
		{
			width = 64;
			height = 9;
			xloop = true;
			_apply();
		}
		~XloopDeformGrid(){}

	private:

		const int type;

		void _createTweakBar()
		{
			if(bar == 0)
			{
				const char* name = "XloopDeformGrid";
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
			case 1:
				{
					ctrl_pts_num = N*2;
					ctrl_pts_x.clear();
					ctrl_pts_y.clear();
					for(int i = 0;i < N;i++)
					{
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(0);
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(height-1);
					}
					ctrl_id = 0;
				}
				break;
			case 2:
				{
					ctrl_pts_num = N;
					ctrl_pts_x.clear();
					ctrl_pts_y.clear();
					for(int i = 0;i < N;i++)
					{
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(0);
					}
					ctrl_id = 0;
				}
				break;
			case 3:
				{
					ctrl_pts_num = N;
					ctrl_pts_x.clear();
					ctrl_pts_y.clear();
					for(int i = 0;i < N;i++)
					{
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(height-1);
					}
					ctrl_id = 0;
				}
				break;
			case 0:
			default:
				{
					ctrl_pts_num = N*2;
					ctrl_pts_x.clear();
					ctrl_pts_y.clear();
					for(int i = 0;i < N;i++)
					{
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(0);
						ctrl_pts_x.push_back(i*width/N);
						ctrl_pts_y.push_back(height-1);
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
			case 1:
				{
					for(int j = 0;j < width;j++)
					{
						bool tmp_fix_flag = false;
						for(int nn = 0;nn < N;nn++)
						{
							if(j == ctrl_pts_x[2*nn])
							{
								tmp_fix_flag = true;
								break;
							}	
						}
						if(tmp_fix_flag)
						{
							fixed_flag.data()[0*width+j] = true;
							fixed_flag.data()[(height-1)*width+j] = true;
						}
					}

				}
				break;
			case 2:
				{
					for(int j = 0;j < width;j++)
					{
						bool tmp_fix_flag = false;
						for(int nn = 0;nn < N;nn++)
						{
							if(j == ctrl_pts_x[nn])
							{
								tmp_fix_flag = true;
								break;
							}	
						}
						if(tmp_fix_flag)
						{
							fixed_flag.data()[0*width+j] = true;
						}
					}
				}
				break;
			case 3: 
				{
					for(int j = 0;j < width;j++)
					{
						bool tmp_fix_flag = false;
						for(int nn = 0;nn < N;nn++)
						{
							if(j == ctrl_pts_x[nn])
							{
								tmp_fix_flag = true;
								break;
							}	
						}
						if(tmp_fix_flag)
						{
							fixed_flag.data()[(height-1)*width+j] = true;
						}
					}
				}
				break;
			case 0:
			default:
				{
					for(int j = 0;j < width;j++)
					{
						bool tmp_fix_flag = false;
						for(int nn = 0;nn < N;nn++)
						{
							if(j == ctrl_pts_x[2*nn])
							{
								tmp_fix_flag = true;
								break;
							}	
						}
						if(tmp_fix_flag)
						{
							for(int i = 0;i < height;i++)
							{
								fixed_flag.data()[i*width+j] = true;
							}
						}
					}
					
				}
				break;
			}

		}

		void _init_coords()
		{
			int w = width;
			int h = height;
			float*& coord_ptr = coords.data();

			const float m_pi = 3.1415926535;
			float radius = w/(2.0*m_pi);
			float shift_r = -0.5*(h-1);
			for(int j = 0;j < w;j++)
			{
				float angle = (float)j/w * 2 * m_pi;

				for(int i = 0;i < h;i++)
				{
					float cur_r = i+shift_r+radius;
					float x = cos(m_pi-angle)*cur_r;
					float y = sin(m_pi-angle)*cur_r;
					coord_ptr[(i*w+j)*2+0] = x;
					coord_ptr[(i*w+j)*2+1] = y;
				}
			}
		}

		void _recompute_fixed()
		{
			switch(type)
			{
			case 1:
				break;
			case 2:
				break;
			case 3:
				break;
			case 0:
			default:
				{
					int W = coords.width();
					int H = coords.height();
					float*& coord_ptr = coords.data();

					for(int i = 0;i < N;i++)
					{
						float pt0[2] = {coord_ptr[(0*W+ctrl_pts_x[i*2])*2+0],coord_ptr[(0*W+ctrl_pts_x[i*2])*2+1]};
						float pt1[2] = {coord_ptr[((H-1)*W+ctrl_pts_x[i*2])*2+0],coord_ptr[((H-1)*W+ctrl_pts_x[i*2])*2+1]};
						for(int j = 1;j < H-1;j++)
						{
							float weight2 = (float)j/(H-1);
							float weight1 = 1.0-weight2;
							coord_ptr[(j*W+ctrl_pts_x[i*2])*2+0] = weight1*pt0[0] + weight2*pt1[0];
							coord_ptr[(j*W+ctrl_pts_x[i*2])*2+1] = weight1*pt0[1] + weight2*pt1[1];
						}
					}
				}
			}
		}
	};
}

#endif