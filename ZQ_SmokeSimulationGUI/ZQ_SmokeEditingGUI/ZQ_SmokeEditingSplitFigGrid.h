#ifndef _ZQ_SMOKE_EDITING_SPLIT_FIGGRID_H_
#define _ZQ_SMOKE_EDITING_SPLIT_FIGGRID_H_

#include "ZQ_SmokeEditingFigGrid.h"

namespace ZQ_SmokeEditing
{

	class SplitDeformGrid: public DeformFigGrid
	{
	public:
		SplitDeformGrid(const bool Uflag = true): flag(Uflag)
		{
			width = 33;
			height = 17;
			_apply();
			xloop = false;
		}

		~SplitDeformGrid(){}

	private:
		const bool flag;


		void _createTweakBar()
		{
			if(bar == 0)
			{
				const char* name = "SplitDeformGrid";
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

		void _addLocalTweakBar()
		{
			if(bar != 0)
			{
				TwAddButton(bar,"centerInc",(TwButtonCallback)OnCenterINC,this,"label='CenterINC'");
				TwAddButton(bar,"centerDec",(TwButtonCallback)OnCenterDEC,this,"label='CenterDEC'");
			}
		}

		static void TW_CALL OnCenterINC(void* clientData)
		{
			SplitDeformGrid* g = (SplitDeformGrid*)clientData;
			g->_centerInc();
		}

		static void TW_CALL OnCenterDEC(void* clientData)
		{
			SplitDeformGrid* g = (SplitDeformGrid*)clientData;
			g->_centerDec();
		}

		void _centerDec()
		{
			if(ctrl_pts_x[4] > 1)
			{
				ctrl_pts_x[4] --;
				_recompute_fixed();
				_reset_deform_iteration();
			}
		}

		void _centerInc()
		{
			int W = coords.width();
			if(ctrl_pts_x[4] < W-2)
			{
				ctrl_pts_x[4] ++;

				_recompute_fixed();
				_reset_deform_iteration();
			}
		}


		void _init_coords()
		{
			int w = width;
			int h = height;
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

		void _init_ctrl_pts()
		{
			ctrl_pts_num = 5;
			ctrl_pts_x.clear();
			ctrl_pts_y.clear();
			ctrl_pts_x.push_back(0);
			ctrl_pts_y.push_back(0);
			ctrl_pts_x.push_back(width-1);
			ctrl_pts_y.push_back(0);
			ctrl_pts_x.push_back(0);
			ctrl_pts_y.push_back(height-1);
			ctrl_pts_x.push_back(width-1);
			ctrl_pts_y.push_back(height-1);
			ctrl_pts_x.push_back(width/2);
			if(flag)
				ctrl_pts_y.push_back(height-1);
			else
				ctrl_pts_y.push_back(0);

			ctrl_id = 0;
		}

		void _init_fixed_pts()
		{
			for(int i = 0;i < width;i++)
			{
				fixed_flag.data()[0*width+i] = true;
				fixed_flag.data()[(height-1)*width+i] = true;
			}
		}


		void _recompute_fixed()
		{
			int W = coords.width();
			int H = coords.height();
			float*& coord_ptr = coords.data();
			int left_pts_space_num = ctrl_pts_x[4];
			int right_pts_space_num = W-1-ctrl_pts_x[4];

			if(flag)
			{
				float pt0[2] = {coord_ptr[((H-1)*W)*2],coord_ptr[((H-1)*W)*2+1]};
				float pt1[2] = {coord_ptr[((H-1)*W+ctrl_pts_x[4])*2],coord_ptr[((H-1)*W+ctrl_pts_x[4])*2+1]};
				float pt2[2] = {coord_ptr[((H-1)*W+W-1)*2],coord_ptr[((H-1)*W+W-1)*2+1]};

				for(int i = 1;i < ctrl_pts_x[4];i++)
				{
					float weight2 = (float)i/left_pts_space_num;
					float weight1 = 1.0-weight2;
					coord_ptr[((H-1)*W+i)*2+0] = weight1*pt0[0]+weight2*pt1[0];
					coord_ptr[((H-1)*W+i)*2+1] = weight1*pt0[1]+weight2*pt1[1];
				}

				for(int i = ctrl_pts_x[4]+1;i < W-1;i++)
				{
					float weight2 = (float)(i-ctrl_pts_x[4])/right_pts_space_num;
					float weight1 = 1.0-weight2;
					coord_ptr[((H-1)*W+i)*2+0] = weight1*pt1[0]+weight2*pt2[0];
					coord_ptr[((H-1)*W+i)*2+1] = weight1*pt1[1]+weight2*pt2[1];
				}

				float pt3[2] = {coord_ptr[0],coord_ptr[1]};
				float pt4[2] = {coord_ptr[(W-1)*2],coord_ptr[(W-1)*2+1]};

				for(int i = 1;i < W-1;i++)
				{
					float weight2 = (float)i/(W-1);
					float weight1 = 1.0-weight2;
					coord_ptr[i*2+0] = weight1*pt3[0]+weight2*pt4[0];
					coord_ptr[i*2+1] = weight1*pt3[1]+weight2*pt4[1];
				}
			}
			else
			{
				float pt0[2] = {coord_ptr[(0)*2],coord_ptr[(0)*2+1]};
				float pt1[2] = {coord_ptr[(ctrl_pts_x[4])*2],coord_ptr[(ctrl_pts_x[4])*2+1]};
				float pt2[2] = {coord_ptr[(W-1)*2],coord_ptr[(W-1)*2+1]};

				for(int i = 1;i < ctrl_pts_x[4];i++)
				{
					float weight2 = (float)i/left_pts_space_num;
					float weight1 = 1.0-weight2;
					coord_ptr[(i)*2+0] = weight1*pt0[0]+weight2*pt1[0];
					coord_ptr[(i)*2+1] = weight1*pt0[1]+weight2*pt1[1];
				}

				for(int i = ctrl_pts_x[4]+1;i < W-1;i++)
				{
					float weight2 = (float)(i-ctrl_pts_x[4])/right_pts_space_num;
					float weight1 = 1.0-weight2;
					coord_ptr[(i)*2+0] = weight1*pt1[0]+weight2*pt2[0];
					coord_ptr[(i)*2+1] = weight1*pt1[1]+weight2*pt2[1];
				}

				float pt3[2] = {coord_ptr[((H-1)*W)*2+0],coord_ptr[((H-1)*W)*2+1]};
				float pt4[2] = {coord_ptr[((H-1)*W+W-1)*2+0],coord_ptr[((H-1)*W+W-1)*2+1]};

				for(int i = 1;i < W-1;i++)
				{
					float weight2 = (float)i/(W-1);
					float weight1 = 1.0-weight2;
					coord_ptr[((H-1)*W+i)*2+0] = weight1*pt3[0]+weight2*pt4[0];
					coord_ptr[((H-1)*W+i)*2+1] = weight1*pt3[1]+weight2*pt4[1];
				}
			}
		}

	};
}

#endif