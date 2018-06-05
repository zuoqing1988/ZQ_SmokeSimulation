#ifndef _ZQ_SMOKE_EDITING_ZIPPER_FIGGRID_H_
#define _ZQ_SMOKE_EDITING_ZIPPER_FIGGRID_H_

namespace ZQ_SmokeEditing
{
	class ZipperDeformGrid: public DeformFigGrid
	{
	public:
		ZipperDeformGrid()
		{
			width = 33;
			height = 33;
			center_width = width/2;
			_apply();
			xloop = false;
		}

		~ZipperDeformGrid(){}

	private:

		int center_width;
		
		void _createTweakBar()
		{
			if(bar == 0)
			{
				const char* name = "ZipperDeformGrid";
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
			ZipperDeformGrid* g = (ZipperDeformGrid*)clientData;
			g->_centerInc();
		}

		static void TW_CALL OnCenterDEC(void* clientData)
		{
			ZipperDeformGrid* g = (ZipperDeformGrid*)clientData;
			g->_centerDec();
		}

		void _centerDec()
		{
			if(ctrl_pts_y[7] > 1)
			{
				nouseful_flag.data()[ctrl_pts_y[7]*width+center_width] = true;
				ctrl_pts_y[7] --;
				_rebuild_matrix();
				_recompute_fixed();
				_reset_deform_iteration();
			}
		}

		void _centerInc()
		{
			int H = coords.height();
			if(ctrl_pts_y[7] < H-1)
			{
				ctrl_pts_y[7] ++;
				nouseful_flag.data()[ctrl_pts_y[7]*width+center_width] = false;
				_rebuild_matrix();
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
			center_width = width/2;
			ctrl_pts_num = 8;
			ctrl_pts_x.clear();
			ctrl_pts_y.clear();
			ctrl_pts_x.push_back(0);
			ctrl_pts_y.push_back(height-1);
			ctrl_pts_x.push_back(center_width-1);
			ctrl_pts_y.push_back(height-1);
			ctrl_pts_x.push_back(center_width+1);
			ctrl_pts_y.push_back(height-1);
			ctrl_pts_x.push_back(width-1);
			ctrl_pts_y.push_back(height-1);
			ctrl_pts_x.push_back(0);
			ctrl_pts_y.push_back(0);
			ctrl_pts_x.push_back(center_width);
			ctrl_pts_y.push_back(0);
			ctrl_pts_x.push_back(width-1);
			ctrl_pts_y.push_back(0);
			ctrl_pts_x.push_back(center_width);
			ctrl_pts_y.push_back(height-1);
			

			ctrl_id = 0;
		}

		void _init_fixed_pts()
		{
			for(int i = 0;i < width;i++)
			{
				fixed_flag.data()[0*width+i] = true;
				fixed_flag.data()[(height-1)*width+i] = true;
			}
			for(int i = 0;i < height;i++)
			{
				fixed_flag.data()[i*width+center_width] = true;
			}
		}


		void _recompute_fixed()
		{
			int W = coords.width();
			int H = coords.height();
			float*& coord_ptr = coords.data();

			float pt0[2] = {coord_ptr[((height-1)*W+0)*2+0],coord_ptr[((height-1)*W+0)*2+1]};
			float pt1[2] = {coord_ptr[((height-1)*W+center_width-1)*2+0],coord_ptr[((height-1)*W+center_width-1)*2+1]};
			for(int i = 1;i < center_width-1;i++)
			{
				float weight2 = (float)i/(center_width-1);
				float weight1 = 1.0 - weight2;
				coord_ptr[((height-1)*W+i)*2+0] = weight1*pt0[0] + weight2*pt1[0];
				coord_ptr[((height-1)*W+i)*2+1] = weight1*pt0[1] + weight2*pt1[1];
			}

			float pt2[2] = {coord_ptr[((height-1)*W+center_width+1)*2+0],coord_ptr[((height-1)*W+center_width+1)*2+1]};
			float pt3[2] = {coord_ptr[((height-1)*W+width-1)*2+0],coord_ptr[((height-1)*W+width-1)*2+1]};
			for(int i = center_width+2;i < width-1;i++)
			{
				float weight2 = (float)(i-center_width-1)/(width-1-center_width-1);
				float weight1 = 1.0 - weight2;
				coord_ptr[((height-1)*W+i)*2+0] = weight1*pt2[0] + weight2*pt3[0];
				coord_ptr[((height-1)*W+i)*2+1] = weight1*pt2[1] + weight2*pt3[1];
			}

			float pt4[2] = {coord_ptr[0],coord_ptr[1]};
			float pt5[2] = {coord_ptr[(center_width)*2+0],coord_ptr[(center_width)*2+1]};
			for(int i = 1;i < center_width;i++)
			{
				float weight2 = (float)i/(center_width);
				float weight1 = 1.0 - weight2;
				coord_ptr[i*2+0] = weight1*pt4[0] + weight2*pt5[0];
				coord_ptr[i*2+1] = weight1*pt4[1] + weight2*pt5[1];
			}

			float pt6[2] = {coord_ptr[(center_width)*2+0],coord_ptr[(center_width)*2+1]};
			float pt7[2] = {coord_ptr[(width-1)*2+0],coord_ptr[(width-1)*2+1]};
			for(int i = center_width+1;i < width-1;i++)
			{
				float weight2 = (float)(i-center_width)/(width-1-center_width);
				float weight1 = 1.0 - weight2;
				coord_ptr[i*2+0] = weight1*pt6[0] + weight2*pt7[0];
				coord_ptr[i*2+1] = weight1*pt6[1] + weight2*pt7[1];
			}

			float pt8[2] = {coord_ptr[(center_width)*2+0],coord_ptr[(center_width)*2+1]};
			float pt9[2] = {coord_ptr[(ctrl_pts_y[7]*W+center_width)*2+0],coord_ptr[(ctrl_pts_y[7]*W+center_width)*2+1]};
			for(int i = 1;i < ctrl_pts_y[7];i++)
			{
				float weight2 = (float)(i)/(ctrl_pts_y[7]);
				float weight1 = 1.0 - weight2;
				coord_ptr[(i*W+center_width)*2+0] = weight1*pt8[0] + weight2*pt9[0];
				coord_ptr[(i*W+center_width)*2+1] = weight1*pt8[1] + weight2*pt9[1];
			}
		}
	};
}

#endif