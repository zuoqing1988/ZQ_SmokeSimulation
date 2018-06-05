#ifndef _ZQ_SMOKE_EDITING_FREE_DEFORM_FIGGRID_H_
#define _ZQ_SMOKE_EDITING_FREE_DEFORM_FIGGRID_H_

#include "ZQ_SmokeEditingFigGrid.h"

namespace ZQ_SmokeEditing
{
	class FreeDeformFigGird : public DeformFigGrid
	{
	public:
		FreeDeformFigGird()
		{
			width = 17;
			height = 17;
			_apply();
			erase_useful_radius = 0.5;
		}
		~FreeDeformFigGird(){}

	private:
		bool mode_erase_useful;
		float erase_useful_radius;
		float erase_useful_coord_x;
		float erase_useful_coord_y;

	public:
		void OnMouseButton(int button_id, int state, int mx, int my)
		{
			switch(button_id)
			{
			case GLUT_LEFT_BUTTON:
				{
					if(edit_mode)
					{
						if(state == GLUT_DOWN)
						{
							if(mode_erase_useful)
							{
								erase_useful_coord_x = mx;
								erase_useful_coord_y = my;
								_erase_useful(mx,my);
							}
							else
								_select_ctrl_pt(mx,my);
						}
					}
				}
				break;
			}
		}

		void OnMouseMove(int mx, int my)
		{
			if(mode_erase_useful)
			{
				erase_useful_coord_x = mx;
				erase_useful_coord_y = my;
				_erase_useful(mx,my);
			}
		}

		void OnMousePassiveMove(int mx, int my)
		{
			if(mode_erase_useful)
			{
				erase_useful_coord_x = mx;
				erase_useful_coord_y = my;
			}
		}

		void OnKeyboard(int key_id, int mx, int my)
		{
			switch(key_id)
			{
			case 'w':case 'W':
				{
					if(edit_mode)
						_move_ctrl_pt(0,move_step);
					else
						_move_all_pts(0,move_step);
				}
				break;
			case 's':case 'S':
				{
					if(edit_mode)
						_move_ctrl_pt(0,-move_step);
					else
						_move_all_pts(0,-move_step);
				}
				break;
			case 'a':case 'A':
				{
					if(edit_mode)
						_move_ctrl_pt(-move_step,0);
					else
						_move_all_pts(-move_step,0);
				}
				break;
			case 'd':case 'D':
				{
					if(edit_mode)
						_move_ctrl_pt(move_step,0);
					else
						_move_all_pts(move_step,0);
				}
				break;
			case 'q':case 'Q':
				{
					_rotate(rot_step);
				}
				break;
			case 'e':case 'E':
				{
					_rotate(-rot_step);
				}
				break;
			case 'f':case 'F':
				{
					_toggle_fixed(mx,my);
				}
				break;
			case 't':case 'T':
				{
					_make_nouseful(mx,my);
				}
				break;
			}
		}

	private:

		void _draw_user_defined(float z_shift)
		{
			if(mode_erase_useful)
			{
				const float m_pi = 3.1415926535;
				float real_radius = grid_scale*erase_useful_radius;
				int circle_N = 32;
				glColor3f(1,1,1);
				glBegin(GL_LINE_LOOP);
				for(int i=0; i < circle_N; i++)
					glVertex3f(erase_useful_coord_x+real_radius*cos(2*m_pi/circle_N*i), erase_useful_coord_y+real_radius*sin(2*m_pi/circle_N*i),z_shift);
				glEnd();
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

		void _init_ctrl_pts()
		{
			ctrl_pts_num = 4;
			ctrl_pts_x.clear();
			ctrl_pts_y.clear();

			ctrl_pts_x.push_back(0);
			ctrl_pts_y.push_back(0);
			ctrl_pts_x.push_back(0);
			ctrl_pts_y.push_back(height-1);
			ctrl_pts_x.push_back(width-1);
			ctrl_pts_y.push_back(0);
			ctrl_pts_x.push_back(width-1);
			ctrl_pts_y.push_back(height-1);

			ctrl_id = 0;
		}

		void _init_fixed_pts()
		{
			fixed_flag.data()[0*width+0] = true;
			fixed_flag.data()[(height-1)*width+0] = true;
			fixed_flag.data()[0*width+width-1] = true;
			fixed_flag.data()[(height-1)*width+width-1] = true;		
		}
		
		void _recompute_fixed(){}

		void _createTweakBar()
		{
			if(bar == 0)
			{
				const char* name = "FreeDeformGrid";
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
			TwAddVarRW(bar,"mode_erase_useful",TW_TYPE_BOOL8,&mode_erase_useful,"label='erase useful'");
			TwAddVarRW(bar,"erase_useful_radius",TW_TYPE_FLOAT,&erase_useful_radius,"label='erase radius' min=0.1 max=1.5 step=0.1");
		}


		bool _select_any_pt(float mx, float my, int& x, int& y)
		{
			int W = coords.width();
			int H = coords.height();
			float*& ptr = coords.data();
			bool*& nouseful_ptr = nouseful_flag.data();
			int cur_find_x = -1;
			int cur_find_y = -1;
			float min_dis2 = 0;
			bool has_find_first = false;
			for(int i = 0;i < H;i++)
			{
				for(int j = 0;j < W;j++)
				{
					if(!nouseful_ptr[i*W+j])
					{
						float cur_dis2 = (mx-ptr[(i*W+j)*2+0]*grid_scale)*(mx-ptr[(i*W+j)*2+0]*grid_scale)
							+ (my-ptr[(i*W+j)*2+1]*grid_scale)*(my-ptr[(i*W+j)*2+1]*grid_scale);
						if(!has_find_first)
						{
							cur_find_x = j;
							cur_find_y = i;
							has_find_first = true;
							min_dis2 = cur_dis2;
						}
						else
						{
							if(cur_dis2 < min_dis2)
							{
								cur_find_x = j;
								cur_find_y = i;
								min_dis2 = cur_dis2;
							}
						}
					}
				}
			}

			if(has_find_first)
			{
				x = cur_find_x;
				y = cur_find_y;
				return true;
			}
			else
				return false;
		}

		void _make_nouseful(float mx, float my)
		{
			int ix,iy;
			if(_select_any_pt(mx,my,ix,iy))
			{
				int W = coords.width();
				bool cur_b = nouseful_flag.data()[iy*W+ix];
				if(!cur_b)
				{
					nouseful_flag.data()[iy*W+ix] = true;
					bool cur_f = fixed_flag.data()[iy*W+ix];
					int del_pp = -1;
					if(cur_f)
					{
						for(int pp = 0;pp < ctrl_pts_num;pp++)
						{
							if(ctrl_pts_x[pp] == ix && ctrl_pts_y[pp] == iy)
							{
								ctrl_pts_x.erase(ctrl_pts_x.begin()+pp);
								ctrl_pts_y.erase(ctrl_pts_y.begin()+pp);
								del_pp = pp;
								break;
							}
						}
						ctrl_pts_num--;
						if(ctrl_pts_num <= 0)
							ctrl_id = -1;
						else
						{
							if(ctrl_id > del_pp)
								ctrl_id --;
							else if(ctrl_id == del_pp)
								ctrl_id = 0;
						}
					}

					_rebuild_matrix();
					_reset_deform_iteration();
				}
			}
		}

		void _toggle_fixed(float mx, float my)
		{
			int ix,iy;
			if(_select_any_pt(mx,my,ix,iy))
			{
				int W = coords.width();
				bool cur_b = nouseful_flag.data()[iy*W+ix];
				if(!cur_b)
				{
					bool cur_f = fixed_flag.data()[iy*W+ix];
					if(cur_f)
					{
						int del_pp = -1;
						for(int pp = 0;pp < ctrl_pts_num;pp++)
						{
							if(ctrl_pts_x[pp] == ix && ctrl_pts_y[pp] == iy)
							{
								ctrl_pts_x.erase(ctrl_pts_x.begin()+pp);
								ctrl_pts_y.erase(ctrl_pts_y.begin()+pp);
								del_pp = pp;
								break;
							}
						}
						ctrl_pts_num--;
						if(ctrl_pts_num <= 0)
							ctrl_id = -1;
						else
						{
							if(ctrl_id > del_pp)
								ctrl_id --;
							else if(ctrl_id == del_pp)
								ctrl_id = 0;
						}
						_reset_deform_iteration();
					}
					else
					{
						ctrl_pts_x.push_back(ix);
						ctrl_pts_y.push_back(iy);
						ctrl_pts_num++;
						ctrl_id = ctrl_pts_num-1;
					}
					fixed_flag.data()[iy*W+ix] = !cur_f;
					_rebuild_matrix();
				}
			}
		}

		void _erase_useful(float mx, float my)
		{
			int W = coords.width();
			int H = coords.height();
			float*& ptr = coords.data();
			bool*& nouseful_ptr = nouseful_flag.data();
			bool*& fixed_ptr = fixed_flag.data();

			bool any_change = false;
			float thresh_dis2 = grid_scale*grid_scale*erase_useful_radius*erase_useful_radius;

			for(int i = 0;i < H;i++)
			{
				for(int j = 0;j < W;j++)
				{
					if(!nouseful_ptr[i*W+j])
					{
						float px = ptr[(i*W+j)*2+0]*grid_scale;
						float py = ptr[(i*W+j)*2+1]*grid_scale;
						float cur_dis2 = (mx-px)*(mx-px)+(my-py)*(my-py);
						if(cur_dis2 < thresh_dis2)
						{
							any_change = true;
							nouseful_flag.data()[i*W+j] = true;
							bool cur_f = fixed_ptr[i*W+j];
							if(cur_f)
							{
								int del_pp = -1;
								for(int pp = 0;pp < ctrl_pts_num;pp++)
								{
									if(ctrl_pts_x[pp] == j && ctrl_pts_y[pp] == i)
									{
										ctrl_pts_x.erase(ctrl_pts_x.begin()+pp);
										ctrl_pts_y.erase(ctrl_pts_y.begin()+pp);
										del_pp = pp;
										break;
									}
								}
								ctrl_pts_num--;
								if(ctrl_pts_num <= 0)
									ctrl_id = -1;
								else
								{
									if(ctrl_id > del_pp)
										ctrl_id --;
									else if(ctrl_id == del_pp)
										ctrl_id = 0;
								}
							}
						}
					}
				}
			}
			if(any_change)
				_rebuild_matrix();
		}
	};
}

#endif