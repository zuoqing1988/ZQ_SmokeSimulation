#include "windows.h"
#include <GL/gl.h>
#include <GL/glut.h>
#include <time.h>
#include "ZQ_SmokeEditingFigGrid3D.h"


using namespace ZQ_SmokeEditing3D;
using namespace ZQ;

BaseFigGrid3D::BaseFigGrid3D()
{
	bar = 0;
	focused = false;
	move_step = 5;
}

BaseFigGrid3D::~BaseFigGrid3D()
{
	if(bar)
	{
		TwDeleteBar(bar);
		bar = 0;
	}
}


void BaseFigGrid3D::SetFocused(bool b)
{
	if(b)
	{
		if(bar == 0)
			InitTweakBar();
	}
	else
	{
		if(bar)
			DeleteTweakBar();
	}
	focused = b;
}

void BaseFigGrid3D::InitTweakBar()
{
	_createTweakBar();

	_addGlobalTweakBar();

	_addLocalTweakBar();
}

void BaseFigGrid3D::_addGlobalTweakBar()
{
	if(bar != 0)
	{	
		TwAddVarRW(bar,"move_step",TW_TYPE_FLOAT,&move_step,"label='move step' min=0.1 max=10 step=0.1");
	}
}


DeformFigGrid3D::DeformFigGrid3D()
{
	const float COLOR_FOR_GRID[3] = {0,0.8,0};
	const float COLOR_FOR_FIXED_POINT[3] = {0.8,0,0.8};
	const float COLOR_FOR_KEY_POINT[3] = {0.8,0,0.8};
	const float COLOR_FOR_CONTROL_POINT[3] = {1,0,0};
	memcpy(color_for_grid,COLOR_FOR_GRID,sizeof(float)*3);
	memcpy(color_for_fixed_pt,COLOR_FOR_FIXED_POINT,sizeof(float)*3);
	memcpy(color_for_key_pt,COLOR_FOR_KEY_POINT,sizeof(float)*3);
	memcpy(color_for_ctrl_pt,COLOR_FOR_CONTROL_POINT,sizeof(float)*3);

	line_width_for_grid = 1;
	point_size_for_fixed_pt = 5;
	point_size_for_key_pt = 10;
	point_size_for_ctrl_pt = 10;
	draw_grid = true;
	border_size = 0;

	edit_mode = false;

	width = 0;
	height = 0;
	line_weight = 0;
	angle_weight = 1;
	distance_weight = 0;
	distance = 1;
	arap = false;
	cur_arap = false;
	neighbor_type = ZQ_GridDeformation3DOptions::NEIGHBOR_6;

	total_iteration = 1000;
	cur_iteration = 0;
	per_frame_iteration = 60;
	with_good_init = false;
	grid_scale = 10;

	move_step = 0.5;
	xloop = false;
	ctrl_pts_num = 0;
	ctrl_id = -1;
}

DeformFigGrid3D::~DeformFigGrid3D()
{

}


bool DeformFigGrid3D::Export(const char* export_fold, const char* script_file, int W, int H, int D)
{
	ZQ_DImage3D<float> tmp_coords(coords);
	float*& coord_ptr = tmp_coords.data();
	int npixels = tmp_coords.npixels();
	tmp_coords.Multiplywith(grid_scale);

	char buf[2000];
	char mapfilename[FILENAME_MAX];
	char objfilename[FILENAME_MAX];

	sprintf(mapfilename,"%s\\map.di3",export_fold);
	if(!tmp_coords.saveImage(mapfilename))
	{
		printf("failed to save file %s\n",mapfilename);
		return false;
	}

	bool has_obj = false;
	for(int i = 0;i < nouseful_flag.npixels();i++)
	{
		if(nouseful_flag.data()[i])
		{
			has_obj = true;
			break;
		}
	}
	if(has_obj)
	{
		ZQ_DImage3D<float> tmp_flag(width,height,depth);
		for(int i = 0;i < nouseful_flag.npixels();i++)
			tmp_flag.data()[i] = nouseful_flag.data()[i] ? 1 : 0;
		sprintf(objfilename,"%s\\obj_mask.di3",export_fold);
		if(!tmp_flag.saveImage(objfilename))
		{
			printf("failed to save file %s\n",objfilename);
			return false;
		}
	}

	
	return true;
}


void DeformFigGrid3D::_move_all_pts(float x, float y, float z)
{
	int npixles = coords.npixels();
	float*& coord_ptr = coords.data();
	for(int i = 0;i < npixles;i++)
	{
		coord_ptr[i*3+0] += x;
		coord_ptr[i*3+1] += y;
		coord_ptr[i*3+2] += z;
	}
}

float DeformFigGrid3D::_cal_dis_point_with_ray(const float* pt, float scale, const float* ray_ori, const float* ray_dir)
{
	float dir_len2 = ray_dir[0]*ray_dir[0]+ray_dir[1]*ray_dir[1]+ray_dir[2]*ray_dir[2];
	if(dir_len2 == 0)
		return 0;
	float v1[3] = {pt[0]*scale-ray_ori[0],pt[1]*scale-ray_ori[1],pt[2]*scale-ray_ori[2]};
	float v1_len2 = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
	if(v1_len2 == 0)
		return 0;
	float cos_theta = (v1[0]*ray_dir[0]+v1[1]*ray_dir[1]+v1[2]*ray_dir[2])/(sqrt(v1_len2*dir_len2));
	float sin_theta = sqrt(1-cos_theta*cos_theta);
	return sqrt(v1_len2)*sin_theta;
}

void DeformFigGrid3D::_reset_deform_iteration()
{
	cur_iteration = total_iteration;
	with_good_init = false;
}


void DeformFigGrid3D::Draw()
{
	if(draw_grid)
	{
		_draw_grid();
	}
	if(focused && draw_grid)
	{
		_draw_control_pts();
	}

	_draw_user_defined();
}

void DeformFigGrid3D::UpdatePerFrame()
{
	if(cur_iteration <= 0)
		return;

	clock_t t1 = clock();
	switch(opt.methodType)
	{
	case ZQ_GridDeformation3DOptions::METHOD_LAPLACIAN:
		{
			ZQ_DImage3D<float> init_coord(coords);
			deform._deformation_laplacian(init_coord.data(),coords.data(),per_frame_iteration);
			with_good_init = true;
			cur_iteration -= per_frame_iteration;
		}
		break;
	case ZQ_GridDeformation3DOptions::METHOD_LAPLACIAN_XLOOP:
		{
			ZQ_DImage3D<float> init_coord(coords);
			deform._deformation_laplacian_XLOOP(init_coord.data(),coords.data(),per_frame_iteration);
			with_good_init = true;
			cur_iteration -= per_frame_iteration;
		}
		break;
	case ZQ_GridDeformation3DOptions::METHOD_ARAP_VERT_AS_CENTER:
		{
			ZQ_DImage3D<float> init_coord(coords);
			deform._deformation_ARAP_VERT(init_coord.data(),coords.data(),1,per_frame_iteration);
			cur_iteration -= per_frame_iteration;
		}
		break;
	case ZQ_GridDeformation3DOptions::METHOD_ARAP_VERT_AS_CENTER_XLOOP:
		{
			ZQ_DImage3D<float> init_coord(coords);
			deform._deformation_ARAP_VERT_XLOOP(init_coord.data(),coords.data(),1,per_frame_iteration);
			cur_iteration -= per_frame_iteration;
		}
		break;
	}
	clock_t t2 = clock();
	printf("solve: %.3f ms\n",(double)(t2-t1));

}

void DeformFigGrid3D::OnMouseButton(int button_id, int state, const float ray_ori[3], const float ray_dir[3])
{
	switch(button_id)
	{
	case GLUT_LEFT_BUTTON:
		{
			if(edit_mode)
			{
				if(state == GLUT_DOWN)
				{
					_select_ctrl_pt(ray_ori,ray_dir);
				}
			}
		}
		break;
	}
}

void DeformFigGrid3D::OnKeyboard(int key_id, const float ray_ori[3], const float ray_dir[3])
{
	switch(key_id)
	{
	case 'w':case 'W':
		{
			if(edit_mode)
				_move_ctrl_pt(0,move_step,0);
			else
				_move_all_pts(0,move_step,0);
		}
		break;
	case 's':case 'S':
		{
			if(edit_mode)
				_move_ctrl_pt(0,-move_step,0);
			else
				_move_all_pts(0,-move_step,0);
		}
		break;
	case 'a':case 'A':
		{
			if(edit_mode)
				_move_ctrl_pt(-move_step,0,0);
			else
				_move_all_pts(-move_step,0,0);
		}
		break;
	case 'd':case 'D':
		{
			if(edit_mode)
				_move_ctrl_pt(move_step,0,0);
			else
				_move_all_pts(move_step,0,0);
		}
		break;
	case 'z':case 'Z':
		{
			if(edit_mode)
				_move_ctrl_pt(0,0,move_step);
			else
				_move_all_pts(0,0,move_step);
		}
		break;
	case 'x':case 'X':
		{
			if(edit_mode)
				_move_ctrl_pt(0,0,-move_step);
			else
				_move_all_pts(0,0,-move_step);
		}
		break;
	}
}


void DeformFigGrid3D::GlobalMove(float x, float y, float z)
{
	int npixels = coords.npixels();
	float*& coord_ptr = coords.data();
	for(int i = 0;i < npixels;i++)
	{
		coord_ptr[i*3+0] += x/grid_scale;
		coord_ptr[i*3+1] += y/grid_scale;
		coord_ptr[i*3+2] += z/grid_scale;
	}

}

void DeformFigGrid3D::OnApply(void* clientData)
{
	DeformFigGrid3D* g = (DeformFigGrid3D*)clientData;
	g->_apply();
}

void DeformFigGrid3D::_apply()
{
	bool rebuild_flag = false;
	if(!coords.matchDimension(width,height,depth,3))
	{
		_init_ctrl_pts();
		coords.allocate(width,height,depth,3);
		nouseful_flag.allocate(width,height,depth);
		fixed_flag.allocate(width,height,depth);

		_init_coords();
		_init_fixed_pts();
		rebuild_flag = true;
	}

	/*if(opt.line_weight != line_weight)
	{
		opt.line_weight = line_weight;
		rebuild_flag = true;
	}
	if(opt.angle_weight != angle_weight)
	{
		opt.angle_weight = angle_weight;
		rebuild_flag = true;
	}
	if(opt.distance_weight != distance_weight)
	{
		opt.distance_weight = distance_weight;
		rebuild_flag = true;
	}*/
	if(opt.distance != distance)
	{
		opt.distance = distance;
		rebuild_flag = true;
	}
	if(cur_arap != arap)
	{
		cur_arap = arap;
		rebuild_flag = true;
	}
	if(opt.neighborType != neighbor_type)
	{
		opt.neighborType = neighbor_type;
		rebuild_flag = true;
	}

	if(rebuild_flag)
	{
		_rebuild_matrix();
		_reset_deform_iteration();
	}
}

void DeformFigGrid3D::_rebuild_matrix()
{
	if(arap)
	{
		opt.methodType = xloop ? ZQ_GridDeformation3DOptions::METHOD_ARAP_VERT_AS_CENTER_XLOOP :
			ZQ_GridDeformation3DOptions::METHOD_ARAP_VERT_AS_CENTER;
	}
	else
	{
		opt.methodType = xloop ? ZQ_GridDeformation3DOptions::METHOD_LAPLACIAN_XLOOP :
			ZQ_GridDeformation3DOptions::METHOD_LAPLACIAN;
	}


	deform.BuildMatrix(width,height,depth,nouseful_flag.data(),fixed_flag.data(),opt);
}

void DeformFigGrid3D::_draw_grid()
{
	int H = coords.height();
	int W = coords.width();
	int D = coords.depth();
	float*& coord_ptr = coords.data();
	bool*& nouseful_ptr = nouseful_flag.data();


	glLineWidth(line_width_for_grid);
	glColor3f(color_for_grid[0],color_for_grid[1],color_for_grid[2]);
	glBegin(GL_LINES);

	if(!xloop)
	{
		for(int k = 0;k < D;k++)
		{
			for(int j = 0;j < H;j++)
			{
				for(int i = 0;i < W-1;i++)
				{
					int offset = k*H*W+j*W+i;
					int offset_2 = k*H*W+j*W+i+1;
					if(nouseful_ptr[offset] || nouseful_ptr[offset_2])
						continue;
					glVertex3f(coord_ptr[offset*3+0]*grid_scale,coord_ptr[offset*3+1]*grid_scale,coord_ptr[offset*3+2]*grid_scale);
					glVertex3f(coord_ptr[offset_2*3+0]*grid_scale,coord_ptr[offset_2*3+1]*grid_scale,coord_ptr[offset_2*3+2]*grid_scale);
				}
			}
		}
	}
	else
	{
		for(int k = 0;k < D;k++)
		{
			for(int j = 0;j < H;j++)
			{
				for(int i = 0;i < W;i++)
				{
					int offset = k*H*W+j*W+i;
					int offset_2 = k*H*W+j*W+(i+1)%W;
					if(nouseful_ptr[offset] || nouseful_ptr[offset_2])
						continue;
					glVertex3f(coord_ptr[offset*3+0]*grid_scale,coord_ptr[offset*3+1]*grid_scale,coord_ptr[offset*3+2]*grid_scale);
					glVertex3f(coord_ptr[offset_2*3+0]*grid_scale,coord_ptr[offset_2*3+1]*grid_scale,coord_ptr[offset_2*3+2]*grid_scale);
				}
			}
		}
	}

	for(int k = 0;k < D;k++)
	{
		for(int j = 0;j < H-1;j++)
		{
			for(int i = 0;i < W;i++)
			{
				int offset = k*H*W+j*W+i;
				int offset_2 = k*H*W+(j+1)*W+i;
				if(nouseful_ptr[offset] || nouseful_ptr[offset_2])
					continue;
				glVertex3f(coord_ptr[offset*3+0]*grid_scale,coord_ptr[offset*3+1]*grid_scale,coord_ptr[offset*3+2]*grid_scale);
				glVertex3f(coord_ptr[offset_2*3+0]*grid_scale,coord_ptr[offset_2*3+1]*grid_scale,coord_ptr[offset_2*3+2]*grid_scale);
			}
		}
	}

	for(int k = 0;k < D-1;k++)
	{
		for(int j = 0;j < H;j++)
		{
			for(int i = 0;i < W;i++)
			{
				int offset = k*H*W+j*W+i;
				int offset_2 = (k+1)*H*W+j*W+i;
				if(nouseful_ptr[offset] || nouseful_ptr[offset_2])
					continue;
				glVertex3f(coord_ptr[offset*3+0]*grid_scale,coord_ptr[offset*3+1]*grid_scale,coord_ptr[offset*3+2]*grid_scale);
				glVertex3f(coord_ptr[offset_2*3+0]*grid_scale,coord_ptr[offset_2*3+1]*grid_scale,coord_ptr[offset_2*3+2]*grid_scale);
			}
		}
	}

	glEnd();

}


void DeformFigGrid3D::_draw_control_pts()
{
	int H = coords.height();
	int W = coords.width();
	int D = coords.depth();
	float*& coord_ptr = coords.data();
	bool*& nouseful_ptr = nouseful_flag.data();

	float scale = grid_scale;

	bool*& fixed_ptr = fixed_flag.data();
	glPointSize(point_size_for_fixed_pt);
	glColor3f(color_for_fixed_pt[0],color_for_fixed_pt[1],color_for_fixed_pt[2]);
	glBegin(GL_POINTS);
	for(int i = 0;i < W*H*D;i++)
	{
		if(!nouseful_ptr[i] && fixed_ptr[i])
			glVertex3f(coord_ptr[i*3+0]*scale,coord_ptr[i*3+1]*scale,coord_ptr[i*3+2]*scale);
	}
	glEnd();

	if(edit_mode)
	{
		//corner points
		glPointSize(point_size_for_key_pt);
		glColor3f(color_for_key_pt[0],color_for_key_pt[1],color_for_key_pt[2]);
		glBegin(GL_POINTS);
		for(int ii = 0;ii < ctrl_pts_num;ii++)
		{
			int offset = ctrl_pts_z[ii]*H*W+ctrl_pts_y[ii]*W+ctrl_pts_x[ii];
			glVertex3f(coord_ptr[offset*3+0]*scale,coord_ptr[offset*3+1]*scale,coord_ptr[offset*3+2]*scale);
		}
		glEnd();

		if(ctrl_id >= 0)
		{
			glPointSize(point_size_for_ctrl_pt);
			glColor3f(color_for_ctrl_pt[0],color_for_ctrl_pt[1],color_for_ctrl_pt[2]);
			glBegin(GL_POINTS);
			{
				int offset = ctrl_pts_z[ctrl_id]*H*W+ctrl_pts_y[ctrl_id]*W+ctrl_pts_x[ctrl_id];
				glVertex3f(coord_ptr[offset*3+0]*scale,coord_ptr[offset*3+1]*scale,coord_ptr[offset*3+2]*scale);
			}
			glEnd();
		}
	}
}

void DeformFigGrid3D::_select_ctrl_pt(const float ray_ori[3], const float ray_dir[3])
{
	if(ctrl_pts_num <= 0)
	{
		ctrl_id = -1;
		return;
	}

	int W = coords.width();
	int H = coords.height();
	int D = coords.depth();
	float*& coord_ptr = coords.data();

	float scale = grid_scale;
	int s_id = 0;
	int offset = ctrl_pts_z[0]*H*W+ctrl_pts_y[0]*W+ctrl_pts_x[0];
	float min_dis = _cal_dis_point_with_ray(coord_ptr+offset*3,scale,ray_ori,ray_dir);
	
	for(int kk = 1;kk < ctrl_pts_num;kk++)
	{
		offset = ctrl_pts_z[kk]*H*W+ctrl_pts_y[kk]*W+ctrl_pts_x[kk];
		float cur_dis = _cal_dis_point_with_ray(coord_ptr+offset*3,scale,ray_ori,ray_dir);
		if(cur_dis < min_dis)
		{
			min_dis = cur_dis;
			s_id = kk;
		}
	}

	ctrl_id = s_id;
}


void DeformFigGrid3D::_move_ctrl_pt(float tran_x, float tran_y, float tran_z)
{
	if(ctrl_id < 0)
		return;

	int W = coords.width();
	int H = coords.height();
	int D = coords.depth();
	float*& coord_ptr = coords.data();

	int offset = ctrl_pts_z[ctrl_id]*H*W + ctrl_pts_y[ctrl_id]*W + ctrl_pts_x[ctrl_id];
	coord_ptr[offset*3+0] += tran_x;
	coord_ptr[offset*3+1] += tran_y;
	coord_ptr[offset*3+2] += tran_z;
	_recompute_fixed();
	_reset_deform_iteration();
}


void DeformFigGrid3D::_addGlobalTweakBar()
{
	const int neighbor_type_num = 2;
	const TwEnumVal neighbor_typeEnum[neighbor_type_num] =
	{
		{ZQ_GridDeformation3DOptions::NEIGHBOR_6, "neighbor 6"},
		{ZQ_GridDeformation3DOptions::NEIGHBOR_26, "neighbor 26"}
	};

	if(bar != 0)
	{
		TwAddVarRW(bar,"width",TW_TYPE_INT32,&width,"label='width' min=4 max=100 group='Group2'");
		TwAddVarRW(bar,"height",TW_TYPE_INT32,&height,"label='height' min=4 max=100 group='Group2'");
		TwAddVarRW(bar,"depth",TW_TYPE_INT32,&depth,"label='depth' min=4 max=100 group='Group2'");
		TwAddVarRW(bar,"distance_scale",TW_TYPE_FLOAT,&distance,"label='scale' min=0 max=100 group='Group2'");
		TwAddVarRW(bar,"ARAP",TW_TYPE_BOOL8,&arap,"label='ARAP' group='Group2'");
		TwType neighborTWType = TwDefineEnum("neighborType", neighbor_typeEnum, neighbor_type_num);
		TwAddVarRW(bar, "neighborType", neighborTWType, &neighbor_type, "label='neighbor type' group='Group2'");
		TwAddButton(bar,"Apply",(TwButtonCallback)OnApply,this,"label='Apply' group='Group2'");

		TwAddVarRW(bar,"edit_mode",TW_TYPE_BOOL8,&edit_mode,"label='edit' group='Group3'");
		TwAddVarRW(bar,"draw_grid",TW_TYPE_BOOL8,&draw_grid,"label='draw grid' group='Group3'");
		TwAddVarRW(bar,"grid_scale",TW_TYPE_FLOAT,&grid_scale,"label='grid_scale' min=1 max=50 group='Group3'");
		TwAddVarRW(bar,"move_step",TW_TYPE_FLOAT,&move_step,"label='move step' min=0.1 max=10 step=0.1 group='Group3'");
		TwAddVarRW(bar,"total_it",TW_TYPE_INT32,&total_iteration,"label='total_it' min=1000 max=100000 step=1000 group='Group3'");
		TwAddVarRW(bar,"perframe_it",TW_TYPE_INT32,&per_frame_iteration,"label='perframe_it' min=20 max=2000 step=20 group='Group3'");	
	}
}

