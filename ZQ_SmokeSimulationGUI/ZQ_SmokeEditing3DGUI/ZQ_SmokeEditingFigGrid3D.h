#ifndef _ZQ_SMOKE_EDITING_FIG_GRID_3D_H_
#define _ZQ_SMOKE_EDITING_FIG_GRID_3D_H_
#pragma once

#include <GL/gl.h>
#include "AntTweakBar.h"
#include "ZQ_DoubleImage3D.h"
#include "ZQ_GridDeformation3D.h"




namespace ZQ_SmokeEditing3D
{
	class BaseFigGrid3D
	{
		enum CONST_VAL{FILENAME_LEN=200};

	public:
		BaseFigGrid3D();
		~BaseFigGrid3D();
	protected:

		TwBar* bar;
		float move_step;
		bool focused;

	public:
		bool IsFocused()const {return focused;}
		void SetFocused(bool b);

		virtual void Draw() = 0;
		virtual void UpdatePerFrame() = 0;
		virtual void OnMouseButton(int button_id, int state, const float ray_ori[3], const float ray_dir[3]) = 0;
		virtual void OnKeyboard(int key_id, const float ray_ori[3], const float ray_dir[3]) = 0;
		virtual void OnMouseMove(const float ray_ori[3], const float ray_dir[3]) = 0;
		virtual void OnMousePassiveMove(const float ray_ori[3], const float ray_dir[3]) = 0;
		virtual bool Export(const char* export_fold, const char* script_file, int width, int height, int depth) = 0;
		virtual int GetWidth() const = 0;
		virtual int GetHeight() const = 0;
		virtual int GetDepth() const = 0;
		virtual void GlobalMove(float x, float y, float z) = 0;

	private:
		void InitTweakBar();
		virtual void _createTweakBar() = 0;
		virtual void _addGlobalTweakBar();
		virtual void _addLocalTweakBar() = 0;
		void DeleteTweakBar(){if(bar) {TwDeleteBar(bar);bar= 0;} }	
	};


	class DeformFigGrid3D: public BaseFigGrid3D
	{
	public:
		DeformFigGrid3D();
		~DeformFigGrid3D();

	protected:
		float color_for_grid[3];
		float color_for_fixed_pt[3];
		float color_for_key_pt[3];
		float color_for_ctrl_pt[3];
		float line_width_for_grid;
		float point_size_for_fixed_pt;
		float point_size_for_key_pt;
		float point_size_for_ctrl_pt;
		bool draw_grid;
		int border_size;


	protected:
		bool xloop;
		int width,height,depth;
		float line_weight;
		float angle_weight;
		float distance_weight;
		float distance;
		bool arap;
		bool cur_arap;
		ZQ::ZQ_GridDeformation3DOptions::NeighborType neighbor_type;
		ZQ::ZQ_DImage3D<float> coords;
		ZQ::ZQ_DImage3D<bool> nouseful_flag;
		ZQ::ZQ_DImage3D<bool> fixed_flag;
		ZQ::ZQ_GridDeformation3DOptions opt;
		ZQ::ZQ_GridDeformation3D<float> deform;
		bool edit_mode;
		int total_iteration;
		int cur_iteration;
		int per_frame_iteration;
		bool with_good_init;
		float grid_scale;

		int ctrl_pts_num;
		std::vector<int> ctrl_pts_x;
		std::vector<int> ctrl_pts_y;
		std::vector<int> ctrl_pts_z;
		int ctrl_id;


		void _reset_deform_iteration();
		static void TW_CALL OnApply(void* clientData);
		void _apply();
		void _rebuild_matrix();

	public:
		void Draw();
		void UpdatePerFrame();
		virtual void OnMouseButton(int button_id, int state,const float ray_ori[3], const float ray_dir[3]);
		virtual void OnKeyboard(int key_id,const float ray_ori[3], const float ray_dir[3]);
		virtual void OnMouseMove(const float ray_ori[3], const float ray_dir[3]){}
		virtual void OnMousePassiveMove(const float ray_ori[3], const float ray_dir[3]){}
		bool Export(const char* export_fold, const char* script_file, int width, int height, int depth);
		int GetWidth() const {return coords.width();}
		int GetHeight() const {return coords.height();}
		int GetDepth() const {return coords.depth();}
		void GlobalMove(float x, float y, float z);

	protected:
		void _select_ctrl_pt(const float ray_ori[3], const float ray_dir[3]);
		void _move_ctrl_pt(float tran_x, float tran_y, float tran_z);
		void _move_all_pts(float x, float y, float z);
		float _cal_dis_point_with_ray(const float* pt, float scale, const float* ray_ori, const float* ray_dir);

	private:

		void _draw_grid();
		void _draw_control_pts();
		virtual void _draw_user_defined(){}


		virtual void _init_coords() = 0;
		virtual void _init_ctrl_pts() = 0;
		virtual void _init_fixed_pts() = 0;
		virtual void _recompute_fixed() = 0;

		virtual void _createTweakBar() = 0;
		virtual void _addGlobalTweakBar();
		virtual void _addLocalTweakBar(){};
	};


}


#endif