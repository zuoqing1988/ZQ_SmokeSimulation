#ifndef _ZQ_SMOKE_EDITING_FIG_GRID_H_
#define _ZQ_SMOKE_EDITING_FIG_GRID_H_

#include <GL/gl.h>
#include "ZQ_DoubleImage.h"
#include "ZQ_MergeImage.h"
#include "ZQ_GridDeformation.h"
#include "AntTweakBar.h"



namespace ZQ_SmokeEditing
{
	class BaseFigGrid
	{
		enum CONST_VAL{FILENAME_LEN=200};
		
	public:
		BaseFigGrid();
		~BaseFigGrid();
	protected:
		char fold[FILENAME_LEN];
		char prefix[FILENAME_LEN];
		char suffix[FILENAME_LEN];
		int total_frame;
		int base_id;

		TwBar* bar;
		float move_step;
		float rot_step;
		GLuint tex_id;
		bool draw_texture;
		bool has_tex;
		char tex_file[FILENAME_LEN];
		bool focused;
		float tex_alpha;

	public:
		bool IsFocused()const {return focused;}
		const char* GetFoldPtr() const {return fold;}
		const char* GetPrefixPtr() const {return prefix;}
		const char* GetSuffixPtr() const {return suffix;}
		int GetTotalFrame() const {return total_frame;}
		int GetBaseId() const {return base_id;}

		void SetFocused(bool b);

		virtual void Draw(float z_shift) = 0;
		virtual void UpdatePerFrame() = 0;
		virtual void OnMouseButton(int button_id, int state, int mx, int my) = 0;
		virtual void OnKeyboard(int key_id, int mx, int my) = 0;
		virtual void OnMouseMove(int mx, int my) = 0;
		virtual void OnMousePassiveMove(int mx, int my) = 0;
		virtual bool Export(const char* export_fold, const char* script_file, int width, int height) = 0;
		virtual int GetWidth() const = 0;
		virtual int GetHeight() const = 0;
		virtual void GlobalMove(float x, float y) = 0;

	protected:
		static void TW_CALL OnLoadTexture(void* clientData);
	private:
		void InitTweakBar();
		virtual void _createTweakBar() = 0;
		virtual void _addGlobalTweakBar();
		virtual void _addLocalTweakBar() = 0;
		void DeleteTweakBar(){if(bar) {TwDeleteBar(bar);bar= 0;} }	
		virtual void _loadTexture(bool _forceUpdate);
	};

	class MaskFigGrid: public BaseFigGrid
	{
	public:
		MaskFigGrid(int w=100, int h=100);
		~MaskFigGrid();

	private:
		int width,height;
		ZQ::ZQ_DImage<bool> mask;

		bool edit_mode;
		int eraser_size;
		int eraser_select_x;
		int eraser_select_y;
		bool eraser_value;
	protected:
		float color_for_true[3];
		float color_for_false[3];
		float color_for_eraser_true[3];
		float color_for_eraser_false[3];
		float alpha_for_true;
		float alpha_for_false;
		
	public:
		void Draw(float z_shift);
		void OnMouseButton(int button_id, int state, int mx, int my);
		void OnMouseMove(int mx, int my);
		void OnMousePassiveMove(int mx, int my);
		void OnKeyboard(int key_id, int mx, int my){}
		void UpdatePerFrame() {}
		bool Export(const char* export_fold, const char* script_file, int width, int height);
		int GetWidth() const {return mask.width();}
		int GetHeight() const {return mask.height();}
		void GlobalMove(float x, float y){}


	private:
		void _createTweakBar();
		void _addGlobalTweakBar();
		void _addLocalTweakBar();

		static void TW_CALL OnResize(void* clientData);
		void _resize();
		void _loadTexture(bool _forceUpdate);

		void _select_eraser_pos(float mx, float my);
		void _update_mask_with_eraser();
	};

	class RigidFigGrid: public BaseFigGrid
	{
		friend class ZQ_SmokeEditingGUI;
	public:
		RigidFigGrid(float w=200, float h=200);
		~RigidFigGrid();

	private:
		float width, height;
		float rot_angle;
		float scale_x,scale_y;
		float tran_x,tran_y;
		float scale_step;
		bool scaleXenable;
		bool scaleYenable;

	protected:
		float line_width_for_border;
		float line_width_for_boarder_focused;
		float color_for_border[3];
		float color_for_border_focused[3];
		float z_pos_for_line;
		float z_pos_for_tex;
		bool draw_border;

	public:
		float GetXsize() const {return width;}
		float GetYsize() const {return height;}
		float GetRotAngle() const {return rot_angle;}
		float GetTranX() const {return tran_x;}
		float GetTranY() const {return tran_y;}

		void Draw(float z_shift);
		void OnKeyboard(int key_id, int mx, int my);
		void UpdatePerFrame(){}
		void OnMouseButton(int button_id, int state, int mx, int my){}
		void OnMouseMove(int mx, int my){}
		void OnMousePassiveMove(int mx, int my){}
		bool Export(const char* export_fold, const char* script_file, int w, int h);
		int GetWidth()const {return width;}
		int GetHeight()const {return height;}
		void GlobalMove(float x, float y){tran_x += x; tran_y += y;}

	private:
		virtual void _createTweakBar();
		virtual void _addGlobalTweakBar();
		virtual void _addLocalTweakBar();

	};

	class DeformFigGrid: public BaseFigGrid
	{
	public:
		DeformFigGrid();
		~DeformFigGrid();

	protected:
		float z_pos_for_grid;
		float z_pos_for_tex;
		float z_pos_for_fixed_pt;
		float z_pos_for_key_pt;
		float z_pos_for_ctrl_pt;
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
		int width,height;
		float line_weight;
		float angle_weight;
		float distance_weight;
		float distance;
		bool arap;
		bool cur_arap;
		ZQ::ZQ_GridDeformationOptions::NeighborType neighbor_type;
		ZQ::ZQ_DImage<float> coords;
		ZQ::ZQ_DImage<bool> nouseful_flag;
		ZQ::ZQ_DImage<bool> fixed_flag;
		ZQ::ZQ_GridDeformationOptions opt;
		ZQ::ZQ_GridDeformation<float> deform;
		bool edit_mode;
		int total_iteration;
		int cur_iteration;
		int per_frame_iteration;
		bool with_good_init;
		float grid_scale;

		int ctrl_pts_num;
		std::vector<int> ctrl_pts_x;
		std::vector<int> ctrl_pts_y;
		int ctrl_id;

		
		void _reset_deform_iteration();
		static void TW_CALL OnApply(void* clientData);
		void _apply();
		void _rebuild_matrix();

	public:
		void Draw(float z_shift);
		void UpdatePerFrame();
		virtual void OnMouseButton(int button_id, int state, int mx, int my);
		virtual void OnKeyboard(int key_id, int mx, int my);
		virtual void OnMouseMove(int mx, int my){}
		virtual void OnMousePassiveMove(int mx, int my){}
		bool Export(const char* export_fold, const char* script_file, int width, int height);
		int GetWidth() const {return coords.width();}
		int GetHeight() const {return coords.height();}
		void GlobalMove(float x, float y);

	protected:
		void _rotate(float angle);
		void _select_ctrl_pt(float mx, float my);
		void _move_ctrl_pt(float tran_x, float tran_y);
		void _move_all_pts(float x, float y);

	private:
		
		void _draw_grid(float z_shift);
		void _draw_tex(float z_shift);
		void _draw_control_pts(float z_shift);
		virtual void _draw_user_defined(float z_shift){}
		
		

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