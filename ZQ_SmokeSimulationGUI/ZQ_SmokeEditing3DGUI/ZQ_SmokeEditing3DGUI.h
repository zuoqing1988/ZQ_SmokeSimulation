#ifndef _ZQ_SMOKE_EDITING_3D_GUI_H_
#define _ZQ_SMOKE_EDITING_3D_GUI_H_
#pragma once

#include <windows.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include "AntTweakBar.h"
#include <vector>
#include "ZQ_SmokeEditingFigGrid3D.h"
#include "ZQ_SmokeEditingCornersFigGrid3D.h"
#include "ZQ_SmokeEditingXloopFigGrid3D.h"

namespace ZQ_SmokeEditing3D
{
	class ZQ_SmokeEditing3DGUI
	{
		enum CONST_VAL{
			FIGGRID_TYPE_MAX_NUM = 30,
			FILENAME_MAXLEN = 260,
			MAX_FIGGRID_NUM = 50
		};

		enum FigGridType{		
			FIGGRID_CORNER8_T0_H,
			FIGGRID_CORNER8_T0_V,
			FIGGRID_CORNER8_T1_H,
			FIGGRID_CORNER8_T1_V,
			FIGGRID_CORNER12_T0_H,
			FIGGRID_CORNER12_T0_V,
			FIGGRID_CORNER12_T1_H,
			FIGGRID_CORNER12_T1_V,
			FIGGRID_CORNER16_T0_H,
			FIGGRID_CORNER16_T0_V,
			FIGGRID_CORNER16_T1_H,
			FIGGRID_CORNER16_T1_V,
			FIGGRID_XLOOP3_T0,
			FIGGRID_XLOOP3_T1,
			FIGGRID_XLOOP4_T0,
			FIGGRID_XLOOP4_T1
		};

	protected:
		ZQ_SmokeEditing3DGUI();
		~ZQ_SmokeEditing3DGUI();

	private:
		static ZQ_SmokeEditing3DGUI* instance;

		static int width;
		static int height;
		static float fovy;
		static float eyepos[3];
		static float sceneZoom;
		static float sceneRotation[4];
		static float sceneRotationStart[4];
		static bool alpha_blend;
		static bool draw_source_group;
		static bool draw_patch_group;

		static TwBar* bar;

		static bool source_or_patch;
		static int cur_focus_fig_id_source;
		static int cur_focus_fig_id_patch;
		static std::vector<BaseFigGrid3D*> source_figGrids;
		static std::vector<BaseFigGrid3D*> patch_figGrids;
		static FigGridType create_type;
		static TwEnumVal typeEnum[FIGGRID_TYPE_MAX_NUM];

		static char export_fold[FILENAME_MAXLEN];
		static char script_filename[FILENAME_MAXLEN];
		static float global_move_step;
	public:
		static ZQ_SmokeEditing3DGUI* GetInstance();
		bool Init(int* argc, char** argv);
		void MainLoop();

	private:
		static void Terminate(void);
		static void Display(void);
		static void Reshape(int width, int height);
		static void  OnMouseMove(int mx, int my);
		static void  OnMousePassiveMove(int mx, int my);
		static void  OnMouseButton(int button, int state, int mx, int my);
		static void  OnKeyBoard(int key, int mx, int my);
		static void  OnSpecial(int key, int mx, int my);

		static void TW_CALL OnSetOriginalView(void* clientData);

		static void OnKeyUp();
		static void OnKeyDown();
		static void OnKeyLeft();
		static void OnKeyRight();

		static void TW_CALL OnCreate(void* clientData);
		static void TW_CALL OnDelete(void* clientData);
		static void TW_CALL OnPrevious(void* clientData);
		static void TW_CALL OnNext(void*clientData);

		static void TW_CALL OnChangeMode(void* clientData);

		static void Draw();
		static void UpdataPerFrame();

		static void _resize_bar(int w,int h);
		static void _mouse_pos_to_ray(int mx, int my, float ray_ori[3], float ray_dir[3]);

		static void TW_CALL OnExport(void* clientData);

		static void SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat);
		static void ConvertQuaternionToMatrix(const float *quat, float *mat);
		static void MultiplyQuaternions(const float *q1, const float *q2, float *qout);
	};
}

#endif