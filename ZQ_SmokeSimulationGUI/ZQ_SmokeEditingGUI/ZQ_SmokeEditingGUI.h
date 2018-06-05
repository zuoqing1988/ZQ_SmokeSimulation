#ifndef _ZQ_SMOKE_EDITING_GUI_H_
#define _ZQ_SMOKE_EDITING_GUI_H_

#include <windows.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include "AntTweakBar.h"
#include <vector>
#include "ZQ_SmokeEditingFigGrid.h"
#include "ZQ_SmokeEditingCornersFigGrid.h"
#include "ZQ_SmokeEditingSplitFigGrid.h"
#include "ZQ_SmokeEditingXloopFigGrid.h"
#include "ZQ_SmokeEditingZipperFigGrid.h"
#include "ZQ_SmokeEditingFreeDeformFigGrid.h"

namespace ZQ_SmokeEditing
{
	class ZQ_SmokeEditingGUI
	{
		enum CONST_VAL{
			FIGGRID_TYPE_MAX_NUM = 30,
			FILENAME_MAXLEN = 260,
			MAX_FIGGRID_NUM = 50
		};

		enum FigGridType{
			FIGGRID_RIGID,
			FIGGRID_CORNER4_H,
			FIGGRID_CORNER4_V,
			FIGGRID_CORNER6_H,
			FIGGRID_CORNER6_V,
			FIGGRID_CORNER8_H,
			FIGGRID_CORNER8_V,
			FIGGRID_CORNER10_H,
			FIGGRID_CORNER10_V,
			FIGGRID_CORNER12_H,
			FIGGRID_CORNER12_V,
			FIGGRID_SPLIT_U,
			FIGGRID_SPLIT_D,
			FIGGRID_XLOOP2_T0,
			FIGGRID_XLOOP3_T0,
			FIGGRID_XLOOP4_T0,
			FIGGRID_XLOOP2_T1,
			FIGGRID_XLOOP3_T1,
			FIGGRID_XLOOP4_T1,
			FIGGRID_XLOOP2_T2,
			FIGGRID_XLOOP3_T2,
			FIGGRID_XLOOP4_T2,
			FIGGRID_XLOOP2_T3,
			FIGGRID_XLOOP3_T3,
			FIGGRID_XLOOP4_T3,
			FIGGRID_ZIPPER,
			FIGGRID_FREEDEFORM
		};

	protected:
		ZQ_SmokeEditingGUI();
		~ZQ_SmokeEditingGUI();

	private:
		static ZQ_SmokeEditingGUI* instance;

		static int width;
		static int height;
		static float center[2];
		static float sceneZoom;
		static bool alpha_blend;
		static bool draw_source_group;
		static bool draw_patch_group;

		static TwBar* bar;

		static bool source_or_patch;
		static bool mask_mode;
		static MaskFigGrid* mask_fig;
		static int cur_focus_fig_id_source;
		static int cur_focus_fig_id_patch;
		static std::vector<BaseFigGrid*> source_figGrids;
		static std::vector<BaseFigGrid*> patch_figGrids;
		static FigGridType create_type;
		static TwEnumVal typeEnum[FIGGRID_TYPE_MAX_NUM];

		static char export_fold[FILENAME_MAXLEN];
		static char script_filename[FILENAME_MAXLEN];
		static int frames;
		static int blur_fsize;
		static float blur_sigma;

		static float global_move_step;

		static bool export_source;
		static bool export_patch;

	public:
		static ZQ_SmokeEditingGUI* GetInstance();
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

		static void TW_CALL OnExport(void* clientData);
	};
}

#endif