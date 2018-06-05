#include "ZQ_SmokeEditingGUI.h"
#include <time.h>

using namespace ZQ_SmokeEditing;


ZQ_SmokeEditingGUI* ZQ_SmokeEditingGUI::instance = 0;
int ZQ_SmokeEditingGUI::width = 1024;
int ZQ_SmokeEditingGUI::height = 768;
float ZQ_SmokeEditingGUI::center[2] = {0,0};
float ZQ_SmokeEditingGUI::sceneZoom = 1.0f;
bool ZQ_SmokeEditingGUI::alpha_blend = true;
bool ZQ_SmokeEditingGUI::draw_source_group = true;
bool ZQ_SmokeEditingGUI::draw_patch_group = true;

TwBar* ZQ_SmokeEditingGUI::bar = 0;

bool ZQ_SmokeEditingGUI::source_or_patch = true;
bool ZQ_SmokeEditingGUI::mask_mode = false;
MaskFigGrid* ZQ_SmokeEditingGUI::mask_fig = 0;
int ZQ_SmokeEditingGUI::cur_focus_fig_id_source = -1;
int ZQ_SmokeEditingGUI::cur_focus_fig_id_patch = -1;
std::vector<BaseFigGrid*> ZQ_SmokeEditingGUI::source_figGrids;
std::vector<BaseFigGrid*> ZQ_SmokeEditingGUI::patch_figGrids;
ZQ_SmokeEditingGUI::FigGridType ZQ_SmokeEditingGUI::create_type = ZQ_SmokeEditingGUI::FIGGRID_RIGID;
TwEnumVal ZQ_SmokeEditingGUI::typeEnum[ZQ_SmokeEditingGUI::FIGGRID_TYPE_MAX_NUM] = 
{ 
	{FIGGRID_RIGID, "rigid"},
	{FIGGRID_CORNER4_H,"4 corners H"},
	{FIGGRID_CORNER4_V,"4 corners V"},
	{FIGGRID_CORNER6_H,"6 corners H"},
	{FIGGRID_CORNER6_V,"6 corners V"},
	{FIGGRID_CORNER8_H,"8 corners H"},
	{FIGGRID_CORNER8_V,"8 corners V"},
	{FIGGRID_CORNER10_H,"10 corners H"},
	{FIGGRID_CORNER10_V,"10 corners V"},
	{FIGGRID_CORNER12_H,"12 corners H"},
	{FIGGRID_CORNER12_V,"12 corners V"},
	{FIGGRID_SPLIT_U, "split U"},
	{FIGGRID_SPLIT_D, "split D"},
	{FIGGRID_XLOOP2_T0,"2 xloop T0"},
	{FIGGRID_XLOOP3_T0,"3 xloop T0"},
	{FIGGRID_XLOOP4_T0,"4 xloop T0"},
	{FIGGRID_XLOOP2_T1,"2 xloop T1"},
	{FIGGRID_XLOOP3_T1,"3 xloop T1"},
	{FIGGRID_XLOOP4_T1,"4 xloop T1"},
	{FIGGRID_XLOOP2_T2,"2 xloop T2"},
	{FIGGRID_XLOOP3_T2,"3 xloop T2"},
	{FIGGRID_XLOOP4_T2,"4 xloop T2"},
	{FIGGRID_XLOOP2_T3,"2 xloop T3"},
	{FIGGRID_XLOOP3_T3,"3 xloop T3"},
	{FIGGRID_XLOOP4_T3,"4 xloop T3"},
	{FIGGRID_ZIPPER,"zipper"},
	{FIGGRID_FREEDEFORM,"free deform"}
};

char ZQ_SmokeEditingGUI::export_fold[FILENAME_MAXLEN] = "export";
char ZQ_SmokeEditingGUI::script_filename[FILENAME_MAXLEN] = "script.bat";
int ZQ_SmokeEditingGUI::frames = 200;
int ZQ_SmokeEditingGUI::blur_fsize = 10;
float ZQ_SmokeEditingGUI::blur_sigma = 20;
bool ZQ_SmokeEditingGUI::export_source = true;
bool ZQ_SmokeEditingGUI::export_patch = true;

float ZQ_SmokeEditingGUI::global_move_step = 10;

ZQ_SmokeEditingGUI::ZQ_SmokeEditingGUI()
{
}

ZQ_SmokeEditingGUI::~ZQ_SmokeEditingGUI()
{
	for(int i = 0;i < source_figGrids.size();i++)
		delete source_figGrids[i];
	for(int i = 0;i < patch_figGrids.size();i++)
		delete patch_figGrids[i];
	delete mask_fig;
}

ZQ_SmokeEditingGUI* ZQ_SmokeEditingGUI::GetInstance()
{
	if(instance == 0)
		instance = new ZQ_SmokeEditingGUI();
	return instance;
}


bool ZQ_SmokeEditingGUI::Init(int* argc, char** argv)
{
	glutInit(argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("Smoke Editing GUI");
	glutCreateMenu(NULL);
	glewInit();

	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	atexit(Terminate);

	TwInit(TW_OPENGL, NULL);

	glutMouseFunc((GLUTmousebuttonfun)OnMouseButton);
	glutMotionFunc((GLUTmousemotionfun)OnMouseMove);
	glutPassiveMotionFunc((GLUTmousemotionfun)OnMousePassiveMove);
	glutKeyboardFunc((GLUTkeyboardfun)OnKeyBoard);
	glutSpecialFunc((GLUTspecialfun)OnSpecial);


	int opened = 0;
	TwGLUTModifiersFunc(glutGetModifiers);

	TwDefine("GLOBAL fontsize=3 fontresizable=false");
	bar = TwNewBar("Menu");
	TwDefine("Menu refresh=0.1 movable=false resizable=false");
	TwDefine("Menu color='0 0 0' alpha='100'");

	TwAddSeparator(bar,"Display Sep","label='Display'");
	TwAddVarRW(bar,"Center X",TW_TYPE_FLOAT,&center[0],"visible=true");
	TwAddVarRW(bar,"Center Y",TW_TYPE_FLOAT,&center[1],"visible=true");
	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &sceneZoom, "visible=true min=0.01 max=100");
	TwAddVarRW(bar,"alpha_blend",TW_TYPE_BOOL8,&alpha_blend,"label='alpha blend'");
	TwAddVarRW(bar,"draw_source_group",TW_TYPE_BOOL8,&draw_source_group,"label='draw source'");
	TwAddVarRW(bar,"draw_source_patch",TW_TYPE_BOOL8,&draw_patch_group,"label='draw patch'");

	TwAddSeparator(bar,"Create Sep","label='Create'");
	TwAddVarRW(bar,"source_or_patch",TW_TYPE_BOOL8,&source_or_patch,"label='mode' readonly=true");
	TwAddButton(bar,"ChangeMode",(TwButtonCallback)OnChangeMode,0,"label='Change Mode'");
	TwType FigGridTWType = TwDefineEnum("FigGridTypeEnum", typeEnum, 27);
	TwAddVarRW(bar, "FigGridType", FigGridTWType, &create_type, "label='type'");
	TwAddButton(bar,"CreateButton",(TwButtonCallback)OnCreate,0,"label='Create'");
	TwAddButton(bar,"DeleteButton",(TwButtonCallback)OnDelete,0,"label='Delete'");

	TwAddSeparator(bar,"Export Sep","label='Export'");
	TwAddVarRW(bar,"blur_fsize",TW_TYPE_INT32,&blur_fsize,"label='fsize'");
	TwAddVarRW(bar,"blur_sigma",TW_TYPE_FLOAT,&blur_sigma,"label='sigma'");
	TwAddVarRW(bar,"frames",TW_TYPE_INT32,&frames,"label='frames'");
	TwAddVarRW(bar,"export_fold",TW_TYPE_CSSTRING(sizeof(export_fold)),&export_fold,"label='fold'");
	TwAddVarRW(bar,"export_source",TW_TYPE_BOOL8,&export_source,"label='export source'");
	TwAddVarRW(bar,"export_patch",TW_TYPE_BOOL8,&export_patch,"label='export patch'");
	TwAddButton(bar,"Export",(TwButtonCallback)OnExport,0,"label='Export'");
	

	mask_fig = new MaskFigGrid(400,400);
	return true;
}

void ZQ_SmokeEditingGUI::MainLoop()
{
	glutMainLoop();
}

void ZQ_SmokeEditingGUI::Terminate()
{
	TwTerminate();
}

void ZQ_SmokeEditingGUI::Display()
{
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
	glPushMatrix(); 

	glScalef(sceneZoom, sceneZoom, sceneZoom);

	glTranslatef(-center[0],-center[1],0);

	if(alpha_blend)
	{
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	}

	clock_t t1 = clock();
	Draw();
	clock_t t2 = clock();

	if(alpha_blend)
	{
		glDisable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);
	}

	glPopMatrix();

	//// Draw tweak bars
	TwDraw();

	// Present frame buffer
	glutSwapBuffers();

	// Recall Display at next frame
	glutPostRedisplay();

	clock_t t3 = clock();
	UpdataPerFrame();
	clock_t t4 = clock();

	//printf("draw:%8.3f update:%8.3f\n",1.0*(t2-t1),1.0*(t4-t3));
	Sleep(30);
}


void ZQ_SmokeEditingGUI::Reshape(int w, int h)
{
	width = w;
	height = h;
	// Set OpenGL viewport and camera
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width/2.0,width/2.0,-height/2.0,height/2.0,-1000,1000);

	// Send the new window size to AntTweakBar
	TwWindowSize(width, height);
	_resize_bar(width,height);
}


void ZQ_SmokeEditingGUI::OnMouseMove(int mx, int my)
{
	float real_mx = (mx - width*0.5)/sceneZoom + center[0];
	float real_my = ((height-1-my) - height*0.5)/sceneZoom + center[1];
	if(!TwEventMouseMotionGLUT(mx,my))
	{
		if(mask_mode)
		{
			mask_fig->OnMouseMove(real_mx,real_my);
		}

		if(source_or_patch)
		{
			if(cur_focus_fig_id_source >=0)
			{
				source_figGrids[cur_focus_fig_id_source]->OnMouseMove(real_mx,real_my);
			}
		}
		else
		{
			if(cur_focus_fig_id_patch >= 0)
			{
				patch_figGrids[cur_focus_fig_id_patch]->OnMouseMove(real_mx,real_my);
			}
		}
	}
}
void ZQ_SmokeEditingGUI::OnMousePassiveMove(int mx, int my)
{
	float real_mx = (mx - width*0.5)/sceneZoom + center[0];
	float real_my = ((height-1-my) - height*0.5)/sceneZoom + center[1];

	if(!TwEventMouseMotionGLUT(mx,my))
	{
		if(mask_mode)
		{
			mask_fig->OnMousePassiveMove(real_mx,real_my);
		}
		if(source_or_patch)
		{
			if(cur_focus_fig_id_source >= 0)
			{
				source_figGrids[cur_focus_fig_id_source]->OnMousePassiveMove(real_mx,real_my);
			}
		}
		else
		{
			if(cur_focus_fig_id_patch >= 0)
			{
				patch_figGrids[cur_focus_fig_id_patch]->OnMousePassiveMove(real_mx,real_my);
			}
		}
		
	}
}

void ZQ_SmokeEditingGUI::OnMouseButton(int button, int state, int mx, int my)
{
	float real_mx = (mx - width*0.5)/sceneZoom + center[0];
	float real_my = ((height-1-my) - height*0.5)/sceneZoom + center[1];

	if(!TwEventMouseButtonGLUT(button,state,mx,my))
	{
		if(mask_mode)
		{
			mask_fig->OnMouseButton(button,state,real_mx,real_my);
		}
		if(source_or_patch)
		{
			if(cur_focus_fig_id_source >= 0)
			{
				source_figGrids[cur_focus_fig_id_source]->OnMouseButton(button,state,real_mx,real_my);
			}
		}
		else
		{
			if(cur_focus_fig_id_patch >= 0)
			{
				patch_figGrids[cur_focus_fig_id_patch]->OnMouseButton(button,state,real_mx,real_my);
			}
		}
	}
}

void ZQ_SmokeEditingGUI::OnKeyBoard(int key, int mx, int my)
{
	if(!TwEventKeyboardGLUT(key,mx,my))
	{	
		switch(key)
		{
		case '+':case '=':
			{
				OnNext(0);
			}
			break;
		case '-':case '_':
			{
				OnPrevious(0);
			}
			break;
		case 'm':case 'M':
			{
				mask_mode = !mask_mode;
				mask_fig->SetFocused(mask_mode);
			}
			break;
		case 'i':case 'I':
			{
				OnKeyUp();
			}
			break;
		case 'k':case 'K':
			{
				OnKeyDown();
			}
			break;
		case 'j':case 'J':
			{
				OnKeyLeft();
			}
			break;
		case 'l':case 'L':
			{
				OnKeyRight();
			}
			break;
		default:
			{
				float real_mx = (mx - width*0.5)/sceneZoom + center[0];
				float real_my = ((height-1-my) - height*0.5)/sceneZoom + center[1];
				if(mask_mode)
					mask_fig->OnKeyboard(key,real_mx,real_my);
				if(source_or_patch)
				{
					if(cur_focus_fig_id_source >= 0)
					{
						source_figGrids[cur_focus_fig_id_source]->OnKeyboard(key,real_mx,real_my);
					}
				}
				else
				{
					if(cur_focus_fig_id_patch >= 0)
					{
						patch_figGrids[cur_focus_fig_id_patch]->OnKeyboard(key,real_mx,real_my);
					}
				}	
			}
		}
	}

	if(key == TW_KEY_ESCAPE) 
	{
		exit(0);
	}
}

void ZQ_SmokeEditingGUI::OnSpecial(int key, int mx, int my)
{
	if(!TwEventSpecialGLUT(key,mx,my))
	{
	}
}

void ZQ_SmokeEditingGUI::OnKeyUp()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(0,global_move_step);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(0,global_move_step);
}

void ZQ_SmokeEditingGUI::OnKeyDown()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(0,-global_move_step);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(0,-global_move_step);
}

void ZQ_SmokeEditingGUI::OnKeyLeft()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(-global_move_step,0);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(-global_move_step,0);
}

void ZQ_SmokeEditingGUI::OnKeyRight()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(global_move_step,0);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(global_move_step,0);
}

void ZQ_SmokeEditingGUI::OnCreate(void* clientData)
{
	BaseFigGrid* fig = 0;
	switch(create_type)
	{
	case FIGGRID_RIGID:
		fig = new RigidFigGrid();
		break;
	case FIGGRID_CORNER4_H:
		fig = new CornersDeformFigGrid<2,true>();
		break;
	case FIGGRID_CORNER4_V:
		fig = new CornersDeformFigGrid<2,false>();
		break;
	case FIGGRID_CORNER6_H:
		fig = new CornersDeformFigGrid<3,true>();
		break;
	case FIGGRID_CORNER6_V:
		fig = new CornersDeformFigGrid<3,false>();
		break;
	case FIGGRID_CORNER8_H:
		fig = new CornersDeformFigGrid<4,true>();
		break;
	case FIGGRID_CORNER8_V:
		fig = new CornersDeformFigGrid<4,false>();
		break;
	case FIGGRID_CORNER10_H:
		fig = new CornersDeformFigGrid<5,true>();
		break;
	case FIGGRID_CORNER10_V:
		fig = new CornersDeformFigGrid<5,false>();
		break;
	case FIGGRID_CORNER12_H:
		fig = new CornersDeformFigGrid<6,true>();
		break;
	case FIGGRID_CORNER12_V:
		fig = new CornersDeformFigGrid<6,false>();
		break;
	case FIGGRID_SPLIT_U:
		fig = new SplitDeformGrid(true);
		break;
	case FIGGRID_SPLIT_D:
		fig = new SplitDeformGrid(false);
		break;
	case FIGGRID_XLOOP2_T0:
		fig = new XloopDeformGrid<2,0>();
		break;
	case FIGGRID_XLOOP3_T0:
		fig = new XloopDeformGrid<3,0>();
		break;
	case FIGGRID_XLOOP4_T0:
		fig = new XloopDeformGrid<4,0>();
		break;
	case FIGGRID_XLOOP2_T1:
		fig = new XloopDeformGrid<2,1>();
		break;
	case FIGGRID_XLOOP3_T1:
		fig = new XloopDeformGrid<3,1>();
		break;
	case FIGGRID_XLOOP4_T1:
		fig = new XloopDeformGrid<4,1>();
		break;
	case FIGGRID_XLOOP2_T2:
		fig = new XloopDeformGrid<2,2>();
		break;
	case FIGGRID_XLOOP3_T2:
		fig = new XloopDeformGrid<3,2>();
		break;
	case FIGGRID_XLOOP4_T2:
		fig = new XloopDeformGrid<4,2>();
		break;
	case FIGGRID_XLOOP2_T3:
		fig = new XloopDeformGrid<2,3>();
		break;
	case FIGGRID_XLOOP3_T3:
		fig = new XloopDeformGrid<3,3>();
		break;
	case FIGGRID_XLOOP4_T3:
		fig = new XloopDeformGrid<4,3>();
		break;
	case FIGGRID_ZIPPER:
		fig = new ZipperDeformGrid();
		break;
	case FIGGRID_FREEDEFORM:
		fig = new FreeDeformFigGird();
		break;
	}

	if(fig != 0)
	{
		if(source_or_patch)
		{
			int cur_fig_num = source_figGrids.size();
			for(int i = 0;i < cur_fig_num;i++)
				source_figGrids[i]->SetFocused(false);
			fig->SetFocused(true);
			source_figGrids.push_back(fig);
			cur_focus_fig_id_source = cur_fig_num;
		}
		else
		{
			int cur_fig_num = patch_figGrids.size();
			for(int i = 0;i < cur_fig_num;i++)
				patch_figGrids[i]->SetFocused(false);
			fig->SetFocused(true);
			patch_figGrids.push_back(fig);
			cur_focus_fig_id_patch = cur_fig_num;
		}
	}
}

void ZQ_SmokeEditingGUI::OnDelete(void* clientData)
{
	if(source_or_patch)
	{
		if(cur_focus_fig_id_source >= 0)
		{
			BaseFigGrid* g = source_figGrids[cur_focus_fig_id_source];
			source_figGrids.erase(source_figGrids.begin()+cur_focus_fig_id_source);
			delete g;
			if(cur_focus_fig_id_source > 0)
				cur_focus_fig_id_source--;
			else
			{
				if(source_figGrids.size() > 0)
				{
					cur_focus_fig_id_source = 0;
				}
				else
				{
					cur_focus_fig_id_source = -1;
				}
			}
			if(cur_focus_fig_id_source >= 0)
				source_figGrids[cur_focus_fig_id_source]->SetFocused(true);
		}
	}
	else
	{
		if(cur_focus_fig_id_patch >= 0)
		{
			BaseFigGrid* g = patch_figGrids[cur_focus_fig_id_patch];
			patch_figGrids.erase(patch_figGrids.begin()+cur_focus_fig_id_patch);
			delete g;
			if(cur_focus_fig_id_patch > 0)
				cur_focus_fig_id_patch--;
			else
			{
				if(patch_figGrids.size() > 0)
				{
					cur_focus_fig_id_patch = 0;
				}
				else
				{
					cur_focus_fig_id_patch = -1;
				}
			}
			if(cur_focus_fig_id_patch >= 0)
				patch_figGrids[cur_focus_fig_id_patch]->SetFocused(true);
		}
	}
	
}


void ZQ_SmokeEditingGUI::OnNext(void*clientData)
{
	if(source_or_patch)
	{
		int num = source_figGrids.size();
		if(num > 1)
		{
			source_figGrids[cur_focus_fig_id_source]->SetFocused(false);
			cur_focus_fig_id_source = (cur_focus_fig_id_source+1)%num;
			source_figGrids[cur_focus_fig_id_source]->SetFocused(true);
		}
	}
	else
	{
		int num = patch_figGrids.size();
		if(num > 1)
		{
			patch_figGrids[cur_focus_fig_id_patch]->SetFocused(false);
			cur_focus_fig_id_patch = (cur_focus_fig_id_patch+1)%num;
			patch_figGrids[cur_focus_fig_id_patch]->SetFocused(true);
		}
	}
	
}

void ZQ_SmokeEditingGUI::OnPrevious(void* clientData)
{
	if(source_or_patch)
	{
		int num = source_figGrids.size();
		if(num > 1)
		{
			source_figGrids[cur_focus_fig_id_source]->SetFocused(false);
			cur_focus_fig_id_source = (cur_focus_fig_id_source-1+num)%num;
			source_figGrids[cur_focus_fig_id_source]->SetFocused(true);
		}
	}
	else
	{
		int num = patch_figGrids.size();
		if(num > 1)
		{
			patch_figGrids[cur_focus_fig_id_patch]->SetFocused(false);
			cur_focus_fig_id_patch = (cur_focus_fig_id_patch-1+num)%num;
			patch_figGrids[cur_focus_fig_id_patch]->SetFocused(true);
		}
	}
}

void ZQ_SmokeEditingGUI::OnChangeMode(void* clientData)
{
	if(source_or_patch)
	{
		if(cur_focus_fig_id_source >= 0)
			source_figGrids[cur_focus_fig_id_source]->SetFocused(false);
	}
	else
	{
		if(cur_focus_fig_id_patch >= 0)
			patch_figGrids[cur_focus_fig_id_patch]->SetFocused(false);
	}

	source_or_patch = !source_or_patch;


	if(source_or_patch)
	{
		if(cur_focus_fig_id_source >= 0)
			source_figGrids[cur_focus_fig_id_source]->SetFocused(true);
	}
	else
	{
		if(cur_focus_fig_id_patch >= 0)
			patch_figGrids[cur_focus_fig_id_patch]->SetFocused(true);
	}
}

void ZQ_SmokeEditingGUI::Draw()
{
	if(draw_source_group)
	{	for(int i = 0;i < source_figGrids.size();i++)
			source_figGrids[i]->Draw(-(float)i);
	}
	if(draw_patch_group)
	{
		for(int i = 0;i < patch_figGrids.size();i++)
			patch_figGrids[i]->Draw(-(float)(source_figGrids.size())-(float)i);
	}
	if(mask_mode)
		mask_fig->Draw(0);
}

void ZQ_SmokeEditingGUI::UpdataPerFrame()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->UpdatePerFrame();
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->UpdatePerFrame();
	if(mask_mode)
		mask_fig->UpdatePerFrame();
}

void ZQ_SmokeEditingGUI::_resize_bar(int w,int h)
{
	char buf[500];
	int bar_width = 200;
	sprintf(buf,"Menu size=\'%d %d\' position=\'%d %d\'",bar_width,500,w-bar_width,0);
	TwDefine(buf);
}

void ZQ_SmokeEditingGUI::OnExport(void* clientData)
{
	char cmdbuf[2000];
	std::string cmdstr;

	char* source_fold = "source";
	char* patch_fold = "patch";
	char* source_patch_fold = "source_patch";
	char* mask_filename = "mask.png";
	char* high_of_source_fold = "source_H";
	char* high_of_patch_fold = "patch_H";
	char* low_of_source_patch_fold = "lowpart";
	char* poisson_edit_high_fold = "highpart";
	char* final_result_fold = "result";

	bool export_all = export_source && export_patch;

	sprintf(cmdbuf,"mkdir %s",export_fold);
	system(cmdbuf);

	if(export_source)
	{
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,source_fold);
		system(cmdbuf);
	}

	if(export_patch)
	{
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,patch_fold);
		system(cmdbuf);
	}
	
	if(export_all)
	{
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,source_patch_fold);
		system(cmdbuf);
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,high_of_source_fold);
		system(cmdbuf);
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,high_of_patch_fold);
		system(cmdbuf);
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,low_of_source_patch_fold);
		system(cmdbuf);
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,poisson_edit_high_fold);
		system(cmdbuf);
		sprintf(cmdbuf,"mkdir %s\\%s",export_fold,final_result_fold);
		system(cmdbuf);
	}
	
	bool wrong_flag = false;

	if(export_source)
	{
		const char* source_bat_filename = "source.bat";
		printf("stage: export source\n");
		if(source_figGrids.size() == 0)
		{
			printf("source has 0 image\n");
			wrong_flag = true;
		}
		else if(source_figGrids.size() == 1)
		{
			std::string cur_fold,cur_script;
			sprintf(cmdbuf,"%s\\%s",export_fold,source_fold);
			cur_fold.append(cmdbuf);
			sprintf(cmdbuf,"mkdir %s",cur_fold.c_str());
			system(cmdbuf);
			cur_script.append(source_bat_filename);
			sprintf(cmdbuf,"type nul>%s",cur_script.c_str());
			system(cmdbuf);
			if(!source_figGrids[0]->Export(cur_fold.c_str(),cur_script.c_str(),mask_fig->GetWidth(),mask_fig->GetHeight()))
			{
				wrong_flag = true;
				printf("failed to export %s \n","source");
			}
		}
		else
		{
			std::string cur_fold,cur_script;
			cur_script.append(source_bat_filename);
			sprintf(cmdbuf,"type nul>%s",cur_script.c_str());
			system(cmdbuf);
			for(int i = 0;i < source_figGrids.size();i++)
			{
				cur_fold.clear();
				sprintf(cmdbuf,"%s\\%s\\%d",export_fold,source_fold,i);
				cur_fold.append(cmdbuf);
				sprintf(cmdbuf,"mkdir %s",cur_fold.c_str());
				system(cmdbuf);
				if(!source_figGrids[i]->Export(cur_fold.c_str(),cur_script.c_str(),mask_fig->GetWidth(),mask_fig->GetHeight()))
				{
					wrong_flag = true;
					printf("failed to export %s \n","source");
				}
			}
		}

		
		//directly merge all sources
		if(wrong_flag)
		{
			printf("not continued\n");
			return;
		}
		FILE* out1 = fopen(source_bat_filename,"a");
		if(source_figGrids.size() > 1)
		{
			for(int fr = 0;fr < frames;fr++)
			{
				cmdstr.clear();
				cmdstr.append("ZQ_MergeImage MethodType MERGE_DIRECTLY OutputFile ");
				sprintf(cmdbuf,"%s\\%s\\%d.png ",export_fold,source_fold,fr);
				cmdstr.append(cmdbuf);
				for(int i = 0;i < source_figGrids.size();i++)
				{
					int base_id = source_figGrids[i]->GetBaseId();
					int total_frame = source_figGrids[i]->GetTotalFrame();
					if(total_frame <= 0)
					{
						wrong_flag = true;
						printf("error: source [%d] total_frame = %d\n",i,total_frame);
					}
					else
					{
						sprintf(cmdbuf,"DirectSource %s\\%s\\%d\\%d.png ",export_fold,source_fold,i,(fr+base_id)%total_frame);
						cmdstr.append(cmdbuf);
					}
				}
				fprintf(out1,"%s\n",cmdstr.c_str());
			}
		}

		fclose(out1);
	}
	

	if(export_patch)
	{
		const char* patch_bat_filename = "patch.bat";
		printf("stage: export patch\n");
		if(patch_figGrids.size() == 0)
		{
			printf("patch has 0 image\n");
			wrong_flag = true;
		}
		else if(patch_figGrids.size() == 1)
		{
			std::string cur_fold,cur_script;
			sprintf(cmdbuf,"%s\\%s",export_fold,patch_fold);
			cur_fold.append(cmdbuf);
			sprintf(cmdbuf,"mkdir %s",cur_fold.c_str());
			system(cmdbuf);
			cur_script.append(patch_bat_filename);
			sprintf(cmdbuf,"type nul>%s",cur_script.c_str());
			system(cmdbuf);
			if(!patch_figGrids[0]->Export(cur_fold.c_str(),cur_script.c_str(),mask_fig->GetWidth(),mask_fig->GetHeight()))
			{
				wrong_flag = true;
				printf("failed to export %s \n","patch");
			}
		}
		else
		{
			std::string cur_fold,cur_script;
			cur_script.append(patch_bat_filename);
			sprintf(cmdbuf,"type nul>%s",cur_script.c_str());
			system(cmdbuf);
			for(int i = 0;i < patch_figGrids.size();i++)
			{
				cur_fold.clear();
				sprintf(cmdbuf,"%s\\%s\\%d",export_fold,patch_fold,i);
				cur_fold.append(cmdbuf);
				sprintf(cmdbuf,"mkdir %s",cur_fold.c_str());
				system(cmdbuf);
				if(!patch_figGrids[i]->Export(cur_fold.c_str(),cur_script.c_str(),mask_fig->GetWidth(),mask_fig->GetHeight()))
				{
					wrong_flag = true;
					printf("failed to export %s \n","patch");
				}
			}
		}


		//directly merge all patches
		if(wrong_flag)
		{
			printf("not continued\n");
			return;
		}
		FILE* out2 = fopen(patch_bat_filename,"a");
		if(patch_figGrids.size() > 1)
		{
			for(int fr = 0;fr < frames;fr++)
			{
				cmdstr.clear();
				cmdstr.append("ZQ_MergeImage MethodType MERGE_DIRECTLY OutputFile ");
				sprintf(cmdbuf,"%s\\%s\\%d.png ",export_fold,patch_fold,fr);
				cmdstr.append(cmdbuf);
				for(int i = 0;i < patch_figGrids.size();i++)
				{
					int base_id = patch_figGrids[i]->GetBaseId();
					int total_frame = patch_figGrids[i]->GetTotalFrame();
					if(total_frame <= 0)
					{
						wrong_flag = true;
						printf("error: patch [%d] total_frame = %d\n",i,total_frame);
					}
					else
					{
						sprintf(cmdbuf,"DirectSource %s\\%s\\%d\\%d.png ",export_fold,patch_fold,i,(fr+base_id)%total_frame);
						cmdstr.append(cmdbuf);
					}
				}
				fprintf(out2,"%s\n",cmdstr.c_str());
			}
		}

		fclose(out2);
	}

	if(export_all)
	{
		//stage3
		printf("stage: get the high part\n");
		FILE* out3 = fopen("highpart.bat","w");
		for(int fr = 0;fr < frames;fr++)
		{
			cmdstr.clear();
			sprintf(cmdbuf,"ZQ_MergeImage MethodType IMAGE_BLUR SourceFile %s\\%s\\%d.png fsize %d sigma %.1f ", export_fold,source_fold,fr,blur_fsize,blur_sigma);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"HighPartFile %s\\%s\\%d.png ", export_fold,high_of_source_fold,fr);
			cmdstr.append(cmdbuf);
			fprintf(out3,"%s\n",cmdstr.c_str());
		}
		
		for(int fr = 0;fr < frames;fr++)
		{
			cmdstr.clear();
			sprintf(cmdbuf,"ZQ_MergeImage MethodType IMAGE_BLUR SourceFile %s\\%s\\%d.png fsize %d sigma %.1f ", export_fold,patch_fold,fr,blur_fsize,blur_sigma);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"HighPartFile %s\\%s\\%d.png ", export_fold,high_of_patch_fold,fr);
			cmdstr.append(cmdbuf);

			fprintf(out3,"%s\n",cmdstr.c_str());
		}

		
		// Poisson editing for high part
		cmdstr.clear();
		sprintf(cmdbuf,"ZQ_PoissonEditing3DCuda %d %s\\%s %s %s\\%s %s %s\\%s %s\\%s",
			frames,
			export_fold,high_of_source_fold,"png",
			export_fold,high_of_patch_fold,"png",
			export_fold,mask_filename,
			export_fold,poisson_edit_high_fold);
		cmdstr.append(cmdbuf);
		fprintf(out3,"%s\n",cmdstr.c_str());
		fclose(out3);

		//stage4
		printf("stage: get the low part\n");
		if(!mask_fig->Export(export_fold,"mask.png",0,0))
		{
			printf("failed to export mask\n");
			wrong_flag = true;
		}
		if(wrong_flag)
		{
			printf("not continued\n");
			return;
		}

		FILE* out4 = fopen("lowpart.bat","w");
		for(int fr = 0;fr < frames;fr++)
		{
			cmdstr.clear();
			sprintf(cmdbuf,"ZQ_MergeImage MethodType MERGE_SOURCE_PATCH OutputFile %s\\%s\\%d.png ",export_fold,source_patch_fold,fr);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"SourceFile %s\\%s\\%d.png ", export_fold,source_fold,fr);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"PatchFile %s\\%s\\%d.png ", export_fold,patch_fold,fr);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"MaskFile %s\\%s ",export_fold,mask_filename);
			cmdstr.append(cmdbuf);

			fprintf(out4,"%s\n",cmdstr.c_str());
		}
		
		// get the low part of source_and_patch
		for(int fr = 0;fr < frames;fr++)
		{
			cmdstr.clear();
			sprintf(cmdbuf,"ZQ_MergeImage MethodType IMAGE_BLUR SourceFile %s\\%s\\%d.png fsize %d sigma %.1f ", export_fold,source_patch_fold,fr,blur_fsize,blur_sigma);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"OutputFile %s\\%s\\%d.png ", export_fold,low_of_source_patch_fold,fr);
			cmdstr.append(cmdbuf);
			fprintf(out4,"%s\n",cmdstr.c_str());
		}
		fclose(out4);



		//stage 5
		printf("stage: merge low and high part\n");

		FILE* out5 = fopen("final.bat","w");
		for(int fr = 0;fr < frames;fr++)
		{
			cmdstr.clear();
			sprintf(cmdbuf,"ZQ_MergeImage MethodType MERGE_LOW_HIGH OutputFile %s\\%s\\%d.png ",export_fold,final_result_fold,fr);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"SourceFile %s\\%s\\%d.png ",export_fold,low_of_source_patch_fold,fr);
			cmdstr.append(cmdbuf);
			sprintf(cmdbuf,"HighPartFile %s\\%s\\%d.png ",export_fold,poisson_edit_high_fold,fr);
			cmdstr.append(cmdbuf);
			fprintf(out5,"%s\n",cmdstr.c_str());
		}
		fclose(out5);
	}
}