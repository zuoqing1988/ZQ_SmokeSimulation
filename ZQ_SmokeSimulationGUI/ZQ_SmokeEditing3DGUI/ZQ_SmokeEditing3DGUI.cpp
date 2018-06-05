#include "ZQ_SmokeEditing3DGUI.h"
#include "ZQ_MathBase.h"
#include <time.h>

using namespace ZQ_SmokeEditing3D;
using namespace ZQ;


ZQ_SmokeEditing3DGUI* ZQ_SmokeEditing3DGUI::instance = 0;
int ZQ_SmokeEditing3DGUI::width = 1024;
int ZQ_SmokeEditing3DGUI::height = 768;
float ZQ_SmokeEditing3DGUI::fovy = 30;
float ZQ_SmokeEditing3DGUI::eyepos[3] = {0,0,200};
float ZQ_SmokeEditing3DGUI::sceneZoom = 1.0f;
float ZQ_SmokeEditing3DGUI::sceneRotation[4] = {0,0,0,1};//{-0.707106781186548,0,0,0.707106781186548};
float ZQ_SmokeEditing3DGUI::sceneRotationStart[4] = {0,0,0,1};//{-0.707106781186548,0,0,0.707106781186548};
bool ZQ_SmokeEditing3DGUI::alpha_blend = true;
bool ZQ_SmokeEditing3DGUI::draw_source_group = true;
bool ZQ_SmokeEditing3DGUI::draw_patch_group = true;

TwBar* ZQ_SmokeEditing3DGUI::bar = 0;
bool ZQ_SmokeEditing3DGUI::source_or_patch = true;
int ZQ_SmokeEditing3DGUI::cur_focus_fig_id_source = -1;
int ZQ_SmokeEditing3DGUI::cur_focus_fig_id_patch = -1;
std::vector<BaseFigGrid3D*> ZQ_SmokeEditing3DGUI::source_figGrids;
std::vector<BaseFigGrid3D*> ZQ_SmokeEditing3DGUI::patch_figGrids;
ZQ_SmokeEditing3DGUI::FigGridType ZQ_SmokeEditing3DGUI::create_type = ZQ_SmokeEditing3DGUI::FIGGRID_CORNER8_T0_H;
TwEnumVal ZQ_SmokeEditing3DGUI::typeEnum[ZQ_SmokeEditing3DGUI::FIGGRID_TYPE_MAX_NUM] = 
{ 
	{FIGGRID_CORNER8_T0_H,"C8 T0 H"},
	{FIGGRID_CORNER8_T0_V,"C8 T0 V"},
	{FIGGRID_CORNER8_T1_H,"C8 T1 H"},
	{FIGGRID_CORNER8_T1_V,"C8 T1 V"},
	{FIGGRID_CORNER12_T0_H,"C12 T0 H"},
	{FIGGRID_CORNER12_T0_V,"C12 T0 V"},
	{FIGGRID_CORNER12_T1_H,"C12 T1 H"},
	{FIGGRID_CORNER12_T1_V,"C12 T1 V"},
	{FIGGRID_CORNER16_T0_H,"C16 T0 H"},
	{FIGGRID_CORNER16_T0_V,"C16 T0 V"},
	{FIGGRID_CORNER16_T1_H,"C16 T1 H"},
	{FIGGRID_CORNER16_T1_V,"C16 T1 V"},
	{FIGGRID_XLOOP3_T0,"xloop3 T0"},
	{FIGGRID_XLOOP3_T1,"xloop3 T1"},
	{FIGGRID_XLOOP4_T0,"xloop4 T0"},
	{FIGGRID_XLOOP4_T1,"xloop4 T1"}
};

char ZQ_SmokeEditing3DGUI::export_fold[FILENAME_MAXLEN] = "export";
char ZQ_SmokeEditing3DGUI::script_filename[FILENAME_MAXLEN] = "script.bat";


float ZQ_SmokeEditing3DGUI::global_move_step = 10;

ZQ_SmokeEditing3DGUI::ZQ_SmokeEditing3DGUI()
{
}

ZQ_SmokeEditing3DGUI::~ZQ_SmokeEditing3DGUI()
{
	for(int i = 0;i < source_figGrids.size();i++)
		delete source_figGrids[i];
	for(int i = 0;i < patch_figGrids.size();i++)
		delete patch_figGrids[i];
}

ZQ_SmokeEditing3DGUI* ZQ_SmokeEditing3DGUI::GetInstance()
{
	if(instance == 0)
		instance = new ZQ_SmokeEditing3DGUI();
	return instance;
}


bool ZQ_SmokeEditing3DGUI::Init(int* argc, char** argv)
{
	glutInit(argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("Smoke Editing 3D GUI");
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
	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &sceneZoom, "visible=true min=0.01 max=100");
	TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &sceneRotation, 
		" label='Object rotation' opened=true help='Change the object orientation.' ");
	TwAddButton(bar,"OriginalView",(TwButtonCallback)OnSetOriginalView,0,"label='Original View'");
	TwAddVarRW(bar,"alpha_blend",TW_TYPE_BOOL8,&alpha_blend,"label='alpha blend'");
	TwAddVarRW(bar,"draw_source_group",TW_TYPE_BOOL8,&draw_source_group,"label='draw source'");
	TwAddVarRW(bar,"draw_source_patch",TW_TYPE_BOOL8,&draw_patch_group,"label='draw patch'");

	TwAddSeparator(bar,"Create Sep","label='Create'");
	TwAddVarRW(bar,"source_or_patch",TW_TYPE_BOOL8,&source_or_patch,"label='mode' readonly=true");
	TwAddButton(bar,"ChangeMode",(TwButtonCallback)OnChangeMode,0,"label='Change Mode'");
	TwType FigGridTWType = TwDefineEnum("FigGridTypeEnum", typeEnum, 16);
	TwAddVarRW(bar, "FigGridType", FigGridTWType, &create_type, "label='type'");
	TwAddButton(bar,"CreateButton",(TwButtonCallback)OnCreate,0,"label='Create'");
	TwAddButton(bar,"DeleteButton",(TwButtonCallback)OnDelete,0,"label='Delete'");

	TwAddSeparator(bar,"Export Sep","label='Export'");
	TwAddVarRW(bar,"export_fold",TW_TYPE_CSSTRING(sizeof(export_fold)),&export_fold,"label='fold'");
	TwAddButton(bar,"Export",(TwButtonCallback)OnExport,0,"label='Export'");

	return true;
}

void ZQ_SmokeEditing3DGUI::MainLoop()
{
	glutMainLoop();
}

void ZQ_SmokeEditing3DGUI::Terminate()
{
	TwTerminate();
}

void ZQ_SmokeEditing3DGUI::Display()
{
	float mat[4*4]; // rotation matrix

	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_NORMALIZE);
	glPushMatrix(); 

	ConvertQuaternionToMatrix(sceneRotation, mat);
	glMultMatrixf(mat);

	glScalef(sceneZoom, sceneZoom, sceneZoom);

	if(alpha_blend)
	{
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE,GL_ONE_MINUS_SRC_ALPHA);
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


void ZQ_SmokeEditing3DGUI::Reshape(int w, int h)
{
	width = w;
	height = h;
	// Set OpenGL viewport and camera
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fovy, (double)width/height, 1, 1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyepos[0],eyepos[1],eyepos[2], 0,0,0, 0,1,0);

	// Send the new window size to AntTweakBar
	TwWindowSize(width, height);
	_resize_bar(width,height);
}


void ZQ_SmokeEditing3DGUI::OnMouseMove(int mx, int my)
{
	
	if(!TwEventMouseMotionGLUT(mx,my))
	{
		float ray_ori[3],ray_dir[3];
		_mouse_pos_to_ray(mx,my,ray_ori,ray_dir);
		if(source_or_patch)
		{
			if(cur_focus_fig_id_source >=0)
			{
				source_figGrids[cur_focus_fig_id_source]->OnMouseMove(ray_ori,ray_dir);
			}
		}
		else
		{
			if(cur_focus_fig_id_patch >= 0)
			{
				patch_figGrids[cur_focus_fig_id_patch]->OnMouseMove(ray_ori,ray_dir);
			}
		}
	}
}
void ZQ_SmokeEditing3DGUI::OnMousePassiveMove(int mx, int my)
{
	if(!TwEventMouseMotionGLUT(mx,my))
	{
		float ray_ori[3],ray_dir[3];
		_mouse_pos_to_ray(mx,my,ray_ori,ray_dir);
		if(source_or_patch)
		{
			if(cur_focus_fig_id_source >= 0)
			{
				source_figGrids[cur_focus_fig_id_source]->OnMousePassiveMove(ray_ori,ray_dir);
			}
		}
		else
		{
			if(cur_focus_fig_id_patch >= 0)
			{
				patch_figGrids[cur_focus_fig_id_patch]->OnMousePassiveMove(ray_ori,ray_dir);
			}
		}

	}
}

void ZQ_SmokeEditing3DGUI::OnMouseButton(int button, int state, int mx, int my)
{
	if(!TwEventMouseButtonGLUT(button,state,mx,my))
	{
		float ray_ori[3],ray_dir[3];
		_mouse_pos_to_ray(mx,my,ray_ori,ray_dir);
		if(source_or_patch)
		{
			if(cur_focus_fig_id_source >= 0)
			{
				source_figGrids[cur_focus_fig_id_source]->OnMouseButton(button,state,ray_ori,ray_dir);
			}
		}
		else
		{
			if(cur_focus_fig_id_patch >= 0)
			{
				patch_figGrids[cur_focus_fig_id_patch]->OnMouseButton(button,state,ray_ori,ray_dir);
			}
		}
	}
}

void ZQ_SmokeEditing3DGUI::OnKeyBoard(int key, int mx, int my)
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
				float ray_ori[3],ray_dir[3];
				_mouse_pos_to_ray(mx,my,ray_ori,ray_dir);
				if(source_or_patch)
				{
					if(cur_focus_fig_id_source >= 0)
					{
						source_figGrids[cur_focus_fig_id_source]->OnKeyboard(key,ray_ori,ray_dir);
					}
				}
				else
				{
					if(cur_focus_fig_id_patch >= 0)
					{
						patch_figGrids[cur_focus_fig_id_patch]->OnKeyboard(key,ray_ori,ray_dir);
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

void ZQ_SmokeEditing3DGUI::OnSpecial(int key, int mx, int my)
{
	if(!TwEventSpecialGLUT(key,mx,my))
	{
	}
}


void ZQ_SmokeEditing3DGUI::OnSetOriginalView(void* clientData)
{
	memcpy(sceneRotation,sceneRotationStart,sizeof(float)*4);
	sceneZoom = 1;

}

void ZQ_SmokeEditing3DGUI::OnKeyUp()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(0,global_move_step,0);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(0,global_move_step,0);
}

void ZQ_SmokeEditing3DGUI::OnKeyDown()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(0,-global_move_step,0);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(0,-global_move_step,0);
}

void ZQ_SmokeEditing3DGUI::OnKeyLeft()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(-global_move_step,0,0);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(-global_move_step,0,0);
}

void ZQ_SmokeEditing3DGUI::OnKeyRight()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->GlobalMove(global_move_step,0,0);
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->GlobalMove(global_move_step,0,0);
}

void ZQ_SmokeEditing3DGUI::OnCreate(void* clientData)
{
	BaseFigGrid3D* fig = 0;
	switch(create_type)
	{
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER8_T0_H:
		{
			fig = new CornersFigGrid3D<2,0,true>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER8_T0_V:
		{
			fig = new CornersFigGrid3D<2,0,false>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER8_T1_H:
		{
			fig = new CornersFigGrid3D<2,1,true>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER8_T1_V:
		{
			fig = new CornersFigGrid3D<2,1,false>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER12_T0_H:
		{
			fig = new CornersFigGrid3D<3,0,true>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER12_T0_V:
		{
			fig = new CornersFigGrid3D<3,0,false>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER12_T1_H:
		{
			fig = new CornersFigGrid3D<3,1,true>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER12_T1_V:
		{
			fig = new CornersFigGrid3D<3,1,false>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER16_T0_H:
		{
			fig = new CornersFigGrid3D<4,0,true>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER16_T0_V:
		{
			fig = new CornersFigGrid3D<4,0,false>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER16_T1_H:
		{
			fig = new CornersFigGrid3D<4,1,true>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_CORNER16_T1_V:
		{
			fig = new CornersFigGrid3D<4,1,false>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_XLOOP3_T0:
		{
			fig = new XloopDeformGrid3D<3,0>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_XLOOP3_T1:
		{
			fig = new XloopDeformGrid3D<3,1>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_XLOOP4_T0:
		{
			fig = new XloopDeformGrid3D<4,0>();
		}
		break;
	case ZQ_SmokeEditing3DGUI::FIGGRID_XLOOP4_T1:
		{
			fig = new XloopDeformGrid3D<4,1>();
		}
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

void ZQ_SmokeEditing3DGUI::OnDelete(void* clientData)
{
	if(source_or_patch)
	{
		if(cur_focus_fig_id_source >= 0)
		{
			BaseFigGrid3D* g = source_figGrids[cur_focus_fig_id_source];
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
			BaseFigGrid3D* g = patch_figGrids[cur_focus_fig_id_patch];
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


void ZQ_SmokeEditing3DGUI::OnNext(void*clientData)
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

void ZQ_SmokeEditing3DGUI::OnPrevious(void* clientData)
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

void ZQ_SmokeEditing3DGUI::OnChangeMode(void* clientData)
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

void ZQ_SmokeEditing3DGUI::Draw()
{
	if(draw_source_group)
	{	
		for(int i = 0;i < source_figGrids.size();i++)
			source_figGrids[i]->Draw();
	}
	if(draw_patch_group)
	{
		for(int i = 0;i < patch_figGrids.size();i++)
			patch_figGrids[i]->Draw();
	}
}

void ZQ_SmokeEditing3DGUI::UpdataPerFrame()
{
	for(int i = 0;i < source_figGrids.size();i++)
		source_figGrids[i]->UpdatePerFrame();
	for(int i = 0;i < patch_figGrids.size();i++)
		patch_figGrids[i]->UpdatePerFrame();
}

void ZQ_SmokeEditing3DGUI::_resize_bar(int w,int h)
{
	char buf[500];
	int bar_width = 200;
	sprintf(buf,"Menu size=\'%d %d\' position=\'%d %d\'",bar_width,500,w-bar_width,0);
	TwDefine(buf);
}

void ZQ_SmokeEditing3DGUI::_mouse_pos_to_ray(int mx, int my, float ray_ori[3], float ray_dir[3])
{
	float cx = width/2;
	float cy = height/2;
	double m_pi = 4*atan(1.0);
	float focal_len = height/2.0/tan(fovy/2.0/180.0*m_pi);
	float view_ray_ori[4] = {0,0,0,1};
	float view_ray_dir[4] = {mx-cx,height-1-my-cy,-focal_len,1};
	double len = sqrt((double)view_ray_dir[0]*view_ray_dir[0]+view_ray_dir[1]*view_ray_dir[1]+view_ray_dir[2]*view_ray_dir[2]);
	view_ray_dir[0] /= len;
	view_ray_dir[1] /= len;
	view_ray_dir[2] /= len;

	float mat_view_to_world[16] = 
	{
		0,0,0,eyepos[0],
		0,0,0,eyepos[1],
		0,0,0,eyepos[2],
		0,0,0,1
	};
	float mat_world_rot[16];
	ConvertQuaternionToMatrix(sceneRotation,mat_world_rot);

	float mat_world_rot_back[16];
	ZQ_MathBase::MatrixInverse(mat_world_rot,4,mat_world_rot_back);

	float total_mat[16];
	ZQ_MathBase::MatrixMul(mat_world_rot,mat_view_to_world,4,4,4,total_mat);

	float ray_ori_4[4],ray_dir_4[4];
	ZQ_MathBase::MatrixMul(total_mat,view_ray_ori,4,4,1,ray_ori_4);
	ZQ_MathBase::MatrixMul(mat_world_rot,view_ray_dir,4,4,1,ray_dir_4);

	memcpy(ray_ori,ray_ori_4,sizeof(float)*3);
	memcpy(ray_dir,ray_dir_4,sizeof(float)*3);

}

void ZQ_SmokeEditing3DGUI::OnExport(void* clientData)
{
	char cmdbuf[2000];
	std::string cmdstr;

	char* source_fold = "source";

	sprintf(cmdbuf,"mkdir %s",export_fold);
	system(cmdbuf);

	
	sprintf(cmdbuf,"mkdir %s\\%s",export_fold,source_fold);
	system(cmdbuf);
	

	bool wrong_flag = false;


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
		if(!source_figGrids[0]->Export(cur_fold.c_str(),cur_script.c_str(),0,0,0))
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
			if(!source_figGrids[i]->Export(cur_fold.c_str(),cur_script.c_str(),0,0,0))
			{
				wrong_flag = true;
				printf("failed to export %s \n","source");
			}
		}
	}
}


void ZQ_SmokeEditing3DGUI::SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat)
{
	float sina2, norm;
	sina2 = (float)sin(0.5f * angle);
	norm = (float)sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
	quat[0] = sina2 * axis[0] / norm;
	quat[1] = sina2 * axis[1] / norm;
	quat[2] = sina2 * axis[2] / norm;
	quat[3] = (float)cos(0.5f * angle);
	float len = sqrt(quat[0]*quat[0] + quat[1]*quat[1] + quat[2]*quat[2] + quat[3]*quat[3]);
	quat[0] /= len;
	quat[1] /= len;
	quat[2] /= len;
	quat[3] /= len;
}

void ZQ_SmokeEditing3DGUI::ConvertQuaternionToMatrix(const float *quat, float *mat)
{
	float yy2 = 2.0f * quat[1] * quat[1];
	float xy2 = 2.0f * quat[0] * quat[1];
	float xz2 = 2.0f * quat[0] * quat[2];
	float yz2 = 2.0f * quat[1] * quat[2];
	float zz2 = 2.0f * quat[2] * quat[2];
	float wz2 = 2.0f * quat[3] * quat[2];
	float wy2 = 2.0f * quat[3] * quat[1];
	float wx2 = 2.0f * quat[3] * quat[0];
	float xx2 = 2.0f * quat[0] * quat[0];
	mat[0*4+0] = - yy2 - zz2 + 1.0f;
	mat[0*4+1] = xy2 + wz2;
	mat[0*4+2] = xz2 - wy2;
	mat[0*4+3] = 0;
	mat[1*4+0] = xy2 - wz2;
	mat[1*4+1] = - xx2 - zz2 + 1.0f;
	mat[1*4+2] = yz2 + wx2;
	mat[1*4+3] = 0;
	mat[2*4+0] = xz2 + wy2;
	mat[2*4+1] = yz2 - wx2;
	mat[2*4+2] = - xx2 - yy2 + 1.0f;
	mat[2*4+3] = 0;
	mat[3*4+0] = mat[3*4+1] = mat[3*4+2] = 0;
	mat[3*4+3] = 1;
}

void ZQ_SmokeEditing3DGUI::MultiplyQuaternions(const float *q1, const float *q2, float *qout)
{
	float qr[4];
	qr[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
	qr[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
	qr[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
	qr[3]  = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
	qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}