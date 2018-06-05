#include "ZQ_SmokeGUI.h"
#include "ZQ_DoubleImage3D.h"
#include "ZQ_MathBase.h"
#include "ZQ_Wavelet.h"
#include "ZQ_CompressedImage.h"

using namespace ZQ;

ZQ_SmokeGUI* ZQ_SmokeGUI::instance = 0;
int ZQ_SmokeGUI::width = 1024;
int ZQ_SmokeGUI::height = 768;
float ZQ_SmokeGUI::fovy = 40;
float ZQ_SmokeGUI::eyepos[3] = {0,0,2};
float ZQ_SmokeGUI::sceneZoom = 1.0f;
float ZQ_SmokeGUI::box_xlen = 1.0f;
float ZQ_SmokeGUI::box_ylen = 1.0f;
float ZQ_SmokeGUI::box_zlen = 1.0f;
float ZQ_SmokeGUI::sceneRotation[4] = {0,0,0,1};//{-0.707106781186548,0,0,0.707106781186548};
float ZQ_SmokeGUI::sceneRotationStart[4] = {0,0,0,1};//{-0.707106781186548,0,0,0.707106781186548};
TwBar* ZQ_SmokeGUI::bar = 0;
TwBar* ZQ_SmokeGUI::bar2 = 0;
bool ZQ_SmokeGUI::drawBox = true;
bool ZQ_SmokeGUI::drawSphere = true;


ZQ_GLSLshader* ZQ_SmokeGUI::glslShaders = 0;
bool ZQ_SmokeGUI::drawVolumeDataFlag = false;
GLuint ZQ_SmokeGUI::volumeTex = 0;
char ZQ_SmokeGUI::volumeDataFile[200] = ".dat";
bool ZQ_SmokeGUI::hasVolumeTex = false;
int ZQ_SmokeGUI::volumeTexWidth = 96;
int ZQ_SmokeGUI::volumeViewPortWidth = ZQ_SmokeGUI::width;
int ZQ_SmokeGUI::volumeViewPortHeight = ZQ_SmokeGUI::height;
int ZQ_SmokeGUI::volumeStep = 300;

GLuint ZQ_SmokeGUI::program = 0;
GLuint ZQ_SmokeGUI::caminfoloc = 0;
GLuint ZQ_SmokeGUI::invViewloc = 0;
GLuint ZQ_SmokeGUI::boxMinloc = 0;
GLuint ZQ_SmokeGUI::boxMaxloc = 0;
GLuint ZQ_SmokeGUI::boxSizeloc = 0;
GLuint ZQ_SmokeGUI::eyePosloc = 0;
GLuint ZQ_SmokeGUI::steploc = 0;
GLuint ZQ_SmokeGUI::texloc = 0;
GLuint ZQ_SmokeGUI::densityScaleloc = 0;
GLuint ZQ_SmokeGUI::redloc = 0;
GLuint ZQ_SmokeGUI::greenloc = 0;
GLuint ZQ_SmokeGUI::blueloc = 0;
float ZQ_SmokeGUI::drawDensityScale = 1.0f;

int ZQ_SmokeGUI::renderRed = 255;
int ZQ_SmokeGUI::renderGreen = 255;
int ZQ_SmokeGUI::renderBlue = 255;


char ZQ_SmokeGUI::sequenceFold[200] = "data0";
char ZQ_SmokeGUI::prefix[200] = "frame";
char ZQ_SmokeGUI::suffix[200] = ".di3";
int ZQ_SmokeGUI::sequenceNum = 100;
int ZQ_SmokeGUI::sequenceBaseID = 0;
int ZQ_SmokeGUI::sequenceDataWidth = 96;
bool ZQ_SmokeGUI::isRenderingSequence = false;
float* ZQ_SmokeGUI::sequenceBuffer[SEQUENCE_BUFFER_NUM] = {0};
int ZQ_SmokeGUI::produceid = 0;
int ZQ_SmokeGUI::consumeid = 0;
bool ZQ_SmokeGUI::producestop = false;
bool ZQ_SmokeGUI::consumestop = false;


const float BOX_LINE_WIDTH = 3.0f;
const float BOX_LINE_COLOR[3] = {0,1,0};

const float LINE_WIDTH = 3.0f;
const float NURBS_WIDTH = 2.0f;
const float CUR_SELECT_POINTSIZE = 8.0f;
const float POINTSIZE = 5.0f;
const float CUR_CREATE_POINTSIZE = 3.0f;
const float CUR_CREAT_POINT_COLOR[3] = {1.0,0.0,0.0};
const float CUR_SELECT_COLOR[3] = {1.0,0.0,0.0};
const float CUR_LINE_COLOR[3] = {0.4,1.0,0.0};
const float SAVED_LINES_COLOR[3] = {1.0,0.4,0.0};
const float NURBS_LINES_COLOR[3] = {1.0,0.0,1.0};


ZQ_SmokeGUI::ZQ_SmokeGUI()
{
	
}

ZQ_SmokeGUI::~ZQ_SmokeGUI()
{
}
ZQ_SmokeGUI* ZQ_SmokeGUI::GetInstance()
{
	if(instance == 0)
		instance = new ZQ_SmokeGUI();
	return instance;
}
bool ZQ_SmokeGUI::Init(int* argc, char** argv)
{


	glutInit(argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(ZQ_SmokeGUI::width, ZQ_SmokeGUI::height);
    glutCreateWindow("Smoke GUI");
    glutCreateMenu(NULL);
	glewInit();

	if(!InitShaders())
		printf("init shader fail\n");


	glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    atexit(Terminate);

	TwInit(TW_OPENGL, NULL);

	glutMouseFunc((GLUTmousebuttonfun)OnMouseButton);
	glutMotionFunc((GLUTmousemotionfun)OnMouseMove);
	glutPassiveMotionFunc((GLUTmousemotionfun)OnMousePassiveMove);
	glutKeyboardFunc((GLUTkeyboardfun)OnKeyBoard);
	glutSpecialFunc((GLUTspecialfun)OnSpecial);

//	float axis[3] = {1,0,0};
//	SetQuaternionFromAxisAngle(axis,-M_PI*90.0/180.0,sceneRotationStart);
//	memcpy(sceneRotation,sceneRotationStart,sizeof(float)*3);

	int opened = 0;
	TwGLUTModifiersFunc(glutGetModifiers);

	bar = TwNewBar("Menu");
	TwDefine("Menu refresh=1");
    TwDefine(" GLOBAL help='Smoke GUI' "); // Message added to the help bar.
    TwDefine(" Menu size='200 600' color='96 96 96' alpha='80' "); // change default tweak bar size and color

	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &sceneZoom, 
          " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' visible=true");
	TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &sceneRotation, 
               " label='Object rotation' opened=true help='Change the object orientation.' ");
	TwAddButton(bar,"OriginalView",(TwButtonCallback)OnSetOriginalView,0,"label='Original View'");
	TwAddVarRW(bar,"DrawBox",TW_TYPE_BOOL8,&drawBox,"label='DrawBox'");
	TwAddVarRW(bar,"DrawSphere",TW_TYPE_BOOL8,&drawSphere,"label='DrawSphere'");
	

	bar2 = TwNewBar("ImportExport");
	TwDefine("ImportExport refresh=1");
    TwDefine(" ImportExport size='200 600' color='96 96 96' alpha='80' position = '800 24'"); // change default tweak bar size and color
	

	TwAddVarRW(bar2,"VolumeDataFile",TW_TYPE_CSSTRING(sizeof(volumeDataFile)),&volumeDataFile,"label='File' group='ImportVolumeData'");
	TwAddVarRW(bar2,"VolumeDataWidth",TW_TYPE_INT32,&volumeTexWidth,"label='Width' min=0 max= 512 group='ImportVolumeData'");
	TwAddButton(bar2,"LoadVolumeData",(TwButtonCallback)OnLoadVolumeData,0,"label='Load' group='ImportVolumeData'");
	opened = 0;
	TwSetParam(bar2,"ImportVolumeData","opened",TW_PARAM_INT32,1,&opened);

	TwAddVarRW(bar2,"DrawVolumeData",TW_TYPE_BOOL8,&drawVolumeDataFlag,"label='DrawVolumeData'");
	TwAddVarRW(bar2,"DensityScale",TW_TYPE_FLOAT,&drawDensityScale,"label='DrawDensityScale' min=0.01 max=100.0 step = 1");
	TwAddVarRW(bar2,"RenderRed",TW_TYPE_INT32,&renderRed,"label='R' min=0 max=255 step=1");
	TwAddVarRW(bar2,"RenderGreen",TW_TYPE_INT32,&renderGreen,"label='G' min=0 max=255 step=1");
	TwAddVarRW(bar2,"RenderBlue",TW_TYPE_INT32,&renderBlue,"label='B' min=0 max=255 step=1");

	TwAddVarRW(bar2,"SeqFold",TW_TYPE_CSSTRING(sizeof(sequenceFold)),&sequenceFold,"label='SeqFold' group='RenderSequence'");
	TwAddVarRW(bar2,"SeqPrefix",TW_TYPE_CSSTRING(sizeof(prefix)),&prefix,"label='Prefix' group='RenderSequence'");
	TwAddVarRW(bar2,"SeqSuffix",TW_TYPE_CSSTRING(sizeof(suffix)),&suffix,"label='Suffix' group='RenderSequence'");
	TwAddVarRW(bar2,"SeqNum",TW_TYPE_INT32,&sequenceNum,"label='SeqNum' min=1 max=10000 group='RenderSequence'");
	TwAddVarRW(bar2,"SeqBaseID",TW_TYPE_INT32,&sequenceBaseID,"label='BaseID' min=0 max=10000 group='RenderSequence'");
	TwAddVarRW(bar2,"SeqDataWidth",TW_TYPE_INT32,&sequenceDataWidth,"label='Width' min=0 max=512 group='RenderSequence'");
	TwAddButton(bar2,"RenderSequenceButton",(TwButtonCallback)OnRenderSequence,0,"label='Render' group='RenderSequence'");
	opened = 0;
	TwSetParam(bar2,"RenderSequence","opened",TW_PARAM_INT32,1,&opened);

	return true;
}

void ZQ_SmokeGUI::MainLoop()
{
	glutMainLoop();
}

void ZQ_SmokeGUI::Terminate()
{
	TwTerminate();
	if(glslShaders)
		delete glslShaders;
	
}

void ZQ_SmokeGUI::Display()
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

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE,GL_SRC_ALPHA);
	
	if(drawVolumeDataFlag)
	{
		DrawVolumeData();
	}
	if(drawSphere)
	{
		DrawSphere();
	}
	
	glDisable(GL_BLEND);
	
	if(drawBox)
		DrawBox();
	
	glEnable(GL_DEPTH_TEST);
	glPopMatrix();

    //// Draw tweak bars
    TwDraw();
	
    // Present frame buffer
    glutSwapBuffers();

    // Recall Display at next frame
    glutPostRedisplay();

	Sleep(10);
}
void ZQ_SmokeGUI::Reshape(int width, int height)
{
	ZQ_SmokeGUI::width = width;
	ZQ_SmokeGUI::height = height;
	ZQ_SmokeGUI::volumeViewPortWidth = width;
	ZQ_SmokeGUI::volumeViewPortHeight = height;
	// Set OpenGL viewport and camera
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluPerspective(fovy, (double)width/height, 0.1, 10);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eyepos[0],eyepos[1],eyepos[2], 0,0,0, 0,1,0);
    
    // Send the new window size to AntTweakBar
    TwWindowSize(width, height);
}

void ZQ_SmokeGUI::OnMouseMove(int mx, int my)
{
	if(!TwEventMouseMotionGLUT(mx,my))
	{
	}
}
void ZQ_SmokeGUI::OnMousePassiveMove(int mx, int my)
{
	if(!TwEventMouseMotionGLUT(mx,my))
	{
	}
}

void ZQ_SmokeGUI::OnMouseButton(int button, int state, int mx, int my)
{
	if(!TwEventMouseButtonGLUT(button,state,mx,my))
	{
	}
}

void ZQ_SmokeGUI::OnKeyBoard(int key, int mx, int my)
{
	if(!TwEventKeyboardGLUT(key,mx,my))
	{
	}
	//printf("key = %d\n", key);
	if(key == TW_KEY_ESCAPE) exit(0);
	//if (key == 'q' || key == 'Q') exit(0);
}
void ZQ_SmokeGUI::OnSpecial(int key, int mx, int my)
{
	if(!TwEventSpecialGLUT(key,mx,my))
	{
	}
}

void ZQ_SmokeGUI::OnSetOriginalView(void* clientData)
{
	memcpy(sceneRotation,sceneRotationStart,sizeof(float)*4);
	sceneZoom = 1;
	
}


void ZQ_SmokeGUI::OnLoadVolumeData(void* clientData)
{
	if(isRenderingSequence)
	{
		printf("rendering sequence... push this button later!\n");
		return ;
	}
	if(!LoadDataToVolumeTex(volumeDataFile,volumeTexWidth))
		printf("load %s fail\n",volumeDataFile);
}


void ZQ_SmokeGUI::OnRenderSequence(void *clientData)
{

	isRenderingSequence = !isRenderingSequence;
	if(isRenderingSequence)
	{
		TwSetParam(bar2,"SeqNum","readonly",TW_PARAM_CSTRING,1,"true");
		TwSetParam(bar2,"SeqFold","readonly",TW_PARAM_CSTRING,1,"true");
		TwSetParam(bar2,"SeqPrefix","readonly",TW_PARAM_CSTRING,1,"true");
		TwSetParam(bar2,"SeqSuffix","readonly",TW_PARAM_CSTRING,1,"true");
		TwSetParam(bar2,"SeqBaseID","readonly",TW_PARAM_CSTRING,1,"true");
		TwSetParam(bar2,"SeqDataWidth","readonly",TW_PARAM_CSTRING,1,"true");
		TwSetParam(bar2,"RenderSequenceButton","label",TW_PARAM_CSTRING,1,"Stop Render");
	}
	else
	{
		TwSetParam(bar2,"SeqNum","readonly",TW_PARAM_CSTRING,1,"false");
		TwSetParam(bar2,"SeqFold","readonly",TW_PARAM_CSTRING,1,"false");
		TwSetParam(bar2,"SeqPrefix","readonly",TW_PARAM_CSTRING,1,"false");
		TwSetParam(bar2,"SeqSuffix","readonly",TW_PARAM_CSTRING,1,"false");
		TwSetParam(bar2,"SeqBaseID","readonly",TW_PARAM_CSTRING,1,"false");
		TwSetParam(bar2,"SeqDataWidth","readonly",TW_PARAM_CSTRING,1,"false");
		TwSetParam(bar2,"RenderSequenceButton","label",TW_PARAM_CSTRING,1,"Render");
	}
}


void ZQ_SmokeGUI::DrawBox()
{
	float pts[8][3] = 
	{
		{-0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen},
		{-0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen},
		{-0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen},
		{-0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen},
		{ 0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen},
		{ 0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen},
		{ 0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen},
		{ 0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen}
	};
	glLineWidth(BOX_LINE_WIDTH);
	glColor3f(BOX_LINE_COLOR[0],BOX_LINE_COLOR[1],BOX_LINE_COLOR[2]);
	glBegin(GL_LINES);
	glVertex3f(pts[0][0],pts[0][1],pts[0][2]);
	glVertex3f(pts[1][0],pts[1][1],pts[1][2]);
	glVertex3f(pts[0][0],pts[0][1],pts[0][2]);
	glVertex3f(pts[2][0],pts[2][1],pts[2][2]);
	glVertex3f(pts[3][0],pts[3][1],pts[3][2]);
	glVertex3f(pts[1][0],pts[1][1],pts[1][2]);
	glVertex3f(pts[3][0],pts[3][1],pts[3][2]);
	glVertex3f(pts[2][0],pts[2][1],pts[2][2]);
	
	glVertex3f(pts[4][0],pts[4][1],pts[4][2]);
	glVertex3f(pts[5][0],pts[5][1],pts[5][2]);
	glVertex3f(pts[4][0],pts[4][1],pts[4][2]);
	glVertex3f(pts[6][0],pts[6][1],pts[6][2]);
	glVertex3f(pts[7][0],pts[7][1],pts[7][2]);
	glVertex3f(pts[5][0],pts[5][1],pts[5][2]);
	glVertex3f(pts[7][0],pts[7][1],pts[7][2]);
	glVertex3f(pts[6][0],pts[6][1],pts[6][2]);

	glVertex3f(pts[0][0],pts[0][1],pts[0][2]);
	glVertex3f(pts[4][0],pts[4][1],pts[4][2]);
	glVertex3f(pts[1][0],pts[1][1],pts[1][2]);
	glVertex3f(pts[5][0],pts[5][1],pts[5][2]);
	glVertex3f(pts[2][0],pts[2][1],pts[2][2]);
	glVertex3f(pts[6][0],pts[6][1],pts[6][2]);
	glVertex3f(pts[3][0],pts[3][1],pts[3][2]);
	glVertex3f(pts[7][0],pts[7][1],pts[7][2]);

	glEnd();
}

void ZQ_SmokeGUI::DrawSphere()
{
	GLUquadricObj* quadratic = gluNewQuadric();
	glColor3f(0.1,0.3,0.1);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	float l_amb[4] = {0.2,0.6,0.2,1.0};
	float l_spe[4] = {0.0,0.0,0.0,1.0};
	float l_dif[4] = {0.1,0.3,0.1,1.0};
	float l_shine = 90;
	float l_pos0[4] = {0,-1,0,0};
	float l_pos1[4] = {0,1,0,0};
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0,GL_AMBIENT,l_amb);
	glLightfv(GL_LIGHT0,GL_SPECULAR,l_spe);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,l_dif);
	glLightf(GL_LIGHT0,GL_SHININESS,l_shine);
	glLightfv(GL_LIGHT0,GL_POSITION,l_pos0);
	glLightfv(GL_LIGHT1,GL_AMBIENT,l_amb);
	glLightfv(GL_LIGHT1,GL_SPECULAR,l_spe);
	glLightfv(GL_LIGHT1,GL_DIFFUSE,l_dif);
	glLightf(GL_LIGHT1,GL_SHININESS,l_shine);
	glLightfv(GL_LIGHT1,GL_POSITION,l_pos1);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	
	/*GLfloat mat_specular[] = { 0.2, 0.6, 0.2, 1.0 };
	GLfloat mat_shininess[] = { 0 };
	GLfloat light_position[] = { 1, -1.0, 1.0, 0.0 };
	glShadeModel (GL_SMOOTH);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);*/

	gluSphere(quadratic, 0.1, 50, 50);
	glDisable(GL_LIGHTING);
	delete quadratic;

}

void ZQ_SmokeGUI::SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat)
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

void ZQ_SmokeGUI::ConvertQuaternionToMatrix(const float *quat, float *mat)
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

void ZQ_SmokeGUI::MultiplyQuaternions(const float *q1, const float *q2, float *qout)
{
    float qr[4];
	qr[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
	qr[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
	qr[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
	qr[3]  = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
    qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}

bool ZQ_SmokeGUI::InitShaders()
{
	glslShaders = new ZQ_GLSLshader();
	if(!glslShaders->CreateShaderFromBytes("raycast",vertexShader,pixelShader))
	{
		delete glslShaders;
		glslShaders = 0;
		return false;
	}
	program = glslShaders->GetShaderHandler("raycast");
	
	caminfoloc = glGetUniformLocation(program,"caminfo");
//	printf("caminfoloc = %d\n",caminfoloc);
	invViewloc = glGetUniformLocation(program,"invView");
//	printf("invViewloc = %d\n",invViewloc);
	boxMinloc = glGetUniformLocation(program,"boxMin");
//	printf("boxMinloc = %d\n",boxMinloc);
	boxMaxloc = glGetUniformLocation(program,"boxMax");
//	printf("boxMaxloc = %d\n",boxMaxloc);
	boxSizeloc = glGetUniformLocation(program,"boxSize");
//	printf("boxSizeloc = %d\n",boxSizeloc);
	eyePosloc = glGetUniformLocation(program,"eyePos");
//	printf("eyePosloc = %d\n",eyePosloc);
	steploc = glGetUniformLocation(program,"step");
//	printf("steploc = %d\n",steploc);
	texloc = glGetUniformLocation(program,"tex");
//	printf("texloc = %d\n",texloc);
	densityScaleloc = glGetUniformLocation(program,"densityScale");
//	printf("densityScale = %d\n",densityScaleloc);

	redloc = glGetUniformLocation(program,"red");
	greenloc = glGetUniformLocation(program,"green");
	blueloc = glGetUniformLocation(program,"blue");


	return true;

}

bool ZQ_SmokeGUI::LoadDataToVolumeTex(const char* file,int N)
{
	ZQ_DImage3D<float> img;
	if(_is_file_extension(file,".di3"))
	{
		if(!img.loadImage(file))
			return false;

		img.collapse();

		int width = img.width();
		int height = img.height();
		int depth = img.depth();

		box_xlen = (float)width/height*box_ylen;
		box_zlen = (float)depth/height*box_ylen;
		
		float*& data = img.data();

		if(hasVolumeTex)
			glDeleteTextures(1,&volumeTex);
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		glGenTextures(1,&volumeTex);

		glBindTexture(GL_TEXTURE_3D,volumeTex);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,width,height,depth,0,GL_RED,GL_FLOAT,data);
		glBindTexture(GL_TEXTURE_3D,0);
		volumeTexWidth = width;
		hasVolumeTex = true;

		return true;

	}
	else if(_is_file_extension(file,".zqci"))
	{
		if(!ZQ_CompressedImage::LoadCompressedImage(file,img))
			return false;

		img.collapse();

		int width = img.width();
		int height = img.height();
		int depth = img.depth();

		box_xlen = (float)width/height*box_ylen;
		box_zlen = (float)depth/height*box_ylen;

		float*& data = img.data();

		if(hasVolumeTex)
			glDeleteTextures(1,&volumeTex);
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		glGenTextures(1,&volumeTex);

		glBindTexture(GL_TEXTURE_3D,volumeTex);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,width,height,depth,0,GL_RED,GL_FLOAT,data);
		glBindTexture(GL_TEXTURE_3D,0);
		volumeTexWidth = width;
		hasVolumeTex = true;

		return true;
	}
	else if(_is_file_extension(file,".zqwv"))
	{
		if(!ZQ_Wavelet<float>::LoadWaveletFile(file,img))
			return false;

		img.collapse();

		int width = img.width();
		int height = img.height();
		int depth = img.depth();

		box_xlen = (float)width/height*box_ylen;
		box_zlen = (float)depth/height*box_ylen;

		float*& data = img.data();

		if(hasVolumeTex)
			glDeleteTextures(1,&volumeTex);
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		glGenTextures(1,&volumeTex);

		glBindTexture(GL_TEXTURE_3D,volumeTex);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,width,height,depth,0,GL_RED,GL_FLOAT,data);
		glBindTexture(GL_TEXTURE_3D,0);
		volumeTexWidth = width;
		hasVolumeTex = true;

		return true;
	}
	else
	{
		FILE* in = fopen(file,"rb");
		if(in == 0)
			return false;
		float* tmpData = new float[N*N*N];
		if(fread(tmpData,sizeof(float),N*N*N,in) != N*N*N)
		{
			delete []tmpData;
			fclose(in);
			return false;
		}
		fclose(in);

		float* tmpdata2 = new float[N*N*N];
		for(int i = 0;i < N;i++)
		{
			for(int j = 0;j < N;j++)
			{
				for(int k = 0;k < N;k++)
				{
					tmpdata2[i*N*N+j*N+k] = tmpData[k*N*N+j*N+i];
				}
			}
		}
		/*printf("N=%d\n",N);*/

		//	EnterCriticalSection(&volumeTexMutex);

		if(hasVolumeTex)
			glDeleteTextures(1,&volumeTex);
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		glGenTextures(1,&volumeTex);

		glBindTexture(GL_TEXTURE_3D,volumeTex);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,N,N,N,0,GL_RED,GL_FLOAT,tmpdata2);
		glBindTexture(GL_TEXTURE_3D,0);
		volumeTexWidth = N;
		hasVolumeTex = true;

		//	LeaveCriticalSection(&volumeTexMutex);
		delete []tmpdata2;
		delete []tmpData;
		return true;
	}
}

void ZQ_SmokeGUI::DrawVolumeData()
{
	if(glslShaders == 0)
		return ;
	static int seq = 0;
	if(isRenderingSequence)
	{
		char buff[200];
		sprintf(buff,"%s\\%s%d%s",sequenceFold,prefix,sequenceBaseID+seq,suffix);
		LoadDataToVolumeTex(buff,sequenceDataWidth);
		seq++;
		if(seq >= sequenceNum)
		{
			isRenderingSequence = false;
			TwSetParam(bar2,"SeqNum","readonly",TW_PARAM_CSTRING,1,"false");
			TwSetParam(bar2,"SeqFold","readonly",TW_PARAM_CSTRING,1,"false");
			TwSetParam(bar2,"SeqPrefix","readonly",TW_PARAM_CSTRING,1,"false");
			TwSetParam(bar2,"SeqSuffix","readonly",TW_PARAM_CSTRING,1,"false");
			TwSetParam(bar2,"SeqFold","readonly",TW_PARAM_CSTRING,1,"false");
			TwSetParam(bar2,"SeqBaseID","readonly",TW_PARAM_CSTRING,1,"false");
			TwSetParam(bar2,"SeqDataWidth","readonly",TW_PARAM_CSTRING,1,"false");
			TwSetParam(bar2,"RenderSequenceButton","label",TW_PARAM_CSTRING,1,"Render");
		}
	}
	else 
	{
		seq = 0;
		if(!hasVolumeTex)
			return ;
	}

	glslShaders->UseShader("raycast");

	float vfocal = volumeViewPortHeight/2.0/atan(fovy/360.0f*M_PI);
	glUniform3f(caminfoloc,float(volumeViewPortWidth),float(volumeViewPortHeight),vfocal);
	float boxmin[3] = {-0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen};
	float boxmax[3] = { 0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen};
	float boxsize[3] = {box_xlen,box_ylen,box_zlen};

	glUniform3f(boxMinloc,boxmin[0],boxmin[1],boxmin[2]);
	glUniform3f(boxMaxloc,boxmax[0],boxmax[1],boxmax[2]);
	glUniform3f(boxSizeloc,boxsize[0],boxsize[1],boxsize[2]);

	glUniform1i(steploc,volumeStep);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_3D,volumeTex);
	glUniform1i(texloc,0);
	glUniform1f(densityScaleloc,drawDensityScale);
	glUniform1f(redloc,renderRed/(255.0f));
	glUniform1f(greenloc,renderGreen/(255.0f));
	glUniform1f(blueloc,renderBlue/(255.0f));

	float tmpmat[16],viewmat[16];
	float invView[16];
	glGetFloatv(GL_MODELVIEW_MATRIX,tmpmat);
	ZQ_MathBase::MatrixTranspose(tmpmat,4,4,viewmat);
	/*printf("view\n");
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
			printf("%12.5f",viewmat[i*4+j]);
		printf("\n");
	}*/
	
	ZQ_MathBase::MatrixInverse(viewmat,4,invView);
	
	float eyeposori[4] = {0,0,0,1};
	float eyeposaft[4];
	ZQ_MathBase::MatrixMul(invView,eyeposori,4,4,1,eyeposaft);

	float invView3[9] = 
	{
		invView[0],invView[1],invView[2],
		invView[4],invView[5],invView[6],
		invView[8],invView[9],invView[10]
	};

	glUniformMatrix3fv(invViewloc,1,GL_TRUE,invView3);
	glUniform3fv(eyePosloc,1,eyeposaft);
	

	
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glBegin(GL_QUADS);
	glVertex3f(-0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen);
	glVertex3f(-0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen);
	glVertex3f( 0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen);
	glVertex3f( 0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen);

	glVertex3f(-0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen);
	glVertex3f(-0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen);

	glVertex3f(-0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen);
	glVertex3f(-0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen);

	glVertex3f(-0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen);
	glVertex3f(-0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen);

	glVertex3f(-0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen);
	glVertex3f(-0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen);
	glVertex3f(-0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen);
	glVertex3f(-0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen);

	glVertex3f( 0.5*box_xlen,-0.5*box_ylen,-0.5*box_zlen);
	glVertex3f( 0.5*box_xlen,-0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen, 0.5*box_ylen, 0.5*box_zlen);
	glVertex3f( 0.5*box_xlen, 0.5*box_ylen,-0.5*box_zlen);

	glEnd();

//	printf("draw volume\n");
	glslShaders->UseShader(0);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	glBindTexture(GL_TEXTURE_3D,0);
}

bool ZQ_SmokeGUI::_is_file_extension(const char* filename, const char* extensionname)
{
	if(filename == 0 || extensionname == 0)
		return false;

	int len1 = strlen(filename);
	int len2 = strlen(extensionname);
	if(len1 < len2)
		return false;

	return _strcmpi(filename+len1-len2,extensionname) == 0;
}
