#ifndef _ZQ_SMOKE_GUI_H_
#define _ZQ_SMOKE_GUI_H_
#pragma once 

#include "windows.h"
#include "ZQ_GLSLShader.h"
#include "AntTweakBar.h"
#include "raycastingshader.h"
#include <GL/glut.h>
#include <math.h>
#include <vector>
#include <iostream>

#include <process.h>

#ifndef M_PI
#define M_PI 3.1415926535
#endif

#define SEQUENCE_BUFFER_NUM 10

class ZQ_SmokeGUI
{
protected:
	ZQ_SmokeGUI();
	~ZQ_SmokeGUI();
private:
	static ZQ_SmokeGUI* instance;

	static int width;
	static int height;
	static float fovy;
	static float eyepos[3];

	static float sceneZoom;
	static float box_xlen;
	static float box_ylen;
	static float box_zlen;
	static float sceneRotation[4];
	static float sceneRotationStart[4];
	static TwBar* bar;
	static TwBar* bar2;
	static bool drawBox;
	static bool drawSphere;

	static ZQ::ZQ_GLSLshader* glslShaders;
	static bool drawVolumeDataFlag;
	static GLuint volumeTex;
	static char volumeDataFile[200];
	static bool hasVolumeTex;
	static int volumeTexWidth;
	static int volumeViewPortWidth;
	static int volumeViewPortHeight;
	static int volumeStep;

	static GLuint program;
	static GLuint caminfoloc;
	static GLuint invViewloc;
	static GLuint boxMinloc;
	static GLuint boxMaxloc;
	static GLuint boxSizeloc;
	static GLuint eyePosloc;
	static GLuint steploc;
	static GLuint texloc;
	static GLuint densityScaleloc;
	static GLuint redloc;
	static GLuint greenloc;
	static GLuint blueloc;
	static float drawDensityScale;
	static int renderRed;
	static int renderGreen;
	static int renderBlue;


	static char sequenceFold[200];
	static char prefix[200];
	static char suffix[200];
	static int sequenceNum;
	static int sequenceBaseID;
	static int sequenceDataWidth;
	static bool isRenderingSequence;
	static float* sequenceBuffer[SEQUENCE_BUFFER_NUM];
	static int produceid;
	static int consumeid;
	static bool producestop;
	static bool consumestop;
	
	static bool exportGuide;

	static char scriptFile[200];
	static char initDensityFile[200];
	static char storeFold[200];


public:
	static ZQ_SmokeGUI* GetInstance();
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

	static void TW_CALL OnLoadVolumeData(void* clientData);

	static void TW_CALL OnRenderSequence(void* clientData);

	static void DrawBox();
	static void DrawSphere();


	static void SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat);
	static void ConvertQuaternionToMatrix(const float *quat, float *mat);
	static void MultiplyQuaternions(const float *q1, const float *q2, float *qout);


	static bool InitShaders();
	static bool LoadDataToVolumeTex(const char* file, int N);
	static void DrawVolumeData();

	static bool _is_file_extension(const char* filename, const char* extensionname);

};

#endif