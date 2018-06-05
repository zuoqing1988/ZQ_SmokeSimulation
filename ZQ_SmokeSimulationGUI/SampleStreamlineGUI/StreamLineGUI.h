#ifndef _ZQ_STREAMLINE_GUI_H_
#define _ZQ_STREAMLINE_GUI_H_

#include <GL/glew.h>

#include "AntTweakBar.h"
#include <GL/glut.h>

#include <math.h>
#include <vector>
#include <iostream>

#include "windows.h"
#include <process.h>


#ifndef M_PI
#define M_PI 3.1415926535
#endif



class StreamLineGUI
{
protected:
	StreamLineGUI();
	~StreamLineGUI();
private:
	static StreamLineGUI* instance;
	
	static int width;
	static int height;
	static float fovy;
	static float eyepos[3];

	static float sceneZoom;
	static float sceneRotation[4];
	static float sceneRotationStart[4];
	
	static TwBar* bar;

	static bool drawBox;
	static bool drawStreamLines;
	static float box_line_width;
	static float box_line_color[3];
	static float line_width;
	static char stream_line_file[200];


	static std::vector<GLfloat*> vertex_lists;
	static std::vector<GLuint*> index_lists;
	static std::vector<GLfloat*> color_lists;
	static std::vector<GLuint> line_nums;

	static unsigned int render_line_num;

	
public:
	static StreamLineGUI* GetInstance();
	bool Init(int* argc, char** argv);
	void MainLoop();

private:
	static void _clearLines();
	static void Terminate(void);
	static void Display(void);
	static void Reshape(int width, int height);
	static void  OnMouseMove(int mx, int my);
	static void  OnMousePassiveMove(int mx, int my);
	static void  OnMouseButton(int button, int state, int mx, int my);
	static void  OnKeyBoard(int key, int mx, int my);
	static void  OnSpecial(int key, int mx, int my);
	static void TW_CALL OnSetOriginalView(void* clientData);
	static void TW_CALL OnLoadFile(void* clientData);
	

	static void DrawBox();
	static void DrawStreamLines();
	

	static void SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat);
	static void ConvertQuaternionToMatrix(const float *quat, float *mat);
	static void MultiplyQuaternions(const float *q1, const float *q2, float *qout);

	

};

void Encode_color_from_vel(float& r, float& g, float& b, float x, float y, float z);

float StreamLineInfoValue(const int vertex_num, const float* points);

void SortStreamLineInfoValue(const int N, int* idx, float* values, bool accend = false);

#endif