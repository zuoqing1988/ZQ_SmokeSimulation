#include "StreamLineGUI.h"

StreamLineGUI* StreamLineGUI::instance = 0;
int StreamLineGUI::width = 1024;
int StreamLineGUI::height = 768;
float StreamLineGUI::fovy = 40;
float StreamLineGUI::eyepos[3] = {0,0,2};
float StreamLineGUI::sceneZoom = 1.0f;
float StreamLineGUI::sceneRotation[4] = {0,0,0,1};//{-0.707106781186548,0,0,0.707106781186548};
float StreamLineGUI::sceneRotationStart[4] = {0,0,0,1};//{-0.707106781186548,0,0,0.707106781186548};

TwBar* StreamLineGUI::bar = 0;
bool StreamLineGUI::drawBox = true;
bool StreamLineGUI::drawStreamLines = false;

std::vector<GLfloat*> StreamLineGUI::vertex_lists;
std::vector<GLuint*> StreamLineGUI::index_lists;
std::vector<GLfloat*> StreamLineGUI::color_lists;
std::vector<GLuint> StreamLineGUI::line_nums;

unsigned int StreamLineGUI::render_line_num = 0;

float StreamLineGUI::box_line_width = 3.0f;
float StreamLineGUI::box_line_color[3] = {0,1,0};
float StreamLineGUI::line_width = 1.0f;
char StreamLineGUI::stream_line_file[200] = "line.txt";


StreamLineGUI::StreamLineGUI()
{
	
}

StreamLineGUI::~StreamLineGUI()
{
}

StreamLineGUI* StreamLineGUI::GetInstance()
{
	if(instance == 0)
		instance = new StreamLineGUI();
	return instance;
}
bool StreamLineGUI::Init(int* argc, char** argv)
{
	glutInit(argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(StreamLineGUI::width, StreamLineGUI::height);
    glutCreateWindow("StreamLine GUI");
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

	bar = TwNewBar("EditLines");
	TwDefine("EditLines refresh=1");
    TwDefine(" GLOBAL help='Sketch smoke' "); // Message added to the help bar.
    TwDefine(" EditLines size='200 600' color='96 96 96' alpha='80' "); // change default tweak bar size and color

	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &sceneZoom, 
          " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' visible=false");
	TwAddVarRW(bar, "ObjRotation", TW_TYPE_QUAT4F, &sceneRotation, 
               " label='Object rotation' opened=true help='Change the object orientation.' ");
	TwAddButton(bar,"OriginalView",(TwButtonCallback)OnSetOriginalView,0,"label='Original View'");
	
	TwAddVarRW(bar,"DrawBox",TW_TYPE_BOOL8,&drawBox,"label='DrawBox'");
	TwAddVarRW(bar,"DrawStreamLines",TW_TYPE_BOOL8,&drawStreamLines,"label='DrawLines'");
	TwAddVarRW(bar,"BoxLineWidth",TW_TYPE_FLOAT,&box_line_width,"min=1 max=3 step=1");
	TwAddVarRW(bar,"LineWidth",TW_TYPE_FLOAT,&line_width,"min=1 max=3 step=1");
	TwAddVarRW(bar,"LineFile",TW_TYPE_CSSTRING(sizeof(stream_line_file)),&stream_line_file,"label='LineFile'");
	TwAddButton(bar,"LoadFileButton",(TwButtonCallback)OnLoadFile,0,"label='LoadFile'");
	TwAddVarRW(bar,"RenderLineNum",TW_TYPE_UINT32,&render_line_num,"min=0 max=9999 step = 10");

	
	
	return true;
}

void StreamLineGUI::MainLoop()
{
	glutMainLoop();
}

void StreamLineGUI::_clearLines()
{
	for(int i = 0;i < vertex_lists.size();i++)
	{
		if(vertex_lists[i])
		{
			delete []vertex_lists[i];
			vertex_lists[i] = 0;
		}
		if(color_lists[i])
		{
			delete []color_lists[i];
			color_lists[i] = 0;
		}
		if(index_lists[i])
		{
			delete []index_lists[i];
			index_lists[i] = 0;
		}
	}
	vertex_lists.clear();
	color_lists.clear();
	index_lists.clear();
	line_nums.clear();
}

void StreamLineGUI::Terminate()
{
	TwTerminate();

	_clearLines();

}

void StreamLineGUI::Display()
{
    float mat[4*4]; // rotation matrix
	

	glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix(); 

	glFrontFace(GL_CCW);
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);


	ConvertQuaternionToMatrix(sceneRotation, mat);
	glMultMatrixf(mat);
	glScalef(sceneZoom, sceneZoom, sceneZoom);
	
	if(drawBox)
		DrawBox();

	if(drawStreamLines)
		DrawStreamLines();
	
	
	glPopMatrix();

    //// Draw tweak bars
    TwDraw();
	
    // Present frame buffer
    glutSwapBuffers();

    // Recall Display at next frame
    glutPostRedisplay();
}
void StreamLineGUI::Reshape(int width, int height)
{
	StreamLineGUI::width = width;
	StreamLineGUI::height = height;

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

void StreamLineGUI::OnMouseMove(int mx, int my)
{
	if(!TwEventMouseMotionGLUT(mx,my))
	{
	}
}
void StreamLineGUI::OnMousePassiveMove(int mx, int my)
{
	if(!TwEventMouseMotionGLUT(mx,my))
	{
		
	}
}

void StreamLineGUI::OnMouseButton(int button, int state, int mx, int my)
{
	if(!TwEventMouseButtonGLUT(button,state,mx,my))
	{

	}
}

void StreamLineGUI::OnKeyBoard(int key, int mx, int my)
{
	if(!TwEventKeyboardGLUT(key,mx,my))
	{
	}

	if(key == TW_KEY_ESCAPE) exit(0);
}
void StreamLineGUI::OnSpecial(int key, int mx, int my)
{
	if(!TwEventSpecialGLUT(key,mx,my))
	{
	}
}

void StreamLineGUI::OnSetOriginalView(void* clientData)
{
	memcpy(sceneRotation,sceneRotationStart,sizeof(float)*4);
	sceneZoom = 1;
	
}

void StreamLineGUI::OnLoadFile(void *clientData)
{
	FILE* in = fopen(stream_line_file,"r");
	if(in == 0)
	{
		printf("fail to load file %s\n",stream_line_file);
		return ;
	}

	std::vector<float*> datas;
	std::vector<float*> vels;
	std::vector<int> nums;
	bool has_vel = false;

	bool flag = true;

	int width = 128, height = 128, depth = 128;

	//fscanf(in,"%d%d%d",&width,&height,&depth);

 	char buf[200] = "";
	char str1[200] = "";
	char str2[200] = "";

	memset(buf,0,200);
	char c = ' ';
	do{
		if(feof(in))
			break;
		c = fgetc(in);
		if(c != ' ' && c != '\t' && c != '\n')
			break;
	}while(true);

	if(feof(in))
	{
		printf("unknown format file %s\n",stream_line_file);
		fclose(in);
		return ;
	}

	buf[0] = c;
	fgets(buf+1,198,in);
	sscanf(buf,"%s%s",str1,str2);

	if(strcmp(str1,"format") == 0)
	{
		if(strcmp(str2,"pos") == 0)
		{
			fscanf(in,"%d%d%d",&width,&height,&depth);
		}
		else if(strcmp(str2,"pos+vel") == 0)
		{
			fscanf(in,"%d%d%d",&width,&height,&depth);
			has_vel = true;
		}
		else
		{
			printf("unknown format file %s\n",stream_line_file);
			fclose(in);
			return ;
		}
	}
	else
	{
		sscanf(buf,"%d%d%d",&width,&height,&depth);
	}


	do{
		memset(buf,0,200);
		char c = ' ';
		do{
			if(feof(in))
				break;
			c = fgetc(in);
			if(c != ' ' && c != '\t' && c != '\n')
				break;
		}while(true);

		if(feof(in))
			break;
		
		buf[0] = c;
		fgets(buf+1,18,in);

		int tmp_num = atoi(buf);
		if(tmp_num <= 1)
		{
			flag = false;
			break;
		}
		if(!has_vel)
		{
			float* tmp_data = new float[3*tmp_num];
			for(int i = 0;i < tmp_num;i++)
			{
				fscanf(in,"%f%f%f",&tmp_data[i*3+0],&tmp_data[i*3+1],&tmp_data[i*3+2]);
			}
			datas.push_back(tmp_data);
			nums.push_back(tmp_num);
		}
		else
		{
			float* tmp_data = new float[3*tmp_num];
			float* tmp_vel = new float[3*tmp_num];
			for(int i = 0;i < tmp_num;i++)
			{
				fscanf(in,"%f%f%f%f%f%f",&tmp_data[i*3+0],&tmp_data[i*3+1],&tmp_data[i*3+2],&tmp_vel[i*3+0],&tmp_vel[i*3+1],&tmp_vel[i*3+2]);
			}
			datas.push_back(tmp_data);
			vels.push_back(tmp_data);
			nums.push_back(tmp_num);
		}
	}while(true);

	if(!flag)
	{
		printf("file %s may be invalid\n",stream_line_file);
		for(int i = 0;i < datas.size();i++)
			delete [](datas[i]);
	}
	else
	{

		int* idx = new int[datas.size()];
		float* values = new float[datas.size()];
		for(int i = 0;i < datas.size();i++)
		{
			idx[i] = i;
			values[i] = StreamLineInfoValue(nums[i],datas[i]);
		}

		SortStreamLineInfoValue(datas.size(),idx,values,false);


		_clearLines();

		int extend_N = __max(width,__max(height,depth));
		double shift_width = (extend_N-width)/(double)extend_N*0.5;
		double shift_height = (extend_N-height)/(double)extend_N*0.5;
		double shift_depth = (extend_N-depth)/(double)extend_N*0.5;

		for(int ll = 0;ll < datas.size();ll++)
		{
			int line_num = nums[idx[ll]]-1;
			GLfloat* vertex_list = new GLfloat[line_num*6];
			GLfloat* color_list = new GLfloat[line_num*6];
			GLuint* index_list = new GLuint[line_num*2];

			float r = rand()%1001/1000.0;
			float g = rand()%1001/1000.0;
			float b = 1.0 - r - g;


			for(int cc = 0; cc < line_num;cc++)
			{
				vertex_list[cc*6+0] = datas[idx[ll]][cc*3+0]/extend_N+shift_width-0.5;
				vertex_list[cc*6+1] = datas[idx[ll]][cc*3+1]/extend_N+shift_height-0.5;
				vertex_list[cc*6+2] = datas[idx[ll]][cc*3+2]/extend_N+shift_depth-0.5;

				vertex_list[cc*6+3] = datas[idx[ll]][cc*3+3]/extend_N+shift_width-0.5;
				vertex_list[cc*6+4] = datas[idx[ll]][cc*3+4]/extend_N+shift_height-0.5;
				vertex_list[cc*6+5] = datas[idx[ll]][cc*3+5]/extend_N+shift_depth-0.5;

				if(has_vel)
					Encode_color_from_vel(r,g,b,vels[idx[ll]][cc*3+0]/width,vels[idx[ll]][cc*3+1]/height,vels[idx[ll]][cc*3+2]/depth);
				color_list[cc*6+0] = r;
				color_list[cc*6+1] = g;
				color_list[cc*6+2] = b;

				if(has_vel)
					Encode_color_from_vel(r,g,b,vels[idx[ll]][cc*3+3]/width,vels[idx[ll]][cc*3+4]/height,vels[idx[ll]][cc*3+5]/depth);
				color_list[cc*6+3] = r;
				color_list[cc*6+4] = g;
				color_list[cc*6+5] = b;

				index_list[cc*2+0] = cc*2+0;
				index_list[cc*2+1] = cc*2+1;
			}

			vertex_lists.push_back(vertex_list);
			color_lists.push_back(color_list);
			index_lists.push_back(index_list);
			line_nums.push_back(line_num);

		}

		render_line_num = datas.size();
		char buf[20];
		sprintf(buf,"%d",datas.size());

		TwSetParam(bar,"RenderLineNum","max",TW_PARAM_CSTRING,1,buf);

		delete []idx;
		delete []values;
	}
}


void StreamLineGUI::DrawBox()
{
	float pts[8][3] = 
	{
		{-0.5,-0.5,-0.5},
		{-0.5,-0.5, 0.5},
		{-0.5, 0.5,-0.5},
		{-0.5, 0.5, 0.5},
		{ 0.5,-0.5,-0.5},
		{ 0.5,-0.5, 0.5},
		{ 0.5, 0.5,-0.5},
		{ 0.5, 0.5, 0.5}
	};
	glLineWidth(box_line_width);
	glColor3f(box_line_color[0],box_line_color[1],box_line_color[2]);
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

void StreamLineGUI::DrawStreamLines()
{
	glLineWidth(line_width);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	for(int i = 0;i < render_line_num;i++)
	{
		glVertexPointer(3,GL_FLOAT,0,vertex_lists[i]);
		glColorPointer(3,GL_FLOAT,0,color_lists[i]);
		glDrawElements(GL_LINES,line_nums[i]*2,GL_UNSIGNED_INT,index_lists[i]);
	}
}

void StreamLineGUI::SetQuaternionFromAxisAngle(const float *axis, float angle, float *quat)
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

void StreamLineGUI::ConvertQuaternionToMatrix(const float *quat, float *mat)
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

void StreamLineGUI::MultiplyQuaternions(const float *q1, const float *q2, float *qout)
{
    float qr[4];
	qr[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
	qr[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
	qr[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
	qr[3]  = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
    qout[0] = qr[0]; qout[1] = qr[1]; qout[2] = qr[2]; qout[3] = qr[3];
}


void Encode_color_from_vel(float& r, float& g, float& b, float x, float y, float z)
{
	float len = sqrt(x*x+y*y+z*z);
	r = fabs(x)*len;
	g = fabs(y)*len;
	b = fabs(z)*len;
	//r = g = b = len;

}

float StreamLineInfoValue(const int vertex_num, const float* points)
{
	if(vertex_num <= 2 || points == 0)
		return 0;

	float val = 0;

	for(int i = 0;i < vertex_num-2;i++)
	{
		float x1 = points[i*3+3] - points[i*3+0];
		float y1 = points[i*3+4] - points[i*3+1];
		float z1 = points[i*3+5] - points[i*3+2];
		float x2 = points[i*3+6] - points[i*3+3];
		float y2 = points[i*3+7] - points[i*3+4];
		float z2 = points[i*3+8] - points[i*3+5];

		float len1 = sqrt(x1*x1+y1*y1+z1*z1);
		float len2 = sqrt(x2*x2+y2*y2+z2*z2);

		if(len1 == 0 || len2 == 0)
			val += 0;
		else
		{
			float c_x = y1*z2-z1*y2;
			float c_y = z1*x2-x1*z2;
			float c_z = x1*y2-y1*x2;

			float len = sqrt(c_x*c_x+c_y*c_y+c_z*c_z);

			val += fabs(len/len1/len2)*(len1+len2);
		}

	}

	return val;

}

void SortStreamLineInfoValue(const int N, int* idx, float* values, bool accend )
{
	if(accend)
	{
		for(int pass = 0;pass < N-1;pass++)
		{
			for(int i = 0;i < N-1-pass;i++)
			{
				if(values[i] > values[i+1])
				{
					int tmp_idx = idx[i];
					idx[i] = idx[i+1];
					idx[i+1] = tmp_idx;
					float tmp_value = values[i];
					values[i] = values[i+1];
					values[i+1] = tmp_value;
				}
			}
		}
	}

	else
	{
		for(int pass = 0; pass < N-1;pass++)
		{
			for(int i = 0;i < N-1-pass;i++)
			{
				if(values[i] < values[i+1])
				{
					int tmp_idx = idx[i];
					idx[i] = idx[i+1];
					idx[i+1] = tmp_idx;
					float tmp_value = values[i];
					values[i] = values[i+1];
					values[i+1] = tmp_value;
				}
			}
		}
	}
}