#include <windows.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include "ZQ_SmokeEditingFigGrid.h"

#include <time.h>

using namespace ZQ_SmokeEditing;
using namespace ZQ;

BaseFigGrid::BaseFigGrid()
{
	tex_id = 0;
	has_tex = false;
	tex_file[0] = '\0';
	fold[0] = '\0';
	prefix[0] = '\0';
	strcpy(suffix,"jpg");
	total_frame = 200;
	base_id = 0;
	bar = 0;
	focused = false;
	move_step = 5;
	draw_texture = true;
	tex_alpha = 1;
}

BaseFigGrid::~BaseFigGrid()
{
	if(bar)
	{
		TwDeleteBar(bar);
		bar = 0;
	}
	if(has_tex)
	{
		glDeleteTextures(1,&tex_id);
		tex_id = 0;
		has_tex = false;
	}
}


void BaseFigGrid::SetFocused(bool b)
{
	if(b)
	{
		if(bar == 0)
			InitTweakBar();
	}
	else
	{
		if(bar)
			DeleteTweakBar();
	}
	focused = b;
}

void BaseFigGrid::InitTweakBar()
{
	_createTweakBar();

	_addGlobalTweakBar();
	
	_addLocalTweakBar();
}

void BaseFigGrid::_addGlobalTweakBar()
{
	if(bar != 0)
	{	
		TwAddVarRW(bar,"tex_file",TW_TYPE_CSSTRING(sizeof(tex_file)),&tex_file,"label='tex_file' group='Group1'");
		TwAddButton(bar,"LoadTexture",(TwButtonCallback)OnLoadTexture,this,"label='Load Tex' group='Group1'");
		TwAddVarRW(bar,"fold",TW_TYPE_CSSTRING(sizeof(fold)),&fold,"label='fold' group='Group1'");
		TwAddVarRW(bar,"prefix",TW_TYPE_CSSTRING(sizeof(prefix)),&prefix,"label='prefix' group='Group1'");
		TwAddVarRW(bar,"suffix",TW_TYPE_CSSTRING(sizeof(suffix)),&suffix,"label='suffix' group='Group1'");
		TwAddVarRW(bar,"total_frame",TW_TYPE_INT32,&total_frame,"label='total frame' group='Group1'");
		TwAddVarRW(bar,"base_id",TW_TYPE_INT32,&base_id,"label='base_id' group='Group1'");

		TwAddVarRW(bar,"move_step",TW_TYPE_FLOAT,&move_step,"label='move step' min=0.1 max=10 step=0.1");
		TwAddVarRW(bar,"rot_step",TW_TYPE_FLOAT,&rot_step,"label='rot_step' min=0.1 max=10 step=0.1");
		TwAddVarRW(bar,"draw_tex",TW_TYPE_BOOL8,&draw_texture,"label='draw tex'");
	}
}


void BaseFigGrid::OnLoadTexture(void* clientData)
{
	BaseFigGrid* g = (BaseFigGrid*)clientData;
	g->_loadTexture(true);
}


void BaseFigGrid::_loadTexture(bool _forceUpdate)
{
	ZQ::ZQ_DImage<float> img;
	if(!ZQ_ImageIO::loadImage(img,tex_file,1))
	{
		printf("failed to load %s\n",tex_file);
		return;
	}
	img.FlipY();

	int im_width = img.width();
	int im_height = img.height();
	int im_nChannels = img.nchannels();
	float*& img_ptr = img.data();

	int nChannels = 3;
	float* data = new float[im_width*im_height*nChannels];
	if(im_nChannels == 1)
	{	
		for(int i = 0;i < im_width*im_height;i++)
		{
			data[i*nChannels+0] = img_ptr[i];
			data[i*nChannels+1] = img_ptr[i];
			data[i*nChannels+2] = img_ptr[i];
			//data[i*nChannels+3] = 0;
		}
	}
	else if(im_nChannels == 3)
	{
		for(int i = 0;i < im_width*im_height;i++)
		{
			data[i*nChannels+0] = img_ptr[i*3+2];
			data[i*nChannels+1] = img_ptr[i*3+1];
			data[i*nChannels+2] = img_ptr[i*3+0];
			//data[i*nChannels+3] = 0;
		}
	}
	else
	{
		delete []data;
		printf("image must be 1 or 3 channels\n");
		return;
	}

	if(has_tex)
	{
		glDeleteTextures(1,&tex_id);
	}
	glPixelStorei(GL_UNPACK_ALIGNMENT,1);

	glGenTextures(1,&tex_id);
	glBindTexture(GL_TEXTURE_2D,tex_id);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,im_width,im_height,0,GL_RGB,GL_FLOAT,data);
	glBindTexture(GL_TEXTURE_2D,0);
	delete []data;
	has_tex = true;
}

MaskFigGrid::MaskFigGrid(int w, int h)
{
	color_for_true[0] = 0.6;
	color_for_true[1] = 0.6;
	color_for_true[2] = 0.6;
	color_for_false[0] = 0;
	color_for_false[1] = 0;
	color_for_false[2] = 0.6;
	alpha_for_true = 0.2;
	alpha_for_false = 0.2;
	color_for_eraser_true[0] = 0.8;
	color_for_eraser_true[1] = 0.8;
	color_for_eraser_true[2] = 0.8;
	color_for_eraser_false[0] = 0;
	color_for_eraser_false[1] = 0;
	color_for_eraser_false[2] = 0.8;

	eraser_value = true;
	eraser_size = 11;
	eraser_select_x = 0;
	eraser_select_y = 0;

	edit_mode = false;

	width = w;
	height = h;
	mask.allocate(width,height);
	_loadTexture(true);
}

MaskFigGrid::~MaskFigGrid()
{

}


void MaskFigGrid::Draw(float z_shift)
{
	int mask_width = mask.width();
	int mask_height = mask.height();
	float offset_x = -mask_width*0.5;
	float offset_y = -mask_height*0.5;

	float z_pos = z_shift;
	float z_pos_for_eraser = 0.001+z_shift;

	float pts[4][3] = 
	{
		{ 0.5*mask_width, 0.5*mask_height,z_pos},
		{-0.5*mask_width, 0.5*mask_height,z_pos},
		{-0.5*mask_width,-0.5*mask_height,z_pos},
		{ 0.5*mask_width,-0.5*mask_height,z_pos}
	};
	float texcoord[4][2] = 
	{
		1,1,
		0,1,
		0,0,
		1,0
	};
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D,tex_id);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glColor4f(1,1,1,tex_alpha);
	glBegin(GL_QUADS);
	glTexCoord2f(texcoord[0][0],texcoord[0][1]);
	glVertex3f(pts[0][0],pts[0][1],pts[0][2]);
	glTexCoord2f(texcoord[1][0],texcoord[1][1]);
	glVertex3f(pts[1][0],pts[1][1],pts[1][2]);
	glTexCoord2f(texcoord[2][0],texcoord[2][1]);
	glVertex3f(pts[2][0],pts[2][1],pts[2][2]);
	glTexCoord2f(texcoord[3][0],texcoord[3][1]);
	glVertex3f(pts[3][0],pts[3][1],pts[3][2]);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	if(edit_mode)
	{
		int mask_width = mask.width();
		int mask_height = mask.height();
		float offset_x = -mask_width*0.5;
		float offset_y = -mask_height*0.5;
		float half_eraser_size = 0.5*eraser_size;
		float pts[4][3] = 
		{
			{ offset_x+eraser_select_x+0.5+half_eraser_size,   offset_y+eraser_select_y+0.5+half_eraser_size, z_pos_for_eraser},
			{ offset_x+eraser_select_x+0.5-half_eraser_size,   offset_y+eraser_select_y+0.5+half_eraser_size, z_pos_for_eraser},
			{ offset_x+eraser_select_x+0.5-half_eraser_size,   offset_y+eraser_select_y+0.5-half_eraser_size, z_pos_for_eraser},
			{ offset_x+eraser_select_x+0.5+half_eraser_size,   offset_y+eraser_select_y+0.5-half_eraser_size, z_pos_for_eraser}
		};

		glLineWidth(1);
		if(eraser_value)
			glColor4f(color_for_eraser_true[0],color_for_eraser_true[1],color_for_eraser_true[2],tex_alpha);
		else
			glColor4f(color_for_eraser_false[0],color_for_eraser_false[1],color_for_eraser_false[2],tex_alpha);
		glBegin(GL_LINES);
		glVertex3f(pts[0][0],pts[0][1],pts[0][2]);
		glVertex3f(pts[1][0],pts[1][1],pts[1][2]);
		glVertex3f(pts[1][0],pts[1][1],pts[1][2]);
		glVertex3f(pts[2][0],pts[2][1],pts[2][2]);
		glVertex3f(pts[2][0],pts[2][1],pts[2][2]);
		glVertex3f(pts[3][0],pts[3][1],pts[3][2]);
		glVertex3f(pts[3][0],pts[3][1],pts[3][2]);
		glVertex3f(pts[0][0],pts[0][1],pts[0][2]);
		glEnd();
	}
}

void MaskFigGrid::_createTweakBar()
{
	if(bar == 0)
	{
		const char* name = "MaskGrid";
		TwBar* tmp_bar = TwGetBarByName(name);
		if(tmp_bar != 0)
			TwDeleteBar(tmp_bar);


		bar = TwNewBar(name);
		TwSetTopBar(bar);
		char buf[200];
		sprintf(buf,"%s refresh=0.1 movable=false resizable=false color='0 0 0' alpha='100' position='0 0' size='200 180'",name);
		TwDefine(buf);
	}
}

void MaskFigGrid::_addGlobalTweakBar()
{

}

void MaskFigGrid::_addLocalTweakBar()
{
	if(bar != 0)
	{
		TwAddVarRW(bar,"tex_alpha",TW_TYPE_FLOAT,&tex_alpha,"label='tex alpha' min=0 max=1 step= 0.1");
		TwAddVarRW(bar,"width",TW_TYPE_INT32,&width,"label='width'");
		TwAddVarRW(bar,"height",TW_TYPE_INT32,&height,"label='height'");
		TwAddButton(bar,"Resize",(TwButtonCallback)OnResize,this,"label='Resize'");
		TwAddVarRW(bar,"eraser_size",TW_TYPE_INT32,&eraser_size,"label='eraser_size' min=1 max=23 step=2");
		TwAddVarRW(bar,"eraser_value",TW_TYPE_BOOL8,&eraser_value,"label='eraser_value'");
		TwAddVarRW(bar,"edit_mode",TW_TYPE_BOOL8,&edit_mode,"label='edit mode'");
	}
}


void MaskFigGrid::OnMouseButton(int button_id, int state, int mx, int my)
{
	if(button_id == GLUT_LEFT_BUTTON)
	{
		if(state == GLUT_DOWN)
		{
			if(edit_mode)
			{
				_select_eraser_pos(mx,my);
				_update_mask_with_eraser();
			}
		}
	}
}

void MaskFigGrid::OnMouseMove(int mx, int my)
{
	if(edit_mode)
	{
		_select_eraser_pos(mx,my);
		_update_mask_with_eraser();
	}
}

void MaskFigGrid::OnMousePassiveMove(int mx, int my)
{
	if(edit_mode)
	{
		_select_eraser_pos(mx,my);
	}
}

bool MaskFigGrid::Export(const char* export_fold, const char* script_file, int width, int height)
{
	char buf[2000];
	sprintf(buf,"%s\\%s",export_fold,script_file);

	ZQ_DImage<bool> tmp_mask(mask);
	tmp_mask.FlipY();
	return ZQ_ImageIO::saveImage(tmp_mask,buf);
}

void MaskFigGrid::OnResize(void* clientData)
{
	MaskFigGrid* g = (MaskFigGrid*)clientData;
	g->_resize();
}



void MaskFigGrid::_resize()
{
	if(!mask.matchDimension(width,height,1))
	{
		ZQ_DImage<bool> old_mask = mask;
		mask.allocate(width,height);

		int old_width = old_mask.width();
		int old_height = old_mask.height();
		bool*& old_mask_data = old_mask.data();
		bool*& mask_data = mask.data();
		int shift_x = (width - old_width)/2;
		int shift_y = (height - old_height)/2;
		for(int i = 0;i < height;i++)
		{
			for(int j = 0;j < width;j++)
			{
				if(i < shift_y || i >= shift_y+old_height || j < shift_x || j >= shift_x+old_width)
				{
					mask_data[i*width+j] = false;
				}
				else
				{
					mask_data[i*width+j] = old_mask_data[(i-shift_y)*old_width+j-shift_x];
				}
			}
		}
		_loadTexture(true);
	}
}

void MaskFigGrid::_loadTexture(bool _forceUpdate)
{

	int mask_width = mask.width();
	int mask_height = mask.height();
	if(mask_width == 0 || mask_height == 0)
		return;

	if(!has_tex || _forceUpdate)
	{
		if(has_tex)
		{
			glDeleteTextures(1,&tex_id);
			has_tex = false;
			tex_id = 0;
		}

		glGenTextures(1,&tex_id);
	}

	int nChannels = 3;
	bool*& mask_data = mask.data();
	float* data = new float[mask_width*mask_height*nChannels];
	for(int i = 0;i < mask_width*mask_height;i++)
	{
		if(mask_data[i])
		{
			data[i*nChannels+0] = color_for_true[0];
			data[i*nChannels+1] = color_for_true[1];
			data[i*nChannels+2] = color_for_true[2];
			//data[i*nChannels+3] = tex_alpha;
		}
		else
		{
			data[i*nChannels+0] = color_for_false[0];
			data[i*nChannels+1] = color_for_false[1];
			data[i*nChannels+2] = color_for_false[2];
			//data[i*nChannels+3] = tex_alpha;
		}
	}
	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
	glBindTexture(GL_TEXTURE_2D,tex_id);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,mask_width,mask_height,0,GL_RGB,GL_FLOAT,data);
	glBindTexture(GL_TEXTURE_2D,0);
	delete []data;
	has_tex = true;

}

void MaskFigGrid::_select_eraser_pos(float real_mx, float real_my)
{
	int mask_width = mask.width();
	int mask_height = mask.height();

	real_mx = real_mx+mask_width*0.5;
	real_my = real_my+mask_height*0.5;

	eraser_select_x = real_mx+0.5;
	eraser_select_y = real_my+0.5;
}

void MaskFigGrid::_update_mask_with_eraser()
{
	int mask_width = mask.width();
	int mask_height = mask.height();
	ZQ_DImage<bool> old_mask = mask;
	bool*& old_mask_data = old_mask.data();
	bool*& mask_data = mask.data();
	int half_eraser_size = eraser_size/2;
	for(int i = eraser_select_y-half_eraser_size;i <= eraser_select_y+half_eraser_size;i++)
	{
		for(int j = eraser_select_x-half_eraser_size;j <= eraser_select_x+half_eraser_size;j++)
		{
			if(i >= 0 && i < mask_height && j >= 0 && j < mask_width)
				mask_data[i*mask_width+j] = eraser_value;
		}
	}

	bool updata_flag = memcmp(mask_data,old_mask_data,mask_width*mask_height*sizeof(bool)) != 0;
	if(updata_flag)
	{
		_loadTexture(false);
	}
}

RigidFigGrid::RigidFigGrid(float w, float h)
{
	width = w;
	height = h;
	rot_angle = 0;
	scale_x = 1;
	scale_y = 1;
	tran_x = 0;
	tran_y = 0;
	rot_step = 1;
	scale_step = 0.1;
	scaleXenable = true;
	scaleYenable = true;

	line_width_for_border = 1;
	line_width_for_boarder_focused = 2;
	color_for_border[0] = 0.5; color_for_border[1] = 0.5; color_for_border[2] = 0.5;
	color_for_border_focused[0] = 1; color_for_border_focused[1] = 0; color_for_border_focused[2] = 0; 
	z_pos_for_line = -0.99;
	z_pos_for_tex = -1.0;

	draw_border = true;
}


void RigidFigGrid::Draw(float z_shift)
{
	double m_pi = 3.1415926535;
	float rot_rad = rot_angle/180.0*m_pi;
	float total_mat[9];
	ZQ::ZQ_MergeImage::MakeMatrix(ZQ::ZQ_Vec2D(scale_x,scale_y),rot_rad,ZQ::ZQ_Vec2D(tran_x,tran_y),total_mat);

	float ori_corners[12] = 
	{
		-width*0.5,   width*0.5,   width*0.5, -width*0.5,
		-height*0.5, -height*0.5,  height*0.5, height*0.5,
		1,           1,           1,          1
	};

	float final_corners[12];
	ZQ_MathBase::MatrixMul(total_mat,ori_corners,3,3,4,final_corners);

	float corners[8];
	corners[0] = final_corners[0];
	corners[1] = final_corners[4];
	corners[2] = final_corners[1];
	corners[3] = final_corners[5];
	corners[4] = final_corners[2];
	corners[5] = final_corners[6];
	corners[6] = final_corners[3];
	corners[7] = final_corners[7];

	float tex_coords[8] = 
	{
		0,0,1,0,1,1,0,1
	};

	/////////

	float z_tex = z_pos_for_tex+z_shift;
	float z_line = z_pos_for_line+z_shift;
	

	if(draw_texture && has_tex)
	{
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D,tex_id);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glColor4f(1,1,1,tex_alpha);
		glBegin(GL_QUADS);
		glTexCoord2f(tex_coords[0],tex_coords[1]);
		glVertex3f(corners[0],corners[1],z_tex);
		glTexCoord2f(tex_coords[2],tex_coords[3]);
		glVertex3f(corners[2],corners[3],z_tex);
		glTexCoord2f(tex_coords[4],tex_coords[5]);
		glVertex3f(corners[4],corners[5],z_tex);
		glTexCoord2f(tex_coords[6],tex_coords[7]);
		glVertex3f(corners[6],corners[7],z_tex);
		glEnd();
		glBindTexture(GL_TEXTURE_2D,0);
		glDisable(GL_TEXTURE_2D);
	}

	if(draw_border)
	{
		if(!focused)
		{
			glLineWidth(line_width_for_border);
			glColor3f(color_for_border[0],color_for_border[1],color_for_border[2]);
		}
		else
		{
			glLineWidth(line_width_for_boarder_focused);
			glColor3f(color_for_border_focused[0],color_for_border_focused[1],color_for_border_focused[2]);
		}

		glBegin(GL_LINES);
		for(int i = 0;i < 4;i++)
		{
			glVertex3f(corners[i*2+0],corners[i*2+1],z_line);
			glVertex3f(corners[(i+1)%4*2+0],corners[(i+1)%4*2+1],z_line);
		}
		glEnd();
	}
}

void RigidFigGrid::_createTweakBar()
{
	if(bar == 0)
	{
		const char* name = "RigidGrid";
		TwBar* tmp_bar = TwGetBarByName(name);
		if(tmp_bar != 0)
			TwDeleteBar(tmp_bar);


		bar = TwNewBar(name);
		TwSetTopBar(bar);
		char buf[200];
		sprintf(buf,"%s refresh=0.1 movable=false resizable=false color='0 0 0' alpha='100' position='0 180' size='200 600'",name);
		TwDefine(buf);
	}
}

void RigidFigGrid::_addGlobalTweakBar()
{
	if(bar != 0)
	{
		TwAddVarRW(bar,"tex_alpha",TW_TYPE_FLOAT,&tex_alpha,"label='tex alpha' min=0 max=1 step= 0.1 group='Group1'");
		TwAddVarRW(bar,"tex_file",TW_TYPE_CSSTRING(sizeof(tex_file)),&tex_file,"label='tex_file' group='Group1'");
		TwAddButton(bar,"LoadTexture",(TwButtonCallback)OnLoadTexture,this,"label='Load Tex' group='Group1'");
		TwAddVarRW(bar,"fold",TW_TYPE_CSSTRING(sizeof(fold)),&fold,"label='fold' group='Group1'");
		TwAddVarRW(bar,"prefix",TW_TYPE_CSSTRING(sizeof(prefix)),&prefix,"label='prefix' group='Group1'");
		TwAddVarRW(bar,"suffix",TW_TYPE_CSSTRING(sizeof(suffix)),&suffix,"label='suffix' group='Group1'");
		TwAddVarRW(bar,"total_frame",TW_TYPE_INT32,&total_frame,"label='total frame' group='Group1'");
		TwAddVarRW(bar,"base_id",TW_TYPE_INT32,&base_id,"label='base_id' group='Group1'");
	}
}

void RigidFigGrid::_addLocalTweakBar()
{
	if(bar != 0)
	{
		TwAddVarRW(bar,"width",TW_TYPE_FLOAT,&width,"label='width' readonly=false group='Group2'");
		TwAddVarRW(bar,"height",TW_TYPE_FLOAT,&height,"label='height' readonly=false group='Group2'");
		TwAddVarRW(bar,"trans_x",TW_TYPE_FLOAT,&tran_x,"label='trans_x' readonly=true group='Group2'");
		TwAddVarRW(bar,"trans_y",TW_TYPE_FLOAT,&tran_y,"label='trans_y' readonly=true group='Group2'");
		TwAddVarRW(bar,"scale_x",TW_TYPE_FLOAT,&scale_x,"label='scale_x' readonly=true group='Group2'");
		TwAddVarRW(bar,"scale_y",TW_TYPE_FLOAT,&scale_y,"label='scale_y' readonly=true group='Group2'");
		TwAddVarRW(bar,"rot_angle",TW_TYPE_FLOAT,&rot_angle,"label='rot' readonly=true group='Group2'");

		TwAddVarRW(bar,"move_step",TW_TYPE_FLOAT,&move_step,"label='move step' min=0.1 max=10 step=0.1 group='Group3'");
		TwAddVarRW(bar,"rot_step",TW_TYPE_FLOAT,&rot_step,"label='rot_step' min=0.1 max=10 step=0.1 group='Group3'");
		TwAddVarRW(bar,"draw_tex",TW_TYPE_BOOL8,&draw_texture,"label='draw tex' group='Group3'");
		TwAddVarRW(bar,"scale_step",TW_TYPE_FLOAT,&scale_step,"label='scale_step' min=0.01 max=0.5 step=0.01 group='Group3'");
		TwAddVarRW(bar,"scale_x_enable",TW_TYPE_BOOL8,&scaleXenable,"label='scale_x' group='Group3'");
		TwAddVarRW(bar,"scale_y_enable",TW_TYPE_BOOL8,&scaleYenable,"label='scale_y' group='Group3'");
		TwAddVarRW(bar,"draw_border",TW_TYPE_BOOL8,&draw_border,"label='draw border' group='Group3'");
	}
}

void RigidFigGrid::OnKeyboard(int key_id, int mx, int my)
{
	switch(key_id)
	{
	case 'w':case 'W':
		{
			tran_y += move_step;
		}
		break;
	case 's':case 'S':
		{
			tran_y -= move_step;
		}
		break;
	case 'a':case 'A':
		{
			tran_x -= move_step;
		}
		break;
	case 'd':case 'D':
		{
			tran_x += move_step;
		}
		break;
	case 'q':case 'Q':
		{
			rot_angle += rot_step;
		}
		break;
	case 'e':case 'E':
		{
			rot_angle -= rot_step;
		}
		break;
	case 'z':case 'Z':
		{
			if(scaleXenable)
				scale_x *= (1+scale_step);
			if(scaleYenable)
				scale_y *= (1+scale_step);
		}
		break;
	case 'x':case 'X':
		{
			if(scaleXenable)
				scale_x *= (1-scale_step);
			if(scaleYenable)
				scale_y *= (1-scale_step);
		}
		break;
	}
}

bool RigidFigGrid::Export(const char* export_fold, const char* script_file, int w, int h)
{
	FILE* out = fopen(script_file,"a");
	if(out == 0)
	{
		printf("failed to open file %s\n",script_file);
		return false;
	}

	char cmdbuf[2000];
	std::string cmdstr;
	for(int i = 0;i < total_frame;i++)
	{
		cmdstr.clear();
		sprintf(cmdbuf,"ZQ_MergeImage MethodType MERGE_DENSITY AxisUp MergeModeBlend OutputFile %s\\%d.png ",export_fold,i);
		cmdstr.append(cmdbuf);
		sprintf(cmdbuf,"BackgroundSize %d %d ",w,h);
		cmdstr.append(cmdbuf);
		sprintf(cmdbuf,"MergeSource %s\\%s%d.%s ",fold,prefix,i,suffix);
		cmdstr.append(cmdbuf);
		sprintf(cmdbuf,"%.1f %.1f %.5f %.5f %.5f ",tran_x,tran_y,rot_angle,width*scale_x,height*scale_y);
		cmdstr.append(cmdbuf);
		fprintf(out,"%s\n",cmdstr.c_str());
	}
	fclose(out);
	return true;
}

DeformFigGrid::DeformFigGrid()
{
	z_pos_for_tex = -1.0;
	z_pos_for_grid = -0.99;
	z_pos_for_fixed_pt = -0.98;
	z_pos_for_key_pt = -0.97;
	z_pos_for_ctrl_pt = -0.96;

	const float COLOR_FOR_GRID[3] = {0,0.8,0};
	const float COLOR_FOR_FIXED_POINT[3] = {0.8,0,0.8};
	const float COLOR_FOR_KEY_POINT[3] = {0.8,0,0.8};
	const float COLOR_FOR_CONTROL_POINT[3] = {1,0,0};
	memcpy(color_for_grid,COLOR_FOR_GRID,sizeof(float)*3);
	memcpy(color_for_fixed_pt,COLOR_FOR_FIXED_POINT,sizeof(float)*3);
	memcpy(color_for_key_pt,COLOR_FOR_KEY_POINT,sizeof(float)*3);
	memcpy(color_for_ctrl_pt,COLOR_FOR_CONTROL_POINT,sizeof(float)*3);

	line_width_for_grid = 1;
	point_size_for_fixed_pt = 5;
	point_size_for_key_pt = 10;
	point_size_for_ctrl_pt = 10;
	draw_grid = true;
	border_size = 0;

	edit_mode = false;

	width = 0;
	height = 0;
	line_weight = 0;
	angle_weight = 1;
	distance_weight = 0;
	distance = 1;
	arap = false;
	cur_arap = false;
	neighbor_type = ZQ_GridDeformationOptions::NEIGHBOR_4;

	total_iteration = 10000;
	cur_iteration = 0;
	per_frame_iteration = 60;
	with_good_init = false;
	grid_scale = 10;

	move_step = 0.5;
	rot_step = 1;
	xloop = false;
	ctrl_pts_num = 0;
	ctrl_id = -1;
}

DeformFigGrid::~DeformFigGrid()
{

}

void DeformFigGrid::_move_all_pts(float x, float y)
{
	int npixles = coords.npixels();
	float*& coord_ptr = coords.data();
	for(int i = 0;i < npixles;i++)
	{
		coord_ptr[i*2+0] += x;
		coord_ptr[i*2+1] += y;
	}
}

void DeformFigGrid::_reset_deform_iteration()
{
	cur_iteration = total_iteration;
	with_good_init = false;
}


void DeformFigGrid::Draw(float z_shift)
{
	if(draw_grid)
	{
		_draw_grid(z_shift);
	}
	if(draw_texture && has_tex)
	{
		_draw_tex(z_shift);
	}
	if(focused && draw_grid)
	{
		_draw_control_pts(z_shift);
	}

	_draw_user_defined(z_shift);
}

void DeformFigGrid::UpdatePerFrame()
{
	if(cur_iteration <= 0)
		return;

	clock_t t1 = clock();
	switch(opt.methodType)
	{
	case ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_DISTANCE_ENERGY:
	case ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_DISTANCE_ENERGY_XLOOP:
		{
			ZQ_DImage<float> init_coord(coords);
			if(!with_good_init)
			{
				int real_it = -1;
				deform._deformation_without_scaling(init_coord.data(),coords.data(),per_frame_iteration,real_it);
				
				with_good_init = true;
				cur_iteration -= per_frame_iteration;
			}

			deform._deformation_with_distance(init_coord.data(),coords.data(),1,per_frame_iteration);
			cur_iteration -= per_frame_iteration;
			
		}
		break;
	case ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_ENERGY:
	case ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_ENERGY_XLOOP:
		{
			ZQ_DImage<float> init_coord(coords);
			int real_it = -1;
			deform._deformation_without_scaling(init_coord.data(),coords.data(),per_frame_iteration,real_it);
			cur_iteration -= per_frame_iteration;
		}
		break;
	case ZQ_GridDeformationOptions::METHOD_ARAP_VERT_AS_CENTER:
		{
			ZQ_DImage<float> init_coord(coords);
			if(!with_good_init)
			{
				int real_it = -1;
				deform._deformation_without_scaling(init_coord.data(),coords.data(),per_frame_iteration,real_it);

				with_good_init = true;
				cur_iteration -= per_frame_iteration;
			}

			deform._deformation_ARAP_VERT(init_coord.data(),coords.data(),1,per_frame_iteration);
			cur_iteration -= per_frame_iteration;
		}
		break;
	case ZQ_GridDeformationOptions::METHOD_ARAP_VERT_AS_CENTER_XLOOP:
		{
			ZQ_DImage<float> init_coord(coords);
			if(!with_good_init)
			{
				int real_it = -1;
				deform._deformation_without_scaling(init_coord.data(),coords.data(),per_frame_iteration,real_it);

				with_good_init = true;
				cur_iteration -= per_frame_iteration;
			}

			deform._deformation_ARAP_VERT_XLOOP(init_coord.data(),coords.data(),1,per_frame_iteration);
			cur_iteration -= per_frame_iteration;
		}
		break;
	}
	clock_t t2 = clock();
	printf("solve: %.3f ms\n",(double)(t2-t1));

}

void DeformFigGrid::OnMouseButton(int button_id, int state, int mx, int my)
{
	switch(button_id)
	{
	case GLUT_LEFT_BUTTON:
		{
			if(edit_mode)
			{
				if(state == GLUT_DOWN)
				{
					_select_ctrl_pt(mx,my);
				}
			}
		}
		break;
	}
}

void DeformFigGrid::OnKeyboard(int key_id, int mx, int my)
{
	switch(key_id)
	{
	case 'w':case 'W':
		{
			if(edit_mode)
				_move_ctrl_pt(0,move_step);
			else
				_move_all_pts(0,move_step);
		}
		break;
	case 's':case 'S':
		{
			if(edit_mode)
				_move_ctrl_pt(0,-move_step);
			else
				_move_all_pts(0,-move_step);
		}
		break;
	case 'a':case 'A':
		{
			if(edit_mode)
				_move_ctrl_pt(-move_step,0);
			else
				_move_all_pts(-move_step,0);
		}
		break;
	case 'd':case 'D':
		{
			if(edit_mode)
				_move_ctrl_pt(move_step,0);
			else
				_move_all_pts(move_step,0);
		}
		break;
	case 'q':case 'Q':
		{
			_rotate(rot_step);
		}
		break;
	case 'e':case 'E':
		{
			_rotate(-rot_step);
		}
		break;
	}
}

bool DeformFigGrid::Export(const char* export_fold, const char* script_file, int width, int height)
{
	ZQ_DImage<float> tmp_coords(coords);
	float*& coord_ptr = tmp_coords.data();
	int npixels = tmp_coords.npixels();
	tmp_coords.FlipY();
	for(int i = 0;i < npixels;i++)
	{
		coord_ptr[i*2+1] = -coord_ptr[i*2+1];
	}
	tmp_coords.Multiplywith(grid_scale);

	float boxmin[2] = {-width/2.0,-height/2.0};
	float boxmax[2] = {width/2.0,height/2.0};

	char buf[2000];
	char mapfilename[FILENAME_MAX];
	char objfilename[FILENAME_MAX];

	sprintf(mapfilename,"%s\\map.di2",export_fold);
	if(!tmp_coords.saveImage(mapfilename))
	{
		printf("failed to save file %s\n",mapfilename);
		return false;
	}
	
	bool has_obj = false;
	for(int i = 0;i < nouseful_flag.npixels();i++)
	{
		if(nouseful_flag.data()[i])
		{
			has_obj = true;
			break;
		}
	}
	if(has_obj)
	{
		ZQ_DImage<bool> tmp_flag = nouseful_flag;
		tmp_flag.FlipY();
		sprintf(objfilename,"%s\\obj_mask.png",export_fold);
		if(!ZQ_ImageIO::saveImage(tmp_flag,objfilename))
		{
			printf("failed to save file %s\n",objfilename);
			return false;
		}
	}

	FILE* out = fopen(script_file,"a");
	if(out == 0)
	{
		printf("failed to open file %s\n",script_file);
		return false;
	}

	std::string cmdstr;
	for(int fr = 0;fr < total_frame;fr++)
	{
		cmdstr.clear();
		cmdstr.append("ZQ_GridDeformationForSmokeWarping render_with_texture map_file ");
		cmdstr.append(mapfilename);
		if(has_obj)
		{
			cmdstr.append(" map_mask_file ");
			cmdstr.append(objfilename);
		}
		cmdstr.append(" tex_file ");
		sprintf(buf,"%s\\%s%d.%s ",fold,prefix,fr,suffix);
		cmdstr.append(buf);
		cmdstr.append(" show_scale 1 show_file ");
		sprintf(buf,"%s\\%d.png ",export_fold,fr);
		cmdstr.append(buf);
		sprintf(buf,"boundingboxmin %.1f %.1f ",boxmin[0],boxmin[1]);
		cmdstr.append(buf);
		sprintf(buf,"boundingboxmax %.1f %.1f ",boxmax[0],boxmax[1]);
		cmdstr.append(buf);
		if(xloop)
			cmdstr.append(" xloop ");
		sprintf(buf," border_size %d ",border_size);
		cmdstr.append(buf);
		fprintf(out,"%s\n",cmdstr.c_str());
	}
	fclose(out);
	return true;
}

void DeformFigGrid::GlobalMove(float x, float y)
{
	int npixels = coords.npixels();
	float*& coord_ptr = coords.data();
	for(int i = 0;i < npixels;i++)
	{
		coord_ptr[i*2] += x/grid_scale;
		coord_ptr[i*2+1] += y/grid_scale;
	}
	
}

void DeformFigGrid::OnApply(void* clientData)
{
	DeformFigGrid* g = (DeformFigGrid*)clientData;
	g->_apply();
}

void DeformFigGrid::_apply()
{
	bool rebuild_flag = false;
	if(!coords.matchDimension(width,height,2))
	{
		_init_ctrl_pts();
		coords.allocate(width,height,2);
		nouseful_flag.allocate(width,height);
		fixed_flag.allocate(width,height);

		_init_coords();
		_init_fixed_pts();
		rebuild_flag = true;
	}

	if(opt.line_weight != line_weight)
	{
		opt.line_weight = line_weight;
		rebuild_flag = true;
	}
	if(opt.angle_weight != angle_weight)
	{
		opt.angle_weight = angle_weight;
		rebuild_flag = true;
	}
	if(opt.distance_weight != distance_weight)
	{
		opt.distance_weight = distance_weight;
		rebuild_flag = true;
	}
	if(opt.distance != distance)
	{
		opt.distance = distance;
		rebuild_flag = true;
	}
	if(cur_arap != arap)
	{
		cur_arap = arap;
		rebuild_flag = true;
	}
	if(opt.neighborType != neighbor_type)
	{
		opt.neighborType = neighbor_type;
		rebuild_flag = true;
	}

	if(rebuild_flag)
	{
		_rebuild_matrix();
		_reset_deform_iteration();
	}
}

void DeformFigGrid::_rebuild_matrix()
{
	if(arap)
	{
		opt.methodType = xloop ? ZQ_GridDeformationOptions::METHOD_ARAP_VERT_AS_CENTER_XLOOP :
			ZQ_GridDeformationOptions::METHOD_ARAP_VERT_AS_CENTER;
	}
	else
	{
		if(distance_weight == 0)
			opt.methodType = xloop ? ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_ENERGY_XLOOP :
			ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_ENERGY;
		else
			opt.methodType = xloop ? ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_DISTANCE_ENERGY_XLOOP :
			ZQ_GridDeformationOptions::METHOD_LINE_ANGLE_DISTANCE_ENERGY;
	}
	

	deform.BuildMatrix(width,height,nouseful_flag.data(),fixed_flag.data(),opt);
}

void DeformFigGrid::_draw_grid(float z_shift)
{
	float z_grid = z_pos_for_grid+z_shift;
	int H = coords.height();
	int W = coords.width();
	float*& coord_ptr = coords.data();
	bool*& nouseful_ptr = nouseful_flag.data();

	if(!xloop)
	{
		glLineWidth(line_width_for_grid);
		glColor3f(color_for_grid[0],color_for_grid[1],color_for_grid[2]);
		glBegin(GL_LINES);
		for(int i = 0;i < H;i++)
		{
			for(int j = 0;j < W-1;j++)
			{
				if(nouseful_ptr[i*W+j] || nouseful_ptr[i*W+j+1])
					continue;
				glVertex3f(coord_ptr[(i*W+j)*2+0]*grid_scale,coord_ptr[(i*W+j)*2+1]*grid_scale,z_grid);
				glVertex3f(coord_ptr[(i*W+j+1)*2+0]*grid_scale,coord_ptr[(i*W+j+1)*2+1]*grid_scale,z_grid);
			}
		}
		for(int i = 0;i < H-1;i++)
		{
			for(int j = 0;j < W;j++)
			{
				if(nouseful_ptr[i*W+j] || nouseful_ptr[(i+1)*W+j])
					continue;
				glVertex3f(coord_ptr[(i*W+j)*2+0]*grid_scale,coord_ptr[(i*W+j)*2+1]*grid_scale,z_grid);
				glVertex3f(coord_ptr[((i+1)*W+j)*2+0]*grid_scale,coord_ptr[((i+1)*W+j)*2+1]*grid_scale,z_grid);
			}
		}
		glEnd();
	}
	else
	{
		glLineWidth(line_width_for_grid);
		glColor3f(color_for_grid[0],color_for_grid[1],color_for_grid[2]);
		glBegin(GL_LINES);
		for(int i = 0;i < H;i++)
		{
			for(int j = 0;j < W;j++)
			{
				int j1 = (j+1)%W;
				if(nouseful_ptr[i*W+j] || nouseful_ptr[i*W+j1])
					continue;
				glVertex3f(coord_ptr[(i*W+j)*2+0]*grid_scale,coord_ptr[(i*W+j)*2+1]*grid_scale,z_grid);
				glVertex3f(coord_ptr[(i*W+j1)*2+0]*grid_scale,coord_ptr[(i*W+j1)*2+1]*grid_scale,z_grid);
			}
		}
		for(int i = 0;i < H-1;i++)
		{
			for(int j = 0;j < W;j++)
			{
				if(nouseful_ptr[i*W+j] || nouseful_ptr[(i+1)*W+j])
					continue;
				glVertex3f(coord_ptr[(i*W+j)*2+0]*grid_scale,coord_ptr[(i*W+j)*2+1]*grid_scale,z_grid);
				glVertex3f(coord_ptr[((i+1)*W+j)*2+0]*grid_scale,coord_ptr[((i+1)*W+j)*2+1]*grid_scale,z_grid);
			}
		}
		glEnd();
	}	
}

void DeformFigGrid::_draw_tex(float z_shift)
{
	float z_tex = z_pos_for_tex+z_shift;

	int H = coords.height();
	int W = coords.width();
	float*& coord_ptr = coords.data();
	bool*& nouseful_ptr = nouseful_flag.data();

	float tex_len_x = xloop ? W : (W-1.0-2*border_size);
	float tex_len_y = H-1.0-2*border_size;
	if(tex_len_x <= 0) 
		tex_len_x = 1;
	else
		tex_len_x = 1.0/tex_len_x;


	if(tex_len_y <= 0) 
		tex_len_y = 1;
	else 
		tex_len_y = 1.0/tex_len_y;

	int i_start = border_size;
	int i_end = H-1-border_size;
	int j_start = xloop ? 0 : border_size;
	int j_end = xloop ? W : W-1-border_size;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D,tex_id);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glColor4f(1,1,1,tex_alpha);
	glBegin(GL_QUADS);

	for(int i = border_size;i < H-1-border_size;i++)
	{
		for(int j = j_start;j < j_end;j++)
		{
			int j1 = (j+1)%W;
			if(nouseful_ptr[i*W+j] || nouseful_ptr[i*W+j1] || nouseful_ptr[(i+1)*W+j1] || nouseful_ptr[(i+1)*W+j])
				continue;
			
			glTexCoord2f((j-j_start)*tex_len_x,(i-i_start)*tex_len_y);
			glVertex3f(coord_ptr[(i*W+j)*2+0]*grid_scale,coord_ptr[(i*W+j)*2+1]*grid_scale,z_tex);

			glTexCoord2f((j+1-j_start)*tex_len_x,(i-i_start)*tex_len_y);
			glVertex3f(coord_ptr[(i*W+j1)*2+0]*grid_scale,coord_ptr[(i*W+j1)*2+1]*grid_scale,z_tex);

			glTexCoord2f((j+1-j_start)*tex_len_x,(i+1-i_start)*tex_len_y);
			glVertex3f(coord_ptr[((i+1)*W+j1)*2+0]*grid_scale,coord_ptr[((i+1)*W+j1)*2+1]*grid_scale,z_tex);

			glTexCoord2f((j-j_start)*tex_len_x,(i+1-i_start)*tex_len_y);
			glVertex3f(coord_ptr[((i+1)*W+j)*2+0]*grid_scale,coord_ptr[((i+1)*W+j)*2+1]*grid_scale,z_tex);
		}
	}
	glEnd();
	glBindTexture(GL_TEXTURE_2D,0);
	glDisable(GL_TEXTURE_2D);
}

void DeformFigGrid::_draw_control_pts(float z_shift)
{

	float z_fixed = z_pos_for_fixed_pt+z_shift;
	float z_key = z_pos_for_key_pt+z_shift;
	float z_ctrl = z_pos_for_ctrl_pt+z_shift;

	int H = coords.height();
	int W = coords.width();
	float*& coord_ptr = coords.data();
	bool*& nouseful_ptr = nouseful_flag.data();

	float scale = grid_scale;

	bool*& fixed_ptr = fixed_flag.data();
	glPointSize(point_size_for_fixed_pt);
	glColor3f(color_for_fixed_pt[0],color_for_fixed_pt[1],color_for_fixed_pt[2]);
	glBegin(GL_POINTS);
	for(int i = 0;i < W*H;i++)
	{
		if(!nouseful_ptr[i] && fixed_ptr[i])
			glVertex3f(coord_ptr[i*2+0]*scale,coord_ptr[i*2+1]*scale,z_fixed);
	}
	glEnd();

	if(edit_mode)
	{
		//corner points
		glPointSize(point_size_for_key_pt);
		glColor3f(color_for_key_pt[0],color_for_key_pt[1],color_for_key_pt[2]);
		glBegin(GL_POINTS);
		for(int ii = 0;ii < ctrl_pts_num;ii++)
			glVertex3f(coord_ptr[(ctrl_pts_y[ii]*W+ctrl_pts_x[ii])*2+0]*scale,coord_ptr[(ctrl_pts_y[ii]*W+ctrl_pts_x[ii])*2+1]*scale,z_key);
		glEnd();

		if(ctrl_id >= 0)
		{
			glPointSize(point_size_for_ctrl_pt);
			glColor3f(color_for_ctrl_pt[0],color_for_ctrl_pt[1],color_for_ctrl_pt[2]);
			glBegin(GL_POINTS);
			glVertex3f(coord_ptr[(ctrl_pts_y[ctrl_id]*W+ctrl_pts_x[ctrl_id])*2+0]*scale,coord_ptr[(ctrl_pts_y[ctrl_id]*W+ctrl_pts_x[ctrl_id])*2+1]*scale,z_ctrl);
			glEnd();
		}
	}
}

void DeformFigGrid::_select_ctrl_pt(float real_mx, float real_my)
{
	if(ctrl_pts_num <= 0)
	{
		ctrl_id = -1;
		return;
	}

	int W = coords.width();
	int H = coords.height();
	float*& coord_ptr = coords.data();

	float scale = grid_scale;
	int s_id = 0;
	float min_dis = sqrt((real_mx-coord_ptr[(ctrl_pts_y[0]*W+ctrl_pts_x[0])*2+0]*scale)*(real_mx-coord_ptr[(ctrl_pts_y[0]*W+ctrl_pts_x[0])*2+0]*scale)
		+(real_my-coord_ptr[(ctrl_pts_y[0]*W+ctrl_pts_x[0])*2+1]*scale)*(real_my-coord_ptr[(ctrl_pts_y[0]*W+ctrl_pts_x[0])*2+1]*scale));
	for(int kk = 1;kk < ctrl_pts_num;kk++)
	{
		float cur_dis = sqrt((real_mx-coord_ptr[(ctrl_pts_y[kk]*W+ctrl_pts_x[kk])*2+0]*scale)*(real_mx-coord_ptr[(ctrl_pts_y[kk]*W+ctrl_pts_x[kk])*2+0]*scale)
			+(real_my-coord_ptr[(ctrl_pts_y[kk]*W+ctrl_pts_x[kk])*2+1]*scale)*(real_my-coord_ptr[(ctrl_pts_y[kk]*W+ctrl_pts_x[kk])*2+1]*scale));
		if(cur_dis < min_dis)
		{
			min_dis = cur_dis;
			s_id = kk;
		}
	}

	ctrl_id = s_id;
}


void DeformFigGrid::_move_ctrl_pt(float tran_x, float tran_y)
{
	if(ctrl_id < 0)
		return;

	int W = coords.width();
	int H = coords.height();
	float*& coord_ptr = coords.data();

	coord_ptr[(ctrl_pts_y[ctrl_id]*W+ctrl_pts_x[ctrl_id])*2+0] += tran_x;
	coord_ptr[(ctrl_pts_y[ctrl_id]*W+ctrl_pts_x[ctrl_id])*2+1] += tran_y;
	_recompute_fixed();
	_reset_deform_iteration();
}

void DeformFigGrid::_rotate(float angle)
{
	int nPixels = coords.npixels();
	if(nPixels == 0)
		return;

	float*& coord_ptr = coords.data();

	double m_pi = 3.1415926535;
	float rot_rad = angle/180.0*2*m_pi;
	float total_mat[9];

	float cx = 0;
	float cy = 0;
	for(int i = 0;i < nPixels;i++)
	{
		cx += coord_ptr[i*2]*grid_scale;
		cy += coord_ptr[i*2+1]*grid_scale;
	}
	cx /= nPixels;
	cy /= nPixels;

	ZQ::ZQ_MergeImage::MakeMatrix(ZQ::ZQ_Vec2D(1,1),rot_rad,ZQ::ZQ_Vec2D(0,0),total_mat);

	for(int i = 0;i < nPixels;i++)
	{
		float ori_pt[3] = {coord_ptr[i*2]*grid_scale-cx,coord_ptr[i*2+1]*grid_scale-cy,1};
		float out_pt[3];
		ZQ_MathBase::MatrixMul(total_mat,ori_pt,3,3,1,out_pt);
		coord_ptr[i*2+0] = (out_pt[0]+cx)/grid_scale;
		coord_ptr[i*2+1] = (out_pt[1]+cy)/grid_scale;
	}
}


void DeformFigGrid::_addGlobalTweakBar()
{
	const int neighbor_type_num = 3;
	const TwEnumVal neighbor_typeEnum[neighbor_type_num] =
	{
		{ZQ_GridDeformationOptions::NEIGHBOR_4, "neighbor 4"},
		{ZQ_GridDeformationOptions::NEIGHBOR_8, "neighbor 8"},
		{ZQ_GridDeformationOptions::NEIGHBOR_12, "neighbor 12"}
	};

	if(bar != 0)
	{
		TwAddVarRW(bar,"tex_alpha",TW_TYPE_FLOAT,&tex_alpha,"label='tex alpha' min=0 max=1 step= 0.1 group='Group1'");
		TwAddVarRW(bar,"border",TW_TYPE_INT32,&border_size,"label='border' min=0 group='Group1'");
		TwAddVarRW(bar,"tex_file",TW_TYPE_CSSTRING(sizeof(tex_file)),&tex_file,"label='tex_file' group='Group1'");
		TwAddButton(bar,"LoadTexture",(TwButtonCallback)OnLoadTexture,this,"label='Load Tex' group='Group1'");
		TwAddVarRW(bar,"fold",TW_TYPE_CSSTRING(sizeof(fold)),&fold,"label='fold' group='Group1'");
		TwAddVarRW(bar,"prefix",TW_TYPE_CSSTRING(sizeof(prefix)),&prefix,"label='prefix' group='Group1'");
		TwAddVarRW(bar,"suffix",TW_TYPE_CSSTRING(sizeof(suffix)),&suffix,"label='suffix' group='Group1'");
		TwAddVarRW(bar,"total_frame",TW_TYPE_INT32,&total_frame,"label='total frame' group='Group1'");
		TwAddVarRW(bar,"base_id",TW_TYPE_INT32,&base_id,"label='base_id' group='Group1'");

		TwAddVarRW(bar,"width",TW_TYPE_INT32,&width,"label='width' min=4 max=100 group='Group2'");
		TwAddVarRW(bar,"height",TW_TYPE_INT32,&height,"label='height' min=4 max=100 group='Group2'");
		TwAddVarRW(bar,"line_weight",TW_TYPE_FLOAT,&line_weight,"label='line' min=0 max=100 group='Group2'");
		TwAddVarRW(bar,"angle_weight",TW_TYPE_FLOAT,&angle_weight,"label='angle' min=0.01 max=100 group='Group2'");
		TwAddVarRW(bar,"distance_weight",TW_TYPE_FLOAT,&distance_weight,"label='distance' min=0 max=100 group='Group2'");
		TwAddVarRW(bar,"distance_scale",TW_TYPE_FLOAT,&distance,"label='scale' min=0 max=100 group='Group2'");
		TwAddVarRW(bar,"ARAP",TW_TYPE_BOOL8,&arap,"label='ARAP' group='Group2'");
		TwType neighborTWType = TwDefineEnum("neighborType", neighbor_typeEnum, neighbor_type_num);
		TwAddVarRW(bar, "neighborType", neighborTWType, &neighbor_type, "label='neighbor type' group='Group2'");
		TwAddButton(bar,"Apply",(TwButtonCallback)OnApply,this,"label='Apply' group='Group2'");

		TwAddVarRW(bar,"edit_mode",TW_TYPE_BOOL8,&edit_mode,"label='edit' group='Group3'");
		TwAddVarRW(bar,"draw_grid",TW_TYPE_BOOL8,&draw_grid,"label='draw grid' group='Group3'");
		TwAddVarRW(bar,"draw_tex",TW_TYPE_BOOL8,&draw_texture,"label='draw tex' group='Group3'");
		TwAddVarRW(bar,"grid_scale",TW_TYPE_FLOAT,&grid_scale,"label='grid_scale' min=1 max=50 group='Group3'");
		TwAddVarRW(bar,"move_step",TW_TYPE_FLOAT,&move_step,"label='move step' min=0.1 max=10 step=0.1 group='Group3'");
		TwAddVarRW(bar,"rot_step",TW_TYPE_FLOAT,&rot_step,"label='rot_step' min=0.1 max=10 step=0.1 group='Group3'");
		TwAddVarRW(bar,"total_it",TW_TYPE_INT32,&total_iteration,"label='total_it' min=1000 max=100000 step=1000 group='Group3'");
		TwAddVarRW(bar,"perframe_it",TW_TYPE_INT32,&per_frame_iteration,"label='perframe_it' min=20 max=2000 step=20 group='Group3'");	
	}
}

