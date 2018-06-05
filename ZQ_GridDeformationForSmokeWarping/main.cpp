#include "ZQ_GridDeformation.h"
#include "ZQ_DoubleImage.h"
#include "ZQ_ImageIO.h"
#include "InsertObjectToGrid.h"
#include "WarpGrid.h"
#include "WarpGrid3D.h"
#include "RenderGridWithTexture.h"
#include "RenderGridWithTexture3D.h"
#include "SampleWithMap.h"
#include "SampleWithMap3D.h"
#include "ResampleMapRBF.h"
#include "ResampleMapRBF3D.h"
#include "Tessellation.h"

using namespace ZQ;

bool main_insert_object(const int argc, const char** argv);
bool main_sample_with_map(const int argc, const char** argv);
bool main_sample_with_map3D(const int argc, const char** argv);
bool main_resample_map_rbf(const int argc, const char** argv);
bool main_resample_map_rbf3D(const int argc, const char** argv);
bool main_warp_grid(const int argc, const char** argv);
bool main_warp_grid3D(const int argc, const char** argv);
bool main_render_with_texture(const int argc, const char** argv);
bool main_render_with_texture3D(const int argc, const char** argv);
bool main_tessellation(const int argc, const char** argv);

void main(int argc, const char** argv)
{
	if(argc < 2)
	{
		printf(".exe methodtype [options]\n");
		return;
	}

	if(_strcmpi(argv[1],"insert_object") == 0)
	{
		main_insert_object(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"sample_with_map") == 0)
	{
		main_sample_with_map(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"sample_with_map3D") == 0)
	{
		main_sample_with_map3D(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"warp_grid") == 0)
	{
		main_warp_grid(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"warp_grid3D") == 0)
	{
		main_warp_grid3D(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"render_with_texture") == 0)
	{
		main_render_with_texture(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"render_with_texture3D") == 0)
	{
		main_render_with_texture3D(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"resample_map_rbf") == 0)
	{
		main_resample_map_rbf(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"resample_map_rbf3D") == 0)
	{
		main_resample_map_rbf3D(argc-2,argv+2);
	}
	else if(_strcmpi(argv[1],"tessellation") == 0)
	{
		main_tessellation(argc-2,argv+2);
	}

}

bool main_insert_object(const int argc, const char** argv)
{
	InsertObjectToGridOptions opt;
	if(!opt.HandleArgs(argc,argv))
		return false;

	if(!opt.has_input_file)
	{
		printf("need input file\n");
		return false;
	}
	const char* mask_file = opt.input_file;
	ZQ_DImage<float> ori_img;
	if(!ZQ_ImageIO::loadImage(ori_img,mask_file,0))
	{
		printf("failed to load %s\n",mask_file);
		return false;
	}
	int width = ori_img.width();
	int height = ori_img.height();
	ZQ_DImage<bool> mask(width,height);
	bool*& objMask = mask.data();
	for(int i = 0;i < width*height;i++)
	{
		objMask[i] = ori_img.data()[i] < 0.5;
	}

	InsertObjectToGrid mApp(width,height);
	mApp.SetObjectMask(objMask);
	mApp.SetMaxIteration(opt.max_iteration);
	mApp.SetFPIteration(opt.FPiteration);
	mApp.SetLineWeight(opt.line_weight);
	mApp.SetAngleWeight(opt.angle_weight);
	mApp.SetDistanceWeight(opt.distance_weight);

	if(!mApp.ShowDialog())
	{
		printf("fail\n");
		return false;
	}
	const float* out_coord = mApp.GetOutCoord();

	if(opt.has_show_file)
	{
		int show_scale = 20;
		int show_shift_x = show_scale/2;
		int show_shift_y = show_scale/2;
		IplImage* show_img = cvCreateImage(cvSize(width*show_scale,height*show_scale),IPL_DEPTH_8U,3);
		cvZero(show_img);

		ZQ_DImage<float> map_img(width,height,2);
		memcpy(map_img.data(),out_coord,sizeof(float)*width*height*2);
		bool retflag = RenderGridWithTexture::render_wireframe(show_img,map_img,mask,show_scale,show_shift_x,show_shift_y,false);
		cvSaveImage(opt.show_file,show_img);
		cvReleaseImage(&show_img);
		if(!retflag)
		{
			printf("failed to render\n");
			return false;
		}
	}

	if(opt.has_map_file)
	{
		ZQ_DImage<float> mapimg(width,height,2);
		float* map_data = mapimg.data();
		memcpy(map_data,out_coord,sizeof(float)*width*height*2);
		if(!mapimg.saveImage(opt.map_file))
		{
			return false;
		}
	}
	
	return true;
}

bool main_sample_with_map(const int argc, const char** argv)
{
	SampleWithMapOptions opt;
	if(!opt.HandleArgs(argc,argv))
	{
		return false;
	}
	if(!SampleWithMap::go(opt))
		return false;
	return true;
}


bool main_sample_with_map3D(const int argc, const char** argv)
{
	SampleWithMap3DOptions opt;
	if(!opt.HandleArgs(argc,argv))
	{
		return false;
	}
	if(!SampleWithMap3D::go(opt))
		return false;
	return true;
}


bool main_warp_grid(const int argc, const char** argv)
{
	WarpGridOptions opt;
	if(!opt.HandleArgs(argc,argv))
	{
		return false;
	}
	if(!WarpGrid::go(opt))
		return false;
	return true;
}

bool main_warp_grid3D(const int argc, const char** argv)
{
	WarpGrid3DOptions opt;
	if(!opt.HandleArgs(argc,argv))
	{
		return false;
	}
	if(!WarpGrid3D::go(opt))
		return false;
	return true;
}

bool main_render_with_texture(const int argc, const char** argv)
{
	RenderGridWithTextureOptions opt;
	if(!opt.HandleArgs(argc,argv))
		return false;

	if(!RenderGridWithTexture::go(opt))
	{
		printf("failed to run render_with_texture");
		return false;
	}
	return true;
}


bool main_render_with_texture3D(const int argc, const char** argv)
{
	/*ZQ_DImage3D<float> tex_img(2,2,2,1);
	float*& tex_data = tex_img.data();
	tex_data[0] = tex_data[1] = tex_data[2] = tex_data[3] = tex_data[4] = tex_data[5] = tex_data[6] = tex_data[7] = 1;
	ZQ_DImage3D<float> map_img(2,2,2,3);
	float coords[24] = {
		-50,-50,-50, 
		-100,-100,100,
		-100,100,-100,
		-100,100,100,
		100,-100,-100,
		100,-100,100,
		100,100,-100,
		50,50,50
	};
	memcpy(map_img.data(),coords,sizeof(float)*24);
	tex_img.saveImage("tex3d.di3");
	map_img.saveImage("map3d.di3");*/

	RenderGridWithTexture3DOptions opt;
	if(!opt.HandleArgs(argc,argv))
		return false;

	if(!RenderGridWithTexture3D::go(opt))
	{
		printf("failed to run render_with_texture3D");
		return false;
	}
	return true;
}

bool main_resample_map_rbf(const int argc, const char** argv)
{
	ResampleMapRBFOptions opt;
	if(!opt.HandleArgs(argc,argv))
		return false;

	if(!ResampleMapRBF::go(opt))
	{
		printf("failed to run resample map\n");
		return false;
	}
	return true;
}


bool main_resample_map_rbf3D(const int argc, const char** argv)
{
	ResampleMapRBF3DOptions opt;
	if(!opt.HandleArgs(argc,argv))
		return false;

	if(!ResampleMapRBF3D::go(opt))
	{
		printf("failed to run resample map 3D\n");
		return false;
	}
	return true;
}

bool main_tessellation(const int argc, const char** argv)
{
	TessellationOptions opt;
	if(!opt.HandleArgs(argc,argv))
		return false;
	if(!Tessellation::Go(opt))
	{
		printf("failed to run tessellation\n");
		return false;
	}
	return true;
}