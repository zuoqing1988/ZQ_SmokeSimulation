#include "ZQ_SmokeMovingObject2D.h"

ZQ_SmokeMovingObject2D::ZQ_SmokeMovingObject2D()
{
	cx = cy = 0;
	posx = posy = 0;
	sizex = sizey = 0;
	curstep = 0;
	occupy = 0;
}

ZQ_SmokeMovingObject2D::~ZQ_SmokeMovingObject2D()
{
	if(occupy)
	{
		delete []occupy;
		occupy = 0;
	}

	vels.clear();
}

void ZQ_SmokeMovingObject2D::GetCenter(int &cx, int &cy) const
{
	cx = this->cx;
	cy = this->cy;
}

void ZQ_SmokeMovingObject2D::GetSize(int &sizex, int &sizey) const
{
	sizex = this->sizex;
	sizey = this->sizey;
}

void ZQ_SmokeMovingObject2D::GetPos(int &posx, int &posy) const
{
	posx = this->posx;
	posy = this->posy;
}

void ZQ_SmokeMovingObject2D::GetCurVel(int &u, int &v) const
{
	if(curstep >= vels.size() || curstep < 0)
		u = v = 0;
	else
	{
		u = vels[curstep].x;
		v = vels[curstep].y;
	}
}

const bool*& ZQ_SmokeMovingObject2D::GetOccupyPtr() const
{
	return (const bool*&)occupy;
}

void ZQ_SmokeMovingObject2D::UpdateOneStep()
{
	if(curstep >= vels.size() || curstep < 0)
		curstep ++;
	else
	{
		posx += vels[curstep].x;
		posy += vels[curstep].y;
		curstep ++;
	}
}

bool LoadMovingObjects2D(const char* file, std::vector<ZQ_SmokeMovingObject2D*>& mvobjs)
{
	FILE * in = fopen(file,"r");
	if(in == 0)
		return false;

	int num_of_mvobj = 0;
	fscanf(in,"%d",&num_of_mvobj);
	for(int i = 0;i < num_of_mvobj;i++)
	{
		int sizex,sizey;
		int cx,cy;
		int posx,posy;
		int step;

		fscanf(in,"%d%d",&sizex,&sizey);
		fscanf(in,"%d%d",&cx,&cy);
		fscanf(in,"%d%d",&posx,&posy);
		fscanf(in,"%d",&step);
		ZQ_SmokeMovingObject2D* obj = new ZQ_SmokeMovingObject2D();
		obj->sizex = sizex;
		obj->sizey = sizey;
		obj->cx = cx;
		obj->cy = cy;
		obj->posx = posx;
		obj->posy = posy;

		for(int st = 0;st < step;st++)
		{
			int u,v;
			fscanf(in,"%d%d",&u,&v);

			obj->vels.push_back(ZQ_SmokeMovingObject2D::Velocity(u,v));
		}

		obj->occupy = new bool[sizex*sizey];
		for(int bb = 0; bb < sizex*sizey;bb++)
		{
			int bval = 0;
			fscanf(in,"%d",&bval);
			obj->occupy[bb] = bool(bval);
		}

		mvobjs.push_back(obj);
	}
	fclose(in);
	return true;
}