#include "ZQ_SmokeMovingObject3D.h"

ZQ_SmokeMovingObject3D::ZQ_SmokeMovingObject3D()
{
	cx = cy = cz = 0;
	posx = posy = posz = 0;
	sizex = sizey = sizez = 0;
	curstep = 0;
	occupy = 0;
}

ZQ_SmokeMovingObject3D::~ZQ_SmokeMovingObject3D()
{
	if(occupy)
	{
		delete []occupy;
		occupy = 0;
	}

	vels.clear();
}

void ZQ_SmokeMovingObject3D::GetCenter(int &cx, int &cy, int &cz) const
{
	cx = this->cx;
	cy = this->cy;
	cz = this->cz;
}

void ZQ_SmokeMovingObject3D::GetSize(int &sizex, int &sizey, int &sizez) const
{
	sizex = this->sizex;
	sizey = this->sizey;
	sizez = this->sizez;
}

void ZQ_SmokeMovingObject3D::GetPos(int &posx, int &posy, int &posz) const
{
	posx = this->posx;
	posy = this->posy;
	posz = this->posz;
}

void ZQ_SmokeMovingObject3D::GetCurVel(int &u, int &v, int &w) const
{
	if(curstep >= vels.size() || curstep < 0)
		u = v = w = 0;
	else
	{
		u = vels[curstep].x;
		v = vels[curstep].y;
		w = vels[curstep].z;
	}
}

const bool*& ZQ_SmokeMovingObject3D::GetOccupyPtr() const
{
	return (const bool*&)occupy;
}

void ZQ_SmokeMovingObject3D::UpdateOneStep()
{
	if(curstep >= vels.size() || curstep < 0)
		curstep ++;
	else
	{
		posx += vels[curstep].x;
		posy += vels[curstep].y;
		posz += vels[curstep].z;
		curstep ++;
	}
}

bool LoadMovingObjects3D(const char* file, std::vector<ZQ_SmokeMovingObject3D*>& mvobjs)
{
	FILE * in = fopen(file,"r");
	if(in == 0)
		return false;

	int num_of_mvobj = 0;
	fscanf(in,"%d",&num_of_mvobj);
	for(int i = 0;i < num_of_mvobj;i++)
	{
		int sizex,sizey,sizez;
		int cx,cy,cz;
		int posx,posy,posz;
		int step;

		fscanf(in,"%d%d%d",&sizex,&sizey,&sizez);
		fscanf(in,"%d%d%d",&cx,&cy,&cz);
		fscanf(in,"%d%d%d",&posx,&posy,&posz);
		fscanf(in,"%d",&step);
		ZQ_SmokeMovingObject3D* obj = new ZQ_SmokeMovingObject3D();
		obj->sizex = sizex;
		obj->sizey = sizey;
		obj->sizez = sizez;
		obj->cx = cx;
		obj->cy = cy;
		obj->cz = cz;
		obj->posx = posx;
		obj->posy = posy;
		obj->posz = posz;

		for(int st = 0;st < step;st++)
		{
			int u,v,w;
			fscanf(in,"%d%d%d",&u,&v,&w);

			obj->vels.push_back(ZQ_SmokeMovingObject3D::Velocity(u,v,w));
		}

		obj->occupy = new bool[sizex*sizey*sizez];
		for(int bb = 0; bb < sizex*sizey*sizez;bb++)
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