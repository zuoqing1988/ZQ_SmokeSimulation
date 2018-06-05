#ifndef _ZQ_SMOKE_MOVING_OBJECT_3D_H_
#define _ZQ_SMOKE_MOVING_OBJECT_3D_H_
#pragma once

#include <vector>
class ZQ_SmokeMovingObject3D
{
public:
	class Velocity
	{
	public:
		int x,y,z;
		Velocity(int u = 0, int v = 0, int w = 0):x(u),y(v),z(w){}
	};
public:
	ZQ_SmokeMovingObject3D();
	~ZQ_SmokeMovingObject3D();

public:
	int cx,cy,cz;//rotation center
	int posx,posy,posz; // position in world
	int sizex, sizey, sizez;
	int curstep;
	std::vector<Velocity> vels;
	bool* occupy;

public:
	bool LoadFromFile(const char* file);
	void GetCenter(int& cx, int& cy, int& cz) const;
	void GetPos(int& posx, int& posy, int& posz) const;
	void GetSize(int& sizex, int& sizey, int& sizez) const;
	void GetCurVel(int& u, int& v, int& w) const;
	const bool*& GetOccupyPtr() const;
	void UpdateOneStep();
};

bool LoadMovingObjects3D(const char* file, std::vector<ZQ_SmokeMovingObject3D*>& mvobjs);

#endif