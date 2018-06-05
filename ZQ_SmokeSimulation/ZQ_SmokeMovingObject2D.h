#ifndef _ZQ_SMOKE_MOVING_OBJECT_2D_H_
#define _ZQ_SMOKE_MOVING_OBJECT_2D_H_
#pragma once

#include <vector>
class ZQ_SmokeMovingObject2D
{
public:
	class Velocity
	{
	public:
		int x,y;
		Velocity(int u = 0, int v = 0):x(u),y(v){}
	};
public:
	ZQ_SmokeMovingObject2D();
	~ZQ_SmokeMovingObject2D();

public:
	int cx,cy;//rotation center
	int posx,posy; // position in world
	int sizex, sizey;
	int curstep;
	std::vector<Velocity> vels;
	bool* occupy;

public:
	bool LoadFromFile(const char* file);
	void GetCenter(int& cx, int& cy) const;
	void GetPos(int& posx, int& posy) const;
	void GetSize(int& sizex, int& sizey) const;
	void GetCurVel(int& u, int& v) const;
	const bool*& GetOccupyPtr() const;
	void UpdateOneStep();
};

bool LoadMovingObjects2D(const char* file, std::vector<ZQ_SmokeMovingObject2D*>& mvobjs);

#endif