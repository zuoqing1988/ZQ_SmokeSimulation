#include "ZQ_SmokeEditing3DGUI.h"


int main(int argc, char** argv)
{
	ZQ_SmokeEditing3D::ZQ_SmokeEditing3DGUI::GetInstance()->Init(&argc,argv);
	ZQ_SmokeEditing3D::ZQ_SmokeEditing3DGUI::GetInstance()->MainLoop();
	return 0;
}