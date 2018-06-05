#include "ZQ_SmokeGUI.h"


int main(int argc, char** argv)
{
	ZQ_SmokeGUI::GetInstance()->Init(&argc,argv);
	ZQ_SmokeGUI::GetInstance()->MainLoop();
	return 0;
}