#include "ZQ_SmokeEditingGUI.h"

#pragma comment(lib,"glew32.lib")

int main(int argc, char** argv)
{
	ZQ_SmokeEditing::ZQ_SmokeEditingGUI::GetInstance()->Init(&argc,argv);
	ZQ_SmokeEditing::ZQ_SmokeEditingGUI::GetInstance()->MainLoop();
	return 0;
}