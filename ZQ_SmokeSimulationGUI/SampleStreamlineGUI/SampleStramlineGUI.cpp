#include "StreamLineGUI.h"

int main(int argc, char** argv)
{
	StreamLineGUI::GetInstance()->Init(&argc,argv);
	StreamLineGUI::GetInstance()->MainLoop();
	return 0;
}