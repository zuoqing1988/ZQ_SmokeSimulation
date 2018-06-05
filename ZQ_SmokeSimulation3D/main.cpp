#include "ZQ_SmokeSimulation3D.h"


void main(int argc, char** argv)
{
	if(argc != 3)
		return;

	char* config = argv[1];
	char* fold = argv[2];

	ZQ_SmokeSimulation3D* simu = new ZQ_SmokeSimulation3D();
	if(!simu->Init(config))
	{
		printf("init fail\n");
		return ;
	}
	simu->Run(fold);

	return ;

}