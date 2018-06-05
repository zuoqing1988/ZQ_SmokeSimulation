#include "ZQ_SmokeSimulation2D.h"



void main(int argc, char** argv)
{
	if(argc != 3)
		return;

	char* config = argv[1];
	char* fold = argv[2];

	ZQ_SmokeSimulation2D* simu = new ZQ_SmokeSimulation2D();
	if(!simu->Init(config))
	{
		printf("init fail\n");
		return ;
	}
	simu->Run(fold);
	return ;

}