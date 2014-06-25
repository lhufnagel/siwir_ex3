#pragma once
#include <string>

typedef unsigned int uint;

typedef struct Input
{

	uint sizex;
	uint sizey;
	uint timesteps;
	float omega;
	std::string vtk_file;
	uint vtk_step;

	Input()
	{
		sizex=50;
		sizey=50;
		timesteps=10000;
		omega=1.9;
		vtk_file="test";
		vtk_step=1000;
	}
} Input;

