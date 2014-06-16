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

/* todo: implementieren, und dann den default konstruktor rauswerfen  */
//bool readInputParams(const char* argv[1],Input* inp);

