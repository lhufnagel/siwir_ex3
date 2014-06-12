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
		sizex=100;
		sizey=100;
		timesteps=8000;
		omega=1.9;
		vtk_file="example";
		vtk_step=150;
	}
} Input;

/* todo: implementieren, und dann den default konstruktor rauswerfen  */
//bool readInputParams(const char* argv[1],Input* inp);

