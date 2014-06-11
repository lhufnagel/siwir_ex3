#include "main.h"
#include "grid.h"
#include "flags.h"
#include "directions.h"
#include "inputParams.h"

using namespace lbm;
int main(int argc, char* argv[])
{
	if (argc != 1/*2*/) {
		cerr<<"error: wrong number of arguments"<<endl;
		cout<<"call ./lbm params.dat "<<endl;
		exit(EXIT_FAILURE);
	}

	double uw[]={0.08,0.}; //x and y

	Input inp;
	/*
		if (!readInputParams(argv[1],&inp))
		{
		cerr << "Failed to open input-file" << endl;
		exit(EXIT_FAILURE);
		}
		*/

	PDF_Field src(inp.sizex+2,inp.sizey+2),
				 dst(inp.sizex+2,inp.sizey+2);
	V_Field velField(inp.sizex+2,inp.sizey+2);
	D_Field densField(inp.sizex+2,inp.sizey+2); //macroscopic
	Flags flagField(inp.sizex+2,inp.sizey+2); //info about wether boundary or fluid cell


	//initialize fields
	//todo probably put this in function
	for (uint x=0; x < inp.sizex+2; x++)
	{
		for (uint y=0; y < inp.sizey+2; y++)
		{
			flagField(x,y)= FLUID;

			velField(x,y,X)=velField(x,y,Y)=0.;

			densField(x,y)=1.;

			for (int i=0;i<stencilSize;i++)
				src(x,y,i) = dst(x,y,i) = weigths[i];
		}
	}

	//initialize boundary field
	for (uint y=1; y < inp.sizey+1; y++)
	{
		flagField(0,y)=flagField(inp.sizex+1,y)=NOSLIP;
	}
	for (uint x=0; x < inp.sizex+2; x++)
	{
		flagField(x,0)=NOSLIP;
		flagField(x,inp.sizey+1)=MOVNOSLIP;
	}


	for (uint tStep = 0; tStep<inp.timesteps; ++tStep)
	{
		//stream it
		for (uint x=0; x < inp.sizex+2; x++)
		{
			for (uint y=0; y < inp.sizey+2; y++)
			{
				if (flagField(x,y) == FLUID)
				{
					for (int i=0;i<stencilSize;i++)
						dst(x + stenNbs[X][i],y + stenNbs[Y][i],i) = src(x,y,i);

				}
			}
		}


	//handle boundary
	//don't use uints here ;)!
	for (int x=0; x < inp.sizex+2; x++)
	{

		for (int y=0; y < inp.sizey+2; y++)
		{
			if (flagField(x,y) == NOSLIP)
			{
				for (uint i=1;i<stencilSize;i++)
				{
					
					if (
						(x + stenNbs[X][i]) > -1 &&
						(x + stenNbs[X][i]) < (inp.sizex+2) &&
						(y + stenNbs[Y][i]) > -1 &&
						(y + stenNbs[Y][i]) < (inp.sizey+2) &&
						(flagField(x + stenNbs[X][i], y + stenNbs[Y][i]) == FLUID))
					{
						dst(x + stenNbs[X][i], y + stenNbs[Y][i],i) = dst(x,y,(i+3)%8+1);
					}
				}

			}
			else if (flagField(x,y) == MOVNOSLIP)
			{
				for (uint i=1;i<stencilSize;i++)
				{
					if (	(x + stenNbs[X][i])> -1 &&
						(x + stenNbs[X][i])< (inp.sizex+2) &&
						(y + stenNbs[Y][i])> -1 &&
						(y + stenNbs[Y][i])< (inp.sizey+2) &&
					(flagField(x + stenNbs[X][i], y + stenNbs[Y][i]) == FLUID))
					{
						                                        									//plus sign to compensate for the scalarproduct
						dst(x + stenNbs[X][i], y + stenNbs[Y][i] ,i) = dst(x,y,(i+3)%8+1) + 6.*weigths[i]*(stenNbs[X][i]*uw[X]+stenNbs[Y][i]*uw[Y]);
					}
				}

			}

		}
	}


		//collide it
		for (uint x=1; x < inp.sizex+1; x++)
		{
			for (uint y=1; y < inp.sizey+1; y++)
			{
				//calculate macroscopic measures
				densField(x,y)=dst(x,y,C);
				for (uint i=1;i<stencilSize;i++)
					densField(x,y)+=dst(x,y,i);

				velField(x,y,X)=
					(dst(x,y,NE)-
					 dst(x,y,W ))+
					(dst(x,y,E )-
					 dst(x,y,SW))+
					(dst(x,y,SE)-
					 dst(x,y,NW));


				velField(x,y,Y)=
					(dst(x,y,N )-
					 dst(x,y,SE))+
					(dst(x,y,NW)-
					 dst(x,y,S ))+
					(dst(x,y,NE)-
					 dst(x,y,SW));

				//seems like we should not divide by the density here...
				//otherwise the algorithm gets instable
				/*
				velField(x,y,X)/=densField(x,y);
				velField(x,y,Y)/=densField(x,y);
				*/

				double rho=densField(x,y);
				double vx=velField(x,y,X);
				double vy=velField(x,y,Y);
				double vsq = 3./2.*(vx*vx+vy*vy);


				for (int i=0;i<stencilSize;i++)
				{
					double scalProd = stenNbs[X][i]*vx+stenNbs[Y][i]*vy;
					dst(x,y,i )-= inp.omega*(dst(x,y,i) - weigths[i]*(rho + 3.*(scalProd) + 9./2*scalProd*scalProd  -vsq)); 
				}
			}
		}

		dst.swap(src);

		if (inp.vtk_step!=0 &&/*  tStep!=0 &&*/ (tStep % inp.vtk_step)==0)
		{
			//todo: put in function

			ofstream file; 
			string s;
			ostringstream outStream;
			outStream << tStep;
			s = outStream.str();
			file.open((inp.vtk_file + s +".vtk").c_str(), ios::out);
			if(!(file.is_open()))
			{ 
				printf("konnte nicht gespeichert werden\n");
				exit(1);
			}
			file << "# vtk DataFile Version 4.0" << endl;
			file << "SiWiRVisFile" << endl;
			file << "ASCII" << endl;
			file << "DATASET STRUCTURED_POINTS" << endl;
			file << "DIMENSIONS " << inp.sizex << " " << inp.sizey << " 1" << endl;
			file << "ORIGIN 0 0 0" << endl;
			file << "SPACING 1 1 1" << endl;
			file << "POINT_DATA "<< inp.sizex*inp.sizey << endl << endl;


			file << "SCALARS flags unsigned_int 1" << endl;
			file << "LOOKUP_TABLE default" << endl;
			for (uint y=1; y < inp.sizey+1; y++)
			{
				for (uint x=1; x < inp.sizex+1; x++)
				{
					file << flagField(x,y) << endl;
				}
			}
			file << endl;


			file << "SCALARS density double 1" << endl;
			file << "LOOKUP_TABLE default" << endl;
			for (uint y=1; y < inp.sizey+1; y++)
			{
				for (uint x=1; x < inp.sizex+1; x++)
				{
					file << densField(x,y) << endl;
				}
			}
			file << endl;

			file << "VECTORS velocity double" << endl;
			for (uint y=1; y < inp.sizey+1; y++)
			{
				for (uint x=1; x < inp.sizex+1; x++)
				{
					file << velField(x,y,X) << " " << velField(x,y,Y) << " 0" << endl;
				}
			}
			file << endl;

			file.close();
		}

	}

	return EXIT_SUCCESS;
}

