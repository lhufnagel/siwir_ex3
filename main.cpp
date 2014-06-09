#include "main.h"
#include "grid.h"
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

	double uw=0.08;

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
	//Flags flagField(inp.sizex+2,inp.sizey+2); //info about wether boundary or fluid cell

	//initializeFlags(flagField);

	for (uint x=0; x < inp.sizex+2; x++)
	{
		for (uint y=0; y < inp.sizey+2; y++)
		{
			velField(x,y,X)=0.;
			velField(x,y,Y)=0.;
			densField(x,y)=1.;

			src(x,y,C )= 4./9.; 

			src(x,y,N )= 1./9.; 
			src(x,y,S )= 1./9.; 
			src(x,y,E )= 1./9.; 
			src(x,y,W )= 1./9.; 

			src(x,y,NE)=  1./36.; 
			src(x,y,SW)=  1./36.; 
			src(x,y,SE)=  1./36.; 
			src(x,y,NW)=  1./36.; 

			dst(x,y,C )= 4./9.; 

			dst(x,y,N )= 1./9.; 
			dst(x,y,S )= 1./9.; 
			dst(x,y,E )= 1./9.; 
			dst(x,y,W )= 1./9.; 

			dst(x,y,NE)=  1./36.; 
			dst(x,y,SW)=  1./36.; 
			dst(x,y,SE)=  1./36.; 
			dst(x,y,NW)=  1./36.; 
		}
	}

/*
	for (uint x=0; x < inp.sizex+2; x++)
	{
		velField(x,inp.sizey+1,X)=uw;
	}
*/


	for (uint tStep = 0; tStep<inp.timesteps; ++tStep)
	{

		//stream it
		for (uint x=1; x < inp.sizex+1; x++)
		{
			for (uint y=1; y < inp.sizey+1; y++)
			{

/*
				cout << "C " <<  dst(x,y,C) << endl;
				cout << "N " <<  dst(x,y,N) << endl;
				cout << "NE "  << dst(x,y,NE) << endl;
				cout << "E " <<  dst(x,y,E) << endl;
				cout << "SE " <<  dst(x,y,SE) << endl;
				cout << "S " <<  dst(x,y,S) << endl;
				cout << "SW "  << dst(x,y,SW) << endl;
				cout << "W " <<  dst(x,y,W) << endl;
				cout << "NW " <<  dst(x,y,NW) << endl;
*/

				dst(x  ,y  ,C )=src(x,y,C );
				dst(x  ,y+1,N )=src(x,y,N );
				dst(x+1,y+1,NE)=src(x,y,NE);
				dst(x+1,y  ,E )=src(x,y,E );
				dst(x+1,y-1,SE)=src(x,y,SE);
				dst(x  ,y-1,S )=src(x,y,S );
				dst(x-1,y-1,SW)=src(x,y,SW);
				dst(x-1,y  ,W )=src(x,y,W );
				dst(x-1,y+1,NW)=src(x,y,NW);

			}
		}

		//handle boundary
		for (uint y=1; y < inp.sizey+1; y++)
		{
			// x=0
			dst(1,y  ,E )=dst(0,y,W );
			dst(1,y+1,NE)=dst(0,y,SW);
			dst(1,y-1,SE)=dst(0,y,NW);

			// x=inp.sizex+1
			dst(inp.sizex,y  ,W )=dst(inp.sizex+1,y,E );
			dst(inp.sizex,y+1,NW)=dst(inp.sizex+1,y,SE);
			dst(inp.sizex,y-1,SW)=dst(inp.sizex+1,y,NE);
		}

		for (uint x=1; x < inp.sizex+1; x++)
		{
			// y=0
			dst(x  , 1,N )=dst(x,0,S );
			dst(x+1,1,NE)=dst(x,0,SW);
			dst(x-1,1,NW)=dst(x,0,SE);

			// y=inp.sizey+1
			dst(x  ,inp.sizey,S )=dst(x,inp.sizey+1,N );//-6./9.*0;
			dst(x+1,inp.sizey,SE)=dst(x,inp.sizey+1,NW)+1./6.*uw;
			dst(x-1,inp.sizey,SW)=dst(x,inp.sizey+1,NE)-1./6.*uw;
		}

		dst(1,1,NE)=dst(0,0,SW);
		dst(inp.sizex,1,NW)=dst(inp.sizex+1,0,SE);
		dst(inp.sizex,inp.sizey,SW)=dst(inp.sizex+1,inp.sizey+1,NE)-1./6.*uw;
		dst(1,inp.sizey,SE)=        dst(0,inp.sizey+1,NW)+1./6.*uw;



		//collide it
		for (uint x=1; x < inp.sizex+1; x++)
		{
			for (uint y=1; y < inp.sizey+1; y++)
			{
				//calculate macroscopic measures
				densField(x,y)=dst(x,y,C)+
					dst(x,y,N )+
					dst(x,y,E )+
					dst(x,y,S )+
					dst(x,y,W )+
					dst(x,y,NE)+
					dst(x,y,SE)+
					dst(x,y,SW)+
					dst(x,y,NW);

				velField(x,y,X)=((dst(x,y,NE)+
							dst(x,y,E )+
							dst(x,y,SE))-
						(dst(x,y,SW)+
						 dst(x,y,W )+
						 dst(x,y,NW)));
					///densField(x,y);
/*

				velField(x,y,X)=(
						(dst(x,y,NE)-
						 dst(x,y,W ))+
						(dst(x,y,E )-
						dst(x,y,SW))+
							(dst(x,y,SE)-
						 dst(x,y,NW)));
					///densField(x,y);

*/
				velField(x,y,Y)=((dst(x,y,N )+
							dst(x,y,NE)+
							dst(x,y,NW))-
						(dst(x,y,SE)+
						 dst(x,y,S )+
						 dst(x,y,SW)));
					///densField(x,y);
/*

				velField(x,y,Y)=(
							(dst(x,y,N )-
							 dst(x,y,SW))+
							(dst(x,y,NW)-
							 dst(x,y,S ))+
							(dst(x,y,NE)-
							dst(x,y,SE)))
					///densField(x,y);

*/

				double rho=densField(x,y);
				double vx=velField(x,y,X);
				double vy=velField(x,y,Y);
				double vsq = 3./2.*(vx*vx+vy*vy);


				dst(x,y,C )-= inp.omega*(dst(x,y,C) - 4./9.*(rho-vsq)); 

				dst(x,y,N )-= inp.omega*(dst(x,y,N) - 1./9.*(rho+3.*vy+9./2.*vy*vy-vsq)); 
				dst(x,y,S )-= inp.omega*(dst(x,y,S) - 1./9.*(rho-3.*vy+9./2.*vy*vy-vsq)); 
				dst(x,y,E )-= inp.omega*(dst(x,y,E) - 1./9.*(rho+3.*vx+9./2.*vx*vx-vsq)); 
				dst(x,y,W )-= inp.omega*(dst(x,y,W) - 1./9.*(rho-3.*vx+9./2.*vx*vx-vsq)); 

				dst(x,y,NE)-= inp.omega*(dst(x,y,NE) - 1./36.*(rho+3.*(vx+vy)+9./2.*(vx+vy)*(vx+vy)-vsq)); 
				dst(x,y,SW)-= inp.omega*(dst(x,y,SW) - 1./36.*(rho-3.*(vx+vy)+9./2.*(vx+vy)*(vx+vy)-vsq)); 
				dst(x,y,SE)-= inp.omega*(dst(x,y,SE) - 1./36.*(rho+3.*(vx-vy)+9./2.*(vx-vy)*(vx-vy)-vsq)); 
				dst(x,y,NW)-= inp.omega*(dst(x,y,NW) - 1./36.*(rho+3.*(vy-vx)+9./2.*(vy-vx)*(vy-vx)-vsq)); 


			}
		}

		if (inp.vtk_step!=0 &&  tStep!=0 && (tStep % inp.vtk_step)==0)
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
			//file << "SCALARS flags unsigned_int 1" << endl;
			//file << "LOOKUP_TABLE default" << endl;
			//file << "1" << endl;
			//file << "..." << endl;
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

		dst.swap(src);
	}

	return EXIT_SUCCESS;
}

