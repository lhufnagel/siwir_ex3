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

	Input inp;
	/*
	 * ENSURE, THAT WE ADD +2 TO sizex AND sizey AFTER READING THEM FROM THE INPUT FILE  ?! (maybe?)
   if (!readInputParams(argv[1],&inp))
   {
      cerr << "Failed to open input-file" << endl;
      exit(EXIT_FAILURE);
   }
	*/

	PDF_Field src(inp.sizex,inp.sizey),dst(inp.sizex,inp.sizey);
	V_Field velField(inp.sizex,inp.sizey);
	D_Field densField(inp.sizex,inp.sizey); //macroscopic!!
	Flags flagField(inp.sizex,inp.sizey); //info about boundary or fluid cell

	//initializeFlags(flagField);
	//initializeVel(velField);

	for (uint tStep = 0; tStep<inp.timesteps; ++tStep)
	{
		
		//stream it
		for (uint x=0; x < inp.sizex; x++)
		{
			for (uint y=0; y < inp.sizey; y++)
			{
				if (flagField(x,y) != 1/*FLUID*/) continue;

				dst(x  ,y  ,C )=src(x,y,C );
				dst(x+1,y  ,N )=src(x,y,N );
				dst(x+1,y+1,NE)=src(x,y,NE);
				dst(x  ,y+1,E )=src(x,y,E );
				dst(x-1,y+1,SE)=src(x,y,SE);
				dst(x-1,y  ,S )=src(x,y,S );
				dst(x-1,y-1,SW)=src(x,y,SW);
				dst(x  ,y-1,W )=src(x,y,W );
				dst(x+1,y-1,NW)=src(x,y,NW);
		
				//calculate macroscopic measures
				densField(x,y)=src(x,y,C)+
				               src(x,y,N )+
				               src(x,y,NE)+
				               src(x,y,E )+
				               src(x,y,SE)+
				               src(x,y,S )+
				               src(x,y,SW)+
				               src(x,y,W )+
				               src(x,y,NW);

				velField(x,y,X)=(src(x,y,N )+
				                src(x,y,NE)+
				                src(x,y,NW)-
				                src(x,y,SE)-
				                src(x,y,S )-
				                src(x,y,SW))
									/densField(x,y);

				velField(x,y,Y)=(src(x,y,NE)+
				                 src(x,y,E )+
				                 src(x,y,SE)-
				                 src(x,y,SW)-
				                 src(x,y,W )-
				                 src(x,y,NW))
									/densField(x,y);
			}
		}


		//collide it
		for (uint x=0; x < inp.sizex; x++)
		{
			for (uint y=0; y < inp.sizey; y++)
			{
				if (flagField(x,y) != 1/*FLUID*/) continue;

				dst(x  ,y  ,C )=src(x,y,C );
				dst(x+1,y  ,N )=src(x,y,N );
				dst(x+1,y+1,NE)=src(x,y,NE);
				dst(x  ,y+1,E )=src(x,y,E );
				dst(x-1,y+1,SE)=src(x,y,SE);
				dst(x-1,y  ,S )=src(x,y,S );
				dst(x-1,y-1,SW)=src(x,y,SW);
				dst(x  ,y-1,W )=src(x,y,W );
				dst(x+1,y-1,NW)=src(x,y,NW);
			}
		}




	}

   return EXIT_SUCCESS;
}

