#include "main.h"
#include "grid.h"
#include "io.h"
#include "flags.h"
#include "directions.h"

using namespace lbm;

int main(int argc, char* argv[])
{
   if (argc != 2) {
      cerr<<"error: wrong number of arguments"<<endl;
      cout<<"call ./lbm params.dat "<<endl;
      exit(EXIT_FAILURE);
   }

   double uw[]={0.08,0.}; //x and y

   FileReader reader;
   reader.readParameters(argv[1]); 
   uint sizex = reader.getParameter<uint>("sizex");
   uint sizey = reader.getParameter<uint>("sizey");
   uint timesteps = reader.getParameter<uint>("timesteps");
   assert(timesteps!=0);
   if (timesteps > 100000) cerr << "Attention: Calculating many timesteps (100000+). Expect Latency" << endl;
   double omega = reader.getParameter<double>("omega");
   if (omega < 0.5 || omega > 1.95) cerr << "Attention: omega not in range of stable simulation (0.5-1.95)" << endl;

   string vtk_file = reader.getParameter<string>("vtk_file");
   assert(!vtk_file.empty());
   uint vtk_step = reader.getParameter<uint>("vtk_step");
   if (0==vtk_step) cerr << "Attention: not writing vtk-output" << endl;

   string geometry = reader.getParameter<string>("geometry");


   gray** geomVals = NULL;
   gray maxVal = 255;
   if (!geometry.empty())
   {
      pgm_init(&argc,argv);
      int rows, cols;
      int format;

      FILE* geomFile;
      if ((geomFile = fopen(geometry.c_str(),"r")) == NULL)
      {
	 cerr << "Failed to open geometry file " << geometry << endl;
	 return EXIT_FAILURE;
      }

      pgm_readpgminit(geomFile,&cols,&rows,&maxVal,&format);

      geomVals = pgm_allocarray( cols,  rows );

      //PGMs horizontal gespiegelt einlesen, 
      //damit wir unseren koordinatenursprung unten links (statt oben rehchts)  haben
      for (int y=rows-1;y>=0;y--){
	 pgm_readpgmrow( geomFile, geomVals[y], cols, maxVal, format );
      }

      sizex = cols;
      sizey = rows;
   }

   assert(sizex < 1000 && sizey < 1000);

   PDF_Field src(sizex+2,sizey+2),
	     dst(sizex+2,sizey+2);
   V_Field velField(sizex+2,sizey+2);
   D_Field densField(sizex+2,sizey+2); //macroscopic
   Flags flagField(sizex+2,sizey+2); //info about wether boundary or fluid cell
   uint cellCount(0); //for MLUP calculation
   uint writeIndex(0);


   //initialize fields
   //todo probably put this in function
   for (uint x=0; x < sizex+2; x++)
   {
      for (uint y=0; y < sizey+2; y++)
      {
	 flagField(x,y)= NOSLIP;

	 if (x>0 && y>0 &&
	       x<sizex+1 && y<sizey+1 &&
	       (geomVals == NULL ||
	       geomVals[y-1][x-1] == maxVal))
	 {
	    flagField(x,y)= FLUID;

	    velField(x,y,X)=velField(x,y,Y)=0.;

	    densField(x,y)=1.;
	    cellCount++;
	 }

	 for (uint i=0;i<stencilSize;i++)
	    src(x,y,i) = dst(x,y,i) = weigths[i];

      }
   }

   //initialize boundary field
   for (uint y=1; y < sizey+1; y++)
   {
      flagField(0,y)=flagField(sizex+1,y)=NOSLIP;
   }
   for (uint x=0; x < sizex+2; x++)
   {
      flagField(x,0)=NOSLIP;
      flagField(x,sizey+1)=MOVNOSLIP;
   }

    // time measurements
    struct timeval start, end;
    long seconds, useconds;
    gettimeofday(&start, NULL); 

//Start lbm here
   for (uint tStep = 0; tStep < /* or <= ? */ timesteps; ++tStep)
   {
      //Streaming
#pragma omp parallel for
      for (uint x=0; x < sizex+2; x++)
      {
	 for (uint y=0; y < sizey+2; y++)
	 {
	    if (flagField(x,y) == FLUID)
	    {
	       for (uint i=0;i<stencilSize;i++)
		  dst(x + stenNbs[X][i],y + stenNbs[Y][i],i) = src(x,y,i);

	    }
	 }
      }


      //Boundary handling
#pragma omp parallel for
      for (uint x=0; x < sizex+2; x++)
      {

	 for (uint y=0; y < sizey+2; y++)
	 {
	    if (flagField(x,y) == NOSLIP)
	    {
	       for (uint i=1;i<stencilSize;i++)
	       {

		  if (
			(x + stenNbs[X][i]) < (sizex+2) && //dont have to check against x+stenNbs[i]>=0 because unsigned int ;)
			(y + stenNbs[Y][i]) < (sizey+2) &&
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
		  if (	(x + stenNbs[X][i])< (sizex+2) &&
			(y + stenNbs[Y][i])< (sizey+2) &&
			(flagField(x + stenNbs[X][i], y + stenNbs[Y][i]) == FLUID))
		  {
		     dst(x + stenNbs[X][i], y + stenNbs[Y][i] ,i) = dst(x,y,(i+3)%8+1) - 6.*weigths[i]*(-stenNbs[X][i]*uw[X]-stenNbs[Y][i]*uw[Y]);
		  }
	       }

	    }

	 }
      }


      
      //Colliding
#pragma omp parallel for
      for (uint x=1; x < sizex+1; x++)
      {
	 for (uint y=1; y < sizey+1; y++)
	 {
	    //calculate macroscopic measures here
	    double rho=dst(x,y,C);
	    for (uint i=1;i<stencilSize;i++)
	       rho += dst(x,y,i);

	    densField(x,y) = rho;


	    double vx=
	       (dst(x,y,E)+
	       dst(x,y,NE )+
	       dst(x,y,SE))-
		(dst(x,y,W)+
		dst(x,y,NW )+
		dst(x,y,SW));

	    velField(x,y,X)=vx;

	    double vy=
	       (dst(x,y,N )+
	       dst(x,y,NW)+
	       dst(x,y,NE))-
		(dst(x,y,S)+
		dst(x,y,SE)+
		dst(x,y,SW));

	    velField(x,y,Y)=vy;

	    //seems like we should not divide by the density here...
	    //otherwise the algorithm gets instable
	    /*
	       velField(x,y,X)/=densField(x,y);
	       velField(x,y,Y)/=densField(x,y);
	       */


	    double vsq = 3./2.*(vx*vx+vy*vy);
	    //actually collide
	    for (uint i=0;i<stencilSize;i++)
	    {
	       double scalProd = stenNbs[X][i]*vx+stenNbs[Y][i]*vy;
	       dst(x,y,i )-= omega*(dst(x,y,i) - weigths[i]*(rho + 3.*(scalProd) + 9./2*scalProd*scalProd  -vsq)); 
	    }
	 }
      }

      dst.swap(src);

      if (vtk_step!=0 && tStep!=0 && (tStep % vtk_step)==0)
      	 writeOutput(velField,densField,flagField,vtk_file, sizex, sizey, writeIndex++); 
   }

    gettimeofday(&end, NULL);
    useconds = end.tv_usec - start.tv_usec;
    seconds = end.tv_sec - start.tv_sec;  
    if(useconds <0){
        useconds+=1000000;
        seconds--;
    }
    double mlups = cellCount*timesteps/(double(seconds+double(useconds/1000000.0))*1000000);
    cout<<"MLUps: "<< mlups <<endl;

   return EXIT_SUCCESS;
}


