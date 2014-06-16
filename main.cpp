#include "main.h"
#include "grid.h"
#include "flags.h"
#include "directions.h"
#include "inputParams.h"

using namespace lbm;

void writeOutput(V_Field& velField,D_Field& densField,Flags& flagField,Input& inp,uint tStep) 
{

   ofstream file; 
   string s;
   ostringstream outStream;
   outStream << tStep;
   s = outStream.str();
   file.open((inp.vtk_file + s +".vtk").c_str(), ios::out);
   if(!(file.is_open()))
   { 
      cerr << "Failed to open " << (inp.vtk_file + s +".vtk") << " for writing" << endl;
      //exit(EXIT_FAILURE);
      return;
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


   pgm_init(&argc,argv);
   int rows, cols;
   gray** refrVals = NULL;
   gray maxValref = 255;
   int format;

   FILE* refrFile;
   if ((refrFile = fopen("./test.pgm","r")) == NULL)
   {
      cerr << "fopen(): " << argv[3] << endl;
      return EXIT_FAILURE;
   }

   pgm_readpgminit(refrFile,&cols,&rows,&maxValref,&format);

   refrVals = pgm_allocarray( cols,  rows );

   //PGMs horizontal gespiegelt einlesen, damit wir unseren koordinatenursprung unten links (statt oben rehchts)  haben
   for (int y=rows-1;y>=0;y--){
      pgm_readpgmrow( refrFile, refrVals[y], cols, maxValref, format );
   }
   inp.sizex = cols;
   inp.sizey = rows;

   PDF_Field src(inp.sizex+2,inp.sizey+2),
	     dst(inp.sizex+2,inp.sizey+2);
   V_Field velField(inp.sizex+2,inp.sizey+2);
   D_Field densField(inp.sizex+2,inp.sizey+2); //macroscopic
   Flags flagField(inp.sizex+2,inp.sizey+2); //info about wether boundary or fluid cell
   uint cellCount(0); //for MLUP calculation


   //initialize fields
   //todo probably put this in function
   for (uint x=0; x < inp.sizex+2; x++)
   {
      for (uint y=0; y < inp.sizey+2; y++)
      {
	 flagField(x,y)= NOSLIP;

	 if (x>0 && y>0 &&
	       x<inp.sizex+1 && y<inp.sizey+1 &&
	       refrVals[y-1][x-1] == maxValref)
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
   for (uint y=1; y < inp.sizey+1; y++)
   {
      flagField(0,y)=flagField(inp.sizex+1,y)=NOSLIP;
   }
   for (uint x=0; x < inp.sizex+2; x++)
   {
      flagField(x,0)=NOSLIP;
      flagField(x,inp.sizey+1)=MOVNOSLIP;
   }

    // time measurements
    struct timeval start, end;
    long seconds, useconds;
    gettimeofday(&start, NULL); 

//Start lbm here
   for (uint tStep = 0; tStep<inp.timesteps; ++tStep)
   {
      //Streaming
#pragma omp parallel for
      for (uint x=0; x < inp.sizex+2; x++)
      {
	 for (uint y=0; y < inp.sizey+2; y++)
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
      for (uint x=0; x < inp.sizex+2; x++)
      {

	 for (uint y=0; y < inp.sizey+2; y++)
	 {
	    if (flagField(x,y) == NOSLIP)
	    {
	       for (uint i=1;i<stencilSize;i++)
	       {

		  if (
			(x + stenNbs[X][i]) < (inp.sizex+2) && //dont have to check against x+stenNbs[i]>=0 bc unsigned int ;)
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
		  if (	(x + stenNbs[X][i])< (inp.sizex+2) &&
			(y + stenNbs[Y][i])< (inp.sizey+2) &&
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
      for (uint x=1; x < inp.sizex+1; x++)
      {
	 for (uint y=1; y < inp.sizey+1; y++)
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
	       dst(x,y,i )-= inp.omega*(dst(x,y,i) - weigths[i]*(rho + 3.*(scalProd) + 9./2*scalProd*scalProd  -vsq)); 
	    }
	 }
      }

      dst.swap(src);

      if (inp.vtk_step!=0 &&  tStep!=0 && (tStep % inp.vtk_step)==0)
	 writeOutput(velField,densField,flagField,inp,tStep);
   }

    gettimeofday(&end, NULL);
    useconds = end.tv_usec - start.tv_usec;
    seconds = end.tv_sec - start.tv_sec;  
    if(useconds <0){
        useconds+=1000000;
        seconds--;
    }
    double mlups = cellCount*inp.timesteps/(double(seconds+double(useconds/1000000.0))*1000000);
    cout<<"MLUps: "<< mlups <<endl;

   return EXIT_SUCCESS;
}


