#include "io.h"

FileReader::FileReader():
   entries()
{
}

FileReader::~FileReader()
{}

void FileReader::readParameters(const char* filename)
{
   ifstream inputFile;
   inputFile.open(filename);
   if (!inputFile.is_open())
   {
      cerr << "Failed to open " << filename << " for reading. End." << endl;
      exit(EXIT_FAILURE);
   }

   while (!inputFile.eof())
   {
      string key,value;
      inputFile >> key >> value;
      
      if (!key.empty() && !value.empty())
	 entries.insert(pair<string,string>(key,value));
   }
   inputFile.close();
}

void writeOutput(V_Field& velField,D_Field& densField,Flags& flagField,string& vtk_file,uint sizex,uint sizey,uint index) 
{

   ofstream file; 
   string s;
   ostringstream outStream;
   outStream << index;
   s = outStream.str();
   file.open((vtk_file + s + ".vtk").c_str(), ios::out);
   if(!(file.is_open()))
   { 
      cerr << "Failed to open " << (vtk_file + s +".vtk") << " for writing" << endl;
      //exit(EXIT_FAILURE);
      return;
   }
   file << "# vtk DataFile Version 4.0" << endl;
   file << "SiwiRVisFile" << endl;
   file << "ASCII" << endl;
   file << "DATASET STRUCTURED_POINTS" << endl;
   file << "DIMENSIONS " << sizex << " " << sizey << " 1" << endl;
   file << "ORIGIN 0 0 0" << endl;
   file << "SPACING 1 1 1" << endl;
   file << "POINT_DATA "<< sizex*sizey << endl << endl;


   file << "SCALARS flags unsigned_int 1" << endl;
   file << "LOOKUP_TABLE default" << endl;
   for (uint y=1; y < sizey+1; y++)
   {
      for (uint x=1; x < sizex+1; x++)
      {
	 file << flagField(x,y) << endl;
      }
   }
   file << endl;


   file << "SCALARS density double 1" << endl;
   file << "LOOKUP_TABLE default" << endl;
   for (uint y=1; y < sizey+1; y++)
   {
      for (uint x=1; x < sizex+1; x++)
      {
	 file << densField(x,y) << endl;
      }
   }
   file << endl;

   file << "VECTORS velocity double" << endl;
   for (uint y=1; y < sizey+1; y++)
   {
      for (uint x=1; x < sizex+1; x++)
      {
	 //avoiding numerical-precision errors. Border in reference files seems to be at 1e-8
	 file << ((fabs(velField(x,y,X)) > 1e-8) ? velField(x,y,X) : 0.) << " " << ((fabs(velField(x,y,Y)) > 1e-8) ? velField(x,y,Y) : 0.) << " 0" << endl;
      }
   }
   file << endl;

   file.close();
}
