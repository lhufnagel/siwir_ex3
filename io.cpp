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
   map<string,string>::iterator it;
   for (it = entries.begin();it!=entries.end();it++)
   {
      cout << "Key: " << it->first <<" v: " << it->second << endl;
   }
}


void writeOutput(V_Field& velField,D_Field& densField,Flags& flagField,Input& inp,uint tStep) 
{

   ofstream file; 
   string s;
   ostringstream outStream;
   outStream << tStep;
   s = outStream.str();
   file.open((inp.vtk_file + s + ".vtk").c_str(), ios::out);
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
