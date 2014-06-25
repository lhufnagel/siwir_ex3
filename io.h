#pragma once
#include "grid.h"
#include "directions.h"
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include <typeinfo>

using namespace std;

using namespace lbm;

typedef unsigned int uint;

class FileReader
{
   public:
      FileReader();
      ~FileReader();

      void readParameters(const char* filename); 
      template<typename T> T getParameter(const std::string& key) const; 
   private:
      std::map<std::string,std::string> entries; 
};


template< typename T> T FileReader::getParameter(const std::string & key) const
{
   if (entries.count(key)==0)
      return T();

   T ret;
   stringstream ss;

   ss << entries.find(key)->second;
   ss >> ret;

   return ret;
}

void writeOutput(V_Field& velField,D_Field& densField,Flags& flagField,string& vtk_file,uint sizex,uint sizey,uint index); 
