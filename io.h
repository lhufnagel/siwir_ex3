#pragma once
#include "grid.h"
#include "inputParams.h"
#include "directions.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
using namespace std;

using namespace lbm;

class FileReader
{
   public:
      FileReader();
      ~FileReader();

      void readParameters(const char* filename); 
      template<typename T>
      T getParameter(const std::string& key); 
   private:
      std::map<std::string,std::string> entries; 
};

void writeOutput(V_Field& velField,D_Field& densField,Flags& flagField,Input& inp,uint tStep);
