#ifndef __VALUES
#define __VALUES

//#include<vector>
#include<string>
using std::string;
#include "NamedObjects.h"
#include "Utils.h"

#ifdef _IN_CPP_CODE

class Values : public NamedObjects {
public:
  int numRows;
  int getsize(){return(numRows);};
  virtual void resize(int nrow)=0;
  
  string buildName;
  string getMVBuildName(){ return buildName;}  ;
  Values() { buildName = "missing"; numRows = 0;};
};

#endif

#endif
