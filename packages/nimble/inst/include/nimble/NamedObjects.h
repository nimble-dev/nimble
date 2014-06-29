#ifndef __NAMEDOBJECTS
#define __NAMEDOBJECTS

#include<map>
#include<string>
using namespace std;
#include "RcppUtils.h"
#include<Rinternals.h>

/*
Class for void pointers to arbitrary objects for R access
*/
class NamedObjects {
public:
  map< string, void * > namedObjects;
  virtual void* getObjectPtr( string &name );
};

extern "C" {
  SEXP getModelObjectPtr(SEXP Sextptr, SEXP Sname); /* should rename to getObjectPtr*/
  SEXP getAvailableNames(SEXP Sextptr);
}




#endif
