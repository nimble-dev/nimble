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
  virtual ~NamedObjects() {};
};

extern "C" {
  SEXP getModelObjectPtr(SEXP Sextptr, SEXP Sname); /* should rename to getObjectPtr*/
  SEXP getAvailableNames(SEXP Sextptr);
}

class NumberedObjects {
public:
	vector<void*> numberedObjects;
	void* getObjectPtr(int index);
	void setObjectPtr(int index, void* newPtr);
	void resize(int size);
	virtual ~NumberedObjects(){};
};

extern "C" {
	SEXP getNumberedObject(SEXP Snp, SEXP index);
	SEXP setNumberedObject(SEXP Snp, SEXP index, SEXP val);
	SEXP resizeNumberedObjects(SEXP Snp, SEXP size);
	SEXP getSizeNumberedObjects(SEXP Snp);
	SEXP newNumberedObjects();
}

void numberedObjects_Finalizer(SEXP Snp);

#endif
