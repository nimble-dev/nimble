#ifndef __NAMEDOBJECTS
#define __NAMEDOBJECTS

#include<map>
#include<string>
using namespace std;
#include "RcppUtils.h"
#include<Rinternals.h>

#ifdef _IN_CPP_CODE

/*
Class for void pointers to arbitrary objects for R access
*/
class NamedObjects {
public:
  map< string, void * > namedObjects;
  virtual void* getObjectPtr( string &name );
  virtual ~NamedObjects() { };
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

void namedObjects_Finalizer(SEXP Sno);


// This is a class which is basically identical to NumberedObjects,
// but it is for items of class <T> which are built on the spot with
// no external pointer. Thus, the finalizer must be specialized for class<T>
template<class T>
class SpecialNumberedObjects : public NumberedObjects{
 public:
  virtual ~SpecialNumberedObjects(){
    int len = numberedObjects.size();
    T* ptr;
    for(int i = 0; i < len; i++){
      ptr = static_cast<T*>(getObjectPtr(i));
      if(ptr != 0){
	delete ptr;	
      }
    }
  }
};

template<class T>
void Special_NumberedObjects_Finalizer(SEXP Snp);


#endif

#endif
