#ifndef __NAMEDOBJECTS
#define __NAMEDOBJECTS

#include<map>
#include<string>
using namespace std;
#include "RcppNimbleUtils.h"
#include<Rinternals.h>

/*
Class for void pointers to arbitrary objects for R access
*/
class NamedObjects {
public:
  map< string, void * > namedObjects;
  map< string, void * > &getNamedObjects() {return(namedObjects);}
  virtual void* getObjectPtr( string &name );
  virtual ~NamedObjects() {//PRINTF("In NamedObjects destructor\n");
  };
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
  SEXP register_namedObjects_Finalizer(SEXP Sno, SEXP Dll, SEXP Slabel);
  SEXP register_numberedObjects_Finalizer(SEXP Sno, SEXP Dll, SEXP Slabel);
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
