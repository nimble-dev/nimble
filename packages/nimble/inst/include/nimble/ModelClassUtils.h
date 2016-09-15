#ifndef __MODELCLASSUTILS
#define __MODELCLASSUTILS

#include "NamedObjects.h"
#include "RcppNimbleUtils.h"
#include "Values.h"

class ModelBase : public NamedObjects{
	public:
	Values* _modelValues;
	virtual Values* getModelValuesPtr(){ return _modelValues; }
	};

extern "C" {
  SEXP getModelValuesPtrFromModel (SEXP rPtr); // gets pointer to modelvalues object in model
  SEXP getModelElementPtr(SEXP Sextptr, SEXP Sname); // Gets the ptr to an element of name Sname from
  // the ModelValues object pointed to by Sextptr
  
  SEXP getMVBuildName(SEXP rPtr);			// gets character string to feed to .Call 
  // that builds a new ModelValue of the same type as
  // is pointed to by rPtr
  
  SEXP derefPtr(SEXP SmultiPtr);
  SEXP setSinglePtrFromSinglePtr(SEXP SPtrTo, SEXP SPtrToSet);
 }


NimArrType** cGetModelElementPtr(SEXP Sextptr, SEXP Sname);	// Gets the ptr to an element of name Sname from
// the ModelValues object pointed to by Sextptr

#endif
