#ifndef __MODELCLASSUTILS
#define __MODELCLASSUTILS

#include "NamedObjects.h"
#include "RcppNimbleUtils.h"
#include "Values.h"

class ModelBase : public NamedObjects {
 public:
  Values* modelValues_;
  virtual Values* getModelValuesPtr() { return modelValues_; }
};

extern "C" {
// gets pointer to modelvalues object in model
SEXP getModelValuesPtrFromModel(SEXP rPtr);
// Gets the ptr to an element of name Sname from
// the ModelValues object pointed to by Sextptr
SEXP getModelElementPtr(SEXP Sextptr, SEXP Sname);

// gets character string to feed to .Call
// that builds a new ModelValue of the same type as
// is pointed to by rPtr
SEXP getMVBuildName(SEXP rPtr);

SEXP derefPtr(SEXP SmultiPtr);
}

// Gets the ptr to an element of name Sname from
// the ModelValues object pointed to by Sextptr
NimArrType** cGetModelElementPtr(SEXP Sextptr, SEXP Sname);

#endif
