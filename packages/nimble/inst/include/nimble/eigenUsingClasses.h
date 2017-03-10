#ifndef __EIGENUSINGCLASSES
#define __EIGENUSINGCLASSES
#include "NimArr.h"
#include "NamedObjects.h"
#include "smartPtrs.h"

class EIGEN_EIGENCLASS : public NamedObjects, public pointedToBase {
public:
  NimArr<1, double> values; //these will be defined in nimble.so, and nimble.a and again in on-the-fly compilation
  NimArr<2, double> vectors;
  SEXP RObjectPointer;
  bool RCopiedFlag;
  
  SEXP  copyToSEXP (   );
  void  createNewSEXP (  );
  void  copyFromSEXP ( SEXP S_nimList_ );
  EIGEN_EIGENCLASS(){	
    namedObjects["values"]=&values;
    namedObjects["vectors"]=&vectors;
    RCopiedFlag = false;
    RObjectPointer = NULL;
  };
};

class EIGEN_SVDCLASS : public NamedObjects, public pointedToBase {
 public:
  NimArr<1, double> d;
  NimArr<2, double> u;
  NimArr<2, double> v;
  SEXP RObjectPointer;
  bool RCopiedFlag;
  
  SEXP  copyToSEXP (   );
  void  createNewSEXP (  );
  void  copyFromSEXP ( SEXP S_nimList_ );
  EIGEN_SVDCLASS (  ) {
    namedObjects["d"]=&d;
    namedObjects["u"]=&u;
    namedObjects["v"]=&v;
    RCopiedFlag = false;
    RObjectPointer = NULL;
  };
};

extern "C" {
SEXP C_nimEigen(SEXP S_x, SEXP S_valuesOnly, SEXP returnList);
SEXP C_nimSvd(SEXP S_x, SEXP S_vectors, SEXP returnList);
}



#endif
