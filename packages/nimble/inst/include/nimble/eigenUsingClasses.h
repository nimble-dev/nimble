#ifndef __EIGENUSINGCLASSES
#define __EIGENUSINGCLASSES
#include "NimArr.h"
#include "NamedObjects.h"
#include "smartPtrs.h"


class EIGEN_EIGENCLASS_R : public pointedToBase {
public:
  NimArr<1, double> values; //these will be defined in nimble.so, and nimble.a and again in on-the-fly compilation
  NimArr<2, double> vectors;
  NimArr<1, double> &getValues() {return(values);}
  NimArr<2, double> &getVectors() {return(vectors);}
  SEXP RObjectPointer;
  
  virtual SEXP copyToSEXP (   );
  void  createNewSEXP (  );
  void  copyFromSEXP ( SEXP S_nimList_ );
  EIGEN_EIGENCLASS_R(){	
    RObjectPointer = NULL;
  };
};


class EIGEN_SVDCLASS_R : public pointedToBase {
 public:
  NimArr<1, double> d;
  NimArr<2, double> u;
  NimArr<2, double> v;
  NimArr<1, double> &getD() {return(d);}
  NimArr<2, double> &getU() {return(u);}
  NimArr<2, double> &getV() {return(v);}
  SEXP RObjectPointer;
  
  virtual SEXP  copyToSEXP (   );
  void  createNewSEXP (  );
  void  copyFromSEXP ( SEXP S_nimList_ );
  EIGEN_SVDCLASS_R (  ) {
    RObjectPointer = NULL;
  };
};

extern "C" {
SEXP C_nimEigen(SEXP S_x, SEXP S_symmetric, SEXP S_valuesOnly, SEXP returnList);
SEXP C_nimSvd(SEXP S_x, SEXP S_vectors, SEXP returnList);
}


#endif
