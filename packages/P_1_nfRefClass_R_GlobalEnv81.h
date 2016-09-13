#ifndef __P_1_NFREFCLASS_R_GLOBALENV81
#define __P_1_NFREFCLASS_R_GLOBALENV81
#include <nimble/NamedObjects.h>
#include <Rinternals.h>
#include <nimble/NimArr.h>
#include <nimble/accessorClasses.h>
#include <nimble/nimDists.h>
#undef eval

class nfRefClass_R_GlobalEnv80 : public NamedObjects {
public:
  double A;
  NimArr<2, double> B;
void  operator() (  );
 nfRefClass_R_GlobalEnv80 (  );
};

extern "C" SEXP  new_nfRefClass_R_GlobalEnv80 (  );

extern "C" SEXP  CALL_nfRefClass_R_GlobalEnv80_operator_oP_cP ( SEXP SextPtrToObject );

class nfRefClass_R_GlobalEnv81 : public NamedObjects {
public:
  nfRefClass_R_GlobalEnv80 * nf1;
void  operator() (  );
 nfRefClass_R_GlobalEnv81 (  );
};

extern "C" SEXP  new_nfRefClass_R_GlobalEnv81 (  );

extern "C" SEXP  CALL_nfRefClass_R_GlobalEnv81_operator_oP_cP ( SEXP SextPtrToObject );
#endif
