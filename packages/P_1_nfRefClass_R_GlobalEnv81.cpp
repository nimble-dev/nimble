#ifndef __P_1_NFREFCLASS_R_GLOBALENV81_CPP
#define __P_1_NFREFCLASS_R_GLOBALENV81_CPP
#include <nimble/EigenTypedefs.h>
#include <iostream>
#include <Rmath.h>
#include <math.h>
#include <nimble/Utils.h>
#include <nimble/accessorClasses.h>
#include <nimble/RcppUtils.h>
#include "P_1_nfRefClass_R_GlobalEnv81.h"
#undef eval

void  nfRefClass_R_GlobalEnv80::operator() (  )  {
_nimble_global_output <<"This is A"<<A<<"\n"<<"\n"; nimble_print_to_R(_nimble_global_output);
}
 nfRefClass_R_GlobalEnv80::nfRefClass_R_GlobalEnv80 (  )  {
namedObjects["A"]=&A;
namedObjects["B"]=&B;
}

SEXP  new_nfRefClass_R_GlobalEnv80 (  )  {
nfRefClass_R_GlobalEnv80 * newObj;
SEXP Sans;
newObj = new nfRefClass_R_GlobalEnv80;
PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
R_RegisterCFinalizerEx(Sans, &namedObjects_Finalizer, FALSE);
UNPROTECT(1);
return(Sans);
}

SEXP  CALL_nfRefClass_R_GlobalEnv80_operator_oP_cP ( SEXP SextPtrToObject )  {
GetRNGstate();
static_cast<nfRefClass_R_GlobalEnv80*>(R_ExternalPtrAddr(SextPtrToObject))->operator()();
PutRNGstate();
return(R_NilValue);
}

void  nfRefClass_R_GlobalEnv81::operator() (  )  {
Map<MatrixXd> Eig_nfVar_oPnf1_B_cP(0,0,0);
_nimble_global_output <<"accessing A:"<<(*nf1).A<<"\n"; nimble_print_to_R(_nimble_global_output);
(*nf1).B(1, 1) = -(1000);
new (&Eig_nfVar_oPnf1_B_cP) Map< MatrixXd >((*nf1).B.getPtr(),(*nf1).B.dim()[0],(*nf1).B.dim()[1]);
_nimble_global_output <<"accessing B:"<<Eig_nfVar_oPnf1_B_cP<<"\n"; nimble_print_to_R(_nimble_global_output);
}
 nfRefClass_R_GlobalEnv81::nfRefClass_R_GlobalEnv81 (  )  {
namedObjects["nf1"]=&nf1;
}

SEXP  new_nfRefClass_R_GlobalEnv81 (  )  {
nfRefClass_R_GlobalEnv81 * newObj;
SEXP Sans;
newObj = new nfRefClass_R_GlobalEnv81;
PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
R_RegisterCFinalizerEx(Sans, &namedObjects_Finalizer, FALSE);
UNPROTECT(1);
return(Sans);
}

SEXP  CALL_nfRefClass_R_GlobalEnv81_operator_oP_cP ( SEXP SextPtrToObject )  {
GetRNGstate();
static_cast<nfRefClass_R_GlobalEnv81*>(R_ExternalPtrAddr(SextPtrToObject))->operator()();
PutRNGstate();
return(R_NilValue);
}
#endif
