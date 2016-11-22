#include "nimble/EigenTypedefs.h"
#include "nimble/dists.h"
#include "nimble/nimDists.h"
#include "nimble/RcppNimbleUtils.h"


void  EIGEN_EIGENCLASS::copyToSEXP ( SEXP S_nimList_ )  {
SEXP S_pxData;
SEXP S_values;
SEXP S_vectors;
PROTECT(S_pxData = allocVector(STRSXP, 1));
SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
PROTECT(S_values = NimArr_2_SEXP<1>(values));
PROTECT(S_vectors = NimArr_2_SEXP<2>(vectors));
defineVar(install("values"), S_values, GET_SLOT(S_nimList_, S_pxData));
defineVar(install("vectors"), S_vectors, GET_SLOT(S_nimList_, S_pxData));
UNPROTECT(3);
}
SEXP  EIGEN_EIGENCLASS::writeToSEXP (  )  {
SEXP S_newNimList;
PROTECT(S_newNimList = makeNewNimbleList());
copyToSEXP(S_newNimList);
UNPROTECT(1);
return(S_newNimList);
}
 EIGEN_EIGENCLASS::EIGEN_EIGENCLASS (  )  {
namedObjects["values"]=&values;
namedObjects["vectors"]=&vectors;
}

SEXP  new_EIGEN_EIGENCLASS (  )  {
EIGEN_EIGENCLASS * newObj;
SEXP Sans;
newObj = new EIGEN_EIGENCLASS;
PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
UNPROTECT(1);
return(Sans);
}
