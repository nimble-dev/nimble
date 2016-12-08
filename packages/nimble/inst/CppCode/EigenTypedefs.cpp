#include "nimble/EigenTypedefs.h"
#include "nimble/dists.h"
#include "nimble/nimDists.h"
#include "nimble/RcppNimbleUtils.h"

/*EIGEN_EIGEN class functions below */
SEXP  EIGEN_EIGENCLASS::copyToSEXP (  )  {
SEXP S_pxData;
SEXP S_values;
SEXP S_vectors;
if (!RCopiedFlag){
PROTECT(S_pxData = allocVector(STRSXP, 1));
SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
PROTECT(S_values = NimArr_2_SEXP<1>(values));
PROTECT(S_vectors = NimArr_2_SEXP<2>(vectors));
defineVar(install("values"), S_values, GET_SLOT(RObjectPointer, S_pxData));
defineVar(install("vectors"), S_vectors, GET_SLOT(RObjectPointer, S_pxData));
RCopiedFlag = true;
UNPROTECT(3);
}
return(RObjectPointer);
}

void  EIGEN_EIGENCLASS::createNewSEXP (  )  {
SEXP S_newNimList;
SEXP S_listName;
PROTECT(S_listName = allocVector(STRSXP, 1));
SET_STRING_ELT(S_listName, 0, mkChar("EIGEN_EIGENCLASS"));
PROTECT(S_newNimList = makeNewNimbleList(S_listName));
RObjectPointer = S_newNimList;
UNPROTECT(2);
}

 EIGEN_EIGENCLASS::EIGEN_EIGENCLASS (  )  {
namedObjects["values"]=&values;
namedObjects["vectors"]=&vectors;
RCopiedFlag = false;
RObjectPointer = NULL;
}

SEXP  new_EIGEN_EIGENCLASS (  )  {
EIGEN_EIGENCLASS * newObj;
SEXP Sans;
newObj = new EIGEN_EIGENCLASS;
PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
UNPROTECT(1);
return(Sans);
}

/*EIGEN_SVD class functions below */

