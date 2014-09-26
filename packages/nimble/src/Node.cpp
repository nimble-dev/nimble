#include "nimble/NimArr.h"
#include "nimble/Node.h"
#include "R.h"
#include "nimble/Utils.h"

SEXP callSimulate(SEXP Sextptr) {
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a a valid external pointer\n");
    return(R_NilValue);
  }
  GetRNGstate();
  static_cast< Node *>(R_ExternalPtrAddr(Sextptr))->simulate();
  PutRNGstate();
  return(R_NilValue);
}

SEXP callCalculate(SEXP Sextptr) {
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a a valid external pointer\n");
    return(R_NilValue);
  }
  double ans = static_cast< Node *>(R_ExternalPtrAddr(Sextptr))->calculate();
  SEXP Sans = allocVector(REALSXP, 1);
  PROTECT(Sans);
  REAL(Sans)[0] = ans;
  UNPROTECT(1);
  return(Sans);
}

SEXP callGetLogProb(SEXP Sextptr) {
    if(!R_ExternalPtrAddr(Sextptr)) {
        PRINTF("Error: Sextptr is not a a valid external pointer\n");
        return(R_NilValue);
    }
    double ans = static_cast< Node *>(R_ExternalPtrAddr(Sextptr))->getLogProb();
    SEXP Sans = allocVector(REALSXP, 1);
    PROTECT(Sans);
    REAL(Sans)[0] = ans;
    UNPROTECT(1);
    return(Sans);
}
