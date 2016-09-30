#include "code.h"

SEXP foo(SEXP Sa) {
  SEXP Sans;
  PROTECT(Sans = allocVector(REALSXP, 1));
  REAL(Sans)[0] = REAL(Sa)[0] + 1;
  UNPROTECT(1);
  return(Sans);
}

