#include "reg.h"

void R_init_code(DllInfo *dll)
{
  value = 1234;
  R_registerRoutines(dll,  NULL, CallEntries, NULL, NULL);
}


SEXP getValue() {
  SEXP Sans;
  PROTECT(Sans = allocVector(REALSXP, 1));
  REAL(Sans)[0] = value;
  UNPROTECT(1);
  return(Sans);
}
