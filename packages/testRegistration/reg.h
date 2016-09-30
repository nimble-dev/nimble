#include "R.h"
#include "Rinternals.h"
#include "code.h"

#include <R_ext/Rdynload.h>

#define FUN(name, numArgs) \
  {#name, (DL_FUNC) &name, numArgs}

#define CFUN(name, numArgs)			\
  {"R_"#name, (DL_FUNC) &name, numArgs}

R_CallMethodDef CallEntries[] = {
  FUN(foo, 1),
  {NULL, NULL, 0}
};

double value(0);

extern "C" SEXP getValue();

extern "C" void R_init_reg(DllInfo *dll);
