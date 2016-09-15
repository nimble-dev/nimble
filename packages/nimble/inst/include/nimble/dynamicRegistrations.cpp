// Code to be included in each on-the-fly (dynamic) nimble compilation
// It will be included by a generated file in each project


#include <nimble/RcppUtils.h>
#include <nimble/ModelClassUtils.h>
#include <nimble/accessorClasses.h>
#include <nimble/dists.h>

#include <R_ext/Rdynload.h>

#define FUN(name, numArgs) \
  {#name, (DL_FUNC) &name, numArgs}

#define CFUN(name, numArgs) \
  {"R_"#name, (DL_FUNC) &name, numArgs}

R_CallMethodDef CallEntries[] = {
  {"getModelValuesPtrFromModel", (DL_FUNC) &getModelValuesPtrFromModel, 1},
}

// Something like this will be generated with each .so/.dll nimble creates
// however it is required to be named R_init_SONAME, so it must be generated for each one.
//  
extern "C"
void
R_init_nimble_routines(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}


