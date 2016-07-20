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
// FUN(setNodeModelPtr, 3),
 FUN(getAvailableNames, 1),
 FUN(getMVElement, 2),
 FUN(setMVElement, 3),
 FUN(resizeManyModelVarAccessor, 2),
 FUN(resizeManyModelValuesAccessor, 2),
 FUN(getModelObjectPtr, 2),
 // FUN(getModelValuesMemberElement, 2),
 FUN(getNRow, 1),
 FUN(addBlankModelValueRows, 2),
 FUN(copyModelValuesElements, 4),
 FUN(C_dwish_chol, 5),
 FUN(C_rwish_chol, 3),
 FUN(C_ddirch, 3),
 FUN(C_rdirch, 1),
 FUN(C_dmulti, 4),
 FUN(C_rmulti, 2),
 FUN(C_dcat, 3),
 FUN(C_rcat, 2),
 FUN(C_dt_nonstandard, 5),
 FUN(C_rt_nonstandard, 4),
 FUN(C_dmnorm_chol, 5),
 FUN(C_rmnorm_chol, 3),
 FUN(C_dinterval, 4),
 FUN(C_rinterval, 3),
 FUN(makeNumericList, 3),
//   The following 4 conflict with names of R functions. So we prefix them with a R_
 CFUN(setPtrVectorOfPtrs, 3),
 CFUN(setOnePtrVectorOfPtrs, 3),
 CFUN(setDoublePtrFromSinglePtr, 2),
// CFUN(setSinglePtrFromSinglePtr, 2),
// FUN(newModelValues, 1),
 {NULL, NULL, 0}
};

/* Arrange to register the routines so that they can be checked and also accessed by 
   their equivalent name as an R symbol, e.g. .Call(setNodeModelPtr, ....)
 */
extern "C"
void
R_init_nimble(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}


