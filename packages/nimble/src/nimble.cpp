#include <Eigen/Dense>
#include <iostream> // must go before other things because R defines a "length" macro
#include <nimble/RcppUtils.h>
#include "nimble/RcppNimbleUtils.h"
//#include <nimble/ModelClassUtils.h>
//#include <nimble/accessorClasses.h>
#include <nimble/nimbleGraph.h>
#include <nimble/EigenTypedefs.h>
#include <nimble/dists.h>

#include <R_ext/Rdynload.h>

#define FUN(name, numArgs) \
  {#name, (DL_FUNC) &name, numArgs}

#define CFUN(name, numArgs) \
  {"R_"#name, (DL_FUNC) &name, numArgs}

R_CallMethodDef CallEntries[] = {
  // {"getModelValuesPtrFromModel", (DL_FUNC) &getModelValuesPtrFromModel, 1},
// FUN(setNodeModelPtr, 3),
// FUN(getAvailableNames, 1),
// FUN(getMVElement, 2),
// FUN(setMVElement, 3),
// FUN(resizeManyModelVarAccessor, 2),
// FUN(resizeManyModelValuesAccessor, 2),
// FUN(getModelObjectPtr, 2),
 // FUN(getModelValuesMemberElement, 2),
 // FUN(getNRow, 1),
 //FUN(addBlankModelValueRows, 2),
 //FUN(copyModelValuesElements, 4),
 FUN(C_dwish_chol, 5),
 FUN(C_rwish_chol, 3),
 FUN(C_dinvwish_chol, 5),
 FUN(C_rinvwish_chol, 3),
 FUN(C_ddirch, 3),
 FUN(C_rdirch, 1),
 FUN(C_dmulti, 4),
 FUN(C_rmulti, 2),
 FUN(C_dcat, 3),
 FUN(C_rcat, 2),
 FUN(C_dt_nonstandard, 5),
 FUN(C_rt_nonstandard, 4),
 FUN(C_qt_nonstandard, 6),
 FUN(C_pt_nonstandard, 6),
 FUN(C_dmnorm_chol, 5),
 FUN(C_rmnorm_chol, 3),
 FUN(C_dinterval, 4),
 FUN(C_rinterval, 3),
 FUN(C_dmvt_chol, 6),
 FUN(C_rmvt_chol, 4),
 FUN(C_dexp_nimble, 3),
 FUN(C_rexp_nimble, 2),
 FUN(C_pexp_nimble, 4),
 FUN(C_qexp_nimble, 4),
 FUN(C_dinvgamma, 4),
 FUN(C_rinvgamma, 3),
 FUN(C_pinvgamma, 5),
 FUN(C_qinvgamma, 5),
 FUN(C_dsqrtinvgamma, 4),
 FUN(C_rsqrtinvgamma, 3),
 FUN(C_dcar_normal, 7),
 FUN(C_nimEigen, 3),
 FUN(C_nimSvd, 3),

 // FUN(makeNumericList, 3),
//   The following 4 conflict with names of R functions. So we prefix them with a R_
 //CFUN(setPtrVectorOfPtrs, 3),
 //CFUN(setOnePtrVectorOfPtrs, 3),
 //CFUN(setDoublePtrFromSinglePtr, 2),
// CFUN(setSinglePtrFromSinglePtr, 2),
// FUN(newModelValues, 1),

 //These don't register well since they are overloaded
// FUN(SEXP_2_double, 3),
 // FUN(double_2_SEXP, 2),
 // FUN(SEXP_2_bool, 3),
 // FUN(bool_2_SEXP, 2),
 // FUN(SEXP_2_int, 3),
 // FUN(int_2_SEXP, 2),
 // FUN(SEXP_2_string, 2),
 // FUN(SEXP_2_stringVector, 2),
 // FUN(string_2_SEXP, 1),
 // FUN(stringVector_2_SEXP, 1)
 FUN(fastMatrixInsert, 4),
 FUN(matrix2ListDouble, 5),
 FUN(matrix2ListInt, 5),
 FUN(C_rankSample, 4),
 FUN(parseVar, 1),

 FUN(setGraph, 7),
 FUN(anyStochDependencies, 1),
 FUN(anyStochParents, 1),
 FUN(getDependencies, 4),
 FUN(getDependencyPathCountOneNode, 2), 
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
    R_useDynamicSymbols(dll, FALSE);
}


