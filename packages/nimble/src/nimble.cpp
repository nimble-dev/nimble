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
 FUN(C_dlkj_corr_cholesky, 4),
 FUN(C_rlkj_corr_cholesky, 2),
 FUN(C_dexp_nimble, 3),
 FUN(C_rexp_nimble, 2),
 FUN(C_pexp_nimble, 4),
 FUN(C_qexp_nimble, 4),
 FUN(C_ddexp, 4),
 FUN(C_rdexp, 3),
 FUN(C_pdexp, 5),
 FUN(C_qdexp, 5),
 FUN(C_dinvgamma, 4),
 FUN(C_rinvgamma, 3),
 FUN(C_pinvgamma, 5),
 FUN(C_qinvgamma, 5),
 FUN(C_dsqrtinvgamma, 4),
 FUN(C_rsqrtinvgamma, 3),
 FUN(C_dcar_normal, 8),
 FUN(C_dcar_proper, 10),
 FUN(C_rcar_proper, 9),
 FUN(C_nimEigen, 4),
 FUN(C_nimSvd, 3),
 FUN(fastMatrixInsert, 4),
 FUN(matrix2ListDouble, 5),
 FUN(matrix2ListInt, 5),
 FUN(C_rankSample, 4),
 FUN(parseVar, 1),
 FUN(makeParsedVarList, 1),
 FUN(C_setGraph, 7),
 FUN(C_anyStochDependencies, 1),
 FUN(C_anyStochParents, 1),
 FUN(C_getDependencies, 4),
 FUN(C_getParents, 5),
 FUN(C_getConditionallyIndependentSets, 7),
 FUN(C_getDependencyPaths, 2), 
 FUN(C_getDependencyPathCountOneNode, 2), 
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
