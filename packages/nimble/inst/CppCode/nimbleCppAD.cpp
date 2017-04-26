/* Definitions only to be included when a nimbleFunction needs CppAD */
#include <nimble/nimbleCppAD.h>

/* nimbleCppADinfoClass is the class to convey information from a nimbleFunction object 
   to generic CppAD driver wrappers like calcGradient.
   Each nimbleFunction enabled for CppAD will have an object of this class. */

/*EIGEN_EIGEN class functions below */
SEXP  NIMBLE_ADCLASS::copyToSEXP (  )  {
	SEXP S_pxData;
	SEXP S_value;
	SEXP S_gradient;
	SEXP S_hessian;
	SEXP S_thirdDerivs;

	PROTECT(S_pxData = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
	PROTECT(S_value = NimArr_2_SEXP<1>(value));
	PROTECT(S_gradient = NimArr_2_SEXP<2>(S_gradient));
	PROTECT(S_hessian = NimArr_2_SEXP<3>(S_hessian));
	PROTECT(S_thirdDerivs = NimArr_2_SEXP<4>(S_thirdDerivs));

	defineVar(install("value"), S_value, GET_SLOT(RObjectPointer, S_pxData));
	defineVar(install("gradient"), S_gradient, GET_SLOT(RObjectPointer, S_pxData));
	defineVar(install("hessian"), S_hessian, GET_SLOT(RObjectPointer, S_pxData));
	defineVar(install("thirdDerivs"), S_thirdDerivs, GET_SLOT(RObjectPointer, S_pxData));

	UNPROTECT(5);
	
	return(RObjectPointer);
}

void  NIMBLE_ADCLASS::createNewSEXP (  )  {
	SEXP S_newNimList;
	SEXP S_listName;
	PROTECT(S_listName = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_listName, 0, mkChar("NIMBLE_ADCLASS"));
	PROTECT(S_newNimList = makeNewNimbleList(S_listName));
	RObjectPointer = S_newNimList;
	UNPROTECT(2);
}

void  NIMBLE_ADCLASS::copyFromSEXP ( SEXP S_nimList_ ) {
	SEXP S_pxData;
	SEXP S_value;
	SEXP S_gradient;
	SEXP S_hessian;
	SEXP S_thirdDerivs;
	RObjectPointer =  S_nimList_;
	PROTECT(S_pxData = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
	PROTECT(S_value = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("value")));
	PROTECT(S_gradient = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("gradient")));
	PROTECT(S_hessian = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("hessian")));
	PROTECT(S_thirdDerivs = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("thirdDerivs")));

	SEXP_2_NimArr<1>(S_value, value);
	SEXP_2_NimArr<2>(S_gradient, gradient);
	SEXP_2_NimArr<3>(S_hessian, hessian);
	SEXP_2_NimArr<4>(S_thirdDerivs, thirdDerivs);

	UNPROTECT(5);
}

void NIMBLE_ADCLASS::resetFlags () {
  RCopiedFlag = false;	
}

NIMBLE_ADCLASS::NIMBLE_ADCLASS(){
  namedObjects["value"]=&value;
  namedObjects["gradient"]=&gradient;
  namedObjects["hessian"]=&hessian;
  namedObjects["thirdDerivs"]=&thirdDerivs;
  RCopiedFlag = false;
}

