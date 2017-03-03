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

void  EIGEN_EIGENCLASS::copyFromSEXP ( SEXP S_nimList_ ) {
	SEXP S_pxData;
	SEXP S_values;
	SEXP S_vectors;
	RObjectPointer =  S_nimList_;
	PROTECT(S_pxData = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
	PROTECT(S_values = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("values")));
	PROTECT(S_vectors = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("vectors")));
	SEXP_2_NimArr<1>(S_values, values);
	SEXP_2_NimArr<2>(S_vectors, vectors);
	UNPROTECT(3);
}

/*EIGEN_SVD class functions below */
SEXP  EIGEN_SVDCLASS::copyToSEXP (  )  {
	SEXP S_pxData;
	SEXP S_d;
	SEXP S_u;
	SEXP S_v;

	if (!RCopiedFlag){
		PROTECT(S_pxData = allocVector(STRSXP, 1));
		SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
		PROTECT(S_d = NimArr_2_SEXP<1>(d));
		PROTECT(S_u = NimArr_2_SEXP<2>(u));
		PROTECT(S_v = NimArr_2_SEXP<2>(v));
		defineVar(install("d"), S_d, GET_SLOT(RObjectPointer, S_pxData));
		defineVar(install("u"), S_u, GET_SLOT(RObjectPointer, S_pxData));
		defineVar(install("v"), S_v, GET_SLOT(RObjectPointer, S_pxData));
		RCopiedFlag = true;
		UNPROTECT(4);
	}
	return(RObjectPointer);
}

void  EIGEN_SVDCLASS::createNewSEXP (  )  {
	SEXP S_newNimList;
	SEXP S_listName;
	PROTECT(S_listName = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_listName, 0, mkChar("EIGEN_SVDCLASS"));
	PROTECT(S_newNimList = makeNewNimbleList(S_listName));
	RObjectPointer = S_newNimList;
	UNPROTECT(2);
}

void  EIGEN_SVDCLASS::copyFromSEXP ( SEXP S_nimList_ ) {
	SEXP S_pxData;
	SEXP S_d;
	SEXP S_v;
	SEXP S_u;
	RObjectPointer =  S_nimList_;
	PROTECT(S_pxData = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
	PROTECT(S_d = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("d")));
	PROTECT(S_v = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("v")));
	PROTECT(S_u = findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("u")));
	SEXP_2_NimArr<1>(S_d, d);
	SEXP_2_NimArr<2>(S_v, v);
	SEXP_2_NimArr<2>(S_u, u);
	UNPROTECT(4);
}

