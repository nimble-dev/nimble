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

EIGEN_EIGENCLASS::EIGEN_EIGENCLASS (  )  {
	namedObjects["values"]=&values;
	namedObjects["vectors"]=&vectors;
	RCopiedFlag = false;
	RObjectPointer = NULL;
}

// SEXP C_nimEigen(SEXP S_x, SEXP S_valuesOnly, SEXP returnList) {
// 	NimArr<2, double> x;
// 	bool valuesOnly;
// 	SEXP_2_NimArr<2>(S_x, x);
// 	valuesOnly = SEXP_2_bool(S_valuesOnly);
// 	nimSmartPtr<EIGEN_EIGENCLASS> C_eigenClass;
// 	Map<MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
// 	C_eigenClass = EIGEN_EIGEN(Eig_x, valuesOnly);
// 	SET_VECTOR_ELT(returnList, 0, NimArr_2_SEXP<1>((*C_eigenClass).values));
// 	SET_VECTOR_ELT(returnList, 1, NimArr_2_SEXP<2>((*C_eigenClass).vectors));
//     return(returnList);
// }

SEXP  new_EIGEN_EIGENCLASS (  )  {
	EIGEN_EIGENCLASS * newObj;
	SEXP Sans;
	newObj = new EIGEN_EIGENCLASS;
	PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
	UNPROTECT(1);
	return(Sans);
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

EIGEN_SVDCLASS::EIGEN_SVDCLASS (  )  {
	namedObjects["d"]=&d;
	namedObjects["u"]=&u;
	namedObjects["v"]=&v;
	RCopiedFlag = false;
	RObjectPointer = NULL;
}

SEXP  new_EIGEN_SVDCLASS (  )  {
	EIGEN_SVDCLASS * newObj;
	SEXP Sans;
	newObj = new EIGEN_SVDCLASS;
	PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
	UNPROTECT(1);
	return(Sans);
}
