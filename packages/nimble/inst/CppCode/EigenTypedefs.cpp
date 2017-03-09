#include "nimble/EigenTypedefs.h"
#include "nimble/dists.h"
#include "nimble/nimDists.h"


vector<int> getSEXPdims(SEXP Sx) {
  if(!isNumeric(Sx)) {PRINTF("Error, getSEXPdims called for something not numeric\n"); return(vector<int>());}
  if(!isVector(Sx)) {PRINTF("Error, getSEXPdims called for something not vector\n"); return(vector<int>());}
  if(!isArray(Sx) & !isMatrix(Sx)) {
    vector<int> ans; 
    ans.resize(1); ans[0] = LENGTH(Sx); return(ans);
  }
  return(SEXP_2_vectorInt(getAttrib(Sx, R_DimSymbol), 0));
}

template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, double> &ans) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_NimArr<1> called for SEXP that is not a numeric or logical!\n");
  int nn = LENGTH(Sn);
  if(ans.size() != 0) PRINTF("Error: trying to reset a NimArr that was already sized\n");
  ans.setSize(nn);
  if(isReal(Sn)) {
     std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());	
  } else {
    if(isInteger(Sn) || isLogical(Sn)) {
      int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
      for(int i = 0; i < nn; ++i) {
	ans(i) = static_cast<double>(iSn[i]);
      }
    } else {
      PRINTF("Error: We could not handle the R input type to SEXP_2_NimArr<1>\n");
    }
  }
}

// Actually this is identical to above so could be done without specialization
template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, int> &ans) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_NimArr<1> called for SEXP that is not a numeric or logical!\n");
  int nn = LENGTH(Sn);
  if(ans.size() != 0) PRINTF("Error: trying to reset a NimArr that was already sized\n");
  ans.setSize(nn);
  if(isReal(Sn)) {
     std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());	
  } else {
    if(isInteger(Sn) || isLogical(Sn)) {
      int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
      for(int i = 0; i < nn; ++i) {
	ans(i) = static_cast<double>(iSn[i]);
      }
    } else {
      PRINTF("Error: We could not handle the R input type to SEXP_2_NimArr<1>\n");
    }
  }
} 

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

extern "C" {
SEXP C_nimEigen(SEXP S_x, SEXP S_valuesOnly, SEXP returnList) {
	NimArr<2, double> x;
	bool valuesOnly;
	SEXP_2_NimArr<2>(S_x, x);
	valuesOnly = SEXP_2_bool(S_valuesOnly);
	Map<MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
    EIGEN_EIGENCLASS C_eigenClass = *EIGEN_EIGEN(Eig_x, valuesOnly);
	C_eigenClass.RObjectPointer = returnList;
	C_eigenClass.copyToSEXP();
    return(returnList);
}

SEXP C_nimSvd(SEXP S_x, SEXP S_vectors, SEXP returnList) {
	NimArr<2, double> x;
	int vectors = SEXP_2_int(S_vectors, 0, 0);
	SEXP_2_NimArr<2>(S_x, x);
	Map<MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
	EIGEN_SVDCLASS C_svdClass = *EIGEN_SVD(Eig_x, vectors);
	C_svdClass.RObjectPointer = returnList;
	C_svdClass.copyToSEXP();
    return(C_svdClass.RObjectPointer);
}

} 


SEXP makeNewNimbleList(SEXP S_listName) {
  SEXP call;
  PROTECT(call = allocVector(LANGSXP, 2));
  SETCAR(call, install("makeNewNimListSEXPRESSIONFromC"));
  SETCADR(call, S_listName);
  UNPROTECT(1);
  return(EVAL(call));
}