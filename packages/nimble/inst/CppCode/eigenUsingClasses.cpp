#include <nimble/EigenTypedefs.h>
#include <nimble/eigenUsingClasses.h>
#include <nimble/RcppNimbleUtils.h>
//#include "nimble/dists.h"
//#include "nimble/nimDists.h"


template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, double> &ans) {
  NIM_ASSERT2(isNumeric(Sn) || isLogical(Sn),
    "SEXP_2_NimArr<1, double> called for SEXP that is not a numeric or logical: actual type %s\n",
    type2str(TYPEOF(Sn)));
  int nn = LENGTH(Sn);
  NIM_ASSERT1(ans.size() == 0, "trying to reset a NimArr that was already sized\n");
  ans.setSize(nn);
  if(isReal(Sn)) {
     std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());
  } else {
    NIM_ASSERT2(isInteger(Sn) || isLogical(Sn),
      "could not handle input of type %s to SEXP_2_NimArr<1, double>\n",
      type2str(TYPEOF(Sn)));
    int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    for(int i = 0; i < nn; ++i) {
      ans(i) = static_cast<double>(iSn[i]);
    }
  }
}

// Actually this is identical to above so could be done without specialization
template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, int> &ans) {
  NIM_ASSERT2(isNumeric(Sn) || isLogical(Sn),
    "SEXP_2_NimArr<1, int> called for SEXP that is not a numeric or logical: actual type %s\n",
    type2str(TYPEOF(Sn)));
  int nn = LENGTH(Sn);
  NIM_ASSERT1(ans.size() == 0, "trying to reset a NimArr that was already sized\n");
  ans.setSize(nn);
  if(isReal(Sn)) {
     std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());
  } else {
    NIM_ASSERT2(isInteger(Sn) || isLogical(Sn),
      "could not handle input of type %s to SEXP_2_NimArr<1, int>\n",
      type2str(TYPEOF(Sn)));
    int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    for(int i = 0; i < nn; ++i) {
      ans(i) = static_cast<double>(iSn[i]);
    }
  }
}

/*EIGEN_EIGEN class functions below */
SEXP  EIGEN_EIGENCLASS_R::copyToSEXP (  )  {
	SEXP S_pxData;
	SEXP S_values;
	SEXP S_vectors;

	PROTECT(S_pxData = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
	PROTECT(S_values = NimArr_2_SEXP<1>(values));
	PROTECT(S_vectors = NimArr_2_SEXP<2>(vectors));
	defineVar(install("values"), S_values, GET_SLOT(RObjectPointer, S_pxData));
	defineVar(install("vectors"), S_vectors, GET_SLOT(RObjectPointer, S_pxData));
	UNPROTECT(3);
	
	return(RObjectPointer);
}

void  EIGEN_EIGENCLASS_R::createNewSEXP (  )  {
	SEXP S_newNimList;
	SEXP S_listName;
	PROTECT(S_listName = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_listName, 0, mkChar("EIGEN_EIGENCLASS"));
	PROTECT(S_newNimList = makeNewNimbleList(S_listName));
	RObjectPointer = S_newNimList;
	UNPROTECT(2);
}

void  EIGEN_EIGENCLASS_R::copyFromSEXP ( SEXP S_nimList_ ) {
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
 SEXP  EIGEN_SVDCLASS_R::copyToSEXP (  )  {
 	SEXP S_pxData;
 	SEXP S_d;
 	SEXP S_u;
 	SEXP S_v;

 	PROTECT(S_pxData = allocVector(STRSXP, 1));
 	SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
 	PROTECT(S_d = NimArr_2_SEXP<1>(d));
 	PROTECT(S_u = NimArr_2_SEXP<2>(u));
 	PROTECT(S_v = NimArr_2_SEXP<2>(v));
 	defineVar(install("d"), S_d, GET_SLOT(RObjectPointer, S_pxData));
 	defineVar(install("u"), S_u, GET_SLOT(RObjectPointer, S_pxData));
 	defineVar(install("v"), S_v, GET_SLOT(RObjectPointer, S_pxData));
 	UNPROTECT(4);
		
 	return(RObjectPointer);
 }

 void  EIGEN_SVDCLASS_R::createNewSEXP (  )  {
 	SEXP S_newNimList;
 	SEXP S_listName;
 	PROTECT(S_listName = allocVector(STRSXP, 1));
 	SET_STRING_ELT(S_listName, 0, mkChar("EIGEN_SVDCLASS"));
 	PROTECT(S_newNimList = makeNewNimbleList(S_listName));
 	RObjectPointer = S_newNimList;
 	UNPROTECT(2);
 }

 void  EIGEN_SVDCLASS_R::copyFromSEXP ( SEXP S_nimList_ ) {
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

SEXP C_nimEigen(SEXP S_x, SEXP S_valuesOnly, SEXP returnList) {
  int* dims = INTEGER(getAttrib(S_x, R_DimSymbol));
  if(!isMatrix(S_x) || dims[0] != dims[1])
    RBREAK("Error (C_nimEigen): 'x' must be a square matrix.\n");
  NimArr<2, double> x;
  bool valuesOnly;
  SEXP_2_NimArr<2>(S_x, x);
  valuesOnly = SEXP_2_bool(S_valuesOnly);
  Eigen::Map<Eigen::MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
  EIGEN_EIGENCLASS_R C_eigenClass = *EIGEN_EIGEN_R(Eig_x, valuesOnly);
  C_eigenClass.RObjectPointer = returnList;
  C_eigenClass.copyToSEXP();
  return(returnList);
}

SEXP C_nimSvd(SEXP S_x, SEXP S_vectors, SEXP returnList) {
  if(!isMatrix(S_x))
    RBREAK("Error (C_nimSvd): 'x' must be a matrix.\n");
  NimArr<2, double> x;
  int vectors = SEXP_2_int(S_vectors, 0, 0);
  SEXP_2_NimArr<2>(S_x, x);
  Eigen::Map<Eigen::MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
  EIGEN_SVDCLASS_R C_svdClass = *EIGEN_SVD_R(Eig_x, vectors);
  C_svdClass.RObjectPointer = returnList;
  C_svdClass.copyToSEXP();
  return(C_svdClass.RObjectPointer);
}
