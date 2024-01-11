/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

#include <nimble/EigenTypedefs.h>
#include <nimble/eigenUsingClasses.h>
#include <nimble/RcppNimbleUtils.h>
//#include "nimble/dists.h"
//#include "nimble/nimDists.h"


template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, double> &ans) {
  NIM_ASSERT2(Rf_isNumeric(Sn) || Rf_isLogical(Sn),
    "SEXP_2_NimArr<1, double> called for SEXP that is not a numeric or logical: actual type %s\n",
    Rf_type2char(TYPEOF(Sn)));
  int nn = LENGTH(Sn);
  NIM_ASSERT1(ans.size() == 0, "trying to reset a NimArr that was already sized\n");
  ans.setSize(nn);
  if(Rf_isReal(Sn)) {
     std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());
  } else {
    NIM_ASSERT2(Rf_isInteger(Sn) || Rf_isLogical(Sn),
      "could not handle input of type %s to SEXP_2_NimArr<1, double>\n",
      Rf_type2char(TYPEOF(Sn)));
    int *iSn = Rf_isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    for(int i = 0; i < nn; ++i) {
      ans(i) = static_cast<double>(iSn[i]);
    }
  }
}

// Actually this is identical to above so could be done without specialization
template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, int> &ans) {
  NIM_ASSERT2(Rf_isNumeric(Sn) || Rf_isLogical(Sn),
    "SEXP_2_NimArr<1, int> called for SEXP that is not a numeric or logical: actual type %s\n",
    Rf_type2char(TYPEOF(Sn)));
  int nn = LENGTH(Sn);
  NIM_ASSERT1(ans.size() == 0, "trying to reset a NimArr that was already sized\n");
  ans.setSize(nn);
  if(Rf_isReal(Sn)) {
     std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());
  } else {
    NIM_ASSERT2(Rf_isInteger(Sn) || Rf_isLogical(Sn),
      "could not handle input of type %s to SEXP_2_NimArr<1, int>\n",
      Rf_type2char(TYPEOF(Sn)));
    int *iSn = Rf_isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
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

	PROTECT(S_pxData = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData"));
	PROTECT(S_values = NimArr_2_SEXP<1>(values));
	PROTECT(S_vectors = NimArr_2_SEXP<2>(vectors));
	Rf_defineVar(Rf_install("values"), S_values, PROTECT(GET_SLOT(RObjectPointer, S_pxData)));
	Rf_defineVar(Rf_install("vectors"), S_vectors, PROTECT(GET_SLOT(RObjectPointer, S_pxData)));
	UNPROTECT(5);
	
	return(RObjectPointer);
}

void  EIGEN_EIGENCLASS_R::createNewSEXP (  )  {
	SEXP S_newNimList;
	SEXP S_listName;
	PROTECT(S_listName = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(S_listName, 0, Rf_mkChar("EIGEN_EIGENCLASS"));
	PROTECT(S_newNimList = makeNewNimbleList(S_listName));
	RObjectPointer = S_newNimList;
	UNPROTECT(2);
}

void  EIGEN_EIGENCLASS_R::copyFromSEXP ( SEXP S_nimList_ ) {
	SEXP S_pxData;
	SEXP S_values;
	SEXP S_vectors;
	RObjectPointer =  S_nimList_;
	PROTECT(S_pxData = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData"));
	PROTECT(S_values = Rf_findVarInFrame(PROTECT(GET_SLOT(S_nimList_, S_pxData)), Rf_install("values")));
	PROTECT(S_vectors = Rf_findVarInFrame(PROTECT(GET_SLOT(S_nimList_, S_pxData)), Rf_install("vectors")));
	SEXP_2_NimArr<1>(S_values, values);
	SEXP_2_NimArr<2>(S_vectors, vectors);
	UNPROTECT(5);
}


/*EIGEN_SVD class functions below */
 SEXP  EIGEN_SVDCLASS_R::copyToSEXP (  )  {
 	SEXP S_pxData;
 	SEXP S_d;
 	SEXP S_u;
 	SEXP S_v;

 	PROTECT(S_pxData = Rf_allocVector(STRSXP, 1));
 	SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData"));
 	PROTECT(S_d = NimArr_2_SEXP<1>(d));
 	PROTECT(S_u = NimArr_2_SEXP<2>(u));
 	PROTECT(S_v = NimArr_2_SEXP<2>(v));
 	Rf_defineVar(Rf_install("d"), S_d, PROTECT(GET_SLOT(RObjectPointer, S_pxData)));
 	Rf_defineVar(Rf_install("u"), S_u, PROTECT(GET_SLOT(RObjectPointer, S_pxData)));
 	Rf_defineVar(Rf_install("v"), S_v, PROTECT(GET_SLOT(RObjectPointer, S_pxData)));
 	UNPROTECT(7);
		
 	return(RObjectPointer);
 }

 void  EIGEN_SVDCLASS_R::createNewSEXP (  )  {
 	SEXP S_newNimList;
 	SEXP S_listName;
 	PROTECT(S_listName = Rf_allocVector(STRSXP, 1));
 	SET_STRING_ELT(S_listName, 0, Rf_mkChar("EIGEN_SVDCLASS"));
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
 	PROTECT(S_pxData = Rf_allocVector(STRSXP, 1));
 	SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData"));
 	PROTECT(S_d = Rf_findVarInFrame(PROTECT(GET_SLOT(S_nimList_, S_pxData)), Rf_install("d")));
 	PROTECT(S_v = Rf_findVarInFrame(PROTECT(GET_SLOT(S_nimList_, S_pxData)), Rf_install("v")));
 	PROTECT(S_u = Rf_findVarInFrame(PROTECT(GET_SLOT(S_nimList_, S_pxData)), Rf_install("u")));
 	SEXP_2_NimArr<1>(S_d, d);
 	SEXP_2_NimArr<2>(S_v, v);
 	SEXP_2_NimArr<2>(S_u, u);
 	UNPROTECT(7);
 }

SEXP C_nimEigen(SEXP S_x, SEXP S_symmetric, SEXP S_valuesOnly, SEXP returnList) {
  int* dims = INTEGER(Rf_getAttrib(S_x, R_DimSymbol));
  if(!Rf_isMatrix(S_x) || dims[0] != dims[1])
    RBREAK("Error (C_nimEigen): 'x' must be a square matrix.\n");
  NimArr<2, double> x;
  bool valuesOnly;
  bool symmetric;
  SEXP_2_NimArr<2>(S_x, x);
  valuesOnly = SEXP_2_bool(S_valuesOnly);
  symmetric = SEXP_2_bool(S_symmetric);
  Eigen::Map<Eigen::MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
  EIGEN_EIGENCLASS_R C_eigenClass = *EIGEN_EIGEN_R(Eig_x, symmetric, valuesOnly);
  C_eigenClass.RObjectPointer = returnList;
  C_eigenClass.copyToSEXP();
  return(returnList);
}

SEXP C_nimSvd(SEXP S_x, SEXP S_vectors, SEXP returnList) {
  if(!Rf_isMatrix(S_x))
    RBREAK("Error (C_nimSvd): 'x' must be a matrix.\n");
  NimArr<2, double> x;
  int vectors = SEXP_2_int(S_vectors, 0);
  SEXP_2_NimArr<2>(S_x, x);
  Eigen::Map<Eigen::MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
  EIGEN_SVDCLASS_R C_svdClass = *EIGEN_SVD_R(Eig_x, vectors);
  C_svdClass.RObjectPointer = returnList;
  C_svdClass.copyToSEXP();
  return(C_svdClass.RObjectPointer);
}
