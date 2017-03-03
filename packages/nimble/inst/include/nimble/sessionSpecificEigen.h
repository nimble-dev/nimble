#ifndef _SESSIONSPECIFICEIGEN
#define _SESSIONSPECIFICEIGEN

#include "EigenTypeDefs.h"
#include "nimArr.h"

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

/* SEXP  new_EIGEN_EIGENCLASS (  )  {
	EIGEN_EIGENCLASS * newObj;
	SEXP Sans;
	newObj = new EIGEN_EIGENCLASS;
	PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
	UNPROTECT(1);
	return(Sans);
} 

SEXP  new_EIGEN_SVDCLASS (  )  {
	EIGEN_SVDCLASS * newObj;
	SEXP Sans;
	newObj = new EIGEN_SVDCLASS;
	PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));
	UNPROTECT(1);
	return(Sans);
}  */
}
#endif
