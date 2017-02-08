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
	nimSmartPtr<EIGEN_EIGENCLASS> C_eigenClass;
	Map<MatrixXd> Eig_x(x.getPtr(), x.dim()[0], x.dim()[1]); 
	C_eigenClass = EIGEN_EIGEN(Eig_x, valuesOnly);
	SET_VECTOR_ELT(returnList, 0, NimArr_2_SEXP<1>((*C_eigenClass).values));
	SET_VECTOR_ELT(returnList, 1, NimArr_2_SEXP<2>((*C_eigenClass).vectors));
    return(returnList);
}
}
#endif
