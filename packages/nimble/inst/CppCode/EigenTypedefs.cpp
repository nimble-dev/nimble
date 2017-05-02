#include <nimble/EigenTypedefs.h>

SEXP  EIGEN_EIGENCLASS::copyToSEXP (  )  {
  if(!RCopiedFlag) {
    EIGEN_EIGENCLASS_R::copyToSEXP();
    RCopiedFlag = true;
  }
  return(RObjectPointer);
}

void EIGEN_EIGENCLASS::resetFlags () {
  RCopiedFlag = false;	
}

EIGEN_EIGENCLASS::EIGEN_EIGENCLASS(){
  //std::cout<<"Constructing EIGEN_EIGENCLASS\n";
  namedObjects["values"]=&values;
  namedObjects["vectors"]=&vectors;
  RCopiedFlag = false;
}

SEXP  EIGEN_SVDCLASS::copyToSEXP (  )  {
  if(!RCopiedFlag) {
    EIGEN_SVDCLASS_R::copyToSEXP();
    RCopiedFlag = true;
  }
  return(RObjectPointer);
}

EIGEN_SVDCLASS::EIGEN_SVDCLASS(){
  //std::cout<<"Constructing EIGEN_SVDCLASS\n";
  namedObjects["d"]=&d;
  namedObjects["u"]=&u;
  namedObjects["v"]=&v;
  RCopiedFlag = false;
}

void EIGEN_SVDCLASS::resetFlags () {
  RCopiedFlag = false;	
}

SEXP  NIMBLE_ADCLASS::copyToSEXP (  )  {
	SEXP S_pxData;
	SEXP S_value;
	SEXP S_gradient;
	SEXP S_hessian;
	SEXP S_thirdDerivs;

	PROTECT(S_pxData = allocVector(STRSXP, 1));
	SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
	PROTECT(S_value = NimArr_2_SEXP<1>(value));
	PROTECT(S_gradient = NimArr_2_SEXP<2>(gradient));
	PROTECT(S_hessian = NimArr_2_SEXP<3>(hessian));
	PROTECT(S_thirdDerivs = NimArr_2_SEXP<4>(thirdDerivs));

	defineVar(install("value"), S_value, GET_SLOT(RObjectPointer, S_pxData));
	defineVar(install("gradient"), S_gradient, GET_SLOT(RObjectPointer, S_pxData));
	defineVar(install("hessian"), S_hessian, GET_SLOT(RObjectPointer, S_pxData));
	defineVar(install("thirdDerivs"), S_thirdDerivs, GET_SLOT(RObjectPointer, S_pxData));

	UNPROTECT(5);
	
	return(RObjectPointer);
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
  RObjectPointer = NULL;
}


void  NIMBLE_ADCLASS::createNewSEXP (  )  {
    printf("creating new ADCLASS\n");
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
