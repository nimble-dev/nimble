#include "nimble/RcppUtils.h"
#include "nimble/Utils.h"
#include<iostream>
#include<sstream>
#include<algorithm>

std::ostringstream _nimble_global_output;

string _NIMBLE_WHITESPACE_UTIL(" \t");
string _NIMBLE_WHITESPACEBRACKET_UTIL(" \t[");

void nimble_print_to_R(std::ostringstream &input) {
  PRINTF("%s", input.str().c_str());
  input.str("");
  input.clear();
}

void multivarTestCall(double *x, int n) {
  _nimble_global_output<<"In multivarTestCall\n";
  for(int iii = 0; iii < n; ++iii) {
    _nimble_global_output<<x[iii]<<" ";
  }
  _nimble_global_output<<"\n";
  nimble_print_to_R(_nimble_global_output);
}

vector<int> getSEXPdims(SEXP Sx) {
  if(!isNumeric(Sx)) {PRINTF("Error, getSEXPdims called for something not numeric\n"); return(vector<int>());}
  if(!isVector(Sx)) {PRINTF("Error, getSEXPdims called for something not vector\n"); return(vector<int>());}
  if(!isArray(Sx) & !isMatrix(Sx)) {
    vector<int> ans; 
    ans.resize(1); ans[0] = LENGTH(Sx); return(ans);
  }
  return(SEXP_2_vectorInt(getAttrib(Sx, R_DimSymbol), 0));
}

string STRSEXP_2_string(SEXP Ss, int i) {
  if(!isString(Ss)) {
    PRINTF("Error: STRSEXP_2_string called for SEXP that is not a string!\n"); 
    return(string(""));
  }
  if(LENGTH(Ss) <= i) {
    PRINTF("Error: STRSEXP_2_string called for (C) element %i of an SEXP that has length %i!\n", i, LENGTH(Ss));
    return(string(""));
  }
  int l = LENGTH(STRING_ELT(Ss, i));
  string ans(CHAR(STRING_ELT(Ss,i)),l);
  return(ans);
}

void STRSEXP_2_vectorString(SEXP Ss, vector<string> &ans) {
  if(!isString(Ss)) {
    PRINTF("Error: STRSEXP_2_vectorString called for SEXP that is not a string!\n"); 
    return;
  }
  int nn = LENGTH(Ss);
  ans.resize(nn);
  for(int i = 0; i < nn; i++) {
    ans[i].assign(CHAR(STRING_ELT(Ss, i)), LENGTH(STRING_ELT(Ss, i)));
  }
}

SEXP string_2_STRSEXP(string v) {
  SEXP Sans;
  PROTECT(Sans = allocVector(STRSXP, 1));
  SET_STRING_ELT(Sans, 0, mkChar(v.c_str()));
  UNPROTECT(1);
  return(Sans);
}

SEXP vectorString_2_STRSEXP(const vector<string> &v) {
  SEXP Sans;
  int nn = v.size();
  PROTECT(Sans = allocVector(STRSXP, nn));
  for(int i = 0; i < nn; i++) {
    SET_STRING_ELT(Sans, i, mkChar(v[i].c_str()));
  }
  UNPROTECT(1);
  return(Sans);
}

vector<double> SEXP_2_vectorDouble( SEXP Sn ) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_vectorDouble called for SEXP that is not a numeric or logica!\n");
  int nn = LENGTH(Sn);
  vector<double> ans(nn);
  if(isReal(Sn)) {
    double *dSn = REAL(Sn);
    copy(dSn, dSn + nn, ans.begin());
  } else {
    if(isInteger(Sn) || isLogical(Sn)) {
      int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
      for(int i = 0; i < nn; ++i) {
	ans[i] = static_cast<double>(iSn[i]);
      }
    } else {
      PRINTF("Error: We could not handle the input type to SEXP_2_vectorDouble\n");
    }
  }
  return(ans);
}

double SEXP_2_double(SEXP Sn, int i) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_double called for SEXP that is not numeric or logical\n");
  if(LENGTH(Sn) <= i) PRINTF("Error: SEXP_2_double called for element %i >= length of %i.\n", i, LENGTH(Sn));
  if(isReal(Sn)) {
    return(REAL(Sn)[i]);
  } 
  if(isInteger(Sn) || isLogical(Sn)) {
    if(isInteger(Sn))
      return(static_cast<double>(INTEGER(Sn)[i]));
    else
      return(static_cast<double>(LOGICAL(Sn)[i]));
  }
  PRINTF("Error: We could not handle the input type to SEXP_2_double\n");
  return(0.);
}

SEXP double_2_SEXP(double v) {
  SEXP Sans;
  PROTECT(Sans = allocVector(REALSXP, 1));
  REAL(Sans)[0] = v;
  UNPROTECT(1);
  return(Sans);
}

SEXP vectorDouble_2_SEXP(const vector<double> &v) {
  SEXP Sans;
  int nn = v.size();
  PROTECT(Sans = allocVector(REALSXP, nn));
  if(nn > 0) {
    copy(v.begin(), v.end(), REAL(Sans));
  }
  UNPROTECT(1);
  return(Sans);
}

SEXP vectorInt_2_SEXP(const vector<int> &v) {
  SEXP Sans;
  int nn = v.size();
  PROTECT(Sans = allocVector(INTSXP, nn));
  if(nn > 0) {
    copy(v.begin(), v.end(), INTEGER(Sans));
  }
  UNPROTECT(1);
  return(Sans);
}


void vectorDouble_2_SEXP(const vector<double> &v, SEXP Sn) {
  int nn = v.size();
  PROTECT(Sn = allocVector(REALSXP, nn));
  if(nn > 0) {
    copy(v.begin(), v.end(), REAL(Sn));
  }
  UNPROTECT(1);
}

void vectorInt_2_SEXP(const vector<int> &v, SEXP Sn) {
  int nn = v.size();
  PROTECT(Sn = allocVector(INTSXP, nn));
  if(nn > 0) {
    copy(v.begin(), v.end(), INTEGER(Sn));
  }
  UNPROTECT(1);
}


struct opIntegerShift {
public:
  int c;
  opIntegerShift(int inputc) : c(inputc) {};
  int operator()(int x) {return(x + c);}
};

SEXP vectorInt_2_SEXP(const vector<int> &v, int offset) {
  SEXP Sans;
  int nn = v.size();
  PROTECT(Sans = allocVector(INTSXP, nn));
  if(nn > 0) {
    if(offset == 0)
      copy(v.begin(), v.end(), INTEGER(Sans));
    else
      std::transform(v.begin(), v.end(), INTEGER(Sans), opIntegerShift(offset));
  }
  UNPROTECT(1);
  return(Sans);
}



vector<int> SEXP_2_vectorInt( SEXP Sn, int offset ) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_vectorInt called for SEXP that is not a numeric or logica!\n");
  int nn = LENGTH(Sn);
  vector<int> ans(nn);
  if(isInteger(Sn) || isLogical(Sn)) {
    int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    if(offset == 0) copy(iSn, iSn + nn, ans.begin());
    else {
      std::transform(iSn, iSn + nn, ans.begin(), opIntegerShift(offset)); 
    }
  } else {
    if(isReal(Sn)) {
      double *dSn = REAL(Sn);
      bool warning = false;
      for(int i = 0; i < nn; ++i) {
	if(dSn[i] != floor(dSn[i])) warning = true;
	ans[i] = static_cast<int>(dSn[i] + offset);
      }
      if(warning) PRINTF("Warning from SEXP_2_vectorInt: some input elements are reals that do not have integer values\n");
    } else {
      PRINTF("Error: We could not handle input type to SEXP_2_vectorInt\n");
    }
  }
  return(ans);
}

int SEXP_2_int(SEXP Sn, int i, int offset ) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_int called for SEXP that is not numeric or logical\n");
  if(LENGTH(Sn) <= i) PRINTF("Error: SEXP_2_int called for element %i% >= length of %i.\n", i, LENGTH(Sn));
  if(isInteger(Sn) || isLogical(Sn)) {
    if(isInteger(Sn))
      return(INTEGER(Sn)[i]);
    else
      return(LOGICAL(Sn)[i]);
  } else {
    if(isReal(Sn)) {
      double ans = REAL(Sn)[i] + offset;
      if(ans != floor(ans)) PRINTF("Warning from SEXP_2_int: input element is a real with a non-integer value\n");
      return(static_cast<int>(ans));
    } else {
      PRINTF("Error: We could not handle input type to  SEXP_2_int\n");
    }
  }
  return(0);
}

SEXP int_2_SEXP(int i) {
  SEXP Sans;
  PROTECT(Sans = allocVector(INTSXP, 1));
  INTEGER(Sans)[0] = i;
  UNPROTECT(1);
  return(Sans);
}

SEXP bool_2_SEXP(bool ind){
	SEXP Sans;
	PROTECT(Sans = allocVector(LGLSXP, 1) );
	LOGICAL(Sans)[0] = ind;
	UNPROTECT(1);
	return(Sans);
}

bool SEXP_2_bool(SEXP Sn, int i) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_bool called for SEXP that is not numeric or logical\n");
  if(LENGTH(Sn) <= i) PRINTF("Error: SEXP_2_bool called for element %i% >= length of %i.\n", i, LENGTH(Sn));
  if(isLogical(Sn)) {
    return(static_cast<bool>(LOGICAL(Sn)[i]));
  } else {
    if(isInteger(Sn)) {
      return(static_cast<bool>(INTEGER(Sn)[i]));
    } else {
      if(isReal(Sn)) {
	return(static_cast<bool>(REAL(Sn)[i]));
      }
    }
  }
  PRINTF("Error: could not handle input type to SEXP_2_bool\n");
  return(false);
}

// bool checkString(SEXP Ss, int len) {
//   if(!isString(Ss)) {
//     PRINTF("Error: something that was supposed to be a string is not\n"); 
//     return(false);
//   }
//   if(LENGTH(Ss) < len) {
//     PRINTF("Error: A string vector that was supposed to have length at least %i does not\n", len);
//     return(false);
//   }
//   return(true);
// }

// bool checkNumeric(SEXP Sval, int len) {
//   if(!isNumeric(Sval)) return(false);
//   if(LENGTH(Sval) != len) return(false);
//   return(true);
// }

// SEXP setVec(SEXP Sextptr, SEXP Svalue) {
//   if(!isNumeric(Svalue)) {
//     PRINTF("Error: Svalue is not a numeric!\n");
//     return(R_NilValue);
//   }
//   if(!R_ExternalPtrAddr(Sextptr)) {
//      RBREAK("Error: Sextptr is not a valid external pointer\n");
//   }
//   //  vector<double> *vecPtr = *static_cast< vector<double>** >(R_ExternalPtrAddr(Sextptr));
//   NimArrBase<double> *vecPtr = *static_cast< NimArrBase<double>** >(R_ExternalPtrAddr(Sextptr));

//  int nn = LENGTH(Svalue);
//   int lengthCpp = vecPtr->size();
//   if(nn != lengthCpp) {
//     PRINTF("Error: length of Svalues does not match C++ vector length.  Using smaller.\n");
//     if(nn > lengthCpp) nn = lengthCpp;
//   }
//   double *value = REAL(Svalue);
//   std::copy(value, value + nn, &(vecPtr->v[0]) );	
//   return(R_NilValue);
// }
  
// SEXP getVec(SEXP Sextptr) {
//   if(!R_ExternalPtrAddr(Sextptr)) {
//     PRINTF("Error: Sextptr is not a valid external pointer\n");
//     return(R_NilValue);
//   }

//   NimArrBase<double> *vecPtr = static_cast< NimArrBase<double> * >(R_ExternalPtrAddr(Sextptr));

//   int len = vecPtr->size();
//   SEXP Sans;  
  
//   PROTECT(Sans = allocVector(REALSXP, len));
//   //std::copy(vecPtr->v.begin(), vecPtr->v.end() , REAL(Sans) );
//   std::copy(vecPtr->v, vecPtr->v + len , REAL(Sans) );
  
//   int numDims = vecPtr->numDims();
//   if(numDims > 1) {
//     SEXP Sdim;
//     PROTECT(Sdim = allocVector(INTSXP, numDims));
//     for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = vecPtr->dimSize(idim);
//     setAttrib(Sans, R_DimSymbol, Sdim);
//     UNPROTECT(2);
//   } else {
//     UNPROTECT(1);
//   } 
//  return(Sans);
// }

// SEXP getVec_Integer(SEXP Sextptr) {
//   if(!R_ExternalPtrAddr(Sextptr)) {
//     PRINTF("Error: Sextptr is not a valid external pointer\n");
//     return(R_NilValue);
//   }
//   NimArrBase<int> *vecPtr = *static_cast< NimArrBase<int> ** >(R_ExternalPtrAddr(Sextptr));
//   SEXP Sans;
//   int len = vecPtr->size();
  
//   PROTECT(Sans = allocVector(REALSXP, len));
//   //  std::copy(vecPtr->v.begin(), vecPtr->v.end(), INTEGER(Sans) );
//   std::copy(vecPtr->v, vecPtr->v + len, INTEGER(Sans) );  
  
//   int numDims = vecPtr->numDims();
//   if(numDims > 1) {
//     SEXP Sdim;
//     PROTECT(Sdim = allocVector(INTSXP, numDims));
//     for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = vecPtr->dimSize(idim);
//     setAttrib(Sans, R_DimSymbol, Sdim);
//     UNPROTECT(2);
//   } else {
//     UNPROTECT(1);
//   } 
//  return(Sans);
// }

// void dontDeleteFinalizer(SEXP ptr){
// 	return;
// }

SEXP SEXP_2_double(SEXP rPtr, SEXP refNum, SEXP rScalar){
    void* vPtr = R_ExternalPtrAddr(rPtr);
    if(vPtr == NULL){
        PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
        return(R_NilValue);
    }
    double* cPtr;
    int cRefNum = INTEGER(refNum)[0];
    if(cRefNum == 1)
        cPtr = static_cast<double*>( vPtr );
    else if(cRefNum == 2)
        cPtr = (*static_cast<double**> ( vPtr ) );
    else return(R_NilValue);
    
    if(isLogical(rScalar) )
      (*cPtr) = LOGICAL(rScalar)[0];
    else if(isInteger(rScalar) )
      (*cPtr) = INTEGER(rScalar)[0];
    else if(isReal(rScalar) )
      (*cPtr) = REAL(rScalar)[0];
    else
        PRINTF("R class not identified. Currently numeric and integers supported\n");
    return(R_NilValue);
}

SEXP SEXP_2_int(SEXP rPtr, SEXP refNum, SEXP rScalar){
    void* vPtr = R_ExternalPtrAddr(rPtr);
    if(vPtr == NULL){
        PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
        return(R_NilValue);
    }
    int* cPtr;
    int cRefNum = INTEGER(refNum)[0];
    if(cRefNum == 1)
        cPtr = static_cast<int*>( vPtr );
    else if(cRefNum == 2)
        cPtr = (*static_cast<int**> ( vPtr ) );
    else return(R_NilValue);
    if(isInteger(rScalar) )
        (*cPtr) = INTEGER(rScalar)[0];
    else if(isReal(rScalar) )
        (*cPtr) = REAL(rScalar)[0];
    else
        PRINTF("R class not identified. Currently numeric and integers supported\n");
    return(R_NilValue);
}		//Note: this is identical to above and one should be removed (ie. just SEXP_2_scalar)
		//But our generated code calls both


  SEXP SEXP_2_bool(SEXP rPtr, SEXP refNum, SEXP rScalar){
    void* vPtr = R_ExternalPtrAddr(rPtr);
    if(vPtr == NULL){
        PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
        return(R_NilValue);
    }
    bool* cPtr;
    int cRefNum = INTEGER(refNum)[0];
    if(cRefNum == 1)
        cPtr = static_cast<bool*>( vPtr );
    else if(cRefNum == 2)
        cPtr = (*static_cast<bool**> ( vPtr ) );
    else return(R_NilValue);
    if(isLogical(rScalar) )
      (*cPtr) = LOGICAL(rScalar)[0];
    else if(isInteger(rScalar) )
        (*cPtr) = INTEGER(rScalar)[0];
    else if(isReal(rScalar) )
        (*cPtr) = REAL(rScalar)[0];
    else
        PRINTF("R class not identified. Currently numeric and integers supported\n");
    return(R_NilValue);
}

SEXP bool_2_SEXP(SEXP rPtr, SEXP refNum){
    void* vPtr = R_ExternalPtrAddr(rPtr);
    if(vPtr == NULL){
        PRINTF("Warning: pointing to NULL in bool_2_SEXP\n");
        return(R_NilValue);
    }
    bool* cPtr;
    int cRefNum = INTEGER(refNum)[0];
    if(cRefNum == 1)
        cPtr = static_cast<bool*>( vPtr );
    else if(cRefNum == 2)
        cPtr = (*static_cast<bool**> ( vPtr ) );
    else {
      PRINTF("Warning: bool_2_SEXP called with reNum != 1 or 2\n");
      return(R_NilValue);
    }
    SEXP Sans;
    PROTECT(Sans = allocVector(LGLSXP, 1));
    LOGICAL(Sans)[0] = (*cPtr);
    UNPROTECT(1);
    return(Sans);
}

SEXP double_2_SEXP(SEXP rPtr, SEXP refNum){
    void* vPtr = R_ExternalPtrAddr(rPtr);
    if(vPtr == NULL){
        PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
        return(R_NilValue);
    }
    double* cPtr;
    int cRefNum = INTEGER(refNum)[0];
    if(cRefNum == 1)
        cPtr = static_cast<double*>( vPtr );
    else if(cRefNum == 2)
        cPtr = (*static_cast<double**> ( vPtr ) );
    else {
      PRINTF("Warning: double_2_SEXP called with reNum != 1 or 2\n");
      return(R_NilValue);
    }
    SEXP Sans;
    PROTECT(Sans = allocVector(REALSXP, 1));
    REAL(Sans)[0] = (*cPtr);    
    UNPROTECT(1);
    return(Sans);
}


SEXP int_2_SEXP(SEXP rPtr, SEXP refNum){
    void* vPtr = R_ExternalPtrAddr(rPtr);
    if(vPtr == NULL){
        PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
        return(R_NilValue);
    }
    int* cPtr;
    int cRefNum = INTEGER(refNum)[0];
    if(cRefNum == 1)
        cPtr = static_cast<int*>( vPtr );
    else if(cRefNum == 2)
        cPtr = (*static_cast<int**> ( vPtr ) );
    else {
        PROBLEM "incorrect value passed to int_2_SEXP"
            ERROR;
    }
    SEXP Sans;
    PROTECT(Sans = allocVector(INTSXP, 1));
    INTEGER(Sans)[0] = (*cPtr);
    UNPROTECT(1);
    return(Sans);
}

SEXP SEXP_2_string(SEXP rPtr, SEXP rString) {
  void* vPtr = R_ExternalPtrAddr(rPtr);
  if(vPtr == NULL){
    PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
    return(R_NilValue);
  }
  *static_cast<string*>(vPtr) = STRSEXP_2_string(rString, 0);
  return(R_NilValue);
}

SEXP SEXP_2_stringVector(SEXP rPtr, SEXP rStringVector) {
  void* vPtr = R_ExternalPtrAddr(rPtr);
  if(vPtr == NULL){
    PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
    return(R_NilValue);
  }
  STRSEXP_2_vectorString(rStringVector, *static_cast<vector<string> *>(vPtr));
  return(R_NilValue);
}

SEXP string_2_SEXP(SEXP rPtr) {
  void* vPtr = R_ExternalPtrAddr(rPtr);
  if(vPtr == NULL){
    PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
    return(R_NilValue);
  }
  return(string_2_STRSEXP(*static_cast<string *>(vPtr)));
}

SEXP stringVector_2_SEXP(SEXP rPtr) {
  void* vPtr = R_ExternalPtrAddr(rPtr);
  if(vPtr == NULL){
    PRINTF("Warning: pointing to NULL in SEXP_2_double\n");
    return(R_NilValue);
  }
  return(vectorString_2_STRSEXP(*static_cast<vector<string> *>(vPtr)));
}

//Next 3 functions are used for CmodelValues objects
SEXP fastMatrixInsert(SEXP matrixInto, SEXP matrix, SEXP rowStart, SEXP colStart){
	SEXP RdimInto = getAttrib(matrixInto, R_DimSymbol);
	PROTECT(RdimInto);
	int Iinto = INTEGER(RdimInto)[0];
	int Jinto = INTEGER(RdimInto)[1];
	
	SEXP Rdim = getAttrib(matrix, R_DimSymbol);
	PROTECT(Rdim);
	int I = INTEGER(Rdim)[0];
	int J = INTEGER(Rdim)[1];
	
	int cRowStart = INTEGER(rowStart)[0] - 1;
	int cColStart = INTEGER(colStart)[0] - 1;
	
	if((I + cRowStart > Iinto) | (J + cColStart > Jinto)){
		UNPROTECT(2);
		PRINTF("Matrix copying not allowed for given indices\n");
		return(R_NilValue);
	}
	for(int i = 0; i < I; i++){
		for(int j = 0; j < J; j++)
			REAL(matrixInto)[ (i + cRowStart) + Iinto * (j + cColStart) ] = REAL(matrix)[i + I * j];
	}
	UNPROTECT(2);
	return(R_NilValue);
}

SEXP matrix2ListDouble(SEXP matrix, SEXP list, SEXP listStartIndex, SEXP RnRows,  SEXP dims){
  //int cStart = INTEGER(listStartIndex)[0] - 1;
	int cNRows = INTEGER(RnRows)[0];
	int len = 1;
	for(int i = 0; i < LENGTH(dims); i++)
		len = len * INTEGER(dims)[i];
	for(int i = 0; i < cNRows; i++){
	SEXP row = PROTECT(allocVector(REALSXP, len) ) ;
	setAttrib(row, R_DimSymbol, dims);
		for(int j = 0; j <  len; j++){		
			REAL(row)[j] = REAL(matrix)[i + cNRows * j];
			}
		SET_VECTOR_ELT(list, i, row);
	UNPROTECT(1);
	}	
	return(R_NilValue);
}

SEXP matrix2ListInt(SEXP matrix, SEXP list, SEXP listStartIndex, SEXP RnRows,  SEXP dims){
  //	int cStart = INTEGER(listStartIndex)[0] - 1;
	int cNRows = INTEGER(RnRows)[0];
	int len = 1;
	for(int i = 0; i < LENGTH(dims); i++)
		len = len * INTEGER(dims)[i];
	for(int i = 0; i < cNRows; i++){
	SEXP row = PROTECT(allocVector(INTSXP, len) ) ;
	setAttrib(row, R_DimSymbol, dims);
		for(int j = 0; j <  len; j++){		
			INTEGER(row)[j] = INTEGER(matrix)[i + cNRows * j];
			}
		SET_VECTOR_ELT(list, i, row);
	UNPROTECT(1);
	}	
	return(R_NilValue);
}



void rawSample(double* p, int c_samps, int N, int* ans, bool unsort, bool silent) {
  vector<double> cdf(N+1);
  cdf[0] = 0;
  bool badVals = false;
  for(int i = 1; i < N+1; i++ ){
    cdf[i] = cdf[i-1] + p[i-1];
    if(!(p[i-1] >= 0)){
      badVals = true;
      if(!silent)
	PRINTF("Warning: negative probability given to rankSample. Returning a uniform sample.\n");
      cdf[N] = 1;
      break;
    }
  }
  double sum = cdf[N];
  if(sum == 0){
    badVals = true;
    if(!silent)
      PRINTF("Warning: sum of weights = 0 in rankSample. Returning a uniform sample.\n");
  }
  if(badVals){
    for(int i = 1; i <= c_samps; i++)
      ans[i-1] = ((i-1) % N) + 1;  // generates a cyclic uniform sample (DT, May 2015)
    return;
  }
  cdf[N] = sum + 1;
  vector<double> sampP(c_samps + 1);
  sampP[0] = 1 - exp( log( unif_rand() ) / c_samps );
	
  sampP[0] = sampP[0] * sum;
  sampP[c_samps] = sum + 1;
	
  for(int i = 1; i < c_samps ; i++)
    sampP[i] = (1 - exp(  log(unif_rand()) / (c_samps - i)   ) )* (sum - sampP[i-1]) + sampP[i-1];
  int curP = 0;
  if(unsort == false){
    for(int i = 1; i <= N; i++){
      while(cdf[i] > (sampP[curP])){
	ans[curP] = i;
	curP++;
      }
    }
    return;	
  }
  // unsort must be true to get here
  vector<double> sortAns(c_samps);	
  for(int i = 1; i <= N; i++){
    while(cdf[i] > (sampP[curP])){
      sortAns[curP] = i;
      curP++;
    }
  }
	
  vector<int> newOrder(c_samps);
  for(int i = 0; i < c_samps;i++)
    newOrder[i] = i;
  int drawIndex;
  for(int i = c_samps-1; i >= 0; i--){
    drawIndex = unif_rand() * i;
    ans[i] = sortAns[newOrder[drawIndex]];
    newOrder[drawIndex] = newOrder[i];
  }
}

SEXP rankSample(SEXP p, SEXP n, SEXP not_used, SEXP s) {
  //PRINTF("in SEXP rankSample\n");
  int N = LENGTH(p);
  int c_samps = INTEGER(n)[0];
  bool silent = LOGICAL(s)[0];
  SEXP output;
  PROTECT(output = allocVector(INTSXP, c_samps ));
  GetRNGstate();
  rawSample(REAL(p), c_samps, N, INTEGER(output), false, silent);
  PutRNGstate();
  UNPROTECT(1);
  return(output);
}

void parseVar(const vector<string> &input, vector<string> &output) {
  int vecSize = input.size();
  std::size_t iEnd, iBegin;
  output.resize( vecSize );
  for(int i = 0; i < vecSize; i++) {
    iBegin = input[i].find_first_not_of(_NIMBLE_WHITESPACE_UTIL);
    iEnd = input[i].find_first_of(_NIMBLE_WHITESPACEBRACKET_UTIL, iBegin);
    if(iBegin < iEnd)
      output[i].assign( input[i].substr(iBegin, iEnd - iBegin) );
    else
      output[i].assign( string("") );
    //    if(iBracket != std::string::npos)

    //    else
    //   output[i].assign( input[i] );
  }
}

SEXP parseVar(SEXP Sinput) {
  vector<string> input, output;
  STRSEXP_2_vectorString(Sinput, input);
  parseVar(input, output);
  return(vectorString_2_STRSEXP(output));
}


SEXP makeNewNimbleList(SEXP S_listName) {
  SEXP SnimbleInternalFunctionsEnv;
  SEXP call;
  PROTECT(SnimbleInternalFunctionsEnv = EVAL(findVar(install("nimbleInternalFunctions"), R_GlobalEnv)));
  PROTECT(call = allocVector(LANGSXP, 2));
  SETCAR(call, install("makeNewNimListSEXPRESSIONFromC"));
  SETCADR(call, S_listName);
  UNPROTECT(2);
  return(eval(call, SnimbleInternalFunctionsEnv));
}