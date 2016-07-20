#include "nimble/NimArrBase.h"
#include "nimble/NamedObjects.h"
#include "nimble/RcppUtils.h"
#include "nimble/Utils.h"
#include<iostream>
#include<sstream>
#include<algorithm>

std::ostringstream _nimble_global_output;

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

SEXP setPtrVectorOfPtrs(SEXP SaccessorPtr, SEXP ScontentsPtr, SEXP Ssize) {
  vectorOfPtrsAccessBase *accessorPtr = static_cast<vectorOfPtrsAccessBase *>(R_ExternalPtrAddr(SaccessorPtr)); 
  void *contentsPtr = static_cast<void *>(R_ExternalPtrAddr(ScontentsPtr));
  int size = INTEGER(Ssize)[0];
  accessorPtr->setTheVec(contentsPtr, size);
  return(R_NilValue);
}

SEXP setOnePtrVectorOfPtrs(SEXP SaccessorPtr, SEXP Si, SEXP ScontentsPtr) {
  vectorOfPtrsAccessBase *accessorPtr = static_cast<vectorOfPtrsAccessBase *>(R_ExternalPtrAddr(SaccessorPtr)); 
  void *contentsPtr = static_cast<void *>(R_ExternalPtrAddr(ScontentsPtr));
  int i = INTEGER(Si)[0];
  accessorPtr->setVecPtr(i, contentsPtr);
  return(R_NilValue);
}

SEXP getOnePtrVectorOfPtrs(SEXP SaccessorPtr, SEXP Si) {
  vectorOfPtrsAccessBase *accessorPtr = static_cast<vectorOfPtrsAccessBase *>(R_ExternalPtrAddr(SaccessorPtr)); 
  int i = INTEGER(Si)[0];
  void *ans = accessorPtr->getVecPtr(i);
  return( R_MakeExternalPtr(ans, R_NilValue, R_NilValue));
}

/*This function was created so we can report the status of different operations
from C++ to R. I found that Rbreak does not stop the R function from which the C++
function was called. This can cause an excessive number (e.g. nrow) of printed 
error reports if the user makes a common mistake. Instead, this function is used to 
return the status to R. If it is false, we use the R function "stop" which will break
for loops, etc. */
 
SEXP returnStatus(bool OK)
	{
	SEXP ans;
    PROTECT(ans = allocVector(LGLSXP, 1) );
    LOGICAL(ans)[0] = OK;
    UNPROTECT(1);
    return(ans);
	}

bool all(bool x) {return(x);}
bool any(bool x) {return(x);}
int t(int x) {return(x);}
double t(double x) {return(x);}
int prod(int x) {return(x);}
double prod(double x) {return(x);}

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
};

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
};


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

vector<int> getSEXPdims(SEXP Sx) {
  if(!isNumeric(Sx)) {PRINTF("Error, getSEXPdims called for something not numeric\n"); return(vector<int>());}
  if(!isVector(Sx)) {PRINTF("Error, getSEXPdims called for something not vector\n"); return(vector<int>());}
  if(!isArray(Sx) & !isMatrix(Sx)) {
    vector<int> ans; 
    ans.resize(1); ans[0] = LENGTH(Sx); return(ans);
  }
  return(SEXP_2_vectorInt(getAttrib(Sx, R_DimSymbol), 0));
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
};

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

bool checkString(SEXP Ss, int len) {
  if(!isString(Ss)) {
    PRINTF("Error: something that was supposed to be a string is not\n"); 
    return(false);
  }
  if(LENGTH(Ss) < len) {
    PRINTF("Error: A string vector that was supposed to have length at least %i does not\n", len);
    return(false);
  }
  return(true);
}

bool checkNumeric(SEXP Sval, int len) {
  if(!isNumeric(Sval)) return(false);
  if(LENGTH(Sval) != len) return(false);
  return(true);
}

SEXP setVec(SEXP Sextptr, SEXP Svalue) {
  if(!isNumeric(Svalue)) {
    PRINTF("Error: Svalue is not a numeric!\n");
    return(R_NilValue);
  }
  if(!R_ExternalPtrAddr(Sextptr)) {
     RBREAK("Error: Sextptr is not a valid external pointer\n");
  }
  //  vector<double> *vecPtr = *static_cast< vector<double>** >(R_ExternalPtrAddr(Sextptr));
  NimArrBase<double> *vecPtr = *static_cast< NimArrBase<double>** >(R_ExternalPtrAddr(Sextptr));

 int nn = LENGTH(Svalue);
  int lengthCpp = vecPtr->size();
  if(nn != lengthCpp) {
    PRINTF("Error: length of Svalues does not match C++ vector length.  Using smaller.\n");
    if(nn > lengthCpp) nn = lengthCpp;
  }
  double *value = REAL(Svalue);
  std::copy(value, value + nn, &(vecPtr->v[0]) );	
  return(R_NilValue);
}
  
SEXP getVec(SEXP Sextptr) {
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a valid external pointer\n");
    return(R_NilValue);
  }

  NimArrBase<double> *vecPtr = static_cast< NimArrBase<double> * >(R_ExternalPtrAddr(Sextptr));

  int len = vecPtr->size();
  SEXP Sans;  
  
  PROTECT(Sans = allocVector(REALSXP, len));
  //std::copy(vecPtr->v.begin(), vecPtr->v.end() , REAL(Sans) );
  std::copy(vecPtr->v, vecPtr->v + len , REAL(Sans) );
  
  int numDims = vecPtr->numDims();
  if(numDims > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(INTSXP, numDims));
    for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = vecPtr->dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  } 
 return(Sans);
}

SEXP getVec_Integer(SEXP Sextptr) {
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a valid external pointer\n");
    return(R_NilValue);
  }
  NimArrBase<int> *vecPtr = *static_cast< NimArrBase<int> ** >(R_ExternalPtrAddr(Sextptr));
  SEXP Sans;
  int len = vecPtr->size();
  
  PROTECT(Sans = allocVector(REALSXP, len));
  //  std::copy(vecPtr->v.begin(), vecPtr->v.end(), INTEGER(Sans) );
  std::copy(vecPtr->v, vecPtr->v + len, INTEGER(Sans) );  
  
  int numDims = vecPtr->numDims();
  if(numDims > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(INTSXP, numDims));
    for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = vecPtr->dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  } 
 return(Sans);
}


/* Cliff's new function for adding blank rows to C model values */
SEXP addBlankModelValueRows(SEXP Sextptr, SEXP numAdded){
    if(!isInteger(numAdded)) {
        PRINTF("Error: numAdded is not an integer!\n");
        return(returnStatus(false) );
    }
    
    if(!R_ExternalPtrAddr(Sextptr)) {
        PRINTF("Error: Sextptr is not a valid external pointer\n");
        return(returnStatus(false) );
    }
    NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
    nimType vecType = (*typePtr).getNimType();
    if(vecType == DOUBLE){
    	VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
	    int new_size = INTEGER(numAdded)[0] + nrowCpp;
    	matPtr->resize(new_size);
    	NimArrBase<double> *thisRow;
    	thisRow = matPtr->getBasePtr(0);
    	int numDims = thisRow->numDims();
	     vector<int> Dims(numDims);
 		for(int i = 0; i < numDims; i++)
      		Dims[i] = thisRow->dimSize(i);
    	for(int i = nrowCpp; i < INTEGER(numAdded)[0] + nrowCpp; i++){
      		thisRow = matPtr->getBasePtr(i); 
      		thisRow->setSize(Dims);
    	}
    return(returnStatus(true) );
    }
    
    else if(vecType == INT){
    	VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
    	  int new_size = INTEGER(numAdded)[0] + nrowCpp;
	  	matPtr->resize(new_size);
    	NimArrBase<int> *thisRow;
    	thisRow = matPtr->getBasePtr(0);
    	int numDims = thisRow->numDims();
  	    vector<int> Dims(numDims);
		for(int i = 0; i < numDims; i++)
      		Dims[i] = thisRow->dimSize(i);
    	for(int i = nrowCpp; i < INTEGER(numAdded)[0] + nrowCpp; i++){
      		thisRow = matPtr->getBasePtr(i); 
      		thisRow->setSize(Dims);
    	}
    return(returnStatus(true) );
    }
 	
 	PRINTF("Data type for VecNimArr not currently supported\n");
 	return(returnStatus(false) ) ;
}

SEXP setNumListRows(SEXP Sextptr, SEXP nRows, SEXP setSize2row1){
    NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
    nimType vecType = (*typePtr).getNimType();
    if(vecType == DOUBLE){
    	VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
	    int new_size = INTEGER(nRows)[0];
    	matPtr->resize(new_size);
    	if(new_size <= nrowCpp)
    		return(returnStatus(true));
    	if(LOGICAL(setSize2row1)[0]	== TRUE){
		    	NimArrBase<double> *thisRow;
    			thisRow = matPtr->getBasePtr(0);
    			int numDims = thisRow->numDims();
		     	vector<int> Dims(numDims);
 				for(int i = 0; i < numDims; i++)
    	  			Dims[i] = thisRow->dimSize(i);
   		   		for(int i = nrowCpp; i < new_size; i++){
   		   			thisRow = matPtr->getBasePtr(i);
    	  			thisRow->setSize(Dims);
    	  		}
    	return(returnStatus(true) );
    	}
    }
    else if(vecType == INT){
    	VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
    	  int new_size = INTEGER(nRows)[0];
	  	matPtr->resize(new_size);
	  	if(new_size <= nrowCpp)
    		return(returnStatus(true));

    	if(LOGICAL(setSize2row1)[0]	== TRUE){
    		NimArrBase<int> *thisRow;
    		thisRow = matPtr->getBasePtr(0);
    		int numDims = thisRow->numDims();
  	    	vector<int> Dims(numDims);
			for(int i = 0; i < numDims; i++)
      			Dims[i] = thisRow->dimSize(i);
      		for(int i = nrowCpp; i < new_size; i++){
      			thisRow = matPtr->getBasePtr(i);
      			thisRow->setSize(Dims);
      		}
      	}
    return(returnStatus(true) );
    }
 	
 	PRINTF("Data type for VecNimArr not currently supported\n");
 	return(returnStatus(false) ) ;
}





/* Cliff's new function for retreiving nrow from CmodelValues object */
SEXP getNRow(SEXP Sextptr){
    if(!R_ExternalPtrAddr(Sextptr)) {
        PRINTF("Error: Sextptr is not a valid external pointer\n");
        return(R_NilValue);
    }
    int nRow = 0;
    NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
    nimType vecType = (*typePtr).getNimType();
    if(vecType == DOUBLE){
    	VecNimArrBase<double>* matPtr = static_cast<VecNimArrBase<double>* > (typePtr);
	    nRow = matPtr->size();
	}
	else if(vecType == INT)	{
    	VecNimArrBase<int>* matPtr = static_cast<VecNimArrBase<int>* > (typePtr);
	    nRow = matPtr->size();
	}
	else{
		PRINTF("Data type of VecNimArr not currently supported\n");
		_nimble_global_output<< "vecType = " << vecType << "\n"; nimble_print_to_R(_nimble_global_output);
//		cout << "vecType DOUBLE = " << DOUBLE << "\n";
		}
	SEXP rNRow;
	PROTECT(rNRow = allocVector(INTSXP, 1) ); 
	INTEGER(rNRow)[0] = nRow;
	UNPROTECT(1);
    return(rNRow);
}


/* Cliff's new function for copying variables of given rows from one CModelValues object to another */
/* This copies the values from SextptrFrom into SextptrTo. It takes values from rowsFrom and copies
into rowsTo */
SEXP copyModelValuesElements(SEXP SextptrFrom, SEXP SextptrTo, SEXP rowsFrom, SEXP rowsTo){
    if(!R_ExternalPtrAddr(SextptrFrom) | !R_ExternalPtrAddr(SextptrFrom)) {
        PRINTF("Error: Sextptr is not a valid external pointer\n");
        return(returnStatus(false));
    }
    if(LENGTH(rowsFrom) != LENGTH(rowsTo) ){
        PRINTF("Error: length of rowsFrom != rowsTo");
        return(returnStatus(false) );
    }
    if(LENGTH(rowsFrom) == 0){
        return(returnStatus(true) );
    }
  
  NimVecType *typePtrFrom = static_cast< NimVecType* >(R_ExternalPtrAddr(SextptrFrom));
  nimType vecTypeFrom = (*typePtrFrom).getNimType();
  NimVecType *typePtrTo = static_cast< NimVecType* >(R_ExternalPtrAddr(SextptrTo));
  nimType vecTypeTo = (*typePtrTo).getNimType();

  if((vecTypeFrom == DOUBLE) & (vecTypeTo == DOUBLE)){
	  VecNimArrBase<double> *matPtrFrom = static_cast< VecNimArrBase<double>* >(typePtrFrom);
	  VecNimArrBase<double> *matPtrTo = static_cast< VecNimArrBase<double>* >(typePtrTo );
	  int k = LENGTH(rowsFrom);
	  int sizeFrom = matPtrFrom->size();
	  int sizeTo = matPtrTo->size();
	  int *indexFrom = INTEGER(rowsFrom);
	  int *indexTo = INTEGER(rowsTo);
	  NimArrBase<double> *thisRowFrom;
	  NimArrBase<double> *thisRowTo;
	  int ncFrom = 0;
	  int ncTo = 0;
    
	  for(int i = 0; i < k; i++){
	  if((indexFrom[i] > sizeFrom) | (indexFrom[i] <= 0))
	    {
	    _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexFrom[i] << " sizeFrom = " << sizeFrom << "\n";
	    nimble_print_to_R(_nimble_global_output);
	    if(i > 0)
	      PRINTF("Warning: partial copy completed before error discovered!\n");
	    return(returnStatus(false) );
    	}
  	if((indexTo[i] > sizeTo) | (indexTo[i] <=0))
  	  {
	    _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexTo[i] << " sizeTo = "<< sizeTo << "\n";
	    nimble_print_to_R(_nimble_global_output);
  	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	thisRowFrom = matPtrFrom->getBasePtr(indexFrom[i] - 1);  
  	thisRowTo = matPtrTo->getBasePtr(indexTo[i] - 1);  
  
  	ncFrom = thisRowFrom->size();
  	ncTo = thisRowTo->size();
  
  	if(ncFrom != ncTo){
  	  PRINTF("Error: ncFrom != ncTo\n");
  	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	for(int j = 0; j < ncFrom; j++)
  		    (*thisRowTo)[j] = (*thisRowFrom)[j];
  	}
  	return(returnStatus(true) );
  }

  if((vecTypeFrom == INT) & (vecTypeTo == INT)){
	  VecNimArrBase<int> *matPtrFrom = static_cast< VecNimArrBase<int>* >(typePtrFrom);
	  VecNimArrBase<int> *matPtrTo = static_cast< VecNimArrBase<int>* >(typePtrTo );
	  int k = LENGTH(rowsFrom);
	  int sizeFrom = matPtrFrom->size();
	  int sizeTo = matPtrTo->size();
	  int *indexFrom = INTEGER(rowsFrom);
	  int *indexTo = INTEGER(rowsTo);
	  NimArrBase<int> *thisRowFrom;
	  NimArrBase<int> *thisRowTo;
	  int ncFrom = 0;
	  int ncTo = 0;
    
	  for(int i = 0; i < k; i++){
	  if((indexFrom[i] > sizeFrom) | (indexFrom[i] <= 0))
	    {
	    _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexFrom[i] << " sizeFrom = " << sizeFrom << "\n";
	    nimble_print_to_R(_nimble_global_output);
	    if(i > 0)
	      PRINTF("Warning: partial copy completed before error discovered!\n");
	    return(returnStatus(false) );
    	}
  	if((indexTo[i] > sizeTo) | (indexTo[i] <=0))
  	  {
  	  _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexTo[i] << " sizeTo = "<< sizeTo << "\n";
	  nimble_print_to_R(_nimble_global_output);
	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	thisRowFrom = matPtrFrom->getBasePtr(indexFrom[i] - 1);  
  	thisRowTo = matPtrTo->getBasePtr(indexTo[i] - 1);  
  
  	ncFrom = thisRowFrom->size();
  	ncTo = thisRowTo->size();
  
  	if(ncFrom != ncTo){
  	  PRINTF("Error: ncFrom != ncTo\n");
  	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	for(int j = 0; j < ncFrom; j++)
  		    (*thisRowTo)[j] = (*thisRowFrom)[j];
  	}
  	return(returnStatus(true) );
  }

  PRINTF("Data type not currently supported for VecNimArr. Copy MV elements failed\n");
  return(returnStatus(false) );

}	

SEXP getMVElement(SEXP Sextptr, SEXP Sindex){
	if(!isInteger(Sindex)) {
    	PRINTF("Error: Sindex is not an integer!\n");
    	return(returnStatus(false) ) ;
	}
  	if(!R_ExternalPtrAddr(Sextptr)) {
  	  PRINTF("Error: Sextptr is not a valid external pointer\n");
  	  return(returnStatus(false) ) ;
  	}
	NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
	nimType vecType = (*typePtr).getNimType();
	int nrowCpp = typePtr->size();
	int index = INTEGER(Sindex)[0];
	if(index > nrowCpp){
		PRINTF("Error: index too large\n");
		return(returnStatus(false) ) ;
	  }
  	if(index < 1){
		PRINTF("Error: index < 1\n");
		return(returnStatus(false) ) ;
  	}
  return(cGetMVElementOneRow(typePtr, vecType, index) ) ;	
}  	
  	

//SEXP cGetMVElementOneRow(NimVecType* typePtr, nimType vecType, int nrowCpp, int index) 	{  	
//ONe
SEXP cGetMVElementOneRow(NimVecType* typePtr, nimType vecType, int index) 	{  	
    if(vecType == DOUBLE){
		VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(typePtr);
//		int nrowCpp = matPtr->size();
		NimArrBase<double> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting R index to C index
		int outputLength = thisRow->size();
  		SEXP Sans;
  		PROTECT(Sans = allocVector(REALSXP, outputLength));
		double *value = REAL(Sans);
  		std::copy(thisRow->getPtr(), thisRow->getPtr() + outputLength, value); /* copy should work now */
  		if(thisRow->numDims() > 1) { // using first row as representative
    		SEXP Sdim;
    		PROTECT(Sdim = allocVector(INTSXP, thisRow->numDims()));
    		for(int idim = 0; idim < thisRow->numDims(); ++idim) INTEGER(Sdim)[idim] = thisRow->dimSize(idim);
    			setAttrib(Sans, R_DimSymbol, Sdim);
    		UNPROTECT(2);
  		}
  		else {
    UNPROTECT(1);
  		}
  	 	return(Sans);
  	}
    else if(vecType == INT){
		VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(typePtr);
//		int nrowCpp = matPtr->size();
		NimArrBase<int> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting R index to C index
		int outputLength = thisRow->size();
  		SEXP Sans;
  		PROTECT(Sans = allocVector(INTSXP, outputLength));
		int *value = INTEGER(Sans);
  		std::copy(thisRow->getPtr(), thisRow->getPtr() + outputLength, value); /* copy should work now */
  		if(thisRow->numDims() > 1) { // using first row as representative
    		SEXP Sdim;
    		PROTECT(Sdim = allocVector(INTSXP, thisRow->numDims()));
    		for(int idim = 0; idim < thisRow->numDims(); ++idim) INTEGER(Sdim)[idim] = thisRow->dimSize(idim);
    			setAttrib(Sans, R_DimSymbol, Sdim);
    		UNPROTECT(2);
  		}
  		else {
    UNPROTECT(1);
  		}
  	  	return(Sans);	
  	}
  	else
  		PRINTF("VecNimArr datatype not supported\n");
	return(R_NilValue);
}

  SEXP setMVElementFromList(SEXP vNimPtr, SEXP rList, SEXP Sindices){
  	// This function needs to get the length of rList
  	// Once it has this length, it marches down the length of rList and sets all the values of
  	// of vNimPtr
  	int listLength = LENGTH(rList);
  	int checkLength = LENGTH(Sindices);
  	if(listLength != checkLength){
  		PRINTF("Number of indices copying to does not match length of list\n");
  		return(R_NilValue);
  	}
	SEXP listElement;
	NimVecType *typePtr = static_cast<NimVecType*>(R_ExternalPtrAddr(vNimPtr) );
	nimType vecType = (*typePtr).getNimType();
	int index;
	for(int i = 0; i < listLength ; i++){
		listElement = VECTOR_ELT(rList, i);
		index = INTEGER(Sindices)[i];
		cSetMVElementSingle(typePtr, vecType, index, listElement);
	}
  	return(R_NilValue);
  }


 void cSetMVElementSingle(NimVecType* typePtr, nimType vecType,  int index, SEXP Svalue) {
    if(vecType == DOUBLE){
		VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(typePtr);
  		NimArrBase<double> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting from R index to C index
  		int nc = thisRow->size();
  		double *value = REAL(Svalue);
  		int new_nc = LENGTH(Svalue);
  		if(new_nc != nc){
  			PRINTF("Error: size of values assigned not the same! If this is during compiling, could be from improperly sized r modelValues variable (see user manual section on modelValues)\n");
		return;
  		}
  for(int i = 0; i < nc; i++)
    (*thisRow)[i] = value[i];
	return;
	}

	else if(vecType == INT){
		VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(typePtr);
  		NimArrBase<int> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting from R index to C index
  		int nc = thisRow->size();
  		double *value = REAL(Svalue);
  		int new_nc = LENGTH(Svalue);
  		if(new_nc != nc){
  			PRINTF("Error: size of values assigned not the same!\n");
		return;
  		}
  	for(int i = 0; i < nc; i++)
    	(*thisRow)[i] = value[i];
	return;
	}
	PRINTF("VecNimArr data type not currently supported\n");
	return;
}



 SEXP getMVElementAsList(SEXP SextPtr, SEXP Sindices){
 	int nIndice = LENGTH(Sindices);
 	SEXP sList;
 	PROTECT(sList = allocVector(VECSXP, nIndice) );
 	NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(SextPtr));
	nimType vecType = (*typePtr).getNimType();
	for(int i = 0; i < nIndice; i++)
		SET_VECTOR_ELT(sList, i, cGetMVElementOneRow(typePtr, vecType, INTEGER(Sindices)[i]) );
	UNPROTECT(1);
	return(sList);
 }


SEXP resizeNumListRow(SEXP Sextptr, SEXP Sindex, SEXP dims){
	int nDims = LENGTH(dims);
	vector<int> cDims(nDims);
	for(int i = 0; i < nDims ; i++)
		cDims[i] = INTEGER(dims)[i];
	NimVecType* listBasePtr = static_cast<NimVecType* >(R_ExternalPtrAddr(Sextptr) ) ;
	int cRow = INTEGER(Sindex)[0]-1;	
	(*listBasePtr).setRowDims(cRow, cDims);	
	return(R_NilValue);
}

SEXP setMVElement(SEXP Sextptr, SEXP Sindex, SEXP Svalue){
    if(!isInteger(Sindex)) {
    PRINTF("Error: Sindex is not an integer!\n");
	return(returnStatus(false) ) ;
  }
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a valid external pointer\n");
	return(returnStatus(false) ) ;
  }
  int index = INTEGER(Sindex)[0];
  if(index < 1){ 
    PRINTF("Error: index < 1\n");   	
	return(returnStatus(false) ) ;
  }

   NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
   nimType vecType = (*typePtr).getNimType();
  int nrowCpp = typePtr->size();
  if(index > nrowCpp){
	PRINTF("Error: index too large\n");
	return(returnStatus(false) ) ;
  }
  cSetMVElementSingle( typePtr, vecType, index, Svalue );
  SEXP SEXP_2_int(SEXP rPtr, SEXP refNum, SEXP rScalar);
  return(returnStatus(true) );
 } 
 
 
 SEXP getMVsize(SEXP Sextptr){
	NimVecType	*nimVecPtr = static_cast<NimVecType*>(R_ExternalPtrAddr(Sextptr) ); 
 	vector<int> cDims = (*nimVecPtr).getRowDims(0);
	int nDims = cDims.size();
 	SEXP Sans;
 	PROTECT(Sans = allocVector(INTSXP, nDims) );
 	for(int i = 0; i < nDims; i++)
 		INTEGER(Sans)[i] = cDims[i];
 	UNPROTECT(1);
 	return(Sans);
 	}


SEXP setVarPointer(SEXP SextptrModelVar, SEXP SextptrStorageVar, SEXP Srownum) {
  // This function is untested in its new (NimArr) version
  if(!R_ExternalPtrAddr(SextptrModelVar)) {
    PRINTF("Error: Sextptr is not a valid external pointer\n");
    return(R_NilValue);
  }
  if(!R_ExternalPtrAddr(SextptrStorageVar)) {
    PRINTF("Error: Sextptr is not a valid external pointer\n");
    return(R_NilValue);
  }
  if(!isInteger(Srownum)) {
    PRINTF("Error: Srownum is not a valid integer\n");
    return(R_NilValue);
  }
  VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(R_ExternalPtrAddr(SextptrStorageVar));
  NimArrBase<double> **vecPtr = static_cast< NimArrBase<double>** >(R_ExternalPtrAddr(SextptrModelVar));
  (*vecPtr) = matPtr->getBasePtr( INTEGER(Srownum)[0] );
  return(R_NilValue);
}



void dontDeleteFinalizer(SEXP ptr){
	return;
}

SEXP NimArrDouble_2_SEXP(NimArrBase<double> &nimArrDbl){
		int len = nimArrDbl.size();
		SEXP Sans;  
		PROTECT(Sans = allocVector(REALSXP, len));
		//		std::copy(nimArrDbl.v.begin(), nimArrDbl.v.end() , REAL(Sans) );
		std::copy(nimArrDbl.v, nimArrDbl.v + len , REAL(Sans) );  
  		int numDims = nimArrDbl.numDims();
  		if(numDims > 1) {
		  SEXP Sdim;
		  PROTECT(Sdim = allocVector(INTSXP, numDims));
		  for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = nimArrDbl.dimSize(idim);
		  setAttrib(Sans, R_DimSymbol, Sdim);
		  UNPROTECT(2);
  		}
  		else {
    		UNPROTECT(1);
  		} 
	return(Sans);
 }

SEXP NimArrInt_2_SEXP(NimArrBase<int> &nimArrInt){
		int len = nimArrInt.size();
		SEXP Sans;  
		PROTECT(Sans = allocVector(INTSXP, len));
		//		std::copy(nimArrInt.v.begin(), nimArrInt.v.end() , INTEGER(Sans) );
		std::copy(nimArrInt.v, nimArrInt.v + len , INTEGER(Sans) );  
  		int numDims = nimArrInt.numDims();
  		if(numDims > 1) {
    		SEXP Sdim;
    		PROTECT(Sdim = allocVector(INTSXP, numDims));
    		for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = nimArrInt.dimSize(idim);
			setAttrib(Sans, R_DimSymbol, Sdim);
		    UNPROTECT(2);
  		}
  		else {
    		UNPROTECT(1);
  		} 
	return(Sans);
 }


NimArrType* getNimTypePtr(SEXP &rPtr, SEXP &refNum)
{
	NimArrType* nimTypePtr;
	int cRefNum = INTEGER(refNum)[0];
	if(cRefNum == 1){
		nimTypePtr = static_cast<NimArrType*>(R_ExternalPtrAddr(rPtr) );
		return(nimTypePtr); 
		}
	if(cRefNum == 2){
		nimTypePtr = (*static_cast<NimArrType**> (R_ExternalPtrAddr(rPtr) ) );
		return(nimTypePtr);
		}
	PRINTF("Warning: number of dereferences not currently supported, returning null pointer\n");
	nimTypePtr = NULL;
	return(nimTypePtr);
}

void SEXP_2_NimArrDouble (SEXP rValues, NimArrBase<double> &NimArrDbl){
  int rLength = LENGTH(rValues);
  if(rLength != NimArrDbl.size() ) {
    PRINTF("Warning: R object of different size than NimArrDouble. R obj has size %i but NimArrDbl has size %i.\n", rLength, NimArrDbl.size());
    return;
  }
  if(isReal(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrDbl[i] = REAL(rValues)[i];
  }
  else if(isInteger(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrDbl[i] = INTEGER(rValues)[i];
  }
  
  else
    PRINTF("WARNING: class of R object not recognized. Should be numeric or integer\n");
  return;
}

void SEXP_2_NimArrInt (SEXP rValues, NimArrBase<int> &NimArrInt){
  int rLength = LENGTH(rValues);
  if(rLength != NimArrInt.size() ) {
    PRINTF("Warning: R object of different size than NimArrInt!\n");
    return;		
  }
  
  if(isInteger(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrInt[i] = INTEGER(rValues)[i];
  }
  else if(isReal(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrInt[i] = REAL(rValues)[i];
  }
  
  else
    PRINTF("WARNING: class of R object not recognized. Should be numeric or integer\n");    
  return;
}


SEXP Nim_2_SEXP(SEXP rPtr, SEXP NumRefers){
	NimArrType* nimTypePtr = getNimTypePtr(rPtr, NumRefers);
	if(!nimTypePtr)
		return(R_NilValue);
	if(	(*nimTypePtr).getNimType() == INT){
		NimArrBase<int>* nimBase = static_cast<NimArrBase<int> *>(nimTypePtr);
		return(NimArrInt_2_SEXP( (*nimBase) ) ) ;
	}
	if(	(*nimTypePtr).getNimType() == DOUBLE){
		NimArrBase<double>* nimBase = static_cast<NimArrBase<double> *>(nimTypePtr);
		return(NimArrDouble_2_SEXP( (*nimBase) ) );
	}
	PRINTF("Datatype of NimArr not found\n");
	return(R_NilValue);
}


 SEXP SEXP_2_Nim(SEXP rPtr, SEXP NumRefers, SEXP rValues, SEXP allowResize){
    vector<int> sexpDims = getSEXPdims( rValues );
    int sexpNumDims = sexpDims.size();
    bool resize = (LOGICAL(allowResize)[0] == TRUE);

	NimArrType* nimTypePtr = getNimTypePtr(rPtr, NumRefers);
	if(!nimTypePtr)
		return(R_NilValue);
	if(	(*nimTypePtr).getNimType() == INT){
		NimArrBase<int>* nimBase = static_cast<NimArrBase<int> *>(nimTypePtr);
        int nimNumDims = (*nimBase).numDims();
        if(nimNumDims != sexpNumDims){
            if((LENGTH(rValues) != (*nimBase).size()) & (sexpNumDims != 1)){
                PRINTF("Incorrect number of dimensions in copying\n");
                return(R_NilValue);
            }
        }
        if(resize == true)
            (*nimBase).setSize(sexpDims);
		SEXP_2_NimArrInt(rValues,  (*nimBase) ) ;
	}
	if(	(*nimTypePtr).getNimType() == DOUBLE){
		
		NimArrBase<double>* nimBase = static_cast<NimArrBase<double> *>(nimTypePtr);
        int nimNumDims = (*nimBase).numDims();
                
        if(nimNumDims != sexpNumDims){
            if((LENGTH(rValues) != (*nimBase).size()) & (sexpNumDims != 1)){
                PRINTF("Incorrect number of dimensions in copying\n");
                return(R_NilValue);
            }
        }
        if((resize == true) & (nimNumDims == sexpNumDims))
            (*nimBase).setSize(sexpDims);
        SEXP_2_NimArrDouble( rValues, (*nimBase) ) ;
	}
	return(R_NilValue);
}

template <class T>
void cNimArr_2_NimArr(NimArrBase<T> &nimFrom, NimArrBase<T> &nimTo){
	if(nimFrom.size() != nimTo.size() ){
	PRINTF("Warning: NimArr's of different sizes!\n");
	return;
	}
	//	copy(nimFrom.v.begin(), nimFrom.v.end(), nimTo.v.begin() );
	copy(nimFrom.v, nimFrom.v + nimFrom.size(), nimTo.v );
	return;
}

template <class T1, class T2>
void cNimArr_2_NimArr_Index(NimArrBase<T1> &nimFrom, int fromBegin, NimArrBase<T2> &nimTo, int toBegin, int length){
	if(nimFrom.size() != nimTo.size() ){
	PRINTF("Warning: NimArr's of different sizes!\n");
	return;
	}
	//	copy(nimFrom.v.begin() + fromBegin, nimFrom.v.end() + length - 1, nimTo.v.begin() + toBegin); // was this a bug?
	copy(nimFrom.v + fromBegin, nimFrom.v + fromBegin + length, nimTo.v + toBegin);
	return;
}

SEXP Nim_2_Nim(SEXP rPtrFrom, SEXP numRefFrom, SEXP rPtrTo, SEXP numRefTo){
	NimArrType* nimTypePtrFrom = getNimTypePtr(rPtrFrom, numRefFrom);
	NimArrType* nimTypePtrTo = getNimTypePtr(rPtrTo, numRefTo);
	if(!nimTypePtrFrom | !nimTypePtrTo)
		return(R_NilValue);
	if(((*nimTypePtrFrom).getNimType() == INT) & ((*nimTypePtrTo).getNimType() == INT))
		{
		NimArrBase<int>* nimArrFrom = static_cast<NimArrBase<int> * >(nimTypePtrFrom);
		NimArrBase<int>* nimArrTo = static_cast<NimArrBase<int> * >(nimTypePtrTo);
		cNimArr_2_NimArr<int>( (*nimArrFrom), (*nimArrTo) ) ;
		}
	else if( ((*nimTypePtrFrom).getNimType() == DOUBLE) & ((*nimTypePtrTo).getNimType() == DOUBLE ))
		{
		NimArrBase<double>* nimArrFrom = static_cast<NimArrBase<double> * >(nimTypePtrFrom);
		NimArrBase<double>* nimArrTo = static_cast<NimArrBase<double> * >(nimTypePtrTo);
		cNimArr_2_NimArr<double>( (*nimArrFrom), (*nimArrTo) ) ;
		}
	else {
		PRINTF("Copying templates not supported!\n");
		}
	return(R_NilValue);
}


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



/*	Our sample class for which we will build our tools on. 
Later, we want this code to be generated by R	*/

class sampleClass: public NamedObjects {
	public:
	NimArr<1,double> x;
	NimArr<2,int> y;
	NimArr<1,int> ind;
	sampleClass(){
		x.setSize(2);
		y.setSize(2, 2);
		ind.setSize(2); 	
		ind[0] = 0;
		x[0] = 0;
		x[1] = 1;
		namedObjects["x"] = &x;
		namedObjects["y"] = &y;
		namedObjects["ind"] = &ind;
	}
  ~sampleClass(){};
};


SEXP makeNumericList(SEXP nDims, SEXP type, SEXP nRows){
	int cDims = INTEGER(nDims)[0];
	int cRows = INTEGER(nRows)[0];
	string cType = STRSEXP_2_string(type, 0);
	if( (cType == "double" || cType == "numeric") && cDims == 1){
		VecNimArr<1, double>* newList = new VecNimArr<1, double>;
//		(*newList).allowResizeRows = true;
		(*newList).resize(cRows);		
		SEXP ptr = R_MakeExternalPtr(newList, R_NilValue, R_NilValue);
		PROTECT(ptr);
		R_RegisterCFinalizerEx(ptr, &VecNimArrFinalizer<1, double>, TRUE);
		UNPROTECT(1);
		return(ptr);
	}
	if( (cType == "double" || cType == "numeric") && cDims == 2){
		VecNimArr<2, double>* newList = new VecNimArr<2, double>;
//		(*newList).allowResizeRows = true;
		(*newList).resize(cRows);
		SEXP ptr = R_MakeExternalPtr(newList, R_NilValue, R_NilValue);
		PROTECT(ptr);
		R_RegisterCFinalizerEx(ptr, &VecNimArrFinalizer<2, double>, TRUE);
		UNPROTECT(1);
		return(ptr);
	}
	if( (cType == "double" || cType == "numeric") && cDims == 3){
		VecNimArr<3, double>* newList = new VecNimArr<3, double>;
//		(*newList).allowResizeRows = true;
		(*newList).resize(cRows);
		SEXP ptr = R_MakeExternalPtr(newList, R_NilValue, R_NilValue);
		PROTECT(ptr);
		R_RegisterCFinalizerEx(ptr, &VecNimArrFinalizer<3, double>, TRUE);
		UNPROTECT(1);
		return(ptr);
	}

	if( (cType == "integer" || cType == "int") && cDims == 1){
		VecNimArr<1, int>* newList = new VecNimArr<1, int>;
//		(*newList).allowResizeRows = true;
		(*newList).resize(cRows);
		SEXP ptr = R_MakeExternalPtr(newList, R_NilValue, R_NilValue);
		PROTECT(ptr);
		R_RegisterCFinalizerEx(ptr, &VecNimArrFinalizer<1, int>, TRUE);
		UNPROTECT(1);
		return(ptr);
	}
	if( (cType == "double" || cType == "numeric") && cDims == 2){
		VecNimArr<2, int>* newList = new VecNimArr<2, int>;
//		(*newList).allowResizeRows = true;
		(*newList).resize(cRows);
		SEXP ptr = R_MakeExternalPtr(newList, R_NilValue, R_NilValue);
		PROTECT(ptr);
		R_RegisterCFinalizerEx(ptr, &VecNimArrFinalizer<2, int>, TRUE);
		UNPROTECT(1);
		return(ptr);
	}
	if( (cType == "double" || cType == "numeric") && cDims == 3){
		VecNimArr<3, int>* newList = new VecNimArr<3, int>;
//		(*newList).allowResizeRows = true;
		(*newList).resize(cRows);
		SEXP ptr = R_MakeExternalPtr(newList, R_NilValue, R_NilValue);
		PROTECT(ptr);
		R_RegisterCFinalizerEx(ptr, &VecNimArrFinalizer<3, int>, TRUE);
		UNPROTECT(1);
		return(ptr);
	}

        PROBLEM "unhandled case in makeNumericList()"
            ERROR;
}

/*	Generates a sampleClass and returns a void pointer to R
Again, will be generated by R later	*/
	
SEXP newSampObject(){
	sampleClass* newObject = new sampleClass; 
	SEXP ptr = R_MakeExternalPtr(newObject, R_NilValue, R_NilValue);
	PROTECT(ptr);
	R_RegisterCFinalizerEx(ptr, &sampleFinalizer, TRUE);
	UNPROTECT(1);
	return(ptr) ;
}


void sampleFinalizer(SEXP ptr){
	sampleClass* nVPtr = static_cast<sampleClass*>(R_ExternalPtrAddr(ptr));
	delete nVPtr;
}

template<int nDim, class T>
void VecNimArrFinalizer(SEXP ptr){
	VecNimArr<nDim, T>* cPtr = static_cast<VecNimArr<nDim, T>*>(R_ExternalPtrAddr(ptr));
	delete cPtr;
}

double nimMod(double a, double b) {
  return(fmod(a, b));
}

bool compareOrderedPair(orderedPair a, orderedPair b) {  //function called for sort 
  return(a.value < b.value);
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

void rankSample(NimArr<1, double> &weights, int &n, NimArr<1, int> &output) {
  bool silent = false;
  rankSample(weights, n, output, silent);
}

void rankSample(NimArr<1, double> &weights, int &n, NimArr<1, int> &output, bool& silent) {
  //PRINTF("in VOID rankSample\n");
  output.setSize(n);
  int N = weights.size();
  //GetRNGstate();
  rawSample(weights.getPtr(), n, N, output.getPtr(), false, silent);
  //PutRNGstate();
}




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

void row2NimArr(SEXP matrix, NimArrBase<double>* nimPtr, int startPoint, int len, int nRows){
	for(int i = 0; i < len; i++)
		(*nimPtr)[i] = REAL(matrix)[startPoint + i * nRows];
}

void row2NimArr(SEXP matrix, NimArrBase<int>* nimPtr, int startPoint, int len, int nRows){
	for(int i = 0; i < len; i++)
		(*nimPtr)[i] = INTEGER(matrix)[startPoint + i * nRows];	
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


SEXP matrix2VecNimArr(SEXP RvecNimPtr, SEXP matrix, SEXP rowStart, SEXP rowEnd){
	int cRowStart = INTEGER(rowStart)[0] - 1;
	int cRowEnd = INTEGER(rowEnd)[0] - 1;
	NimVecType* vecTypePtr = static_cast<NimVecType*>(R_ExternalPtrAddr(RvecNimPtr) );
	vector<int> rowDims = vecTypePtr->getRowDims(0);
	int len = 0;
	for(unsigned int i = 0; i < rowDims.size(); i++)
		len = len + rowDims[i];		
	int nRows = LENGTH(matrix) / len;
	nimType varType = vecTypePtr->getNimType();
	if(varType == DOUBLE){
		VecNimArrBase<double>* vecPtr = static_cast<VecNimArrBase<double>*>(vecTypePtr);
		for(int i = cRowStart; i <= cRowEnd; i++)
			row2NimArr(matrix, static_cast<NimArrBase<double>*>(vecPtr->getRowTypePtr(i) ), cRowStart + i , len, nRows);
	}
	if(varType == INT){
		VecNimArrBase<int>* vecPtr = static_cast<VecNimArrBase<int>*>(vecTypePtr);
		for(int i = cRowStart; i <= cRowEnd; i++)
			row2NimArr(matrix, static_cast<NimArrBase<int>*>(vecPtr->getRowTypePtr(i) ), cRowStart + i , len, nRows);
	}
	return(R_NilValue);
}


SEXP getEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sIndex){
	SEXP ans;
	int cIndex = INTEGER(sIndex)[0] - 1;
	ans = findVar(install(CHAR(STRING_ELT(sString, cIndex ))), sEnv);
	PROTECT(ans);
	UNPROTECT(1);
	return(ans);
	}

SEXP getEnvVar(SEXP sString, SEXP sEnv){
  	return(getEnvVar_Sindex(sString, sEnv, ScalarInteger(1)));
}

SEXP setEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sVal, SEXP sIndex){
	int cIndex = INTEGER(sIndex)[0] - 1;
	setVar(install(CHAR(STRING_ELT(sString, cIndex ))), sVal, sEnv);
	return(R_NilValue);
}

 SEXP setEnvVar(SEXP sString, SEXP sEnv, SEXP sVal){						
  	return(setEnvVar_Sindex(sString, sEnv, sVal, ScalarInteger(1)));
  }
  
  
/* tools for R's optim	*/


//   FRAMEWORK FOR USING R'S OPTIM FUNCTIONS IN NIMBLE




void bareBonesOptim(NimArr<1, double> initPar, optimfn objFxn, void* nfPtr, int nargs,  ...){
	int n1 = initPar.size();
	NM_OptimControl* nmControl = new NM_OptimControl(1.0, 0.5, 2.0,0.0001, 0.0001, 500, 0);
	OptimAns* optAns = new OptimAns(n1);
	optAns->Fmin = 100;
	vector<void*>* otherArgs = new vector<void*>(nargs);
	va_list vl;
	va_start(vl, nargs);
	for(int i = 0; i < nargs; i++)
		(*otherArgs)[i] = va_arg(vl, void*);
	va_end(vl);
	
    void* vp_otherArgs = static_cast<void*>(otherArgs);
    
	nimble_optim(nfPtr, static_cast<OptimControl*>(nmControl), optAns,
					 initPar, vp_otherArgs, objFxn);
	
	Rprintf("Called bareBonesOptim, final value = %f\n", optAns->Fmin);
	
    delete nmControl;
    delete optAns;
    delete otherArgs;
	
}	
	
// NEW CORE GENERIC FUNCTIONS

//We hope this function will be generic (i.e. not requiring generated function)
//objFxn is the double (int, double*, void*) that we want to optimize
//objFxn has to be built dynamically 
//Also, not 100% sure what to do about the nimbleFunction it self. Right now I'm 
//assuming it will be passed as a void pointer
void nimble_optim(void* nimFun, OptimControl* control, OptimAns* ans,
				 	NimArr<1, double> par, void* otherArgs,
				 	optimfn objFxn){

	//Steps required:	
	//Unpack parts of control (things that are required for all optimizers)
	//Decide optimizer {
	//Based on optimizer, unpack specialized part of optimControl and call appropriate optimizer	
	//  unpack results of optimizer into ans
	//	}
	// ?return ans (depends if we want to allow return of nimble functions)

	// Generic Unpacking	
	int n = par.size();
	double* xin = &(par[0]);	//IF we want to leave par untouched by optim, we could copy these values instead of pointing at them
	double* x = &(ans->par[0]);
	double* Fmin = &(ans->Fmin);
	int* fail = &(ans->fail);
	int* fncount = &(ans->fncount);
	
	vector<void*>* ex = new vector<void*> (2);
	(*ex)[0] = nimFun;
	(*ex)[1] = otherArgs;

	
	//Choosing and applying optimizer
	//Just using Nelder Mead as example here
	if( control->optimType == 1 ){
		NM_OptimControl* nmControl = static_cast<NM_OptimControl*>(control);	
			nmmin(n, xin, x, Fmin, objFxn, fail, nmControl->abstol, nmControl->intol,
			static_cast<void*> (ex), nmControl->alpha, nmControl->beta, nmControl->gamma, nmControl->trace,
			fncount, nmControl->maxit);
	
	}		
	
	
		delete ex;
		
	//actually, nothing else to do at this point! 
	//Could return ans if we decided to build it on the fly here instead of providing it as argument
}	



void nimble_optim_withVarArgs(void* nimFun, OptimControl* control, OptimAns* ans,
				 	NimArr<1, double> par, optimfn objFxn,
				 	int numOtherArgs, ...){

	//Steps required:	
	//Unpack parts of control (things that are required for all optimizers)
	//Decide optimizer {
	//Based on optimizer, unpack specialized part of optimControl and call appropriate optimizer	
	//  unpack results of optimizer into ans
	//	}
	// ?return ans (depends if we want to allow return of nimble functions)

	// Generic Unpacking	
	int n = par.size();
	double* xin = &(par[0]);	//IF we want to leave par untouched by optim, we could copy these values instead of pointing at them
	double* x = &(ans->par[0]);
	double* Fmin = &(ans->Fmin);
	int* fail = &(ans->fail);
	int* fncount = &(ans->fncount);
	vector<void*>* otherArgs = new vector<void*>(numOtherArgs);
	va_list vl;
	va_start(vl, numOtherArgs);
	for(int i = 0; i < numOtherArgs; i++)
		(*otherArgs)[i] = va_arg(vl, void*);
	va_end(vl);
	vector<void*>* ex = new vector<void*> (2);
	(*ex)[0] = nimFun;
	(*ex)[1] = otherArgs;
	//Choosing and applying optimizer
	//Just using Nelder Mead as example here
	if( control->optimType == 1 ){
		NM_OptimControl* nmControl = static_cast<NM_OptimControl*>(control);	
			nmmin(n, xin, x, Fmin, objFxn, fail, nmControl->abstol, nmControl->intol,
			static_cast<void*> (ex), nmControl->alpha, nmControl->beta, nmControl->gamma, nmControl->trace,
			fncount, nmControl->maxit);
	
	}		
		delete ex;
		delete otherArgs;
	//actually, nothing else to do at this point! 
	//Could return ans if we decided to build it on the fly here instead of providing it as argument
}	



/* forwardsolve, backsolve */
//
//NimArr<1, double> forwardsolve(NimArr<2, double> A, NimArr<1, double> b) {
//  NimArr<1, double> ans = NimArr<1, double>();
//  return(ans);
//}
// 
//NimArr<2, double> forwardsolve(NimArr<2, double> A, NimArr<2, double> B) {
//  NimArr<2, double> ans = NimArr<2, double>();
//  return(ans);
//}






