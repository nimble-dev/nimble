#ifndef __RCPPNIMBLEUTILS
#define __RCPPNIMBLEUTILS

#include "NimArrBase.h"
#include "NimArr.h"
#include "RcppUtils.h"

// all of this is in RcppUtils.h
/* #include "R.h" */
/* #include "Utils.h" */
/* #include <string> */
/* #include <vector> */
/* #include<iostream> */
/* #include<sstream> */

/* #include <Rinternals.h> */

/* #include <R_ext/Applic.h>	/\* this is required for optim *\/ */
/* #include <stdarg.h> 		/\* this is required for variable number of arguments *\/ */

using namespace std;

bool all(bool x);
bool any(bool x);
int t(int x);
double t(double x);
int prod(int x);
double prod(double x);

SEXP cGetMVElementOneRow(NimVecType* typePtr, nimType vecType, int index);
//SEXP cGetMVElementOneRow(NimVecType* typePtr, nimType vecType, int nrowCpp, int index);
void cSetMVElementSingle(NimVecType* typePtr, nimType vecType,  int index, SEXP Svalue);
 
//bool checkString(SEXP Ss, int len);
//bool checkNumeric(SEXP Sval, int len);

vector<int> getSEXPdims(SEXP Sx);

extern "C" {
  SEXP setDoublePtrFromSinglePtr(SEXP SdoublePtr, SEXP SsinglePtr);

  //  SEXP setVec(SEXP Sextptr, SEXP Svalue);
  //  SEXP getVec(SEXP Sextptr);
  //  SEXP getVec_Integer(SEXP Sextptr);
  
  SEXP addBlankModelValueRows(SEXP Sextptr, SEXP numAdded);
  SEXP getNRow(SEXP Sextptr);
  SEXP copyModelValuesElements(SEXP SextptrFrom, SEXP SextptrTo, SEXP rowsFrom, SEXP rowsTo);
  SEXP getMVElement(SEXP Sextptr, SEXP Sindex);
  SEXP getMVsize(SEXP Sextptr);

  SEXP getMVElementAsList(SEXP SextPtr, SEXP Sindices);
  SEXP setMVElementFromList(SEXP vNimPtr, SEXP rList, SEXP Sindices);

  SEXP matrix2VecNimArr(SEXP RvecNimPtr, SEXP matrix, SEXP rowStart, SEXP rowEnd);

  //  SEXP printMVElement(SEXP Sextptr, SEXP Sindex);
  SEXP setMVElement(SEXP Sextptr, SEXP Sindex, SEXP Svalue);

  SEXP resizeNumListRow(SEXP Sextptr, SEXP Sindex, SEXP dims); 	// resizes a particular row of a numericlist

//  SEXP setNumList(SEXP Sextptr, SEXP Sindex, SEXP Svalue);   automatically resizes. Might want to use later
   SEXP setNumListRows(SEXP Sextptr, SEXP nRows, SEXP setSize2row1);		// this sets the number of rows in a numericList (really, any VecNimArr)


  //  SEXP setVarPointer(SEXP SextptrModelVar, SEXP SextptrStorageVar, SEXP Srownum);
  SEXP makeNumericList(SEXP nDims, SEXP type, SEXP nRows);

  SEXP newSampObject();	//  Creates our new object from sampleClass (will be generated automatically later)
							//	Just for use in demos
							// 	To get a pointer to an element from sampleClass, use
							//	getModelObjectPtr (from the NamedObjects.cpp file)
					
  SEXP Nim_2_SEXP(SEXP rPtr, SEXP NumRefers);	//	Returns SEXP object with correct data type and dimensions. NumRefers
												//  should be an integer with the number of dereferencing required for rPtr
												//	So if rPtr is a pointer to a NimArr, NumRefers = 1
												//	If rPtr is a pointer to a pointer to a NimArr, NumRefers = 2
												//	Currently only NumRefers = 1 and 2 are allowed, but easily updated
												//	by extending "getNimTypePtr" function
	
  SEXP SEXP_2_Nim(SEXP rPtr, SEXP NumRefers, SEXP rValues, SEXP allowResize); //	Copies values from rValues to NimArr. Same behavior
															  // 	with NumRefers as above. Also, type checking is done
															  // 	by R.internals functions INTEGER and REAL
															  
  //  SEXP Nim_2_Nim(SEXP rPtrFrom, SEXP numRefFrom, SEXP rPtrTo, SEXP numRefTo);	
												//  Copies from one NimArr to another. Type checks
												//	For now, both NimArr's must be either double or int. We can add other
												//  types or allow conversion by extending Nim_2_Nim and the cNim_2_Nim options 

  SEXP setPtrVectorOfPtrs(SEXP SaccessorPtr, SEXP ScontentsPtr, SEXP Ssize);
  SEXP setOnePtrVectorOfPtrs(SEXP SaccessorPtr, SEXP Si, SEXP ScontentsPtr);
  //SEXP getOnePtrVectorOfPtrs(SEXP SaccessorPtr, SEXP Si);
  
  
  
  SEXP getEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sIndex);// This is a utility for looking up a field of an environment
  														 // sString is a character vector with the field name we want
  														 // sEnv is the environment
  														 // sIndex is the index of the sString that contains the name
  														 // we actually want to use. 
  														 //	Look up by sString[sIndex] is done to allow for easy looping
  														 //Important Note: sIndex = 1 looks up the first name (i.e. use R indexing, not C) 
  SEXP getEnvVar(SEXP sString, SEXP sEnv);	// Same as above, but uses sIndex = 1 (i.e. sString is a single character string)
  														 
  SEXP setEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sVal, SEXP sIndex);	//Same as getEnvVar_Sindex, but this function sets rather than gets														
  SEXP setEnvVar(SEXP sString, SEXP sEnv, SEXP sVal);						//Same as above but uses sIndex = 1

  SEXP register_VecNimArr_Finalizer(SEXP Sp, SEXP Dll);
}

void VecNimArr_Finalizer(SEXP Sp);

class vectorOfPtrsAccessBase {
 public:
  virtual void setTheVec(void *tV, int size)=0;
  virtual void setVecPtr(int i, void *v)=0;
  virtual void *getVecPtr(int i)=0;
};

template<class T>
class vectorOfPtrsAccess : public vectorOfPtrsAccessBase {
 public:
  vector<T*> *theVec;
  void setTheVec(void* v, int size) {theVec = static_cast<vector<T*>*>(v); theVec->resize(size);}
  void setVecPtr(int i, void *v) {
    (*theVec)[i] = static_cast<T*>(v);
  }
  void *getVecPtr(int i) {return(static_cast<void *>( (*theVec)[i] ) ); }
};


/*
  Apparently partial specialization of function templates is not allowed.
  So these are witten for doubles, and when we get to integers and logicals we can 
  use overlaoding or different names.
 */
template<int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, double> &ans );
template<int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, int> &ans );

template<int ndim>
SEXP NimArr_2_SEXP(const NimArr<ndim, double> &val);
template<int ndim>
SEXP NimArr_2_SEXP(const NimArr<ndim, int> &val);

template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, double> &ans); 
template<>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, int> &ans); 

template<int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, double> &ans) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_NimArr<ndim> called for SEXP that is not a numeric or logica!\n");
  vector<int> inputDims(getSEXPdims(Sn));
  if(inputDims.size() != ndim) PRINTF("Error: Wrong number of input dimensions in SEXP_2_NimArr<ndim, double> called for SEXP that is not a numeric!\n");
  // if(ans.size() != 0) PRINTF("Error: trying to reset a NimArr that was already sized\n");
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if(isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr() );
  } else {
    if(isInteger(Sn) || isLogical(Sn)) {
      int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
      std::copy(iSn, iSn + nn, ans.getPtr()); //v);
    } else {
      PRINTF("Error: could not handle input type to SEXP_2_NimArr\n");
    }
  }
}

// ACTUALLY THIS IS IDENTICAL CODE TO ABOVE, SO THEY COULD BE COMBINED WITHOUT TEMPLATE SPECIALIZATION
template<int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, int> &ans) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_NimArr<ndim> called for SEXP that is not a numeric or logica!\n");
  vector<int> inputDims(getSEXPdims(Sn));
  if(inputDims.size() != ndim) PRINTF("Error: Wrong number of input dimensions in SEXP_2_NimArr<ndim, double> called for SEXP that is not a numeric!\n");
  // if(ans.size() != 0) PRINTF("Error: trying to reset a NimArr that was already sized\n");
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if(isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr() );
  } else {
    if(isInteger(Sn) || isLogical(Sn)) {
      int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
      std::copy(iSn, iSn + nn, ans.getPtr()); //v);
    } else {
      PRINTF("Error: could not handle input type to SEXP_2_NimArr\n");
    }
  }
}

template<int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, bool> &ans) {
  if(!(isNumeric(Sn) || isLogical(Sn))) PRINTF("Error: SEXP_2_NimArr<ndim> called for SEXP that is not a numeric or logica!\n");
  vector<int> inputDims(getSEXPdims(Sn));
  if(inputDims.size() != ndim) PRINTF("Error: Wrong number of input dimensions in SEXP_2_NimArr<ndim, double> called for SEXP that is not a numeric!\n");
  // if(ans.size() != 0) PRINTF("Error: trying to reset a NimArr that was already sized\n");
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if(isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr() );
  } else {
    if(isInteger(Sn) || isLogical(Sn)) {
      int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
      std::copy(iSn, iSn + nn, ans.getPtr()); //v);
    } else {
      PRINTF("Error: could not handle input type to SEXP_2_NimArr\n");
    }
  }
}


template<int ndim>
SEXP NimArr_2_SEXP(NimArr<ndim, double> &val) {
  SEXP Sans;
  int outputLength = val.size();
  PROTECT(Sans = allocVector(REALSXP, outputLength));
  double *ans = REAL(Sans);

  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if(val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(INTSXP, val.numDims() ) );
    for(int idim = 0; idim < val.numDims(); ++idim) INTEGER(Sdim)[idim] = val.dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return(Sans);
}

template<int ndim>
SEXP NimArr_2_SEXP(NimArr<ndim, int> &val) {
  SEXP Sans;
  int outputLength = val.size();
  PROTECT(Sans = allocVector(INTSXP, outputLength));
  int *ans = INTEGER(Sans);

  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if(val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(INTSXP, val.numDims() ) );
    for(int idim = 0; idim < val.numDims(); ++idim) INTEGER(Sdim)[idim] = val.dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return(Sans);
}

template<int ndim>
SEXP NimArr_2_SEXP(NimArr<ndim, bool> &val) {
  SEXP Sans;
  int outputLength = val.size();
  PROTECT(Sans = allocVector(LGLSXP, outputLength));
  int *ans = LOGICAL(Sans);

  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if(val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(LGLSXP, val.numDims() ) );
    for(int idim = 0; idim < val.numDims(); ++idim) LOGICAL(Sdim)[idim] = val.dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return(Sans);
}



void row2NimArr(SEXP &matrix, NimArrBase<double> &nimPtr, int startPoint, int len, int nRows);
void row2NimArr(SEXP &matrix, NimArrBase<int> &nimPtr, int startPoint, int len, int nRows);

template <class T>
void cNimArr_2_NimArr(NimArrBase<T> &nimFrom, NimArrBase<T> &nimTo);

template <class T1, class T2>
void cNimArr_2_NimArr_Index(NimArrBase<T1> &nimFrom, int fromBegin, NimArrBase<T2> &nimTo, int toBegin, int length);

//void trashNamedObjectsPtr(SEXP rPtr);
//void trashElementPtr(SEXP rPtr);
void sampleFinalizer(SEXP ptr);
template<int nDim, class T>
void VecNimArrFinalizer(SEXP ptr);

double nimMod(double a, double b);

template<class T>
int length(vector<T> vec)
	{
	return(vec.size());
	}
	
/* class orderedPair	//simple class which is used to be sorted by value, but remember what the original order was. used in rawSample */
/* 	{ */
/* 	public: */
/* 	double value; */
/* 	int rank; */
/* 	}; */

/* bool compareOrderedPair(orderedPair a, orderedPair b);	 //function called for sort  */


void rankSample(NimArr<1, double>& weights, int& n, NimArr<1, int>& output);
void rankSample(NimArr<1, double>& weights, int& n, NimArr<1, int>& output, bool& silent);


/*	optim tools	*/


//	 NEW CLASSES (may be classes for which nimble functions can inherit from to allow for easy 


//This is a class that will be used to store outcome of a call to optim
//By giving the optim functions pointers to these elements, 
//there is actually no need to "unpack" the results
//The following two classes COULD be a nimbleFunction
class OptimAns{
	public: 
	NimArr<1, double> par;
	double Fmin;
	int fail, fncount;
	OptimAns(int n){
		par = NimArr<1, double> (n);
	};
};		


//This is a class that will be used to store inputs to various optims
//Will further specialize classes for each version of optim!


class OptimControl{
	public:
	int optimType, maxit, trace;
	//Choices of optimType: 1 = Nelder Mead, 2 = BFG, 3 = BFG with Box Constraints, 4 = Conjugate Gradient, 5 = Simulated Annealing
};		
					
					
//	This is a specialized control specifically for Nelder Mead optimizer						
class NM_OptimControl : public OptimControl{
	public:
	double alpha, beta, gamma, abstol, intol;
	NM_OptimControl(double ialpha, double ibeta, double igamma, double iabstol, double iintol, int imaxit, int itrace){
		alpha = ialpha; beta = ibeta; gamma = igamma; abstol = iabstol; intol = iintol;
		maxit = imaxit; trace = itrace;
		optimType = 1;
	};
};


void bareBonesOptim(NimArr<1, double> initPar, optimfn objFxn, void* nfPtr, int nargs,  ...);

void nimble_optim(void* nimFun, OptimControl* control, OptimAns* ans,
				 	NimArr<1, double> par, void* otherArgs,
				 	optimfn objFxn);

void nimble_optim_withVarArgs(void* nimFun, OptimControl* control, OptimAns* ans,
				 	NimArr<1, double> par, optimfn objFxn,
				 	int numOtherArgs, ...);
					
SEXP getClassElement(SEXP Sobject, const char *name);
void setClassElement(SEXP Sobject, const char *name, SEXP setObject);
SEXP makeNewNimbleList();

#endif


