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

#ifndef __RCPPNIMBLEUTILS
#define __RCPPNIMBLEUTILS
#include "NimArrBase.h"
#include "NimArr.h"
#include "RcppUtils.h"

/* The following two macros are for use by copyFromRobject methods 
   in compiled nimbleFunctions. */
#define SETUP_S_xData \
  SEXP S_string_xData; \
  SEXP S_xData; \
  PROTECT(S_string_xData = Rf_allocVector(STRSXP, 1)); \
  SET_STRING_ELT(S_string_xData, 0, Rf_mkChar(".xData")); \
  PROTECT(S_xData = GET_SLOT(Robject, S_string_xData));

#define COPY_NUMERIC_VECTOR_FROM_R_OBJECT(varName) \
  { \
  std::string svarName(varName); \
  SEXP_2_Nim_for_copyFromRobject(getObjectPtr(svarName),\
				 PROTECT(Rf_findVarInFrame(S_xData,	\
							   Rf_install(varName)))); \
  }

#define COPY_VALUE_MAP_ACCESSORS_FROM_NODE_NAMES(varName, derivsEnabled)	\
  { \
  std::string svarName(varName); \
  std::string svarName_AD_(svarName + "_AD_"); \
  populateValueMapAccessorsFromNodeNames_copyFromRobject(getObjectPtr(svarName), \
							 PROTECT(Rf_findVarInFrame(S_xData, \
										   Rf_install(varName))), \
							 derivsEnabled, \
							 getObjectPtr(svarName_AD_, false) \
							 );		\
  }

#define COPY_NODE_FXN_VECTOR_FROM_R_OBJECT(varName) \
  { \
  std::string svarName(varName); \
  populateNodeFxnVectorNew_copyFromRobject(getObjectPtr(svarName),\
					   PROTECT(Rf_findVarInFrame(S_xData, \
								     Rf_install(varName)))); \
                    }

#define COPY_NODE_FXN_VECTOR_DERIVS_FROM_R_OBJECT(varName) \
  { \
  std::string svarName(varName); \
  populateNodeFxnVectorNew_copyFromRobject_forDerivs(getObjectPtr(svarName),\
					   PROTECT(Rf_findVarInFrame(S_xData, \
								     Rf_install(varName)))); \
                    }

#define COPY_NIMBLE_FXN_FROM_R_OBJECT(varName) \
  { \
  std::string svarName(varName); \
  setNimbleFxnPtr_copyFromRobject(getObjectPtr(svarName),\
					    PROTECT(Rf_findVarInFrame(S_xData, \
								      Rf_install(varName)))); \
  }

#define COPY_DOUBLE_SCALAR_FROM_R_OBJECT(varName) \
  { \
  std::string svarName(varName); \
  populate_SEXP_2_double_for_copyFromRobject(getObjectPtr(svarName),\
					     PROTECT(Rf_findVarInFrame(S_xData, \
								     Rf_install(varName)))); \
  }

#define COPY_INTEGER_SCALAR_FROM_R_OBJECT(varName) \
  { \
  std::string svarName(varName); \
  populate_SEXP_2_int_for_copyFromRobject(getObjectPtr(svarName),\
					     PROTECT(Rf_findVarInFrame(S_xData, \
								     Rf_install(varName)))); \
  }

#define COPY_LOGICAL_SCALAR_FROM_R_OBJECT(varName) \
  { \
  std::string svarName(varName); \
  populate_SEXP_2_bool_for_copyFromRobject(getObjectPtr(svarName),\
					     PROTECT(Rf_findVarInFrame(S_xData, \
								     Rf_install(varName)))); \
  }


using namespace std;

bool all(bool x);
bool any(bool x);
int t(int x);
double t(double x);
int prod(int x);
double prod(double x);

SEXP cGetMVElementOneRow(NimVecType* typePtr, nimType vecType, int index);
void cSetMVElementSingle(NimVecType* typePtr, nimType vecType,  int index, SEXP Svalue);

void SEXP_2_Nim_internal(NimArrType* nimTypePtr,
			 SEXP rValues,
			 bool resize);
void SEXP_2_Nim_for_copyFromRobject(void *NimArrPtr, SEXP rValues);

void setNimbleFxnPtr_copyFromRobject(void *nf_to, SEXP S_NF_from);

extern "C" {
  SEXP setDoublePtrFromSinglePtr(SEXP SdoublePtr, SEXP SsinglePtr);
  SEXP setSmartPtrFromSinglePtr(SEXP SdoublePtr, SEXP SsinglePtr);
  SEXP setSmartPtrFromDoublePtr(SEXP SdoublePtr, SEXP SsinglePtr);
  
  SEXP setVecNimArrRows(SEXP Sextptr, SEXP nRows, SEXP setSize2row1);
  SEXP addBlankModelValueRows(SEXP Sextptr, SEXP numAdded);
  SEXP getNRow(SEXP Sextptr);
  SEXP copyModelValuesElements(SEXP SextptrFrom, SEXP SextptrTo, SEXP rowsFrom, SEXP rowsTo);
  SEXP getMVElement(SEXP Sextptr, SEXP Sindex);
  SEXP getMVsize(SEXP Sextptr);

  SEXP getMVElementAsList(SEXP SextPtr, SEXP Sindices);
  SEXP setMVElementFromList(SEXP vNimPtr, SEXP rList, SEXP Sindices);

  SEXP matrix2VecNimArr(SEXP RvecNimPtr, SEXP matrix, SEXP rowStart, SEXP rowEnd);

  SEXP setMVElement(SEXP Sextptr, SEXP Sindex, SEXP Svalue);

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

NimArr<1, double> vectorDouble_2_NimArr(vector<double> input);

/*
  Apparently partial specialization of function templates is not allowed.
  So these are witten for doubles, and when we get to integers and logicals we can 
  use overlaoding or different names.
 */
/* Try overloading these for all needed copy operations */ 
void SEXP_2_NimArr(SEXP Sn, double &x);
void SEXP_2_NimArr(SEXP Sn, int &x);
void SEXP_2_NimArr(SEXP Sn, bool &x);
void SEXP_2_NimArr(SEXP Sn, std::string &x);
void SEXP_2_NimArr(SEXP Sn, std::vector<std::string> &x);

SEXP NimArr_2_SEXP(double x);
SEXP NimArr_2_SEXP(int x);
SEXP NimArr_2_SEXP(bool x);
SEXP NimArr_2_SEXP(std::string &x);
SEXP NimArr_2_SEXP(const std::vector<std::string> &x);

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
  NIM_ASSERT3(Rf_isNumeric(Sn) || Rf_isLogical(Sn),
    "SEXP_2_NimArr<%d, double> called for SEXP that is not a numeric or logical: actual type %s\n",
    ndim, Rf_type2char(TYPEOF(Sn)));
  vector<int> inputDims(getSEXPdims(Sn));
  NIM_ASSERT4(inputDims.size() == ndim,
    "Wrong number of input dimensions in SEXP_2_NimArr<%d, double> called for SEXP that is not a numeric: expected %d, actual %d\n",
    ndim, ndim, static_cast<int>(inputDims.size()));
  // NIM_ASSERT(ans.size() == 0, "trying to reset a NimArr that was already sized\n");
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if(Rf_isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr() );
  } else {
    NIM_ASSERT3(Rf_isInteger(Sn) || Rf_isLogical(Sn),
      "could not handle input of type %s to SEXP_2_NimArr<%d, double>\n",
      Rf_type2char(TYPEOF(Sn)), ndim);
    int *iSn = Rf_isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    std::copy(iSn, iSn + nn, ans.getPtr()); //v);
  }
}

// ACTUALLY THIS IS IDENTICAL CODE TO ABOVE, SO THEY COULD BE COMBINED WITHOUT TEMPLATE SPECIALIZATION
template<int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, int> &ans) {
  NIM_ASSERT3(Rf_isNumeric(Sn) || Rf_isLogical(Sn),
    "SEXP_2_NimArr<%d, int> called for SEXP that is not a numeric or logical: actual type %s\n",
    ndim, Rf_type2char(TYPEOF(Sn)));
  vector<int> inputDims(getSEXPdims(Sn));
  NIM_ASSERT4(inputDims.size() == ndim,
    "Wrong number of input dimensions in SEXP_2_NimArr<%d, int> called for SEXP that is not a numeric: expected %d, actual %d\n",
    ndim, ndim, inputDims.size());
  // NIM_ASSERT(ans.size() == 0, "trying to reset a NimArr that was already sized\n");
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if(Rf_isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr() );
  } else {
    NIM_ASSERT3(Rf_isInteger(Sn) || Rf_isLogical(Sn),
      "could not handle input type %s to SEXP_2_NimArr<%d, int>\n",
      Rf_type2char(TYPEOF(Sn)), ndim);
    int *iSn = Rf_isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    std::copy(iSn, iSn + nn, ans.getPtr()); //v);
  }
}

template<int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, bool> &ans) {
  NIM_ASSERT3(Rf_isNumeric(Sn) || Rf_isLogical(Sn),
    "SEXP_2_NimArr<%d, bool> called for SEXP that is not a numeric or logical: actual type %s\n",
    ndim, Rf_type2char(TYPEOF(Sn)));
  vector<int> inputDims(getSEXPdims(Sn));
  NIM_ASSERT4(inputDims.size() == ndim,
    "Wrong number of input dimensions in SEXP_2_NimArr<%d, bool> called for SEXP that is not a numeric: expected %d, actual %d\n",
    ndim, ndim, inputDims.size());
  // NIM_ASSERT(ans.size() == 0, "trying to reset a NimArr that was already sized\n");
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if(Rf_isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr() );
  } else {
    NIM_ASSERT3(Rf_isInteger(Sn) || Rf_isLogical(Sn),
      "could not handle input type %s to SEXP_2_NimArr<%d, bool>\n",
      Rf_type2char(TYPEOF(Sn)), ndim);
    int *iSn = Rf_isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    std::copy(iSn, iSn + nn, ans.getPtr()); //v);
  }
}


template<int ndim>
SEXP NimArr_2_SEXP(NimArr<ndim, double> &val) {
  SEXP Sans;
  int outputLength = val.size();
  PROTECT(Sans = Rf_allocVector(REALSXP, outputLength));
  double *ans = REAL(Sans);

  NimArr_map_2_allocatedMemory(val, &ans, outputLength);
  if(val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = Rf_allocVector(INTSXP, val.numDims() ) );
    for(int idim = 0; idim < val.numDims(); ++idim) INTEGER(Sdim)[idim] = val.dimSize(idim);
    Rf_setAttrib(Sans, R_DimSymbol, Sdim);
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
  PROTECT(Sans = Rf_allocVector(INTSXP, outputLength));
  int *ans = INTEGER(Sans);

  NimArr_map_2_allocatedMemory(val, &ans, outputLength);
  //  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if(val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = Rf_allocVector(INTSXP, val.numDims() ) );
    for(int idim = 0; idim < val.numDims(); ++idim) INTEGER(Sdim)[idim] = val.dimSize(idim);
    Rf_setAttrib(Sans, R_DimSymbol, Sdim);
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
  PROTECT(Sans = Rf_allocVector(LGLSXP, outputLength));
  int *ans = LOGICAL(Sans);

  NimArr_map_2_allocatedMemory(val, &ans, outputLength);
  //  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if(val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = Rf_allocVector(LGLSXP, val.numDims() ) );
    for(int idim = 0; idim < val.numDims(); ++idim) LOGICAL(Sdim)[idim] = val.dimSize(idim);
    Rf_setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return(Sans);
}

template<int ndim>
void SEXP_list_2_NimArr_double_vec(SEXP Sn, vector<NimArr<ndim, double> > &ans) {
			int numListElements = Rf_length(Sn);
			ans.resize(numListElements);
			for(int i = 0; i < numListElements; i++){
				SEXP_2_NimArr<ndim>(VECTOR_ELT(Sn, i), ans[i]);
			}
}

template<int ndim>
void SEXP_list_2_NimArr_int_vec(SEXP Sn, vector<NimArr<ndim, int> > &ans) {
			int numListElements = Rf_length(Sn);
			ans.resize(numListElements);
			for(int i = 0; i < numListElements; i++){
				SEXP_2_NimArr<ndim>(VECTOR_ELT(Sn, i), ans[i]);
			}
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

void rankSample(NimArr<1, double>& weights, int n, NimArr<1, int>& output);
void rankSample(NimArr<1, double>& weights, int n, NimArr<1, int>& output, bool silent);

#endif
