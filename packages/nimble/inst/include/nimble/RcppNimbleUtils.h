#ifndef __RCPPNIMBLEUTILS
#define __RCPPNIMBLEUTILS
#include "NimArr.h"
#include "NimArrBase.h"
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
/* #include <stdarg.h> 		/\* this is required for variable number of
 * arguments *\/ */

using namespace std;

bool all(bool x);
bool any(bool x);
int t(int x);
double t(double x);
int prod(int x);
double prod(double x);

SEXP cGetMVElementOneRow(NimVecType *typePtr, nimType vecType, int index);
void cSetMVElementSingle(NimVecType *typePtr, nimType vecType, int index,
                         SEXP Svalue);

extern "C" {
SEXP setDoublePtrFromSinglePtr(SEXP SdoublePtr, SEXP SsinglePtr);
SEXP setSmartPtrFromSinglePtr(SEXP SdoublePtr, SEXP SsinglePtr);
SEXP setSmartPtrFromDoublePtr(SEXP SdoublePtr, SEXP SsinglePtr);

SEXP setVecNimArrRows(SEXP Sextptr, SEXP nRows, SEXP setSize2row1);
SEXP addBlankModelValueRows(SEXP Sextptr, SEXP numAdded);
SEXP getNRow(SEXP Sextptr);
SEXP copyModelValuesElements(SEXP SextptrFrom, SEXP SextptrTo, SEXP rowsFrom,
                             SEXP rowsTo);
SEXP getMVElement(SEXP Sextptr, SEXP Sindex);
SEXP getMVsize(SEXP Sextptr);

SEXP getMVElementAsList(SEXP SextPtr, SEXP Sindices);
SEXP setMVElementFromList(SEXP vNimPtr, SEXP rList, SEXP Sindices);

SEXP matrix2VecNimArr(SEXP RvecNimPtr, SEXP matrix, SEXP rowStart, SEXP rowEnd);

SEXP setMVElement(SEXP Sextptr, SEXP Sindex, SEXP Svalue);

// Creates our new object from sampleClass (will be generated automatically
// later).
// Just for use in demos.
// To get a pointer to an element from sampleClass, use
// getModelObjectPtr (from the NamedObjects.cpp file).
SEXP newSampObject();

// Returns SEXP object with correct data type and dimensions.  NumRefers
// should be an integer with the number of dereferencing required for rPtr
// So if rPtr is a pointer to a NimArr, NumRefers = 1
// If rPtr is a pointer to a pointer to a NimArr, NumRefers = 2
// Currently only NumRefers = 1 and 2 are allowed, but easily updated
// by extending "getNimTypePtr" function
SEXP Nim_2_SEXP(SEXP rPtr, SEXP NumRefers);

//	Copies values from rValues to NimArr.
// Same behavior with NumRefers as above. Also, type checking is done by
// R.internals functions INTEGER and REAL
SEXP SEXP_2_Nim(SEXP rPtr, SEXP NumRefers, SEXP rValues, SEXP allowResize);

SEXP setPtrVectorOfPtrs(SEXP SaccessorPtr, SEXP ScontentsPtr, SEXP Ssize);
SEXP setOnePtrVectorOfPtrs(SEXP SaccessorPtr, SEXP Si, SEXP ScontentsPtr);

// This is a utility for looking up a field of an environment sString is a
// character vector with the field name we want sEnv is the environment sIndex
// is the index of the sString that contains the name we actually want to use.
// Look up by sString[sIndex] is done to allow for easy looping
// Important Note: sIndex = 1 looks up the first name (i.e. use R indexing, not
// C).
SEXP getEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sIndex);

// Same as above, but uses sIndex = 1 (i.e. sString is a single character
// string)
SEXP getEnvVar(SEXP sString, SEXP sEnv);

// Same as getEnvVar_Sindex, but this function sets rather than gets
SEXP setEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sVal, SEXP sIndex);

// Same as above but uses sIndex = 1
SEXP setEnvVar(SEXP sString, SEXP sEnv, SEXP sVal);

SEXP register_VecNimArr_Finalizer(SEXP Sp, SEXP Dll);
}

void VecNimArr_Finalizer(SEXP Sp);

class vectorOfPtrsAccessBase {
 public:
  virtual void setTheVec(void *tV, int size) = 0;
  virtual void setVecPtr(int i, void *v) = 0;
  virtual void *getVecPtr(int i) = 0;
};

template <class T>
class vectorOfPtrsAccess : public vectorOfPtrsAccessBase {
 public:
  vector<T *> *theVec;
  void setTheVec(void *v, int size) {
    theVec = static_cast<vector<T *> *>(v);
    theVec->resize(size);
  }
  void setVecPtr(int i, void *v) { (*theVec)[i] = static_cast<T *>(v); }
  void *getVecPtr(int i) { return (static_cast<void *>((*theVec)[i])); }
};

NimArr<1, double> vectorDouble_2_NimArr(vector<double> input);

/*
  Apparently partial specialization of function templates is not allowed.
  So these are witten for doubles, and when we get to integers and logicals we
  can use overlaoding or different names.
 */
template <int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, double> &ans);
template <int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, int> &ans);

template <int ndim>
SEXP NimArr_2_SEXP(const NimArr<ndim, double> &val);
template <int ndim>
SEXP NimArr_2_SEXP(const NimArr<ndim, int> &val);

template <>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, double> &ans);
template <>
void SEXP_2_NimArr<1>(SEXP Sn, NimArr<1, int> &ans);

template <int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, double> &ans) {
  NIM_ASSERT3(isNumeric(Sn) || isLogical(Sn),
              "SEXP_2_NimArr<%d, double> called for SEXP that is not a numeric "
              "or logical: actual type %s\n",
              ndim, type2str(TYPEOF(Sn)));
  vector<int> inputDims(getSEXPdims(Sn));
  NIM_ASSERT4(inputDims.size() == ndim,
              "Wrong number of input dimensions in SEXP_2_NimArr<%d, double> "
              "called for SEXP that is not a numeric: expected %d, actual %d\n",
              ndim, ndim, inputDims.size());
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if (isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());
  } else {
    NIM_ASSERT3(
        isInteger(Sn) || isLogical(Sn),
        "could not handle input of type %s to SEXP_2_NimArr<%d, double>\n",
        type2str(TYPEOF(Sn)), ndim);
    int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    std::copy(iSn, iSn + nn, ans.getPtr());
  }
}

// ACTUALLY THIS IS IDENTICAL CODE TO ABOVE, SO THEY COULD BE COMBINED WITHOUT
// TEMPLATE SPECIALIZATION
template <int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, int> &ans) {
  NIM_ASSERT3(isNumeric(Sn) || isLogical(Sn),
              "SEXP_2_NimArr<%d, int> called for SEXP that is not a numeric or "
              "logical: actual type %s\n",
              ndim, type2str(TYPEOF(Sn)));
  vector<int> inputDims(getSEXPdims(Sn));
  NIM_ASSERT4(inputDims.size() == ndim,
              "Wrong number of input dimensions in SEXP_2_NimArr<%d, int> "
              "called for SEXP that is not a numeric: expected %d, actual %d\n",
              ndim, ndim, inputDims.size());
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if (isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());
  } else {
    NIM_ASSERT3(isInteger(Sn) || isLogical(Sn),
                "could not handle input type %s to SEXP_2_NimArr<%d, int>\n",
                type2str(TYPEOF(Sn)), ndim);
    int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    std::copy(iSn, iSn + nn, ans.getPtr());
  }
}

template <int ndim>
void SEXP_2_NimArr(SEXP Sn, NimArr<ndim, bool> &ans) {
  NIM_ASSERT3(isNumeric(Sn) || isLogical(Sn),
              "SEXP_2_NimArr<%d, bool> called for SEXP that is not a numeric "
              "or logical: actual type %s\n",
              ndim, type2str(TYPEOF(Sn)));
  vector<int> inputDims(getSEXPdims(Sn));
  NIM_ASSERT4(inputDims.size() == ndim,
              "Wrong number of input dimensions in SEXP_2_NimArr<%d, bool> "
              "called for SEXP that is not a numeric: expected %d, actual %d\n",
              ndim, ndim, inputDims.size());
  ans.setSize(inputDims);
  int nn = LENGTH(Sn);
  if (isReal(Sn)) {
    std::copy(REAL(Sn), REAL(Sn) + nn, ans.getPtr());
  } else {
    NIM_ASSERT3(isInteger(Sn) || isLogical(Sn),
                "could not handle input type %s to SEXP_2_NimArr<%d, bool>\n",
                type2str(TYPEOF(Sn)), ndim);
    int *iSn = isInteger(Sn) ? INTEGER(Sn) : LOGICAL(Sn);
    std::copy(iSn, iSn + nn, ans.getPtr());
  }
}

template <int ndim>
SEXP NimArr_2_SEXP(NimArr<ndim, double> &val) {
  SEXP Sans;
  int outputLength = val.size();
  PROTECT(Sans = allocVector(REALSXP, outputLength));
  double *ans = REAL(Sans);

  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if (val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(INTSXP, val.numDims()));
    for (int idim = 0; idim < val.numDims(); ++idim)
      INTEGER(Sdim)[idim] = val.dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return (Sans);
}

template <int ndim>
SEXP NimArr_2_SEXP(NimArr<ndim, int> &val) {
  SEXP Sans;
  int outputLength = val.size();
  PROTECT(Sans = allocVector(INTSXP, outputLength));
  int *ans = INTEGER(Sans);

  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if (val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(INTSXP, val.numDims()));
    for (int idim = 0; idim < val.numDims(); ++idim)
      INTEGER(Sdim)[idim] = val.dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return (Sans);
}

template <int ndim>
SEXP NimArr_2_SEXP(NimArr<ndim, bool> &val) {
  SEXP Sans;
  int outputLength = val.size();
  PROTECT(Sans = allocVector(LGLSXP, outputLength));
  int *ans = LOGICAL(Sans);

  std::copy(val.getPtr(), val.getPtr() + outputLength, ans);
  if (val.numDims() > 1) {
    SEXP Sdim;
    PROTECT(Sdim = allocVector(LGLSXP, val.numDims()));
    for (int idim = 0; idim < val.numDims(); ++idim)
      LOGICAL(Sdim)[idim] = val.dimSize(idim);
    setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  } else {
    UNPROTECT(1);
  }
  return (Sans);
}

void row2NimArr(SEXP &matrix, NimArrBase<double> &nimPtr, int startPoint,
                int len, int nRows);
void row2NimArr(SEXP &matrix, NimArrBase<int> &nimPtr, int startPoint, int len,
                int nRows);

template <class T>
void cNimArr_2_NimArr(NimArrBase<T> &nimFrom, NimArrBase<T> &nimTo);

template <class T1, class T2>
void cNimArr_2_NimArr_Index(NimArrBase<T1> &nimFrom, int fromBegin,
                            NimArrBase<T2> &nimTo, int toBegin, int length);

void sampleFinalizer(SEXP ptr);
template <int nDim, class T>
void VecNimArrFinalizer(SEXP ptr);

double nimMod(double a, double b);

template <class T>
int length(vector<T> vec) {
  return (vec.size());
}

void rankSample(NimArr<1, double> &weights, int &n, NimArr<1, int> &output);
void rankSample(NimArr<1, double> &weights, int &n, NimArr<1, int> &output,
                bool &silent);

#endif
