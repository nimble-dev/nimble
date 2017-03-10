#ifndef __NIMBLE_EIGEN_TYPEDEFS
#define __NIMBLE_EIGEN_TYPEDEFS
#include <Eigen/Dense>
#include <iostream> // must go before other things because R defines a "length" macro
#include <nimble/NamedObjects.h>
#include <nimble/NimArr.h>
#include <nimble/smartPtrs.h>
#include <nimble/RcppNimbleUtils.h>

using namespace Eigen;
#include "nimbleEigen.h"


typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;
typedef Array<bool, Dynamic, Dynamic> ArrayXXb;

typedef Stride<Dynamic, Dynamic> EigStrDyn;
typedef Map<MatrixXd, Unaligned, EigStrDyn > EigenMapStrd;
typedef Map<MatrixXi, Unaligned, EigStrDyn > EigenMapStri;
typedef Map<MatrixXb, Unaligned, EigStrDyn > EigenMapStrb;

vector<int> getSEXPdims(SEXP Sx);

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

SEXP makeNewNimbleList(SEXP S_listName);

//#define EIGEN_FS(x,y)       (x).triangularView<Eigen::Lower>().solve(y)
//#define EIGEN_BS(x,y)       (x).triangularView<Eigen::Upper>().solve(y)
//#define EIGEN_SOLVE(x,y)    (x).lu().solve(y)

/* MatrixXd EIGEN_FS(const Eigen::Map<MatrixXd, Unaligned, Stride<Dynamic, Dynamic> >& x, const Eigen::Map<MatrixXd, Unaligned, Stride<Dynamic, Dynamic> >& y) { */
/*   MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble */
/*   MatrixXd ans = x.triangularView<Eigen::Lower>().solve(ycopy); */
/*   return(ans); */
/* } */

// Diagonal: Tried using Eigen::DiagonalMatrix<> but one can't do much with this.  It doesn't have same flexibility and member functions as a regular Eigen::MatrixXd
// So we'll use a NullaryExpr approach instead


/* extern "C" { */
/*   SEXP C_nimEigen(SEXP, SEXP, SEXP); */
/* } */

template<class derived1, class derived2>
MatrixXd EIGEN_FS(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd ans = x.template triangularView<Eigen::Lower>().solve(ycopy);
  return(ans);
}

template<class derived1, class derived2>
MatrixXd EIGEN_BS(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd ans = x.template triangularView<Eigen::Upper>().solve(ycopy);
  return(ans);
}

template<class derived1, class derived2>
MatrixXd EIGEN_SOLVE(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd ans = x.lu().solve(ycopy);
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const ArrayBase<Derived1>& v1, const ArrayBase<Derived2>& v2) { 
  double ans = (v1 * v2).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const MatrixBase<Derived1>& v1, const MatrixBase<Derived2>& v2) {
  double ans = (v1.cwiseProduct(v2)).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const MatrixBase<Derived1>& v1, const ArrayBase<Derived2>& v2) {
  double ans = (v1.cwiseProduct(v2.matrix())).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const ArrayBase<Derived1>& v1, const MatrixBase<Derived2>& v2) {
  double ans = (v2.cwiseProduct(v1.matrix())).sum();
  return(ans);
}

template <typename Derived>
double sd(const ArrayBase<Derived>& v) { // array only for scalar -
  double mean = v.mean();
  double n = v.size();
  double ans = sqrt((v-mean).matrix().squaredNorm() / (n-1.));
  return(ans);
}

template <typename Derived>
double var(const ArrayBase<Derived>& v) { // array only for scalar -
  double mean = v.mean();
  double n = v.size();
  double ans = (v-mean).matrix().squaredNorm() / (n-1.);
  return(ans);
}

template <typename Derived>
double logdet(const MatrixBase<Derived>& v) {
  return(log(v.determinant()));
}



SEXP makeNewNimbleList(SEXP S_listName);


#endif
