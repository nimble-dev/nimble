#ifndef __NIMBLE_EIGEN_TYPEDEFS
#define __NIMBLE_EIGEN_TYPEDEFS
#include <Eigen/Dense>
#include <iostream> // must go before other things because R defines a "length" macro
#include <nimble/NamedObjects.h>
#include <nimble/NimArr.h>
#include <nimble/smartPtrs.h>
#include <nimble/RcppNimbleUtils.h>
#include <nimble/eigenUsingClasses.h>

using namespace Eigen;
#include "nimbleEigen.h"


typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;
typedef Array<bool, Dynamic, Dynamic> ArrayXXb;

typedef Stride<Dynamic, Dynamic> EigStrDyn;
typedef Map<MatrixXd, Unaligned, EigStrDyn > EigenMapStrd;
typedef Map<MatrixXi, Unaligned, EigStrDyn > EigenMapStri;
typedef Map<MatrixXb, Unaligned, EigStrDyn > EigenMapStrb;

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

class EIGEN_EIGENCLASS : public EIGEN_EIGENCLASS_R, public NamedObjects {
public:
  bool RCopiedFlag;
  SEXP  copyToSEXP();
  EIGEN_EIGENCLASS();
};

template<class Derived>
nimSmartPtr<EIGEN_EIGENCLASS>   EIGEN_EIGEN(const Eigen::MatrixBase<Derived> &x, bool valuesOnly) {
    nimSmartPtr<EIGEN_EIGENCLASS> returnClass = new EIGEN_EIGENCLASS;
    EIGEN_EIGEN_INTERNAL(x, valuesOnly, returnClass.getRealPtr());
    return(returnClass);
}

template<class Derived>
nimSmartPtr<EIGEN_EIGENCLASS_R>   EIGEN_EIGEN_R(const Eigen::MatrixBase<Derived> &x, bool valuesOnly) {
    nimSmartPtr<EIGEN_EIGENCLASS_R> returnClass = new EIGEN_EIGENCLASS_R;
    EIGEN_EIGEN_INTERNAL(x, valuesOnly, returnClass.getRealPtr());
    return(returnClass);
}

template<class Derived>
void EIGEN_EIGEN_INTERNAL(const Eigen::MatrixBase<Derived> &x, bool valuesOnly, EIGEN_EIGENCLASS_R *returnClass) {
  // nimSmartPtr<EIGEN_EIGENCLASS> returnClass = new EIGEN_EIGENCLASS;
  returnClass->getValues().initialize(0, 0, x.rows());
	Map<VectorXd> Eig_eigVals(returnClass->getValues().getPtr(),x.rows());
	if(x.rows() != x.cols()) {
	 _nimble_global_output <<"Run-time size error: expected matrix argument to eigen() to be square."<<"\n"; nimble_print_to_R(_nimble_global_output);
	}
    Eigen::DecompositionOptions eigOpts = valuesOnly ? EigenvaluesOnly : ComputeEigenvectors;
	SelfAdjointEigenSolver<MatrixXd> solver(x, eigOpts); // The MatrixXd here doesn't seem generic, but I couldn't get it to work otherwise and it would be odd to do an Eigen decomposition on anything else. -Perry
	Eig_eigVals = solver.eigenvalues();
	if(!valuesOnly){
	  returnClass->getVectors().initialize(0, 0, x.rows(), x.cols());
	  Map<MatrixXd> Eig_eigVecs(returnClass->getVectors().getPtr(),x.rows(),x.cols());
	  Eig_eigVecs = solver.eigenvectors().rowwise().reverse();	
	}
	//	return(returnClass);
};

/* template<class Derived> */
/* nimSmartPtr<EIGEN_SVDCLASS>   EIGEN_SVD(const Eigen::MatrixBase<Derived> &x, int vectors) { */
/*     nimSmartPtr<EIGEN_SVDCLASS> returnClass = new EIGEN_SVDCLASS; */
/* 	int nu = min(x.rows(), x.cols()); */
/* 	(*returnClass).d.initialize(0, 0, nu); */
/* 	Map<VectorXd> Svd_d((*returnClass).d.getPtr(), nu); */
/* 	JacobiSVD<MatrixXd> svd; */

/* 	/\* note: if nu > 16, bidiagonialization algo. is recommended on eigen website.  not currently available w/ nimble's version of eigen, but may be in future. *\/ */
/* 	if(vectors == 0){ */
/* 	  svd.compute(x); */
/* 	} */
/* 	else{ */
/* 	  int leftSVs = nu; */
/* 	  int rightSVs = nu; */
/* 	  if(vectors == 1){ */
/* 	  	svd.compute(x, ComputeThinU | ComputeThinV); */
/* 	  } */
/* 	  if(vectors == 2){ */
/* 	    leftSVs = x.rows(); */
/* 		rightSVs = x.cols(); */
/* 		svd.compute(x, ComputeFullU | ComputeFullV); */
/* 	  } */
/* 	  (*returnClass).u.initialize(0, 0, x.rows(), leftSVs); */
/* 	  (*returnClass).v.initialize(0, 0, x.cols(), rightSVs); */
/* 	  Map<MatrixXd> Svd_u((*returnClass).u.getPtr(), x.rows(), leftSVs); */
/* 	  Map<MatrixXd> Svd_v((*returnClass).v.getPtr(), x.cols(), rightSVs); */
/* 	  Svd_u = svd.matrixU(); */
/* 	  Svd_v = svd.matrixV(); */
/* 	} */
/* 	Svd_d = svd.singularValues();  */
/* 	return(returnClass); */
/* }; */


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
