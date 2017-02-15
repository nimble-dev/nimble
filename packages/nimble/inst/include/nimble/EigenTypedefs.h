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

class EIGEN_EIGENCLASS : public NamedObjects, public pointedToBase {
public:
  NimArr<1, double> values;
  NimArr<2, double> vectors;
  SEXP RObjectPointer;
  bool RCopiedFlag;
  
SEXP  copyToSEXP (   );
void  createNewSEXP (  );
void  copyFromSEXP ( SEXP S_nimList_ );
 EIGEN_EIGENCLASS(){	
    namedObjects["values"]=&values;
	namedObjects["vectors"]=&vectors;
	RCopiedFlag = false;
	RObjectPointer = NULL;
 };
};


template<class Derived>
nimSmartPtr<EIGEN_EIGENCLASS>   EIGEN_EIGEN(const Eigen::MatrixBase<Derived> &x, bool valuesOnly) {
    nimSmartPtr<EIGEN_EIGENCLASS> returnClass = new EIGEN_EIGENCLASS;
	(*returnClass).values.initialize(0, 0, x.rows());
	Map<VectorXd> Eig_eigVals((*returnClass).values.getPtr(),x.rows());
	if(x.rows() != x.cols()) {
	 _nimble_global_output <<"Run-time size error: expected matrix argument to eigen() to be square."<<"\n"; nimble_print_to_R(_nimble_global_output);
	}
    Eigen::DecompositionOptions eigOpts = valuesOnly ? EigenvaluesOnly : ComputeEigenvectors;
	SelfAdjointEigenSolver<MatrixXd> solver(x, eigOpts); // The MatrixXd here doesn't seem generic, but I couldn't get it to work otherwise and it would be odd to do an Eigen decomposition on anything else. -Perry
	Eig_eigVals = solver.eigenvalues().reverse();
	if(!valuesOnly){
	  (*returnClass).vectors.initialize(0, 0, x.rows(), x.cols());
	  Map<MatrixXd> Eig_eigVecs((*returnClass).vectors.getPtr(),x.rows(),x.cols());
	  Eig_eigVecs = solver.eigenvectors().reverse();	
	}
	return(returnClass);
};



class EIGEN_SVDCLASS : public NamedObjects, public pointedToBase {
public:
  NimArr<1, double> d;
  NimArr<2, double> u;
  NimArr<2, double> v;
  SEXP RObjectPointer;
  bool RCopiedFlag;
  
SEXP  copyToSEXP (   );
void  createNewSEXP (  );
 EIGEN_SVDCLASS (  ) {
	namedObjects["d"]=&d;
	namedObjects["u"]=&u;
	namedObjects["v"]=&v;
	RCopiedFlag = false;
	RObjectPointer = NULL;
 };
};

template<class Derived>
nimSmartPtr<EIGEN_SVDCLASS>   EIGEN_SVD(const Eigen::MatrixBase<Derived> &x) {
    nimSmartPtr<EIGEN_SVDCLASS> returnClass = new EIGEN_SVDCLASS;
	int nu = min(x.rows(), x.cols());
	(*returnClass).d.initialize(0, 0, nu);
	(*returnClass).u.initialize(0, 0, x.rows(), nu);
	(*returnClass).v.initialize(0, 0, x.cols(), nu);

	Map<VectorXd> Svd_d((*returnClass).d.getPtr(), nu);
	Map<MatrixXd> Svd_u((*returnClass).u.getPtr(), x.rows(), nu);
	Map<MatrixXd> Svd_v((*returnClass).v.getPtr(), x.cols(), nu);

	/* if(nu < 16){    not currently available in nimble's version of EIGEN library, but may be available after upgrade */		
	JacobiSVD<MatrixXd> svd(x, ComputeThinU | ComputeThinV);
	/* } */
	/* else{ // if minimum dimension length > 16, use bidiagonialization, as recommended on eigen website
		Eigen::BDCSVD<MatrixXd> svd(x, ComputeThinU | ComputeThinV);
	} */
	Svd_d = svd.singularValues();
	Svd_u = svd.matrixU();
	Svd_v = svd.matrixV();
	return(returnClass);
};

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

#endif
