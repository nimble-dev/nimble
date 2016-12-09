#ifndef __NIMBLE_EIGEN_TYPEDEFS
#define __NIMBLE_EIGEN_TYPEDEFS
#include <Eigen/Dense>
#include <iostream> // must go before other things because R defines a "length" macro
#include <nimble/NamedObjects.h>
#include <nimble/NimArr.h>
#include <nimble/smartPtrs.h>
#include <nimble/RcppNimbleUtils.h>

using namespace Eigen;

typedef Stride<Dynamic, Dynamic> EigStrDyn;
typedef Map<MatrixXd, Unaligned, EigStrDyn > EigenMapStr;

//#define EIGEN_FS(x,y)       (x).triangularView<Eigen::Lower>().solve(y)
//#define EIGEN_BS(x,y)       (x).triangularView<Eigen::Upper>().solve(y)
//#define EIGEN_SOLVE(x,y)    (x).lu().solve(y)

/* MatrixXd EIGEN_FS(const Eigen::Map<MatrixXd, Unaligned, Stride<Dynamic, Dynamic> >& x, const Eigen::Map<MatrixXd, Unaligned, Stride<Dynamic, Dynamic> >& y) { */
/*   MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble */
/*   MatrixXd ans = x.triangularView<Eigen::Lower>().solve(ycopy); */
/*   return(ans); */
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
 EIGEN_EIGENCLASS (  );
};

template<class T>
nimSmartPtr<EIGEN_EIGENCLASS>   EIGEN_EIGEN(NimArr<2, T> &x, bool valuesOnly) {
    nimSmartPtr<EIGEN_EIGENCLASS> returnClass = new EIGEN_EIGENCLASS;
	(*returnClass).values.initialize(0, 0, x.dim()[0]);
	Map<VectorXd> Eig_eigVals(0,0);
	new (&Eig_eigVals) Map< VectorXd >((*returnClass).values.getPtr(),x.dim()[0]);
	Map<MatrixXd> Eig_x(x.getPtr(),x.dim()[0],x.dim()[1]);
	if(x.dim()[0] != x.dim()[1]) {
	 _nimble_global_output <<"Run-time size error: expected matrix argument to eigen() to be square."<<"\n"; nimble_print_to_R(_nimble_global_output);
	}	
	EigenSolver<MatrixXd> solver(Eig_x, !valuesOnly);
	Eig_eigVals = solver.eigenvalues().real();
	if(!valuesOnly){
	    (*returnClass).vectors.initialize(0, 0, x.dim()[0], x.dim()[0]);
		Map<MatrixXd> Eig_eigVecs(0,0,0);
		new (&Eig_eigVecs) Map< MatrixXd >((*returnClass).vectors.getPtr(),x.dim()[0],x.dim()[0]);
		Eig_eigVecs = solver.eigenvectors().real();	
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
EIGEN_SVDCLASS (  );
};

template<class T>
nimSmartPtr<EIGEN_SVDCLASS>   EIGEN_SVD(NimArr<2, T> &x) {
    nimSmartPtr<EIGEN_SVDCLASS> returnClass = new EIGEN_SVDCLASS;
	int nu = min(x.dim()[0], x.dim()[1]);
	(*returnClass).d.initialize(0, 0, nu);
	(*returnClass).u.initialize(0, 0, x.dim()[0], nu);
	(*returnClass).v.initialize(0, 0, x.dim()[1], nu);

	Map<VectorXd> Svd_d(0,0);
	new (&Svd_d) Map< VectorXd >((*returnClass).d.getPtr(), nu);
	Map<MatrixXd> Svd_u(0,0,0);
	new (&Svd_u) Map< VectorXd >((*returnClass).u.getPtr(), x.dim()[0], nu);
	Map<MatrixXd> Svd_v(0,0,0);
	new (&Svd_v) Map< VectorXd >((*returnClass).v.getPtr(), x.dim()[1], nu);
	
	Map<MatrixXd> Svd_x(x.getPtr(),x.dim()[0],x.dim()[1]);
	/* if(nu < 16){   */
		JacobiSVD<MatrixXd> svd(Svd_x, ComputeThinU | ComputeThinV);
	/* } */
	/* else{ // if minimum dimension length > 16, use bidiagonialization, as recommended on eigen website
		Eigen::BDCSVD<MatrixXd> svd(Svd_x, ComputeThinU | ComputeThinV);
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
