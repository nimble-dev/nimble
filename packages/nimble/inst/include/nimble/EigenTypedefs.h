#ifndef __NIMBLE_EIGEN_TYPEDEFS
#define __NIMBLE_EIGEN_TYPEDEFS

#include <Eigen/Dense>
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

template<class T>
class EIGEN_EIGENCLASS : public NamedObjects, public pointedToBase {
public:
  double nlScalar;
  NimArr<1, double> values;
  NimArr<2, double> vectors;
void  copyFromSEXP ( SEXP S_nimList_ );
void  copyToSEXP ( SEXP S_nimList_ );
SEXP  writeToSEXP (  );
 EIGEN_EIGENCLASS (  );
};


template<class derived1>
nimSmartPtr<EIGEN_EIGENCLASS>   EIGEN_EIGEN(const MatrixBase<derived1> &x, bool valuesOnly) {
  MatrixXd xcopy = x;
  nimSmartPtr<EIGEN_EIGENCLASS> returnClass = new EIGEN_EIGENCLASS;
  Map<MatrixXd> Eig_eigVals(0,0,0);
  new (&Eig_eigVals) Map< MatrixXd >((*returnClass).values.getPtr(),x.dim()[0],1);
  EigenSolver<MatrixXd> solver = EigenSolver<MatrixXd>(xcopy, valuesOnly);
  Eig_eigVals = solver.eigenvectors().real();
  if(!valuesOnly){
	Map<MatrixXd> Eig_eigVecs(0,0,0);
	new (&Eig_eigVecs) Map< MatrixXd >((*returnClass).vectors.getPtr(),x.dim()[0],x.dim()[0]);
	Eig_eigVecs = solver.eigenvalues().real();
  }
  return(returnClass);
};

/* template<class derived1>
MatrixXd EIGEN_EIGENVECS(const MatrixBase<derived1> &x) {
  MatrixXd xcopy = x;	
  MatrixXd ans = EigenSolver<MatrixXd>(xcopy, true).eigenvectors().real();
  return(ans); 
}

template<class derived1>
VectorXd EIGEN_EIGENVALS(const MatrixBase<derived1> &x) {
  MatrixXd xcopy = x;	
  VectorXd ans = EigenSolver<MatrixXd>(xcopy, false).eigenvalues().real();
  return(ans); 
} */

template<class derived1>
VectorXd EIGEN_SVDD(const MatrixBase<derived1> &x) {
  MatrixXd xcopy = x;	  
  VectorXd ans = JacobiSVD<MatrixXd>(xcopy).singularValues();
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

#endif
