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

template <typename Type>
struct EigenTemplateTypes {
  typedef Matrix<Type, Dynamic, Dynamic> typeMatrixXd;
  typedef Map<typeMatrixXd, Unaligned, EigStrDyn > typeEigenMapStrd;
  typedef Eigen::Map<const typeMatrixXd, Unaligned, EigStrDyn > typeEigenConstMapStrd;
};

typedef typename EigenTemplateTypes<CppAD::AD<double> >::typeMatrixXd MatrixXd_CppAD;

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
  void  resetFlags();
  EIGEN_EIGENCLASS();
};

class EIGEN_SVDCLASS : public EIGEN_SVDCLASS_R, public NamedObjects {
public:
  bool RCopiedFlag;
  SEXP copyToSEXP();
  void  resetFlags();
  EIGEN_SVDCLASS();
};


template<class Derived>
nimSmartPtr<EIGEN_EIGENCLASS>   EIGEN_EIGEN(const Eigen::MatrixBase<Derived> &x,  bool symmetric, bool valuesOnly) {
    nimSmartPtr<EIGEN_EIGENCLASS> returnClass = new EIGEN_EIGENCLASS;
    EIGEN_EIGEN_INTERNAL(x, symmetric, valuesOnly, returnClass.getRealPtr());
    return(returnClass);
}

template<class Derived>
nimSmartPtr<EIGEN_SVDCLASS>   EIGEN_SVD(const Eigen::MatrixBase<Derived> &x, int vectors) {
    nimSmartPtr<EIGEN_SVDCLASS> returnClass = new EIGEN_SVDCLASS;
    EIGEN_SVD_INTERNAL(x, vectors, returnClass.getRealPtr());
    return(returnClass);
}


template<class Derived>
nimSmartPtr<EIGEN_EIGENCLASS_R>   EIGEN_EIGEN_R(const Eigen::MatrixBase<Derived> &x,  bool symmetric, bool valuesOnly) {
    nimSmartPtr<EIGEN_EIGENCLASS_R> returnClass = new EIGEN_EIGENCLASS_R;
    EIGEN_EIGEN_INTERNAL(x, symmetric, valuesOnly, returnClass.getRealPtr());
    return(returnClass);
}

template<class Derived>
nimSmartPtr<EIGEN_SVDCLASS_R>   EIGEN_SVD_R(const Eigen::MatrixBase<Derived> &x, int vectors) {
    nimSmartPtr<EIGEN_SVDCLASS_R> returnClass = new EIGEN_SVDCLASS_R;
    EIGEN_SVD_INTERNAL(x, vectors, returnClass.getRealPtr());
    return(returnClass);
}

template<class Derived>
bool EIGEN_CHECKSYMMETRY(const Eigen::MatrixBase<Derived> &x) {
	for(int i = 0; i<x.rows(); i++){
		for(int j = i + 1; j<x.rows(); j++){
			if(x(i,j) != x(j,i)){
				return(false);
			}
		}
	}
	return(true);
}


template<class Derived>
void EIGEN_EIGEN_INTERNAL(const Eigen::MatrixBase<Derived> &x,  bool symmetric, bool valuesOnly, EIGEN_EIGENCLASS_R *returnClass) {
    returnClass->getValues().initialize(0, 0, x.rows());
	Map<VectorXd> Eig_eigVals(returnClass->getValues().getPtr(),x.rows());
	if(!valuesOnly){
		returnClass->getVectors().initialize(0, 0, x.rows(), x.cols());
	}
	if(x.rows() != x.cols()) {
	 _nimble_global_output <<"Run-time size error: expected matrix argument to eigen() to be square."<<"\n"; nimble_print_to_R(_nimble_global_output);
	}
    Eigen::DecompositionOptions eigOpts = valuesOnly ? EigenvaluesOnly : ComputeEigenvectors;
	if(!symmetric){
		symmetric = EIGEN_CHECKSYMMETRY(x);
	}
	if(symmetric){ 
		SelfAdjointEigenSolver<MatrixXd> solver1(x, eigOpts); // The MatrixXd here doesn't seem generic, but I couldn't get it to work otherwise and it would be odd to do an Eigen decomposition on anything else. -Perry
		Eig_eigVals = solver1.eigenvalues().reverse();
		if(!valuesOnly){
			Map<MatrixXd> Eig_eigVecs(returnClass->getVectors().getPtr(),x.rows(),x.cols());
			Eig_eigVecs = solver1.eigenvectors().rowwise().reverse();	
		}
	}
	else{ 
		EigenSolver<MatrixXd> solver2(x, eigOpts); 
		vector<pair<double,int> > sortIndices;
		for(int i = 0; i<x.rows(); i++){
			sortIndices.push_back(make_pair(abs(solver2.eigenvalues().real()(i)),i));
		}
		std::sort(sortIndices.begin(),sortIndices.end());
		std::reverse(sortIndices.begin(), sortIndices.end());
		for(int i = 0; i < x.rows() ; ++i){
			if(solver2.eigenvalues().imag()(sortIndices[i].second) != 0){
				_nimble_global_output <<"Run-time warning: matrix used in call to nimEigen() has a complex eigenvalue."<<"\n"; nimble_print_to_R(_nimble_global_output);
				Eig_eigVals(i) = NAN;
			}
			else{
				Eig_eigVals(i) = solver2.eigenvalues().real()(sortIndices[i].second);
			}
		}
		if(!valuesOnly){
			Map<MatrixXd> Eig_eigVecs(returnClass->getVectors().getPtr(),x.rows(),x.cols());
			MatrixXd sorted_eigVecs(x.rows(), x.rows());
			for(int i = 0; i<x.rows(); i++){
				sorted_eigVecs.col(i) = solver2.eigenvectors().real().col(sortIndices[i].second);
				for(int j = 0; j<x.rows(); j++){
					if(solver2.eigenvectors().imag()(j, sortIndices[i].second) != 0){
						_nimble_global_output <<"Run-time warning: matrix matrix used in call to nimEigen() has a complex valued eigenvector."<<"\n"; nimble_print_to_R(_nimble_global_output);
						sorted_eigVecs(j, i) = NAN;
					}
				}
			}
			Eig_eigVecs = sorted_eigVecs;
		}
	} 
}

 template<class Derived> 
void EIGEN_SVD_INTERNAL(const Eigen::MatrixBase<Derived> &x, int vectors, EIGEN_SVDCLASS_R *returnClass) { 
    // nimSmartPtr<EIGEN_SVDCLASS> returnClass = new EIGEN_SVDCLASS; 
 	int nu = min(x.rows(), x.cols()); 
 	(*returnClass).d.initialize(0, 0, nu); 
 	Map<VectorXd> Svd_d((*returnClass).d.getPtr(), nu); 
 	JacobiSVD<MatrixXd> svd; 
 	/* note: if nu > 16, bidiagonialization algo. is recommended on eigen website.  not currently available w/ nimble's version of eigen, but may be in future. */
 	if(vectors == 0){ 
	  svd.compute(x); 
 	} 
 	else{ 
 	  int leftSVs = nu; 
 	  int rightSVs = nu; 
 	  if(vectors == 1){ 
 	  	svd.compute(x, ComputeThinU | ComputeThinV); 
 	  } 
 	  if(vectors == 2){ 
 	    leftSVs = x.rows(); 
 		rightSVs = x.cols(); 
 		svd.compute(x, ComputeFullU | ComputeFullV); 
 	  } 
 	  (*returnClass).u.initialize(0, 0, x.rows(), leftSVs); 
 	  (*returnClass).v.initialize(0, 0, x.cols(), rightSVs); 
 	  Map<MatrixXd> Svd_u((*returnClass).u.getPtr(), x.rows(), leftSVs); 
 	  Map<MatrixXd> Svd_v((*returnClass).v.getPtr(), x.cols(), rightSVs); 
 	  Svd_u = svd.matrixU(); 
 	  Svd_v = svd.matrixV(); 
 	} 
	Svd_d = svd.singularValues();  
 	//return(returnClass); 
 }


template<class derived1, class derived2>
MatrixXd EIGEN_FS(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd ans = x.template triangularView<Eigen::Lower>().solve(ycopy);
  return(ans);
}

template<class derived1, class derived2>
MatrixXd_CppAD nimDerivs_EIGEN_FS_no_atomic(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd_CppAD ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd_CppAD ans = x.template triangularView<Eigen::Lower>().solve(ycopy);
  return(ans);
}

template<class derived1, class derived2>
MatrixXd EIGEN_BS(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd ans = x.template triangularView<Eigen::Upper>().solve(ycopy);
  return(ans);
}

template<class derived1, class derived2>
MatrixXd_CppAD nimDerivs_EIGEN_BS_no_atomic(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd_CppAD ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd_CppAD ans = x.template triangularView<Eigen::Upper>().solve(ycopy);
  return(ans);
}

template<class derived1, class derived2>
MatrixXd EIGEN_SOLVE(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd ans = x.lu().solve(ycopy);
  return(ans);
}

template<class derived1, class derived2>
MatrixXd_CppAD nimDerivs_EIGEN_SOLVE(const MatrixBase<derived1> &x, const MatrixBase<derived2> &y) {
  MatrixXd_CppAD ycopy = y; // in case y was a map, which it will always be from nimble
  MatrixXd_CppAD ans = x.lu().solve(ycopy);
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const ArrayBase<Derived1>& v1, const ArrayBase<Derived2>& v2) { 
  double ans = (v1 * v2).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
typename Derived2::Scalar nimDerivs_eigenInprod(const ArrayBase<Derived1>& v1, const ArrayBase<Derived2>& v2) { 
  typename Derived2::Scalar ans = (v1 * v2).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const MatrixBase<Derived1>& v1, const MatrixBase<Derived2>& v2) {
  double ans = (v1.cwiseProduct(v2)).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
typename Derived2::Scalar nimDerivs_eigenInprod(const MatrixBase<Derived1>& v1, const MatrixBase<Derived2>& v2) {
  typename Derived2::Scalar ans = (v1.cwiseProduct(v2)).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const MatrixBase<Derived1>& v1, const ArrayBase<Derived2>& v2) {
  double ans = (v1.cwiseProduct(v2.matrix())).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
typename Derived2::Scalar nimDerivs_eigenInprod(const MatrixBase<Derived1>& v1, const ArrayBase<Derived2>& v2) {
  typename Derived2::Scalar ans = (v1.cwiseProduct(v2.matrix())).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
double eigenInprod(const ArrayBase<Derived1>& v1, const MatrixBase<Derived2>& v2) {
  double ans = (v2.cwiseProduct(v1.matrix())).sum();
  return(ans);
}

template <typename Derived1, typename Derived2>
  typename Derived2::Scalar nimDerivs_eigenInprod(const ArrayBase<Derived1>& v1, const MatrixBase<Derived2>& v2) {
  typename Derived2::Scalar ans = (v2.cwiseProduct(v1.matrix())).sum();
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

template<class Derived>
typename Derived::Scalar nimDerivs_var(const ArrayBase<Derived>& v) {
  typedef typename Derived::Scalar T;
  T mean = v.mean();
  int n = v.size();
  T ans = (v-mean).matrix().squaredNorm() / T(n-1);
  return ans;
}

template<class Derived>
typename Derived::Scalar nimDerivs_sd(const ArrayBase<Derived>& v) {
  return CppAD::sqrt(nimDerivs_var(v));
}

template <typename Derived>
typename Derived::Scalar det(const MatrixBase<Derived>& v) {
  return(v.determinant());
}

template <typename Derived>
typename Derived::Scalar logdet(const MatrixBase<Derived>& v) {
  return(log(v.determinant()));
}

#endif
