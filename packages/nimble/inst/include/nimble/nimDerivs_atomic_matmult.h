#ifndef _NIMDERIVS_ATOMIC_MATMULT
#define _NIMDERIVS_ATOMIC_MATMULT


// See end of this file for global object definitions and functions to call
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimDerivs_vecmat_utils.h"

/*
Atomic class for matrix multiplication.
This benefits from TMB's implementation but is updated to CppAD's atomic_three system.
We have Y = X1 %*% X2.
X1 is n1-x-n2.
X2 is n2-x-n3.
Y is n1-x-n3
The input X is c(n2, X1, X2), giving the view that
 Y = f(X), where f: R^n -> R^m, 
 where m = n1*n3 and n = 1 + n1*n2 + n2*n3 = 1 + n2*(n1 + n3)
(The "1+" is for the element whose values will be n1, which is placed first.  Y does not actually depend on this element.
(The m and n are CppAD's notation.)
We need one (not two) of the three sizes as an input to deduce all dimensions.
We choose to have that be n1, for convenience.
We can then deduce the dimensions from inputs to forward() and reverse() as follows.
nrow = (q+1) = (order_up+1), where q = order_up are CppAD notations.
taylor_y is nrow-x-m, so m = taylor_y.size()/nrow.
taylor_x is nrow-x-n, so n = taylor_x.size()/nrow.
n3 = m/n1
n2 = (n-1) / (n1 + n3)

(Another strategy would be to have n2 be a class member and have different 
objects for each usage.  However it is not clear how we would manage those objects.
Hence we stick with a single global (static) object and include n2 in the
input (X) vector.)

-----
Value:
Y = X1 %*% X2
(More fully: Y(t) = H(X(t)) = G(F(X(t)))
-----
Forward first order
dY  = dX1 %*% X2 + X1 %*% dX2

----
Reverse first order:
According to Giles(2008):
X1adjoint = Yadjoint %*% t(X2)
X2adjoint = t(X1) %*% Yadjoint

These derive from
dS = <X1adjoint, dX1> + <X2adjoint, dX2> = <Yadjoint, dY> = <Yadjoint, dX1 %*% X2 + X1 %*% dX2> = <Yadjoint,  dX1 %*% X2> + <Yadjoint,  X1 %*% dX2>
   = <Yadjoint %*% t(X2), dX1> + <X1^T %*% Yadjoint, dX2>
Rules
<A, B> = tr(A^T B) = tr(B A^T) = <B, A>
<A, BC> = tr(A^T BC) = tr(C A^T B) = <A C^T, B>
<A, BC> = tr(A^T BC) = tr(A^TB C) = <B^T A, C>
<A, BCD> = <A D^T, BC> = <B^T A D^T, C>
----
Reverse second order:
Fdot =  X1dot %*% X2 + X1 %*% X2dot
dYdot = dFdot = dX1dot %*% X2 + X1dot %*% dX2 + dX1 %*% X2dot + X1 %*% dX2dot
S0 = H = G(F(X1, X2), Fdot(X1, X2, X1dot, X2dot)) = G( X1 %*% X2, X1dot %*% X2 + X1 %*% X2dot ) ## F() is same as Y()
X1adjoint = Yadjoint %*% t(X2) + Ydot_adjoint %*% t(X2dot)
X2adjoint = t(X1) %*% Yadjoint + t(X1dot) %*% Ydot_adjoint
X1dot_adjoint = Ydot_adjoint %*% t(X2)
X2dot_adjoint = t(X1) %*% Ydot_adjoint

These derive from:
dS = <Yadjoint, dY> + <Ydot_adjoint, dYdot>
X1adjoint: first term is the same as from reverse first order
second term: <Ydot_adjoint, dX1dot %*% X2 + X1dot %*% dX2 + dX1 %*% X2dot + X1 %*% dX2dot>
            = <Ydot_adjoint,  dX1dot %*% X2> + <Ydot_adjoint, X1dot %*% dX2> +  <Ydot_adjoint, dX1 %*% X2dot> + <Ydot_adjoint, X1 %*% dX2dot>
            = <Ydot_adjoint %*% X2^T, dX1dot> + <X1dot^T %*% Ydot_adjoint, dX2> + <Ydot_adjoint %*% X2dot^T, dX1> + <X1^T %*% Ydot_adjoint, dX2dot>
            = <X1dot_adjoint term, dX1dot> +         <X2adjoint_term, dX2> +...

*/

void atomic_matmult(const MatrixXd_CppAD &x1,
		    const MatrixXd_CppAD &x2,
		    MatrixXd_CppAD &y);

#define USE_NEW_DYNAMIC // Now required, for the matrix_category handling

enum matrix_category {lower_diagonal, upper_diagonal, square_full, non_square, unknown};

MatrixXd_CppAD nimDerivs_matmult(const MatrixXd_CppAD &x1,
				 const MatrixXd_CppAD &x2);

class atomic_matmult_class :  public CppAD::atomic_three< double > {
 private:
  int n1_;
  std::vector<double> X_stored;
  std::vector<CppAD::AD<double> > X_AD_stored;
  matrix_category X1cat_, X2cat_;
  bool x1_is_constant_, x2_is_constant_;
  
 public:
  double * get_X_stored_ptr() {return &X_stored[0];}
  CppAD::AD<double> * get_X_AD_stored_ptr() {return &X_AD_stored[0];}
  matrix_category transpose_matrix_cat(const matrix_category &Xcat) const {
    if(Xcat == lower_diagonal) return upper_diagonal;
    if(Xcat == upper_diagonal) return lower_diagonal;
    return Xcat;
  };
  void set_X_stored(const MatrixXd_CppAD &X);
  matrix_category &X1cat() {return X1cat_;}
  matrix_category &X2cat() {return X2cat_;}
  matrix_category const &X1cat() const {return X1cat_;}
  matrix_category const &X2cat() const {return X2cat_;}
  matrix_category transpose_X1cat() {return transpose_matrix_cat(X1cat_);}
  matrix_category transpose_X2cat() {return transpose_matrix_cat(X2cat_);}
  matrix_category transpose_X1cat() const {return transpose_matrix_cat(X1cat_);}
  matrix_category transpose_X2cat() const {return transpose_matrix_cat(X2cat_);}

  bool &X1constant() {return x1_is_constant_;}
  bool &X2constant() {return x2_is_constant_;}
  bool const &X1constant() const {return x1_is_constant_;}
  bool const &X2constant() const {return x2_is_constant_;}
  
  atomic_matmult_class(const std::string& name);
  void set_n1(int n1__) {n1_ = n1__;}
  int get_n1() {return n1_;}

#ifdef USE_NEW_DYNAMIC
 private:
  // suffix C for const
  EigenTemplateTypes<double>::typeEigenConstMapStrd X1mapC, X2mapC, dX1mapC, dX2mapC;
  EigenTemplateTypes<double>::typeEigenConstMapStrd Yadjoint_mapC, X1dot_mapC, X2dot_mapC, Ydot_adjoint_mapC;
  
  EigenTemplateTypes<double>::typeEigenMapStrd Ymap, dY_map;
  EigenTemplateTypes<double>::typeEigenMapStrd X1adjoint_map, X2adjoint_map, X1dot_adjoint_map,  X2dot_adjoint_map;

  // prefix m for meta
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mX1mapC, mX2mapC, mdX1mapC, mdX2mapC;
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd mYadjoint_mapC, mX1dot_mapC, mX2dot_mapC, mYdot_adjoint_mapC;

  EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd mYmap, mdY_map;
  EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd mX1adjoint_map, mX2adjoint_map, mX1dot_adjoint_map,  mX2dot_adjoint_map;

  EigenTemplateTypes<CppAD::AD<double>>::typeMatrixXd mTerm1, mTerm2;
#endif

  
 private:
 // parameter_x is a vector of: value for constants; value for dynamics (from last new_dynamic); unspecified for variables
 // type_x is a vector of: CppAD::constant_enum, CppAD::dynamic_enum, or CppAD::variable_enum
 // type_y is return vector of same options as type_x
 virtual bool for_type(
		       const CppAD::vector<double>&               parameter_x ,
		       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
		       CppAD::vector<CppAD::ad_type_enum>&        type_y      );

 virtual bool rev_depend(
			 const CppAD::vector<double>&          parameter_x ,
			 const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
			 CppAD::vector<bool>&                depend_x    ,
			 const CppAD::vector<bool>&          depend_y
			 );
 
 virtual bool forward(
		      const CppAD::vector<double>&               parameter_x  ,
		      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
		      size_t                              need_y       ,
		      size_t                              order_low    ,
		      size_t                              order_up     ,
		      const CppAD::vector<double>&               taylor_x     ,
		      CppAD::vector<double>&                     taylor_y     );
 
 virtual bool forward(
		      const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
		      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
		      size_t                              need_y       ,
		      size_t                              order_low    ,
		      size_t                              order_up     ,
		      const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
		      CppAD::vector<CppAD::AD<double> >&                     taylor_y     );
 
 virtual bool reverse(
		      const CppAD::vector<double>&               parameter_x ,
		      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
		      size_t                              order_up    ,
		      const CppAD::vector<double>&               taylor_x    ,
		      const CppAD::vector<double>&               taylor_y    ,
		      CppAD::vector<double>&                     partial_x   ,
		      const CppAD::vector<double>&               partial_y   );

 virtual bool reverse(
		      const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
		      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
		      size_t                              order_up    ,
		      const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
		      const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
		      CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
		      const CppAD::vector<CppAD::AD<double> >&               partial_y   );
};

#endif
