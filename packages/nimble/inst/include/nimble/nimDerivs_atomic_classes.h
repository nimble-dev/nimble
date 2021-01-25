#ifndef _NIMDERIVS_ATOMIC_CLASSES
#define _NIMDERIVS_ATOMIC_CLASSES

// See end of this file for global object definitions and functions to call
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "nimDerivs_atomic_matmult.h"
#include "nimDerivs_atomic_matinverse.h"
#include "nimDerivs_atomic_backsolve.h"
#include "nimDerivs_atomic_forwardsolve.h"
#include "nimDerivs_atomic_cholesky.h"

template<class T>
class unary_atomic_class : public CppAD::atomic_three<T> {
  // This layer in the class hierarchy simply provides the same for_type and rev_depend
  // for all atomic classes representing unary functions below.
 public:
 unary_atomic_class(const std::string& name) : CppAD::atomic_three<T>(name) {};
 private:
   // for_type is essential.
  // This determines which elements of y are constants, dynamic parameters, or variables
  //   depending on types of x.
  // If omitted, calls to forward (and possibly reverse) will not happen.
  virtual bool for_type(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
  { 
    type_y[0] = type_x[0];
    return true;
  }
  // rev_depend is used when the optimize() method is called for a tape (an ADFun).
  virtual bool rev_depend(
			  const CppAD::vector<double>&          parameter_x ,
			  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
			  CppAD::vector<bool>&                depend_x    ,
			  const CppAD::vector<bool>&          depend_y
			  ) {
    depend_x[0] = depend_y[0];
    return true;
  }
};

class atomic_lgamma_class : public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  atomic_lgamma_class(const std::string& name, int baseOrder_);

private:
    int baseOrder;
  virtual bool forward(
      const CppAD::vector<double>&               parameter_x  ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
      size_t                              need_y       ,
      size_t                              order_low    ,
      size_t                              order_up     ,
      const CppAD::vector<double>&               taylor_x     ,
      CppAD::vector<double>&                     taylor_y     );
  virtual bool forward(
		       const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
		       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
		       size_t                              need_y       ,
		       size_t                              order_low    ,
		       size_t                              order_up     ,
		       const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
		       CppAD::vector< CppAD::AD<double> >&                     taylor_y     );
  virtual bool reverse(
      const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
      CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
      const CppAD::vector< CppAD::AD<double> >&               partial_y   );
  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   );
};

/******************************/

class atomic_pnorm1_class :public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  atomic_pnorm1_class(const std::string& name) : 
    unary_atomic_class<double>(name)        
    { }

private:
   virtual bool forward(
      const CppAD::vector<double>&               parameter_x  ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
      size_t                              need_y       ,
      size_t                              order_low    ,
      size_t                              order_up     ,
      const CppAD::vector<double>&               taylor_x     ,
      CppAD::vector<double>&                     taylor_y     )
  {
    // This can do more with checking sizes, need_y, orders, and
    // what is a parameter.
    // need_y determines whether all elements need calculation,
    // or only those for parameters, or only those for variables.
    if(order_low == 0) {
      taylor_y[0] = Rf_pnorm5(taylor_x[0], 0, 1, 1, 0); // (x, mu, sigma, lower_tail, log)
      // F(x)
    }
    double Fprime(0.);
    if(order_low <= 1 && order_up >= 1) {
      Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
      taylor_y[1] = Fprime * taylor_x[1];
      // F'(x) (x')
    }
    if(order_low <= 2 && order_up >= 2) {
      if(Fprime == 0) {     // only calculate Fprime if it is really needed
	if(taylor_x[2] != 0)
	  Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
      }
      taylor_y[2] = 0.5 * (- taylor_x[0]) * taylor_y[1] * taylor_x[1] +
	Fprime*taylor_x[2]; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
    }
    return true;
  }

  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   )
  {
    double Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
    partial_x[0] = 0;
    if(order_up >= 1) partial_x[1] = 0;
    if(order_up >= 2) return false; // not implemented
    if(order_up >= 1) { // Needed for 2nd derivs
      partial_x[0] += partial_y[1] * (- taylor_x[0]) * Fprime * taylor_x[1];
      partial_x[1] += partial_y[1] * Fprime;
    } 
    // Needed for 1st derivs
    partial_x[0] += partial_y[0] * Fprime;
    return true;
  }
};

/*******************************************/

class atomic_qnorm1_class : public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  atomic_qnorm1_class(const std::string& name) : 
    unary_atomic_class<double>(name)        
    { }

private:
  virtual bool forward(
      const CppAD::vector<double>&               parameter_x  ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
      size_t                              need_y       ,
      size_t                              order_low    ,
      size_t                              order_up     ,
      const CppAD::vector<double>&               taylor_x     ,
      CppAD::vector<double>&                     taylor_y     )
  {
    // This can do more with checking sizes, need_y, orders, and
    // what is a parameter.
    // need_y determines whether all elements need calculation,
    // or only those for parameters, or only those for variables.
    if(order_low == 0) {
      taylor_y[0] = Rf_qnorm5(taylor_x[0], 0, 1, 1, 0); // (p, mu, sigma, lower_tail, log)
      // F(x)
    }
    double z = taylor_y[0];
    double invFprime(0.);
    if(order_low <= 1 && order_up >= 1) {
      invFprime = Rf_dnorm4(z, 0, 1, 0);
      taylor_y[1] = taylor_x[1] / invFprime;
      // F'(x) (x')
    }
    if(order_low <= 2 && order_up >= 2) {
      taylor_y[2] = 0.5 * z * taylor_y[1] * taylor_y[1];
      if(taylor_x[2] != 0) {
        if(invFprime == 0)
          invFprime = Rf_dnorm4(z, 0, 1, 0);
        taylor_y[2] += taylor_x[2] / invFprime; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
      }
    }
    return true;
  }

  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   )
  {
    double z = taylor_y[0];
    double invFprime = Rf_dnorm4(z, 0, 1, 0);
    partial_x[0] = 0;
    if(order_up >= 1) partial_x[1] = 0;
    if(order_up >= 2) return false; // not implemented
    if(order_up >= 1) { // Needed for 2nd derivs
      partial_x[0] += partial_y[1] * (z/(invFprime*invFprime)) * taylor_x[1];
      partial_x[1] += partial_y[1] / invFprime;
    } 
    // Needed for 1st derivs
    partial_x[0] += partial_y[0] / invFprime;
    return true;
  }
};

template<class T>
T nimDerivs_lgammafn(T x, int baseOrder) {
  static atomic_lgamma_class static_atomic_lgamma0("nimDerivs_lgamma", 0);
  static atomic_lgamma_class static_atomic_lgamma1("nimDerivs_lgamma", 1);
  static atomic_lgamma_class static_atomic_lgamma2("nimDerivs_lgamma", 2);
  static atomic_lgamma_class static_atomic_lgamma3("nimDerivs_lgamma", 3);
  static atomic_lgamma_class static_atomic_lgamma4("nimDerivs_lgamma", 4);
  CppAD::vector<T> in(1);
  CppAD::vector<T> out(1);
  in[0] = x;
  switch(baseOrder) {
  case 0:
    static_atomic_lgamma0(in, out);
    break;
  case 1:
    static_atomic_lgamma1(in, out);
    break;
  case 2:
    static_atomic_lgamma2(in, out);
    break;
  case 3:
    static_atomic_lgamma3(in, out);
    break;
  case 4:
    static_atomic_lgamma4(in, out);
    break;
  default:
    std::cout<<"Error: attempting lgamma derivative beyond order 4."<<std::endl;
  }
  return out[0];
}

template<class T>
T nimDerivs_lgammafn(T x) {
  return nimDerivs_lgammafn<T>(x, 0); // simply putting a default value of 0 above leads to failure in std::ptr_fun<CppAD::AD<double>, CppAD::AD<double>>(nimDerivs_lgammafn) in Eigen-style vectorization
}

template<class T>
T nimDerivs_pnorm1(T x) {
  static atomic_pnorm1_class static_atomic_pnorm1("nimDerivs_pnorm1");
  CppAD::vector<T> in(1);
  CppAD::vector<T> out(1);
  in[0] = x;
  static_atomic_pnorm1(in, out);
  return out[0];
}

template<class T>
T nimDerivs_qnorm1(T x) {
  static atomic_qnorm1_class static_atomic_qnorm1("nimDerivs_qnorm1");
  CppAD::vector<T> in(1);
  CppAD::vector<T> out(1);
  in[0] = x;
  static_atomic_qnorm1(in, out);
  return out[0];
}

#endif
