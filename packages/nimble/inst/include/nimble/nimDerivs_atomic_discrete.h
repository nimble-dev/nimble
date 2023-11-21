#ifndef _NIMDERIVS_ATOMIC_DISCRETE
#define _NIMDERIVS_ATOMIC_DISCRETE

#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimbleCppAD.h"

template<class ftor>
  class atomic_discrete_class : public unary_atomic_class<double> {
 public:
 atomic_discrete_class(const std::string& name) :
  unary_atomic_class<double>(name) {};
 private:
  ftor op;
  
  bool forward(
	       const CppAD::vector<double>&               parameter_x  ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
	       size_t                              need_y       ,
	       size_t                              order_low    ,
	       size_t                              order_up     ,
	       const CppAD::vector<double>&               taylor_x     ,
	       CppAD::vector<double>&                     taylor_y     ) {
    if((order_low <= 0) && (order_up >= 0)) {
      taylor_y[0] = op(taylor_x[0]);
    }
    if((order_low <= 1) && (order_up >= 1)) {
      taylor_y[1] = 0.;
    }
    return true;
  }
  bool forward(
	       const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
	       size_t                              need_y       ,
	       size_t                              order_low    ,
	       size_t                              order_up     ,
	       const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
	       CppAD::vector< CppAD::AD<double> >&                     taylor_y     ) {
    if((order_low <= 0) && (order_up >= 0)) {
      taylor_y[0] = op(taylor_x[0]);
    }
    if((order_low <= 1) && (order_up >= 1)) {
      taylor_y[1] = CppAD::AD<double>(0.);
    }
    return true;
  }
  bool reverse(
	       const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
	       size_t                              order_up    ,
	       const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
	       const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
	       CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
	       const CppAD::vector< CppAD::AD<double> >&               partial_y   ) {
    if(order_up >= 0) {
      partial_x[0] = CppAD::AD<double>(0.);
    }
    if(order_up >= 1) {
      partial_x[1] = CppAD::AD<double>(0.);
    }
    return true;
  }
  bool reverse(
	       const CppAD::vector<double>&               parameter_x ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
	       size_t                              order_up    ,
	       const CppAD::vector<double>&               taylor_x    ,
	       const CppAD::vector<double>&               taylor_y    ,
	       CppAD::vector<double>&                     partial_x   ,
	       const CppAD::vector<double>&               partial_y   ) {
    if(order_up >= 0) {
      partial_x[0] = 0.;
    }
    if(order_up >= 1) {
      partial_x[1] = 0.;
    }
    return true;
  }
};

#endif
