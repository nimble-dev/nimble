#ifndef _NIMDERIVS_ATOMIC_CLASSES
#define _NIMDERIVS_ATOMIC_CLASSES

// See end of this file for global object definitions and functions to call
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "nimbleCppAD.h"
#include "nimDerivs_atomic_pow_int.h"
#include "nimDerivs_atomic_zround.h"
#include "nimDerivs_atomic_log_pow_int.h"
#include "nimDerivs_atomic_matmult.h"
#include "nimDerivs_atomic_matinverse.h"
#include "nimDerivs_atomic_backsolve.h"
#include "nimDerivs_atomic_forwardsolve.h"
#include "nimDerivs_atomic_cholesky.h"
#include "nimDerivs_atomic_cache.h"
#include "nimDerivs_atomic_probit.h"

CppAD::AD<double> nimDerivs_lgammafn_base(CppAD::AD<double> x, int baseOrder, bool verbose = FALSE);
CppAD::AD<double> nimDerivs_lgammafn(CppAD::AD<double> x);

class atomic_lgamma_class : public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  atomic_lgamma_class(const std::string& name, int baseOrder_);
  atomic_lgamma_class(const std::string& name, int baseOrder_, bool verbose_);
  ~atomic_lgamma_class() {}; //std::cout<<"destructing atomic_lgamma_class"<<std::endl;
private:
    int baseOrder;
    bool verbose;
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

template<class T>
T nimDerivs_lgammafn_verbose(T x, int baseOrder ) {
  static atomic_lgamma_class static_atomic_lgamma0("nimDerivs_lgamma", 0, true);
  static atomic_lgamma_class static_atomic_lgamma1("nimDerivs_lgamma", 1, true);
  static atomic_lgamma_class static_atomic_lgamma2("nimDerivs_lgamma", 2, true);
  static atomic_lgamma_class static_atomic_lgamma3("nimDerivs_lgamma", 3, true);
  static atomic_lgamma_class static_atomic_lgamma4("nimDerivs_lgamma", 4, true);
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

atomic_lgamma_class* new_atomic_lgamma(void* tape_mgr, const std::string& name, int bO);
void delete_atomic_lgamma(void* tape_mgr, atomic_lgamma_class *atomic_lgamma);

template<class T>
T nimDerivs_lgammafn_verbose(T x) {
  return nimDerivs_lgammafn_verbose<T>(x, 0);
}

/***********************************************/

CppAD::AD<double> nimDerivs_gammafn(CppAD::AD<double> x);
CppAD::AD<double> nimDerivs_gammafn(CppAD::AD<double> x);

class atomic_gammafn_class : public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  atomic_gammafn_class(const std::string& name);
  ~atomic_gammafn_class() {}; //std::cout<<"destructing atomic_gammafn_class"<<std::endl;
private:
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

atomic_gammafn_class* new_atomic_gamma(void* tape_mgr, const std::string& name);
void delete_atomic_gamma(void* tape_mgr, atomic_gammafn_class *atomic_gamma);



#endif
