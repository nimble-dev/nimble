#ifndef _NIMDERIVS_ATOMIC_PROBIT
#define _NIMDERIVS_ATOMIC_PROBIT

#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimbleCppAD.h"

CppAD::AD<double> nimDerivs_probit(const CppAD::AD<double> x);
CppAD::AD<double> nimDerivs_iprobit(const CppAD::AD<double> x);

atomic_probit_class *track_atomic_probit(void* tape_mgr_ptr,
					 std::vector<CppAD::local::atomic_index_info>* vec_ptr);
atomic_iprobit_class *track_atomic_iprobit(void* tape_mgr_ptr,
					 std::vector<CppAD::local::atomic_index_info>* vec_ptr);

/******************************/

class atomic_iprobit_class :public unary_atomic_class<double> {
 public:
  // iprobit(x) = probit^{-1}(x) = pnorm(x, 0, 1)
  // From atomic_three_get_started
  atomic_iprobit_class(const std::string& name);

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
		       const CppAD::vector< CppAD::AD<double> >&  parameter_x  ,
		       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
		       size_t                                     need_y       ,
		       size_t                                     order_low    ,
		       size_t                                     order_up     ,
		       const CppAD::vector< CppAD::AD<double> >& taylor_x     ,
		       CppAD::vector< CppAD::AD<double> >&       taylor_y     );

  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   );
  
  virtual bool reverse(
      const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
      CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
      const CppAD::vector< CppAD::AD<double> >&               partial_y   );
};

 
/*******************************************/

class atomic_probit_class : public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  // probit(x) is qnorm(x, 0, 1)
  atomic_probit_class(const std::string& name);
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
		       const CppAD::vector< CppAD::AD<double> >&  parameter_x  ,
		       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
		       size_t                                     need_y       ,
		       size_t                                     order_low    ,
		       size_t                                     order_up     ,
		       const CppAD::vector< CppAD::AD<double> >& taylor_x     ,
		       CppAD::vector< CppAD::AD<double> >&       taylor_y     );
  
  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   );

    virtual bool reverse(
      const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
      CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
      const CppAD::vector< CppAD::AD<double> >&               partial_y   );


};

#endif
