#ifndef _NIMDERIVS_ATOMIC_POW_INT
#define _NIMDERIVS_ATOMIC_POW_INT

// See end of this file for global object definitions and functions to call
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimbleCppAD.h"

CppAD::AD<double> nimDerivs_pow_int(const CppAD::AD<double> &a,
				    const CppAD::AD<double> &b);

atomic_pow_int_class *track_atomic_pow_int(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr);

class atomic_pow_int_class : public CppAD::atomic_three<double>, public nimble_atomic_base {
public:
  atomic_pow_int_class(const std::string& name);
 private:
  bool for_type(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      CppAD::vector<CppAD::ad_type_enum>&        type_y      );
  // Not sure this is ever needed.
  bool for_type(
		const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
		const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
		CppAD::vector<CppAD::ad_type_enum>&        type_y      );
  // rev_depend is used when the optimize() method is called for a tape (an ADFun).
  bool rev_depend(
		  const CppAD::vector<double>&          parameter_x ,
		  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
		  CppAD::vector<bool>&                depend_x    ,
		  const CppAD::vector<bool>&          depend_y
		  );
  bool forward(
	       const CppAD::vector<double>&               parameter_x  ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
	       size_t                              need_y       ,
	       size_t                              order_low    ,
	       size_t                              order_up     ,
	       const CppAD::vector<double>&               taylor_x     ,
	       CppAD::vector<double>&                     taylor_y     );
  bool forward(
	       const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
	       size_t                              need_y       ,
	       size_t                              order_low    ,
	       size_t                              order_up     ,
	       const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
	       CppAD::vector< CppAD::AD<double> >&                     taylor_y     );
  bool reverse(
	       const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
	       size_t                              order_up    ,
	       const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
	       const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
	       CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
	       const CppAD::vector< CppAD::AD<double> >&               partial_y   );
  bool reverse(
	       const CppAD::vector<double>&               parameter_x ,
	       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
	       size_t                              order_up    ,
	       const CppAD::vector<double>&               taylor_x    ,
	       const CppAD::vector<double>&               taylor_y    ,
	       CppAD::vector<double>&                     partial_x   ,
	       const CppAD::vector<double>&               partial_y   );
};

#endif
