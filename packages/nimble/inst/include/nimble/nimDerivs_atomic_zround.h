#ifndef _NIMDERIVS_ATOMIC_ZROUND
#define _NIMDERIVS_ATOMIC_ZROUND

#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimbleCppAD.h"
#include "nimDerivs_atomic_discrete.h"

CppAD::AD<double> nimDerivs_zround(const CppAD::AD<double> x);

atomic_zround_class *track_atomic_zround(void* tape_mgr_ptr,
					 std::vector<CppAD::local::atomic_index_info>* vec_ptr);

class atomic_zround_class : public CppAD::atomic_three<double>, public nimble_atomic_base {
public:
  atomic_zround_class(const std::string& name);
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

/*************************/

CppAD::AD<double> nimDerivs_floor(const CppAD::AD<double> x);

atomic_floor_class *track_atomic_floor(void* tape_mgr_ptr,
				       std::vector<CppAD::local::atomic_index_info>* vec_ptr);

class nimDerivs_floor_class {
 public:
  double operator()(double x);
  double operator()(int x);
  CppAD::AD<double> operator()(const CppAD::AD<double> &x);
};

class atomic_floor_class : public atomic_discrete_class<nimDerivs_floor_class> {
public:
  atomic_floor_class(const std::string& name);
};

/***********/

CppAD::AD<double> nimDerivs_ceil(const CppAD::AD<double> x);

atomic_ceil_class *track_atomic_ceil(void* tape_mgr_ptr,
				       std::vector<CppAD::local::atomic_index_info>* vec_ptr);

class nimDerivs_ceil_class {
 public:
  double operator()(double x);
  double operator()(int x);
  CppAD::AD<double> operator()(const CppAD::AD<double> &x);
};

class atomic_ceil_class : public atomic_discrete_class<nimDerivs_ceil_class> {
public:
  atomic_ceil_class(const std::string& name);
};

/***********/

CppAD::AD<double> nimDerivs_ftrunc(const CppAD::AD<double> x);

atomic_ftrunc_class *track_atomic_ftrunc(void* tape_mgr_ptr,
				       std::vector<CppAD::local::atomic_index_info>* vec_ptr);

class nimDerivs_ftrunc_class {
 public:
  double operator()(double x);
  double operator()(int x);
  CppAD::AD<double> operator()(const CppAD::AD<double> &x);
};

class atomic_ftrunc_class : public atomic_discrete_class<nimDerivs_ftrunc_class> {
public:
  atomic_ftrunc_class(const std::string& name);
};

/***********/

CppAD::AD<double> nimDerivs_nimRound(const CppAD::AD<double> x);

atomic_nimRound_class *track_atomic_nimRound(void* tape_mgr_ptr,
				       std::vector<CppAD::local::atomic_index_info>* vec_ptr);

class nimDerivs_nimRound_class {
 public:
  double operator()(double x);
  double operator()(int x);
  CppAD::AD<double> operator()(const CppAD::AD<double> &x);
};

class atomic_nimRound_class : public atomic_discrete_class<nimDerivs_nimRound_class> {
public:
  atomic_nimRound_class(const std::string& name);
};

/***********/

CppAD::AD<double> nimDerivs_nimStep(const CppAD::AD<double> x);

atomic_nimStep_class *track_atomic_nimStep(void* tape_mgr_ptr,
				       std::vector<CppAD::local::atomic_index_info>* vec_ptr);

class nimDerivs_nimStep_class {
 public:
  double operator()(double x);
  double operator()(int x);
  CppAD::AD<double> operator()(const CppAD::AD<double> &x);
};

class atomic_nimStep_class : public atomic_discrete_class<nimDerivs_nimStep_class> {
public:
  atomic_nimStep_class(const std::string& name);
};

#define TTT_ CppAD::AD<double>

inline TTT_ nimDerivs_pairmax(TTT_ x1, TTT_ x2) {
  TTT_ cond = nimDerivs_nimStep(x1-x2);
  return CppAD::azmul(cond, x1) + CppAD::azmul(TTT_(1) - cond, x2);
//  return(CondExpGt(x1, x2, x1, x2));
}

inline TTT_ nimDerivs_pairmin(TTT_ x1, TTT_ x2) {
  TTT_ cond = nimDerivs_nimStep(x1-x2);
  return CppAD::azmul((TTT_(1)-cond), x1) + CppAD::azmul(cond, x2);
//  return(CondExpLt(x1, x2, x1, x2));
}

inline TTT_ nimDerivs_nimEquals(TTT_ x1, TTT_ x2){
  TTT_ cond = nimDerivs_nimStep(x1-x2) * nimDerivs_nimStep(x2-x1);
  return cond;
//  return(CondExpEq(x1, x2, T(1), T(0)));
}

inline TTT_ nimDerivs_CondExpEq(TTT_ x1, TTT_ x2, TTT_ y1, TTT_ y2) {
  TTT_ cond = nimDerivs_nimEquals(x1, x2);
  return CppAD::azmul(cond, y1) + CppAD::azmul((TTT_(1) - cond), y2);
//  return(CondExpEq(x1, x2, T(1), T(0)));
}

inline TTT_ nimDerivs_CondExpGe(TTT_ x1, TTT_ x2, TTT_ y1, TTT_ y2) {
  TTT_ cond = nimDerivs_nimStep(x1-x2);
  return CppAD::azmul(cond, y1) + CppAD::azmul((TTT_(1) - cond), y2);
}

inline TTT_ nimDerivs_CondExpGt(TTT_ x1, TTT_ x2, TTT_ y1, TTT_ y2) {
  TTT_ cond = 1-nimDerivs_nimStep(x2-x1);
  return CppAD::azmul(cond, y1) + CppAD::azmul((TTT_(1) - cond), y2);
}

inline TTT_ nimDerivs_CondExpLe(TTT_ x1, TTT_ x2, TTT_ y1, TTT_ y2) {
  TTT_ cond = nimDerivs_nimStep(x2-x1);
  return CppAD::azmul(cond, y1) + CppAD::azmul((TTT_(1) - cond), y2);
}

inline TTT_ nimDerivs_CondExpLt(TTT_ x1, TTT_ x2, TTT_ y1, TTT_ y2) {
  TTT_ cond = 1-nimDerivs_nimStep(x1-x2);
  return CppAD::azmul(cond, y1) + CppAD::azmul((TTT_(1) - cond), y2);
}

#undef TTT_

#endif

