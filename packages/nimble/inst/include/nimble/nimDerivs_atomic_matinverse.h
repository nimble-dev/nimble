#ifndef _NIMDERIVS_ATOMIC_MATINVERSE
#define _NIMDERIVS_ATOMIC_MATINVERSE

#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimDerivs_vecmat_utils.h"
#include "nimDerivs_atomic_matmult.h"
#include "nimDerivs_atomic_cache.h"

void atomic_matinverse(const MatrixXd_CppAD &x,
		       MatrixXd_CppAD &y);

MatrixXd_CppAD nimDerivs_matinverse(const MatrixXd_CppAD &x);

class atomic_matinverse_class :  public CppAD::atomic_three< double >, public nimble_atomic_base {
 public:
  atomic_matinverse_class(const std::string& name);
 private:
  friend class atomic_cache_class<double>;
  friend class atomic_cache_class<CppAD::AD<double> >;
  atomic_cache_class<double> double_cache;
  atomic_cache_class< CppAD::AD<double> > CppADdouble_cache;

  typedef EigenTemplateTypes<double>::typeEigenConstMapStrd EigenConstMap;
  typedef EigenTemplateTypes<double>::typeEigenMapStrd EigenMap;
  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd metaEigenConstMap;
  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd metaEigenMap;
  typedef EigenTemplateTypes<CppAD::AD<double>>::typeMatrixXd metaEigenMatrixXd;
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
		       const CppAD::vector<CppAD::AD<double> >&  parameter_x  ,
		       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
		       size_t                              need_y       ,
		       size_t                              order_low    ,
		       size_t                              order_up     ,
		       const CppAD::vector<CppAD::AD<double> >&    taylor_x     ,
		       CppAD::vector<CppAD::AD<double> >&          taylor_y     );
  
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

atomic_matinverse_class* new_atomic_matinverse(void* tape_mgr, const std::string& name);
void delete_atomic_matinverse(void* tape_mgr, atomic_matinverse_class *atomic_matinverse);

#endif
