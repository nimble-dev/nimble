#ifndef _NIMDERIVS_ATOMIC_BACKSOLVE
#define _NIMDERIVS_ATOMIC_BACKSOLVE

#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimDerivs_vecmat_utils.h"
#include "nimDerivs_atomic_matmult.h"

void atomic_backsolve(const MatrixXd_CppAD &A,
		      const MatrixXd_CppAD &B,
		      MatrixXd_CppAD &y);

MatrixXd_CppAD nimDerivs_EIGEN_BS(const MatrixXd_CppAD &A,
				   const MatrixXd_CppAD &B);

class atomic_backsolve_class :  public CppAD::atomic_three< double > {
 public:
  atomic_backsolve_class(const std::string& name);
 private:
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


#endif

