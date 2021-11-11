#ifndef _NIMDERIVS_ATOMIC_SOLVE_BASE
#define _NIMDERIVS_ATOMIC_SOLVE_BASE

#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "nimDerivs_vecmat_utils.h"
#include "nimDerivs_atomic_matmult.h"
#include "nimDerivs_atomic_cache.h"

#define USE_NEW_DYNAMIC_SOLVE_BASE

class atomic_solve_base_class {    
 public:
  friend class atomic_cache_class<double>;
  friend class atomic_cache_class<CppAD::AD<double> >;
  atomic_cache_class<double> double_cache;
  atomic_cache_class< CppAD::AD<double> > CppADdouble_cache;

  std::vector<double> X_stored;
  std::vector<CppAD::AD<double> > X_AD_stored;
  bool A_is_constant_, B_is_constant_;
  bool A_is_variable_, B_is_variable_;

  double * get_X_stored_ptr() {return &X_stored[0];}
  CppAD::AD<double> * get_X_AD_stored_ptr() {return &X_AD_stored[0];}
  void fill_X_AD_stored() {
    X_AD_stored.resize( X_stored.size() );
    for(int i = 0; i < X_stored.size(); ++i) X_AD_stored[i] = X_stored[i];
  };
  void clear_X_AD_stored() {X_AD_stored.clear();}
  std::vector<double> &get_X_stored() {return X_stored;}
 public:
  void set_X_stored(const MatrixXd_CppAD &X) {
    int n1 = X.rows();
    int n2 = X.cols();
    X_stored.resize(n1 * n2);
    mat2vec_v(X, X_stored, 0);
  }

  bool &Aconstant() {return A_is_constant_;}
  bool &Bconstant() {return B_is_constant_;}
  bool const &Aconstant() const {return A_is_constant_;}
  bool const &Bconstant() const {return B_is_constant_;}

  bool &Avariable() {return A_is_variable_;}
  bool &Bvariable() {return B_is_variable_;}
  bool const &Avariable() const {return A_is_variable_;}
  bool const &Bvariable() const {return B_is_variable_;}

  typedef EigenTemplateTypes<double>::typeEigenConstMapStrd EigenConstMap;
  typedef EigenTemplateTypes<double>::typeEigenMapStrd EigenMap;
  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenConstMapStrd metaEigenConstMap;
  typedef EigenTemplateTypes<CppAD::AD<double> >::typeEigenMapStrd metaEigenMap;
  typedef EigenTemplateTypes<CppAD::AD<double>>::typeMatrixXd metaEigenMatrixXd;

#ifdef USE_NEW_DYNAMIC_SOLVE_BASE
 public:
  EigenConstMap Amap, Bmap, dA_map, dB_map;
  EigenMap      Ymap, dY_map;
  EigenConstMap YmapC, Yadjoint_map, Ydot_adjoint_map, Ydot_map, Adot_map;
  EigenMap      Aadjoint_map, Badjoint_map, Adot_adjoint_map, Bdot_adjoint_map;

  metaEigenConstMap mAmap, mBmap, mdA_map, mdB_map;
  metaEigenMap      mYmap, mdY_map;
  metaEigenConstMap mYmapC, mYadjoint_map, mYdot_adjoint_map, mYdot_map, mAdot_map;
  metaEigenMap      mAadjoint_map, mBadjoint_map, mAdot_adjoint_map, mBdot_adjoint_map;
  // Badjoint calcs are used for Aadjoint even if B is constant.
  // Here is a little scheme to provide local memory instead of using CppAD-allocated
  // memory for Badjoint in the case that B is constant.
  std::vector<double> Badjoint_memory_if_needed;
  double *get_Badjoint_memory(int s) {Badjoint_memory_if_needed.resize(s);
    return &Badjoint_memory_if_needed[0]; }
  std::vector<double> Bdot_adjoint_memory_if_needed;
  double *get_Bdot_adjoint_memory(int s) {Bdot_adjoint_memory_if_needed.resize(s);
    return &Bdot_adjoint_memory_if_needed[0]; }
  std::vector<CppAD::AD<double>> Badjoint_AD_memory_if_needed;
  CppAD::AD<double> *get_Badjoint_AD_memory(int s) {Badjoint_AD_memory_if_needed.resize(s);
    return &Badjoint_AD_memory_if_needed[0]; }
  void clear_Badjoint_AD_memory() {Badjoint_AD_memory_if_needed.resize(0);}
  std::vector<CppAD::AD<double>> Bdot_adjoint_AD_memory_if_needed;
  CppAD::AD<double> *get_Bdot_adjoint_AD_memory(int s) {Bdot_adjoint_AD_memory_if_needed.resize(s);
    return &Bdot_adjoint_AD_memory_if_needed[0]; }
  void clear_Bdot_adjoint_AD_memory() {Bdot_adjoint_AD_memory_if_needed.resize(0);}
 public:
#endif
  
  /* virtual bool for_type( */
  /* 			const CppAD::vector<double>&               parameter_x , */
  /* 			const CppAD::vector<CppAD::ad_type_enum>&  type_x      , */
  /* 			CppAD::vector<CppAD::ad_type_enum>&        type_y      ); */

  /* virtual bool rev_depend( */
  /* 			  const CppAD::vector<double>&          parameter_x , */
  /* 			  const CppAD::vector<CppAD::ad_type_enum>&  type_x      , */
  /* 			  CppAD::vector<bool>&                depend_x    , */
  /* 			  const CppAD::vector<bool>&          depend_y */
  /* 			  ); */

  /* virtual bool forward( */
  /* 		       const CppAD::vector<double>&               parameter_x  , */
  /* 		       const CppAD::vector<CppAD::ad_type_enum>&  type_x       , */
  /* 		       size_t                              need_y       , */
  /* 		       size_t                              order_low    , */
  /* 		       size_t                              order_up     , */
  /* 		       const CppAD::vector<double>&               taylor_x     , */
  /* 		       CppAD::vector<double>&                     taylor_y     ); */

  /* virtual bool forward( */
  /* 		       const CppAD::vector<CppAD::AD<double> >&  parameter_x  , */
  /* 		       const CppAD::vector<CppAD::ad_type_enum>&  type_x       , */
  /* 		       size_t                              need_y       , */
  /* 		       size_t                              order_low    , */
  /* 		       size_t                              order_up     , */
  /* 		       const CppAD::vector<CppAD::AD<double> >&    taylor_x     , */
  /* 		       CppAD::vector<CppAD::AD<double> >&          taylor_y     ); */
  
  /* virtual bool reverse( */
  /* 		       const CppAD::vector<double>&               parameter_x , */
  /* 		       const CppAD::vector<CppAD::ad_type_enum>&  type_x      , */
  /* 		       size_t                              order_up    , */
  /* 		       const CppAD::vector<double>&               taylor_x    , */
  /* 		       const CppAD::vector<double>&               taylor_y    , */
  /* 		       CppAD::vector<double>&                     partial_x   , */
  /* 		       const CppAD::vector<double>&               partial_y   ); */
  /* virtual bool reverse( */
  /* 		       const CppAD::vector<CppAD::AD<double> >&               parameter_x , */
  /* 		       const CppAD::vector<CppAD::ad_type_enum>&  type_x      , */
  /* 		       size_t                              order_up    , */
  /* 		       const CppAD::vector<CppAD::AD<double> >&               taylor_x    , */
  /* 		       const CppAD::vector<CppAD::AD<double> >&               taylor_y    , */
  /* 		       CppAD::vector<CppAD::AD<double> >&                     partial_x   , */
  /* 		       const CppAD::vector<CppAD::AD<double> >&               partial_y   ); */

  atomic_solve_base_class();
};


#endif

