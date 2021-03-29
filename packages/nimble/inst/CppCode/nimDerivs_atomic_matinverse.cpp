#include <nimble/nimDerivs_atomic_matinverse.h>
#include "nimDerivs_atomic_cache.cpp"
/*
Atomic class for matrix inverse.
This benefits from TMB's implementation but is updated to CppAD's atomic_three system.

We follow Giles' (2008; arXiv) useful results for first-order forward and reverse.
Matrix multiplication will be * instead of %*%
-----
Value:
Y = X^{-1} = solve(X).
X is n-x-n.
Y is n-x-n.

-----
Forward first order
YX = I
dY * X + Y * dX = 0
dY = -Y * dX * X^-1 = -Y * dX * Y

-----
Reverse first order
According to Giles(2008):
Xadjoint = -t(Y) %*% Yadjoint %*% t(Y);

This comes from (my notation)
<Yadjoint, dY> = <Yadjoint,  -Y * dX * Y>
               = <Yadjoint * Y^T, -Y * dX > // See rules in matmult
               = < -Y^T * Yadjoint * Y^T , dX>
               = < Xadjoint, dX >
-----
Reverse second order
dS = <Yadjoint, dY> + <Ydot_adjoint, dYdot>

Xadjoint =  -(Y^T * Yadjoint * Y^T) +  Y^T * Ydot_adjoint * (Y * Xdot * Y)^T + (Y * Xdot * Y)^T * Ydot_adjoint * Y^T
Xdot_adjoint = -Y^T * Ydot_adjoint * Y^T

Thus comes from:
(First term same as from reverse first order)
Second term:
Fdot = -Y * Xdot * Y = -X^-1 * Xdot * X^-1
dFdot = -(dY * Xdot * Y) -(Y * dXdot * Y) - (Y * Xdot * dY),
    where dY = d(X^-1) = -Y * dX * Y
dFdot = (Y * dX * Y * Xdot * Y) - (Y * dXdot * Y) + (Y * dXdot * Y * dX * Y)
<Ydot_adjoint, dYdot> = <Ydot_adjoint,  Y * dX * Y * Xdot * Y> + <Ydot_adjoint, -Y * dXdot * Y> + <Ydot_adjoint, Y * Xdot * Y * dX * Y>
                      = <Y^T * Ydot_adjoint * (Y * Xdot * Y)^T, dX> + <-Y^T * Ydot_adjoint * Y^T , dXdot > + <(Y * Xdot * Y)^T * Ydot_adjoint * Y^T  , dX>
                      = <Y^T * Ydot_adjoint * (Y * Xdot * Y)^T + (Y * Xdot * Y)^T * Ydot_adjoint * Y^T, dX> +  <-Y^T * Ydot_adjoint * Y^T , dXdot >
                      = <Xadjoint_term, dX> + <Xdot_adjoint_term, dXdot>
*/

atomic_matinverse_class::atomic_matinverse_class(const std::string& name) :
  CppAD::atomic_three<double>(name) {};

bool atomic_matinverse_class::for_type(
				       const CppAD::vector<double>&               parameter_x ,
				       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				       CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  // printf("In matinverse for_type\n");
  size_t n = type_y.size();
  // All types must be the same because matrix inversion is an all-to-all mapping 
  CppAD::ad_type_enum final_type(CppAD::constant_enum);
  CppAD::ad_type_enum this_type;
  for(size_t i = 0; i < n; ++i) {
    this_type = type_x[i];
    if(this_type == CppAD::dynamic_enum)
      final_type = CppAD::dynamic_enum;
    if(this_type == CppAD::variable_enum) {
      final_type = CppAD::variable_enum;
      break;
    }
  }
  for(size_t i = 0; i < n; ++i) type_y[i] = final_type;
  return true;
}

bool atomic_matinverse_class::rev_depend(
					 const CppAD::vector<double>&          parameter_x ,
					 const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
					 CppAD::vector<bool>&                depend_x    ,
					 const CppAD::vector<bool>&          depend_y
					 ) {
  // printf("In matinverse reverse_depend\n");
  bool any_depend_y(false);
  size_t n = depend_y.size();
  for(size_t i = 0; i < n; ++i) {
    if(depend_y[i]) {
      any_depend_y = true;
      break;
    }
  }
  if(depend_x.size() != n) {
    printf("In matinverse rev_depend, somehow size of depend_x does not match size of depend_y.  That should never happen.\n");
  }
  n = depend_x.size();
  for(size_t i = 0; i < n; ++i) {
    depend_x[i] = any_depend_y;
  }
  return false;
}

bool atomic_matinverse_class::forward(
				      const CppAD::vector<double>&               parameter_x  ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				      size_t                              need_y       ,
				      size_t                              order_low    ,
				      size_t                              order_up     ,
				      const CppAD::vector<double>&               taylor_x     ,
				      CppAD::vector<double>&                     taylor_y     ) {
  //forward mode
  // printf("In matinverse forward\n");
  int nrow = order_up + 1;
  //  std::cout<<"in matinverse forward with order_low = "<<order_low<<" and order_up = "<<order_up<<std::endl;
    
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
  EigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  //  std::cout<<"Xmap inputted for forward 0\n"<<Xmap<<std::endl;
  
  // std::cout<<"Xmap\n"<<Xmap<<std::endl;
  if(order_low <= 0 & order_up >= 0) { // value
    // We could compile different cases depending on need for strides or not.  
    EigenMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    Ymap = Xmap.inverse().eval(); // This eval is necessary if nrow > 1
    // std::cout<<"Ymap calculated for forward 0\n"<<Ymap<<std::endl;
    double_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
    
  }
  if(order_low <= 1 & order_up >= 1) {
    // printf("In forward >1\n");
    double_cache.check_and_set_cache(this,
				     parameter_x,
				     type_x,
				     0,
				     order_up,
				     taylor_x,
				     taylor_y.size());
    int cache_nrow = double_cache.nrow();
    EigenMap Ymap(double_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );

    //    std::cout<<"cached Ymap "<<Ymap<<std::endl;
    //    std::cout<<"cache_nrow = "<<cache_nrow<<" n = "<<n<<std::endl;
    //    double_cache.show_taylor_y();

    //    std::cout<<"taylor Y in cache "<<std::endl;
    //double_cache.show_taylor_y();
    
    // EigenMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    // std::cout<<"inputted Ymap "<<Ymap<<std::endl;
    EigenMap dYmap(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    EigenConstMap dXmap(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow));
    dYmap = -Ymap * dXmap * Ymap;
    // double_cache.set_cache( 1, 1, order_up, taylor_x, taylor_y ); // This would be the right way to store 1st order, but 1st order is not used by any reverse mode orders, so it is not needed.
  }
  return true;
}

bool atomic_matinverse_class::forward(
				      const CppAD::vector<CppAD::AD<double> >&  parameter_x  ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				      size_t                              need_y       ,
				      size_t                              order_low    ,
				      size_t                              order_up     ,
				      const CppAD::vector<CppAD::AD<double> >&    taylor_x     ,
				      CppAD::vector<CppAD::AD<double> >&          taylor_y     ) {
  //forward mode
  // printf("In matinverse meta-forward\n");
  int nrow = order_up + 1;
  // std::cout<<"nrow = "<<nrow<<std::endl;
    
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
  metaEigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  if(order_low <= 0 & order_up >= 0) { // value
    // We could compile different cases depending on need for strides or not.  
    metaEigenMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    Ymap = nimDerivs_matinverse(Xmap);
    CppADdouble_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
  }
  if(order_low <= 1 & order_up >= 1) {
    // printf("In forward >1\n");
    CppADdouble_cache.check_and_set_cache(this,
					  parameter_x,
					  type_x,
					  0,
					  order_up,
					  taylor_x,
					  taylor_y.size());
    int cache_nrow = CppADdouble_cache.nrow();
    metaEigenMap Ymap(CppADdouble_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
    metaEigenMap dYmap(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    metaEigenConstMap dXmap(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow));
    dYmap = nimDerivs_matmult(-Ymap, nimDerivs_matmult( dXmap,  Ymap ) );
    // CppADdouble_cache.set_cache( 1, 1, order_up, taylor_x, taylor_y ); // This would be the right way to store 1st order, but 1st order is not used by any reverse mode orders, so it is not needed.
  }
  return true;
}
  
bool atomic_matinverse_class::reverse(
				      const CppAD::vector<double>&               parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      size_t                              order_up    ,
				      const CppAD::vector<double>&               taylor_x    ,
				      const CppAD::vector<double>&               taylor_y    ,
				      CppAD::vector<double>&                     partial_x   ,
				      const CppAD::vector<double>&               partial_y   )
{
  //reverse mode
  //  printf("In matinverse reverse\n");
  //  std::cout<<"In matinverse reverse"<<std::endl;
  int nrow = order_up + 1;
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));

  double_cache.check_and_set_cache(this,
				   parameter_x,
				   type_x,
				   0, // only use cached values up to order 0
				   order_up,
				   taylor_x,
				   taylor_y.size());
  int cache_nrow = double_cache.nrow();  
  //  std::cout<<"cache_nrow = "<<cache_nrow<<std::endl;
  EigenConstMap Ymap(double_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
  //  std::cout<<"Ymap\n"<<Ymap<<std::endl;

  EigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  // std::cout<<"Xmap\n"<<Xmap<<std::endl;
  // Eigen::MatrixXd Y = Xmap.inverse();
  // std::cout<<"Y (recomputed) \n"<<Y<<std::endl;
  
  EigenMap Xadjoint_map(&partial_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  if(order_up >= 0) {
    EigenConstMap Yadjoint_map(&partial_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    // std::cout<<"Yadjoint_map\n"<<Yadjoint_map<<std::endl;
    // Xadjoint_map = -Y.transpose() * Yadjoint_map *  Y.transpose();
    // std::cout<<"Xadjoint_map (from the recomputed Y)\n"<<Xadjoint_map<<std::endl;
    Xadjoint_map = -Ymap.transpose() * Yadjoint_map *  Ymap.transpose();
    // std::cout<<"Xadjoint_map (as is)\n"<<Xadjoint_map<<std::endl;
  }

  if(order_up >= 1) {
    EigenConstMap Xdot_map(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow ) );     
    EigenConstMap Ydot_adjoint_map(&partial_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    EigenMap Xdot_adjoint_map(&partial_x[1], n, n, EigStrDyn(nrow*n, nrow) );
    Eigen::MatrixXd Y_Xdot_Y_transpose = (Ymap * Xdot_map * Ymap).transpose();
    Eigen::MatrixXd Ydot_adjoint_Ytranspose = Ydot_adjoint_map * Ymap.transpose();
    Xadjoint_map += Ymap.transpose() * Ydot_adjoint_map *  Y_Xdot_Y_transpose +
      Y_Xdot_Y_transpose * Ydot_adjoint_Ytranspose;
    Xdot_adjoint_map = -Ymap.transpose() *  Ydot_adjoint_Ytranspose;
  }
  return true;
}

bool atomic_matinverse_class::reverse(
				      const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      size_t                              order_up    ,
				      const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
				      const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
				      CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
				      const CppAD::vector<CppAD::AD<double> >&               partial_y   )
{
  //reverse mode
  // printf("In matinverse reverse\n");
  int nrow = order_up + 1;
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));

  CppADdouble_cache.check_and_set_cache(this,
					parameter_x,
					type_x,
					0, // only use cached values up to order 0
					order_up,
					taylor_x,
					taylor_y.size());
  int cache_nrow = CppADdouble_cache.nrow();
  metaEigenConstMap Ymap(CppADdouble_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
  metaEigenMap Xadjoint_map(&partial_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  if(order_up >= 0) {
    metaEigenConstMap Yadjoint_map(&partial_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    Xadjoint_map = nimDerivs_matmult(-Ymap.transpose(), nimDerivs_matmult( Yadjoint_map,  Ymap.transpose() ) );
  }

  if(order_up >= 1) {
    metaEigenConstMap Xdot_map(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow ) );     
    metaEigenConstMap Ydot_adjoint_map(&partial_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    metaEigenMap Xdot_adjoint_map(&partial_x[1], n, n, EigStrDyn(nrow*n, nrow) );
    metaEigenMatrixXd Y_Xdot_Y_transpose = nimDerivs_matmult(Ymap, nimDerivs_matmult( Xdot_map, Ymap)).transpose();
    metaEigenMatrixXd Ydot_adjoint_Ytranspose = nimDerivs_matmult(Ydot_adjoint_map, Ymap.transpose());
    Xadjoint_map += nimDerivs_matmult(Ymap.transpose(), nimDerivs_matmult( Ydot_adjoint_map,  Y_Xdot_Y_transpose)) +
      nimDerivs_matmult(Y_Xdot_Y_transpose, Ydot_adjoint_Ytranspose);
    Xdot_adjoint_map = nimDerivs_matmult(-Ymap.transpose(), Ydot_adjoint_Ytranspose);
  }
  return true;
}

void atomic_matinverse(const MatrixXd_CppAD &x, // This (non-template) type forces any incoming expression to be evaluated
			 MatrixXd_CppAD &y) {
  // static atomic_matinverse_class atomic_matinverse("atomic_matinverse"); // this has no state information so the same object can be used for all cases
  atomic_matinverse_class *atomic_matinverse;
  int n = x.rows();
  std::vector<CppAD::AD<double> > xVec(n*n);
  mat2vec(x, xVec);
  std::vector<CppAD::AD<double> > yVec(n*n);
  atomic_matinverse = new atomic_matinverse_class("atomic_matinverse");
  (*atomic_matinverse)(xVec, yVec);
  y.resize(n, n);
  vec2mat(yVec, y);  
}

MatrixXd_CppAD nimDerivs_matinverse(const MatrixXd_CppAD &x) {
  MatrixXd_CppAD ans;
  atomic_matinverse(x, ans);
  return ans;

}
