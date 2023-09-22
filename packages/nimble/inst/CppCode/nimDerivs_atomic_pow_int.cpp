/*********************************/
/*****      pow_int           ****/
/*********************************/
/*
 y = a^b
 b must be integer
 b cannot be a CppAD variable.  That is, there can be no derivatives wrt b.
 
 dy = b a^(b-1) dx
 d2y = b (b-1)
 
 a is x[0]
 b is x[1]
 y is y[0]

*/
#include <nimble/nimDerivs_atomic_pow_int.h>
#include <nimble/nimDerivs_atomic_zround.h>

atomic_pow_int_class::atomic_pow_int_class(const std::string& name) : 
  CppAD::atomic_three<double>(name)
{ }

bool atomic_pow_int_class::for_type(
	      const CppAD::vector<double>&               parameter_x ,
	      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
	      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = max(type_x[0], type_x[1]); // type of y is type of a. type of b is ignored.
  return true;
}
// Not sure this is ever needed.
bool atomic_pow_int_class::for_type(
	      const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
	      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
	      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = max(type_x[0], type_x[1]);
  return true;
}

bool atomic_pow_int_class::rev_depend(
				      const CppAD::vector<double>&          parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      CppAD::vector<bool>&                depend_x    ,
				      const CppAD::vector<bool>&          depend_y
				      ) {
  depend_x[0] = depend_y[0];
  depend_x[1] = depend_y[0]; // y is also calculated from b, even though b should never be CppAD::variable.
  return true;
}

// This makes a version of round for CppAD with derivatives set to 0.
// This is needed for meta-taping b * nimDerivs_pow_int(a, b-1).
// the "b * <stuff>" will create a derivative in b, but using
// CppADround(b) * <stuff> sets that to 0.
//  Elsewhere I encountered difficulties with this approach, so have replaced it with an atomic
// double CppADround(const double &x) {return round(x);}
// CPPAD_DISCRETE_FUNCTION(double, CppADround)
  
bool atomic_pow_int_class::forward(
     const CppAD::vector<double>&               parameter_x  ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
     size_t                              need_y       ,
     size_t                              order_low    ,
     size_t                              order_up     ,
     const CppAD::vector<double>&               taylor_x     ,
     CppAD::vector<double>&                     taylor_y     )
{
  int nrow(order_up + 1);
  double a(taylor_x[0]);
  int b(round(taylor_x[0 + 1*nrow])); // b must be integer.  We could error-trap if needed.
  bool a_zero(a==0);
  double log_fabs_a(log(fabs(a)));
  // sign of y=a^b is: 1 if a >= 0 (a == 0 will become 0 anyway)
  //                     even(b) if a < 0, where even(b) gives 1 for b even and -1 for b odd.
  double sign_a = a >= 0 ? 1 : -1;
  double sign_y = a >= 0 ? 1 : ( b % 2 == 0 ? 1 : -1); 
  if((order_low <= 0) && (order_up >= 0)) {
    // Cases with a == 0 take special care:
    // 0^(+ive) = 0 by special case.
    // 0^0 becomes exp(azmul(0, -Inf)) = exp(0) = 1, correct
    // 0^(-ive) becomes exp(azmul(-ive, -Inf)) = exp(Inf) = Inf, correct
    //  Cases with a != 0 are fine.
    taylor_y[0] = (a_zero && (b > 0)) ? 0 : sign_y * exp( CppAD::azmul(b, log_fabs_a) ); // a^b = exp(b * log(a))
  }
  if((order_low <= 1) && (order_up >= 1)) {
    // dy = b a^(b-1) da = b * exp( (b-1) * log(a) ) * da
    double pow_a_bm1 = (a_zero && (b > 1)) ? 0 : exp( CppAD::azmul(b-1, log_fabs_a) );
    taylor_y[1] = sign_a * sign_y * CppAD::azmul(b, pow_a_bm1) * taylor_x[1];
  }
  return true;
}

bool atomic_pow_int_class::forward(
				   const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                              need_y       ,
				   size_t                              order_low    ,
				   size_t                              order_up     ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
				   CppAD::vector< CppAD::AD<double> >&                     taylor_y     ) {
  int nrow(order_up + 1);
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = nimDerivs_pow_int(taylor_x[0], taylor_x[0 + 1*nrow]);
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = CppAD::azmul(nimDerivs_nimRound(taylor_x[0 + 1*nrow]),
			       nimDerivs_pow_int(taylor_x[0], taylor_x[0 + 1*nrow] - 1)) * taylor_x[1];
  }
  return true;
}

bool atomic_pow_int_class::reverse(
     const CppAD::vector<double>&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector<double>&               taylor_x    ,
     const CppAD::vector<double>&               taylor_y    ,
     CppAD::vector<double>&                     partial_x   ,
     const CppAD::vector<double>&               partial_y   )
{
  // <y_adj, dy> = <y_adj, b a^(b-1) da> = <a_adj, da>
  // giving a_adj = b a^(b-1) y_adj.
  // This uses that b_adj = 0 because pow_int is a function of a with b an integer.
  //
  // y_dot = b a^(b-1) a_dot
  // <y_dot_adj, dy_dot> = < y_dot_adj, b (b-1) a^(b-2) a_dot da + b a^(b-1) da_dot>
  //                     = < y_dot_adj b (b-1) a^(b-2) a_dot, da> + <y_dot_adj b a^(b-1), da_dot>
  // giving a_adj += y_dot_adj b (b-1) a^(b-2) a_dot
  // and a_dot_adj = y_dot_adj b a^(b-1)
  //
  // and all b adjoint terms should be 0
  
  int nrow(order_up + 1);
  // std::cout<<"reverse "<<nrow<<" "<<partial_x.size()<<" "<<partial_y.size()<<std::endl;
  double a(taylor_x[0]);
  int b(round(taylor_x[0 + 1*nrow])); // b must be integer.  We could error-trap if needed.
  bool a_zero(a==0);
  double log_fabs_a(log(fabs(a)));
  double sign_a = a >= 0 ? 1 : -1;
  double sign_y = a >= 0 ? 1 : ( b % 2 == 0 ? 1 : -1); 
  //  double dy_da(a_zero ? 0 : sign_a * sign_y * b * exp( (b-1) * log_fabs_a) );
  double pow_a_bm1( (a_zero && (b > 1)) ? 0 : exp( CppAD::azmul(b-1, log_fabs_a) ) );
  double dy_da(sign_a * sign_y * CppAD::azmul(b, pow_a_bm1) );
  partial_x[0] = 0;  // a_adjoint
  partial_x[0 + 1*nrow] = 0; // b_adjoint
  if(order_up >= 0) {
    partial_x[0] += dy_da * partial_y[0];
  }
  if(order_up >= 1) {
    partial_x[1] = dy_da * partial_y[1]; // a_dot_adjoint
    double pow_a_bm2( (a_zero && (b > 2)) ? 0 : exp( CppAD::azmul(b-2, log_fabs_a) ) );
    partial_x[0] += sign_y * CppAD::azmul(b * (b-1), pow_a_bm2) * taylor_x[1] * partial_y[1];
    partial_x[1 + 1*nrow] = 0; // b_dot_adjoint
  }
  return true;
}

bool atomic_pow_int_class::reverse(
				   const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   size_t                              order_up    ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
				   CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
				   const CppAD::vector< CppAD::AD<double> >&               partial_y   ) {
  int nrow(order_up + 1);
  // std::cout<<"meta-reverse "<<nrow<<" "<<partial_x.size()<<" "<<partial_y.size()<<std::endl;
  CppAD::AD<double> rb = nimDerivs_nimRound(taylor_x[0 + 1*nrow]);
  CppAD::AD<double> dy_da = CppAD::azmul(rb, nimDerivs_pow_int(taylor_x[0], taylor_x[0 + 1*nrow] - 1));
  partial_x[0] = 0;  // a_adjoint
  partial_x[0 + 1*nrow] = 0; // b_adjoint
  if(order_up >= 0) {
    partial_x[0] += dy_da * partial_y[0];
  }
  if(order_up >= 1) {
    partial_x[1] = dy_da * partial_y[1]; // a_dot_adjoint
    partial_x[0] += CppAD::azmul(rb * (rb-1),
				 nimDerivs_pow_int(taylor_x[0], taylor_x[0 + 1*nrow] - 2)) * taylor_x[1] * partial_y[1];
    partial_x[1 + 1*nrow] = 0; // b_dot_adjoint
  }
  return true;
};

CppAD::AD<double> nimDerivs_pow_int(const CppAD::AD<double> &a,
				    const CppAD::AD<double> &b) {
  // We can't use statics for reasons of multiple DLLs and bad interaction with CppAD statics.
  //  static atomic_pow_int_class atomic_pow_int("nimDerivs_pow_int");
  atomic_pow_int_class* atomic_pow_int;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_pow_int = new atomic_pow_int_class("atomic_pow_int");
  } else {
    atomic_pow_int = track_atomic_pow_int(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(2);
  in[0] = a;
  in[1] = b;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_pow_int)(in, out);
  if(!recording) {
    delete atomic_pow_int;
  }
  return out[0];
}
