/*********************************/
/*****     log_pow_int        ****/
/*********************************/
/*
 y = f(a, b) = log(pow_int(a, b)) = round(b) * log(a), only for a >= 0
 The key is to handle log(0^0) correctly, with derivatives in b always 0.

 We thought it should work to use azmul(CPPAD_DISCRETE_FUNCTION(double, nimround)(b), log(a)), 
 where nimround is a function that returns round(x).

 However, a CPPAD_DISCRETE_FUNCTION fails through double-taping.  Furthermore,
 when we wrote our own atomic for zround (round(x) with all derivs set to 0),
 we still had the problem that double-taped reverse order failed because we would have
 adjoint of log(a) = 0 but then adjoint of a becomes Inf or NaN.

 So we are writing this atomic to handle the entire log(pow_int(a, b)) case.
 This also uses the atomic log_a_over_b, which enforces a/b = 0/0 = 0
 and derivatives wrt a = 0.

 df/da = round(b) / a, or 0 if round(b)/a = 0/0, i.e. if round(b) = 0
 df/db = 0 by definition
 d2f/da2 = -round(b) / a^2, or 0 if round(b) = 0
 d2f/dadb = 0
 d2f/db2 = 0

 dy = (df/da) da + (df/db) db

 <yadj, dy> = <aadj, da> + <badj, db>
 yadj (df/da) da + yadj (df/db) db = aadj da + badj db
 aadj = yadj (df/da) = yadj (b/a), or 0 if b/a = 0/0
 badj = (df/db) db = 0 by definition

dS = <ydotadj, dydot> + <yadj, dj>
ydot = (df/da) adot
dydot = d(df/da) adot + (df/da) dadot
dydot = (d2f/da2) da adot + (d2f/dadb) db adot + (df/da) dadot
<ydotadj, dydot> = ydotadj (d2f/da2) da adot + ydotadj (d2f/dadb) db adot + ydotadj (df/da) dadot
This identifies terms
aadj += ydotadj (d2f/da2) adot
badj += 0 b/c d2f/dadb = 0
adotadj = ydotadj (df/da)
bdotadj = 0
*/
#include <nimble/nimDerivs_atomic_log_pow_int.h>


atomic_log_pow_int_class::atomic_log_pow_int_class(const std::string& name) : 
  CppAD::atomic_three<double>(name)
{ }

bool atomic_log_pow_int_class::for_type(
				   const CppAD::vector<double>&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = max(type_x[0], type_x[1]);
  return true;
}
// Not sure this is ever needed.
bool atomic_log_pow_int_class::for_type(
				   const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = max(type_x[0], type_x[1]);
  return true;
}

bool atomic_log_pow_int_class::rev_depend(
				     const CppAD::vector<double>&          parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     CppAD::vector<bool>&                depend_x    ,
				     const CppAD::vector<bool>&          depend_y
				     ) {
  depend_x[0] = depend_x[1] = depend_y[0];
  return true;
}

bool atomic_log_pow_int_class::forward(
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
  
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = b == 0 ? 0 : b * log(a);
    
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = b == 0 ? 0 : (b/a) * taylor_x[1];
  }
  return true;
}

bool atomic_log_pow_int_class::forward(
				   const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                              need_y       ,
				   size_t                              order_low    ,
				   size_t                              order_up     ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
				   CppAD::vector< CppAD::AD<double> >&                     taylor_y     ) {
 int nrow(order_up + 1);
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = nimDerivs_log_pow_int(taylor_x[0], taylor_x[0 + 1*nrow]);
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = nimDerivs_zb_over_a(taylor_x[0], taylor_x[0 + 1*nrow]) * taylor_x[1];
  }
  return true;
}

bool atomic_log_pow_int_class::reverse(
     const CppAD::vector<double>&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector<double>&               taylor_x    ,
     const CppAD::vector<double>&               taylor_y    ,
     CppAD::vector<double>&                     partial_x   ,
     const CppAD::vector<double>&               partial_y   )
{
  int nrow(order_up + 1);
  double a(taylor_x[0]);
  int b(round(taylor_x[0 + 1*nrow]));

  if(order_up >= 0) {
    partial_x[0] = b == 0 ? 0 : (b/a) * partial_y[0]; //a_adjoint
    partial_x[0 + 1*nrow] = 0.; // b_adjoint
  }
  if(order_up >= 1) {
    partial_x[0] += (b == 0 ? 0 : (-b/(a*a)) * taylor_x[1] * partial_y[1]); // a_adj
    partial_x[1]  = (b == 0 ? 0 : (b/a) * partial_y[1]); // adot_adjoint
    partial_x[1 + 1*nrow] = 0; // bdot_adjoint
  }
  return true;
}

bool atomic_log_pow_int_class::reverse(
				       const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
				       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				       size_t                              order_up    ,
				       const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
				       const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
				       CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
				       const CppAD::vector< CppAD::AD<double> >&               partial_y   ) {
  int nrow(order_up + 1);
  CppAD::AD<double> a(taylor_x[0]);
  CppAD::AD<double> b(taylor_x[0 + 1*nrow]);
  if(order_up >= 0) {
    partial_x[0] = nimDerivs_zb_over_a(a, b) * partial_y[0]; //a_adjoint
    partial_x[0 + 1*nrow] = CppAD::AD<double>(0.); // b_adjoint
  }
  if(order_up >= 1) {
    partial_x[0] += -nimDerivs_zb_over_a(a*a, b) * taylor_x[1] * partial_y[1];
    partial_x[1]  = nimDerivs_zb_over_a(a, b) * partial_y[1]; // adot_adjoint
    partial_x[1 + 1*nrow] = CppAD::AD<double>(0.); // bdot_adjoint
  }
  return true;
};

CppAD::AD<double> nimDerivs_log_pow_int(const CppAD::AD<double> &a,
					const CppAD::AD<double> &b) {
  atomic_log_pow_int_class* atomic_log_pow_int;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_log_pow_int = new atomic_log_pow_int_class("atomic_log_pow_int");
  } else {
    atomic_log_pow_int = track_atomic_log_pow_int(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
						  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(2);
  in[0] = a;
  in[1] = b;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_log_pow_int)(in, out);
  if(!recording) {
    delete atomic_log_pow_int;
  }
  return out[0];
}

/*******************************************/
/*             zb_over_a                   */
/*******************************************/

/*
  f(a, b) = round(b) / a, with 0/0 = 0
 
  df/da = -round(b) / (a*a), or 0 if round(b) = 0
  df/db = 0
  d2f/da2 = 2 round(b) / (a*a*a), or 0 if round(b) = 0
  d2f/dadb = 0
  d2f/db2 = 0

  dy = (df/da) da
  ydot = (df/da) adot
  dydot = (d2f/da2) da adot + 0 + 0 + (df/da) dadot
  
  dS = <yadj, dy> = yadj (df/da) da
  giving
  aadj = yadj (df/da)
  

  dS += <ydotadj, dydot> 
  giving
  aadj += ydotadj (d2f/da2) adot
  adotadj = ydotadj (df/da)
*/

atomic_zb_over_a_class::atomic_zb_over_a_class(const std::string& name) : 
  CppAD::atomic_three<double>(name)
{ }

bool atomic_zb_over_a_class::for_type(
				   const CppAD::vector<double>&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = max(type_x[0], type_x[1]);
  return true;
}
// Not sure this is ever needed.
bool atomic_zb_over_a_class::for_type(
				   const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  type_y[0] = max(type_x[0], type_x[1]);
  return true;
}

bool atomic_zb_over_a_class::rev_depend(
				     const CppAD::vector<double>&          parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     CppAD::vector<bool>&                depend_x    ,
				     const CppAD::vector<bool>&          depend_y
				     ) {
  depend_x[0] = depend_x[1] = depend_y[0];
  return true;
}

bool atomic_zb_over_a_class::forward(
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
  
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = b == 0 ? 0 : b / a;
    
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = b == 0 ? 0 : (-b/(a*a)) * taylor_x[1];
  }
  return true;
}

bool atomic_zb_over_a_class::forward(
				   const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                              need_y       ,
				   size_t                              order_low    ,
				   size_t                              order_up     ,
				   const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
				   CppAD::vector< CppAD::AD<double> >&                     taylor_y     ) {
 int nrow(order_up + 1);
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = nimDerivs_zb_over_a(taylor_x[0], taylor_x[0 + 1*nrow]);
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = -nimDerivs_zb_over_a(taylor_x[0]*taylor_x[0], taylor_x[0 + 1*nrow]) * taylor_x[1];
  }
  return true;
}

bool atomic_zb_over_a_class::reverse(
     const CppAD::vector<double>&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector<double>&               taylor_x    ,
     const CppAD::vector<double>&               taylor_y    ,
     CppAD::vector<double>&                     partial_x   ,
     const CppAD::vector<double>&               partial_y   )
{
  int nrow(order_up + 1);
  double a(taylor_x[0]);
  int b(round(taylor_x[0 + 1*nrow]));

  if(order_up >= 0) {
    partial_x[0] = b == 0 ? 0 : (-b/(a*a)) * partial_y[0]; //a_adjoint
    partial_x[0 + 1*nrow] = 0.; // b_adjoint
  }
  if(order_up >= 1) {
    partial_x[0] += (b == 0 ? 0 : (2*b/(a*a*a)) * taylor_x[1] * partial_y[1]); // a_adj
    partial_x[1]  = (b == 0 ? 0 : (-b/(a*a)) * partial_y[1]); // adot_adjoint
    partial_x[1 + 1*nrow] = 0; // bdot_adjoint
  }
  return true;
}

bool atomic_zb_over_a_class::reverse(
				       const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
				       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				       size_t                              order_up    ,
				       const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
				       const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
				       CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
				       const CppAD::vector< CppAD::AD<double> >&               partial_y   ) {
  int nrow(order_up + 1);
  CppAD::AD<double> a(taylor_x[0]);
  CppAD::AD<double> b(taylor_x[0 + 1*nrow]);
  if(order_up >= 0) {
    partial_x[0] = -nimDerivs_zb_over_a(a*a, b) * partial_y[0]; //a_adjoint
    partial_x[0 + 1*nrow] = 0.; // b_adjoint
  }
  if(order_up >= 1) {
    partial_x[0] += 2 * nimDerivs_zb_over_a(a*a*a, b) * taylor_x[1] * partial_y[1];
    partial_x[1]  = -nimDerivs_zb_over_a(a*a, b) * partial_y[1]; // adot_adjoint
    partial_x[1 + 1*nrow] = 0; // bdot_adjoint
  }
  return true;
};

CppAD::AD<double> nimDerivs_zb_over_a(const CppAD::AD<double> &a,
					const CppAD::AD<double> &b) {
  atomic_zb_over_a_class* atomic_zb_over_a;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_zb_over_a = new atomic_zb_over_a_class("atomic_zb_over_a");
  } else {
    atomic_zb_over_a = track_atomic_zb_over_a(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
						  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(2);
  in[0] = a;
  in[1] = b;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_zb_over_a)(in, out);
  if(!recording) {
    delete atomic_zb_over_a;
  }
  return out[0];
}
