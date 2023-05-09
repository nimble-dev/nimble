#include <nimble/nimDerivs_atomic_probit.h>

/**********************/
/* probit and iprobit */
/**********************/

atomic_iprobit_class::atomic_iprobit_class(const std::string& name) : 
  unary_atomic_class<double>(name)        
{ }

bool atomic_iprobit_class::forward(
				   const CppAD::vector<double>&               parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                              need_y       ,
				   size_t                              order_low    ,
				   size_t                              order_up     ,
				   const CppAD::vector<double>&               taylor_x     ,
				   CppAD::vector<double>&                     taylor_y     )
{
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = Rf_pnorm5(taylor_x[0], 0, 1, 1, 0); // (x, mu, sigma, lower_tail, log)
    // F(x)
  }
  double Fprime(0.);
  if((order_low <= 1) && (order_up >= 1)) {
    Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
    taylor_y[1] = Fprime * taylor_x[1];
    // F'(x) (x')
  }
  // Actually nimble will never call Forward above order 1, so following code is untested.
  if((order_low <= 2) && (order_up >= 2)) {
    if(Fprime == 0) {     // only calculate Fprime if it is really needed
      if(taylor_x[2] != 0)
	Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
    }
    taylor_y[2] = 0.5 * (- taylor_x[0]) * taylor_y[1] * taylor_x[1] +
      Fprime*taylor_x[2]; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
  }
  return true;
}

bool atomic_iprobit_class::forward(
				   const CppAD::vector< CppAD::AD<double> >&  parameter_x  ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				   size_t                                     need_y       ,
				   size_t                                     order_low    ,
				   size_t                                     order_up     ,
				   const CppAD::vector< CppAD::AD<double> >& taylor_x     ,
				   CppAD::vector< CppAD::AD<double> >&       taylor_y     ) {
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = nimDerivs_iprobit(taylor_x[0]);
    // F(x)
  }
  if(order_up == 0) return true;
  CppAD::AD<double> Fprime =
    exp(-(CppAD::AD<double>(0.5)) * (taylor_x[0] * taylor_x[0]) - CppAD::AD<double>(M_LN_SQRT_2PI));
  
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = Fprime * taylor_x[1];
    // F'(x) (x')
  }
  // Actually nimble will never call Forward above order 1, so following code is untested.
  if((order_low <= 2) && (order_up >= 2)) {
    taylor_y[2] = 0.5 * (- taylor_x[0]) * taylor_y[1] * taylor_x[1] +
      Fprime*taylor_x[2]; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
  }
  return true;
}

bool atomic_iprobit_class::reverse(
				   const CppAD::vector<double>&               parameter_x ,
				   const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				   size_t                              order_up    ,
				   const CppAD::vector<double>&               taylor_x    ,
				   const CppAD::vector<double>&               taylor_y    ,
				   CppAD::vector<double>&                     partial_x   ,
				   const CppAD::vector<double>&               partial_y   )
{
  double Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
  partial_x[0] = 0;
  if(order_up >= 1) partial_x[1] = 0;
  if(order_up >= 2) return false; // not implemented
  if(order_up >= 1) { // Needed for 2nd derivs
    partial_x[0] += partial_y[1] * (- taylor_x[0]) * Fprime * taylor_x[1];
    partial_x[1] += partial_y[1] * Fprime;
  } 
  // Needed for 1st derivs
  partial_x[0] += partial_y[0] * Fprime;
  return true;
}

bool atomic_iprobit_class::reverse(
      const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
      const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
      CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
      const CppAD::vector< CppAD::AD<double> >&               partial_y   ) {
  CppAD::AD<double> Fprime = exp(-(CppAD::AD<double>(0.5)) * (taylor_x[0] * taylor_x[0]) - CppAD::AD<double>(M_LN_SQRT_2PI));
  partial_x[0] = 0;
  if(order_up >= 1) partial_x[1] = 0;
  if(order_up >= 2) return false; // not implemented
  if(order_up >= 1) {             // Needed for 2nd derivs
    partial_x[0] += partial_y[1] * (- taylor_x[0]) * Fprime * taylor_x[1];
    partial_x[1] += partial_y[1] * Fprime;
  } 
  // Needed for 1st derivs
  partial_x[0] += partial_y[0] * Fprime;
  return true;
}

CppAD::AD<double> nimDerivs_iprobit(const CppAD::AD<double> x) {
  atomic_iprobit_class* atomic_iprobit;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_iprobit = new atomic_iprobit_class("atomic_iprobit");
  } else {
    atomic_iprobit = track_atomic_iprobit(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(1);
  in[0] = x;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_iprobit)(in, out);
  if(!recording) {
    delete atomic_iprobit;
  }
  return out[0];
}

/****************/

atomic_probit_class::atomic_probit_class(const std::string& name) : 
  unary_atomic_class<double>(name)        
{ }

bool atomic_probit_class::forward(
				  const CppAD::vector<double>&               parameter_x  ,
				  const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				  size_t                              need_y       ,
				  size_t                              order_low    ,
				  size_t                              order_up     ,
				  const CppAD::vector<double>&               taylor_x     ,
				  CppAD::vector<double>&                     taylor_y     )
{
  // This can do more with checking sizes, need_y, orders, and
  // what is a parameter.
  // need_y determines whether all elements need calculation,
  // or only those for parameters, or only those for variables.
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = Rf_qnorm5(taylor_x[0], 0, 1, 1, 0); // (p, mu, sigma, lower_tail, log)
    // F(x)
  }
  double z = taylor_y[0];
  double invFprime(0.);
  if(order_low <= 1 && order_up >= 1) {
    invFprime = Rf_dnorm4(z, 0, 1, 0);
    taylor_y[1] = taylor_x[1] / invFprime;
    // F'(x) (x')
  }
  // Following will not be called by nimble, so it is untested.
  if((order_low <= 2) && (order_up >= 2)) {
    taylor_y[2] = 0.5 * z * taylor_y[1] * taylor_y[1];
    if(taylor_x[2] != 0) {
      if(invFprime == 0)
	invFprime = Rf_dnorm4(z, 0, 1, 0);
      taylor_y[2] += taylor_x[2] / invFprime; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
    }
  }
  return true;
}

bool atomic_probit_class::forward(
				  const CppAD::vector< CppAD::AD<double> >&  parameter_x  ,
				  const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				  size_t                                     need_y       ,
				  size_t                                     order_low    ,
				  size_t                                     order_up     ,
				  const CppAD::vector< CppAD::AD<double> >& taylor_x     ,
				  CppAD::vector< CppAD::AD<double> >&       taylor_y     ) {
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = nimDerivs_probit(taylor_x[0]);
    // F(x)
  }
  if(order_up == 0) return true;
  CppAD::AD<double> z = taylor_y[0];
  CppAD::AD<double> invFprime = exp(-(CppAD::AD<double>(0.5)) * z * z - CppAD::AD<double>(M_LN_SQRT_2PI));
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = taylor_x[1] / invFprime;
    // F'(x) (x')
  }
  // Following will not be called by nimble, so it is untested.
  if((order_low <= 2) && (order_up >= 2)) {
    taylor_y[2] = 0.5 * z * taylor_y[1] * taylor_y[1];
    taylor_y[2] += taylor_x[2] / invFprime; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
  }
  return true;
}

bool atomic_probit_class::reverse(
				  const CppAD::vector<double>&               parameter_x ,
				  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				  size_t                              order_up    ,
				  const CppAD::vector<double>&               taylor_x    ,
				  const CppAD::vector<double>&               taylor_y    ,
				  CppAD::vector<double>&                     partial_x   ,
				  const CppAD::vector<double>&               partial_y   )
{
  double z = taylor_y[0];
  double invFprime = Rf_dnorm4(z, 0, 1, 0);
  partial_x[0] = 0;
  if(order_up >= 1) partial_x[1] = 0;
  if(order_up >= 2) return false; // not implemented
  if(order_up >= 1) { // Needed for 2nd derivs
    partial_x[0] += partial_y[1] * (z/(invFprime*invFprime)) * taylor_x[1];
    partial_x[1] += partial_y[1] / invFprime;
  } 
  // Needed for 1st derivs
  partial_x[0] += partial_y[0] / invFprime;
  return true;
};

bool atomic_probit_class::reverse(
				  const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
				  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				  size_t                              order_up    ,
				  const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
				  const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
				  CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
				  const CppAD::vector< CppAD::AD<double> >&               partial_y   ) {
  CppAD::AD<double> z = taylor_y[0];
  CppAD::AD<double> invFprime = exp(-(CppAD::AD<double>(0.5)) * z * z - CppAD::AD<double>(M_LN_SQRT_2PI));
  partial_x[0] = 0;
  if(order_up >= 1) partial_x[1] = 0;
  if(order_up >= 2) return false; // not implemented
  if(order_up >= 1) { // Needed for 2nd derivs
    partial_x[0] += partial_y[1] * (z/(invFprime*invFprime)) * taylor_x[1];
    partial_x[1] += partial_y[1] / invFprime;
  } 
  // Needed for 1st derivs
  partial_x[0] += partial_y[0] / invFprime;
  return true;
}

CppAD::AD<double> nimDerivs_probit(const CppAD::AD<double> x) {
  atomic_probit_class* atomic_probit;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_probit = new atomic_probit_class("atomic_probit");
  } else {
    atomic_probit = track_atomic_probit(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  CppAD::vector< CppAD::AD<double> > in(1);
  in[0] = x;
  CppAD::vector< CppAD::AD<double> > out(1);
  (*atomic_probit)(in, out);
  if(!recording) {
    delete atomic_probit;
  }
  return out[0];
}
