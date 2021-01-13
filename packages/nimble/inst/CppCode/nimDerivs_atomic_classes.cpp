#include <nimble/nimDerivs_atomic_classes.h>

atomic_lgamma_class::atomic_lgamma_class(const std::string& name, int baseOrder_) : 
  unary_atomic_class<double>(name),
  baseOrder(baseOrder_)
{ }

bool atomic_lgamma_class::forward(
     const CppAD::vector<double>&               parameter_x  ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
     size_t                              need_y       ,
     size_t                              order_low    ,
     size_t                              order_up     ,
     const CppAD::vector<double>&               taylor_x     ,
     CppAD::vector<double>&                     taylor_y     )
{
  // std::cout<<"lgamma forward baseOrder = "<<baseOrder<<std::endl;
  // std::cout<<"forward "<<order_low<<" "<<order_up<<" "<<taylor_x[0]<<std::endl;
     if(order_low <= 0 & order_up >= 0) {
	  if(baseOrder == 0)
	       taylor_y[0] = Rf_lgammafn(taylor_x[0]);
	  else
	       taylor_y[0] = Rf_psigamma(taylor_x[0], baseOrder-1);
	  //  std::cout<<"taylor_y[0] "<<taylor_y[0]<<std::endl;
     }
     double fprime;
     if(order_low <= 2 & order_up >= 1)
	  fprime = Rf_psigamma(taylor_x[0], baseOrder);
     if(order_low <= 1 & order_up >= 1) {
	  taylor_y[1] = fprime * taylor_x[1];
	  // std::cout<<"taylor_y[1] "<<taylor_y[1]<<std::endl;
	  // f'(x) (x')
     }
     if(order_low <= 2 & order_up >= 2) {
	  taylor_y[2] = 0.5 * (Rf_psigamma(taylor_x[0], baseOrder+1) * taylor_x[1] * taylor_x[1] + 
			       fprime * 2 * taylor_x[2]);
	  // 0.5 * ((f''(x)) (x')^2 + 2*f'(x) (x''))
	  // Note x'' is a taylor coeff so it is 0.5*2nd deriv
     }
     return true;
}
bool atomic_lgamma_class::reverse(
     const CppAD::vector<double>&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector<double>&               taylor_x    ,
     const CppAD::vector<double>&               taylor_y    ,
     CppAD::vector<double>&                     partial_x   ,
     const CppAD::vector<double>&               partial_y   )
{
  // std::cout<<"lgamma reverse baseOrder = "<<baseOrder<<std::endl;
  // std::cout<<"reverse "<<order_up<<" "<<taylor_x[0]<<" "<<partial_y[0]<<std::endl;
     partial_x[0] = 0;
     if(order_up >= 1) partial_x[1] = 0;
     double fprime = Rf_psigamma(taylor_x[0], baseOrder);
     if(order_up >= 1) {
	  partial_x[1] += partial_y[1] * fprime;
	  partial_x[0] += partial_y[1] * Rf_psigamma(taylor_x[0], baseOrder+1) * taylor_x[1];
	  // std::cout<<partial_x[1]<<" ";
     }
     partial_x[0] += partial_y[0] * fprime;
     // std::cout<<"fprime = "<<fprime<<" partial_y[0] = "<<partial_y[0]<<" partial_x[0] = "<<partial_x[0]<<" partial_y[1] = "<<partial_y[1]<<" partial_x[1] = "<<partial_x[1]<<"\n";
     // std::cout<<partial_x[0]<<std::endl;
     // dG/dx = dG/dy  dy/dx 
     return true;
}

bool atomic_lgamma_class::forward(
     const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
     size_t                              need_y       ,
     size_t                              order_low    ,
     size_t                              order_up     ,
     const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
     CppAD::vector< CppAD::AD<double> >&                     taylor_y     )
{
  // printf("In lgamma meta-forward for orders %lu %lu\n", order_low, order_up);
  //std::cout<<"lgamma meta-forward baseOrder = "<<baseOrder<<std::endl;
  // std::cout<<"tape_id and handle:"<< CppAD::AD<double>::get_tape_id_nimble()<<" "<< CppAD::AD<double>::get_tape_handle_nimble()<<"\n";
  // std::cout<<"atomic info:"<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<"\n";
     if(order_low <= 0 & order_up >= 0) {
	  taylor_y[0] = nimDerivs_lgammafn(taylor_x[0], baseOrder); // This puts it in the new tape being recorded
     }
     // std::cout<<"Constant(taylor_y[0]) = "<<CppAD::Constant(taylor_y[0])<<"\n";
     // std::cout<<"Value = "<<CppAD::Value(taylor_y[0])<<" which should equal "<<Rf_lgammafn(CppAD::Value(taylor_x[0]))<<"\n";
     
     CppAD::AD<double> fprime;
     if(order_low <= 2 & order_up >= 1)
	  fprime = nimDerivs_lgammafn(taylor_x[0], baseOrder+1);
     if(order_low <= 1 & order_up >= 1) {
	  taylor_y[1] = fprime * taylor_x[1];
     }
     if(order_low <= 2 & order_up >= 2) {
	  taylor_y[2] = 0.5 * (nimDerivs_lgammafn(taylor_x[0], baseOrder+2) * taylor_x[1] * taylor_x[1] + 
			       fprime * 2 * taylor_x[2]);
     }
     return true;
}

bool atomic_lgamma_class::reverse(
     const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
     const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
     CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
     const CppAD::vector< CppAD::AD<double> >&               partial_y   ) 
{
  // std::cout<<"lgamma meta-reverse baseOrder = "<<baseOrder<<std::endl;
  //  printf("In lgamma meta-reverse for order_up %lu\n", order_up);
  partial_x[0] = CppAD::AD<double>(0);
     if(order_up >= 1) partial_x[1] = CppAD::AD<double>(0);
     CppAD::AD<double> fprime = nimDerivs_lgammafn(taylor_x[0], baseOrder+1);
     if(order_up >= 1) {
	  partial_x[1] += partial_y[1] * fprime;
	  partial_x[0] += partial_y[1] * nimDerivs_lgammafn(taylor_x[0], baseOrder+2) * taylor_x[1];
     }
     partial_x[0] += partial_y[0] * fprime;
     //     std::cout<<"fprime = "<<CppAD::Value(fprime)<<" partial_y[0] = "<<CppAD::Value(partial_y[0])<<" partial_x[0] = "<<CppAD::Value(partial_x[0])<<" partial_y[1] = "<<CppAD::Value(partial_y[1])<<" partial_x[1] = "<<CppAD::Value(partial_x[1])<<"\n";

     return true;
}
