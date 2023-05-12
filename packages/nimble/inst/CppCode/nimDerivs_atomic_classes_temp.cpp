#include <nimble/nimDerivs_atomic_classes.h>

/*
  The following tricky/awkward scheme for creating and deleting
  our custom atomic objects is set up to ensure they are created and deleted
  in the same DLL.
*/

#define ATOMIC_NEW_DELETE_(NAME) \
  atomic_##NAME##_class* new_atomic_##NAME(void *tape_mgr_ptr, const std::string& name) {\
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->new_atomic_##NAME(name);\
  }									\
  void delete_atomic_##NAME(void* tape_mgr_ptr, atomic_##NAME##_class *atomic_obj) { \
    reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->delete_atomic_##NAME(atomic_obj);	\
  }									\
  atomic_##NAME##_class* nimble_CppAD_tape_mgr::new_atomic_##NAME(const std::string& name) { \
    atomic_##NAME##_class* obj = new atomic_##NAME##_class(name);	\
    return obj;								\
   }									\
  void nimble_CppAD_tape_mgr::delete_atomic_##NAME(atomic_##NAME##_class *atomic_obj) { \
    delete atomic_obj;							\
  }\

ATOMIC_NEW_DELETE_(gammafn)
ATOMIC_NEW_DELETE_(backsolve)
ATOMIC_NEW_DELETE_(cholesky)
ATOMIC_NEW_DELETE_(forwardsolve)
ATOMIC_NEW_DELETE_(matinverse)
ATOMIC_NEW_DELETE_(matmult)
ATOMIC_NEW_DELETE_(pow_int)
ATOMIC_NEW_DELETE_(zround)
ATOMIC_NEW_DELETE_(floor)
ATOMIC_NEW_DELETE_(ceil)
ATOMIC_NEW_DELETE_(ftrunc)
ATOMIC_NEW_DELETE_(nimRound)
ATOMIC_NEW_DELETE_(log_pow_int)
ATOMIC_NEW_DELETE_(zb_over_a)
ATOMIC_NEW_DELETE_(probit)
ATOMIC_NEW_DELETE_(iprobit)

atomic_lgamma_class* new_atomic_lgamma(void* tape_mgr_ptr, const std::string& name, int bO) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->new_atomic_lgamma(name, bO);
}

void delete_atomic_lgamma(void* tape_mgr_ptr, atomic_lgamma_class *atomic_lgamma) {
  reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->delete_atomic_lgamma(atomic_lgamma);
}

atomic_lgamma_class* nimble_CppAD_tape_mgr::new_atomic_lgamma(const std::string& name, int bO) {
  atomic_lgamma_class* obj = new atomic_lgamma_class(name, bO);
  //  std::cout<<"creating new atomic_lgamma "<<obj<<" in "<<this<<std::endl;
  return obj;
}
void nimble_CppAD_tape_mgr::delete_atomic_lgamma(atomic_lgamma_class *atomic_lgamma) {
  //  std::cout<<"deleting atomic_lgamma "<<atomic_lgamma<<" in "<<this<<std::endl;
  delete atomic_lgamma;
}

atomic_gammafn_class *track_atomic_gammafn(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_gammafn(vec_ptr);
}
atomic_gammafn_class *nimble_CppAD_tape_mgr::get_atomic_gammafn(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!gammafn_exists) {
    gammafn_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_gammafn("atomic_gammafn_managed"), vec_ptr) );
    gammafn_exists = true;
  }
  return dynamic_cast<atomic_gammafn_class*>(atomic_ptrs[gammafn_index].first);
}

atomic_pow_int_class *track_atomic_pow_int(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_pow_int(vec_ptr);
}
atomic_pow_int_class *nimble_CppAD_tape_mgr::get_atomic_pow_int(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!pow_int_exists) {
    pow_int_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_pow_int("atomic_pow_int_managed"), vec_ptr) );
    pow_int_exists = true;
  }
  return dynamic_cast<atomic_pow_int_class*>(atomic_ptrs[pow_int_index].first);
}

atomic_zround_class *track_atomic_zround(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_zround(vec_ptr);
}
atomic_zround_class *nimble_CppAD_tape_mgr::get_atomic_zround(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!zround_exists) {
    zround_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_zround("atomic_zround_managed"), vec_ptr) );
    zround_exists = true;
  }
  return dynamic_cast<atomic_zround_class*>(atomic_ptrs[zround_index].first);
}

atomic_floor_class *track_atomic_floor(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_floor(vec_ptr);
}
atomic_floor_class *nimble_CppAD_tape_mgr::get_atomic_floor(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!floor_exists) {
    floor_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_floor("atomic_floor_managed"), vec_ptr) );
    floor_exists = true;
  }
  return dynamic_cast<atomic_floor_class*>(atomic_ptrs[floor_index].first);
}

atomic_ceil_class *track_atomic_ceil(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_ceil(vec_ptr);
}
atomic_ceil_class *nimble_CppAD_tape_mgr::get_atomic_ceil(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!ceil_exists) {
    ceil_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_ceil("atomic_ceil_managed"), vec_ptr) );
    ceil_exists = true;
  }
  return dynamic_cast<atomic_ceil_class*>(atomic_ptrs[ceil_index].first);
}

atomic_ftrunc_class *track_atomic_ftrunc(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_ftrunc(vec_ptr);
}
atomic_ftrunc_class *nimble_CppAD_tape_mgr::get_atomic_ftrunc(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!ftrunc_exists) {
    ftrunc_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_ftrunc("atomic_ftrunc_managed"), vec_ptr) );
    ftrunc_exists = true;
  }
  return dynamic_cast<atomic_ftrunc_class*>(atomic_ptrs[ftrunc_index].first);
}

atomic_nimRound_class *track_atomic_nimRound(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_nimRound(vec_ptr);
}
atomic_nimRound_class *nimble_CppAD_tape_mgr::get_atomic_nimRound(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!nimRound_exists) {
    nimRound_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_nimRound("atomic_nimRound_managed"), vec_ptr) );
    nimRound_exists = true;
  }
  return dynamic_cast<atomic_nimRound_class*>(atomic_ptrs[nimRound_index].first);
}

atomic_probit_class *track_atomic_probit(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_probit(vec_ptr);
}
atomic_probit_class *nimble_CppAD_tape_mgr::get_atomic_probit(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!probit_exists) {
    probit_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_probit("atomic_probit_managed"), vec_ptr) );
    probit_exists = true;
  }
  return dynamic_cast<atomic_probit_class*>(atomic_ptrs[probit_index].first);
}

atomic_iprobit_class *track_atomic_iprobit(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_iprobit(vec_ptr);
}
atomic_iprobit_class *nimble_CppAD_tape_mgr::get_atomic_iprobit(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!iprobit_exists) {
    iprobit_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_iprobit("atomic_iprobit_managed"), vec_ptr) );
    iprobit_exists = true;
  }
  return dynamic_cast<atomic_iprobit_class*>(atomic_ptrs[iprobit_index].first);
}


atomic_log_pow_int_class *track_atomic_log_pow_int(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_log_pow_int(vec_ptr);
}
atomic_log_pow_int_class *nimble_CppAD_tape_mgr::get_atomic_log_pow_int(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!log_pow_int_exists) {
    log_pow_int_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_log_pow_int("atomic_log_pow_int_managed"), vec_ptr) );
    log_pow_int_exists = true;
  }
  return dynamic_cast<atomic_log_pow_int_class*>(atomic_ptrs[log_pow_int_index].first);
}

atomic_zb_over_a_class *track_atomic_zb_over_a(void* tape_mgr_ptr,
					   std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_zb_over_a(vec_ptr);
}
atomic_zb_over_a_class *nimble_CppAD_tape_mgr::get_atomic_zb_over_a(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(!zb_over_a_exists) {
    zb_over_a_index = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_zb_over_a("atomic_zb_over_a_managed"), vec_ptr) );
    zb_over_a_exists = true;
  }
  return dynamic_cast<atomic_zb_over_a_class*>(atomic_ptrs[zb_over_a_index].first);
}

/***********************************/

atomic_lgamma_class *track_atomic_lgamma(int baseOrder,
					 void* tape_mgr_ptr,
					 std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  return reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->get_atomic_lgamma(baseOrder, vec_ptr);
}

atomic_lgamma_class *nimble_CppAD_tape_mgr::get_atomic_lgamma(int baseOrder,
							      std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  if(baseOrder > 4) baseOrder = 4;
  if(!lgamma_exists[baseOrder]) {
    lgamma_index[baseOrder] = atomic_ptrs.size();
    atomic_ptrs.push_back(atomic_pair(new_atomic_lgamma("atomic_lgamma_managed", baseOrder), vec_ptr) );
    lgamma_exists[baseOrder] = true;
  }
  return dynamic_cast<atomic_lgamma_class*>(atomic_ptrs[lgamma_index[baseOrder] ].first);
}

/***********************************/

atomic_lgamma_class::atomic_lgamma_class(const std::string& name, int baseOrder_) : 
  unary_atomic_class<double>(name),
  baseOrder(baseOrder_),
  verbose(false)
{ }

atomic_lgamma_class::atomic_lgamma_class(const std::string& name, int baseOrder_, bool verbose_) : 
  unary_atomic_class<double>(name),
  baseOrder(baseOrder_),
  verbose(verbose_)
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
  if(verbose) {
    std::cout<<"lgamma forward baseOrder = "<<baseOrder<<" low="<<order_low<<" up="<<order_up<<" tx[0]="<<taylor_x[0]<<" type_x[0]="<<type_x[0]<<" need_y="<<need_y<<std::endl;
    std::cout<<"tape_id and handle:"<< CppAD::AD<double>::get_tape_id_nimble()<<" "<< CppAD::AD<double>::get_tape_handle_nimble()<<"\n";
    std::cout<<"atomic info:"<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<"\n";
  }
  if((order_low <= 0) && (order_up >= 0)) {
    if(baseOrder == 0)
      taylor_y[0] = Rf_lgammafn(taylor_x[0]);
    else
      taylor_y[0] = Rf_psigamma(taylor_x[0], baseOrder-1);
    if(verbose) {
      std::cout<<"taylor_y[0] "<<taylor_y[0]<<" ";
    }
  }
  double fprime;
  if((order_low <= 2) && (order_up >= 1)) {
    fprime = Rf_psigamma(taylor_x[0], baseOrder);
    if(verbose) std::cout<<"fprime "<<fprime<<" ";
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = fprime * taylor_x[1];
    if(verbose) std::cout<<"taylor_x[1] "<<taylor_x[1]<<" taylor_y[1] "<<taylor_y[1]<<" ";
	  // f'(x) (x')
  }
  if((order_low <= 2) && (order_up >= 2)) {
    taylor_y[2] = 0.5 * (Rf_psigamma(taylor_x[0], baseOrder+1) * taylor_x[1] * taylor_x[1] + 
			 fprime * 2 * taylor_x[2]);
    if(verbose) std::cout<<"taylor_x[2] "<<taylor_x[2]<<" taylor_y[2] "<<taylor_y[2]<<" ";
    // 0.5 * ((f''(x)) (x')^2 + 2*f'(x) (x''))
    // Note x'' is a taylor coeff so it is 0.5*2nd deriv
  }
  if(verbose) std::cout<<std::endl;
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
  if(verbose) {
    std::cout<<"lgamma reverse baseOrder = "<<baseOrder<<" up="<<order_up<<" tx[0]="<<taylor_x[0]<<" ty[0]="<<taylor_y[0]<<" py[0]="<<partial_y[0]<<" type_x[0]="<<type_x[0]<<std::endl;
    std::cout<<"tape_id and handle:"<< CppAD::AD<double>::get_tape_id_nimble()<<" "<< CppAD::AD<double>::get_tape_handle_nimble()<<"\n";
    std::cout<<"atomic info:"<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<"\n";
  }
  // std::cout<<"lgamma reverse baseOrder = "<<baseOrder<<std::endl;
  // std::cout<<"reverse "<<order_up<<" "<<taylor_x[0]<<" "<<partial_y[0]<<std::endl;
     partial_x[0] = 0;
     if(order_up >= 1) partial_x[1] = 0;
     double fprime = Rf_psigamma(taylor_x[0], baseOrder);
     if(verbose) std::cout<<"fprime "<<fprime<<" ";
     if(order_up >= 1) {
	  partial_x[1] += partial_y[1] * fprime;
	  partial_x[0] += partial_y[1] * Rf_psigamma(taylor_x[0], baseOrder+1) * taylor_x[1];
	  if(verbose) std::cout<<"partial_x[1] "<<partial_x[1]<<" first step of partial_x[0] "<<partial_x[0]<<" ";
     }
     partial_x[0] += partial_y[0] * fprime;
     if(verbose) std::cout<<"partial_x[0] "<<partial_x[0]<<std::endl;
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
  if(verbose) {
    std::cout<<"lgamma meta-forward baseOrder = "<<baseOrder<<" low="<<order_low<<" up="<<order_up<<" tx[0]="<<CppAD::Value(taylor_x[0])<<" type_x[0]="<<type_x[0]<<" need_y="<<need_y<<std::endl;
    std::cout<<"tape_id and handle:"<< CppAD::AD<double>::get_tape_id_nimble()<<" "<< CppAD::AD<double>::get_tape_handle_nimble()<<"\n";
    std::cout<<"atomic info:"<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<"\n";
  }
     if((order_low <= 0) && (order_up >= 0)) {
       taylor_y[0] = nimDerivs_lgammafn(taylor_x[0], baseOrder, verbose); // This puts it in the new tape being recorded
	  if(verbose) {
	    std::cout<<"taylor_y[0] "<<CppAD::Value(taylor_y[0])<<" ";
	  }
     }
     // std::cout<<"Constant(taylor_y[0]) = "<<CppAD::Constant(taylor_y[0])<<"\n";
     // std::cout<<"Value = "<<CppAD::Value(taylor_y[0])<<" which should equal "<<Rf_lgammafn(CppAD::Value(taylor_x[0]))<<"\n";
     
     CppAD::AD<double> fprime;
     if((order_low <= 2) && (order_up >= 1)) {
       fprime = nimDerivs_lgammafn(taylor_x[0], baseOrder+1, verbose);
	  if(verbose) std::cout<<"fprime "<<CppAD::Value(fprime)<<" ";
     }
     if((order_low <= 1) && (order_up >= 1)) {
	  taylor_y[1] = fprime * taylor_x[1];
	  if(verbose) std::cout<<"taylor_x[1] "<<CppAD::Value(taylor_x[1])<<" taylor_y[1] "<<CppAD::Value(taylor_y[1])<<" ";
     }
     if((order_low <= 2) && (order_up >= 2)) {
       taylor_y[2] = 0.5 * (nimDerivs_lgammafn(taylor_x[0], baseOrder+2, verbose) * taylor_x[1] * taylor_x[1] + 
			       fprime * 2 * taylor_x[2]);
	  if(verbose) std::cout<<"taylor_x[2] "<<CppAD::Value(taylor_x[2])<<" taylor_y[2] "<<CppAD::Value(taylor_y[2])<<" ";
     }
     if(verbose) std::cout<<std::endl;
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
  if(verbose) {
    std::cout<<"lgamma meta-reverse baseOrder = "<<baseOrder<<" up="<<order_up<<" tx[0]="<<taylor_x[0]<<" ty[0]="<<taylor_y[0]<<" py[0]="<<partial_y[0]<<" type_x[0]="<<type_x[0]<<std::endl;
    std::cout<<"tape_id and handle:"<< CppAD::AD<double>::get_tape_id_nimble()<<" "<< CppAD::AD<double>::get_tape_handle_nimble()<<"\n";
    std::cout<<"atomic info:"<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<"\n";
  }

  partial_x[0] = CppAD::AD<double>(0);
     if(order_up >= 1) partial_x[1] = CppAD::AD<double>(0);
     CppAD::AD<double> fprime = nimDerivs_lgammafn(taylor_x[0], baseOrder+1, verbose);
     if(verbose) std::cout<<"fprime "<<CppAD::Value(fprime)<<" ";
     if(order_up >= 1) {
	  partial_x[1] += partial_y[1] * fprime;
	  partial_x[0] += partial_y[1] * nimDerivs_lgammafn(taylor_x[0], baseOrder+2, verbose) * taylor_x[1];
	  if(verbose) std::cout<<"partial_x[1] "<<CppAD::Value(partial_x[1])<<" first step of partial_x[0] "<<CppAD::Value(partial_x[0])<<" ";
     }
     partial_x[0] += partial_y[0] * fprime;
     if(verbose) std::cout<<"partial_x[0] "<<CppAD::Value(partial_x[0])<<std::endl;

     //     std::cout<<"fprime = "<<CppAD::Value(fprime)<<" partial_y[0] = "<<CppAD::Value(partial_y[0])<<" partial_x[0] = "<<CppAD::Value(partial_x[0])<<" partial_y[1] = "<<CppAD::Value(partial_y[1])<<" partial_x[1] = "<<CppAD::Value(partial_x[1])<<"\n";

     return true;
}

CppAD::AD<double> nimDerivs_lgammafn(CppAD::AD<double> x, int baseOrder, bool verbose) {
  if(verbose) {
    return nimDerivs_lgammafn_verbose(x, baseOrder);
  }
  atomic_lgamma_class *atomic_lgamma;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_lgamma = new atomic_lgamma_class("nimDerivs_lgamma", baseOrder);
  } else {
    if(baseOrder > 4) {
      std::cout<<"Error: lgamma derivatives requested for higher order than supported. "<<std::endl;
      baseOrder = 4;
    }
    atomic_lgamma = track_atomic_lgamma(baseOrder,
					CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  //  void *tape_mgr = CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr();

  // We can't use statics for reasons of multiple DLLs and bad interaction with CppAD statics.
 
  CppAD::vector<CppAD::AD<double>> in(1);
  CppAD::vector<CppAD::AD<double>> out(1);
  in[0] = x;
  (*atomic_lgamma)(in, out);
   if(!recording) {
    delete atomic_lgamma;
  }
  /* if(CppAD::AD<double>::get_tape_handle_nimble() == nullptr) { */
  /*   delete_atomic_lgamma(tape_mgr, atomic_lgamma); */
  /* } else { */
  /*   track_nimble_atomic(atomic_lgamma, */
  /* 			CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(), */
  /* 			CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() ); */
  /* } */
  return out[0];
}

CppAD::AD<double> nimDerivs_lgammafn(CppAD::AD<double> x) {
  return nimDerivs_lgammafn(x, 0); // a relic of a problem with default value when writing this previously using templates.
}

/***************************************/
/* gammafn                             */
/* We can't just use exp(lgammafn(x))  */
/* because lgamma is defined as log of */
/* ABSOLUTE VALUE of gamma(x)          */
/***************************************/

atomic_gammafn_class::atomic_gammafn_class(const std::string& name) : 
  unary_atomic_class<double>(name)
{ }

bool atomic_gammafn_class::forward(
     const CppAD::vector<double>&               parameter_x  ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
     size_t                              need_y       ,
     size_t                              order_low    ,
     size_t                              order_up     ,
     const CppAD::vector<double>&               taylor_x     ,
     CppAD::vector<double>&                     taylor_y     )
{
  bool verbose = false;
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = Rf_gammafn(taylor_x[0]);
  }
  double lgamma_prime;
  if((order_low <= 2) && (order_up >= 1)) {
    lgamma_prime = Rf_psigamma(taylor_x[0], 0) ;
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = lgamma_prime * taylor_y[0] * taylor_x[1];
    if(verbose) std::cout<<"taylor_x[1] "<<taylor_x[1]<<" taylor_y[1] "<<taylor_y[1]<<" ";
	  // f'(x) (x')
  }
  // This is never used from nimble
  if((order_low <= 2) && (order_up >= 2)) {
    taylor_y[2] = 0.5 * taylor_y[0] *
      ((Rf_psigamma(taylor_x[0], 1) + lgamma_prime*lgamma_prime) * taylor_x[1] * taylor_x[1] + 
       lgamma_prime * 2 * taylor_x[2]);
    if(verbose) std::cout<<"taylor_x[2] "<<taylor_x[2]<<" taylor_y[2] "<<taylor_y[2]<<" ";
    // 0.5 * ((f''(x)) (x')^2 + 2*f'(x) (x''))
    // Note x'' is a taylor coeff so it is 0.5*2nd deriv
  }
  if(verbose) std::cout<<std::endl;
  return true;
}
bool atomic_gammafn_class::reverse(
     const CppAD::vector<double>&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector<double>&               taylor_x    ,
     const CppAD::vector<double>&               taylor_y    ,
     CppAD::vector<double>&                     partial_x   ,
     const CppAD::vector<double>&               partial_y   )
{
  bool verbose = false;
  partial_x[0] = 0;
  if(order_up >= 1) partial_x[1] = 0;
  double lgamma_prime = Rf_psigamma(taylor_x[0], 0);
  double fprime = lgamma_prime * taylor_y[0];
  if(order_up >= 1) {
    partial_x[1] += partial_y[1] * fprime;
    double fprimeprime = (Rf_psigamma(taylor_x[0], 1) + lgamma_prime * lgamma_prime) * taylor_y[0];
    partial_x[0] += partial_y[1] * fprimeprime * taylor_x[1];
    if(verbose) std::cout<<"partial_x[1] "<<partial_x[1]<<" first step of partial_x[0] "<<partial_x[0]<<" ";
  }
  partial_x[0] += partial_y[0] * fprime;
  if(verbose) std::cout<<"partial_x[0] "<<partial_x[0]<<std::endl;
  return true;
}

bool atomic_gammafn_class::forward(
     const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
     size_t                              need_y       ,
     size_t                              order_low    ,
     size_t                              order_up     ,
     const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
     CppAD::vector< CppAD::AD<double> >&                     taylor_y     )
{
  bool verbose = false;
  if((order_low <= 0) && (order_up >= 0)) {
    taylor_y[0] = nimDerivs_gammafn(taylor_x[0]); // This puts it in the new tape being recorded
    if(verbose) {
      std::cout<<"taylor_y[0] "<<CppAD::Value(taylor_y[0])<<" ";
    }
  }
  CppAD::AD<double> lgamma_prime;
  if((order_low <= 2) && (order_up >= 1)) {
    lgamma_prime = nimDerivs_lgammafn(taylor_x[0], 1);
  }
  if((order_low <= 1) && (order_up >= 1)) {
    taylor_y[1] = lgamma_prime * taylor_y[0] * taylor_x[1];
    if(verbose) std::cout<<"taylor_x[1] "<<CppAD::Value(taylor_x[1])<<" taylor_y[1] "<<CppAD::Value(taylor_y[1])<<" ";
  }
  if((order_low <= 2) && (order_up >= 2)) {
    taylor_y[2] = 0.5 * taylor_y[0] *
      ((nimDerivs_lgammafn(taylor_x[0], 2) + lgamma_prime*lgamma_prime) * taylor_x[1] * taylor_x[1] + 
       lgamma_prime * 2 * taylor_x[2]);
    if(verbose) std::cout<<"taylor_x[2] "<<CppAD::Value(taylor_x[2])<<" taylor_y[2] "<<CppAD::Value(taylor_y[2])<<" ";
  }
  if(verbose) std::cout<<std::endl;
  return true;
}

bool atomic_gammafn_class::reverse(
     const CppAD::vector< CppAD::AD<double> >&               parameter_x ,
     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
     size_t                              order_up    ,
     const CppAD::vector< CppAD::AD<double> >&               taylor_x    ,
     const CppAD::vector< CppAD::AD<double> >&               taylor_y    ,
     CppAD::vector< CppAD::AD<double> >&                     partial_x   ,
     const CppAD::vector< CppAD::AD<double> >&               partial_y   ) 
{
  bool verbose = false;
  partial_x[0] = CppAD::AD<double>(0);
  if(order_up >= 1) partial_x[1] = CppAD::AD<double>(0);
  CppAD::AD<double> lgamma_prime = nimDerivs_lgammafn(taylor_x[0], 1);
  CppAD::AD<double>  fprime = lgamma_prime * taylor_y[0];
  if(verbose) std::cout<<"fprime "<<CppAD::Value(fprime)<<" ";
  if(order_up >= 1) {
    partial_x[1] += partial_y[1] * fprime;
    CppAD::AD<double> fprimeprime = (nimDerivs_lgammafn(taylor_x[0], 2) + lgamma_prime * lgamma_prime) * taylor_y[0];
    partial_x[0] += partial_y[1] * fprimeprime * taylor_x[1];
    if(verbose) std::cout<<"partial_x[1] "<<CppAD::Value(partial_x[1])<<" first step of partial_x[0] "<<CppAD::Value(partial_x[0])<<" ";
  }
  partial_x[0] += partial_y[0] * fprime;
  if(verbose) std::cout<<"partial_x[0] "<<CppAD::Value(partial_x[0])<<std::endl;
  return true;
}

CppAD::AD<double> nimDerivs_gammafn(CppAD::AD<double> x) {
  atomic_gammafn_class *atomic_gammafn;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_gammafn = new atomic_gammafn_class("nimDerivs_gammafn");
  } else {
    atomic_gammafn = track_atomic_gammafn(CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
					  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  } 
  CppAD::vector<CppAD::AD<double>> in(1);
  CppAD::vector<CppAD::AD<double>> out(1);
  in[0] = x;
  (*atomic_gammafn)(in, out);
  if(!recording) {
    delete atomic_gammafn;
  }
  return out[0];
}

