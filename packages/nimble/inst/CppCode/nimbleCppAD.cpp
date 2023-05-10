#include <nimble/nimbleCppAD.h>

#define USE_CPPAD_OPTIMIZE_FOR_MODEL_TAPES // comment this out to turn off atomics for nimDerivs(model$calculate(...),...)

void copy_CppADdouble_to_double(CppAD::AD<double> *first, CppAD::AD<double> *last, double *output) {
  CppAD::AD<double> *orig;
  double *result = output;
  for(orig = first; orig != last ; )
    *result++ = CppAD::Value(*orig++);
}


void copy_CppADdouble_to_double(CppAD::AD<double> &from, double &to) {
  to = CppAD::Value(from);
}

void copy_CppADdouble_to_double(NimArrBase< CppAD::AD<double> > &from, NimArrBase< double > &to) {
  to.setSize(from.getSizeVec(), false, false);
  copy_CppADdouble_to_double(from.getPtr(), from.getPtr() + from.size(), to.getPtr());
}

void setSizeGeneric( ) {

}


#ifdef _TIME_AD_GENERAL // need also _TIME_AD in nimbleCppAD.h
ad_timer derivs_main_timer("derivs_main");
ad_timer derivs_calc_timer("derivs_calc");
ad_timer derivs_node_iteration_timer("derivs_node_iteration");
ad_timer derivs_getDerivs_timer("derivs_getDerivs");
ad_timer derivs_run_tape_timer("derivs_run_tape");
int id(0);
// Include the following line to see node by node timing timing
// results within use of getDerivs and use of CppAD tapes.  This
// generates a lot of output.
//#define _SHOW_NODE_BY_NODE 
void derivs_getDerivs_timer_start() {derivs_getDerivs_timer.start(false);}
void derivs_run_tape_timer_start() {derivs_run_tape_timer.start(false);}
#ifdef _SHOW_NODE_BY_NODE 
void derivs_getDerivs_timer_stop() {derivs_getDerivs_timer.stop(true);}
void derivs_run_tape_timer_stop() {derivs_run_tape_timer.stop(true);}
void derivs_show_id() {std::cout<<"(id "<<id<<")"<<std::endl;}
#else
void derivs_getDerivs_timer_stop() {derivs_getDerivs_timer.stop(false);}
void derivs_run_tape_timer_stop() {derivs_run_tape_timer.stop(false);}
void derivs_show_id() {}
#endif
  void derivs_tick_id() {++id;}

SEXP report_AD_timers() {
  derivs_main_timer.show_report();
  derivs_calc_timer.show_report();
  derivs_node_iteration_timer.show_report();
  derivs_getDerivs_timer.show_report();
  derivs_run_tape_timer.show_report();
  return(R_NilValue);
}

SEXP reset_AD_timers(SEXP SreportInterval) {
  derivs_main_timer.reset();
  derivs_calc_timer.reset();
  derivs_node_iteration_timer.reset();
  derivs_getDerivs_timer.reset();
  derivs_run_tape_timer.reset();
  derivs_main_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_calc_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_node_iteration_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_getDerivs_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_run_tape_timer.set_interval(INTEGER(SreportInterval)[0]);
  return(R_NilValue);
}

#endif

#define _DERIVS_FULLTAPE
#ifdef  _DERIVS_FULLTAPE

nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
						 NodeVectorClassNew_derivs &nodes,
						 const NimArr<1, double> &derivOrders) {
  // std::cout<<"need to propagate const-ness"<<std::endl;
  std::cout<<"CALLING A FUNCTION THAT WE THOUGHT COULD BE DEPRECATED.  PLEASE REPORT TO PERRY."<<std::endl;
  
  if(!nodes.tapeRecorded()) nodes.recordTape();
#ifdef _TIME_AD_GENERAL
  derivs_main_timer.start();
#endif

  nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
#ifdef _TIME_AD_GENERAL
  derivs_getDerivs_timer.start();
#endif
  std::vector<double> independentVars;
  nodes.runTape_setIndependent(independentVars);
  std::vector<double> dependentVars;

#ifdef _TIME_AD_GENERAL
  derivs_run_tape_timer.start();
#endif
  
  nodes.runTape_runTape(independentVars, dependentVars,
			derivOrders, ansList);
#ifdef _TIME_AD_GENERAL
  derivs_run_tape_timer.stop();
#endif

#ifdef _TIME_AD_GENERAL
  derivs_getDerivs_timer.stop();
#endif

#ifdef _TIME_AD_GENERAL
  derivs_main_timer.stop();
#endif

  return ansList;
}

#endif

nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
    NodeVectorClassNew_derivs &nodes, const double derivOrders) {
  std::cout<<"CALLING A FUNCTION THAT WE THOUGHT COULD BE DEPRECATED.  PLEASE REPORT TO PERRY."<<std::endl;
  NimArr<1, double> orders(1);
  orders[0] = derivOrders;
  return (nimDerivs_calculate(nodes, orders));
}



nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
    NodeVectorClassNew_derivs &nodes, int iNodeFunction,
    NimArr<1, double> &derivOrders) {
  std::cout<<"CALLING A FUNCTION THAT WE THOUGHT COULD BE DEPRECATED.  PLEASE REPORT TO PERRY."<<std::endl;
  nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
  return (ADlist);
}

void setOrdersFound(const NimArr<1, double> &derivOrders,
		    bool *ordersFound,
		    int &maxOrder) {
  /* ordersFound must have length (at least) 3.*/
  int orderSize = derivOrders.size();
  double const* array_derivOrders = derivOrders.getConstPtr();
  maxOrder = 0;
  maxOrder =
    *std::max_element(array_derivOrders, array_derivOrders + orderSize);
  std::fill(ordersFound, ordersFound + 3, false);
  for (int i = 0; i < orderSize; i++) {
    if ((array_derivOrders[i] > 2) | (array_derivOrders[i] < 0)) {
      printf("Error: Derivative orders must be between 0 and 2.\n");
    }
    ordersFound[static_cast<int>(array_derivOrders[i])] = true;
  }
}


void update_dynamicVars(NodeVectorClassNew_derivs &NV,
			nimbleCppADinfoClass &ADinfo) {
  // This version may become deprecated
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  NimArr<1, double> NimArrValues;
  if(length_extraInput > 0) {
    NimArrValues.setSize(length_extraInput);
    ADinfo.dynamicVars.resize(length_extraInput);
    getValues(NimArrValues, NV.model_extraInput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_extraInput,
	       ADinfo.dynamicVars.begin());
    std::cout<<"calling new_dynamic"<<std::endl;
    ADinfo.ADtape()->new_dynamic(ADinfo.dynamicVars);
    std::cout<<"done in new_dynamic"<<std::endl;
  }
}

void update_dynamicVars(nimbleCppADinfoClass &ADinfo) {
  if(!ADinfo.updaterNV()) return;
  // copy from model --> dynamicVars and call new_dynamic
  NodeVectorClassNew_derivs &NV = *(ADinfo.updaterNV());
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  NimArr<1, double> NimArrValues;
  if(length_extraInput > 0) {
    NimArrValues.setSize(length_extraInput);
    ADinfo.dynamicVars.resize(length_extraInput);
    getValues(NimArrValues, NV.model_extraInput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_extraInput,
	       ADinfo.dynamicVars.begin());
    //    std::cout<<"calling new_dynamic"<<std::endl;
    ADinfo.ADtape()->new_dynamic(ADinfo.dynamicVars);
    //    std::cout<<"done in new_dynamic"<<std::endl;
  }
}

void update_dynamicVars_meta(nimbleCppADinfoClass &ADinfo) {
  if(!ADinfo.updaterNV()) return;
  // std::cout<<"in update_dynamicVars_meta"<<std::endl;
  // copy from model --> dynamicVars and call new_dynamic
  NodeVectorClassNew_derivs &NV = *(ADinfo.updaterNV());
  int length_extraInput = NV.model_AD_extraInput_accessor.getTotalLength();
  NimArr<1, CppAD::AD<double> > NimArrValuesAD;
  if(length_extraInput > 0) {
    NimArrValuesAD.setSize(length_extraInput);
    ADinfo.dynamicVars_meta.resize(length_extraInput);
    getValues_AD_AD(NimArrValuesAD, NV.model_AD_extraInput_accessor);
    std::copy( NimArrValuesAD.getPtr(),
	       NimArrValuesAD.getPtr() + length_extraInput,
	       ADinfo.dynamicVars_meta.begin());
  }  
  // call to new_dynamic is inside getDerivs_meta after base2ad.
}

void update_dynamicVars_meta(NodeVectorClassNew_derivs &NV,
			     nimbleCppADinfoClass &ADinfo) {
  // This version may become deprecated
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  NimArr<1, double> NimArrValues;
  if(length_extraInput > 0) {
    NimArrValues.setSize(length_extraInput);
    ADinfo.dynamicVars.resize(length_extraInput);
    getValues(NimArrValues, NV.model_extraInput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_extraInput,
	       ADinfo.dynamicVars.begin());
  }
}

bool check_inf_nan_gdi(double v) {
  if(((v == -std::numeric_limits<double>::infinity()) ||
      (v == std::numeric_limits<double>::infinity())) ||
     (std::isnan(v))) {
    return true;
  }
  return false;
}

bool check_inf_nan_gdi(CppAD::AD<double> v) {
  return false;
}

template<typename BASE, class TAPETYPE, class ADCLASS>
void getDerivs_internal(vector<BASE> &independentVars,			
			TAPETYPE *ADtape,
			const NimArr<1, double> &derivOrders,
			const NimArr<1, double> &wrtVector,
			nimSmartPtr<ADCLASS> &ansList) {
  // std::cout<<"entering getDerivs_internal"<<std::endl;
  // std::cout<<"independentVars: ";
  // for(size_t ijk = 0; ijk < independentVars.size(); ++ijk)
  //   std::cout<<independentVars[ijk]<<" ";
  // std::cout<<std::endl;
  // std::cout<<"derivOrders: ";
  // for(size_t ijk = 0; ijk < derivOrders.size(); ++ijk)
  //   std::cout<<derivOrders[ijk]<<" ";
  // std::cout<<std::endl;
  // std::cout<<"wrtVector: ";
  // for(size_t ijk = 0; ijk < wrtVector.size(); ++ijk)
  //   std::cout<<wrtVector[ijk]<<" ";
  // std::cout<<std::endl;
  
#ifdef _TIME_AD_GENERAL
  derivs_getDerivs_timer_start();
  derivs_tick_id();
  derivs_show_id();  
#endif
  std::size_t n = independentVars.size();  // dim of independent vars

  std::size_t wrt_n = wrtVector.size();            // dim of wrt vars
  if(wrt_n == 2){
    if(wrtVector[1] == -1){ // 2nd element -1 is a filler to ensure there is a vector out of compilation
      wrt_n = 1;
    }
  }
  bool wrtAll = wrtVector[0] == -1; // 1st element -1 is a flag to behave as if wrtVector has all elements
  if(wrtAll) wrt_n = n;

  int maxOrder;
  bool ordersFound[3];
  setOrdersFound(derivOrders, ordersFound, maxOrder);
  //  std::cout<<"orders: "<<ordersFound[0]<<" "<<ordersFound[1]<<" "<<ordersFound[2]<<" "<<maxOrder<<std::endl;
  // std::cout<<"maxOrder = "<<maxOrder<<std::endl;
  vector<BASE> value_ans;
  #ifdef _TIME_AD_GENERAL
  derivs_run_tape_timer_start();
#endif
  value_ans = ADtape->Forward(0, independentVars);
  //  std::cout<<"value_ans.size() = "<<value_ans.size()<<std::endl;
#ifdef _TIME_AD_GENERAL
  derivs_run_tape_timer_stop();
#endif
  if (ordersFound[0]) {
    ansList->value.setSize(value_ans.size(), false, false);
    std::copy(value_ans.begin(), value_ans.end(), ansList->value.getPtr());
  }
  if(maxOrder > 0){
    std::size_t q = value_ans.size();
    vector<bool> infIndicators(q, false); // default values will be false 
    for(size_t inf_ind = 0; inf_ind < q; inf_ind++){
      if(check_inf_nan_gdi(value_ans[inf_ind])) {
	infIndicators[inf_ind] = true;
      }
    }
    if (ordersFound[1]) {
      ansList->jacobian.setSize(q, wrt_n, false, false); 
    }
    if (ordersFound[2]) {
      ansList->hessian.setSize(wrt_n, wrt_n, q, false, false);
    }
    vector<BASE> cppad_derivOut;
    std::vector<BASE> w(q, 0);
    for (size_t dy_ind = 0; dy_ind < q; dy_ind++) {
      w[dy_ind] = 1;
      if (maxOrder == 1) {   
	if(!infIndicators[dy_ind]){
#ifdef _TIME_AD_GENERAL
	  derivs_run_tape_timer_start();
#endif
	  cppad_derivOut = ADtape->Reverse(1, w);
#ifdef _TIME_AD_GENERAL
	  derivs_run_tape_timer_stop();
#endif
	}
      } else {
	//	std::cout<<"wrtAll = "<<wrtAll<<std::endl;
	// if(!wrtAll) {
	//   for(size_t ijk = 0; ijk < wrt_n; ijk++) {
	//     std::cout<<wrtVector[ijk]<<" ";
	//   }
	//   std::cout<<std::endl;
	// }
	for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
	  if(!infIndicators[dy_ind]){
	    int dx1_ind = wrtAll ? vec_ind : wrtVector[vec_ind] - 1;
	    std::vector<BASE> x1(n, 0);  // vector specifying first derivatives.
	    // first specify coeffs for first dim
	    // of s across all directions r, then
	    // second dim, ...
	    x1[dx1_ind] = 1;
#ifdef _TIME_AD_GENERAL
	    derivs_run_tape_timer_start();
#endif
	    // std::cout<<"Forward 1 x1: ";
	    // for(int ijk = 0; ijk < x1.size(); ++ijk)
	    //   std::cout<<x1[ijk]<<" ";
	    // std::cout<<std::endl;	    
	    // vector<BASE> forwardOut;
	    // forwardOut = ADtape->Forward(1, x1);
	    ADtape->Forward(1, x1);
	    //	    std::cout<<"forwardOut 1 result (dx1_ind = "<< dx1_ind << ", forwardOut.size() = "<< forwardOut.size() <<"): ";
	    // for(int ijk = 0; ijk < forwardOut.size(); ++ijk)
	    //   std::cout<<forwardOut[ijk]<<" ";
	    // std::cout<<std::endl;
	    cppad_derivOut = ADtape->Reverse(2, w);
	    // std::cout<<"reverse 2 result: ";
	    // for(int ijk = 0; ijk < cppad_derivOut.size(); ++ijk)
	    //   std::cout<<cppad_derivOut[ijk]<<" ";
	    // std::cout<<std::endl;

#ifdef _TIME_AD_GENERAL
	    derivs_run_tape_timer_stop();
#endif
	  }
	  for (size_t vec_ind2 = 0; vec_ind2 < wrt_n; vec_ind2++) {
	    if(!infIndicators[dy_ind]){
	      int dx2_ind = wrtAll ? vec_ind2 : wrtVector[vec_ind2] - 1;
	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] =
		cppad_derivOut[dx2_ind * 2 + 1];
	    }
	    else{
	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] = 
		CppAD::numeric_limits<BASE>::quiet_NaN();
	    }
	  }
	}
      }
      if (ordersFound[1]) {
	BASE *LHS = ansList->jacobian.getPtr() + dy_ind;
	if(!infIndicators[dy_ind]){
	  if(wrtAll) {
	    for (size_t vec_ind3 = 0; vec_ind3 < wrt_n; ++vec_ind3, LHS += q) {
	      *LHS = cppad_derivOut[vec_ind3 * maxOrder];
	    }
	  } else {
	    double const *wrtVector_p = wrtVector.getConstPtr();
	    double const *wrtVector_p_end = wrtVector_p + wrt_n;
	    for(; wrtVector_p != wrtVector_p_end; LHS += q ) {
	      *LHS = cppad_derivOut[(static_cast<int>(*wrtVector_p++) - 1) * maxOrder];
	    }
	  }
	} else {
	  for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
	    *LHS = CppAD::numeric_limits<BASE>::quiet_NaN();
	    LHS += q;
	  }
	}

	// for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
	//   if(!infIndicators[dy_ind]){
	//     int dx1_ind = wrtVector[vec_ind] - 1;
	//     ansList->jacobian[vec_ind * q + dy_ind] =
	//       cppad_derivOut[dx1_ind * maxOrder + 0];
	//   }
	//   else{
	//     ansList->jacobian[vec_ind * q + dy_ind] =
	//       CppAD::numeric_limits<double>::quiet_NaN();
	//   }     
	// }
	
      }
      w[dy_ind] = 0;
    }
  }
#ifdef _TIME_AD_GENERAL
  derivs_getDerivs_timer_stop();
#endif
};

void nimbleFunctionCppADbase::getDerivs_meta(nimbleCppADinfoClass &ADinfo,
					     const NimArr<1, double> &derivOrders,
					     const NimArr<1, double> &wrtVector,
					     const nimbleCppADrecordingInfoClass &nimRecInfo,
					     nimSmartPtr<NIMBLE_ADCLASS_META> &ansList) {
  //  std::cout<<"Entering getDerivs_meta"<<std::endl;
  //  std::cout<<"ADinfo is at :"<< &ADinfo <<"\n";

  //  if(!nimRecInfo.recording_cp()) return;
  
  bool orderIncludesZero(false);
  for(int i = 0; i < derivOrders.size(); ++i) {
    orderIncludesZero |= (derivOrders[i] == 0);
  }
  // std::cout << "orderIncludesZero = " << orderIncludesZero << std::endl;
  bool oldUpdateModel = ADinfo.updateModel();
  ADinfo.updateModel() = orderIncludesZero;

  // We need to use the tricks to have CppAD statics in (potentially) two compilation units
  // that were not linked together (one for model, one for the current algorithm)
  // be matching.  This makes it so the CppAD tape handle and atomic information are shared.
  // This is necessary even in this double taping step because of our atomic classes,
  // which might reside in the other compilation unit.  During double taping, an atomic
  // such as lgamma puts itself or other statics onto the new tape, and returns
  // CppAD::AD variables created in the other compilation unit.
  set_CppAD_tape_info_for_model my_tape_info_RAII_; // must be at function scope, not declared inside if(){}
  
  if(ADinfo.nodeFunPtrSet()) {
    //    std::cout<<"tape_id and handle:"<<  nimRecInfo.tape_id_cp() <<" "<< nimRecInfo.tape_handle_cp() <<"\n";
    //   std::cout<<"atomic info:"<<nimRecInfo.atomic_vec_ptr_cp()<<"\n";
    my_tape_info_RAII_.set_from_nodeFunPtr(ADinfo.nodeFunPtr(),
					   nimRecInfo.tape_id_cp(), //CppAD::AD<double>::get_tape_id_nimble(),
					   nimRecInfo.tape_handle_cp());//CppAD::AD<double>::get_tape_handle_nimble());
    set_CppAD_atomic_info_for_model(ADinfo.nodeFunPtr(),
				    nimRecInfo.atomic_vec_ptr_cp());
    //    std::cout<<"done setting nodeFunPtr\n";
  }

  CppAD::ADFun< CppAD::AD<double>, double > innerTape;
  innerTape = ADinfo.ADtape()->base2ad();
  innerTape.new_dynamic(ADinfo.dynamicVars_meta);

  //  std::cout<<" after making inner tape\n";
  //  std::cout<<"tape_id and handle:"<< CppAD::AD<double>::get_tape_id_nimble()<<" "<< CppAD::AD<double>::get_tape_handle_nimble()<<"\n";
  //  std::cout<<"atomic info:"<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<"\n";
  getDerivs_internal< CppAD::AD<double>,
		      CppAD::ADFun< CppAD::AD<double>, double >,
		      NIMBLE_ADCLASS_META>(ADinfo.independentVars_meta,
					   &innerTape,
					   derivOrders,
					   wrtVector,
					   ansList);
  ADinfo.updateModel() = oldUpdateModel;
  //  std::cout<<"Exiting getDerivs_meta"<<std::endl;
}
  
void nimbleFunctionCppADbase::getDerivs(nimbleCppADinfoClass &ADinfo,
                                        const NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
  // std::cout<<"Entering getDerivs"<<std::endl;
  bool orderIncludesZero(false);
  for(int i = 0; i < derivOrders.size(); ++ i) {
    orderIncludesZero |= (derivOrders[i] == 0);
  }
  //  std::cout << "orderIncludesZero = " << orderIncludesZero << std::endl;
  bool oldUpdateModel = ADinfo.updateModel();
  ADinfo.updateModel() = orderIncludesZero;
  getDerivs_internal<double,
		     CppAD::ADFun<double>,
		     NIMBLE_ADCLASS>(ADinfo.independentVars,
				     ADinfo.ADtape(),
				     derivOrders,
				     wrtVector,
				     ansList);
  ADinfo.updateModel() = oldUpdateModel;
  //  std::cout<<"Exiting getDerivs"<<std::endl;
}


CppAD::ADFun<double>* calculate_recordTape(NodeVectorClassNew_derivs &NV,
					   bool includeExtraOutputs,
					   nimbleCppADinfoClass &ADinfo) {
  vector< CppAD::AD<double> > dependentVars(1);
  NimArr<1, double> NimArrValues;
  NimArr<1, CppAD::AD<double> > NimArrValues_AD;
  
  // 1. Copy all constantNodes values from model -> model_AD
  int length_constant = NV.model_constant_accessor.getTotalLength();
  if(length_constant > 0) {
    NimArr<1, double> NimArrValues;
    NimArr<1, CppAD::AD<double> > NimArrValues_AD;
    NimArrValues.setSize(length_constant);
    NimArrValues_AD.setSize(length_constant);
    getValues(NimArrValues, NV.model_constant_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_constant,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_constant_accessor);
  }

  // 2. Copy all wrtNodes values from model -> model_AD, AND
  // 3. Copy all wrtNodes values from model -> independentVars, AND
  // 4. [Deleted]
  int length_wrt = NV.model_wrt_accessor.getTotalLength();
  int length_independent = length_wrt;
  // std::cout<<"recording with length "<<length_independent<<std::endl;
  // std::cout<<" length_wrt = "<<length_wrt<<std::endl;
  vector< CppAD::AD<double> > independentVars(length_independent);
  if(length_wrt > 0) {
    NimArrValues.setSize(length_wrt);
    getValues(NimArrValues, NV.model_wrt_accessor);
    // 2
    NimArrValues_AD.setSize(length_wrt);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_wrt,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
    // 3
    std::copy(  NimArrValues.getPtr(),
		NimArrValues.getPtr() + length_wrt,
		independentVars.begin() );
  }
  // 5a. Copy all extraInputNodes values from model -> model_AD (ditto, may be redundant)
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  if(length_extraInput > 0) {
    NimArrValues.setSize(length_extraInput);
    NimArrValues_AD.setSize(length_extraInput);
    getValues(NimArrValues, NV.model_extraInput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_extraInput,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_extraInput_accessor);
  }
  // 5b. Copy all extraInputNodes into dynamicVars
  // std::cout<<"Don't forget to set the CppAD statics as needed"<<std::endl;
  vector< CppAD::AD<double> > dynamicVars;
  dynamicVars.resize(length_extraInput);
  if(length_extraInput > 0) {
    std::copy( NimArrValues_AD.getPtr(),
	       NimArrValues_AD.getPtr() + length_extraInput,
	       dynamicVars.begin() );
  }
  
  // 6. Start taping
  size_t abort_op_index = 0;    // per CppAD examples, these match CppAD default values
  bool   record_compare = true; // but must be provided explicitly to get to the dynamic parameter (4th) argument
  // std::cout<<"recording with "<<dynamicVars.size()<<std::endl;
  // std::cout<<"Before independent: tape handle address = "<< CppAD::AD<double>::get_handle_address_nimble() <<std::endl;
  CppAD::Independent(independentVars, abort_op_index, record_compare, dynamicVars);
  ADinfo.set_internal_tape(CppAD::AD<double>::get_tape_handle_nimble());
  //  std::cout<<"After independent: tape handle address = "<< CppAD::AD<double>::get_handle_address_nimble() <<std::endl;
  {
    set_CppAD_tape_info_for_model my_tape_info_RAII_(NV,
						     CppAD::AD<double>::get_tape_id_nimble(),
						     CppAD::AD<double>::get_tape_handle_nimble());
    set_CppAD_atomic_info_for_model(NV, CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage());
    
    // 7. [deleted]
    // 8. [deleted]
    
    // 9. Copy all wrtNodes AD objects from independentVars -> model_AD.
    if(length_extraInput > 0) {
      NimArrValues_AD.setSize(length_extraInput);    
      std::copy(dynamicVars.begin(),
		dynamicVars.begin() + length_extraInput,
		NimArrValues_AD.getPtr());
      setValues_AD_AD(NimArrValues_AD, NV.model_AD_extraInput_accessor);
    }
    if(length_wrt > 0) {
      NimArrValues_AD.setSize(length_wrt);
      std::copy(independentVars.begin(),
		independentVars.begin() + length_wrt,
		NimArrValues_AD.getPtr());
      setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
    }
    // 10. call calculate.  This also sets up the extraOutput step
    nimbleCppADrecordingInfoClass recordingInfo(true, &ADinfo);
    CppAD::AD<double> logProb = calculate_ADproxyModel(NV,
						       includeExtraOutputs, // if true, model will be updated from tape.
						       recordingInfo);
    dependentVars[0] = logProb;
    // 13. Finish taping, AND
    // 14. Call tape->optimize()
    // make it locally to get the right globals during recording and playback
    // DO NOT USE THE CONSTRUCTOR VERSION BECAUSE IT ALWAYS DOES .Forward(0)
    // INSTEAD MAKE THE BLANK OBJECT AND USE .Dependent(...)
    // TRY USING CppAD's vector type
  } // These {} ensure that the destructor for the my_tape_info_RAII_ is called before Dependent, which is necessary in some cases (depending on libnimble.a vs libnimble.so)
  CppAD::ADFun<double>* ansTape = new CppAD::ADFun<double>;
  ADinfo.sum_dummyOutputs_to_dependentVars(dependentVars);
  ansTape->Dependent(independentVars, dependentVars);
  //  std::cout<<"about to call optimize"<<std::endl;
  //  std::cout<<"tape handle address = "<< CppAD::AD<double>::get_handle_address_nimble() <<std::endl;
#ifdef USE_CPPAD_OPTIMIZE_FOR_MODEL_TAPES
  ansTape->optimize(); //("no_compare_op") makes almost no difference;
#endif
  // std::cout<<"done with optimize"<<std::endl;
  return ansTape;
}

void nimbleFunctionCppADbase::getDerivs_calculate_internal(nimbleCppADinfoClass &ADinfo,
							   // CppAD::ADFun<double>* &tapePtr,
							   NodeVectorClassNew_derivs &nodes,
							   const NimArr<1, double> &derivOrders,
							   const NimArr<1, double> &wrtVector,
                 bool do_update,
							   bool reset,
							   nimSmartPtr<NIMBLE_ADCLASS> ansList) {
  // If use_meta_tape is true, then double-recording will be used.
  // This means that a first tape will be recorded, then a second tape will be recorded of obtaining 1st order derivatives from the first tape.
  // When derivatives are requested, 0th order will be obtained from the regular model (not AD tape),
  // 1st order will be 0th order of the second tape, and 2nd order will be 1st order of the second tape.
  // If use_meta_tape is false, then the (single, first) tape will be used "directly" for 0th, 1st or 2nd order.
  using std::cout;
  using std::endl;
  bool use_meta_tape = true;
  //  cout<<"in getDerivs_calculate_internal"<<endl;
  // Record tape(s) if this is the first time or if reset is true.
  if(ADinfo.ADtape_empty() || reset) {
    // Delete previous tape if it exists.
    if(!ADinfo.ADtape_empty())
      ADinfo.ADtape_reset();
    if(!use_meta_tape) {
      ADinfo.ADtape() = calculate_recordTape(nodes, true, ADinfo); // sets internal tape for atomic tracking
    } else {
      CppAD::ADFun< double > *firstTape;
      firstTape = calculate_recordTape(nodes, false, ADinfo); // sets internal tape for atomic tracking
      CppAD::ADFun< CppAD::AD<double>, double > innerTape;
      // Make original tape use CppAD::AD<double> instead of double
      set_CppAD_atomic_info_for_model(nodes, CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage());
      innerTape = firstTape->base2ad();
      int length_wrt = nodes.model_wrt_accessor.getTotalLength();
      int length_independent = length_wrt;
      int length_extraInput = nodes.model_extraInput_accessor.getTotalLength();
      vector< CppAD::AD<double> > dependentVars(length_wrt); // This will be the jacobian from the first tape, i.e. value of the second tape
      vector< CppAD::AD<double> > independentVars(length_independent);
      vector< CppAD::AD<double> > dynamicVars(length_extraInput);
      size_t abort_op_index = 0;
      bool   record_compare = true;
      NimArr<1, double > NimArrVars;
      NimArrVars.setSize(length_wrt);
      // Initialize values of independentVars, before recording.
      getValues(NimArrVars, nodes.model_wrt_accessor);
      for(int iii = 0; iii < length_wrt; ++iii)
	independentVars[iii] = NimArrVars[iii];
      // std::copy(NimArrVars.getPtr(),
      // 		NimArrVars.getPtr() + length_wrt,
      // 		independentVars.begin());
      //
      // Initialize values of dynamicVars, before recording.
      if(length_extraInput > 0) {
	NimArr<1, double> NimArr_dynamicVars;
	NimArr_dynamicVars.setSize(length_extraInput);
	getValues(NimArr_dynamicVars, nodes.model_extraInput_accessor);
	for(int iii = 0; iii < length_extraInput; ++iii)
	  dynamicVars[iii] = NimArr_dynamicVars[iii];
	// std::copy( NimArr_dynamicVars.getPtr(),
	// 	   NimArr_dynamicVars.getPtr() + length_extraInput,
	// 	   dynamicVars.begin() );
      }
      nimSmartPtr<NIMBLE_ADCLASS_META> ansList_meta = new NIMBLE_ADCLASS_META;
      // start recording new (second) tape
      CppAD::Independent(independentVars, abort_op_index, record_compare, dynamicVars);
      ADinfo.set_internal_tape(CppAD::AD<double>::get_tape_handle_nimble());
      // Trick CppAD statics to work across nimble compilation units
      {
	set_CppAD_tape_info_for_model my_tape_info_RAII_(nodes,
							 CppAD::AD<double>::get_tape_id_nimble(),
							 CppAD::AD<double>::get_tape_handle_nimble());
	set_CppAD_atomic_info_for_model(nodes, CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage());
	// Set up inputs to first tape (recorded in second tape)
	ADinfo.independentVars_meta.resize(length_wrt);
	ADinfo.dynamicVars_meta.resize(length_extraInput);
	if(length_extraInput > 0) {
	  for(int iii = 0; iii < length_extraInput; ++iii)
	    ADinfo.dynamicVars_meta[iii] = dynamicVars[iii]; // std::copy does not seem to work for CppAD recording
	  // std::copy(dynamicVars.begin(),
	  // 	  dynamicVars.end(),
	  // 	  ADinfo.dynamicVars_meta.begin());
	}
	for(int iii = 0; iii < length_wrt; ++iii)
	  ADinfo.independentVars_meta[iii] = independentVars[iii];
	// std::copy(independentVars.begin(),
	// 		independentVars.end(),
	// 		ADinfo.independentVars_meta.begin());      
	
	NimArr<1, double> derivOrders_meta;
	derivOrders_meta.setSize(1);
	derivOrders_meta[0] = 1;
	// std::cout<<ADinfo.dynamicVars_meta.size()<<std::endl;
	innerTape.new_dynamic(ADinfo.dynamicVars_meta);
	
	getDerivs_internal< CppAD::AD<double>,
			    CppAD::ADFun< CppAD::AD<double>, double >,
			    NIMBLE_ADCLASS_META>(ADinfo.independentVars_meta,
						 &innerTape,
						 derivOrders_meta,
						 wrtVector,
						 ansList_meta);
	for(int iii = 0; iii < length_wrt; ++iii)
	  dependentVars[iii] = ansList_meta->jacobian[iii];
	ADinfo.ADtape() = new CppAD::ADFun<double>;
      } // These {} ensure the RAII object's destructor is called before Dependent, which is important on OS's (linux) with libnimble.so instead of libnimble.a
      ADinfo.ADtape()->Dependent(dependentVars);
#ifdef USE_CPPAD_OPTIMIZE_FOR_MODEL_TAPES
      ADinfo.ADtape()->optimize();
#endif
      delete firstTape;
    }
  }

  // Recording, if needed, is done.
  // From here on is use of the tape(s).  This may be used much more often than recording section above.
  
  // std::cout<<"getDerivs_calculate_internal A"<<std::endl;
  // Copy values from the model into the independentVars
  int length_wrt = nodes.model_wrt_accessor.getTotalLength();
  int length_independent = length_wrt;
  ADinfo.independentVars.resize(length_independent);
    
  NimArr<1, double > NimArrVars;
  NimArrVars.setSize(length_wrt);
  getValues(NimArrVars, nodes.model_wrt_accessor);
  
  std::copy(NimArrVars.getPtr(),
	    NimArrVars.getPtr() + length_wrt,
	    ADinfo.independentVars.begin());
  /* set dynamic */
  // Copy extraInput (CppAD "dynamic") values from the model into the dynamicVars
  // *and* set them in the tape.
  size_t length_extraNodes_accessor = nodes.model_extraInput_accessor.getTotalLength();
  if(length_extraNodes_accessor > 0) {
    NimArr<1, double> NimArr_dynamicVars;
    NimArr_dynamicVars.setSize(length_extraNodes_accessor);
    getValues(NimArr_dynamicVars, nodes.model_extraInput_accessor);
    std::vector<double> dynamicVars(length_extraNodes_accessor);
    std::copy( NimArr_dynamicVars.getPtr(),
	       NimArr_dynamicVars.getPtr() + length_extraNodes_accessor,
	       dynamicVars.begin() );
    //std::cout<<"Setting new_dynamic to"<<std::endl;
    //for(int ijk = 0; ijk < length_extraNodes_accessor; ijk++)
    //  std::cout<<dynamicVars[ijk]<<" ";
    //std::cout<<std::endl;
    ADinfo.ADtape()->new_dynamic(dynamicVars);
  }

  if(use_meta_tape) {
    // manage orders and use regular calculate for value
    // and tape for jacobian or hessian
    int maxOrder;
    bool ordersFound[3];
    setOrdersFound(derivOrders, ordersFound, maxOrder);
    if(ordersFound[0]) {
      ansList->value.setSize(1, false, false);
      ansList->value[0] = calculate(nodes);
    }
    NimArr<1, double> derivOrders_nested;
    int higherOrders = 0;
    if(ordersFound[1]) ++higherOrders;      
    if(ordersFound[2]) ++higherOrders;
    if(higherOrders) {
      derivOrders_nested.setSize(higherOrders, false, false);
      higherOrders = 0;
      if(ordersFound[1]) derivOrders_nested[higherOrders++] = 0; // If Jacobian was requested, get value of meta tape
      if(ordersFound[2]) derivOrders_nested[higherOrders] = 1; // If Hessian was requested, get Jacobian of meta tape
      nimSmartPtr<NIMBLE_ADCLASS> ansList_nested = new NIMBLE_ADCLASS;
      //std::cout<<"about to call getDerivs_internal"<<std::endl;
      getDerivs_internal<double,
			 CppAD::ADFun<double>,
			 NIMBLE_ADCLASS>(ADinfo.independentVars,
					 ADinfo.ADtape(),
					 derivOrders_nested,
					 wrtVector, // NOTE: This will not behave fully correctly in non-default use without further thought.
					 ansList_nested);
      if(ordersFound[1]) {
	ansList->jacobian.setSize(1, length_wrt, false, false);
	for(int ii = 0; ii < length_wrt; ++ii) //We could rewrite this with better pointer magic
	  ansList->jacobian[ ii ] = ansList_nested->value[ ii ];
      }
      if(ordersFound[2]) {
	ansList->hessian.setSize(length_wrt, length_wrt, 1, false, false);
	for(int ii = 0; ii < length_wrt; ++ii)
	  for(int jj = 0; jj < length_wrt; ++jj)
	    ansList->hessian[jj + ii*length_wrt ] = ansList_nested->jacobian[jj + ii*length_wrt]; //orientation shouldn't matter due to symmetry
      }
    }
  } else {
  /* run tape */
    //std::cout<<"running tape"<<std::endl;
    getDerivs_internal<double,
		       CppAD::ADFun<double>,
		       NIMBLE_ADCLASS>(ADinfo.independentVars,
				       ADinfo.ADtape(),
				       derivOrders,
				       wrtVector,
				       ansList);
  }
}

NimArr<1, double> make_vector_if_necessary(int a){
      NimArr<1, double> intArray;
      intArray.setSize(1, false, false);
      intArray[0] = a;
      return(intArray);
}

NimArr<1, double> make_vector_if_necessary(double a){
      NimArr<1, double> intArray;
      intArray.setSize(1, false, false);
      intArray[0] = a;
      return(intArray);
}

NimArr<1, double> make_vector_if_necessary(NimArr<1, double> a){
      return(a);
}

NimArr<1, double> make_vector_if_necessary(NimArr<1, int> a){
      NimArr<1, double> intArray;
      intArray.setSize(a.size(), false, false);
      std::copy(a.getPtr(),
		a.getPtr() + a.size(),
		intArray.getPtr());
      return(intArray);
}

void setValues_AD_AD_taping(NimArr<1, CppAD::AD<double> > &v,
			    ManyVariablesMapAccessor &MVA_AD,
			    ManyVariablesMapAccessor &MVA_orig,
			    nimbleCppADrecordingInfoClass &recordingInfo){
  size_t totalLength = MVA_orig.getTotalLength();
  if(!recordingInfo.recording()) {
    // This could cause confusion:
    // There are always two calls to functions with CppAD types.
    // The first call occurs with recording *off*.  This ensures that any nested taping occurs.  (The nested call with itself do both of its calls, resulting in a tape, as part of the first call of an outer [meta] taping process.)

    // At one stage of development, we always pulled the values out of the tape
    // and put them in the model during the non-recording (first) call.  The below code block did that.
    // Now we are turning that off.
    
    // We decided to never invoke updates of calculated values to the model during recording.
    // However, if order 0 is not involved, then the model should not be modified.

    // A potential concern is if updating the model ever plays a role in ensuring that it has valid
    // values during the taping call.  I am not sure if that really comes into play or not.
    // If it does, the user probably wants to run the function (or an appropriate model$calculate()
    // prior to taping (which happens via nimDerivs())

    // Code that previously ensured model updates from values()<- during the non-taping AD call.
    /*
    NimArr<1, double> dv;
    dv.setSize(totalLength);
    for(size_t ii = 0; ii < totalLength; ++ii) {
      dv[ii] = CppAD::Value(v[ii]);
    }
    setValues(dv, MVA_orig);
    */
  } else {
    // Make a new extraOutputObject, which during tape execution will copy to the real (non-AD) model
    std::vector< CppAD::AD<double> > extraOutputDummyResult(1);
    extraOutputDummyResult[0] = 0;
    std::vector< CppAD::AD<double> > extraOutputs(totalLength);
    for(size_t ii = 0; ii < totalLength; ++ii) {
      extraOutputs[ii] = v[ii];
    }
    // std::cout << "remember to get the extraOutputObject destructed." << std::endl;
    bool oldUpdateModel = recordingInfo.ADinfoPtr()->updateModel();
    recordingInfo.ADinfoPtr()->updateModel() = false;

    atomic_extraOutputObject* localExtraOutputObject =
      new atomic_extraOutputObject("copying-extraOutputObject",
				   &MVA_orig,
				   recordingInfo.ADinfoPtr()); // These objects are stable members of nimbleFunction classes, so the pointer should be stable.

    // Operate the object so it is recorded in the tape
    (*localExtraOutputObject)(extraOutputs, extraOutputDummyResult);
    recordingInfo.ADinfoPtr()->updateModel() = oldUpdateModel;
    // It is unclear whether we then need to use the output in a way that forces CppAD to keep it as part of the calculation graph.
    //   The concern is that otherwise CppAD might optimize it away by determining that nothing really depends on it.
    //   The following line is essentially a no-operation for this purpose (extraOutputDummyResult[0] will always be 0 in value).
    // We do this for all v[i] defensively, to reduce the change that the atomic step is optimized out of the model.
    //    for(size_t i = 0; i < v.size(); ++i)
    //  v[i] += extraOutputDummyResult[0];
    recordingInfo.ADinfoPtr()->add_dummyOutput(extraOutputDummyResult[0]);
  }
  setValues_AD_AD(v, MVA_AD); // record copying on the tape
}

bool atomic_extraOutputObject::forward(
				      const CppAD::vector<double>&               parameter_x  ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				      size_t                              need_y       ,
				      size_t                              order_low    ,
				      size_t                              order_up     ,
				      const CppAD::vector<double>&               taylor_x     ,
				      CppAD::vector<double>&                     taylor_y     ) {
  // std::cout<<"Entering atomic_extraOutputObject named "<<objName<<" with p = "<<p<<" and q = "<<q<<std::endl;
  // std::cout<<"handle address: "<<CppAD::AD<double>::get_handle_address_nimble()<<std::endl;
  // std::cout<<"tx.size() = "<<tx.size()<<" and ty.size() = "<<ty.size()<<std::endl;
  // return flag
  bool ok = true;
    
  size_t length_modelOutput_accessor = MVMA_->getTotalLength(); //NV_->model_modelOutput_accessor.getTotalLength();
  // std::cout << "length_modelOutput_accessor =" << length_modelOutput_accessor << std::endl;
  // if(length_modelOutput_accessor < tx.size()) {
  //   std::cout<<"Problem: length_modelOutput_accessor < vx.size()"<<std::endl;
  // }
  // 0th-order
  size_t q = order_up;
  if( (order_low <= 0) & (ADinfoPtr_->updateModel()) ) {
    // Put values in model.
    //   std::cout<<"putting values in model\n";
    NimArr<1, double> NimArr_tx;
    NimArr_tx.setSize(length_modelOutput_accessor);
    for(size_t i = 0; i < length_modelOutput_accessor; ++i) {
      NimArr_tx[i] = taylor_x[i*(q+1) + 0];
    }
    setValues(NimArr_tx, *MVMA_); //NV_->model_modelOutput_accessor);
    // std::cout<<"done putting values in model\n";
  }

  for(unsigned int i = 0; i < taylor_y.size(); ++i) {
      taylor_y[i] = 0.;
  }

  return ok;
}

bool atomic_extraOutputObject::forward(
				       const CppAD::vector< CppAD::AD<double> >&               parameter_x  ,
				       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				       size_t                              need_y       ,
				       size_t                              order_low    ,
				       size_t                              order_up     ,
				       const CppAD::vector< CppAD::AD<double> >&               taylor_x     ,
				       CppAD::vector< CppAD::AD<double> >&                     taylor_y     ) {
  // This overload should only be used in a case of double-taping.
  // What we do here is what gets recorded on the outer tape when this, the inner tape, is played.
  // We could put a copy of the same, but the decision right now is to put nothing and simply record zeros.
  for(unsigned int i = 0; i < taylor_y.size(); ++i) {
    taylor_y[i] = CppAD::AD<double>(0.);
  }
  return true;
}

CppAD::AD<double> calculate_ADproxyModel(NodeVectorClassNew_derivs &nodes,
					 bool includeExtraOutputStep,
					 nimbleCppADrecordingInfoClass &recordingInfo) {
  //   std::cout <<"entering calculate_ADproxyModel"<< std::endl;
  //  std::cout<<"handle address: "<<CppAD::AD<double>::get_handle_address_nimble()<<std::endl;
  
  CppAD::AD<double> ans = 0;
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  // std::cout<<"starting node calcs in calculate_ADproxyModel"<<std::endl;

  // Set the nodeFunPtr in ADinfo to a single nodeFun (it could be any).
  // This is used during getDerivs_meta to do the tricks to get CppAD
  // to work across our two compilation units.  This is necessary for our atomic classes.
  // Currently, it is thought that this is only necessary if there are
  // nodeFunctions called (i.e. it is not the trivial case of no nodes).
  // It is thought that this is not necessary if the only use of the model
  // compilation unit is copying into and out of the ADproxyModel variables.
  //
  // Note that, at the time this is called, the compilation unit tricks
  // have already been invoked for the current taping from generated code.
  // Stashing it here is for future double taping of this tape.
  if(instructions.size() > 0) {
    nodeFun *nodeFunInModelDLL;
    nodeFunInModelDLL = instructions[0].nodeFunPtr;
    recordingInfo.ADinfoPtr()->set_nodeFunPtr(nodeFunInModelDLL);
  }
  
  for(; iNode != iNodeEnd; iNode++)
    ans += iNode->nodeFunPtr->calculateBlock_ADproxyModel(iNode->operand);
  if(includeExtraOutputStep && recordingInfo.recording()) {
    //   std::cout<<"starting extraOutputStep"<<std::endl;
    if(instructions.begin() != iNodeEnd) {
      //    std::cout<<"will do extraOutputStep"<<std::endl;
      // It is arbitrary to call this for the first node,
      // but it is important to have it done by a nodeFun
      // because that will be in the right compilation unit
      // to access the right globals (statics) from CppAD.
      bool oldUpdateModel = recordingInfo.ADinfoPtr()->updateModel();
      recordingInfo.ADinfoPtr()->updateModel() = false;
      
      instructions.begin()->nodeFunPtr->setup_extraOutput_step( nodes,
      								ans,
								recordingInfo.ADinfoPtr() );
      recordingInfo.ADinfoPtr()->updateModel() = oldUpdateModel;
    }
    //    std::cout<<"done with extraOutputStep"<<std::endl;
  }
  return(ans);
}

void nimble_CppAD_tape_mgr::add_atomic_ptr(nimble_atomic_base *new_atomic_ptr, std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  // std::cout<<"adding atomic_ptr "<<new_atomic_ptr<<" to "<<this<<" with atomic vec_ptr "<<vec_ptr<<std::endl;
  atomic_ptrs.push_back(atomic_pair(new_atomic_ptr, vec_ptr) );
};

void nimble_CppAD_tape_mgr::reset() {
  //  std::cout<<"Doing nimble_CppAD_tape_mgr::reset "<<this<<" "<<atomic_ptrs.size()<<std::endl;
  if(ADtape_) {
    delete ADtape_;
    ADtape_ = nullptr;
  }
  // std::cout<<"Doing nimble_CppAD_tape_mgr::reset 2"<<std::endl;
  //  std::vector<CppAD::local::atomic_index_info>* temp_vec_ptr;
  // temp_vec_ptr = CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage();
  for(size_t i = 0; i < atomic_ptrs.size(); ++i) {
    //    std::cout<<i<<std::endl;
    CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage(1, atomic_ptrs[i].second);
    /* This sets the atomic vec mgr back to what it was when the atomic was constructed.
       This should be redundant and unnecessary because nimble_atomic_base handles the same 
       issue.  So the delete line will call ~nimble_atomic_base before ~atomic_three<>,
       and the formder should also reset the info vec manager. We leave this here
       as a failsafe and because the behavior this solves has been difficult to follow. */
    atomic_ptrs[i].first->set_CppAD_atomic_info_vec_manager(atomic_ptrs[i].second);
    delete atomic_ptrs[i].first;
  }
  //  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage(1, temp_vec_ptr);
  atomic_ptrs.resize(0);
  gammafn_exists = false;
  gammafn_index = 0;
  pow_int_exists = false;
  pow_int_index = 0;
  zround_exists = false;
  zround_index = 0;
  floor_exists = false;
  floor_index = 0;
  ceil_exists = false;
  ceil_index = 0;
  ftrunc_exists = false;
  ftrunc_index = 0;
  nimRound_exists = false;
  nimRound_index = 0;
  log_pow_int_exists = false;
  log_pow_int_index = 0;
  zb_over_a_exists = false;
  zb_over_a_index = 0;
  probit_exists = false;
  probit_index = 0;
  iprobit_exists = false;
  iprobit_index = 0;
  for(int i = 0; i < 5; i++) {
    lgamma_exists[i] = false;
    lgamma_index[i] = 0;
  }
  dummyOutputs.clear();
  //  std::cout<<"Done with nimble_CppAD_tape_mgr::reset"<<std::endl;
}

void nimble_CppAD_tape_mgr::add_dummyOutput(CppAD::AD<double> &dummy) {
  dummyOutputs.push_back(dummy);
}

void add_dummyOutput(nimbleCppADinfoClass *ADinfoPtr, CppAD::AD<double> &dummy) {
  ADinfoPtr->add_dummyOutput(dummy);
}

void nimble_CppAD_tape_mgr::sum_dummyOutputs_to_dependentVars(std::vector<CppAD::AD<double> > &depVars) {
  if(dummyOutputs.size() == 0) return;
  CppAD::AD<double> firstSum = dummyOutputs[0];
  if(dummyOutputs.size() > 1) {
    for(size_t i = 1; i < dummyOutputs.size(); ++i) {
      firstSum += dummyOutputs[i];
    }
  }
  for(size_t j = 0; j < depVars.size(); ++j) {
    depVars[j] += firstSum;
  }
}

void nimble_CppAD_tape_mgr::set_internal_tape(CppAD::local::ADTape<double>* internal_tape_ptr) {
  internal_tape_ptr_ = internal_tape_ptr;
  internal_tape_ptr->nimble_CppAD_tape_mgr_ptr() = static_cast<void*>(this);
};
nimble_CppAD_tape_mgr::nimble_CppAD_tape_mgr() :
  ADtape_(0),
  internal_tape_ptr_(0),
  gammafn_index(0),
  gammafn_exists(false),
  pow_int_index(0),
  pow_int_exists(false),
  zround_index(0),
  zround_exists(false),
  floor_index(0),
  floor_exists(false),
  ceil_index(0),
  ceil_exists(false),
  ftrunc_index(0),
  ftrunc_exists(false),
  nimRound_index(0),
  nimRound_exists(false),
  log_pow_int_index(0),
  log_pow_int_exists(false),
  zb_over_a_index(0),
  zb_over_a_exists(false),
  probit_index(0),
  probit_exists(false),
  iprobit_index(0),
  iprobit_exists(false)
{
  for(int i = 0; i < 5; i++) {
    lgamma_exists[i] = false;
    lgamma_index[i] = 0;
  }
};
nimble_CppAD_tape_mgr::~nimble_CppAD_tape_mgr() {
  // std::cout<<"Doing ~nimble_CppAD_tape_mgr "<<this<<std::endl;
  reset();
};


nimble_atomic_base::nimble_atomic_base() {
  vec_ptr_where_constructed = CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage();
}

nimble_atomic_base::~nimble_atomic_base() {
  //  std::cout<<"trying to delete nimble_atomic_base"<<std::endl;
  set_CppAD_atomic_info_vec_manager(vec_ptr_where_constructed);
}

void nimble_atomic_base::set_CppAD_atomic_info_vec_manager( std::vector<CppAD::local::atomic_index_info>* vec_ptr ) {
  // std::cout<<"Setting atomic vec_ptr from value of "<<CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage()<<" to "<<vec_ptr<<std::endl;
  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage(1, vec_ptr);
}

void track_nimble_atomic(nimble_atomic_base *obj, void *tape_mgr_ptr, std::vector<CppAD::local::atomic_index_info>* vec_ptr ) {
  // std::cout<<"track_nimble_atomic for "<<obj<<" in "<<tape_mgr_ptr<<" with atomic vec_ptr "<< vec_ptr<<std::endl;
  reinterpret_cast<nimble_CppAD_tape_mgr*>(tape_mgr_ptr)->add_atomic_ptr(obj, vec_ptr);
}
