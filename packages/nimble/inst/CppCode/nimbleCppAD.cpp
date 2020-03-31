#include <nimble/nimbleCppAD.h>

void copy_CppADdouble_to_double(CppAD::AD<double> *first, CppAD::AD<double> *last, double *output) {
  CppAD::AD<double> *orig;
  double *result = output;
  for(orig = first; orig != last ; )
    *result++ = CppAD::Value(*orig++);
}

#ifdef _TIME_AD
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

  if(!nodes.tapeRecorded()) nodes.recordTape();
#ifdef _TIME_AD
  derivs_main_timer.start();
#endif

  nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
#ifdef _TIME_AD
  derivs_getDerivs_timer.start();
#endif
  std::vector<double> independentVars;
  nodes.runTape_setIndependent(independentVars);
  std::vector<double> dependentVars;

#ifdef _TIME_AD
  derivs_run_tape_timer.start();
#endif
  
  nodes.runTape_runTape(independentVars, dependentVars,
			derivOrders, ansList);
#ifdef _TIME_AD
  derivs_run_tape_timer.stop();
#endif

#ifdef _TIME_AD
  derivs_getDerivs_timer.stop();
#endif

#ifdef _TIME_AD
  derivs_main_timer.stop();
#endif

  return ansList;
}

#endif

nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
    NodeVectorClassNew_derivs &nodes, const double derivOrders) {
  NimArr<1, double> orders(1);
  orders[0] = derivOrders;
  return (nimDerivs_calculate(nodes, orders));
}



nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
    NodeVectorClassNew_derivs &nodes, int iNodeFunction,
    NimArr<1, double> &derivOrders) {
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
//std::vector<double> &dynamicVars,
//			CppAD::ADFun<double>* &tapePtr) {
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  NimArr<1, double> NimArrValues;
  if(length_extraInput > 0) {
    NimArrValues.setSize(length_extraInput);
    ADinfo.dynamicVars.resize(length_extraInput);
    getValues(NimArrValues, NV.model_extraInput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_extraInput,
	       ADinfo.dynamicVars.begin());
    ADinfo.ADtape->new_dynamic(ADinfo.dynamicVars);
  }
}

void update_dynamicVars_meta(NodeVectorClassNew_derivs &NV,
			     nimbleCppADinfoClass &ADinfo) {
  //std::vector< CppAD::AD<double> > &dynamicVars,
  //			     CppAD::ADFun<double>* &tapePtr) {
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
  std::cout<<"Use of dynamicVars in meta-taping is ambiguous."<<std::endl;
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
  
#ifdef _TIME_AD
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
  #ifdef _TIME_AD
  derivs_run_tape_timer_start();
#endif
  value_ans = ADtape->Forward(0, independentVars);
  //  std::cout<<"value_ans.size() = "<<value_ans.size()<<std::endl;
#ifdef _TIME_AD
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
      std::cout<<"Fix the inf and nan checking for CppAD::AD<double> case"<<std::endl;
      // if(((value_ans[inf_ind] == -std::numeric_limits<double>::infinity()) |
      //     (value_ans[inf_ind] == std::numeric_limits<double>::infinity())) | 
      // 	 (std::isnan(value_ans[inf_ind]))){
      // 	infIndicators[inf_ind] = true;
      // }
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
#ifdef _TIME_AD
	  derivs_run_tape_timer_start();
#endif
	  cppad_derivOut = ADtape->Reverse(1, w);
#ifdef _TIME_AD
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
#ifdef _TIME_AD
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

#ifdef _TIME_AD
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
#ifdef _TIME_AD
  derivs_getDerivs_timer_stop();
#endif
};

void nimbleFunctionCppADbase::getDerivs_meta(nimbleCppADinfoClass &ADinfo,
					     const NimArr<1, double> &derivOrders,
					     const NimArr<1, double> &wrtVector,
					     nimSmartPtr<NIMBLE_ADCLASS_META> &ansList) {
  // std::cout<<"Entering getDerivs_meta"<<std::endl;
  CppAD::ADFun< CppAD::AD<double>, double > innerTape;
  innerTape = ADinfo.ADtape->base2ad();
  getDerivs_internal< CppAD::AD<double>,
		      CppAD::ADFun< CppAD::AD<double>, double >,
		      NIMBLE_ADCLASS_META>(ADinfo.independentVars_meta,
					   &innerTape,
					   derivOrders,
					   wrtVector,
					   ansList);
  // std::cout<<"Exiting getDerivs_meta"<<std::endl;
}
  
void nimbleFunctionCppADbase::getDerivs(nimbleCppADinfoClass &ADinfo,
                                        const NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
  // std::cout<<"Entering getDerivs"<<std::endl;
  getDerivs_internal<double,
		     CppAD::ADFun<double>,
		     NIMBLE_ADCLASS>(ADinfo.independentVars,
		       ADinfo.ADtape,
		       derivOrders,
		       wrtVector,
		       ansList);
  // std::cout<<"Exiting getDerivs"<<std::endl;
    
// #ifdef _TIME_AD
//   derivs_getDerivs_timer_start();
//   derivs_tick_id();
//   derivs_show_id();  
// #endif
//   std::size_t n = ADinfo.independentVars.size();  // dim of independent vars

//   std::size_t wrt_n = wrtVector.size();            // dim of wrt vars
//   if(wrt_n == 2){
//     if(wrtVector[1] == -1){
//       wrt_n = 1;
//     }
//   }
//   int orderSize = derivOrders.size();
//   double const* array_derivOrders = derivOrders.getConstPtr();

//   int maxOrder =
//     *std::max_element(array_derivOrders, array_derivOrders + orderSize);
//   bool ordersFound[3] = {false};

//   for (int i = 0; i < orderSize; i++) {
//     if ((array_derivOrders[i] > 2) | (array_derivOrders[i] < 0)) {
//       printf("Error: Derivative orders must be between 0 and 2.\n");
//     }
//     ordersFound[static_cast<int>(array_derivOrders[i])] = true;
//   }
//   vector<double> value_ans;
// #ifdef _TIME_AD
//   derivs_run_tape_timer_start();
// #endif
//   value_ans = ADinfo.ADtape->Forward(0, ADinfo.independentVars);
// #ifdef _TIME_AD
//   derivs_run_tape_timer_stop();
// #endif
//   if (ordersFound[0]) {
//     ansList->value = vectorDouble_2_NimArr(value_ans);
//   }
//   if(maxOrder > 0){
//     std::size_t q = value_ans.size();
//     vector<bool> infIndicators(q); // default values will be false 
//     for(size_t inf_ind = 0; inf_ind < q; inf_ind++){
//       if(((value_ans[inf_ind] == -std::numeric_limits<double>::infinity()) |
//           (value_ans[inf_ind] == std::numeric_limits<double>::infinity())) | 
// 	 (std::isnan(value_ans[inf_ind]))){
// 	infIndicators[inf_ind] = true;
//       }
//     }
//     if (ordersFound[1]) {
//       ansList->jacobian.setSize(q, wrt_n, false, false); // setSize may be costly.  Possible to setSize outside of fxn, within chain rule algo, and only resize when necessary?
//     }
//     if (ordersFound[2]) {
//       ansList->hessian.setSize(wrt_n, wrt_n, q, false, false);
//     }
//     vector<double> cppad_derivOut;
//     std::vector<double> w(q, 0);
//     for (size_t dy_ind = 0; dy_ind < q; dy_ind++) {
//       //      std::vector<double> w(q, 0);
//       w[dy_ind] = 1;
//       if (maxOrder == 1) {   
// 	if(!infIndicators[dy_ind]){
// #ifdef _TIME_AD
// 	  derivs_run_tape_timer_start();
// #endif
// 	  cppad_derivOut = ADinfo.ADtape->Reverse(1, w);
// #ifdef _TIME_AD
// 	  derivs_run_tape_timer_stop();
// #endif
// 	}
//       } else {
// 	for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
// 	  if(!infIndicators[dy_ind]){
// 	    int dx1_ind = wrtVector[vec_ind] - 1;
// 	    std::vector<double> x1(n, 0);  // vector specifying first derivatives.
// 	    // first specify coeffs for first dim
// 	    // of s across all directions r, then
// 	    // second dim, ...
// 	    x1[dx1_ind] = 1;
// #ifdef _TIME_AD
// 	    derivs_run_tape_timer_start();
// #endif
// 	    ADinfo.ADtape->Forward(1, x1);
// 	    cppad_derivOut = ADinfo.ADtape->Reverse(2, w);
// #ifdef _TIME_AD
// 	    derivs_run_tape_timer_stop();
// #endif
// 	  }
// 	  for (size_t vec_ind2 = 0; vec_ind2 < wrt_n; vec_ind2++) {
// 	    if(!infIndicators[dy_ind]){
// 	      int dx2_ind = wrtVector[vec_ind2] - 1;
// 	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] =
// 		cppad_derivOut[dx2_ind * 2 + 1];
// 	    }
// 	    else{
// 	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] = 
// 		CppAD::numeric_limits<double>::quiet_NaN();
// 	    }
// 	  }
// 	}
//       }
//       if (ordersFound[1]) {
// 	double *LHS = ansList->jacobian.getPtr() + dy_ind;
// 	if(!infIndicators[dy_ind]){
// 	  double const *wrtVector_p = wrtVector.getConstPtr();
// 	  double const *wrtVector_p_end = wrtVector_p + wrt_n;
// 	    for(; wrtVector_p != wrtVector_p_end; LHS += q ) {
// 	      *LHS = cppad_derivOut[(static_cast<int>(*wrtVector_p++) - 1) * maxOrder];
// 	    }
// 	} else {
// 	  for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
// 	    *LHS = CppAD::numeric_limits<double>::quiet_NaN();
// 	    LHS += q;
// 	  }
// 	}

// 	// for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
// 	//   if(!infIndicators[dy_ind]){
// 	//     int dx1_ind = wrtVector[vec_ind] - 1;
// 	//     ansList->jacobian[vec_ind * q + dy_ind] =
// 	//       cppad_derivOut[dx1_ind * maxOrder + 0];
// 	//   }
// 	//   else{
// 	//     ansList->jacobian[vec_ind * q + dy_ind] =
// 	//       CppAD::numeric_limits<double>::quiet_NaN();
// 	//   }     
// 	// }
	
//       }
//       w[dy_ind] = 0;
//     }
//   }
// #ifdef _TIME_AD
//   derivs_getDerivs_timer_stop();
// #endif
}


CppAD::ADFun<double>* calculate_recordTape(NodeVectorClassNew_derivs &NV) {
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
  CppAD::Independent(independentVars, abort_op_index, record_compare, dynamicVars);

  set_CppAD_tape_info_for_model my_tape_info_RAII_(NV,
						   CppAD::AD<double>::get_tape_id_nimble(),
						   CppAD::AD<double>::get_tape_handle_nimble());
  set_CppAD_atomic_info_for_model(NV, CppAD::local::atomic_index_info_vec_manager<double>::manage());

  // 7. [deleted]
  // 8. [deleted]

  // 9. Copy all wrtNodes AD objects from independentVars -> model_AD.
  if(length_wrt > 0) {
    NimArrValues_AD.setSize(length_wrt);
    std::copy(independentVars.begin(),
	      independentVars.begin() + length_wrt,
	      NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
  }
  // 10. call calculate.  This also sets up the extraOutput step
  CppAD::AD<double> logProb = calculate_ADproxyModel(NV,
						     true,
						     true);
  dependentVars[0] = logProb;
  // 13. Finish taping, AND
  // 14. Call tape->optimize()
  // make it locally to get the right globals during recording and playback
  // DO NOT USE THE CONSTRUCTOR VERSION BECAUSE IT ALWAYS DOES .Forward(0)
  // INSTEAD MAKE THE BLANK OBJECT AND USE .Dependent(...)
  // TRY USING CppAD's vector type
  CppAD::ADFun<double>* ansTape = new CppAD::ADFun<double>;
  ansTape->Dependent(independentVars, dependentVars);
  ansTape->optimize(); //("no_compare_op") makes almost no difference;
  return ansTape;
}

void nimbleFunctionCppADbase::getDerivs_calculate_internal(nimbleCppADinfoClass &ADinfo,
							   CppAD::ADFun<double>* &tapePtr,
							   NodeVectorClassNew_derivs &nodes,
							   const NimArr<1, double> &derivOrders,
							   const NimArr<1, double> &wrtVector,
							   nimSmartPtr<NIMBLE_ADCLASS> ansList) {
  bool use_meta_tape = true;

  if(!tapePtr) {
    if(!use_meta_tape) {
      tapePtr = calculate_recordTape(nodes);
    } else {
      CppAD::ADFun< double > *firstTape;
      firstTape = calculate_recordTape(nodes);
      CppAD::ADFun< CppAD::AD<double>, double > innerTape;
      // Make original tape use CppAD::AD<double> instead of double
      set_CppAD_atomic_info_for_model(nodes, CppAD::local::atomic_index_info_vec_manager<double>::manage());
      innerTape = firstTape->base2ad();
      int length_wrt = nodes.model_wrt_accessor.getTotalLength();
      int length_independent = length_wrt;
      int length_extraInput = nodes.model_extraInput_accessor.getTotalLength();
      vector< CppAD::AD<double> > dependentVars(length_wrt); // This will be the gradient
      vector< CppAD::AD<double> > independentVars(length_independent);
      vector< CppAD::AD<double> > dynamicVars(length_extraInput);
      size_t abort_op_index = 0; 
      bool   record_compare = true;
      NimArr<1, double > NimArrVars;
      NimArrVars.setSize(length_wrt);
      getValues(NimArrVars, nodes.model_wrt_accessor);
      std::copy(NimArrVars.getPtr(),
		NimArrVars.getPtr() + length_wrt,
		independentVars.begin());
      if(length_extraInput > 0) {
	NimArr<1, double> NimArr_dynamicVars;
	NimArr_dynamicVars.setSize(length_extraInput);
	getValues(NimArr_dynamicVars, nodes.model_extraInput_accessor);
	std::copy( NimArr_dynamicVars.getPtr(),
		   NimArr_dynamicVars.getPtr() + length_extraInput,
		   dynamicVars.begin() );
      }
      nimSmartPtr<NIMBLE_ADCLASS_META> ansList_meta = new NIMBLE_ADCLASS_META;
      CppAD::Independent(independentVars, abort_op_index, record_compare, dynamicVars); // start recording new tape
      set_CppAD_tape_info_for_model my_tape_info_RAII_(nodes,
						       CppAD::AD<double>::get_tape_id_nimble(),
						       CppAD::AD<double>::get_tape_handle_nimble());
      set_CppAD_atomic_info_for_model(nodes, CppAD::local::atomic_index_info_vec_manager<double>::manage());
      ADinfo.independentVars_meta.resize(length_wrt);
      std::copy(independentVars.begin(),
		independentVars.end(),
		ADinfo.independentVars_meta.begin());
      NimArr<1, double> derivOrders_meta;
      derivOrders_meta.setSize(1);
      derivOrders_meta[0] = 1;
      getDerivs_internal< CppAD::AD<double>,
			  CppAD::ADFun< CppAD::AD<double>, double >,
			  NIMBLE_ADCLASS_META>(ADinfo.independentVars_meta,
					       &innerTape,
					       derivOrders_meta,
					       wrtVector,
					       ansList_meta);
      for(size_t iii = 0; iii < length_wrt; ++iii)
	dependentVars[iii] = ansList_meta->jacobian[iii];
      tapePtr = new CppAD::ADFun<double>;
      tapePtr->Dependent(dependentVars);
      tapePtr->optimize();
      delete firstTape;
    }
  }

  // std::cout<<"getDerivs_calculate_internal A"<<std::endl;
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
  size_t length_extraNodes_accessor = nodes.model_extraInput_accessor.getTotalLength();
  if(length_extraNodes_accessor > 0) {
    NimArr<1, double> NimArr_dynamicVars;
    NimArr_dynamicVars.setSize(length_extraNodes_accessor);
    getValues(NimArr_dynamicVars, nodes.model_extraInput_accessor);
    std::vector<double> dynamicVars(length_extraNodes_accessor);
    std::copy( NimArr_dynamicVars.getPtr(),
	       NimArr_dynamicVars.getPtr() + length_extraNodes_accessor,
	       dynamicVars.begin() );
    tapePtr->new_dynamic(dynamicVars);
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
      getDerivs_internal<double,
			 CppAD::ADFun<double>,
			 NIMBLE_ADCLASS>(ADinfo.independentVars,
					 tapePtr,
					 derivOrders_nested,
					 wrtVector, // NOTE: This will not behave fully correctly in non-default use without further thought.
					 ansList_nested);
      if(ordersFound[1]) {
	ansList->jacobian.setSize(1, length_wrt, false, false);
	for(size_t ii = 0; ii < length_wrt; ++ii) //We could rewrite this with better pointer magic
	  ansList->jacobian[ ii ] = ansList_nested->value[ ii ];
      }
      if(ordersFound[2]) {
	ansList->hessian.setSize(length_wrt, length_wrt, 1, false, false);
	for(size_t ii = 0; ii < length_wrt; ++ii)
	  for(size_t jj = 0; jj < length_wrt; ++jj)
	    ansList->hessian[jj + ii*length_wrt ] = ansList_nested->jacobian[jj + ii*length_wrt]; //orientation shouldn't matter due to symmetry
      }
    }
  } else {
  /* run tape */
    getDerivs_internal<double,
		       CppAD::ADFun<double>,
		       NIMBLE_ADCLASS>(ADinfo.independentVars,
				       tapePtr,
				       derivOrders,
				       wrtVector,
				       ansList);
  }
}

NimArr<1, double> make_vector_if_necessary(int a){
      NimArr<1, double> intArray;
      intArray.setSize(1);
      intArray[0] = a;
      return(intArray);
}

NimArr<1, double> make_vector_if_necessary(double a){
      NimArr<1, double> intArray;
      intArray.setSize(1);
      intArray[0] = a;
      return(intArray);
}

NimArr<1, double> make_vector_if_necessary(NimArr<1, double> a){
      return(a);
}
