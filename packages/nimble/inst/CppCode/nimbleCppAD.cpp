#include <nimble/nimbleCppAD.h>

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

void nimbleFunctionCppADbase::getDerivs(nimbleCppADinfoClass &ADinfo,
                                        const NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
#ifdef _TIME_AD
  derivs_getDerivs_timer_start();
  derivs_tick_id();
  derivs_show_id();  
#endif
  std::size_t n = ADinfo.independentVars.size();  // dim of independent vars

  std::size_t wrt_n = wrtVector.size();            // dim of wrt vars
  if(wrt_n == 2){
    if(wrtVector[1] == -1){
      wrt_n = 1;
    }
  }
  int orderSize = derivOrders.size();
  double const* array_derivOrders = derivOrders.getConstPtr();

  int maxOrder =
    *std::max_element(array_derivOrders, array_derivOrders + orderSize);
  bool ordersFound[3] = {false};

  for (int i = 0; i < orderSize; i++) {
    if ((array_derivOrders[i] > 2) | (array_derivOrders[i] < 0)) {
      printf("Error: Derivative orders must be between 0 and 2.\n");
    }
    ordersFound[static_cast<int>(array_derivOrders[i])] = true;
  }
  vector<double> value_ans;
#ifdef _TIME_AD
  derivs_run_tape_timer_start();
#endif
  value_ans = ADinfo.ADtape->Forward(0, ADinfo.independentVars);
#ifdef _TIME_AD
  derivs_run_tape_timer_stop();
#endif
  if (ordersFound[0]) {
    ansList->value = vectorDouble_2_NimArr(value_ans);
  }
  if(maxOrder > 0){
    std::size_t q = value_ans.size();
    vector<bool> infIndicators(q); // default values will be false 
    for(size_t inf_ind = 0; inf_ind < q; inf_ind++){
      if(((value_ans[inf_ind] == -std::numeric_limits<double>::infinity()) |
          (value_ans[inf_ind] == std::numeric_limits<double>::infinity())) | 
	 (std::isnan(value_ans[inf_ind]))){
	infIndicators[inf_ind] = true;
      }
    }
    if (ordersFound[1]) {
      ansList->jacobian.setSize(q, wrt_n, false, false); // setSize may be costly.  Possible to setSize outside of fxn, within chain rule algo, and only resize when necessary?
    }
    if (ordersFound[2]) {
      ansList->hessian.setSize(wrt_n, wrt_n, q, false, false);
    }
    vector<double> cppad_derivOut;
    std::vector<double> w(q, 0);
    for (size_t dy_ind = 0; dy_ind < q; dy_ind++) {
      //      std::vector<double> w(q, 0);
      w[dy_ind] = 1;
      if (maxOrder == 1) {   
	if(!infIndicators[dy_ind]){
#ifdef _TIME_AD
	  derivs_run_tape_timer_start();
#endif
	  cppad_derivOut = ADinfo.ADtape->Reverse(1, w);
#ifdef _TIME_AD
	  derivs_run_tape_timer_stop();
#endif
	}
      } else {
	for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
	  if(!infIndicators[dy_ind]){
	    int dx1_ind = wrtVector[vec_ind] - 1;
	    std::vector<double> x1(n, 0);  // vector specifying first derivatives.
	    // first specify coeffs for first dim
	    // of s across all directions r, then
	    // second dim, ...
	    x1[dx1_ind] = 1;
#ifdef _TIME_AD
	    derivs_run_tape_timer_start();
#endif
	    ADinfo.ADtape->Forward(1, x1);
	    cppad_derivOut = ADinfo.ADtape->Reverse(2, w);
#ifdef _TIME_AD
	    derivs_run_tape_timer_stop();
#endif
	  }
	  for (size_t vec_ind2 = 0; vec_ind2 < wrt_n; vec_ind2++) {
	    if(!infIndicators[dy_ind]){
	      int dx2_ind = wrtVector[vec_ind2] - 1;
	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] =
		cppad_derivOut[dx2_ind * 2 + 1];
	    }
	    else{
	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] = 
		CppAD::numeric_limits<double>::quiet_NaN();
	    }
	  }
	}
      }
      if (ordersFound[1]) {
	double *LHS = ansList->jacobian.getPtr() + dy_ind;
	if(!infIndicators[dy_ind]){
	  double const *wrtVector_p = wrtVector.getConstPtr();
	  double const *wrtVector_p_end = wrtVector_p + wrt_n;
	    for(; wrtVector_p != wrtVector_p_end; LHS += q ) {
	      *LHS = cppad_derivOut[(static_cast<int>(*wrtVector_p++) - 1) * maxOrder];
	    }
	} else {
	  for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
	    *LHS = CppAD::numeric_limits<double>::quiet_NaN();
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
