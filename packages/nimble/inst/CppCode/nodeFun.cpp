#include "nimble/nodeFun.h"
#include "nimble/accessorClasses.h"


// This should be called before CppAD::Independent
void nodeFun::initialize_AD_model_before_recording(NodeVectorClassNew_derivs &NV) {
  
  NimArr<1, double> NimArrValues;
  NimArr<1, CppAD::AD<double> > NimArrValues_AD;
  
  // 1. Copy all constantNodes values from model -> model_AD
  int length_constant = NV.model_constant_accessor.getTotalLength();
  if(length_constant > 0) {
    NimArrValues.setSize(length_constant);
    NimArrValues_AD.setSize(length_constant);
    getValues(NimArrValues, NV.model_constant_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_constant,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_constant_accessor);
  }
  
  // Do the same for modelOutput nodes -- somewhat defensive -- code might use getLogProb before calculate but in that case it might not record correctly anyway
  int length_modelOutput = NV.model_modelOutput_accessor.getTotalLength();
  if(length_modelOutput > 0) {
    NimArrValues.setSize(length_modelOutput);
    NimArrValues_AD.setSize(length_modelOutput);
    getValues(NimArrValues, NV.model_modelOutput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_modelOutput,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_modelOutput_accessor);
  }

  // And the same for with-respect-to nodes
  int length_wrt = NV.model_wrt_accessor.getTotalLength();
  if(length_wrt > 0) {
    NimArrValues.setSize(length_wrt);
    NimArrValues_AD.setSize(length_wrt);
    getValues(NimArrValues, NV.model_wrt_accessor);

    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_wrt,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
  }
}  

void nodeFun::set_atomic_info_from_nodeFun(std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage(1, vec_ptr);
}

void nodeFun::set_tape_ptr_from_nodeFun(CppAD::tape_id_t tape_id,
					CppAD::local::ADTape<double>* tape_handle_,
					bool recover) {
  CppAD::AD<double>::set_tape_info_nimble(tape_id, tape_handle_, recover);
}

void nodeFun::setup_extraOutput_step(NodeVectorClassNew_derivs &NV,
				     CppAD::AD<double> &logProb,
				     nimbleCppADinfoClass* ADinfoPtr) {
  // Call right after logProb = calculate_ADproxyModel.
  NV.extraOutputObject = 
    runExtraOutputObject(NV,
			 // extraOutputs,
			 //extraOutputDummyResult,
			 logProb,
			 ADinfoPtr);
}

void nodeFun::setTapeIndependent(std::vector< CppAD::AD<double> > &independentVars) {
  CppAD::Independent(independentVars);
}
void nodeFun::finishADFun(CppAD::ADFun< double > &ADtape,
			  std::vector< CppAD::AD<double> > &independentVars,
			  std::vector< CppAD::AD<double> > &dependentVars) {
  // make it locally to get the right globals during recording and playback
  // DO NOT USE THE CONSTRUCTOR VERSION BECAUSE IT ALWAYS DOES .Forward(0)
  // INSTEAD MAKE THE BLANK OBJECT AND USE .Dependent(...)
  // TRY USING CppAD's vector type
  CppAD::ADFun<double> localADtape(independentVars, dependentVars);
  //  CppAD::ADFun<double> localADtape;
  //  localADtape.Dependent(independentVars, dependentVars);
  localADtape.optimize(); //("no_compare_op") makes almost no difference;
  // then copy it back to object from another DLL
  ADtape = localADtape;
  //  ADtape.capacity_order(2); // makes only a tiny improvement
}

void nodeFun::runTape(CppAD::ADFun< double > &ADtape,
		      std::vector< double > &independentVars,
		      std::vector< double > &dependentVars,
		      const NimArr<1, double> &derivOrders,
		      nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
  // dependentVars = ADtape.Jacobian(independentVars);
  std::size_t n = independentVars.size();
  // Always calculate value, regardless of valueFlag.
  // Is this the right behavior?
  dependentVars = ADtape.Forward(0, independentVars);
  ansList->value.setSize(1, false, false);
  ansList->value[0] = dependentVars[0];
  
  bool hessianFlag = false;   // Are second order derivs requested?
  bool jacobianFlag = false;  // Are first order derivs requested?

  for (int i = 0; i < derivOrders.dimSize(0); i++) {
    if (derivOrders[i] == 1) {
      jacobianFlag = true;
    }
    if (derivOrders[i] == 2) {
      hessianFlag = true;
      jacobianFlag = true;  // If second order is asked for, first order must be
                            // calculated as well (although we could technically
                            // skip copying this into output)
    }
  }
  if(jacobianFlag) {
    std::vector<double> rev;
    std::vector<double> w(1);
    w[0] = 1;
    rev = ADtape.Reverse(1, w);
    size_t jacSize = rev.size()-1;
    ansList->jacobian.setSize( 1, jacSize, false, false);
    std::copy(rev.begin(), rev.begin() + jacSize, ansList->jacobian.getPtr());
    if(hessianFlag) {
      ansList->hessian.setSize( jacSize, jacSize, 1, false, false);
      for (size_t dx_ind = 0; dx_ind < jacSize; dx_ind++) {
        std::vector<double> x1(n, 0);
        x1[dx_ind] = 1;
        ADtape.Forward(1, x1);
        rev = ADtape.Reverse(2, w);
	for (size_t dx_ind2 = 0; dx_ind2 < jacSize; dx_ind2++) {
	  ansList->hessian[jacSize * dx_ind + dx_ind2] =
	    rev[dx_ind2 * 2 + 1];
        }
      }
    }
  }
}


atomic_extraOutputObject*
nodeFun::runExtraOutputObject(NodeVectorClassNew_derivs &NV,
			      CppAD::AD<double> &logProb,
			      nimbleCppADinfoClass* ADinfoPtr) {
  size_t length_modelOutput = NV.model_modelOutput_accessor.getTotalLength();
  // std::cout<<"runExtraOutputObject length_modelOutput = "<<length_modelOutput<<std::endl;
  NimArr<1, CppAD::AD<double> > NAV_AD;
  NAV_AD.setSize(length_modelOutput);
  std::vector< CppAD::AD<double> > extraOutputs(length_modelOutput);
  std::vector< CppAD::AD<double> > extraOutputDummyResult(1);
  extraOutputDummyResult[0] = 0;
  // std::cout<<"runExtraOutputObject A"<<std::endl;
  getValues_AD_AD( NAV_AD, NV.model_AD_modelOutput_accessor);
  
  for(size_t ii = 0; ii < length_modelOutput; ++ii) {
    extraOutputs[ii] = NAV_AD[ii];
  }
  // std::cout<<"runExtraOutputObject B"<<std::endl;
  
  // Calling the constructor with a local object will use statics in the correct DLL
  // std::cout<<"runExtraOutputObject C"<<std::endl;
  atomic_extraOutputObject* localExtraOutputObject =
    new atomic_extraOutputObject("extraOutputObject",
				 &(NV.model_modelOutput_accessor),
				 ADinfoPtr);
				 //&NV);
  // And then calling operator() with the local object will again use statics in the correct DLL
  // std::cout<<"runExtraOutputObject D"<<std::endl;
  (*localExtraOutputObject)(extraOutputs, extraOutputDummyResult);
  // Should be ok to copy.  No further uses after taping involve statics,
  // so there should be no further issue with which DLL we are in.
  // std::cout<<"runExtraOutputObject E"<<std::endl;
  logProb += extraOutputDummyResult[0];
  // std::cout<<"runExtraOutputObject E"<<std::endl;
  return localExtraOutputObject;
}

void
nodeFun::delete_extraOutputObject(NodeVectorClassNew_derivs &NV) {
  // std::cout<<"About to delete extraOutputObject"<<std::endl;
  delete NV.extraOutputObject;
  // std::cout<<"Done deleting extraOutputObject"<<std::endl;
}

// CppAD::AD<double> nodeFun::call_calculate_ADproxyModel(NodeVectorClassNew_derivs &NV) {
//   CppAD::AD<double> logProb = ::calculate_ADproxyModel( NV );
//   return logProb;
// }


// This makes a proxy "y = f(x)"
// where y are a set of variables we want to get into the
// calculations without introducing them via CppAD::Independent.
// This doesn't make much difference for 0-order forward performance
// but should simplify use of higher-order forward and reverse
// calls.
// By definition any derivatives are zero for this proxy f() 

////////
// This makes a proxy "y = f(x)"
// where y is a dead-end variable (nothing else should depend on it,
// to avoid unnecessary tape operations) and x is a vector of
// variables whose values we want to copy back into the 
// to the regular model.  This will typically include log probability
// values and deterministic results.
//
// We use it at the last step, like
// dummy = f(logProbs_and_determ).
//
// By definition the value and any derivatives for this proxy f()
// are zero.
atomic_extraOutputObject::atomic_extraOutputObject(const std::string& name,
						   ManyVariablesMapAccessor* MVMA,
						   nimbleCppADinfoClass* ADinfoPtr) :
  CppAD::atomic_three<double>(name),
  ADinfoPtr_(ADinfoPtr),
  MVMA_(MVMA),
  objName(name)
{
  // std::cout<< "Constructing atomic_extraOutputObject named "<< name <<std::endl;
  // std::cout<<"handle address: "<<CppAD::AD<double>::get_handle_address_nimble()<<std::endl;
}

// reverse
bool
atomic_extraOutputObject::reverse(
				  const CppAD::vector<double>&               parameter_x ,
				  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				  size_t                              order_up    ,
				  const CppAD::vector<double>&               taylor_x    ,
				  const CppAD::vector<double>&               taylor_y    ,
				  CppAD::vector<double>&                     partial_x   ,
				  const CppAD::vector<double>&               partial_y   ) {
  size_t n = partial_x.size();
  //  std::cout<<"in atomic_extraOutputObject::reverse"<<std::endl;
  for(size_t i = 0; i < n; ++i)
    partial_x[i] = 0;
  return true;
}

bool
atomic_extraOutputObject::reverse(
				  const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				  size_t                              order_up    ,
				  const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
				  const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
				  CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
				  const CppAD::vector<CppAD::AD<double> >&               partial_y   ) {
  size_t n = partial_x.size();
  //  std::cout<<"in atomic_extraOutputObject::reverse"<<std::endl;
  for(size_t i = 0; i < n; ++i)
    partial_x[i] = CppAD::AD<double>(0);
  return true;
}
