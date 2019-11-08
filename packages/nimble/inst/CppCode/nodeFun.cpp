#include "nimble/nodeFun.h"
#include "nimble/accessorClasses.h"


// useADreconfigure version ("2")
// Need to call this at the start of a function that will use a model in a tape.
// Need to think about nested cases.
// Need to think about extraInputVar and extraOutputVar
void nodeFun::prepare_to_record_model_tape(NodeVectorClassNew_derivs &NV, CppAD::AD<double> &extraInputDummy) { // move extraInputDummy in NV? Careful with compilation units.
  //  vector< CppAD::AD<double> > dependentVars(1);
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
  
  // Do the same for modelOutput nodes -- somewhat defensive -- code might use getLogProb before calculate but in that case it might not record correctly anyway
  int length_modelOutput = NV.model_modelOutput_accessor.getTotalLength();
  if(length_modelOutput > 0) {
    NimArr<1, double> NimArrValues;
    NimArr<1, CppAD::AD<double> > NimArrValues_AD;
    NimArrValues.setSize(length_modelOutput);
    NimArrValues_AD.setSize(length_modelOutput);
    getValues(NimArrValues, NV.model_modelOutput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_modelOutput,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_modelOutput_accessor);
  }
  
  int length_wrt = NV.model_wrt_accessor.getTotalLength();
  if(length_wrt > 0) {
    NimArrValues.setSize(length_wrt);
    getValues(NimArrValues, NV.model_wrt_accessor);
    // 2
    NimArrValues_AD.setSize(length_wrt);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_wrt,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
  }
  
  // 5. Copy all extraInputNodes values from model -> model_AD (ditto, may be redundant)  
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

  // record extraInputObject
  if(length_extraInput > 0) {
    std::vector< CppAD::AD<double> > extraInputDummyInput(1);
    extraInputDummyInput[0] = extraInputOutputVars[ 0 ];
    std::vector< CppAD::AD<double> > extraInputResults(length_extraInput);
      // After this, the tape treats extraInputResults as a function of
      // extraInputDummy.  This means during tape use, the extraInputObject
      // functor will be called, and the extraInputResults will contain
      // node values from the model.  By copying those during taping
      // into the model_AD, those nodes play the right roles in the tape.
      // 7.

    NV.extraInputObject = 
      runExtraInputObject(NV,
    			  extraInputDummyInput,
    			  extraInputResults);

    /* TO-DO: Make the runExtraInputObject actually copy the values, so that valid values go into the model for taping.*/
    // Actually I think this could go with above steps.
    
    // During recording this will not put values in extraInputResults
    // but I don't think that will be a problem.  We are recording via
    // ADtape.Dependent(X, Y), which does not call Forward(0) for values.
    // 8.
    
    NimArrValues_AD.setSize(length_extraInput);
    std::copy(extraInputResults.begin(),
	      extraInputResults.end(),
	      NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_extraInput_accessor);
  }
}

void nodeFun::update_model_after_taping(NodeVectorClassNew_derivs &NV) {
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();  
 

    NV.extraOutputObject = 
      runExtraOutputObject(NV,
    			   // extraOutputs,
    			   //extraOutputDummyResult,
    			   logProb);
    dependentVars[0] = logProb;
        
}

// Original version
void nodeFun::recordTape(NodeVectorClassNew_derivs &NV) {
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
  
  // This does not solve the nan bug when tape optimization has been done
  // Do the same for modelOutput nodes -- experimental for nan bug
  int length_modelOutput = NV.model_modelOutput_accessor.getTotalLength();
  if(length_modelOutput > 0) {
    NimArr<1, double> NimArrValues;
    NimArr<1, CppAD::AD<double> > NimArrValues_AD;
    NimArrValues.setSize(length_modelOutput);
    NimArrValues_AD.setSize(length_modelOutput);
    getValues(NimArrValues, NV.model_modelOutput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_modelOutput,
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_modelOutput_accessor);
  }
  
  // 2. Copy all wrtNodes values from model -> model_AD, AND
  // 3. Copy all wrtNodes values from model -> independentVars, AND
  // 4. independentVars should have an extra dummy element for getExtraInputs and setExtraInputs
  int length_wrt = NV.model_wrt_accessor.getTotalLength();
  int length_independent = length_wrt + 1; // extra element is dummy
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
  independentVars[ length_independent - 1 ] = 0;
  
  // 5. Copy all extraInputNodes values from model -> model_AD (ditto, may be redundant)
  
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
  
  // 6. Start taping    
  // Steps done via nodeFunInModelDLL-> need to happen in the model DLL,
  // so that the same set of CppAD globals will be used as in taping.
  setTapeIndependent(independentVars);
  
  // 7. Instantiate and call getExtraInputs atomic object, AND
  // 8. Copy arr extraInputNodes AD objects from extraInputResult -> model_AD
  if(length_extraInput > 0) {
    std::vector< CppAD::AD<double> > extraInputDummyInput(1);
    extraInputDummyInput[0] = independentVars[ length_independent - 1 ];
    std::vector< CppAD::AD<double> > extraInputResults(length_extraInput);
      // After this, the tape treats extraInputResults as a function of
      // extraInputDummy.  This means during tape use, the extraInputObject
      // functor will be called, and the extraInputResults will contain
      // node values from the model.  By copying those during taping
      // into the model_AD, those nodes play the right roles in the tape.
      // 7.

    NV.extraInputObject = 
      runExtraInputObject(NV,
    			  extraInputDummyInput,
    			  extraInputResults);

    // During recording this will no put values in extraInputResults
    // but I don't think that will be a problem.  We are recording via
    // ADtape.Dependent(X, Y), which does not call Forward(0) for values.
    // 8.
    
    NimArrValues_AD.setSize(length_extraInput);
    std::copy(extraInputResults.begin(),
	      extraInputResults.end(),
	      NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_extraInput_accessor);
  }
  // 9. Copy all wrtNodes AD objects from independentVars -> model_AD.
    if(length_wrt > 0) {
      NimArrValues_AD.setSize(length_wrt);
      std::copy(independentVars.begin(),
		independentVars.begin() + length_wrt,
		NimArrValues_AD.getPtr());
      setValues_AD_AD(NimArrValues_AD, NV.model_AD_wrt_accessor);
    }
    // 10. call calculate
    CppAD::AD<double> logProb = call_calculate_ADproxyModel( NV );
    // 11. Copy logProb to dependentVars[0]
    
    // 12. Call setModelOutputs to put AD modelOutputs into model
    /* int length_modelOutput = model_modelOutput_accessor.getTotalLength(); */
    /* NimArr<1, CppAD::AD<double> > NAV_AD; */
    /* NAV_AD.setSize(length_modelOutput); */
    /* std::vector< CppAD::AD<double> > extraOutputs(length_modelOutput); */
    /* std::vector< CppAD::AD<double> > extraOutputDummyResult(1); */

    /* getValues_AD_AD( NAV_AD, model_AD_modelOutput_accessor); */
    /* for(size_t ii = 0; ii < length_modelOutput; ++ii) { */
    /*   extraOutputs[ii] = NAV_AD[ii]; */
    /* } */
    /* std::copy(NAV_AD.getPtr(), */
    /* 	      NAV_AD.getPtr() + length_modelOutput, */
    /* 	      extraOutputs.begin()); */

    NV.extraOutputObject = 
      runExtraOutputObject(NV,
    			   // extraOutputs,
    			   //extraOutputDummyResult,
    			   logProb);
    dependentVars[0] = logProb;
    // 13. Finish taping, AND
    // 14. Call tape->optimize()
    
    finishADFun(NV.ADtape,
		independentVars,
		dependentVars);
    
        
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

atomic_extraInputObject*
nodeFun::runExtraInputObject(NodeVectorClassNew_derivs &NV,
			     std::vector< CppAD::AD<double> > &extraInputDummyInput,
			     std::vector< CppAD::AD<double> > &extraInputResult) {
  // Calling the constructor with a local object will use statics in the correct DLL
  atomic_extraInputObject* localExtraInputObject =
    new atomic_extraInputObject("extraInputObject",
				&NV);
  // And then calling operator() with the local object will again use statics in the correct DLL
  (*localExtraInputObject)(extraInputDummyInput, extraInputResult);
  // Should be ok to copy.  No further uses after taping involve statics,
  // so there should be no further issue with which DLL we are in.
  return localExtraInputObject;
}

void
nodeFun::delete_extraInputObject(NodeVectorClassNew_derivs &NV) {
  delete NV.extraInputObject;
}

atomic_extraOutputObject*
nodeFun::runExtraOutputObject(NodeVectorClassNew_derivs &NV,
			      CppAD::AD<double> &logProb) {
  size_t length_modelOutput = NV.model_modelOutput_accessor.getTotalLength();
  NimArr<1, CppAD::AD<double> > NAV_AD;
  NAV_AD.setSize(length_modelOutput);
  std::vector< CppAD::AD<double> > extraOutputs(length_modelOutput);
  std::vector< CppAD::AD<double> > extraOutputDummyResult(1);
  
  getValues_AD_AD( NAV_AD, NV.model_AD_modelOutput_accessor);
  
  for(size_t ii = 0; ii < length_modelOutput; ++ii) {
    extraOutputs[ii] = NAV_AD[ii];
  }
  
  // Calling the constructor with a local object will use statics in the correct DLL
  atomic_extraOutputObject* localExtraOutputObject =
    new atomic_extraOutputObject("extraOutputObject",
				&NV);
  // And then calling operator() with the local object will again use statics in the correct DLL
  (*localExtraOutputObject)(extraOutputs, extraOutputDummyResult);
  // Should be ok to copy.  No further uses after taping involve statics,
  // so there should be no further issue with which DLL we are in.
  logProb += extraOutputDummyResult[0];
  return localExtraOutputObject;
}

void
nodeFun::delete_extraOutputObject(NodeVectorClassNew_derivs &NV) {
  delete NV.extraOutputObject;
}

CppAD::AD<double> nodeFun::call_calculate_ADproxyModel(NodeVectorClassNew_derivs &NV) {
  CppAD::AD<double> logProb = ::calculate_ADproxyModel( NV );
  return logProb;
}


// This makes a proxy "y = f(x)"
// where y are a set of variables we want to get into the
// calculations without introducing them via CppAD::Independent.
// This doesn't make much difference for 0-order forward performance
// but should simplify use of higher-order forward and reverse
// calls.
// By definition any derivatives are zero for this proxy f() 
atomic_extraInputObject::atomic_extraInputObject(const std::string& name,
						 NodeVectorClassNew_derivs* NV) :
  CppAD::atomic_base<double>(name, bool_sparsity_enum),
    NV_(NV)
    { }

bool atomic_extraInputObject::forward(
				      size_t                    p ,
				      size_t                    q ,
				      const ADvector<bool>&      vx ,
				      ADvector<bool>&      vy ,
				      const ADvector<double>&    tx ,
				      ADvector<double>&    ty
				      )
{
    
  // std::cout<<"Entering atomic_extraInputObject with p = "<<p<<" and q = "<<q<<std::endl;
  // std::cout<<"vx.size() = "<<vx.size()<<" and vy.size() = "<<vy.size()<<std::endl;
      
  // return flag
  bool ok = q <= 1; // max order currently implemented
  if( ! ok )
    return ok;
    
  if(vx.size() > 0) {// only true for Forward(0)
    for(size_t i = 0; i < vx.size(); ++i)
      std::cout<<vx[i]<<" ";
    std::cout<<std::endl;
    for(size_t i = 0; i < vy.size(); ++i)
      vy[i] = true;
  }

  size_t lengthExtraInput = vy.size();
  size_t length_extraNodes_accessor = NV_->model_extraInput_accessor.getTotalLength();
  if(length_extraNodes_accessor < lengthExtraInput) {
    // We should really check that the length is the right multiple of q
    std::cout<<"Problem: length_extraNodes_accessor < vy.size()"<<std::endl;
    return false;
  }

  // 0th-order
  if( p <= 0 ) {
    // Get values from model.
    //    std::cout<<"copying extra values"<<std::endl;
    NimArr<1, double> NimArr_ty;
    NimArr_ty.setSize(length_extraNodes_accessor);
    getValues(NimArr_ty, NV_->model_extraInput_accessor);
    for(size_t i = 0; i < length_extraNodes_accessor; ++i) {
      ty[i*(q+1) + 0] = NimArr_ty[i];
    }
  }
  // If 0 is the max order, return
  if( q <= 0 )
    return ok;
    
  // std::cout<<"got to order 1"<<std::endl;
    
  if( p <= 1 ) {
    for(size_t i = 0; i < length_extraNodes_accessor; ++i) {
      ty[i*(q+1) + 1] = 0.;
    }
  }    
  //  std::cout<<"finished order 1"<<std::endl;
    
    
  if( q <= 1 )
    return ok;
    
  std::cout<<"uh oh, got to order >1"<<std::endl;
    
  return false; 
}

// reverse
bool
atomic_extraInputObject::reverse(
				 size_t                    q ,
				 const ADvector<double>&    tx ,
				 const ADvector<double>&    ty ,
				 ADvector<double>&    px ,
				 const ADvector<double>&    py
				 )
{
  size_t n = px.size();
  //std::cout<<"in atomic_extraInputObject::reverse"<<std::endl;
  for(size_t i = 0; i < n; ++i)
    px[i] = 0;
  return true;
}

bool
atomic_extraInputObject::for_sparse_jac(
					size_t                     q ,
					const ADvector<bool>&   r ,
					ADvector<bool>&         s ,
					const ADvector<double>&      x )
{
  // std::cout<<"in for_sparse_jac r.size() = "<<r.size()<<" s.size() = "<<s.size()<<std::endl;
  // for(size_t j = 0; j < r.size(); ++j)
  //   std::cout<<r[j]<<" ";
  // std::cout<<std::endl;
  // sparsity for first row of S(x) = f'(x) * R
  for(size_t j = 0; j < s.size(); j++)
    s[ j ] = true;
  
  return true;
}

bool
atomic_extraInputObject::rev_sparse_jac(
					size_t                     q ,
					const ADvector<bool>&   rt ,
					ADvector<bool>&         st ,
					const ADvector<double>&      x )
{ 
  // std::cout<<"in rev_sparse_jac rt.size() = "<<rt.size()<<" st.size() = "<<st.size()<<std::endl;
  // for(size_t j = 0; j < rt.size(); ++j)
  //   std::cout<<rt[j]<<" ";
  // std::cout<<std::endl;

  // sparsity for first row of S(x) = f'(x) * R
  for(size_t j = 0; j < st.size(); j++)
    st[ j ] = true;
  
  return true;
}



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
						   NodeVectorClassNew_derivs* NV) :
  CppAD::atomic_base<double>(name, bool_sparsity_enum),
  NV_(NV)
{ }


bool atomic_extraOutputObject::forward(
				      size_t                    p ,
				      size_t                    q ,
				      const ADvector<bool>&      vx ,
				      ADvector<bool>&      vy ,
				      const ADvector<double>&    tx ,
				      ADvector<double>&    ty
				      )
{
  // std::cout<<"Entering atomic_extraOutputObject with p = "<<p<<" and q = "<<q<<std::endl;
  // std::cout<<"tx.size() = "<<tx.size()<<" and ty.size() = "<<ty.size()<<std::endl;
      
  // return flag
  bool ok = true;
    
  if(vx.size() > 0) {// only true for Forward(0)
    for(unsigned int i = 0; i < vy.size(); ++i) // should have size 1, but we'll handle anything.
      vy[i] = true;
  }
      
  size_t length_modelOutput_accessor = NV_->model_modelOutput_accessor.getTotalLength();
  // std::cout << "length_modelOutput_accessor =" << length_modelOutput_accessor << std::endl;
  // if(length_modelOutput_accessor < tx.size()) {
  //   std::cout<<"Problem: length_modelOutput_accessor < vx.size()"<<std::endl;
  // }
  
  // 0th-order
  if( p <= 0 ) {
    // Put values in model.
    //    std::cout<<"copying modelOutputs to model"<<std::endl;
    NimArr<1, double> NimArr_tx;
    NimArr_tx.setSize(length_modelOutput_accessor);
    for(size_t i = 0; i < length_modelOutput_accessor; ++i) {
      NimArr_tx[i] = tx[i*(q+1) + 0];
      // std::cout<<NimArr_tx[i]<<" ";
    }
    // std::cout<<std::endl;
    setValues(NimArr_tx, NV_->model_modelOutput_accessor);
  }

  for(unsigned int i = 0; i < ty.size(); ++i) {
      ty[i] = 0.;
  }

  return ok;
}

// reverse
bool
atomic_extraOutputObject::reverse(
				 size_t                    q ,
				 const ADvector<double>&    tx ,
				 const ADvector<double>&    ty ,
				 ADvector<double>&    px ,
				 const ADvector<double>&    py
				 )
{
  size_t n = px.size();
  //  std::cout<<"in atomic_extraOutputObject::reverse"<<std::endl;
  for(size_t i = 0; i < n; ++i)
    px[i] = 0;
  return true;
}


bool
atomic_extraOutputObject::for_sparse_jac(
				       size_t                     q ,
				       const ADvector<bool>&   r ,
				       ADvector<bool>&         s ,
				       const ADvector<double>&      x )
{ 
  // sparsity for first row of S(x) = f'(x) * R
  for(size_t j = 0; j < s.size(); j++)
    s[ j ] = true;
  
  return true;
}

bool
atomic_extraOutputObject::rev_sparse_jac(
					 size_t                     q ,
					 const ADvector<bool>&   rt ,
					 ADvector<bool>&         st ,
					 const ADvector<double>&      x )
{ 
  // sparsity for first row of S(x) = f'(x) * R
  for(size_t j = 0; j < st.size(); j++)
    st[ j ] = true;
  
  return true;
}
