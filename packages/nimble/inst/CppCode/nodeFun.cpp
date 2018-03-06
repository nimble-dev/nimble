#include "nimble/nodeFun.h"
#include "nimble/accessorClasses.h"

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
  CppAD::ADFun<double> localADtape;
  localADtape.Dependent(independentVars, dependentVars);
  std::cout<<"Tape memory size before optimize: "<< localADtape.size_op_seq() <<std::endl;
  localADtape.optimize(); //("no_compare_op") makes almost no difference;
  std::cout<<"Tape memory size after optimize: "<< localADtape.size_op_seq() <<std::endl;
  // then copy it back to object from another DLL
  ADtape = localADtape;
  ADtape.capacity_order(2); // makes only a tiny improvement
}

void nodeFun::runTape(CppAD::ADFun< double > &ADtape,
		      std::vector< double > &independentVars,
		      std::vector< double > &dependentVars) {
  dependentVars = ADtape.Forward(0, independentVars);
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
  std::cout<<"NEED TO ADD CLEAR DESTRUCTION OF extraInputObject"<<std::endl;
  return localExtraInputObject;
}

CppAD::AD<double> nodeFun::call_calculate_ADproxyModel(NodeVectorClassNew_derivs &NV) {
  CppAD::AD<double> logProb = ::calculate_ADproxyModel( NV );
  return logProb;
}



bool atomic_extraInputObject::forward(
				      size_t                    p ,
				      size_t                    q ,
				      const ADvector<bool>&      vx ,
				      ADvector<bool>&      vy ,
				      const ADvector<double>&    tx ,
				      ADvector<double>&    ty
				      )
{
    
  std::cout<<"Entering atomic with p = "<<p<<" and q = "<<q<<std::endl;
  std::cout<<"tx.size() = "<<tx.size()<<" and ty.size() = "<<ty.size()<<std::endl;

      
  // return flag
  bool ok = q <= 1; // max order currently implement
  if( ! ok )
    return ok;
    
  if(vx.size() > 0) {// only true for Forward(0)
    for(unsigned int i = 0; i < vy.size(); ++i)
      vy[i] = true;
  }
      
  int lengthExtraInput = vy.size();
  int length_extraNodes_accessor = NV_->model_extraInput_accessor.getTotalLength();
  if(length_extraNodes_accessor < lengthExtraInput) {
    // We should really check that the length is the right multiple of q
    std::cout<<"Problem: length_extraNodes_accessor < vy.size()"<<std::endl;
    return false;
  }

  // 0th-order
  if( p <= 0 ) {
    // Get values from model.
    NimArr<1, double> NimArr_ty;
    NimArr_ty.setSize(length_extraNodes_accessor);
    getValues(NimArr_ty, NV_->model_extraInput_accessor);
    for(unsigned int i = 0; i < length_extraNodes_accessor; ++i) {
      ty[i*(q+1) + 0] = NimArr_ty[i];
    }
  }
  // If 0 is the max order, return
  if( q <= 0 )
    return ok;
    
  std::cout<<"got to order 1"<<std::endl;
    
  if( p <= 1 ) {
    for(unsigned int i = 0; i < length_extraNodes_accessor; ++i) {
      ty[i*(q+1) + 1] = 0.;
    }
  }    
  std::cout<<"finished order 1"<<std::endl;
    
    
  if( q <= 1 )
    return ok;
    
  std::cout<<"uh oh, got to order >1"<<std::endl;
    
  return false; 
}
