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
  CppAD::ADFun<double> localADtape(independentVars, dependentVars);
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

CppAD::AD<double> nodeFun::call_calculate_ADproxyModel(NodeVectorClassNew_derivs &NV) {
  CppAD::AD<double> logProb = ::calculate_ADproxyModel( NV );
  return logProb;
}
