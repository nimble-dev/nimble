#include "nimble/nodeFun.h"
#include "nimble/accessorClasses.h"

void nodeFun::setTapeIndependent(std::vector< CppAD::AD<double> > &independentVars) {
  CppAD::Independent(independentVars);
}
void nodeFun::finishADFun(CppAD::ADFun< double > &ADtape,
			  std::vector< CppAD::AD<double> > &independentVars,
			  std::vector< CppAD::AD<double> > &dependentVars) {
  // make it locally to get the right globals during recording and playback
  CppAD::ADFun<double> localADtape(independentVars, dependentVars);
  std::cout<<"Tape memory size before optimize: "<< localADtape.size_op_seq() <<std::endl;
  localADtape.optimize();
  std::cout<<"Tape memory size after optimize: "<< localADtape.size_op_seq() <<std::endl;
  // then copy it back to object from another DLL
  ADtape = localADtape;
}

CppAD::AD<double> nodeFun::call_calculate_ADproxyModel(NodeVectorClassNew_derivs &NV) {
  CppAD::AD<double> logProb = ::calculate_ADproxyModel( NV );
  return logProb;
}
