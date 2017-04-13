#ifndef __ADtypes
#define __ADtypes
/* THIS WAS A DRAFT OF CONCEPTS THAT IS NOT USED.  
   I AM LEAVING IT FOR NOW BUT IT SHOULD BE CLEANED UP UNLESS 
   IT TURNS OUT TO BE USEFUL.  */

// this shuold be included before EigenTypesDefs or R.h
#include <cppad/cppad.hpp> // need to get a path into Makevars for this

using CppAD::AD;

class ADtapes { //permanent
public:
  vector<CppAD::ADFun<double> *> ADtapePtrVec;
  CppAD::ADFun<double> *getADtapePtr(int i) {return ADtapePtrVec[i];}
  //  int NumIndependentVars, NumResponseVars;
  //ADtapes(int niv, int nrv) : NumIndependentVars(niv), NumResponseVars(nrv) {};
  //int &numIndependentVars() {return(NumIndependentVars);}
  //int &numResponseVars() {return(NumResponseVars);}
  void recordTape(int NumIndependentVars, int NumResponseVars);
  virtual void callFunction(vector< AD<double> > &ADindependentVars, vector< AD<double> > &ADresponseVars)=0;
};

void ADtapes::recordTape() { //permanent
  vector< AD<double> > ADindependentVars;
  vector< AD<double> > ADresponseVars;
  //ADindependentVars.resize(NumIndependentVars); // let callFunction resize these, since it knows the relevant sizes
  //ADresponseVars.resize(NumResponseVars);
  CppAD::ADFun<double> *ADtapePtr = new CppAD::ADFun<double>;
  callFunction(ADindependentVars, ADresponseVars); // this can expand to call multiple class methods via an integer argument
  //  ADtapePtr = new CppAD::ADFun<double>(ADindependentVars, ADresponseVars);
  ADtapePtr->Dependent(ADindependentVars, ADresponseVars);
  cout<<"size before optimize = "<< ADtapePtr->size_var() <<"\n";
  ADtapePtr->optimize();
  cout<<"size after optimize = "<< ADtapePtr->size_var() <<"\n";
  ADtapePtrVec.push_back(ADtapePtr);
}

#endif
