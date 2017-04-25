/* Definitions only to be included when a nimbleFunction needs CppAD */
#include <cppad/cppad.hpp>
#include <nimble/RcppNimbleUtils.h>
#include<vector>
#include<cstdio>
/* nimbleCppADinfoClass is the class to convey information from a nimbleFunction object 
   to generic CppAD driver wrappers like calcGradient.
   Each nimbleFunction enabled for CppAD will have an object of this class. */
class nimbleCppADinfoClass {
 public:
  std::vector<double> independentVars;
  CppAD::ADFun<double> *ADtape;
};


/* nimbleFunctionCppADbase is a base class to be inherited by all
   CppAD-enabled nimbleFunctions. Some of these functions might 
   make more sense as stand-alone functions.  Let's see. */
class nimbleFunctionCppADbase {
 public:
  vector<double> getGradient(nimbleCppADinfoClass &ADinfo) {
    printf("inside getGradient\n");
    vector<double> ans = ADinfo.ADtape->RevOne(ADinfo.independentVars, 0);
    return(ans);
  }
  
  NimArr<2, double> getHessian(nimbleCppADinfoClass &ADinfo) {
	NimArr<2, double> nimAnsMat;
	std::size_t p = length(ADinfo.independentVars);
    std::vector<int> i (p);
    std::vector<int> j (p);
	
	nimAnsMat.initialize(0, 1, 1, 1, p, p);
	for(std::size_t independentVar=0; independentVar < p; independentVar++ ) {
	  i[independentVar] = 0; 
	  j[independentVar] = independentVar; 
	}
	printf("inside getHessian\n");
	vector<double> ans = ADinfo.ADtape->RevTwo(ADinfo.independentVars, i, j);
	std::copy(ans.begin(), ans.end(), nimAnsMat.getPtr());  //could this line ever cause error?
	return(nimAnsMat);
  }
};
