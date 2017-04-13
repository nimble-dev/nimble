/* Definitions only to be included when a nimbleFunction needs CppAD */

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
};
