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
  NimArr<1, double> getGradient(nimbleCppADinfoClass &ADinfo) {
    printf("inside getGradient\n");
    NimArr<1, double> ans = vectorDouble_2_NimArr(ADinfo.ADtape->RevOne(ADinfo.independentVars, 0));
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

class NIMBLE_ADCLASS : public NamedObjects, public pointedToBase {
public:
  NimArr<1, double> value;
  NimArr<2, double> gradient;
  NimArr<3, double> hessian;
  NimArr<4, double> thirdDerivs;

  SEXP RObjectPointer;
  bool RCopiedFlag;
void  copyFromSEXP ( SEXP S_nimList_ );
SEXP  copyToSEXP (  );
void  createNewSEXP (  );
void  resetFlags (  );
 NIMBLE_ADCLASS (  );
};
// use template for integer dim of nimArr?
// seems like we'll always want to use combos of 0 and 1 to specify directions for derivatives