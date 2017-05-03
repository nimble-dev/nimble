/* Definitions only to be included when a nimbleFunction needs CppAD */
#include <cppad/cppad.hpp>
#include <nimble/RcppNimbleUtils.h>
#include <nimble/EigenTypedefs.h>
//#include <nimble/RcppNimbleUtils.h>
//#include <nimble/NimArr.h>
//#include <nimble/smartPtrs.h>
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
   nimSmartPtr<NIMBLE_ADCLASS> getDerivs(nimbleCppADinfoClass &ADinfo) {
		nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
		std::size_t n = length(ADinfo.independentVars); // dim of independent vars

		ADlist->value = vectorDouble_2_NimArr(ADinfo.ADtape->Forward(0, ADinfo.independentVars));
		std::size_t q = ADlist->value.size();
		ADlist->gradient.initialize(0, false, n, q);
		ADlist->hessian.initialize(0, false, n, n, q);
		//ADlist->thirdDerivs.initialize(0, false, n, n, n, q);
		
		vector<double> gradient_ans (n*q, -1);
		vector<double> hessian_ans (n*n*q, -1);
		//vector<double> thirdDeriv_ans (n*n*n*q, -1);
		
		vector<double> cppad_hessOut;
		//vector<double> cppad_thirdDerivOut;

		for(size_t dy_ind = 0; dy_ind < q; dy_ind++){
			std::vector<double> w (q, 0);
			w[dy_ind] = 1;
			for(size_t dx1_ind = 0; dx1_ind < n; dx1_ind++){
				std::vector<double> x1 (n, 0); // vector specifying first derivatives.  first specify coeffs for first dim of s across all directions r, then second dim, ...
				x1[dx1_ind] = 1;

				ADinfo.ADtape->Forward(1, x1); // may want separate case for r=1?	
				cppad_hessOut = ADinfo.ADtape->Reverse(2, w);
			
				// printf("thirdDerivOut length: %d \n", length(cppad_thirdDerivOut));
				// for(size_t i = 0; i < length(cppad_thirdDerivOut); i++){
					// printf("thirdDeriv element %d is %f \n", i+1, cppad_thirdDerivOut[i]);
				// }
				
				for(size_t i = 0; i < n; i++){
					gradient_ans[n*dy_ind + i] = cppad_hessOut[i*2 + 0];
					hessian_ans[n*n*dy_ind + n*dx1_ind + i]  = cppad_hessOut[i*2 + 1];
				}
			}
		}
		
		
		std::copy(gradient_ans.begin(), gradient_ans.end(), ADlist->gradient.getPtr());
		std::copy(hessian_ans.begin(), hessian_ans.end(), ADlist->hessian.getPtr());
		return(ADlist);
  } 
/* NimArr<1, double> getGradient(nimbleCppADinfoClass &ADinfo) {
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
  */

};
