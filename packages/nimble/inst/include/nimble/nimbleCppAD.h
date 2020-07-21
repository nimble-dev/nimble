/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

/* Definitions only to be included when a nimbleFunction needs CppAD */
#include <cppad/cppad.hpp>
#include <nimble/RcppNimbleUtils.h>
#include <nimble/predefinedNimbleLists.h>
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

template<class Type>
Type nimDerivs_dnorm(Type x, Type mean, Type sd, int give_log=0)
{
  Type logres;
  logres=-log(Type(sqrt(2*M_PI))*sd)-Type(.5)*pow((x-mean)/sd,2);
  if(give_log)return logres; else return exp(logres);
}


/* nimbleFunctionCppADbase is a base class to be inherited by all
   CppAD-enabled nimbleFunctions. Some of these functions might 
   make more sense as stand-alone functions.  Let's see. */
class nimbleFunctionCppADbase {
 public:
   nimSmartPtr<NIMBLE_ADCLASS> getDerivs(nimbleCppADinfoClass &ADinfo, NimArr<1, double> &derivOrders) {
		nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
		std::size_t n = length(ADinfo.independentVars); // dim of independent vars
		
		int orderSize = derivOrders.size();
		double array_derivOrders[orderSize];
		std::memcpy(array_derivOrders, derivOrders.getPtr(), orderSize*sizeof(double));
		int maxOrder = *std::max_element(array_derivOrders, array_derivOrders + orderSize);
		bool ordersFound[3] = {false};

		for(int i = 0; i < orderSize; i++){
			if( (array_derivOrders[i] > 2) | (array_derivOrders[i] < 0)){ 
				printf("Error: Derivative orders must be between 0 and 2.\n");
				return(ADlist);
			}
			ordersFound[static_cast<int>(array_derivOrders[i])] = true;
		}
		
		vector<double> value_ans = ADinfo.ADtape->Forward(0, ADinfo.independentVars);
		
		if(ordersFound[0] == true){
			ADlist->value = vectorDouble_2_NimArr(value_ans);
			if(maxOrder == 0){
				return(ADlist);
			}
		}
		
		std::size_t q = length(value_ans);

		if(ordersFound[1] == true){
			ADlist->gradient.initialize(0, false, n, q);
		}
		if(ordersFound[2] == true){
			ADlist->hessian.initialize(0, false, n, n, q);
		}
		
		vector<double> cppad_derivOut;
		vector<double> hessian_ans (n*n*q, -1);
		vector<double> gradient_ans (n*q, -1);

		for(size_t dy_ind = 0; dy_ind < q; dy_ind++){
			std::vector<double> w (q, 0);
			w[dy_ind] = 1;
			if(maxOrder == 1){
				cppad_derivOut = ADinfo.ADtape->Reverse(1, w);
			}
			else{
				for(size_t dx1_ind = 0; dx1_ind < n; dx1_ind++){				
					std::vector<double> x1 (n, 0); // vector specifying first derivatives.  first specify coeffs for first dim of s across all directions r, then second dim, ...
					x1[dx1_ind] = 1;

					ADinfo.ADtape->Forward(1, x1); // may want separate case for r=1?	
					cppad_derivOut = ADinfo.ADtape->Reverse(2, w);
				
					for(size_t i = 0; i < n; i++){
						hessian_ans[n*n*dy_ind + n*dx1_ind + i]  = cppad_derivOut[i*2 + 1];
					}
				}
			}
			if(ordersFound[1] == true){
				for(size_t i = 0; i < n; i++){
					gradient_ans[n*dy_ind + i] = cppad_derivOut[i*maxOrder + 0];
				}
			}
		}
		
		if(ordersFound[1] == true){
			std::copy(gradient_ans.begin(), gradient_ans.end(), ADlist->gradient.getPtr());
		}
		if(ordersFound[2] == true){
			std::copy(hessian_ans.begin(), hessian_ans.end(), ADlist->hessian.getPtr());
		}

		return(ADlist);
  } 

};
