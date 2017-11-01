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
//#include <TMB/tmbDists.h>
#include <nimble/accessorClasses.h>
#include <nimble/EigenTypedefs.h>
#include <nimble/nodeFun.h>
//#include <nimble/RcppNimbleUtils.h>
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


/* nimbleFunctionCppADbase is a base class to be inherited by all
   CppAD-enabled nimbleFunctions. Some of these functions might 
   make more sense as stand-alone functions.  Let's see. */
class nimbleFunctionCppADbase {
 public:
   nimSmartPtr<NIMBLE_ADCLASS> getDerivs(nimbleCppADinfoClass &ADinfo, NimArr<1, double> &derivOrders, NimArr<1, double> &wrtVector) {
		nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
		std::size_t n = length(ADinfo.independentVars); // dim of independent vars
		std::size_t wrt_n = wrtVector.size(); // dim of wrt vars
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
			ADlist->gradient.initialize(0, false, q, wrt_n);
		}
		if(ordersFound[2] == true){
			ADlist->hessian.initialize(0, false, wrt_n, wrt_n, q);
		}

		vector<double> cppad_derivOut;
		vector<double> hessian_ans (wrt_n*wrt_n*q, -1);
		vector<double> gradient_ans (wrt_n*q, -1);
		for(size_t dy_ind = 0; dy_ind < q; dy_ind++){
			std::vector<double> w (q, 0);
			w[dy_ind] = 1;
			if(maxOrder == 1){
				cppad_derivOut = ADinfo.ADtape->Reverse(1, w);
			}
			else{
 				for(size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++){
					int dx1_ind = wrtVector[vec_ind] - 1;	
 					std::vector<double> x1 (n, 0); // vector specifying first derivatives.  first specify coeffs for first dim of s across all directions r, then second dim, ...
					x1[dx1_ind] = 1;
					ADinfo.ADtape->Forward(1, x1); // may want separate case for r=1?	
					cppad_derivOut = ADinfo.ADtape->Reverse(2, w);
  					for(size_t vec_ind2 = 0; vec_ind2 < wrt_n; vec_ind2++){
						int dx2_ind = wrtVector[vec_ind2] - 1;	
						hessian_ans[wrt_n*wrt_n*dy_ind + wrt_n*vec_ind + vec_ind2]  = cppad_derivOut[dx2_ind*2 + 1];
					}   
  				} 
			}
 			if(ordersFound[1] == true){
 				for(size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++){
					int dx1_ind = wrtVector[vec_ind] - 1;	
					gradient_ans[vec_ind*q + dy_ind] = cppad_derivOut[dx1_ind*maxOrder + 0];
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

nimSmartPtr<NIMBLE_ADCLASS>  NIM_DERIVS_CALCULATE(NodeVectorClassNew_derivs &nodes, NimArr<1, double> &derivOrders);
nimSmartPtr<NIMBLE_ADCLASS>  NIM_DERIVS_CALCULATE(NodeVectorClassNew_derivs &nodes,  int iNodeFunction, NimArr<1, double> &derivOrders);


template<class Type>
Type nimDerivs_dnorm(Type x, Type mean, Type sd, int give_log)
{
  Type logres;
  logres=-log(Type(sqrt(2*M_PI))*sd)-Type(.5)*pow((x-mean)/sd,2);
  if(give_log)return logres; else return exp(logres);
}

template<class Type> 
Type nimDerivs_dexp(Type x, Type rate, int give_log)
{
	if(!give_log)
		return CppAD::CondExpGe(x,Type(0),rate*exp(-rate*x),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(rate)-rate*x,Type(-INFINITY));
}

template<class Type> 
Type nimDerivs_dweib(Type x, Type shape, Type scale, int give_log)
{
	if(!give_log)
		return CppAD::CondExpGe(x,Type(0),shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape)),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(shape) - log(scale) + (shape-1)*(log(x)-log(scale)) - pow(x/scale,shape),Type(-INFINITY));
}

template<class Type> 
Type nimDerivs_dbinom(Type k, Type size, Type prob, int give_log)
{
	Type logres = lgamma(size+1)-lgamma(k+1)-lgamma(size-k+1)+k*log(prob)+(size-k)*log(1-prob);
	if(!give_log) return exp(logres);
	else return logres;
}

template <class Type>
Type nimDerivs_dbeta(Type x, Type shape1, Type shape2, int give_log)
{
	Type res = exp(lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1);
	if(!give_log) 
		return res;
	else 
		return CppAD::CondExpEq(x,Type(0),log(res),lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2) + (shape1-1)*log(x) + (shape2-1)*log(1-x));
}

template <class Type>
Type nimDerivs_dlogis(Type x, Type location, Type scale, int give_log)
{
	Type logres = -(x-location)/scale - log(scale) - 2*log(1+exp(-(x-location)/scale));
	if(!give_log) return exp(logres);
	else return logres;
}

template <class Type>
Type nimDerivs_dt(Type x, Type df, int give_log)
{
	Type logres = lgamma((df+1)/2) - Type(1)/2*log(df*M_PI) -lgamma(df/2) - (df+1)/2*log(1+x*x/df);
	if(!give_log) return exp(logres);
	else return logres;
}

template <class Type>
Type nimDerivs_dmulti(vector<Type> x, vector<Type> p, int give_log=0)
{
	vector<Type> xp1 = x+Type(1);
	Type logres = lgamma(x.sum() + Type(1)) - lgamma(xp1).sum() + (x*log(p)).sum();
	if(give_log) return logres;
	else return exp(logres);
}

// template<class Type>
// Type lgamma(Type x){
  // CppAD::vector<Type> tx(2);
  // tx[0] = x;
  // tx[1] = Type(0);
  // return atomic::D_lgamma(tx)[0];
// }

// template<class Type>
// Type lfactorial(Type x){
  // CppAD::vector<Type> tx(2);
  // tx[0] = x + Type(1);
  // tx[1] = Type(0);
  // return atomic::D_lgamma(tx)[0];
// }

// template<class Type>
// inline Type nimDerivs_dnbinom(const Type &x, const Type &size, const Type &prob,
		    // int give_log)
// {
  // Type n=size;
  // Type p=prob;
  // Type logres = lgamma(x+n)-lgamma(n)-lgamma(x+Type(1))+
    // n*log(p)+x*log(Type(1)-p);
  // if (give_log) return logres; else return exp(logres);
// }

// template<class Type>
// inline Type nimDerivs_dpois(const Type &x, const Type &lambda, int give_log)
// {
  // Type logres = -lambda + x*log(lambda) - lgamma(x+Type(1));
  // if (give_log) return logres; else return exp(logres);
// }

// template<class Type>
// Type nimDerivs_dgamma(Type y, Type shape, Type scale, int give_log)
// {
  // Type logres=-lgamma(shape)+(shape-Type(1.0))*log(y)-y/scale-shape*log(scale);
  // if(give_log)return logres; else return exp(logres);
// }
