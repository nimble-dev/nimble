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
#include <cppad/utility/nan.hpp>
//#include <TMB/tmbDists.h>
#include <nimble/EigenTypedefs.h>
#include <nimble/accessorClasses.h>
#include <nimble/nodeFun.h>
//#include <nimble/RcppNimbleUtils.h>
#include <nimble/predefinedNimbleLists.h>
//#include <nimble/RcppNimbleUtils.h>
//#include <nimble/NimArr.h>
//#include <nimble/smartPtrs.h>
#include <cstdio>
#include <vector>


/* nimbleCppADinfoClass is the class to convey information from a nimbleFunction
   object
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
  void                        getDerivs(nimbleCppADinfoClass &ADinfo,
                                        NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
    std::size_t n = ADinfo.independentVars.size();  // dim of independent vars

    std::size_t wrt_n = wrtVector.size();            // dim of wrt vars
    int orderSize = derivOrders.size();
    double array_derivOrders[orderSize];

    std::memcpy(array_derivOrders, derivOrders.getPtr(),
                orderSize * sizeof(double));

    int maxOrder =
        *std::max_element(array_derivOrders, array_derivOrders + orderSize);
    bool ordersFound[3] = {false};

    for (int i = 0; i < orderSize; i++) {
      if ((array_derivOrders[i] > 2) | (array_derivOrders[i] < 0)) {
        printf("Error: Derivative orders must be between 0 and 2.\n");
      }
      ordersFound[static_cast<int>(array_derivOrders[i])] = true;
    }
    vector<double> value_ans;
    value_ans = ADinfo.ADtape->Forward(0, ADinfo.independentVars);
    if (ordersFound[0] == true) {
      ansList->value = vectorDouble_2_NimArr(value_ans);
    }
    if(maxOrder > 0){
      std::size_t q = value_ans.size();
      vector<bool> infIndicators(q); 
      for(size_t inf_ind = 0; inf_ind < q; inf_ind++){
        if(((value_ans[inf_ind] == -std::numeric_limits<double>::infinity()) |
          (value_ans[inf_ind] == std::numeric_limits<double>::infinity())) | 
          (isnan(value_ans[inf_ind]))){
          infIndicators[inf_ind] = true;
        }
        else{
          infIndicators[inf_ind] = false;
        }
      }
      if (ordersFound[1] == true) {
        ansList->gradient.setSize(q, wrt_n, false, false); // setSize may be costly.  Possible to setSize outside of fxn, within chain rule algo, and only resize when necessary?
      }
      if (ordersFound[2] == true) {
        ansList->hessian.setSize(wrt_n, wrt_n, q, false, false);
      }
        vector<double> cppad_derivOut;
        for (size_t dy_ind = 0; dy_ind < q; dy_ind++) {
          std::vector<double> w(q, 0);
          w[dy_ind] = 1;
          if (maxOrder == 1) {   
            if(infIndicators[dy_ind] == false){
              cppad_derivOut = ADinfo.ADtape->Reverse(1, w);
            }
          } else {
            for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
              if(infIndicators[dy_ind] == false){
                int dx1_ind = wrtVector[vec_ind] - 1;
                std::vector<double> x1(n, 0);  // vector specifying first derivatives.
                                              // first specify coeffs for first dim
                                              // of s across all directions r, then
                                              // second dim, ...
                x1[dx1_ind] = 1;
                ADinfo.ADtape->Forward(1, x1);
                cppad_derivOut = ADinfo.ADtape->Reverse(2, w);
              }
              for (size_t vec_ind2 = 0; vec_ind2 < wrt_n; vec_ind2++) {
                if(infIndicators[dy_ind] == false){
                  int dx2_ind = wrtVector[vec_ind2] - 1;
                  ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] =
                      cppad_derivOut[dx2_ind * 2 + 1];
                }
                else{
                  ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] = 
                  CppAD::numeric_limits<double>::quiet_NaN();
                }
              }
            }
          }
          if (ordersFound[1] == true) {
            for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
              if(infIndicators[dy_ind] == false){
                int dx1_ind = wrtVector[vec_ind] - 1;
                ansList->gradient[vec_ind * q + dy_ind] =
                    cppad_derivOut[dx1_ind * maxOrder + 0];
              }
              else{
                ansList->gradient[vec_ind * q + dy_ind] =
                CppAD::numeric_limits<double>::quiet_NaN();
              }     
            }
          }
      }
    }
  };

   nimSmartPtr<NIMBLE_ADCLASS> getDerivs_wrapper(nimbleCppADinfoClass &ADinfo,
                                        NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector){
            nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
            getDerivs(ADinfo, derivOrders, wrtVector, ansList);
            return(ansList);
    }
};

nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    const NodeVectorClassNew_derivs &nodes, const NimArr<1, double> &derivOrders);
nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    const NodeVectorClassNew_derivs &nodes, const double derivOrders);
nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    NodeVectorClassNew_derivs &nodes, int iNodeFunction,
    NimArr<1, double> &derivOrders);
