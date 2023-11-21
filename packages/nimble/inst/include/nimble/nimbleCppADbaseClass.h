#ifndef NIMBLECPPAD_BASECLASS_H_
#define NIMBLECPPAD_BASECLASS_H_

/* Definitions only to be included when a nimbleFunction needs CppAD */
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <nimble/EigenTypedefs.h>
#include <nimble/accessorClasses.h>
#include <nimble/nodeFun.h>
#include <nimble/predefinedNimbleLists.h>
#include <cstdio>
#include <vector>
#include <algorithm>

class nimbleFunctionCppADbase {
public:
  void getDerivs(nimbleCppADinfoClass &ADinfo,
                 const NimArr<1, double> &derivOrders,
                 const NimArr<1, double> &wrtVector,
                 nimSmartPtr<NIMBLE_ADCLASS> &ansList);

  nimSmartPtr<NIMBLE_ADCLASS> getDerivs_wrapper(nimbleCppADinfoClass &ADinfo,
                                                const NimArr<1, double> &derivOrders,
                                                const NimArr<1, double> &wrtVector){
    nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
    getDerivs(ADinfo, derivOrders, wrtVector, ansList);
    return(ansList);
  }

  void getDerivs_meta(nimbleCppADinfoClass &ADinfo,
                      const NimArr<1, double> &derivOrders,
                      const NimArr<1, double> &wrtVector,
                      const nimbleCppADrecordingInfoClass &nimRecInfo,
                      nimSmartPtr<NIMBLE_ADCLASS_META> &ansList);

  nimSmartPtr<NIMBLE_ADCLASS_META> getDerivs_wrapper_meta(nimbleCppADinfoClass &ADinfo,
                                                          const NimArr<1, double> &derivOrders,
                                                          const NimArr<1, double> &wrtVector,
                                                          const nimbleCppADrecordingInfoClass &nimRecInfo){
    nimSmartPtr<NIMBLE_ADCLASS_META> ansList = new NIMBLE_ADCLASS_META;
    getDerivs_meta(ADinfo, derivOrders, wrtVector, nimRecInfo, ansList);
    return(ansList);
  }

  void getDerivs_calculate_internal(nimbleCppADinfoClass &ADinfo,
                                    //CppAD::ADFun<double>* &tapePtr,
                                    NodeVectorClassNew_derivs &nodes,
                                    const NimArr<1, double> &derivOrders,
                                    const NimArr<1, double> &wrtVector,
                                    bool do_update,
                                    bool reset,
                                    nimSmartPtr<NIMBLE_ADCLASS> ansList);
  /* This form is not actually generated in code at the time of this writing:*/
  nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(nimbleCppADinfoClass &ADinfo,
                                                  //CppAD::ADFun<double>* &tapePtr,
                                                  NodeVectorClassNew_derivs &nodes,
                                                  const NimArr<1, double> &derivOrders,
                                                  const NimArr<1, double> &wrtVector,
                                                  bool do_update,
                                                  bool reset){
    nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
    getDerivs_calculate_internal(ADinfo,// tapePtr,
                                 nodes, derivOrders, wrtVector, do_update, reset, ansList);
    return(ansList);
  }
  /* This is the form that would be generated in code, with no wrtVector*/
  nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(nimbleCppADinfoClass &ADinfo,
                                                  //CppAD::ADFun<double>* &tapePtr,
                                                  NodeVectorClassNew_derivs &nodes,
                                                  const NimArr<1, double> &derivOrders,
                                                  bool do_update,
                                                  bool reset) {
    NimArr<1, double> wrtVector; // with new default functionality, this could be set to simply length 1 with value -1
    int totlen = nodes.model_wrt_accessor.getTotalLength();
    wrtVector.setSize(totlen,
                      false,
                      false);
    for(int ii = 0; ii < totlen; ++ii) {
      wrtVector[ii] = ii + 1; // This gets -1 in use, as if it were from R.
    }
    nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
    getDerivs_calculate_internal(ADinfo, //tapePtr,
                                 nodes, derivOrders, wrtVector, do_update, reset, ansList);
    return(ansList);
  }
};


#endif // NIMBLECPPAD_BASECLASS_H_
