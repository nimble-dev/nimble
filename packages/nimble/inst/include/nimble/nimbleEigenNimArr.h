#ifndef __NIMBLE_EIGEN_NIMARR
#define __NIMBLE_EIGEN_NIMARR

#include <iostream>
#include <vector>
#include <cstdlib>
#include<Rmath.h>
#include "nimble/accessorClasses.h"

template<typename Derived>
double calculate(NodeVectorClassNew &nodes, const Derived &indices) {
  double ans(0);
  const vector<oneNodeUseInfo> &useInfoVec = nodes.getUseInfoVec();
  int len = indices.size();
  //  vector<oneNodeUseInfo>::const_iterator iNode(useInfoVec.begin());
  //  vector<oneNodeUseInfo>::const_iterator iNodeEnd(useInfoVec.end());
  int thisIndex;
  for(int i = 0; i < len; ++i) {
    thisIndex = indices(i)-1; // indices is R-based
    ans += useInfoVec[ thisIndex ].nodeFunPtr->calculateBlock(useInfoVec[ thisIndex ].useInfo );
  }
    //    ans += iNode->nodeFunPtr->calculateBlock(iNode->useInfo);
  return(ans);
}


template<typename NimArrOutput, typename VectorInput>
void assignVectorToNimArr(NimArrOutput &output, const VectorInput &input) {
  if(output.isMap()) {
    if(output.size() != input.size()) {
      std::cout<<"Error: mismatched sizes in assignment from which()\n";
      return;
    }
    for(int i = 0; i < input.size(); i++) output[i] = input[i];
  } else {
    output.setSize(input.size(), false, false); // would need same scalar types before trying memcpy.
    for(int i = 0; i < input.size(); i++) output.valueNoMap(i) = input[i];
  }
  return;
}

template<typename NimArrOutput, typename DerivedBool>
  void setWhich(NimArrOutput &output, const DerivedBool &BoolArg) {
  std::vector<int> ans;
  ans.reserve(BoolArg.size());
  bool nextBool;
  for(unsigned int i = 0; i < BoolArg.size(); i++ ) {
    nextBool = nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedBool>::Flags & LinearAccessBit, int, DerivedBool, typename Eigen::internal::traits<DerivedBool>::Index >::getCoeff(BoolArg, i);
    if(nextBool) ans.push_back(i + 1); // That 1 makes it one-based indexing, which will then be adjusted back to zero-based when used for indexing something else
  }
  assignVectorToNimArr<NimArrOutput, std::vector<int> >(output, ans);
  /* if(output.isMap()) { */
  /*   if(output.size() != ans.size()) { */
  /*     std::cout<<"Error: mismatched sizes in assignment from which()\n"; */
  /*     return; */
  /*   } */
  /* } else { */
  /*   output.setSize(ans.size(), false, false); // would need same scalar types before trying memcpy. */
  /* } */
  /* for(int i = 0; i < ans.size(); i++) output[i] = ans[i]; */
  return;
}

template<typename NimArrOutput, typename DerivedX, typename DerivedTimes>
  void setRepVectorTimes(NimArrOutput &output, const DerivedX &EigenX, const DerivedTimes &EigenTimes) { // uses for rep(when the "times" argument is a vector, so we pre-compute the result.  No need to worry about length_out because if it is there only times[1] is used, and we won't be in this function at all.
  typedef typename Eigen::internal::traits<DerivedX>::Scalar Scalar;
  std::vector<Scalar> ans;
  if(!(EigenX.size() == EigenTimes.size())) std::cout<<"Error: x and times have different lengths\n";
  unsigned int i, j, thisTimes;
  Scalar thisVal;
  for(i = 0; i < EigenX.size(); i++ ) {
    thisTimes = nimble_eigen_coeff_impl< nimble_eigen_traits<DerivedTimes>::nimbleUseLinearAccess, typename Eigen::internal::traits<DerivedTimes>::Scalar, DerivedTimes, unsigned int >::getCoeff(EigenTimes, i);
    thisVal = nimble_eigen_coeff_impl< nimble_eigen_traits<DerivedX>::nimbleUseLinearAccess, Scalar, DerivedX, unsigned int >::getCoeff(EigenX, i);
    for(j = 0; j < thisTimes; j++) 
      ans.push_back(thisVal);
  }
  assignVectorToNimArr<NimArrOutput, std::vector<Scalar> >(output, ans);
  return;
}

#endif
