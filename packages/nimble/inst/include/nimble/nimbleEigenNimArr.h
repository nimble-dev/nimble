#ifndef __NIMBLE_EIGEN_NIMARR
#define __NIMBLE_EIGEN_NIMARR

#include <iostream>
#include <vector>
#include <cstdlib>
#include<Rmath.h>

template<typename NimArrOutput, typename DerivedBool>
  void setWhich(NimArrOutput &output, const DerivedBool &BoolArg) {
  std::vector<int> ans;
  ans.reserve(BoolArg.size());
  bool nextBool;
  for(unsigned int i = 0; i < BoolArg.size(); i++ ) {
    nextBool = nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedBool>::Flags & LinearAccessBit, int, DerivedBool, typename Eigen::internal::traits<DerivedBool>::Index >::getCoeff(BoolArg, i);
    if(nextBool) ans.push_back(i + 1); // That 1 makes it one-based indexing, which will then be adjusted back to zero-based when used for indexing something else
  }
  if(output.isMap()) {
    if(output.size() != ans.size()) {
      std::cout<<"Error: mismatched sizes in assignment from which()\n";
      return;
    }
  } else {
    output.setSize(ans.size()); // would need same scalar types before trying memcpy.
  }
  for(int i = 0; i < ans.size(); i++) output[i] = ans[i];
  return;
}


#endif
