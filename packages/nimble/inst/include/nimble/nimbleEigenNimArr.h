#ifndef __NIMBLE_EIGEN_NIMARR
#define __NIMBLE_EIGEN_NIMARR

#include <iostream>
#include <vector>
#include <cstdlib>
#include<Rmath.h>

template<typename NimArrOutput, typename VectorInput>
void assignVectorToNimArr(NimArrOutput &output, const VectorInput &input) {
  if(output.isMap()) {
    if(output.size() != input.size()) {
      std::cout<<"Error: mismatched sizes in assignment from which()\n";
      return;
    }
  } else {
    output.setSize(input.size()); // would need same scalar types before trying memcpy.
  }
  for(int i = 0; i < input.size(); i++) output[i] = input[i];
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

template<typename NimArrOutput, typename DerivedX, typename DerivedTimes>
  void setRepVectorTimes(NimArrOutput &output, const DerivedX &EigenX, const DerivedTimes &EigenTimes) { // uses for rep(when the "times" argument is a vector, so we pre-compute the result.  No need to worry about length_out because if it is there only times[1] is used, and we won't be in this function at all.
  typedef typename Eigen::internal::traits<DerivedX>::Scalar Scalar;
  std::vector<Scalar> ans;
  if(!(EigenX.size() == EigenTimes.size())) std::cout<<"Error: x and times have different lengths\n";
  unsigned int i, j, thisTimes;
  Scalar thisVal;
  for(i = 0; i < EigenX.size(); i++ ) {
    thisTimes = nimble_eigen_coeff_impl< nimble_eigen_traits<DerivedTimes>::LinearAccessBit, typename Eigen::internal::traits<DerivedTimes>::Scalar, DerivedTimes, unsigned int >::getCoeff(EigenTimes, i);
    thisVal = nimble_eigen_coeff_impl< nimble_eigen_traits<DerivedX>::LinearAccessBit, Scalar, DerivedX, unsigned int >::getCoeff(EigenX, i);
    for(j = 0; j < thisTimes; j++) 
      ans.push_back(thisVal);
  }
  assignVectorToNimArr<NimArrOutput, std::vector<Scalar> >(output, ans);
  return;
}

// split above into two pieces
/* template<typename DerivedBool> */
/* std::vector<int> which2vector(const DerivedBool &BoolArg) { */
/*   std::vector<int> ans; */
/*   ans.reserve(BoolArg.size()); */
/*   bool nextBool; */
/*   for(unsigned int i = 0; i < BoolArg.size(); i++ ) { */
/*     nextBool = nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedBool>::Flags & LinearAccessBit, int, DerivedBool, typename Eigen::internal::traits<DerivedBool>::Index >::getCoeff(BoolArg, i); */
/*     if(nextBool) ans.push_back(i + 1); // That 1 makes it one-based indexing, which will then be adjusted back to zero-based when used for indexing something else */
/*   } */
/*   return(ans); */
/* } */

// instead of setWhich(ans, eigen_expression)
// we'll need assignVectorToNimArr(ans, which2vector(eigen_expression))
// that means that from ans <- which(expression)
// we need to create ans <- assignVectorToNimArr(which2vector(eigen_expression))
// But I'm not sure that's any better than before.
// Maybe instead stick with setWhich, setRepVectorTimes, etc. but inside those use assignVectorToNimArr.


#endif
