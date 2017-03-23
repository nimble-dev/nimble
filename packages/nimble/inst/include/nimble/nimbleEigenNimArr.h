#ifndef __NIMBLE_EIGEN_NIMARR
#define __NIMBLE_EIGEN_NIMARR

#include <iostream>
#include <vector>
#include <cstdlib>
#include<Rmath.h>
#include "nimble/accessorClasses.h"

template<typename Derived>
double calculate(NodeVectorClassNew &nodes, const Derived &indices, bool logical=false) {
  double ans(0);
  const vector<oneNodeUseInfo> &useInfoVec = nodes.getUseInfoVec();
  int len = indices.size();
  if(!logical) {
    int thisIndex;
    for(int i = 0; i < len; ++i) {
      thisIndex = indices(i)-1; // indices is R-based
      ans += useInfoVec[ thisIndex ].nodeFunPtr->calculateBlock(useInfoVec[ thisIndex ].useInfo );
    }
    return(ans);
  }
  // logical = true. treat indices as a logical vector
  for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      ans += useInfoVec[ i ].nodeFunPtr->calculateBlock(useInfoVec[ i ].useInfo );
    }
  }
  return(ans);
}
  
template<typename Derived>  
  double calculateDiff(NodeVectorClassNew &nodes, const Derived &indices, bool logical=false) {
  double ans(0);
  const vector<oneNodeUseInfo> &useInfoVec = nodes.getUseInfoVec();
  int len = indices.size();
  if(!logical) {
    int thisIndex;
    for(int i = 0; i < len; ++i) {
      thisIndex = indices(i)-1; // indices is R-based
      ans += useInfoVec[ thisIndex ].nodeFunPtr->calculateDiffBlock(useInfoVec[ thisIndex ].useInfo );
    }
    return(ans);
  }
  // logical = true. treat indices as a logical vector
  for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      ans += useInfoVec[ i ].nodeFunPtr->calculateDiffBlock(useInfoVec[ i ].useInfo );
    }
  }
  return(ans);
}

template<typename Derived>  
  double getLogProb(NodeVectorClassNew &nodes, const Derived &indices, bool logical=false) {
  double ans(0);
  const vector<oneNodeUseInfo> &useInfoVec = nodes.getUseInfoVec();
  int len = indices.size();
  if(!logical) {
    int thisIndex;
    for(int i = 0; i < len; ++i) {
      thisIndex = indices(i)-1; // indices is R-based
      ans += useInfoVec[ thisIndex ].nodeFunPtr->getLogProbBlock(useInfoVec[ thisIndex ].useInfo );
    }
    return(ans);
  }
  // logical = true. treat indices as a logical vector
  for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      ans += useInfoVec[ i ].nodeFunPtr->getLogProbBlock(useInfoVec[ i ].useInfo );
    }
  }
  return(ans);
}

template<typename Derived>  
  void simulate(NodeVectorClassNew &nodes, const Derived &indices, bool logical=false) {
  const vector<oneNodeUseInfo> &useInfoVec = nodes.getUseInfoVec();
  int len = indices.size();
  if(!logical) {
   int thisIndex;
   for(int i = 0; i < len; ++i) {
    thisIndex = indices(i)-1; // indices is R-based
    useInfoVec[ thisIndex ].nodeFunPtr->simulateBlock(useInfoVec[ thisIndex ].useInfo );
   }
  } else {
  // logical = true. treat indices as a logical vector
   for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      useInfoVec[ i ].nodeFunPtr->simulateBlock(useInfoVec[ i ].useInfo );
    }
   }
  }
}


/* template<typename NimArrOutput, typename NimArrInput> */
/*   void setSizeNimArrToNimArr(NimArrOutput &target, NimArrInput &input, bool copyValues = true, bool fillZeros = true ) { */
/*   if(target.isMap()) { */
/*     _nimble_global_output << "Error from C++: using setSizeNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output); */
/*     return; */
/*   }  */
/*   if(input.numDims() != 1) { */
/*     _nimble_global_output << "Error from C++: using setSizeNimArrToNimArr with input that is not 1-dimensional.\n"; nimble_print_to_R( _nimble_global_output); */
/*   } */
/*   if(input.size() != target.numDims()) { */
/*     _nimble_global_output << "Error from C++: using setSizeNimArrToNimArr with wrong number of dimensions provided.\n"; nimble_print_to_R( _nimble_global_output); */
/*   } */
/*   switch(target.numDims()) { */
/*   case 1: */
/*     target.setSize(input[0], copyValues, fillZeros); */
/*     break; */
/*   case 2: */
/*     target.setSize(input[0], input[1], copyValues, fillZeros); */
/*     break; */
/*   case 3: */
/*     target.setSize(input[0], input[1], input[2], copyValues, fillZeros); */
/*     break; */
/*   case 4: */
/*     target.setSize(input[0], input[1], input[2], input[3], copyValues, fillZeros); */
/*     break; */
/*   default: */
/*     _nimble_global_output << "Error from C++: using setSizeNimArrToNimArr with invalid number of dimensions.\n"; nimble_print_to_R( _nimble_global_output); */
/*     break; */
/*   } */

/*   return; */
/* } */


template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, int size1) {
  if(!init) return;
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, false, false); // would need same scalar types before trying memcpy.
    if(input.size() < size1) {
      _nimble_global_output << "Warning from C++: not enough initialization values.\n"; nimble_print_to_R( _nimble_global_output);
    }
    if(input.isMap() | (output.getNimType() != input.getNimType()))
      for(int i = 0; i < size1; i++) output.valueNoMap(i) = input[i];
    else {
      memcpy(output.getPtr(), input.getPtr(), size1 * output.element_size() );
    }
  }
  return;
}

template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, int size1, int size2) {
  if(!init) return;
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, size2, false, false); // would need same scalar types before trying memcpy.
    int totsize = size1 * size2;
    if(input.size() < totsize) {
      _nimble_global_output << "Warning from C++: not enough initialization values.\n"; nimble_print_to_R( _nimble_global_output);
    }
    if(input.isMap() | (output.getNimType() != input.getNimType()))
      for(int i = 0; i < totsize; i++) output.valueNoMap(i) = input[i];
    else {
      memcpy(output.getPtr(), input.getPtr(), totsize * output.element_size() );
    }
   }
  return;
}

template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, int size1, int size2, int size3) {
  if(!init) return;
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, size2, size3, false, false); // would need same scalar types before trying memcpy.
    int totsize = size1 * size2 * size3;
    if(input.size() < totsize) {
      _nimble_global_output << "Warning from C++: not enough initialization values.\n"; nimble_print_to_R( _nimble_global_output);
    }
    if(input.isMap() | (output.getNimType() != input.getNimType()))
      for(int i = 0; i < totsize; i++) output.valueNoMap(i) = input[i];
    else {
      memcpy(output.getPtr(), input.getPtr(), totsize * output.element_size() );
    }
   }
  return;
}


template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, int size1, int size2, int size3, int size4) {
  if(!init) return;
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, size2, size3, size4, false, false); // would need same scalar types before trying memcpy.
    int totsize = size1 * size2 * size3 * size4;
    if(input.size() < totsize) {
      _nimble_global_output << "Warning from C++: not enough initialization values.\n"; nimble_print_to_R( _nimble_global_output);
    }
    if(input.isMap() | (output.getNimType() != input.getNimType()))
      for(int i = 0; i < totsize; i++) output.valueNoMap(i) = input[i];
    else {
      memcpy(output.getPtr(), input.getPtr(), totsize * output.element_size() );
    }
   }
  return;
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
