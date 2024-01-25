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

#ifndef __NIMBLE_EIGEN_NIMARR
#define __NIMBLE_EIGEN_NIMARR

#include <iostream>
#include <vector>
#include <cstdlib>
#include<Rmath.h>
#include "nimble/accessorClasses.h"

template<typename Derived>
int getNodesLength_Indices(ManyVariablesMapAccessor &MMVAPtr, const Derived &indices) {
  // vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int k = indices.size();
  //SingleVariableMapAccessBase* curSingleAccess;
  int totLength = 0;
  int thisIndex;
  for(int i = 0; i < k ; i++) {
    thisIndex = indices(i);
    //curSingleAccess = (*SMVA_Vec)[thisIndex];
    totLength += MMVAPtr.getNodeLength(thisIndex); // subtracts 1 internally - should be changed
  }
  return(totLength);
}


template<typename Derived, class T>
void ManyModelAccessIndexRange_2_nimArr(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr, const Derived &indices) {
  vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = indices.size();
  int nextNumVals, thisIndex;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++) {
    thisIndex = indices(i)-1;
    curSingleAccess = (*SMVA_Vec)[thisIndex];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd)
      PRINTF("Warning: in ManyModelAccessIndex_2_nimArr (range), accessor larger than NimArr!\n");
    SingleModelAccess_2_nimArr<T>(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
    nimCurrent += nextNumVals;
    nimCurrentOffset += nextNumVals * nimArrStride;
  }
}

// This could just be getValues, overloading other getValues prototypes
// It became getValuesIndexRange for clarity and possibility of separate handling in size processing
// But ultimately I think it is not different.
template<class Derived, class T>
void getValuesIndexRange(NimArr<1, T> &nimArr, ManyVariablesMapAccessor &MVA, const Derived &indices){
  //  int sizeFromIndices = getNodesLength_Indices(MVA, indices);
  //if(!(nimArr.isMap())) nimArr.setSize(sizeFromIndices, false, false);
  ManyModelAccessIndexRange_2_nimArr(MVA, nimArr, indices);
}


template<typename Derived, class T>
void nimArr_2_ManyModelAccessIndexRange(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr, const Derived &indices) {
  vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = indices.size();
  int nextNumVals, thisIndex;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++) {
    thisIndex = indices(i)-1;
    curSingleAccess = (*SMVA_Vec)[thisIndex];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd)
      PRINTF("Warning: in nimArr_2_ManyModelAccessIndexRange (range), accessor larger than NimArr!\n");
    nimArr_2_SingleModelAccess<T>(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
    nimCurrent += nextNumVals;
    nimCurrentOffset += nextNumVals * nimArrStride;
  }
}


template<class Derived, class T>
void setValuesIndexRange(NimArrBase<T> &nimArr, ManyVariablesMapAccessor &MVA, const Derived &indices){
  nimArr_2_ManyModelAccessIndexRange(MVA, nimArr, indices);
}


template<typename Derived>
double calculate(NodeVectorClassNew &nodes, const Derived &indices, bool logical) {
  double ans(0);
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  int len = indices.size();
  if(!logical) {
    int thisIndex;
    for(int i = 0; i < len; ++i) {
      thisIndex = indices(i)-1; // indices is R-based
      ans += instructions[ thisIndex ].nodeFunPtr->calculateBlock(instructions[ thisIndex ].operand);
    }
    return(ans);
  }
  // logical = true. treat indices as a logical vector
  for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      ans += instructions[ i ].nodeFunPtr->calculateBlock(instructions[ i ].operand);
    }
  }
  return(ans);
}
  
template<typename Derived>  
  double calculateDiff(NodeVectorClassNew &nodes, const Derived &indices, bool logical) {
  double ans(0);
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  int len = indices.size();
  if(!logical) {
    int thisIndex;
    for(int i = 0; i < len; ++i) {
      thisIndex = indices(i)-1; // indices is R-based
      ans += instructions[ thisIndex ].nodeFunPtr->calculateDiffBlock(instructions[ thisIndex ].operand);
    }
    return(ans);
  }
  // logical = true. treat indices as a logical vector
  for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      ans += instructions[ i ].nodeFunPtr->calculateDiffBlock(instructions[ i ].operand);
    }
  }
  return(ans);
}

template<typename Derived>  
  double getLogProb(NodeVectorClassNew &nodes, const Derived &indices, bool logical) {
  double ans(0);
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  int len = indices.size();
  if(!logical) {
    int thisIndex;
    for(int i = 0; i < len; ++i) {
      thisIndex = indices(i)-1; // indices is R-based
      ans += instructions[ thisIndex ].nodeFunPtr->getLogProbBlock(instructions[ thisIndex ].operand);
    }
    return(ans);
  }
  // logical = true. treat indices as a logical vector
  for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      ans += instructions[ i ].nodeFunPtr->getLogProbBlock(instructions[ i ].operand);
    }
  }
  return(ans);
}

template<typename Derived>  
  void simulate(NodeVectorClassNew &nodes, const Derived &indices, bool logical) {
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  int len = indices.size();
  if(!logical) {
   int thisIndex;
   for(int i = 0; i < len; ++i) {
    thisIndex = indices(i)-1; // indices is R-based
    instructions[ thisIndex ].nodeFunPtr->simulateBlock(instructions[ thisIndex ].operand);
   }
  } else {
  // logical = true. treat indices as a logical vector
   for(int i = 0; i < len; ++i) {
    if(indices(i)) {
      instructions[ i ].nodeFunPtr->simulateBlock(instructions[ i ].operand);
    }
   }
  }
}


template<typename NimArrOutput, typename NimArrInput>
  void setSizeNimArrToNimArr(NimArrOutput &target, NimArrInput &input, bool copyValues = true, bool fillZeros = true ) {
  if(target.isMap()) {
    _nimble_global_output << "Warning from C++: using setSize with a map.\n"; nimble_print_to_R( _nimble_global_output);
    return;
  }
  if(input.numDims() != 1) {
    _nimble_global_output << "Warning from C++: using setSize with input that is not 1-dimensional.\n"; nimble_print_to_R( _nimble_global_output);
  }
  if(input.size() != target.numDims()) {
    _nimble_global_output << "Warning from C++: using setSize with wrong number of dimensions provided.\n"; nimble_print_to_R( _nimble_global_output);
  }
  vector<int> newSizes;
  for(int i = 0; i < target.numDims(); ++i) newSizes.push_back(input[i]); // input could be a map with strides
  target.setSize(newSizes, copyValues, fillZeros);
  return;
}

// This is called from one of the assignNimArrToNimArr functions, which are generated from numeric / integer / logical / matrix / array with NON-SCALAR initialization
template<typename NimArrOutput, typename NimArrInput>
  void copyNimArrToNimArrInternal(NimArrOutput &output, NimArrInput &input, int totSize, bool fillZeros, bool recycle) {
  int sizeToCopy;
  if(input.size() < totSize) {
    //     _nimble_global_output << "Warning from C++: not enough initialization values.\n"; nimble_print_to_R( _nimble_global_output);
      sizeToCopy = input.size();
  } else {
    if(input.size() > totSize) {
      //   _nimble_global_output << "Warning from C++: too many initialization values.\n"; nimble_print_to_R( _nimble_global_output);
    }
    sizeToCopy = totSize;
  }
  if(input.isMap() || (output.getNimType() != input.getNimType())) { 
    for(int i = 0; i < sizeToCopy; i++) output.valueNoMap(i) = input[i];
    if(sizeToCopy < totSize) {
      if(recycle) {
	int iStart = sizeToCopy;
	while(iStart + sizeToCopy <= totSize) {
	  for(int i = 0; i < sizeToCopy; i++) output.valueNoMap(iStart + i) = input[i];
	  iStart += sizeToCopy;
	}
	if(iStart < totSize)
	  for(int i = 0; i < (totSize - iStart); i++) output.valueNoMap(iStart + i) = input[i];
      } else {
	if(fillZeros) 
	  for(int i = sizeToCopy; i < totSize; i++) output.valueNoMap(i) = 0;
      }
    }
  } else {
    memcpy(output.getPtr(), input.getPtr(), sizeToCopy * output.element_size() );
    if(sizeToCopy < totSize) {
      if(recycle) {
	int iStart = sizeToCopy;
	while(iStart + sizeToCopy <= totSize) {
	  memcpy(output.getPtr() + iStart, input.getPtr(), sizeToCopy * output.element_size() );
	  iStart += sizeToCopy;
	}
	if(iStart < totSize)
	  memcpy(output.getPtr() + iStart, input.getPtr(), (totSize - iStart) * output.element_size() );
      } else {
	if(fillZeros)
	  std::fill(output.getPtr() + sizeToCopy, output.getPtr() + totSize, 0);
      }
    }
  }
  return;
}

template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, bool fillZeros, bool recycle, int size1) {
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, false, false); // would need same scalar types before trying memcpy.
    if(!init) return;
    copyNimArrToNimArrInternal(output, input, size1, fillZeros, recycle);
  }
  return;
}

template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, bool fillZeros, bool recycle, int size1, int size2) {
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, size2, false, false); // would need same scalar types before trying memcpy.
    if(!init) return;
    int totsize = size1 * size2;
    copyNimArrToNimArrInternal(output, input, totsize, fillZeros, recycle);
  }
  return;
}

template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, bool fillZeros, bool recycle, int size1, int size2, int size3) {
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, size2, size3, false, false); // would need same scalar types before trying memcpy.
    if(!init) return;
    int totsize = size1 * size2 * size3;
    copyNimArrToNimArrInternal(output, input, totsize, fillZeros, recycle);
  }
  return;
}


template<typename NimArrOutput, typename NimArrInput>
  void assignNimArrToNimArr(NimArrOutput &output, NimArrInput &input, bool init, bool fillZeros, bool recycle, int size1, int size2, int size3, int size4) {
  if(output.isMap()) {
    _nimble_global_output << "Error from C++: using assignNimArrToNimArr with a map.\n"; nimble_print_to_R( _nimble_global_output);
  } else {
    output.setSize(size1, size2, size3, size4, false, false); // would need same scalar types before trying memcpy.
    if(!init) return;
    int totsize = size1 * size2 * size3 * size4;
    copyNimArrToNimArrInternal(output, input, totsize, fillZeros, recycle);
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
    nextBool = nimble_eigen_coeff_impl< Eigen::internal::traits<DerivedBool>::Flags & LinearAccessBit, int, DerivedBool, Eigen::Index >::getCoeff(BoolArg, i);
    if(nextBool) ans.push_back(i + 1); // That 1 makes it one-based indexing, which will then be adjusted back to zero-based when used for indexing something else
  }
  assignVectorToNimArr<NimArrOutput, std::vector<int> >(output, ans);
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
