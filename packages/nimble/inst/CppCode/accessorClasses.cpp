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

#include "nimble/accessorClasses.h"
#include "nimble/RcppUtils.h"
#include<sstream>
using std::istringstream;

NimArr<1, double> getParam_1D_double(int paramID, const NodeInstruction &useInfo, int iNodeFunction) {
  if(iNodeFunction == 0) paramID += 0;
  return(useInfo.nodeFunPtr->getParam_1D_double_block(paramID, useInfo.operand));
}

NimArr<2, double> getParam_2D_double(int paramID, const NodeInstruction &useInfo, int iNodeFunction) {
  if(iNodeFunction == 0) paramID += 0;
  return(useInfo.nodeFunPtr->getParam_2D_double_block(paramID, useInfo.operand));
}

// code for getBound copied over from getParam; none of this currently used as we only have scalar bounds for multivariate nodes
NimArr<1, double> getBound_1D_double(int boundID, const NodeInstruction &useInfo, int iNodeFunction) {
  if(iNodeFunction == 0) boundID += 0;
  return(useInfo.nodeFunPtr->getBound_1D_double_block(boundID, useInfo.operand));
}

NimArr<2, double> getBound_2D_double(int boundID, const NodeInstruction &useInfo, int iNodeFunction) {
  if(iNodeFunction == 0) boundID += 0;
  return(useInfo.nodeFunPtr->getBound_2D_double_block(boundID, useInfo.operand));
}

// see include/nimble/nimbleEigenNimArr.h for some templated versions of calculate, simulate, calculateDiff and getLogProb used for arbitrary index vectors
double calculate(NodeVectorClassNew &nodes) {
  double ans(0);
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  for(; iNode != iNodeEnd; iNode++)
    ans += iNode->nodeFunPtr->calculateBlock(iNode->operand);
  return(ans);
}

/* This is a function because it should be done once and not undone. */
void set_CppAD_atomic_info_for_model(NodeVectorClassNew_derivs &nodes,
				     std::vector<CppAD::local::atomic_index_info>* vec_ptr) {
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  if(instructions.size() > 0) {
    nodeFun *nodeFunInModelDLL;
    nodeFunInModelDLL = instructions[0].nodeFunPtr;
    nodeFunInModelDLL->set_atomic_info_from_nodeFun(vec_ptr);
  }
}

/* This is a class so it works by RAII */
set_CppAD_tape_info_for_model::set_CppAD_tape_info_for_model(NodeVectorClassNew_derivs &nodes,
				   CppAD::tape_id_t tape_id,
				   CppAD::local::ADTape<double>* tape_handle_) { // tape handle
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  not_empty  = instructions.size() > 0;
  if(not_empty) {
    nodeFunInModelDLL = instructions[0].nodeFunPtr;
    nodeFunInModelDLL->set_tape_ptr_from_nodeFun(tape_id, tape_handle_, false);
  }
}

set_CppAD_tape_info_for_model::~set_CppAD_tape_info_for_model() {
  if(not_empty) {
    nodeFunInModelDLL->set_tape_ptr_from_nodeFun(0, 0, true);
  }
}

CppAD::AD<double> calculate_ADproxyModel(NodeVectorClassNew_derivs &nodes,
					 bool includeExtraOutputStep,
					 bool recordingInfo__isRecording) {
  std::cout <<"entering calculate_ADproxyModel"<< std::endl;//"entering calculate_ADproxyModel" << "\n";
  std::cout<<"handle address: "<<CppAD::AD<double>::get_handle_address_nimble()<<std::endl;

  CppAD::AD<double> ans = 0;
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  std::cout<<"starting node calcs in calculate_ADproxyModel"<<std::endl;
  
  for(; iNode != iNodeEnd; iNode++)
    ans += iNode->nodeFunPtr->calculateBlock_ADproxyModel(iNode->operand);
  if(includeExtraOutputStep && recordingInfo__isRecording) {
    std::cout<<"starting extraOutputStep"<<std::endl;
    if(instructions.begin() != iNodeEnd) {
      std::cout<<"will do extraOutputStep"<<std::endl;
      // It is arbitrary to call this for the first node,
      // but it is important to have it done by a nodeFun
      // because that will be in the right compilation unit
      // to access the right globals (statics) from CppAD.
      instructions.begin()->nodeFunPtr->setup_extraOutput_step( nodes,
      								ans );
    }
    std::cout<<"done with extraOutputStep"<<std::endl;
  }

  return(ans);
}

// void assign_extraInputDummy(NodeVectorClassNew_derivs &nodes,
// 			    CppAD::AD<double> &extraInputDummy) {
//   nodes.extraInputDummy = extraInputDummy;
// }



// void setup_extraInput_step(NodeVectorClassNew_derivs &nodes) {
//   const vector<NodeInstruction> &instructions = nodes.getInstructions();
//   if(instructions.size() == 0) {
//       printf("No nodes for initialize_AD_model_before_recording\n");
//       return;
//   }
//   nodeFun* nodeFunInModelDLL = instructions[0].nodeFunPtr;
//   nodeFunInModelDLL->setup_extraInput_step(nodes);
// }

void initialize_AD_model_before_recording(NodeVectorClassNew_derivs &nodes) {
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  if(instructions.size() == 0) {
      printf("No nodes for initialize_AD_model_before_recording\n");
      return;
  }
  nodeFun* nodeFunInModelDLL = instructions[0].nodeFunPtr;
  nodeFunInModelDLL->initialize_AD_model_before_recording(nodes);
}

void init_dynamicVars(NodeVectorClassNew_derivs &NV,
		      std::vector<CppAD::AD<double> > &dynamicVars) {
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  NimArr<1, double> NimArrValues;

  if(length_extraInput > 0) {
    NimArrValues.setSize(length_extraInput);
    dynamicVars.resize(length_extraInput);
    getValues(NimArrValues, NV.model_extraInput_accessor);
    std::copy( NimArrValues.getPtr(),
	       NimArrValues.getPtr() + length_extraInput,
	       dynamicVars.begin());
  }
}

void copy_dynamicVars_to_model(NodeVectorClassNew_derivs &NV,
			       std::vector<CppAD::AD<double> > &dynamicVars) {
  int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
  NimArr<1, CppAD::AD<double> > NimArrValues_AD;

  if(length_extraInput > 0) {
    NimArrValues_AD.setSize(length_extraInput);
    std::copy( dynamicVars.begin(),
	       dynamicVars.end(),
	       NimArrValues_AD.getPtr());
    setValues_AD_AD(NimArrValues_AD, NV.model_AD_extraInput_accessor);
  }
}

// void update_dynamicVars(NodeVectorClassNew_derivs &NV,
// 			std::vector<double> &dynamicVars) {
//   int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
//   NimArr<1, double> NimArrValues;
//   if(length_extraInput > 0) {
//     NimArrValues.setSize(length_extraInput);
//     dynamicVars.resize(length_extraInput);
//     getValues(NimArrValues, NV.model_extraInput_accessor);
//     std::copy( NimArrValues.getPtr(),
// 	       NimArrValues.getPtr() + length_extraInput,
// 	       dynamicVars.begin());
//   }
// }

// void update_dynamicVars_meta(NodeVectorClassNew_derivs &NV,
// 			     std::vector< CppAD::AD<double> > &dynamicVars) {
//   int length_extraInput = NV.model_extraInput_accessor.getTotalLength();
//   NimArr<1, double> NimArrValues;
//   if(length_extraInput > 0) {
//     NimArrValues.setSize(length_extraInput);
//     dynamicVars.resize(length_extraInput);
//     getValues(NimArrValues, NV.model_extraInput_accessor);
//     std::copy( NimArrValues.getPtr(),
// 	       NimArrValues.getPtr() + length_extraInput,
// 	       dynamicVars.begin());
//   }
// }

double calculate(NodeVectorClassNew &nodes, int iNodeFunction) {
  if(nodes.getInstructions().size() < static_cast<unsigned int>(iNodeFunction)) {
    PRINTF("Warning in calculate: index of requested set of nodes is too large\n");
    return(0);
  }
  const NodeInstruction &oneUseInfo = nodes.getInstructions()[iNodeFunction-1];
  return(oneUseInfo.nodeFunPtr->calculateBlock(oneUseInfo.operand));
}

double calculateDiff(NodeVectorClassNew &nodes) {
  double ans(0);
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  for(; iNode != iNodeEnd; iNode++)
    ans += iNode->nodeFunPtr->calculateDiffBlock(iNode->operand);
  return(ans);
}

double calculateDiff(NodeVectorClassNew &nodes, int iNodeFunction) {
  if(nodes.getInstructions().size() < static_cast<unsigned int>(iNodeFunction)) {
    PRINTF("Warning in calculateDiff: index of requested set of nodes is too large\n");
    return(0);
  }
  const NodeInstruction &oneUseInfo = nodes.getInstructions()[iNodeFunction-1];
  return(oneUseInfo.nodeFunPtr->calculateDiffBlock(oneUseInfo.operand));
}

double getLogProb(NodeVectorClassNew &nodes) {
  double ans(0);
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  for(; iNode != iNodeEnd; iNode++)
    ans += iNode->nodeFunPtr->getLogProbBlock(iNode->operand);
  return(ans);
}

double getLogProb(NodeVectorClassNew &nodes, int iNodeFunction) {
  if(nodes.getInstructions().size() < static_cast<unsigned int>(iNodeFunction)) {
    PRINTF("Warning in getLogProb: index of requested set of nodes is too large\n");
    return(0);
  }
  const NodeInstruction &oneUseInfo = nodes.getInstructions()[iNodeFunction-1];
  return(oneUseInfo.nodeFunPtr->getLogProbBlock(oneUseInfo.operand));
}

void simulate(NodeVectorClassNew &nodes) {
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  for(; iNode != iNodeEnd; iNode++)
    iNode->nodeFunPtr->simulateBlock(iNode->operand);
}

void simulate(NodeVectorClassNew &nodes, int iNodeFunction) {
  if(nodes.getInstructions().size() < static_cast<unsigned int>(iNodeFunction)) {
    PRINTF("Warning in simulate: index of requested set of nodes is too large\n");
    return;
  }
  const NodeInstruction &oneUseInfo = nodes.getInstructions()[iNodeFunction-1];
  oneUseInfo.nodeFunPtr->simulateBlock(oneUseInfo.operand);
}

SingleVariableMapAccessBase::~SingleVariableMapAccessBase(){}

ManyVariablesMapAccessor::~ManyVariablesMapAccessor(){
    for(unsigned int i = 0; i < varAccessors.size(); ++i)
      delete static_cast<SingleVariableMapAccess*>(varAccessors[i]);
}

#ifdef __NIMBLE_DEBUG_ACCESSORS
void ManyVariablesMapAccessor::check() {
  if(varAccessors.size() == 0) PRINTF("Run-time error: using nimCopy on map accessor with length 0\n");
}

void ManyVariablesMapAccessor::check(int i) {
  if(varAccessors.size() == 0) PRINTF("Run-time error: using nimCopy on map accessor with length 0\n");
}
#endif

ManyModelValuesMapAccessor::~ManyModelValuesMapAccessor() {
    for(unsigned int i = 0; i < varAccessors.size(); ++i)
      delete static_cast<SingleModelValuesMapAccess*>(varAccessors[i]);
}

#ifdef __NIMBLE_DEBUG_ACCESSORS
void ManyModelValuesMapAccessor::check() {
  if(varAccessors.size() == 0) PRINTF("Run-time error: using nimCopy on map accessor with length 0\n");
}

void ManyModelValuesMapAccessor::check(int i) {
  if(varAccessors.size() == 0) PRINTF("Run-time error: using nimCopy on map accessor with length 0\n");
  if(i < 0 || i >= varAccessors.size()) PRINTF("Run-time error: using nimCopy on map accessor with an invalid row\n");
}
#endif

void ManyModelValuesMapAccessor::setRow(int i) {
  if(i != currentRow) {
    currentRow = i;
    vector<SingleVariableMapAccessBase *>::iterator iVar, iEnd;
    iEnd = varAccessors.end();
    for(iVar = varAccessors.begin(); iVar != iEnd; ++iVar) {
      static_cast<SingleModelValuesMapAccess*>(*iVar)->setRow(i);
    }
  }
//  return(varAccessors);
}

void nimArr_2_SingleModelAccess_AD_AD(SingleVariableMapAccessBase* SMVAPtr,
				      NimArrBase< CppAD::AD<double> > &nimArr,
				      int nimBegin,
				      int nimStride) {
  NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
  if(SMVAPtr->getSingleton()) {
    (*static_cast<NimArrBase< CppAD::AD<double> >* >(SMA_NimTypePtr))[SMVAPtr->offset] = (*nimArr.getVptr())[nimBegin];
  } else {
    dynamicMapCopyFlatToDim< CppAD::AD<double> , CppAD::AD<double> >(static_cast<NimArrBase< CppAD::AD<double> >* >(SMA_NimTypePtr),
								     SMVAPtr->getOffset(),
								     SMVAPtr->getStrides(),
								     SMVAPtr->getSizes(),
								     &nimArr,
								     nimBegin,
								     nimStride);
  }
}


void nimArr_2_ManyModelAccess_AD_AD(ManyVariablesMapAccessor &MMVAPtr, NimArrBase< CppAD::AD<double> > &nimArr){
  vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = SMVA_Vec->size();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = (*SMVA_Vec)[i];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd){
      PRINTF("Warning: in nimArr_2_ManyModelAccess, accessor larger than NimArr!\n");
      break;
    }
    nimArr_2_SingleModelAccess_AD_AD(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
    nimCurrent += nextNumVals;
    nimCurrentOffset += nextNumVals * nimArrStride;
  }
  if(nimCurrent != nimEnd)
    PRINTF("Warning: after completing nimArr_2_ManyModelAccess, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
}

template<class T>
void nimArr_2_ManyModelAccess(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr){
  vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = SMVA_Vec->size();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = (*SMVA_Vec)[i];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd){
      PRINTF("Warning: in nimArr_2_ManyModelAccess, accessor larger than NimArr!\n");
      break;
    }
    nimArr_2_SingleModelAccess<T>(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
    nimCurrent += nextNumVals;
    nimCurrentOffset += nextNumVals * nimArrStride;
  }
  if(nimCurrent != nimEnd)
    PRINTF("Warning: after completing nimArr_2_ManyModelAccess, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
}

template<class T>
void nimArr_2_ManyModelAccessIndex(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr, int index){
  vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  curSingleAccess = (*SMVA_Vec)[index];
  nextNumVals = (*curSingleAccess).getLength();
  if(nextNumVals + nimCurrent > nimEnd)
    PRINTF("Warning: in nimArr_2_ManyModelAccessIndex, accessor larger than NimArr!\n");
  nimArr_2_SingleModelAccess<T>(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
  nimCurrent += nextNumVals;
  if(nimCurrent != nimEnd)
    PRINTF("Warning: after completing nimArr_2_ManyModelAccessIndex, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
}
  
///////////////
// [accessors]_2_nimArr
// nimArr is "to". SMVAPtr is "from"
void SingleModelAccess_2_nimArr_AD_AD(SingleVariableMapAccessBase* SMVAPtr,
				      NimArrBase< CppAD::AD<double> > &nimArr,
				      int nimBegin,
				      int nimStride) {
  NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
  if(SMVAPtr->getSingleton()) {
    (*nimArr.getVptr())[nimBegin] = (*static_cast<NimArrBase< CppAD::AD<double> >* >(SMA_NimTypePtr))[SMVAPtr->offset];
  } else {
    dynamicMapCopyDimToFlat< CppAD::AD<double>, CppAD::AD<double> >(&nimArr,
								  nimBegin,
								  nimStride,
								  static_cast<NimArrBase< CppAD::AD<double> >* >(SMA_NimTypePtr),
								  SMVAPtr->getOffset(),
								  SMVAPtr->getStrides(),
								  SMVAPtr->getSizes() );
  }
}

void ManyModelAccess_2_nimArr_AD_AD(ManyVariablesMapAccessor &MMVAPtr, NimArrBase< CppAD::AD<double> > &nimArr){
  const vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = SMVA_Vec->size();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = (*SMVA_Vec)[i];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd){
      PRINTF("Warning: in nimArr_2_ManyModelAccess, accessor larger than NimArr!\n");
      break;
    }
    SingleModelAccess_2_nimArr_AD_AD(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
    nimCurrent += nextNumVals;
    nimCurrentOffset += nextNumVals * nimArrStride;
  }
  if(nimCurrent != nimEnd)
    PRINTF("Warning: after completing ManyModelAccess_2_nimArr, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
}

template<class T>
void ManyModelAccess_2_nimArr(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr){
  const vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = SMVA_Vec->size();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = (*SMVA_Vec)[i];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd){
      PRINTF("Warning: in nimArr_2_ManyModelAccess, accessor larger than NimArr!\n");
      break;
    }
    SingleModelAccess_2_nimArr<T>(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
    nimCurrent += nextNumVals;
    nimCurrentOffset += nextNumVals * nimArrStride;
  }
  if(nimCurrent != nimEnd)
    PRINTF("Warning: after completing ManyModelAccess_2_nimArr, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
}


template<class T>
void ManyModelAccessIndex_2_nimArr(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr, int index) {
  vector<SingleVariableMapAccessBase*> *SMVA_Vec = &(MMVAPtr.getMapAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  curSingleAccess = (*SMVA_Vec)[index];
  nextNumVals = (*curSingleAccess).getLength();
  if(nextNumVals + nimCurrent > nimEnd)
    PRINTF("Warning: in ManyModelAccessIndex_2_nimArr, accessor larger than NimArr!\n");
  SingleModelAccess_2_nimArr<T>(curSingleAccess, nimArr, nimCurrentOffset, nimArrStride);
}


//////////
void setValues(NimArrBase<double> &nimArr, ManyVariablesMapAccessor &MVA){
  nimArr_2_ManyModelAccess<double>(MVA, nimArr);
}
void setValues(NimArrBase<int> &nimArr, ManyVariablesMapAccessor &MVA){
  nimArr_2_ManyModelAccess<int>(MVA, nimArr);
}
void setValues_AD_AD(NimArrBase< CppAD::AD<double> > &nimArr,  ManyVariablesMapAccessor &MVA){
  nimArr_2_ManyModelAccess_AD_AD(MVA, nimArr);
}

void setValues_AD_AD_taping(NimArr<1, CppAD::AD<double> > &v,
			    ManyVariablesMapAccessor &MVA_AD,
			    ManyVariablesMapAccessor &MVA_orig,
			    bool recording){
  size_t totalLength = MVA_orig.getTotalLength();  
  if(!recording) {
    NimArr<1, double> dv;
    dv.setSize(totalLength);
    for(size_t ii = 0; ii < totalLength; ++ii) {
      dv[ii] = CppAD::Value(v[ii]);
    }
    setValues(dv, MVA_orig);
  } else {
    // Make a new extraOutputObject, which during tape execution will copy to the real (non-AD) model
    std::vector< CppAD::AD<double> > extraOutputDummyResult(1);
    std::vector< CppAD::AD<double> > extraOutputs(totalLength);
    for(size_t ii = 0; ii < totalLength; ++ii) {
      extraOutputs[ii] = v[ii];
    }
    std::cout << "remember to get the extraOutputObject destructed." << std::endl;
    atomic_extraOutputObject* localExtraOutputObject =
      new atomic_extraOutputObject("copying-extraOutputObject",
				   &MVA_orig); // These objects are stable members of nimbleFunction classes, so the pointer should be stable.
  // Operate the object so it is recorded in the tape
    (*localExtraOutputObject)(extraOutputs, extraOutputDummyResult);
    // It is unclear whether we then need to use the output in a way that forces CppAD to keep it as part of the calculation graph.
    //   The concern is that otherwise CppAD might optimize it away by determining that nothing really depends on it.
    //   The following line is essentially a no-operation for this purpose (extraOutputDummyResult[0] will always be 0 in value).
    v[0] += extraOutputDummyResult[0];
  }
  setValues_AD_AD(v, MVA_AD); // record copying on the tape
}


void setValues(NimArrBase<double> &nimArr, ManyVariablesMapAccessor &MVA, int index){
  nimArr_2_ManyModelAccessIndex<double>(MVA, nimArr, index-1);
}
void setValues(NimArrBase<int> &nimArr, ManyVariablesMapAccessor &MVA, int index){
  nimArr_2_ManyModelAccessIndex<int>(MVA, nimArr, index-1);
}
void getValues(NimArr<1, double> &nimArr, ManyVariablesMapAccessor &MVA){
  ManyModelAccess_2_nimArr<double>(MVA, nimArr);
}
void getValues(NimArr<1, int> &nimArr, ManyVariablesMapAccessor &MVA){
  ManyModelAccess_2_nimArr<int>(MVA, nimArr);
}
void getValues_AD_AD(NimArr<1, CppAD::AD<double> > &nimArr, ManyVariablesMapAccessor &MVA){
  ManyModelAccess_2_nimArr_AD_AD(MVA, nimArr);
}

void getValues(NimArr<1, double> &nimArr, ManyVariablesMapAccessor &MVA, int index){
  ManyModelAccessIndex_2_nimArr<double>(MVA, nimArr, index-1);
} 

void getValues(NimArr<1, int> &nimArr, ManyVariablesMapAccessor &MVA, int index){
  ManyModelAccessIndex_2_nimArr<int>(MVA, nimArr, index-1);
}



void nimCopy(const copierVectorClass &copiers) {
  vector<copierClass*>::const_iterator iCopy;
  for(iCopy = copiers.copyVector.begin(); iCopy != copiers.copyVector.end(); iCopy++) {
    (*iCopy)->copy(copiers.rowInfo);
  }
}

// A bit cheesy here: we add an unused argument only to trigger the right overloaded type in a simple way

void nimCopy(copierVectorClass &copiers, int rowFrom) {
  copiers.rowFrom() = rowFrom-1;
  nimCopy(copiers);
}

void nimCopy(copierVectorClass &copiers, int rowFrom, int rowTo) { // only ever called with rowFrom = 0 for right overloading
  copiers.rowTo() = rowTo-1;
  nimCopy(copiers);
}

void nimCopy(copierVectorClass &copiers, int rowFrom, int rowTo, int unused) { // only ever called when rowFrom and rowTo both needed
  copiers.rowFrom() = rowFrom-1;
  copiers.rowTo() = rowTo-1;
  nimCopy(copiers);
}


void copierVectorClass::setup(ManyVariablesMapAccessorBase *from, ManyVariablesMapAccessorBase *to, int isFromMV, int isToMV) {
  // Imitates old version of nimCopy but populates copyVector of correct choices of derived copierClass objects
  vector<SingleVariableMapAccessBase *> *fromAccessors = &(from->getMapAccessVector());
  vector<SingleVariableMapAccessBase *> *toAccessors = &(to->getMapAccessVector());

#ifdef __NIMBLE_DEBUG_ACCESSORS
  PRINTF("Entering nimCopy\n");
  from->check();
  to->check();
#endif

  if(fromAccessors->size() != toAccessors->size()) {
    _nimble_global_output<<"Error in setting up a copierVector: from and to access vectors have sizes "<<fromAccessors->size() << " and " << toAccessors->size() << "\n";
    nimble_print_to_R(_nimble_global_output);
  }
  copyVector.resize( fromAccessors->size() );
  vector<SingleVariableMapAccessBase *>::iterator iFrom, iTo, iFromEnd;
  iFromEnd = fromAccessors->end();
  iTo =  toAccessors->begin();
  int i = 0;
  for(iFrom = fromAccessors->begin(); iFrom != iFromEnd; iFrom++) {
    copyVector[i] = makeOneCopyClass(*iFrom, *iTo, isFromMV, isToMV); 
    iTo++;
    i++;
  }
}

copierVectorClass::copierVectorClass() {}

copierVectorClass::~copierVectorClass() {
  vector<copierClass*>::iterator iCopy;
  for(iCopy = copyVector.begin(); iCopy != copyVector.end(); iCopy++) {
    delete (*iCopy);
  }
}


void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to) { 

  vector<SingleVariableMapAccessBase *> *fromAccessors = &(from.getMapAccessVector());
  vector<SingleVariableMapAccessBase *> *toAccessors = &(to.getMapAccessVector());

#ifdef __NIMBLE_DEBUG_ACCESSORS
  PRINTF("Entering nimCopy\n");
  from->check();
  to->check();
#endif

  if(fromAccessors->size() != toAccessors->size()) {
    _nimble_global_output<<"Error in nimCopy: from and to access vectors have sizes "<<fromAccessors->size() << " and " << toAccessors->size() << "\n";
    nimble_print_to_R(_nimble_global_output);
  }

  vector<SingleVariableMapAccessBase *>::iterator iFrom, iTo, iFromEnd;
  iFrom = fromAccessors->begin();
  iFromEnd = fromAccessors->end();
  iTo =  toAccessors->begin();
  for( ; iFrom != iFromEnd; ++iFrom) {
    nimCopyOne(*iFrom, *iTo);
    ++iTo;
  }
}


void nimCopy(ManyVariablesMapAccessorBase &from, int rowFrom, ManyVariablesMapAccessorBase &to) {
#ifdef __NIMBLE_DEBUG_ACCESSORS
  PRINTF("Entering nimCopy with rowFrom\n");
  from.check(rowFrom-1);
#endif

  from.setRow(rowFrom - 1);
  nimCopy(from, to);
}

void nimCopy(ManyVariablesMapAccessorBase &from, int rowFrom, ManyVariablesMapAccessorBase &to, int rowTo) {
#ifdef __NIMBLE_DEBUG_ACCESSORS
  PRINTF("Entering nimCopy with rowFrom and rowTo\n");
  from.check(rowFrom-1);
  to.check(rowTo-1);
#endif
  to.setRow(rowTo - 1);
  from.setRow(rowFrom - 1);
  nimCopy(from, to);
}

void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to, int rowTo) {
#ifdef __NIMBLE_DEBUG_ACCESSORS
  PRINTF("Entering nimCopy with rowTo\n");
  to.check(rowTo-1);
#endif

  to.setRow(rowTo - 1);
  nimCopy(from, to);
}

copierClassBuilderCase< singletonCopierClass_M2M<double, double>, singletonCopierClass_M2M<double, int>, singletonCopierClass_M2M<int, int>, singletonCopierClass_M2M<int, double> > globalCopierBuilder_singleton_M2M;

copierClassBuilderCase< singletonCopierClass_M2MV<double, double>, singletonCopierClass_M2MV<double, int>, singletonCopierClass_M2MV<int, int>, singletonCopierClass_M2MV<int, double> > globalCopierBuilder_singleton_M2MV;

copierClassBuilderCase< singletonCopierClass_MV2M<double, double>, singletonCopierClass_MV2M<double, int>, singletonCopierClass_MV2M<int, int>, singletonCopierClass_MV2M<int, double> > globalCopierBuilder_singleton_MV2M;

copierClassBuilderCase< singletonCopierClass_MV2MV<double, double>, singletonCopierClass_MV2MV<double, int>, singletonCopierClass_MV2MV<int, int>, singletonCopierClass_MV2MV<int, double> > globalCopierBuilder_singleton_MV2MV;


copierClassBuilderCase< blockCopierClass<double, double, 1>, blockCopierClass<double, int, 1>, blockCopierClass<int, int, 1>, blockCopierClass<int, double, 1> > globalCopierClassBuilderBlock1;

copierClassBuilderCase< blockCopierClass<double, double, 2>, blockCopierClass<double, int, 2>, blockCopierClass<int, int, 2>, blockCopierClass<int, double, 2> > globalCopierClassBuilderBlock2;

copierClassBuilderCase< blockCopierClass<double, double, 3>, blockCopierClass<double, int, 3>, blockCopierClass<int, int, 3>, blockCopierClass<int, double, 3> > globalCopierClassBuilderBlock3;

copierClassBuilderCase< blockCopierClass<double, double, 4>, blockCopierClass<double, int, 4>, blockCopierClass<int, int, 4>, blockCopierClass<int, double, 4> > globalCopierClassBuilderBlock4;

copierClassBuilderCase< blockCopierClass<double, double, 5>, blockCopierClass<double, int, 5>, blockCopierClass<int, int, 5>, blockCopierClass<int, double, 5> > globalCopierClassBuilderBlock5;

copierClassBuilderCase< blockCopierClass<double, double, 6>, blockCopierClass<double, int, 6>, blockCopierClass<int, int, 6>, blockCopierClass<int, double, 6> > globalCopierClassBuilderBlock6;

copierClass* makeOneCopyClass(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int isFromMV, int isToMV) { // like nimCopyOne but it returns an appropriate derived copierClass object
  copierClassBuilderClass *copierClassBuilder;

  if(to->getSingleton()) {
#ifdef __NIMBLE_DEBUG_ACCESSORS
    if(!from->getSingleton()) PRINTF("Run-time error: to is a singleton but from is not a singleton\n");
    singletonCopyCheck(fromNimArr, from->getOffset());
    singletonCopyCheck(toNimArr, to->getOffset());
#endif
    switch(isFromMV) {
    case 0:
      switch(isToMV) {
      case 0:
	copierClassBuilder = &globalCopierBuilder_singleton_M2M;
	break;
      case 1:
	copierClassBuilder = &globalCopierBuilder_singleton_M2MV;
	break;
      default:
	NIMERROR("problem in makeOneCopyClass");
	return 0;
	break;
      };
      break;
    case 1:
      switch(isToMV) {
      case 0:
	copierClassBuilder = &globalCopierBuilder_singleton_MV2M;
	break;
      case 1:
	copierClassBuilder = &globalCopierBuilder_singleton_MV2MV;
	break;
      default:
	NIMERROR("problem in makeOneCopyClass");
	return 0;
	break;
      };
      break;
    default:
      NIMERROR("problem in makeOneCopyClass");
      return 0;
      break;
    }
    return copierClassBuilder->build(from, to, isFromMV, isToMV);
  }

  int mapDim = to->getStrides().size();
  switch(mapDim) {
  case 1:
    copierClassBuilder = &globalCopierClassBuilderBlock1;
    break;
  case 2:
    copierClassBuilder = &globalCopierClassBuilderBlock2;
    break;
  case 3:
    copierClassBuilder = &globalCopierClassBuilderBlock3;
    break;
  case 4:
    copierClassBuilder = &globalCopierClassBuilderBlock4;
    break;
  case 5:
    copierClassBuilder = &globalCopierClassBuilderBlock5;
    break;
  case 6:
    copierClassBuilder = &globalCopierClassBuilderBlock6;
    break;
  default:
    NIMERROR("problem in makeOneCopyClass");
    return 0;
    break;
  }
  return copierClassBuilder->build(from, to, isFromMV, isToMV);
}

void nimCopyOne(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to) { // map version
  nimType fromType, toType;
  NimArrType *fromNimArr, *toNimArr;
  fromNimArr = from->getNimArrPtr();
  toNimArr = to->getNimArrPtr();
  fromType = fromNimArr->getNimType();
  toType = toNimArr->getNimType();
  if(to->getSingleton()) {
#ifdef __NIMBLE_DEBUG_ACCESSORS
    if(!from->getSingleton()) PRINTF("Run-time error: to is a singleton but from is not a singleton\n");
    singletonCopyCheck(fromNimArr, from->getOffset());
    singletonCopyCheck(toNimArr, to->getOffset());
#endif
    switch(fromType) {
    case nimType::DOUBLE:
      switch(toType) {
      case nimType::DOUBLE:
	(*static_cast<NimArrBase<double> *>(toNimArr))[to->getOffset()] = (*static_cast<NimArrBase<double> *>(fromNimArr))[from->getOffset()];
	break;
    case nimType::INT:
	(*static_cast<NimArrBase<int> *>(toNimArr))[to->getOffset()] = (*static_cast<NimArrBase<double> *>(fromNimArr))[from->getOffset()];
      break;
    default:
      _nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
      nimble_print_to_R(_nimble_global_output);
    }
    break;
  case nimType::INT:
    switch(toType) {
    case nimType::DOUBLE:
      (*static_cast<NimArrBase<double> *>(toNimArr))[to->getOffset()] = (*static_cast<NimArrBase<int> *>(fromNimArr))[from->getOffset()];
      break;
    case nimType::INT:
	(*static_cast<NimArrBase<int> *>(toNimArr))[to->getOffset()] = (*static_cast<NimArrBase<int> *>(fromNimArr))[from->getOffset()];
      break;
    default:
      _nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
      nimble_print_to_R(_nimble_global_output);
    }
    break;
  default:
    _nimble_global_output<<"Error in nimCopyOne: unknown type for source\n";
    nimble_print_to_R(_nimble_global_output);
    }
  } else {
#ifdef __NIMBLE_DEBUG_ACCESSORS
    dynamicMapCopyCheck(toNimArr, to->getOffset(), to->getStrides(), to->getSizes());
    dynamicMapCopyCheck(fromNimArr, from->getOffset(), from->getStrides(), from->getSizes());
#endif
    switch(fromType) {
    case nimType::DOUBLE:
      switch(toType) {
      case nimType::DOUBLE:
	dynamicMapCopy<double, double>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      case nimType::INT:
	dynamicMapCopy<double, int>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      default:
	_nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
	nimble_print_to_R(_nimble_global_output);
      }
      break;
    case nimType::INT:
      switch(toType) {
      case nimType::DOUBLE:
	dynamicMapCopy<int, double>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      case nimType::INT:
	dynamicMapCopy<int, int>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      default:
	_nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
      }
      break;
    default:
      _nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
    }
  }
}

void singletonCopyCheck(NimArrType *NAT, int offset) {
  nimType NATtype = NAT->getNimType();
  int NATsize;
  switch(NATtype) {
  case nimType::INT:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->size(); //getVptr()->size();
    break;
  case nimType::DOUBLE:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->size(); //getVptr()->size();
    break;
  default:
    PRINTF("Error with a NimArrType type\n");
    return;
    break;
  }
  if(offset < 0 || offset >= NATsize) PRINTF("Run-time error: bad singleton offset\n");
}

void dynamicMapCopyCheck(NimArrType *NAT, int offset, vector<int> &strides, vector<int> &sizes) {
  nimType NATtype = NAT->getNimType();
  int NATsize;
  switch(NATtype) {
  case nimType::INT:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->size(); //getVptr()->size();
    break;
  case nimType::DOUBLE:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->size(); //getVptr()->size();
    break;
  default:
    PRINTF("Error with a NimArrType type\n");
    return;
    break;
  }
  if(offset < 0 || offset >= NATsize) PRINTF("Run-time error: bad offset\n");
  int lastOffset = offset;
  for(unsigned int i = 0; i < strides.size(); ++i) {
    lastOffset += sizes[i] * strides[i];
  }
  if(lastOffset < 0 || lastOffset >= NATsize) PRINTF("Run-time error: bad lastOffset\n");
}

SEXP getListElement(SEXP list, const char *str){
	SEXP ans = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);
	PROTECT(ans);
	PROTECT(names);
	for(int i = 0; i < LENGTH(list); i++){
		if(strcmp(CHAR(STRING_ELT(names, i) ), str) == 0){
			ans = VECTOR_ELT(list, i);
			break;
		}
	}
	UNPROTECT(2);
	return(ans);
}

void populateNodeFxnVectorNew_internal_forDerivs(NodeVectorClassNew_derivs* nfv, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS, SEXP SderivInfo) {
  (*nfv).populateDerivsInfo(SderivInfo);
  int len = LENGTH(S_ROWINDS);
  if(len == 0) return;
  int* gids = INTEGER(S_GIDs);
  int* rowinds = INTEGER(S_ROWINDS);
  int index;
  NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumberedObj));
  int nextRowInd;
  (*nfv).instructions.clear();
  for(int i = 0; i < len; i++){
    index = gids[i] - 1;
    nextRowInd = rowinds[i]-1;
    if(nextRowInd == -1) { // should only happen from a scalar, so there is one dummy indexedNodeInfo
      nextRowInd = 0;
    }
    (*nfv).instructions.push_back(NodeInstruction(static_cast<nodeFun*>(numObj->getObjectPtr(index)), nextRowInd));
  }
}

void populateNodeFxnVectorNew_copyFromRobject_forDerivs(void *nodeFxnVec_to, SEXP S_nodeFxnVec_from ) {
   SEXP S_indexingInfo;
   SEXP S_pxData;
   PROTECT(S_pxData = Rf_allocVector(STRSXP, 1));
   SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData"));
  PROTECT(S_indexingInfo = VECTOR_ELT(S_nodeFxnVec_from, 1));
  SEXP S_declIDs;
  PROTECT(S_declIDs = VECTOR_ELT(S_indexingInfo, 0));
  SEXP S_rowIndices;
  PROTECT(S_rowIndices = VECTOR_ELT(S_indexingInfo, 1));
  SEXP S_numberedPtrs;
  PROTECT(S_numberedPtrs = Rf_findVarInFrame(PROTECT(GET_SLOT(
							      Rf_findVarInFrame(PROTECT(GET_SLOT(
												 Rf_findVarInFrame(PROTECT(GET_SLOT(
																    VECTOR_ELT(S_nodeFxnVec_from,
																	       2
																	       ),
																    S_pxData)),
														   Rf_install("CobjectInterface")
														   ),
												 S_pxData)),
										Rf_install(".nodeFxnPointers_byDeclID")),
							      S_pxData)),
					     Rf_install(".ptr")
					     )
  );
  SEXP SderivInfo;
  PROTECT(SderivInfo = VECTOR_ELT(S_nodeFxnVec_from, 3));
  NodeVectorClassNew_derivs* nfv_derivs = static_cast<NodeVectorClassNew_derivs*>(nodeFxnVec_to);
  populateNodeFxnVectorNew_internal_forDerivs(nfv_derivs, S_declIDs, S_numberedPtrs, S_rowIndices, SderivInfo);
  UNPROTECT(9);
}

SEXP populateNodeFxnVectorNew_byDeclID_forDerivs(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS, SEXP SderivInfo){
  NodeVectorClassNew_derivs* nfv_derivs = static_cast<NodeVectorClassNew_derivs*>(R_ExternalPtrAddr(SnodeFxnVec)) ;
  populateNodeFxnVectorNew_internal_forDerivs(nfv_derivs, S_GIDs, SnumberedObj, S_ROWINDS, SderivInfo);
  return(R_NilValue);
}


void populateNodeFxnVectorNew_internal(NodeVectorClassNew* nfv, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS ) {
  int len = LENGTH(S_ROWINDS);
  if(len == 0) return;
  int* gids = INTEGER(S_GIDs);
  int* rowinds = INTEGER(S_ROWINDS);
  int index;
  NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumberedObj));
  int nextRowInd;
  (*nfv).instructions.clear();
  for(int i = 0; i < len; i++){
    index = gids[i] - 1;
    nextRowInd = rowinds[i]-1;
    if(nextRowInd == -1) { // should only happen from a scalar, so there is one dummy indexedNodeInfo
      nextRowInd = 0;
    }
    (*nfv).instructions.push_back(NodeInstruction(static_cast<nodeFun*>(numObj->getObjectPtr(index)), nextRowInd));
  }
}

void populateNodeFxnVectorNew_copyFromRobject(void *nodeFxnVec_to, SEXP S_nodeFxnVec_from ) {
  SEXP S_indexingInfo;
   SEXP S_pxData;
   S_pxData = PROTECT(Rf_allocVector(STRSXP, 1));
   SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData"));
   S_indexingInfo = PROTECT(VECTOR_ELT(S_nodeFxnVec_from, 1));
  SEXP S_declIDs;
  S_declIDs = PROTECT(VECTOR_ELT(S_indexingInfo, 0));
  SEXP S_rowIndices;
  S_rowIndices = PROTECT(VECTOR_ELT(S_indexingInfo, 1));
  SEXP S_numberedPtrs;
 // equivalent to S_nodeFxnVec_from[["model"]]$CobjectInterface$.nodeFxnPointers_byDeclID$.ptr
  // implemented by S_nodeFxnVec_from[[2]]@.xData[["CobjectInterface"]]@.xData[[".nodeFxnPointers_byDeclID"]]@.xData[[".ptr"]]
  S_numberedPtrs = PROTECT(Rf_findVarInFrame(PROTECT(GET_SLOT(
							     PROTECT(Rf_findVarInFrame(PROTECT(GET_SLOT(
													PROTECT(Rf_findVarInFrame(PROTECT(GET_SLOT(
																    VECTOR_ELT(S_nodeFxnVec_from,
																	       2
																	       ),
																    S_pxData)),
														   Rf_install("CobjectInterface")
														   )),
												 S_pxData)),
										       Rf_install(".nodeFxnPointers_byDeclID"))),
							      S_pxData)),
					     Rf_install(".ptr")
					     )
	  );
  NodeVectorClassNew* nfv = static_cast<NodeVectorClassNew*>(nodeFxnVec_to);
  populateNodeFxnVectorNew_internal(nfv, S_declIDs, S_numberedPtrs, S_rowIndices);
  UNPROTECT(10);
}

SEXP populateNodeFxnVectorNew_byDeclID(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS){
  NodeVectorClassNew* nfv = static_cast<NodeVectorClassNew*>(R_ExternalPtrAddr(SnodeFxnVec) ) ;
  populateNodeFxnVectorNew_internal(nfv, S_GIDs, SnumberedObj, S_ROWINDS);
  return(R_NilValue);
}

SEXP populateIndexedNodeInfoTable(SEXP StablePtr, SEXP StableContents) {
  SEXP Sdim;
  PROTECT(Sdim = Rf_getAttrib(StableContents, R_DimSymbol));
  if(LENGTH(Sdim) != 2) {
    PRINTF("Warning from populateIndexedNodeInfoTable: LENGTH(Sdim) != 2");
    UNPROTECT(1);
    return R_NilValue;
  }
  int nrow = INTEGER(Sdim)[0];
  int ncol = INTEGER(Sdim)[1];
  vector<indexedNodeInfo> *tablePtr = static_cast<vector<indexedNodeInfo> *>(R_ExternalPtrAddr(StablePtr));
  if(nrow == 0) {
    void *vptr=0;
    tablePtr->push_back(indexedNodeInfo(static_cast<int *>(vptr), 0, 0));
    if(ncol != 0) {
      PRINTF("Warning from populateIndexedNodeInfoTable: nrow == 0 but ncol != 0.");
    }
  } else {
    if(!Rf_isNumeric(StableContents)) {
      PRINTF("Warning from populateIndexedNodeInfoTable: StableContents is not numeric");
      UNPROTECT(1);
      return R_NilValue;
    }
    if(Rf_isInteger(StableContents)) {
      int *contentsPtr = INTEGER(StableContents);
      tablePtr->reserve(nrow);
      for(int i = 0; i < nrow; i++) {
	tablePtr->push_back(indexedNodeInfo(contentsPtr + i, ncol, nrow));
      }
    } else {
      double *contentsPtrd = REAL(StableContents);
      tablePtr->reserve(nrow);
      for(int i = 0; i < nrow; i++) {
	tablePtr->push_back(indexedNodeInfo(contentsPtrd + i, ncol, nrow));
      }
    }
  }
  UNPROTECT(1);
  return R_NilValue;
}

SEXP populateCopierVector(SEXP ScopierVector, SEXP SfromPtr, SEXP StoPtr, SEXP SintIsFromMV, SEXP SintIsToMV) {
  copierVectorClass *copierVector = static_cast<copierVectorClass*>(R_ExternalPtrAddr(ScopierVector));
  ManyVariablesMapAccessorBase* fromValuesAccessor = static_cast<ManyVariablesMapAccessorBase*>(R_ExternalPtrAddr(SfromPtr) );
  ManyVariablesMapAccessorBase* toValuesAccessor = static_cast<ManyVariablesMapAccessorBase*>(R_ExternalPtrAddr(StoPtr) );
  int isFromMV = INTEGER(SintIsFromMV)[0];
  int isToMV   = INTEGER(SintIsToMV)[0];
  copierVector->setup(fromValuesAccessor, toValuesAccessor, isFromMV, isToMV);
  return(R_NilValue);
}

class mapInfoClass {
public:
  int offset;
  vector<int> sizes;
  vector<int> strides;
};

SEXP varAndIndices2Rlist(const varAndIndicesClass &input) {
  SEXP Soutput, Sindices;
  PROTECT(Soutput = Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Soutput, 0, string_2_STRSEXP(input.varName));
  int numinds = input.indices.size();
  PROTECT(Sindices = Rf_allocVector(VECSXP, numinds));
  for(int i = 0; i < numinds; i++) {
    SET_VECTOR_ELT(Sindices, i, vectorInt_2_SEXP(input.indices[i]));
  }
  SET_VECTOR_ELT(Soutput, 1, Sindices);

  vector<string> newNames(2);
  newNames[0].assign("varName");
  newNames[1].assign("indices");
  SEXP SnewNames;
  PROTECT(SnewNames = vectorString_2_STRSEXP(newNames));
  Rf_setAttrib(Soutput, R_NamesSymbol, SnewNames);

  UNPROTECT(3);
  return(Soutput);
}

SEXP getVarAndIndices(SEXP Sstring) {
  string input(STRSEXP_2_string(Sstring, 0));
  varAndIndicesClass output;
  parseVarAndInds(input, output);
  return(varAndIndices2Rlist(output));
}

void varAndIndices2mapParts(const varAndIndicesClass &varAndInds, int snDim, const vector<int> &sizes, mapInfoClass &output) {
  unsigned int nDim = static_cast<unsigned int>(snDim);
  output.sizes.resize(0);
  output.strides.resize(0);
  //  bool sizeOne(sizes.size() == 0);
  int Rindexing(1); // assume indexing comes in R form (Starting at 1).  output does not depend on indexing.
  int offset = 0;
  int currentStride = 1;
  if((nDim > 0) & (varAndInds.indices.size() == 0)) {
    if(sizes.size() == 0) output.sizes.push_back(1); else output.sizes = sizes;
    output.strides.push_back(1);
    if(nDim > 1) {
      for(unsigned int i = 1; i < nDim; i++) output.strides.push_back( output.strides[i-1] * output.sizes[i-1] );
    }
  } else {
    //    vector<bool> blockBool(nDim, false);
    int thisSize;
    if(nDim != sizes.size()) {
      _nimble_global_output<<"Confused in varAndInds2MapParts: nDim != sizes.size()\n";
      nimble_print_to_R(_nimble_global_output);
    }
    if(nDim != varAndInds.indices.size()) {
      _nimble_global_output<<"Confused in varAndInds2MapParts: nDim != varAndInds.indices.size()\n";
      nimble_print_to_R(_nimble_global_output);
    }
    for(unsigned int i = 0; i < nDim; ++i) {
      thisSize = varAndInds.indices[i].size();
      switch(thisSize) {
      case 0: // index is blank
	//	blockBool[i] = true;
	output.sizes.push_back( sizes[i] );
	output.strides.push_back( currentStride );
	break;
      case 1: // index is single number
	offset += (varAndInds.indices[i][0] - Rindexing) * currentStride;
	break;
      case 2:  // index is a block range
	offset += (varAndInds.indices[i][0] - Rindexing) * currentStride;
	output.sizes.push_back(varAndInds.indices[i][1] - varAndInds.indices[i][0] + 1);
	output.strides.push_back( currentStride );
	break;
      default:
	_nimble_global_output<<"Confused in varAndInds2MapParts: an index content has length > 2\n";
	nimble_print_to_R(_nimble_global_output);
	break;
      }
      currentStride *= sizes[i];
    }
  }
  output.offset = offset;
}

SEXP mapInfo2Rlist(const mapInfoClass &input) {
  SEXP Soutput;
  PROTECT(Soutput = Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(Soutput, 0, int_2_SEXP(input.offset));
  SET_VECTOR_ELT(Soutput, 1, vectorInt_2_SEXP(input.sizes));
  SET_VECTOR_ELT(Soutput, 2, vectorInt_2_SEXP(input.strides));
  vector<string> newNames(3);
  newNames[0].assign("offset");
  newNames[1].assign("sizes");
  newNames[2].assign("strides");
  SEXP SnewNames;
  PROTECT(SnewNames = vectorString_2_STRSEXP(newNames));
  Rf_setAttrib(Soutput, R_NamesSymbol, SnewNames);
  UNPROTECT(2);
  return(Soutput);
}

SEXP varAndIndices2mapParts(SEXP SvarAndIndicesExtPtr, SEXP Ssizes, SEXP SnDim) {
  varAndIndicesClass *varAndIndicesPtr = static_cast<varAndIndicesClass *>(R_ExternalPtrAddr(SvarAndIndicesExtPtr));
  vector<int> sizes(SEXP_2_vectorInt(Ssizes, 0));
  int nDim(SEXP_2_int(SnDim, 0, 0));
  mapInfoClass output;
  varAndIndices2mapParts(*varAndIndicesPtr, nDim, sizes, output);
  return(mapInfo2Rlist(output));
}

SEXP var2mapParts(SEXP Sinput, SEXP Ssizes, SEXP SnDim) {
  string input(STRSEXP_2_string(Sinput, 0));
  varAndIndicesClass varAndIndices;
  parseVarAndInds(input, varAndIndices);
  vector<int> sizes(SEXP_2_vectorInt(Ssizes, 0));
  int nDim(SEXP_2_int(SnDim, 0, 0));
  mapInfoClass output;
  varAndIndices2mapParts(varAndIndices, nDim, sizes, output);
  return(mapInfo2Rlist(output));
}

//#define _DEBUG_POPULATE_MAP_ACCESSORS


void populateValueMapAccessorsFromNodeNames_internal(ManyVariablesMapAccessorBase* valuesAccessor,
						     SEXP SnodeNames,
						     SEXP SsizesAndNdims,
						     SEXP SModelOrModelValuesPtr) {
  vector<string> nodeNames;
  STRSEXP_2_vectorString(SnodeNames, nodeNames);
  NamedObjects *sourceNamedObject = static_cast<NamedObjects*>(R_ExternalPtrAddr(SModelOrModelValuesPtr));
  int numNames = nodeNames.size();
  valuesAccessor->resize(numNames);
  vector<SingleVariableMapAccessBase *> *singleAccessors = &(valuesAccessor->getMapAccessVector());
  varAndIndicesClass varAndIndices;
  int nDim;
  vector<int> sizes;
  SEXP SoneSizesAndNdims;
  mapInfoClass mapInfo;
  int totalLength = 0;
  
#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
  _nimble_global_output<<"New: "<<numNames<<"\n";
  nimble_print_to_R(_nimble_global_output);
#endif

  for(int i = 0; i < numNames; i++) {
    PROTECT(SoneSizesAndNdims = VECTOR_ELT(SsizesAndNdims, i));
    sizes = SEXP_2_vectorInt(VECTOR_ELT(SoneSizesAndNdims, 0));
    nDim = SEXP_2_int(VECTOR_ELT(SoneSizesAndNdims, 1));
#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
    _nimble_global_output<<nodeNames[i]<<"\n";
    nimble_print_to_R(_nimble_global_output);
#endif
    parseVarAndInds(nodeNames[i], varAndIndices);
#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
    _nimble_global_output<<varAndIndices.varName<<" parsed: (";
    for(int j = 0; j < varAndIndices.indices.size(); j++) {
      _nimble_global_output<<"(";
      for(int k = 0; k < varAndIndices.indices[j].size(); k++) _nimble_global_output<<varAndIndices.indices[j][k]<<" ";
      _nimble_global_output<<") ";
    }
    _nimble_global_output<<"\n";
    _nimble_global_output<<"input nDim "<< nDim << "(";
    for(int j = 0; j < sizes.size(); j++) _nimble_global_output<<sizes[j]<<" ";
    _nimble_global_output<<")\n";
    nimble_print_to_R(_nimble_global_output);
#endif
    varAndIndices2mapParts(varAndIndices, nDim, sizes, mapInfo);
    (*singleAccessors)[i]->getOffset() = mapInfo.offset;
    (*singleAccessors)[i]->getSizes() = mapInfo.sizes;
    (*singleAccessors)[i]->calculateLength();
    totalLength += (*singleAccessors)[i]->getLength();
    (*singleAccessors)[i]->getStrides() = mapInfo.strides;
    (*singleAccessors)[i]->getSingleton() = (*singleAccessors)[i]->getLength() == 1;
    (*singleAccessors)[i]->setObject( sourceNamedObject->getObjectPtr(varAndIndices.varName) );

#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
    _nimble_global_output<< varAndIndices.varName <<" "<<mapInfo.offset<<", (";
    for(int j = 0; j < mapInfo.sizes.size(); j++) _nimble_global_output<<mapInfo.sizes[j]<<",";
    _nimble_global_output<<"), "<<(*singleAccessors)[i]->getLength()<<", (";
    for(int j = 0; j < mapInfo.strides.size(); j++) _nimble_global_output<<mapInfo.strides[j]<<",";
    _nimble_global_output<<"), "<<(*singleAccessors)[i]->getSingleton()<<".\n";
    nimble_print_to_R(_nimble_global_output);
#endif


    UNPROTECT(1);
  }
  valuesAccessor->getTotalLength() = totalLength;
#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
  _nimble_global_output<<totalLength<<"\n";
  nimble_print_to_R(_nimble_global_output);
#endif
}

void populateValueMapAccessorsFromNodeNames_copyFromRobject(void *VvaluesAccessor,
							    SEXP Sargs) {
  //  _nimble_global_output<<"In new copy system\n";
  //  nimble_print_to_R(_nimble_global_output);
  ManyVariablesMapAccessorBase* valuesAccessor = static_cast<ManyVariablesMapAccessorBase*>(VvaluesAccessor);
  SEXP SnodeNames, SsizesAndNdims, SModelOrModelValuesPtr;
  PROTECT(SnodeNames = VECTOR_ELT(Sargs, 0));
  PROTECT(SsizesAndNdims = VECTOR_ELT(Sargs, 1));
  PROTECT(SModelOrModelValuesPtr = VECTOR_ELT(Sargs, 2));
  populateValueMapAccessorsFromNodeNames_internal(valuesAccessor,
						  SnodeNames,
						  SsizesAndNdims,
						  SModelOrModelValuesPtr);
  UNPROTECT(3);
}

// call both with a bunch of output generated...
SEXP populateValueMapAccessorsFromNodeNames(SEXP StargetPtr, SEXP SnodeNames, SEXP SsizesAndNdims, SEXP SModelOrModelValuesPtr ) {
  ManyVariablesMapAccessorBase* valuesAccessor = static_cast<ManyVariablesMapAccessorBase*>(R_ExternalPtrAddr(StargetPtr) );
  populateValueMapAccessorsFromNodeNames_internal(valuesAccessor, SnodeNames, SsizesAndNdims, SModelOrModelValuesPtr);
  return(R_NilValue);
}

SEXP populateValueMapAccessors(SEXP StargetPtr, SEXP SsourceList, SEXP SModelOrModelValuesPtr ) {
  // typeCode = 1 for model, 2 for modelValues
  ManyVariablesMapAccessorBase* valuesAccessor = static_cast<ManyVariablesMapAccessorBase*>(R_ExternalPtrAddr(StargetPtr) );
  int numAccessors = LENGTH(SsourceList);
  valuesAccessor->resize(numAccessors);
  vector<SingleVariableMapAccessBase *> *singleAccessors = &(valuesAccessor->getMapAccessVector());

  NamedObjects *sourceNamedObject = static_cast<NamedObjects*>(R_ExternalPtrAddr(SModelOrModelValuesPtr));
  SEXP SoneSource;
  string varName;
  int totalLength = 0;

#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
  _nimble_global_output<<"Old: "<<numAccessors<<"\n";
  nimble_print_to_R(_nimble_global_output);
#endif


  for(int i = 0; i < numAccessors; ++i) {
    PROTECT(SoneSource = VECTOR_ELT(SsourceList, i));
    (*singleAccessors)[i]->getOffset() = SEXP_2_int(VECTOR_ELT(SoneSource, 0));
    (*singleAccessors)[i]->getSizes() = SEXP_2_vectorInt(VECTOR_ELT(SoneSource, 1));
    (*singleAccessors)[i]->calculateLength();
    totalLength += (*singleAccessors)[i]->getLength();
    (*singleAccessors)[i]->getStrides() = SEXP_2_vectorInt(VECTOR_ELT(SoneSource, 2));
    (*singleAccessors)[i]->getSingleton() = static_cast<bool>(SEXP_2_int(VECTOR_ELT(SoneSource, 4)));
    varName = STRSEXP_2_string(VECTOR_ELT(SoneSource, 3), 0);
    (*singleAccessors)[i]->setObject( sourceNamedObject->getObjectPtr(varName) );
#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
    _nimble_global_output<< varName <<" "<<(*singleAccessors)[i]->getOffset()<<", (";
    for(int j = 0; j < (*singleAccessors)[i]->getSizes().size(); j++) _nimble_global_output<<(*singleAccessors)[i]->getSizes()[j]<<",";
    _nimble_global_output<<"), "<<(*singleAccessors)[i]->getLength()<<", (";
    for(int j = 0; j < (*singleAccessors)[i]->getStrides().size(); j++) _nimble_global_output<<(*singleAccessors)[i]->getStrides()[j]<<",";
    _nimble_global_output<<"), "<<(*singleAccessors)[i]->getSingleton()<<".\n";
    nimble_print_to_R(_nimble_global_output);
#endif

    UNPROTECT(1);
  }
  valuesAccessor->getTotalLength() = totalLength;
#ifdef _DEBUG_POPULATE_MAP_ACCESSORS
  _nimble_global_output<<totalLength<<"\n";
  nimble_print_to_R(_nimble_global_output);
#endif

  return(R_NilValue);
}

///

template<class T>
void cRemoveAccessor(T* aPtr, int index, bool removeAll){
	if(removeAll == TRUE){
		(*aPtr).varAccessors.erase( (*aPtr).varAccessors.begin(), (*aPtr).varAccessors.end() );
		return;
	}
	if(index < 0 || index >= static_cast<signed int>((*aPtr).varAccessors.size()) ){
		PRINTF("Warning: attempting to delete invalid index value of ManyAccessor\n");
		return;
	}
	(*aPtr).varAccessors.erase( (*aPtr).varAccessors.begin() + index);
	return;
}

template<class Many, class Single>
void cAddAccessor(Many* mPtr, Single* sPtr, bool addAtEnd, int index){
	int size = (*mPtr).varAccessors.size();
	if(addAtEnd == TRUE ) {
		(*mPtr).varAccessors.push_back(sPtr);
		return;
	}
	if((index >= size) | (index < 0)){
		PRINTF("Invalid index passed to addAccessor\n");
		return;
	}
	(*mPtr).varAccessors[index] = sPtr;
	return;
}

SEXP manualSetNRows(SEXP Sextptr, SEXP nRows){
  	Values* vPtr = static_cast<Values*>(R_ExternalPtrAddr(Sextptr) );
  	(*vPtr).numRows = INTEGER(nRows)[0];
  return(R_NilValue);
  }

indexedNodeInfo generateDummyIndexedNodeInfo(){
	vector<double> dummyVec;
	dummyVec.push_back(0);
	indexedNodeInfo dummyInfo = dummyVec;
	return(dummyInfo);
};	
