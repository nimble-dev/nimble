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


// 2. We will only delete the singleVariableAccessors when we delete the ManyVariablesAccesor that they are contained in
		// Cliff's note: actually, we have changed the paradigm. We now make one variable accessor per node, which can appear in many different
		// ManyVariableAccessors. All unique singleAccessors will be held in a single SpecializedNumberedObjects, which takes care of the finalizer.
ManyVariablesAccessor::~ManyVariablesAccessor(){
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

// 3.
void ManyModelValuesAccessor::setRow(int i) {
  if(i != currentRow) {
    currentRow = i;
    vector<SingleVariableAccessBase *>::iterator iVar, iEnd;
    iEnd = varAccessors.end();
    for(iVar = varAccessors.begin(); iVar != iEnd; ++iVar) {
      static_cast<SingleModelValuesAccess*>(*iVar)->setRow(i);
    }
  }
//  return(varAccessors);
}

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

/////////////////////
// new copying functions and getValues and setValues
// The nimArr may be a map, but it must be 1D (whether a map or not)

/////////
// nimArr_2_[accessors]
// nimArr is "from".  SMVAPtr is "to"

// moved to accessorClasses.h so it can be seen from nimbleEigenNimArr.h
// template<class T>
// void nimArr_2_SingleModelAccess(SingleVariableMapAccessBase* SMVAPtr, NimArrBase<T> &nimArr, int nimBegin, int nimStride){


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

// Moved to accessorClasses.h so nimbleEigenNimArr.h can see it
//template<class T>
//void SingleModelAccess_2_nimArr(SingleVariableMapAccessBase* SMVAPtr, NimArrBase<T> &nimArr, int nimBegin, int nimStride){


template<class T>
void ManyModelAccess_2_nimArr(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr){
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
  //SingleModelAccess_2_nimArr<T>(curSingleAccess, nimArr, 0, 0 ); 
}


//////////
//
void setValues(NimArrBase<double> &nimArr, ManyVariablesMapAccessor &MVA){
	nimArr_2_ManyModelAccess<double>(MVA, nimArr);
}

void setValues(NimArrBase<int> &nimArr, ManyVariablesMapAccessor &MVA){
	nimArr_2_ManyModelAccess<int>(MVA, nimArr);
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

void getValues(NimArr<1, double> &nimArr, ManyVariablesMapAccessor &MVA, int index){
  ManyModelAccessIndex_2_nimArr<double>(MVA, nimArr, index-1);
  } 

void getValues(NimArr<1, int> &nimArr, ManyVariablesMapAccessor &MVA, int index){
  ManyModelAccessIndex_2_nimArr<int>(MVA, nimArr, index-1);
}

////////////////////////
// Old versions
// Functions for copying

//template<int D, class T>
//void nimArr_2_SingleModelAccess(SingleVariableAccessBase* SMVAPtr, NimArr<D, T> &nimArr, int nimBegin){
template<class T>
void nimArr_2_SingleModelAccess(SingleVariableAccessBase* SMVAPtr, NimArrBase<T> &nimArr, int nimBegin){
	NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
	nimType SMA_Type = (*SMA_NimTypePtr).getNimType();
	int SMA_length = (*SMVAPtr).getLength();
	if(SMA_Type == DOUBLE){
		NimArrBase<double>* SMA_NimArrPtr = static_cast<NimArrBase<double>*>(SMA_NimTypePtr);
		std::copy(nimArr.getPtr() + nimBegin,
		 nimArr.getPtr() + SMA_length + nimBegin,
		 SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() );
	}
	else if(SMA_Type == INT){
		NimArrBase<int>* SMA_NimArrPtr = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
		std::copy(nimArr.getPtr() + nimBegin,
		 nimArr.getPtr() + SMA_length +nimBegin,
		 SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() );
	}
	else {
		PRINTF("Copying type for nimArr_2_SingleModelAccess not supported\n");
	}
}


//template<int D, class T>
//void nimArr_2_ManyModelAccess(ManyVariablesAccessor &MMVAPtr, NimArr<D, T> &nimArr){
template<class T>
void nimArr_2_ManyModelAccess(ManyVariablesAccessor &MMVAPtr, NimArrBase<T> &nimArr) {
  vector<SingleVariableAccessBase*> *SMVA_Vec = &(MMVAPtr.getAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int k = SMVA_Vec->size();
  int nextNumVals;
  SingleVariableAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = (*SMVA_Vec)[i];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd){
      PRINTF("Warning: in nimArr_2_ManyModelAccess, accessor larger than NimArr!\n");
      break;
    }
    nimArr_2_SingleModelAccess<T>(curSingleAccess, nimArr, nimCurrent);
    nimCurrent = nimCurrent + nextNumVals;
  }
  if(nimCurrent != nimEnd)
    PRINTF("Warning: after completing nimArr_2_ManyModelAccess, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
}

template<int D, class T>
void SingleModelAccess_2_nimArr(SingleVariableAccessBase* SMVAPtr, NimArr<D, T> &nimArr, int nimBegin) {
  NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
  nimType SMA_Type = (*SMA_NimTypePtr).getNimType();
  int SMA_length = (*SMVAPtr).getLength();
  if(SMA_Type == DOUBLE){
    NimArrBase<double>* SMA_NimArrPtr = static_cast<NimArrBase<double>*>(SMA_NimTypePtr);
    std::copy(SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart(),
	      SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() + SMA_length ,
	      nimArr.getPtr()  + nimBegin);
  }
  else if(SMA_Type == INT){
    NimArrBase<int>* SMA_NimArrPtr = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
    std::copy(SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart(),
	      SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() + SMA_length ,
	      nimArr.getPtr()  + nimBegin);
  }
  else {
    PRINTF("Copying type for nimArr_2_SingleModelAccess not supported\n");
  }
}


template<int D, class T>
void ManyModelAccess_2_nimArr(ManyVariablesAccessor &MMVAPtr, NimArr<D, T> &nimArr) {
  vector<SingleVariableAccessBase*> *SMVA_Vec = &(MMVAPtr.getAccessVector());
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int k = SMVA_Vec->size();
  int nextNumVals;
  SingleVariableAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = (*SMVA_Vec)[i];
    nextNumVals = (*curSingleAccess).getLength();
    if(nextNumVals + nimCurrent > nimEnd){
      PRINTF("Warning: in ManyModelAccess_2_nimArr, accessor larger than NimArr!\n");
      break;
    }
    SingleModelAccess_2_nimArr<D, T>(curSingleAccess, nimArr, nimCurrent);
    nimCurrent = nimCurrent + nextNumVals;
  }
  if(nimCurrent != nimEnd)
    PRINTF("Warning: after completing ManyModelAccess_2_nimArr, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
}

void setValues(NimArrBase<double> &nimArr, ManyVariablesAccessor &MVA){
	nimArr_2_ManyModelAccess<double>(MVA, nimArr);
}
void setValues(NimArrBase<int> &nimArr, ManyVariablesAccessor &MVA){
	nimArr_2_ManyModelAccess<int>(MVA, nimArr);
}

void getValues(NimArr<1, double> &nimArr, ManyVariablesAccessor &MVA){
	ManyModelAccess_2_nimArr<1, double>(MVA, nimArr);
}
void getValues(NimArr<1, int> &nimArr, ManyVariablesAccessor &MVA){
	ManyModelAccess_2_nimArr<1, int>(MVA, nimArr);
}


// new copierClass versions
// remember to look at calculate() too to avoid copies every time.
void nimCopy(const copierVectorClass &copiers) {
  vector<copierClass*>::const_iterator iCopy;
  //  PRINTF("iterating over %i\n", copiers.copyVector.size());
  // int i = 0;
  for(iCopy = copiers.copyVector.begin(); iCopy != copiers.copyVector.end(); iCopy++) {
    //   PRINTF("Starting copier %i\n", i);
    (*iCopy)->copy(copiers.rowInfo);
    //    PRINTF("Done with copier %i\n", i);
    //   i++;
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
  //PRINTF("Ready to set up length %i\n", copyVector.size());
  vector<SingleVariableMapAccessBase *>::iterator iFrom, iTo, iFromEnd;
  iFromEnd = fromAccessors->end();
  iTo =  toAccessors->begin();
  int i = 0;
  for(iFrom = fromAccessors->begin(); iFrom != iFromEnd; iFrom++) {
    //PRINTF("setting up %i\n", i);
    copyVector[i] = makeOneCopyClass(*iFrom, *iTo, isFromMV, isToMV); // switched from isFromMV and isToMV
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

// new MapAccessor versions
void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to) { //map version

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
  //  dynamicMapCopy<double, double>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
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
    case DOUBLE:
      switch(toType) {
      case DOUBLE:
	(*static_cast<NimArrBase<double> *>(toNimArr))[to->getOffset()] = (*static_cast<NimArrBase<double> *>(fromNimArr))[from->getOffset()];
	break;
    case INT:
	(*static_cast<NimArrBase<int> *>(toNimArr))[to->getOffset()] = (*static_cast<NimArrBase<double> *>(fromNimArr))[from->getOffset()];
      break;
    default:
      _nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
      nimble_print_to_R(_nimble_global_output);
    }
    break;
  case INT:
    switch(toType) {
    case DOUBLE:
      (*static_cast<NimArrBase<double> *>(toNimArr))[to->getOffset()] = (*static_cast<NimArrBase<int> *>(fromNimArr))[from->getOffset()];
      break;
    case INT:
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
    case DOUBLE:
      switch(toType) {
      case DOUBLE:
	dynamicMapCopy<double, double>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      case INT:
	dynamicMapCopy<double, int>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      default:
	_nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
	nimble_print_to_R(_nimble_global_output);
      }
      break;
    case INT:
      switch(toType) {
      case DOUBLE:
	dynamicMapCopy<int, double>(toNimArr, to->getOffset(), to->getStrides(), to->getSizes(), fromNimArr, from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      case INT:
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
  case INT:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->size(); //getVptr()->size();
    break;
  case DOUBLE:
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
  case INT:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->size(); //getVptr()->size();
    break;
  case DOUBLE:
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

// old versions:
void nimCopy(ManyVariablesAccessorBase &from, ManyVariablesAccessorBase &to) {

  vector<SingleVariableAccessBase *> *fromAccessors = &(from.getAccessVector());
  vector<SingleVariableAccessBase *> *toAccessors = &(to.getAccessVector());

  if(fromAccessors->size() != toAccessors->size()) {
    _nimble_global_output<<"Error in nimCopy: from and to access vectors have sizes "<<fromAccessors->size() << " and " << toAccessors->size() << "\n";
    nimble_print_to_R(_nimble_global_output);
  }

  vector<SingleVariableAccessBase *>::iterator iFrom, iTo, iFromEnd;
  iFrom = fromAccessors->begin();
  iFromEnd = fromAccessors->end();
  iTo =  toAccessors->begin();
  for( ; iFrom != iFromEnd; ++iFrom) {
    nimCopyOne(*iFrom, *iTo);
    ++iTo;
  }
}


void nimCopy(ManyVariablesAccessorBase &from, int rowFrom, ManyVariablesAccessorBase &to) {
	from.setRow(rowFrom - 1);
	nimCopy(from, to);
}

void nimCopy(ManyVariablesAccessorBase &from, int rowFrom, ManyVariablesAccessorBase &to, int rowTo) {
	to.setRow(rowTo - 1);
	from.setRow(rowFrom - 1);
	nimCopy(from, to);
}

void nimCopy(ManyVariablesAccessorBase &from, ManyVariablesAccessorBase &to, int rowTo) {
	to.setRow(rowTo - 1);
	nimCopy(from, to);
}




void nimCopyOne(SingleVariableAccessBase *from, SingleVariableAccessBase *to) {
  nimType fromType, toType;
  NimArrType *fromNimArr, *toNimArr;
  fromNimArr = from->getNimArrPtr();
  toNimArr = to->getNimArrPtr();
  fromType = fromNimArr->getNimType();
  toType = toNimArr->getNimType();
  switch(fromType) {
  case DOUBLE:
    switch(toType) {
    case DOUBLE:
      nimCopyOneTyped<double, double>(from, to);
      break;
    case INT:
      nimCopyOneTyped<double, int>(from, to);
      break;
    default:
      _nimble_global_output<<"Error in nimCopyOne: unknown type for destination\n";
      nimble_print_to_R(_nimble_global_output);
    }
    break;
  case INT:
    switch(toType) {
//    int *ifromPtr =  static_cast<NimArrBase<int> *>( (from)->getPtr() );
    case DOUBLE:
      nimCopyOneTyped<int, double>(from, to);
      break;
    case INT:
      nimCopyOneTyped<int, int>(from, to);
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
}


SingleModelValuesAccess* cMakeSingleModelValuesAccessor(NimVecType* varPtr, int beginIndex, int endIndex, int row){
	SingleModelValuesAccess* singleMValuesPtr = new SingleModelValuesAccess;
	(*singleMValuesPtr).flatIndexStart = beginIndex;
	(*singleMValuesPtr).flatIndexEnd = endIndex;
	(*singleMValuesPtr).length = endIndex - beginIndex + 1;
	(*singleMValuesPtr).currentRow = row;
	(*singleMValuesPtr).pVVar = varPtr;
	return(singleMValuesPtr);
}



SingleVariableAccess* cMakeSingleVariableAccessor(NimArrType** varPtr, int beginIndex, int endIndex){
	SingleVariableAccess* sVAPtr = new SingleVariableAccess;
	sVAPtr->flatIndexStart = beginIndex;
	sVAPtr->flatIndexEnd = endIndex;
	sVAPtr->length = endIndex - beginIndex + 1;
	sVAPtr->ppVar = varPtr;
	return(sVAPtr);
}

SEXP getMVAccessorValues(SEXP accessor){
	void* vPtr = R_ExternalPtrAddr(accessor);
	if(vPtr == NULL)
		return(R_NilValue);
	SingleModelValuesAccess* SVAptr = static_cast<SingleModelValuesAccess*>(vPtr);
	NimArrType* nimTypePtr = (*SVAptr).getNimArrPtr();
	nimType arrType = (*nimTypePtr).getNimType();
	SEXP rOutput;
	if(arrType == DOUBLE){
		PROTECT(rOutput = allocVector(REALSXP, (*SVAptr).getLength() ) );
		NimArrBase<double>* NimArrPtr = static_cast<NimArrBase<double>*>(nimTypePtr) ;
//		std::copy(	NimArrPtr->getPtr() + SVAptr->getIndexStart(),
//			 	NimArrPtr->getPtr() + SVAptr->getIndexEnd(),
//			 	REAL(rOutput) );			// Not sure why this isn't working...
		int begin = SVAptr->getIndexStart();
		int length = SVAptr->getLength();
		for(int i = 0; i <length; i++)
			REAL(rOutput)[i] = (*NimArrPtr)[begin + i];
		UNPROTECT(1);
		return(rOutput);
	}
	if(arrType == INT){
		PROTECT(rOutput = allocVector(INTSXP, (*SVAptr).getLength() ) );
		NimArrBase<int>* NimArrPtr = static_cast<NimArrBase<int>*>(nimTypePtr) ;
		int begin = SVAptr->getIndexStart();
		int length = SVAptr->getLength();
		for(int i = 0; i <length; i++)
			INTEGER(rOutput)[i] = (*NimArrPtr)[begin + i];
		UNPROTECT(1);
		return(rOutput);
	}
	return(R_NilValue);
}

SEXP getModelAccessorValues(SEXP accessor){
	void* vPtr = R_ExternalPtrAddr(accessor);
	if(vPtr == NULL)
		return(R_NilValue);
	SingleVariableAccess* SVAptr = static_cast<SingleVariableAccess*>(vPtr);
	NimArrType* nimTypePtr = (*SVAptr).getNimArrPtr();
	nimType arrType = (*nimTypePtr).getNimType();
	SEXP rOutput;
	if(arrType == DOUBLE){
		PROTECT(rOutput = allocVector(REALSXP, (*SVAptr).getLength() ) );
		NimArrBase<double>* NimArrPtr = static_cast<NimArrBase<double>*>(nimTypePtr) ;
//		std::copy(	NimArrPtr->getPtr() + SVAptr->getIndexStart(),
//			 	NimArrPtr->getPtr() + SVAptr->getIndexEnd(),
//			 	REAL(rOutput) );			// Not sure why this isn't working...
		int begin = SVAptr->getIndexStart();
		int length = SVAptr->getLength();
		for(int i = 0; i <length; i++)
			REAL(rOutput)[i] = (*NimArrPtr)[begin + i];
		UNPROTECT(1);
		return(rOutput);
	}
	if(arrType == INT){
		PROTECT(rOutput = allocVector(INTSXP, (*SVAptr).getLength() ) );
		NimArrBase<int>* NimArrPtr = static_cast<NimArrBase<int>*>(nimTypePtr) ;
		int begin = SVAptr->getIndexStart();
		int length = SVAptr->getLength();
		for(int i = 0; i <length; i++)
			INTEGER(rOutput)[i] = (*NimArrPtr)[begin + i];
		UNPROTECT(1);
		return(rOutput);
	}

	return(R_NilValue);
}

SEXP getListElement(SEXP list, const char *str){
	SEXP ans = R_NilValue, names = getAttrib(list, R_NamesSymbol);
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


SEXP populateNumberedObject_withSingleModelValuesAccessors(SEXP mvPtr, SEXP varName, SEXP GIDs, SEXP curRow, SEXP SnumbObj){
//	int cIndex = INTEGER(beginIndex)[0] - 1;
	int len = LENGTH(GIDs);
	int cRow = INTEGER(curRow)[0] - 1;
	int* cGIDs = INTEGER(GIDs);
	string vName = STRSEXP_2_string(varName, 0);
	Values* modelValuesPtr = static_cast<Values*>(R_ExternalPtrAddr(mvPtr));
	NimVecType* nimPtr = static_cast<NimVecType*> (modelValuesPtr->getObjectPtr(vName));
	SingleModelValuesAccess* smva;
	NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumbObj) );
	void* vPtr;
	for(int i = 0; i < len; i++){
		smva = cMakeSingleModelValuesAccessor(nimPtr, i, i, cRow);
		vPtr = static_cast<void*>(smva);
		numObj->numberedObjects[cGIDs[i] - 1] = vPtr;
	}
	return(R_NilValue);
}

SEXP populateNumberedObject_withSingleModelVariablesAccessors(SEXP modelPtr, SEXP varName, SEXP sGIDS, SEXP SvalidIndices, SEXP SnumbObj){
	int len = LENGTH(sGIDS);
	int* validIndices = INTEGER(SvalidIndices);
	int* gids = INTEGER(sGIDS);
	string vName = STRSEXP_2_string(varName, 0);
	ModelBase* modelValuesPtr = static_cast<ModelBase*>(R_ExternalPtrAddr(modelPtr));
	NimArrType** nimPtr = static_cast<NimArrType**> (modelValuesPtr->getObjectPtr(vName));
	SingleVariableAccess* smva;
	NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumbObj) );
	void* vPtr;
	for(int i = 0; i < len; i++){
		smva = cMakeSingleVariableAccessor(nimPtr, validIndices[i] - 1, validIndices[i] - 1);
		vPtr = static_cast<void*>(smva);
		numObj->numberedObjects[gids[i] - 1] = vPtr;
	}
	return(R_NilValue);
}


SEXP populateNodeFxnVectorNew_byDeclID(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS){
  //std::cout<<"in populateNodeFxnVectorNew_byDeclID\n";
  int len = LENGTH(S_ROWINDS);
  if(len == 0) return(R_NilValue);
  //std::cout<<"len = "<<len<<"\n";
  int* gids = INTEGER(S_GIDs);
  int* rowinds = INTEGER(S_ROWINDS);
  int index;
  NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumberedObj));
  NodeVectorClassNew* nfv = static_cast<NodeVectorClassNew*>(R_ExternalPtrAddr(SnodeFxnVec) ) ;
  //  (*nfv).instructions.resize(len);
  //int previousIndex = -1;
  int nextRowInd;
  for(int i = 0; i < len; i++){
    index = gids[i] - 1;
    //    std::cout<<"index "<<index<<" i "<<i<<" rowinds[i]-1 "<<rowinds[i]-1<<"\n";
    nextRowInd = rowinds[i]-1;
    if(nextRowInd == -1) { // should only happen from a scalar, so there is one dummy indexedNodeInfo
      nextRowInd = 0;
    }
    (*nfv).instructions.push_back(NodeInstruction(static_cast<nodeFun*>(numObj->getObjectPtr(index)), nextRowInd));
  }
  //  std::cout<<"done with "<<(*nfv).instructions.size()<<"\n";
  return(R_NilValue);
}

SEXP populateIndexedNodeInfoTable(SEXP StablePtr, SEXP StableContents) {
  SEXP Sdim;
  //  std::cout<<"in populateIndexedNodeInfoTable\n";
  PROTECT(Sdim = getAttrib(StableContents, R_DimSymbol));
  if(LENGTH(Sdim) != 2) {PRINTF("Warning from populateIndexedNodeInfoTable: LENGTH(Sdim) != 2"); return(R_NilValue);}
  int nrow = INTEGER(Sdim)[0];
  int ncol = INTEGER(Sdim)[1];
  //std::cout<<"nrow "<<nrow<<" ncol "<<ncol<<"\n";
  vector<indexedNodeInfo> *tablePtr = static_cast<vector<indexedNodeInfo> *>(R_ExternalPtrAddr(StablePtr));
  if(nrow == 0) {
    void *vptr=0;
    tablePtr->push_back(indexedNodeInfo(static_cast<int *>(vptr), 0, 0));
    if(ncol != 0) {PRINTF("Warning from populateIndexedNodeInfoTable: nrow == 0 but ncol != 0.");}
  } else {

    if(!isNumeric(StableContents)) {PRINTF("Warning from populateIndexedNodeInfoTable: StableContents is not numeric"); return(R_NilValue);}
    if(isInteger(StableContents)) {
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
  //  std::cout<<"done with size "<<tablePtr->size()<<"\n";
  UNPROTECT(1);
  return(R_NilValue);
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

///NEW

string _NIMBLE_WHITESPACE(" \t");
string _NIMBLE_WHITESPACEBRACKET(" \t[");
string _NIMBLE_NUMERICS("0123456789.");
string _NIMBLE_SPACECOMMABRACKET(" ,]");

// std::stoi not consistent across C++ and not worth the portability worries, so we made our own using istringstream
int nimble_stoi(const string &input) {
  istringstream converter;
  std::size_t iStart(input.find_first_not_of(_NIMBLE_WHITESPACE));
  std::size_t iEnd(input.find_first_not_of(_NIMBLE_NUMERICS, iStart));
  if(iEnd != std::string::npos && iEnd > iStart) iEnd--;
  converter.str(input.substr(iStart, iEnd - iStart + 1));
  int ans;
  converter >> ans;
  return ans;
}

class varAndIndicesClass {
public:
  string varName;
  vector< vector<int> > indices;
};

class mapInfoClass {
public:
  int offset;
  vector<int> sizes;
  vector<int> strides;
};

void parseVarAndInds(const string &input, varAndIndicesClass &output) { //string &varName, vector<vector<int> > &inds) {
  output.indices.resize(0);
  std::size_t iBracket = input.find_first_of('[');
  //std::cout<<iBracket<<"\n";
  if(iBracket == std::string::npos) { // no bracket
    output.varName = input;
    //    std::cout<<output.varName<<"\n";
    return;
  }
  output.varName = input.substr(0, iBracket);
  //  std::cout<<output.varName<<"\n";

  string restOfInput = input.substr(iBracket+1);
  //  std::cout<<restOfInput<<"\n";
  bool done(false);

  vector<int> oneAns;
  std::size_t iColon, iComma;
  int firstNum, secondNum;
  std::size_t iNextStart, iNonBlank;
  iBracket = restOfInput.find_first_of(']');
  if(iBracket == std::string::npos) {
    _nimble_global_output<<"problem in parseVarAndInds: there is no closing bracket\n";
    nimble_print_to_R(_nimble_global_output);
  }
  while(!done) {
    iColon   = restOfInput.find_first_of(':');
    iComma   = restOfInput.find_first_of(',');
    if((iColon < iBracket) & (iColon < iComma)) { // next is a colon expr like 2:5
      firstNum = nimble_stoi(restOfInput);
      // test x[11 :4]
      iNextStart = iColon + 1;
      restOfInput = restOfInput.substr(iNextStart);
      iComma   = restOfInput.find_first_of(',');
      iBracket = restOfInput.find_first_of(']');
      //std::cout<<restOfInput<<"\n";
      secondNum   = nimble_stoi(restOfInput);
      if(iComma < iBracket) iNextStart = iComma + 1; else {iNextStart = iBracket; done = true;}
      restOfInput = restOfInput.substr(iNextStart);
      //    std::cout<<"got "<<firstNum<<" : "<<secondNum<<"\n";
      //    std::cout<<restOfInput<<"\n";
      oneAns.push_back(firstNum);
      oneAns.push_back(secondNum);
      output.indices.push_back(oneAns);
      oneAns.clear();
    } else {
      // test for blanks
      // this bit ends in either a comma or the bracket
      if(iComma >= iBracket) {iComma = iBracket; done = true;} // now iComma is the ending index after here
      iNonBlank = restOfInput.find_first_not_of(_NIMBLE_SPACECOMMABRACKET);
      if(iNonBlank < iComma) { // there is a number
	firstNum = nimble_stoi(restOfInput);
	if(iComma < iBracket) iNextStart = iComma + 1; else iNextStart = iBracket;
	//	  if(iEndOfNum < iColon) iEndOfNum = iColon;
	restOfInput = restOfInput.substr(iNextStart);
	//	std::cout<<"got "<<firstNum<<"\n";
	//	std::cout<<restOfInput<<"\n";
	oneAns.push_back(firstNum);
	output.indices.push_back(oneAns);
	oneAns.clear();
      } else { // there is a blank
	output.indices.push_back(oneAns);
	if(iComma < iBracket) iNextStart = iComma + 1; else iNextStart = iBracket;
	restOfInput = restOfInput.substr(iNextStart);
	//	std::cout<<"got blank\n";
	//	std::cout<<restOfInput<<"\n";
      }
    }
    iBracket = restOfInput.find_first_of(']');
  }
}

SEXP varAndIndices2Rlist(const varAndIndicesClass &input) {
  SEXP Soutput, Sindices;
  PROTECT(Soutput = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Soutput, 0, string_2_STRSEXP(input.varName));
  int numinds = input.indices.size();
  PROTECT(Sindices = allocVector(VECSXP, numinds));
  for(int i = 0; i < numinds; i++) {
    SET_VECTOR_ELT(Sindices, i, vectorInt_2_SEXP(input.indices[i]));
  }
  SET_VECTOR_ELT(Soutput, 1, Sindices);

  vector<string> newNames(2);
  newNames[0].assign("varName");
  newNames[1].assign("indices");
  SEXP SnewNames;
  PROTECT(SnewNames = vectorString_2_STRSEXP(newNames));
  setAttrib(Soutput, R_NamesSymbol, SnewNames);

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
  PROTECT(Soutput = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(Soutput, 0, int_2_SEXP(input.offset));
  SET_VECTOR_ELT(Soutput, 1, vectorInt_2_SEXP(input.sizes));
  SET_VECTOR_ELT(Soutput, 2, vectorInt_2_SEXP(input.strides));
  vector<string> newNames(3);
  newNames[0].assign("offset");
  newNames[1].assign("sizes");
  newNames[2].assign("strides");
  SEXP SnewNames;
  PROTECT(SnewNames = vectorString_2_STRSEXP(newNames));
  setAttrib(Soutput, R_NamesSymbol, SnewNames);
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

// call both with a bunch of output generated...
SEXP populateValueMapAccessorsFromNodeNames(SEXP StargetPtr, SEXP SnodeNames, SEXP SsizesAndNdims, SEXP SModelOrModelValuesPtr ) {
  vector<string> nodeNames;
  STRSEXP_2_vectorString(SnodeNames, nodeNames);
  NamedObjects *sourceNamedObject = static_cast<NamedObjects*>(R_ExternalPtrAddr(SModelOrModelValuesPtr));
  ManyVariablesMapAccessorBase* valuesAccessor = static_cast<ManyVariablesMapAccessorBase*>(R_ExternalPtrAddr(StargetPtr) );
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

SEXP populateModelVariablesAccessors_byGID(SEXP SmodelVariableAccessorVector, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_LP_GIDs, SEXP S_LP_numberedObj){
	int len = LENGTH(S_GIDs);
	int* gids = INTEGER(S_GIDs);
	int len_LP = LENGTH(S_LP_GIDs);
	int* LP_gids = INTEGER(S_LP_GIDs);
	int index;
	NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumberedObj));
	NumberedObjects* LP_numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(S_LP_numberedObj));
	ManyVariablesAccessor* accessVector = static_cast<ManyVariablesAccessor*>(R_ExternalPtrAddr(SmodelVariableAccessorVector) );
	accessVector->varAccessors.resize(len + len_LP);
	for(int i = 0; i < len; i++){
		index = gids[i] - 1;
		accessVector->varAccessors[i] = static_cast<SingleVariableAccessBase*>(numObj->getObjectPtr(index));
	//	Rprintf("skipping the initialization of modelVariable accessor for debugging purposes");
	}
	for(int i = 0; i < len_LP; i++){
		index = LP_gids[i] - 1;
		accessVector->varAccessors[len + i] = static_cast<SingleVariableAccessBase*>(LP_numObj->getObjectPtr(index) ) ;
	//	Rprintf("skipping the initialization of logProbs for debugging purposes\n");
	}
	return(R_NilValue);
}

SEXP removeModelVariableAccessor(SEXP rPtr, SEXP index, SEXP removeAll){
	int cIndex = INTEGER(index)[0] - 1;
	bool cRemoveAll = LOGICAL(removeAll)[0];
	void* vPtr = R_ExternalPtrAddr(rPtr);
	if(vPtr == NULL){
		PRINTF("Warning: pointer to null passed to removeNodeFun\n");
		return(R_NilValue);
	}
	ManyVariablesAccessor* mVAPtr = static_cast<ManyVariablesAccessor*>(vPtr);
	cRemoveAccessor<ManyVariablesAccessor>(mVAPtr, cIndex, cRemoveAll);
	return(R_NilValue);
}

SEXP removeModelValuesAccessor(SEXP rPtr, SEXP index, SEXP removeAll){
	int cIndex = INTEGER(index)[0] - 1;
	bool cRemoveAll = LOGICAL(removeAll)[0];
	void* vPtr = R_ExternalPtrAddr(rPtr);
	if(vPtr == NULL){
		PRINTF("Warning: pointer to null passed to removeNodeFun\n");
		return(R_NilValue);
	}
	ManyModelValuesAccessor* mMVAPtr = static_cast<ManyModelValuesAccessor*>(vPtr);
	cRemoveAccessor<ManyModelValuesAccessor>(mMVAPtr, cIndex, cRemoveAll);
	return(R_NilValue);
}


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

SEXP resizeManyModelVarAccessor(SEXP manyModelVarPtr, SEXP size){
	int cSize = INTEGER(size)[0];
	ManyVariablesAccessor* MMVVec = static_cast<ManyVariablesAccessor*>(R_ExternalPtrAddr(manyModelVarPtr) ) ;
	(*MMVVec).varAccessors.resize(cSize);
	return(R_NilValue);
}

SEXP resizeManyModelValuesAccessor(SEXP manyModelValuesPtr, SEXP size){
	int cSize = INTEGER(size)[0];
	ManyModelValuesAccessor* MMVVec = static_cast<ManyModelValuesAccessor*>(R_ExternalPtrAddr(manyModelValuesPtr) ) ;
	(*MMVVec).varAccessors.resize(cSize);
	return(R_NilValue);
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

SEXP addSingleVariableAccessor(SEXP MVAPtr, SEXP SVAPtr, SEXP addAtEnd, SEXP index){
	void* vMVAPtr = R_ExternalPtrAddr(MVAPtr);
	void* vSVAPtr = R_ExternalPtrAddr(SVAPtr);
	if((vMVAPtr == NULL) | (vSVAPtr == NULL)){
		PRINTF("Warning: pointer to null passed to addNodeFun\n");
		return(R_NilValue);
	}
	int cIndex = - 1;
	bool cAddAtEnd = LOGICAL(addAtEnd)[0];
	if(cAddAtEnd == FALSE)
		cIndex = INTEGER(index)[0] - 1;
	ManyVariablesAccessor* cMVAPtr = static_cast<ManyVariablesAccessor*>(vMVAPtr);
	SingleVariableAccess* cSVAPtr = static_cast<SingleVariableAccess*>(vSVAPtr);
	cAddAccessor<ManyVariablesAccessor, SingleVariableAccess>(cMVAPtr, cSVAPtr, cAddAtEnd, cIndex);
	return(R_NilValue);
}

SEXP addSingleModelValuesAccessor(SEXP MMVAPtr, SEXP SMVAPtr, SEXP addAtEnd, SEXP index){
	void* vMMVAPtr = R_ExternalPtrAddr(MMVAPtr);
	void* vSMVAPtr = R_ExternalPtrAddr(SMVAPtr);
	if((vMMVAPtr == NULL) | (vSMVAPtr == NULL)){
		PRINTF("Warning: pointer to null passed to addNodeFun\n");
		return(R_NilValue);
	}
	int cIndex = INTEGER(index)[0] - 1;
	bool cAddAtEnd = LOGICAL(addAtEnd)[0];
	ManyModelValuesAccessor* cMMVAPtr = static_cast<ManyModelValuesAccessor*>(vMMVAPtr);
	SingleModelValuesAccess* cSMVAPtr = static_cast<SingleModelValuesAccess*>(vSMVAPtr);
	cAddAccessor<ManyModelValuesAccessor, SingleModelValuesAccess>(cMMVAPtr, cSMVAPtr, cAddAtEnd, cIndex);
	return(R_NilValue);
}

SEXP manualSetNRows(SEXP Sextptr, SEXP nRows){
  	Values* vPtr = static_cast<Values*>(R_ExternalPtrAddr(Sextptr) );
  	(*vPtr).numRows = INTEGER(nRows)[0];
  return(R_NilValue);
  }

