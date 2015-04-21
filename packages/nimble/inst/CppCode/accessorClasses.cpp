#include "nimble/accessorClasses.h"
#include "nimble/RcppUtils.h"

// 1. NodeVectors
double calculate(NodeVectorClass &nodes) {
  double ans(0);
  vector<nodeFun *> nodeFunPtrs = nodes.getNodeFunctionPtrs();
  int vecSize = nodeFunPtrs.size();
  for(int i = 0; i < vecSize; i++)
  	 	ans +=	nodeFunPtrs[i]->calculate();
  return(ans);
}

double getLogProb(NodeVectorClass &nodes) {
  double ans(0);
  vector<nodeFun *> nodeFunPtrs = nodes.getNodeFunctionPtrs();
  vector<nodeFun *>::iterator endNode(nodeFunPtrs.end());
  for( vector<nodeFun *>::iterator iNodes = nodeFunPtrs.begin();
       iNodes != endNode;
       ++iNodes) {
    ans += (*iNodes)->getLogProb();
  }
  return(ans);
}

void simulate(NodeVectorClass &nodes) {
  vector<nodeFun *> nodeFunPtrs = nodes.getNodeFunctionPtrs();
  vector<nodeFun *>::iterator endNode(nodeFunPtrs.end());
  for( vector<nodeFun *>::iterator iNodes = nodeFunPtrs.begin();
       iNodes != endNode;
       ++iNodes) {
     (*iNodes)->simulate();
  }
}





// 2. We will only delete the singleVariableAccessors when we delete the ManyVariablesAccesor that they are contained in
		// Cliff's note: actually, we have changed the paradigm. We now make one variable accessor per node, which can appear in many different
		// ManyVariableAccessors. All unique singleAccessors will be held in a single SpecializedNumberedObjects, which takes care of the finalizer.
ManyVariablesAccessor::~ManyVariablesAccessor(){
//	int k = varAccessors.size();
//	Rprintf("skipping the deletion of singleVariableAccesses for debugging purposes\n");
//	for(int i = 0; i < k; i++)
//		delete static_cast<SingleVariableAccess*>(varAccessors[i]);
}


SingleVariableMapAccessBase::~SingleVariableMapAccessBase(){}

ManyVariablesMapAccessor::~ManyVariablesMapAccessor(){
    for(int i = 0; i < varAccessors.size(); ++i)
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

ManyModelValuesMapAccessor::ManyModelValuesMapAccessor() : currentRow(0) {
  
}


ManyModelValuesMapAccessor::~ManyModelValuesMapAccessor() {
    for(int i = 0; i < varAccessors.size(); ++i)
      delete static_cast<SingleModelValuesMapAccess*>(varAccessors[i]);
};

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
template<class T>
void nimArr_2_SingleModelAccess(SingleVariableMapAccessBase* SMVAPtr, NimArrBase<T> &nimArr, int nimBegin, int nimStride){
  NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
  nimType SMA_Type = (*SMA_NimTypePtr).getNimType();
  NimArrBase<double>* SMA_NimArrPtrD;
  NimArrBase<int>* SMA_NimArrPtrI;
  
  if(SMVAPtr->getSingleton()) {
    switch(SMA_Type) {
    case DOUBLE:
      (*static_cast<NimArrBase<double> *>(SMA_NimTypePtr))[SMVAPtr->offset] = (*nimArr.getVptr())[nimBegin];
      break;
    case INT:
      (*static_cast<NimArrBase<double> *>(SMA_NimTypePtr))[SMVAPtr->offset] = (*nimArr.getVptr())[nimBegin];
      break;
    default:
      PRINTF("Copying type for nimArr_2_SingleModelAccess not supported\n");
      break;
    }
  } else {
    switch(SMA_Type) {
    case DOUBLE:
      SMA_NimArrPtrD = static_cast<NimArrBase<double>*>(SMA_NimTypePtr);
      dynamicMapCopyFlatToDim<T, double>(SMA_NimArrPtrD, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes(), &nimArr, nimBegin, nimStride);
      // std::copy(nimArr.getPtr() + nimBegin,
      // 	      nimArr.getPtr() + SMA_length + nimBegin, 
      // 	      SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() );
      break;
    case INT:
      SMA_NimArrPtrI = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
      dynamicMapCopyFlatToDim<T, int>(SMA_NimArrPtrI, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes(), &nimArr, nimBegin, nimStride);
      // std::copy(nimArr.getPtr() + nimBegin,
      // 	      nimArr.getPtr() + SMA_length +nimBegin, 
      // 	      SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() );
      break;
    default:
      PRINTF("Copying type for nimArr_2_SingleModelAccess not supported\n");
      break;
    }
  }
}

template<class T>
void nimArr_2_ManyModelAccess(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr){
  vector<SingleVariableMapAccessBase*> SMVA_Vec = MMVAPtr.getMapAccessVector();
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = SMVA_Vec.size();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = SMVA_Vec[i];
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

///////////////
// [accessors]_2_nimArr
// nimArr is "to". SMVAPtr is "from"
template<class T>
void SingleModelAccess_2_nimArr(SingleVariableMapAccessBase* SMVAPtr, NimArrBase<T> &nimArr, int nimBegin, int nimStride){
  NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
  nimType SMA_Type = (*SMA_NimTypePtr).getNimType();
  NimArrBase<double>* SMA_NimArrPtrD;
  NimArrBase<int>* SMA_NimArrPtrI;
  if(SMVAPtr->getSingleton()) {
    switch(SMA_Type) {
    case DOUBLE:
      (*nimArr.getVptr())[nimBegin] = (*static_cast<NimArrBase<double> *>(SMA_NimTypePtr))[SMVAPtr->offset];
      break;
    case INT:
      (*nimArr.getVptr())[nimBegin] = (*static_cast<NimArrBase<int> *>(SMA_NimTypePtr))[SMVAPtr->offset];
      break;
    default:
      PRINTF("Copying type for SingleModelAccess_2_nimArr not supported\n");
      break;
    }
  } else {
    switch(SMA_Type) {
    case DOUBLE:
      SMA_NimArrPtrD = static_cast<NimArrBase<double>*>(SMA_NimTypePtr);
      dynamicMapCopyDimToFlat<double, T>(&nimArr, nimBegin, nimStride, SMA_NimArrPtrD, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes() );
      
      // 	std::copy(SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart(),
      // 	SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() + SMA_length ,
      // 	 nimArr.getPtr()  + nimBegin);

      break;
    case INT:
      SMA_NimArrPtrI = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
      dynamicMapCopyDimToFlat<int, T>(&nimArr, nimBegin, nimStride, SMA_NimArrPtrI, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes());
      
      // 	std::copy(SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart(),
      // 	SMA_NimArrPtr->getPtr() + SMVAPtr->getIndexStart() + SMA_length ,
      // 	 nimArr.getPtr()  + nimBegin);
      break;
    default:
      PRINTF("Copying type for SingleModelAccess_2_nimArr not supported\n");
      break;
    }
  }
}


template<class T>
void ManyModelAccess_2_nimArr(ManyVariablesMapAccessor &MMVAPtr, NimArrBase<T> &nimArr){
  vector<SingleVariableMapAccessBase*> SMVA_Vec = MMVAPtr.getMapAccessVector();
  int nimCurrent = 0;
  int nimEnd = nimArr.size();
  int nimArrStride = nimArr.strides()[0];
  int nimCurrentOffset = nimArr.getOffset();
  int k = SMVA_Vec.size();
  int nextNumVals;
  SingleVariableMapAccessBase* curSingleAccess;
  for(int i = 0; i < k ; i++){
    curSingleAccess = SMVA_Vec[i];
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

//////////
//
void setValues(NimArrBase<double> &nimArr, ManyVariablesMapAccessor &MVA){
	nimArr_2_ManyModelAccess<double>(MVA, nimArr);
}

void setValues(NimArrBase<int> &nimArr, ManyVariablesMapAccessor &MVA){
	nimArr_2_ManyModelAccess<int>(MVA, nimArr);
}

void getValues(NimArr<1, double> &nimArr, ManyVariablesMapAccessor &MVA){
	ManyModelAccess_2_nimArr<double>(MVA, nimArr);
}
void getValues(NimArr<1, int> &nimArr, ManyVariablesMapAccessor &MVA){
	ManyModelAccess_2_nimArr<int>(MVA, nimArr);
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
void nimArr_2_ManyModelAccess(ManyVariablesAccessor &MMVAPtr, NimArrBase<T> &nimArr){
	vector<SingleVariableAccessBase*> SMVA_Vec = MMVAPtr.getAccessVector();
	int nimCurrent = 0;
	int nimEnd = nimArr.size();
	int k = SMVA_Vec.size();
	int nextNumVals;
	SingleVariableAccessBase* curSingleAccess;
	for(int i = 0; i < k ; i++){
		curSingleAccess = SMVA_Vec[i];
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
void SingleModelAccess_2_nimArr(SingleVariableAccessBase* SMVAPtr, NimArr<D, T> &nimArr, int nimBegin){
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
void ManyModelAccess_2_nimArr(ManyVariablesAccessor &MMVAPtr, NimArr<D, T> &nimArr){
	vector<SingleVariableAccessBase*> SMVA_Vec = MMVAPtr.getAccessVector();
	int nimCurrent = 0;
	int nimEnd = nimArr.size();
	int k = SMVA_Vec.size();
	int nextNumVals;
	SingleVariableAccessBase* curSingleAccess;
	for(int i = 0; i < k ; i++){
		curSingleAccess = SMVA_Vec[i];
		nextNumVals = (*curSingleAccess).getLength();
		if(nextNumVals + nimCurrent > nimEnd){
			PRINTF("Warning: in nimArr_2_ManyModelAccess, accessor larger than NimArr!\n");
			break;
			}
		SingleModelAccess_2_nimArr<D, T>(curSingleAccess, nimArr, nimCurrent);
		nimCurrent = nimCurrent + nextNumVals;
				
	}
	if(nimCurrent != nimEnd)
		PRINTF("Warning: after completing nimArr_2_ManyModelAccess, nimCurrent != nimEnd. Perhaps the NimArr was longer than the accessor?\n");
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


// new MapAccessor versions
void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to) { //map version

  vector<SingleVariableMapAccessBase *> fromAccessors = from.getMapAccessVector();
  vector<SingleVariableMapAccessBase *> toAccessors = to.getMapAccessVector();

#ifdef __NIMBLE_DEBUG_ACCESSORS
  PRINTF("Entering nimCopy\n");
  from.check();
  to.check();
#endif
  
  if(fromAccessors.size() != toAccessors.size()) {
    std::cout<<"Error in nimCopy: from and to access vectors have sizes "<<fromAccessors.size() << " and " << toAccessors.size() << "\n";
  }

  vector<SingleVariableMapAccessBase *>::iterator iFrom, iTo, iFromEnd;
  iFrom = fromAccessors.begin();
  iFromEnd = fromAccessors.end();
  iTo =  toAccessors.begin();
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
};


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
      cout<<"Error in nimCopyOne: unknown type for destination\n";
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
      cout<<"Error in nimCopyOne: unknown type for destination\n";
    }
    break;
  default:
    cout<<"Error in nimCopyOne: unknown type for source\n";
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
	cout<<"Error in nimCopyOne: unknown type for destination\n";
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
	cout<<"Error in nimCopyOne: unknown type for destination\n";
      }
      break;
    default:
      cout<<"Error in nimCopyOne: unknown type for destination\n";
    }
  }
}

void singletonCopyCheck(NimArrType *NAT, int offset) {
  nimType NATtype = NAT->getNimType();
  int NATsize;
  switch(NATtype) {
  case INT:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->getVptr()->size();
    break;
  case DOUBLE:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->getVptr()->size();
    break;
  default:
    PRINTF("Error with a NimArrType type\n");
    break;
  }
  if(offset < 0 || offset >= NATsize) PRINTF("Run-time error: bad singleton offset\n");
}

void dynamicMapCopyCheck(NimArrType *NAT, int offset, vector<int> &strides, vector<int> &sizes) {
  nimType NATtype = NAT->getNimType();
  int NATsize;
  switch(NATtype) {
  case INT:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->getVptr()->size();
    break;
  case DOUBLE:
    NATsize = static_cast<NimArrBase<int>*>(NAT)->getVptr()->size();
    break;
  default:
    PRINTF("Error with a NimArrType type\n");
    break;
  }
  if(offset < 0 || offset >= NATsize) PRINTF("Run-time error: bad offset\n");
  int lastOffset = offset;
  for(int i = 0; i < strides.size(); ++i) {
    lastOffset += sizes[i] * strides[i];
  }
  if(lastOffset < 0 || lastOffset >= NATsize) PRINTF("Run-time error: bad lastOffset\n");
}

// old versions:
void nimCopy(ManyVariablesAccessorBase &from, ManyVariablesAccessorBase &to) {

  vector<SingleVariableAccessBase *> fromAccessors = from.getAccessVector();
  vector<SingleVariableAccessBase *> toAccessors = to.getAccessVector();

  if(fromAccessors.size() != toAccessors.size()) {
    std::cout<<"Error in nimCopy: from and to access vectors have sizes "<<fromAccessors.size() << " and " << toAccessors.size() << "\n";
  }

  vector<SingleVariableAccessBase *>::iterator iFrom, iTo, iFromEnd;
  iFrom = fromAccessors.begin();
  iFromEnd = fromAccessors.end();
  iTo =  toAccessors.begin();
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
};




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
      cout<<"Error in nimCopyOne: unknown type for destination\n";
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
      cout<<"Error in nimCopyOne: unknown type for destination\n";
    }
    break;
  default:
    cout<<"Error in nimCopyOne: unknown type for source\n";
  }
}





SEXP makeSingleVariableAccessor(SEXP rModelPtr, SEXP elementName,  SEXP beginIndex, SEXP endIndex){
	SingleVariableAccess* sVAPtr = NULL;
	void* vPtr = R_ExternalPtrAddr(rModelPtr);
	if(vPtr != NULL){
		ModelBase* mPtr = static_cast<ModelBase*>(vPtr);
		string eName = STRSEXP_2_string(elementName, 0);
		sVAPtr = new SingleVariableAccess;
		(*sVAPtr).flatIndexStart = INTEGER(beginIndex)[0] - 1;
		(*sVAPtr).flatIndexEnd = INTEGER(endIndex)[0] - 1;
		(*sVAPtr).length = (*sVAPtr).flatIndexEnd - (*sVAPtr).flatIndexStart + 1;
		(*sVAPtr).ppVar = static_cast<NimArrType**> (mPtr->getObjectPtr(eName) );
		}
	
	SEXP rPtr;
	PROTECT(rPtr = R_MakeExternalPtr(sVAPtr, R_NilValue, R_NilValue) );
	R_RegisterCFinalizerEx(rPtr, &dontDeleteFinalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
}


SEXP makeSingleModelValuesAccessor(SEXP rModelValuesPtr, SEXP elementName,  SEXP curRow, SEXP beginIndex, SEXP endIndex){
	SingleModelValuesAccess* sMVAPtr = NULL;
	void* vPtr = R_ExternalPtrAddr(rModelValuesPtr);
	if(vPtr != NULL){
		Values* MVPtr = static_cast<Values*>(vPtr);
		int cRow = INTEGER(curRow)[0] - 1;
		string eName = STRSEXP_2_string(elementName, 0);
		sMVAPtr = new SingleModelValuesAccess;
		(*sMVAPtr).flatIndexStart = INTEGER(beginIndex)[0] - 1;
		(*sMVAPtr).flatIndexEnd = INTEGER(endIndex)[0] - 1;
		(*sMVAPtr).length = (*sMVAPtr).flatIndexEnd - (*sMVAPtr).flatIndexStart + 1;
		(*sMVAPtr).currentRow = cRow;
		(*sMVAPtr).pVVar = static_cast<NimVecType*> (MVPtr->getObjectPtr(eName) );
		}
	
	SEXP rPtr;
	PROTECT(rPtr = R_MakeExternalPtr(sMVAPtr, R_NilValue, R_NilValue) );
	R_RegisterCFinalizerEx(rPtr, &dontDeleteFinalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
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

SEXP newNodeFxnVector(SEXP size){
	NodeVectorClass* nVPtr = new NodeVectorClass;
	int cSize = INTEGER(size)[0];
	(*nVPtr).nodeFunPtrs.resize(cSize);
	SEXP rPtr = R_MakeExternalPtr(nVPtr, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	R_RegisterCFinalizerEx(rPtr, &NodeVector_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
	}

SEXP resizeNodeFxnVector(SEXP nodeFxnVecPtr, SEXP size){
	int cSize = INTEGER(size)[0];
	NodeVectorClass* nodeVec = static_cast<NodeVectorClass*>(R_ExternalPtrAddr(nodeFxnVecPtr) ) ;
	(*nodeVec).nodeFunPtrs.resize(cSize);
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


/*
SEXP populateNumberedObject_withSingleVariableAccessors(SEXP modelPtr, SEXP varName, SEXP Sgids, SEXP SnumbObj){
	int len = LENGTH(Sgids);
	int* gids = INTEGER(Sgids);
	string vName = STRSEXP_2_string(varName, 0);
	ModelBase* cModPtr = static_cast<ModelBase*>(R_ExternalPtrAddr(modelPtr) );
	NimArrType** varPtr = static_cast<NimArrType**>(cModPtr->getObjectPtr(vName) );
	SingleVariableAccess* smva;
	NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumbObj));
	void* vPtr;
	for(int i = 0; i < len; i++){
		smva = cMakeSingleVariableAccessor(varPtr, i, i);
		vPtr = static_cast<void*>(smva);
		numObj->numberedObjects[gids[i] - 1] = vPtr;
	}
	return(R_NilValue);
}
*/



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



/*			This is no longer used
SEXP populateNodeFxnVector(SEXP nodeFxnVec, SEXP nodeNames, SEXP nodeEnv){
	int numNodes = LENGTH(nodeNames);
	SEXP thisPtr;
	PROTECT(thisPtr);
	SEXP thisName = PROTECT(allocVector(STRSXP, 1) );
	SEXP falseObj = PROTECT( ScalarLogical(FALSE) );
	SEXP indexObj = PROTECT( ScalarInteger(0) );
	for(int i = 0; i < numNodes; i++){
		INTEGER(indexObj)[0]++;
		thisPtr = getEnvVar_Sindex(nodeNames, nodeEnv, indexObj );
		if(TYPEOF(thisPtr) != EXTPTRSXP)
			error("trying to copy a non-existant pointer");
		thisPtr = addNodeFun(nodeFxnVec, thisPtr, falseObj, indexObj );
	}
	UNPROTECT(4);
	return(R_NilValue);
}
*/

SEXP populateNodeFxnVector_byGID(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj){
	int len = LENGTH(S_GIDs);
	int* gids = INTEGER(S_GIDs);
	int index;
	NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumberedObj));
	NodeVectorClass* nfv = static_cast<NodeVectorClass*>(R_ExternalPtrAddr(SnodeFxnVec) ) ;
	(*nfv).nodeFunPtrs.resize(len);
	for(int i = 0; i < len; i++){
		index = gids[i] - 1;
		(*nfv).nodeFunPtrs[i] = static_cast<nodeFun*>(numObj->getObjectPtr(index));
		}
	return(R_NilValue);
}

SEXP populateModelValuesAccessors_byGID(SEXP SmodelValuesAccessorVector, SEXP S_GIDs, SEXP SnumberedObj){
	int len = LENGTH(S_GIDs);
	int* gids = INTEGER(S_GIDs);
	int index;
	NumberedObjects* numObj = static_cast<NumberedObjects*>(R_ExternalPtrAddr(SnumberedObj));
	ManyModelValuesAccessor* accessVector = static_cast<ManyModelValuesAccessor*>(R_ExternalPtrAddr(SmodelValuesAccessorVector) );
	(*accessVector).varAccessors.resize(len);
	for(int i = 0; i < len; i++){
		index = gids[i] - 1;
		(*accessVector).varAccessors[i] = static_cast<SingleModelValuesAccess*>(numObj->getObjectPtr(index));
		}
	return(R_NilValue);
}

///NEW

SEXP populateValueMapAccessors(SEXP StargetPtr, SEXP SsourceList, SEXP SModelOrModelValuesPtr ) {
  // typeCode = 1 for model, 2 for modelValues
  ManyVariablesMapAccessorBase* valuesAccessor = static_cast<ManyVariablesMapAccessorBase*>(R_ExternalPtrAddr(StargetPtr) );
  int numAccessors = LENGTH(SsourceList);
  valuesAccessor->resize(numAccessors);
  vector<SingleVariableMapAccessBase *> *singleAccessors = &(valuesAccessor->getMapAccessVector());
  
  NamedObjects *sourceNamedObject = static_cast<NamedObjects*>(R_ExternalPtrAddr(SModelOrModelValuesPtr));
  SEXP SoneSource;
  string varName;
  for(int i = 0; i < numAccessors; ++i) {
    PROTECT(SoneSource = VECTOR_ELT(SsourceList, i));
    (*singleAccessors)[i]->getOffset() = SEXP_2_int(VECTOR_ELT(SoneSource, 0));
    (*singleAccessors)[i]->getSizes() = SEXP_2_vectorInt(VECTOR_ELT(SoneSource, 1));
    (*singleAccessors)[i]->calculateLength();
    (*singleAccessors)[i]->getStrides() = SEXP_2_vectorInt(VECTOR_ELT(SoneSource, 2));
    (*singleAccessors)[i]->getSingleton() = static_cast<bool>(SEXP_2_int(VECTOR_ELT(SoneSource, 4)));
    varName = STRSEXP_2_string(VECTOR_ELT(SoneSource, 3), 0);
    (*singleAccessors)[i]->setObject( sourceNamedObject->getObjectPtr(varName) );
    UNPROTECT(1);
  }
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


void cAddNodeFun(NodeVectorClass* nVPtr, nodeFun* nFPtr, bool addAtEnd, int index){
	int size = (*nVPtr).nodeFunPtrs.size();
	if(addAtEnd == TRUE) {	
		(*nVPtr).nodeFunPtrs.push_back(nFPtr);
		return;
	}
	if((index >= size) | (index < 0)){
		PRINTF("Invalid index passed to addNodeFun\n");
		return;
	}
	(*nVPtr).nodeFunPtrs[index] = nFPtr;
	return;	
}

void cRemoveNodeFun(NodeVectorClass* nVPtr, int index, bool removeAll){
	if(removeAll == TRUE){
		(*nVPtr).nodeFunPtrs.erase( (*nVPtr).nodeFunPtrs.begin(), (*nVPtr).nodeFunPtrs.end() );
		return;
	}
	if((index >= static_cast<signed int>((*nVPtr).nodeFunPtrs.size())) | (index < 0)){
		PRINTF("Warning: attempted to delete nodeFunction from nodeFunctionVector with invalid index\n");
		return;
	}
	(*nVPtr).nodeFunPtrs.erase( (*nVPtr).nodeFunPtrs.begin() + index);
	return; 
}

SEXP removeNodeFun(SEXP rPtr, SEXP index, SEXP removeAll){
	int cIndex = INTEGER(index)[0] - 1;
	bool cRemoveAll = LOGICAL(removeAll)[0];
	void* vPtr = R_ExternalPtrAddr(rPtr);
	if(vPtr == NULL){
		PRINTF("Warning: pointer to null passed to removeNodeFun\n");
		return(R_NilValue);
	} 
	NodeVectorClass* nVPtr = static_cast<NodeVectorClass*>(vPtr);
	cRemoveNodeFun(nVPtr, cIndex, cRemoveAll);
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

/*
SEXP setNodeModelPtr(SEXP nodeFxnPtr, SEXP modelElementPtr, SEXP nodeElementName){
	SEXP nodeElementPtr = getModelObjectPtr(nodeFxnPtr, nodeElementName);
	PROTECT(nodeElementPtr);
	NimArrType** nodeNimArrTypePtr = static_cast<NimArrType**>(R_ExternalPtrAddr(nodeElementPtr) ) ;
	NimArrType** cModelElementPtr = static_cast<NimArrType**>(R_ExternalPtrAddr(modelElementPtr)); 
	nodeNimArrTypePtr = cModelElementPtr;
	UNPROTECT(1);
	return(R_NilValue);
}
*/

SEXP addNodeFun(SEXP nVPtr, SEXP nFPtr, SEXP addAtEnd, SEXP index){
	void* vNVPtr = R_ExternalPtrAddr(nVPtr);
	void* vNFPtr = R_ExternalPtrAddr(nFPtr);
	int cIndex = - 1;
	bool cAddAtEnd = LOGICAL(addAtEnd)[0];
	if(cAddAtEnd == FALSE)
		cIndex = INTEGER(index)[0] - 1;
	NodeVectorClass* cNVPtr = static_cast<NodeVectorClass*>(vNVPtr);
	nodeFun* cNFPtr = static_cast<nodeFun*>(vNFPtr);
	cAddNodeFun(cNVPtr, cNFPtr, cAddAtEnd, cIndex);
	return(R_NilValue);
}


SEXP newManyVariableAccessor(SEXP size){
	ManyVariablesAccessor* nMVAPtr = new ManyVariablesAccessor;
	int cSize = INTEGER(size)[0];
	(*nMVAPtr).varAccessors.resize(cSize);
	SEXP rPtr = R_MakeExternalPtr(nMVAPtr, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	R_RegisterCFinalizerEx(rPtr, &ManyVariable_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
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


SEXP newManyModelValuesAccessor(SEXP size){
	ManyModelValuesAccessor* nMVAPtr = new ManyModelValuesAccessor;
	int cSize = INTEGER(size)[0];
	(*nMVAPtr).varAccessors.resize(cSize);
	SEXP rPtr = R_MakeExternalPtr(nMVAPtr, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	R_RegisterCFinalizerEx(rPtr, &ManyMV_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
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


void SingleModelValuesAccessor_NumberedObjects_Finalizer(SEXP Snp){
	SpecialNumberedObjects<SingleModelValuesAccess>* np 
	= static_cast<SpecialNumberedObjects<SingleModelValuesAccess>*>(R_ExternalPtrAddr(Snp));
	delete np;		
}

SEXP new_SingleModelValuesAccessor_NumberedObjects(){
	SpecialNumberedObjects<SingleModelValuesAccess>* np = new SpecialNumberedObjects<SingleModelValuesAccess>;
	SEXP rPtr = R_MakeExternalPtr(np, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	R_RegisterCFinalizerEx(rPtr, &SingleModelValuesAccessor_NumberedObjects_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
	}

void SingleVariableAccessBase_NumberedObjects_Finalizer(SEXP Snp){
	SpecialNumberedObjects<SingleVariableAccessBase>* np 
	= static_cast<SpecialNumberedObjects<SingleVariableAccessBase>*>(R_ExternalPtrAddr(Snp));
	delete np;	
}

SEXP new_SingleModelVariablesAccessor_NumberedObjects(){
	SpecialNumberedObjects<SingleVariableAccessBase>* np = new SpecialNumberedObjects<SingleVariableAccessBase>;
	SEXP rPtr = R_MakeExternalPtr(np, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	R_RegisterCFinalizerEx(rPtr, &SingleVariableAccessBase_NumberedObjects_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
	}

void  SingleMVA_Finalizer ( SEXP Sv ) {
	SingleModelValuesAccess* oldObj;
	oldObj = static_cast<SingleModelValuesAccess *>(R_ExternalPtrAddr(Sv));
	delete oldObj;
}


void  SingleVA_Finalizer ( SEXP Sv ) {
	SingleVariableAccess* oldObj;
	oldObj = static_cast<SingleVariableAccess *>(R_ExternalPtrAddr(Sv));
	delete oldObj;
}

void NodeVector_Finalizer( SEXP Sv) {
	NodeVectorClass* nVPtr = static_cast<NodeVectorClass *>(R_ExternalPtrAddr(Sv));
	delete nVPtr;
}

void ManyVariable_Finalizer(SEXP Sv){
	ManyVariablesAccessor* mVAPtr = static_cast<ManyVariablesAccessor*>(R_ExternalPtrAddr(Sv));
	delete mVAPtr;
}

void ManyMV_Finalizer(SEXP Sv){
	ManyModelValuesAccessor* mVPtr = static_cast<ManyModelValuesAccessor*>(R_ExternalPtrAddr(Sv));
	delete mVPtr;
}
