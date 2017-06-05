// NOTE: THESE FUNCTIONS WILL ASSUME MODELS ARE OF CLASS
// "ModelBase". THIS CLASS CAN BE FOUND IN "ModelClassUtils.h" 

// This file contains new components for accessing and copying from groups of node Functions or subsets of variables in models or modelValues
#ifndef __ACCESSORCLASSES
#define __ACCESSORCLASSES

#include <iostream>

#include "NimArrBase.h"
#include "NimArr.h"			
#include "ModelClassUtils.h"
#include "RcppNimbleUtils.h"
#include <Rinternals.h>
#include "R.h"


using std::cout;

#include "nodeFun.h" 

//#define __NIMBLE_DEBUG_ACCESSORS

/////////////////////////////////
// 1. NodeVectors:
/////////////////////////////////
class oneNodeUseInfo {
 public:
  nodeFun *nodeFunPtr;
  useInfoForIndexedNodeInfo useInfo;
  oneNodeUseInfo(nodeFun *nFP, int firstIndex) {
    nodeFunPtr = nFP;
    useInfo.indicesForIndexedNodeInfo.push_back(firstIndex);
  }
};

class NodeVectorClassNew {
 public:
  vector<oneNodeUseInfo> useInfoVec;
  vector<oneNodeUseInfo> &getUseInfoVec() {return(useInfoVec);}
};

///// Using NodeVectors:
// utilities for calling node functions from a vector of node pointers
// see .cpp file for definitions
double calculate(NodeVectorClassNew &nodes);
double calculate(NodeVectorClassNew &nodes, int iNodeFunction);
double calculateDiff(NodeVectorClassNew &nodes);
double calculateDiff(NodeVectorClassNew &nodes, int iNodeFunction);
double getLogProb(NodeVectorClassNew &nodes);
double getLogProb(NodeVectorClassNew &nodes, int iNodeFunction);
void simulate(NodeVectorClassNew &nodes);
void simulate(NodeVectorClassNew &nodes, int iNodeFunction);

// ideas on efficiency
/* ## could propagate const-ness through getParam_0D_double_block etc. - that helped */
/* ## could pull getParam implementations back to .h for inlining. */
/*  ## could make a const version of operator[] and operator() for NimArray's */
/*  ## the getParam_..._block is not needed. -skipping this didn't help   */
/*  ## could do this all on a new branch of newNodeFxn to be able to compare */
/*  ## instead of the extra argument with default 0, use overloaded versions - might have helped - didn't compare carefully */
/*  ## try making all calculate, simulate etc. const */
 

// for these, if there is use of iNodeFunction, it is generated directly from cppOutputGetParam
//getParam_0D
inline double getParam_0D_double(int paramID, const oneNodeUseInfo &useInfo) {
  return(useInfo.nodeFunPtr->getParam_0D_double_block(paramID, useInfo.useInfo));
}
inline double getParam_0D_double(int paramID, const oneNodeUseInfo &useInfo, int iNodeFunction) { 
  /* iNodeFunction sometimes needs to be generated in a call even if not needed */
  /* but we want to avoid compiled warnings about an unused argument */
  /* the following line of code tries to make the compiler think iNodeFunction will be used */
  if(iNodeFunction) paramID += 0;
  return(useInfo.nodeFunPtr->getParam_0D_double_block(paramID, useInfo.useInfo));
}
template<typename paramIDtype>
inline double getParam_0D_double(const paramIDtype &paramID, const oneNodeUseInfo &useInfo, int iNodeFunction) {
return(useInfo.nodeFunPtr->getParam_0D_double_block(paramID[iNodeFunction], useInfo.useInfo));
}

//getParam_1D
NimArr<1, double> getParam_1D_double(int paramID, const oneNodeUseInfo &useInfo, int iNodeFunction = 0);

template<class paramIDtype>
NimArr<1, double> getParam_1D_double(const paramIDtype &paramID, const oneNodeUseInfo &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getParam_1D_double_block(paramID[iNodeFunction], useInfo.useInfo));
}

//getParam_2D
NimArr<2, double> getParam_2D_double(int paramID, const oneNodeUseInfo &useInfo, int iNodeFunction = 0);
/* template<typename paramIDtype> */
template<class paramIDtype>
NimArr<2, double> getParam_2D_double(const paramIDtype &paramID, const oneNodeUseInfo &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getParam_2D_double_block(paramID[iNodeFunction], useInfo.useInfo));
}

// code for getBound copied over from getParam; only 0D currently used as we set bounds same for all values in multivariate nodes
//getBound_0D
inline double getBound_0D_double(int boundID, const oneNodeUseInfo &useInfo) {
  return(useInfo.nodeFunPtr->getBound_0D_double_block(boundID, useInfo.useInfo));
}
inline double getBound_0D_double(int boundID, const oneNodeUseInfo &useInfo, int iNodeFunction) { 
  /* iNodeFunction sometimes needs to be generated in a call even if not needed */
  /* but we want to avoid compiled warnings about an unused argument */
  /* the following line of code tries to make the compiler think iNodeFunction will be used */
  if(iNodeFunction) boundID += 0;
  return(useInfo.nodeFunPtr->getBound_0D_double_block(boundID, useInfo.useInfo));
}
template<typename boundIDtype>
inline double getBound_0D_double(const boundIDtype &boundID, const oneNodeUseInfo &useInfo, int iNodeFunction) {
return(useInfo.nodeFunPtr->getBound_0D_double_block(boundID[iNodeFunction], useInfo.useInfo));
}

//getBound_1D
NimArr<1, double> getBound_1D_double(int boundID, const oneNodeUseInfo &useInfo, int iNodeFunction = 0);


template<class boundIDtype>
NimArr<1, double> getBound_1D_double(const boundIDtype &boundID, const oneNodeUseInfo &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getBound_1D_double_block(boundID[iNodeFunction], useInfo.useInfo));
}


//getBound_2D
NimArr<2, double> getBound_2D_double(int boundID, const oneNodeUseInfo &useInfo, int iNodeFunction = 0);
/* template<typename boundIDtype> */
template<class boundIDtype>
NimArr<2, double> getBound_2D_double(const boundIDtype &boundID, const oneNodeUseInfo &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getBound_2D_double_block(boundID[iNodeFunction], useInfo.useInfo));
}


/////////////////////
// new version of variable accessors using maps (offset and strided windows into multivariate objects (NimArr<>s) )
/////////////////////

// single and many base classes
class SingleVariableMapAccessBase {
 public:
  int offset, length;
  bool singleton;
  vector<int> sizes, strides;
 SingleVariableMapAccessBase() : offset(0), singleton(false) {length = 0;} /* length(0) triggers the Rf_length macro mixup */
  virtual ~SingleVariableMapAccessBase();
  virtual NimArrType *getNimArrPtr()=0;
  void calculateLength() {
    length = 1;
    for(unsigned int i = 0 ; i < sizes.size(); ++i) {length *= sizes[i];}
  }
  int &getLength() {return(length);}
  int &getOffset() {return(offset);}
  vector<int> &getSizes() {return(sizes);}
  vector<int> &getStrides() {return(strides);}
  bool &getSingleton() {return(singleton);}
  virtual void setObject(void *object)=0;
  virtual void *getObject()=0;
};


class ManyVariablesMapAccessorBase {
 public:
  int totalLength;
 ManyVariablesMapAccessorBase() : totalLength(0) {}
  int &getTotalLength() {return(totalLength);}
  virtual vector<SingleVariableMapAccessBase *> &getMapAccessVector()=0;
  virtual void  setRow(int i) = 0;
  virtual ~ManyVariablesMapAccessorBase() {};
  virtual void resize(int n) = 0;
#ifdef __NIMBLE_DEBUG_ACCESSORS
  virtual void check(int i) = 0;
  virtual void check() = 0;
#endif
};

// single and many variable classes

class SingleVariableMapAccess : public SingleVariableMapAccessBase {
 public:
  NimArrType **ppVar; // I think we have to do some casting when populating this
  virtual NimArrType *getNimArrPtr() {return(*ppVar);}
  ~SingleVariableMapAccess() {};
  void setObject(void *object) {ppVar = static_cast<NimArrType**>(object);}
  void *getObject() {return(ppVar);}
};

class ManyVariablesMapAccessor : public ManyVariablesMapAccessorBase {
 public:
  vector<SingleVariableMapAccessBase *> varAccessors;
  virtual vector<SingleVariableMapAccessBase *> &getMapAccessVector() {return(varAccessors);}
  int &getNodeLength(int index) {return(varAccessors[index-1]->getLength());}
  ~ManyVariablesMapAccessor();
  void setRow(int i){PRINTF("Bug detected in code: attempting to setRow for model. Can only setRow for modelValues\n");}
  void resize(int n){
    // this is a destructive resize only intended to be used once at setup
    if(varAccessors.size() != 0) PRINTF("Run-time Warning: resizing a ManyVariablesMapAccessor that was not empty.\n");
    varAccessors.resize(n); for(int i = 0; i < n; ++i) varAccessors[i] = new SingleVariableMapAccess;
  }
#ifdef __NIMBLE_DEBUG_ACCESSORS
  void check();
  void check(int i);
#endif
};

// single and many modelValues classes
class SingleModelValuesMapAccess : public SingleVariableMapAccessBase {
 public:
  NimVecType *pVVar;
  int currentRow;
  virtual NimArrType *getNimArrPtr() {return(pVVar->getRowTypePtr(currentRow));}
 SingleModelValuesMapAccess() : currentRow(0) {};
  ~SingleModelValuesMapAccess() {};
  void setRow(int i) {currentRow = i;}
  int getRow() {return(currentRow);}
  void setObject(void *object) {pVVar = static_cast<NimVecType*>(object);}
  void *getObject() {return(pVVar);}
};

class ManyModelValuesMapAccessor : public ManyVariablesMapAccessorBase {
  public:
  int currentRow;
 ManyModelValuesMapAccessor() : currentRow(0) {}
  vector<SingleVariableMapAccessBase *> varAccessors;
  virtual vector<SingleVariableMapAccessBase *> &getMapAccessVector() {return(varAccessors);}
  virtual void setRow(int i);// see .cpp
  ~ManyModelValuesMapAccessor();
  void resize(int n){
    // this is a destructive resize only intended to be used once at setup
    if(varAccessors.size() != 0) PRINTF("Run-time Warning: resizing a ManyVariablesMapAccessor that was not empty.\n");
    varAccessors.resize(n); for(int i = 0; i < n; ++i) varAccessors[i] = new SingleModelValuesMapAccess;
    currentRow = 0; // constructor for the singles sets their currentRows to 0
  }
#ifdef __NIMBLE_DEBUG_ACCESSORS
  void check();
  void check(int i);
#endif
};

void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to);
void nimCopyOne(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to);
void nimCopy(ManyVariablesMapAccessorBase &from, int rowFrom, ManyVariablesMapAccessorBase &to);
void nimCopy(ManyVariablesMapAccessorBase &from, int rowFrom, ManyVariablesMapAccessorBase &to, int rowTo);
void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to, int rowTo);

void setValues(NimArrBase<double> &nimArr, ManyVariablesMapAccessor &MVA);
void setValues(NimArrBase<int> &nimArr, ManyVariablesMapAccessor &MVA);
void setValues(NimArrBase<double> &nimArr, ManyVariablesMapAccessor &MVA, int index);
void setValues(NimArrBase<int> &nimArr, ManyVariablesMapAccessor &MVA, int index);
void getValues(NimArr<1, double> &nimArr, ManyVariablesMapAccessor &MVA);
void getValues(NimArr<1, int> &nimArr, ManyVariablesMapAccessor &MVA);
void getValues(NimArr<1, double> &nimArr, ManyVariablesMapAccessor &MVA, int index);
void getValues(NimArr<1, int> &nimArr, ManyVariablesMapAccessor &MVA, int index);

///////
// copierClass
///////

class rowInfoClass {
 public:
  int rowFrom, rowTo;
 rowInfoClass() : rowFrom(0), rowTo(0) {}
};

class copierClass { // virtual base class
 public:
  virtual void copy(const rowInfoClass &rowInfo) const=0;
  virtual ~copierClass() {};
};

class copierVectorClass {
 public:
  rowInfoClass rowInfo;
  int &rowTo() {return rowInfo.rowTo;}
  int &rowFrom() {return rowInfo.rowFrom;}
  vector<copierClass*> copyVector;
  copierVectorClass(); // objects with setup content are not ready at instantiation.
  void setup(ManyVariablesMapAccessorBase *from, ManyVariablesMapAccessorBase *to, int isFromMV, int isToMV);// constructor that creates new entries
  ~copierVectorClass(); //destructor that deletes them
};


copierClass* makeOneCopyClass(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int isFromMV, int isToMV);

template<class Tfrom, class Tto>
  class singletonCopierClass_M2MV : public copierClass {
 public:
  int toOffset, fromOffset;
 singletonCopierClass_M2MV() : toOffset(0), fromOffset(0) {}
  singletonCopierClass_M2MV(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int notUsed1, int notUsed2) {
    fromNimArr =  static_cast<NimArrType**>(from->getObject());
    toVecNimArr = static_cast<NimVecType*>(to->getObject());
    toOffset = to->getOffset();
    fromOffset = from->getOffset();
  }
  NimArrType **fromNimArr; 
  NimVecType *toVecNimArr; 
  void copy(const rowInfoClass &rowInfo) const {
    (*static_cast<VecNimArrBase<Tto> *>(toVecNimArr)->getBasePtr(rowInfo.rowTo)->getVptr())[toOffset] = (*static_cast<NimArrBase<Tfrom> *>(*fromNimArr)->getVptr())[fromOffset];
  }
  ~singletonCopierClass_M2MV() {};
};

template<class Tfrom, class Tto>
  class singletonCopierClass_M2M : public copierClass {
 public:
  int toOffset, fromOffset;
 singletonCopierClass_M2M() : toOffset(0), fromOffset(0) {}
  singletonCopierClass_M2M(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int notUsed1, int notUsed2) {
    fromNimArr =  static_cast<NimArrType**>(from->getObject());
    toNimArr = static_cast<NimArrType**>(to->getObject());
    toOffset = to->getOffset();
    fromOffset = from->getOffset();
  }
  NimArrType **fromNimArr; 
  NimArrType **toNimArr;
  void copy(const rowInfoClass &rowInfo) const {
    (*static_cast<NimArrBase<Tto> *>(*toNimArr)->getVptr())[toOffset] = (*static_cast<NimArrBase<Tfrom> *>(*fromNimArr)->getVptr())[fromOffset];
  }
  ~singletonCopierClass_M2M() {};
};

template<class Tfrom, class Tto>
  class singletonCopierClass_MV2M : public copierClass {
 public:
  int toOffset, fromOffset;
 singletonCopierClass_MV2M() : toOffset(0), fromOffset(0) {}
  singletonCopierClass_MV2M(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int notUsed1, int notUsed2) {
    fromVecNimArr =  static_cast<NimVecType*>(from->getObject());
    toNimArr = static_cast<NimArrType**>(to->getObject());
    toOffset = to->getOffset();
    fromOffset = from->getOffset();
  }
  NimVecType *fromVecNimArr; 
  NimArrType **toNimArr;
  void copy(const rowInfoClass &rowInfo) const {
    (*static_cast<NimArrBase<Tto> *>(*toNimArr)->getVptr())[toOffset] = (*static_cast<VecNimArrBase<Tfrom> *>(fromVecNimArr)->getBasePtr(rowInfo.rowFrom)->getVptr())[fromOffset];
  }
  ~singletonCopierClass_MV2M() {};
};

template<class Tfrom, class Tto>
  class singletonCopierClass_MV2MV : public copierClass {
 public:
  int toOffset, fromOffset;
 singletonCopierClass_MV2MV() : toOffset(0), fromOffset(0) {}
  singletonCopierClass_MV2MV(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int notUsed1, int notUsed2) {
    fromVecNimArr =  static_cast<NimVecType*>(from->getObject());
    toVecNimArr = static_cast<NimVecType*>(to->getObject());
    toOffset = to->getOffset();
    fromOffset = from->getOffset();
  }
  NimVecType *fromVecNimArr; 
  NimVecType *toVecNimArr;
  void copy(const rowInfoClass &rowInfo) const {
    (*static_cast<VecNimArrBase<Tto> *>(toVecNimArr)->getBasePtr(rowInfo.rowTo)->getVptr())[toOffset] = (*static_cast<VecNimArrBase<Tfrom> *>(fromVecNimArr)->getBasePtr(rowInfo.rowFrom)->getVptr())[fromOffset];
  }
  ~singletonCopierClass_MV2MV() {};
};

class blockCopierClassBase : public copierClass {
 public:
  //  int toOffset, fromOffset;
  bool isFromMV, isToMV;
 blockCopierClassBase() : isFromMV(false), isToMV(false) {}
  // put common initializations here
  // for now leave strides and sizes dynamically accessed
  blockCopierClassBase(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int isFromMV1, int isToMV1) {
    //  toOffset = to->getOffset();
    //  fromOffset = from->getOffset();
    isFromMV = isFromMV1;
    isToMV = isToMV1;
  }
};

template<class Tfrom, class Tto, int mapDim>
  class blockCopierClass : public blockCopierClassBase {
 public:
  SingleVariableMapAccessBase *fromPtr, *toPtr;
 blockCopierClass(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int isFromMV1, int isToMV1)
   : blockCopierClassBase(from, to, isFromMV1, isToMV1)
  {
    fromPtr = from;
    toPtr = to;
  }
  void copy(const rowInfoClass &rowInfo) const {
    // set rows if needed here.
    if(isFromMV) {
      static_cast<SingleModelValuesMapAccess *>(fromPtr)->setRow(rowInfo.rowFrom);
    }
    if(isToMV) {
      static_cast<SingleModelValuesMapAccess *>(toPtr)->setRow(rowInfo.rowTo);
    }
    dynamicMapCopyDim<Tfrom, Tto, mapDim>(toPtr->getNimArrPtr(), toPtr->getOffset(), toPtr->getStrides(), toPtr->getSizes(), fromPtr->getNimArrPtr(), fromPtr->getOffset(), fromPtr->getStrides(), fromPtr->getSizes() );
  }
};

class copierClassBuilderClass {
 public:
  virtual copierClass *build(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int isFromMV, int isToMV)=0;
};

template<class TDD, class TDI, class TID, class TII>
class copierClassBuilderCase : public copierClassBuilderClass {
 public:
  copierClassBuilderCase() {}
  copierClass *build(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to, int isFromMV, int isToMV) { // branch on types, singletons
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
	  return new TDD(from, to, isFromMV, isToMV);
	  break;
	case INT:
	  return new TDI(from, to, isFromMV, isToMV);
	  break;
	default:
	  NIMERROR("problem in copierClassBuilderCase");
	  return 0;
	  break;
	};
      case INT:
	switch(toType) {
	case DOUBLE:
	  return new TID(from, to, isFromMV, isToMV);
	  break;
	case INT:
	  return new TII(from, to, isFromMV, isToMV);
	  break;
	default:
	  NIMERROR("problem in copierClassBuilderCase");
	  return 0;
	  break;
	};
      default:
	NIMERROR("problem in copierClassBuilderCase");
	return 0;
	break;
      }
  }
};

extern copierClassBuilderCase< singletonCopierClass_M2M<double, double>, singletonCopierClass_M2M<double, int>, singletonCopierClass_M2M<int, int>, singletonCopierClass_M2M<int, double> > globalCopierBuilder_singleton_M2M;

extern copierClassBuilderCase< singletonCopierClass_M2MV<double, double>, singletonCopierClass_M2MV<double, int>, singletonCopierClass_M2MV<int, int>, singletonCopierClass_M2MV<int, double> > globalCopierBuilder_singleton_M2MV;

extern copierClassBuilderCase< singletonCopierClass_MV2M<double, double>, singletonCopierClass_MV2M<double, int>, singletonCopierClass_MV2M<int, int>, singletonCopierClass_MV2M<int, double> > globalCopierBuilder_singleton_MV2M;

extern copierClassBuilderCase< singletonCopierClass_MV2MV<double, double>, singletonCopierClass_MV2MV<double, int>, singletonCopierClass_MV2MV<int, int>, singletonCopierClass_MV2MV<int, double> > globalCopierBuilder_singleton_MV2MV;


extern copierClassBuilderCase< blockCopierClass<double, double, 1>, blockCopierClass<double, int, 1>, blockCopierClass<int, int, 1>, blockCopierClass<int, double, 1> > globalCopierClassBuilderBlock1;

extern copierClassBuilderCase< blockCopierClass<double, double, 2>, blockCopierClass<double, int, 2>, blockCopierClass<int, int, 2>, blockCopierClass<int, double, 2> > globalCopierClassBuilderBlock2;

extern copierClassBuilderCase< blockCopierClass<double, double, 3>, blockCopierClass<double, int, 3>, blockCopierClass<int, int, 3>, blockCopierClass<int, double, 3> > globalCopierClassBuilderBlock3;

extern copierClassBuilderCase< blockCopierClass<double, double, 4>, blockCopierClass<double, int, 4>, blockCopierClass<int, int, 4>, blockCopierClass<int, double, 4> > globalCopierClassBuilderBlock4;

/////////////////////////////////
// nimCopy function
/////////////////////////////////
void nimCopy(const copierVectorClass &copiers);
void nimCopy(copierVectorClass &copiers, int rowFrom);
void nimCopy(copierVectorClass &copiers, int rowFrom, int rowTo);
void nimCopy(copierVectorClass &copiers, int rowFrom, int rowTo, int unused);
	
void dynamicMapCopyCheck(NimArrType *NAT, int offset, vector<int> &strides, vector<int> &sizes);
void singletonCopyCheck(NimArrType *NAT, int offset);

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
      break;
    case INT:
      SMA_NimArrPtrI = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
      dynamicMapCopyDimToFlat<int, T>(&nimArr, nimBegin, nimStride, SMA_NimArrPtrI, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes());
      break;
    default:
      PRINTF("Copying type for SingleModelAccess_2_nimArr not supported\n");
      break;
    }
  }
}

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
      break;
    case INT:
      SMA_NimArrPtrI = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
      dynamicMapCopyFlatToDim<T, int>(SMA_NimArrPtrI, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes(), &nimArr, nimBegin, nimStride);
      break;
    default:
      PRINTF("Copying type for nimArr_2_SingleModelAccess not supported\n");
      break;
    }
  }
}

extern "C" {
  SEXP resizeManyModelVarAccessor(SEXP manyModelVarPtr, SEXP size);
  SEXP resizeManyModelValuesAccessor(SEXP manyModelValuesPtr, SEXP size);	 
  SEXP manualSetNRows(SEXP Sextptr, SEXP nRows);

  SEXP getVarAndIndices(SEXP Sstring);
  SEXP varAndIndices2mapParts(SEXP SvarAndIndicesExtPtr, SEXP Ssizes, SEXP SnDim);
  SEXP var2mapParts(SEXP Sinput, SEXP Ssizes, SEXP SnDim);
  
  SEXP populateNodeFxnVectorNew_byDeclID(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS);
  SEXP populateIndexedNodeInfoTable(SEXP StablePtr, SEXP StableContents);
  SEXP populateValueMapAccessorsFromNodeNames(SEXP StargetPtr, SEXP SnodeNames, SEXP SsizesAndNdims, SEXP SModelOrModelValuesPtr );
  SEXP populateValueMapAccessors(SEXP StargetPtr, SEXP SsourceList, SEXP SModelOrModelValuesPtr );

  SEXP populateNumberedObject_withSingleModelValuesAccessors(SEXP mvPtr, SEXP varName, SEXP GIDs, SEXP curRow, SEXP SnumbObj);

  SEXP populateCopierVector(SEXP ScopierVector, SEXP SfromPtr, SEXP StoPtr, SEXP SintIsFromMV, SEXP SintIsToMV);
  
  SEXP populateNumberedObject_withSingleModelVariablesAccessors(SEXP modelPtr, SEXP varName, SEXP sGIDS, SEXP SvalidIndices, SEXP SnumbObj);
  SEXP populateModelVariablesAccessors_byGID(SEXP SmodelVariableAccessorVector, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_LP_GIDs, SEXP S_LP_numberedObj);
}
void NodeVector_Finalizer( SEXP Sv);
void ManyVariable_Finalizer(SEXP Sv);
void ManyMV_Finalizer(SEXP Sv);

#endif
