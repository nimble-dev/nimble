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
//#include <nimble/nimbleCppAD.h>

using std::cout;

#include "nodeFun.h" 

//#define __NIMBLE_DEBUG_ACCESSORS

/////////////////////////////////
// 1. NodeVectors:
/////////////////////////////////

// This is a single instruction applying a nodeFun to a single node represented by operand.
class NodeInstruction {
 public:
  nodeFun *nodeFunPtr;
  int operand;  // An index into a nodeFun's indexedNodeInfoTable field.
  NodeInstruction(nodeFun *nFP, int op) : nodeFunPtr(nFP), operand(op) {}
};

// This is a collection of instructions denoting a sort of "program".
class NodeVectorClassNew {
 public:
  vector<NodeInstruction> instructions;
  vector<NodeInstruction> &getInstructions() { return instructions; }
};

#define _DERIVS_FULLTAPE
#ifdef _DERIVS_FULLTAPE


#else

// This is a collection of instructions denoting a sort of "program".
/* class NodeVectorClassNew_derivs : public NodeVectorClassNew { */
/* public: */
/*   vector<vector<NimArr<1, int> > > parentIndicesList; */
/*   vector<vector<vector<NimArr<1, int> > > > topLevelWrtDeps; */
/*   NimArr<1, int> stochNodeIndicators; */
/*   NimArr<1, int> isAddedScalarNode; */
/* 	NimArr<1, int> calcNodeIndicators; */
/* 	vector<NimArr<1, double> > cppWrtArgIndices; */
/* 	NimArr<1, int> wrtLineNums; */
/* 	NimArr<1, int> nodeLengths; */
/* 	vector<NimArr<1, int> > wrtToIndices; */
/* 	vector<NimArr<1, int> > wrtFromIndices; */
/* 	vector<NimArr<1, int> > wrtLineIndices; */
/*   vector<NimArr<1, int> > lineWrtArgSizeInfo; */
/*   vector<NimArr<1, int> > allNeededWRTCopyVars; */
/*   vector<NimArr<2, double> > thisAddedNodeJacobianList; */
/* 	int totalOutWrtSize; */
/* 	int totalWrtSize; */
/*   NimArr<1, int> cumulativeWrtLineNums; */
/*   NimArr<1, int> wrtLineSize; */
/* 	const vector<NodeInstruction> &getConstInstructions() const {  */
/*     return instructions; } */

/* 	void populateDerivsInfo(SEXP SderivsInfo) { */
/* 		SEXP S_pxData; */
/* 		SEXP S_parentInds; */
/*     SEXP S_thisList; */
/*     SEXP S_thisListI; */
/* 		SEXP S_stochNodeIndicators; */
/* 		SEXP S_calcNodeIndicators; */
/* 		SEXP S_cppWrtArgIndices; */
/* 		SEXP S_wrtLineNums; */
/* 		SEXP S_wrtToIndices; */
/* 		SEXP S_wrtFromIndices; */
/* 		SEXP S_wrtLineIndices; */
/* 		SEXP S_lineWrtArgSizeInfo; */
/*     SEXP S_nodeLengths; */
/*     SEXP S_topLevelWrtDeps; */
/*     SEXP S_allNeededWRTCopyVars; */
/*     SEXP S_isAddedScalarNode; */
/*     SEXP S_thisAddedNodeJacobianList; */
/* 		int numNodes; */
/*     int numNodesI; */
/* 		PROTECT(S_pxData = Rf_allocVector(STRSXP, 1)); */
/* 		SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData")); */
/* 		PROTECT(S_parentInds = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												Rf_install("parentIndicesList"))); */
/* 		numNodes = Rf_length(S_parentInds); */
/*     parentIndicesList.resize(numNodes); */
/* 		for(int i = 0; i < numNodes; i++){ */
/* 			PROTECT(S_thisList =  VECTOR_ELT(S_parentInds, i)); */
/* 			SEXP_list_2_NimArr_int_vec(S_thisList, parentIndicesList[i]); */
/*     } */

/*     PROTECT(S_thisAddedNodeJacobianList = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/*                         Rf_install("thisAddedNodeJacobianList"))); */
/*     SEXP_list_2_NimArr_double_vec(S_thisAddedNodeJacobianList, thisAddedNodeJacobianList); */

/*     PROTECT(S_topLevelWrtDeps = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												Rf_install("topLevelWrtDeps"))); */
/*     topLevelWrtDeps.resize(numNodes); */
/*     int sumNodesI = 0; */
/* 		for(int i = 0; i < numNodes; i++){ */
/*       PROTECT(S_thisList =  VECTOR_ELT(S_topLevelWrtDeps, i)); */
/*       numNodesI  = Rf_length(S_thisList); */
/*       sumNodesI += numNodesI; */
/*       topLevelWrtDeps[i].resize(numNodesI); */
/* 		  for(int j = 0; j < numNodesI; j++){ */
/* 			  PROTECT(S_thisListI =  VECTOR_ELT(S_thisList, j)); */
/* 		  	SEXP_list_2_NimArr_int_vec(S_thisListI, topLevelWrtDeps[i][j]); */
/*       } */
/*     } */

/*     PROTECT(S_isAddedScalarNode = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 														  Rf_install("isAddedScalarNode"))); */
/* 		SEXP_2_NimArr(S_isAddedScalarNode, isAddedScalarNode); */
/* 		PROTECT(S_stochNodeIndicators = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 														  Rf_install("stochNodeIndicators"))); */
/* 		SEXP_2_NimArr(S_stochNodeIndicators, stochNodeIndicators); */
/* 		PROTECT(S_calcNodeIndicators = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("calcNodeIndicators"))); */
/*   		SEXP_2_NimArr(S_calcNodeIndicators, calcNodeIndicators); */
/* 		PROTECT(S_wrtLineNums = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("wrtLineNums"))); */
/*   		SEXP_2_NimArr(S_wrtLineNums, wrtLineNums); */
/* 		PROTECT(S_nodeLengths = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("nodeLengths"))); */
/*   		SEXP_2_NimArr(S_nodeLengths, nodeLengths); */
/* 		PROTECT(S_cppWrtArgIndices = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("cppWrtArgIndices"))); */
/* 		SEXP_list_2_NimArr_double_vec(S_cppWrtArgIndices, cppWrtArgIndices); */
/* 		PROTECT(S_wrtToIndices = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("wrtToIndices"))); */
/* 		SEXP_list_2_NimArr_int_vec(S_wrtToIndices, wrtToIndices); */
/* 		PROTECT(S_wrtFromIndices = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("wrtFromIndices"))); */
/* 		SEXP_list_2_NimArr_int_vec(S_wrtFromIndices, wrtFromIndices); */
/* 		PROTECT(S_wrtLineIndices = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("wrtLineIndices"))); */
/* 		SEXP_list_2_NimArr_int_vec(S_wrtLineIndices, wrtLineIndices); */
/* 		PROTECT(S_lineWrtArgSizeInfo = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("lineWrtArgSizeInfo"))); */
/* 		SEXP_list_2_NimArr_int_vec(S_lineWrtArgSizeInfo, lineWrtArgSizeInfo); */
/* 		PROTECT(S_allNeededWRTCopyVars = Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, S_pxData)), */
/* 												         Rf_install("allNeededWRTCopyVars"))); */
/* 		SEXP_list_2_NimArr_int_vec(S_allNeededWRTCopyVars, allNeededWRTCopyVars); */

/* 		UNPROTECT(29 + 2*numNodes + sumNodesI); */
		
/* 		totalOutWrtSize = 0; */
/* 		for(int i = 0; i < length(wrtToIndices); i++){ */
/* 			totalOutWrtSize += wrtToIndices[i].dimSize(0); */
/* 		} */
		
/* 		cumulativeWrtLineNums.initialize(-1, 1, numNodes); */
/* 		wrtLineSize.setSize(wrtLineNums.dimSize(0)); */
/* 		totalWrtSize = 0; */
/* 		for(int i = 0; i < wrtLineNums.dimSize(0); i++){ */
/* 			cumulativeWrtLineNums[wrtLineNums[i] - 1] = i; */
/* 			wrtLineSize[i] = nodeLengths[wrtLineNums[i] -1]; */
/* 			totalWrtSize += wrtLineSize[i]; */
/* 		} */
/* 	} */
/* }; */
#endif


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
inline double getParam_0D_double(int paramID, const NodeInstruction &useInfo) {
  return(useInfo.nodeFunPtr->getParam_0D_double_block(paramID, useInfo.operand));
}
inline double getParam_0D_double(int paramID, const NodeInstruction &useInfo, int iNodeFunction) { 
  /* iNodeFunction sometimes needs to be generated in a call even if not needed */
  /* but we want to avoid compiled warnings about an unused argument */
  /* the following line of code tries to make the compiler think iNodeFunction will be used */
  if(iNodeFunction) paramID += 0;
  return(useInfo.nodeFunPtr->getParam_0D_double_block(paramID, useInfo.operand));
}
template<typename paramIDtype>
inline double getParam_0D_double(const paramIDtype &paramID, const NodeInstruction &useInfo, int iNodeFunction) {
return(useInfo.nodeFunPtr->getParam_0D_double_block(paramID[iNodeFunction], useInfo.operand));
}

//getParam_1D
NimArr<1, double> getParam_1D_double(int paramID, const NodeInstruction &useInfo, int iNodeFunction = 0);

template<class paramIDtype>
NimArr<1, double> getParam_1D_double(const paramIDtype &paramID, const NodeInstruction &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getParam_1D_double_block(paramID[iNodeFunction], useInfo.operand));
}

//getParam_2D
NimArr<2, double> getParam_2D_double(int paramID, const NodeInstruction &useInfo, int iNodeFunction = 0);
template<class paramIDtype>
NimArr<2, double> getParam_2D_double(const paramIDtype &paramID, const NodeInstruction &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getParam_2D_double_block(paramID[iNodeFunction], useInfo.operand));
}

// code for getBound copied over from getParam; only 0D currently used as we set bounds same for all values in multivariate nodes
//getBound_0D
inline double getBound_0D_double(int boundID, const NodeInstruction &useInfo) {
  return(useInfo.nodeFunPtr->getBound_0D_double_block(boundID, useInfo.operand));
}
inline double getBound_0D_double(int boundID, const NodeInstruction &useInfo, int iNodeFunction) { 
  /* iNodeFunction sometimes needs to be generated in a call even if not needed */
  /* but we want to avoid compiled warnings about an unused argument */
  /* the following line of code tries to make the compiler think iNodeFunction will be used */
  if(iNodeFunction) boundID += 0;
  return(useInfo.nodeFunPtr->getBound_0D_double_block(boundID, useInfo.operand));
}
template<typename boundIDtype>
inline double getBound_0D_double(const boundIDtype &boundID, const NodeInstruction &useInfo, int iNodeFunction) {
return(useInfo.nodeFunPtr->getBound_0D_double_block(boundID[iNodeFunction], useInfo.operand));
}

//getBound_1D
NimArr<1, double> getBound_1D_double(int boundID, const NodeInstruction &useInfo, int iNodeFunction = 0);


template<class boundIDtype>
NimArr<1, double> getBound_1D_double(const boundIDtype &boundID, const NodeInstruction &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getBound_1D_double_block(boundID[iNodeFunction], useInfo.operand));
}


//getBound_2D
NimArr<2, double> getBound_2D_double(int boundID, const NodeInstruction &useInfo, int iNodeFunction = 0);
template<class boundIDtype>
NimArr<2, double> getBound_2D_double(const boundIDtype &boundID, const NodeInstruction &useInfo, int iNodeFunction) {
  return(useInfo.nodeFunPtr->getBound_2D_double_block(boundID[iNodeFunction], useInfo.operand));
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
  virtual const vector<SingleVariableMapAccessBase *> &getConstMapAccessVector() const =0;
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
  virtual const vector<SingleVariableMapAccessBase *> &getConstMapAccessVector() const {return(varAccessors);}
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
  virtual const vector<SingleVariableMapAccessBase *> &getConstMapAccessVector() const {return(varAccessors);}
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
void setValues_AD_AD(NimArrBase< CppAD::AD<double> > &nimArr, ManyVariablesMapAccessor &MVA);
void setValues(NimArrBase<double> &nimArr, ManyVariablesMapAccessor &MVA, int index);
void setValues(NimArrBase<int> &nimArr, ManyVariablesMapAccessor &MVA, int index);
void getValues(NimArr<1, double> &nimArr, ManyVariablesMapAccessor &MVA);
void getValues(NimArr<1, int> &nimArr, ManyVariablesMapAccessor &MVA);
void getValues_AD_AD(NimArr<1, CppAD::AD<double> > &nimArr, ManyVariablesMapAccessor &MVA);
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
      case nimType::DOUBLE:
	switch(toType) {
	case nimType::DOUBLE:
	  return new TDD(from, to, isFromMV, isToMV);
	  break;
	case nimType::INT:
	  return new TDI(from, to, isFromMV, isToMV);
	  break;
	default:
	  NIMERROR("problem in copierClassBuilderCase");
	  return 0;
	  break;
	};
      case nimType::INT:
	switch(toType) {
	case nimType::DOUBLE:
	  return new TID(from, to, isFromMV, isToMV);
	  break;
	case nimType::INT:
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

//
void SingleModelAccess_2_nimArr_AD_AD(SingleVariableMapAccessBase* SMVAPtr,
				      NimArrBase< CppAD::AD<double> > &nimArr,
				      int nimBegin,
				      int nimStride);

template<class T>
void SingleModelAccess_2_nimArr(SingleVariableMapAccessBase* SMVAPtr, NimArrBase<T> &nimArr, int nimBegin, int nimStride){
  NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
  nimType SMA_Type = (*SMA_NimTypePtr).getNimType();
  NimArrBase<double>* SMA_NimArrPtrD;
  NimArrBase<int>* SMA_NimArrPtrI;
  if(SMVAPtr->getSingleton()) {
    switch(SMA_Type) {
    case nimType::DOUBLE:
      (*nimArr.getVptr())[nimBegin] = (*static_cast<NimArrBase<double> *>(SMA_NimTypePtr))[SMVAPtr->offset];
      break;
    case nimType::INT:
      (*nimArr.getVptr())[nimBegin] = (*static_cast<NimArrBase<int> *>(SMA_NimTypePtr))[SMVAPtr->offset];
      break;
    default:
      PRINTF("Copying type for SingleModelAccess_2_nimArr not supported\n");
      break;
    }
  } else {
    switch(SMA_Type) {
    case nimType::DOUBLE:
      SMA_NimArrPtrD = static_cast<NimArrBase<double>*>(SMA_NimTypePtr);
      dynamicMapCopyDimToFlat<double, T>(&nimArr, nimBegin, nimStride, SMA_NimArrPtrD, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes() );
      break;
    case nimType::INT:
      SMA_NimArrPtrI = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
      dynamicMapCopyDimToFlat<int, T>(&nimArr, nimBegin, nimStride, SMA_NimArrPtrI, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes());
      break;
    default:
      PRINTF("Copying type for SingleModelAccess_2_nimArr not supported\n");
      break;
    }
  }
}

void nimArr_2_SingleModelAccess_AD_AD(SingleVariableMapAccessBase* SMVAPtr, NimArrBase< CppAD::AD<double> > &nimArr, int nimBegin, int nimStride);

template<class T>
void nimArr_2_SingleModelAccess(SingleVariableMapAccessBase* SMVAPtr, NimArrBase<T> &nimArr, int nimBegin, int nimStride){

  NimArrType* SMA_NimTypePtr = (*SMVAPtr).getNimArrPtr();
  nimType SMA_Type = (*SMA_NimTypePtr).getNimType();
  NimArrBase<double>* SMA_NimArrPtrD;
  NimArrBase<int>* SMA_NimArrPtrI;

  if(SMVAPtr->getSingleton()) {
    switch(SMA_Type) {
    case nimType::DOUBLE:
      (*static_cast<NimArrBase<double> *>(SMA_NimTypePtr))[SMVAPtr->offset] = (*nimArr.getVptr())[nimBegin];
      break;
    case nimType::INT:
      (*static_cast<NimArrBase<double> *>(SMA_NimTypePtr))[SMVAPtr->offset] = (*nimArr.getVptr())[nimBegin];
      break;
    default:
      PRINTF("Copying type for nimArr_2_SingleModelAccess not supported\n");
      break;
    }
  } else {
    switch(SMA_Type) {
    case nimType::DOUBLE:
      SMA_NimArrPtrD = static_cast<NimArrBase<double>*>(SMA_NimTypePtr);
      dynamicMapCopyFlatToDim<T, double>(SMA_NimArrPtrD, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes(), &nimArr, nimBegin, nimStride);
      break;
    case nimType::INT:
      SMA_NimArrPtrI = static_cast<NimArrBase<int>*>(SMA_NimTypePtr);
      dynamicMapCopyFlatToDim<T, int>(SMA_NimArrPtrI, SMVAPtr->getOffset(), SMVAPtr->getStrides(), SMVAPtr->getSizes(), &nimArr, nimBegin, nimStride);
      break;
    default:
      PRINTF("Copying type for nimArr_2_SingleModelAccess not supported\n");
      break;
    }
  }
}

void populateNodeFxnVectorNew_internal_forDerivs(NodeVectorClassNew* nfv, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS, SEXP SderivInfo );
void populateNodeFxnVectorNew_copyFromRobject_forDerivs(void *nodeFxnVec_to, SEXP S_nodeFxnVec_from );
void populateNodeFxnVectorNew_internal(NodeVectorClassNew* nfv, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS );
void populateNodeFxnVectorNew_copyFromRobject(void *nodeFxnVec_to, SEXP S_nodeFxnVec_from );
void populateValueMapAccessorsFromNodeNames_internal(ManyVariablesMapAccessorBase* valuesAccessor,
						     SEXP SnodeNames,
						     SEXP SsizesAndNdims,
						     SEXP SModelOrModelValuesPtr);
void populateValueMapAccessorsFromNodeNames_copyFromRobject(void *VvaluesAccessor,
							    SEXP Sargs);
extern "C" {
  SEXP resizeManyModelVarAccessor(SEXP manyModelVarPtr, SEXP size);
  SEXP resizeManyModelValuesAccessor(SEXP manyModelValuesPtr, SEXP size);	 
  SEXP manualSetNRows(SEXP Sextptr, SEXP nRows);

  SEXP getVarAndIndices(SEXP Sstring);
  SEXP varAndIndices2mapParts(SEXP SvarAndIndicesExtPtr, SEXP Ssizes, SEXP SnDim);
  SEXP var2mapParts(SEXP Sinput, SEXP Ssizes, SEXP SnDim);
  
  SEXP populateNodeFxnVectorNew_byDeclID(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS);
  SEXP populateNodeFxnVectorNew_byDeclID_forDerivs(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_ROWINDS, SEXP SderivInfo);
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

indexedNodeInfo generateDummyIndexedNodeInfo();


#ifdef _DERIVS_FULLTAPE

class NodeVectorClassNew_derivs;

CppAD::AD<double> calculate_ADproxyModel(NodeVectorClassNew_derivs &nodes);

class NodeVectorClassNew_derivs : public NodeVectorClassNew {
 public:
  // Groups of nodes are labeled as:
  // wrt = "with respect to", the nodes with repsect to which we want derivatives
  // neededParents = "needed parent", addiitonal nodes that must be treated as CppAD independent variables
  //       These are nodes that are parents of anything in the calculations that are not themselves in the calculation anyway.
  //       Example z <- f(x, y), and wrt includes 'x'.  Then y is a needed parent node.
  // data = data nodes that must be included as indepedent variables
  // independentNodes = c(wrt, neededParents, data)
  // outputNodes = deterministic outputs and logProb vaules to copy back into the model
  //
  // To run the tape: (these steps are split into different functions for timing and control)
  // runTape_setIndependent:
  // 1. Copy all wrtNodes values from model -> independentVars;
  // runTape_runTape
  // 2. Run tape
  // runTape_unpackDependent
  // 3. Copy logProb value (scalar)
  // 4. Collect derivatives
  // 5. Arrange derivatives into return object
  //
  // To record the tape:
  // 1. Copy all constantNodes values from model -> model_AD
  // 2. Copy all wrtNodes values from model -> model_AD (may be redundant now but may be useful later).
  // 3. Copy all wrtNodes values from model -> independentVars
  // 4. independentVars should have an extra dummy element for getExtraInputs and setExtraInputs
  // 5. Copy all extraInputNodes values from model -> model_AD (ditto, may be redundant)
  // 6. Start taping
  // 7. Instantiate and call getExtraInputs atomic object
  // 8. Copy all extraInputNodes AD objects from extraInputResult -> model_AD
  // 9. Copy all wrtNodes AD objects from independentVars -> model_AD.
  // 10. Call calculate()
  // 11. Copy logProb to dependentVars[0]
  // 12. Call setModelOutputs to put AD modelOutputs into model
  // 13. Finish taping
  // 14. Call tape->optimize() 
  ManyVariablesMapAccessor model_wrt_accessor;
  ManyVariablesMapAccessor& get_model_wrt_accessor() {return model_wrt_accessor;}
  ManyVariablesMapAccessor model_AD_wrt_accessor;
  ManyVariablesMapAccessor model_extraInput_accessor;
  ManyVariablesMapAccessor model_AD_extraInput_accessor;
  ManyVariablesMapAccessor model_modelOutput_accessor;
  ManyVariablesMapAccessor model_AD_modelOutput_accessor;
  ManyVariablesMapAccessor model_constant_accessor;
  ManyVariablesMapAccessor model_AD_constant_accessor;
  CppAD::ADFun< double > ADtape;
  atomic_extraInputObject *extraInputObject;
  atomic_extraOutputObject *extraOutputObject;
  bool tapeRecorded_;
 NodeVectorClassNew_derivs() : tapeRecorded_(false) {}
  ~NodeVectorClassNew_derivs() {
    if(instructions.size() == 0) return;
    nodeFun* nodeFunInModelDLL = instructions[0].nodeFunPtr;
    if(extraInputObject) nodeFunInModelDLL->delete_extraInputObject(*this);
    if(extraOutputObject) nodeFunInModelDLL->delete_extraOutputObject(*this);
  }
  bool tapeRecorded() {return(tapeRecorded_);}
  void recordTape() {
    if(instructions.size() == 0) {
      printf("No nodes for calculation\n");
      return;
    }
    nodeFun* nodeFunInModelDLL = instructions[0].nodeFunPtr;
    nodeFunInModelDLL->recordTape(*this);
    tapeRecorded_ = true;
  }
  void runTape_setIndependent(std::vector<double> &independentVars) {
    //   std::cout<<"runTape_setInd"<<std::endl;
    // 1. Copy all independentNodes from model -> independentVars;
    int length_wrt = model_wrt_accessor.getTotalLength();
    int length_independent = length_wrt + 1; // + 1 for a dummy for extraInputNodes
    // std::cout<<"length_independent ="<<length_independent<<std::endl;
    //std::vector< double > independentVars(length_independent);

    independentVars.resize(length_independent);
    
    NimArr<1, double > NimArrVars;
    NimArrVars.setSize(length_wrt);
    getValues(NimArrVars, model_wrt_accessor);

    std::copy(NimArrVars.getPtr(),
	      NimArrVars.getPtr() + length_wrt,
	      independentVars.begin());

    //std::cout<<"done runTape_setInd"<<std::endl;
    // 2. Run tape
    //    std::vector<double> dependentVars;
    //dependentVars = ADtape.Forward(0, independentVars);
    //return(dependentVars[0]);
    // 3. Copy logProb value (scalar)
    // 4. Copy outputNodes values from dependentVars -> model
    // 5. Collect derivatives
    // 6. Sift derivatives into return object
  }
  void runTape_runTape(std::vector<double> &independentVars,
		       std::vector<double> &dependentVars,
		       const NimArr<1, double> &derivOrders,
		       nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
    nodeFun* nodeFunInModelDLL = instructions[0].nodeFunPtr;
    //std::cout<<"runTape_runTape"<<std::endl;
    nodeFunInModelDLL->runTape(ADtape, independentVars, dependentVars,
			       derivOrders, ansList); //ADtape.Forward(0, independentVars);
    // std::cout<<"runTape_runTape before reverse "<<q<<std::endl;
    //    cppad_derivOut = ADtape.Reverse(1, w);
    //std::cout<<"done runTape_runTape"<<std::endl;
  }
  double runTape_unpackDependent(std::vector<double> &dependentVars) {
    //  std::cout<<"runTape_unpackDep"<<std::endl;
    return(dependentVars[0]);
  }
  /* double runTape( const NimArr<1, double> &derivOrders) { */
  /*   // 1. Copy all independentNodes from model -> independentVars; */
  /*   int length_independent = model_independent_accessor.getTotalLength(); */
  /*   std::vector< double > independentVars(length_independent); */
  /*   NimArr<1, double > NimArrIndependentVars; */
  /*   NimArrIndependentVars.setSize(length_independent); */
  /*   getValues(NimArrIndependentVars, model_independent_accessor); */
  /*   std::copy(NimArrIndependentVars.getPtr(), */
  /* 	      NimArrIndependentVars.getPtr() + length_independent, */
  /* 	      independentVars.begin()); */
  /*   // 2. Run tape */
  /*   std::vector<double> dependentVars; */
  /*   dependentVars = ADtape.Forward(0, independentVars); */
  /*   return(dependentVars[0]); */
  /*   // 3. Copy logProb value (scalar) */
  /*   // 4. Copy outputNodes values from dependentVars -> model */
  /*   // 5. Collect derivatives */
  /*   // 6. Sift derivatives into return object */
  /* } */
  void populateDerivsInfo(SEXP SderivsInfo) {
    SEXP SpxData;
    SEXP Smodel, SCobjInt, SbasePtr, SADptrs, SbasePtrAD;

    PROTECT(SpxData = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(SpxData, 0, Rf_mkChar(".xData"));
    
    //Smodel <- SderivsInfo$model
    PROTECT(Smodel =
	    Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, SpxData)),
			      Rf_install("model")));
    //SCobjInt <- Smodel$CobjectInterface
    PROTECT(SCobjInt = 
	    Rf_findVarInFrame(PROTECT(GET_SLOT(Smodel, SpxData)),
			      Rf_install("CobjectInterface")));
    // SbasePtr <- SCobjInt$.basePtr
    PROTECT(SbasePtr = 
	    Rf_findVarInFrame(PROTECT(GET_SLOT(SCobjInt, SpxData)),
			      Rf_install(".basePtr")));
    // SADptrs <- SCobjInt$.ADptrs
    PROTECT(SADptrs = 
	    Rf_findVarInFrame(PROTECT(GET_SLOT(SCobjInt, SpxData)),
			      Rf_install(".ADptrs")));
    // SbasePtrAD <- SADptrs[[".ADptrs"]]
    PROTECT(SbasePtrAD = 
	    VECTOR_ELT(SADptrs, 0));


    //Swrt <- SderivsInfo$wrtMapInfo
    SEXP Swrt;
    SEXP SwrtNodeNames, SwrtSizesAndNdims;
    PROTECT(Swrt =
	    Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, SpxData)),
			      Rf_install("wrtMapInfo")));
    // SwrtNodeNames = Swrt[[1]]
    PROTECT(SwrtNodeNames = VECTOR_ELT(Swrt, 0));
    // SwrtNizesAndNdims = Swrt[[2]]
    PROTECT(SwrtSizesAndNdims = VECTOR_ELT(Swrt, 1));
    
    populateValueMapAccessorsFromNodeNames_internal(&model_wrt_accessor,
						    SwrtNodeNames,
						    SwrtSizesAndNdims,
						    SbasePtr);

    populateValueMapAccessorsFromNodeNames_internal(&model_AD_wrt_accessor,
						    SwrtNodeNames,
						    SwrtSizesAndNdims,
						    SbasePtrAD);

    //SextraInput <- SderivsInfo$extraInputMapInfo
    SEXP SextraInput;
    SEXP SextraInputNodeNames, SextraInputSizesAndNdims;
    PROTECT(SextraInput =
	    Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, SpxData)),
			      Rf_install("extraInputMapInfo")));
    // SextraInputNodeNames = SextraInput[[1]]
    PROTECT(SextraInputNodeNames = VECTOR_ELT(SextraInput, 0));
    // SextraInputNizesAndNdims = SextraInput[[2]]
    PROTECT(SextraInputSizesAndNdims = VECTOR_ELT(SextraInput, 1));
    
    populateValueMapAccessorsFromNodeNames_internal(&model_extraInput_accessor,
						    SextraInputNodeNames,
						    SextraInputSizesAndNdims,
						    SbasePtr);

    populateValueMapAccessorsFromNodeNames_internal(&model_AD_extraInput_accessor,
						    SextraInputNodeNames,
						    SextraInputSizesAndNdims,
						    SbasePtrAD);

    //SmodelOutput <- SderivsInfo$modelOutputMapInfo
    SEXP SmodelOutput;
    SEXP SmodelOutputNodeNames, SmodelOutputSizesAndNdims;
    PROTECT(SmodelOutput =
	    Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, SpxData)),
			      Rf_install("modelOutputMapInfo")));
    // SmodelOutputNodeNames = SmodelOutput[[1]]
    PROTECT(SmodelOutputNodeNames = VECTOR_ELT(SmodelOutput, 0));
    // SmodelOutputNizesAndNdims = SmodelOutput[[2]]
    PROTECT(SmodelOutputSizesAndNdims = VECTOR_ELT(SmodelOutput, 1));
    
    populateValueMapAccessorsFromNodeNames_internal(&model_modelOutput_accessor,
						    SmodelOutputNodeNames,
						    SmodelOutputSizesAndNdims,
						    SbasePtr);

    populateValueMapAccessorsFromNodeNames_internal(&model_AD_modelOutput_accessor,
						    SmodelOutputNodeNames,
						    SmodelOutputSizesAndNdims,
						    SbasePtrAD);

    //Sconstant <- SderivsInfo$constantMapInfo
    SEXP Sconstant;
    SEXP SconstantNodeNames, SconstantSizesAndNdims;
    PROTECT(Sconstant =
	    Rf_findVarInFrame(PROTECT(GET_SLOT(SderivsInfo, SpxData)),
			      Rf_install("constantMapInfo")));
    // SconstantNodeNames = Sconstant[[1]]
    PROTECT(SconstantNodeNames = VECTOR_ELT(Sconstant, 0));
    // SconstantNizesAndNdims = Sconstant[[2]]
    PROTECT(SconstantSizesAndNdims = VECTOR_ELT(Sconstant, 1));
    
    populateValueMapAccessorsFromNodeNames_internal(&model_constant_accessor,
						    SconstantNodeNames,
						    SconstantSizesAndNdims,
						    SbasePtr);

    populateValueMapAccessorsFromNodeNames_internal(&model_AD_constant_accessor,
						    SconstantNodeNames,
						    SconstantSizesAndNdims,
						    SbasePtrAD);

    
    UNPROTECT(26);
  }
};



#endif

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

  
#endif
