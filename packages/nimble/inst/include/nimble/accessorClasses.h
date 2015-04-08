// VERY IMPORTANT NOTE: THESE FUNCTIONS WILL ASSUME MODELS ARE OF CLASS
// "ModelBase". THIS CLASS CAN BE FOUND IN "ModelClassUtils.h" 



// This file contains new components for accessing and copying from groups of node Functions or subsets of variables in models or modelValues
#ifndef __ACCESSORCLASSES
#define __ACCESSORCLASSES

#include <iostream>

#include "NimArrBase.h"
#include "NimArr.h"			
#include "ModelClassUtils.h"
#include "RcppUtils.h"
#include <Rinternals.h>
#include "R.h"


using std::cout;

#include "nodeFun.h" 

/////////////////////////////////
// 1. NodeVectors:
/////////////////////////////////
class NodeVectorClass {
 public:
  vector<nodeFun *> nodeFunPtrs;
  virtual vector<nodeFun *> &getNodeFunctionPtrs() {return(nodeFunPtrs);}
  // to be inherited and implemented differently when we have dynamic dependencies
  virtual ~NodeVectorClass() {};
};

///// Using NodeVectors:
// utilities for calling node functions from a vector of node pointers
// see .cpp file for definitions
double calculate(NodeVectorClass &nodes);
double getLogProb(NodeVectorClass &nodes);
void simulate(NodeVectorClass &nodes);


/////////////////////
// new version of variable accessors using maps (offset and strided windows into multivariate objects (NimArr<>s) )
/////////////////////

// single and many base classes
class SingleVariableMapAccessBase {
 public:
  int offset;
  bool singleton;
  vector<int> sizes, strides; 
  virtual ~SingleVariableMapAccessBase() {};
  virtual NimArrType *getNimArrPtr()=0;
  int &getOffset() {return(offset);}
  vector<int> &getSizes() {return(sizes);}
  vector<int> &getStrides() {return(strides);}
  bool &getSingleton() {return(singleton);}
  virtual void setObject(void *object);
};


class ManyVariablesMapAccessorBase {
 public:
  virtual vector<SingleVariableMapAccessBase *> &getMapAccessVector()=0;
  virtual void  setRow(int i) = 0;
  virtual ~ManyVariablesMapAccessorBase() {};
};

// single and many variable classes

class SingleVariableMapAccess : public SingleVariableMapAccessBase {
  NimArrType **ppVar; // I think we have to do some casting when populating this
  virtual NimArrType *getNimArrPtr() {return(*ppVar);}
  ~SingleVariableMapAccess() {};
  void setObject(void *object) {ppVar = static_cast<NimArrType**>(object);}
};

class ManyVariablesMapAccessor : public ManyVariablesMapAccessorBase {
 public:
  vector<SingleVariableMapAccessBase *> varAccessors;
  virtual vector<SingleVariableMapAccessBase *> &getMapAccessVector() {return(varAccessors);}
  ~ManyVariablesMapAccessor();
  void setRow(int i){PRINTF("Bug detected in code: attempting to setRow for model. Can only setRow for modelValues\n");}
};

// single and many modelValues classes
class SingleModelValuesMapAccess : public SingleVariableMapAccessBase {
 public:
  NimVecType *pVVar;   
  int currentRow;
  virtual NimArrType *getNimArrPtr() {return(pVVar->getRowTypePtr(currentRow));} 
  ~SingleModelValuesMapAccess() {};
  void setRow(int i) {currentRow = i;}
  int getRow() {return(currentRow);}
  void setObject(void *object) {pVVar = static_cast<NimVecType*>(object);}
};

class ManyModelValuesMapAccessor : public ManyVariablesMapAccessorBase {
  public:
  int currentRow;
  vector<SingleVariableMapAccessBase *> varAccessors;
  virtual vector<SingleVariableMapAccessBase *> &getMapAccessVector() {return(varAccessors);}
  virtual void setRow(int i);// see .cpp
  ~ManyModelValuesMapAccessor() {};
};

void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to);
void nimCopyOne(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to);

void nimCopy(ManyVariablesMapAccessorBase &from, ManyVariablesMapAccessorBase &to) { //map version

  vector<SingleVariableMapAccessBase *> fromAccessors = from.getMapAccessVector();
  vector<SingleVariableMapAccessBase *> toAccessors = to.getMapAccessVector();

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

void nimCopyOne(SingleVariableMapAccessBase *from, SingleVariableMapAccessBase *to) { // map version
  nimType fromType, toType;
  NimArrType *fromNimArr, *toNimArr;
  fromNimArr = from->getNimArrPtr();
  toNimArr = to->getNimArrPtr();

  fromType = fromNimArr->getNimType();
  toType = toNimArr->getNimType();

  if(to->singleton) {
    switch(fromType) {
    case DOUBLE:
      switch(toType) {
      case DOUBLE:
	(*static_cast<NimArrBase<double> *>(toNimArr))[from->offset] = (*static_cast<NimArrBase<double> *>(fromNimArr))[from->offset];
      break;
    case INT:
	(*static_cast<NimArrBase<int> *>(toNimArr))[from->offset] = (*static_cast<NimArrBase<double> *>(fromNimArr))[from->offset];
      break;
    default:
      cout<<"Error in nimCopyOne: unknown type for destination\n";
    }
    break;
  case INT:
    switch(toType) {
    case DOUBLE:
      (*static_cast<NimArrBase<double> *>(toNimArr))[from->offset] = (*static_cast<NimArrBase<int> *>(fromNimArr))[from->offset];
      break;
    case INT:
	(*static_cast<NimArrBase<int> *>(toNimArr))[from->offset] = (*static_cast<NimArrBase<int> *>(fromNimArr))[from->offset];
      break;
    default:
      cout<<"Error in nimCopyOne: unknown type for destination\n";
    }
    break;
  default:
    cout<<"Error in nimCopyOne: unknown type for source\n";
    }
  } else {
    
    switch(fromType) {
    case DOUBLE:
      switch(toType) {
      case DOUBLE:
	static_cast<NimArrBase<double> *>(toNimArr)->genericMapCopy<double>(to->getOffset(), to->getStrides(), to->getSizes(), static_cast<NimArrBase<double> *>(fromNimArr), from->getOffset(), from->getStrides(), from->getSizes() );
      break;
      case INT:
	static_cast<NimArrBase<int> *>(toNimArr)->genericMapCopy<double>(to->getOffset(), to->getStrides(), to->getSizes(), static_cast<NimArrBase<double> *>(fromNimArr), from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      default:
	cout<<"Error in nimCopyOne: unknown type for destination\n";
      }
      
    case INT:
      switch(toType) {
      case DOUBLE:
	static_cast<NimArrBase<double> *>(toNimArr)->genericMapCopy<int>(to->getOffset(), to->getStrides(), to->getSizes(), static_cast<NimArrBase<int> *>(fromNimArr), from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      case INT:
	static_cast<NimArrBase<int> *>(toNimArr)->genericMapCopy<int>(to->getOffset(), to->getStrides(), to->getSizes(), static_cast<NimArrBase<int> *>(fromNimArr), from->getOffset(), from->getStrides(), from->getSizes() );
	break;
      default:
	cout<<"Error in nimCopyOne: unknown type for destination\n";
      }
      
    default:
      cout<<"Error in nimCopyOne: unknown type for destination\n";
    }
  }
}

///////////////////////////////
// 2. Variable accessors
//
// These are trickier because there are different types inside of NimArr<>s.
// We now have a type system so we can extract the type in the copy function.
// We *may* need to templatize types here, but I think we can avoid it.
///////////////////////////////

// Base class for access to one NimArr<>
class SingleVariableAccessBase {
 public:
  int flatIndexStart, flatIndexEnd; // For 1-3: In R these should be 1,3; In C++: 0, 3
  int length; // copying function can check for singletons if it wants to
  int getIndexStart() {return(flatIndexStart);}  
  int getIndexEnd() {return(flatIndexEnd);}
  int getLength() {return(length);}
  virtual NimArrType *getNimArrPtr()=0; //
  virtual ~SingleVariableAccessBase() {};
};

// Derived class for access to one NimArr<> in a model
class SingleVariableAccess : public SingleVariableAccessBase {
 public:
  NimArrType **ppVar; // I think we have to do some casting when populating this
  virtual NimArrType *getNimArrPtr() {return(*ppVar);}
  ~SingleVariableAccess() {};
};

// Base class for vector of single variables accessors
class ManyVariablesAccessorBase {
 public:
  virtual vector<SingleVariableAccessBase *> &getAccessVector()=0;
  virtual void  setRow(int i) = 0;
  virtual ~ManyVariablesAccessorBase() {};
};

// Derived class for vector of single variable accessors to NimArr<>s in a model
class ManyVariablesAccessor : public ManyVariablesAccessorBase {
 public:
  vector<SingleVariableAccessBase *> varAccessors;
  virtual vector<SingleVariableAccessBase *> &getAccessVector() {return(varAccessors);}
  ~ManyVariablesAccessor();
  void setRow(int i){PRINTF("Bug detected in code: attempting to setRow for model. Can only setRow for modelValues\n");}
};

/////////////////////////////////
// 3. modelValues accessors
// 
// These are like model variable accessors but need a row index
/////////////////////////////////

// Derived class for access to one NimArr<> in a modelValues
// This uses some things not yet written: a vecNimArrType base class with a getRowTypePtr().
class SingleModelValuesAccess : public SingleVariableAccessBase {
 public:
  NimVecType *pVVar;   // Cliff and I talked about making a vecNimArrType base class
  int currentRow;
  virtual NimArrType *getNimArrPtr() {return(pVVar->getRowTypePtr(currentRow));} // Need to put a function like this in vecNimArrType base class
  ~SingleModelValuesAccess() {};
  void setRow(int i) {currentRow = i;}
  int getRow() {return(currentRow);}
};

// Derived class for a vector of single variable accessors to NimArr<>s in a modelValues
// The row(i) member function sets the currentRow of all the single accessors and returns the vector of them
class ManyModelValuesAccessor : public ManyVariablesAccessorBase {
  public:
  int currentRow;
  vector<SingleVariableAccessBase *> varAccessors;
  virtual vector<SingleVariableAccessBase *> &getAccessVector() {return(varAccessors);}
  virtual void setRow(int i);// see .cpp
  ~ManyModelValuesAccessor() {};
};

/////////////////////////////////
// nimCopy function
//
// Now we are ready to write a fairly general copy function
// I am calling it nimCopy to avoid name conflicts with std::copy or others.
/////////////////////////////////

void nimCopy(ManyVariablesAccessorBase &from, ManyVariablesAccessorBase &to);
void nimCopyOne(SingleVariableAccessBase *from, SingleVariableAccessBase *to);

// This templated piece is given in the .h file
template<class Tfrom, class Tto>
void nimCopyOneTyped(SingleVariableAccessBase *fromSVA, SingleVariableAccessBase *toSVA) {
  NimArrBase<Tfrom> *fromNimPtr = static_cast<NimArrBase<Tfrom> *> ( (fromSVA)->getNimArrPtr() ); //	I don't believe static casting should be necessary
  NimArrBase<Tto>  *toNimPtr = static_cast<NimArrBase<Tto> *> ( (toSVA)->getNimArrPtr() );		//	Same
  if(fromSVA->getLength() != toSVA->getLength()) {
    cout<<"Error in nimCopyOneTyped: lengths do not match.\n";
    cout << "FromLength = " << fromSVA->getLength() << " ToLength = "<< toSVA->getLength() << "\n";
  	return;
  }
  if(fromSVA->getLength() == 1) {
    (*toNimPtr)[toSVA->getIndexStart()] = (*fromNimPtr)[fromSVA->getIndexStart()];
    return;
  }
  
  
  std::copy( fromNimPtr->getPtr() + fromSVA->getIndexStart(),
	    fromNimPtr->getPtr() + fromSVA->getIndexEnd() + 1,
	    toNimPtr->getPtr() + toSVA->getIndexStart());
}

void nimCopy(ManyVariablesAccessorBase &from, int rowFrom, ManyVariablesAccessorBase &to);

void nimCopy(ManyVariablesAccessorBase &from, int rowFrom, ManyVariablesAccessorBase &to, int rowTo);

void nimCopy(ManyVariablesAccessorBase &from, ManyVariablesAccessorBase &to, int rowTo);
	


/* template<int D, class T> */
/* void nimArr_2_SingleModelAccess(SingleVariableAccess* SMVAPtr, NimArr<D, T>* nimArrPtr, int nimBegin); */
/* template<int D, class T> */
/* void nimArr_2_ManyModelAccess(ManyVariablesAccessor &MMVAPtr, NimArr<D, T>* nimArrPtr); */

template<class T>
void nimArr_2_SingleModelAccess(SingleVariableAccess* SMVAPtr, NimArrBase<T>* nimArrPtr, int nimBegin);
template<class T>
void nimArr_2_ManyModelAccess(ManyVariablesAccessor &MMVAPtr, NimArrBase<T>* nimArrPtr);

template<int D, class T>
void SingleModelAccess_2_nimArr(SingleVariableAccess* SMVAPtr, NimArr<D, T>* nimArrPtr, int nimBegin);
template<int D, class T>
void ManyModelAccess_2_nimArr(ManyVariablesAccessor &MMVAPtr, NimArr<D, T>* nimArrPtr);

/* void setValues(NimArr<1, double> &nimArr, ManyVariablesAccessor &MVA); */
/* void setValues(NimArr<1, int> &nimArr, ManyVariablesAccessor &MVA); */

void setValues(NimArrBase<double> &nimArr, ManyVariablesAccessor &MVA);
void setValues(NimArrBase<int> &nimArr, ManyVariablesAccessor &MVA);

void getValues(NimArr<1, double> &nimArr, ManyVariablesAccessor &MVA);
void getValues(NimArr<1, int> &nimArr, ManyVariablesAccessor &MVA);


double calculate(NodeVectorClass &nodes);
double getLogProb(NodeVectorClass &nodes);
void simulate(NodeVectorClass &nodes);

void cAddNodeFun(NodeVectorClass nVObj, nodeFun* nFPtr);
void cAddVariableAccessor(ManyVariablesAccessor* mVAPtr, SingleVariableAccess* sVAPtr, int index);
template<class T>
void cRemoveAccessor(T* aPtr, int index, bool removeAll);
void setModelValuesAccessorRow(ManyVariablesAccessorBase &access);


SingleModelValuesAccess* cMakeSingleModelValuesAccessor(NimVecType* varPtr, int beginIndex, int endIndex, int row);
SingleVariableAccess* cMakeSingleVariableAccessor(NimArrType** varPtr, int beginIndex, int endIndex);

extern "C" {
	SEXP makeSingleVariableAccessor(SEXP rModelPtr, SEXP elementName,  SEXP beginIndex, SEXP endIndex);
	SEXP makeSingleModelValuesAccessor(SEXP rModelValuesPtr, SEXP elementName,  SEXP curRow, SEXP beginIndex, SEXP endIndex);

	SEXP getModelAccessorValues(SEXP accessor);
	SEXP getMVAccessorValues(SEXP accessor);

	SEXP newNodeFxnVector(SEXP size);
	SEXP setNodeModelPtr(SEXP nodeFxnPtr, SEXP modelElementPtr, SEXP nodeElementName);
	SEXP resizeNodeFxnVector(SEXP nodeFxnVecPtr, SEXP size);
	SEXP addNodeFun(SEXP nVPtr, SEXP nFPtr, SEXP addAtEnd, SEXP index);
	SEXP removeNodeFun(SEXP rPtr, SEXP index, SEXP removeAll);
	
	SEXP newManyVariableAccessor(SEXP size);
	SEXP addSingleVariableAccessor(SEXP MVAPtr, SEXP SVAPtr, SEXP addAtEnd, SEXP index);
	SEXP resizeManyModelVarAccessor(SEXP manyModelVarPtr, SEXP size);
	SEXP removeModelVariableAccessor(SEXP rPtr, SEXP index, SEXP removeAll);
	
	SEXP newManyModelValuesAccessor(SEXP size);
	SEXP resizeManyModelValuesAccessor(SEXP manyModelValuesPtr, SEXP size);
	SEXP addSingleModelValuesAccessor(SEXP MVAPtr, SEXP SVAPtr, SEXP addAtEnd, SEXP index);
	SEXP removeModelValuesAccessor(SEXP rPtr, SEXP index, SEXP removeAll);
	 
	SEXP manualSetNRows(SEXP Sextptr, SEXP nRows);

	SEXP populateNodeFxnVector(SEXP nodeFxnVec, SEXP nodeNames, SEXP );
    SEXP populateNodeFxnVector_byGID(SEXP SnodeFxnVec, SEXP S_GIDs, SEXP SnumberedObj);
	
//	SEXP populateNumberedObject_withSingleModelValuesAccessors(SEXP mvPtr, SEXP varName, SEXP beginIndex, SEXP varLength, SEXP curRow, SEXP SnumbObj);
	SEXP populateNumberedObject_withSingleModelValuesAccessors(SEXP mvPtr, SEXP varName, SEXP GIDs, SEXP curRow, SEXP SnumbObj);
	SEXP populateModelValuesAccessors_byGID(SEXP SmodelValuesAccessorVector, SEXP S_GIDs, SEXP SnumberedObj);
	
	SEXP populateNumberedObject_withSingleModelVariablesAccessors(SEXP modelPtr, SEXP varName, SEXP sGIDS, SEXP SvalidIndices, SEXP SnumbObj);
	SEXP populateModelVariablesAccessors_byGID(SEXP SmodelVariableAccessorVector, SEXP S_GIDs, SEXP SnumberedObj, SEXP S_LP_GIDs, SEXP S_LP_numberedObj);

	SEXP new_SingleModelValuesAccessor_NumberedObjects();
	SEXP new_SingleModelVariablesAccessor_NumberedObjects();
}
void  SingleVA_Finalizer ( SEXP Sv );
void  SingleMVA_Finalizer ( SEXP Sv );
void NodeVector_Finalizer( SEXP Sv);
void ManyVariable_Finalizer(SEXP Sv);
void ManyMV_Finalizer(SEXP Sv);
#endif 			
