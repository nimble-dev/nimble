#ifndef __NIMARRBASE
#define __NIMARRBASE

/* fix to avoid warnings exemplified by edison.nersc.gov SUSE Linux - Github issue #214 */
#if defined __GNUC__ && __GNUC__>=6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif


#include <vector>
#include <string>
#include <cstring>
#include <R.h>
#include <typeinfo>
#include <iostream>

/* #ifdef _WIN32 */
/* #define _WIN3264 */
/* #endif */

/* #ifdef _WIN64 */
/* #define _WIN3264 */
/* #endif */

/* #ifdef _WIN3264 */
/* #pragma GCC diagnostic ignored "-Wmaybe-uninitialized" */
/* #endif */

using std::vector;

enum nimType {INT = 1, DOUBLE = 2, BOOL = 3, UNDEFINED = -1};

 class NimArrType{
	public:
	nimType myType;
	virtual nimType getNimType() const {return(myType);};
	virtual ~NimArrType(){//Rprintf("In NimArrType destructor\n");
	};
 };


  class NimVecType{
	public:
	nimType myType;
	virtual nimType getNimType() const {return(myType);};
	virtual NimArrType* getRowTypePtr(int row)=0;
	virtual int size() =0;
	virtual void setRowDims(int row, vector<int> dims) = 0;
	virtual vector<int> getRowDims(int row) = 0;
	virtual ~NimVecType(){};
  };

template<class T>
class NimArrBase: public NimArrType {
 public:
  //vector<T> v;
  T *v;
  //vector<T> *vPtr;
  T **vPtr;
  std::size_t element_size() {return(sizeof(T));}
  void setVptr() {vPtr = &v;}
  //  vector<T> *getVptr() const {return(vPtr);}
  T **getVptr() const{return(vPtr);}
  bool own_v;
  //vector<int> NAdims;
  //  int *NAdims;
  int NAdims[4];
  //  const vector<int> &dim() const {return(NAdims);}
  const int* dim() const {return(NAdims);}
    //vector<int> NAstrides;
  //  int *NAstrides;
  int NAstrides[4];
  int stride1, offset; // everyone has a stride1, and the flat [] operator needs it, so it is here.
  int getOffset() {return(offset);}
  bool boolMap;
  bool isMap() const {return(boolMap);}
  //const vector<int> &strides() const {return(NAstrides);}
  const int* strides() const {return(NAstrides);}
  int NAlength; // name length can cause problems if R headers have been #include'd, as they have #define length Rf_length
  int size() const {return(NAlength);}
  virtual int numDims() const = 0;
  virtual int dimSize(int i) const = 0;
  T &operator[](int i) const {return((*vPtr)[offset + i * stride1]);} // generic for nDim > 1, overloaded for other dimensions
  T &valueNoMap(int i) const {return(*(v + i));} // only to be used if not a map 
  virtual int calculateIndex(vector<int> &i) const =0;
  T *getPtr() {return(&((*vPtr)[0]));}
  virtual void setSize(vector<int> sizeVec, bool copyValues = true, bool fillZeros = true)=0;
  void setLength(int l, bool copyValues = true, bool fillZeros = true) {
    if(NAlength==l) {
      if((!copyValues) & fillZeros) fillAllValues(static_cast<T>(0));
      return;
    }
    T *new_v = new T[l];
    if(own_v) {
      if(copyValues) {
	if(l < NAlength) std::copy(v, v + l, new_v);
	else {
	  std::copy(v, v + NAlength, new_v);
	  if(fillZeros) {
	    std::fill(new_v + NAlength, new_v + l, static_cast<T>(0));
	  }
	}
      } else {
	if(fillZeros)
	  std::fill(new_v, new_v + l, static_cast<T>(0));
      }
      delete[] v;
    }
    NAlength = l;
    //v.resize(l);
    v = new_v;
    own_v = true;
  } // Warning, this does not make sense if vPtr is pointing to someone else's vMemory.
  void fillAllValues(T value) { std::fill(v, v + NAlength, value); }
  void fillAllValues(T value, bool fillZeros, bool recycle) {
    if(recycle) {
      std::fill(v, v + NAlength, value);
    } else {
      if(NAlength > 0) v[0] = value;
      if(NAlength > 1)
	if(fillZeros)
	  std::fill(v + 1, v + NAlength, static_cast<T>(0));
    }
  }
  void setMyType() {
    myType = UNDEFINED;
    if(typeid(T) == typeid(int) )
      myType = INT;
    if(typeid(T) == typeid(double) )
      myType = DOUBLE;
    if(typeid(T) == typeid(bool) )
      myType = BOOL;

  }
  virtual ~NimArrBase(){
    //delete[] NAdims;
    //delete[] NAstrides;
    //Rprintf("In NimArrBase destructor\n");
    if(own_v) delete [] v;
  };
 NimArrBase(const NimArrBase<T> &other) : // do we ever use this case?
  //NAdims(other.dim()),
    own_v(false), // isn't a map but we'll only set to true when giving it values.
      offset(0),
      boolMap(false),
    NAlength(other.size())
      {
	//	NAdims = new int[other.numDims()];
	std::memcpy(NAdims, other.dim(), other.numDims()*sizeof(int));
	// std::cout<<"Using copy constructor for a NimArrBase<T>\n";
   	myType = other.getNimType();
      };

 NimArrBase() : v(), vPtr(&v), own_v(false), offset(0), boolMap(false), NAlength(0) {
    //  std::cout<<"Creating a NimArr with &v = "<<&v<<" and "<<" vPtr = "<<vPtr<<"\n";
    setMyType();
  }
 NimArrBase(const T* &vm, int off) : vPtr(&vm), own_v(false), offset(off), boolMap(true) {
    setMyType();
  }
  template<class Tfrom>
    void genericMapCopy(int offset, vector<int> &str, vector<int> &is, NimArrBase<Tfrom> *from, int fromOffset, vector<int> &fromStr, vector<int> &fromIs);
};

template<class T>
class VecNimArrBase : public NimVecType {
 public:
  virtual void resize(int i)=0;
  virtual NimArrBase<T>* getBasePtr(int i)=0;
  NimArrType* getRowTypePtr(int row){
   	return(static_cast<NimArrType *> (getBasePtr(row) )  );
   }

  VecNimArrBase() {
    myType = UNDEFINED;
    if(typeid(T) == typeid(int) )
      myType = INT;
    if(typeid(T) == typeid(double) )
      myType = DOUBLE;
    if(typeid(T) == typeid(bool) )
      myType = BOOL;

  }

  ~VecNimArrBase(){};
};


#endif
