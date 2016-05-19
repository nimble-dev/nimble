#ifndef __NIMARRBASE
#define __NIMARRBASE
#include <vector>
#include <string>
#include <R.h>
#include <typeinfo>
#include <iostream>

using std::vector;

enum nimType {INT = 1, DOUBLE = 2, UNDEFINED = -1};

 class NimArrType{
	public:
	nimType myType;
	virtual nimType getNimType() const {return(myType);};
	virtual ~NimArrType(){};
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
  vector<T> v;
  vector<T> *vPtr;
  void setVptr() {vPtr = &v;}
  vector<T> *getVptr() const {return(vPtr);}
  vector<int> NAdims;
  const vector<int> &dim() const {return(NAdims);}
  vector<int> NAstrides;
  int stride1, offset; // everyone has a stride1, and the flat [] operator needs it, so it is here. 
  int getOffset() {return(offset);}
  bool boolMap;
  bool isMap() const {return(boolMap);}
  const vector<int> &strides() const {return(NAstrides);}
  int NAlength; // name length can cause problems if R headers have been #include'd, as they have #define length Rf_length
  int size() const {return(NAlength);}
  virtual int numDims() const = 0;
  virtual int dimSize(int i) const = 0;
  T &operator[](int i) {return((*vPtr)[offset + i * stride1]);} // could be misused for nDim > 1
  virtual int calculateIndex(vector<int> &i)=0;
  T *getPtr() {return(&((*vPtr)[0]));}
  virtual void setSize(vector<int> sizeVec)=0;
  void fillAllValues(T value) { std::fill(v.begin(), v.end(), value); }
  void setLength(int l) {NAlength = l; v.resize(l); } // Warning, this does not make sense if vPtr is pointing to someone else's vMemory. 
  void setMyType() {
    myType = UNDEFINED;
    if(typeid(T) == typeid(int) )
      myType = INT;
    if(typeid(T) == typeid(double) )
      myType = DOUBLE;
  }
  ~NimArrBase(){};
 NimArrBase(const NimArrBase<T> &other) :
  NAdims(other.dim()),
    offset(0),
    boolMap(false),
    NAlength(other.size())
      {
	// std::cout<<"Using copy constructor for a NimArrBase<T>\n";
   	myType = other.getNimType();
      };

 NimArrBase() : v(), vPtr(&v), offset(0), boolMap(false), NAlength(0) {
    //  std::cout<<"Creating a NimArr with &v = "<<&v<<" and "<<" vPtr = "<<vPtr<<"\n";
    setMyType();
  }
 NimArrBase(const vector<T> &vm, int off) : vPtr(&vm), offset(off), boolMap(true) {
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
  }
  
  ~VecNimArrBase(){};
};


#endif
