#ifndef __NIMARRBASE
#define __NIMARRBASE

/* fix to avoid warnings exemplified by edison.nersc.gov SUSE Linux - Github
 * issue #214 */
#if defined __GNUC__ && __GNUC__ >= 6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <R.h>
#include <cstring>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

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

enum nimType { INT = 1, DOUBLE = 2, BOOL = 3, UNDEFINED = -1 };

class NimArrType {
 public:
  nimType myType;
  virtual nimType getNimType() const { return myType; }
  virtual ~NimArrType() {}
};

class NimVecType {
 public:
  nimType myType;
  virtual nimType getNimType() const { return myType; }
  virtual NimArrType *getRowTypePtr(int row) = 0;
  virtual int size() = 0;
  virtual void setRowDims(int row, vector<int> dims) = 0;
  virtual vector<int> getRowDims(int row) = 0;
  virtual ~NimVecType() {}
};

template <class T>
class NimArrBase : public NimArrType {
 public:
  T *v;
  T **vPtr;
  std::size_t element_size() { return (sizeof(T)); }
  void setVptr() { vPtr = &v; }
  T **getVptr() const { return vPtr; }
  bool own_v;
  int NAdims[4];
  const int *dim() const { return NAdims; }
  int NAstrides[4];
  // Everyone has a stride1, and the flat [] operator needs it, so it is here.
  int stride1, offset;
  int getOffset() { return offset; }
  bool boolMap;
  bool isMap() const { return boolMap; }
  const int *strides() const { return NAstrides; }
  // Name length can cause problems if R headers have been #include'd, as they
  // have #define length Rf_length.
  int NAlength;
  int size() const { return NAlength; }
  virtual int numDims() const = 0;
  virtual int dimSize(int i) const = 0;
  // Generic for nDim > 1, overloaded for other dimensions.
  T &operator[](int i) const { return ((*vPtr)[offset + i * stride1]); }
  // Only to be used if not a map .
  T &valueNoMap(int i) const { return (*(v + i)); }
  virtual int calculateIndex(vector<int> &i) const = 0;
  T *getPtr() { return (&((*vPtr)[0])); }
  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) = 0;
  // Warning, this does not make sense if vPtr is pointing to someone else's
  // vMemory.
  void setLength(int l, bool copyValues = true, bool fillZeros = true) {
    if (NAlength == l) {
      if ((!copyValues) & fillZeros) fillAllValues(static_cast<T>(0));
      return;
    }
    T *new_v = new T[l];
    if (own_v) {
      if (copyValues) {
        if (l < NAlength) {
          std::copy(v, v + l, new_v);
        } else {
          std::copy(v, v + NAlength, new_v);
          if (fillZeros) {
            std::fill(new_v + NAlength, new_v + l, static_cast<T>(0));
          }
        }
      } else {
        if (fillZeros) std::fill(new_v, new_v + l, static_cast<T>(0));
      }
      delete[] v;
    }
    NAlength = l;
    v = new_v;
    own_v = true;
  }
  void fillAllValues(T value) { std::fill(v, v + NAlength, value); }
  void fillAllValues(T value, bool fillZeros, bool recycle) {
    if (recycle) {
      std::fill(v, v + NAlength, value);
    } else {
      if (NAlength > 0) v[0] = value;
      if (NAlength > 1)
        if (fillZeros) std::fill(v + 1, v + NAlength, static_cast<T>(0));
    }
  }
  void setMyType() {
    myType = UNDEFINED;
    if (typeid(T) == typeid(int)) {
      myType = INT;
    } else if (typeid(T) == typeid(double)) {
      myType = DOUBLE;
    } else if (typeid(T) == typeid(bool)) {
      myType = BOOL;
    }
  }
  virtual ~NimArrBase() {
    if (own_v) delete[] v;
  }
  // Do we ever use this case?
  NimArrBase(const NimArrBase<T> &other)
      :  // own_v isn't a map but we'll only set to true when giving it values.
        own_v(false),
        offset(0),
        boolMap(false),
        NAlength(other.size()) {
    std::memcpy(NAdims, other.dim(), other.numDims() * sizeof(int));
    myType = other.getNimType();
  }

  NimArrBase()
      : v(), vPtr(&v), own_v(false), offset(0), boolMap(false), NAlength(0) {
    setMyType();
  }
  NimArrBase(const T *&vm, int off)
      : vPtr(&vm), own_v(false), offset(off), boolMap(true) {
    setMyType();
  }
  template <class Tfrom>
  void genericMapCopy(int offset, vector<int> &str, vector<int> &is,
                      NimArrBase<Tfrom> *from, int fromOffset,
                      vector<int> &fromStr, vector<int> &fromIs);
};

template <class T>
class VecNimArrBase : public NimVecType {
 public:
  virtual void resize(int i) = 0;
  virtual NimArrBase<T> *getBasePtr(int i) = 0;
  NimArrType *getRowTypePtr(int row) {
    return (static_cast<NimArrType *>(getBasePtr(row)));
  }

  VecNimArrBase() {
    myType = UNDEFINED;
    if (typeid(T) == typeid(int)) {
      myType = INT;
    } else if (typeid(T) == typeid(double)) {
      myType = DOUBLE;
    } else if (typeid(T) == typeid(bool)) {
      myType = BOOL;
    }
  }

  ~VecNimArrBase() {}
};

#endif
