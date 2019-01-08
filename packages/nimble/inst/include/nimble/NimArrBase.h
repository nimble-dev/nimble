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

#ifndef __NIMARRBASE
#define __NIMARRBASE

/* fix to avoid warnings exemplified by edison.nersc.gov SUSE Linux - Github
 * issue #214 */
/* based on CRAN complaint, remove this as of 0.6-9 */
/* #if defined __GNUC__ && __GNUC__ >= 6 */
/* #pragma GCC diagnostic ignored "-Wignored-attributes" */
/* #endif */

#include <R.h>
#include <Eigen/Core>
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

// This allows us to experiment with different allocators. In theory aligned
// memory should improve Eigen performance, but microbenchmarks have failed to
// substantiate this.
#if !defined(NIMBLE_ALIGNED_MALLOC)
#define NIMBLE_ALIGNED_MALLOC 0
#endif  // !defined(NIMBLE_ALIGNED_MALLOC)

// These wrappers help us use Eigen::internal components, while minimizing the
// cost of updating if those components change in a future Eigen release.
template <class T>
inline T *nimble_malloc(size_t size) {
#if NIMBLE_ALIGNED_MALLOC
  return static_cast<T *>(Eigen::internal::aligned_malloc(size * sizeof(T)));
#else   // NIMBLE_ALIGNED_MALLOC
  return new T[size];
#endif  // NIMBLE_ALIGNED_MALLOC
}

template <class T>
inline void nimble_free(T *ptr) {
#if NIMBLE_ALIGNED_MALLOC
  Eigen::internal::aligned_free(ptr);
#else   // NIMBLE_ALIGNED_MALLOC
  delete[] ptr;
#endif  // NIMBLE_ALIGNED_MALLOC
}

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
  T** &getVptrRef() {return vPtr;}
  T **getVptr() const { return vPtr; }
  bool own_v;
  int NAdims[5];
  const int *dim() const { return NAdims; }
  int NAstrides[5];
  // Everyone has a stride1, and the flat [] operator needs it, so it is here.
  int stride1, offset;
  int getOffset() const { return offset; }
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
  const T *getConstPtr() const { return (&((*vPtr)[0])); }
  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) = 0;
  // Warning, this does not make sense if vPtr is pointing to someone else's
  // vMemory.
  void setLength(int l, bool copyValues = true, bool fillZeros = true) {
    if (NAlength == l) {
      if ((!copyValues) & fillZeros) fillAllValues(static_cast<T>(0));
      return;
    }
    T *new_v = nimble_malloc<T>(l);
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
      nimble_free(v);
    } else {
      if (fillZeros) std::fill(new_v, new_v + l, static_cast<T>(0));
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
    if (own_v) nimble_free(v);
  }
  // Base class copy constructor:
  // This will be used by derived class copy constructors,
  // which will be called when a VecNimArr resizes its vector of NimArr<>s.
  NimArrBase(const NimArrBase<T> &other)
    :   own_v(false), // this may get set to true in a derived copy constructor
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
