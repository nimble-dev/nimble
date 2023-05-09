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

#ifndef __NIMARR
#define __NIMARR

#include <cstdlib>
#include "NimArrBase.h"
#include "Utils.h"

template <int ndim, class T>
class NimArr;

// needed for externalCalls
template<int nDim, class T>
T* nimArrPtr_copyIfNeeded(NimArr<nDim, T> &orig, NimArr<nDim, T> &possibleCopy) {
  if(orig.isMap()) {
    // if(!isMapEntire(orig)) { // not ready for this yet.
      possibleCopy = orig;
      return(possibleCopy.getPtr());
      // }
  }
  return(orig.getPtr());
}

template<int nDim, class T>
  void nimArrPtr_copyBackIfNeeded(T* tptr, NimArr<nDim, T> &orig, NimArr<nDim, T> &possibleCopy) {
  if(tptr == orig.getPtr()) return;
  if(tptr != possibleCopy.getPtr()) {
    NIMERROR("Problem in unconverting from an external call\n");
  }
  orig.mapCopy(possibleCopy);
}

// only case where T and T2 need to be different is
// T = bool and
// T2 = int, since R represents bools as ints
template<int nDim, class T, class T2>
  void NimArr_map_2_allocatedMemory(const NimArr<nDim, T> &val, T2 **ans, int length) {
  if(val.isMap()) {
    NimArr<nDim, T2> target;
    vector<int> sizes(nDim);
    vector<int> strides(nDim);
    strides[0] = 1;
    for(int iii = 0; iii < nDim; ++iii) {
      sizes[iii] = val.dim()[iii];
      if(iii > 0)
	strides[iii] = strides[iii-1]*sizes[iii-1];
    }
    // just use dummy to be able to call setMap
    NimArr<nDim, T2> dummy;
    target.setMap(dummy, 0, strides, sizes);
    // then reset the Vptr to ans
    target.getVptrRef() = ans;
    target.mapCopy(val);
  } else {
    std::copy(val.getConstPtr(), val.getConstPtr() + length, *ans);
  }
}

template<int nDim, class T>
  bool isMapEntire(const NimArr<nDim, T> &v) {
  const int *s = v.strides();
  if(s[0] != 1) {return false;}
  if(v.getOffset() != 0) {return false;} // This could be relaxed, but we'd have to fix how functions handle this case.
  if(nDim == 1) {return true;}
  const int *d = v.dim();
  int nextStride = 1;
  for(unsigned int iii = 1; iii < nDim; ++iii) {
    nextStride *= d[iii-1];
    if(s[iii] != nextStride) {return false;}
  }
  return true;
}

// needed for some dists and some deriv functions
// Was previously in nimDists.cpp
template<int nDim, class T>
NimArr<nDim, T> &nimArrCopyIfNeeded(NimArr<nDim, T> &orig, NimArr<nDim, T> &possibleCopy) {
  if(orig.isMap()) {
    if(!isMapEntire<nDim, T>(orig)) {
      possibleCopy = orig;
      return(possibleCopy);
    }
  }
  return(orig);
}

// Here is the specialization for 1 dimensions (for any type, T = double, int or
// bool).
template <class T>
class NimArr<1, T> : public NimArrBase<T> {
 public:
  int size1;
  int calculateIndex(int i) const {
    return NimArrBase<T>::offset + NimArrBase<T>::stride1 * i;
  }
  int calculateIndex(vector<int> &i) const { return calculateIndex(i[0]); }
  T &operator()(int i) const {
    // could add asserts here
    return (*NimArrBase<T>::vPtr)[calculateIndex(i)];
  }

  ~NimArr<1, T>() {}

  template <class Tother>
  NimArr<1, T> &mapCopy(const NimArr<1, Tother> &other) {
    if (size1 != other.size1) {
      PRINTF("Error in mapCopy.  Sizes don't match: %i != %i \n", size1,
             other.size1);
    }
    T *to(*NimArrBase<T>::vPtr + NimArrBase<T>::offset);
    Tother *from(*other.vPtr + other.offset);
    int otherStride = other.stride1;
    for (int iii = 0; iii < size1; iii++) {
      *to = *from;
      to += NimArrBase<T>::stride1;
      from += otherStride;
    }
    return *this;
  }

  // This makes the copied-to object contiguous-memory.
  // Use mapCopy to copy into an existing map.
  template <class Tother>
  NimArr<1, T> &operator=(const NimArr<1, Tother> &other) {
    if (NimArrBase<T>::isMap()) {
      return mapCopy(other);
    }

    NimArrBase<T>::NAdims[0] = other.dim()[0];
    size1 = NimArrBase<T>::NAdims[0];

    NimArrBase<T>::NAlength = other.size();
    NimArrBase<T>::setMyType();

    NimArrBase<T>::boolMap = false;
    NimArrBase<T>::offset = 0;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    if (other.boolMap) {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(size1);
      NimArrBase<T>::own_v = true;
      T *to(NimArrBase<T>::v);
      T *toEnd(to + size1);
      Tother *from(*other.getVptr() + other.offset);

      int otherStride = other.stride1;
      for (; to != toEnd; to++) {
        *to = *from;
        from += otherStride;
      }
    } else {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(size1);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + size1, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
    return *this;
  }

  NimArr<1, T> &operator=(const NimArr<1, T> &other) {
    if (NimArrBase<T>::isMap()) {
      return (mapCopy(other));
    }

    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 1 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];

    NimArrBase<T>::NAlength = other.size();
    NimArrBase<T>::myType = other.getNimType();

    NimArrBase<T>::boolMap = false;
    NimArrBase<T>::offset = 0;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    if (other.boolMap) {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);

      int otherStride1 = other.stride1;
      for (int iRow = 0; iRow < size1; iRow++) {
        *to = *from;
        to++;
        from += otherStride1;
      }
    } else {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
    return (*this);
  }

  NimArr<1, T>(const NimArr<1, T> &other) : NimArrBase<T>(other) {
    NimArrBase<T>::NAdims[0] = other.dim()[0]; // redundant with base class copy constructor?
    size1 = NimArrBase<T>::NAdims[0];

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    if (other.boolMap) {
      NimArrBase<T>::v = nimble_malloc<T>(size1);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *toEnd(to + size1);
      T *from(*other.vPtr + other.offset);

      int otherStride = other.stride1;
      for (; to != toEnd; to++) {
        *to = *from;
        from += otherStride;
      }
    } else {
      NimArrBase<T>::v = nimble_malloc<T>(size1);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + size1, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
  }

  NimArr<1, T>() : NimArrBase<T>() { setSize(0); }

  void setMap(NimArrBase<T> &source, int off, int str1, int is1) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = NimArrBase<T>::NAlength = size1 = is1;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
  }

  void setMap(NimArrBase<T> &source, int off, vector<int> &str,
              vector<int> &is) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = NimArrBase<T>::NAlength = size1 = is[0];
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str[0];
  }

  template <class Tfrom>
  void dynamicMapCopy(int offset, vector<int> &str, vector<int> &is,
                      NimArrBase<Tfrom> *from, int fromOffset,
                      vector<int> &fromStr, vector<int> &fromIs) {
    if (NimArrBase<T>::isMap() || from->isMap()) {
      PRINTF("Error, dynamicMapCopy is not set up for nested maps\n");
    }
    NimArr<1, T> mapTo;
    mapTo.setMap(*this, offset, str, is);
    NimArr<1, Tfrom> mapFrom;
    mapFrom.setMap(*from, fromOffset, fromStr, fromIs);
    mapTo.mapCopy(mapFrom);
  }

  NimArr<1, T>(vector<T> &vm, int off, int str1, int is1)
      : NimArrBase<T>(vm, off) {
    NimArrBase<T>::NAdims[0] = NimArrBase<T>::NAlength = size1 = is1;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
  }

  NimArr<1, T>(int is1) : NimArrBase<T>() { setSize(is1); }

  void initialize(T value, bool init, int is1) {
    setSize(is1, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value);
    }
  }

  void initialize(T value, bool init, bool fillZeros, bool recycle, int is1) {
    setSize(is1, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value, fillZeros, recycle);
    }
  }

  void setSize(int is1, bool copyValues = true, bool fillZeros = true) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::setLength(size1, copyValues, fillZeros);
  }
  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) {
    setSize(sizeVec[0], copyValues, fillZeros);
  }
  virtual int numDims() const { return 1; }
  virtual int dimSize(int i) const {
    switch (i) {
      case 0:
        return size1;
        break;
      default:
        PRINTF("Error, incorrect dimension given to dimSize\n");
        return 0;
    }
  }
};

template <class T>
void dimNimArr(NimArr<1, int> &output, NimArrBase<T> &input) {
  output.setSize(input.numDims(), false, false);
  std::copy(input.dim(), input.dim() + input.numDims(), output.getPtr());
}

// Here is the specialization for 2 dimensions.

template <class T>
class NimArr<2, T> : public NimArrBase<T> {
 public:
  int size1, size2, stride2;
  int calculateIndex(int i, int j) const {
    return NimArrBase<T>::offset + NimArrBase<T>::stride1 * i + stride2 * j;
  }
  int calculateIndex(vector<int> &i) const {
    return calculateIndex(i[0], i[1]);
  }
  T &operator()(int i, int j) const {
    // could add asserts here
    return (*NimArrBase<T>::vPtr)[calculateIndex(i, j)];
  }

  T &operator[](int i) const {
    std::div_t divRes = std::div(i, size1);
    return (*NimArrBase<T>::vPtr)[calculateIndex(divRes.rem, divRes.quot)];
  }

  ~NimArr<2, T>() {}

  template <class Tother>
  NimArr<2, T> &mapCopy(const NimArr<2, Tother> &other) {
    if (size1 != other.size1) {
      PRINTF("Error in mapCopy.  Sizes 1 don't match: %i != %i \n", size1,
             other.size1);
    }
    if (size2 != other.size2) {
      PRINTF("Error in mapCopy.  Sizes 2 don't match: %i != %i \n", size2,
             other.size2);
    }

    T *to(*NimArrBase<T>::vPtr + NimArrBase<T>::offset);
    Tother *from(*other.vPtr + other.offset);

    int otherStride1 = other.stride1;
    int otherStride2 = other.stride2;
    for (int iCol = 0; iCol < size2; iCol++) {
      for (int iRow = 0; iRow < size1; iRow++) {
        *to = *from;
        to += NimArrBase<T>::stride1;
        from += otherStride1;
      }
      from += (-size1 * otherStride1) + otherStride2;
      to += (-size1 * NimArrBase<T>::stride1) + stride2;
    }
    return *this;
  }

  NimArr<2, T> &operator=(const NimArr<2, T> &other) {
    if (NimArrBase<T>::isMap()) {
      return mapCopy(other);
    }

    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 2 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];

    NimArrBase<T>::NAlength = other.size();
    NimArrBase<T>::myType = other.getNimType();

    NimArrBase<T>::boolMap = false;
    NimArrBase<T>::offset = 0;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    if (other.boolMap) {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);

      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      for (int iCol = 0; iCol < size2; iCol++) {
        for (int iRow = 0; iRow < size1; iRow++) {
          *to = *from;
          to++;
          from += otherStride1;
        }
        from += (-size1 * otherStride1) + otherStride2;
      }
    } else {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
    return *this;
  }

  NimArr<2, T>(const NimArr<2, T> &other) : NimArrBase<T>(other) {
    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 2 * sizeof(int));

    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    if (other.boolMap) {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);

      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      for (int iCol = 0; iCol < size2; iCol++) {
        for (int iRow = 0; iRow < size1; iRow++) {
          *to = *from;
          to++;
          from += otherStride1;
        }
        from += (-size1 * otherStride1) + otherStride2;
      }
    } else {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
  }

  NimArr<2, T>() : NimArrBase<T>() { setSize(0, 0); }

  void setMap(NimArrBase<T> &source, int off, int str1, int str2, int is1,
              int is2) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAlength = size1 * size2;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
  }

  void setMap(NimArrBase<T> &source, int off, vector<int> &str,
              vector<int> &is) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is[0];
    NimArrBase<T>::NAdims[1] = size2 = is[1];
    NimArrBase<T>::NAlength = size1 * size2;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str[0];
    NimArrBase<T>::NAstrides[1] = stride2 = str[1];
  }

  template <class Tfrom>
  void dynamicMapCopy(int offset, vector<int> &str, vector<int> &is,
                      NimArrBase<Tfrom> *from, int fromOffset,
                      vector<int> &fromStr, vector<int> &fromIs) {
    if (NimArrBase<T>::isMap() || from->isMap()) {
      PRINTF("Error, dynamicMapCopy is not set up for nested maps\n");
    }
    NimArr<2, T> mapTo;
    mapTo.setMap(*this, offset, str, is);
    NimArr<2, Tfrom> mapFrom;
    mapFrom.setMap(*from, fromOffset, fromStr, fromIs);
    mapTo.mapCopy(mapFrom);
  }

  NimArr<2, T>(vector<T> &vm, int off, int str1, int str2, int is1, int is2)
      : NimArrBase<T>(vm, off) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAlength = size1 * size2;
    // not setSize because it uses the allocated vm
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
  }

  NimArr<2, T>(int is1, int is2) : NimArrBase<T>() { setSize(is1, is2); }

  void initialize(T value, bool init, int is1, int is2) {
    setSize(is1, is2, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value);
    }
  }

  void initialize(T value, bool init, bool fillZeros, bool recycle, int is1,
                  int is2) {
    setSize(is1, is2, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value, fillZeros, recycle);
    }
  }

  void setSize(int is1, int is2, bool copyValues = true,
               bool fillZeros = true) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::setLength(size1 * size2, copyValues, fillZeros);
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = is1;
  }
  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) {
    setSize(sizeVec[0], sizeVec[1], copyValues, fillZeros);
  }
  virtual int numDims() const { return 2; }
  virtual int dimSize(int i) const {
    switch (i) {
      case 0:
        return size1;
        break;
      case 1:
        return size2;
        break;
      default:
        PRINTF("Error, incorrect dimension given to dimSize\n");
        return 0;
    }
  }
};

// Here is the specialization for 3 dimensions.

template <class T>
class NimArr<3, T> : public NimArrBase<T> {
 public:
  int size1, size2, size3, stride2, stride3;
  int calculateIndex(int i, int j, int k) const {
    return NimArrBase<T>::offset + NimArrBase<T>::stride1 * i + stride2 * j +
           stride3 * k;
  }
  int calculateIndex(vector<int> &i) const {
    return calculateIndex(i[0], i[1], i[2]);
  }
  T &operator()(int i, int j, int k) const {
    // could add asserts here
    return (*NimArrBase<T>::vPtr)[calculateIndex(i, j, k)];
  }

  T &operator[](int i) const {
    std::div_t divRes1 = std::div(i, size1);
    std::div_t divRes2 = std::div(divRes1.quot, size2);
    return (*NimArrBase<T>::vPtr)[calculateIndex(divRes1.rem, divRes2.rem,
                                                 divRes2.quot)];
  }

  ~NimArr<3, T>() {}

  template <class Tother>
  NimArr<3, T> &mapCopy(const NimArr<3, Tother> &other) {
    if (size1 != other.size1) {
      PRINTF("Error in mapCopy.  Sizes 1 don't match: %i != %i \n", size1,
             other.size1);
    }
    if (size2 != other.size2) {
      PRINTF("Error in mapCopy.  Sizes 2 don't match: %i != %i \n", size2,
             other.size2);
    }
    if (size3 != other.size3) {
      PRINTF("Error in mapCopy.  Sizes 3 don't match: %i != %i \n", size3,
             other.size3);
    }

    T *to(*NimArrBase<T>::vPtr + NimArrBase<T>::offset);
    Tother *from(*other.vPtr + other.offset);

    int otherStride1 = other.stride1;
    int otherStride2 = other.stride2;
    int otherStride3 = other.stride3;
    for (int i3rd = 0; i3rd < size3; i3rd++) {
      for (int iCol = 0; iCol < size2; iCol++) {
        for (int iRow = 0; iRow < size1; iRow++) {
          *to = *from;
          to += NimArrBase<T>::stride1;
          from += otherStride1;
        }
        from += (-size1 * otherStride1) + otherStride2;
        to += (-size1 * NimArrBase<T>::stride1) + stride2;
      }
      from += (-size2 * otherStride2) + otherStride3;
      to += (-size2 * stride2 + stride3);
    }
    return *this;
  }

  NimArr<3, T> &operator=(const NimArr<3, T> &other) {
    if (NimArrBase<T>::isMap()) {
      return mapCopy(other);
    }

    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 3 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];

    NimArrBase<T>::NAlength = other.size();
    NimArrBase<T>::myType = other.getNimType();

    NimArrBase<T>::boolMap = false;
    NimArrBase<T>::offset = 0;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    if (other.boolMap) {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);
      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      for (int i3rd = 0; i3rd < size3; i3rd++) {
        for (int iCol = 0; iCol < size2; iCol++) {
          for (int iRow = 0; iRow < size1; iRow++) {
            *to = *from;
            to++;
            from += otherStride1;
          }
          from += (-size1 * otherStride1) + otherStride2;
        }
        from += (-size2 * otherStride2) + otherStride3;
      }
    } else {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
    return *this;
  }

  NimArr<3, T>(const NimArr<3, T> &other) : NimArrBase<T>(other) {
    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 3 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    if (other.boolMap) {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);
      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      for (int i3rd = 0; i3rd < size3; i3rd++) {
        for (int iCol = 0; iCol < size2; iCol++) {
          for (int iRow = 0; iRow < size1; iRow++) {
            *to = *from;
            to++;
            from += otherStride1;
          }
          from += (-size1 * otherStride1) + otherStride2;
        }
        from += (-size2 * otherStride2) + otherStride3;
      }
    } else {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
  }

  NimArr<3, T>() : NimArrBase<T>() { setSize(0, 0, 0); }

  void setMap(NimArrBase<T> &source, int off, int str1, int str2, int str3,
              int is1, int is2, int is3) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;

    NimArrBase<T>::NAlength = size1 * size2 * size3;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;
  }

  void setMap(NimArrBase<T> &source, int off, vector<int> &str,
              vector<int> &is) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is[0];
    NimArrBase<T>::NAdims[1] = size2 = is[1];
    NimArrBase<T>::NAdims[2] = size3 = is[2];
    NimArrBase<T>::NAlength = size1 * size2 * size3;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str[0];
    NimArrBase<T>::NAstrides[1] = stride2 = str[1];
    NimArrBase<T>::NAstrides[2] = stride3 = str[2];
  }

  template <class Tfrom>
  void dynamicMapCopy(int offset, vector<int> &str, vector<int> &is,
                      NimArrBase<Tfrom> *from, int fromOffset,
                      vector<int> &fromStr, vector<int> &fromIs) {
    if (NimArrBase<T>::isMap() || from->isMap()) {
      PRINTF("Error, dynamicMapCopy is not set up for nested maps\n");
    }
    NimArr<3, T> mapTo;
    mapTo.setMap(*this, offset, str, is);
    NimArr<3, Tfrom> mapFrom;
    mapFrom.setMap(*from, fromOffset, fromStr, fromIs);
    mapTo.mapCopy(mapFrom);
  }

  NimArr<3, T>(vector<T> &vm, int off, int str1, int str2, int str3, int is1,
               int is2, int is3)
      : NimArrBase<T>(vm, off) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;

    NimArrBase<T>::NAlength = size1 * size2 * size3;
  }

  NimArr<3, T>(int is1, int is2, int is3) : NimArrBase<T>() {
    setSize(is1, is2, is3);
  }

  void initialize(T value, bool init, int is1, int is2, int is3) {
    setSize(is1, is2, is3, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value);
    }
  }

  void initialize(T value, bool init, bool fillZeros, bool recycle, int is1,
                  int is2, int is3) {
    setSize(is1, is2, is3, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value, fillZeros, recycle);
    }
  }

  void setSize(int is1, int is2, int is3, bool copyValues = true,
               bool fillZeros = true) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = is1;
    NimArrBase<T>::NAstrides[2] = stride3 = is1 * is2;
    NimArrBase<T>::setLength(stride3 * size3, copyValues, fillZeros);
  }

  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) {
    setSize(sizeVec[0], sizeVec[1], sizeVec[2], copyValues, fillZeros);
  }
  virtual int numDims() const { return 3; }
  virtual int dimSize(int i) const {
    switch (i) {
      case 0:
        return size1;
        break;
      case 1:
        return size2;
        break;
      case 2:
        return size3;
        break;
      default:
        PRINTF("Error, incorrect dimension given to dimSize\n");
        return 0;
    }
  }
};

// Here is the specialization for 4 dimensions.

template <class T>
class NimArr<4, T> : public NimArrBase<T> {
 public:
  int size1, size2, size3, size4, stride2, stride3, stride4;
  int calculateIndex(int i, int j, int k, int l) const {
    return NimArrBase<T>::offset + NimArrBase<T>::stride1 * i + stride2 * j +
           stride3 * k + stride4 * l;
  }
  int calculateIndex(vector<int> &i) const {
    return calculateIndex(i[0], i[1], i[2], i[3]);
  }
  T &operator()(int i, int j, int k, int l) const {
    // could add asserts here
    return (*NimArrBase<T>::vPtr)[calculateIndex(i, j, k, l)];
  }

  T &operator[](int i) const {
    std::div_t divRes1 = std::div(i, size1);
    std::div_t divRes2 = std::div(divRes1.quot, size2);
    std::div_t divRes3 = std::div(divRes2.quot, size3);
    return (*NimArrBase<T>::vPtr)[calculateIndex(divRes1.rem, divRes2.rem,
                                                 divRes3.rem, divRes3.quot)];
  }

  ~NimArr<4, T>() {}

  template <class Tother>
  NimArr<4, T> &mapCopy(const NimArr<4, Tother> &other) {
    if (size1 != other.size1) {
      PRINTF("Error in mapCopy.  Sizes 1 don't match: %i != %i \n", size1,
             other.size1);
    }
    if (size2 != other.size2) {
      PRINTF("Error in mapCopy.  Sizes 2 don't match: %i != %i \n", size2,
             other.size2);
    }
    if (size3 != other.size3) {
      PRINTF("Error in mapCopy.  Sizes 3 don't match: %i != %i \n", size3,
             other.size3);
    }
    if (size4 != other.size4) {
      PRINTF("Error in mapCopy.  Sizes 4 don't match: %i != %i \n", size4,
             other.size4);
    }
    T *to(*NimArrBase<T>::vPtr + NimArrBase<T>::offset);
    Tother *from(*other.vPtr + other.offset);

    int otherStride1 = other.stride1;
    int otherStride2 = other.stride2;
    int otherStride3 = other.stride3;
    int otherStride4 = other.stride4;
    for (int i4th = 0; i4th < size4; i4th++) {
      for (int i3rd = 0; i3rd < size3; i3rd++) {
        for (int iCol = 0; iCol < size2; iCol++) {
          for (int iRow = 0; iRow < size1; iRow++) {
            *to = *from;
            to += NimArrBase<T>::stride1;
            from += otherStride1;
          }
          from += (-size1 * otherStride1) + otherStride2;
          to += (-size1 * NimArrBase<T>::stride1) + stride2;
        }
        from += (-size2 * otherStride2) + otherStride3;
        to += (-size2 * stride2 + stride3);
      }
      from += (-size3 * otherStride3) + otherStride4;
      to += (-size3 * stride3 + stride4);
    }
    return *this;
  }

  NimArr<4, T> &operator=(const NimArr<4, T> &other) {
    if (NimArrBase<T>::isMap()) {
      return mapCopy(other);
    }

    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 4 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];
    size4 = NimArrBase<T>::NAdims[3];

    NimArrBase<T>::NAlength = other.size();
    NimArrBase<T>::myType = other.getNimType();

    NimArrBase<T>::boolMap = false;
    NimArrBase<T>::offset = 0;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    NimArrBase<T>::NAstrides[3] = stride4 = size1 * size2 * size3;
    if (other.boolMap) {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);

      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      int otherStride4 = other.stride4;
      for (int i4th = 0; i4th < size4; i4th++) {
        for (int i3rd = 0; i3rd < size3; i3rd++) {
          for (int iCol = 0; iCol < size2; iCol++) {
            for (int iRow = 0; iRow < size1; iRow++) {
              *to = *from;
              to++;
              from += otherStride1;
            }
            from += (-size1 * otherStride1) + otherStride2;
          }
          from += (-size2 * otherStride2) + otherStride3;
        }
        from += (-size3 * otherStride3) + otherStride4;
      }
    } else {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
    return *this;
  }

  NimArr<4, T>(const NimArr<4, T> &other) : NimArrBase<T>(other) {
    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 4 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];
    size4 = NimArrBase<T>::NAdims[3];

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    NimArrBase<T>::NAstrides[3] = stride4 = size1 * size2 * size3;
    if (other.boolMap) {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);
      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      int otherStride4 = other.stride4;
      for (int i4th = 0; i4th < size4; i4th++) {
        for (int i3rd = 0; i3rd < size3; i3rd++) {
          for (int iCol = 0; iCol < size2; iCol++) {
            for (int iRow = 0; iRow < size1; iRow++) {
              *to = *from;
              to++;
              from += otherStride1;
            }
            from += (-size1 * otherStride1) + otherStride2;
          }
          from += (-size2 * otherStride2) + otherStride3;
        }
        from += (-size3 * otherStride3) + otherStride4;
      }

    } else {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
  }

  NimArr<4, T>() : NimArrBase<T>() { setSize(0, 0, 0, 0); }

  void setMap(NimArrBase<T> &source, int off, int str1, int str2, int str3,
              int str4, int is1, int is2, int is3, int is4) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;
    NimArrBase<T>::NAstrides[3] = stride4 = str4;
  }

  void setMap(NimArrBase<T> &source, int off, vector<int> &str,
              vector<int> &is) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is[0];
    NimArrBase<T>::NAdims[1] = size2 = is[1];
    NimArrBase<T>::NAdims[2] = size3 = is[2];
    NimArrBase<T>::NAdims[3] = size4 = is[3];

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str[0];
    NimArrBase<T>::NAstrides[1] = stride2 = str[1];
    NimArrBase<T>::NAstrides[2] = stride3 = str[2];
    NimArrBase<T>::NAstrides[3] = stride4 = str[3];
  }

  template <class Tfrom>
  void dynamicMapCopy(int offset, vector<int> &str, vector<int> &is,
                      NimArrBase<Tfrom> *from, int fromOffset,
                      vector<int> &fromStr, vector<int> &fromIs) {
    if (NimArrBase<T>::isMap() || from->isMap()) {
      PRINTF("Error, dynamicMapCopy is not set up for nested maps\n");
    }
    NimArr<4, T> mapTo;
    mapTo.setMap(*this, offset, str, is);
    NimArr<4, Tfrom> mapFrom;
    mapFrom.setMap(*from, fromOffset, fromStr, fromIs);
    mapTo.mapCopy(mapFrom);
  }

  NimArr<4, T>(vector<T> &vm, int off, int str1, int str2, int str3, int str4,
               int is1, int is2, int is3, int is4)
      : NimArrBase<T>(vm, off) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;
    NimArrBase<T>::NAstrides[3] = stride4 = str4;

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4;
  }

  NimArr<4, T>(int is1, int is2, int is3, int is4) : NimArrBase<T>() {
    setSize(is1, is2, is3, is4);
  }

  void initialize(T value, bool init, int is1, int is2, int is3, int is4) {
    setSize(is1, is2, is3, is4, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value);
    }
  }

  void initialize(T value, bool init, bool fillZeros, bool recycle, int is1,
                  int is2, int is3, int is4) {
    setSize(is1, is2, is3, is4, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value, fillZeros, recycle);
    }
  }

  void setSize(int is1, int is2, int is3, int is4, bool copyValues = true,
               bool fillZeros = true) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = is1;
    NimArrBase<T>::NAstrides[2] = stride3 = is1 * is2;
    NimArrBase<T>::NAstrides[3] = stride4 = is1 * is2 * is3;
    NimArrBase<T>::setLength(stride4 * size4, copyValues, fillZeros);
  }

  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) {
    setSize(sizeVec[0], sizeVec[1], sizeVec[2], sizeVec[3], copyValues,
            fillZeros);
  }
  virtual int numDims() const { return 4; }
  virtual int dimSize(int i) const {
    switch (i) {
      case 0:
        return size1;
        break;
      case 1:
        return size2;
        break;
      case 2:
        return size3;
        break;
      case 3:
        return size4;
        break;
      default:
        PRINTF("Error, incorrect dimension given to dimSize\n");
        return 0;
    }
  }
};


// Here is the specialization for 5 dimensions.

template <class T>
class NimArr<5, T> : public NimArrBase<T> {
 public:
  int size1, size2, size3, size4, size5, stride2, stride3, stride4, stride5;
  int calculateIndex(int i, int j, int k, int l, int m) const {
    return NimArrBase<T>::offset + NimArrBase<T>::stride1 * i + stride2 * j +
           stride3 * k + stride4 * l + stride5 * m;
  }
  int calculateIndex(vector<int> &i) const {
    return calculateIndex(i[0], i[1], i[2], i[3], i[4]);
  }
  T &operator()(int i, int j, int k, int l, int m) const {
    // could add asserts here
    return (*NimArrBase<T>::vPtr)[calculateIndex(i, j, k, l, m)];
  }

  T &operator[](int i) const {
    std::div_t divRes1 = std::div(i, size1);
    std::div_t divRes2 = std::div(divRes1.quot, size2);
    std::div_t divRes3 = std::div(divRes2.quot, size3);
    std::div_t divRes4 = std::div(divRes3.quot, size4);
    return (*NimArrBase<T>::vPtr)[calculateIndex(divRes1.rem, divRes2.rem,
                                                 divRes3.rem, divRes4.rem,
						 divRes4.quot)];
  }

  ~NimArr<5, T>() {}

  template <class Tother>
  NimArr<5, T> &mapCopy(const NimArr<5, Tother> &other) {
    if (size1 != other.size1) {
      PRINTF("Error in mapCopy.  Sizes 1 don't match: %i != %i \n", size1,
             other.size1);
    }
    if (size2 != other.size2) {
      PRINTF("Error in mapCopy.  Sizes 2 don't match: %i != %i \n", size2,
             other.size2);
    }
    if (size3 != other.size3) {
      PRINTF("Error in mapCopy.  Sizes 3 don't match: %i != %i \n", size3,
             other.size3);
    }
    if (size4 != other.size4) {
      PRINTF("Error in mapCopy.  Sizes 4 don't match: %i != %i \n", size4,
             other.size4);
    }
    if (size5 != other.size5) {
      PRINTF("Error in mapCopy.  Sizes 5 don't match: %i != %i \n", size5,
             other.size5);
    }
    T *to(*NimArrBase<T>::vPtr + NimArrBase<T>::offset);
    Tother *from(*other.vPtr + other.offset);

    int otherStride1 = other.stride1;
    int otherStride2 = other.stride2;
    int otherStride3 = other.stride3;
    int otherStride4 = other.stride4;
    int otherStride5 = other.stride5;
    for (int i5th = 0; i5th < size5; i5th++) {
      for (int i4th = 0; i4th < size4; i4th++) {
        for (int i3rd = 0; i3rd < size3; i3rd++) {
          for (int iCol = 0; iCol < size2; iCol++) {
            for (int iRow = 0; iRow < size1; iRow++) {
              *to = *from;
              to += NimArrBase<T>::stride1;
              from += otherStride1;
            }
            from += (-size1 * otherStride1) + otherStride2;
            to += (-size1 * NimArrBase<T>::stride1) + stride2;
          }
          from += (-size2 * otherStride2) + otherStride3;
          to += (-size2 * stride2 + stride3);
        }
        from += (-size3 * otherStride3) + otherStride4;
        to += (-size3 * stride3 + stride4);
      }
      from += (-size4 * otherStride4) + otherStride5;
      to += (-size4 * stride4 + stride5);
    }
    return *this;
  }

  NimArr<5, T> &operator=(const NimArr<5, T> &other) {
    if (NimArrBase<T>::isMap()) {
      return mapCopy(other);
    }

    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 5 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];
    size4 = NimArrBase<T>::NAdims[3];
    size5 = NimArrBase<T>::NAdims[4];

    NimArrBase<T>::NAlength = other.size();
    NimArrBase<T>::myType = other.getNimType();

    NimArrBase<T>::boolMap = false;
    NimArrBase<T>::offset = 0;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    NimArrBase<T>::NAstrides[3] = stride4 = size1 * size2 * size3;
    NimArrBase<T>::NAstrides[4] = stride5 = size1 * size2 * size3 * size4;
    if (other.boolMap) {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);

      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      int otherStride4 = other.stride4;
      int otherStride5 = other.stride5;
      for (int i5th = 0; i5th < size5; i5th++) {
        for (int i4th = 0; i4th < size4; i4th++) {
          for (int i3rd = 0; i3rd < size3; i3rd++) {
            for (int iCol = 0; iCol < size2; iCol++) {
              for (int iRow = 0; iRow < size1; iRow++) {
                *to = *from;
                to++;
                from += otherStride1;
              }
              from += (-size1 * otherStride1) + otherStride2;
            }
            from += (-size2 * otherStride2) + otherStride3;
          }
          from += (-size3 * otherStride3) + otherStride4;
        }
        from += (-size4 * otherStride4) + otherStride5;
      }
    } else {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
    return *this;
  }

  NimArr<5, T>(const NimArr<5, T> &other) : NimArrBase<T>(other) {
    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 5 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];
    size4 = NimArrBase<T>::NAdims[3];
    size5 = NimArrBase<T>::NAdims[4];

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    NimArrBase<T>::NAstrides[3] = stride4 = size1 * size2 * size3;
    NimArrBase<T>::NAstrides[4] = stride5 = size1 * size2 * size3 * size4;
    if (other.boolMap) {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);
      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      int otherStride4 = other.stride4;
      int otherStride5 = other.stride5;
      for (int i5th = 0; i5th < size5; i5th++) {
        for (int i4th = 0; i4th < size4; i4th++) {
          for (int i3rd = 0; i3rd < size3; i3rd++) {
            for (int iCol = 0; iCol < size2; iCol++) {
              for (int iRow = 0; iRow < size1; iRow++) {
                *to = *from;
                to++;
                from += otherStride1;
              }
              from += (-size1 * otherStride1) + otherStride2;
            }
            from += (-size2 * otherStride2) + otherStride3;
          }
          from += (-size3 * otherStride3) + otherStride4;
        }
	from += (-size4 * otherStride4) + otherStride5;
      }

    } else {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
  }

  NimArr<5, T>() : NimArrBase<T>() { setSize(0, 0, 0, 0, 0); }

  void setMap(NimArrBase<T> &source, int off, int str1, int str2, int str3,
              int str4, int str5, int is1, int is2, int is3, int is4, int is5) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;
    NimArrBase<T>::NAdims[4] = size5 = is5;

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4 * size5;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;
    NimArrBase<T>::NAstrides[3] = stride4 = str4;
    NimArrBase<T>::NAstrides[4] = stride5 = str5;
  }

  void setMap(NimArrBase<T> &source, int off, vector<int> &str,
              vector<int> &is) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is[0];
    NimArrBase<T>::NAdims[1] = size2 = is[1];
    NimArrBase<T>::NAdims[2] = size3 = is[2];
    NimArrBase<T>::NAdims[3] = size4 = is[3];
    NimArrBase<T>::NAdims[4] = size5 = is[4];

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4 * size5;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str[0];
    NimArrBase<T>::NAstrides[1] = stride2 = str[1];
    NimArrBase<T>::NAstrides[2] = stride3 = str[2];
    NimArrBase<T>::NAstrides[3] = stride4 = str[3];
    NimArrBase<T>::NAstrides[4] = stride5 = str[4];
  }

  template <class Tfrom>
  void dynamicMapCopy(int offset, vector<int> &str, vector<int> &is,
                      NimArrBase<Tfrom> *from, int fromOffset,
                      vector<int> &fromStr, vector<int> &fromIs) {
    if (NimArrBase<T>::isMap() || from->isMap()) {
      PRINTF("Error, dynamicMapCopy is not set up for nested maps\n");
    }
    NimArr<5, T> mapTo;
    mapTo.setMap(*this, offset, str, is);
    NimArr<5, Tfrom> mapFrom;
    mapFrom.setMap(*from, fromOffset, fromStr, fromIs);
    mapTo.mapCopy(mapFrom);
  }

  NimArr<5, T>(vector<T> &vm, int off, int str1, int str2, int str3, int str4, int str5,
               int is1, int is2, int is3, int is4, int is5)
      : NimArrBase<T>(vm, off) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;
    NimArrBase<T>::NAdims[4] = size5 = is5;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;
    NimArrBase<T>::NAstrides[3] = stride4 = str4;
    NimArrBase<T>::NAstrides[4] = stride5 = str5;

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4 * size5;
  }

  NimArr<5, T>(int is1, int is2, int is3, int is4, int is5) : NimArrBase<T>() {
    setSize(is1, is2, is3, is4, is5);
  }

  void initialize(T value, bool init, int is1, int is2, int is3, int is4, int is5) {
    setSize(is1, is2, is3, is4, is5, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value);
    }
  }

  void initialize(T value, bool init, bool fillZeros, bool recycle, int is1,
                  int is2, int is3, int is4, int is5) {
    setSize(is1, is2, is3, is4, is5, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value, fillZeros, recycle);
    }
  }

  void setSize(int is1, int is2, int is3, int is4, int is5, bool copyValues = true,
               bool fillZeros = true) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;
    NimArrBase<T>::NAdims[4] = size5 = is5;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = is1;
    NimArrBase<T>::NAstrides[2] = stride3 = is1 * is2;
    NimArrBase<T>::NAstrides[3] = stride4 = is1 * is2 * is3;
    NimArrBase<T>::NAstrides[4] = stride5 = is1 * is2 * is3 * is4;
    NimArrBase<T>::setLength(stride5 * size5, copyValues, fillZeros);
  }

  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) {
    setSize(sizeVec[0], sizeVec[1], sizeVec[2], sizeVec[3], sizeVec[4], copyValues,
            fillZeros);
  }
  virtual int numDims() const { return 5; }
  virtual int dimSize(int i) const {
    switch (i) {
      case 0:
        return size1;
        break;
      case 1:
        return size2;
        break;
      case 2:
        return size3;
        break;
      case 3:
        return size4;
        break;
      case 4:
        return size5;
        break;
      default:
        PRINTF("Error, incorrect dimension given to dimSize\n");
        return 0;
    }
  }
};




// Here is the specialization for 6 dimensions.

template <class T>
class NimArr<6, T> : public NimArrBase<T> {
 public:
  int size1, size2, size3, size4, size5, size6, stride2, stride3, stride4, stride5, stride6;
  int calculateIndex(int i, int j, int k, int l, int m, int n) const {
    return NimArrBase<T>::offset + NimArrBase<T>::stride1 * i + stride2 * j +
           stride3 * k + stride4 * l + stride5 * m + stride6 * n;
  }
  int calculateIndex(vector<int> &i) const {
    return calculateIndex(i[0], i[1], i[2], i[3], i[4], i[5]);
  }
  T &operator()(int i, int j, int k, int l, int m, int n) const {
    // could add asserts here
    return (*NimArrBase<T>::vPtr)[calculateIndex(i, j, k, l, m, n)];
  }

  T &operator[](int i) const {
    std::div_t divRes1 = std::div(i, size1);
    std::div_t divRes2 = std::div(divRes1.quot, size2);
    std::div_t divRes3 = std::div(divRes2.quot, size3);
    std::div_t divRes4 = std::div(divRes3.quot, size4);
    std::div_t divRes5 = std::div(divRes4.quot, size5);
    return (*NimArrBase<T>::vPtr)[calculateIndex(divRes1.rem, divRes2.rem,
                                                 divRes3.rem, divRes4.rem,
						 divRes5.rem, divRes5.quot)];
  }

  ~NimArr<6, T>() {}

  template <class Tother>
  NimArr<6, T> &mapCopy(const NimArr<6, Tother> &other) {
    if (size1 != other.size1) {
      PRINTF("Error in mapCopy.  Sizes 1 don't match: %i != %i \n", size1,
             other.size1);
    }
    if (size2 != other.size2) {
      PRINTF("Error in mapCopy.  Sizes 2 don't match: %i != %i \n", size2,
             other.size2);
    }
    if (size3 != other.size3) {
      PRINTF("Error in mapCopy.  Sizes 3 don't match: %i != %i \n", size3,
             other.size3);
    }
    if (size4 != other.size4) {
      PRINTF("Error in mapCopy.  Sizes 4 don't match: %i != %i \n", size4,
             other.size4);
    }
    if (size5 != other.size5) {
      PRINTF("Error in mapCopy.  Sizes 5 don't match: %i != %i \n", size5,
             other.size5);
    }
    if (size6 != other.size6) {
      PRINTF("Error in mapCopy.  Sizes 6 don't match: %i != %i \n", size6,
             other.size6);
    }
    T *to(*NimArrBase<T>::vPtr + NimArrBase<T>::offset);
    Tother *from(*other.vPtr + other.offset);

    int otherStride1 = other.stride1;
    int otherStride2 = other.stride2;
    int otherStride3 = other.stride3;
    int otherStride4 = other.stride4;
    int otherStride5 = other.stride5;
    int otherStride6 = other.stride6;
    for (int i6th = 0; i6th < size6; i6th++) {
      for (int i5th = 0; i5th < size5; i5th++) {
        for (int i4th = 0; i4th < size4; i4th++) {
          for (int i3rd = 0; i3rd < size3; i3rd++) {
            for (int iCol = 0; iCol < size2; iCol++) {
              for (int iRow = 0; iRow < size1; iRow++) {
                *to = *from;
                to += NimArrBase<T>::stride1;
                from += otherStride1;
              }
              from += (-size1 * otherStride1) + otherStride2;
              to += (-size1 * NimArrBase<T>::stride1) + stride2;
            }
            from += (-size2 * otherStride2) + otherStride3;
            to += (-size2 * stride2 + stride3);
          }
          from += (-size3 * otherStride3) + otherStride4;
          to += (-size3 * stride3 + stride4);
        }
        from += (-size4 * otherStride4) + otherStride5;
        to += (-size4 * stride4 + stride5);
      }
      from += (-size5 * otherStride5) + otherStride6;
      to += (-size5 * stride5 + stride6);
    }
    return *this;
  }

  NimArr<6, T> &operator=(const NimArr<6, T> &other) {
    if (NimArrBase<T>::isMap()) {
      return mapCopy(other);
    }

    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 6 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];
    size4 = NimArrBase<T>::NAdims[3];
    size5 = NimArrBase<T>::NAdims[4];
    size6 = NimArrBase<T>::NAdims[5];

    NimArrBase<T>::NAlength = other.size();
    NimArrBase<T>::myType = other.getNimType();

    NimArrBase<T>::boolMap = false;
    NimArrBase<T>::offset = 0;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    NimArrBase<T>::NAstrides[3] = stride4 = size1 * size2 * size3;
    NimArrBase<T>::NAstrides[4] = stride5 = size1 * size2 * size3 * size4;
    NimArrBase<T>::NAstrides[5] = stride6 = size1 * size2 * size3 * size4 * size5;
    if (other.boolMap) {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);

      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      int otherStride4 = other.stride4;
      int otherStride5 = other.stride5;
      int otherStride6 = other.stride6;
      for (int i6th = 0; i6th < size6; i6th++) {
        for (int i5th = 0; i5th < size5; i5th++) {
          for (int i4th = 0; i4th < size4; i4th++) {
            for (int i3rd = 0; i3rd < size3; i3rd++) {
              for (int iCol = 0; iCol < size2; iCol++) {
                for (int iRow = 0; iRow < size1; iRow++) {
                  *to = *from;
                  to++;
                  from += otherStride1;
                }
                from += (-size1 * otherStride1) + otherStride2;
              }
              from += (-size2 * otherStride2) + otherStride3;
            }
            from += (-size3 * otherStride3) + otherStride4;
          }
          from += (-size4 * otherStride4) + otherStride5;
        }
        from += (-size5 * otherStride5) + otherStride6;
      }
    } else {
      if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
    return *this;
  }

  NimArr<6, T>(const NimArr<6, T> &other) : NimArrBase<T>(other) {
    std::memcpy(NimArrBase<T>::NAdims, other.dim(), 6 * sizeof(int));
    size1 = NimArrBase<T>::NAdims[0];
    size2 = NimArrBase<T>::NAdims[1];
    size3 = NimArrBase<T>::NAdims[2];
    size4 = NimArrBase<T>::NAdims[3];
    size5 = NimArrBase<T>::NAdims[4];
    size6 = NimArrBase<T>::NAdims[5];

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = size1;
    NimArrBase<T>::NAstrides[2] = stride3 = size1 * size2;
    NimArrBase<T>::NAstrides[3] = stride4 = size1 * size2 * size3;
    NimArrBase<T>::NAstrides[4] = stride5 = size1 * size2 * size3 * size4;
    NimArrBase<T>::NAstrides[5] = stride6 = size1 * size2 * size3 * size4 * size5;
    if (other.boolMap) {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;

      T *to(NimArrBase<T>::v);
      T *from(*other.vPtr + other.offset);
      int otherStride1 = other.stride1;
      int otherStride2 = other.stride2;
      int otherStride3 = other.stride3;
      int otherStride4 = other.stride4;
      int otherStride5 = other.stride5;
      int otherStride6 = other.stride6;
      for (int i6th = 0; i6th < size6; i6th++) {
        for (int i5th = 0; i5th < size5; i5th++) {
          for (int i4th = 0; i4th < size4; i4th++) {
            for (int i3rd = 0; i3rd < size3; i3rd++) {
              for (int iCol = 0; iCol < size2; iCol++) {
                for (int iRow = 0; iRow < size1; iRow++) {
                  *to = *from;
                  to++;
                  from += otherStride1;
                }
                from += (-size1 * otherStride1) + otherStride2;
              }
              from += (-size2 * otherStride2) + otherStride3;
            }
            from += (-size3 * otherStride3) + otherStride4;
          }
          from += (-size4 * otherStride4) + otherStride5;
        }
        from += (-size5 * otherStride5) + otherStride6;
      }

    } else {
      NimArrBase<T>::v = nimble_malloc<T>(NimArrBase<T>::NAlength);
      NimArrBase<T>::own_v = true;
      std::copy(other.v, other.v + NimArrBase<T>::NAlength, NimArrBase<T>::v);
    }
    NimArrBase<T>::setVptr();
  }

  NimArr<6, T>() : NimArrBase<T>() { setSize(0, 0, 0, 0, 0, 0); }

  void setMap(NimArrBase<T> &source, int off,
	      int str1, int str2, int str3, int str4, int str5, int str6,
	      int is1, int is2, int is3, int is4, int is5, int is6) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;
    NimArrBase<T>::NAdims[4] = size5 = is5;
    NimArrBase<T>::NAdims[5] = size6 = is6;

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4 * size5 * size6;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;
    NimArrBase<T>::NAstrides[3] = stride4 = str4;
    NimArrBase<T>::NAstrides[4] = stride5 = str5;
    NimArrBase<T>::NAstrides[5] = stride6 = str6;
  }

  void setMap(NimArrBase<T> &source, int off, vector<int> &str,
              vector<int> &is) {
    if (NimArrBase<T>::own_v) nimble_free(NimArrBase<T>::v);
    NimArrBase<T>::boolMap = true;
    NimArrBase<T>::offset = off;
    NimArrBase<T>::vPtr = source.getVptr();
    NimArrBase<T>::own_v = false;
    NimArrBase<T>::NAdims[0] = size1 = is[0];
    NimArrBase<T>::NAdims[1] = size2 = is[1];
    NimArrBase<T>::NAdims[2] = size3 = is[2];
    NimArrBase<T>::NAdims[3] = size4 = is[3];
    NimArrBase<T>::NAdims[4] = size5 = is[4];
    NimArrBase<T>::NAdims[5] = size6 = is[5];

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4 * size5 * size6;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str[0];
    NimArrBase<T>::NAstrides[1] = stride2 = str[1];
    NimArrBase<T>::NAstrides[2] = stride3 = str[2];
    NimArrBase<T>::NAstrides[3] = stride4 = str[3];
    NimArrBase<T>::NAstrides[4] = stride5 = str[4];
    NimArrBase<T>::NAstrides[5] = stride6 = str[5];
  }

  template <class Tfrom>
  void dynamicMapCopy(int offset, vector<int> &str, vector<int> &is,
                      NimArrBase<Tfrom> *from, int fromOffset,
                      vector<int> &fromStr, vector<int> &fromIs) {
    if (NimArrBase<T>::isMap() || from->isMap()) {
      PRINTF("Error, dynamicMapCopy is not set up for nested maps\n");
    }
    NimArr<6, T> mapTo;
    mapTo.setMap(*this, offset, str, is);
    NimArr<6, Tfrom> mapFrom;
    mapFrom.setMap(*from, fromOffset, fromStr, fromIs);
    mapTo.mapCopy(mapFrom);
  }

  NimArr<6, T>(vector<T> &vm, int off,
	       int str1, int str2, int str3, int str4, int str5, int str6,
               int is1, int is2, int is3, int is4, int is5, int is6)
      : NimArrBase<T>(vm, off) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;
    NimArrBase<T>::NAdims[4] = size5 = is5;
    NimArrBase<T>::NAdims[5] = size6 = is6;

    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = str1;
    NimArrBase<T>::NAstrides[1] = stride2 = str2;
    NimArrBase<T>::NAstrides[2] = stride3 = str3;
    NimArrBase<T>::NAstrides[3] = stride4 = str4;
    NimArrBase<T>::NAstrides[4] = stride5 = str5;
    NimArrBase<T>::NAstrides[5] = stride6 = str6;

    NimArrBase<T>::NAlength = size1 * size2 * size3 * size4 * size5 * size6;
  }

  NimArr<6, T>(int is1, int is2, int is3, int is4, int is5, int is6) : NimArrBase<T>() {
    setSize(is1, is2, is3, is4, is5, is6);
  }

  void initialize(T value, bool init, int is1, int is2, int is3, int is4, int is5, int is6) {
    setSize(is1, is2, is3, is4, is5, is6, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value);
    }
  }

  void initialize(T value, bool init, bool fillZeros, bool recycle, int is1,
                  int is2, int is3, int is4, int is5, int is6) {
    setSize(is1, is2, is3, is4, is5, is6, false, false);
    if (init) {
      NimArrBase<T>::fillAllValues(value, fillZeros, recycle);
    }
  }

  void setSize(int is1, int is2, int is3, int is4, int is5, int is6,
	       bool copyValues = true, bool fillZeros = true) {
    NimArrBase<T>::NAdims[0] = size1 = is1;
    NimArrBase<T>::NAdims[1] = size2 = is2;
    NimArrBase<T>::NAdims[2] = size3 = is3;
    NimArrBase<T>::NAdims[3] = size4 = is4;
    NimArrBase<T>::NAdims[4] = size5 = is5;
    NimArrBase<T>::NAdims[5] = size6 = is6;
    NimArrBase<T>::NAstrides[0] = NimArrBase<T>::stride1 = 1;
    NimArrBase<T>::NAstrides[1] = stride2 = is1;
    NimArrBase<T>::NAstrides[2] = stride3 = is1 * is2;
    NimArrBase<T>::NAstrides[3] = stride4 = is1 * is2 * is3;
    NimArrBase<T>::NAstrides[4] = stride5 = is1 * is2 * is3 * is4;
    NimArrBase<T>::NAstrides[5] = stride6 = is1 * is2 * is3 * is4 * is5;
    NimArrBase<T>::setLength(stride6 * size6, copyValues, fillZeros);
  }

  virtual void setSize(vector<int> sizeVec, bool copyValues = true,
                       bool fillZeros = true) {
    setSize(sizeVec[0], sizeVec[1], sizeVec[2], sizeVec[3], sizeVec[4], sizeVec[5],
	    copyValues, fillZeros);
  }
  virtual int numDims() const { return 6; }
  virtual int dimSize(int i) const {
    switch (i) {
      case 0:
        return size1;
        break;
      case 1:
        return size2;
        break;
      case 2:
        return size3;
        break;
      case 3:
        return size4;
        break;
      case 4:
        return size5;
        break;
      case 5:
        return size6;
        break;
      default:
        PRINTF("Error, incorrect dimension given to dimSize\n");
        return 0;
    }
  }
};






////////////////////////////////////
// VecNimArr
///////////////////////////////////

template <int ndim, class T>
class VecNimArr : public VecNimArrBase<T> {
 public:
  ~VecNimArr<ndim, T>() {}
  std::vector<NimArr<ndim, T> > values;
  NimArr<ndim, T> &operator[](unsigned int i) {
    if (i >= values.size()) {
      PRINTF(
          "Error accessing a VecNimArr element: requested element %i but "
          "values.size() is only %i\n",
          i, values.size());
      PRINTF(
          "Returning the first element if possible to avoid a segfault "
          "crash\n");
      return values[0];
    }
    return values[i];
  }
  virtual void resize(int i) { values.resize(i); }
  void resizeNoPtr(int i) { values.resize(i); }
  int getsizeNoPtr() { return values.size(); }
  virtual NimArrBase<T> *getBasePtr(int i) { return &(values[i]); }
  virtual int size() { return values.size(); }

  virtual void setRowDims(int row, vector<int> dims) {
    if (dims.size() != ndim) {
      PRINTF(
          "Error: number of dimensions incorrect in resize of numericList\n");
      return;
    }
    NimArrBase<T> *nimBasePtr = VecNimArr<ndim, T>::getBasePtr(row);
    (*nimBasePtr).setSize(dims);
    return;
  }

  virtual vector<int> getRowDims(int row) {
    if (row >= size()) return vector<int>(0);
    NimArrBase<T> *nimBasePtr = VecNimArr<ndim, T>::getBasePtr(row);
    int nRowDims = (*nimBasePtr).numDims();
    vector<int> rowDims(nRowDims);
    for (int i = 0; i < nRowDims; i++) rowDims[i] = (*nimBasePtr).dimSize(i);
    return rowDims;
  }

  void setSize(int row, int d1) {
    if (ndim != 1) {
      PRINTF(
          "Error: number of dimensions incorrect in resize of numericList\n");
      return;
    }
    vector<int> dims;
    dims.resize(1);
    dims[0] = d1;
    NimArrBase<T> *nimBasePtr = VecNimArr<ndim, T>::getBasePtr(row);
    (*nimBasePtr).setSize(dims);
    return;
  }

  void setSize(int row, int d1, int d2) {
    if (ndim != 2) {
      PRINTF(
          "Error: number of dimensions incorrect in resize of numericList\n");
      return;
    }
    vector<int> dims;
    dims.resize(2);
    dims[0] = d1;
    dims[1] = d2;
    NimArrBase<T> *nimBasePtr = VecNimArr<ndim, T>::getBasePtr(row);
    (*nimBasePtr).setSize(dims);
    return;
  }

  void setSize(int row, int d1, int d2, int d3) {
    if (ndim != 3) {
      PRINTF(
          "Error: number of dimensions incorrect in resize of numericList\n");
      return;
    }
    vector<int> dims;
    dims.resize(3);
    dims[0] = d1;
    dims[1] = d2;
    dims[2] = d3;
    NimArrBase<T> *nimBasePtr = VecNimArr<ndim, T>::getBasePtr(row);
    (*nimBasePtr).setSize(dims);
    return;
  }
};

template <class Tfrom, class Tto, int mapDim>
void dynamicMapCopyDim(NimArrType *toNimArr, int toOffset, vector<int> &toStr,
                       vector<int> &toIs, NimArrType *fromNimArr,
                       int fromOffset, vector<int> &fromStr,
                       vector<int> &fromIs) {
  NimArr<mapDim, Tfrom> mapFrom;
  mapFrom.setMap(*static_cast<NimArrBase<Tfrom> *>(fromNimArr), fromOffset,
                 fromStr, fromIs);
  NimArr<mapDim, Tto> mapTo;
  mapTo.setMap(*static_cast<NimArrBase<Tto> *>(toNimArr), toOffset, toStr,
               toIs);
  mapTo.mapCopy(mapFrom);
}

template <class Tfrom, class Tto>
void dynamicMapCopy(NimArrType *toNimArr, int toOffset, vector<int> &toStr,
                    vector<int> &toIs, NimArrType *fromNimArr, int fromOffset,
                    vector<int> &fromStr, vector<int> &fromIs) {
  int mapDim = toStr.size();
  // must be the same as fromStr.sizes();
  if (static_cast<NimArrBase<Tfrom> *>(fromNimArr)->isMap() ||
      static_cast<NimArrBase<Tto> *>(toNimArr)->isMap()) {
    PRINTF("Error, dynamicMapCopy is not set up for nested maps\n");
  }
  switch (mapDim) {
    case 1:
      dynamicMapCopyDim<Tfrom, Tto, 1>(toNimArr, toOffset, toStr, toIs,
                                       fromNimArr, fromOffset, fromStr, fromIs);
      break;
    case 2:
      dynamicMapCopyDim<Tfrom, Tto, 2>(toNimArr, toOffset, toStr, toIs,
                                       fromNimArr, fromOffset, fromStr, fromIs);
      break;
    case 3:
      dynamicMapCopyDim<Tfrom, Tto, 3>(toNimArr, toOffset, toStr, toIs,
                                       fromNimArr, fromOffset, fromStr, fromIs);
      break;
    case 4:
      dynamicMapCopyDim<Tfrom, Tto, 4>(toNimArr, toOffset, toStr, toIs,
                                       fromNimArr, fromOffset, fromStr, fromIs);
      break;
    default:
      PRINTF(
          "Error in copying (dynamicMapCopy): more than 4 dimensions not "
          "supported yet\n");
  }
}

template <class Tfrom, class Tto, int mapDim>
void dynamicMapCopyFlatToDimFixed(NimArrBase<Tto> *toNimArr, int toOffset,
                                  vector<int> &toStr, vector<int> &toIs,
                                  NimArrBase<Tfrom> *fromNimArr, int fromOffset,
                                  int fromStr) {
  NimArr<mapDim, Tfrom> mapFrom;
  vector<int> fromStrVec(mapDim);
  fromStrVec[0] = fromStr;
  for (int i = 1; i < mapDim; i++) {
    fromStrVec[i] = toIs[i - 1] * fromStrVec[i - 1];
  }
  mapFrom.setMap(*fromNimArr, fromOffset, fromStrVec, toIs);

  NimArr<mapDim, Tto> mapTo;
  mapTo.setMap(*toNimArr, toOffset, toStr, toIs);
  mapTo.mapCopy(mapFrom);
}

// the from must be 1D
// the from can be a map.
// the to cannot be a map
template <class Tfrom, class Tto>
void dynamicMapCopyFlatToDim(NimArrBase<Tto> *toNimArr, int toOffset,
                             vector<int> &toStr, vector<int> &toIs,
                             NimArrBase<Tfrom> *fromNimArr, int fromOffset,
                             int fromStr) {
  int mapDim = toStr.size();
  // must be the same as fromStr.sizes();
  if (toNimArr->isMap()) {
    PRINTF("Error, dynamicMapCopyFlatToDim is not set up for nested maps\n");
  }
  switch (mapDim) {
    case 1:
      dynamicMapCopyFlatToDimFixed<Tfrom, Tto, 1>(
          toNimArr, toOffset, toStr, toIs, fromNimArr, fromOffset, fromStr);
      break;
    case 2:
      dynamicMapCopyFlatToDimFixed<Tfrom, Tto, 2>(
          toNimArr, toOffset, toStr, toIs, fromNimArr, fromOffset, fromStr);
      break;
    case 3:
      dynamicMapCopyFlatToDimFixed<Tfrom, Tto, 3>(
          toNimArr, toOffset, toStr, toIs, fromNimArr, fromOffset, fromStr);
      break;
    case 4:
      dynamicMapCopyFlatToDimFixed<Tfrom, Tto, 4>(
          toNimArr, toOffset, toStr, toIs, fromNimArr, fromOffset, fromStr);
      break;
    default:
      PRINTF(
          "Error in copying (dynamicMapCopyFlatToDim): more than 4 dimensions "
          "not supported yet\n");
  }
}

template <class Tfrom, class Tto, int mapDim>
void dynamicMapCopyDimToFlatFixed(NimArrBase<Tto> *toNimArr, int toOffset,
                                  int toStr, NimArrBase<Tfrom> *fromNimArr,
                                  int fromOffset, vector<int> &fromStr,
                                  vector<int> &fromIs) {
  NimArr<mapDim, Tto> mapTo;
  vector<int> toStrVec(mapDim);
  toStrVec[0] = toStr;
  for (int i = 1; i < mapDim; i++) {
    toStrVec[i] = fromIs[i - 1] * toStrVec[i - 1];
  }
  mapTo.setMap(*toNimArr, toOffset, toStrVec, fromIs);

  NimArr<mapDim, Tfrom> mapFrom;
  mapFrom.setMap(*fromNimArr, fromOffset, fromStr, fromIs);
  mapTo.mapCopy(mapFrom);
}

// the to must be 1D
// the to can be a map.
// the from cannot be a map
template <class Tfrom, class Tto>
void dynamicMapCopyDimToFlat(NimArrBase<Tto> *toNimArr, int toOffset, int toStr,
                             NimArrBase<Tfrom> *fromNimArr, int fromOffset,
                             vector<int> &fromStr, vector<int> &fromIs) {
  int mapDim = fromStr.size();
  // must be the same as fromStr.sizes();
  if (fromNimArr->isMap()) {
    PRINTF("Error, dynamicMapCopyFlatToDim is not set up for nested maps\n");
  }
  switch (mapDim) {
    case 1:
      dynamicMapCopyDimToFlatFixed<Tfrom, Tto, 1>(
          toNimArr, toOffset, toStr, fromNimArr, fromOffset, fromStr, fromIs);
      break;
    case 2:
      dynamicMapCopyDimToFlatFixed<Tfrom, Tto, 2>(
          toNimArr, toOffset, toStr, fromNimArr, fromOffset, fromStr, fromIs);
      break;
    case 3:
      dynamicMapCopyDimToFlatFixed<Tfrom, Tto, 3>(
          toNimArr, toOffset, toStr, fromNimArr, fromOffset, fromStr, fromIs);
      break;
    case 4:
      dynamicMapCopyDimToFlatFixed<Tfrom, Tto, 4>(
          toNimArr, toOffset, toStr, fromNimArr, fromOffset, fromStr, fromIs);
      break;
    default:
      PRINTF(
          "Error in copying (dynamicMapCopyDimToFlat): more than 4 dimensions "
          "not supported yet\n");
  }
}

#endif
