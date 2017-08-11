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

#ifndef __NIMBLE_SMART_POINTERS
#define __NIMBLE_SMART_POINTERS

#include "Utils.h"

//#define _DEBUG_SMARTPTRS

class nimSmartPtrBase {
 public:
  virtual void setPtrFromVoidPtr(void* & inputPtr)=0;
  virtual ~nimSmartPtrBase() {
#ifdef _DEBUG_SMARTPTRS
    PRINTF("smartPtrBase destructing\n");
#endif
  };
  virtual void* getVoidPtrToRealPtr()=0;
};

template<typename T>
class nimSmartPtr : public nimSmartPtrBase {
 public:
  T* realPtr;
  T* getRealPtr() {return(realPtr);}
  T& operator*() {return *realPtr;}
  T* operator->() {return realPtr;}
  void setPtr(const nimSmartPtr & input) {
    if(realPtr == input.realPtr) return;
    if(realPtr) realPtr->removeWatcher();
    realPtr = input.realPtr;
    realPtr->newWatcher();
  };
  void setPtrFromT(T* & inputPtr) {
    if(realPtr == inputPtr) return;
    if(realPtr) realPtr->removeWatcher();
    realPtr = inputPtr;
    realPtr->newWatcher();
  };
  void setPtrFromVoidPtr(void* & inputPtr) {
    T* tempPtr = static_cast<T*>(inputPtr);
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("setting pointer at %p to realPtr = %p (cast to T* from void* %p)\n", &realPtr, tempPtr, inputPtr);
#endif
    setPtrFromT( tempPtr ); 
  };
  void* getVoidPtrToRealPtr() {
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("getting void pointer %p to realPtr = %p\n", static_cast<void*>(&realPtr), realPtr);
#endif
    return(static_cast<void*>(&realPtr));
  }
  bool equalsPtr(const nimSmartPtr & otherPtr) {
	  return(realPtr == otherPtr.realPtr);
  }

  nimSmartPtr & operator=(const nimSmartPtr & rhs) {
    if(this == &rhs)
      return *this;
    setPtr(rhs);
    return *this;
  }
  
  nimSmartPtr & operator=(const T* & rhs) {
    if(realPtr == rhs)
      return *this;
    setPtrFromT(rhs);
    return *this;
  }

  nimSmartPtr() : realPtr(0) {
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("smartPtr constructing %p\n", &realPtr);
#endif
  };
  nimSmartPtr(const nimSmartPtr &rhs) {
    realPtr = rhs.realPtr;
    realPtr->newWatcher();
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("smartPtr constructing (from another smartPtr) with ptrToPtr = %p, and realPtr = %p\n", &realPtr, realPtr); 
#endif
  }

  nimSmartPtr(T* rhs) {
    realPtr = rhs;
    realPtr->newWatcher();
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("smartPtr constructing (from a T*) with ptrToPtr = %p, and realPtr = %p\n", &realPtr, realPtr); 
#endif
  }

  ~nimSmartPtr() {
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("smartPtr destructing with ptrToPtr = %p, and realPtr = %p\n", &realPtr, realPtr);
#endif
    if(realPtr != 0)
      realPtr->removeWatcher();
  };
};

// could be dangerous if used inside a vector<> due to the delete this

class pointedToBase {
 public:
  int watcherCount;
 pointedToBase() : watcherCount(0) {};
  void newWatcher() {
    watcherCount++;
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("Adding watcher to %p (now has %i watchers).\n", this, watcherCount);
#endif
  }
  void removeWatcher() {
    watcherCount--;
#ifdef _DEBUG_SMARTPTRS 
    PRINTF("Removing watcher to %p (now has %i watchers).\n", this, watcherCount);
#endif
    if(watcherCount <= 0) {
      if(watcherCount < 0) {
	PRINTF("Error, watcherCount went below 0.\n");
      }
#ifdef _DEBUG_SMARTPTRS 
      PRINTF("pointedToBase self-destructing\n");
#endif
      delete this;
    }
  }
  virtual ~pointedToBase() {};
};

extern "C" {
  SEXP register_pointedToBase_Finalizer(SEXP Snp, SEXP Dll, SEXP Slabel);
  SEXP register_smartPtrBase_Finalizer(SEXP Snp, SEXP Dll, SEXP Slabel);
}

// example
/* class pointedToDerived : public pointedToBase { */
/*  public: */
/*   pointedToDerived() {}; */
/*   void hw() {PRINTF("hello world\n");} */
/*   ~pointedToDerived() {PRINTF("In derived destructor\n");} */
/* }; */


#endif
