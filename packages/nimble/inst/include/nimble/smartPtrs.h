#ifndef __NIMBLE_SMART_POINTERS
#define __NIMBLE_SMART_POINTERS

#include "Utils.h"

class nimSmartPtrBase {
 public:
  virtual void setPtrFromVoidPtr(void* & inputPtr)=0;
  virtual ~nimSmartPtrBase() {PRINTF("smartPtrBase destructing\n"); };
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
    PRINTF("setting pointer at %p to realPtr = %p (cast to T* from void* %p)\n", &realPtr, tempPtr, inputPtr);
    setPtrFromT( tempPtr ); 
  };
  void* getVoidPtrToRealPtr() {
    PRINTF("getting void pointer %p to realPtr = %p\n", static_cast<void*>(&realPtr), realPtr);
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

  nimSmartPtr() : realPtr(0) {PRINTF("smartPtr constructing %p\n", &realPtr); };
  nimSmartPtr(const nimSmartPtr &rhs) {
    realPtr = rhs.realPtr;
    realPtr->newWatcher();
    PRINTF("smartPtr constructing (from another smartPtr) with ptrToPtr = %p, and realPtr = %p\n", &realPtr, realPtr); 
  }

  nimSmartPtr(T* rhs) {
    realPtr = rhs;
    realPtr->newWatcher();
    PRINTF("smartPtr constructing (from a T*) with ptrToPtr = %p, and realPtr = %p\n", &realPtr, realPtr); 
  }

  ~nimSmartPtr() {
    PRINTF("smartPtr destructing with ptrToPtr = %p, and realPtr = %p\n", &realPtr, realPtr);
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
    PRINTF("Adding watcher to %p (now has %i watchers).\n", this, watcherCount);
  }
  void removeWatcher() {
    watcherCount--;
    PRINTF("Removing watcher to %p (now has %i watchers).\n", this, watcherCount);
    if(watcherCount <= 0) {
      if(watcherCount < 0) {
	PRINTF("Error, watcherCount went below 0.\n");
      }
      PRINTF("pointedToBase self-destructing\n");
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
