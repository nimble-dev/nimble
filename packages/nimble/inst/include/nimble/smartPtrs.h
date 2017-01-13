#ifndef __NIMBLE_SMART_POINTERS
#define __NIMBLE_SMART_POINTERS

#include "Utils.h"

class nimSmartPtrBase {
	public:
	virtual void setPtrFromVoidPtr(void* & inputPtr)=0;
};

template<typename T>
class nimSmartPtr : public nimSmartPtrBase {
 public:
  T* realPtr;
  T& operator*() {return *realPtr;}
  T* operator->() {return realPtr;}
  void setPtr(const nimSmartPtr & input) {
    if(realPtr) realPtr->removeWatcher();
    realPtr = input.realPtr;
    realPtr->newWatcher();
  };
  void setPtrFromT(T* & inputPtr) {
    if(realPtr) realPtr->removeWatcher();
    realPtr = inputPtr;
	PRINTF("STEP 4a \n");
    realPtr->newWatcher();
  };
  void setPtrFromVoidPtr(void* & inputPtr) {
	    PRINTF("STEP 1");
		T* tempPtr = static_cast<T*>(inputPtr);
	    setPtrFromT( tempPtr ); 
  PRINTF("STEP 2");
  };
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

  nimSmartPtr() : realPtr(0) {};
  nimSmartPtr(const nimSmartPtr &rhs) {
    realPtr = rhs.realPtr;
    realPtr->newWatcher();
  }

  nimSmartPtr(T* rhs) {
    realPtr = rhs;
    realPtr->newWatcher();
  }

  ~nimSmartPtr() {
    realPtr->removeWatcher();
  };
};

// could be dangerous if used inside a vector<> due to the delete this

class pointedToBase {
 public:
  int watcherCount;
 pointedToBase() : watcherCount(0) {};
  void newWatcher() {
	  printf("%d \n", watcherCount);
	  PRINTF("Watcher ++ \n");
watcherCount++;
	  PRINTF("watcherCount is %d \n", watcherCount);
}
  void removeWatcher() {
	  	PRINTF("Watcher -- \n");

    watcherCount--;
	PRINTF("watcherCount is %d \n", watcherCount);

    if(watcherCount <= 0) {
      if(watcherCount < 0) {
	PRINTF("Error, a watcherCount went below 0. \n");
      }
      delete this;
    }
  }
  virtual ~pointedToBase() {};
};

// example
/* class pointedToDerived : public pointedToBase { */
/*  public: */
/*   pointedToDerived() {}; */
/*   void hw() {PRINTF("hello world\n");} */
/*   ~pointedToDerived() {PRINTF("In derived destructor\n");} */
/* }; */


#endif
