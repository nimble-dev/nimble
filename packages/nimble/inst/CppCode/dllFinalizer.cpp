#include <unordered_map>
#include <utility>
#include "nimble/dllFinalizer.h"

typedef std::pair<SEXP, R_CFinalizer_t> DllAndFinalizer;

std::unordered_map<SEXP, DllAndFinalizer> RnimblePtrs;
typedef std::unordered_map<SEXP, DllAndFinalizer>::iterator RnimblePtrsIterator;

void finalizeOneObject(RnimblePtrsIterator RNPiter) {
    R_CFinalizer_t cfun;
    cfun = RNPiter->second.second;
    if(cfun)
      cfun(RNPiter->first);
    R_ClearExternalPtr(RNPiter->first);
    RnimblePtrs.erase(RNPiter);
}

void RNimble_PtrFinalizer(SEXP obj) {
  // add check that obj is an external pointer
  RnimblePtrsIterator value = RnimblePtrs.find(obj);
  if(value == RnimblePtrs.end()) {
    PRINTF("Trying to finalize a pointer whose object was already cleared\n");
    return;
  }
  finalizeOneObject(value);
}

void RegisterNimblePointer(SEXP ptr, SEXP Dll, R_CFinalizer_t finalizer) { // same as RegisterNimbleFinalizer
  RnimblePtrsIterator value = RnimblePtrs.find(ptr);
  // add check that ptr is an external pointer
  if(value != RnimblePtrs.end()) {
    PRINTF("Error: trying to register a finalizer for an external pointer that has already has one\n");
    return;
  }
  RnimblePtrs[ptr] = DllAndFinalizer(Dll, finalizer);
  R_RegisterCFinalizerEx(ptr, RNimble_PtrFinalizer, TRUE);
}

SEXP RNimble_Ptr_ManualFinalizer(SEXP obj) {
  RNimble_PtrFinalizer(obj);
  return(R_NilValue);
}

SEXP RNimble_Ptr_CheckAndRunAllDllFinalizers(SEXP Dll) {
  RnimblePtrsIterator RNPiter;
  int objectsFound(0);
  for(RNPiter = RnimblePtrs.begin();
      RNPiter != RnimblePtrs.end();
      ) {
    if(RNPiter->second.first == Dll) {
      objectsFound++;
      finalizeOneObject(RNPiter++); // pass by copy a current iterator value and increment, since finalizerOneObject will use .erase()
    } else {
      ++RNPiter;
    }
  }
  if(objectsFound > 0) {
    PRINTF("Warning: %i objects were cleared from a DLL\n", objectsFound);
  }
  return(R_NilValue);
}
