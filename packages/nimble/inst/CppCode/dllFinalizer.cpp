#include <unordered_map>
#include <utility>
#include "nimble/dll.h"

typedef std::pair<SEXP, R_CFinalizer_t> DllAndFinalizer;

void RegisterNimblePointer(SEXP ptr, SEXP Dll, R_CFinalizer_t finalizer) {
  RnimblePtrs[ptr] = DllAndFinalizer(Dll, finalizer);
  R_RegisterCFinalizerEx(ptr, RNimble_PtrFinalizer, TRUE);
}

void RNimble_PtrFinalizer(SEXP obj) {
  R_CFinalizer_t cfun;
  DllAndFinalizer::iterator value = RnimblePtrs.find(obj);
  if(value == RnimblePtrs.end()) {
    std::cout<<"Trying to finalize a pointer whose object was already cleared\n";
    return;
  }
}
