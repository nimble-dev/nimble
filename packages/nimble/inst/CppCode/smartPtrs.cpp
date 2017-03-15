#include <iostream>
#include <nimble/smartPtrs.h>
#include <nimble/dllFinalizer.h>

void pointedToBase_Finalizer(SEXP Snp){
  // std::cout<< "In namedObjects_Finalizer\n";
  pointedToBase* np = static_cast<pointedToBase*>(R_ExternalPtrAddr(Snp));
  if(np) {
    np->removeWatcher(); /* object will naturally self-destruct if watcher count goes to 0*/
  }
  R_ClearExternalPtr(Snp);
}

SEXP register_pointedToBase_Finalizer(SEXP Snp, SEXP Dll, SEXP Slabel) {
  std::cout<< "In register_pointedToBase_Finalizer\n";
  RegisterNimbleFinalizer(Snp, Dll, &pointedToBase_Finalizer, Slabel);
  return(Snp);
}

