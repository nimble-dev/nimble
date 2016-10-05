#ifndef __DLLFINALIZER
#define __DLLFINALIZER

// modeled after DTL ideas in dll_auto_unload branch

#include "Utils.h"
#include <Rinternals.h>

void RegisterNimblePointer(SEXP ptr, SEXP Dll, R_CFinalizer_t finalizer, SEXP Slabel);
extern "C" {
  SEXP RNimble_Ptr_ManualFinalizer(SEXP obj);
  SEXP RNimble_Ptr_CheckAndRunAllDllFinalizers(SEXP Dll, SEXP Sforce);
  SEXP CountDllObjects();
}

#define RegisterNimbleFinalizer RegisterNimblePointer

#endif
