#ifndef __DLLFINALIZER
#define __DLLFINALIZER

// modeled after DTL ideas in dll_auto_unload branch

#include <Rinternals.h>

void RegisterNimblePointer(SEXP ptr, SEXP Dll, R_CFinalizer_t finalizer);

#define RegisterNimbleFinalizer RegisterNimblePointer

#endif
