#ifndef __NIMBLEFUNCTIONUTILS
#define __NIMBLEFUNCTIONUTILS

#include "NamedObjects.h"
#include "RcppUtils.h"
#include "NimArr.h"


extern "C" {
	SEXP setNFPointer (SEXP RnfPtr, SEXP RmodelElementPtr);		// This function takes an a pointer from an NimbleFunction
																// object and points it at the pointer of a modelValues object
}






#endif