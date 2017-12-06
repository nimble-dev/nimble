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

#include "nimble/ModelClassUtils.h"

SEXP getModelValuesPtrFromModel (SEXP rPtr){
	void* vPtr = R_ExternalPtrAddr(rPtr);
	if(vPtr == 0){
		PRINTF("Warning: rPtr points to null!\n");
		return(R_NilValue);
	}
	ModelBase* modelBasePtr = static_cast<ModelBase*> (vPtr);
	void* cModelValuesPtr = static_cast<void *>((*modelBasePtr).getModelValuesPtr());
	SEXP rModelValuesPtr;
 	PROTECT(rModelValuesPtr = R_MakeExternalPtr(cModelValuesPtr, R_NilValue, R_NilValue));
	UNPROTECT(1);
	return(rModelValuesPtr) ;
}

NimArrType** cGetModelElementPtr(SEXP Sextptr, SEXP Sname) {
  if(!Rf_isString(Sname)) {
    PRINTF("Error: Sname is not character!\n");
    return(NULL);
  }
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a a valid external pointer\n");
    return(NULL);
  }
  string name = STRSEXP_2_string(Sname, 0);
  ModelBase *m = static_cast< ModelBase *>(R_ExternalPtrAddr(Sextptr));
  return( static_cast <NimArrType**> (m->getObjectPtr(name) ) );// I think this has wrong pointer-depth -Perry.  getObjectPtr always returns an address, so it is a **NimArr<>
}

SEXP getModelElementPtr(SEXP Sextptr, SEXP Sname){
	NimArrType** vPtr = cGetModelElementPtr(Sextptr, Sname) ;
	SEXP rPtr = R_MakeExternalPtr(vPtr, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	UNPROTECT(1);
	return(rPtr);	
}


// This is for the new method of building modelValues. 
// Returns a character string rName s.t. in R, .Call(rName)
// will build our modelValues object and return a pointer
SEXP getMVBuildName (SEXP rPtr){
	void* vPtr = R_ExternalPtrAddr(rPtr);
	if(vPtr == 0){
		PRINTF("Warning: rPtr points to null!\n");
		return(R_NilValue);
	}
	Values* modelValuesPtr = static_cast<Values*> (vPtr);
	string mVBuildName = (*modelValuesPtr).getMVBuildName();
	if (mVBuildName == "missing")
		PRINTF("Warning: buildName missing in modelValues! \nConstructor must assign the string buildName to the name which call a SEXP that builds modelValues object\n");
	SEXP rName = Rf_allocVector(STRSXP, 1);
	PROTECT(rName);
	SEXP mvbn = Rf_mkChar(mVBuildName.c_str() );
	PROTECT(mvbn);
	SET_STRING_ELT(rName, 0, mvbn );
	UNPROTECT(2);
	return(rName);
}

SEXP derefPtr(SEXP SmultiPtr) {
  void **doublePtr = static_cast<void **>(R_ExternalPtrAddr(SmultiPtr));
  return(R_MakeExternalPtr( static_cast<void *>(*doublePtr), R_NilValue, R_NilValue) );
}

