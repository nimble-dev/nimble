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

#include "nimble/NamedObjects.h"
#include "nimble/Utils.h"
#include "nimble/Model.h"
#include "nimble/dllFinalizer.h"
#include "R.h"

//#define _DEBUG_NAMEDOBJECTS

void NamedObjects::NO_hw(  ) {
  PRINTF("hello world from NamedObjects\n");
}

SEXP copyFromRobject(SEXP Sextptr, SEXP Robject) {
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a a valid external pointer\n");
  }
  NamedObjects *m;
  m = static_cast< NamedObjects *>(R_ExternalPtrAddr(Sextptr));
  m->copyFromRobject(Robject);
  return(R_NilValue);
}

void* NamedObjects::getObjectPtr( string &name, bool warn_not_found ) {
  //cout<<name<<"\n";
  map<string, void *>::iterator iMO;
  iMO= namedObjects.find(name);		
  if(iMO == namedObjects.end()) {
    if(warn_not_found) {
      //std::cout<<"Error, could not find "<<name<<"\n";
      PRINTF("Error, could not find name\n");
      //    cout << "Name = " << name << "\n";
      _nimble_global_output << "Name = " << name << "\n"; nimble_print_to_R( _nimble_global_output);
      iMO = namedObjects.begin();
      _nimble_global_output << "Available Name 1 = " << iMO->first << "\n"; nimble_print_to_R( _nimble_global_output);
    }
    return(0);
  }
#ifdef _DEBUG_NAMEDOBJECTS
  PRINTF("Getting pointer to %s: %p\n", name.c_str(), iMO->second);
#endif
  return(	(iMO->second) ) ;
}

//GlobalObjects globalObjects;
SEXP getModelObjectPtr(SEXP Sextptr, SEXP Sname) {
  if(!Rf_isString(Sname)) {
    PRINTF("Error: Sname is not character!\n");
    return(R_NilValue);
  }
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a a valid external pointer\n");
    return(R_NilValue);
  }
  string name = STRSEXP_2_string(Sname, 0);
  NamedObjects *m = static_cast< NamedObjects *>(R_ExternalPtrAddr(Sextptr));
  void* cPtr = m->getObjectPtr(name);
  if(cPtr != 0){
	  SEXP output = R_MakeExternalPtr( cPtr, R_NilValue, R_NilValue);
	  PROTECT(output);
	  UNPROTECT(1);
	  return(output);
	}
	return(R_NilValue);
}

SEXP getAvailableNames(SEXP Sextptr) {
  ///_nimble_global_output << "In getAvailableNames\n"; nimble_print_to_R( _nimble_global_output);
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a a valid external pointer\n");
    return(R_NilValue);
  }
  NamedObjects *m;
  m = static_cast< NamedObjects *>(R_ExternalPtrAddr(Sextptr));
  
  SEXP Sans;
  int numNames = m->namedObjects.size();
  //  _nimble_global_output << "numNames = "<<numNames<<"\n"; nimble_print_to_R( _nimble_global_output);
  PROTECT(Sans = Rf_allocVector(STRSXP, numNames));
  //  m->hw();
  map<string, void *>::iterator iNO = m->getNamedObjects().begin();
  for(int i = 0; i < numNames; ++i, ++iNO) {
    // _nimble_global_output << "starting "<<i<<"\n"; nimble_print_to_R( _nimble_global_output);
    //_nimble_global_output << iNO->first.c_str() <<" \n";
    //nimble_print_to_R( _nimble_global_output);
    SET_STRING_ELT(Sans, i, PROTECT(Rf_mkChar(iNO->first.c_str())));
    //_nimble_global_output << "done with "<<i<<" "<<iNO->first<<" \n"; nimble_print_to_R( _nimble_global_output);
  }
  UNPROTECT(numNames + 1);
  return(Sans);
}

void* NumberedObjects::getObjectPtr(int index){
	return(numberedObjects[index]);
}

void NumberedObjects::setObjectPtr(int index, void* newPtr){
	numberedObjects[index] = newPtr;
}


void NumberedObjects::resize(int size){
	numberedObjects.resize(size);
}

SEXP getNumberedObject(SEXP Snp, SEXP index){
	NumberedObjects* np = static_cast<NumberedObjects*>(R_ExternalPtrAddr(Snp));
	void* vp = np->getObjectPtr(INTEGER(index)[0] - 1);
	SEXP ans = R_MakeExternalPtr(vp, R_NilValue, R_NilValue);
	PROTECT(ans);
	UNPROTECT(1);
	return(ans);
}

SEXP setNumberedObject(SEXP Snp, SEXP index, SEXP val){
	NumberedObjects* np = static_cast<NumberedObjects*>(R_ExternalPtrAddr(Snp));
	void* vp = static_cast<void*>(R_ExternalPtrAddr(val));
	np->setObjectPtr(INTEGER(index)[0] - 1, vp);
	return(R_NilValue);
}

SEXP resizeNumberedObjects(SEXP Snp, SEXP size){
	NumberedObjects* np = static_cast<NumberedObjects*>(R_ExternalPtrAddr(Snp));
	np->resize(INTEGER(size)[0]);
	return(R_NilValue);
}

SEXP getSizeNumberedObjects(SEXP Snp){
	NumberedObjects* np = static_cast<NumberedObjects*>(R_ExternalPtrAddr(Snp));
	SEXP ans = Rf_ScalarInteger(np->numberedObjects.size());
	PROTECT(ans);
	UNPROTECT(1);
	return(ans);
}

SEXP register_numberedObjects_Finalizer(SEXP Snp, SEXP Dll, SEXP Slabel) {
  //std::cout<< "In register_numberedObjects_Finalizer\n";
  //  R_RegisterCFinalizerEx(Snp, &numberedObjects_Finalizer, TRUE);
  RegisterNimbleFinalizer(Snp, Dll, &numberedObjects_Finalizer, Slabel);
  return(Snp);
}

void numberedObjects_Finalizer(SEXP Snp){
  //std::cout<< "In numberedObjects_Finalizer\n";
  NumberedObjects* np = static_cast<NumberedObjects*>(R_ExternalPtrAddr(Snp));
  if(np) delete np;
  R_ClearExternalPtr(Snp);
}

SEXP register_namedObjects_Finalizer(SEXP Snp, SEXP Dll, SEXP Slabel) {
  //  std::cout<< "In register_namedObjects_Finalizer\n";
  //R_RegisterCFinalizerEx(Snp, &namedObjects_Finalizer, TRUE);
  RegisterNimbleFinalizer(Snp, Dll, &namedObjects_Finalizer, Slabel);
  return(Snp);
}

void namedObjects_Finalizer(SEXP Snp){
  //  std::cout<< "In namedObjects_Finalizer\n";
  NamedObjects* np = static_cast<NamedObjects*>(R_ExternalPtrAddr(Snp));
  if(np) delete np;
  R_ClearExternalPtr(Snp);
}


SEXP newNumberedObjects(){
	NumberedObjects* np = new NumberedObjects;
	SEXP rPtr = R_MakeExternalPtr(np, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	// register from R to ensure finalizer will never be unloaded with this dll
	//R_RegisterCFinalizerEx(rPtr, &numberedObjects_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
}

