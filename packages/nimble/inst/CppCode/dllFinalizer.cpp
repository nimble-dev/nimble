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

#include <map>
#include <utility>
#include <string>
#include <vector>
//#include "nimble/RcppUtils.h"
#include "nimble/dllFinalizer.h"

//#define _DEBUG_DLLFINALIZER
//typedef std::pair<SEXP, R_CFinalizer_t> DllAndFinalizer;

struct DllAndFinalizer {
public:
  SEXP Dll;
  R_CFinalizer_t Finalizer;
  string Label;
};

std::map<SEXP, DllAndFinalizer> RnimblePtrs;
typedef std::map<SEXP, DllAndFinalizer>::iterator RnimblePtrsIterator;

void finalizeOneObject(RnimblePtrsIterator RNPiter) {
    R_CFinalizer_t cfun;
    cfun = RNPiter->second.Finalizer;
    if(cfun) {
      std::cout<<"finalizing "<<RNPiter->first<<std::endl;
      cfun(RNPiter->first);
      std::cout<<"done finalizing "<<RNPiter->first<<std::endl;
    }
    R_ClearExternalPtr(RNPiter->first);
    RnimblePtrs.erase(RNPiter);
}

void RNimble_PtrFinalizer(SEXP obj) {
  // add check that obj is an external pointer
  RnimblePtrsIterator value = RnimblePtrs.find(obj);
  if(value == RnimblePtrs.end()) {
    //PRINTF("Trying to finalize a pointer whose object was already cleared\n");
    return;
  }
  finalizeOneObject(value);
}

string local_STRSEXP_2_string(SEXP Ss, int i) {
  if(!Rf_isString(Ss)) {
    PRINTF("Error: STRSEXP_2_string called for SEXP that is not a string!\n"); 
    return(string(""));
  }
  if(LENGTH(Ss) <= i) {
    PRINTF("Error: STRSEXP_2_string called for (C) element %i of an SEXP that has length %i!\n", i, LENGTH(Ss));
    return(string(""));
  }
  int l = LENGTH(STRING_ELT(Ss, i));
  string ans(CHAR(STRING_ELT(Ss,i)),l);
  return(ans);
}


void RegisterNimblePointer(SEXP ptr, SEXP Dll, R_CFinalizer_t finalizer, SEXP Slabel) { // same as RegisterNimbleFinalizer
  RnimblePtrsIterator value = RnimblePtrs.find(ptr);
  // add check that ptr is an external pointer
  string label;
  if(Slabel != R_NilValue) {
    label = local_STRSEXP_2_string(Slabel, 0);
  } else {
    label = string("");
  }

  if(value != RnimblePtrs.end()) {
    PRINTF("Error: trying to register a finalizer for an external pointer (for %s) that has already has one\n", label.c_str());
    return;
  }
  DllAndFinalizer newDLLAndFinalizer;
  newDLLAndFinalizer.Dll = Dll;
  newDLLAndFinalizer.Finalizer = finalizer;
  newDLLAndFinalizer.Label = label;
#ifdef _DEBUG_DLLFINALIZER
   PRINTF("Adding label %s\n",  newDLLAndFinalizer.Label.c_str());
#endif
   RnimblePtrs[ptr] = newDLLAndFinalizer;
   R_RegisterCFinalizerEx(ptr, RNimble_PtrFinalizer, TRUE);
}

SEXP CountDllObjects() {
  SEXP Sans;
  PROTECT(Sans = Rf_allocVector(INTSXP, 1));
  INTEGER(Sans)[0] = RnimblePtrs.size();
  UNPROTECT(1);
  return(Sans);
}

SEXP RNimble_Ptr_ManualFinalizer(SEXP obj) {
  RNimble_PtrFinalizer(obj);
  return(R_NilValue);
}


SEXP local_vectorString_2_STRSEXP(const std::vector<string> &v) {
  SEXP Sans;
  int nn = v.size();
  PROTECT(Sans = Rf_allocVector(STRSXP, nn));
  for(int i = 0; i < nn; i++) {
    SET_STRING_ELT(Sans, i, Rf_mkChar(v[i].c_str()));
  }
  UNPROTECT(1);
  return(Sans);
}

SEXP RNimble_Ptr_CheckAndRunAllDllFinalizers(SEXP Dll, SEXP Sforce) {
  RnimblePtrsIterator RNPiter;
  int objectsFound(0);
  bool force = LOGICAL(Sforce)[0];
  std::vector<string> objectsFoundLabels;
  for(RNPiter = RnimblePtrs.begin();
      RNPiter != RnimblePtrs.end();
      ) {
    if(RNPiter->second.Dll == Dll) {
      objectsFound++;
#ifdef _DEBUG_DLLFINALIZER
      PRINTF("Found label %s\n",  RNPiter->second.Label.c_str());
#endif	
      objectsFoundLabels.push_back(RNPiter->second.Label);
      if(force) {
	finalizeOneObject(RNPiter++);// pass by copy a current iterator value and increment, since finalizerOneObject will use .erase()
      } else
	++RNPiter; 
    } else {
      ++RNPiter;
    }
  }
  if(objectsFound > 0) {
    if(force) 
      PRINTF("Warning: %i objects were force-cleared from a DLL\n", objectsFound);
    else
      PRINTF("Warning: %i objects were found from a DLL\n", objectsFound);
  }
  //  SEXP Sans;
  //  PROTECT(Sans = Rf_allocVector(INTSXP, 1));
  //  INTEGER(Sans)[0] = objectsFound;
  //  UNPROTECT(1);
  //  return(Sans);

  return(local_vectorString_2_STRSEXP(objectsFoundLabels));
}
