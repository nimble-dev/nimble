#include <map>
#include <utility>
#include <string>
#include <vector>
//#include "nimble/RcppUtils.h"
#include "nimble/dllFinalizer.h"


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
    if(cfun)
      cfun(RNPiter->first);
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
  if(!isString(Ss)) {
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
  PRINTF("Adding label %s\n",  newDLLAndFinalizer.Label.c_str());
  RnimblePtrs[ptr] = newDLLAndFinalizer;
  R_RegisterCFinalizerEx(ptr, RNimble_PtrFinalizer, TRUE);
}

SEXP CountDllObjects() {
  SEXP Sans;
  PROTECT(Sans = allocVector(INTSXP, 1));
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
  PROTECT(Sans = allocVector(STRSXP, nn));
  for(int i = 0; i < nn; i++) {
    SET_STRING_ELT(Sans, i, mkChar(v[i].c_str()));
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
      PRINTF("Found label %s\n",  RNPiter->second.Label.c_str());
	
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
  //  PROTECT(Sans = allocVector(INTSXP, 1));
  //  INTEGER(Sans)[0] = objectsFound;
  //  UNPROTECT(1);
  //  return(Sans);

  return(local_vectorString_2_STRSEXP(objectsFoundLabels));
}
