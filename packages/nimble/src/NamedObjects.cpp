#include "nimble/NamedObjects.h"
#include "nimble/Utils.h"
#include "nimble/Model.h"
#include "R.h"


void* NamedObjects::getObjectPtr( string &name ) {
  //  cout<<name<<"\n";
  map<string, void *>::iterator iMO;
  iMO= namedObjects.find(name);		
  if(iMO == namedObjects.end()) {
    //std::cout<<"Error, could not find "<<name<<"\n";
    PRINTF("Error, could not find name\n");
    cout << "Name = " << name << "\n";
    iMO = namedObjects.begin();
    cout << "Available Name 1 = " << iMO->first << "\n";
    return(0);
  }
  return(	(iMO->second) ) ;
}

//GlobalObjects globalObjects;
SEXP getModelObjectPtr(SEXP Sextptr, SEXP Sname) {
  if(!isString(Sname)) {
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
};

SEXP getAvailableNames(SEXP Sextptr) {
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a a valid external pointer\n");
    return(R_NilValue);
  }
  NamedObjects *m = static_cast< NamedObjects *>(R_ExternalPtrAddr(Sextptr));

  SEXP Sans;
  int numNames = m->namedObjects.size();
  PROTECT(Sans = allocVector(STRSXP, numNames));
  map<string, void *>::iterator iNO = m->namedObjects.begin();
  for(int i = 0; i < numNames; ++i, ++iNO) {
    SET_STRING_ELT(Sans, i, mkChar(iNO->first.c_str()));
  }
  UNPROTECT(1);
  return(Sans);
};






