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
}

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
	SEXP ans = ScalarInteger(np->numberedObjects.size());
	PROTECT(ans);
	UNPROTECT(1);
	return(ans);
}


void numberedObjects_Finalizer(SEXP Snp){
  NumberedObjects* np = static_cast<NumberedObjects*>(R_ExternalPtrAddr(Snp));
  if(np) delete np;
  R_ClearExternalPtr(Snp);
}

void namedObjects_Finalizer(SEXP Snp){
  NamedObjects* np = static_cast<NamedObjects*>(R_ExternalPtrAddr(Snp));
  if(np) delete np;
  R_ClearExternalPtr(Snp);
}


SEXP newNumberedObjects(){
	NumberedObjects* np = new NumberedObjects;
	SEXP rPtr = R_MakeExternalPtr(np, R_NilValue, R_NilValue);
	PROTECT(rPtr);
	R_RegisterCFinalizerEx(rPtr, &numberedObjects_Finalizer, TRUE);
	UNPROTECT(1);
	return(rPtr);
}

/*		This is an example of a finalizer for the templated numbered objects
		They must be written by hand
		for each template
void SingleModelValuesAccessor_NumberedObjects_Finalizer(SEXP Snp){
	SpecialNumberedObjects<SingleModelValuesAccess>* np 
	= static_cast<SpecialNumberedObjects<SingleModelValuesAccess>*>(R_ExternalPtrAddr(Snp));
	delete np;
}
*/
