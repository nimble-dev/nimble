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

#include "nimble/NimArrBase.h"
#include "nimble/NamedObjects.h"
#include "nimble/RcppNimbleUtils.h"
#include "nimble/dllFinalizer.h"
//#include "nimble/RcppUtils.h"
#include "nimble/Utils.h"
#include "nimble/smartPtrs.h"
#include<iostream>
#include<sstream>
#include<algorithm>

// Work around R's unqualified use of `error` in the `ERROR` macro
// (see $R_SOURCE/src/include/R_ext/RS.h).
#define error Rf_error

void SEXP_2_NimArr(SEXP Sn, double &x) {x = SEXP_2_double(Sn);}
void SEXP_2_NimArr(SEXP Sn, int &x) {x = SEXP_2_int(Sn);}
void SEXP_2_NimArr(SEXP Sn, bool &x) {x = SEXP_2_bool(Sn);}
void SEXP_2_NimArr(SEXP Sn, std::string &x) {x = STRSEXP_2_string(Sn);}
void SEXP_2_NimArr(SEXP Sn, std::vector<std::string> &x) {STRSEXP_2_vectorString(Sn, x);}

SEXP NimArr_2_SEXP(double x) {return(double_2_SEXP(x));}
SEXP NimArr_2_SEXP(int x) {return(int_2_SEXP(x));}
SEXP NimArr_2_SEXP(bool x) {return(bool_2_SEXP(x));}
SEXP NimArr_2_SEXP(std::string &x) {return(string_2_STRSEXP(x));}
SEXP NimArr_2_SEXP(const std::vector<std::string> &x) {return(vectorString_2_STRSEXP(x));}

// This is used when we have a NimArr<>* in a model and a NimArr<>** that needs to point to it.
// We assume we have an extptr to each
SEXP setDoublePtrFromSinglePtr(SEXP SdoublePtr, SEXP SsinglePtr) {
  void *singlePtr = R_ExternalPtrAddr(SsinglePtr); // this is really a **
  void **doublePtr = static_cast<void **>(R_ExternalPtrAddr(SdoublePtr)); // this is really a ***.  
  *doublePtr = singlePtr;
  return(R_NilValue);
}

void setNimbleFxnPtr_copyFromRobject(void *nf_to, SEXP S_NF_from) {
  void **doublePtr = static_cast<void **>(nf_to);
  SEXP Scnf, SsinglePtr;
  SEXP S_pxData;
   PROTECT(S_pxData = Rf_allocVector(STRSXP, 1));
   SET_STRING_ELT(S_pxData, 0, Rf_mkChar(".xData"));
  // environment(modelVar)$.CobjectInterface
   PROTECT(Scnf = Rf_findVarInFrame(PROTECT(GET_SLOT(
						       S_NF_from,
						       S_pxData)),
				    Rf_install(".CobjectInterface"))
	   );
   int unprotectCount = 3;
   if(Rf_isNewList(Scnf)) {
    // multi-interface
    //Cnf[[1]]$basePtrList[[ Cnf[[2]] ]]
    SEXP Sindex;
    PROTECT(Sindex = VECTOR_ELT(Scnf, 1));
    int index = (Rf_isInteger(Sindex) ? INTEGER(Sindex)[0] : REAL(Sindex)[0]); 
    index--; // From 1-based to 0-based indexing 
    PROTECT(SsinglePtr = VECTOR_ELT(
				    Rf_findVarInFrame(PROTECT(GET_SLOT(
								       VECTOR_ELT(Scnf,
										  0),
								       S_pxData)),
						      Rf_install("basePtrList")),
				    index
				    )
	    );
    unprotectCount += 3;
  } else {
    // printf("in non-list\n");
    // full interface
    // Cnf$.basePtr
    PROTECT(SsinglePtr =  Rf_findVarInFrame(PROTECT(GET_SLOT(
							     Scnf,
							     S_pxData)),
					    Rf_install(".basePtr")));
    unprotectCount += 2;
  }
  void *singlePtr = R_ExternalPtrAddr(SsinglePtr);
  *doublePtr = singlePtr;
  UNPROTECT(unprotectCount);
}

// probably deprecated
SEXP setSmartPtrFromSinglePtr(SEXP StoPtr, SEXP SfromPtr) {
  void *fromPtr = R_ExternalPtrAddr(SfromPtr); 
  nimSmartPtrBase *toSmartPtr = static_cast<nimSmartPtrBase*>(R_ExternalPtrAddr(StoPtr));  
  toSmartPtr->setPtrFromVoidPtr(fromPtr);
  return(R_NilValue);
}

SEXP setSmartPtrFromDoublePtr(SEXP StoPtr, SEXP SfromPtr) {
  void *fromPtr = *static_cast<void**>(R_ExternalPtrAddr(SfromPtr)); 
  nimSmartPtrBase *toSmartPtr = static_cast<nimSmartPtrBase*>(R_ExternalPtrAddr(StoPtr));  
  toSmartPtr->setPtrFromVoidPtr(fromPtr);
  return(R_NilValue);
}

SEXP setPtrVectorOfPtrs(SEXP SaccessorPtr, SEXP ScontentsPtr, SEXP Ssize) {
  vectorOfPtrsAccessBase *accessorPtr = static_cast<vectorOfPtrsAccessBase *>(R_ExternalPtrAddr(SaccessorPtr)); 
  void *contentsPtr = static_cast<void *>(R_ExternalPtrAddr(ScontentsPtr));
  int size = INTEGER(Ssize)[0];
  accessorPtr->setTheVec(contentsPtr, size);
  return(R_NilValue);
}

SEXP setOnePtrVectorOfPtrs(SEXP SaccessorPtr, SEXP Si, SEXP ScontentsPtr) {
  vectorOfPtrsAccessBase *accessorPtr = static_cast<vectorOfPtrsAccessBase *>(R_ExternalPtrAddr(SaccessorPtr)); 
  void *contentsPtr = static_cast<void *>(R_ExternalPtrAddr(ScontentsPtr));
  int i = INTEGER(Si)[0];
  accessorPtr->setVecPtr(i, contentsPtr);
  return(R_NilValue);
}

/*This function was created so we can report the status of different operations
from C++ to R. I found that Rbreak does not stop the R function from which the C++
function was called. This can cause an excessive number (e.g. nrow) of printed 
error reports if the user makes a common mistake. Instead, this function is used to 
return the status to R. If it is false, we use the R function "stop" which will break
for loops, etc. */
 
SEXP returnStatus(bool OK)
	{
	SEXP ans;
    PROTECT(ans = Rf_allocVector(LGLSXP, 1) );
    LOGICAL(ans)[0] = OK;
    UNPROTECT(1);
    return(ans);
	}

bool all(bool x) {return(x);}
bool any(bool x) {return(x);}
int t(int x) {return(x);}
double t(double x) {return(x);}
int prod(int x) {return(x);}
double prod(double x) {return(x);}

NimArr<1, double> vectorDouble_2_NimArr(vector<double> input) {
  NimArr<1, double> output;
  output.setSize(input.size(), false, false);
  std::copy(input.begin(), input.end(), output.getPtr());
  return(output);
}

SEXP setVecNimArrRows(SEXP Sextptr, SEXP nRows, SEXP setSize2row1){
    NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
    nimType vecType = (*typePtr).getNimType();
    if(vecType == nimType::DOUBLE){
    	VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
	    int new_size = INTEGER(nRows)[0];
    	matPtr->resize(new_size);
    	if(new_size <= nrowCpp)
    		return(returnStatus(true));
    	if(LOGICAL(setSize2row1)[0]	== TRUE){
		    	NimArrBase<double> *thisRow;
    			thisRow = matPtr->getBasePtr(0);
    			int numDims = thisRow->numDims();
		     	vector<int> Dims(numDims);
 				for(int i = 0; i < numDims; i++)
    	  			Dims[i] = thisRow->dimSize(i);
   		   		for(int i = nrowCpp; i < new_size; i++){
   		   			thisRow = matPtr->getBasePtr(i);
    	  			thisRow->setSize(Dims);
    	  		}
    	return(returnStatus(true) );
    	}
    }
    else if(vecType == nimType::INT){
    	VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
    	  int new_size = INTEGER(nRows)[0];
	  	matPtr->resize(new_size);
	  	if(new_size <= nrowCpp)
    		return(returnStatus(true));

    	if(LOGICAL(setSize2row1)[0]	== TRUE){
    		NimArrBase<int> *thisRow;
    		thisRow = matPtr->getBasePtr(0);
    		int numDims = thisRow->numDims();
  	    	vector<int> Dims(numDims);
			for(int i = 0; i < numDims; i++)
      			Dims[i] = thisRow->dimSize(i);
      		for(int i = nrowCpp; i < new_size; i++){
      			thisRow = matPtr->getBasePtr(i);
      			thisRow->setSize(Dims);
      		}
      	}
    return(returnStatus(true) );
    }
 	
 	PRINTF("Data type for VecNimArr not currently supported\n");
 	return(returnStatus(false) ) ;
}

/* Add blank rows to C model values */
SEXP addBlankModelValueRows(SEXP Sextptr, SEXP numAdded){
    if(!Rf_isInteger(numAdded)) {
        PRINTF("Error: numAdded is not an integer!\n");
        return(returnStatus(false) );
    }
    
    if(!R_ExternalPtrAddr(Sextptr)) {
        PRINTF("Error: Sextptr is not a valid external pointer\n");
        return(returnStatus(false) );
    }
    NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
    nimType vecType = (*typePtr).getNimType();
    if(vecType == nimType::DOUBLE){
    	VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
	    int new_size = INTEGER(numAdded)[0] + nrowCpp;
    	matPtr->resize(new_size);
    	NimArrBase<double> *thisRow;
    	thisRow = matPtr->getBasePtr(0);
    	int numDims = thisRow->numDims();
	     vector<int> Dims(numDims);
 		for(int i = 0; i < numDims; i++)
      		Dims[i] = thisRow->dimSize(i);
    	for(int i = nrowCpp; i < INTEGER(numAdded)[0] + nrowCpp; i++){
      		thisRow = matPtr->getBasePtr(i); 
      		thisRow->setSize(Dims);
    	}
    return(returnStatus(true) );
    }
    
    else if(vecType == nimType::INT){
    	VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(R_ExternalPtrAddr(Sextptr));
    	int nrowCpp = matPtr->size();
    	  int new_size = INTEGER(numAdded)[0] + nrowCpp;
	  	matPtr->resize(new_size);
    	NimArrBase<int> *thisRow;
    	thisRow = matPtr->getBasePtr(0);
    	int numDims = thisRow->numDims();
  	    vector<int> Dims(numDims);
		for(int i = 0; i < numDims; i++)
      		Dims[i] = thisRow->dimSize(i);
    	for(int i = nrowCpp; i < INTEGER(numAdded)[0] + nrowCpp; i++){
      		thisRow = matPtr->getBasePtr(i); 
      		thisRow->setSize(Dims);
    	}
    return(returnStatus(true) );
    }
 	
 	PRINTF("Data type for VecNimArr not currently supported\n");
 	return(returnStatus(false) ) ;
}

/* Cliff's new function for retreiving nrow from CmodelValues object */
SEXP getNRow(SEXP Sextptr){
    if(!R_ExternalPtrAddr(Sextptr)) {
        PRINTF("Error: Sextptr is not a valid external pointer\n");
        return(R_NilValue);
    }
    int nRow = 0;
    NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
    nimType vecType = (*typePtr).getNimType();
    if(vecType == nimType::DOUBLE){
    	VecNimArrBase<double>* matPtr = static_cast<VecNimArrBase<double>* > (typePtr);
	    nRow = matPtr->size();
	}
	else if(vecType == nimType::INT)	{
    	VecNimArrBase<int>* matPtr = static_cast<VecNimArrBase<int>* > (typePtr);
	    nRow = matPtr->size();
	}
	else{
		PRINTF("Data type of VecNimArr not currently supported\n");
		//_nimble_global_output<< "vecType = " << vecType << "\n"; nimble_print_to_R(_nimble_global_output);
//		cout << "vecType DOUBLE = " << DOUBLE << "\n";
		}
	SEXP rNRow;
	PROTECT(rNRow = Rf_allocVector(INTSXP, 1) ); 
	INTEGER(rNRow)[0] = nRow;
	UNPROTECT(1);
    return(rNRow);
}


/* Cliff's new function for copying variables of given rows from one CModelValues object to another */
/* This copies the values from SextptrFrom into SextptrTo. It takes values from rowsFrom and copies
into rowsTo */
SEXP copyModelValuesElements(SEXP SextptrFrom, SEXP SextptrTo, SEXP rowsFrom, SEXP rowsTo){
    if(!R_ExternalPtrAddr(SextptrFrom) | !R_ExternalPtrAddr(SextptrFrom)) {
        PRINTF("Error: Sextptr is not a valid external pointer\n");
        return(returnStatus(false));
    }
    if(LENGTH(rowsFrom) != LENGTH(rowsTo) ){
        PRINTF("Error: length of rowsFrom != rowsTo");
        return(returnStatus(false) );
    }
    if(LENGTH(rowsFrom) == 0){
        return(returnStatus(true) );
    }
  
  NimVecType *typePtrFrom = static_cast< NimVecType* >(R_ExternalPtrAddr(SextptrFrom));
  nimType vecTypeFrom = (*typePtrFrom).getNimType();
  NimVecType *typePtrTo = static_cast< NimVecType* >(R_ExternalPtrAddr(SextptrTo));
  nimType vecTypeTo = (*typePtrTo).getNimType();

  if((vecTypeFrom == nimType::DOUBLE) & (vecTypeTo == nimType::DOUBLE)){
	  VecNimArrBase<double> *matPtrFrom = static_cast< VecNimArrBase<double>* >(typePtrFrom);
	  VecNimArrBase<double> *matPtrTo = static_cast< VecNimArrBase<double>* >(typePtrTo );
	  int k = LENGTH(rowsFrom);
	  int sizeFrom = matPtrFrom->size();
	  int sizeTo = matPtrTo->size();
	  int *indexFrom = INTEGER(rowsFrom);
	  int *indexTo = INTEGER(rowsTo);
	  NimArrBase<double> *thisRowFrom;
	  NimArrBase<double> *thisRowTo;
	  int ncFrom = 0;
	  int ncTo = 0;
    
	  for(int i = 0; i < k; i++){
	  if((indexFrom[i] > sizeFrom) | (indexFrom[i] <= 0))
	    {
	    _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexFrom[i] << " sizeFrom = " << sizeFrom << "\n";
	    nimble_print_to_R(_nimble_global_output);
	    if(i > 0)
	      PRINTF("Warning: partial copy completed before error discovered!\n");
	    return(returnStatus(false) );
    	}
  	if((indexTo[i] > sizeTo) | (indexTo[i] <=0))
  	  {
	    _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexTo[i] << " sizeTo = "<< sizeTo << "\n";
	    nimble_print_to_R(_nimble_global_output);
  	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	thisRowFrom = matPtrFrom->getBasePtr(indexFrom[i] - 1);  
  	thisRowTo = matPtrTo->getBasePtr(indexTo[i] - 1);  
  
  	ncFrom = thisRowFrom->size();
  	ncTo = thisRowTo->size();
  
  	if(ncFrom != ncTo){
  	  PRINTF("Error: ncFrom != ncTo\n");
  	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	for(int j = 0; j < ncFrom; j++)
  		    (*thisRowTo)[j] = (*thisRowFrom)[j];
  	}
  	return(returnStatus(true) );
  }

  if((vecTypeFrom == nimType::INT) & (vecTypeTo == nimType::INT)){
	  VecNimArrBase<int> *matPtrFrom = static_cast< VecNimArrBase<int>* >(typePtrFrom);
	  VecNimArrBase<int> *matPtrTo = static_cast< VecNimArrBase<int>* >(typePtrTo );
	  int k = LENGTH(rowsFrom);
	  int sizeFrom = matPtrFrom->size();
	  int sizeTo = matPtrTo->size();
	  int *indexFrom = INTEGER(rowsFrom);
	  int *indexTo = INTEGER(rowsTo);
	  NimArrBase<int> *thisRowFrom;
	  NimArrBase<int> *thisRowTo;
	  int ncFrom = 0;
	  int ncTo = 0;
    
	  for(int i = 0; i < k; i++){
	  if((indexFrom[i] > sizeFrom) | (indexFrom[i] <= 0))
	    {
	    _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexFrom[i] << " sizeFrom = " << sizeFrom << "\n";
	    nimble_print_to_R(_nimble_global_output);
	    if(i > 0)
	      PRINTF("Warning: partial copy completed before error discovered!\n");
	    return(returnStatus(false) );
    	}
  	if((indexTo[i] > sizeTo) | (indexTo[i] <=0))
  	  {
  	  _nimble_global_output<<"Warning: invalid index to copy from. Index = " << indexTo[i] << " sizeTo = "<< sizeTo << "\n";
	  nimble_print_to_R(_nimble_global_output);
	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	thisRowFrom = matPtrFrom->getBasePtr(indexFrom[i] - 1);  
  	thisRowTo = matPtrTo->getBasePtr(indexTo[i] - 1);  
  
  	ncFrom = thisRowFrom->size();
  	ncTo = thisRowTo->size();
  
  	if(ncFrom != ncTo){
  	  PRINTF("Error: ncFrom != ncTo\n");
  	  if(i > 0)
  	    PRINTF("Warning: partial copy completed before error discovered!\n");
  	  return(returnStatus(false) );
  	  }
  	for(int j = 0; j < ncFrom; j++)
  		    (*thisRowTo)[j] = (*thisRowFrom)[j];
  	}
  	return(returnStatus(true) );
  }

  PRINTF("Data type not currently supported for VecNimArr. Copy MV elements failed\n");
  return(returnStatus(false) );

}	

SEXP getMVElement(SEXP Sextptr, SEXP Sindex){
	if(!Rf_isInteger(Sindex)) {
    	PRINTF("Error: Sindex is not an integer!\n");
    	return(returnStatus(false) ) ;
	}
  	if(!R_ExternalPtrAddr(Sextptr)) {
  	  PRINTF("Error: Sextptr is not a valid external pointer\n");
  	  return(returnStatus(false) ) ;
  	}
	NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
	nimType vecType = (*typePtr).getNimType();
	int nrowCpp = typePtr->size();
	int index = INTEGER(Sindex)[0];
	if(index > nrowCpp){
		PRINTF("Error: index too large\n");
		return(returnStatus(false) ) ;
	  }
  	if(index < 1){
		PRINTF("Error: index < 1\n");
		return(returnStatus(false) ) ;
  	}
  return(cGetMVElementOneRow(typePtr, vecType, index) ) ;	
}  	

// This is not called directly from R.  It is called from getMVElement
SEXP cGetMVElementOneRow(NimVecType* typePtr, nimType vecType, int index) 	{  	
    if(vecType == nimType::DOUBLE){
		VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(typePtr);
		NimArrBase<double> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting R index to C index
		int outputLength = thisRow->size();
  		SEXP Sans;
  		PROTECT(Sans = Rf_allocVector(REALSXP, outputLength));
		double *value = REAL(Sans);
  		std::copy(thisRow->getPtr(), thisRow->getPtr() + outputLength, value); /* copy should work now */
  		if(thisRow->numDims() > 1) { // using first row as representative
    		SEXP Sdim;
    		PROTECT(Sdim = Rf_allocVector(INTSXP, thisRow->numDims()));
    		for(int idim = 0; idim < thisRow->numDims(); ++idim) INTEGER(Sdim)[idim] = thisRow->dimSize(idim);
    			Rf_setAttrib(Sans, R_DimSymbol, Sdim);
    		UNPROTECT(2);
  		}
  		else {
    UNPROTECT(1);
  		}
  	 	return(Sans);
  	}
    else if(vecType == nimType::INT){
		VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(typePtr);
//		int nrowCpp = matPtr->size();
		NimArrBase<int> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting R index to C index
		int outputLength = thisRow->size();
  		SEXP Sans;
  		PROTECT(Sans = Rf_allocVector(INTSXP, outputLength));
		int *value = INTEGER(Sans);
  		std::copy(thisRow->getPtr(), thisRow->getPtr() + outputLength, value); /* copy should work now */
  		if(thisRow->numDims() > 1) { // using first row as representative
    		SEXP Sdim;
    		PROTECT(Sdim = Rf_allocVector(INTSXP, thisRow->numDims()));
    		for(int idim = 0; idim < thisRow->numDims(); ++idim) INTEGER(Sdim)[idim] = thisRow->dimSize(idim);
    			Rf_setAttrib(Sans, R_DimSymbol, Sdim);
    		UNPROTECT(2);
  		}
  		else {
    UNPROTECT(1);
  		}
  	  	return(Sans);	
  	}
  	else
  		PRINTF("VecNimArr datatype not supported\n");
	return(R_NilValue);
}

SEXP setMVElementFromList(SEXP vNimPtr, SEXP rList, SEXP Sindices){
  // This function needs to get the length of rList
  // Once it has this length, it marches down the length of rList and sets all the values of
  // of vNimPtr
  int listLength = LENGTH(rList);
  int checkLength = LENGTH(Sindices);
  if(listLength != checkLength){
    PRINTF("Number of indices copying to does not match length of list\n");
    return(R_NilValue);
  }
  SEXP listElement;
  NimVecType *typePtr = static_cast<NimVecType*>(R_ExternalPtrAddr(vNimPtr) );
  nimType vecType = (*typePtr).getNimType();
  int index;
  for(int i = 0; i < listLength ; i++){
    listElement = VECTOR_ELT(rList, i);
    index = INTEGER(Sindices)[i];
    cSetMVElementSingle(typePtr, vecType, index, listElement);
  }
  return(R_NilValue);
}


 void cSetMVElementSingle(NimVecType* typePtr, nimType vecType,  int index, SEXP Svalue) {
    if(vecType == nimType::DOUBLE){
		VecNimArrBase<double> *matPtr = static_cast< VecNimArrBase<double>* >(typePtr);
  		NimArrBase<double> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting from R index to C index
  		int nc = thisRow->size();
  		double *value = REAL(Svalue);
  		int new_nc = LENGTH(Svalue);
  		if(new_nc != nc){
  			PRINTF("Error: size of values assigned not the same! If this is during compiling, could be from improperly sized r modelValues variable (see user manual section on modelValues)\n");
		return;
  		}
  for(int i = 0; i < nc; i++)
    (*thisRow)[i] = value[i];
	return;
	}

	else if(vecType == nimType::INT){
		VecNimArrBase<int> *matPtr = static_cast< VecNimArrBase<int>* >(typePtr);
  		NimArrBase<int> *thisRow;
  		thisRow = matPtr->getBasePtr(index - 1);		//Converting from R index to C index
  		int nc = thisRow->size();
  		double *value = REAL(Svalue);
  		int new_nc = LENGTH(Svalue);
  		if(new_nc != nc){
  			PRINTF("Error: size of values assigned not the same!\n");
		return;
  		}
  	for(int i = 0; i < nc; i++)
    	(*thisRow)[i] = value[i];
	return;
	}
	PRINTF("VecNimArr data type not currently supported\n");
	return;
}



 SEXP getMVElementAsList(SEXP SextPtr, SEXP Sindices){
 	int nIndice = LENGTH(Sindices);
 	SEXP sList;
 	PROTECT(sList = Rf_allocVector(VECSXP, nIndice) );
 	NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(SextPtr));
	nimType vecType = (*typePtr).getNimType();
	for(int i = 0; i < nIndice; i++)
		SET_VECTOR_ELT(sList, i, cGetMVElementOneRow(typePtr, vecType, INTEGER(Sindices)[i]) );
	UNPROTECT(1);
	return(sList);
 }

SEXP setMVElement(SEXP Sextptr, SEXP Sindex, SEXP Svalue){
    if(!Rf_isInteger(Sindex)) {
    PRINTF("Error: Sindex is not an integer!\n");
	return(returnStatus(false) ) ;
  }
  if(!R_ExternalPtrAddr(Sextptr)) {
    PRINTF("Error: Sextptr is not a valid external pointer\n");
	return(returnStatus(false) ) ;
  }
  int index = INTEGER(Sindex)[0];
  if(index < 1){ 
    PRINTF("Error: index < 1\n");   	
	return(returnStatus(false) ) ;
  }

   NimVecType *typePtr = static_cast< NimVecType* >(R_ExternalPtrAddr(Sextptr));
   nimType vecType = (*typePtr).getNimType();
  int nrowCpp = typePtr->size();
  if(index > nrowCpp){
	PRINTF("Error: index too large\n");
	return(returnStatus(false) ) ;
  }
  cSetMVElementSingle( typePtr, vecType, index, Svalue );
  //SEXP SEXP_2_int(SEXP rPtr, SEXP refNum, SEXP rScalar);
  return(returnStatus(true) );
 } 
 
 
 SEXP getMVsize(SEXP Sextptr){
	NimVecType	*nimVecPtr = static_cast<NimVecType*>(R_ExternalPtrAddr(Sextptr) ); 
 	vector<int> cDims = (*nimVecPtr).getRowDims(0);
	int nDims = cDims.size();
 	SEXP Sans;
 	PROTECT(Sans = Rf_allocVector(INTSXP, nDims) );
 	for(int i = 0; i < nDims; i++)
 		INTEGER(Sans)[i] = cDims[i];
 	UNPROTECT(1);
 	return(Sans);
 	}

SEXP NimArrDouble_2_SEXP(NimArrBase<double> &nimArrDbl){
		int len = nimArrDbl.size();
		SEXP Sans;  
		PROTECT(Sans = Rf_allocVector(REALSXP, len));
		//		std::copy(nimArrDbl.v.begin(), nimArrDbl.v.end() , REAL(Sans) );
		std::copy(nimArrDbl.v, nimArrDbl.v + len , REAL(Sans) );  
  		int numDims = nimArrDbl.numDims();
  		if(numDims > 1) {
		  SEXP Sdim;
		  PROTECT(Sdim = Rf_allocVector(INTSXP, numDims));
		  for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = nimArrDbl.dimSize(idim);
		  Rf_setAttrib(Sans, R_DimSymbol, Sdim);
		  UNPROTECT(2);
  		}
  		else {
    		UNPROTECT(1);
  		} 
	return(Sans);
 }

SEXP NimArrInt_2_SEXP(NimArrBase<int> &nimArrInt){
		int len = nimArrInt.size();
		SEXP Sans;  
		PROTECT(Sans = Rf_allocVector(INTSXP, len));
		//		std::copy(nimArrInt.v.begin(), nimArrInt.v.end() , INTEGER(Sans) );
		std::copy(nimArrInt.v, nimArrInt.v + len , INTEGER(Sans) );  
  		int numDims = nimArrInt.numDims();
  		if(numDims > 1) {
    		SEXP Sdim;
    		PROTECT(Sdim = Rf_allocVector(INTSXP, numDims));
    		for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = nimArrInt.dimSize(idim);
			Rf_setAttrib(Sans, R_DimSymbol, Sdim);
		    UNPROTECT(2);
  		}
  		else {
    		UNPROTECT(1);
  		} 
	return(Sans);
 }

SEXP NimArrBool_2_SEXP(NimArrBase<bool> &nimArrBl){
  int len = nimArrBl.size();
  SEXP Sans;  
  PROTECT(Sans = Rf_allocVector(LGLSXP, len));
  //		std::copy(nimArrInt.v.begin(), nimArrInt.v.end() , INTEGER(Sans) );
  std::copy(nimArrBl.v, nimArrBl.v + len , LOGICAL(Sans) );  
  int numDims = nimArrBl.numDims();
  if(numDims > 1) {
    SEXP Sdim;
    PROTECT(Sdim = Rf_allocVector(INTSXP, numDims));
    for(int idim = 0; idim < numDims; ++idim) INTEGER(Sdim)[idim] = nimArrBl.dimSize(idim);
    Rf_setAttrib(Sans, R_DimSymbol, Sdim);
    UNPROTECT(2);
  }
  else {
    UNPROTECT(1);
  } 
  return(Sans);
 }


NimArrType* getNimTypePtr(SEXP &rPtr, SEXP &refNum)
{
	NimArrType* nimTypePtr;
	int cRefNum = INTEGER(refNum)[0];
	if(cRefNum == 1){
		nimTypePtr = static_cast<NimArrType*>(R_ExternalPtrAddr(rPtr) );
		return(nimTypePtr); 
		}
	if(cRefNum == 2){
		nimTypePtr = (*static_cast<NimArrType**> (R_ExternalPtrAddr(rPtr) ) );
		return(nimTypePtr);
		}
	PRINTF("Warning: number of dereferences not currently supported, returning null pointer\n");
	nimTypePtr = NULL;
	return(nimTypePtr);
}

void SEXP_2_NimArrDouble (SEXP rValues, NimArrBase<double> &NimArrDbl){
  int rLength = LENGTH(rValues);
  if(rLength != NimArrDbl.size() ) {
    PRINTF("Warning: R object of different size than NimArrDouble. R obj has size %i but NimArrDbl has size %i.\n", rLength, NimArrDbl.size());
    return;
  }
  if(Rf_isReal(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrDbl[i] = REAL(rValues)[i];
  }
  else if(Rf_isInteger(rValues) | Rf_isLogical(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrDbl[i] = INTEGER(rValues)[i];
  }
  
  else
    PRINTF("WARNING: class of R object not recognized. Should be numeric or integer\n");
  return;
}

void SEXP_2_NimArrInt (SEXP rValues, NimArrBase<int> &NimArrInt){
  int rLength = LENGTH(rValues);
  if(rLength != NimArrInt.size() ) {
    PRINTF("Warning: R object of different size than NimArrInt!\n");
    return;		
  }
  
  if(Rf_isInteger(rValues) | Rf_isLogical(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrInt[i] = INTEGER(rValues)[i];
  }
  else if(Rf_isReal(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrInt[i] = REAL(rValues)[i];
  }
  
  else
    PRINTF("WARNING: class of R object not recognized. Should be numeric or integer\n");    
  return;
}

void SEXP_2_NimArrBool (SEXP rValues, NimArrBase<bool> &NimArrBl){
  int rLength = LENGTH(rValues);
  if(rLength != NimArrBl.size() ) {
    PRINTF("Warning: R object of different size than NimArrBl!\n");
    return;		
  }

  // In R, Logical is represented as integer
  if(Rf_isInteger(rValues) | Rf_isLogical(rValues)) {
    for(int i = 0; i < rLength; i++)
      NimArrBl[i] = INTEGER(rValues)[i];
  }
  else if(Rf_isReal(rValues) ) {
    for(int i = 0; i < rLength; i++)
      NimArrBl[i] = REAL(rValues)[i];
  }
  
  else
    PRINTF("WARNING: class of R object not recognized. Should be numeric or integer\n");    
  return;
}


SEXP Nim_2_SEXP(SEXP rPtr, SEXP NumRefers){
	NimArrType* nimTypePtr = getNimTypePtr(rPtr, NumRefers);
	if(!nimTypePtr)
		return(R_NilValue);
	if(	(*nimTypePtr).getNimType() == nimType::INT){
		NimArrBase<int>* nimBase = static_cast<NimArrBase<int> *>(nimTypePtr);
		return(NimArrInt_2_SEXP( (*nimBase) ) ) ;
	}
	if(	(*nimTypePtr).getNimType() == nimType::DOUBLE){
		NimArrBase<double>* nimBase = static_cast<NimArrBase<double> *>(nimTypePtr);
		return(NimArrDouble_2_SEXP( (*nimBase) ) );
	}
	if(	(*nimTypePtr).getNimType() == nimType::BOOL){
	  NimArrBase<bool>* nimBase = static_cast<NimArrBase<bool> *>(nimTypePtr);
	  return(NimArrBool_2_SEXP( (*nimBase) ) );
	}

	PRINTF("Datatype of NimArr not found\n");
	return(R_NilValue);
}

void SEXP_2_Nim_internal(NimArrType* nimTypePtr,
			 SEXP rValues,
			 bool resize = true) {
    vector<int> sexpDims = getSEXPdims( rValues );
    int sexpNumDims = sexpDims.size();
    
    if(!nimTypePtr)
      return;
    
    if((*nimTypePtr).getNimType() == nimType::INT){
      NimArrBase<int>* nimBase = static_cast<NimArrBase<int> *>(nimTypePtr);
      int nimNumDims = (*nimBase).numDims();
      if(nimNumDims != sexpNumDims) {
	if((LENGTH(rValues) != (*nimBase).size()) & 
	   (sexpNumDims != 1)){
	  PRINTF("Incorrect number of dimensions in copying\n");
	  return;
	}
      }
      if(resize)
	(*nimBase).setSize(sexpDims);
      SEXP_2_NimArrInt(rValues,  (*nimBase) ) ;
    }
    if((*nimTypePtr).getNimType() == nimType::DOUBLE){
      
      NimArrBase<double>* nimBase = static_cast<NimArrBase<double> *>(nimTypePtr);
      int nimNumDims = (*nimBase).numDims();
      
      if(nimNumDims != sexpNumDims){
	if((LENGTH(rValues) != (*nimBase).size()) &
	   (sexpNumDims != 1)){
	  PRINTF("Incorrect number of dimensions in copying\n");
	  return;
	}
      }
      if((resize) & (nimNumDims == sexpNumDims))
	(*nimBase).setSize(sexpDims);
      SEXP_2_NimArrDouble( rValues, (*nimBase) ) ;
    }
    if((*nimTypePtr).getNimType() == nimType::BOOL){
      NimArrBase<bool>* nimBase = static_cast<NimArrBase<bool> *>(nimTypePtr);
      int nimNumDims = (*nimBase).numDims();
      if(nimNumDims != sexpNumDims){
	if((LENGTH(rValues) != (*nimBase).size()) &
	   (sexpNumDims != 1)){
	  PRINTF("Incorrect number of dimensions in copying\n");
	  return;
	}
      }
      if(resize)
	(*nimBase).setSize(sexpDims);
      SEXP_2_NimArrBool(rValues,  (*nimBase) ) ;
    }
}

void SEXP_2_Nim_for_copyFromRobject(void *NimArrPtr, SEXP rValues) {
  NimArrType* nimTypePtr = static_cast<NimArrType*>(NimArrPtr);
  SEXP_2_Nim_internal(nimTypePtr, rValues);
}

SEXP SEXP_2_Nim(SEXP rPtr, SEXP NumRefers, SEXP rValues, SEXP allowResize){
    bool resize = (LOGICAL(allowResize)[0] == TRUE);
    NimArrType* nimTypePtr = getNimTypePtr(rPtr, NumRefers);
    SEXP_2_Nim_internal(nimTypePtr, rValues, resize);
    return(R_NilValue);
}

void VecNimArr_Finalizer(SEXP Sp) {
  //  std::cout<< "In VecNimArr_Finalizer\n";
  NimVecType* np = static_cast<NimVecType*>(R_ExternalPtrAddr(Sp));
  if(np) delete np;
  R_ClearExternalPtr(Sp);

}

SEXP register_VecNimArr_Finalizer(SEXP Sp, SEXP Dll) {
  //  std::cout<< "In register_VecNimArr_Finalizer\n";
  //  R_RegisterCFinalizerEx(Sp, &VecNimArr_Finalizer, TRUE);
  RegisterNimbleFinalizer(Sp, Dll, &VecNimArr_Finalizer, R_NilValue);
  return(Sp);
}

double nimMod(double a, double b) {
  return(fmod(a, b));
}

// bool compareOrderedPair(orderedPair a, orderedPair b) {  //function called for sort 
//   return(a.value < b.value);
// }

void rankSample(NimArr<1, double> &weights, int n, NimArr<1, int> &output) {
  bool silent = false;
  rankSample(weights, n, output, silent);
}

void rankSample(NimArr<1, double> &weights, int n, NimArr<1, int> &output, bool silent) {
  //PRINTF("in VOID rankSample\n");
  output.setSize(n);
  int N = weights.size();
  //GetRNGstate();
  rawSample(weights.getPtr(), n, N, output.getPtr(), false, silent);
  //PutRNGstate();
}




void row2NimArr(SEXP matrix, NimArrBase<double>* nimPtr, int startPoint, int len, int nRows){
	for(int i = 0; i < len; i++)
		(*nimPtr)[i] = REAL(matrix)[startPoint + i * nRows];
}

void row2NimArr(SEXP matrix, NimArrBase<int>* nimPtr, int startPoint, int len, int nRows){
  for(int i = 0; i < len; i++)
		(*nimPtr)[i] = INTEGER(matrix)[startPoint + i * nRows];	
}

SEXP matrix2VecNimArr(SEXP RvecNimPtr, SEXP matrix, SEXP rowStart, SEXP rowEnd){
	int cRowStart = INTEGER(rowStart)[0] - 1;
	int cRowEnd = INTEGER(rowEnd)[0] - 1;
	NimVecType* vecTypePtr = static_cast<NimVecType*>(R_ExternalPtrAddr(RvecNimPtr) );

	// Look at dimensions of input matrix...
	SEXP RmatrixDim;
	PROTECT(RmatrixDim = Rf_getAttrib(matrix, R_DimSymbol));
	// ... and of the NimArr
	vector<int> rowDims = vecTypePtr->getRowDims(0);
	int len = 1;
	for(unsigned int i = 0; i < rowDims.size(); i++)
	  len = len * rowDims[i];
	// len is the product of the dimensions of the NimArr
	int nRows;
	// A variety of cases were written below, but actually only
	// the matrix case is supported.  This provides an inverse to as.matrix(mv).
	// copying of individual variables can be done as cmv['v'] <- v
	if(RmatrixDim == R_NilValue) {
	  NIMERROR("In matrix2VecNimArr: matrix argument must be a matrix, but a vector was provided\n");
	  // // If matrix is a vector, then use the NimArr to determine the dimensions.
	  // double dnRows = ((double) Rf_length(matrix)) / ((double) len);
	  // if(dnRows != floor(dnRows)) {
	  //   NIMERROR("In matrix2VecNimArr: Length of matrix is not congruent with dimensions of modelValues variable.\n");
	  // }
	  // nRows = Rf_length(matrix) / len;	
	} else if(Rf_length(RmatrixDim) == 1) {
	  NIMERROR("In matrix2VecNimArr: matrix argument must be a matrix, but a one-dimensional array was provided\n");
	  // // If matrix is not a vector, then if it is a 1D array, handle it like a vector
	  // double dnRows = ((double) Rf_length(matrix)) / ((double) len);
	  // if(dnRows != floor(dnRows)) {
	  //   NIMERROR("In matrix2VecNimArr: Length of matrix is not congruent with dimensions of modelValues variable.\n");
	  // }
	  // nRows = LENGTH(matrix) / len;
	} else {
	  // If matrix has dimension >= 2, then check that the dimensions are congruent with those of the NimArr.
	  // Assume that the first index goes over "rows", so the remaining indices should match
	  // either individually, or should match the product of the NimArr dims.
	  if(Rf_length(RmatrixDim) == 2) {
	    // It is two-dimensional, then the second dimension should match the flattened size
	    // (i.e. product of dimensions) of the NimArr
	    if( INTEGER(RmatrixDim)[1] != len ) {
	      NIMERROR("In matrix2VecNimArr: Length of matrix is not congruent with dimensions of modelValues.  Second dimension should match the size of the modelValues variable.\n");
	    }
	    nRows = INTEGER(RmatrixDim)[0];
	  } else {
	    // if(Rf_length(RmatrixDim) != rowDims.size()-1) {
	    //   NIMERROR("In matrix2VecNimArr: Length of matrix is not congruent with dimensions of modelValues.  If matrix is > 2 dimensions, the number of dimenions of the matrix should be one greater than the number of dimensions of the modelValues variable.\n");
	    // }
	    // for(int iii = 1; iii < Rf_length(RmatrixDim); iii++) {
	    //   if(INTEGER(RmatrixDim)[iii] != rowDims[iii-1]) {
	    // 	NIMERROR("In matrix2VecNimArr: Length of matrix is not congruent with dimensions of modelValues.  If matrix is > 2 dimensions, dimensions past the first should match those of the modelValues variable.\n");
	    //   }
	    // }
	    NIMERROR("In matrix2VecNimArr: matrix argument must be a matrix, but an array with >2 dimensions^ was provided\n");
	  }
	}

	nimType varType = vecTypePtr->getNimType();
	if(varType == nimType::DOUBLE){
		VecNimArrBase<double>* vecPtr = static_cast<VecNimArrBase<double>*>(vecTypePtr);
		for(int i = cRowStart; i <= cRowEnd; i++)
			row2NimArr(matrix, static_cast<NimArrBase<double>*>(vecPtr->getRowTypePtr(i) ), cRowStart + i , len, nRows);
	}
	if(varType == nimType::INT){
		VecNimArrBase<int>* vecPtr = static_cast<VecNimArrBase<int>*>(vecTypePtr);
		for(int i = cRowStart; i <= cRowEnd; i++)
			row2NimArr(matrix, static_cast<NimArrBase<int>*>(vecPtr->getRowTypePtr(i) ), cRowStart + i , len, nRows);
	}
	UNPROTECT(1);
	return(R_NilValue);
}


SEXP getEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sIndex){
	SEXP ans;
	int cIndex = INTEGER(sIndex)[0] - 1;
	ans = Rf_findVar(Rf_install(CHAR(STRING_ELT(sString, cIndex ))), sEnv);
	PROTECT(ans);
	UNPROTECT(1);
	return(ans);
	}

SEXP getEnvVar(SEXP sString, SEXP sEnv){
  	return(getEnvVar_Sindex(sString, sEnv, Rf_ScalarInteger(1)));
}

SEXP setEnvVar_Sindex(SEXP sString, SEXP sEnv, SEXP sVal, SEXP sIndex){
	int cIndex = INTEGER(sIndex)[0] - 1;
	Rf_setVar(Rf_install(CHAR(STRING_ELT(sString, cIndex ))), sVal, sEnv);
	return(R_NilValue);
}

 SEXP setEnvVar(SEXP sString, SEXP sEnv, SEXP sVal){						
  	return(setEnvVar_Sindex(sString, sEnv, sVal, Rf_ScalarInteger(1)));
  }
  
