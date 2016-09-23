// Code to be included in each on-the-fly (dynamic) nimble compilation
// It will be included by a generated file in each project
#ifndef __DYNAMIC_REGISTRATIONS
#define __DYNAMIC_REGISTRATIONS

#include <nimble/RcppUtils.h>
#include <nimble/RcppNimbleUtils.h>
#include <nimble/ModelClassUtils.h>
#include <nimble/accessorClasses.h>
#include <nimble/dists.h>
#include <nimble/NamedObjects.h>
#include <nimble/dllFinalizer.h>

#include <R_ext/Rdynload.h>

#define FUN(name, numArgs) \
  {#name, (DL_FUNC) &name, numArgs}

#define CFUN(name, numArgs) \
  {"R_"#name, (DL_FUNC) &name, numArgs}

R_CallMethodDef CallEntries[] = {
  // from ModelClassUtils.h
  {"getModelValuesPtrFromModel", (DL_FUNC) &getModelValuesPtrFromModel, 1},
  FUN(getModelElementPtr, 2),
  FUN(getMVBuildName, 1),
  FUN(derefPtr, 1),
  //  FUN(setSinglePtrFromSinglePtr, 2),

  // accessorClasses
  //FUN(makeSingleVariableAccessor, 4),
  //  FUN(makeSingleModelValuesAccessor, 5),
  FUN(getModelAccessorValues, 1),
  FUN(getMVAccessorValues, 1),
  //  FUN(setNodeModelPtr, 3),
  //FUN(newManyVariableAccessor, 1),
  FUN(addSingleVariableAccessor, 4),
  FUN(resizeManyModelVarAccessor, 2),
  FUN(removeModelVariableAccessor, 3),
  //FUN(newManyModelValuesAccessor, 1),
  FUN(resizeManyModelValuesAccessor, 2),
  FUN(addSingleModelValuesAccessor, 4),
  FUN(removeModelValuesAccessor, 3),
  FUN(manualSetNRows, 2),
  //FUN(getVarAndIndicesExtPtr, 2),
  FUN(getVarAndIndices, 1),
  FUN(varAndIndices2mapParts, 3),
  FUN(var2mapParts, 3),
  FUN(populateNodeFxnVector_byGID, 3),
  FUN(populateNodeFxnVectorNew_byDeclID, 4),
  FUN(populateIndexedNodeInfoTable, 2),
  FUN(populateValueMapAccessorsFromNodeNames, 4),
  FUN(populateValueMapAccessors, 3),
  FUN(populateNumberedObject_withSingleModelValuesAccessors, 5),
  FUN(populateCopierVector, 5),
  FUN(populateNumberedObject_withSingleModelVariablesAccessors, 5),
  FUN(populateModelVariablesAccessors_byGID, 5),
  //  FUN(new_SingleModelValuesAccessor_NumberedObjects, 0),
  //FUN(new_SingleModelVariablesAccessor_NumberedObjects, 0),

  //dists
  // these don't need to be linked into an on-the-fly dll because they will be in nimble.so

  //RcppUtils
  // ditto, these will be in nimble.so

  //RcppNimbleUtils
  FUN(setDoublePtrFromSinglePtr, 2),
  FUN(addBlankModelValueRows, 2),
  FUN(getNRow, 1),
  FUN(copyModelValuesElements, 4),
  FUN(getMVElement, 2),
  FUN(getMVsize, 1),
  FUN(getMVElementAsList, 2),
  FUN(setMVElementFromList, 3),
  FUN(matrix2VecNimArr, 4),
  //  FUN(printMVElement, 2),
  FUN(setMVElement, 3),
  FUN(resizeNumListRow, 3),
  FUN(setNumListRows, 3),
  //FUN(setVarPointer, 3),
  FUN(makeNumericList, 3),
  FUN(Nim_2_SEXP, 2),
  FUN(SEXP_2_Nim, 4),
  FUN(setPtrVectorOfPtrs, 3),
  FUN(setOnePtrVectorOfPtrs, 3),
  //FUN(getOnePtrVectorOfPtrs, 2),
  FUN(getEnvVar_Sindex, 3),
  FUN(getEnvVar, 2),
  FUN(setEnvVar_Sindex, 4),
  FUN(setEnvVar, 3),

  // from NamedObjects
  FUN(getModelObjectPtr, 2),
  FUN(getAvailableNames, 1),
  FUN(getNumberedObject, 2),
  FUN(setNumberedObject, 3),
  FUN(resizeNumberedObjects, 2),
  FUN(getSizeNumberedObjects, 1),
  FUN(newNumberedObjects, 0),
  FUN(register_namedObjects_Finalizer, 2),
  FUN(register_numberedObjects_Finalizer, 2),
  FUN(register_VecNimArr_Finalizer, 2),

  FUN(RNimble_Ptr_ManualFinalizer, 1),
  FUN(RNimble_Ptr_CheckAndRunAllDllFinalizers, 1),
  
  {NULL, NULL, 0}
};


/* R_CMethodDef CEntries[] = { */
/*   FUN(RegisterNimblePointer, 3), */
/*   {NULL, NULL, 0} */
/* } */

// Something like this will be generated with each .so/.dll nimble creates
// however it is required to be named R_init_SONAME, so it must be generated for each one.
//
// extern "C"
// void
// R_init_nimble_on_the_fly(DllInfo *dll)
// {
//     R_registerRoutines(dll, CEntries OR NULL, CallEntries, NULL, NULL);
// }

#endif
