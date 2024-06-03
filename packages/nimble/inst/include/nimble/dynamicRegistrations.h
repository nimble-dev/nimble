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
#include <nimble/smartPtrs.h>
#include <nimble/predefinedNimbleLists.h>
#include <nimble/nimOptim.h>

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

  // nimOptim
  FUN(CALL_NimOptimProblem_fn, 2),
  FUN(CALL_NimOptimProblem_gr, 2),
  FUN(CALL_NimOptimProblem_he, 2),

  // accessorClasses
  FUN(manualSetNRows, 2),
  FUN(getVarAndIndices, 1),
  FUN(varAndIndices2mapParts, 3),
  FUN(var2mapParts, 3),
  FUN(populateNodeFxnVectorNew_byDeclID, 4),
  FUN(populateNodeFxnVectorNew_byDeclID_forDerivs, 5),
  FUN(populateIndexedNodeInfoTable, 2),
  FUN(populateValueMapAccessorsFromNodeNames, 4),
  FUN(populateValueMapAccessors, 3),
  FUN(populateCopierVector, 5),

  //dists
  // these don't need to be linked into an on-the-fly dll because they will be in nimble.so

  //RcppUtils
  // ditto, these will be in nimble.so

  //RcppNimbleUtils
  FUN(setDoublePtrFromSinglePtr, 2),
  FUN(setSmartPtrFromSinglePtr, 2),
  FUN(setSmartPtrFromDoublePtr, 2),
  FUN(addBlankModelValueRows, 2),
  FUN(getNRow, 1),
  FUN(copyModelValuesElements, 4),
  FUN(getMVElement, 2),
  FUN(getMVsize, 1),
  FUN(getMVElementAsList, 2),
  FUN(setMVElementFromList, 3),
  FUN(matrix2VecNimArr, 4),
  FUN(setMVElement, 3),
  FUN(Nim_2_SEXP, 2),
  FUN(SEXP_2_Nim, 4),
  FUN(extract_double_2_SEXP, 2),
  FUN(populate_SEXP_2_double, 3),
  FUN(extract_int_2_SEXP, 2),
  FUN(populate_SEXP_2_int, 3),
  FUN(extract_bool_2_SEXP, 2),
  FUN(populate_SEXP_2_bool, 3),
  
  FUN(populate_SEXP_2_string, 2),
  FUN(extract_string_2_SEXP, 1),
  FUN(populate_SEXP_2_stringVector, 2),
  FUN(extract_stringVector_2_SEXP, 1),

  FUN(setVecNimArrRows, 3),
  
  FUN(setPtrVectorOfPtrs, 3),
  FUN(setOnePtrVectorOfPtrs, 3),
  FUN(getEnvVar_Sindex, 3),
  FUN(getEnvVar, 2),
  FUN(setEnvVar_Sindex, 4),
  FUN(setEnvVar, 3),

  // from NamedObjects
  FUN(getModelObjectPtr, 2),
  FUN(getAvailableNames, 1),
  FUN(copyFromRobject, 2),
  FUN(getNumberedObject, 2),
  FUN(setNumberedObject, 3),
  FUN(resizeNumberedObjects, 2),
  FUN(getSizeNumberedObjects, 1),
  FUN(newNumberedObjects, 0),
  FUN(register_namedObjects_Finalizer, 3),
  FUN(register_numberedObjects_Finalizer, 3),
  FUN(register_VecNimArr_Finalizer, 2),
  FUN(register_pointedToBase_Finalizer, 3),
  FUN(register_smartPtrBase_Finalizer, 3),
  FUN(RNimble_Ptr_ManualFinalizer, 1),
  FUN(RNimble_Ptr_CheckAndRunAllDllFinalizers, 2),
  FUN(CountDllObjects, 1),

  // predefinedNimbleList
  FUN(new_EIGEN_EIGENCLASS, 0),
  FUN(EIGEN_EIGENCLASS_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(EIGEN_EIGENCLASS_castDerivedPtrPtrToPairOfPtrsSEXP, 1),

  FUN(new_EIGEN_SVDCLASS, 0),
  FUN(EIGEN_SVDCLASS_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(EIGEN_SVDCLASS_castDerivedPtrPtrToPairOfPtrsSEXP, 1),

  FUN(new_OptimResultNimbleList, 0),
  FUN(OptimResultNimbleList_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(OptimResultNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP, 1),
  
  FUN(new_OptimControlNimbleList, 0),
  FUN(OptimControlNimbleList_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(OptimControlNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP, 1),

  FUN(new_NIMBLE_ADCLASS, 0),
  FUN(NIMBLE_ADCLASS_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(NIMBLE_ADCLASS_castDerivedPtrPtrToPairOfPtrsSEXP, 1),

  FUN(new_waicNimbleList, 0),
  FUN(waicNimbleList_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(waicNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP, 1),
  
  FUN(new_waicDetailsNimbleList, 0),
  FUN(waicDetailsNimbleList_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(waicDetailsNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP, 1),

  FUN(new_AGHQuad_params, 0),
  FUN(AGHQuad_params_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(AGHQuad_params_castDerivedPtrPtrToPairOfPtrsSEXP, 1),

  FUN(new_AGHQuad_summary, 0),
  FUN(AGHQuad_summary_castPtrPtrToNamedObjectsPtrSEXP, 1),
  FUN(AGHQuad_summary_castDerivedPtrPtrToPairOfPtrsSEXP, 1),

  {NULL, NULL, 0}
};

/* Reference prototype: */
/* R_CMethodDef CEntries[] = { */
/*   FUN(RegisterNimblePointer, 3), */
/*   {NULL, NULL, 0} */
/* } */

// Something like the following will be generated with each .so/.dll nimble creates.
// It is required to be named R_init_SONAME, so it must be generated for each one.
//
// extern "C"
// void
// R_init_nimble_on_the_fly(DllInfo *dll)
// {
//     R_registerRoutines(dll, CEntries OR NULL, CallEntries, NULL, NULL);
// }

#endif
