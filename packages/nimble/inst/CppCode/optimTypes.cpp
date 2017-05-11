#include <Rmath.h>
#include <math.h>
#include <nimble/EigenTypedefs.h>
#include <nimble/RcppUtils.h>
#include <nimble/Utils.h>
#include <nimble/accessorClasses.h>
#include <nimble/optimTypes.h>
#include <nimble/smartPtrs.h>
#include <iostream>
#undef eval

void OptimResultNimbleList::copyFromSEXP(SEXP S_nimList_) {
    SEXP S_pxData;
    SEXP S_par;
    SEXP S_value;
    SEXP S_counts;
    SEXP S_convergence;
    SEXP S_message;
    SEXP S_hessian;
    RObjectPointer = S_nimList_;
    PROTECT(S_pxData = allocVector(STRSXP, 1));
    SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
    PROTECT(S_par =
                findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("par")));
    PROTECT(S_value = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("value")));
    PROTECT(S_counts = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                      install("counts")));
    PROTECT(S_convergence = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                           install("convergence")));
    PROTECT(S_message = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                       install("message")));
    PROTECT(S_hessian = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                       install("hessian")));
    SEXP_2_NimArr<1>(S_par, par);
    value = SEXP_2_double(S_value);
    SEXP_2_NimArr<1>(S_counts, counts);
    convergence = SEXP_2_int(S_convergence);
    message = STRSEXP_2_string(S_message);
    SEXP_2_NimArr<2>(S_hessian, hessian);
    UNPROTECT(7);
}
SEXP OptimResultNimbleList::copyToSEXP() {
    SEXP S_pxData;
    SEXP S_par;
    SEXP S_value;
    SEXP S_counts;
    SEXP S_convergence;
    SEXP S_message;
    SEXP S_hessian;
    if (!RCopiedFlag) {
        PROTECT(S_pxData = allocVector(STRSXP, 1));
        SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
        PROTECT(S_par = NimArr_2_SEXP<1>(par));
        PROTECT(S_value = double_2_SEXP(value));
        PROTECT(S_counts = NimArr_2_SEXP<1>(counts));
        PROTECT(S_convergence = int_2_SEXP(convergence));
        PROTECT(S_message = string_2_STRSEXP(message));
        PROTECT(S_hessian = NimArr_2_SEXP<2>(hessian));
        defineVar(install("par"), S_par, GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("value"), S_value,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("counts"), S_counts,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("convergence"), S_convergence,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("message"), S_message,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("hessian"), S_hessian,
                  GET_SLOT(RObjectPointer, S_pxData));
        RCopiedFlag = true;
        UNPROTECT(7);
    }
    return (RObjectPointer);
}
void OptimResultNimbleList::createNewSEXP() {
    SEXP S_newNimList;
    SEXP S_listName;
    PROTECT(S_listName = allocVector(STRSXP, 1));
    SET_STRING_ELT(S_listName, 0, mkChar("OptimResultNimbleList"));
    PROTECT(S_newNimList = makeNewNimbleList(S_listName));
    RObjectPointer = S_newNimList;
    UNPROTECT(2);
}
void OptimResultNimbleList::resetFlags() { RCopiedFlag = false; }
OptimResultNimbleList::OptimResultNimbleList() {
    RCopiedFlag = false;
    RObjectPointer = NULL;
    namedObjects["par"] = &par;
    namedObjects["value"] = &value;
    namedObjects["counts"] = &counts;
    namedObjects["convergence"] = &convergence;
    namedObjects["message"] = &message;
    namedObjects["hessian"] = &hessian;
    namedObjects["RObjectPointer"] = &RObjectPointer;
    namedObjects["RCopiedFlag"] = &RCopiedFlag;
}

SEXP new_OptimResultNimbleList() {
    nimSmartPtr<OptimResultNimbleList> *ptrToSmartPtr;
    OptimResultNimbleList *newObj;
    SEXP SptrToSmartPtr;
    newObj = new OptimResultNimbleList;
    ptrToSmartPtr = new nimSmartPtr<OptimResultNimbleList>;
    ptrToSmartPtr->setPtrFromT(newObj);
    PROTECT(SptrToSmartPtr =
                R_MakeExternalPtr(ptrToSmartPtr, R_NilValue, R_NilValue));
    UNPROTECT(1);
    return (OptimResultNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP(
        SptrToSmartPtr));
}

SEXP OptimResultNimbleList_castPtrPtrToNamedObjectsPtrSEXP(SEXP input) {
    return (R_MakeExternalPtr(
        dynamic_cast<NamedObjects *>(reinterpret_cast<OptimResultNimbleList *>(
            *static_cast<void **>(R_ExternalPtrAddr(input)))),
        R_NilValue, R_NilValue));
}

SEXP OptimResultNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP(SEXP input) {
    nimSmartPtrBase *ptrToSmartPtrBase;
    nimSmartPtr<OptimResultNimbleList> *ptrToSmartPtr;
    void *ptrToPtr;
    SEXP SptrToSmartPtrBase;
    SEXP SptrToPtr;
    SEXP Sans;
    ptrToSmartPtr = static_cast<nimSmartPtr<OptimResultNimbleList> *>(
        R_ExternalPtrAddr(input));
    ptrToSmartPtrBase = dynamic_cast<nimSmartPtrBase *>(ptrToSmartPtr);
    ptrToPtr = ptrToSmartPtr->getVoidPtrToRealPtr();
    PROTECT(SptrToSmartPtrBase =
                R_MakeExternalPtr(ptrToSmartPtrBase, R_NilValue, R_NilValue));
    PROTECT(SptrToPtr = R_MakeExternalPtr(ptrToPtr, R_NilValue, R_NilValue));
    PROTECT(Sans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(Sans, 0, SptrToSmartPtrBase);
    SET_VECTOR_ELT(Sans, 1, SptrToPtr);
    UNPROTECT(3);
    return (Sans);
}
