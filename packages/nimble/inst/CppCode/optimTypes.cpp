// DO NOT EDIT BY HAND.
// This file was automatically generated by nimble/packages/generateStaticCode.R

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

void OptimControlNimbleList::copyFromSEXP(SEXP S_nimList_) {
    SEXP S_pxData;
    SEXP S_trace;
    SEXP S_parscale;
    SEXP S_ndeps;
    SEXP S_maxIt;
    SEXP S_abstol;
    SEXP S_reltol;
    SEXP S_alpha;
    SEXP S_beta;
    SEXP S_gamma;
    SEXP S_REPORT;
    SEXP S_type;
    SEXP S_lmm;
    SEXP S_factr;
    SEXP S_pgtol;
    SEXP S_temp;
    SEXP S_tmax;
    RObjectPointer = S_nimList_;
    PROTECT(S_pxData = allocVector(STRSXP, 1));
    SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
    PROTECT(S_trace = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("trace")));
    PROTECT(S_parscale = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                        install("parscale")));
    PROTECT(S_ndeps = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("ndeps")));
    PROTECT(S_maxIt = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("maxit")));
    PROTECT(S_abstol = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                      install("abstol")));
    PROTECT(S_reltol = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                      install("reltol")));
    PROTECT(S_alpha = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("alpha")));
    PROTECT(S_beta = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                    install("beta")));
    PROTECT(S_gamma = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("gamma")));
    PROTECT(S_REPORT = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                      install("REPORT")));
    PROTECT(S_type = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                    install("type")));
    PROTECT(S_lmm =
                findVarInFrame(GET_SLOT(S_nimList_, S_pxData), install("lmm")));
    PROTECT(S_factr = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("factr")));
    PROTECT(S_pgtol = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                     install("pgtol")));
    PROTECT(S_temp = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                    install("temp")));
    PROTECT(S_tmax = findVarInFrame(GET_SLOT(S_nimList_, S_pxData),
                                    install("tmax")));
    trace = SEXP_2_int(S_trace);
    SEXP_2_NimArr<1>(S_parscale, parscale);
    SEXP_2_NimArr<1>(S_ndeps, ndeps);
    maxit = SEXP_2_int(S_maxIt);
    abstol = SEXP_2_double(S_abstol);
    reltol = SEXP_2_double(S_reltol);
    alpha = SEXP_2_double(S_alpha);
    beta = SEXP_2_double(S_beta);
    gamma = SEXP_2_double(S_gamma);
    REPORT = SEXP_2_int(S_REPORT);
    type = SEXP_2_int(S_type);
    lmm = SEXP_2_int(S_lmm);
    factr = SEXP_2_double(S_factr);
    pgtol = SEXP_2_double(S_pgtol);
    temp = SEXP_2_double(S_temp);
    tmax = SEXP_2_int(S_tmax);
    UNPROTECT(17);
}
SEXP OptimControlNimbleList::copyToSEXP() {
    SEXP S_pxData;
    SEXP S_trace;
    SEXP S_parscale;
    SEXP S_ndeps;
    SEXP S_maxIt;
    SEXP S_abstol;
    SEXP S_reltol;
    SEXP S_alpha;
    SEXP S_beta;
    SEXP S_gamma;
    SEXP S_REPORT;
    SEXP S_type;
    SEXP S_lmm;
    SEXP S_factr;
    SEXP S_pgtol;
    SEXP S_temp;
    SEXP S_tmax;
    if (!RCopiedFlag) {
        PROTECT(S_pxData = allocVector(STRSXP, 1));
        SET_STRING_ELT(S_pxData, 0, mkChar(".xData"));
        PROTECT(S_trace = int_2_SEXP(trace));
        PROTECT(S_parscale = NimArr_2_SEXP<1>(parscale));
        PROTECT(S_ndeps = NimArr_2_SEXP<1>(ndeps));
        PROTECT(S_maxIt = int_2_SEXP(maxit));
        PROTECT(S_abstol = double_2_SEXP(abstol));
        PROTECT(S_reltol = double_2_SEXP(reltol));
        PROTECT(S_alpha = double_2_SEXP(alpha));
        PROTECT(S_beta = double_2_SEXP(beta));
        PROTECT(S_gamma = double_2_SEXP(gamma));
        PROTECT(S_REPORT = int_2_SEXP(REPORT));
        PROTECT(S_type = int_2_SEXP(type));
        PROTECT(S_lmm = int_2_SEXP(lmm));
        PROTECT(S_factr = double_2_SEXP(factr));
        PROTECT(S_pgtol = double_2_SEXP(pgtol));
        PROTECT(S_temp = double_2_SEXP(temp));
        PROTECT(S_tmax = int_2_SEXP(tmax));
        defineVar(install("trace"), S_trace,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("parscale"), S_parscale,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("ndeps"), S_ndeps,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("maxit"), S_maxIt,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("abstol"), S_abstol,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("reltol"), S_reltol,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("alpha"), S_alpha,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("beta"), S_beta, GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("gamma"), S_gamma,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("REPORT"), S_REPORT,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("type"), S_type, GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("lmm"), S_lmm, GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("factr"), S_factr,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("pgtol"), S_pgtol,
                  GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("temp"), S_temp, GET_SLOT(RObjectPointer, S_pxData));
        defineVar(install("tmax"), S_tmax, GET_SLOT(RObjectPointer, S_pxData));
        RCopiedFlag = true;
        UNPROTECT(17);
    }
    return (RObjectPointer);
}
void OptimControlNimbleList::createNewSEXP() {
    SEXP S_newNimList;
    SEXP S_listName;
    PROTECT(S_listName = allocVector(STRSXP, 1));
    SET_STRING_ELT(S_listName, 0, mkChar("OptimControlNimbleList"));
    PROTECT(S_newNimList = makeNewNimbleList(S_listName));
    RObjectPointer = S_newNimList;
    UNPROTECT(2);
}
void OptimControlNimbleList::resetFlags() { RCopiedFlag = false; }
OptimControlNimbleList::OptimControlNimbleList() {
    RCopiedFlag = false;
    RObjectPointer = NULL;
    namedObjects["trace"] = &trace;
    namedObjects["parscale"] = &parscale;
    namedObjects["ndeps"] = &ndeps;
    namedObjects["maxit"] = &maxit;
    namedObjects["abstol"] = &abstol;
    namedObjects["reltol"] = &reltol;
    namedObjects["alpha"] = &alpha;
    namedObjects["beta"] = &beta;
    namedObjects["gamma"] = &gamma;
    namedObjects["REPORT"] = &REPORT;
    namedObjects["type"] = &type;
    namedObjects["lmm"] = &lmm;
    namedObjects["factr"] = &factr;
    namedObjects["pgtol"] = &pgtol;
    namedObjects["temp"] = &temp;
    namedObjects["tmax"] = &tmax;
    namedObjects["RObjectPointer"] = &RObjectPointer;
    namedObjects["RCopiedFlag"] = &RCopiedFlag;
}

SEXP new_OptimControlNimbleList() {
    nimSmartPtr<OptimControlNimbleList> *ptrToSmartPtr;
    OptimControlNimbleList *newObj;
    SEXP SptrToSmartPtr;
    newObj = new OptimControlNimbleList;
    ptrToSmartPtr = new nimSmartPtr<OptimControlNimbleList>;
    ptrToSmartPtr->setPtrFromT(newObj);
    PROTECT(SptrToSmartPtr =
                R_MakeExternalPtr(ptrToSmartPtr, R_NilValue, R_NilValue));
    UNPROTECT(1);
    return (OptimControlNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP(
        SptrToSmartPtr));
}

SEXP OptimControlNimbleList_castPtrPtrToNamedObjectsPtrSEXP(SEXP input) {
    return (R_MakeExternalPtr(
        dynamic_cast<NamedObjects *>(reinterpret_cast<OptimControlNimbleList *>(
            *static_cast<void **>(R_ExternalPtrAddr(input)))),
        R_NilValue, R_NilValue));
}

SEXP OptimControlNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP(SEXP input) {
    nimSmartPtrBase *ptrToSmartPtrBase;
    nimSmartPtr<OptimControlNimbleList> *ptrToSmartPtr;
    void *ptrToPtr;
    SEXP SptrToSmartPtrBase;
    SEXP SptrToPtr;
    SEXP Sans;
    ptrToSmartPtr = static_cast<nimSmartPtr<OptimControlNimbleList> *>(
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
