// DO NOT EDIT BY HAND.
// This file was automatically generated by nimble/packages/generateStaticCode.R

#ifndef __NIMBLE_OPTIMTYPES_H
#define __NIMBLE_OPTIMTYPES_H

#include <Rinternals.h>
#include <nimble/NamedObjects.h>
#include <nimble/NimArr.h>
#include <nimble/accessorClasses.h>
#include <nimble/nimDists.h>
#include <nimble/smartPtrs.h>
#undef eval

class OptimResultNimbleList : public NamedObjects, public pointedToBase {
   public:
    NimArr<1, double> par;
    double value;
    NimArr<1, int> counts;
    int convergence;
    std::string message;
    NimArr<2, double> hessian;
    SEXP RObjectPointer;
    bool RCopiedFlag;
    void copyFromSEXP(SEXP S_nimList_);
    SEXP copyToSEXP();
    void createNewSEXP();
    void resetFlags();
    OptimResultNimbleList();
};

extern "C" SEXP new_OptimResultNimbleList();

extern "C" SEXP OptimResultNimbleList_castPtrPtrToNamedObjectsPtrSEXP(
    SEXP input);

extern "C" SEXP OptimResultNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP(
    SEXP input);

class OptimControlNimbleList : public NamedObjects, public pointedToBase {
   public:
    int trace;
    double fnscale;
    NimArr<1, double> parscale;
    NimArr<1, double> ndeps;
    int maxit;
    double abstol;
    double reltol;
    double alpha;
    double beta;
    double gamma;
    int REPORT;
    int type;
    int lmm;
    double factr;
    double pgtol;
    double temp;
    int tmax;
    SEXP RObjectPointer;
    bool RCopiedFlag;
    void copyFromSEXP(SEXP S_nimList_);
    SEXP copyToSEXP();
    void createNewSEXP();
    void resetFlags();
    OptimControlNimbleList();
};

extern "C" SEXP new_OptimControlNimbleList();

extern "C" SEXP OptimControlNimbleList_castPtrPtrToNamedObjectsPtrSEXP(
    SEXP input);

extern "C" SEXP OptimControlNimbleList_castDerivedPtrPtrToPairOfPtrsSEXP(
    SEXP input);

#endif  // __NIMBLE_OPTIMTYPES_H
