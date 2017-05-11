// DO NOT EDIT - This file was semiautomatically generated.
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

#endif  // __NIMBLE_OPTIMTYPES_H
