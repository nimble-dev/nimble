#include <R_ext/Applic.h>
#include <nimble/nimOptim.h>

// This is defined in R_ext/Applic.h:
// typedef double optimfn(int, double *, void *);

nimSmartPtr<OptimResultNimbleList> nimFakeOptim(NimArr<1, double>& par,
                                                NimObjectiveFn fn) {
    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;

    // This is fake.
    result->par = par;
    result->value = fn(par);
    result->counts.setSize(0);
    result->convergence = 0;
    result->message = "";
    result->hessian.setSize(0, 0);

    return result;
}
