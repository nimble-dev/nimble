#include <nimble/nimOptim.h>

nimSmartPtr<OptimResultNimbleList> nimFakeOptim(const NimArr<1, double>& par) {
    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;

    // This is fake.
    result->par = par;
    result->value = 0.0;
    result->counts.setSize(0);
    result->convergence = 0;
    result->message = "";
    result->hessian.setSize(0, 0);

    return result;
}
