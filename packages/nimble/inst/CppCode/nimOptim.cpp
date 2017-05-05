#include <R_ext/Applic.h>
#include <cmath>
#include <limits>
#include <nimble/nimOptim.h>
#include <string.h>

// This is defined in R_ext/Applic.h:
// typedef double optimfn(int, double *, void *);

// This wraps a NimObjectiveFn as an R optimfn.
static double optimfn_wrapper(int n, double* par, void* fn) {
    NimArr<1, double> par_copy;
    par_copy.setSize(n);
    memcpy(par_copy.getPtr(), par, n * sizeof(double));
    return ((NimObjectiveFn*)fn)(par_copy);
}

// This attempts to match the behavior of optim() defined in the documentation
//   https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
// and in the reference implementation
//   https://svn.r-project.org/R/trunk/src/library/stats/R/optim.R
//   https://svn.r-project.org/R/trunk/src/library/stats/src/optim.c
//   https://svn.r-project.org/R/trunk/src/include/R_ext/Applic.h
nimSmartPtr<OptimResultNimbleList> nimOptim(NimArr<1, double>& par,
                                            NimObjectiveFn fn) {
    const int n = par.dimSize(0);

    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;
    result->par = par;
    result->counts.setSize(2);
    result->hessian.setSize(n, n);

    // Use Nelder-Mead by default.
    double* Bvec = par.getPtr();
    double* X = result->par.getPtr();
    double* Fmin = &(result->value);
    int* fail = &(result->convergence);
    double abstol = -INFINITY;
    double intol = std::sqrt(std::numeric_limits<double>::epsilon());
    void* ex = (void*)(fn);
    double alpha = 1.0;
    double bet = 0.5;
    double gamm = 2.0;
    int trace = 0;
    int* fncount = result->counts.getPtr();
    int maxit = 100;
    nmmin(n, Bvec, X, Fmin, optimfn_wrapper, fail, abstol, intol,
          ex, alpha, bet, gamm, trace, fncount, maxit);

    return result;
}

// DEPRECATED
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
