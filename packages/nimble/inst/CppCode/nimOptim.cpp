#include <R_ext/Applic.h>
#include <nimble/nimOptim.h>
#include <cmath>
#include <limits>

// This wraps a NimObjectiveFn as an R `optimfn` type
//   typedef double optimfn(int, double *, void *);
// defined in
//   https://svn.r-project.org/R/trunk/src/include/R_ext/Applic.h
static double optimfn_wrapper(int n, double* par, void* fn) {
    NimArr<1, double> nim_par;
    nim_par.setSize(n, false, false);
    std::copy(par, par + n, nim_par.getPtr());
    return ((NimObjectiveFn*)fn)(nim_par);
}

nimSmartPtr<OptimResultNimbleList> nimOptim(NimArr<1, double>& par,
                                            NimObjectiveFn fn) {
    return nimOptim_internal(par, optimfn_wrapper, (void*)(fn));
}

// This attempts to match the behavior of optim() defined in the documentation
//   https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
// and in the reference implementation
//   https://svn.r-project.org/R/trunk/src/library/stats/R/optim.R
//   https://svn.r-project.org/R/trunk/src/library/stats/src/optim.c
//   https://svn.r-project.org/R/trunk/src/include/R_ext/Applic.h
nimSmartPtr<OptimResultNimbleList> nimOptim_internal(NimArr<1, double>& par,
                                                     optimfn fn, void* ex) {
    const int n = par.dimSize(0);
    NimArr<1, double> par_nomap = par;

    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;
    result->par = par;
    result->counts.initialize(NA_INTEGER, true, 2);
    // result->hessian is not set.

    // Use Nelder-Mead by default.
    double* Bvec = par_nomap.getPtr();
    double* X = result->par.getPtr();
    double* Fmin = &(result->value);
    int* fail = &(result->convergence);
    double abstol = -INFINITY;
    double reltol = std::sqrt(std::numeric_limits<double>::epsilon());
    double alpha = 1.0;
    double bet = 0.5;
    double gamm = 2.0;
    int trace = 0;
    int* fncount = result->counts.getPtr();
    int maxit = 100;
    nmmin(n, Bvec, X, Fmin, fn, fail, abstol, reltol, ex, alpha, bet, gamm,
          trace, fncount, maxit);

    return result;
}
