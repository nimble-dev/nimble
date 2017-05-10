#include <R_ext/Applic.h>
#include <nimble/nimOptim.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <limits>

double NimOptimProblem::fn(int n, double* par, void* ex) {
    NimOptimProblem* problem = static_cast<NimOptimProblem*>(ex);
    problem->par_.setSize(n, false, false);
    std::copy(par, par + n, problem->par_.getPtr());
    return problem->function();
}

void NimOptimProblem::gr(int n, double* par, double* ans, void* ex) {
    NimOptimProblem* problem = static_cast<NimOptimProblem*>(ex);
    problem->par_.setSize(n, false, false);
    std::copy(par, par + n, problem->par_.getPtr());
    problem->ans_.setSize(n, false, false);
    problem->gradient();
    std::copy(problem->ans_.getPtr(), problem->ans_.getPtr() + n, ans);
}

// This attempts to match the behavior of optim() defined in the documentation
//   https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
// and in the reference implementation
//   https://svn.r-project.org/R/trunk/src/library/stats/R/optim.R
//   https://svn.r-project.org/R/trunk/src/library/stats/src/optim.c
//   https://svn.r-project.org/R/trunk/src/include/R_ext/Applic.h
nimSmartPtr<OptimResultNimbleList> NimOptimProblem::solve(
    NimArr<1, double>& par, const char* method) {
    const int n = par.dimSize(0);
    NimArr<1, double> par_nomap = par;

    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;
    result->par = par;
    result->counts.initialize(NA_INTEGER, true, 2);
    // result->hessian is not set.

    if (strcmp(method, "Nelder-Mead") == 0) {
        double* Bvec = par_nomap.getPtr();
        double* X = result->par.getPtr();
        double* Fmin = &(result->value);
        int* fail = &(result->convergence);
        double abstol = -INFINITY;
        double reltol = std::sqrt(std::numeric_limits<double>::epsilon());
        void* ex = this;
        double alpha = 1.0;
        double bet = 0.5;
        double gamm = 2.0;
        int trace = 0;
        int* fncount = result->counts.getPtr();
        int maxit = 100;
        nmmin(n, Bvec, X, Fmin, NimOptimProblem::fn, fail, abstol, reltol, ex,
              alpha, bet, gamm, trace, fncount, maxit);
    } else {
        NIMERROR("Unknown method: %s", method);
    }

    return result;
}
