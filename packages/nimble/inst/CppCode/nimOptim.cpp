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
    NimArr<1, double>& par, const std::string& method) {
    NIM_ASSERT(!par.isMap(), "Internal error: failed to handle mapped NimArr");

    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;
    result->par = par;
    result->counts.initialize(NA_INTEGER, true, 2);
    // result->hessian is not set.

    // Parameters common to all methods.
    const int n = par.dimSize(0);
    double *Fmin = &(result->value);
    int* fail = &(result->convergence);
    double abstol = -INFINITY;
    double reltol = std::sqrt(std::numeric_limits<double>::epsilon());
    void* ex = this;
    int trace = 0;
    int* fncount = &(result->counts[0]);
    int* grcount = &(result->counts[0]);
    int maxit = 100;

    if (method == "Nelder-Mead") {
        double* Bvec = par.getPtr();
        double* X = result->par.getPtr();
        double alpha = 1.0;
        double bet = 0.5;
        double gamm = 2.0;
        nmmin(n, Bvec, X, Fmin, NimOptimProblem::fn, fail, abstol, reltol, ex,
              alpha, bet, gamm, trace, fncount, maxit);
    } else if (method == "BFGS") {
        double* b = par.getPtr();
        std::vector<int> mask(n, 1);
        int nREPORT = 10;
        vmmin(n, b, Fmin, NimOptimProblem::fn, NimOptimProblem::gr, maxit,
              trace, mask.data(), abstol, reltol, nREPORT, ex, fncount, grcount,
              fail);
    } else if (method == "CG") {
        double* Bvec = par.getPtr();
        double* X = result->par.getPtr();
        int type = 1;
        cgmin(n, Bvec, X, Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail,
              abstol, reltol, ex, type, trace, fncount, grcount, maxit);
    } else {
        NIMERROR("Unknown method: %s", method.c_str());
    }

    return result;
}
