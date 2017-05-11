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

    // Default control parameters.
    const int trace = 0;
    const int maxit = 100;
    const double abstol = -INFINITY;
    const double reltol = std::sqrt(std::numeric_limits<double>::epsilon());
    const double alpha = 1.0;
    const double beta = 0.5;
    const double gamma = 2.0;
    const int REPORT = 10;
    const int type = 1;
    const int lmm = 5;
    const double factr = 1e7;
    const double pgtol = 0;

    // Parameters common to all methods.
    const int n = par.dimSize(0);
    double* Fmin = &(result->value);
    int* fail = &(result->convergence);
    void* ex = this;
    int* fncount = &(result->counts[0]);
    int* grcount = &(result->counts[0]);

    if (method == "Nelder-Mead") {
        double* Bvec = par.getPtr();
        double* X = result->par.getPtr();
        nmmin(n, Bvec, X, Fmin, NimOptimProblem::fn, fail, abstol, reltol, ex,
              alpha, beta, gamma, trace, fncount, maxit);
    } else if (method == "BFGS") {
        double* b = par.getPtr();
        std::vector<int> mask(n, 1);
        vmmin(n, b, Fmin, NimOptimProblem::fn, NimOptimProblem::gr, maxit,
              trace, mask.data(), abstol, reltol, REPORT, ex, fncount, grcount,
              fail);
    } else if (method == "CG") {
        double* Bvec = par.getPtr();
        double* X = result->par.getPtr();
        int type = 1;
        cgmin(n, Bvec, X, Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail,
              abstol, reltol, ex, type, trace, fncount, grcount, maxit);
    } else if (method == "L-BFGS-B") {
        double* x = result->par.getPtr();
        std::vector<double> lower(n, -INFINITY);
        std::vector<double> upper(n, INFINITY);
        std::vector<int> nbd(n, 0);  // 0 means no active constraints.
        char msg[60];
        lbfgsb(n, lmm, x, lower.data(), upper.data(), nbd.data(), Fmin,
               NimOptimProblem::fn, NimOptimProblem::gr, fail, ex, factr, pgtol,
               fncount, grcount, maxit, msg, trace, REPORT);
        result->message = msg;
    } else {
        NIMERROR("Unknown method: %s", method.c_str());
    }

    return result;
}
