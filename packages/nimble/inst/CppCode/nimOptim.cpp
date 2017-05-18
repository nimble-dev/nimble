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

nimSmartPtr<OptimControlNimbleList> nimOptimDefaultControl() {
    nimSmartPtr<OptimControlNimbleList> control = new OptimControlNimbleList;
    control->trace = 0;
    control->parscale.initialize(1.0, true, 1);
    control->ndeps.initialize(1e-3, true, 1);
    control->maxit = NA_INTEGER;  // Context-dependent.
    control->abstol = -INFINITY;
    control->reltol = std::sqrt(std::numeric_limits<double>::epsilon());
    control->alpha = 1.0;
    control->beta = 0.5;
    control->gamma = 2.0;
    control->REPORT = NA_INTEGER;  // Context-dependent.
    control->type = 1;
    control->lmm = 5;
    control->factr = 1e7;
    control->pgtol = 0;
    control->tmax = 10;
    control->temp = 10.0;
    return control;
}

// This attempts to match the behavior of optim() defined in the documentation
//   https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
// and in the reference implementation
//   https://svn.r-project.org/R/trunk/src/library/stats/R/optim.R
//   https://svn.r-project.org/R/trunk/src/library/stats/src/optim.c
//   https://svn.r-project.org/R/trunk/src/include/R_ext/Applic.h
nimSmartPtr<OptimResultNimbleList> NimOptimProblem::solve(
    NimArr<1, double>& par, const std::string& method, NimArr<1, double>& lower,
    NimArr<1, double>& upper, nimSmartPtr<OptimControlNimbleList> control,
    bool hessian) {
    NIM_ASSERT(!par.isMap(), "Internal error: failed to handle mapped NimArr");
    const int n = par.dimSize(0);

    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;
    result->par = par;
    result->counts.initialize(NA_INTEGER, true, 2);
    if (hessian) {
        result->hessian.initialize(NA_REAL, true, n, n);
    }

    // Set context-dependent default control values.
    if (control->maxit == NA_INTEGER) {
        if (method == "Nelder-Mead") {
            control->maxit = 500;
        } else if (method == "SANN") {
            control->maxit = 10000;
        } else {
            control->maxit = 100;
        }
    }
    if (control->REPORT == NA_INTEGER) {
        if (method == "SANN") {
            control->REPORT = 100;
        } else {
            control->REPORT = 10;
        }
    }

    // Parameters common to all methods.
    double* dpar = par.getPtr();
    double* X = result->par.getPtr();
    double* Fmin = &(result->value);
    int* fail = &(result->convergence);
    void* ex = this;
    int* fncount = &(result->counts[0]);
    int* grcount = &(result->counts[1]);

    if (method == "Nelder-Mead") {
        nmmin(n, dpar, X, Fmin, NimOptimProblem::fn, fail, control->abstol,
              control->reltol, ex, control->alpha, control->beta,
              control->gamma, control->trace, fncount, control->maxit);
    } else if (method == "BFGS") {
        std::vector<int> mask(n, 1);
        vmmin(n, dpar, Fmin, NimOptimProblem::fn, NimOptimProblem::gr,
              control->maxit, control->trace, mask.data(), control->abstol,
              control->reltol, control->REPORT, ex, fncount, grcount, fail);
        result->par = par;
    } else if (method == "CG") {
        cgmin(n, dpar, X, Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail,
              control->abstol, control->reltol, ex, control->type,
              control->trace, fncount, grcount, control->maxit);
    } else if (method == "L-BFGS-B") {
        if (lower.dimSize(0) == 1) lower.initialize(lower[0], true, n);
        if (upper.dimSize(0) == 1) upper.initialize(upper[0], true, n);
        std::vector<int> nbd(n, 0);  // 0 means no active constraints.
        char msg[60];
        lbfgsb(n, control->lmm, X, lower.getPtr(), upper.getPtr(), nbd.data(),
               Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail, ex,
               control->factr, control->pgtol, fncount, grcount, control->maxit,
               msg, control->trace, control->REPORT);
        result->message = msg;
    } else if (method == "SANN") {
        const int trace = control->trace ? control->REPORT : 0;
        samin(n, dpar, Fmin, NimOptimProblem::fn, control->maxit, control->tmax,
              control->temp, trace, ex);
        result->par = par;
        result->counts[0] = control->maxit;
    } else {
        NIMERROR("Unknown method: %s", method.c_str());
    }

    // Compute Hessian.
    if (hessian) {
        Rf_warning("Hessian computation is not implemented");  // TODO
    }

    return result;
}
