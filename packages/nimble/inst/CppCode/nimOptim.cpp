/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

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
    return problem->function() / problem->control_->fnscale;
}

void NimOptimProblem::gr(int n, double* par, double* ans, void* ex) {
    NimOptimProblem* problem = static_cast<NimOptimProblem*>(ex);
    problem->par_.setSize(n, false, false);
    std::copy(par, par + n, problem->par_.getPtr());
    problem->ans_.setSize(n, false, false);
    problem->gradient();
    for (int i = 0; i < n; ++i) {
        ans[i] = problem->ans_[i] / problem->control_->fnscale;
    }
}

nimSmartPtr<OptimControlNimbleList> nimOptimDefaultControl() {
    nimSmartPtr<OptimControlNimbleList> control = new OptimControlNimbleList;
    control->trace = 0;
    control->fnscale = 1;
    control->parscale.initialize(1.0, true, 1);
    control->ndeps.initialize(1e-3, true, 1);
    control->maxit = NA_INTEGER;  // Context-dependent.
    control->abstol = -INFINITY;
    control->reltol = std::sqrt(std::numeric_limits<double>::epsilon());
    control->alpha = 1.0;
    control->beta = 0.5;
    control->gamma = 2.0;
    control->REPORT = 10;
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
    NimArr<1, double>& par) {
    NIM_ASSERT1(!par.isMap(), "Internal error: failed to handle mapped NimArr");
    const int n = par.dimSize(0);

    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;
    result->par = par;
    result->counts.initialize(NA_INTEGER, true, 2);
    if (hessian_) {
        result->hessian.initialize(NA_REAL, true, n, n);
    }

    // Set context-dependent default control_ values.
    if (control_->maxit == NA_INTEGER) {
        if (method_ == "Nelder-Mead") {
            control_->maxit = 500;
        } else {
            control_->maxit = 100;
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

    if (method_ == "Nelder-Mead") {
        nmmin(n, dpar, X, Fmin, NimOptimProblem::fn, fail, control_->abstol,
              control_->reltol, ex, control_->alpha, control_->beta,
              control_->gamma, control_->trace, fncount, control_->maxit);
    } else if (method_ == "BFGS") {
        std::vector<int> mask(n, 1);
        vmmin(n, dpar, Fmin, NimOptimProblem::fn, NimOptimProblem::gr,
              control_->maxit, control_->trace, mask.data(), control_->abstol,
              control_->reltol, control_->REPORT, ex, fncount, grcount, fail);
        result->par = par;
    } else if (method_ == "CG") {
        cgmin(n, dpar, X, Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail,
              control_->abstol, control_->reltol, ex, control_->type,
              control_->trace, fncount, grcount, control_->maxit);
    } else if (method_ == "L-BFGS-B") {
        if (lower_.dimSize(0) == 1) lower_.initialize(lower_[0], true, n);
        if (upper_.dimSize(0) == 1) upper_.initialize(upper_[0], true, n);
        NIM_ASSERT_SIZE(lower_, n);
        NIM_ASSERT_SIZE(upper_, n);
        std::vector<int> nbd(n, 0);
        for (int i = 0; i < n; ++i) {
            if (std::isfinite(lower_[i])) nbd[i] |= 1;
            if (std::isfinite(upper_[i])) nbd[i] |= 2;
        }
        char msg[60];
        lbfgsb(n, control_->lmm, X, lower_.getPtr(), upper_.getPtr(),
               nbd.data(), Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail,
               ex, control_->factr, control_->pgtol, fncount, grcount,
               control_->maxit, msg, control_->trace, control_->REPORT);
        result->message = msg;
    } else {
        NIMERROR("Unknown method_: %s", method_.c_str());
    }
    result->value *= control_->fnscale;

    // Compute Hessian.
    if (hessian_) {
        Rf_warning("Hessian computation is not implemented");  // TODO
    }

    return result;
}
