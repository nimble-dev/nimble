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
    // std::cout<<"fn ";
    // int n5 = n > 5 ? 5 : n;
    // for(int i = 0; i < n5; ++i) std::cout<<par[i]<<"\t"; 
    double ans = problem->function() / problem->control_->fnscale;
    // std::cout<<"ans = "<<ans<<"\t";
    if(isnan(ans)) ans = std::numeric_limits<double>::infinity();
    // std::cout<<"returning "<<ans<<std::endl;
    return ans;
}

void NimOptimProblem::gr(int n, double* par, double* ans, void* ex) {
    NimOptimProblem* problem = static_cast<NimOptimProblem*>(ex);
    problem->par_.setSize(n, false, false);
    std::copy(par, par + n, problem->par_.getPtr());
    // std::cout<<"gr ";
    // int n5 = n > 5 ? 5 : n;
    // for(int i = 0; i < n5; ++i) std::cout<<par[i]<<"\t"; 
    problem->ans_.setSize(n, false, false);
    problem->gradient();
    for (int i = 0; i < n; ++i) {
        ans[i] = problem->ans_[i] / problem->control_->fnscale;
    }
    // std::cout<<"ans: ";
    // for(int i = 0; i < n5; ++i) std::cout<<ans[i]<<"\t"; 
    // std::cout<<std::endl;
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
	/* Shouldn't dpar be copied to X and X passed to the function? */
	/* It looks like this method replaces X with the arg max result.*/
	/* Oh, I see that instead, after vmmin, we have result->par = par,*/
	/* This will place the last evaluated parameters in the result.*/
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
	NimOptimProblem::calc_hessian(result->par, result->hessian);
    }
    return result;
}

void NimOptimProblem::calc_hessian(NimArr<1, double> par,
				   NimArr<2, double> &hessian) {
  // Notice par is copied but hessian is by reference.
  double ndeps = 0.001; //  This should be obtained from control list but is hard-wired for now.
  double epsilon = ndeps;
  int n = par.dimSize(0);
  void *ex = this;
  double* dpar = par.getPtr();
  NimArr<1, double> ansUpper;
  NimArr<1, double> ansLower;
  ansUpper.setSize(n, false, false);
  ansLower.setSize(n, false, false);
  hessian.setSize(n, n, false, false);
  int i, j;
  for(i = 0; i < n; ++i) {
    dpar[i] += epsilon;
    gr(n, dpar, ansUpper.getPtr(), ex);
    dpar[i] -= 2*epsilon;
    gr(n, dpar, ansLower.getPtr(), ex);
    for(j = 0; j < n; ++j) {
      // Following R's optimhess, we want to multiply by fnscale here to return answer to original scale.
      hessian(i, j) = control_->fnscale * (ansUpper[j] - ansLower[j]) / (2.*epsilon);
    }
    dpar[i] += epsilon;
  }
  // Follow the "symmetrize" step of R's optimhess in stats, i.e. average the two relevant elements.
  for(i = 0; i < n; ++i) {
    for(j = 0; j < i; ++j) {
      double tmp = 0.5* (hessian(i,j) + hessian(j,i));
      hessian(i, j) = hessian(j, i) = tmp;
    }
  }
}

/*
Notes on R's implementation of Hessian.
---------------------------------------
Source code of optim (type "optim" in R) uses .External2(C_optim, ...) followed by .External2(C_optimhess, ...).
What each of these calls can be seen by stats:::C_optim and stats:::C_optimhess.
They call "optim" and "optimhess" respectively.  These are in base package stats.
Within R source code, these files are in src/library/stats/src.

The .External2 format passes the R call in a standard four arguments, of which the third is a list of actual arguments.

(We see that parscale is applied internally, which means if we want to imitate it's behavior, we'll have to implement it directly.)

The actual optimization algorithms such as nmmin and vmmin can be found in src/appl/optim.c.

optimhess calculates the Hessian.  We see that it uses finite element method for the gradient at points p + eps and p - eps, where eps comes from the
 control list argument ndeps, which acccording to help(optim) can be user supplied or defaults to 0.001.  

The gradient is evaluated via fmingr, which either uses the supplied gradient function or uses finite element differences of +/- eps.

  For the finite element case, this means in effect that the function (fminfn) is evaluated at p + 2*eps, p [twice, once in each call to fmingr], and p - 2*eps.

The exception is for a case with bounds (L-BFGS-B), in which case, inside fmingr, the p + eps and p-eps are reduced to upper boundary or lower boundary, 
respectively, and the corresponding epsilons are adjusted. 

fminfn and fmingr are defined in the same file as optim and optimhess.  (Our gr function above mimics the behavior of fmingr.)

These functions wrap calls to the R evaluator for the provided 
objective and gradient functions, respectively.  Within these functions, parscale and fnscale are applied. (Our fn and gr use fnscale but not parscale.)
*/
