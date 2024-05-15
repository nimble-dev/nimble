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
#include <nimble/accessorClasses.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <limits>

double NimOptimProblem::fn(int n, double* par, void* ex) {
    NimOptimProblem* problem = static_cast<NimOptimProblem*>(ex);
    problem->par_.setSize(n, false, false);
    double* problem_par = problem->par_.getPtr();
    double* problem_parscale = problem->working_parscale.getPtr();
    for(int i = 0; i < n; ++i)
        *problem_par++ = par[i] * problem_parscale[i];
    double ans = problem->function() / problem->control_->fnscale;
    if(isnan(ans)) ans = std::numeric_limits<double>::infinity();
    return ans;
}

SEXP CALL_NimOptimProblem_fn(SEXP Sx, SEXP Sextptr) {
    NimArr<1, double> x;
    SEXP_2_NimArr<1>(Sx, x);
    void *ex = R_ExternalPtrAddr(Sextptr);
    double ans = NimOptimProblem::fn(x.size(), x.getPtr(), ex);
    SEXP S_returnValue;
    PROTECT(S_returnValue = double_2_SEXP(ans)); // CHECK what is needed
    UNPROTECT(1);
    return S_returnValue;
}

void NimOptimProblem::gr(int n, double* par, double* ans, void* ex) {
    NimOptimProblem* problem = static_cast<NimOptimProblem*>(ex);
    problem->par_.setSize(n, false, false);
    double* problem_par = problem->par_.getPtr();
    double* problem_parscale = problem->working_parscale.getPtr();
    for(int i = 0; i < n; ++i)
        *problem_par++ = par[i] * problem_parscale[i];
    problem->ans_.setSize(n, false, false);
    problem->gradient();
    for (int i = 0; i < n; ++i) {
        ans[i] = problem->ans_[i] / problem->control_->fnscale;
    }
}

SEXP CALL_NimOptimProblem_gr(SEXP Sx, SEXP Sextptr) {
    NimArr<1, double> x;
    SEXP_2_NimArr<1>(Sx, x);
    void *ex = R_ExternalPtrAddr(Sextptr);
    SEXP S_returnValue = PROTECT(Rf_allocVector(REALSXP, x.size()));
    NimOptimProblem::gr(x.size(), x.getPtr(), REAL(S_returnValue), ex);
    UNPROTECT(1);
    return S_returnValue;
}

void NimOptimProblem::he(int n, double* par, double* ans, void* ex) {
    NimOptimProblem* problem = static_cast<NimOptimProblem*>(ex);
    problem->par_.setSize(n, false, false);
    double* problem_par = problem->par_.getPtr();
    double* problem_parscale = problem->working_parscale.getPtr();
    for(int i = 0; i < n; ++i)
        *problem_par++ = par[i] * problem_parscale[i];
    problem->ans_hessian_.setSize(n, n, false, false);
    problem->hessian_callback();
    for (int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            ans[i + j*n] = problem->ans_hessian_(i, j) / problem->control_->fnscale;
        }
    }
}

SEXP CALL_NimOptimProblem_he(SEXP Sx, SEXP Sextptr) {
    NimArr<1, double> x;
    SEXP_2_NimArr<1>(Sx, x);
    void *ex = R_ExternalPtrAddr(Sextptr);
    size_t n = x.size();
    SEXP S_returnValue = PROTECT(Rf_allocVector(REALSXP, n*n));
    NimOptimProblem::he(x.size(), x.getPtr(), REAL(S_returnValue), ex);
    SEXP Sdim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(Sdim)[0] = INTEGER(Sdim)[1] = n;
    Rf_setAttrib(S_returnValue, R_DimSymbol, Sdim);
    UNPROTECT(2);
    return S_returnValue;
}

nimSmartPtr<OptimControlNimbleList> nimOptimDefaultControl() {
    nimSmartPtr<OptimControlNimbleList> control = new OptimControlNimbleList;
    control->trace = 0;
    control->fnscale = 1;
    control->parscale.initialize(NA_REAL, true, 1);
    control->ndeps.initialize(NA_REAL, true, 1);
    control->maxit = NA_INTEGER;  // Context-dependent.
    control->abstol = -INFINITY;
    control->reltol = std::sqrt(std::numeric_limits<double>::epsilon());
    control->alpha = 1.0;
    control->beta = 0.5;
    control->gamma = 2.0;
    control->REPORT = NA_INTEGER;
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
    NimArr<1, double>& init_par) {
    NIM_ASSERT1(!init_par.isMap(), "Internal error: failed to handle mapped NimArr");
    const int n = init_par.dimSize(0);
    // Set context-dependent default control_ values.
    int working_maxit = control_->maxit;
    if (working_maxit == NA_INTEGER) {
        if (method_ == "Nelder-Mead") {
            working_maxit = 500;
        } else {
            working_maxit = 100;
        }
    }
    if((control_->parscale.size() == 1) && (R_IsNA(control_->parscale[0]))) {
        working_parscale.initialize(1.0, true, n);
    } else {
        working_parscale = control_->parscale;
    }
    bool parscale_error = (working_parscale.size() != n);
    if((control_->ndeps.size() == 1) && (R_IsNA(control_->ndeps[0]))) {
        working_ndeps.initialize(1e-3, true, n);
    } else {
        working_ndeps = control_->ndeps;
    }

    if(working_ndeps.size() != n) {
        if(parscale_error)
            NIMERROR("In compiled optim (aka nimOptim) call: lengths for control parameters parscale and ndeps must equal length(par).");
        else
            NIMERROR("In compiled optim (aka nimOptim) call: length for control parameter ndeps must equal length(par).");
    } else {
        if(parscale_error)
            NIMERROR("In compiled optim (aka nimOptim) call: length for control parameters parscale must equal length(par).");
    }
    int working_REPORT = control_->REPORT;
    if(working_REPORT == NA_INTEGER) {
        if (method_ == "BFGS" || method_ == "L-BFGS-B")
            working_REPORT = 10;
        if (method_ == "SANN") // SANN is not even supported, but I'm including this so we don't forget if we ever do support SANN
            working_REPORT = 100;
    }

    NimArr<1, double> par = init_par;
    for(int i = 0; i < n; ++i)
        par[i] /= working_parscale[i];

    nimSmartPtr<OptimResultNimbleList> result = new OptimResultNimbleList;
    result->par = par;
    result->counts.initialize(NA_INTEGER, true, 2);
    if (hessian_) {
        result->hessian.initialize(NA_REAL, true, n, n);
    }

    // Parameters common to all methods.
    double* dpar = par.getPtr();
    double* X = result->par.getPtr();
    double* Fmin = &(result->value);
    int* fail = &(result->convergence);
    void* ex = this;
    int* fncount = &(result->counts[0]);
    int* grcount = &(result->counts[1]);

    bool calc_hessian_after = hessian_;

    if (method_ == "Nelder-Mead") {
        nmmin(n, dpar, X, Fmin, NimOptimProblem::fn, fail, control_->abstol,
              control_->reltol, ex, control_->alpha, control_->beta,
              control_->gamma, control_->trace, fncount, working_maxit);
    } else if (method_ == "BFGS") {
        std::vector<int> mask(n, 1);
        for(size_t iii = 0; iii < n; ++iii) X[iii] = dpar[iii];
        vmmin(n, X, Fmin, NimOptimProblem::fn, NimOptimProblem::gr,
              working_maxit, control_->trace, mask.data(), control_->abstol,
              control_->reltol, working_REPORT, ex, fncount, grcount, fail);
    } else if (method_ == "CG") {
        cgmin(n, dpar, X, Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail,
              control_->abstol, control_->reltol, ex, control_->type,
              control_->trace, fncount, grcount, working_maxit);
    } else {
        // from here on, all methods need lower_ and upper_ set up.
        if (lower_.dimSize(0) == 1) lower_.initialize(lower_[0], true, n);
        if (upper_.dimSize(0) == 1) upper_.initialize(upper_[0], true, n);
        NIM_ASSERT_SIZE(lower_, n);
        NIM_ASSERT_SIZE(upper_, n);
        if (method_ == "L-BFGS-B") {
            std::vector<int> nbd(n, 0);
            for (int i = 0; i < n; ++i) {
                if (std::isfinite(lower_[i])) nbd[i] |= 1;
                if (std::isfinite(upper_[i])) nbd[i] |= 2;
            }
            char msg[60];
            lbfgsb(n, control_->lmm, X, lower_.getPtr(), upper_.getPtr(),
                   nbd.data(), Fmin, NimOptimProblem::fn, NimOptimProblem::gr, fail,
                   ex, control_->factr, control_->pgtol, fncount, grcount,
                   working_maxit, msg, control_->trace, working_REPORT);
            result->message = msg;
        } else { //"nvm" or other
            SEXP SLANG;
            SEXP SLANGpiece;
            SEXP SANS;
            SEXP Spar;
            SEXP Smethod;
            Spar = PROTECT(NimArr_2_SEXP(par));
            SEXP Slower = PROTECT(NimArr_2_SEXP(lower_));
            SEXP Supper = PROTECT(NimArr_2_SEXP(upper_));
            SEXP Sgr_provided = PROTECT(bool_2_SEXP(gr_provided_));
            SEXP She_provided = PROTECT(bool_2_SEXP(he_provided_));
            SEXP Sreturn_hessian = PROTECT(bool_2_SEXP(hessian_));

            // if (!(*control_).RObjectPointer) control_->createNewSEXP();
            // SEXP Scontrol = PROTECT(Scontrol = control_->copyToSEXP());
            // control_->resetFlags();

            SEXP Scontrol = PROTECT(Rf_allocVector(VECSXP, 6));
            SET_VECTOR_ELT(Scontrol, 0, PROTECT(double_2_SEXP(control_->abstol)));
            SET_VECTOR_ELT(Scontrol, 1, PROTECT(double_2_SEXP(control_->reltol)));
            SET_VECTOR_ELT(Scontrol, 2, PROTECT(int_2_SEXP(control_->maxit)));
            SET_VECTOR_ELT(Scontrol, 3, PROTECT(NimArr_2_SEXP<1>(control_->parscale)));
            SET_VECTOR_ELT(Scontrol, 4, PROTECT(double_2_SEXP(control_->fnscale)));
            SET_VECTOR_ELT(Scontrol, 5, PROTECT(int_2_SEXP(control_->trace)));
            vector<string> newNames(6);
            newNames[0].assign("abstol");
            newNames[1].assign("reltol");
            newNames[2].assign("maxit");
            newNames[3].assign("parscale");
            newNames[4].assign("fnscale");
            newNames[5].assign("trace");
            SEXP SnewNames;
            PROTECT(SnewNames = vectorString_2_STRSEXP(newNames));
            Rf_setAttrib(Scontrol, R_NamesSymbol, SnewNames);

            Smethod = PROTECT(string_2_STRSEXP(method_));
            SLANGpiece = (SLANG = PROTECT(Rf_allocVector(LANGSXP, 10)));
            SET_NEXT_LANG_ARG(SLANGpiece, Rf_install("custom_optim"));
            SET_TAG(SLANGpiece, Rf_install("method"));
            SET_NEXT_LANG_ARG(SLANGpiece, Smethod);
            SET_TAG(SLANGpiece, Rf_install("par"));
            SET_NEXT_LANG_ARG(SLANGpiece, Spar);
            SET_TAG(SLANGpiece, Rf_install("lower"));
            SET_NEXT_LANG_ARG(SLANGpiece, Slower);
            SET_TAG(SLANGpiece, Rf_install("upper"));
            SET_NEXT_LANG_ARG(SLANGpiece, Supper);
            SET_TAG(SLANGpiece, Rf_install("control"));
            SET_NEXT_LANG_ARG(SLANGpiece, Scontrol);

            SET_TAG(SLANGpiece, Rf_install("hessian"));
            SET_NEXT_LANG_ARG(SLANGpiece, Sreturn_hessian);

            SET_TAG(SLANGpiece, Rf_install("use_gr"));
            SET_NEXT_LANG_ARG(SLANGpiece, Sgr_provided);
            SET_TAG(SLANGpiece, Rf_install("use_he"));
            SET_NEXT_LANG_ARG(SLANGpiece, She_provided);

            SET_TAG(SLANGpiece, Rf_install("extptr"));
            SET_NEXT_LANG_ARG(SLANGpiece,
                              PROTECT(R_MakeExternalPtr(this, R_NilValue, R_NilValue)));
            SEXP SnimbleInternalFunctionsEnv =
                PROTECT(Rf_eval(PROTECT(Rf_findVar(Rf_install("nimbleInternalFunctions"),
                                                   R_GlobalEnv)),
                                R_GlobalEnv));

            PutRNGstate();
            SANS = PROTECT(Rf_eval(SLANG, SnimbleInternalFunctionsEnv));
            GetRNGstate();
            SEXP Sresult_par = PROTECT(getListElement(SANS, "par"));
            SEXP Sresult_value = PROTECT(getListElement(SANS, "value"));
            SEXP Sresult_msg = PROTECT(getListElement(SANS, "message"));
            SEXP_2_NimArr<1>(Sresult_par, result->par);
            result->value = SEXP_2_double(Sresult_value);
            if(!Rf_isString(Sresult_msg)) result->message = std::string("");
            else result->message = STRSEXP_2_string(Sresult_msg, 0);
            result->convergence = SEXP_2_int(PROTECT(getListElement(SANS, "convergence")));
            SEXP_2_NimArr<1>(PROTECT(getListElement(SANS, "counts")), result->counts);
            if(hessian_) {
                SEXP Sresult_hessian = PROTECT(getListElement(SANS, "hessian"));
                if( Sresult_hessian != R_NilValue ) {
                    SEXP_2_NimArr<2>(Sresult_hessian, result->hessian);
                    calc_hessian_after = false;
                }
                UNPROTECT(1);
            }
            UNPROTECT(25);
            //        return result;
        }
    }
    //else {
        //  NIMERROR("Unknown method_: %s", method_.c_str());
        //}
    result->value *= control_->fnscale;

    // Compute Hessian.
    // Parameters are still on the optimization scale,
    // i.e. divided by parscale
    if (calc_hessian_after) {
        NimOptimProblem::calc_hessian(result->par, result->hessian);
    }

    for(int i = 0; i < n; ++i)
        result->par[i] *= working_parscale[i];
    return result;
}

void NimOptimProblem::calc_hessian(NimArr<1, double> par,
                                   NimArr<2, double> &hessian) {
    Rprintf("entering calc_hessian\n");
    for(int i = 0; i < par.size(); ++i) std::cout<<par[i]<<" ";
    std::cout<<std::endl;
    // Notice par is copied but hessian is by reference.
    double *ndeps = working_ndeps.getPtr();
    double *parscale = working_parscale.getPtr();
    //  double ndeps = 0.001; //  This should be obtained from control list but is hard-wired for now.
  double epsilon;
  int n = par.dimSize(0);
  void *ex = this;
  double* dpar = par.getPtr(); // This is already divided by parscale
  NimArr<1, double> ansUpper;
  NimArr<1, double> ansLower;
  ansUpper.setSize(n, false, false);
  ansLower.setSize(n, false, false);
  hessian.setSize(n, n, false, false);
  int i, j;
  for(i = 0; i < n; ++i) {
    std::cout<<"i = "<<i<<std::endl;
      // It is strange to divide ndeps by parscale, but that's
      // exactly what R's C code for optimhess does
      epsilon = ndeps[i] / parscale[i];
    std::cout<<"epsilon = "<<epsilon<<std::endl;
    dpar[i] += epsilon;
    gr(n, dpar, ansUpper.getPtr(), ex);
    dpar[i] -= 2*epsilon;
    gr(n, dpar, ansLower.getPtr(), ex);
    for(j = 0; j < n; ++j) {
      // Following R's optimhess, we want to multiply by fnscale here to return answer to original scale.
      hessian(i, j) = control_->fnscale * (ansUpper[j] - ansLower[j]) / (2.*epsilon*parscale[i]*parscale[j]);
      std::cout<<ansUpper[j]<<" "<<ansLower[j]<<" "<<parscale[i]<<" "<<parscale[j]<<" "<<hessian(i,j)<<std::endl;
    }
    dpar[i] += epsilon;
  }
  // Follow the "symmetrize" step of R's optimhess in stats, i.e. average the two relevant elements.
  for(i = 0; i < n; ++i) {
    for(j = 0; j < i; ++j) {
      double tmp = 0.5* (hessian(i,j) + hessian(j,i));
      hessian(i, j) = hessian(j, i) = tmp;
      std::cout<<hessian(i,j)<<" ";
    }
    std::cout<<std::endl;
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
