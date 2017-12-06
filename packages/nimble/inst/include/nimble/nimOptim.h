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

#ifndef __NIMBLE_NIMOPTIM_H
#define __NIMBLE_NIMOPTIM_H

#include <nimble/NimArr.h>
#include <nimble/predefinedNimbleLists.h>
#include <nimble/smartPtrs.h>
#include <algorithm>
#include <string>

// ---------------------------------------------------------------------------
// This class and wrapper are equivalent to std::bind(-,-) in C++11, so that:
//   object->method(par) == NimBoundMethod<T>(&T::method, object)(par);
//                       == NimBind(&T::method, object)(par);
//                       == std::bind(&T::method, object)(par);
template <class T>
class NimBoundMethod {
   public:
    NimBoundMethod(double (T::*method)(NimArr<1, double>&), T* object)
        : method_(method), object_(object) {}
    double operator()(NimArr<1, double>& par) {
        return (object_->*method_)(par);
    }

   private:
    double (T::*method_)(NimArr<1, double>&);
    T* object_;
};

template <class T>
inline NimBoundMethod<T> NimBind(double (T::*method)(NimArr<1, double>&),
                                 T* object) {
    return NimBoundMethod<T>(method, object);
}

// ---------------------------------------------------------------------------
// NimOptimProblem class hierarchy

class NimOptimProblem {
   public:
    NimOptimProblem(const std::string& method, NimArr<1, double>& lower,
                    NimArr<1, double>& upper,
                    nimSmartPtr<OptimControlNimbleList> control, bool hessian)
        : method_(method),
          lower_(lower),
          upper_(upper),
          control_(control),
          hessian_(hessian) {}
    nimSmartPtr<OptimResultNimbleList> solve(NimArr<1, double>& par);

   private:
    // These function and gradient callbacks for R's optim() where `this` is
    // passed in as the final argument `void * ex`.
    static double fn(int, double*, void*);
    static void gr(int, double*, double*, void*);

   protected:
    // These are callbacks used internally by fn() and gr().
    virtual double function() = 0;
    virtual void gradient() = 0;

    // Problem parameters.
    const std::string& method_;
    NimArr<1, double>& lower_;
    NimArr<1, double>& upper_;
    nimSmartPtr<OptimControlNimbleList> control_;
    const bool hessian_;

    // Temporaries.
    NimArr<1, double> par_;  // Argument for fn() and gr().
    NimArr<1, double> ans_;  // Result of gradient.
};

template <class Fn>
class NimOptimProblem_Fun : public NimOptimProblem {
   public:
    NimOptimProblem_Fun(Fn fn, const std::string& method,
                        NimArr<1, double>& lower, NimArr<1, double>& upper,
                        nimSmartPtr<OptimControlNimbleList> control,
                        bool hessian)
        : NimOptimProblem(method, lower, upper, control, hessian), fn_(fn) {}

   protected:
    virtual double function() { return fn_(par_); }
    virtual void gradient();  // Uses finite difference approximation.

   private:
    Fn fn_;
};

template <class Fn>
void NimOptimProblem_Fun<Fn>::gradient() {
    const int n = par_.dimSize(0);
    NIM_ASSERT_SIZE(ans_, n);
    NimArr<1, double>& ndeps = control_->ndeps;
    if (ndeps.dimSize(0) == 1) ndeps.initialize(ndeps[0], true, n);
    NIM_ASSERT_SIZE(ndeps, n);
    NimArr<1, double> par_h = par_;
    if (method_ == "L-BFGS-B") {
        // Constrained optimization.
        for (int i = 0; i < n; ++i) {
            const double h_pos = std::min(ndeps[i], upper_[i] - par_[i]);
            const double h_neg = std::min(ndeps[i], par_[i] - lower_[i]);
            par_h[i] = par_[i] + h_pos;
            const double pos = fn_(par_h);
            par_h[i] = par_[i] - h_neg;
            const double neg = fn_(par_h);
            par_h[i] = par_[i];
            ans_[i] = (pos - neg) / (h_pos + h_neg);
        }
    } else {
        // Unconstrained optimization.
        for (int i = 0; i < n; ++i) {
            const double h = ndeps[i];
            par_h[i] = par_[i] + h;
            const double pos = fn_(par_h);
            par_h[i] = par_[i] - h;
            const double neg = fn_(par_h);
            par_h[i] = par_[i];
            ans_[i] = (pos - neg) / (2 * h);
        }
    }
}

template <class Fn, class Gr>
class NimOptimProblem_Fun_Grad : public NimOptimProblem {
   public:
    NimOptimProblem_Fun_Grad(Fn fn, Gr gr, const std::string& method,
                             NimArr<1, double>& lower, NimArr<1, double>& upper,
                             nimSmartPtr<OptimControlNimbleList> control,
                             bool hessian)
        : NimOptimProblem(method, lower, upper, control, hessian),
          fn_(fn),
          gr_(gr) {}

   protected:
    virtual double function() { return fn_(par_); }
    virtual void gradient() { ans_ = gr_(par_); }

   private:
    Fn fn_;
    Gr gr_;
};

// ---------------------------------------------------------------------------
// These nimOptim() and nimOptimDefaultControl() will appear in generated code

nimSmartPtr<OptimControlNimbleList> nimOptimDefaultControl();

template <class Fn, class Gr>
inline nimSmartPtr<OptimResultNimbleList> nimOptim(
    NimArr<1, double>& par, Fn fn, Gr gr, const std::string& method,
    NimArr<1, double>& lower, NimArr<1, double>& upper,
    nimSmartPtr<OptimControlNimbleList> control, bool hessian) {
    return NimOptimProblem_Fun_Grad<Fn, Gr>(fn, gr, method, lower, upper,
                                            control, hessian)
        .solve(par);
}

// This handles the special case where gr is not specified (i.e. is "NULL").
template <class Fn>
inline nimSmartPtr<OptimResultNimbleList> nimOptim(
    NimArr<1, double>& par, Fn fn, const char* gr, const std::string& method,
    NimArr<1, double>& lower, NimArr<1, double>& upper,
    nimSmartPtr<OptimControlNimbleList> control, bool hessian) {
    NIM_ASSERT1(
        strcmp(gr, "NULL") == 0,
        "Internal error: failed to handle gradient argument type in optim()");
    return NimOptimProblem_Fun<Fn>(fn, method, lower, upper, control, hessian)
        .solve(par);
}

// TODO Handle double -> vector conversion in R instead of here.
template <class Fn, class Gr>
inline nimSmartPtr<OptimResultNimbleList> nimOptim(
    NimArr<1, double>& par, Fn fn, Gr gr, const std::string& method,
    double lower, double upper, nimSmartPtr<OptimControlNimbleList> control,
    bool hessian) {
    NimArr<1, double> lower_vector;
    lower_vector.initialize(lower, true, 1);
    NimArr<1, double> upper_vector;
    upper_vector.initialize(upper, true, 1);
    return NimOptimProblem_Fun_Grad<Fn, Gr>(fn, gr, method, lower_vector,
                                            upper_vector, control, hessian)
        .solve(par);
}

// TODO Handle double -> vector conversion in R instead of here.
// This handles the special case where gr is not specified (i.e. is "NULL").
template <class Fn>
inline nimSmartPtr<OptimResultNimbleList> nimOptim(
    NimArr<1, double>& par, Fn fn, const char* gr, const std::string& method,
    double lower, double upper, nimSmartPtr<OptimControlNimbleList> control,
    bool hessian) {
    NIM_ASSERT1(
        strcmp(gr, "NULL") == 0,
        "Internal error: failed to handle gradient argument type in optim()");
    NimArr<1, double> lower_vector;
    lower_vector.initialize(lower, true, 1);
    NimArr<1, double> upper_vector;
    upper_vector.initialize(upper, true, 1);
    return NimOptimProblem_Fun<Fn>(fn, method, lower_vector, upper_vector,
                                   control, hessian)
        .solve(par);
}

#endif  // __NIMBLE_NIMOPTIM_H
