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
#include <nimble/nimbleCppAD.h>
#include <algorithm>
#include <string>

// ---------------------------------------------------------------------------
// This class and wrapper are equivalent to std::bind(-,-) in C++11, so that:
//   object->method(par) == NimBoundMethod<T>(&T::method, object)(par);
//                       == NimBind(&T::method, object)(par);
//                       == std::bind(&T::method, object)(par);
template <class RT, class T>
class NimBoundMethod {
   public:
    NimBoundMethod(RT (T::*method)(NimArr<1, double>&), T* object)
        : method_(method), object_(object) {}
    RT operator()(NimArr<1, double>& par) {
        return (object_->*method_)(par);
    }

   private:
    RT (T::*method_)(NimArr<1, double>&);
    T* object_;
};

template <class RT, class T>
inline NimBoundMethod<RT, T> NimBind(RT (T::*method)(NimArr<1, double>&),
                                 T* object) {
  return NimBoundMethod<RT, T>(method, object);
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
    void calc_hessian(NimArr<1, double> par,
		      NimArr<2, double> &hessian);
   protected:
    // These are callbacks used internally by fn() and gr().
    virtual double function() = 0;
    virtual void gradient() = 0;

    // Problem parameters.
    const std::string& method_;
    NimArr<1, double>& lower_;
    NimArr<1, double>& upper_;
    nimSmartPtr<OptimControlNimbleList> control_;
    bool hessian_;

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

class NimOptimProblem_model : public NimOptimProblem {
 private:
  int length_wrt;
  bool use_gr;
  NodeVectorClassNew_derivs &nodes;
  NimArr<1, double> lastP;
  double lastValue;
  NimArr<1, double> lastGradient;
  NimArr<1, double> derivOrders;
  NimArr<1, double> get_wrt() {
    NimArr<1, double > wrt;
    wrt.setSize(length_wrt);
    getValues(wrt, nodes.get_model_wrt_accessor());
    return wrt;
  }

 public:
  NimOptimProblem_model(
			NodeVectorClassNew_derivs &nodes_,
			bool use_gr_,
			const std::string& method,
			NimArr<1, double>& lower, NimArr<1, double>& upper,
			nimSmartPtr<OptimControlNimbleList> control,
			bool hessian
			)
    : NimOptimProblem(method, lower, upper, control, hessian),
    use_gr(use_gr_),
    nodes(nodes_) {
      derivOrders.setSize(2);
      derivOrders[0] = 0;
      derivOrders[1] = 1;
      length_wrt = nodes.get_model_wrt_accessor().getTotalLength();
      lastP.setSize(length_wrt, false, false);
      lastGradient.setSize(length_wrt, false, false);
    }

  nimSmartPtr<OptimResultNimbleList> solve() {
    NimArr<1, double> initP = get_wrt();
    bool hessian = NimOptimProblem::hessian_;
    NimOptimProblem::hessian_ = false;
    nimSmartPtr<OptimResultNimbleList> result = NimOptimProblem::solve(initP);
    NimOptimProblem::hessian_ = hessian;
    if(hessian) {
      /* TO-DO: We should be able to re-use last-run tape calculation */
      /* and thereby not need to re-do 0-order forward (value) calculation, */
      /* but for now we do re-do it.*/
      derivOrders[1] = 2;
      setValues(par_, nodes.get_model_wrt_accessor());
      nimSmartPtr<NIMBLE_ADCLASS> ADresult = nimDerivs_calculate(nodes, derivOrders);
      NimArr<2, double> hessianMap;
      int n = initP.size();
      // offset = 0; stride1 = 1; stride2 = n; size1 = n; size2 = n;
      hessianMap.setMap(ADresult->hessian, 0, 1, n, n, n);
      result->hessian = hessianMap;
      derivOrders[1] = 1; // will probably never be used again, but just in case, reset this.
    }
    return result;
  }
  
 private:
  void run_calculations() {
    setValues(par_, nodes.get_model_wrt_accessor());
    nimSmartPtr<NIMBLE_ADCLASS> ADresult = nimDerivs_calculate(nodes, derivOrders);
    if(par_.size() != length_wrt) {
      std::cout<<"Error in C++: wrong length par_ in nimOptim_model"<<std::endl;
    }
    std::memcpy(lastP.getPtr(), par_.getPtr(), length_wrt * sizeof(double));
    lastValue = ADresult->value[0];
    if(ADresult->jacobian.dim()[0] != 1) {
      std::cout<<"Error in C++: jacobian number of rows != 1 in nimOptim_model"<<std::endl;
    } 
    std::memcpy(lastGradient.getPtr(), ADresult->jacobian.getPtr(), length_wrt * sizeof(double));
  }

    
 protected:
  virtual double function() {
    run_calculations();
    return lastValue;
  };
  virtual void gradient() {
    bool sameAsLastCall = true;
    if(par_.size() != length_wrt) {
      std::cout<<"Error in C++: wrong length par_ in nimOptim_model"<<std::endl;
    }
    for(size_t i = 0; i < length_wrt; ++i) {
      if(lastP[i] != par_[i]) {
	sameAsLastCall = false;
	break;
      }
    }
    if(!sameAsLastCall) run_calculations();
    if(ans_.size() != length_wrt) {
      std::cout<<"Error in C++: wrong length ans_ in nimOptim_model"<<std::endl;
    }
    std::memcpy(ans_.getPtr(), lastGradient.getPtr(), length_wrt * sizeof(double));
  };
};

// ---------------------------------------------------------------------------
// These nimOptim() and nimOptimDefaultControl() will appear in generated code

nimSmartPtr<OptimControlNimbleList> nimOptimDefaultControl();

inline nimSmartPtr<OptimResultNimbleList> nimOptim_model(
							 NodeVectorClassNew_derivs &nodes,
							 bool use_gr,
							 const std::string& method,
							 double lower,
							 double upper,
							 nimSmartPtr<OptimControlNimbleList> control,
							 bool hessian = false
							 ) {
  NimArr<1, double> lower_vector;
  lower_vector.initialize(lower, true, 1);
  NimArr<1, double> upper_vector;
  upper_vector.initialize(upper, true, 1);
  control->fnscale = -1.;
  return NimOptimProblem_model(nodes, use_gr, method, lower_vector, upper_vector, control, hessian ).solve( );
}

inline nimSmartPtr<OptimResultNimbleList> nimOptim_model(
							 NodeVectorClassNew_derivs &nodes,
							 bool use_gr,
							 const std::string& method,
							 NimArr<1, double> &lower_vector,
							 NimArr<1, double> &upper_vector,
							 nimSmartPtr<OptimControlNimbleList> control,
							 bool hessian = false
							 ) {
  control->fnscale = -1.;
  return NimOptimProblem_model(nodes, use_gr, method, lower_vector, upper_vector, control, hessian ).solve( );
}


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
