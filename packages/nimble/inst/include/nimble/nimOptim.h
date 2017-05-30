#ifndef __NIMBLE_NIMOPTIM_H
#define __NIMBLE_NIMOPTIM_H

#include <nimble/NimArr.h>
#include <nimble/predefinedNimbleLists.h>
#include <nimble/smartPtrs.h>
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
    nimSmartPtr<OptimResultNimbleList> solve(
        NimArr<1, double>& par, const std::string& method,
        NimArr<1, double>& lower, NimArr<1, double>& upper,
        nimSmartPtr<OptimControlNimbleList> control, bool hessian);

   private:
    // These are callbacks for R's optim() where this is passed in as the final
    // argument `void * ex`.
    static double fn(int, double*, void*);
    static void gr(int, double*, double*, void*);

   protected:
    // These are callbacks used internally by fn() and gr().
    virtual double function() = 0;
    virtual void gradient() { NIMERROR("Gradient is not defined"); }

    // These are used as temporaries for C <-> NimArr conversion.
    NimArr<1, double> par_;
    NimArr<1, double> ans_;
};

template <class Fn>
class NimOptimProblem_Fun : public NimOptimProblem {
   public:
    NimOptimProblem_Fun(Fn fn) : fn_(fn) {}

   protected:
    virtual double function() { return fn_(par_); }

   private:
    Fn fn_;
};

template <class Fn, class Gr>
class NimOptimProblem_Fun_Grad : public NimOptimProblem {
   public:
    NimOptimProblem_Fun_Grad(Fn fn, Gr gr) : fn_(fn), gr_(gr) {}

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
    return NimOptimProblem_Fun_Grad<Fn, Gr>(fn, gr).solve(
        par, method, lower, upper, control, hessian);
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
    return NimOptimProblem_Fun<Fn>(fn).solve(par, method, lower, upper, control,
                                             hessian);
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
    return NimOptimProblem_Fun_Grad<Fn, Gr>(fn, gr).solve(
        par, method, lower_vector, upper_vector, control, hessian);
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
    return NimOptimProblem_Fun<Fn>(fn).solve(par, method, lower_vector,
                                             upper_vector, control, hessian);
}

#endif  // __NIMBLE_NIMOPTIM_H
