#ifndef __NIMBLE_NIMOPTIM_H
#define __NIMBLE_NIMOPTIM_H

#include <nimble/NimArr.h>
#include <nimble/optimTypes.h>
#include <nimble/smartPtrs.h>
#include <algorithm>

// These match the interface to R's optim().
typedef double optimfn(int, double *, void *);
typedef void optimgr(int, double *, double *, void *);

nimSmartPtr<OptimResultNimbleList> nimOptim_internal(NimArr<1, double> &par,
                                                     optimfn fn, void *ex,
                                                     const char *method);

// ---------------------------------------------------------------------------
// Interface for RCfunction functions.

typedef double NimObjectiveFn(NimArr<1, double> &par);
typedef NimArr<1, double> NimObjectiveGr(NimArr<1, double> &par);

// This behaves like R's native optim() function on RCfunction functions.
nimSmartPtr<OptimResultNimbleList> nimOptim(NimArr<1, double> &par,
                                            NimObjectiveFn fn, void *gr,
                                            const char *method);

// ---------------------------------------------------------------------------
// Interface for nimbleFunction classes.

// This class and wrapper are equivalent to std::bind(-,-) in C++11, so that:
//   object->method(par) == NimBoundMethod<T>(&T::method, object)(par);
//                       == NimBind(&T::method, object)(par);
//                       == std::bind(&T::method, object)(par);
template <class T>
class NimBoundMethod {
   public:
    NimBoundMethod(double (T::*method)(NimArr<1, double> &), T *object)
        : method_(method), object_(object) {}
    double operator()(NimArr<1, double> &par) {
        return (object_->*method_)(par);
    }

   private:
    double (T::*method_)(NimArr<1, double> &);
    T *object_;
};

template <class T>
inline NimBoundMethod<T> NimBind(double (T::*method)(NimArr<1, double> &),
                                 T *object) {
    return NimBoundMethod<T>(method, object);
}

template <class T>
double optimfn_method_wrapper(int n, double *par, void *fn) {
    NimArr<1, double> nim_par;
    nim_par.setSize(n, false, false);
    std::copy(par, par + n, nim_par.getPtr());
    return static_cast<NimBoundMethod<T> *>(fn)->operator()(nim_par);
}

// This behaves like R's native optim() function on nimbleFunction classes.
template <class T>
nimSmartPtr<OptimResultNimbleList> nimOptim(NimArr<1, double> &par,
                                            NimBoundMethod<T> fn, void *gr,
                                            const char *method) {
    return nimOptim_internal(par, &optimfn_method_wrapper<T>, &fn, method);
}

#endif  // __NIMBLE_NIMOPTIM_H
