#ifndef __UTILS
#define __UTILS

#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>
#include <limits>
#include <string>
#include "R.h"
using std::string;

// A utility function that will return floor(x) unless x is within numerical
// imprecision of an integer, in which case it will return round(x)
int floorOrEquivalent(double x);

int rFunLength(int Arg);
int rFunLength(double Arg);
int rFunLength(bool Arg);

class nimbleTimerClass_ {
 public:
  clock_t t_start;
  void startNimbleTimer() { t_start = clock(); }
  double endNimbleTimer() {
    return (((double)(clock() - t_start)) / CLOCKS_PER_SEC);
  }
};

#if defined(__GNUG__) || defined(__clang__)
#define NIM_LIKELY(x) __builtin_expect(!!(x), 1)
#define NIM_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define NIM_LIKELY(x) (x)
#define NIM_UNLIKELY(x) (x)
#endif

#define PRINTF Rprintf
#define NIMERROR error
#define RBREAK(msg)      \
  {                      \
    PRINTF(msg);         \
    return (R_NilValue); \
  }

#define NIM_ASSERT1(cond, msg)   \
  {                              \
    if (NIM_UNLIKELY(!(cond))) { \
      NIMERROR("Error: ", msg);  \
    }                            \
  }
#define NIM_ASSERT2(cond, msg, msgArg1)  \
  {                                      \
    if (NIM_UNLIKELY(!(cond))) {         \
      NIMERROR("Error: ", msg, msgArg1); \
    }                                    \
  }
#define NIM_ASSERT3(cond, msg, msgArg1, msgArg2)  \
  {                                               \
    if (NIM_UNLIKELY(!(cond))) {                  \
      NIMERROR("Error: ", msg, msgArg1, msgArg2); \
    }                                             \
  }
#define NIM_ASSERT4(cond, msg, msgArg1, msgArg2, msgArg3)  \
  {                                                        \
    if (NIM_UNLIKELY(!(cond))) {                           \
      NIMERROR("Error: ", msg, msgArg1, msgArg2, msgArg3); \
    }                                                      \
  }

#define NIM_ASSERT_SIZE(my_array, n)                                  \
  NIM_ASSERT3(my_array.dimSize(0) == n,                               \
              #my_array " has wrong size: expected %d, actual %d", n, \
              my_array.dimSize(0));

// Code copied from nmath.h - useful utilities.
#define MATHLIB_ERROR(fmt, x) error(fmt, x);
#define MATHLIB_WARNING(fmt, x) warning(fmt, x)
#define MATHLIB_WARNING2(fmt, x, x2) warning(fmt, x, x2)
#define MATHLIB_WARNING3(fmt, x, x2, x3) warning(fmt, x, x2, x3)
#define MATHLIB_WARNING4(fmt, x, x2, x3, x4) warning(fmt, x, x2, x3, x4)

#define ML_POSINF R_PosInf
#define ML_NEGINF R_NegInf
#define ML_NAN R_NaN

#define _(String) String

#define ML_VALID(x) (!ISNAN(x))

#define ME_NONE 0
/*	no error */
#define ME_DOMAIN 1
/*	argument out of domain */
#define ME_RANGE 2
/*	value out of range */
#define ME_NOCONV 4
/*	process did not converge */
#define ME_PRECISION 8
/*	does not have "full" precision */
#define ME_UNDERFLOW 16
/*	and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN    \
  {                          \
    ML_ERROR(ME_DOMAIN, ""); \
    return ML_NAN;           \
  }

#define ML_ERROR(x, s)                                                    \
  {                                                                       \
    if (x > ME_DOMAIN) {                                                  \
      const char *msg = "";                                               \
      switch (x) {                                                        \
        case ME_DOMAIN:                                                   \
          msg = _("argument out of domain in '%s'\n");                    \
          break;                                                          \
        case ME_RANGE:                                                    \
          msg = _("value out of range in '%s'\n");                        \
          break;                                                          \
        case ME_NOCONV:                                                   \
          msg = _("convergence failed in '%s'\n");                        \
          break;                                                          \
        case ME_PRECISION:                                                \
          msg = _("full precision may not have been achieved in '%s'\n"); \
          break;                                                          \
        case ME_UNDERFLOW:                                                \
          msg = _("underflow occurred in '%s'\n");                        \
          break;                                                          \
      }                                                                   \
      MATHLIB_WARNING(msg, s);                                            \
    }                                                                     \
  }

// Code copied from dpq.h - useful utilities for return values of dist
// functions.
/* "DEFAULT" */
/* --------- */
#define R_D__0 (log_p ? ML_NEGINF : 0.)       /* 0 */
#define R_D__1 (log_p ? 0. : 1.)              /* 1 */
#define R_DT_0 (lower_tail ? R_D__0 : R_D__1) /* 0 */
#define R_DT_1 (lower_tail ? R_D__1 : R_D__0) /* 1 */
#define R_D_val(x) (log_p ? log(x) : (x))     /*  x  in pF(x,..) */

#define give_log log_p

#define R_D_forceint(x) floor((x) + 0.5)
#define R_D_nonint(x) (fabs((x)-floor((x) + 0.5)) > 1e-7)
/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))

#define R_D_nonint_check(x)                   \
  if (R_D_nonint(x)) {                        \
    MATHLIB_WARNING("non-integer x = %f", x); \
    return R_D__0;                            \
  }

#define EIGEN_CHOL(x) (x).selfadjointView<Eigen::Upper>().llt().matrixU()

bool decide(double lMHr);

void nimStop(string msg);
void nimStop();

bool nimNot(bool x);

// needed for link functions
double ilogit(double x);
double icloglog(double x);
double iprobit(double x);
double probit(double x);
double cloglog(double x);
int nimEquals(double x1, double x2);
double nimbleIfElse(bool condition, double x1, double x2);
double lfactorial(double x);
double factorial(double x);
double logit(double x);
double nimRound(double x);
double pairmax(double x1, double x2);
double pairmin(double x1, double x2);
int nimStep(double x);
double cube(double x);
double inprod(double v1, double v2);

inline double nimble_NaN() {
  return std::numeric_limits<double>::has_quiet_NaN
             ? std::numeric_limits<double>::quiet_NaN()
             : (0. / 0.);
}
#endif
