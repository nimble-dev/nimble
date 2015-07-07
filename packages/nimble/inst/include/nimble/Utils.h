#ifndef __UTILS
#define __UTILS

//#include <vector>
#include "R.h"
#include<Rinternals.h>
#include<Rmath.h>

//using namespace std;

#define PRINTF Rprintf
#define NIMERROR error
#define RBREAK(msg) {PRINTF(msg); return(R_NilValue);}

// code copied from nmath.h - useful utilities 
# define MATHLIB_ERROR(fmt,x)		error(fmt,x);
# define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)

#define ML_POSINF	R_PosInf
#define ML_NEGINF	R_NegInf
#define ML_NAN		R_NaN

#define _(String) String

#define ML_VALID(x)	(!ISNAN(x))

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }

#define ML_ERROR(x, s) { \
   if(x > ME_DOMAIN) { \
       const char *msg = ""; \
       switch(x) { \
       case ME_DOMAIN: \
	   msg = _("argument out of domain in '%s'\n");	\
	   break; \
       case ME_RANGE: \
	   msg = _("value out of range in '%s'\n");	\
	   break; \
       case ME_NOCONV: \
	   msg = _("convergence failed in '%s'\n");	\
	   break; \
       case ME_PRECISION: \
	   msg = _("full precision may not have been achieved in '%s'\n"); \
	   break; \
       case ME_UNDERFLOW: \
	   msg = _("underflow occurred in '%s'\n");	\
	   break; \
       } \
       MATHLIB_WARNING(msg, s); \
   } \
}


// code copied from dpq.h - useful utilities for return values of dist functions
                                                        /* "DEFAULT" */
                                                        /* --------- */
#define R_D__0  (log_p ? ML_NEGINF : 0.)                /* 0 */
#define R_D__1  (log_p ? 0. : 1.)                       /* 1 */
#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)          /* 0 */
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0)          /* 1 */
#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */

#define give_log log_p

#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x) 	  (fabs((x) - floor((x)+0.5)) > 1e-7)
/* [neg]ative or [non int]eger : */
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))

#define R_D_nonint_check(x) 				\
   if(R_D_nonint(x)) {					\
	MATHLIB_WARNING("non-integer x = %f", x);	\
	return R_D__0;					\
   }

bool decide(double lMHr);
//void allocate(vector< vector <double> > *vv, int sampleSize, int variableSize);

// needed for link functions
double ilogit(double x);
double icloglog(double x);
double iprobit(double x);
double probit(double x);
//double abs(double x);
double cloglog(double x);
int nimbleEquals(double x1, double x2);
double nimbleIfElse(bool condition, double x1, double x2);
double lfactorial(double x);
double factorial(double x);
//double loggam(double x);
double logit(double x);
double nimbleRound(double x);
double pairmax(double x1, double x2);
double pairmin(double x1, double x2);
//double phi(double x);
int nimbleStep(double x); 
double cube(double x);
double inprod(double v1, double v2);

#endif
