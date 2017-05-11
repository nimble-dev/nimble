#ifndef __DISTS
#define __DISTS

#include "Utils.h"

bool R_IsNA(double*, int);
bool R_isnancpp(double*, int);
bool R_FINITE_VEC(double*, int);

extern "C" {
// BLAS/LAPACK routines
  // not needed when #include "R_ext/lapack.h"
  /*
  int dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*); 
  int dtrmv_(char*, char*, char*, int*, double*, int*, double*, int*);
  int dtrsv_(char*, char*, char*, int*, double*, int*, double*, int*);
  int dtrsm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
  int dtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
  int dpotrf_(char*, int*, double*, int*, int*);
  */

  // NIMBLE C wrappers called from R
  SEXP C_dmnorm_chol(SEXP, SEXP, SEXP, SEXP, SEXP); 
  SEXP C_rmnorm_chol(SEXP, SEXP, SEXP);
  SEXP C_dmvt_chol(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
  SEXP C_rmvt_chol(SEXP, SEXP, SEXP, SEXP); 
  SEXP C_dwish_chol(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_rwish_chol(SEXP, SEXP, SEXP);
  SEXP C_dcat(SEXP, SEXP, SEXP);
  SEXP C_rcat(SEXP, SEXP);
  SEXP C_dmulti(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rmulti(SEXP, SEXP);
  SEXP C_ddirch(SEXP, SEXP, SEXP);
  SEXP C_rdirch(SEXP);
  SEXP C_dinterval(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rinterval(SEXP, SEXP, SEXP);
  SEXP C_dexp_nimble(SEXP, SEXP, SEXP);
  SEXP C_rexp_nimble(SEXP, SEXP);
  SEXP C_pexp_nimble(SEXP, SEXP, SEXP, SEXP);
  SEXP C_qexp_nimble(SEXP, SEXP, SEXP, SEXP);
  SEXP C_dinvgamma(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rinvgamma(SEXP, SEXP, SEXP);
  SEXP C_pinvgamma(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_qinvgamma(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_dsqrtinvgamma(SEXP, SEXP, SEXP, SEXP);
  SEXP C_rsqrtinvgamma(SEXP, SEXP, SEXP);

}

// NOTE: R CMD SHLIB seems to handle C++ code without using wrapping the functions in 'extern "C"'; note that some of these functions have a bit of C++ syntax

// core scalar d/r functions provided by NIMBLE to extend R
double dcat(double, double*, int, int);
double rcat(double*, int);
double dmulti(double*, double, double*, int, int);
void rmulti(int*, double, double*, int);
double ddirch(double*, double*, int, int);
void rdirch(double*, double*, int);

double dmnorm_chol(double*, double*, double*, int, double, int, int);
void rmnorm_chol(double *, double*, double*, int, double);
double dmvt_chol(double*, double*, double*, double, int, double, int, int);
void rmvt_chol(double *, double*, double*, double, int, double);
double dwish_chol(double*, double*, double, int, double, int, int);
void rwish_chol(double*, double*, double, int, double, int);

double dinterval(double, double, double*, int, int);
double rinterval(double, double*, int);

// SHOULD BE IN nimDists.h 
// Chris comment on above line: I don't think this is the case...

extern "C" {
  SEXP C_rt_nonstandard(SEXP, SEXP, SEXP, SEXP);
  SEXP C_dt_nonstandard(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_pt_nonstandard(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP C_qt_nonstandard(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
}

double dt_nonstandard(double, double, double, double, int);
double rt_nonstandard(double, double, double);
double pt_nonstandard(double, double, double, double, int, int);
double qt_nonstandard(double, double, double, double, int, int);

double dconstraint(double, double, int);
double rconstraint(double);

double rexp_nimble(double);
double dexp_nimble(double, double, int);
double pexp_nimble(double, double, int, int);
double qexp_nimble(double, double, int, int);

double rinvgamma(double, double);
double dinvgamma(double, double, double, int);
double pinvgamma(double, double, double, int, int);
double qinvgamma(double, double, double, int, int);

double rsqrtinvgamma(double, double);
double dsqrtinvgamma(double, double, double, int);

double rflat();
double dflat(double, int);
double rhalfflat();
double dhalfflat(double, int);
#endif
