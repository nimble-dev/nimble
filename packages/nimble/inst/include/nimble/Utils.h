#ifndef __UTILS
#define __UTILS

//#include <vector>
#include "R.h"
#include<Rinternals.h>
#include<Rmath.h>

//using namespace std;

#define PRINTF Rprintf
#define RBREAK(msg) {PRINTF(msg); return(R_NilValue);}


// code copied from dpq.h - useful utilities for return values of dist functions
                                                        /* "DEFAULT" */
                                                        /* --------- */
#define R_D__0  (log_p ? ML_NEGINF : 0.)                /* 0 */
#define R_D__1  (log_p ? 0. : 1.)                       /* 1 */
#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)          /* 0 */
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0)          /* 1 */


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
