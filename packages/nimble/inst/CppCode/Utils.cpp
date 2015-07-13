#include "nimble/Utils.h"
#define PRINTF Rprintf
#define NIMERROR error
#define RBREAK(msg) {PRINTF(msg); return(R_NilValue);}

bool decide(double lMHr) { // simple function accept or reject based on log Metropolis-Hastings ratio
  if(isnan(lMHr)) return(false);
  if(lMHr > 0) return(true);
  if(runif(0,1) < exp(lMHr)) return(true);
  return(false);
}

void nimStop(string msg) {NIMERROR(msg.c_str());}
void nimStop() {NIMERROR("");}

double ilogit(double x) {return(1./(1. + exp(-x)));}

double icloglog(double x) {return(1.-exp(-exp(x)));}

double iprobit(double x) {return(pnorm(x, 0., 1., 1, 0));}

double probit(double x) {return(qnorm(x, 0., 1., 1, 0));}
//double abs(double x) {return(fabs(x));}
double cloglog(double x) {return(log(-log(1.-x)));}
int nimEquals(double x1, double x2) {return(x1 == x2 ? 1 : 0);}
double nimbleIfElse(bool condition, double x1, double x2) {return(condition ? x1 : x2);}
double lfactorial(double x) {return(lgamma1p(x));}
double factorial(double x) {return(gammafn(1+x));}
//double loggam(double x) {return(lgammafn(x));}
double logit(double x) {return(log(x / (1.-x)));}
double nimbleRound(double x) {return(fround(x, 0.));}
double pairmax(double x1, double x2) {return(x1 > x2 ? x1 : x2);}
double pairmin(double x1, double x2) {return(x1 < x2 ? x1 : x2);}
//double phi(double x) {return(pnorm(x, 0., 1., 1, 0));}
int nimStep(double x) { return(x >= 0 ? 1 : 0);}
double cube(double x) {return(x*x*x);}
double inprod(double v1, double v2) {return(v1*v2);}
