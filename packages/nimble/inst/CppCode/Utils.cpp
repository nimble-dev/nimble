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

#include<limits>
#include "nimble/Utils.h"

// A utility function that will return floor(x) unless x is within numerical imprecision of an integer, in which case it will return round(x)
int floorOrEquivalent(double x) {
  double roundX = round(x);
  double sqrtEpsilon = sqrt(std::numeric_limits<double>::epsilon());
  bool shouldBeExactInteger(false);
  if(fabs(x) > sqrtEpsilon) { // This algorithm for numerical equivalence imitates what happens in all.equal.numeric in R
    if(fabs(x - roundX) / fabs(x) < sqrtEpsilon)
      shouldBeExactInteger = true;
  } else {
    if(fabs(x - roundX) < sqrtEpsilon) // in the present context, this would mean length should be zero anyway
      shouldBeExactInteger = true;
  }
  if(shouldBeExactInteger)
    return(static_cast<int>(roundX));
  return(floor(x));
}

int rFunLength(int Arg) {
  return Arg;
}

int rFunLength(double Arg) {
  return Arg;
}

int rFunLength(bool Arg) {
  return Arg;
}

bool decide(double lMHr) { // simple function accept or reject based on log Metropolis-Hastings ratio
  if(ISNAN(lMHr)) {
    // double junk = runif(0,1); // purely for debugging to ensure
    // same sequence of RNG uses for comparison purposes.
    return(false);
  }
  if(lMHr > 0) return(true);
  if(runif(0,1) < exp(lMHr)) return(true);
  return(false);
}

void nimStop(string msg) {NIMERROR(msg.c_str());}
void nimStop() {NIMERROR("");}

bool nimNot(bool x) {return(!x);}

double pow_int(double a, double b) {
  return pow(a, fround(b, 0.));
}

double ilogit(double x) {return(1./(1. + exp(-x)));}

double icloglog(double x) {return(1.-exp(-exp(x)));}

double iprobit(double x) {return(pnorm(x, 0., 1., 1, 0));}

double probit(double x) {return(qnorm(x, 0., 1., 1, 0));}
//double abs(double x) {return(fabs(x));}
double cloglog(double x) {return(log(-log(1.-x)));}
int nimEquals(double x1, double x2) {return(x1 == x2 ? 1 : 0);}
double nimbleIfElse(bool condition, double x1, double x2) {return(condition ? x1 : x2);}
double lfactorial(double x) {return(lgammafn(x + 1));} //{return(lgamma1p(x));} not numerically equivalent to lgamma(x+1), which is what is done from R. lgamma seems to clean up numerical zeros
double factorial(double x) {return(gammafn(1+x));}
//double loggam(double x) {return(lgammafn(x));}
double logit(double x) {return(log(x / (1.-x)));}
double nimRound(double x) {return(fround(x, 0.));}
double pairmax(double x1, double x2) {return(x1 > x2 ? x1 : x2);}
double pairmin(double x1, double x2) {return(x1 < x2 ? x1 : x2);}
//double phi(double x) {return(pnorm(x, 0., 1., 1, 0));}
int nimStep(double x) { return(x >= 0 ? 1 : 0);}
double cube(double x) {return(x*x*x);}
double inprod(double v1, double v2) {return(v1*v2);}
