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

// define this to include timing code
#define _TIME_AD

/* Definitions only to be included when a nimbleFunction needs CppAD */
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
//#include <TMB/tmbDists.h>
#include <nimble/EigenTypedefs.h>
#include <nimble/accessorClasses.h>
#include <nimble/nodeFun.h>
//#include <nimble/RcppNimbleUtils.h>
#include <nimble/predefinedNimbleLists.h>
//#include <nimble/RcppNimbleUtils.h>
//#include <nimble/NimArr.h>
//#include <nimble/smartPtrs.h>
#include <cstdio>
#include <vector>

#ifdef _TIME_AD
extern "C" {
  SEXP reset_AD_timers(SEXP SreportInterval);
}
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>
typedef std::chrono::high_resolution_clock::time_point timetype;
class ad_timer {
private:
  timetype t1, t2;
  double totaltime;
  unsigned int ticks, totalticks;
  unsigned int report_interval;
  bool touched;
public:
  std::string name;
  ad_timer(std::string myname) {
    touched = false;
    name = myname;
    reset();
    set_interval(100);
    std::cout << std::chrono::high_resolution_clock::duration::period::den << std::endl;
  }
  void reset() {
    totaltime = 0;
    ticks = 0;
    totalticks = 0;
  }
  void set_interval(int ri) {
    report_interval = static_cast<unsigned int>(ri);
  }
  void start(bool verbose = false) {
    touched = true;
      t1 = std::chrono::high_resolution_clock::now();
    if(verbose) {
      printf("start %g\n", static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1.time_since_epoch()).count()));
    }
  }
  void tick() {
    ++ticks;
  }
  void stop(bool verbose = false) {
    t2 = std::chrono::high_resolution_clock::now();
    double oldtotaltime = totaltime;
    if(verbose) {
      printf("stop %g (%g) (%g)",
	     static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t2.time_since_epoch()).count()),
	     static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()),
	     totaltime);
    }
    totaltime += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count());
    if(verbose) {
      printf(" (%g) (%g)\n", totaltime, oldtotaltime-totaltime);
    }
  }
  void report() {
    tick();
    if(ticks >= report_interval) {
      totalticks += ticks;
      printf("Reporting time for %s (%i): %g (%i)\n",
	     name.c_str(),
	     totalticks,
	     totaltime,
	     static_cast<int>(touched));
      ticks = 0;
    }
  }
};
void derivs_getDerivs_timer_start();
void derivs_getDerivs_timer_stop();
void derivs_run_tape_timer_start();
void derivs_run_tape_timer_stop();
#endif

/* nimbleCppADinfoClass is the class to convey information from a nimbleFunction
   object
   to generic CppAD driver wrappers like calcjacobian.
   Each nimbleFunction enabled for CppAD will have an object of this class. */

class nimbleCppADinfoClass {
 public:
  std::vector<double> independentVars;
  CppAD::ADFun<double> *ADtape;
};

/* nimbleFunctionCppADbase is a base class to be inherited by all
   CppAD-enabled nimbleFunctions. Some of these functions might
   make more sense as stand-alone functions.  Let's see. */
class nimbleFunctionCppADbase {
public:
  void                        getDerivs(nimbleCppADinfoClass &ADinfo,
                                        NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList);
  
   nimSmartPtr<NIMBLE_ADCLASS> getDerivs_wrapper(nimbleCppADinfoClass &ADinfo,
                                        NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector){
            nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
            getDerivs(ADinfo, derivOrders, wrtVector, ansList);
            return(ansList);
    }
};

nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    const NodeVectorClassNew_derivs &nodes, const NimArr<1, double> &derivOrders);
nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    const NodeVectorClassNew_derivs &nodes, const double derivOrders);
nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    NodeVectorClassNew_derivs &nodes, int iNodeFunction,
    NimArr<1, double> &derivOrders);
