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
// #define _TIME_AD
// To see all timing components, use:
//.Call(getNativeSymbolInfo("report_AD_timers", compiled_model$dll))
//.Call(getNativeSymbolInfo("report_AD_timers", compiled_nf_using_derivs$dll))
// Calling both is necessary because NIMBLE will create one DLL for the model
// and one for nimbleFunctions using the model, and the two DLLs won't
// see the same timer objects.  The zeros in each indicate which timers
// were not used by stuff in that DLL.

#ifndef _NIMBLE_CPPAD
#define _NIMBLE_CPPAD

/* Definitions only to be included when a nimbleFunction needs CppAD */
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <nimble/EigenTypedefs.h>
#include <nimble/accessorClasses.h>
#include <nimble/nodeFun.h>
#include <nimble/predefinedNimbleLists.h>
#include <cstdio>
#include <vector>

#ifdef _TIME_AD
extern "C" {
  SEXP reset_AD_timers(SEXP SreportInterval);
  SEXP report_AD_timers();
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
    std::cout << "setting up "<<name<<std::endl;
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
    if(verbose) {
      std::cout<<name<<" increment "<<static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count())<<std::endl;
    }
    totaltime += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count());
  }
  void report() {
    tick();
    if(ticks >= report_interval) {
      totalticks += ticks;
      show_report();
      ticks = 0;
    }
  }
  void show_report() {
      printf("Reporting time for %s (%i): %g (%i)\n",
	     name.c_str(),
	     totalticks,
	     totaltime,
	     static_cast<int>(touched));
  }
};

void derivs_getDerivs_timer_start();
void derivs_getDerivs_timer_stop();
void derivs_run_tape_timer_start();
void derivs_run_tape_timer_stop();
void derivs_tick_id();
void show_tick_id();
#endif

/* nimbleCppADinfoClass is the class to convey information from a nimbleFunction
   object
   to generic CppAD driver wrappers like calcjacobian.
   Each nimbleFunction enabled for CppAD will have an object of this class. */

class nimbleCppADrecordingInfoClass {
 private:
  bool recording_;
  CppAD::tape_id_t tape_id_;
  CppAD::local::ADTape<double>* tape_handle_;
 public:
  bool& recording() {return recording_;}
  CppAD::tape_id_t& tape_id() {return tape_id_;}
  CppAD::local::ADTape<double>* &tape_handle() {return tape_handle_;}
 nimbleCppADrecordingInfoClass(bool r_, CppAD::tape_id_t tid_, CppAD::local::ADTape<double>* th_) :
  recording_(r_),
    tape_id_(tid_),
    tape_handle_(th_) {}
 nimbleCppADrecordingInfoClass(CppAD::tape_id_t tid_, CppAD::local::ADTape<double>* th_) :
  recording_(false),
    tape_id_(tid_),
    tape_handle_(th_) {}
 nimbleCppADrecordingInfoClass() : recording_(false) {}
};

class nimbleCppADinfoClass {
 public:
  std::vector<double> independentVars;
  std::vector<double> dynamicVars;
  CppAD::ADFun<double> *ADtape;
};

/* nimbleFunctionCppADbase is a base class to be inherited by all
   CppAD-enabled nimbleFunctions. Some of these functions might
   make more sense as stand-alone functions.  Let's see. */
class nimbleFunctionCppADbase {
public:
  void                        getDerivs(nimbleCppADinfoClass &ADinfo,
                                        const NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList);
  
   nimSmartPtr<NIMBLE_ADCLASS> getDerivs_wrapper(nimbleCppADinfoClass &ADinfo,
                                        const NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector){
            nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
            getDerivs(ADinfo, derivOrders, wrtVector, ansList);
            return(ansList);
    }
};

nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
    NodeVectorClassNew_derivs &nodes, const NimArr<1, double> &derivOrders);
nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
    NodeVectorClassNew_derivs &nodes, const double derivOrders);
nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(
    NodeVectorClassNew_derivs &nodes, int iNodeFunction,
    NimArr<1, double> &derivOrders);

NimArr<1, double> make_vector_if_necessary(int);
NimArr<1, double> make_vector_if_necessary(double);
NimArr<1, double> make_vector_if_necessary(NimArr<1, double>);

#endif
