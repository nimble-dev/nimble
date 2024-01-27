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
#include <algorithm>

/*
  nimble_atomic_base has two purposes:
1. It allows a base calss pointer to any nimble atomic
2. Its constructor records the atomic vec manager (a static in CppAD), and
   its destructor resets that to its original location.
   This prevents crashes when the atomic_three<> destructor is called.

   ***It is necessary to declare atomic_three<> inheritance before nimble_atomic_base
   inheritance to guarantee the correct order of destructor calls.***
 */
class nimble_atomic_base {
 public:
  nimble_atomic_base();
  virtual ~nimble_atomic_base();
  std::vector<CppAD::local::atomic_index_info>* vec_ptr_where_constructed;
  /* This needs to be virtual to avoid any compiler inlining it and then
     crossing boundaries of the static variables used in CppAD.
     By making it virtual, the impementation will be looked up at
     run-time and will be in the correct DLL. */
  virtual void set_CppAD_atomic_info_vec_manager( std::vector<CppAD::local::atomic_index_info>* vec_ptr );
};

template<class T>
class unary_atomic_class : public CppAD::atomic_three<T>, public nimble_atomic_base {
  // This layer in the class hierarchy simply provides the same for_type and rev_depend
  // for all atomic classes representing unary functions below.
 public:
 unary_atomic_class(const std::string& name) : CppAD::atomic_three<T>(name) {};
  virtual ~unary_atomic_class() {};
 private:
   // for_type is essential.
  // This determines which elements of y are constants, dynamic parameters, or variables
  //   depending on types of x.
  // If omitted, calls to forward (and possibly reverse) will not happen.
  virtual bool for_type(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
  {
    //    std::cout<<"in for_type\n";
    type_y[0] = type_x[0];
    return true;
  }
  // Not sure this is ever needed.
  virtual bool for_type(
			const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
  {
    //    std::cout<<"in meta for_type\n";
    type_y[0] = type_x[0];
    return true;
  }
  // rev_depend is used when the optimize() method is called for a tape (an ADFun).
  virtual bool rev_depend(
			  const CppAD::vector<double>&          parameter_x ,
			  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
			  CppAD::vector<bool>&                depend_x    ,
			  const CppAD::vector<bool>&          depend_y
			  ) {
    //  std::cout<<"in rev_depend\n";
    depend_x[0] = depend_y[0];
    return true;
  }
  // Not sure this is ever needed
  virtual bool rev_depend(
			  const CppAD::vector<CppAD::AD<double> >&          parameter_x ,
			  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
			  CppAD::vector<bool>&                depend_x    ,
			  const CppAD::vector<bool>&          depend_y
			  ) {
    //    std::cout<<"in meta rev_depend\n";
    depend_x[0] = depend_y[0];
    return true;
  }
};

void track_nimble_atomic(nimble_atomic_base *obj, void *tape_mgr_ptr, std::vector<CppAD::local::atomic_index_info>* vec_ptr );

class atomic_lgamma_class;
class atomic_gammafn_class;
class atomic_pow_int_class;
class atomic_backsolve_class;
class atomic_forwardsolve_class;
class atomic_cholesky_class;
class atomic_matmult_class;
class atomic_matinverse_class;
class atomic_zround_class;
class nimDerivs_floor_class;
template<class ftor> class atomic_discrete_class;
class atomic_floor_class;
class atomic_ceil_class;
class atomic_ftrunc_class;
class atomic_nimRound_class;
class atomic_log_pow_int_class;
class atomic_zb_over_a_class;
class atomic_probit_class;
class atomic_iprobit_class;
class atomic_dyn_ind_get_class;
class atomic_dyn_ind_set_class;

void copy_CppADdouble_to_double(CppAD::AD<double> *first, CppAD::AD<double> *last, double *output);
void copy_CppADdouble_to_double(NimArrBase< CppAD::AD<double> > &from, NimArrBase< double > &to);
void copy_CppADdouble_to_double(CppAD::AD<double> &from, double &to);

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

/* A place to manage atomic objects associated with a tape */
/* The system of atomic_pair and the "new_" and "delete_" functions
   came from efforts to manage the atomic vec manager (a static inside CppAD) used
   when any atomic is created or destructed. This was necessary in some
   cases because nimble pass objects between DLLs, in a way that we hope to
   clean up in the future but has always worked but has caused some headaches
   when interacting with CppAD statics.  At this point the atomic_pair
   system may be entirely unnecessary because nimble_atomic_base does the job
   of recording the atomic vec manager on construction and ensuring
   the same one is used during destruction.  However we are leaving both
   systems in place in further needs are revealed. */
class nimble_CppAD_tape_mgr {
 public:
  typedef std::pair<nimble_atomic_base *, std::vector<CppAD::local::atomic_index_info>* > atomic_pair;
  std::vector<atomic_pair> atomic_ptrs;
  void add_atomic_ptr(nimble_atomic_base *new_atomic_ptr, std::vector<CppAD::local::atomic_index_info>* vec_ptr);
  CppAD::ADFun<double> *ADtape_;
  CppAD::ADFun<double>* &ADtape() {return ADtape_;};
  CppAD::local::ADTape<double>* internal_tape_ptr_;
  void reset();
  void set_internal_tape(CppAD::local::ADTape<double>* internal_tape_ptr);
  nimble_CppAD_tape_mgr();
  ~nimble_CppAD_tape_mgr();

  int lgamma_index[5];
  bool lgamma_exists[5];
  atomic_lgamma_class* new_atomic_lgamma(const std::string& name, int bO);
  void delete_atomic_lgamma(atomic_lgamma_class *atomic_lgamma);
  atomic_lgamma_class *get_atomic_lgamma(int baseOrder,
					 std::vector<CppAD::local::atomic_index_info>* vec_ptr);
  int gammafn_index;
  bool gammafn_exists;
  atomic_gammafn_class* new_atomic_gammafn(const std::string& name);
  void delete_atomic_gammafn(atomic_gammafn_class *atomic_gammafn);
  atomic_gammafn_class *get_atomic_gammafn(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int pow_int_index;
  bool pow_int_exists;
  atomic_pow_int_class* new_atomic_pow_int(const std::string& name);
  void delete_atomic_pow_int(atomic_pow_int_class *atomic_pow_int);
  atomic_pow_int_class *get_atomic_pow_int(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int zround_index;
  bool zround_exists;
  atomic_zround_class* new_atomic_zround(const std::string& name);
  void delete_atomic_zround(atomic_zround_class *atomic_zround);
  atomic_zround_class *get_atomic_zround(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int floor_index;
  bool floor_exists;
  atomic_floor_class* new_atomic_floor(const std::string& name);
  void delete_atomic_floor(atomic_floor_class *atomic_floor);
  atomic_floor_class *get_atomic_floor(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int ceil_index;
  bool ceil_exists;
  atomic_ceil_class* new_atomic_ceil(const std::string& name);
  void delete_atomic_ceil(atomic_ceil_class *atomic_ceil);
  atomic_ceil_class *get_atomic_ceil(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int ftrunc_index;
  bool ftrunc_exists;
  atomic_ftrunc_class* new_atomic_ftrunc(const std::string& name);
  void delete_atomic_ftrunc(atomic_ftrunc_class *atomic_ftrunc);
  atomic_ftrunc_class *get_atomic_ftrunc(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int nimRound_index;
  bool nimRound_exists;
  atomic_nimRound_class* new_atomic_nimRound(const std::string& name);
  void delete_atomic_nimRound(atomic_nimRound_class *atomic_nimRound);
  atomic_nimRound_class *get_atomic_nimRound(std::vector<CppAD::local::atomic_index_info>* vec_ptr);
  
  int log_pow_int_index;
  bool log_pow_int_exists;
  atomic_log_pow_int_class* new_atomic_log_pow_int(const std::string& name);
  void delete_atomic_log_pow_int(atomic_log_pow_int_class *atomic_log_pow_int);
  atomic_log_pow_int_class *get_atomic_log_pow_int(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int zb_over_a_index;
  bool zb_over_a_exists;
  atomic_zb_over_a_class* new_atomic_zb_over_a(const std::string& name);
  void delete_atomic_zb_over_a(atomic_zb_over_a_class *atomic_zb_over_a);
  atomic_zb_over_a_class *get_atomic_zb_over_a(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int probit_index;
  bool probit_exists;
  atomic_probit_class* new_atomic_probit(const std::string& name);
  void delete_atomic_probit(atomic_probit_class *atomic_probit);
  atomic_probit_class *get_atomic_probit(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int iprobit_index;
  bool iprobit_exists;
  atomic_iprobit_class* new_atomic_iprobit(const std::string& name);
  void delete_atomic_iprobit(atomic_iprobit_class *atomic_iprobit);
  atomic_iprobit_class *get_atomic_iprobit(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int dyn_ind_get_index;
  bool dyn_ind_get_exists;
  atomic_dyn_ind_get_class* new_atomic_dyn_ind_get(const std::string& name);
  void delete_atomic_dyn_ind_get(atomic_dyn_ind_get_class *atomic_dyn_ind_get);
  atomic_dyn_ind_get_class *get_atomic_dyn_ind_get(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  int dyn_ind_set_index;
  bool dyn_ind_set_exists;
  atomic_dyn_ind_set_class* new_atomic_dyn_ind_set(const std::string& name);
  void delete_atomic_dyn_ind_set(atomic_dyn_ind_set_class *atomic_dyn_ind_set);
  atomic_dyn_ind_set_class *get_atomic_dyn_ind_set(std::vector<CppAD::local::atomic_index_info>* vec_ptr);

  atomic_backsolve_class* new_atomic_backsolve(const std::string& name);
  void delete_atomic_backsolve(atomic_backsolve_class *atomic_backsolve);
  atomic_forwardsolve_class* new_atomic_forwardsolve(const std::string& name);
  void delete_atomic_forwardsolve(atomic_forwardsolve_class *atomic_forwardsolve);
  atomic_cholesky_class* new_atomic_cholesky(const std::string& name);
  void delete_atomic_cholesky(atomic_cholesky_class *atomic_cholesky);
  atomic_matmult_class* new_atomic_matmult(const std::string& name);
  void delete_atomic_matmult(atomic_matmult_class *atomic_matmult);
  atomic_matinverse_class* new_atomic_matinverse(const std::string& name);
  void delete_atomic_matinverse(atomic_matinverse_class *atomic_matinverse);

  std::vector<CppAD::AD<double> > dummyOutputs;
  void add_dummyOutput(CppAD::AD<double> &dummy);
  void sum_dummyOutputs_to_dependentVars(std::vector<CppAD::AD<double> > &depVars);
};

/* nimbleCppADinfoClass is the class to convey information from a nimbleFunction
   object
   to generic CppAD driver wrappers like calcjacobian.
   Each nimbleFunction enabled for CppAD will have an object of this class. */
class nimbleCppADinfoClass {
 public:
  std::vector<double> independentVars;
  std::vector<double> dynamicVars;
  std::vector< CppAD::AD<double> > independentVars_meta;
  std::vector< CppAD::AD<double> > dynamicVars_meta;
  bool metaFlag;
  nimble_CppAD_tape_mgr ADtape_mgr_;
  bool ADtape_empty() {return ADtape_mgr_.ADtape() == 0;}
  void ADtape_reset() {ADtape_mgr_.reset();}
  void set_internal_tape(CppAD::local::ADTape<double>* internal_tape_ptr) {
    ADtape_mgr_.set_internal_tape(internal_tape_ptr);}
  void add_dummyOutput(CppAD::AD<double> &dummy) {ADtape_mgr_.add_dummyOutput(dummy); };
  void sum_dummyOutputs_to_dependentVars(std::vector<CppAD::AD<double> > &depVars) {ADtape_mgr_.sum_dummyOutputs_to_dependentVars(depVars); }
  CppAD::ADFun<double>* &ADtape() {return ADtape_mgr_.ADtape(); }
  //  CppAD::ADFun<double> *ADtape;
  NodeVectorClassNew_derivs *updaterNV_;
  NodeVectorClassNew_derivs *updaterNV() {return updaterNV_;}
  nimbleCppADinfoClass& setUpdaterNV(NodeVectorClassNew_derivs &UNV) {
    updaterNV_ = &UNV;
    return *this;
  }
  bool updateModel_;
  bool &updateModel() {return updateModel_;}
  bool nodeFunPtrSet_;
  nodeFun *nodeFunPtr_;
  bool nodeFunPtrSet() {return nodeFunPtrSet_;}
  nodeFun *nodeFunPtr() {return nodeFunPtr_;}
  void set_nodeFunPtr(nodeFun *nfp) {
    nodeFunPtr_ = nfp;
    nodeFunPtrSet_ = true;
  }
  void clear_nodeFunPtr() {
    nodeFunPtr_ = 0;
    nodeFunPtrSet_ = false;
  }
 nimbleCppADinfoClass() :
  metaFlag(false),
    //    ADtape(0),
    updaterNV_(0),
    updateModel_(true),
    nodeFunPtrSet_(false),
    nodeFunPtr_(0)
      {}
  ~nimbleCppADinfoClass() {
    /* if(ADtape) { */
    /*   delete ADtape; */
    /*   ADtape = 0; */
    /* } */
  }
};

void add_dummyOutput(nimbleCppADinfoClass *ADinfoPtr, CppAD::AD<double> &dummy);

class nimbleCppADrecordingInfoClass {
 private:
  bool recording_;
  CppAD::tape_id_t tape_id_;
  CppAD::local::ADTape<double>* tape_handle_;
  nimbleCppADinfoClass *ADinfoPtr_;
  std::vector<CppAD::local::atomic_index_info>* atomic_vec_ptr_;
 public:
  bool& recording() {return recording_;}
  bool recording_cp() const {return recording_;}
  CppAD::tape_id_t& tape_id() {return tape_id_;}
  CppAD::local::ADTape<double>* &tape_handle() {return tape_handle_;}
  CppAD::tape_id_t tape_id_cp() const {return tape_id_;}
  CppAD::local::ADTape<double>* tape_handle_cp() const {return tape_handle_;}
  std::vector<CppAD::local::atomic_index_info>* &atomic_vec_ptr() {return atomic_vec_ptr_;}
  std::vector<CppAD::local::atomic_index_info>* atomic_vec_ptr_cp() const {return atomic_vec_ptr_;}

  nimbleCppADinfoClass* &ADinfoPtr() {return ADinfoPtr_;}
 nimbleCppADrecordingInfoClass(bool r_, CppAD::tape_id_t tid_, CppAD::local::ADTape<double>* th_) :
  recording_(r_),
    tape_id_(tid_),
    tape_handle_(th_),
    ADinfoPtr_(0),
    atomic_vec_ptr_(0) {}
 nimbleCppADrecordingInfoClass(CppAD::tape_id_t tid_, CppAD::local::ADTape<double>* th_, nimbleCppADinfoClass *ADinfoPtr) :
  recording_(false),
    tape_id_(tid_),
    tape_handle_(th_),
    ADinfoPtr_(ADinfoPtr),
    atomic_vec_ptr_(0)
    {}
 nimbleCppADrecordingInfoClass(CppAD::tape_id_t tid_, CppAD::local::ADTape<double>* th_, std::vector<CppAD::local::atomic_index_info>* avp_, nimbleCppADinfoClass *ADinfoPtr) :
  recording_(false),
    tape_id_(tid_),
    tape_handle_(th_),
    ADinfoPtr_(ADinfoPtr),
    atomic_vec_ptr_(avp_)
    {}
 nimbleCppADrecordingInfoClass(bool r_,  nimbleCppADinfoClass *ADinfoPtr) :
  recording_(r_),
    tape_id_(0),
    tape_handle_(0),
    ADinfoPtr_(ADinfoPtr),
    atomic_vec_ptr_(0)
    {}
 nimbleCppADrecordingInfoClass() : recording_(false) {}
};

void setValues_AD_AD_taping(NimArr<1, CppAD::AD<double> > &v, ManyVariablesMapAccessor &MVA_AD, ManyVariablesMapAccessor &MVA_orig, nimbleCppADrecordingInfoClass &recordingInfo);

void update_dynamicVars(NodeVectorClassNew_derivs &NV,
                        nimbleCppADinfoClass &ADinfo);
void update_dynamicVars(nimbleCppADinfoClass &ADinfo);
//std::vector<double> &dynamicVars,
//CppAD::ADFun<double>* &tapePtr);
void update_dynamicVars_meta(NodeVectorClassNew_derivs &NV,
                             nimbleCppADinfoClass &ADinfo);
void update_dynamicVars_meta(nimbleCppADinfoClass &ADinfo);
//std::vector< CppAD::AD<double> > &dynamicVars,
//CppAD::ADFun<double>* &tapePtr);

/* nimbleFunctionCppADbase is a base class to be inherited by all
   CppAD-enabled nimbleFunctions. Some of these functions might
   make more sense as stand-alone functions.  Let's see. */
// class nimbleFunctionCppADbase {
// public:
//   void getDerivs(nimbleCppADinfoClass &ADinfo,
//                  const NimArr<1, double> &derivOrders,
//                  const NimArr<1, double> &wrtVector,
//                  nimSmartPtr<NIMBLE_ADCLASS> &ansList);

//   nimSmartPtr<NIMBLE_ADCLASS> getDerivs_wrapper(nimbleCppADinfoClass &ADinfo,
//                                                 const NimArr<1, double> &derivOrders,
//                                                 const NimArr<1, double> &wrtVector){
//     nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
//     getDerivs(ADinfo, derivOrders, wrtVector, ansList);
//     return(ansList);
//   }

//   void getDerivs_meta(nimbleCppADinfoClass &ADinfo,
//                       const NimArr<1, double> &derivOrders,
//                       const NimArr<1, double> &wrtVector,
//                       const nimbleCppADrecordingInfoClass &nimRecInfo,
//                       nimSmartPtr<NIMBLE_ADCLASS_META> &ansList);

//   nimSmartPtr<NIMBLE_ADCLASS_META> getDerivs_wrapper_meta(nimbleCppADinfoClass &ADinfo,
//                                                           const NimArr<1, double> &derivOrders,
//                                                           const NimArr<1, double> &wrtVector,
//                                                           const nimbleCppADrecordingInfoClass &nimRecInfo){
//     nimSmartPtr<NIMBLE_ADCLASS_META> ansList = new NIMBLE_ADCLASS_META;
//     getDerivs_meta(ADinfo, derivOrders, wrtVector, nimRecInfo, ansList);
//     return(ansList);
//   }

//   void getDerivs_calculate_internal(nimbleCppADinfoClass &ADinfo,
//                                     //CppAD::ADFun<double>* &tapePtr,
//                                     NodeVectorClassNew_derivs &nodes,
//                                     const NimArr<1, double> &derivOrders,
//                                     const NimArr<1, double> &wrtVector,
//                                     bool do_update,
//                                     bool reset,
//                                     nimSmartPtr<NIMBLE_ADCLASS> ansList);
//   /* This form is not actually generated in code at the time of this writing:*/
//   nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(nimbleCppADinfoClass &ADinfo,
//                                                   //CppAD::ADFun<double>* &tapePtr,
//                                                   NodeVectorClassNew_derivs &nodes,
//                                                   const NimArr<1, double> &derivOrders,
//                                                   const NimArr<1, double> &wrtVector,
//                                                   bool do_update,
//                                                   bool reset){
//     nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
//     getDerivs_calculate_internal(ADinfo,// tapePtr,
//                                  nodes, derivOrders, wrtVector, do_update, reset, ansList);
//     return(ansList);
//   }
//   /* This is the form that would be generated in code, with no wrtVector*/
//   nimSmartPtr<NIMBLE_ADCLASS> nimDerivs_calculate(nimbleCppADinfoClass &ADinfo,
//                                                   //CppAD::ADFun<double>* &tapePtr,
//                                                   NodeVectorClassNew_derivs &nodes,
//                                                   const NimArr<1, double> &derivOrders,
//                                                   bool do_update,
//                                                   bool reset) {
//     NimArr<1, double> wrtVector; // with new default functionality, this could be set to simply length 1 with value -1
//     int totlen = nodes.model_wrt_accessor.getTotalLength();
//     wrtVector.setSize(totlen,
//                       false,
//                       false);
//     for(int ii = 0; ii < totlen; ++ii) {
//       wrtVector[ii] = ii + 1; // This gets -1 in use, as if it were from R.
//     }
//     nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
//     getDerivs_calculate_internal(ADinfo, //tapePtr,
//                                  nodes, derivOrders, wrtVector, do_update, reset, ansList);
//     return(ansList);
//   }
// };

inline nimbleCppADinfoClass& set_tape_ptr(nimbleCppADinfoClass &ADtapeSetup,
                                          CppAD::ADFun<double>* &ADtapePtr,
                                          bool do_this) {
  if(!ADtapePtr) ADtapePtr = new CppAD::ADFun<double>;
  if(do_this) ADtapeSetup.ADtape() = ADtapePtr;
  return ADtapeSetup;
}

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
NimArr<1, double> make_vector_if_necessary(NimArr<1, int>);

template <typename T>
NimArr<1, double> make_vector_if_necessary(const T &x) { // This case catches eigen ops
  typedef typename Eigen::internal::traits<T>::Scalar Scalar;
  Eigen::Matrix<Scalar, Dynamic, Dynamic> materialized_x = x;
  size_t len = x.size();
  NimArr<1, Scalar> NimArr_x;
  NimArr_x.setSize(len, 0, 0);
  Map< MatrixXd > Eig_NimArr_x(NimArr_x.getPtr(),len,1);
  Eig_NimArr_x = materialized_x;
  return make_vector_if_necessary(NimArr_x);
}


#endif
