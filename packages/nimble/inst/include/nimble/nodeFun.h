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

#ifndef __NODEFUN
#define __NODEFUN
#include "NimArr.h"
#include "smartPtrs.h"
#include "NamedObjects.h"
#include <cppad/cppad.hpp>

#define ADvector CppAD::vector

class NIMBLE_ADCLASS : public NamedObjects, public pointedToBase {
 public:
  NimArr<1, double> value;
  NimArr<2, double> jacobian;
  NimArr<3, double> hessian;
  NimArr<4, double> thirdDerivs;
  SEXP RObjectPointer;
  bool RCopiedFlag;
  void copyFromSEXP(SEXP S_nimList_);
  SEXP copyToSEXP();
  void createNewSEXP();
  void resetFlags();
  void copyFromRobject(SEXP Robject);
  NIMBLE_ADCLASS();
};

class NIMBLE_ADCLASS_META : public pointedToBase {
 public:
  NimArr<1, CppAD::AD<double> > value;
  NimArr<2, CppAD::AD<double> > jacobian;
  NimArr<3, CppAD::AD<double> > hessian;
};

class NodeVectorClassNew_derivs;

class atomic_extraInputObject : public CppAD::atomic_base<double> {
  public:
 atomic_extraInputObject(const std::string& name,
			 NodeVectorClassNew_derivs* NV);
 private:
  NodeVectorClassNew_derivs* NV_;//for access to model_extraInput_accessor;
  
  virtual bool forward(
		       size_t                    p ,
		       size_t                    q ,
		       const ADvector<bool>&      vx ,
		       ADvector<bool>&      vy ,
		       const ADvector<double>&    tx ,
		       ADvector<double>&    ty
		       );
  virtual
    bool reverse(
		 size_t                    q ,
		 const ADvector<double>&    tx ,
		 const ADvector<double>&    ty ,
		 ADvector<double>&    px ,
		 const ADvector<double>&    py
		 );
  virtual
    bool for_sparse_jac(
			size_t                     q ,
			const ADvector<bool>&   rt ,
			ADvector<bool>&         st ,
			const ADvector<double>&      x );
  virtual
    bool rev_sparse_jac(
			size_t                     q ,
			const ADvector<bool>&   rt ,
			ADvector<bool>&         st ,
			const ADvector<double>&      x );
};


class atomic_extraOutputObject : public CppAD::atomic_base<double> {
  public:
 atomic_extraOutputObject(const std::string& name,
			  NodeVectorClassNew_derivs* NV);
 private:
  NodeVectorClassNew_derivs* NV_;//for access to model_modelOutput_accessor;
  
  virtual bool forward(
		       size_t                    p ,
		       size_t                    q ,
		       const ADvector<bool>&      vx ,
		       ADvector<bool>&      vy ,
		       const ADvector<double>&    tx ,
		       ADvector<double>&    ty
		       );
  virtual
    bool reverse(
		 size_t                    q ,
		 const ADvector<double>&    tx ,
		 const ADvector<double>&    ty ,
		 ADvector<double>&    px ,
		 const ADvector<double>&    py
		 );
  virtual
    bool for_sparse_jac(
			size_t                     q ,
			const ADvector<bool>&   r ,
			ADvector<bool>&         s ,
			const ADvector<double>&      x );
  virtual
    bool rev_sparse_jac(
			size_t                     q ,
			const ADvector<bool>&   rt ,
			ADvector<bool>&         st ,
			const ADvector<double>&      x );
};


// This contains the indexed information -- often indices themselves but also any partially evaluated values --
// to calculate/simulate/getLogProb for one node withina new node function.
// Typically we'll have a vector of these in a node function or some such way of packaging information.
class indexedNodeInfo {
 public:
  vector<double> info; // although the main purposes of this is for indices, it will sometimes hold constants (sometimes from partially evaluated expressions)
  indexedNodeInfo() {}
  template<typename itertype>
    indexedNodeInfo(itertype iStart, int ncol, int stride = 1) {
    info.reserve(ncol);
    for(int i = 0; i < ncol; i++, iStart += stride) {
      info.push_back(*iStart);
    }
  }
  indexedNodeInfo(const vector<double>& newinfo) : info(newinfo) {}
  operator const vector<double>& () const { return info; }
};

// Derived classes can set up indexedNodeInfo needs in different ways.
// The default is to have a vector<indexedNodeInfo> and then pass along index vectors.
class nodeFun : public NamedObjects { 
 public:  // etc. put iNI into all cases
  vector<indexedNodeInfo> indexedNodeInfoTable;
  vector<indexedNodeInfo> *getIndexedNodeInfoTablePtr() {return(&indexedNodeInfoTable);}
  nodeFun() {
    namedObjects["indexedNodeInfoTable"] = getIndexedNodeInfoTablePtr();
  }

  virtual double calculate(const indexedNodeInfo &iNI) const =0;
  virtual CppAD::AD<double> calculate_ADproxyModel(const indexedNodeInfo &iNI) const {
    printf("Error in C++: Dummy calculate_ADproxyModel is being used\n.");
    return(CppAD::AD<double>(0)); // need a default in case generating deriv functions is turned off.
  };
  CppAD::AD<double> calculateBlock_ADproxyModel(int operand) const { return calculate_ADproxyModel(indexedNodeInfoTable[operand]); }

  //  virtual void calculateWithArgs_deriv(const indexedNodeInfo &iNI, const NimArr<1, double> & ARG2_nimDerivsOrders_, const NimArr<1, double> & ARG3_wrtVector_, nimSmartPtr<NIMBLE_ADCLASS> ansList) = 0;
  virtual double calculateDiff(const indexedNodeInfo &iNI) const =0;
  // virtual CppAD::AD<double> calculateDiff_ADproxyModel(const indexedNodeInfo &iNI) const =0;
  virtual void simulate(const indexedNodeInfo &iNI) const =0;
  virtual double getLogProb(const indexedNodeInfo &iNI) const =0;

  virtual double getParam_0D_double(int paramID, const indexedNodeInfo &iNI) const {return(0./0.);} 
  virtual NimArr<1, double> getParam_1D_double(int paramID, const indexedNodeInfo &iNI) const {NimArr<1, double> ans; return(ans);}
  virtual NimArr<2, double> getParam_2D_double(int paramID, const indexedNodeInfo &iNI) const {NimArr<2, double> ans; return(ans);}

  virtual double getBound_0D_double(int boundID, const indexedNodeInfo &iNI) const {return(0./0.);} 
  virtual NimArr<1, double> getBound_1D_double(int boundID, const indexedNodeInfo &iNI) const {NimArr<1, double> ans; return(ans);}
  virtual NimArr<2, double> getBound_2D_double(int boundID, const indexedNodeInfo &iNI) const {NimArr<2, double> ans; return(ans);}

  double calculateBlock(int operand) const { return calculate(indexedNodeInfoTable[operand]); }
  /* void calculateWithArgs_derivBlock(int operand, NimArr<1, double> &derivOrders, const NimArr<1, double> &wrtVector, nimSmartPtr<NIMBLE_ADCLASS> ansList) { */
  /* 	  return(calculateWithArgs_deriv(indexedNodeInfoTable[operand], derivOrders, wrtVector, ansList)); */
  /* } */
  double calculateDiffBlock(int operand) const { return calculateDiff(indexedNodeInfoTable[operand]); }
  double getLogProbBlock(int operand) const { return getLogProb(indexedNodeInfoTable[operand]); }
  void simulateBlock(int operand) const { simulate(indexedNodeInfoTable[operand]); }

  double getParam_0D_double_block(int paramID, int operand) const {
    return(getParam_0D_double(paramID, indexedNodeInfoTable[operand]));
  }
  NimArr<1, double> getParam_1D_double_block(int paramID, int operand) const {
    return(getParam_1D_double(paramID, indexedNodeInfoTable[operand]));
  }
  NimArr<2, double> getParam_2D_double_block(int paramID, int operand) const {
    return(getParam_2D_double(paramID, indexedNodeInfoTable[operand]));
  }
  double getBound_0D_double_block(int boundID, int operand) const {
    return(getBound_0D_double(boundID, indexedNodeInfoTable[operand]));
  }
  NimArr<1, double> getBound_1D_double_block(int boundID, int operand) const {
    return(getBound_1D_double(boundID, indexedNodeInfoTable[operand]));
  }
  NimArr<2, double> getBound_2D_double_block(int boundID, int operand) const {
    return(getBound_2D_double(boundID, indexedNodeInfoTable[operand]));
  }

  // useADreconfigure functions (version "2")
  void initialize_AD_model_before_recording(NodeVectorClassNew_derivs &NV);
  virtual void set_atomic_info_from_nodeFun(std::vector<CppAD::local::atomic_index_info>* vec_ptr);
  virtual void set_tape_ptr_from_nodeFun(CppAD::tape_id_t tape_id,
					 CppAD::local::ADTape<double>* tape_handle_,
					 bool recover);
  // virtual void setup_extraInput_step(NodeVectorClassNew_derivs &NV);
  virtual void setup_extraOutput_step(NodeVectorClassNew_derivs &NV,
				      CppAD::AD<double> &logProb);

  // Next 3 functions are virtual to ensure the code in the model DLL
  // will be used so that the correct CppAD globals / statics will be found.
  void recordTape(NodeVectorClassNew_derivs &NV);
  virtual void setTapeIndependent(std::vector< CppAD::AD<double> > &independentVars);
  virtual void finishADFun(CppAD::ADFun< double > &ADtape,
			   std::vector< CppAD::AD<double> > &independentVars,
			   std::vector< CppAD::AD<double> > &dependentVars);
  virtual void runTape(CppAD::ADFun< double > &ADtape,
		       std::vector< double > &independentVars,
		       std::vector< double > &dependentVars,
		       const NimArr<1, double> &derivOrders,
		       nimSmartPtr<NIMBLE_ADCLASS> &ansList);
  virtual atomic_extraInputObject*
    runExtraInputObject(NodeVectorClassNew_derivs &NV,
			std::vector< CppAD::AD<double> > &extraInputDummyInput,
			std::vector< CppAD::AD<double> > &extraInputResult);
  void
    delete_extraInputObject(NodeVectorClassNew_derivs &NV);
  virtual atomic_extraOutputObject*
    runExtraOutputObject(NodeVectorClassNew_derivs &NV,
			 CppAD::AD<double> &logProb);
  void
    delete_extraOutputObject(NodeVectorClassNew_derivs &NV);
  virtual CppAD::AD<double> call_calculate_ADproxyModel(NodeVectorClassNew_derivs &NV);
};

#endif
