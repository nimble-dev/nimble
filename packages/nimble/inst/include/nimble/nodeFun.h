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
  virtual double calculateDiff(const indexedNodeInfo &iNI) const =0;
  virtual void simulate(const indexedNodeInfo &iNI) const =0;
  virtual double getLogProb(const indexedNodeInfo &iNI) const =0;

  virtual double getParam_0D_double(int paramID, const indexedNodeInfo &iNI) const {return(0./0.);} 
  virtual NimArr<1, double> getParam_1D_double(int paramID, const indexedNodeInfo &iNI) const {NimArr<1, double> ans; return(ans);}
  virtual NimArr<2, double> getParam_2D_double(int paramID, const indexedNodeInfo &iNI) const {NimArr<2, double> ans; return(ans);}

  virtual double getBound_0D_double(int boundID, const indexedNodeInfo &iNI) const {return(0./0.);} 
  virtual NimArr<1, double> getBound_1D_double(int boundID, const indexedNodeInfo &iNI) const {NimArr<1, double> ans; return(ans);}
  virtual NimArr<2, double> getBound_2D_double(int boundID, const indexedNodeInfo &iNI) const {NimArr<2, double> ans; return(ans);}

  double calculateBlock(int operand) const { return calculate(indexedNodeInfoTable[operand]); }
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
};

#endif
