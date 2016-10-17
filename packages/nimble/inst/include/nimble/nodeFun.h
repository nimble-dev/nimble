#ifndef __NODEFUN
#define __NODEFUN
#include "NimArr.h"

// this contains the indexed information -- often indices themselves but also any partially evaluated values -- to calculate/simulate/getLogProb for one node withina new node function
// typically we'll have a vector of these in a node function or some such way of packaging information 
class indexedNodeInfo {
 public:
  vector<double> info; // although the main purposes of this is for indices, it will sometimes hold constants (sometimes from partially evaluated expressions)
  indexedNodeInfo() {};
  template<typename itertype>
    indexedNodeInfo(itertype iStart, int ncol, int stride = 1) {
    //  std::cout<<"in templated constructor for indexedNodeInfo iStart "<<iStart<<" ncol "<<ncol<<" stride "<<stride<< "\n";
    info.reserve(ncol);
    for(int i = 0; i < ncol; i++, iStart += stride) {
      //      std::cout<<"pushing back "<<*iStart<<"\n";
      info.push_back(*iStart);
    }
  }
  indexedNodeInfo(vector<double> newinfo) {info = newinfo;};
  operator vector<double>() const {return(info);}
};

// this will be the information for a block of indexedNodeInfo to use within one call to a new node function to trigger operations for one or more nodes 
// first simple version will be a vector of integers that index a vector<indexedNodeInfo> in the node function
// but future modes of operation can be added.
class useInfoForIndexedNodeInfo {
 public:
  vector<int> indicesForIndexedNodeInfo;
};

/* class nodeFun_old : public NamedObjects {  */
/*   public:  */
/*    virtual double calculate()=0;  */
/*    virtual double calculateDiff()=0; */
/*    virtual void simulate()=0; */
/*    virtual double getLogProb()=0;  */
/*    virtual double getParam_0D_double(int paramID) {return(0./0.);} */
/*    virtual NimArr<1, double> getParam_1D_double(int paramID) {NimArr<1, double> ans; return(ans);} */
/*    virtual NimArr<2, double> getParam_2D_double(int paramID) {NimArr<2, double> ans; return(ans);} */
/* }; */

// derived classes can set up indexedNodeInfo needs in different ways
// a default will be to have a vector<indexedNodeInfo> and then pass along index vectors
class nodeFun : public NamedObjects { 
 public:  // etc. put iNI into all cases
  vector<indexedNodeInfo> indexedNodeInfoTable;
  vector<indexedNodeInfo> *getIndexedNodeInfoTablePtr() {return(&indexedNodeInfoTable);}
  // virtual void shout() {PRINTF("shouting\n");}
  // carry these here to allow compilation for now -- can the old and new systems coexist?
   /* virtual double calculate()=0;  */
   /* virtual double calculateDiff()=0; */
   /* virtual void simulate()=0; */
   /* virtual double getLogProb()=0;  */
   /* virtual double getParam_0D_double(int paramID) {return(0./0.);} */
   /* virtual NimArr<1, double> getParam_1D_double(int paramID) {NimArr<1, double> ans; return(ans);} */
   /* virtual NimArr<2, double> getParam_2D_double(int paramID) {NimArr<2, double> ans; return(ans);} */
  
  nodeFun() {
    namedObjects["indexedNodeInfoTable"] = getIndexedNodeInfoTablePtr();
  };
  
  virtual double calculate(const indexedNodeInfo &iNI) const =0;
  virtual double calculateDiff(const indexedNodeInfo &iNI) const =0;
  virtual void simulate(const indexedNodeInfo &iNI) const =0;
  virtual double getLogProb(const indexedNodeInfo &iNI) const =0;

  virtual double getParam_0D_double(int paramID, const indexedNodeInfo &iNI) const {return(0./0.);} 
  virtual NimArr<1, double> getParam_1D_double(int paramID, const indexedNodeInfo &iNI) const {NimArr<1, double> ans; return(ans);}
  virtual NimArr<2, double> getParam_2D_double(int paramID, const indexedNodeInfo &iNI) const {NimArr<2, double> ans; return(ans);}

  double calculateBlock(const useInfoForIndexedNodeInfo &biNI) const {
    double ans(0);
    vector<int>::const_iterator iIndex(biNI.indicesForIndexedNodeInfo.begin());
    vector<int>::const_iterator iIndexEnd(biNI.indicesForIndexedNodeInfo.end());
    for(; iIndex != iIndexEnd; iIndex++) {
      ans += calculate(indexedNodeInfoTable[ *iIndex ]);
    }
    return(ans);
  }; 
  double calculateDiffBlock(const useInfoForIndexedNodeInfo &biNI) const {
    double ans(0);
    vector<int>::const_iterator iIndex(biNI.indicesForIndexedNodeInfo.begin());
    vector<int>::const_iterator iIndexEnd(biNI.indicesForIndexedNodeInfo.end());
    for(; iIndex != iIndexEnd; iIndex++) {
      ans += calculateDiff(indexedNodeInfoTable[ *iIndex ]);
    }
    return(ans);
  }; 
  double getLogProbBlock(const useInfoForIndexedNodeInfo &biNI) const {
    double ans(0);
    vector<int>::const_iterator iIndex(biNI.indicesForIndexedNodeInfo.begin());
    vector<int>::const_iterator iIndexEnd(biNI.indicesForIndexedNodeInfo.end());
    for(; iIndex != iIndexEnd; iIndex++) {
      ans += getLogProb(indexedNodeInfoTable[ *iIndex ]);
    }
    return(ans);
  }; 
  void simulateBlock(const useInfoForIndexedNodeInfo &biNI) const {
    vector<int>::const_iterator iIndex(biNI.indicesForIndexedNodeInfo.begin());
    vector<int>::const_iterator iIndexEnd(biNI.indicesForIndexedNodeInfo.end());
    for(; iIndex != iIndexEnd; iIndex++) {
      simulate(indexedNodeInfoTable[ *iIndex ]);
    }
  };
  double getParam_0D_double_block(int paramID, const useInfoForIndexedNodeInfo &biNI) const {
    return(getParam_0D_double(paramID, indexedNodeInfoTable[ biNI.indicesForIndexedNodeInfo[0] ]));
  }
  NimArr<1, double> getParam_1D_double_block(int paramID, const useInfoForIndexedNodeInfo &biNI) const {
    return(getParam_1D_double(paramID, indexedNodeInfoTable[ biNI.indicesForIndexedNodeInfo[0] ]));
  }
  NimArr<2, double> getParam_2D_double_block(int paramID, const useInfoForIndexedNodeInfo &biNI) const {
    return(getParam_2D_double(paramID, indexedNodeInfoTable[ biNI.indicesForIndexedNodeInfo[0] ]));
  }
};

#endif
