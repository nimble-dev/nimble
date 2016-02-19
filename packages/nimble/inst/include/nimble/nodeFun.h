#ifndef __NODEFUN
#define __NODEFUN
#include "NimArr.h"

// this contains the indexed information -- often indices themselves but also any partially evaluated values -- to calculate/simulate/getLogProb for one node withina new node function
// typically we'll have a vector of these in a node function or some such way of packaging information 
class indexedNodeInfo {
 public:
  vector<int> info;
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

  // carry these here to allow compilation for now -- can the old and new systems coexist?
   virtual double calculate()=0; 
   virtual double calculateDiff()=0;
   virtual void simulate()=0;
   virtual double getLogProb()=0; 
   virtual double getParam_0D_double(int paramID) {return(0./0.);}
   virtual NimArr<1, double> getParam_1D_double(int paramID) {NimArr<1, double> ans; return(ans);}
   virtual NimArr<2, double> getParam_2D_double(int paramID) {NimArr<2, double> ans; return(ans);}
  
  /* virtual double getParam_0D_double(int paramID, const indexedNodeInfo &iNI) {return(0./0.);} */
  /* virtual NimArr<1, double> getParam_1D_double(int paramID) {NimArr<1, double> ans; return(ans);} */
  /* virtual NimArr<2, double> getParam_2D_double(int paramID) {NimArr<2, double> ans; return(ans);} */
  
  virtual double calculate_indexedNodeInfo(const indexedNodeInfo &iNI)=0;
   // same for calculateDiff, simulate and getLogProb
  double calculateBlock(const useInfoForIndexedNodeInfo &biNI); // same implementation should work for all derived classes
};


#endif
