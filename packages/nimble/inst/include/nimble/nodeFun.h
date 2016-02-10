#ifndef __NODEFUN
#define __NODEFUN
#include "NimArr.h"

class nodeFun : public NamedObjects { 
  public: 
   virtual double calculate()=0; 
   virtual double calculateDiff()=0;
   virtual void simulate()=0;
   virtual double getLogProb()=0; 
   virtual double getParam_0D_double(int paramID) {return(0./0.);}
   virtual NimArr<1, double> getParam_1D_double(int paramID) {NimArr<1, double> ans; return(ans);}
   virtual NimArr<2, double> getParam_2D_double(int paramID) {NimArr<2, double> ans; return(ans);}
};
 
#endif
