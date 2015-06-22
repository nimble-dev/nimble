#ifndef __NODEFUN
#define __NODEFUN


class nodeFun : public NamedObjects { 
  public: 
   virtual double calculate()=0; 
   virtual double calculateDiff()=0;
   virtual void simulate()=0;
   virtual double getLogProb()=0; 
 };
 
#endif
