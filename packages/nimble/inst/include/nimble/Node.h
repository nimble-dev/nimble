#ifndef __NODE
#define __NODE

//#include "NimArr.h"
#include "Values.h"
#include <Rinternals.h>

class Node {
public:
  //  CalcStates *calcStates;
  virtual void simulate()=0;
  virtual double calculate()=0;
  virtual double getLogProb()=0;
  virtual double calculateDiff()=0;
  virtual double get_param_0D_double(int paramID) {return(0./0.);}
};

class DetermNode : public Node{
public:
  virtual void simulate()=0;
  virtual double calculate() { simulate(); return(0.); } // default behavior
  virtual double calculateDiff() { simulate(); return(0.); }
  virtual double getLogProb() { return(0.); } // default behavior
};

class StochNode : public Node {
public:
  virtual void simulate()=0;
  virtual double calculate()=0;
  virtual double calculateDiff()=0;
  virtual double getLogProb()=0;
};

#define DEFINE_DETERMNODE_CLASS(name) \
class name : public DetermNode {	\
public:				\
virtual void simulate();		\
};

// If the default getLogProb is ok, only simulate and calculate need specialization
#define DEFINE_STOCHNODE_SCALAR_CLASS(name) \
  class name : public StochNode {	\
  public:				\
  virtual void simulate();		\
  virtual double calculate();		\
  virtual double calculateDiff(); \
  virtual double getLogProb();      \
  };

// If getLogProb also needs specialization:
#define DEFINE_STOCHNODE_GENERAL_CLASS(name) \
  class name : public StochNode {	\
  public:				\
  virtual void simulate();		\
  virtual double calculate();         \
  virtual double calculateDiff();     \
  virtual double getLogProb();        \
  };

extern "C" {
  SEXP callSimulate(SEXP Sextptr);
  SEXP callCalculate(SEXP Sextptr);
  SEXP callCalculateDiff(SEXP Sextptr);
  SEXP callGetLogProb(SEXP Sextptr);
}

#endif
