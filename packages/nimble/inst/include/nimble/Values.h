#ifndef __VALUES
#define __VALUES

//#include<vector>
#include<string>
// #include<tr1/unordered_set>			//Cliff debug
//#include<map>
using std::string;
#include "NamedObjects.h"
#include "Utils.h"




class Values : public NamedObjects {
public:
  /* typedef vector< double > valueType; */
  /* typedef vector< vector< double > > valueVecType; */
  /* typedef valueVecType::iterator valueVecIterator; */

 // int SampleSize; // deprecate this
  int numRows;
  int getsize(){return(numRows);};
  virtual void resize(int nrow)=0;
  
  string buildName;
  string getMVBuildName(){ return buildName;}  ;
  Values() { buildName = "missing"; numRows = 0;};
};




//
//
//
//}
/* typedef Values* (*ValuesGenerator)(); */

/* class ValuesFactory { */
/*  public: */
/*   map<string, ValuesGenerator> ValuesGeneratorMap; */
/* //  tr1::unordered_set<Values *> ValuesCollection;   //Cliff debug */
/*   Values *makeNew(string &typelabel); */
/*   void registerGenerator(string typelabel, ValuesGenerator vg); */
/*   void removeValues(Values *empty); */
/*   ValuesFactory(){}; */
/*   ~ValuesFactory(); */
/* }; */



/* extern "C" { */
/*   SEXP newModelValues(SEXP typelabel); */
  
/* } */

/* extern ValuesFactory valuesFactory; */

/* class ValuesRegisterer { */
/*  public: */
/*   ValuesRegisterer(string typelabel, ValuesGenerator vg) { */
/*     PRINTF("Trying to register\n"); */
/*     valuesFactory.registerGenerator(typelabel, vg); */
/*   } */
/* }; */


//void valuesFinalizer(SEXP Sv);

    //void allocate(vector< vector <double> > *vv, int sampleSize, int variableSize);

#endif
