#include <nimble/EigenTypedefs.h>

SEXP  EIGEN_EIGENCLASS::copyToSEXP (  )  {
  if(!RCopiedFlag) {
    EIGEN_EIGENCLASS_R::copyToSEXP();
    RCopiedFlag = true;
  }
  return(RObjectPointer);
}

void EIGEN_EIGENCLASS::resetFlags () {
  RCopiedFlag = false;	
}

EIGEN_EIGENCLASS::EIGEN_EIGENCLASS(){
  //std::cout<<"Constructing EIGEN_EIGENCLASS\n";
  namedObjects["values"]=&values;
  namedObjects["vectors"]=&vectors;
  RCopiedFlag = false;
}

SEXP  EIGEN_SVDCLASS::copyToSEXP (  )  {
  if(!RCopiedFlag) {
    EIGEN_SVDCLASS_R::copyToSEXP();
    RCopiedFlag = true;
  }
  return(RObjectPointer);
}

EIGEN_SVDCLASS::EIGEN_SVDCLASS(){
  //std::cout<<"Constructing EIGEN_SVDCLASS\n";
  namedObjects["d"]=&d;
  namedObjects["u"]=&u;
  namedObjects["v"]=&v;
  RCopiedFlag = false;
}

void EIGEN_SVDCLASS::resetFlags () {
  RCopiedFlag = false;	
}
