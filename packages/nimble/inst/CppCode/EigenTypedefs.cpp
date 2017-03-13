#include <nimble/EigenTypedefs.h>

SEXP  EIGEN_EIGENCLASS::copyToSEXP (  )  {
  if(!RCopiedFlag) {
    EIGEN_EIGENCLASS_R::copyToSEXP();
    RCopiedFlag = true;
  }
  return(RObjectPointer);
}

EIGEN_EIGENCLASS::EIGEN_EIGENCLASS(){
  //std::cout<<"Constructing EIGEN_EIGENCLASS\n";
  namedObjects["values"]=&values;
  namedObjects["vectors"]=&vectors;
  RCopiedFlag = false;
}
