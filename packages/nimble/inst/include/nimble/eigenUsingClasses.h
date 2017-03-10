#ifndef __EIGENUSINGCLASSES
#define __EIGENUSINGCLASSES

class EIGEN_EIGENCLASS : public NamedObjects, public pointedToBase {
public:
  NimArr<1, double> values;
  NimArr<2, double> vectors;
  SEXP RObjectPointer;
  bool RCopiedFlag;
  
SEXP  copyToSEXP (   );
void  createNewSEXP (  );
void  copyFromSEXP ( SEXP S_nimList_ );
 EIGEN_EIGENCLASS(){	
    namedObjects["values"]=&values;
	namedObjects["vectors"]=&vectors;
	RCopiedFlag = false;
	RObjectPointer = NULL;
 };
};


template<class Derived>
nimSmartPtr<EIGEN_EIGENCLASS>   EIGEN_EIGEN(const Eigen::MatrixBase<Derived> &x, bool valuesOnly) {
    nimSmartPtr<EIGEN_EIGENCLASS> returnClass = new EIGEN_EIGENCLASS;
	(*returnClass).values.initialize(0, 0, x.rows());
	Map<VectorXd> Eig_eigVals((*returnClass).values.getPtr(),x.rows());
	if(x.rows() != x.cols()) {
	 _nimble_global_output <<"Run-time size error: expected matrix argument to eigen() to be square."<<"\n"; nimble_print_to_R(_nimble_global_output);
	}
    Eigen::DecompositionOptions eigOpts = valuesOnly ? EigenvaluesOnly : ComputeEigenvectors;
	SelfAdjointEigenSolver<MatrixXd> solver(x, eigOpts); // The MatrixXd here doesn't seem generic, but I couldn't get it to work otherwise and it would be odd to do an Eigen decomposition on anything else. -Perry
	Eig_eigVals = solver.eigenvalues().reverse();
	if(!valuesOnly){
	  (*returnClass).vectors.initialize(0, 0, x.rows(), x.cols());
	  Map<MatrixXd> Eig_eigVecs((*returnClass).vectors.getPtr(),x.rows(),x.cols());
	  Eig_eigVecs = solver.eigenvectors().reverse();	
	}
	return(returnClass);
};



class EIGEN_SVDCLASS : public NamedObjects, public pointedToBase {
public:
  NimArr<1, double> d;
  NimArr<2, double> u;
  NimArr<2, double> v;
  SEXP RObjectPointer;
  bool RCopiedFlag;
  
SEXP  copyToSEXP (   );
void  createNewSEXP (  );
void  copyFromSEXP ( SEXP S_nimList_ );
 EIGEN_SVDCLASS (  ) {
	namedObjects["d"]=&d;
	namedObjects["u"]=&u;
	namedObjects["v"]=&v;
	RCopiedFlag = false;
	RObjectPointer = NULL;
 };
};



template<class Derived>
nimSmartPtr<EIGEN_SVDCLASS>   EIGEN_SVD(const Eigen::MatrixBase<Derived> &x, int vectors) {
    nimSmartPtr<EIGEN_SVDCLASS> returnClass = new EIGEN_SVDCLASS;
	int nu = min(x.rows(), x.cols());
	(*returnClass).d.initialize(0, 0, nu);
	Map<VectorXd> Svd_d((*returnClass).d.getPtr(), nu);
	JacobiSVD<MatrixXd> svd;

	/* note: if nu > 16, bidiagonialization algo. is recommended on eigen website.  not currently available w/ nimble's version of eigen, but may be in future. */
	if(vectors == 0){
	  svd.compute(x);
	}
	else{
	  int leftSVs = nu;
	  int rightSVs = nu;
	  if(vectors == 1){
	  	svd.compute(x, ComputeThinU | ComputeThinV);
	  }
	  if(vectors == 2){
	    leftSVs = x.rows();
		rightSVs = x.cols();
		svd.compute(x, ComputeFullU | ComputeFullV);
	  }
	  (*returnClass).u.initialize(0, 0, x.rows(), leftSVs);
	  (*returnClass).v.initialize(0, 0, x.cols(), rightSVs);
	  Map<MatrixXd> Svd_u((*returnClass).u.getPtr(), x.rows(), leftSVs);
	  Map<MatrixXd> Svd_v((*returnClass).v.getPtr(), x.cols(), rightSVs);
	  Svd_u = svd.matrixU();
	  Svd_v = svd.matrixV();
	}
	Svd_d = svd.singularValues(); 
	return(returnClass);
};

extern "C" {
SEXP C_nimEigen(SEXP S_x, SEXP S_valuesOnly, SEXP returnList);
SEXP C_nimSvd(SEXP S_x, SEXP S_vectors, SEXP returnList);
}

#endif
