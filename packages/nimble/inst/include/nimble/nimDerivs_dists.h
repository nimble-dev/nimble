#ifndef _NIMDERIVS_DISTS__
#define _NIMDERIVS_DISTS__

#define _USE_ATOMICS_IN_DMNORM // remove this line to turn off atomics in dmnorm

// The code here benefited from studying source code of Template Model Builder(TMB).
// In the end, there is little if any code directly used verbatim from TMB under open-source licensing.
// Where possible we have noted pieces of code modified from TMB.

#include <cppad/cppad.hpp>
#include <R.h>
#include <Rinternals.h>
#include "nimDerivs_atomic_classes.h"
//#include <cppad/utility/nan.hpp>

#include "NimArr.h" // This includes Rmath.h via Utils.h.

/* discrete-round is four lines of simple code from TMB, but I don't think we even call it. */
/* because there was a bug using CPPAD_DISCRETE_FUNCTION from CppAD with base2ad. */
/* That bug was fixed, but after we worked around it anyway. */
inline double discrete_round(const double &x)
{     
  double out_x = round(x);
  return(out_x);
}
CPPAD_DISCRETE_FUNCTION(double, discrete_round)

// Functions here connect from nimble-generated C++ to CppAD.
// Note that these functions are called rarely, typically once, because
// they are only used when recording a CppAD tape of operations.  This tape
// is then the object used to obtain derivatives of the recorded operations.
// For this reason, function call overhead is not a major concern.
//
// For clarify, we have named functions uniquely rather than relying on
// overloading.  Calls to these functions are typically code-generated.

/* dnorm: normal distribution */
/* This case is modified from TMB code so that give_log can be different
   when the tape is used from when it is recorded. */
template<class Type>
Type nimDerivs_dnorm(Type x, Type mean, Type sd, Type give_log)
{
  Type res = -log(Type(sqrt(2*M_PI))*sd)-Type(.5)*pow((x-mean)/sd,2);
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_dnorm_logFixed(Type x, Type mean, Type sd, int give_log)
{
  Type res = -log(Type(sqrt(2*M_PI))*sd)-Type(.5)*pow((x-mean)/sd,2);
  if(!give_log) res = exp(res);
  return(res);
}


/* dmnorm: Multivariate normal distribution */
template<class Type>
Type nimDerivs_nimArr_dmnorm_chol(NimArr<1, Type> &x, NimArr<1, Type> &mean, NimArr<2, Type> &chol, Type prec_param, Type give_log, Type overwrite_inputs) { 
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  if(!CppAD::Constant(prec_param))
    std::cout<<"Warning: In dmnorm, prec_param with value = "<<CppAD::Value(prec_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make prec_param a constant to avoid this warning."<<std::endl;

  int n = x.dimSize(0);
  int i;
  Type dens = Type(-n * M_LN_SQRT_2PI);
  Type sumDens = Type(0.);
  for(i = 0; i < n*n; i += n + 1) 
	  sumDens += log(chol[i]);

  dens += CppAD::CondExpEq(prec_param, Type(1), sumDens, -sumDens);

  MatrixXt xCopy(n, 1);
  for(i = 0; i < n; i++)
    xCopy(i, 0) = x[i] - mean[i];

  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
#ifdef _USE_ATOMICS_IN_DMNORM
  if(CppAD::Value(prec_param) == 1) {
    xCopy = nimDerivs_matmult(mapChol.template triangularView<Eigen::Upper>(), xCopy);
  } else {
    xCopy = nimDerivs_EIGEN_FS(mapChol.transpose(), xCopy);
  }
#else
  if(CppAD::Value(prec_param) == 1) {
    xCopy = mapChol.template triangularView<Eigen::Upper>()*xCopy;
  } else {
    xCopy = mapChol.transpose().template triangularView<Eigen::Lower>().solve(xCopy);
  }
#endif
  /* xCopy = CppAD::CondExpEq(prec_param, Type(1), */
  /*                          mapChol.template triangularView<Eigen::Upper>()*xCopy, */
  /*                          mapChol.template triangularView<Eigen::Upper>().transpose().solve(xCopy) ); */
  // Note that with solve(), transpose of U appears slightly less costly than if input 'chol' were L,
  // presumably because of column-major order interacting well with the solve.
  xCopy = xCopy.array()*xCopy.array();
  dens += -Type(0.5)*(xCopy.sum());
  dens = CppAD::CondExpEq(give_log, Type(1), dens, exp(dens));
  return(dens);
}

template<class Type>
Type nimDerivs_nimArr_dmnorm_chol_logFixed(NimArr<1, Type> &x, NimArr<1, Type> &mean, NimArr<2, Type> &chol, Type prec_param, int give_log, Type overwrite_inputs) {
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;

  if(!CppAD::Constant(prec_param))
    std::cout<<"Warning: In dmnorm, prec_param with value = "<<CppAD::Value(prec_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make prec_param a constant to avoid this warning."<<std::endl;
  
  int n = x.dimSize(0);
  int i;
  Type dens = Type(-n * M_LN_SQRT_2PI);
  Type sumDens = Type(0.);
  for(i = 0; i < n*n; i += n + 1) 
	  sumDens += log(chol[i]);

  dens += CppAD::CondExpEq(prec_param, Type(1), sumDens, -sumDens);

  MatrixXt xCopy(n, 1);
  for(i = 0; i < n; i++)
    xCopy(i, 0) = x[i] - mean[i];

  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
#ifdef _USE_ATOMICS_IN_DMNORM
  if(CppAD::Value(prec_param) == 1) {
    xCopy = nimDerivs_matmult(mapChol.template triangularView<Eigen::Upper>(), xCopy);
  } else {
    xCopy = nimDerivs_EIGEN_FS(mapChol.transpose(), xCopy);
  }
#else
  if(CppAD::Value(prec_param) == 1) {
    xCopy = mapChol.template triangularView<Eigen::Upper>()*xCopy;
  } else {
    xCopy = mapChol.transpose().template triangularView<Eigen::Lower>().solve(xCopy);
    //    xCopy = mapChol.template triangularView<Eigen::Upper>().transpose().solve(xCopy);
  }
#endif
  /*  */
  /* xCopy = CppAD::CondExpEq(prec_param, Type(1), */
  /*                          mapChol.template triangularView<Eigen::Upper>()*xCopy, */
  /*                          mapChol.template triangularView<Eigen::Upper>().transpose().solve(xCopy)); */
  xCopy = xCopy.array()*xCopy.array();
  dens += -Type(0.5)*(xCopy.sum());
  if(!give_log){
    dens = exp(dens);
  }
  return(dens);
}

/* dmvt: Multivariate t distribution */
template<class Type>
Type nimDerivs_nimArr_dmvt_chol(NimArr<1, Type> &x, NimArr<1, Type> &mu, NimArr<2, Type> &chol, Type df, Type prec_param, Type give_log, Type overwrite_inputs) { 
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;

  if(!CppAD::Constant(prec_param))
    std::cout<<"Warning: In dmvt, prec_param with value = "<<CppAD::Value(prec_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make prec_param a constant to avoid this warning."<<std::endl;
  
  int n = x.dimSize(0);
  int i;
  Type dens = nimDerivs_lgammafn((df + Type(n)) / Type(2)) - nimDerivs_lgammafn(df / Type(2)) - Type(n) * Type(M_LN_SQRT_PI) - Type(n) * log(df) / Type(2); 
  Type sumDens = Type(0.);
  for(i = 0; i < n*n; i += n + 1) 
	  sumDens += log(chol[i]);

  dens += CppAD::CondExpEq(prec_param, Type(1), sumDens, -sumDens);

  MatrixXt xCopy(n, 1);
  for(i = 0; i < n; i++)
    xCopy(i, 0) = x[i] - mu[i];

  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  // Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
  if(CppAD::Value(prec_param) == 1) {
    xCopy = mapChol.template triangularView<Eigen::Upper>()*xCopy;
  } else {
    xCopy = mapChol.template triangularView<Eigen::Upper>().transpose().solve(xCopy);
  }
  /* xCopy = CppAD::CondExpEq(prec_param, Type(1), */
  /*                          mapChol.template triangularView<Eigen::Upper>()*xCopy, */
  /*                          mapChol.template triangularView<Eigen::Upper>().transpose().solve(xCopy)); */
  xCopy = xCopy.array()*xCopy.array();
  
  dens += -Type(0.5)*(df + Type(n)) * log(Type(1.) + xCopy.sum() / df);
  
  dens = CppAD::CondExpEq(give_log, Type(1), dens, exp(dens));
  return(dens);
}

template<class Type>
Type nimDerivs_nimArr_dmvt_chol_logFixed(NimArr<1, Type> &x, NimArr<1, Type> &mu, NimArr<2, Type> &chol, Type df, Type prec_param, int give_log, Type overwrite_inputs) { 
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;

  if(!CppAD::Constant(prec_param))
    std::cout<<"Warning: In dmvt, prec_param with value = "<<CppAD::Value(prec_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make prec_param a constant to avoid this warning."<<std::endl;

  int n = x.dimSize(0);
  int i;
  Type dens = nimDerivs_lgammafn((df + Type(n)) / Type(2.)) - nimDerivs_lgammafn(df / Type(2.)) - Type(n) * Type(M_LN_SQRT_PI) - Type(n) * log(df) / Type(2.); 
  Type sumDens = Type(0.);
  for(i = 0; i < n*n; i += n + 1) 
	  sumDens += log(chol[i]);

  dens += CppAD::CondExpEq(prec_param, Type(1), sumDens, -sumDens);

  MatrixXt xCopy(n, 1);
  for(i = 0; i < n; i++)
    xCopy(i, 0) = x[i] - mu[i];

  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  // Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
  if(CppAD::Value(prec_param) == 1) {
    xCopy = mapChol.template triangularView<Eigen::Upper>()*xCopy;
  } else {
    xCopy = mapChol.template triangularView<Eigen::Upper>().transpose().solve(xCopy);
  }

  /* xCopy = CppAD::CondExpEq(prec_param, Type(1), */
  /*                          mapChol.template triangularView<Eigen::Upper>()*xCopy, */
  /*                          mapChol.template triangularView<Eigen::Upper>().transpose().solve(xCopy)); */
  xCopy = xCopy.array()*xCopy.array();
  
  dens += -Type(0.5)*(df + Type(n)) * log(Type(1.) + xCopy.sum() / df);
  
  if(!give_log){
    dens = exp(dens);
  }
  return(dens);
}

/* dlkj: LKJ correlation cholesky factor */
template<class Type>
Type nimDerivs_nimArr_dlkj_corr_cholesky(NimArr<2, Type> &x, Type eta, Type p, Type give_log) {
  // Should p be of type Type or double? If double/int no matching call as passed a CppAD::AD<double>
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  typedef Eigen::Array<Type, Eigen::Dynamic, 1> ArrayXt;
  int n = x.dimSize(0);   // otherwise not sure how to cast `p` to be acceptable to mapX

  NimArr<2, Type> x_possible_copy;
  Eigen::Map<MatrixXt > mapX(nimArrCopyIfNeeded<2, Type>(x, x_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapX(x.getPtr(), n, n);

  // not sure how to use Eigen:seq in arithmetic or to initialize an array
  ArrayXt pseq(n);
  for(int i = 0; i < n; i++)
    pseq(i) = i+1;

  /*
  Type dens = Type(0.);
  int counter = 0;
  for(int k = 2; k <= n; k++) {
    counter += p+1;
    dens += log(x[counter]) * (Type(p) - k + Type(2.0)*eta - Type(2.0));
  }
  */
  Type dens = (mapX.diagonal().array().log() *
               (Type(p) - pseq + Type(2.0)*eta - Type(2.0))).sum();

  dens = CppAD::CondExpEq(give_log, Type(1), dens, exp(dens));
  return(dens);
}

/*  Type sumDens = Type(0.);
  for(i = 0; i < n*n; i += n + 1) 
  sumDens += log(chol[i]);*/


template<class Type>
Type nimDerivs_nimArr_dlkj_corr_cholesky_logFixed(NimArr<2, Type> &x, Type eta, Type p, int give_log) {

  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  typedef Eigen::Array<Type, Eigen::Dynamic, 1> ArrayXt;
  int n = x.dimSize(0);   // otherwise not sure how to cast `p` to be acceptable to mapX
  NimArr<2, Type> x_possible_copy;
  Eigen::Map<MatrixXt > mapX(nimArrCopyIfNeeded<2, Type>(x, x_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapX(x.getPtr(), n, n);

  // not sure how to use Eigen:seq in arithmetic
  ArrayXt pseq(n);
  for(int i = 0; i < n; i++)
    pseq(i) = i+1;

  Type dens = (mapX.diagonal().array().log() *
               (Type(p) - pseq + Type(2.0)*eta - Type(2.0))).sum();

  if(!give_log){
    dens = exp(dens);
  }
  return(dens);
}
  
/* dwish: Wishart distribution */
template<class Type>
Type nimDerivs_nimArr_dwish_chol(NimArr<2, Type> &x, NimArr<2, Type> &chol, Type df, Type scale_param, Type give_log, Type overwrite_inputs) { 

  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  if(!CppAD::Constant(scale_param))
    std::cout<<"Warning: In dwish, scale_param with value = "<<CppAD::Value(scale_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make scale_param a constant to avoid this warning."<<std::endl;

  int n = x.dimSize(0);
  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  //Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
  NimArr<2, Type> x_possible_copy;
  Eigen::Map<MatrixXt > mapX(nimArrCopyIfNeeded<2, Type>(x, x_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapX(x.getPtr(), n, n);
  int p = x.dim()[0];
  Type dens = (df * mapChol.diagonal().array().log()).sum();
  dens = CppAD::CondExpEq(scale_param, Type(1), -dens, dens);

  dens += -(df*p/Type(2.) * Type(M_LN2) + p*(p-Type(1.))*Type(M_LN_SQRT_PI/2.));
  for(int i = 0; i < p; i++)
    dens -= nimDerivs_lgammafn((df - Type(i)) / Type(2.));
 
  /* Lower may be more efficient: https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html#details */
  MatrixXt Lx = mapX.template selfadjointView<Eigen::Lower>().llt().matrixL();
  dens += (df - p - Type(1.)) * Lx.diagonal().array().log().sum();
  
  if(CppAD::Value(scale_param) == 1) {
    dens -= Type(0.5) * (mapChol.template triangularView<Eigen::Upper>().transpose().solve(Lx)).squaredNorm();
  } else {
    dens -= Type(0.5) * (mapChol.template triangularView<Eigen::Upper>()*Lx).squaredNorm();
  }

  /* dens -= Type(0.5) * CppAD::CondExpEq(scale_param, Type(1), */
  /*                          mapChol.template triangularView<Eigen::Upper>().transpose().solve(Lx).squaredNorm(), */
  /*                          (mapChol.template triangularView<Eigen::Upper>()*Lx).squaredNorm()); */

  dens = CppAD::CondExpEq(give_log, Type(1), dens, exp(dens));
  return(dens);
}


template<class Type>
Type nimDerivs_nimArr_dwish_chol_logFixed(NimArr<2, Type> &x, NimArr<2, Type> &chol, Type df, Type scale_param, int give_log, Type overwrite_inputs) { 

  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  if(!CppAD::Constant(scale_param))
    std::cout<<"Warning: In dwish, scale_param with value = "<<CppAD::Value(scale_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make scale_param a constant to avoid this warning."<<std::endl;
  
  int n = x.dimSize(0);
  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  //Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
  NimArr<2, Type> x_possible_copy;
  Eigen::Map<MatrixXt > mapX(nimArrCopyIfNeeded<2, Type>(x, x_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapX(x.getPtr(), n, n);
  int p = x.dim()[0];
  Type dens = (df * mapChol.diagonal().array().log()).sum();
  dens = CppAD::CondExpEq(scale_param, Type(1), -dens, dens);

  dens += -(df*p/Type(2.) * Type(M_LN2) + p*(p-Type(1.))*Type(M_LN_SQRT_PI/2.));
  for(int i = 0; i < p; i++)
    dens -= nimDerivs_lgammafn((df - Type(i)) / Type(2.));
 
  /* Lower may be more efficient: https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html#details */
  MatrixXt Lx = mapX.template selfadjointView<Eigen::Lower>().llt().matrixL();
  dens += (df - p - Type(1.)) * Lx.diagonal().array().log().sum();

  if(CppAD::Value(scale_param) == 1) {
    dens -= Type(0.5) * (mapChol.template triangularView<Eigen::Upper>().transpose().solve(Lx)).squaredNorm();
  } else {
    dens -= Type(0.5) * (mapChol.template triangularView<Eigen::Upper>()*Lx).squaredNorm();
  }

  /* dens -= Type(0.5) * CppAD::CondExpEq(scale_param, Type(1), */
  /*                          mapChol.template triangularView<Eigen::Upper>().transpose().solve(Lx).squaredNorm(), */
  /*                          (mapChol.template triangularView<Eigen::Upper>()*Lx).squaredNorm()); */

  if(!give_log){
    dens = exp(dens);
  }
  return(dens);
}

/* dinvwish: Inverse Wishart distribution */
template<class Type>
Type nimDerivs_nimArr_dinvwish_chol(NimArr<2, Type> &x, NimArr<2, Type> &chol, Type df, Type scale_param, Type give_log, Type overwrite_inputs) { 

  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  if(!CppAD::Constant(scale_param))
    std::cout<<"Warning: In dinvwish, scale_param with value = "<<CppAD::Value(scale_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make scale_param a constant to avoid this warning."<<std::endl;
  int n = x.dimSize(0);
  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  //Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
  NimArr<2, Type> x_possible_copy;
  Eigen::Map<MatrixXt > mapX(nimArrCopyIfNeeded<2, Type>(x, x_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapX(x.getPtr(), n, n);
  int p = x.dim()[0];

  Type dens = (df * mapChol.diagonal().array().log()).sum();
  dens = CppAD::CondExpEq(scale_param, Type(1), dens, -dens);

  dens += -(df*p/Type(2.) * Type(M_LN2) + p*(p-Type(1.))*Type(M_LN_SQRT_PI/2.));
  for(int i = 0; i < p; i++)
    dens -= nimDerivs_lgammafn((df - Type(i)) / Type(2.));
 
  /* Lower may be more efficient: https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html#details */
  MatrixXt Lx = mapX.template selfadjointView<Eigen::Lower>().llt().matrixL();
  dens -= (df + p + Type(1.)) * Lx.diagonal().array().log().sum();

  /* in basic tests, use of .solve(Identity) was 2-3x faster than .inverse();
  I don't see a way to use .inverse that recognizes the triangular input */
  // std::cout<<"Baking in scale_param = "<<scale_param<<" in dinvwish."<<std::endl;
  if(CppAD::Value(scale_param) == 1) {
    dens -= Type(0.5) * (Lx.template triangularView<Eigen::Lower>().solve(mapChol.transpose())).squaredNorm();
  } else {
    dens -= Type(0.5) * (Lx.template triangularView<Eigen::Lower>().solve(mapChol.template triangularView<Eigen::Upper>().solve(MatrixXt::Identity(n,n)))).squaredNorm();
  }

  /* dens -= Type(0.5) * CppAD::CondExpEq(scale_param, Type(1), */
  /*                          (Lx.template triangularView<Eigen::Lower>().solve(mapChol.transpose())).squaredNorm(), */
  /*                          (Lx.template triangularView<Eigen::Lower>().solve(mapChol.template triangularView<Eigen::Upper>().solve(MatrixXd::Identity(n, n)))).squaredNorm()); */

  dens = CppAD::CondExpEq(give_log, Type(1), dens, exp(dens));
  return(dens);
}

template<class Type>
Type nimDerivs_nimArr_dinvwish_chol_logFixed(NimArr<2, Type> &x, NimArr<2, Type> &chol, Type df, Type scale_param, int give_log, Type overwrite_inputs) { 

  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  if(!CppAD::Constant(scale_param))
    std::cout<<"Warning: In dinvwish, scale_param with value = "<<CppAD::Value(scale_param)<<" is a variable but will be fixed in the tape until the tape is reset.  Make scale_param a constant to avoid this warning."<<std::endl;
  int n = x.dimSize(0);
  NimArr<2, Type> chol_possible_copy;
  Eigen::Map<MatrixXt > mapChol(nimArrCopyIfNeeded<2, Type>(chol, chol_possible_copy).getPtr(), n, n);
  //Eigen::Map<MatrixXt > mapChol(chol.getPtr(), n, n);
  NimArr<2, Type> x_possible_copy;
  Eigen::Map<MatrixXt > mapX(nimArrCopyIfNeeded<2, Type>(x, x_possible_copy).getPtr(), n, n);
  //  Eigen::Map<MatrixXt > mapX(x.getPtr(), n, n);
  int p = x.dim()[0];

  Type dens = (df * mapChol.diagonal().array().log()).sum();
  dens = CppAD::CondExpEq(scale_param, Type(1), dens, -dens);

  dens += -(df*p/Type(2.) * Type(M_LN2) + p*(p-Type(1.))*Type(M_LN_SQRT_PI/2.));
  for(int i = 0; i < p; i++)
    dens -= nimDerivs_lgammafn((df - Type(i)) / Type(2.));

  /* Lower may be more efficient: https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html#details */
  MatrixXt Lx = mapX.template selfadjointView<Eigen::Lower>().llt().matrixL();
  dens -= (df + p + Type(1.)) * Lx.diagonal().array().log().sum();

  /* in basic tests, use of .solve(Identity) was 2-3x faster than .inverse();
  I don't see a way to use .inverse that recognizes the triangular input */
  // std::cout<<"Baking in scale_param = "<<scale_param<<" in dinvwish."<<std::endl;

  if(CppAD::Value(scale_param) == 1) {
    dens -= Type(0.5) * (Lx.template triangularView<Eigen::Lower>().solve(mapChol.transpose())).squaredNorm();
  } else {
    dens -= Type(0.5) * (Lx.template triangularView<Eigen::Lower>().solve(mapChol.template triangularView<Eigen::Upper>().solve(MatrixXt::Identity(n,n)))).squaredNorm();
  }

  if(!give_log){
    dens = exp(dens);
  }
  return(dens);

}

/* Multinomial */
template<class Type>
Type nimDerivs_nimArr_dmulti(const NimArr<1, Type> &x,
			     const Type &size,
			     const NimArr<1, Type> &prob,
			     const Type &give_log) 
{
  Type sumProb = 0.0;
  Type sumX = 0.0;

  Type logProb = nimDerivs_lgammafn(size + 1);
  int K = x.dimSize(0);
  for(int i = 0; i < K; i++) {
    sumProb += prob[i];
    sumX += x[i];
  }

  for(int i = 0; i < K; i++) {
    logProb += nimDerivs_log_pow_int(prob[i]/sumProb, x[i]) - nimDerivs_lgammafn(x[i] + Type(1.));
  }
  logProb = CppAD::CondExpEq(sumX, size,
			     logProb,
			     -Type(std::numeric_limits<double>::infinity()));
  Type ans = CppAD::CondExpEq(give_log, Type(1),
			      logProb,
			      exp(logProb));
  return ans;
}

template<class Type>
Type nimDerivs_nimArr_dmulti_logFixed(const NimArr<1, Type> &x,
				      const Type &size,
				      const NimArr<1, Type> &prob,
				      int give_log) 
{
  Type sumProb(0.0);
  Type sumX(0.0);

  Type logProb = nimDerivs_lgammafn(size + 1);
  int K = x.size();
  for(int i = 0; i < K; i++) {
    sumProb += prob[i];
    sumX += x[i];
  }

  for(int i = 0; i < K; i++) {
    logProb += nimDerivs_log_pow_int(prob[i]/sumProb, x[i]) - nimDerivs_lgammafn(x[i] + Type(1.));
  }
  logProb = CppAD::CondExpEq(sumX, size,
			     logProb,
			     -Type(std::numeric_limits<double>::infinity()));
  if(give_log) {
    return logProb;
  } else {
    return exp(logProb);
  }
}
/* Categorial */
template<class Type>
Type nimDerivs_nimArr_dcat(const Type &x,
		    const NimArr<1, Type> &prob,
		    const Type &give_log)
{
  CppAD::VecAD<double> vecProb(prob.size());
  Type sumProb = 0.;
  CppAD::AD<double> i; // Per CppAD example, index assigning into VecAD must be a CppAD::AD<>
  for(i = 0; size_t(Integer(i) ) < prob.size(); i += 1.) {
    vecProb[i] = prob[ size_t(Integer(i) ) ];
  }
  for(int j = 0; j < prob.size(); ++j) {
    sumProb += prob[ j ];
  }
  Type ansProb = vecProb[x - 1] / sumProb;
  Type ans = CppAD::CondExpEq(give_log, Type(1),
			      log(ansProb),
			      ansProb);
  return(ans);
}

template<class Type>
Type nimDerivs_nimArr_dcat_logFixed(const Type &x,
			     const NimArr<1, Type> &prob,
			     int give_log)
{
  CppAD::VecAD<double> vecProb(prob.size());
  Type sumProb = 0.;
  CppAD::AD<double> i;
  for(i = 0; size_t(Integer(i) ) < prob.size(); i += 1.) {
    vecProb[i] = prob[ size_t(Integer(i) ) ];
  }
  for(int j = 0; j < prob.size(); ++j) {
    sumProb += prob[ j ];
  }
  Type ansProb = vecProb[x - 1] / sumProb;
  if(give_log) {
    return log(ansProb);
  } else {
    return ansProb;
  }
}

  // Some of the following functions seem to be missing Type() casting of constants. --CP 2020-04-22
  
/* dt: Student's t distribution */
/* This case is modified from TMB code so that give_log can be different
   when the tape is used from when it is recorded.
   It also uses nimble's atomic for lgammafn*/
template <class Type>
Type nimDerivs_dt(Type x, Type df, Type give_log)
{
  Type res = nimDerivs_lgammafn((df+1)/2) - Type(1)/2*log(df*M_PI) -nimDerivs_lgammafn(df/2) - (df+1)/2*log(1+x*x/df);
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template <class Type>
Type nimDerivs_dt_logFixed(Type x, Type df, int give_log)
{
  Type res = nimDerivs_lgammafn((df+1)/2) - Type(1)/2*log(df*M_PI) -nimDerivs_lgammafn(df/2) - (df+1)/2*log(1+x*x/df);
  if(!give_log) res = exp(res);
  return res;
}


/* Nonstandard t-distribution */
/* This is not in TMB. */
template <class Type>
Type nimDerivs_dt_nonstandard(Type x, Type df, Type mu, Type sigma, Type give_log)
{
  Type res =  nimDerivs_dt_logFixed( (x - mu)/sigma, df, 1) - log(sigma);
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template <class Type>
Type nimDerivs_dt_nonstandard_logFixed(Type x, Type df, Type mu, Type sigma, int give_log)
{
  Type res =  nimDerivs_dt_logFixed( (x - mu)/sigma, df, 1) - log(sigma);
  if(!give_log){
	res = exp(res);
  }
  return(res);
}

/* dgamma: gamma distribution */
/* Code modified from TMB. */
template<class Type>
Type nimDerivs_dgamma(Type y, Type shape, Type scale, Type give_log)
{
  Type res = -nimDerivs_lgammafn(shape)+(shape-Type(1.0))*log(y)-y/scale-shape*log(scale);
	res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
	return(res);
}

template<class Type>
Type nimDerivs_dgamma_logFixed(Type y, Type shape, Type scale, int give_log)
{
  Type res = -nimDerivs_lgammafn(shape)+(shape-Type(1.0))*log(y)-y/scale-shape*log(scale);
  if(!give_log){
    res = exp(res);
  }
  return(res);
}

/* inverse gamma */
/* Code modified from TMB. */
template<class Type>
Type nimDerivs_dinvgamma(Type x, Type shape, Type rate, Type give_log)
{
  Type xinv = Type(1.0)/x;
  Type res = nimDerivs_dgamma_logFixed(xinv, shape, rate, 1) - 2*log(x);
  res = CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_dinvgamma_logFixed(Type x, Type shape, Type rate, int give_log)
{
  Type xinv = Type(1.0)/x;
  Type res = nimDerivs_dgamma_logFixed(xinv, shape, rate, 1) - 2*log(x);
  if(!give_log){
    res = exp(res); //nimDerivs_dgamma_logFixed(xinv, shape, rate, give_log) * xinv * xinv;
  }
  return(res);
}

template<class Type>
Type nimDerivs_dsqrtinvgamma(Type x, Type shape, Type rate, Type give_log)
{
  Type res = nimDerivs_dinvgamma_logFixed(x*x, shape, rate, 1) + log(2*x);
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_dsqrtinvgamma_logFixed(Type x, Type shape, Type rate, int give_log)
{
  Type res = nimDerivs_dinvgamma_logFixed(x*x, shape, rate, 1) + log(2*x);
  if(!give_log){
    res = exp(res);
  }
  return(res);
}


/* Negative binomial */
/* Code modified from TMB */
template<class Type>
inline Type nimDerivs_dnbinom(const Type &x, const Type &size, const Type &prob,
		    Type give_log)
{
  Type n=size;
  Type p=prob;
  // p > 0, so  n*log(p) is ok.
  // but x is integer and p may be == 1, so nimDerivs_log_pow_int(Type(1.)-p, x)
  //    is needed instead of x * log(1.-p)
  Type res = nimDerivs_lgammafn(x+n)-nimDerivs_lgammafn(n)-nimDerivs_lgammafn(x+Type(1))+
    n*log(p) + nimDerivs_log_pow_int(Type(1.)-p, x); 
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
inline Type nimDerivs_dnbinom_logFixed(const Type &x, const Type &size, const Type &prob,
		    int give_log)
{
  Type n=size;
  Type p=prob;
  Type res = nimDerivs_lgammafn(x+n)-nimDerivs_lgammafn(n)-nimDerivs_lgammafn(x+Type(1))+
    n*log(p) + nimDerivs_log_pow_int(Type(1.)-p, x); 
  if(!give_log){
    res = exp(res);
  }
	return(res);
}

/* Poisson */
/* Code modified from TMB */
template<class Type>
inline Type nimDerivs_dpois(const Type &x, const Type &lambda, Type give_log)
{
  Type res = -lambda + nimDerivs_log_pow_int(lambda, x) - nimDerivs_lgammafn(x+Type(1));
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
inline Type nimDerivs_dpois_logFixed(const Type &x, const Type &lambda, int give_log)
{
  Type res = -lambda +  nimDerivs_log_pow_int(lambda, x) - nimDerivs_lgammafn(x+Type(1));
  if(!give_log){
    res = exp(res);
  }
  return(res);
}

/* dexp: exponential distribution */
/* Parameterization available in nimble. Code modified from TMB*/
template<class Type> 
Type nimDerivs_dexp_nimble(Type x, Type rate, Type give_log)
{
  Type res = log(rate)-rate*x;
  res = CppAD::CondExpEq(give_log, Type(1),
			 res,
			 exp(res));
	/* Type res = CppAD::CondExpEq(give_log, Type(0), */
	/* 			    rate*exp(-rate*x), */
	/* 			    log(rate)-rate*x); */
	return(res);
}

template<class Type> 
Type nimDerivs_dexp_nimble_logFixed(Type x, Type rate, int give_log)
{
	Type res;
	if(!give_log){
		res = rate*exp(-rate*x);
	} 
	else{
		res = log(rate)-rate*x;
	}
	return(res);
}

/* Parameterization available in R and TMB */
/* This case is modified from TMB code so that give_log can be different
   when the tape is used from when it is recorded. */
template<class Type>
Type nimDerivs_dexp(Type x, Type scale, Type give_log)
{
  Type res = -log(scale)-x/scale;
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
    //	Type res = CppAD::CondExpEq(give_log, Type(0), exp(-x/scale)/scale, -log(scale)-x/scale);
	return(res);
}

/* This case is modified from TMB code so that we do not use CondExpGe for x > 0.*/
template<class Type>
Type nimDerivs_dexp_logFixed(Type x, Type scale, int give_log)
{
	Type res;
	if(!give_log){
		res = exp(-x/scale)/scale;
	}
	else{
		res = -log(scale)-x/scale;
	}
	return(res);
}

template<class Type>
Type nimDerivs_dflat(Type x, Type give_log)
{
  Type res = CppAD::CondExpEq(give_log, Type(1), Type(0), Type(1));
  return(res);
}
template<class Type>
Type nimDerivs_dflat_logFixed(Type x, int give_log)
{
  if(give_log) return(Type(0));
  return(Type(1));
}
template<class Type>
Type nimDerivs_dhalfflat(Type x, Type give_log)
{
  Type res = CppAD::CondExpGe(x, Type(0), Type(0), -CppAD::numeric_limits<Type>::max());
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}
template<class Type>
Type nimDerivs_dhalfflat_logFixed(Type x, int give_log)
{
  Type res;
  if(give_log) {
    res = CppAD::CondExpGe(x, Type(0), Type(0), -CppAD::numeric_limits<Type>::max());
  } else {
    res = CppAD::CondExpGe(x, Type(0), Type(1), Type(0));
  }
  return(res);
}

/* dunif: uniform distribution */
/* TMB does not have these cases. */
template<class Type> 
Type nimDerivs_dunif(Type x, Type a, Type b, Type give_log)
{   
  Type res = Type(1.0)/(b-a);
  res = CppAD::CondExpEq(give_log, Type(1), log(res), res);
  Type ans = CppAD::CondExpLe(x, b, res, -CppAD::numeric_limits<Type>::max());
  ans = CppAD::CondExpGe(x, a, ans, -CppAD::numeric_limits<Type>::max());
  return(ans);
}

template<class Type> 
Type nimDerivs_dunif_logFixed(Type x, Type a, Type b, int give_log)
{   
  Type res = Type(1.)/(b-a);
  if(give_log){
    res = log(res);
  }
  Type ans = CppAD::CondExpLe(x, b, res, -CppAD::numeric_limits<Type>::max());
  ans = CppAD::CondExpGe(x, a, ans, -CppAD::numeric_limits<Type>::max());
  return(ans);
}

/* dweibull: Weibull distribution */
/* This case is modified from TMB code so that give_log can be different
   when the tape is used from when it is recorded. */
/* Also the conditional check on x > 0 has been removed. */
template<class Type> 
Type nimDerivs_dweibull(Type x, Type shape, Type scale, Type give_log)
{
	Type res = shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape));
	res = CppAD::CondExpEq(give_log, Type(0), res, log(res));
	return(res);
}

/* This case is modified from TMB to remove the conditional check on x > 0. */
template<class Type> 
Type nimDerivs_dweibull_logFixed(Type x, Type shape, Type scale, int give_log)
{
	Type res = shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape));
	if(give_log){
		res = log(res);
	}
	return(res);
}

/* dbinom: Binomial distribution */
/* This case is modified from TMB code so that give_log can be different
   when the tape is used from when it is recorded. */
/* Both cases are modified from TMB to avoid the use of CondExpGt to decide
   which terms to include.  We have observed CondExpGt to give slow performance,
   but we have not benchmarked these cases specifically.*/
template<class Type>
Type nimDerivs_dbinom(Type k, Type size, Type prob, Type give_log)
{
  Type res = nimDerivs_lgammafn(size+1)-nimDerivs_lgammafn(k+1)-nimDerivs_lgammafn(size-k+1)+
    CppAD::azmul(k, log(prob))+ CppAD::azmul(size-k, log(1-prob));
  res = CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_dbinom_logFixed(Type k, Type size, Type prob, int give_log)
{
  Type res = nimDerivs_lgammafn(size+1)-nimDerivs_lgammafn(k+1)-nimDerivs_lgammafn(size-k+1)+
    CppAD::azmul(k, log(prob))+ CppAD::azmul(size-k, log(1-prob));
  if(!give_log){
    res = exp(res);
  }
  return(res);
}

/* dchisq: chi-squared distribution */
/* TMB does not appear to have this distribution. */
template<class Type>
Type nimDerivs_dchisq(Type x, Type df, Type give_log)
{	
  Type res = (df/Type(2.0) - Type(1.0))*log(x) - (x/Type(2.0)) - (df/Type(2.0))*log(Type(2.0)) - nimDerivs_lgammafn(df/Type(2.0));
  res = CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_dchisq_logFixed(Type x, Type df, int give_log)
{	
  Type res = (df/Type(2.0) - Type(1.0))*log(x) - (x/Type(2.0)) - (df/Type(2.0))*log(Type(2.0)) - nimDerivs_lgammafn(df/Type(2.0));
  if(!give_log){
	res = exp(res);
  }
  return(res);
}

/* dbeta: beta distribution */
/* These are modified from TMB to skip the check of x==0, but we may want the behavior for that case implemented in TMB. */
template <class Type>
Type nimDerivs_dbeta(Type x, Type shape1, Type shape2, Type give_log)
{
  Type res = nimDerivs_lgammafn(shape1+shape2) - nimDerivs_lgammafn(shape1) - nimDerivs_lgammafn(shape2) + 
    (shape1-1)*log(x) + (shape2-1)*log(1-x);
  res = CondExpEq(give_log, Type(1), res, exp(res));

  /* The following code, when used with tape optimize()ation, gives a crash when the wrong bracnch of CondExpEq is triggered */
  /* Type res = CondExpEq(give_log, Type(0),  */
  /* 		       exp(nimDerivs_lgammafn(shape1+shape2) - nimDerivs_lgammafn(shape1) -  */
  /* 			   nimDerivs_lgammafn(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1), */
  /* 		       nimDerivs_lgammafn(shape1+shape2) - nimDerivs_lgammafn(shape1) - nimDerivs_lgammafn(shape2) +  */
  /* 		       (shape1-1)*log(x) + (shape2-1)*log(1-x)); */
  return(res);
}

template <class Type>
Type nimDerivs_dbeta_logFixed(Type x, Type shape1, Type shape2, int give_log)
{
  Type res;
  if(give_log){
    res = nimDerivs_lgammafn(shape1+shape2) - nimDerivs_lgammafn(shape1) - nimDerivs_lgammafn(shape2) + 
      (shape1-1)*log(x) + (shape2-1)*log(1-x);
  }
  else{
    res = exp(nimDerivs_lgammafn(shape1+shape2) - nimDerivs_lgammafn(shape1) - 
	      nimDerivs_lgammafn(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1);
  } 
  return(res);
}

/* dlogis: logistic distribution */
/* This case is modified from TMB code so that give_log can be different
   when the tape is used from when it is recorded. */
template <class Type>
Type nimDerivs_dlogis(Type x, Type location, Type scale, Type give_log)
{
  Type res = -(x-location)/scale - log(scale) - 2*log(1+exp(-(x-location)/scale));
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template <class Type>
Type nimDerivs_dlogis_logFixed(Type x, Type location, Type scale, int give_log)
{
  Type res = -(x-location)/scale - log(scale) - 2*log(1+exp(-(x-location)/scale));
  if(!give_log) res = exp(res);
  return res;
}

/* dlnorm: log-normal distribtuion */
/* This does not seem to be implemented in TMB. */
template<class Type>
Type nimDerivs_dlnorm(Type x, Type mean, Type sd, Type give_log)
{
  Type res =  -log(x) - log(sd) - Type(0.5)*log(Type(2.0*M_PI)) - Type(.5)*pow((log(x)-mean)/sd, 2); 
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_dlnorm_logFixed(Type x, Type mean, Type sd, int give_log)
{
  Type res =  -log(x) - log(sd) - Type(0.5)*log(Type(2.0*M_PI)) - Type(.5)*pow((log(x)-mean)/sd, 2); 
  if(!give_log){
	  res = exp(res);
  }
  return(res);
}

/* ddexp: double exponential (Lapalace) distribution */
template<class Type>
Type nimDerivs_ddexp(Type x, Type location, Type scale, Type give_log)
{
  Type res = log(Type(0.5)) - log(scale) - abs(x-location)/scale;
  //  Type res = (Type(1.0/2.0) / scale) * exp(-abs(x - location) / scale);
  //  res = CppAD::CondExpEq(give_log, Type(1), log(res), res);
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_ddexp_logFixed(Type x, Type location, Type scale, int give_log)
{
  Type res = log(Type(0.5)) - log(scale) - abs(x-location)/scale;
  //  Type res = Type(1.0/2.0) / scale * exp(-abs(x - location) / scale);
  if(!give_log){
	  res = exp(res);
  }
  return(res);
}

/* ddirch: Dirichlet distribution */
/* This is not in TMB */
template <class Type>
Type nimDerivs_nimArr_ddirch(NimArr<1, Type> &x, NimArr<1, Type> &alpha, Type give_log)
{
  Type K = alpha.size();
  Type n = x.size();
  Type logres = CppAD::CondExpLt(K, Type(1), CppAD::numeric_limits<Type>::quiet_NaN(), Type(0));
  logres = CppAD::CondExpEq(K, n, logres, CppAD::numeric_limits<Type>::quiet_NaN());
  
  Type sumAlpha = Type(0.0);
  Type sumX = Type(0.0);
  Type dens = Type(0.0);
  
  for(int i = 0; i < K; i++) {
    logres += CppAD::CondExpGt(alpha[i], Type(0), Type(0),   CppAD::numeric_limits<Type>::quiet_NaN());
    logres += CppAD::CondExpGe(x[i], Type(0), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += CppAD::CondExpLe(x[i], Type(1), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += (alpha[i]-Type(1)) * log(x[i]) - nimDerivs_lgammafn(alpha[i]) ;
    sumAlpha += alpha[i];
    sumX += x[i];
  }
  logres += CppAD::CondExpLe(sumX, Type(1.0 + 10.0*std::numeric_limits<double>::epsilon()), 
			     CppAD::CondExpGe(sumX, Type(1.0 - 10.0*std::numeric_limits<double>::epsilon()), Type(0), -Type(std::numeric_limits<double>::infinity())),
			     -Type(std::numeric_limits<double>::infinity()));
  logres += nimDerivs_lgammafn(sumAlpha);
  logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
  return(logres);
}

template <class Type>
Type nimDerivs_nimArr_ddirch_logFixed(NimArr<1, Type> &x, NimArr<1, Type> &alpha, int give_log)
{
  Type K = alpha.size();
  Type n = x.size();
  Type logres = CppAD::CondExpLt(K, Type(1), CppAD::numeric_limits<Type>::quiet_NaN(), Type(0));
  logres = CppAD::CondExpEq(K, n, logres, CppAD::numeric_limits<Type>::quiet_NaN());
  
  Type sumAlpha = Type(0.0);
  Type sumX = Type(0.0);
  Type dens = Type(0.0);
  
  for(int i = 0; i < K; i++) {
    logres += CppAD::CondExpGt(alpha[i], Type(0), Type(0), CppAD::numeric_limits<Type>::quiet_NaN());
    logres += CppAD::CondExpGe(x[i], Type(0), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += CppAD::CondExpLe(x[i], Type(1), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += (alpha[i]-Type(1)) * log(x[i]) - nimDerivs_lgammafn(alpha[i]) ;
    sumAlpha += alpha[i];
    sumX += x[i];
  }
  logres += CppAD::CondExpLe(sumX, Type(1.0 + 10.0*std::numeric_limits<double>::epsilon()), 
			     CppAD::CondExpGe(sumX, Type(1.0 - 10.0*std::numeric_limits<double>::epsilon()), Type(0), -Type(std::numeric_limits<double>::infinity())),
			     -Type(std::numeric_limits<double>::infinity()));
  logres += nimDerivs_lgammafn(sumAlpha);
  if(!give_log){
    logres = exp(logres);
  }
  return(logres);
}

/*************/
/* Functions */


/* pow */
/* If x > 0, do simple x^y.*/
/* Otherwise, if y is not an integer, return NaN */
/* The CondExp expressions may slow down CppAD tapes.*/
//template<class Type>
#define Type CppAD::AD<double>
inline Type nimDerivs_pow(Type x, Type y) {

  /* We experimented with CppAD conditionals but still had cases that 
     crashed or had constants baked into tapes when they shouldn't be
     As a result we implemented nimDerivs_pow_int for cases where 
     y will only have integer values and never need its derivatives taken.
  */

  /* Type outVal = CppAD::CondExpGt(x, Type(0), */
  /* 			       CppAD::pow(x, y), */
  /* 			       CppAD::CondExpEq(y, discrete_round(y), */
  /* 						CppAD::pow(x, CppAD::Integer(y)), */
  /* 						Type(CppAD::numeric_limits<Type>::quiet_NaN()))); */

  /* There is still a need to identify hard-coded cases with integer y and 
     divert them to nimDerivs_pow_int.  Currently nimble-generated code
     does not cast the second argument to Type() for pow, so one of the 
     over-loaded cases below (for y being int or double) will be called.*/
  Type outVal = CppAD::pow(x, y);
  return(outVal);
}

// This may not longer ever be used from nimble-generated code.
//template<class Type>
inline Type nimDerivs_pow(Type x, int y) {
  Type outVal = nimDerivs_pow_int(x, y);
  return(outVal);
}

//template<class Type>
inline Type nimDerivs_pow(Type x, double y) {
  Type outVal;
  if(fabs(y - round(y)) < 1e-8) // allows binary rounding error on computed indices
    outVal = nimDerivs_pow_int(x, round(y));
  else
    outVal = nimDerivs_pow(x, Type(y));
  return(outVal);
}

/* probit */
/* This is modified from TMB's qnorm, to avoid using the mean and sd arguments unnecessarily. */
/* template<class T> */
/* T nimDerivs_probit(T p){ */
/*   return nimDerivs_qnorm1<T>(p); */
/* } */

/* /\* iprobit *\/ */
/* template<class T> */
/* T nimDerivs_iprobit(T q){ */
/*   return nimDerivs_pnorm1<T>(q); */
/* } */

/* This seems necessary for use in Eigen unaryExpr via ptr_fun */
/* which is necessary because otherwise Eigen circumvents regular call to fabs. */
#define T CppAD::AD<double>
//template<class T>
inline T nimDerivs_fabs(T x) {
  return fabs(x);
}

/* lfactorial */
//template<class T>
inline T nimDerivs_lfactorial(T x) {
  return nimDerivs_lgammafn(x + 1.);
}

/* factorial */
//template<class Type>
inline Type nimDerivs_factorial(Type x) {
  return nimDerivs_gammafn(x + 1.); // Even R does not restrict to x >= 0
}
#undef Type
#undef T

/* /\* gammafn *\/ */
/* template<class Type> */
/* Type nimDerivs_gammafn(Type x) { */
/*   return exp(nimDerivs_lgammafn(x)); */
/* } */

/* AD recycling rule. See nimbleEigen.h for the non-AD equivalents. */
template<>
struct nimble_eigen_traits<CppAD::AD<double> > {
  enum {nimbleUseLinearAccess = int(1)};
};

template<>
struct nimble_size_impl<CppAD::AD<double> > {
  static unsigned int getSize(const CppAD::AD<double> &arg) {return 1;}
  static unsigned int getRows(const CppAD::AD<double> &arg) {return 1;}
  static unsigned int getCols(const CppAD::AD<double> &arg) {return 1;}
  //  static unsigned int getDiagRows(const CppAD::AD<double> &arg) {return floor(arg);}
  //  static unsigned int getDiagCols(const CppAD::AD<double> &arg) {return floor(arg);}
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_impl<true, result_type, CppAD::AD<double>, Index> {
  static result_type getCoeff(const CppAD::AD<double> &Arg, Index i) {return Arg;}
  //  static result_type getDiagCoeff(const CppAD::AD<double> &Arg, Index i) {return 1.0;}
};

template<typename result_type, typename Index>
struct nimble_eigen_coeff_mod_impl<true, result_type, CppAD::AD<double>, Index> {
  static result_type getCoeff(const CppAD::AD<double> Arg, Index i, unsigned int size) {return Arg;}
};

MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dbinom, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS2_1scalar(nimDerivs_dexp_nimble, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS2_1scalar(nimDerivs_dexp, CppAD::AD<double>) // broken
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dnbinom, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS2_1scalar(nimDerivs_dpois, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS2_1scalar(nimDerivs_dchisq, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dbeta, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dnorm, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dgamma, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dinvgamma, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dsqrtinvgamma, CppAD::AD<double>) // broken
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_ddexp, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dlnorm, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dlogis, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dunif, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS3_1scalar(nimDerivs_dweibull, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS4_1scalar(nimDerivs_dt_nonstandard, CppAD::AD<double>)
MAKE_RECYCLING_RULE_CLASS2_1scalar(nimDerivs_dt, CppAD::AD<double>)

MAKE_RECYCLING_RULE_CLASS1_1scalar(nimDerivs_pow_int, CppAD::AD<double>)

#define nimDerivs_nimNewMatrixD newMatrix_impl<MatrixXd_CppAD, CppAD::AD<double> >::newMatrix
#define nimDerivs_nimDiagonalD diagonal_impl<MatrixXd_CppAD, CppAD::AD<double> >::diagonal
#define nimDerivs_nimCd concatenate_impl<MatrixXd_CppAD, CppAD::AD<double> >::concatenate
#define nimDerivs_nimNonseqIndexedd nonseqIndexed_impl<MatrixXd_CppAD, CppAD::AD<double> >::nonseqIndexed
#define nimDerivs_nimRepd rep_impl<MatrixXd_CppAD, CppAD::AD<double> >::rep

#endif
