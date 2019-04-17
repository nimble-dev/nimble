#ifndef _NIMDERIVS_TMB__
#define _NIMDERIVS_TMB__

// Some of the following is extracted from TMB.hpp.
// We don't need all of what that includes, so we have
// extracted pieces that allow the subset of TMB that we
// use to work.
#define TMB_EXTERN
#define CSKIP(x) x
#define IF_TMB_PRECOMPILE(x)
#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include <TMB/atomic_math.hpp>
// #include <TMB/atomic_macro.hpp> // loaded by atomic_math
#include <TMB/dnorm.hpp>
#include <TMB/lgamma.hpp>
#include <TMB/distributions_R.hpp>

// Functions here connect from nimble-generated C++ to TMB code that uses CppAD.
// TMB provides many nice distribution functions.
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

/* This wraps a call to TMB's dnorm. */
template<class Type>
Type nimDerivs_dnorm_logFixed(Type x, Type mean, Type sd, int give_log)
{
  return dnorm<Type>(x, mean, sd, give_log);
}

/* pnorm: normal distribution */
/* This wraps a call to TMB's pnorm. */ 
template<class Type>
Type nimDerivs_pnorm(Type q, Type mean = 0., Type sd = 1.){
  return pnorm<Type>(q, mean, sd);
}

/* anorm: normal distribution */
/* This wraps a call to TMB's qnorm. */ 
template<class Type>
Type nimDerivs_qnorm(Type p, Type mean = 0., Type sd = 1.){
  return qnorm<Type>(p, mean, sd);
}

/* dmnorm: Multivariate normal distribution */
/* TMB uses specialized atomic classes.  nimble currently has a "naive" and likely less efficiency implementation: */
/* TO-DO: support by prec_param = 0 and 1*/
template<class Type>
Type nimDerivs_nimArr_dmnorm_chol(NimArr<1, Type> &x, NimArr<1, Type> &mean, NimArr<2, Type> &chol, Type prec_param, Type give_log, Type overwrite_inputs) { 
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  int n = x.dimSize(0);
  int i;
  Type dens = Type(-n * M_LN_SQRT_2PI);
  Type sumDens = Type(0);
  for(i = 0; i < n*n; i += n + 1) 
	  sumDens += log(chol[i]);

  dens += CppAD::CondExpEq(prec_param, Type(1), sumDens, -sumDens);

  MatrixXt xCopy(n, 1);
  for(i = 0; i < n; i++)
    xCopy(i, 0) = x[i] - mean[i];

  Eigen::Map<MatrixXt > eigenChol(chol.getPtr(), n, n);
  xCopy = eigenChol*xCopy;
  xCopy = xCopy.array()*xCopy.array();
  dens += -Type(0.5)*xCopy.sum();
  dens = CppAD::CondExpEq(give_log, Type(1), dens, exp(dens));
  return(dens);
}

template<class Type>
Type nimDerivs_nimArr_dmnorm_chol_logFixed(NimArr<1, Type> &x, NimArr<1, Type> &mean, NimArr<2, Type> &chol, Type prec_param, int give_log, Type overwrite_inputs) { 
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  int n = x.dimSize(0);
  int i;
  Type dens = Type(-n * M_LN_SQRT_2PI);
  Type sumDens = Type(0);
  for(i = 0; i < n*n; i += n + 1) 
	  sumDens += log(chol[i]);

  dens += CppAD::CondExpEq(prec_param, Type(1), sumDens, -sumDens);

  MatrixXt xCopy(n, 1);
  for(i = 0; i < n; i++)
    xCopy(i, 0) = x[i] - mean[i];

  Eigen::Map<MatrixXt > eigenChol(chol.getPtr(), n, n);
  xCopy = eigenChol*xCopy;
  xCopy = xCopy.array()*xCopy.array();
  dens += -Type(0.5)*xCopy.sum();
  if(!give_log){
    dens = exp(dens);
  }
  return(dens);
}

/* dt: Student's t distribution */
/* This case is modified from TMB code so that give_log can be different
   when the tape is used from when it is recorded. */
template <class Type>
Type nimDerivs_dt(Type x, Type df, Type give_log)
{
  Type res = lgamma((df+1)/2) - Type(1)/2*log(df*M_PI) -lgamma(df/2) - (df+1)/2*log(1+x*x/df);
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

/* This wraps a call to TMB's dt. */
template <class Type>
Type nimDerivs_dt_logFixed(Type x, Type df, int give_log)
{
  return dt<Type>(x, df, give_log);
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


/* dexp: exponential distribution */
/* Parameterization available in nimble. Code modified from TMB*/
template<class Type> 
Type nimDerivs_dexp_nimble(Type x, Type rate, Type give_log)
{
	Type res = CppAD::CondExpEq(give_log, Type(0),
				    rate*exp(-rate*x),
				    log(rate)-rate*x);
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
	Type res = CppAD::CondExpEq(give_log, Type(0), exp(-x/scale)/scale, -log(scale)-x/scale);
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

/* dunif: uniform distribution */
/* TMB does not have these cases. */
template<class Type> 
Type nimDerivs_dunif(Type x, Type a, Type b, Type give_log)
{   
	Type res = 1/(b-a);
	res = CppAD::CondExpEq(give_log, Type(1), log(res), res);
	return(res);
}

template<class Type> 
Type nimDerivs_dunif_logFixed(Type x, Type a, Type b, int give_log)
{   
	Type res = 1/(b-a);
	if(give_log){
		res = log(res);
	}
	return(res);
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
  Type res = lgamma(size+1)-lgamma(k+1)-lgamma(size-k+1)+k*log(prob)+(size-k)*log(1-prob);
  res = CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template<class Type>
Type nimDerivs_dbinom_logFixed(Type k, Type size, Type prob, int give_log)
{
  Type res = lgamma(size+1)-lgamma(k+1)-lgamma(size-k+1)+k*log(prob)+(size-k)*log(1-prob);
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
  Type res = (df/Type(2.0) - Type(1.0))*log(x) - (x/Type(2.0)) - (df/Type(2.0))*log(Type(2.0)) - lgamma(df/Type(2.0));
  res = CondExpEq(give_log, Type(0), exp(res), res);
  return(res);
}

template<class Type>
Type nimDerivs_dchisq_logFixed(Type x, Type df, int give_log)
{	
  Type res = (df/Type(2.0) - Type(1.0))*log(x) - (x/Type(2.0)) - (df/Type(2.0))*log(Type(2.0)) - lgamma(df/Type(2.0));
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
  Type res = CondExpEq(give_log, Type(0), 
		       exp(lgamma(shape1+shape2) - lgamma(shape1) - 
			   lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1),
		       lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2) + 
		       (shape1-1)*log(x) + (shape2-1)*log(1-x));
  return(res);
}

template <class Type>
Type nimDerivs_dbeta_logFixed(Type x, Type shape1, Type shape2, int give_log)
{
  Type res;
  if(give_log){
    res = lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2) + 
      (shape1-1)*log(x) + (shape2-1)*log(1-x);
  }
  else{
    res = exp(lgamma(shape1+shape2) - lgamma(shape1) - 
	      lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1);
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
  return dlogis<Type>(x, location, scale, give_log);
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


/* ddirch: Dirichlet distribution */
/* This is not in TMB */
template <class Type>
Type nimDerivs_nimArr_ddirch(NimArr<1, Type> &x, NimArr<1, Type> &alpha, Type give_log)
{
  Type K = alpha.size();
  Type n = x.size();
  Type logres = CppAD::CondExpLt(K, Type(1), Type(CppAD::numeric_limits<Type>::quiet_NaN()), Type(0));
  logres = CppAD::CondExpEq(K, n, logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  
  Type sumAlpha = Type(0.0);
  Type sumX = Type(0.0);
  Type dens = Type(0.0);
  
  for(int i = 0; i < K; i++) {
    logres += CppAD::CondExpGt(alpha[i], Type(0), Type(0),   Type(CppAD::numeric_limits<Type>::quiet_NaN()));
    logres += CppAD::CondExpGe(x[i], Type(0), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += CppAD::CondExpLe(x[i], Type(1), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += (alpha[i]-Type(1)) * log(x[i]) - lgamma(alpha[i]) ;
    sumAlpha += alpha[i];
    sumX += x[i];
  }
  logres += CppAD::CondExpLe(sumX, Type(1.0 + 100.0*numeric_limits<Type>::epsilon()), 
			     CppAD::CondExpGe(sumX, Type(1.0 - 100.0*numeric_limits<Type>::epsilon()), Type(0), -Type(std::numeric_limits<double>::infinity())),
			     -Type(std::numeric_limits<double>::infinity()));
  logres += lgamma(sumAlpha);
  logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
  return(logres);
}

template <class Type>
Type nimDerivs_nimArr_ddirch_logFixed(NimArr<1, Type> &x, NimArr<1, Type> &alpha, int give_log)
{
  Type K = alpha.size();
  Type n = x.size();
  Type logres = CppAD::CondExpLt(K, Type(1), Type(CppAD::numeric_limits<Type>::quiet_NaN()), Type(0));
  logres = CppAD::CondExpEq(K, n, logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  
  Type sumAlpha = Type(0.0);
  Type sumX = Type(0.0);
  Type dens = Type(0.0);
  
  for(int i = 0; i < K; i++) {
    logres += CppAD::CondExpGt(alpha[i], Type(0), Type(0),   Type(CppAD::numeric_limits<Type>::quiet_NaN()));
    logres += CppAD::CondExpGe(x[i], Type(0), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += CppAD::CondExpLe(x[i], Type(1), Type(0), -Type(std::numeric_limits<double>::infinity()));
    logres += (alpha[i]-Type(1)) * log(x[i]) - lgamma(alpha[i]) ;
    sumAlpha += alpha[i];
    sumX += x[i];
  }
  logres += CppAD::CondExpLe(sumX, Type(1.0 + 100.0*numeric_limits<Type>::epsilon()), 
			     CppAD::CondExpGe(sumX, Type(1.0 - 100.0*numeric_limits<Type>::epsilon()), Type(0), -Type(std::numeric_limits<double>::infinity())),
			     -Type(std::numeric_limits<double>::infinity()));
  logres += lgamma(sumAlpha);
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
template<class Type> 
Type nimDerivs_pow(Type x, Type y) {
	Type outVal = CppAD::CondExpGt(x, Type(0),
				       CppAD::pow(x, y),
				       CppAD::CondExpEq(y, discrete_round(y),
							CppAD::pow(x, CppAD::Integer(y)),
							Type(CppAD::numeric_limits<Type>::quiet_NaN())));
	return(outVal);
}

template<class Type> 
Type nimDerivs_pow(int x, Type y) {
	Type outVal = CppAD::CondExpGt(x, Type(0),
				       CppAD::pow(x, y),
				       CppAD::CondExpEq(y, discrete_round(y),
							CppAD::pow(x, CppAD::Integer(y)),
							Type(CppAD::numeric_limits<Type>::quiet_NaN())));
	return(outVal);
}

template<class Type> 
Type nimDerivs_pow(Type x, int y) {
	Type outVal = CppAD::pow(x, y);
	return(outVal);
}


/* probit */
/* This is modified from TMB's qnorm, to avoid using the mean and sd arguments unnecessarily. */
template<class T>
T nimDerivs_probit(T p){
  CppAD::vector<T> tx(1);
  tx[0] = p;
  return tmb_atomic::qnorm1(tx)[0];
}

/* iprobit */
template<class T>
T nimDerivs_iprobit(T q){
  CppAD::vector<T> tx(1);
  tx[0] = q;
  return tmb_atomic::pnorm1(tx)[0];
}

/* lfactorial */
template<class Type> 
Type nimDerivs_lfactorial(Type x) {
  return lfactorial<Type>(x);
} 

/* factorial */
template<class Type> 
Type nimDerivs_factorial(Type x) {
  return exp(lfactorial<Type>(x));
}

#endif
