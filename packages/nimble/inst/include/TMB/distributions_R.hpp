// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

// Copyright (C) 2018 the NIMBLE authors 
// License: GPL (>=2)

/**	\file
	\brief Probability distribution functions.
	*/
	#ifdef WITH_LIBTMB
	#define CSKIP(x) ;
	#define TMB_EXTERN extern
	#else
	#define CSKIP(x) x
	#define TMB_EXTERN
	#endif
	#define IF_TMB_PRECOMPILE(x)
	#include <Eigen/Dense>
	#include <Eigen/Sparse>
    #include <TMB/lgamma.hpp>
  	#include <nimble/NimArr.h>

/** \brief Distribution function of the normal distribution (following R argument convention).
    \ingroup R_style_distribution
*/
#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

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

template<class Type> 
Type nimDerivs_pow(Type x, Type y) {
	Type outVal = CppAD::CondExpGt(x, Type(0), CppAD::pow(x, y), CppAD::CondExpEq(y, discrete_round(y), CppAD::pow(x, CppAD::Integer(y)), Type(CppAD::numeric_limits<Type>::quiet_NaN())));
	return(outVal);
}

template<class Type> 
Type nimDerivs_pow(int x, Type y) {
	Type outVal = CppAD::CondExpGt(x, Type(0), CppAD::pow(x, y), CppAD::CondExpEq(y, discrete_round(y), CppAD::pow(x, CppAD::Integer(y)), Type(CppAD::numeric_limits<Type>::quiet_NaN())));
	return(outVal);
}

template<class Type> 
Type nimDerivs_pow(Type x, int y) {
	Type outVal = CppAD::pow(x, y);
	return(outVal);
}

template<class T>
T nimDerivs_gammafn(T x) {
  return(exp(lgamma(x)));
}

template<class T>
T nimDerivs_lgammafn(T x) {
  return(lgamma(x));
}

template<class T>
T nimDerivs_atan(T x) {
  return(CppAD::atan(x));
}

template<class T>
T nimDerivs_cosh(T x) {
  return(CppAD::cosh(x));
}

template<class T>
T nimDerivs_sinh(T x) {
  return(CppAD::sinh(x));
}

template<class T>
T nimDerivs_tanh(T x) {
  return(CppAD::tanh(x));
}

template<class T>
T nimDerivs_acosh(T x) {
  return(CppAD::acosh(x));
}

template<class T>
T nimDerivs_asinh(T x) {
  return(CppAD::asinh(x));
}

template<class T>
T nimDerivs_atanh(T x) {
  return(CppAD::atanh(x));
}

template<class T>
T nimDerivs_log1p(T x) {
  return(CppAD::log1p(x));
}

// template<class Type>
// Type nimDerivs_nimArr_dwish_chol(NimArr<2, Type> &xNimArr, NimArr<2, Type> &cholNimArr,
// 	 Type df, Type scale_param, Type give_log, Type overwrite_inputs){
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
//   int p = xNimArr.dim()[0];
//   int i, j;

//   Type dens = -(df*p/2 * Type(M_LN2) + p*(p-1)*M_LN_SQRT_PI/2);
//   for(i = 0; i < p; i++)
//     dens -= lgamma((df - Type(i)) / Type(2.0));

//   Type sumLogCholNimArr;
//   for(i = 0; i < p; i++){ 
// 	sumLogCholNimArr += df * log(cholNimArr(i, i));
//   }
//   dens += CppAD::CondExpEq(scale_param, Type(1), -sumLogCholNimArr, sumLogCholNimArr);

// //   determinant of x using Cholesky:
//   Eigen::Map<MatrixXt > eigenX(xNimArr.getPtr(), p, p);

//    Eigen::LLT<MatrixXt > llt(eigenX);
//    MatrixXt    eigenXChol = llt.matrixL();
 
// eigenXChol.template triangularView<Eigen::Upper>() =
// 		eigenXChol.transpose().template triangularView<Eigen::Upper>();

// for(i = 0; i < p; i++){ 
//   	dens += (df - p - 1) *  log(eigenXChol(i, i));
// }

//  Eigen::Map<MatrixXt > eigenChol(cholNimArr.getPtr(), p, p); // may need to create new eigen matrix instead of mapping here
//  Type tmp_dens = 0.0;
//  MatrixXt eigenSolved = eigenChol.transpose().lu().solve(eigenXChol.transpose());
//  MatrixXt cholMultX = eigenChol*eigenX;
// 	for(j = 0; j < p; j++){ 
// 		for(i = 0; i <= j; i++){ 
// 			 tmp_dens +=  CppAD::CondExpEq(scale_param, Type(1), eigenSolved(j,i)*eigenSolved(j,i), 
// 			   cholMultX(i, j) * eigenChol(i, j));
// 		}
// 	}
//   dens += -0.5 * tmp_dens;
//   return CppAD::CondExpEq(give_log, Type(1.0), dens, exp(dens));
// }

template<class Type>
Type nimDerivs_pnorm(Type q, Type mean = 0., Type sd = 1.){
  CppAD::vector<Type> tx(1);
  tx[0] = (q - mean) / sd;
  return ::atomic::pnorm1(tx)[0];
}

// template<class Type>
// Type nimDerivs_nimArr_dwish_chol_logFixed(NimArr<2, Type> &xNimArr, NimArr<2, Type> &cholNimArr,
// 	 Type df, Type scale_param, int give_log, Type overwrite_inputs){
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
//   int p = xNimArr.dim()[0];
//   int i, j;

//   Type dens = -(df*p/2 * Type(M_LN2) + p*(p-1)*M_LN_SQRT_PI/2);
//   for(i = 0; i < p; i++)
//     dens -= lgamma((df - Type(i)) / Type(2.0));

//   Type sumLogCholNimArr;
//   for(i = 0; i < p; i++){ 
// 	sumLogCholNimArr += df * log(cholNimArr(i, i));
//   }
//   dens += CppAD::CondExpEq(scale_param, Type(1), -sumLogCholNimArr, sumLogCholNimArr);

// //   determinant of x using Cholesky:
//   Eigen::Map<MatrixXt > eigenX(xNimArr.getPtr(), p, p);

//    Eigen::LLT<MatrixXt > llt(eigenX);
//    MatrixXt    eigenXChol = llt.matrixL();
 
// eigenXChol.template triangularView<Eigen::Upper>() =
// 		eigenXChol.transpose().template triangularView<Eigen::Upper>();

// for(i = 0; i < p; i++){ 
//   	dens += (df - p - 1) *  log(eigenXChol(i, i));
// }

//  Eigen::Map<MatrixXt > eigenChol(cholNimArr.getPtr(), p, p); // may need to create new eigen matrix instead of mapping here
//  Type tmp_dens = 0.0;
//  MatrixXt eigenSolved = eigenChol.transpose().lu().solve(eigenXChol.transpose());
//  MatrixXt cholMultX = eigenChol*eigenX;
// 	for(j = 0; j < p; j++){ 
// 		for(i = 0; i <= j; i++){ 
// 			 tmp_dens +=  CppAD::CondExpEq(scale_param, Type(1), eigenSolved(j,i)*eigenSolved(j,i), 
// 			   cholMultX(i, j) * eigenChol(i, j));
// 		}
// 	}
//   dens += -0.5 * tmp_dens;
//   if(!give_log){
//     dens = exp(dens);
//   }
//   return dens;
// }

/** \brief Quantile function of the normal distribution (following R argument convention).
    \ingroup R_style_distribution
*/
template<class Type>
Type nimDerivs_qnorm(Type p, Type mean = 0., Type sd = 1.){
  CppAD::vector<Type> tx(1);
  tx[0] = p;
  return sd*::atomic::qnorm1(tx)[0] + mean;
}

template<class T>
T nimDerivs_probit(T x){
  return(nimDerivs_qnorm(x, 0., 1., 1, 0));
}

template<class T>
T nimDerivs_iprobit(T x){
  return(nimDerivs_pnorm(x, 0., 1., 1, 0));
}
// VECTORIZE3_ttt(qnorm)
// VECTORIZE1_t(qnorm)

/** \brief Distribution function of the gamma distribution (following R argument convention).
    \ingroup R_style_distribution
*/
// template<class Type>
// Type pgamma(Type q, Type shape, Type scale = 1.){
//   CppAD::vector<Type> tx(4);
//   tx[0] = q/scale;
//   tx[1] = shape;
//   tx[2] = Type(0);        // 0'order deriv
//   tx[3] = -lgamma(shape); // normalize
//   return atomic::D_incpl_gamma_shape(tx)[0];
// }
// VECTORIZE3_ttt(pgamma)

// /** \brief Quantile function of the gamma distribution (following R argument convention).
//     \ingroup R_style_distribution
// */
// template<class Type>
// Type qgamma(Type q, Type shape, Type scale = 1.){
//   CppAD::vector<Type> tx(3);
//   tx[0] = q;
//   tx[1] = shape;
//   tx[2] = -lgamma(shape); // normalize
//   return atomic::inv_incpl_gamma(tx)[0] * scale;
// }
// VECTORIZE3_ttt(qgamma)

// /** \brief Distribution function of the poisson distribution (following R argument convention).
//     \ingroup R_style_distribution
// */
// template<class Type>
// Type ppois(Type q, Type lambda){
//   CppAD::vector<Type> tx(2);
//   tx[0] = q;
//   tx[1] = lambda;
//   return atomic::ppois(tx)[0];
// }
// VECTORIZE2_tt(ppois)

// /** 	@name Exponential distribution.
// 	Functions relative to the exponential distribution.
// 	*/
// /**@{*/
// /**	\brief Cumulative distribution function of the exponential distribution.
// 	\ingroup R_style_distribution
// 	\param rate Rate parameter. Must be strictly positive.
// 	*/
// template<class Type> 
// Type pexp(Type x, Type rate)
// {
// 	return CppAD::CondExpGe(x,Type(0),1-exp(-rate*x),Type(0));
// }

// // Vectorize pexp
// VECTORIZE2_tt(pexp)

// /**	\brief Probability density function of the exponential distribution.
// 	\ingroup R_style_distribution
// 	\param rate Rate parameter. Must be strictly positive.
// 	\param give_log true if one wants the log-probability, false otherwise.
// 	*/
template<class Type> 
Type nimDerivs_dexp_nimble(Type x, Type rate, Type give_log)
{
	Type res = CppAD::CondExpEq(give_log, Type(0) ,rate*exp(-rate*x), log(rate)-rate*x);
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
	
template<class Type>
Type nimDerivs_dexp(Type x, Type scale, Type give_log)
{
	Type res = CppAD::CondExpEq(give_log, Type(0), exp(-x/scale)/scale, -log(scale)-x/scale);
	return(res);
}

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
// Vectorize dexp
// VECTORIZE3_tti(dexp)

/**	\brief Inverse cumulative distribution function of the exponential distribution.
	\ingroup R_style_distribution
	\param rate Rate parameter. Must be strictly positive.
	*/
// template <class Type>
// Type qexp(Type p, Type rate)
// {
// 	return -log(1-p)/rate;
// }

// // Vectorize qexp.
// VECTORIZE2_tt(qexp)
// /**@}*/


// /**	@name Weibull distribution.
// 	Functions relative to the Weibull distribution.
// 	*/
// /**@{*/
// /** 	\brief Cumulative distribution function of the Weibull distribution.
// 	\ingroup R_style_distribution
// 	\param shape Shape parameter. Must be strictly positive.
// 	\param scale Scale parameter. Must be strictly positive.
// 	*/
// // template<class Type> 
// // Type pweibull(Type x, Type shape, Type scale)
// // {
// // 	return CppAD::CondExpGe(x,Type(0),1-exp(-pow(x/scale,shape)),Type(0));
// // }

// // // Vectorize pweibull
// // VECTORIZE3_ttt(pweibull)

// /** 	\brief Probability density function of the Weibull distribution.
// 	\ingroup R_style_distribution
// 	\param shape Shape parameter. Must be strictly positive.
// 	\param scale Scale parameter. Must be strictly positive.
// 	\param give_log true if one wants the log-probability, false otherwise.
// 	*/
template<class Type> 
Type nimDerivs_dweibull(Type x, Type shape, Type scale, Type give_log)
{
	Type res = shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape));
	res = CppAD::CondExpEq(give_log, Type(0), res, log(res));
	return(res);
}

template<class Type> 
Type nimDerivs_dweibull_logFixed(Type x, Type shape, Type scale, int give_log)
{
	Type res = shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape));
	if(give_log){
		res = log(res);
	}
	return(res);
}

// Vectorize dweibull
// VECTORIZE4_ttti(dweibull)

/**	\brief Inverse cumulative distribution function of the Weibull distribution.
	\ingroup R_style_distribution
	\param p Probability ; must be between 0 and 1.
	\param shape Shape parameter. Must be strictly positive.
	\param scale Scale parameter. Must be strictly positive.
	*/
// template<class Type> 
// Type qweibull(Type p, Type shape, Type scale)
// {
// 	Type res = scale * pow( (-log(1-p)) , 1/shape );
// 	res = CppAD::CondExpLt(p,Type(0),Type(0),res);
// 	res = CppAD::CondExpGt(p,Type(1),Type(0),res);
// 	return res;
// }

// // Vectorize qweibull
// VECTORIZE3_ttt(qweibull)
// /**@}*/

// /**	\brief Probability mass function of the binomial distribution.
// 	\ingroup R_style_distribution
// 	\param k Number of successes.
// 	\param size Number of trials.
// 	\param prob Probability of success.
// 	\param give_log true if one wants the log-probability, false otherwise.
// 	*/

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

template<class Type> 
Type nimDerivs_lfactorial(Type x) {
	return(lgamma(x + 1));
} 

template<class Type> 
Type nimDerivs_factorial(Type x) {
	return(exp(lgamma(1+x)));
}



// Vectorize dbinom
// VECTORIZE4_ttti(dbinom)

/** \brief Density of binomial distribution parameterized via logit(prob)

    This version should be preferred when working on the logit scale
    as it is numerically stable for probabilities close to 0 or 1.

    \ingroup R_style_distribution
*/
// template<class Type>
// Type dbinom_robust(Type k, Type size, Type logit_p, int give_log=0)
// {
//   CppAD::vector<Type> tx(4);
//   tx[0] = k;
//   tx[1] = size;
//   tx[2] = logit_p;
//   tx[3] = 0;
//   Type ans = atomic::log_dbinom_robust(tx)[0]; /* without norm. constant */
//   if (size > 1) {
//     ans += lgamma(size+1.) - lgamma(k+1.) - lgamma(size-k+1.);
//   }
//   return ( give_log ? ans : exp(ans) );
// }
// VECTORIZE4_ttti(dbinom_robust)

// /**	\brief Probability density function of the beta distribution.
// 	\ingroup R_style_distribution
// 	\param shape1 First shape parameter. Must be strictly positive.
// 	\param shape2 Second shape parameter. Must be strictly positive.
// 	\param give_log true if one wants the log-probability, false otherwise.
// 	*/
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
// Vectorize dbeta
// VECTORIZE4_ttti(dbeta)

/**	\brief Probability density function of the Fisher distribution.
	\ingroup R_style_distribution
	\param df1 Degrees of freedom 1.
	\param df2 Degrees of freedom 2.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
// template <class Type>
// Type df(Type x, Type df1, Type df2, int give_log)
// {
// 	Type logres = lgamma((df1+df2)/2.) - lgamma(df1/2.) - lgamma(df2/2.) + df1/2.*log(Type(df1)/df2) + (df1/2.-1)*log(x) - (df1+df2)/2.*log(1+Type(df1)/df2*x);
// 	if(!give_log) return exp(logres);
// 	else return logres;
// }

// //Vectorize df
// VECTORIZE4_ttti(df)

// /**	\brief Probability density function of the logistic distribution.
// 	\ingroup R_style_distribution
// 	\param location Location parameter.
// 	\param scale Scale parameter. Must be strictly positive.
// 	\param give_log true if one wants the log-probability, false otherwise.
// 	*/
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
	if(!give_log){
		res = exp(res);
	}
	return(res);
}

// Vectorize dlogis
// VECTORIZE4_ttti(dlogis)

/**	\brief Probability density function of the skew-normal distribution.
	\ingroup R_style_distribution
	\param alpha Slant parameter.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
// template <class Type>
// Type dsn(Type x, Type alpha, int give_log=0)
// {
	
// 	if(!give_log) return 2 * dnorm(x,Type(0),Type(1),0) * pnorm(alpha*x);
// 	else return log(2.0) + log(dnorm(x,Type(0),Type(1),0)) + log(pnorm(alpha*x));
// }

// // Vectorize dsn
// VECTORIZE3_tti(dsn)

// /** 	\brief Probability density function of the Student t-distribution.
// 	\ingroup R_style_distribution
// 	\param df Degree of freedom.
// 	\param give_log true if one wants the log-probability, false otherwise.
// 	*/	

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
  if(!give_log){
	  res = exp(res);
  }
  return(res);
}

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

template <class Type>
Type nimDerivs_dt(Type x, Type df, Type give_log)
{
  Type res = lgamma((df+1)/2) - Type(1)/2*log(df*M_PI) -lgamma(df/2) - (df+1)/2*log(1+x*x/df);
  res = CppAD::CondExpEq(give_log, Type(1), res, exp(res));
  return(res);
}

template <class Type>
Type nimDerivs_dt_logFixed(Type x, Type df, int give_log)
{
  Type res = lgamma((df+1)/2) - Type(1)/2*log(df*M_PI) -lgamma(df/2) - (df+1)/2*log(1+x*x/df);
  if(!give_log){
	res = exp(res);
  }
  return(res);
}


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

// template <class Type>
// Type nimDerivs_nimArr_dmvt_chol(NimArr<1, Type> &x, NimArr<1, Type> &mu, NimArr<2, Type> &chol, Type df, Type prec_param, Type give_log, Type overwrite_inputs) { 
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VectorXt;

//   int n = x.size();
//   Type logres = CppAD::CondExpEq(Type(mu.size()), Type(n), Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN()));
//   logres += CppAD::CondExpEq(Type(chol.dim()[0]), Type(n),  CppAD::CondExpEq(Type(chol.dim()[1]), Type(n), Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN())), Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  
//   logres += lgamma((df + n) / Type(2)) - lgamma(df / Type(2)) - n * Type(M_LN_SQRT_PI) - n * log(df) / Type(2);
//   int i;
//   Type logCholSum;
//   for(i = 0; i < n*n; i += n + 1){
// 	logCholSum += log(chol[i]);
//   }
//   logres += CppAD::CondExpEq(prec_param, Type(1), logCholSum, -logCholSum);	
//   VectorXt eigenXcopy(n);	
//   for(i = 0; i < n; i++)
// 	eigenXcopy(i,0) = x[i] - mu[i];

//   Eigen::Map<MatrixXt > eigenChol(chol.getPtr(), n, n); 
//   MatrixXt eigenXcopyCov = eigenChol.template triangularView<Eigen::Upper>()*eigenXcopy;
//   MatrixXt eigenXcopyPrec = eigenChol.template triangularView<Eigen::Upper>().solve(eigenXcopy).transpose();

//   if(Integer(prec_param) == 0){
// 	eigenXcopy = eigenChol.template triangularView<Eigen::Upper>()*eigenXcopy;
//   }
//   else{
// 	 eigenXcopy = eigenChol.template triangularView<Eigen::Upper>().solve(eigenXcopy).transpose();
//   }
//   // sum of squares to calculate quadratic form
//   Type tmp = Type(0.0);
//   for(i = 0; i < n; i++)
// 	tmp += CppAD::CondExpEq(prec_param, Type(0),  eigenXcopyCov(i,0) * eigenXcopyCov(i,0), eigenXcopyPrec(i,0) * eigenXcopyPrec(i,0));

//   logres += Type(-0.5) * (df + n) * log(Type(1) + tmp / df);
// 	logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
// 	return(logres);
// }

// template <class Type>
// Type nimDerivs_nimArr_dmvt_chol_logFixed(NimArr<1, Type> &x, NimArr<1, Type> &mu, NimArr<2, Type> &chol, Type df, Type prec_param, int give_log, Type overwrite_inputs) { 
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VectorXt;

//   int n = x.size();
//   Type logres = CppAD::CondExpEq(Type(mu.size()), Type(n), Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN()));
//   logres += CppAD::CondExpEq(Type(chol.dim()[0]), Type(n),  CppAD::CondExpEq(Type(chol.dim()[1]), Type(n), Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN())), Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  
//   logres += lgamma((df + n) / Type(2)) - lgamma(df / Type(2)) - n * Type(M_LN_SQRT_PI) - n * log(df) / Type(2);
//   int i;
//   Type logCholSum;
//   for(i = 0; i < n*n; i += n + 1){
// 	logCholSum += log(chol[i]);
//   }
//   logres += CppAD::CondExpEq(prec_param, Type(1), logCholSum, -logCholSum);	
//   VectorXt eigenXcopy(n);	
//   for(i = 0; i < n; i++)
// 	eigenXcopy(i,0) = x[i] - mu[i];

//   Eigen::Map<MatrixXt > eigenChol(chol.getPtr(), n, n); 
//   MatrixXt eigenXcopyCov = eigenChol.template triangularView<Eigen::Upper>()*eigenXcopy;
//   MatrixXt eigenXcopyPrec = eigenChol.template triangularView<Eigen::Upper>().solve(eigenXcopy).transpose();

//   if(Integer(prec_param) == 0){
// 	eigenXcopy = eigenChol.template triangularView<Eigen::Upper>()*eigenXcopy;
//   }
//   else{
// 	 eigenXcopy = eigenChol.template triangularView<Eigen::Upper>().solve(eigenXcopy).transpose();
//   }
//   // sum of squares to calculate quadratic form
//   Type tmp = Type(0.0);
//   for(i = 0; i < n; i++)
// 	tmp += CppAD::CondExpEq(prec_param, Type(0),  eigenXcopyCov(i,0) * eigenXcopyCov(i,0), eigenXcopyPrec(i,0) * eigenXcopyPrec(i,0));

//   logres += Type(-0.5) * (df + n) * log(Type(1) + tmp / df);
//   if(!give_log){
//     logres = exp(logres);
//   }
// 	return(logres);
// }




// template<class Type>
// CppAD::AD<Type> nimDerivs_nimArr_dcat(CppAD::AD<Type> x, CppAD::VecAD<Type> &prob, int give_log){
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VectorXt;
//   int n = prob.size();
// //   Eigen::Map<VectorXt > eigenProb(prob.getPtr(), prob.size(), 1);
// //   Type roundX = discrete_round(x);
// //   Type u;
// //   for(u = 0; u < n; u += 1.){
// // 	adProb[u] = prob[ Integer(u) ];
// //   }
//   CppAD::AD<Type> thisProb = prob[x - Type(1)];  
//   Type probSum = Type(12);
//   CppAD::AD<Type> res = thisProb/probSum;

// //   Type resMult = CppAD::CondExpGt(roundX, Type(0), Type(1), Type(0));
// //   resMult = CppAD::CondExpLe(roundX, Type(n), resMult, Type(0));
// //   resMult = CppAD::CondExpEq(x, roundX,  resMult, Type(0));
// //   CppAD::AD<Type> outRes = res*resMult;
// //   outRes  = CppAD::CondExpEq(Type(give_log), Type(0), outRes, log(outRes));
//   return(res);
// };

// Vectorize dt
// VECTORIZE3_tti(dt)

/** 	\brief Probability mass function of the multinomial distribution.
	\ingroup R_style_distribution
	\param x Vector of length K of integers.
        \param p Vector of length K, specifying the probability for the K classes (note, unlike in R these must sum to 1).
	\param give_log true if one wants the log-probability, false otherwise.
	*/

// template <class Type>
// Type nimDerivs_nimArr_dmulti(NimArr<1, Type> &x, Type size, NimArr<1, Type> &prob, int give_log)
// {
//   typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VectorXt;
//   Eigen::Map<VectorXt > eigenX(x.getPtr(), x.size(), 1);
//   Eigen::Map<VectorXt > eigenP(prob.getPtr(), prob.size(), 1);
//   VectorXt oneVec = VectorXt::Ones(x.size(), 1);
//   VectorXt eigenXp1 = eigenX + oneVec;
//   Type logres = lgamma(eigenX.sum() + Type(1)) - lgamma(eigenXp1.sum()) + Type((eigenX*eigenP.array().log().matrix()).sum());
// 	if(give_log) return logres;
// 	else return exp(logres);
// }



