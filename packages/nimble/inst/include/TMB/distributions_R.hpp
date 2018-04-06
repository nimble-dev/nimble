// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

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
Type nimDerivs_pow(Type x, Type y) {
	Type outVal = CppAD::CondExpGt(x, Type(0), CppAD::pow(x, y), CppAD::CondExpEq(y, discrete_round(y), CppAD::pow(x, CppAD::Integer(y)), Type(CppAD::numeric_limits<Type>::quiet_NaN())));
	return(outVal);
}

template<class Type> 
Type nimDerivs_pow(Type x, int y) {
	Type outVal = CppAD::pow(x, y);
	return(outVal);
}

template<class Type>
Type nimDerivs_nimArr_dwish_chol(NimArr<2, Type> &xNimArr, NimArr<2, Type> &cholNimArr,
	 Type df, Type scale_param, Type give_log, Type overwrite_inputs){
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  int p = xNimArr.dim()[0];
  int i, j;

//   if (R_IsNA(x, p*p) || R_IsNA(chol, p*p) || R_IsNA(df) || R_IsNA(scale_param))
//     return NA_REAL;
// #ifdef IEEE_754
//   if (R_isnancpp(x, p*p) || R_isnancpp(chol, p*p) || R_isnancpp(df) || R_isnancpp(scale_param))
//     return R_NaN;
// #endif

  // also covers df < 0
//   if(df < (double) p) ML_ERR_return_NAN;

//   if(!R_FINITE_VEC(x, p*p) || !R_FINITE_VEC(chol, p*p)) return R_D__0;
  Type dens = -(df*p/2 * Type(M_LN2) + p*(p-1)*M_LN_SQRT_PI/2);
  for(i = 0; i < p; i++)
    dens -= lgamma((df - Type(i)) / Type(2.0));

  Type sumLogCholNimArr;
  for(i = 0; i < p; i++){ 
	sumLogCholNimArr += df * log(cholNimArr(i, i));
  }
  dens += CppAD::CondExpEq(scale_param, Type(1), -sumLogCholNimArr, sumLogCholNimArr);

//   determinant of x using Cholesky:
  Eigen::Map<MatrixXt > eigenX(xNimArr.getPtr(), p, p);

   Eigen::LLT<MatrixXt > llt(eigenX);
   MatrixXt    eigenXChol = llt.matrixL();
 
eigenXChol.template triangularView<Eigen::Upper>() =
		eigenXChol.transpose().template triangularView<Eigen::Upper>();

for(i = 0; i < p; i++){ 
  	dens += (df - p - 1) *  log(eigenXChol(i, i));
}

 Eigen::Map<MatrixXt > eigenChol(cholNimArr.getPtr(), p, p); // may need to create new eigen matrix instead of mapping here
 Type tmp_dens = 0.0;
 MatrixXt eigenSolved = eigenChol.transpose().lu().solve(eigenXChol.transpose());
 MatrixXt cholMultX = eigenChol*eigenX;
	for(j = 0; j < p; j++){ 
		for(i = 0; i <= j; i++){ 
			 tmp_dens +=  CppAD::CondExpEq(scale_param, Type(1), eigenSolved(j,i)*eigenSolved(j,i), 
			   cholMultX(i, j) * eigenChol(i, j));
		}
	}
  dens += -0.5 * tmp_dens;
  return CppAD::CondExpEq(give_log, Type(1.0), dens, exp(dens));
}

template<class Type>
Type nimDerivs_pnorm(Type q, Type mean = 0., Type sd = 1.){
  CppAD::vector<Type> tx(1);
  tx[0] = (q - mean) / sd;
  return ::atomic::pnorm1(tx)[0];
}
// VECTORIZE3_ttt(pnorm)
// VECTORIZE1_t(pnorm)

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
	Type res = CppAD::CondExpEq(give_log, Type(0), CppAD::CondExpGe(rate, Type(0), CppAD::CondExpGe(x,Type(0),rate*exp(-rate*x),Type(0)), Type(CppAD::numeric_limits<Type>::quiet_NaN())),
					CppAD::CondExpGe(rate, Type(0), CppAD::CondExpGe(x,Type(0),log(rate)-rate*x,Type(-INFINITY)), Type(CppAD::numeric_limits<Type>::quiet_NaN())));
	// if(!give_log){
	// 	Type res = CppAD::CondExpGe(x,Type(0),rate*exp(-rate*x),Type(0));
	// 	return  CppAD::CondExpGe(rate, Type(0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
	// }
	// else{
	// 	Type res = CppAD::CondExpGe(x,Type(0),log(rate)-rate*x,Type(-INFINITY));
	// 	return  CppAD::CondExpGe(rate, Type(0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
	// }
	return(res);
}
	
template<class Type>
Type nimDerivs_dexp(Type x, Type scale, Type give_log)
{
	Type res = CppAD::CondExpEq(give_log, Type(0),  CppAD::CondExpGe(scale, Type(0),  CppAD::CondExpGe(x,Type(0),exp(-x/scale)/scale,Type(0)), Type(CppAD::numeric_limits<Type>::quiet_NaN())),
					CppAD::CondExpGe(scale, Type(0),  CppAD::CondExpGe(x,Type(0),-log(scale)-x/scale,Type(-INFINITY)), Type(CppAD::numeric_limits<Type>::quiet_NaN())));
	// if(!give_log){
	// 	Type res = CppAD::CondExpGe(x,Type(0),exp(-x/scale)/scale,Type(0));
	// 	return  CppAD::CondExpGe(scale, Type(0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
	// }
	// else{
	// 	Type res = CppAD::CondExpGe(x,Type(0),-log(scale)-x/scale,Type(-INFINITY));
	// 	return CppAD::CondExpGe(scale, Type(0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
	// }
	return(res);
}

	
template<class Type> 
Type nimDerivs_dunif(Type x, Type a, Type b, Type give_log)
{   
	Type res = CppAD::CondExpGt(b, a, Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN())); 
	res += CppAD::CondExpGe(x, a, CppAD::CondExpLe(x, b, 1/(b-a), Type(0)), Type(0));
	res = CppAD::CondExpEq(give_log, Type(1), log(res), res);
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
	Type res = CppAD::CondExpGt(shape, Type(0), CppAD::CondExpGt(scale, Type(0), Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN())), Type(CppAD::numeric_limits<Type>::quiet_NaN()));
	res += CppAD::CondExpGe(x, Type(0), shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape)), Type(0));
	res = CppAD::CondExpEq(give_log, Type(0), res, log(res));
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
  Type disc_k = discrete_round(k);
  Type disc_size = discrete_round(size);
  Type logres = CondExpEq(k, disc_k,  lgamma(disc_size+1)-lgamma(disc_k+1)-lgamma(disc_size-disc_k+1)+disc_k*log(prob)+(disc_size-disc_k)*log(1-prob),
   -Type( std::numeric_limits<double>::infinity()));
  logres = CondExpEq(size, disc_size, logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres = CondExpGe(size, Type(0.0), logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres = CondExpGe(prob, Type(0.0), logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres = CondExpLe(prob, Type(1.0), logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres = CondExpEq(give_log, Type(1), logres, exp(logres));
  return(logres);
}

template<class Type>
Type nimDerivs_dchisq(Type x, Type df, Type give_log)
{	
  Type logres = CondExpGt(x, Type(0.0), (df/Type(2.0) - Type(1.0))*log(x) - (x/Type(2.0)) - (df/Type(2.0))*log(Type(2.0)) - lgamma(df/Type(2.0)), 
	  -Type( std::numeric_limits<double>::infinity()));	
  logres =  CondExpEq(x, Type(0.0), Type( std::numeric_limits<double>::infinity()), logres); // to match R: dchisq(0, df) is infinity in R.
  logres = CondExpGt(df, Type(0.0), logres, -Type( std::numeric_limits<double>::infinity()));	
  logres =  CondExpLt(df, Type(0.0), Type(CppAD::numeric_limits<Type>::quiet_NaN()), logres); 
  logres = CondExpEq(give_log, Type(0), CondExpEq(df, Type(0.0), zero_NaNderiv(df), exp(logres)), logres);
  return(logres);
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
							CondExpGe(x, Type(0.0), CondExpLe(x, Type(1.0), 
							exp(lgamma(shape1+shape2) - lgamma(shape1) - 
							lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1),
							Type(0.0)), Type(0.0)),
							CondExpLe(x, Type(1.0), CondExpGe(x, Type(0.0), CondExpEq(x, Type(0.0), 
								log(exp(lgamma(shape1+shape2) - lgamma(shape1) - 
								lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1)),
								lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2) + 
								(shape1-1)*log(x) + (shape2-1)*log(1-x)), -Type(std::numeric_limits<double>::infinity() )), -Type(std::numeric_limits<double>::infinity() )));
		res = CondExpGt(shape1, Type(0.0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN())) ;
		res = CondExpGt(shape2, Type(0.0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN())) ;

	// Type res;
	// if(!give_log){
	// 	res = CondExpLe(x, Type(1.0), 
	// 						exp(lgamma(shape1+shape2) - lgamma(shape1) - 
	// 						lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1),
	// 						Type(0.0));
	// 	res = CondExpGe(x, Type(0.0), res, Type(0.0)) ;
	// 	res = CondExpGt(shape1, Type(0.0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN())) ;
	// 	res = CondExpGt(shape2, Type(0.0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN())) ;
	// }
	// if(give_log){
	// 	res = CondExpEq(x, Type(0.0), 
	// 							log(exp(lgamma(shape1+shape2) - lgamma(shape1) - 
	// 							lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1)),
	// 							lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2) + 
	// 							(shape1-1)*log(x) + (shape2-1)*log(1-x));
	// 	res = CondExpGe(x, Type(0.0), res, -Type(std::numeric_limits<double>::infinity() ));
	// 	res = CondExpLe(x, Type(1.0), res, -Type(std::numeric_limits<double>::infinity() ));
	// 	res = CondExpGt(shape1, Type(0.0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN())) ;
	// 	res = CondExpGt(shape2, Type(0.0), res, Type(CppAD::numeric_limits<Type>::quiet_NaN())) ;
	// } 
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
	Type logres = -(x-location)/scale - log(scale) - 2*log(1+exp(-(x-location)/scale));
	logres = CppAD::CondExpGt(scale, Type(0.0), logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
	logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
	return(logres);
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
  Type logres;
  logres=-log(Type(sqrt(2*M_PI))*sd)-Type(.5)*pow((x-mean)/sd,2);
  logres = CppAD::CondExpGe(sd, Type(0.0), logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres = CppAD::CondExpEq(sd, Type(0.0), -Type(std::numeric_limits<double>::infinity()) ,  logres);
  logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
  return(logres);
}

template<class Type>
Type nimDerivs_dlnorm(Type x, Type mean, Type sd, Type give_log)
{
  Type logres = CppAD::CondExpGe(x, Type(0.0), -log(x) - log(sd) - Type(0.5)*log(Type(2.0*M_PI)) - Type(.5)*pow((log(x)-mean)/sd, 2), -Type(std::numeric_limits<double>::infinity())); 
  logres = CppAD::CondExpGe(sd, Type(0.0), logres, Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres = CppAD::CondExpEq(sd, Type(0.0), -Type(std::numeric_limits<double>::infinity()) ,  logres);
  logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
  return(logres);
}


template <class Type>
Type nimDerivs_dt(Type x, Type df, Type give_log)
{
  Type logres = lgamma((df+1)/2) - Type(1)/2*log(df*M_PI) -lgamma(df/2) - (df+1)/2*log(1+x*x/df);
  logres =  CondExpLt(df, Type(0.0), Type(CppAD::numeric_limits<Type>::quiet_NaN()), logres); 
  logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
  return(logres);
}



template <class Type>
Type nimDerivs_dt_nonstandard(Type x, Type df, Type mu, Type sigma, Type give_log)
{
  Type logres = CondExpGe(sigma, Type(0), CondExpEq(sigma, Type(0), CondExpEq(x, mu, Type(std::numeric_limits<double>::infinity()), -Type(std::numeric_limits<double>::infinity())), Type(0)),
	Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres +=  nimDerivs_dt( (x - mu)/sigma, df, 1) - log(sigma);
  logres =  CondExpLt(df, Type(0.0), Type(CppAD::numeric_limits<Type>::quiet_NaN()), logres); 
  logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
  return(logres);
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
Type nimDerivs_nimArr_dmvt_chol(NimArr<1, Type> &x, NimArr<1, Type> &mu, NimArr<2, Type> &chol, Type df, Type prec_param, Type give_log, Type overwrite_inputs) { 
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXt;
  typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VectorXt;

  int n = x.size();
  Type logres = CppAD::CondExpEq(Type(mu.size()), Type(n), Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  logres += CppAD::CondExpEq(Type(chol.dim()[0]), Type(n),  CppAD::CondExpEq(Type(chol.dim()[1]), Type(n), Type(0), Type(CppAD::numeric_limits<Type>::quiet_NaN())), Type(CppAD::numeric_limits<Type>::quiet_NaN()));
  
  logres += lgamma((df + n) / Type(2)) - lgamma(df / Type(2)) - n * Type(M_LN_SQRT_PI) - n * log(df) / Type(2);
  int i;
  Type logCholSum;
  for(i = 0; i < n*n; i += n + 1){
	logCholSum += log(chol[i]);
  }
  logres += CppAD::CondExpEq(prec_param, Type(1), logCholSum, -logCholSum);	
  VectorXt eigenXcopy(n);	
  for(i = 0; i < n; i++)
	eigenXcopy(i,0) = x[i] - mu[i];

  Eigen::Map<MatrixXt > eigenChol(chol.getPtr(), n, n); 
  MatrixXt eigenXcopyCov = eigenChol.template triangularView<Eigen::Upper>()*eigenXcopy;
  MatrixXt eigenXcopyPrec = eigenChol.template triangularView<Eigen::Upper>().solve(eigenXcopy).transpose();

  if(Integer(prec_param) == 0){
	eigenXcopy = eigenChol.template triangularView<Eigen::Upper>()*eigenXcopy;
  }
  else{
	 eigenXcopy = eigenChol.template triangularView<Eigen::Upper>().solve(eigenXcopy).transpose();
  }
  // sum of squares to calculate quadratic form
  Type tmp = Type(0.0);
  for(i = 0; i < n; i++)
	tmp += CppAD::CondExpEq(prec_param, Type(0),  eigenXcopyCov(i,0) * eigenXcopyCov(i,0), eigenXcopyPrec(i,0) * eigenXcopyPrec(i,0));

  logres += Type(-0.5) * (df + n) * log(Type(1) + tmp / df);
	logres = CppAD::CondExpEq(give_log, Type(1), logres, exp(logres));
	return(logres);
}




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


/** 	@name Sinh-asinh distribution.
  	Functions relative to the sinh-asinh distribution.
		*/
/**@{*/
/**	\brief Probability density function of the sinh-asinh distribution.
  	\ingroup R_style_distribution
	\param mu Location.
	\param sigma Scale.
	\param nu Skewness.
	\param tau Kurtosis.
	\param give_log true if one wants the log-probability, false otherwise.
					
	Notation adopted from R package "gamlss.dist".
				
	Probability density given in (2) in __Jones and Pewsey (2009) Biometrika (2009) 96 (4): 761-780__.
								
	It is not possible to call this function with nu a vector or tau a vector.
*/
// template <class Type>
// Type dSHASHo(Type x, Type mu, Type sigma, Type nu, Type tau, int give_log = 0)
// {
// 	// TODO : Replace log(x+sqrt(x^2+1)) by a better approximation for asinh(x).
		
// 	Type z = (x-mu)/sigma;
//    	Type c = cosh(tau*log(z+sqrt(z*z+1))-nu);
//    	Type r = sinh(tau*log(z+sqrt(z*z+1))-nu);
//    	Type logres = -log(sigma) + log(tau) -0.5*log(2*M_PI) -0.5*log(1+(z*z)) +log(c) -0.5*(r*r);
					   	
//   	if(!give_log) return exp(logres);
//    	else return logres;
// }

// // Vectorize dSHASHo
// VECTORIZE6_ttttti(dSHASHo)

// /**	\brief Cumulative distribution function of the sinh-asinh distribution.
//   	\ingroup R_style_distribution
// 	\param mu Location.
// 	\param sigma Scale.
// 	\param nu Skewness.
// 	\param tau Kurtosis.
// 	\param give_log true if one wants the log-probability, false otherwise.
		
// 	Notation adopted from R package "gamlss.dist".
	
// 	It is not possible to call this function with nu a vector or tau a vector.
// */
// template <class Type>
// Type pSHASHo(Type q,Type mu,Type sigma,Type nu,Type tau,int give_log=0)
// {
// 	// TODO : Replace log(x+sqrt(x^2+1)) by a better approximation for asinh(x).

// 	Type z = (q-mu)/sigma;
// 	Type r = sinh(tau * log(z+sqrt(z*z+1)) - nu);
// 	Type p = pnorm(r);
				  	
// 	if (!give_log) return p;
// 	else return log(p);
// }

// // Vectorize pSHASHo
// VECTORIZE6_ttttti(pSHASHo)

// /**	\brief Quantile function of the sinh-asinh distribution.
// 	\ingroup R_style_distribution
// 	\param mu Location.
// 	\param sigma Scale.
// 	\param nu Skewness.
// 	\param tau Kurtosis.
// 	\param log_p true if p is log-probability, false otherwise.
	
// 	Notation adopted from R package "gamlss.dist".
	
// 	It is not possible to call this function with nu a vector or tau a vector.
// 	*/
// template <class Type>
// Type qSHASHo(Type p, Type mu, Type sigma, Type nu, Type tau, int log_p = 0)
// {
// 	// TODO : Replace log(x+sqrt(x^2+1)) by a better approximation for asinh(x).

//    	if(!log_p) return mu + sigma*sinh((1/tau)* log(qnorm(p)+sqrt(qnorm(p)*qnorm(p)+1)) + (nu/tau));
//    	else return mu + sigma*sinh((1/tau)*log(qnorm(exp(p))+sqrt(qnorm(exp(p))*qnorm(exp(p))+1))+(nu/tau));
// }

// // Vectorize qSHASHo
// VECTORIZE6_ttttti(qSHASHo)

// /**	\brief Transforms a normal variable into a sinh-asinh variable.
// 	\param mu Location parameter of the result sinh-asinh distribution.
// 	\param sigma Scale parameter of the result sinh-asinh distribution.
// 	\param nu Skewness parameter of the result sinh-asinh distribution.
// 	\param tau Kurtosis parameter of the result sinh-asinh distribution.
// 	\param log_p true if p is log-probability, false otherwise.
	
// 	It is not possible to call this function with nu a vector or tau a vector.
// 	*/
// template <class Type>
// Type norm2SHASHo(Type x, Type mu, Type sigma, Type nu, Type tau, int log_p = 0)
// {

// 	return qSHASHo(pnorm(x),mu,sigma,nu,tau,log_p);
// }

// // Vectorize norm2SHASHo
// VECTORIZE6_ttttti(norm2SHASHo)
// /**@}*/

// /** \brief Distribution function of the beta distribution (following R
//     argument convention).
//     \note Non-centrality parameter (ncp) not implemented.
//     \ingroup R_style_distribution
// */
// template<class Type>
// Type pbeta(Type q, Type shape1, Type shape2){
//   CppAD::vector<Type> tx(4);
//   tx[0] = q;
//   tx[1] = shape1;
//   tx[2] = shape2;
//   tx[3] = 0; // order
//   Type ans = atomic::pbeta(tx)[0];
//   return ans;
// }
// VECTORIZE3_ttt(pbeta)

// /** \brief Quantile function of the beta distribution (following R
//     argument convention).
//     \note Non-centrality parameter (ncp) not implemented.
//     \ingroup R_style_distribution
// */
// template<class Type>
// Type qbeta(Type p, Type shape1, Type shape2){
//   CppAD::vector<Type> tx(3);
//   tx[0] = p;
//   tx[1] = shape1;
//   tx[2] = shape2;
//   Type ans = atomic::qbeta(tx)[0];
//   return ans;
// }
// VECTORIZE3_ttt(qbeta)

// /** \brief besselK function (same as besselK from R).
//     \note Derivatives wrt. both arguments are implemented
//     \ingroup special_functions
// */
// template<class Type>
// Type besselK(Type x, Type nu){
//   Type ans;
//   if(CppAD::Variable(nu)) {
//     CppAD::vector<Type> tx(3);
//     tx[0] = x;
//     tx[1] = nu;
//     tx[2] = 0;
//     ans = atomic::bessel_k(tx)[0];
//   } else {
//     CppAD::vector<Type> tx(2);
//     tx[0] = x;
//     tx[1] = nu;
//     ans = atomic::bessel_k_10(tx)[0];
//   }
//   return ans;
// }
// VECTORIZE2_tt(besselK)

// /** \brief besselI function (same as besselI from R).
//     \note Derivatives wrt. both arguments are implemented
//     \ingroup special_functions
// */
// template<class Type>
// Type besselI(Type x, Type nu){
//   Type ans;
//   if(CppAD::Variable(nu)) {
//     CppAD::vector<Type> tx(3);
//     tx[0] = x;
//     tx[1] = nu;
//     tx[2] = 0;
//     ans = atomic::bessel_i(tx)[0];
//   } else {
//     CppAD::vector<Type> tx(2);
//     tx[0] = x;
//     tx[1] = nu;
//     ans = atomic::bessel_i_10(tx)[0];
//   }
//   return ans;
// }
// VECTORIZE2_tt(besselI)

// /** \brief besselJ function (same as besselJ from R).
//     \note Derivatives wrt. both arguments are implemented
//     \ingroup special_functions
// */
// template<class Type>
// Type besselJ(Type x, Type nu){
//   CppAD::vector<Type> tx(3);
//   tx[0] = x;
//   tx[1] = nu;
//   tx[2] = 0;
//   Type ans = atomic::bessel_j(tx)[0];
//   return ans;
// }
// VECTORIZE2_tt(besselJ)

// /** \brief besselY function (same as besselY from R).
//     \note Derivatives wrt. both arguments are implemented
//     \ingroup special_functions
// */
// template<class Type>
// Type besselY(Type x, Type nu){
//   CppAD::vector<Type> tx(3);
//   tx[0] = x;
//   tx[1] = nu;
//   tx[2] = 0;
//   Type ans = atomic::bessel_y(tx)[0];
//   return ans;
// }
// VECTORIZE2_tt(besselY)

// /** \brief dtweedie function (same as dtweedie.series from R package
//     'tweedie').

//     Silently returns NaN if not within the valid parameter range:
//     \f[ (0 \leq y) \land (0 < \mu) \land (0 < \phi) \land (1 < p) \land (p < 2) \f] .

//     \note Parameter order differs from the R version.

//     \warning The derivative wrt. the y argument is disabled
//     (zero). Hence the tweedie distribution can only be used for *data*
//     (not random effects).

//     \ingroup R_style_distribution
// */
// template<class Type>
// Type dtweedie(Type y, Type mu, Type phi, Type p, int give_log = 0) {
//   CppAD::vector<Type> tx(5);
//   tx[0] = y;
//   tx[1] = mu;
//   tx[2] = phi;
//   tx[3] = p;
//   tx[4] = 0;
//   Type ans = atomic::log_dtweedie(tx)[0];
//   return ( give_log ? ans : exp(ans) );
// }


// /** \brief Conway-Maxwell-Poisson log normalizing constant.

//     \f[ Z(\lambda, \nu) = \sum_{i=0}^{\infty} \frac{\lambda^i}{(i!)^\nu} \f] .

//     \param loglambda \f$ \log(\lambda) \f$
//     \param nu \f$ \nu \f$

//     \return \f$ \log Z(\lambda, \nu) \f$
// */
// template<class Type>
// Type compois_calc_logZ(Type loglambda, Type nu) {
//   CppAD::vector<Type> tx(3);
//   tx[0] = loglambda;
//   tx[1] = nu;
//   tx[2] = 0;
//   return atomic::compois_calc_logZ(tx)[0];
// }
// VECTORIZE2_tt(compois_calc_logZ)

// /** \brief Conway-Maxwell-Poisson. Calculate log(lambda) from
//     log(mean).

//     \param logmean \f$ \log(E[X]) \f$
//     \param nu \f$ \nu \f$

//     \return \f$ \log \lambda \f$
// */
// template<class Type>
// Type compois_calc_loglambda(Type logmean, Type nu) {
//   CppAD::vector<Type> tx(3);
//   tx[0] = logmean;
//   tx[1] = nu;
//   tx[2] = 0;
//   return atomic::compois_calc_loglambda(tx)[0];
// }
// VECTORIZE2_tt(compois_calc_loglambda)

// /** \brief Conway-Maxwell-Poisson. Calculate density.

//     \f[ P(X=x) \propto \frac{\lambda^x}{(x!)^\nu}\:,x=0,1,\ldots \f]

//     Silently returns NaN if not within the valid parameter range:
//     \f[ (0 \leq x) \land (0 < \lambda) \land (0 < \nu) \f] .

//     \param x Observation
//     \param mode Approximate mode \f$ \lambda^\nu \f$
//     \param nu   \f$ \nu \f$

//     \ingroup R_style_distribution
// */
// template<class T1, class T2, class T3>
// T1 dcompois(T1 x, T2 mode, T3 nu, int give_log = 0) {
//   T2 loglambda = nu * log(mode);
//   T1 ans = x * loglambda - nu * lfactorial(x);
//   ans -= compois_calc_logZ(loglambda, nu);
//   return ( give_log ? ans : exp(ans) );
// }

// /** \brief Conway-Maxwell-Poisson. Calculate density parameterized via
//     the mean.

//     Silently returns NaN if not within the valid parameter range:
//     \f[ (0 \leq x) \land (0 < E[X]) \land (0 < \nu) \f] .

//     \param x Observation
//     \param mean \f$ E[X] \f$
//     \param nu   \f$ \nu \f$

//     \ingroup R_style_distribution
// */
// template<class T1, class T2, class T3>
// T1 dcompois2(T1 x, T2 mean, T3 nu, int give_log = 0) {
//   T2 logmean = log(mean);
//   T2 loglambda = compois_calc_loglambda(logmean, nu);
//   T1 ans = x * loglambda - nu * lfactorial(x);
//   ans -= compois_calc_logZ(loglambda, nu);
//   return ( give_log ? ans : exp(ans) );
// }

// /********************************************************************/
// /* SIMULATON CODE                                                   */
// /********************************************************************/

// extern "C" {
//   double Rf_rnorm(double mu, double sigma);
// }
// /** \brief Simulate from a normal distribution  */
// template<class Type>
// Type rnorm(Type mu, Type sigma)
// {
//   return Rf_rnorm(asDouble(mu), asDouble(sigma));
// }
// VECTORIZE2_tt(rnorm)
// VECTORIZE2_n(rnorm)

// extern "C" {
//   double Rf_rpois(double mu);
// }
// /** \brief Simulate from a Poisson distribution  */
// template<class Type>
// Type rpois(Type mu)
// {
//   return Rf_rpois(asDouble(mu));
// }
// VECTORIZE1_t(rpois)
// VECTORIZE1_n(rpois)

// extern "C" {
//   double Rf_runif(double a, double b);
// }
// /** \brief Simulate from a uniform distribution  */
// template<class Type>
// Type runif(Type a, Type b)
// {
//   return Rf_runif(asDouble(a), asDouble(b));
// }
// VECTORIZE2_tt(runif)
// VECTORIZE2_n(runif)

// extern "C" {
//   double Rf_rbinom(double size, double prob);
// }
// /** \brief Simulate from a binomial distribution  */
// template<class Type>
// Type rbinom(Type size, Type prob)
// {
//   return Rf_rbinom(asDouble(size), asDouble(prob));
// }
// VECTORIZE2_tt(rbinom)
// VECTORIZE2_n(rbinom)

// extern "C" {
//   double Rf_rgamma(double shape, double scale);
// }
// /** \brief Simulate from a gamma distribution  */
// template<class Type>
// Type rgamma(Type shape, Type scale)
// {
//   return Rf_rgamma(asDouble(shape), asDouble(scale));
// }
// VECTORIZE2_tt(rgamma)
// VECTORIZE2_n(rgamma)

// extern "C" {
//   double Rf_rexp(double rate);
// }
// /** \brief Simulate from an exponential distribution */
// template<class Type>
// Type rexp(Type rate)
// {
//   return Rf_rexp(asDouble(rate));
// }

// VECTORIZE1_t(rexp)
// VECTORIZE1_n(rexp)

// extern "C" {
// 	double Rf_rbeta(double shape1, double shape2);
// }
// /** \brief Simulate from a beta distribution */
// template<class Type>
// Type rbeta(Type shape1, Type shape2)
// {
// 	return Rf_rbeta(asDouble(shape1), asDouble(shape2));
// }

// VECTORIZE2_tt(rbeta)
// VECTORIZE2_n(rbeta)

// extern "C" {
// 	double Rf_rf(double df1, double df2);
// }
// /** \brief Simulate from an F distribution */
// template<class Type>
// Type rf(Type df1, Type df2)
// {
// 	return Rf_rf(asDouble(df1), asDouble(df2));
// }

// VECTORIZE2_tt(rf)
// VECTORIZE2_n(rf)

// extern "C" {
// 	double Rf_rlogis(double location, double scale);
// }
// /** \brief Simulate from a logistic distribution */
// template<class Type>
// Type rlogis(Type location, Type scale)
// {
// 	return Rf_rlogis(asDouble(location), asDouble(scale));
// }

// VECTORIZE2_tt(rlogis)
// VECTORIZE2_n(rlogis)

// extern "C" {
// 	double Rf_rt(double df);
// }
// /** \brief Simulate from a Student's t distribution */
// template<class Type>
// Type rt(Type df)
// {
// 	return Rf_rt(asDouble(df));
// }

// VECTORIZE1_t(rt)
// VECTORIZE1_n(rt)

// extern "C" {
// 	double Rf_rweibull(double shape, double scale);
// }
// /** \brief Simulate from a Weibull distribution */
// template<class Type>
// Type rweibull(Type shape, Type scale)
// {
// 	return Rf_rweibull(asDouble(shape), asDouble(scale));
// }

// VECTORIZE2_tt(rweibull)
// VECTORIZE2_n(rweibull)

// /** \brief Simulate from a Conway-Maxwell-Poisson distribution  */
// template<class Type>
// Type rcompois(Type mode, Type nu)
// {
//   Type loglambda = nu * log(mode);
//   return atomic::compois_utils::simulate(asDouble(loglambda), asDouble(nu));
// }
// VECTORIZE2_tt(rcompois)
// VECTORIZE2_n(rcompois)

// /** \brief Simulate from a Conway-Maxwell-Poisson distribution  */
// template<class Type>
// Type rcompois2(Type mean, Type nu)
// {
//   Type logmean = log(mean);
//   Type loglambda = compois_calc_loglambda(logmean, nu);
//   return atomic::compois_utils::simulate(asDouble(loglambda), asDouble(nu));
// }
// VECTORIZE2_tt(rcompois2)

// // Note: Vectorize manually to avoid many identical calls to
// // 'calc_loglambda'.
// template<class Type>
// vector<Type> rcompois2(int n, Type mean, Type nu)
// {
//   Type logmean = log(mean);
//   Type loglambda = compois_calc_loglambda(logmean, nu);
//   Type mode = exp(loglambda / nu);
//   return rcompois(n, mode, nu);
// }




