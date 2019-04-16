#ifndef _NIMDERIVS_TMB__
#define _NIMDERIVS_TMB__

#include <TMB/dnorm.hpp>

// Functions here connect from nimble-generated C++ to TMB code that uses CppAD.
// TMB provides many nice distribution functions.
// Note that these functions are called rarely, typically once, because
// they are only used when recording a CppAD tape of operations.  This tape
// is then the object used to obtain derivatives of the recorded operations.
// For this reason, function call overhead is not a major concern.
//
// For clarify, we have named functions uniquely rather than relying on
// overloading.  Calls to these functions are typically code-generated.

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

/* TO-DO: dlnorm: Is it in TMB? */

#endif _NIMDERIVS_TMB__
