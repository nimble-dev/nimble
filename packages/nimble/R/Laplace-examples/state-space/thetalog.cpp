#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  PARAMETER_VECTOR(u);  
  PARAMETER(logr0);       Type r0 = exp(logr0);
  PARAMETER(logpsi);      Type psi = exp(logpsi);
  PARAMETER(logK);        Type K = exp(logK);
  PARAMETER(logQ);        Type Q = exp(logQ); 
  PARAMETER(logR);        Type R = exp(logR);
  
  int n = y.size();
  Type f = 0;
  for(int t = 1; t < n; t++){
    // Type mean = u[t-1] + r0 * (1 - pow(exp(u[t-1])/K, psi));
    // Type mean = u[t-1] + r0 * (1 - psi* exp(u[t-1]) / K);
    Type mean = u[t-1] + r0 * psi * (1 - exp(u[t-1]) / K);
    f -= dnorm(u[t], mean, sqrt(Q), true);
  }
  
  for(int t = 0; t < n; t++){
    f -= dnorm(y[t], u[t], sqrt(R), true);
  }
  
  REPORT(f);
  REPORT(u);
  return f;
  
}
