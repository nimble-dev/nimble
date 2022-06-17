#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(X);
  DATA_VECTOR(t);
  PARAMETER_VECTOR(u);
  PARAMETER(alpha_transform);
  PARAMETER(beta_transform);
  
  int N = X.size();
  Type ans=0;
  vector<Type> lambda(N);
  // exp transformation for original random effects 
  vector<Type> theta = exp(u);
  Type alpha = exp(alpha_transform);
  Type beta = exp(beta_transform);
  
  // Negative log-likelihood
  for(int i = 0; i < N; i++){
    ans -= dgamma(theta[i], alpha, 1/beta, true);
    lambda[i] = theta[i] * t[i];
    ans -= dpois(X[i], lambda[i], true) + u[i];
  }

  REPORT(u);
  REPORT(ans);

  return ans;
}
