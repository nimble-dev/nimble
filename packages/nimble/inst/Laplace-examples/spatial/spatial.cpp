// Spatial poisson GLMM on a grid, with exponentially decaying correlation function
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(dd);
  PARAMETER_VECTOR(b);
  //PARAMETER(a);
  PARAMETER(log_a);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(u);

  using namespace density;
  int i,j;
  Type res=0;
  Type a = exp(log_a);

  vector<Type> eta(n); 
  eta = X*b + exp(log_sigma)*u;

  // 
  matrix<Type> cov(n,n); 
  for (i=0;i<n;i++)
  {
    cov(i,i)=Type(1);
    for (j=0;j<i;j++)
    {
      cov(i,j)=exp(-a*dd(i,j));			// Exponentially decaying correlation
      cov(j,i)=cov(i,j);
    }
  }

  MVNORM_t<Type> neg_log_density(cov);
  res+=neg_log_density(u);

  // logdpois = N log lam - lam (I added the constant '-lgamma(y[i]+1)' below)
  for(i=0;i<n;i++) res -= y[i]*eta[i]-exp(eta[i])-lgamma(y[i]+1);
  // for(i=0;i<n;i++) res -= eta[i]-y[i]*exp(eta[i]);
  
  REPORT(res); // Negative log-likelihood value
  REPORT(eta)
  REPORT(u);
  // ADREPORT(a);
  return res;

}
