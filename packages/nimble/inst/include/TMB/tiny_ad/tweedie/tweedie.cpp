/************************************************************/
/*               Tweedie density evaluation                 */
/*              Author:  Wayne Zhang                        */
/*            actuary_zhang@hotmail.com                     */
/************************************************************/

/**
 * @file tweedie.c
 * @brief Approximating the density of the Tweedie compound Poisson 
 * distribution using the series evaluation method.
 * @author Wayne Zhang                         
 */

/** the threshold used in finding the bounds of the series */
#define TWEEDIE_DROP 37.0
/** the loop increment used in finding the bounds of the series */
#define TWEEDIE_INCRE 5
/** the max number of terms allowed in the finite sum approximation*/
#define TWEEDIE_NTERM 20000

/*
 * Compute the log density for the tweedie compound Poisson distribution.
 * This is based on the dtweedie.series function in R, but bounds
 * are not determined using all observations because this could result in 
 * dramatically slower computation in certain circumstances. 
 *
 * @param n the number of observations
 * @param y  the vector of observations
 * @param mu the vector of means
 * @param phi scalar: the dispersion parameter
 * @param p scalar: the index parameter
 * @param wts the optional vector of weights (set to NULL if no weight is supplied)
 * @param ans the vector that stores the computed log density
 *
 */

template<class Float>
void dtweedie(int n, Float *y, Float *mu, Float phi, Float p,
                       Float *wts, Float *ans){
    
  int np = 0, pos = 0 ;  
  Float p1 = p - 1.0, p2 = 2.0 - p;
  Float a = - p2 / p1, a1 = 1.0 / p1;
  Float cc, j, w, phiw, sum_ww = 0.0, ww_max ;

  /* # positive values */
  for (int i = 0; i < n; i++)
      if (y[i] != 0) np++ ;
  /* all zeros in the data (probably nonsense) */
  if (np == 0) {
    for (int i = 0; i < n; i++){
      phiw = wts ? (phi / wts[i]) : phi ; 
      ans[i] = -pow(mu[i], p2) / (phiw * p2) ;
    }
    return ;
  }
  /* only need the lower bound and the # terms to be stored */
  int jh, *jl = Calloc(np, int), *jd = Calloc(np, int); 
  Float *jmax = Calloc(np, Float), *logz = Calloc(np, Float);

  /* compute jmax for each y > 0*/
  cc = a * log(p1) - log(p2);
  for (int i = 0; i < n; i++){
      if (y[i] != 0) {
          phiw = wts ? (phi / wts[i]) : phi ;
          jmax[pos] = fmax2(1.0, pow(y[i], p2) / (phiw * p2));
          logz[pos] = - a * log(y[i]) - a1 * log(phiw) + cc;
          pos++ ;
      }
  }

  /* find bounds in the summation */
  for (int i = 0; i < np; i++){
    /* locate upper bound */
    cc = logz[i] + a1 + a * log(-a);
    j = jmax[i] ;
    w = a1 * j ;
    while (1) {
      j += TWEEDIE_INCRE ;
      if (j * (cc - a1 * log(j)) < (w - TWEEDIE_DROP))
        break ;
    }
    jh = ceil(j);  
    /* locate lower bound */
    j = jmax[i];
    while (1) {
      j -= TWEEDIE_INCRE ;
      if (j < 1 || j * (cc - a1 * log(j)) < w - TWEEDIE_DROP)
        break ;
    }
    jl[i] = imax2(1, floor(j)) ; 
    jd[i] = jh - jl[i] + 1;
  }

  /* set limit for # terms in the sum */
  int nterms = imin2(imax(jd, np), TWEEDIE_NTERM), iterm ;
  Float *ww = Calloc(nterms, Float) ;
  /* evaluate series using the finite sum*/
  pos = 0;  
  for (int i = 0; i < n; i++){    
      phiw = wts ? (phi / wts[i]) : phi ; 
      /* y == 0 */
      ans[i] = -pow(mu[i], p2) / (phiw * p2) ;
      /* y > 0 */
      if (y[i] != 0) {
          sum_ww = 0.0 ;
	  iterm = imin2(jd[pos], nterms) ;    /* avoid stepping out of bounds */
          for (int k = 0; k < iterm; k++){
              j = k + jl[pos] ;
              ww[k] = j * logz[pos] - lgammafn(1 + j) - lgammafn(-a * j);
          }
          ww_max = dmax(ww, iterm) ;
          for (int k = 0; k < iterm; k++)
              sum_ww += exp(ww[k] - ww_max);
          ans[i] += -y[i] / (phiw * p1 * pow(mu[i], p1)) - log(y[i]) +
              log(sum_ww) + ww_max  ;
          pos++;
      }
  }
  Free(jmax); 
  Free(logz);
  Free(jl); 
  Free(jd);
  Free(ww);
}
