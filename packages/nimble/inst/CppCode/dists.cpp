// various additional distributions functions needed by NIMBLE
// Author: Chris Paciorek 
// Date: Initial development in February 2014
// uses various BLAS routines and constants from R's C API
// compile as "R CMD SHLIB dists.cpp"

//#include "Utils.h" // moved to dists.h
#include "nimble/dists.h"
#include <R_ext/Lapack.h>

bool R_IsNA(double* P, int s) {
  for(int i = 0; i < s; ++i) if(R_IsNA(P[i])) return(true);
  return(false);
}

bool R_isnancpp(double* P, int s) {
  for(int i = 0; i < s; ++i) if(R_isnancpp(P[i])) return(true);
  return(false);
}

bool R_FINITE_VEC(double* P, int s) {
  for(int i = 0; i < s; ++i) if(!R_FINITE(P[i])) return(false);
  return(true);
}

double dwish_chol(double* x, double* chol, double df, int p, double scale_param, int give_log) {
  char uplo('U');
  char side('L');
  char diag('N');
  char transT('T');
  char transN('N');
  int info(0);
  double alpha(1.0);

  int i;

  if (R_IsNA(x, p*p) || R_IsNA(chol, p*p) || R_IsNA(df) || R_IsNA(scale_param))
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(x, p*p) || R_isnancpp(chol, p*p) || R_isnancpp(df) || R_isnancpp(scale_param))
    return R_NaN;
#endif

  // also covers df < 0
  if(df < (double) p) ML_ERR_return_NAN;

  if(!R_FINITE_VEC(x, p*p) || !R_FINITE_VEC(chol, p*p)) return R_D__0;

  double dens = -(df*p/2 * M_LN2 + p*(p-1)*M_LN_SQRT_PI/2);
  for(i = 0; i < p; i++)
    dens -= lgammafn((df - i) / 2);

  if(scale_param) {
    for(i = 0; i < p*p; i += p + 1) 
      dens -= df * log(chol[i]);
  } else {
    for(i = 0; i < p*p; i += p + 1) 
      dens += df * log(chol[i]);
  }

  // determinant of x using Cholesky:
  // dpotrf overwrites x so need to copy first since the trace calc below also overwrites
  double* tmp = new double[p*p];
  for(i = 0; i < p*p; i++) 
    tmp[i] = x[i];
  F77_CALL(dpotrf)(&uplo, &p, tmp, &p, &info);
  for(i = 0; i < p*p; i += p + 1) 
    dens += (df - p - 1) * log(tmp[i]);
  delete [] tmp;

  // do upper-triangular solves for scale parameterization
  // or upper-triangular multiplies for rate parameterization
  // dtr{m,s}m is a BLAS level-3 function
  double tmp_dens = 0.0;
  if(scale_param) {
    F77_CALL(dtrsm)(&side, &uplo, &transT, &diag, &p, &p, &alpha, 
           chol, &p, x, &p);
    F77_CALL(dtrsm)(&side, &uplo, &transN, &diag, &p, &p, &alpha, 
           chol, &p, x, &p);
    for(i = 0; i < p*p; i += p + 1) 
      tmp_dens += x[i];
    dens += -0.5 * tmp_dens;
  } else {
    // this could be improved by doing efficient U^T U multiply followed by direct product multiply
    F77_CALL(dtrmm)(&side, &uplo, &transN, &diag, &p, &p, &alpha, 
           chol, &p, x, &p);
    // at this point don't need to do all the multiplications so don't call dtrmm again
    // direct product of upper triangles
    for(int j = 0; j < p; j++) {
      for(i = 0; i <= j; i++) {
        tmp_dens += x[j*p+i] * chol[j*p+i];
      }
    }
    dens += -0.5 * tmp_dens;
  }

  return give_log ? dens : exp(dens);
}


SEXP C_dwish_chol(SEXP x, SEXP chol, SEXP df, SEXP scale_param, SEXP return_log) 
// calculates Wishart density given Cholesky of scale or rate matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!isReal(x) || !isReal(chol) || !isReal(df) || !isReal(scale_param) || !isLogical(return_log))
    RBREAK("Error (C_dwish_chol): invalid input type for one of the arguments.\n");
  int p = pow(LENGTH(chol), 0.5);
  int give_log = (int) LOGICAL(return_log)[0];
  double scale = REAL(scale_param)[0];

  double* c_x = REAL(x);
  double* c_chol = REAL(chol);
  double c_df = REAL(df)[0];
 
  if(c_df < p)
    RBREAK("Error (C_dwish_chol): inconsistent degrees of freedom and dimension.\n");
  
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, 1));  
  REAL(ans)[0] = dwish_chol(c_x, c_chol, c_df, p, scale, give_log);
  UNPROTECT(1);
  return ans;
}


void rwish_chol(double *Z, double* chol, double df, int p, double scale_param) {
  // Ok that this returns the array as return value? NIMBLE C code will need to free the memory
  char uplo('U');
  char sideL('L');
  //  char sideR('R');
  char diag('N');
  char transT('T');
  char transN('N');
  //  int info(0);
  double alpha(1.0);
  double beta(0.0);

  int i, j, uind, lind;
  
#ifdef IEEE_754
  if (R_isnancpp(chol, p*p) || R_isnancpp(df) || R_isnancpp(scale_param)) {
    for(j = 0; j < p*p; j++) 
      Z[j] = R_NaN;
    return;
  }
#endif
    
  // also covers df < 0
  if(df < (double) p) {
    for(j = 0; j < p*p; j++) 
      Z[j] = R_NaN;
    return;
  }
  
  // fill diags with sqrts of chi-squares and upper triangle (for scale_param) with std normals
  for(j = 0; j < p; j++) {
    // double *Z_j = &Z[j*p];
    //Z_j[j] = sqrt(rchisq(df - (double) j)); 
    Z[j*p + j] = sqrt(rchisq(df - (double) j)); 
    for(i = 0; i < j; i++) {
      uind = i + j * p, /* upper triangle index */
      lind = j + i * p; /* lower triangle index */
      Z[(scale_param ? uind : lind)] = norm_rand();
      Z[(scale_param ? lind : uind)] = 0;
    }
    /*
    for(i = 0; i < j; i++) 
      Z_j[i] = norm_rand();
    for (i = j + 1; i < p; i++)
      Z_j[i] = 0;
    */
  }
 
  // multiply Z*chol, both upper triangular or solve(chol, Z^T)
  // would be more efficient if make use of fact that right-most matrix is triangular
  if(scale_param) F77_CALL(dtrmm)(&sideL, &uplo, &transN, &diag, &p, &p, &alpha, Z, &p, chol, &p);
  else F77_CALL(dtrsm)(&sideL, &uplo, &transN, &diag, &p, &p, &alpha, chol, &p, Z, &p);

  // cp result to Z or chol so can be used as matrix to multiply against and overwrite
  if(scale_param) {
    for(j = 0; j < p*p; j++) 
      Z[j] = chol[j];
  } else {
    for(j = 0; j < p*p; j++) 
      chol[j] = Z[j]; // FIXME: check this case.  We don't want to overwite an input argument.
  }

  // do crossprod of result; again this would be more efficient use fact that t(input) is lower-tri
  if(scale_param) F77_CALL(dtrmm)(&sideL, &uplo, &transT, &diag, &p, &p, &alpha, chol, &p, Z, &p);
  else F77_CALL(dgemm)(&transN, &transT, &p, &p, &p, &alpha, chol, &p, chol, &p, &beta, Z, &p); 

}

SEXP C_rwish_chol(SEXP chol, SEXP df, SEXP scale_param) 
// generates single Wishart draw given Cholesky of scale or rate matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!isReal(chol) || !isReal(df) || !isReal(scale_param))
    RBREAK("Error (C_rwish_chol): invalid input type for one of the arguments.\n");
  int n_chol = LENGTH(chol);
  int p = pow(n_chol, 0.5);
  double scale = REAL(scale_param)[0];

  double* c_chol = REAL(chol);
  double c_df = REAL(df)[0];

  if(c_df < p)
    RBREAK("Error (C_rwish_chol): inconsistent degrees of freedom and dimension.\n");

  GetRNGstate(); 

  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, n_chol));  
  rwish_chol(REAL(ans), c_chol, c_df, p, scale);
  
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}



double ddirch(double* x, double* alpha, int K, int give_log) 
// scalar function that can be called directly by NIMBLE with same name as in R
{
  double sumAlpha(0.0);
  double sumX(0.0);
  double dens(0.0);

  if (R_IsNA(x, K) || R_IsNA(alpha, K))
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(x, K) || R_isnancpp(alpha, K))
    return R_NaN;
#endif
  
  for(int i = 0; i < K; i++) {
    if(alpha[i] <= 0.0) ML_ERR_return_NAN;
    if(x[i] < 0.0 || x[i] > 1.0) return R_D__0;
    dens += (alpha[i]-1) * log(x[i]) - lgammafn(alpha[i]) ;
    sumAlpha += alpha[i];
    sumX += x[i];
  }
  if(sumX > 1.0 + 10*DBL_EPSILON || sumX < 1.0 - 10*DBL_EPSILON) {
    return R_D__0;
  }

  dens += lgammafn(sumAlpha);
  return give_log ? dens : exp(dens);
}

void rdirch(double *ans, double* alpha, int K) 
// scalar function that can be called directly by NIMBLE with same name as in R
{
  int i, j;
  //  double* ans = new double[K];
#ifdef IEEE_754
  if (R_isnancpp(alpha, K)) {
    for(j = 0; j < K; j++) 
      ans[j] = R_NaN;
    return;
  }
#endif

  double sum(0.0);
  for(i = 0; i < K; i++) {
    if(alpha[i] <= 0.0) {
      for(j = 0; j < K; j++) 
        ans[j] = R_NaN;
      return;
    }
    ans[i] = rgamma(alpha[i], 1);
    sum += ans[i];
  }
  for(i = 0; i < K; i++) {
    ans[i] /= sum;
  }
}

SEXP C_ddirch(SEXP x, SEXP alpha, SEXP return_log) 
{
  if(!isReal(x) || !isReal(alpha) || !isLogical(return_log)) 
    RBREAK("Error (C_ddirch): invalid input type for one of the arguments.\n");
  int K = LENGTH(alpha);
  if(LENGTH(x) != K)
    RBREAK("Error (C_ddirch): length of x must equal length of alpha.\n")
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(K == 0) {
    return alpha;
  }

  double* c_x = REAL(x);
  double* c_alpha = REAL(alpha);

  PROTECT(ans = allocVector(REALSXP, 1));  
  REAL(ans)[0] = ddirch(c_x, c_alpha, K, give_log);

  UNPROTECT(1);
  return ans;
}
 

SEXP C_rdirch(SEXP alpha) {
  if(!isReal(alpha))
    RBREAK("Error (C_rdirch): invalid input type for the argument.\n");
  int K = LENGTH(alpha);

  SEXP ans;

  if(K == 0) {
    PROTECT(ans = allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }

  double* c_alpha = REAL(alpha);

  GetRNGstate(); 

  PROTECT(ans = allocVector(REALSXP, K));  
  rdirch(REAL(ans), c_alpha, K);
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}



double dmulti(double* x, double size, double* prob, int K, int give_log) // Calling functions need to copy first arg to int if needed
// scalar function that can be called directly by NIMBLE with same name as in R
{
  double sumProb(0.0);
  double sumX(0.0);

  if (R_IsNA(x, K) || R_IsNA(prob, K) || R_IsNA(size))
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(x, K) || R_isnancpp(prob, K) || R_isnancpp(size))
    return R_NaN;
#endif

  if(R_D_negInonint(size))
        ML_ERR_return_NAN;
  size = R_D_forceint(size);

  double dens = lgammafn(size + 1);
  for(int i = 0; i < K; i++) {
    if (prob[i] < 0 || prob[i] > 1) ML_ERR_return_NAN;
    R_D_nonint_check(x[i]);
    if (x[i] < 0 || !R_FINITE(x[i])) return R_D__0;

    x[i] = R_D_forceint(x[i]);
    sumProb += prob[i];
    sumX += x[i];

    if(!(x[i] == 0.0 && prob[i] == 0.0))
      dens += x[i]*log(prob[i]) - lgammafn(x[i] + 1);
  }

  if(sumX > size + 10*DBL_EPSILON || sumX < size - 10*DBL_EPSILON) {
    return R_D__0;
  }
  if(sumProb > 1.0 + 10*DBL_EPSILON || sumProb < 1.0 - 10*DBL_EPSILON) {
    ML_ERR_return_NAN;
  }

  return give_log ? dens : exp(dens);
}




void rmulti(int *ans, double size, double* prob, int K) // Calling functions need to copy first arg back and forth to double if needed
// scalar function that can be called directly by NIMBLE with same name as in R
// just call Rmath's rmultinom, which passes result by pointer
// IMPORTANT: have ans and size as int when sent to rmultinom as Rmath rmultinom has these types
// Nimble does a copy in nimArr_rmulti
{
  rmultinom((int) size, prob, K, ans);
}

SEXP C_dmulti(SEXP x, SEXP size, SEXP prob, SEXP return_log) 
{
  if(!isReal(x) || !isReal(size) || !isReal(prob) || !isLogical(return_log)) 
    RBREAK("Error (C_dmulti): invalid input type for one of the arguments.\n");
  int K = LENGTH(prob);
  if(LENGTH(x) != K)
    RBREAK("Error (C_dmulti): length of x must equal size.\n")
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  int i;

  if(K == 0) {
    return prob;
  }

  double* c_x = REAL(x);
  double* c_prob = REAL(prob);
  double c_size = REAL(size)[0];

  double sum = 0.0;
  for(i = 0; i < K; i++) 
    sum += c_prob[i];
  if(sum > 1.0 + 10*DBL_EPSILON || sum < 1.0 - 10*DBL_EPSILON)
    RBREAK("Error (C_dmulti): sum of probabilities is not equal to 1.\n");

  sum = 0.0;
  for(i = 0; i < K; i++) 
    sum += c_x[i];
  if(c_size > sum + 10*DBL_EPSILON || c_size < sum - 10*DBL_EPSILON)
    RBREAK("Error (C_dmulti): sum of values is not equal to size.\n");

  PROTECT(ans = allocVector(REALSXP, 1));  
  REAL(ans)[0] = dmulti(c_x, c_size, c_prob, K, give_log);

  UNPROTECT(1);
  return ans;
}
 

SEXP C_rmulti(SEXP size, SEXP prob) {
  if(!isReal(size) || !isReal(prob))
    RBREAK("Error (C_rmulti): invalid input type for one of the arguments.\n");
  int K = LENGTH(prob);

  SEXP ans;
  int i;

  if(K == 0) {
    PROTECT(ans = allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }

  double* c_prob = REAL(prob);
  double c_size = REAL(size)[0];

  double sum = 0.0;
  for(i = 0; i < K; i++) 
    sum += c_prob[i];
  if(sum > 1.0 + 10*DBL_EPSILON || sum < 1.0 - 10*DBL_EPSILON)
    RBREAK("Error (C_rmulti): sum of probabilities is not equal to 1.\n");

  GetRNGstate(); 

  PROTECT(ans = allocVector(INTSXP, K));  
  rmulti(INTEGER(ans), c_size, c_prob, K);
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}


double dcat(double x, double* prob, int K, int give_log)
// scalar function that can be called directly by NIMBLE with same name as in R
{
  if (R_IsNA(x) || R_IsNA(prob, K))
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(x) || R_isnancpp(prob, K))
    return R_NaN;
#endif

  double sumProb(0.0);
  for(int i = 0; i < K; i++) 
    sumProb += prob[i];
  if(sumProb > 1.0 + 10*DBL_EPSILON || sumProb < 1.0 - 10*DBL_EPSILON) {
    ML_ERR_return_NAN;
  }

  R_D_nonint_check(x);
  x = R_D_forceint(x);

  if(x > K || x < 1) return R_D__0;
  return give_log ? log(prob[(int) x - 1]) : prob[(int) x - 1];
}

double rcat(double* prob, int K)
// scalar function that can be called directly by NIMBLE with same name as in R
// problem is no apparent way to call this w/o also passing the number of categories
// relying on sum to 1 risks accessing beyond storage
// we'll need to figure out how to inject the number of categories w/in NIMBLE
{
  double u = unif_rand();
  double prob_cum = prob[0];
  int value = 1;

  double sumProb(0.0);
  for(int i = 0; i < K; i++) 
    sumProb += prob[i];
  if(sumProb > 1.0 + 10*DBL_EPSILON || sumProb < 1.0 - 10*DBL_EPSILON) {
    ML_ERR_return_NAN;
  }

  while(u > prob_cum && value < K) {
    prob_cum += prob[value];
    value++;
  }

  return (double) value;
}
 
SEXP C_dcat(SEXP x, SEXP prob, SEXP return_log) 
{
  // this will call NIMBLE's dcat() for computation on scalars
  // prob must be a single vector of probs adding to one, but x can be a vector

  if(!isReal(x) || !isReal(prob) || !isLogical(return_log)) 
    RBREAK("Error (C_dcat): invalid input type for one of the arguments.\n");
  int n_x = LENGTH(x);
  int K = LENGTH(prob);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  int i;

  if(n_x == 0) {
    return x;
  }

  double* c_x = REAL(x);
  double* c_prob = REAL(prob);

  double sum = 0.0;
  for(i = 0; i < K; i++) 
    sum += c_prob[i];
  if(sum > 1.0 + 10*DBL_EPSILON || sum < 1.0 - 10*DBL_EPSILON)
    RBREAK("Error (C_dcat): sum of probabilities is not equal to 1.\n");

  PROTECT(ans = allocVector(REALSXP, n_x));  
  for(i = 0; i < n_x; i++) {
    REAL(ans)[i] = dcat(c_x[i], c_prob, K, give_log);
  }

  UNPROTECT(1);
  return ans;
}
  
SEXP C_rcat(SEXP n, SEXP prob) {
  if(!isInteger(n) || !isReal(prob))
    RBREAK("Error (C_rcat): invalid input type for one of the arguments.\n");
  int n_values = INTEGER(n)[0];
  int K = LENGTH(prob);

  SEXP ans;
  int i;

  if(n_values == 0) {
    PROTECT(ans = allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    RBREAK("Error (C_rcat): n must be non-negative.\n");

  double* c_prob = REAL(prob);

  double sum = 0.0;
  for(i = 0; i < K; i++) 
    sum += c_prob[i];
  if(sum > 1.0 + 10*DBL_EPSILON || sum < 1.0 - 10*DBL_EPSILON)
    RBREAK("Error (C_rcat): sum of probabilities is not equal to 1.\n");

  GetRNGstate(); 
  PROTECT(ans = allocVector(INTSXP, n_values));  

  for(i = 0; i < n_values; i++) 
    INTEGER(ans)[i] = (int) rcat(c_prob, K);

  PutRNGstate();
  UNPROTECT(1);
  return ans;

}
  
double dmnorm_chol(double* x, double* mean, double* chol, int n, double prec_param, int give_log) {
  char uplo('U');
  char transPrec('N');
  char transCov('T');
  char diag('N');
  int lda(n);
  int incx(1);

  double dens = -n * M_LN_SQRT_2PI;
  int i;
  // add diagonals of Cholesky

  if (R_IsNA(x, n) || R_IsNA(mean, n) || R_IsNA(chol, n*n) || R_IsNA(prec_param))
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(x, n) || R_isnancpp(mean, n) || R_isnancpp(chol, n*n) || R_isnancpp(prec_param))
    return R_NaN;
#endif

  if(!R_FINITE_VEC(x, n) || !R_FINITE_VEC(mean, n) || !R_FINITE_VEC(chol, n*n)) return R_D__0;


  if(prec_param) {
    for(i = 0; i < n*n; i += n + 1) 
      dens += log(chol[i]);
  } else {
    for(i = 0; i < n*n; i += n + 1) 
      dens -= log(chol[i]);
  }
  for(i = 0; i < n; i++) 
    x[i] -= mean[i];

  // do matrix-vector multiply with upper-triangular matrix stored column-wise as full n x n matrix (prec parameterization)
  // or upper-triangular (transpose) solve (cov parameterization)
  // dtr{m,s}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrmv)(&uplo, &transPrec, &diag, &n, chol, &lda, x, &incx);
  else F77_CALL(dtrsv)(&uplo, &transCov, &diag, &n, chol, &lda, x, &incx);

  // sum of squares to calculate quadratic form
  double tmp = 0.0;
  for(i = 0; i < n; i++)
    tmp += x[i] * x[i];

  dens += -0.5 * tmp;

  return give_log ? dens : exp(dens);
}

SEXP C_dmnorm_chol(SEXP x, SEXP mean, SEXP chol, SEXP prec_param, SEXP return_log) 
// calculates mv normal density given Cholesky of precision matrix or covariance matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!isReal(x) || !isReal(mean) || !isReal(chol) || !isReal(prec_param) || !isLogical(return_log))
    RBREAK("Error (C_dmnorm_chol): invalid input type for one of the arguments.\n");
  int n_x = LENGTH(x);
  int n_mean = LENGTH(mean);
  int give_log = (int) LOGICAL(return_log)[0];
  double prec = REAL(prec_param)[0];

  double* c_x = REAL(x);
  double* c_mean = REAL(mean);
  double* c_chol = REAL(chol);

  double* xcopy = new double[n_x];
  for(int i = 0; i < n_x; ++i) 
    xcopy[i] = c_x[i];

  double* full_mean;
  if(n_mean < n_x) {
    full_mean = new double[n_x];
    int i_mean = 0;
    for(int i = 0; i < n_x; i++) {
      full_mean[i] = c_mean[i_mean++];
      if(i_mean == n_mean) i_mean = 0;
    }
  } else full_mean = c_mean;
  
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, 1));  
  REAL(ans)[0] = dmnorm_chol(xcopy, full_mean, c_chol, n_x, prec, give_log);
  if(n_mean < n_x)
    delete [] full_mean;
  delete [] xcopy;
  UNPROTECT(1);
  return ans;
}

void rmnorm_chol(double *ans, double* mean, double* chol, int n, double prec_param) {
  // Ok that this returns the array as return value? NIMBLE C code will need to free the memory
  char uplo('U');
  char transPrec('N');
  char transCov('T');
  char diag('N');
  int lda(n);
  int incx(1);
  
  int i, j;

#ifdef IEEE_754
  if (R_isnancpp(mean, n) || R_isnancpp(chol, n*n) || R_isnancpp(prec_param)) {
    for(j = 0; j < n; j++) 
      ans[j] = R_NaN;
    return;
  }
#endif

  if(!R_FINITE_VEC(chol, n*n)) { 
    for(j = 0; j < n; j++) 
      ans[j] = R_NaN;
    return;
  }

  double* devs = new double[n];
  //  double* ans = new double[n];

  for(i = 0; i < n; i++) 
    devs[i] = norm_rand();


  // do upper-triangular solve or (transpose) multiply
  // dtr{s,m}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrsv)(&uplo, &transPrec, &diag, &n, chol, &lda, devs, &incx);
  else F77_CALL(dtrmv)(&uplo, &transCov, &diag, &n, chol, &lda, devs, &incx);

  for(i = 0; i < n; i++) 
    ans[i] = mean[i] + devs[i];

  delete [] devs;
}

SEXP C_rmnorm_chol(SEXP mean, SEXP chol, SEXP prec_param) 
// generates single mv normal draw given Cholesky of precision matrix or covariance matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!isReal(mean) || !isReal(chol) || !isReal(prec_param))
    RBREAK("Error (C_rmnorm_chol): invalid input type for one of the arguments.\n");
  int n_mean = LENGTH(mean);
  int n_chol = LENGTH(chol);
  int n_values = pow(n_chol, 0.5);
  double prec = REAL(prec_param)[0];

  int i;

  double* c_mean = REAL(mean);
  double* c_chol = REAL(chol);
  double* full_mean; 

  if(n_mean < n_values) {
    full_mean = new double[n_values];
    int i_mean = 0;
    for(i = 0; i < n_values; i++) {
      full_mean[i] = c_mean[i_mean++];
      if(i_mean == n_mean) i_mean = 0;
    }
  } else full_mean = c_mean;

  GetRNGstate(); 

  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, n_values));  
  rmnorm_chol(REAL(ans), full_mean, c_chol, n_values, prec);

  PutRNGstate();
  if(n_mean < n_values) 
    delete [] full_mean;
  UNPROTECT(1);
  return ans;
}

// Begin multivariate t

double dmvt_chol(double* x, double* mean, double* chol, double df, int n, double prec_param, int give_log) {
  char uplo('U');
  char transPrec('N');
  char transCov('T');
  char diag('N');
  int lda(n);
  int incx(1);
  
  double dens = lgammafn((df + n) / 2) - lgammafn(df / 2) - n * M_LN_SQRT_PI - n * log(df) / 2;
  int i;
  // add diagonals of Cholesky
  
  if (R_IsNA(x, n) || R_IsNA(mean, n) || R_IsNA(chol, n*n) || R_IsNA(df) || R_IsNA(prec_param))
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(x, n) || R_isnancpp(mean, n) || R_isnancpp(chol, n*n) || R_IsNA(df) || R_isnancpp(prec_param))
    return R_NaN;
#endif
  
  if(!R_FINITE_VEC(x, n) || !R_FINITE_VEC(mean, n) || !R_FINITE_VEC(chol, n*n)) return R_D__0;
  
  
  if(prec_param) {
    for(i = 0; i < n*n; i += n + 1) 
      dens += log(chol[i]);
  } else {
    for(i = 0; i < n*n; i += n + 1) 
      dens -= log(chol[i]);
  }
  for(i = 0; i < n; i++) 
    x[i] -= mean[i];
  
  // do matrix-vector multiply with upper-triangular matrix stored column-wise as full n x n matrix (prec parameterization)
  // or upper-triangular (transpose) solve (cov parameterization)
  // dtr{m,s}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrmv)(&uplo, &transPrec, &diag, &n, chol, &lda, x, &incx);
  else F77_CALL(dtrsv)(&uplo, &transCov, &diag, &n, chol, &lda, x, &incx);
  
  // sum of squares to calculate quadratic form
  double tmp = 0.0;
  for(i = 0; i < n; i++)
    tmp += x[i] * x[i];
  
  dens += -0.5 * (df + n) * log(1 + tmp / df);
  
  return give_log ? dens : exp(dens);
}

SEXP C_dmvt_chol(SEXP x, SEXP mean, SEXP chol, SEXP df, SEXP prec_param, SEXP return_log) 
  // calculates mv normal density given Cholesky of precision matrix or covariance matrix
  // Cholesky matrix should be given as a numeric vector in column-major order
  //   including all n x n elements; lower-triangular elements are ignored
{
  if(!isReal(x) || !isReal(mean) || !isReal(chol) || !isReal(df) || !isReal(prec_param) || !isLogical(return_log))
    RBREAK("Error (C_dmvt_chol): invalid input type for one of the arguments.\n");
  int n_x = LENGTH(x);
  int n_mean = LENGTH(mean);
  int give_log = (int) LOGICAL(return_log)[0];
  double c_df = REAL(df)[0];
  double prec = REAL(prec_param)[0];
  
  double* c_x = REAL(x);
  double* c_mean = REAL(mean);
  double* c_chol = REAL(chol);
  
  double* xcopy = new double[n_x];
  for(int i = 0; i < n_x; ++i) 
    xcopy[i] = c_x[i];
  
  double* full_mean;
  if(n_mean < n_x) {
    full_mean = new double[n_x];
    int i_mean = 0;
    for(int i = 0; i < n_x; i++) {
      full_mean[i] = c_mean[i_mean++];
      if(i_mean == n_mean) i_mean = 0;
    }
  } else full_mean = c_mean;
  
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, 1));  
  REAL(ans)[0] = dmvt_chol(xcopy, full_mean, c_chol, c_df, n_x, prec, give_log);
  if(n_mean < n_x)
    delete [] full_mean;
  delete [] xcopy;
  UNPROTECT(1);
  return ans;
}

void rmvt_chol(double *ans, double* mean, double* chol, double df, int n, double prec_param) {
  // Ok that this returns the array as return value? NIMBLE C code will need to free the memory
  char uplo('U');
  char transPrec('N');
  char transCov('T');
  char diag('N');
  int lda(n);
  int incx(1);
  
  int i, j;
  
#ifdef IEEE_754
  if (R_isnancpp(mean, n) || R_isnancpp(chol, n*n) || R_isnancpp(df) || R_isnancpp(prec_param)) {
    for(j = 0; j < n; j++) 
      ans[j] = R_NaN;
    return;
  }
#endif
  
  if(!R_FINITE_VEC(chol, n*n)) { 
    for(j = 0; j < n; j++) 
      ans[j] = R_NaN;
    return;
  }
  
  double* devs = new double[n];
  //  double* ans = new double[n];
  
  for(i = 0; i < n; i++) 
    devs[i] = norm_rand();
  
  // sample from chi-squared and calculate scaling factor
  double scaling = sqrt(df / rchisq(df));
  
  // do upper-triangular solve or (transpose) multiply
  // dtr{s,m}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrsv)(&uplo, &transPrec, &diag, &n, chol, &lda, devs, &incx);
  else F77_CALL(dtrmv)(&uplo, &transCov, &diag, &n, chol, &lda, devs, &incx);
  
  for(i = 0; i < n; i++) 
    ans[i] = mean[i] + devs[i] * scaling;
  
  delete [] devs;
}

SEXP C_rmvt_chol(SEXP mean, SEXP chol, SEXP df, SEXP prec_param) 
  // generates single mv normal draw given Cholesky of precision matrix or covariance matrix
  // Cholesky matrix should be given as a numeric vector in column-major order
  //   including all n x n elements; lower-triangular elements are ignored
{
  if(!isReal(mean) || !isReal(chol) || !isReal(df) || !isReal(prec_param))
    RBREAK("Error (C_rmvt_chol): invalid input type for one of the arguments.\n");
  int n_mean = LENGTH(mean);
  int n_chol = LENGTH(chol);
  int n_values = pow(n_chol, 0.5);
  double c_df = REAL(df)[0];
  double prec = REAL(prec_param)[0];
  
  int i;
  
  double* c_mean = REAL(mean);
  double* c_chol = REAL(chol);
  double* full_mean; 
  
  if(n_mean < n_values) {
    full_mean = new double[n_values];
    int i_mean = 0;
    for(i = 0; i < n_values; i++) {
      full_mean[i] = c_mean[i_mean++];
      if(i_mean == n_mean) i_mean = 0;
    }
  } else full_mean = c_mean;
  
  GetRNGstate(); 
  
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, n_values));  
  rmvt_chol(REAL(ans), full_mean, c_chol, c_df, n_values, prec);
  
  PutRNGstate();
  if(n_mean < n_values) 
    delete [] full_mean;
  UNPROTECT(1);
  return ans;
}


double dt_nonstandard(double x, double df, double mu, double sigma, int give_log)
// scalar function that can be called directly by NIMBLE with same name as in R
// 'n' is name of 'df' in R's C dt
{
  // standardize and multiply by Jacobian of transformation
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(df))
    return x + mu + sigma + df;
#endif
  if(!R_FINITE(sigma)) return R_D__0;
  if(sigma <= 0.0) {
    if(sigma < 0.0) ML_ERR_return_NAN;
    return (x == mu) ? ML_POSINF : R_D__0;
  }         
  if(give_log) return dt( (x - mu)/sigma, df, give_log) - log(sigma);
  else return dt( (x - mu)/sigma, df, give_log) / sigma;
}

double rt_nonstandard(double df, double mu, double sigma)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(df))
    ML_ERR_return_NAN;
#endif
  if (!R_FINITE(sigma) || sigma < 0.0) ML_ERR_return_NAN;
    
  return mu + sigma * rt(df);
}

double pt_nonstandard(double q, double df, double mu, double sigma, int lower_tail, int log_p)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
    if(ISNAN(q) || ISNAN(mu) || ISNAN(sigma) || ISNAN(df))
        return q + mu + sigma + df;
#endif

  if(!R_FINITE(q) && mu == q) return ML_NAN; /* q-mu is NaN */
  if(sigma <= 0.0) {
    if(sigma < 0.0) ML_ERR_return_NAN;
    return (q < mu) ? R_DT_0 : R_DT_1;
  }   
   
  return( pt( (q - mu) / sigma, df, lower_tail, log_p) );
}

double qt_nonstandard(double p, double df, double mu, double sigma, int lower_tail, int log_p)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma) || ISNAN(df))
    return p + mu + sigma + df;
#endif

  if(sigma < 0.0) ML_ERR_return_NAN;
  if(sigma == 0.0) return mu;

  return( mu + sigma * qt( p, df, lower_tail, log_p ) );
}


SEXP C_dt_nonstandard(SEXP x, SEXP df, SEXP mu, SEXP sigma, SEXP return_log) {
  if(!isReal(x) || !isReal(df) || !isReal(mu) || !isReal(sigma) || !isLogical(return_log)) 
    RBREAK("Error (C_dt_nonstandard): invalid input type for one of the arguments.");
  int n_x = LENGTH(x);
  int n_mu = LENGTH(mu);
  int n_sigma = LENGTH(sigma);
  int n_df = LENGTH(df);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(n_x == 0) {
    return x;
  }
    
  PROTECT(ans = allocVector(REALSXP, n_x));  
  double* c_x = REAL(x);
  double* c_mu = REAL(mu);
  double* c_sigma = REAL(sigma);
  double* c_df = REAL(df);

  // FIXME: abstract the recycling as a function
  if(n_mu == 1 && n_sigma == 1 && n_df == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_x; i++) 
      REAL(ans)[i] = dt_nonstandard(c_x[i], *c_df, *c_mu, *c_sigma, give_log);
  } else {
    int i_mu = 0;
    int i_sigma = 0;
    int i_df = 0;
    for(int i = 0; i < n_x; i++) {
      REAL(ans)[i] = dt_nonstandard(c_x[i], c_df[i_df++], c_mu[i_mu++], c_sigma[i_sigma++], give_log);
      //c_mu[i_mu++] + c_sigma[i_sigma++] * rt(c_df[i_df++]);
      // implement recycling:
      if(i_mu == n_mu) i_mu = 0;
      if(i_sigma == n_sigma) i_sigma = 0;
      if(i_df == n_df) i_df = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
  
SEXP C_rt_nonstandard(SEXP n, SEXP df, SEXP mu, SEXP sigma) {
  // this will call R's rt() for computation on scalars
  if(!isInteger(n) || !isReal(df) || !isReal(mu) || !isReal(sigma))
    RBREAK("Error (C_rt_nonstandard): invalid input type for one of the arguments.");
  int n_mu = LENGTH(mu);
  int n_sigma = LENGTH(sigma);
  int n_df = LENGTH(df);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = allocVector(REALSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rt_nonstandard): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = allocVector(REALSXP, n_values));  
  double* c_mu = REAL(mu);
  double* c_sigma = REAL(sigma);
  double* c_df = REAL(df);
  if(n_mu == 1 && n_sigma == 1 && n_df == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_values; i++) 
      REAL(ans)[i] = rt_nonstandard(*c_df, *c_mu, *c_sigma);
  } else {
    int i_mu = 0;
    int i_sigma = 0;
    int i_df = 0;
    for(int i = 0; i < n_values; i++) {
      REAL(ans)[i] = rt_nonstandard(c_df[i_df++], c_mu[i_mu++], c_sigma[i_sigma++]);
      //c_mu[i_mu++] + c_sigma[i_sigma++] * rt(c_df[i_df++]);
      // implement recycling:
      if(i_mu == n_mu) i_mu = 0;
      if(i_sigma == n_sigma) i_sigma = 0;
      if(i_df == n_df) i_df = 0;
    }
  }
    
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}
  
SEXP C_pt_nonstandard(SEXP q, SEXP df, SEXP mu, SEXP sigma, SEXP lower_tail, SEXP log_p) {
  if(!isReal(q) || !isReal(df) || !isReal(mu) || !isReal(sigma) || !isLogical(lower_tail) || !isLogical(log_p))
    RBREAK("Error (C_pt_nonstandard): invalid input type for one of the arguments.");
  int n_q = LENGTH(q);
  int n_mu = LENGTH(mu);
  int n_sigma = LENGTH(sigma);
  int n_df = LENGTH(df);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_q == 0) {
    return q;
  }
    
  PROTECT(ans = allocVector(REALSXP, n_q));  
  double* c_q = REAL(q);
  double* c_mu = REAL(mu);
  double* c_sigma = REAL(sigma);
  double* c_df = REAL(df);

  // FIXME: abstract the recycling as a function
  if(n_mu == 1 && n_sigma == 1 && n_df == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_q; i++) 
      REAL(ans)[i] = pt_nonstandard(c_q[i], *c_df, *c_mu, *c_sigma, c_lower_tail, c_log_p);
  } else {
    int i_mu = 0;
    int i_sigma = 0;
    int i_df = 0;
    for(int i = 0; i < n_q; i++) {
      REAL(ans)[i] = pt_nonstandard(c_q[i], c_df[i_df++], c_mu[i_mu++], c_sigma[i_sigma++], c_lower_tail, c_log_p);
      //c_mu[i_mu++] + c_sigma[i_sigma++] * rt(c_df[i_df++]);
      // implement recycling:
      if(i_mu == n_mu) i_mu = 0;
      if(i_sigma == n_sigma) i_sigma = 0;
      if(i_df == n_df) i_df = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
 
SEXP C_qt_nonstandard(SEXP p, SEXP df, SEXP mu, SEXP sigma, SEXP lower_tail, SEXP log_p) {
  if(!isReal(p) || !isReal(df) || !isReal(mu) || !isReal(sigma) || !isLogical(lower_tail) || !isLogical(log_p))
    RBREAK("Error (C_qt_nonstandard): invalid input type for one of the arguments.");
  int n_p = LENGTH(p);
  int n_mu = LENGTH(mu);
  int n_sigma = LENGTH(sigma);
  int n_df = LENGTH(df);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_p == 0) {
    return p;
  }
    
  PROTECT(ans = allocVector(REALSXP, n_p));  
  double* c_p = REAL(p);
  double* c_mu = REAL(mu);
  double* c_sigma = REAL(sigma);
  double* c_df = REAL(df);

  // FIXME: abstract the recycling as a function
  if(n_mu == 1 && n_sigma == 1 && n_df == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_p; i++) 
      REAL(ans)[i] = qt_nonstandard(c_p[i], *c_df, *c_mu, *c_sigma, c_lower_tail, c_log_p);
  } else {
    int i_mu = 0;
    int i_sigma = 0;
    int i_df = 0;
    for(int i = 0; i < n_p; i++) {
      REAL(ans)[i] = qt_nonstandard(c_p[i], c_df[i_df++], c_mu[i_mu++], c_sigma[i_sigma++], c_lower_tail, c_log_p);
      //c_mu[i_mu++] + c_sigma[i_sigma++] * rt(c_df[i_df++]);
      // implement recycling:
      if(i_mu == n_mu) i_mu = 0;
      if(i_sigma == n_sigma) i_sigma = 0;
      if(i_df == n_df) i_df = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}


double dinterval(double x, double t, double* c, int K, int give_log)
// scalar function that can be called directly by NIMBLE with same name as in R
{
  if (R_IsNA(c, K) || R_IsNA(x) || R_IsNA(t)) 
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(c, K) || R_isnancpp(x) || R_isnancpp(t))
    return R_NaN;
#endif

  R_D_nonint_check(x);
  x = R_D_forceint(x);

  // we do not check that c is in increasing order, to save time
  int int_x = (int) x;
  if(int_x < 0 || int_x > K) return R_D__0; 
  if(int_x == 0 && t <= c[int_x]) return R_D__1;
  if(int_x == K && t > c[int_x-1]) return R_D__1;
  if(t <= c[int_x] && t > c[int_x - 1]) return R_D__1;
  return R_D__0;
}


double rinterval(double t, double* c, int K)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (R_isnancpp(c, K) || R_isnancpp(t) )
    ML_ERR_return_NAN;
#endif

  // we do not check that c is in increasing order, to save time
  for(int i = 0; i < K; i++) {
    if(t <= c[i]) return (double) i;
  }
  return (double) K;
}


SEXP C_dinterval(SEXP x, SEXP t, SEXP c, SEXP return_log) {
  // this will call NIMBLE's dinterval() for computation on scalars
  // c must be a single vector of cutpoints; x and t can be vectors
  if(!isReal(x) || !isReal(t) || !isReal(c) || !isLogical(return_log)) 
    RBREAK("Error (C_dinterval): invalid input type for one of the arguments.");
  int n_x = LENGTH(x);
  int n_t = LENGTH(t);
  int n_c = LENGTH(c);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(n_x == 0) {
    return x;
  }
    
  PROTECT(ans = allocVector(REALSXP, n_x));  
  double* c_x = REAL(x);
  double* c_t = REAL(t);
  double* c_c = REAL(c);

  // FIXME: abstract the recycling as a function
  if(n_t == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_x; i++) 
      REAL(ans)[i] = dinterval(c_x[i], *c_t, c_c, n_c, give_log);
  } else {
    int i_t = 0;
    for(int i = 0; i < n_x; i++) {
      REAL(ans)[i] = dinterval(c_x[i], c_t[i_t++], c_c, n_c, give_log);
      // implement recycling:
      if(i_t == n_t) i_t = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
  
SEXP C_rinterval(SEXP n, SEXP t, SEXP c) {
  if(!isInteger(n) || !isReal(t) || !isReal(c))
    RBREAK("Error (C_rinterval): invalid input type for one of the arguments.");
  int n_t = LENGTH(t);
  int K = LENGTH(c);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rinterval): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = allocVector(INTSXP, n_values));  
  double* c_t = REAL(t);
  double* c_c = REAL(c);
  if(n_t == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_values; i++) 
      INTEGER(ans)[i] = (int) rinterval(*c_t, c_c, K);
  } else {
    int i_t = 0;
    for(int i = 0; i < n_values; i++) {
      INTEGER(ans)[i] = (int) rinterval(c_t[i_t++], c_c, K);
      // implement recycling:
      if(i_t == n_t) i_t = 0;
    }
  }
    
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}


double dconstraint(double x, double cond, int give_log)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(cond))
    return x + cond;
#endif
  // any non-zero will be treated as true
  if(x == cond || x == 0.0) return R_D__1;
  else return R_D__0;
}


double rconstraint(double cond)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(cond) )
    ML_ERR_return_NAN;
#endif
  return cond;
}

// we need our own exp implementation because R dexp uses rate and C exp uses scale
  
double rexp_nimble(double rate)
{
  return rexp( 1/rate );
} 

double dexp_nimble(double x, double rate, int give_log)
{
  return dexp(x, 1/rate, give_log); 
} 

double pexp_nimble(double q, double rate, int lower_tail, int log_p)
{
  return pexp(q, 1/rate, lower_tail, log_p);
} 

double qexp_nimble(double p, double rate, int lower_tail, int log_p)
{
  return qexp(p, 1 / rate, lower_tail, log_p);
}

SEXP C_dexp_nimble(SEXP x, SEXP rate, SEXP return_log) {
  if(!isReal(x) || !isReal(rate) || !isLogical(return_log))
    RBREAK("Error (C_dexp_nimble): invalid input type for one of the arguments.");
  int n_x = LENGTH(x);
  int n_rate = LENGTH(rate);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(n_x == 0) {
    return x;
  }
    
  PROTECT(ans = allocVector(REALSXP, n_x));  
  double* c_x = REAL(x);
  double* c_rate = REAL(rate);

  // FIXME: abstract the recycling as a function
  if(n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_x; i++) 
      REAL(ans)[i] = dexp_nimble(c_x[i], *c_rate, give_log);
  } else {
    int i_rate = 0;
    for(int i = 0; i < n_x; i++) {
      REAL(ans)[i] = dexp_nimble(c_x[i], c_rate[i_rate++], give_log);
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
  
SEXP C_rexp_nimble(SEXP n, SEXP rate) {
  // this will call rexp_nimble for computation on scalars
  if(!isInteger(n) || !isReal(rate) )
    RBREAK("Error (C_rexp_nimble): invalid input type for one of the arguments.");
  int n_rate = LENGTH(rate);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = allocVector(REALSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rexp_nimble): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = allocVector(REALSXP, n_values));  
  double* c_rate = REAL(rate);
  if(n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_values; i++) 
      REAL(ans)[i] = rexp_nimble(*c_rate);
  } else {
    int i_rate = 0;
    for(int i = 0; i < n_values; i++) {
      REAL(ans)[i] = rexp_nimble(c_rate[i_rate++]);
      // implement recycling:
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}
  
SEXP C_pexp_nimble(SEXP q, SEXP rate, SEXP lower_tail, SEXP log_p) {
  if(!isReal(q) || !isReal(rate) || !isLogical(lower_tail) || !isLogical(log_p)) 
    RBREAK("Error (C_pexp_nimble): invalid input type for one of the arguments.");
  int n_q = LENGTH(q);
  int n_rate = LENGTH(rate);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_q == 0) {
    return q;
  }
    
  PROTECT(ans = allocVector(REALSXP, n_q));  
  double* c_q = REAL(q);
  double* c_rate = REAL(rate);

  // FIXME: abstract the recycling as a function
  if(n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_q; i++) 
      REAL(ans)[i] = pexp_nimble(c_q[i], *c_rate, c_lower_tail, c_log_p);
  } else {
    int i_rate = 0;
    for(int i = 0; i < n_q; i++) {
      REAL(ans)[i] = pexp_nimble(c_q[i], c_rate[i_rate++], c_lower_tail, c_log_p);
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
  
SEXP C_qexp_nimble(SEXP p, SEXP rate, SEXP lower_tail, SEXP log_p) {
  if(!isReal(p) || !isReal(rate) || !isLogical(lower_tail) || !isLogical(log_p)) 
    RBREAK("Error (C_qexp_nimble): invalid input type for one of the arguments.");
  int n_p = LENGTH(p);
  int n_rate = LENGTH(rate);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_p == 0) {
    return p;
  }
    
  PROTECT(ans = allocVector(REALSXP, n_p));  
  double* c_p = REAL(p);
  double* c_rate = REAL(rate);

  // FIXME: abstract the recycling as a function
  if(n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_p; i++) 
      REAL(ans)[i] = qexp_nimble(c_p[i], *c_rate, c_lower_tail, c_log_p);
  } else {
    int i_rate = 0;
    for(int i = 0; i < n_p; i++) {
      REAL(ans)[i] = qexp_nimble(c_p[i], c_rate[i_rate++], c_lower_tail, c_log_p);
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
