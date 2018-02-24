/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

// various additional distributions functions needed by NIMBLE
// Author: Christopher Paciorek 
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

double dwish_chol(double* x, double* chol, double df, int p, double scale_param, int give_log, int overwrite_inputs) {
  char uplo('U');
  char sideL('L');
  char sideR('R');
  char diag('N');
  char transN('N');
  int info(0);
  double alpha(1.0);
  double* xChol;

  int i, j;

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
  if(overwrite_inputs & (int) scale_param)  // if !scale_param we need x below
    xChol = x;
  else {
    xChol = new double[p*p];
    for(i = 0; i < p*p; i++) 
      xChol[i] = x[i];
  }
  F77_CALL(dpotrf)(&uplo, &p, xChol, &p, &info);
  for(i = 0; i < p*p; i += p + 1) 
    dens += (df - p - 1) * log(xChol[i]);
  
  // R %*% x = t(chol) %*% chol %*% x (could also do with chol(x) but no more efficient
  // solve(S, x) = crossproduct( chol(x) %*% inverse(chol) )

  // dtr{m,s}m is a BLAS level-3 function
  double tmp_dens = 0.0;
  if(scale_param) {
    // chol(x) %*% inverse(chol)
    F77_CALL(dtrsm)(&sideR, &uplo, &transN, &diag, &p, &p, &alpha, 
                    chol, &p, xChol, &p);
    // trace of crossproduct of result is sum of squares of elements
    for(j = 0; j < p; j++) 
      for(i = 0; i <= j; i++) 
        tmp_dens += xChol[j*p+i]*xChol[j*p+i];
  } else {
    double* xCopy;
    if(overwrite_inputs) 
      xCopy = x;
    else {
      xCopy = new double[p*p]; 
      for(i = 0; i < p*p; i++) 
        xCopy[i] = x[i];
    }
    // chol %*% x
    F77_CALL(dtrmm)(&sideL, &uplo, &transN, &diag, &p, &p, &alpha, 
           chol, &p, xCopy, &p);
    // trace crossproduct of t(chol) with result is sum of product of upper-triangular elements
    for(j = 0; j < p; j++) {
      for(i = 0; i <= j; i++) {
        tmp_dens += xCopy[j*p+i] * chol[j*p+i];
      }
    }
    if(!overwrite_inputs)  
      delete [] xCopy;
  }

  if(!(overwrite_inputs & (int) scale_param))
    delete [] xChol;

    // attempt to improve above calcs by doing efficient U^T U multiply followed by direct product multiply, however this would not make use of threading provided by BLAS and even with one thread seems to be no faster
    // U^T*U directly followed by direct product with x
    /* 
    double tmp_summand;
    int minij;
    for(j = 0; j < p; j++) 
      for(i = 0; i < p; i++) {
        tmp_summand = 0.0;
        minij = i <= j ? i : j;
        for(int k = 0; k < minij; k++)
          tmp_summand += chol[j*p+k]*chol[i*p+k]; // U^T U 
        // double if not on diagonal to account for direct product of lower triangle too
        if(i != j) tmp_summand *= 2;
        tmp_dens += xCopy[j*p+i] * tmp_summand;
      }
    */
    
  dens += -0.5 * tmp_dens;

  return give_log ? dens : exp(dens);
}


SEXP C_dwish_chol(SEXP x, SEXP chol, SEXP df, SEXP scale_param, SEXP return_log) 
// calculates Wishart density given Cholesky of scale or rate matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isMatrix(x) || !Rf_isMatrix(chol) || !Rf_isReal(x) || !Rf_isReal(chol))
    RBREAK("Error (C_dwish_chol): 'x' and 'chol' must be real matrices.\n")
  if(!Rf_isReal(df) || !Rf_isReal(scale_param) || !Rf_isLogical(return_log))
    RBREAK("Error (C_dwish_chol): invalid input type for one of the arguments.\n");
  int* dims = INTEGER(Rf_getAttrib(x, R_DimSymbol));
  if(dims[0] != dims[1])
    RBREAK("Error (C_dwish_chol): 'x' must be a square matrix.\n");
  int p = dims[0];
  dims = INTEGER(Rf_getAttrib(chol, R_DimSymbol));
  if(dims[0] != dims[1] || dims[0] != p)
   RBREAK("Error (C_dwish_chol): 'chol' must be a square matrix with the same dimensions as 'x'.\n");

  int give_log = (int) LOGICAL(return_log)[0];
  double scale = REAL(scale_param)[0];

  double* c_x = REAL(x);
  double* c_chol = REAL(chol);
  double c_df = REAL(df)[0];
 
  if(c_df < p)
    RBREAK("Error (C_dwish_chol): inconsistent degrees of freedom and dimension.\n");
  
  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, 1));  
  REAL(ans)[0] = dwish_chol(c_x, c_chol, c_df, p, scale, give_log, 0); 
  UNPROTECT(1);
  return ans;
}


void rwish_chol(double *Z, double* chol, double df, int p, double scale_param, int overwrite_inputs) {
  char uplo('U');
  char sideL('L');
  char diag('N');
  char transT('T');
  char transN('N');
  double alpha(1.0);
  double beta(0.0);

  double* cholCopy;
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
  
  // fill diags with sqrts of chi-squares and upper triangle (for scale_param) with std normals - crossproduct of result is standardized Wishart; based on rWishart in stats package
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
  // would be more efficient if make use of fact that right-most matrix is triangular, but no available BLAS routine and hand-coding would eliminate use of threading and might well not be faster
  if(overwrite_inputs)
    cholCopy = chol;
  else {
    cholCopy = new double[p*p];
    if(scale_param)
      for(i = 0; i < p*p; i++) 
        cholCopy[i] = chol[i];
  }
  if(scale_param) F77_CALL(dtrmm)(&sideL, &uplo, &transN, &diag, &p, &p, &alpha, Z, &p, cholCopy, &p);
  else F77_CALL(dtrsm)(&sideL, &uplo, &transN, &diag, &p, &p, &alpha, chol, &p, Z, &p);

  // cp result to Z or chol so can be used as matrix to multiply against and overwrite
  if(scale_param) {
    for(j = 0; j < p*p; j++) 
      Z[j] = cholCopy[j];
  } else {
    for(j = 0; j < p*p; j++) 
      cholCopy[j] = Z[j]; 
  }

  // do crossprod of result
  // for dtrmm call, again this would be more efficient if use fact that RHS upper triangular, but no available BLAS routine and hand-coding would eliminate use of threading and might well not be faster
  if(scale_param) F77_CALL(dtrmm)(&sideL, &uplo, &transT, &diag, &p, &p, &alpha, cholCopy, &p, Z, &p);
  else F77_CALL(dgemm)(&transN, &transT, &p, &p, &p, &alpha, cholCopy, &p, cholCopy, &p, &beta, Z, &p); 
  if(!overwrite_inputs)
    delete [] cholCopy;
}

SEXP C_rwish_chol(SEXP chol, SEXP df, SEXP scale_param) 
// generates single Wishart draw given Cholesky of scale or rate matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isMatrix(chol) || !Rf_isReal(chol))
    RBREAK("Error (Crwish_chol): 'chol' must be a real matrix.\n")
  if(!Rf_isReal(df) || !Rf_isReal(scale_param))
    RBREAK("Error (C_rwish_chol): invalid input type for one of the arguments.\n");
  int* dims = INTEGER(Rf_getAttrib(chol, R_DimSymbol));
  if(dims[0] != dims[1])
    RBREAK("Error (C_rwish_chol): 'chol' must be a square matrix.\n");
  int p = dims[0];

  int n_chol = LENGTH(chol);
  double scale = REAL(scale_param)[0];

  double* c_chol = REAL(chol);
  double c_df = REAL(df)[0];

  if(c_df < p)
    RBREAK("Error (C_rwish_chol): inconsistent degrees of freedom and dimension.\n");

  GetRNGstate(); 

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, n_chol));  
  rwish_chol(REAL(ans), c_chol, c_df, p, scale, 0);
  
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}


double dinvwish_chol(double* x, double* chol, double df, int p, double scale_param, int give_log, int overwrite_inputs) {
  char uplo('U');
  char sideL('L');
  char sideR('R');
  char diag('N');
  char transT('T');
  char transN('N');
  int info(0);
  double alpha(1.0);
  double* xChol;

  int i, j;

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

  // determinant of S or R
  if(scale_param) {
    for(i = 0; i < p*p; i += p + 1) 
      dens += df * log(chol[i]);
  } else {
    for(i = 0; i < p*p; i += p + 1) 
      dens -= df * log(chol[i]);
  }

  // determinant of x using Cholesky:
  if(overwrite_inputs) 
    xChol = x;
  else {
    xChol = new double[p*p];
    for(i = 0; i < p*p; i++) 
      xChol[i] = x[i];
  }
  F77_CALL(dpotrf)(&uplo, &p, xChol, &p, &info);
  for(i = 0; i < p*p; i += p + 1) 
    dens -= (df + p + 1) * log(xChol[i]);
  
  double tmp_dens = 0.0;

  // do upper-triangular solve for scale parameterization
  // or upper-triangular multiplies for rate parameterization

  // S %*% inverse(x) = crossprod(chol %*% inverse(chol(x)) 
  // inverse(R) %*% inverse(x) = crossproduct(solve(t(chol(x)), inverse(chol)))
  // dtr{m,s}m is a BLAS level-3 function
  if(scale_param) {
    double* cholCopy;
    if(overwrite_inputs)
      cholCopy = chol;
    else {
      cholCopy = new double[p*p];
      for(i = 0; i < p*p; i++) 
        cholCopy[i] = chol[i];
    }
    // chol %*% inverse(chol(x))
    F77_CALL(dtrsm)(&sideR, &uplo, &transN, &diag, &p, &p, &alpha, 
           xChol, &p, cholCopy, &p);
    // trace of crossproduct is sum of elements squared
    for(j = 0; j < p; j++) {
      for(i = 0; i <= j; i++) {
        tmp_dens += cholCopy[j*p+i] * cholCopy[j*p+i];
      }
    }
    if(!overwrite_inputs)
      delete [] cholCopy;
  } else {
    double* iden = new double[p*p];
    for(j = 0; j < p; j++)
      for(i = 0; i < p; i++)
        if(i == j) iden[j*p+i] = 1.0; else iden[j*p+i] = 0.0;
    // inverse of chol
    F77_CALL(dtrsm)(&sideL, &uplo, &transN, &diag, &p, &p, &alpha, 
                    chol, &p, iden, &p);
    // solve(t(chol(x)), result)
    F77_CALL(dtrsm)(&sideL, &uplo, &transT, &diag, &p, &p, &alpha, 
                    xChol, &p, iden, &p);
    // trace of crossproduct is sum of elements squared
    for(j = 0; j < p; j++) {
      for(i = 0; i < p; i++) {
        tmp_dens += iden[j*p+i] * iden[j*p+i];
      }
    }
    delete [] iden;
  }
  dens += -0.5 * tmp_dens;
  return give_log ? dens : exp(dens);
}


SEXP C_dinvwish_chol(SEXP x, SEXP chol, SEXP df, SEXP scale_param, SEXP return_log) 
// calculates Inverse Wishart density given Cholesky of scale or rate matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isReal(x) || !Rf_isReal(chol) || !Rf_isReal(df) || !Rf_isReal(scale_param) || !Rf_isLogical(return_log))
    RBREAK("Error (C_dinvwish_chol): invalid input type for one of the arguments.\n");
  int p = pow(LENGTH(chol), 0.5);
  int give_log = (int) LOGICAL(return_log)[0];
  double scale = REAL(scale_param)[0];

  double* c_x = REAL(x);
  double* c_chol = REAL(chol);
  double c_df = REAL(df)[0];
 
  if(c_df < p)
    RBREAK("Error (C_dinvwish_chol): inconsistent degrees of freedom and dimension.\n");
  
  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, 1));  
  REAL(ans)[0] = dinvwish_chol(c_x, c_chol, c_df, p, scale, give_log, 0); 
  UNPROTECT(1);
  return ans;
}


void rinvwish_chol(double *Z, double* chol, double df, int p, double scale_param, int overwrite_inputs) {
  char uploU('U');
  char uploL('L');
  char sideL('L');
  char diag('N');
  char transT('T');
  char transN('N');
  double alpha(1.0);
  double beta(0.0);

  // calculations here generate W ~ Wi(S) for scale_param = 1 and W ~ Wi(S^-1) for scale_param = 0,
  // and then invert the result, giving InvWi(S^-1) for scale_param = 1 and InvWi(S) for scale_param = 0,
  // because scale for Wishart is inverse of scale for inverse-Wishart
  // therefore generate under opposite parameterization
  scale_param = 1 - scale_param;  


  double* cholCopy;
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
  
  // fill diags with sqrts of chi-squares and upper triangle (for scale_param) with std normals - crossproduct of result is standardized Wishart; based on rWishart in stats package
  for(j = 0; j < p; j++) {
    Z[j*p + j] = sqrt(rchisq(df - (double) j)); 
    for(i = 0; i < j; i++) {
      uind = i + j * p, /* upper triangle index */
      lind = j + i * p; /* lower triangle index */
      Z[(scale_param ? uind : lind)] = norm_rand();
      Z[(scale_param ? lind : uind)] = 0;
    }
  }
 
  if(overwrite_inputs)
    cholCopy = chol;
  else {
    cholCopy = new double[p*p];
    for(i = 0; i < p*p; i++) 
      cholCopy[i] = chol[i];
  }

  // with S parameterization: t(chol)%*%t(Z)%*%Z%*%chol is Wishart-distributed, so take inverse based on pieces
  // with R parameterization: inverse(chol)%*%Z%*%t(Z)%*%inverse(t(chol)) is Wishart-distributed, again inverse of pieces
  if(scale_param) {
    // Z %*% chol
    F77_CALL(dtrmm)(&sideL, &uploU, &transN, &diag, &p, &p, &alpha, Z, &p, cholCopy, &p);
    double* iden = new double[p*p];
    for(j = 0; j < p; j++)
      for(i = 0; i < p; i++)
        if(i == j) iden[j*p+i] = 1.0; else iden[j*p+i] = 0.0;
    // inverse(Z %*% chol)
    F77_CALL(dtrsm)(&sideL, &uploU, &transN, &diag, &p, &p, &alpha, cholCopy, &p, iden, &p);
    // crossproduct of result
    F77_CALL(dgemm)(&transN, &transT, &p, &p, &p, &alpha, iden, &p, iden, &p, &beta, Z, &p);    
    delete [] iden;
  } else {    
    // solve(Z, chol)
    F77_CALL(dtrsm)(&sideL, &uploL, &transN, &diag, &p, &p, &alpha, Z, &p, cholCopy, &p);
    // crossproduct of result
    F77_CALL(dgemm)(&transT, &transN, &p, &p, &p, &alpha, cholCopy, &p, cholCopy, &p, &beta, Z, &p);
  }

  if(!overwrite_inputs)
    delete [] cholCopy;
}

SEXP C_rinvwish_chol(SEXP chol, SEXP df, SEXP scale_param) 
// generates single Inverse Wishart draw given Cholesky of scale or rate matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isReal(chol) || !Rf_isReal(df) || !Rf_isReal(scale_param))
    RBREAK("Error (C_rinvwish_chol): invalid input type for one of the arguments.\n");
  int n_chol = LENGTH(chol);
  int p = pow(n_chol, 0.5);
  double scale = REAL(scale_param)[0];

  double* c_chol = REAL(chol);
  double c_df = REAL(df)[0];

  if(c_df < p)
    RBREAK("Error (C_rinvwish_chol): inconsistent degrees of freedom and dimension.\n");

  GetRNGstate(); 

  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, n_chol));  
  rinvwish_chol(REAL(ans), c_chol, c_df, p, scale, 0); 
  
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
  if(!Rf_isReal(x) || !Rf_isReal(alpha) || !Rf_isLogical(return_log)) 
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

  PROTECT(ans = Rf_allocVector(REALSXP, 1));  
  REAL(ans)[0] = ddirch(c_x, c_alpha, K, give_log);

  UNPROTECT(1);
  return ans;
}
 

SEXP C_rdirch(SEXP alpha) {
  if(!Rf_isReal(alpha))
    RBREAK("Error (C_rdirch): invalid input type for the argument.\n");
  int K = LENGTH(alpha);

  SEXP ans;

  if(K == 0) {
    PROTECT(ans = Rf_allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }

  double* c_alpha = REAL(alpha);

  GetRNGstate(); 

  PROTECT(ans = Rf_allocVector(REALSXP, K));  
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
  double logSumProb;

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
    if (prob[i] < 0) ML_ERR_return_NAN;
    R_D_nonint_check(x[i]);
    if (x[i] < 0 || !R_FINITE(x[i])) return R_D__0;

    x[i] = R_D_forceint(x[i]);
    sumProb += prob[i];
    sumX += x[i];
  }
  logSumProb = log(sumProb);

  for(int i = 0; i < K; i++) {  
    if(!(x[i] == 0.0 && prob[i] == 0.0))
      dens += x[i]*(log(prob[i]) - logSumProb) - lgammafn(x[i] + 1);
  }

  if(sumX > size + 10*DBL_EPSILON || sumX < size - 10*DBL_EPSILON) {
    return R_D__0;
  }

  return give_log ? dens : exp(dens);
}




void rmulti(int *ans, double size, double* prob, int K) // Calling functions need to copy first arg back and forth to double if needed
// scalar function that can be called directly by NIMBLE with same name as in R
// just call Rmath's rmultinom, which passes result by pointer
// IMPORTANT: have ans and size as int when sent to rmultinom as Rmath rmultinom has these types
// Nimble does a copy in nimArr_rmulti
{
  /* rmultinom requires normalized probs (in R this is done in Rmultinom interface
     function via FixupProb before passing to rmultinom) */
  double sumProb = 0.0;
  int i;

  for(i = 0; i < K; i++) 
    sumProb += prob[i];
  if (sumProb <= 0.0) {
    for(i = 0; i < K; i++) 
      ans[i] = R_NaN;
    return;
  }
  for(i = 0; i < K; i++) 
    prob[i] /= sumProb;
  rmultinom((int) size, prob, K, ans);
}

SEXP C_dmulti(SEXP x, SEXP size, SEXP prob, SEXP return_log) 
{
  if(!Rf_isReal(x) || !Rf_isReal(size) || !Rf_isReal(prob) || !Rf_isLogical(return_log)) 
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

  sum = 0.0;
  for(i = 0; i < K; i++) 
    sum += c_x[i];
  if(c_size > sum + 10*DBL_EPSILON || c_size < sum - 10*DBL_EPSILON)
    RBREAK("Error (C_dmulti): sum of values is not equal to size.\n");

  PROTECT(ans = Rf_allocVector(REALSXP, 1));  
  REAL(ans)[0] = dmulti(c_x, c_size, c_prob, K, give_log);

  UNPROTECT(1);
  return ans;
}
 

SEXP C_rmulti(SEXP size, SEXP prob) {
  if(!Rf_isReal(size) || !Rf_isReal(prob))
    RBREAK("Error (C_rmulti): invalid input type for one of the arguments.\n");
  int K = LENGTH(prob);

  SEXP ans;

  if(K == 0) {
    PROTECT(ans = Rf_allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }

  double* c_prob = REAL(prob);
  double c_size = REAL(size)[0];

  GetRNGstate(); 

  PROTECT(ans = Rf_allocVector(INTSXP, K));  
  // note that if NaN set in rmulti, the INTEGER() casts it to NA
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

  R_D_nonint_check(x);
  x = R_D_forceint(x);

  if(x > K || x < 1) return R_D__0;
  return give_log ? log(prob[(int) x - 1]) - log(sumProb) : prob[(int) x - 1]/sumProb;
}

double rcat(double* prob, int K)
// scalar function that can be called directly by NIMBLE with same name as in R
// problem is no apparent way to call this w/o also passing the number of categories
// relying on sum to 1 risks accessing beyond storage
// we'll need to figure out how to inject the number of categories w/in NIMBLE
{
  double prob_cum = prob[0];
  int value = 1;

  double sumProb(0.0);
  for(int i = 0; i < K; i++) 
    sumProb += prob[i];

  double u = unif_rand() * sumProb;
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

  if(!Rf_isReal(x) || !Rf_isReal(prob) || !Rf_isLogical(return_log)) 
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

  PROTECT(ans = Rf_allocVector(REALSXP, n_x));  
  for(i = 0; i < n_x; i++) {
    REAL(ans)[i] = dcat(c_x[i], c_prob, K, give_log);
  }

  UNPROTECT(1);
  return ans;
}
  
SEXP C_rcat(SEXP n, SEXP prob) {
  if(!Rf_isInteger(n) || !Rf_isReal(prob))
    RBREAK("Error (C_rcat): invalid input type for one of the arguments.\n");
  int n_values = INTEGER(n)[0];
  int K = LENGTH(prob);

  SEXP ans;
  int i;

  if(n_values == 0) {
    PROTECT(ans = Rf_allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    RBREAK("Error (C_rcat): n must be non-negative.\n");

  double* c_prob = REAL(prob);

  GetRNGstate(); 
  PROTECT(ans = Rf_allocVector(INTSXP, n_values));  

  for(i = 0; i < n_values; i++) 
    INTEGER(ans)[i] = (int) rcat(c_prob, K);

  PutRNGstate();
  UNPROTECT(1);
  return ans;

}
  


double dmnorm_chol(double* x, double* mean, double* chol, int n, double prec_param, int give_log, int overwrite_inputs) {
  char uplo('U');
  char transPrec('N');
  char transCov('T');
  char diag('N');
  int lda(n);
  int incx(1);
  double* xCopy;

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

  if(overwrite_inputs) {
    xCopy = x;
    for(i = 0; i < n; i++) 
      xCopy[i] -= mean[i];
  } else {
    xCopy = new double[n];
    for(i = 0; i < n; i++)
      xCopy[i] = x[i] - mean[i];
  }

  // do matrix-vector multiply with upper-triangular matrix stored column-wise as full n x n matrix (prec parameterization)
  // or upper-triangular (transpose) solve (cov parameterization)
  // dtr{m,s}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrmv)(&uplo, &transPrec, &diag, &n, chol, &lda, xCopy, &incx);
  else F77_CALL(dtrsv)(&uplo, &transCov, &diag, &n, chol, &lda, xCopy, &incx);

  // sum of squares to calculate quadratic form
  double tmp = 0.0;
  for(i = 0; i < n; i++)
    tmp += xCopy[i] * xCopy[i];

  dens += -0.5 * tmp;

  if(!overwrite_inputs)
    delete [] xCopy;

  return give_log ? dens : exp(dens);
}

SEXP C_dmnorm_chol(SEXP x, SEXP mean, SEXP chol, SEXP prec_param, SEXP return_log) 
// calculates mv normal density given Cholesky of precision matrix or covariance matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isMatrix(chol) || !Rf_isReal(chol))
    RBREAK("Error (C_dmnorm_chol): 'chol' must be a real matrix.\n");
  if(!Rf_isReal(x) || !Rf_isReal(mean))
    RBREAK("Error (C_dmnorm_chol): 'x' and 'mean' should be real valued.\n");
  if(!Rf_isReal(prec_param) || !Rf_isLogical(return_log))
    RBREAK("Error (C_dmnorm_chol): invalid input type for one of the arguments.\n");

  int* dims = INTEGER(Rf_getAttrib(chol, R_DimSymbol));
  if(dims[0] != dims[1])
    RBREAK("Error (C_dmnorm_chol): 'chol' must be a square matrix.\n");
  int p = dims[0];

  int n_x = LENGTH(x);
  if(n_x != p)
    RBREAK("Error (C_dmnorm_chol): 'x' and 'chol' are not of compatible sizes.\n");
  int n_mean = LENGTH(mean);
  int give_log = (int) LOGICAL(return_log)[0];
  double prec = REAL(prec_param)[0];

  double* c_x = REAL(x);
  double* c_mean = REAL(mean);
  double* c_chol = REAL(chol);

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
  PROTECT(ans = Rf_allocVector(REALSXP, 1));  
  REAL(ans)[0] = dmnorm_chol(c_x, full_mean, c_chol, n_x, prec, give_log, 0);
  if(n_mean < n_x)
    delete [] full_mean;
  UNPROTECT(1);
  return ans;
}

void rmnorm_chol(double *ans, double* mean, double* chol, int n, double prec_param) {
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

  for(i = 0; i < n; i++) 
    ans[i] = norm_rand();

  // do upper-triangular solve or (transpose) multiply
  // dtr{s,m}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrsv)(&uplo, &transPrec, &diag, &n, chol, &lda, ans, &incx);
  else F77_CALL(dtrmv)(&uplo, &transCov, &diag, &n, chol, &lda, ans, &incx);

  for(i = 0; i < n; i++) 
    ans[i] += mean[i];
}

SEXP C_rmnorm_chol(SEXP mean, SEXP chol, SEXP prec_param) 
// generates single mv normal draw given Cholesky of precision matrix or covariance matrix
// Cholesky matrix should be given as a numeric vector in column-major order
//   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isMatrix(chol) || !Rf_isReal(chol))
    RBREAK("Error (C_rmnorm_chol): 'chol' should be a real matrix.\n");
  if(!Rf_isReal(mean))
    RBREAK("Error (C_rmnorm_chol): 'mean' should be real-valued\n");
  if(!Rf_isReal(prec_param))
    RBREAK("Error (C_rmnorm_chol): invalid input type for one of the arguments.\n");

  int* dims = INTEGER(Rf_getAttrib(chol, R_DimSymbol));
  if(dims[0] != dims[1])
    RBREAK("Error (C_dmnorm_chol): 'chol' must be a square matrix.\n");
  int n_values = dims[0];

  int n_mean = LENGTH(mean);
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
  PROTECT(ans = Rf_allocVector(REALSXP, n_values));  
  rmnorm_chol(REAL(ans), full_mean, c_chol, n_values, prec);

  PutRNGstate();
  if(n_mean < n_values) 
    delete [] full_mean;
  UNPROTECT(1);
  return ans;
}

// Begin multivariate t

double dmvt_chol(double* x, double* mu, double* chol, double df, int n, double prec_param, int give_log, int overwrite_inputs) {
  char uplo('U');
  char transPrec('N');
  char transCov('T');
  char diag('N');
  int lda(n);
  int incx(1);
  double* xCopy;
  
  double dens = lgammafn((df + n) / 2) - lgammafn(df / 2) - n * M_LN_SQRT_PI - n * log(df) / 2;
  int i;
  
  if (R_IsNA(x, n) || R_IsNA(mu, n) || R_IsNA(chol, n*n) || R_IsNA(df) || R_IsNA(prec_param))
    return NA_REAL;
#ifdef IEEE_754
  if (R_isnancpp(x, n) || R_isnancpp(mu, n) || R_isnancpp(chol, n*n) || R_IsNA(df) || R_isnancpp(prec_param))
    return R_NaN;
#endif
  
  if(!R_FINITE_VEC(x, n) || !R_FINITE_VEC(mu, n) || !R_FINITE_VEC(chol, n*n)) return R_D__0;
  
  // add diagonals of Cholesky
  if(prec_param) {
    for(i = 0; i < n*n; i += n + 1) 
      dens += log(chol[i]);
  } else {
    for(i = 0; i < n*n; i += n + 1) 
      dens -= log(chol[i]);
  }

  if(overwrite_inputs) {
    xCopy = x;
    for(i = 0; i < n; i++) 
      xCopy[i] -= mu[i];
  } else {
    xCopy = new double[n];
    for(i = 0; i < n; i++)
      xCopy[i] = x[i] - mu[i];
  }
  
  // do matrix-vector multiply with upper-triangular matrix stored column-wise as full n x n matrix (prec parameterization)
  // or upper-triangular (transpose) solve (cov parameterization)
  // dtr{m,s}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrmv)(&uplo, &transPrec, &diag, &n, chol, &lda, xCopy, &incx);
  else F77_CALL(dtrsv)(&uplo, &transCov, &diag, &n, chol, &lda, xCopy, &incx);
  
  // sum of squares to calculate quadratic form
  double tmp = 0.0;
  for(i = 0; i < n; i++)
    tmp += xCopy[i] * xCopy[i];
  
  dens += -0.5 * (df + n) * log(1 + tmp / df);
  
  if(!overwrite_inputs)
    delete [] xCopy;

  return give_log ? dens : exp(dens);
}

SEXP C_dmvt_chol(SEXP x, SEXP mu, SEXP chol, SEXP df, SEXP prec_param, SEXP return_log) 
  // calculates mv normal density given Cholesky of precision matrix or covariance matrix
  // Cholesky matrix should be given as a numeric vector in column-major order
  //   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isMatrix(chol) || !Rf_isReal(chol))
    RBREAK("Error (C_dmvt_chol): 'chol' must be a real matrix.\n");
  if(!Rf_isReal(x) || !Rf_isReal(mu))
    RBREAK("Error (C_dmvt_chol): 'x' and 'mu' should be real valued.\n");
  if(!Rf_isReal(df) || !Rf_isReal(prec_param) || !Rf_isLogical(return_log))
    RBREAK("Error (C_dmvt_chol): invalid input type for one of the arguments.\n");

  int* dims = INTEGER(Rf_getAttrib(chol, R_DimSymbol));
  if(dims[0] != dims[1])
    RBREAK("Error (C_dmvt_chol): 'chol' must be a square matrix.\n");
  int p = dims[0];

  int n_x = LENGTH(x);
  if(n_x != p)
    RBREAK("Error (C_dmvt_chol): 'x' and 'chol' are not of compatible sizes.\n");
  int n_mu = LENGTH(mu);
  int give_log = (int) LOGICAL(return_log)[0];
  double c_df = REAL(df)[0];
  double prec = REAL(prec_param)[0];
  
  double* c_x = REAL(x);
  double* c_mu = REAL(mu);
  double* c_chol = REAL(chol);
  
  double* full_mu;
  if(n_mu < n_x) {
    full_mu = new double[n_x];
    int i_mu = 0;
    for(int i = 0; i < n_x; i++) {
      full_mu[i] = c_mu[i_mu++];
      if(i_mu == n_mu) i_mu = 0;
    }
  } else full_mu = c_mu;
  
  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, 1));  
  REAL(ans)[0] = dmvt_chol(c_x, full_mu, c_chol, c_df, n_x, prec, give_log, 0); 
  if(n_mu < n_x)
    delete [] full_mu;
  UNPROTECT(1);
  return ans;
}

void rmvt_chol(double *ans, double* mu, double* chol, double df, int n, double prec_param) {
  char uplo('U');
  char transPrec('N');
  char transCov('T');
  char diag('N');
  int lda(n);
  int incx(1);
  
  int i, j;
  
#ifdef IEEE_754
  if (R_isnancpp(mu, n) || R_isnancpp(chol, n*n) || R_isnancpp(df) || R_isnancpp(prec_param)) {
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
    
  for(i = 0; i < n; i++) 
    ans[i] = norm_rand();
  
  // sample from chi-squared and calculate scaling factor
  double scaling = sqrt(df / rchisq(df));
  
  // do upper-triangular solve or (transpose) multiply
  // dtr{s,m}v is a BLAS level-2 function
  if(prec_param) F77_CALL(dtrsv)(&uplo, &transPrec, &diag, &n, chol, &lda, ans, &incx);
  else F77_CALL(dtrmv)(&uplo, &transCov, &diag, &n, chol, &lda, ans, &incx);
  
  for(i = 0; i < n; i++) 
    ans[i] = mu[i] + ans[i] * scaling;
}

SEXP C_rmvt_chol(SEXP mu, SEXP chol, SEXP df, SEXP prec_param) 
  // generates single mv normal draw given Cholesky of precision matrix or covariance matrix
  // Cholesky matrix should be given as a numeric vector in column-major order
  //   including all n x n elements; lower-triangular elements are ignored
{
  if(!Rf_isMatrix(chol) || !Rf_isReal(chol))
    RBREAK("Error (C_rmvt_chol): 'chol' should be a real matrix.\n");
  if(!Rf_isReal(mu))
    RBREAK("Error (C_rmvt_chol): 'mu' should be real-valued\n");
  if(!Rf_isReal(prec_param))
    RBREAK("Error (C_rmvt_chol): invalid input type for one of the arguments.\n");

  int* dims = INTEGER(Rf_getAttrib(chol, R_DimSymbol));
  if(dims[0] != dims[1])
    RBREAK("Error (C_rmvt_chol): 'chol' must be a square matrix.\n");
  int n_values = dims[0];

  int n_mu = LENGTH(mu);
  double c_df = REAL(df)[0];
  double prec = REAL(prec_param)[0];
  
  int i;
  
  double* c_mu = REAL(mu);
  double* c_chol = REAL(chol);
  double* full_mu; 
  
  if(n_mu < n_values) {
    full_mu = new double[n_values];
    int i_mu = 0;
    for(i = 0; i < n_values; i++) {
      full_mu[i] = c_mu[i_mu++];
      if(i_mu == n_mu) i_mu = 0;
    }
  } else full_mu = c_mu;
  
  GetRNGstate(); 
  
  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, n_values));  
  rmvt_chol(REAL(ans), full_mu, c_chol, c_df, n_values, prec);
  
  PutRNGstate();
  if(n_mu < n_values) 
    delete [] full_mu;
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
  if(!Rf_isReal(x) || !Rf_isReal(df) || !Rf_isReal(mu) || !Rf_isReal(sigma) || !Rf_isLogical(return_log)) 
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
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_x));  
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
  if(!Rf_isInteger(n) || !Rf_isReal(df) || !Rf_isReal(mu) || !Rf_isReal(sigma))
    RBREAK("Error (C_rt_nonstandard): invalid input type for one of the arguments.");
  int n_mu = LENGTH(mu);
  int n_sigma = LENGTH(sigma);
  int n_df = LENGTH(df);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = Rf_allocVector(REALSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rt_nonstandard): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_values));  
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
  if(!Rf_isReal(q) || !Rf_isReal(df) || !Rf_isReal(mu) || !Rf_isReal(sigma) || !Rf_isLogical(lower_tail) || !Rf_isLogical(log_p))
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
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_q));  
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
  if(!Rf_isReal(p) || !Rf_isReal(df) || !Rf_isReal(mu) || !Rf_isReal(sigma) || !Rf_isLogical(lower_tail) || !Rf_isLogical(log_p))
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
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_p));  
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
  if(!Rf_isReal(x) || !Rf_isReal(t) || !Rf_isReal(c) || !Rf_isLogical(return_log)) 
    RBREAK("Error (C_dinterval): invalid input type for one of the arguments.");
  int n_x = LENGTH(x);
  int n_t = LENGTH(t);
  int n_c = LENGTH(c);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(n_x == 0) {
    return x;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_x));  
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
  if(!Rf_isInteger(n) || !Rf_isReal(t) || !Rf_isReal(c))
    RBREAK("Error (C_rinterval): invalid input type for one of the arguments.");
  int n_t = LENGTH(t);
  int K = LENGTH(c);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = Rf_allocVector(INTSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rinterval): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = Rf_allocVector(INTSXP, n_values));  
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
  if(!Rf_isReal(x) || !Rf_isReal(rate) || !Rf_isLogical(return_log))
    RBREAK("Error (C_dexp_nimble): invalid input type for one of the arguments.");
  int n_x = LENGTH(x);
  int n_rate = LENGTH(rate);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(n_x == 0) {
    return x;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_x));  
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
  if(!Rf_isInteger(n) || !Rf_isReal(rate) )
    RBREAK("Error (C_rexp_nimble): invalid input type for one of the arguments.");
  int n_rate = LENGTH(rate);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = Rf_allocVector(REALSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rexp_nimble): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_values));  
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
  if(!Rf_isReal(q) || !Rf_isReal(rate) || !Rf_isLogical(lower_tail) || !Rf_isLogical(log_p)) 
    RBREAK("Error (C_pexp_nimble): invalid input type for one of the arguments.");
  int n_q = LENGTH(q);
  int n_rate = LENGTH(rate);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_q == 0) {
    return q;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_q));  
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
  if(!Rf_isReal(p) || !Rf_isReal(rate) || !Rf_isLogical(lower_tail) || !Rf_isLogical(log_p)) 
    RBREAK("Error (C_qexp_nimble): invalid input type for one of the arguments.");
  int n_p = LENGTH(p);
  int n_rate = LENGTH(rate);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_p == 0) {
    return p;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_p));  
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


// used solely in conjugacy for dhalfflat (so dsqrtinvgamma really just a placeholder)
double dsqrtinvgamma(double x, double shape, double rate, int give_log)
{
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(shape) || ISNAN(rate))
    return x + shape + rate;
#endif
  double out = dinvgamma(x*x, shape, rate, 1) + log(2*x);
  if(give_log) return(out);
  else return(exp(out));
}

double rsqrtinvgamma(double shape, double rate)
{
#ifdef IEEE_754
  if (ISNAN(shape) || ISNAN(rate))
    ML_ERR_return_NAN;
#endif
  return(pow(rinvgamma(shape, rate), 0.5));
}

SEXP C_dsqrtinvgamma(SEXP x, SEXP shape, SEXP rate, SEXP return_log) {
  if(!Rf_isReal(x) || !Rf_isReal(shape) || !Rf_isReal(rate) || !Rf_isLogical(return_log)) 
    RBREAK("Error (C_dsqrtinvgamma): invalid input type for one of the arguments.");
  int n_x = LENGTH(x);
  int n_shape = LENGTH(shape);
  int n_rate = LENGTH(rate);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(n_x == 0) {
    return x;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_x));  
  double* c_x = REAL(x);
  double* c_shape = REAL(shape);
  double* c_rate = REAL(rate);

  // FIXME: abstract the recycling as a function
  if(n_rate == 1 && n_shape == 1 && n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_x; i++) 
      REAL(ans)[i] = dsqrtinvgamma(c_x[i], *c_shape, *c_rate, give_log);
  } else {
    int i_shape = 0;
    int i_rate = 0;
    for(int i = 0; i < n_x; i++) {
      REAL(ans)[i] = dsqrtinvgamma(c_x[i], c_shape[i_shape++], c_rate[i_rate++], give_log);
      // implement recycling:
      if(i_shape == n_shape) i_shape = 0;
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
  
SEXP C_rsqrtinvgamma(SEXP n, SEXP shape, SEXP rate) {
  if(!Rf_isInteger(n) || !Rf_isReal(shape) || !Rf_isReal(rate))
    RBREAK("Error (C_rsqrtinvgamma): invalid input type for one of the arguments.");
  int n_shape = LENGTH(shape);
  int n_rate = LENGTH(rate);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = Rf_allocVector(REALSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rsqrtinvgamma): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_values));  
  double* c_shape = REAL(shape);
  double* c_rate = REAL(rate);
  if(n_rate == 1 && n_shape == 1 && n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_values; i++) 
      REAL(ans)[i] = rsqrtinvgamma(*c_shape, *c_rate);
  } else {
    int i_shape = 0;
    int i_rate = 0;
    for(int i = 0; i < n_values; i++) {
      REAL(ans)[i] = rsqrtinvgamma(c_shape[i_shape++], c_rate[i_rate++]);
      // implement recycling:
      if(i_shape == n_shape) i_shape = 0;
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}

double dflat(double x, int give_log)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(x))
    return x;
#endif
  if(give_log) return 0;
  else return 1;
}

double rflat()
// scalar function that can be called directly by NIMBLE with same name as in R
{
  return R_NaN ;
}

double dhalfflat(double x, int give_log)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(x))
    return x;
#endif
  if(x >= 0) x = 0; else x = ML_NEGINF;
  if(!give_log) x = exp(x);
  return x;
}

double rhalfflat()
// scalar function that can be called directly by NIMBLE with same name as in R
{
  return R_NaN;
}

// NIMBLE's C parameterization of invgamma is (shape,rate) because
// when passing to R's C gamma, which uses (shape,scale), 
// the gamma scale parameter is the invgamma rate parameter
double dinvgamma(double x, double shape, double rate, int give_log)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(shape) || ISNAN(rate))
    return x + shape + rate;
#endif
  double xinv = 1/x;
  if(give_log) return(dgamma(xinv, shape, rate, give_log) - 2*log(x));
  else return(dgamma(xinv, shape, rate, give_log) * xinv * xinv);
}

double rinvgamma(double shape, double rate)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(shape) || ISNAN(rate))
    ML_ERR_return_NAN;
#endif
  return(1 / rgamma(shape, rate));
}

double pinvgamma(double q, double shape, double rate, int lower_tail, int log_p)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if(ISNAN(q) || ISNAN(shape) || ISNAN(rate))
    return q + shape + rate;
#endif
  return(pgamma(1/q, shape, rate, !lower_tail, log_p));
}

double qinvgamma(double p, double shape, double rate, int lower_tail, int log_p)
// scalar function that can be called directly by NIMBLE with same name as in R
{
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(shape) || ISNAN(rate))
    return p + shape + rate;
#endif
  return(1 / qgamma(p, shape, rate, !lower_tail, log_p));
}


SEXP C_dinvgamma(SEXP x, SEXP shape, SEXP rate, SEXP return_log) {
  if(!Rf_isReal(x) || !Rf_isReal(shape) || !Rf_isReal(rate) || !Rf_isLogical(return_log)) 
    RBREAK("Error (C_dinvgamma): invalid input type for one of the arguments.");
  int n_x = LENGTH(x);
  int n_shape = LENGTH(shape);
  int n_rate = LENGTH(rate);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
    
  if(n_x == 0) {
    return x;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_x));  
  double* c_x = REAL(x);
  double* c_shape = REAL(shape);
  double* c_rate = REAL(rate);

  // FIXME: abstract the recycling as a function
  if(n_rate == 1 && n_shape == 1 && n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_x; i++) 
      REAL(ans)[i] = dinvgamma(c_x[i], *c_shape, *c_rate, give_log);
  } else {
    int i_shape = 0;
    int i_rate = 0;
    for(int i = 0; i < n_x; i++) {
      REAL(ans)[i] = dinvgamma(c_x[i], c_shape[i_shape++], c_rate[i_rate++], give_log);
      // implement recycling:
      if(i_shape == n_shape) i_shape = 0;
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
  
SEXP C_rinvgamma(SEXP n, SEXP shape, SEXP rate) {
  if(!Rf_isInteger(n) || !Rf_isReal(shape) || !Rf_isReal(rate))
    RBREAK("Error (C_rinvgamma): invalid input type for one of the arguments.");
  int n_shape = LENGTH(shape);
  int n_rate = LENGTH(rate);
  int n_values = INTEGER(n)[0];
  SEXP ans;
    
  if(n_values == 0) {
    PROTECT(ans = Rf_allocVector(REALSXP, 0));
    UNPROTECT(1);
    return ans;
  }
  if(n_values < 0)
    // should formalize using R's C error-handling API
    RBREAK("Error (C_rinvgamma): n must be non-negative.\n");
    
  GetRNGstate(); 
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_values));  
  double* c_shape = REAL(shape);
  double* c_rate = REAL(rate);
  if(n_rate == 1 && n_shape == 1 && n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_values; i++) 
      REAL(ans)[i] = rinvgamma(*c_shape, *c_rate);
  } else {
    int i_shape = 0;
    int i_rate = 0;
    for(int i = 0; i < n_values; i++) {
      REAL(ans)[i] = rinvgamma(c_shape[i_shape++], c_rate[i_rate++]);
      // implement recycling:
      if(i_shape == n_shape) i_shape = 0;
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}
  
SEXP C_pinvgamma(SEXP q, SEXP shape, SEXP rate, SEXP lower_tail, SEXP log_p) {
  if(!Rf_isReal(q) || !Rf_isReal(shape) || !Rf_isReal(rate) || !Rf_isLogical(lower_tail) || !Rf_isLogical(log_p))
    RBREAK("Error (C_pinvgamma): invalid input type for one of the arguments.");
  int n_q = LENGTH(q);
  int n_shape = LENGTH(shape);
  int n_rate = LENGTH(rate);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_q == 0) {
    return q;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_q));  
  double* c_q = REAL(q);
  double* c_shape = REAL(shape);
  double* c_rate = REAL(rate);

  // FIXME: abstract the recycling as a function
  if(n_rate == 1 && n_shape == 1 && n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_q; i++) 
      REAL(ans)[i] = pinvgamma(c_q[i], *c_shape, *c_rate, c_lower_tail, c_log_p);
  } else {
    int i_shape = 0;
    int i_rate = 0;
    for(int i = 0; i < n_q; i++) {
      REAL(ans)[i] = pinvgamma(c_q[i], c_shape[i_shape++], c_rate[i_rate++], c_lower_tail, c_log_p);
      // implement recycling:
      if(i_shape == n_shape) i_shape = 0;
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}
 
SEXP C_qinvgamma(SEXP p, SEXP shape, SEXP rate, SEXP lower_tail, SEXP log_p) {
  if(!Rf_isReal(p) || !Rf_isReal(shape) || !Rf_isReal(rate) || !Rf_isLogical(lower_tail) || !Rf_isLogical(log_p))
    RBREAK("Error (C_qinvgamma): invalid input type for one of the arguments.");
  int n_p = LENGTH(p);
  int n_shape = LENGTH(shape);
  int n_rate = LENGTH(rate);
  int c_lower_tail = (int) LOGICAL(lower_tail)[0];
  int c_log_p = (int) LOGICAL(log_p)[0];
  SEXP ans;
    
  if(n_p == 0) {
    return p;
  }
    
  PROTECT(ans = Rf_allocVector(REALSXP, n_p));  
  double* c_p = REAL(p);
  double* c_shape = REAL(shape);
  double* c_rate = REAL(rate);

  // FIXME: abstract the recycling as a function
  if(n_rate == 1 && n_shape == 1 && n_rate == 1) {
    // if no parameter vectors, more efficient not to deal with multiple indices
    for(int i = 0; i < n_p; i++) 
      REAL(ans)[i] = qinvgamma(c_p[i], *c_shape, *c_rate, c_lower_tail, c_log_p);
  } else {
    int i_shape = 0;
    int i_rate = 0;
    for(int i = 0; i < n_p; i++) {
      REAL(ans)[i] = qinvgamma(c_p[i], c_shape[i_shape++], c_rate[i_rate++], c_lower_tail, c_log_p);
      // implement recycling:
      if(i_shape == n_shape) i_shape = 0;
      if(i_rate == n_rate) i_rate = 0;
    }
  }
    
  UNPROTECT(1);
  return ans;
}

SEXP C_dcar_normal(SEXP x, SEXP adj, SEXP weights, SEXP num, SEXP tau, SEXP c, SEXP zero_mean, SEXP return_log) {
  if(!Rf_isReal(x) || !Rf_isReal(adj) || !Rf_isReal(weights) || !Rf_isReal(num) || !Rf_isReal(tau) || !Rf_isReal(c) || !Rf_isReal(zero_mean) || !Rf_isLogical(return_log))
    RBREAK("Error (C_dcar_normal): invalid input type for one of the arguments.");
  
  int N = LENGTH(x);
  int L = LENGTH(adj);
  
  double* c_x = REAL(x);
  double* c_adj = REAL(adj);
  double* c_weights = REAL(weights);
  double* c_num = REAL(num);
  double c_tau = REAL(tau)[0];
  int c_c = (int) REAL(c)[0];
  int c_zero_mean = (int) REAL(zero_mean)[0];
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
  
  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  REAL(ans)[0] = dcar_normal(c_x, c_adj, c_weights, c_num, c_tau, c_c, c_zero_mean, N, L, give_log);
  
  UNPROTECT(1);
  return ans;
}


double dcar_normal(double* x, double* adj, double* weights, double* num, double tau, int c, int zero_mean, int N, int L, int give_log) {
  // This method implements the following density calculation:
  // p(x1, ..., xn, tau) = (tau/2/pi)^((N-c)/2) * exp(-tau/2 * sum_{i != j) w_ij (xi-xj)^2 ),
  // where tau is precision, c = (number of islands), and N is the length of x.
  // This is the density on an N-c dimensional space, and improper on the remaining c dimensions.
  if(tau < 0) {
    return R_NaN;
  }
  //PRINTF("c is equal to %d\n", c);
  double lp = 0;
  int count = 0;
  double xi, xj;
  for(int i = 0; i < N; i++) {
    xi = x[i];
    for(int j = 0; j < num[i]; j++) {
      xj = x[ (int) adj[count] - 1 ];
      lp += weights[count] * pow(xi-xj,2);
      count++;
    }
  }
  if(count != L) {
    ML_ERR_return_NAN;
  }
  lp *= 1/2.0;     // accounts for double-summing over all (xi,xj) pairs
  lp *= (-1/2.0) * tau;
  lp += (N-c)/2.0 * (log(tau) - M_LN_2PI);
  if(give_log) {
    return(lp);
  }
  return(exp(lp));
}

SEXP C_dcar_proper(SEXP x, SEXP mu, SEXP C, SEXP adj, SEXP num, SEXP M, SEXP tau, SEXP gamma, SEXP evs, SEXP return_log) {
  if(!Rf_isReal(x) || !Rf_isReal(mu) || !Rf_isReal(C) || !Rf_isReal(adj) || !Rf_isReal(num) || !Rf_isReal(M) || !Rf_isReal(tau) || !Rf_isReal(gamma) || !Rf_isReal(evs) || !Rf_isLogical(return_log))
    RBREAK("Error (C_dcar_proper): invalid input type for one of the arguments.");
  
  int N = LENGTH(x);
  int L = LENGTH(adj);
  
  double* c_x = REAL(x);
  double* c_mu = REAL(mu);
  double* c_C = REAL(C);
  double* c_adj = REAL(adj);
  double* c_num = REAL(num);
  double* c_M = REAL(M);
  double c_tau = REAL(tau)[0];
  double c_gamma = REAL(gamma)[0];
  double* c_evs = REAL(evs);
  int give_log = (int) LOGICAL(return_log)[0];
  SEXP ans;
  
  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  REAL(ans)[0] = dcar_proper(c_x, c_mu, c_C, c_adj, c_num, c_M, c_tau, c_gamma, c_evs, N, L, give_log);
  
  UNPROTECT(1);
  return ans;
}

double dcar_proper(double* x, double* mu, double* C, double* adj, double* num, double* M, double tau, double gamma, double* evs, int N, int L, int give_log) {
  // This method implements the following density calculation:
  // x ~ MVN( mean = mu,   cov  = (I-gamma*Cmatrix)^-1* %*% Mmatrix / tau )
  // note that the scalar 'tau' is a precision term.
  if(tau < 0) {
    return R_NaN;
  }
  double lp = 0;
  double xi, xj;
  int i, j, index;
  int count = 0;
  for(i = 0; i < N; i++) {
    xi = x[i] - mu[i];
    lp += pow(xi, 2) / M[i];      // (x-mu)' M^-1 (x-mu)
    //PRINTF("lp is equal to %f\n", lp);
    for(j = 0; j < num[i]; j++) {
      index = (int) adj[count] - 1;
      xj = x[index] - mu[index];
      lp -= gamma * xi * xj * C[count] / M[i];  // -gamma (x-mu)' M^-1 C (x-mu)
      count++;
    }
  }
  if(count != L) {
    ML_ERR_return_NAN;
  }
  lp *= (-1/2.0) * tau;
  // now add -1/2*log(|2*pi*Sigma|) to lp:
  // det(Sigma) = det((I-gamma*Cmatrix)^-1 %*% M / tau) = prod(M) / prod(1 - gamma*evs) / tau^N
  // where evs = eigen(Cmatrix)$values, and tau is precision
  for(i = 0; i < N; i++) {
    lp += (log(1 - gamma*evs[i]) - log(M[i])) / 2.0;
  }
  lp += N / 2.0 * (log(tau) - M_LN_2PI);
  if(give_log) {
    return(lp);
  }
  return(exp(lp));
}


SEXP C_rcar_proper(SEXP n, SEXP mu, SEXP C, SEXP adj, SEXP num, SEXP M, SEXP tau, SEXP gamma, SEXP evs) {
  if(!Rf_isInteger(n) || !Rf_isReal(mu) || !Rf_isReal(C) || !Rf_isReal(adj) || !Rf_isReal(num) || !Rf_isReal(M) || !Rf_isReal(tau) || !Rf_isReal(gamma) || !Rf_isReal(evs))
    RBREAK("Error (C_rcar_proper): invalid input type for one of the arguments.");
  
  int N = LENGTH(mu);
  int L = LENGTH(adj);
  
  double* c_mu = REAL(mu);
  double* c_C = REAL(C);
  double* c_adj = REAL(adj);
  double* c_num = REAL(num);
  double* c_M = REAL(M);
  double c_tau = REAL(tau)[0];
  double c_gamma = REAL(gamma)[0];
  double* c_evs = REAL(evs);
  
  SEXP ans;
  
  GetRNGstate();
  
  PROTECT(ans = Rf_allocVector(REALSXP, N));
  
  rcar_proper(REAL(ans), c_mu, c_C, c_adj, c_num, c_M, c_tau, c_gamma, c_evs, N, L);
  
  PutRNGstate();
  UNPROTECT(1);
  return ans;
}


void rcar_proper(double* ans, double* mu, double* C, double* adj, double* num, double* M, double tau, double gamma, double* evs, int N, int L) {
  
#ifdef IEEE_754
  if (R_isnancpp(mu, N) || R_isnancpp(C, L) || R_isnancpp(adj, L) || R_isnancpp(num, N) || R_isnancpp(M, N) || R_isnancpp(tau) || R_isnancpp(gamma) || R_isnancpp(evs, N)) {
    for(int i = 0; i < N; i++) {
      ans[i] = R_NaN;
    }
    return;
  }
#endif
  
  int i, j;
  
  double* Qchol = new double[N*N];    // precison Q = tau * M^-1 %*% (I - gamma*C)
  for(i = 0; i < N*N; i++) {
    Qchol[i] = 0;
  }
  int count = 0;
  for(i = 0; i < N; i++) {
    Qchol[ i*(N+1) ] = tau / M[i];    // main diagonal
    for(j = 0; j < num[i]; j++) {
      Qchol[  i  +  N*((int) adj[count] - 1)  ] = -tau * gamma * C[count] / M[i];
      count++;
    }
  }
  
  // now take chol(precision), and write that into Qchol
  // F77_CALL(dpotrf)( UPLO, N, A, LDA, INFO ) does Choleskey factorization of a symmetric pos. def. matrix
  // http://www.netlib.org/lapack/double/dpotrf.f
  char uplo('U');
  int info(0);
  F77_CALL(dpotrf)(&uplo, &N, Qchol, &N, &info);
  
  if(!R_FINITE_VEC(Qchol, N*N)) {
    for(i = 0; i < N; i++) {
      ans[i] = R_NaN;
    }
    delete [] Qchol;
    return;
  }
  
  for(i = 0; i < N; i++) {
    ans[i] = norm_rand();
  }
  // upper-triangular solve or (transpose) multiply
  char transPrec('N');
  char diag('N');
  int lda(N);
  int incx(1);
  F77_CALL(dtrsv)(&uplo, &transPrec, &diag, &N, Qchol, &lda, ans, &incx);
  for(i = 0; i < N; i++) {
    ans[i] += mu[i];
  }
  delete [] Qchol;
}

