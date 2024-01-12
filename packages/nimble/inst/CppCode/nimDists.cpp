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

#include <iostream> // must go before other things because R defines a "length" macro
#include "nimble/dists.h"
#include "nimble/nimDists.h"
#include "nimble/RcppUtils.h"

// These are used for compilation of DSL code and presumably belong elsewhere. Where?
bool R_IsNA_ANY(NimArr<1, double> &P) {
  int s = P.size();
  for(int i = 0; i < s; ++i) if(R_IsNA(P[i])) return(true);
  return(false);
}

bool R_isnancpp_ANY(NimArr<1, double> &P) {
  int s = P.size();
  for(int i = 0; i < s; ++i) if(R_isnancpp(P[i])) return(true);
  return(false);
}


// template<int nDim, class T>
// NimArr<nDim, T> &nimArrCopyIfNeeded(NimArr<nDim, T> &orig, NimArr<nDim, T> &possibleCopy) {
//   if(orig.isMap()) {
//     //    printf("It is a map\n");
//     if(!isMapEntire<nDim, T>(orig)) {
//       // NOTE that a map that has offset != 0 always returns false.
//       // This could be improved.
//       // E.g. mu[1:3, i] could always be handled without copying
//       // currently isMapEntire returns FALSE if i > 0 because then offset != 0
//       // We would need more care with checking when to copy back at end of functions
//       // below.
//       //      printf("It is not an entire map\n");
//       possibleCopy = orig;
//       return(possibleCopy);
//     }
//   }
//   return(orig);
// }

double nimArr_dmulti(NimArr<1, double> &x, double size, NimArr<1, double> &prob, int give_log) {
  int K = prob.size();
  if(K == 0) return(0.);
  if(x.size() != K) {
    _nimble_global_output<<"Error in nimArr_dmulti: incompatible sizes for x ("<<x.size()<<") and prob("<<K<<").\n";
    nimble_print_to_R(_nimble_global_output);
  }
  double *xptr; 
  double *probptr;
  NimArr<1, double> xCopy;
  NimArr<1, double> probCopy;
  xCopy = x; // copy from double to int
  xptr = xCopy.getPtr();
  probptr = nimArrCopyIfNeeded<1, double>(prob, probCopy).getPtr();
  double ans = dmulti(xptr, size, probptr, K, give_log);
  return(ans);
}

void nimArr_rmulti(NimArr<1, double> &ans, double size, NimArr<1, double> &prob) {
  int K = prob.size();
  if(K == 0) return;
  int *ansptr;
  double *probptr;
  NimArr<1, int> ansCopy;
  NimArr<1, double> probCopy;

  if(ans.isMap()) {
    if(ans.size() != K) {
      _nimble_global_output<<"Error in nimArr_rmulti: ans size does not match prob.\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }
  ansCopy.setSize(K);
  ansptr = ansCopy.getPtr();
  probptr = nimArrCopyIfNeeded<1, double>(prob, probCopy).getPtr();
  rmulti(ansptr, size, probptr, K);
  ans = ansCopy; // copy from int to double
}

double nimArr_dcat(double x, NimArr<1, double> &prob, int give_log) {
  int K = prob.size();
  double *probptr;
  NimArr<1, double> probCopy;
  probptr = nimArrCopyIfNeeded<1, double>(prob, probCopy).getPtr();
  double ans = dcat(x, probptr, K, give_log);
  return(ans);
}

double nimArr_rcat(NimArr<1, double> &prob) {
  int K = prob.size();
  double *probptr;
  NimArr<1, double> probCopy;
  probptr = nimArrCopyIfNeeded<1, double>(prob, probCopy).getPtr();
  double ans = rcat(probptr, K);
  return(ans);
}

double nimArr_ddirch(NimArr<1, double> &x, NimArr<1, double> &alpha, int give_log) {
  double *xptr, *alphaptr;
  NimArr<1, double> xCopy, alphaCopy;

  int K = alpha.size();
  if(K == 0) return(0.);
  if(x.size() != K) {
    _nimble_global_output<<"Error in nimArr_ddirch: length of x must equal length of alpha.\n";
    nimble_print_to_R(_nimble_global_output);
  }

  xptr = nimArrCopyIfNeeded<1, double>(x, xCopy).getPtr();
  alphaptr = nimArrCopyIfNeeded<1, double>(alpha, alphaCopy).getPtr();
  double ans = ddirch(xptr, alphaptr, K, give_log);
  return(ans);
}

void nimArr_rdirch(NimArr<1, double> &ans, NimArr<1, double> &alpha) {
  double *ansptr, *alphaptr;
  NimArr<1, double> ansCopy, alphaCopy;

  int K = alpha.size();
  if(K == 0) return;
  if(!ans.isMap()) {
    ans.setSize(K);
  } else {
    if(ans.size() != K) {
      _nimble_global_output<<"Error in nimArr_rdirch: ans size does not match alpha.\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }
  ansptr = nimArrCopyIfNeeded<1, double>(ans, ansCopy).getPtr();
  alphaptr = nimArrCopyIfNeeded<1, double>(alpha, alphaCopy).getPtr();
  
  rdirch(ansptr, alphaptr, K);
  if(ansptr != ans.getPtr()) {ans = ansCopy;} 
  //  if(ans.isMap()) {ans = ansCopy;}
}

double nimArr_dwish_chol(NimArr<2, double> &x, NimArr<2, double> &chol, double df, double scale_param, int give_log, int overwrite_inputs) {
  double *xptr, *cholptr;
  NimArr<2, double> xCopy, cholCopy;
  int p = x.dim()[0];
  if((x.dim()[1] != p) || (chol.dim()[0] != p) || (chol.dim()[1] != p)) {
    _nimble_global_output<<"Error in nimArr_dwish_chol: some dimensions are not right\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(df < p) {
    _nimble_global_output<<"Error in nimArr_dwish_chol: inconsistent degrees of freedom and dimension.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  xptr = nimArrCopyIfNeeded<2, double>(x, xCopy).getPtr();
  cholptr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();
  double ans = dwish_chol(xptr, cholptr, df, p, scale_param, give_log, overwrite_inputs);
  return(ans);
}


void nimArr_rwish_chol(NimArr<2, double> &ans, NimArr<2, double> &chol, double df, double scale_param, int overwrite_inputs) {
  double *ansptr, *cholptr;
  NimArr<2, double> ansCopy, cholCopy;
  int p = chol.dim()[0];
  if(chol.dim()[1] != p) {
    _nimble_global_output<<"Error in nimArr_rwish_chol: chol is not square\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(df < p) {
    _nimble_global_output<<"Error in nimArr_rwish_chol: inconsistent degrees of freedom and dimension.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(!ans.isMap()) {
    ans.setSize(p, p);
  } else {
    if((ans.dim()[0] != p) || (ans.dim()[1] != p)) {
      _nimble_global_output<<"Error in nimArr_rwish_chol: ans sizes do not match chol.\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }
  ansptr = nimArrCopyIfNeeded<2, double>(ans, ansCopy).getPtr();
  cholptr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();

  rwish_chol(ansptr, cholptr, df, p, scale_param, overwrite_inputs);
  if(ansptr != ans.getPtr()) {ans = ansCopy;} 
  //  if(ans.isMap()) {ans = ansCopy;}
}

double nimArr_dinvwish_chol(NimArr<2, double> &x, NimArr<2, double> &chol, double df, double scale_param, int give_log, int overwrite_inputs) {
  double *xptr, *cholptr;
  NimArr<2, double> xCopy, cholCopy;
  int p = x.dim()[0];
  if((x.dim()[1] != p) || (chol.dim()[0] != p) || (chol.dim()[1] != p)) {
    _nimble_global_output<<"Error in nimArr_dinvwish_chol: some dimensions are not right\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(df < p) {
    _nimble_global_output<<"Error in nimArr_dinvwish_chol: inconsistent degrees of freedom and dimension.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  xptr = nimArrCopyIfNeeded<2, double>(x, xCopy).getPtr();
  cholptr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();
  double ans = dinvwish_chol(xptr, cholptr, df, p, scale_param, give_log, overwrite_inputs);
  return(ans);
}


void nimArr_rinvwish_chol(NimArr<2, double> &ans, NimArr<2, double> &chol, double df, double scale_param, int overwrite_inputs) {
  double *ansptr, *cholptr;
  NimArr<2, double> ansCopy, cholCopy;
  int p = chol.dim()[0];
  if(chol.dim()[1] != p) {
    _nimble_global_output<<"Error in nimArr_rinvwish_chol: chol is not square\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(df < p) {
    _nimble_global_output<<"Error in nimArr_rinvwish_chol: inconsistent degrees of freedom and dimension.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(!ans.isMap()) {
    ans.setSize(p, p);
  } else {
    if((ans.dim()[0] != p) || (ans.dim()[1] != p)) {
      _nimble_global_output<<"Error in nimArr_rinvwish_chol: ans sizes do not match chol.\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }
  ansptr = nimArrCopyIfNeeded<2, double>(ans, ansCopy).getPtr();
  cholptr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();

  rinvwish_chol(ansptr, cholptr, df, p, scale_param, overwrite_inputs);
  if(ansptr != ans.getPtr()) {ans = ansCopy;} 
  //if(ans.isMap()) {ans = ansCopy;}
}


double nimArr_dmnorm_chol(NimArr<1, double> &x, NimArr<1, double> &mean, NimArr<2, double> &chol, double prec_param, int give_log, int overwrite_inputs) { 
  double *xptr, *meanptr, *cholptr;
  NimArr<1, double> xCopy, meanCopy;
  NimArr<2, double> cholCopy;
  xptr = nimArrCopyIfNeeded<1, double>(x, xCopy).getPtr();
  int n = x.size();
  meanptr = nimArrCopyIfNeeded<1, double>(mean, meanCopy).getPtr();
  if(mean.size() != n) {
    _nimble_global_output<<"Error in nimArr_dmnorm_chol: mean and x are different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  cholptr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();
  if((chol.dim()[0] != n) || (chol.dim()[1] != n)) {
    _nimble_global_output<<"Error in nimArr_dmnorm_chol: chol does not match size size of x.\n";
    nimble_print_to_R(_nimble_global_output);
  }

  double ans;
  ans = dmnorm_chol(xptr, meanptr, cholptr, n, prec_param, give_log, overwrite_inputs);
  return(ans);
}


void nimArr_rmnorm_chol(NimArr<1, double> &ans, NimArr<1, double> &mean, NimArr<2, double> &chol, double prec_param) {

  NimArr<1, double> ansCopy, meanCopy;
  NimArr<2, double> cholCopy;
  double *ansPtr, *meanPtr, *cholPtr;

  int n = mean.size();
  if(!ans.isMap()) {
    ans.setSize(n);
  } else {
    if(ans.size() != n) {
      _nimble_global_output<<"Error in nimArr_rmnorm_chol: answer size ("<< ans.size() <<") does not match mean size ("<<n<<").\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }
  ansPtr = nimArrCopyIfNeeded<1, double>(ans, ansCopy).getPtr();
  meanPtr = nimArrCopyIfNeeded<1, double>(mean, meanCopy).getPtr();
  cholPtr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();
  rmnorm_chol(ansPtr, meanPtr, cholPtr, n, prec_param);

  if(ansPtr != ans.getPtr()) {ans = ansCopy;} 
  // if(ans.isMap()) {
  //   ans = ansCopy;
  // }
}

// Begin multivariate t

double nimArr_dmvt_chol(NimArr<1, double> &x, NimArr<1, double> &mu, NimArr<2, double> &chol, double df, double prec_param, int give_log, int overwrite_inputs) { 
  
  double *xptr, *muptr, *cholptr;
  NimArr<1, double> xCopy, muCopy;
  NimArr<2, double> cholCopy;
  xptr = nimArrCopyIfNeeded<1, double>(x, xCopy).getPtr();
  int n = x.size();
  muptr = nimArrCopyIfNeeded<1, double>(mu, muCopy).getPtr();
  if(mu.size() != n) {
    _nimble_global_output<<"Error in nimArr_dmvt_chol: mu and x and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  cholptr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();
  if((chol.dim()[0] != n) || (chol.dim()[1] != n)) {
    _nimble_global_output<<"Error in nimArr_dmvt_chol: chol does not match size size of x.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  
  double ans;
  ans = dmvt_chol(xptr, muptr, cholptr, df, n, prec_param, give_log, overwrite_inputs);
  
  return(ans);
}


void nimArr_rmvt_chol(NimArr<1, double> &ans, NimArr<1, double> &mu, NimArr<2, double> &chol, double df, double prec_param) {
  
  NimArr<1, double> ansCopy, muCopy;
  NimArr<2, double> cholCopy;
  double *ansPtr, *muPtr, *cholPtr;
  
  int n = mu.size();
  if(!ans.isMap()) {
    ans.setSize(n);
  } else {
    if(ans.size() != n) {
      _nimble_global_output<<"Error in nimArr_rmvt_chol: answer size ("<< ans.size() <<") does not match mu size ("<<n<<").\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }
  ansPtr = nimArrCopyIfNeeded<1, double>(ans, ansCopy).getPtr();
  muPtr = nimArrCopyIfNeeded<1, double>(mu, muCopy).getPtr();
  cholPtr = nimArrCopyIfNeeded<2, double>(chol, cholCopy).getPtr();
  rmvt_chol(ansPtr, muPtr, cholPtr, df, n, prec_param);

  if(ansPtr != ans.getPtr()) {ans = ansCopy;} 
  // if(ans.isMap()) {
  //   ans = ansCopy;
  // }
}

// Begin LKJ

double nimArr_dlkj_corr_cholesky(NimArr<2, double> &x, double eta, int p, int give_log) { 
  
  double *xptr;
  NimArr<2, double> xCopy;
  xptr = nimArrCopyIfNeeded<2, double>(x, xCopy).getPtr();

  if((x.dim()[0] != p) || (x.dim()[1] != p)) {
    _nimble_global_output<<"Error in nimArr_dlkj_corr_cholesky: some dimensions are not right\n";
    nimble_print_to_R(_nimble_global_output);
  }
  
  double ans;
  ans = dlkj_corr_cholesky(xptr, eta, p, give_log);
  
  return(ans);
}

void nimArr_rlkj_corr_cholesky(NimArr<2, double> &ans, double eta, int p) {
  
  NimArr<2, double> ansCopy;
  double *ansPtr;
  
  if(!ans.isMap()) {
    ans.setSize(p, p);
  } else {
    if((ans.dim()[0] != p) || (ans.dim()[1] != p)) {
      _nimble_global_output<<"Error in nimArr_rlkj_corr_cholesky: ans sizes do not match p.\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }

  ansPtr = nimArrCopyIfNeeded<2, double>(ans, ansCopy).getPtr();
  rlkj_corr_cholesky(ansPtr, eta, p);

  if(ansPtr != ans.getPtr()) {ans = ansCopy;} 
}

// the next two handle when 'c' is a vector (e.g., interval censoring)
//  and the following two when 'c' is a scalar (e.g., left- and right-censoring)
double nimArr_dinterval(double x, double t, NimArr<1, double> &c, int give_log) {
  int K = c.size();
  double *cptr;
  NimArr<1, double> cCopy;
  cptr = nimArrCopyIfNeeded<1, double>(c, cCopy).getPtr();
  double ans = dinterval(x, t, cptr, K, give_log);
  return(ans);
}

int nimArr_rinterval(double t, NimArr<1, double> &c) {
  int K = c.size();
  double *cptr;
  NimArr<1, double> cCopy;
  cptr = nimArrCopyIfNeeded<1, double>(c, cCopy).getPtr();
  int ans = rinterval(t, cptr, K);
  return(ans);
}


double nimArr_dinterval(double x, double t, double c, int give_log) {
  double ans = dinterval(x, t, &c, 1, give_log);
  return(ans);
}

int nimArr_rinterval(double t, double c) {
  int ans = rinterval(t, &c, 1);
  return(ans);
}

double nimArr_dcar_normal(NimArr<1, double> &x, NimArr<1, double> &adj, NimArr<1, double> &weights, NimArr<1, double> &num, double tau, int c, int zero_mean, int give_log) {
  
  //return give_log ? ML_POSINF: R_D__0;
  
  double *xptr, *adjptr, *weightsptr, *numptr;
  NimArr<1, double> xCopy, adjCopy, weightsCopy, numCopy;
  
  xptr = nimArrCopyIfNeeded<1, double>(x, xCopy).getPtr();
  adjptr = nimArrCopyIfNeeded<1, double>(adj, adjCopy).getPtr();
  weightsptr = nimArrCopyIfNeeded<1, double>(weights, weightsCopy).getPtr();
  numptr = nimArrCopyIfNeeded<1, double>(num, numCopy).getPtr();
  
  int N = x.size();
  if(num.size() != N) {
    _nimble_global_output<<"Error in nimArr_dcar_normal: x and num and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  int L = adj.size();
  if(weights.size() != L) {
    _nimble_global_output<<"Error in nimArr_dcar_normal: adj and weights and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }

  double ans = dcar_normal(xptr, adjptr, weightsptr, numptr, tau, c, zero_mean, N, L, give_log);
  return(ans);
}


void nimArr_rcar_normal(NimArr<1, double> &ans, NimArr<1, double> &adj, NimArr<1, double> &weights, NimArr<1, double> &num, double tau, int c, int zero_mean) {
  //
  // it's important that simulation via rcar_normal() does *not* set all values to NA (or NaN),
  // since initializeModel() will call this simulate method if there are any NA's present,
  // (which is allowed for island components), which over-writes all the other valid initial values.
  //
}


double nimArr_dcar_proper(NimArr<1, double> &x, NimArr<1, double> &mu, NimArr<1, double> &C, NimArr<1, double> &adj, NimArr<1, double> &num, NimArr<1, double> &M, double tau, double gamma, NimArr<1, double> &evs, int give_log) {
  
  //return give_log ? ML_POSINF: R_D__0;
  
  double *xptr, *muptr, *Cptr, *adjptr, *numptr, *Mptr, *evsptr;
  NimArr<1, double> xCopy, muCopy, CCopy, adjCopy, numCopy, MCopy, evsCopy;
  
  xptr = nimArrCopyIfNeeded<1, double>(x, xCopy).getPtr();
  muptr = nimArrCopyIfNeeded<1, double>(mu, muCopy).getPtr();
  Cptr = nimArrCopyIfNeeded<1, double>(C, CCopy).getPtr();
  adjptr = nimArrCopyIfNeeded<1, double>(adj, adjCopy).getPtr();
  numptr = nimArrCopyIfNeeded<1, double>(num, numCopy).getPtr();
  Mptr = nimArrCopyIfNeeded<1, double>(M, MCopy).getPtr();
  evsptr = nimArrCopyIfNeeded<1, double>(evs, evsCopy).getPtr();
  
  int N = x.size();
  if(mu.size() != N) {
    _nimble_global_output<<"Error in nimArr_dcar_proper: x and mu and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(num.size() != N) {
    _nimble_global_output<<"Error in nimArr_dcar_proper: x and num and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(M.size() != N) {
    _nimble_global_output<<"Error in nimArr_dcar_proper: x and M and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(evs.size() != N) {
    _nimble_global_output<<"Error in nimArr_dcar_proper: x and evs and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }
  int L = adj.size();
  if(C.size() != L) {
    _nimble_global_output<<"Error in nimArr_dcar_proper: adj and C and different sizes.\n";
    nimble_print_to_R(_nimble_global_output);
  }

  double ans = dcar_proper(xptr, muptr, Cptr, adjptr, numptr, Mptr, tau, gamma, evsptr, N, L, give_log);
  return(ans);
}


void nimArr_rcar_proper(NimArr<1, double> &ans, NimArr<1, double> &mu, NimArr<1, double> &C, NimArr<1, double> &adj, NimArr<1, double> &num, NimArr<1, double> &M, double tau, double gamma, NimArr<1, double> &evs) {

  double *ansptr, *muptr, *Cptr, *adjptr, *numptr, *Mptr, *evsptr;
  NimArr<1, double> ansCopy, muCopy, CCopy, adjCopy, numCopy, MCopy, evsCopy;
  
  int N = num.size();
  int L = adj.size();
  
  if(!ans.isMap()) {
    ans.setSize(N);
  } else {
    if(ans.size() != N) {
      _nimble_global_output<<"Error in nimArr_rcar_proper: answer size ("<< ans.size() <<") does not match num size ("<<N<<").\n";
      nimble_print_to_R(_nimble_global_output);
    }
  }
  if(mu.size() != N) {
    _nimble_global_output<<"Error in nimArr_rcar_proper: mu size ("<< mu.size() <<") does not match num size ("<<N<<").\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(C.size() != L) {
    _nimble_global_output<<"Error in nimArr_rcar_proper: C size ("<< C.size() <<") does not match adj size ("<<L<<").\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(M.size() != N) {
    _nimble_global_output<<"Error in nimArr_rcar_proper: M size ("<< M.size() <<") does not match num size ("<<N<<").\n";
    nimble_print_to_R(_nimble_global_output);
  }
  if(evs.size() != N) {
    _nimble_global_output<<"Error in nimArr_rcar_proper: evs size ("<< evs.size() <<") does not match num size ("<<N<<").\n";
    nimble_print_to_R(_nimble_global_output);
  }
  
  ansptr = nimArrCopyIfNeeded<1, double>(ans, ansCopy).getPtr();
  muptr = nimArrCopyIfNeeded<1, double>(mu, muCopy).getPtr();
  Cptr = nimArrCopyIfNeeded<1, double>(C, CCopy).getPtr();
  adjptr = nimArrCopyIfNeeded<1, double>(adj, adjCopy).getPtr();
  numptr = nimArrCopyIfNeeded<1, double>(num, numCopy).getPtr();
  Mptr = nimArrCopyIfNeeded<1, double>(M, MCopy).getPtr();
  evsptr = nimArrCopyIfNeeded<1, double>(evs, evsCopy).getPtr();
  
  rcar_proper(ansptr, muptr, Cptr, adjptr, numptr, Mptr, tau, gamma, evsptr, N, L);

  if(ansptr != ans.getPtr()) {ans = ansCopy;} 
}

