#ifndef TINY_AD_TWEEDIE_H
#define TINY_AD_TWEEDIE_H

#include "../gamma/gamma.hpp"

namespace tweedie_utils {
using gamma_utils::lgammafn;

template<class Float>
Float dmax (Float *x, int n){
  Float s = x[0] ;
  for (int i = 1; i < n; i++)
      if (x[i] > s) s = x[i] ; 
  return s ;
}
#define imax dmax

template<class S, class T>
T fmax2(S x, T y) {return (x < y) ? y : x;}

template<class S>
int imax2(S x, int y) {return (x < y) ? y : x;}
template<class S>
int imin2(S x, int y) {return (x < y) ? x : y;}

#ifndef Calloc
#define Calloc(n,type) (type*)calloc(n, sizeof(type))
#endif
#ifndef Free
#define Free free
#endif
#include "tweedie.cpp"
#undef TWEEDIE_DROP
#undef TWEEDIE_INCRE
#undef TWEEDIE_NTERM
#undef imax

template<class Float>
Float dtweedie(Float y, Float mu, Float phi, Float p, int give_log = 0) {
  Float ans;
  bool ok = (0 <= y) & (0 < mu) & (0 < phi) & (1 < p) & (p < 2);
  if (!ok) return NAN;
  dtweedie<Float>(1, &y, &mu, phi, p, NULL, &ans);
  return ( give_log ? ans : exp(ans) );
}

} // End namespace tweedie_utils

#endif
