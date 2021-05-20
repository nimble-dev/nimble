#ifndef _NIMDERIVS_MATMULT_UTILS
#define _NIMDERIVS_MATMULT_UTILS

// copied from include/EigenTypeDefs.h, to be found there in full use.
#include "EigenTypedefs.h"
typedef typename EigenTemplateTypes<CppAD::AD<double> >::typeMatrixXd MatrixXd_CppAD;

template<typename MT, typename VT>
void mat2vec(const MT &m, VT &v, size_t offset = 0) {
  size_t nrow = m.rows();
  size_t ncol = m.cols();
  for(size_t i=0; i < nrow; ++i) {
    for(size_t j=0; j < ncol; ++j) {
      v[i + j*nrow + offset] = m(i,j);
    }
  }
}

template<typename MT, typename VT>
void mat2vec_v(const MT &m, VT &v, size_t offset = 0) {
  size_t nrow = m.rows();
  size_t ncol = m.cols();
  for(size_t i=0; i < nrow; ++i) {
    for(size_t j=0; j < ncol; ++j) {
      v[i + j*nrow + offset] = CppAD::Value(m(i,j));
    }
  }
}

template<typename VT, typename MT>
void vec2mat(const VT &v, MT &m, size_t offset = 0) {
  size_t nrow = m.rows();
  size_t ncol = m.cols();
  for(size_t i=0; i < nrow; ++i) {
    for(size_t j=0; j < ncol; ++j) {
      m(i,j) = v[i + j*nrow + offset];
    }
  }
}


#endif
