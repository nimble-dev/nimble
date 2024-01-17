#ifndef _NIMDERIVS_ATOMIC_DYN_IND
#define _NIMDERIVS_ATOMIC_DYN_IND

#include <cppad/cppad.hpp>
#include <cppad/utility/nan.hpp>
#include "nimbleCppAD.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

atomic_dyn_ind_get_class* new_atomic_dyn_ind_get(void* tape_mgr, const std::string& name);
void delete_atomic_dyn_ind_get(void* tape_mgr, atomic_dyn_ind_get_class *atomic_dyn_ind_get);
atomic_dyn_ind_set_class* new_atomic_dyn_ind_set(void* tape_mgr, const std::string& name);
void delete_atomic_dyn_ind_set(void* tape_mgr, atomic_dyn_ind_set_class *atomic_dyn_ind_set);

CppAD::AD<double> dyn_ind_get(const CppAD::vector<CppAD::AD<double> > &x,
                              const CppAD::AD<double> &index,
                              const size_t offset=0,
                              const size_t nrow=1, const int nx_=-1);
void dyn_ind_set(CppAD::vector<CppAD::AD<double> > &y,
                 const CppAD::AD<double> &index,
                 const CppAD::AD<double> &x,
                 const size_t offset=0,
                 const size_t nrow=1, const int ny_=-1);

/////////////
// Compile-time type processing for dyn_ind_get and dyn_ind_set templates
template<class I_>
struct IsCppAD {
  // This case catches all NON-CppAD
  static constexpr bool no=true;
  // omit yes
  static constexpr bool ans=false;
};

template<class S_>
struct IsCppAD<CppAD::AD<S_> > {
  // This catches all yes CppAD
  // omit no
  static constexpr bool yes=true;
  static constexpr bool ans=true;
};

template<class S_>
struct IsCppAD<const CppAD::AD<S_> > {
  static constexpr bool yes=true;
  static constexpr bool ans=true;
};

template<bool TF>
struct dyn_ind_truefalse {
  static constexpr bool value = true;
};
template<>
struct dyn_ind_truefalse<false> {
  static constexpr bool value = false;
};


template<class I1, class I2>
struct AnyCppAD2D {
  typedef CppAD::AD<double> type;
  // make viable if I1 or I2 is a CppAD type
  static constexpr bool value =
    std::conditional<IsCppAD<I1>::ans,
                     dyn_ind_truefalse<true>,
                     typename std::conditional<IsCppAD<I2>::ans,
                                               dyn_ind_truefalse<true>,
                                               dyn_ind_truefalse<false> >::type >::type::value;
};


template<class I1, class I2, class I3>
struct AnyCppAD3D {
  typedef CppAD::AD<double> type;
  // make viable if I1 or I2 is a CppAD type
  static constexpr bool value =
    std::conditional<IsCppAD<I1>::ans,
                     dyn_ind_truefalse<true>,
                     typename std::conditional<IsCppAD<I2>::ans,
                                               dyn_ind_truefalse<true>,
                                               typename std::conditional<IsCppAD<I3>::ans,
                                                                         dyn_ind_truefalse<true>,
                                                                         dyn_ind_truefalse<false> >::type >::type>::type::value;
};

// Extract a value whether from a CppAD or non-CppAD type
template<class T>
struct CppADvalue {
  static size_t v(const T &x){return static_cast<size_t>(x);}
};

template<class T>
struct CppADvalue<CppAD::AD<T> > {
  static size_t v(const CppAD::AD<T> &x){return static_cast<size_t>(CppAD::Value(x));}
};

///////////////////
// dyn_ind_get
//////////////////
// 1D get
// Non-CppAD case for t[i] when t is a pointer. Neither t nor i are CppAD.
template<class T_, class I_ >
inline typename std::conditional<IsCppAD<I_>::no, T_, void>::type
stoch_ind_get(const T_ *t, const I_ &i) {
  //std::cout<< "pointer get case"<<std::endl;
  return t[i];
}

// 1D get
// index is not CppAD type
template<class I_>
typename std::conditional<IsCppAD<I_>::no,
                          CppAD::AD<double>,
                          void>::type
stoch_ind_get(const NimArr<1, CppAD::AD<double> > &x, const I_ &i ) {
  //std::cout<<"in template 1D stoch_ind_get"<<std::endl;
  return x[i];
}

// 1D get
// index is a CppAD type
// (This could be directly coded with I_ set to CppAD::AD<double>,
//  but we put it in template form to get that right building to higher dimensional cases.)
template<class I_>
inline typename std::conditional<IsCppAD<I_>::yes,
                                 CppAD::AD<double>,
                                 void>::type
stoch_ind_get(const NimArr<1, CppAD::AD<double> > &x,
              const I_ &i ) {
  //std::cout<<"in CppAD 1D stoch_ind_get"<<std::endl;
  if(CppAD::Constant(i)) return x[CppADvalue<I_>::v(i)];
  CppAD::vector< CppAD::AD<double> > x_(x.size());
  // x is 1D for this case.
  for(size_t i = 0; i < x.size(); ++i) x_[i] = x[i];
  return dyn_ind_get(x_, i); // dummy version: x[CppAD::Value(i)];
}

// 2D get
// i1 and i2 are not CppAD
template<class I1_, class I2_,
         typename std::enable_if<!AnyCppAD2D<I1_, I2_>::value, int>::type = 1>
CppAD::AD<double> stoch_ind_get(const NimArr<2, CppAD::AD<double> > &x,
                                const I1_ &i1,
                                const I2_ &i2) {
  //std::cout<<"in template 2D stoch_ind_get"<<std::endl;
  return x(i1, i2);
}

// 2D get
// at least one of i1 or i2 is CppAD
template<class I1_, class I2_,
         typename std::enable_if<AnyCppAD2D<I1_, I2_>::value, int>::type = 1>
CppAD::AD<double> stoch_ind_get(const NimArr<2, CppAD::AD<double> > &x,
                                const I1_ &i1,
                                const I2_ &i2) {
  //std::cout<<"in CppAD 2D stoch_ind_get"<<std::endl;
  //if(x.isMap()) std::cout<<"have not implemented NimArr map case yet."<<std::endl;
  const int *dim = x.dim();
  CppAD::AD<double> flat_i = i1 + i2*dim[0];
  if(CppAD::Constant(flat_i)) return x(CppADvalue<I1_>::v(i1), CppADvalue<I2_>::v(i2));
  CppAD::vector< CppAD::AD<double> > x_(x.size());
  const int *xdim = x.dim();
  size_t new_i = 0;
  for(size_t j2 = 0; j2 < xdim[1]; ++j2 ) {
    for(size_t j1 = 0; j1 < xdim[0]; ++j1) {
      x_[new_i++] = x(j1, j2);
    }
  }
  return dyn_ind_get(x_, flat_i);
}

// 3D get
// i1, i2, and i3 and not CppAD
template<class I1_, class I2_, class I3_,
         typename std::enable_if<!AnyCppAD3D<I1_, I2_, I3_>::value, int>::type = 1>
CppAD::AD<double> stoch_ind_get(const NimArr<3, CppAD::AD<double> > &x,
                                const I1_ &i1,
                                const I2_ &i2,
                                const I3_ &i3) {
  //std::cout<<"in template 2D stoch_ind_get"<<std::endl;
  return x(i1, i2, i3);
}


// 3D get
// at least one index is CppAD
template<class I1_, class I2_, class I3_,
         typename std::enable_if<AnyCppAD3D<I1_, I2_, I3_>::value, int>::type = 1>
CppAD::AD<double> stoch_ind_get(const NimArr<3, CppAD::AD<double> > &x,
                                const I1_ &i1,
                                const I2_ &i2,
                                const I3_ &i3) {
  //std::cout<<"in CppAD 2D stoch_ind_get"<<std::endl;
  if(x.isMap()) std::cout<<"have not implemented NimArr map case yet."<<std::endl;
  const int *strides = x.strides();
  CppAD::AD<double> flat_i = i1*strides[0] + i2*strides[1] + i3*strides[2]; // ignore x.offset because map case needs special handling anyway
  if(CppAD::Constant(flat_i)) return x(CppADvalue<I1_>::v(i1), CppADvalue<I2_>::v(i2), CppADvalue<I3_>::v(i3));
  CppAD::vector< CppAD::AD<double> > x_(x.size());
  const int *xdim = x.dim();
  size_t new_i = 0;
  for(size_t j3 = 0; j3 < xdim[2]; ++j3 ) {
    for(size_t j2 = 0; j2 < xdim[1]; ++j2 ) {
      for(size_t j1 = 0; j1 < xdim[0]; ++j1) {
        x_[new_i++] = x(j1, j2, j3);
      }
    }
  }
  return dyn_ind_get(x_, flat_i);
}

///////////////////
// stoch_ind_set

// 1D set
// base class
template<class I_>
class stoch_ind_set1_base_c {
  public:
  typedef NimArr<1, CppAD::AD<double> > ArrADd;
  ArrADd *x_ptr;
  stoch_ind_set1_base_c(ArrADd &x_, const I_ &i_) : x_ptr(&x_), i(i_) {}
  I_ i;
};

// 1D set
// index is not CppAD
template<class I_>
class stoch_ind_set1_c : public stoch_ind_set1_base_c<I_> {
  using stoch_ind_set1_base_c<I_>::stoch_ind_set1_base_c;
  using stoch_ind_set1_base_c<I_>::x_ptr;
  using stoch_ind_set1_base_c<I_>::i;
  public:
  CppAD::AD<double> operator=(const CppAD::AD<double> &v) {
    //std::cout<<"In operator= for template 1D stoch_ind_set_c"<<std::endl;
    (*x_ptr)[i] = v;
    return v;
  }
};

// 1D set
// index is CppAD
template<>
class stoch_ind_set1_c<CppAD::AD<double> > : public stoch_ind_set1_base_c<CppAD::AD<double> > {
  using stoch_ind_set1_base_c<CppAD::AD<double> >::stoch_ind_set1_base_c;
  public:
  CppAD::AD<double> operator=(const CppAD::AD<double> &v) {
    //std::cout<<"In operator= for CppAD 1D stoch_ind_set_c"<<std::endl;
    // kluge: assume x is 1D
    if(CppAD::Constant(i)) {
      (*x_ptr)[CppAD::Value(i)] = v;
    } else {
      CppAD::vector< CppAD::AD<double> > x_((*x_ptr).size());
      for(size_t j = 0; j < x_.size(); ++j) x_[j] = (*x_ptr)[j];
      dyn_ind_set(x_, i, v);
      for(size_t j = 0; j < x_.size(); ++j) (*x_ptr)[j] = x_[j];
    }
    return v;
  }
};

// 1D set
// function to catch types and dispatch to classes
template<class I_>
stoch_ind_set1_c<I_> stoch_ind_set(NimArr<1, CppAD::AD<double> > &x, const I_ &i) {
  return stoch_ind_set1_c<I_>(x, i);
}

// 2D set
// base class
template<class I1_, class I2_>
class stoch_ind_set2_base_c {
  public:
  typedef NimArr<2, CppAD::AD<double> > ArrADd;
  ArrADd *x_ptr;
  stoch_ind_set2_base_c(ArrADd &x_, const I1_ &i1_, const I2_ &i2_) : x_ptr(&x_), i1(i1_), i2(i2_) {}
  I1_ i1; I2_ i2;
};

// 2D set
// template base class for specialization
template<class I1_, class I2_, bool UseCppAD = false>
class stoch_ind_set2_c;

// 2D set
// no indicies are CppAD
template<class I1_, class I2_>
class stoch_ind_set2_c<I1_, I2_, false> : public stoch_ind_set2_base_c<I1_, I2_> {
  using stoch_ind_set2_base_c<I1_, I2_>::stoch_ind_set2_base_c;
  using stoch_ind_set2_base_c<I1_, I2_>::x_ptr;
  using stoch_ind_set2_base_c<I1_, I2_>::i1;
  using stoch_ind_set2_base_c<I1_, I2_>::i2;
  public:
  CppAD::AD<double> operator=(const CppAD::AD<double> &v) {
    //std::cout<<"In operator= for template 2D stoch_ind_set_c"<<std::endl;
    (*x_ptr)(i1, i2) = v;
    return v;
  }
};

// 2D set
// at least one index is CppAD
template<class I1_, class I2_>
class stoch_ind_set2_c<I1_, I2_, true> : public stoch_ind_set2_base_c<I1_, I2_> {
  using stoch_ind_set2_base_c<I1_, I2_>::stoch_ind_set2_base_c;
  using stoch_ind_set2_base_c<I1_, I2_>::x_ptr;
  using stoch_ind_set2_base_c<I1_, I2_>::i1;
  using stoch_ind_set2_base_c<I1_, I2_>::i2;
  public:
  CppAD::AD<double> operator=(const CppAD::AD<double> &v) {
    //std::cout<<"In operator= for CppAD 2D stoch_ind_set_c"<<std::endl;
    //if((*x_ptr).isMap()) std::cout<<"have not implemented NimArr map case yet."<<std::endl;
    const int *dim = (*x_ptr).dim();
    CppAD::AD<double> flat_i = i1 + i2*dim[0]; // ignore x.offset because map case needs special handling anyway
    if(CppAD::Constant(flat_i)) {
      (*x_ptr)(CppADvalue<I1_>::v(i1), CppADvalue<I2_>::v(i2)) = v;
    } else {
      const int *xdim = (*x_ptr).dim();
      size_t new_i = 0;
      CppAD::vector< CppAD::AD<double> > x_((*x_ptr).size());
      for(size_t j2 = 0; j2 < xdim[1]; ++j2 ) {
        for(size_t j1 = 0; j1 < xdim[0]; ++j1) {
          x_[new_i++] = (*x_ptr)(j1, j2);
        }
      }
      dyn_ind_set(x_, flat_i, v);
      new_i = 0;
      for(size_t j2 = 0; j2 < xdim[1]; ++j2 ) {
        for(size_t j1 = 0; j1 < xdim[0]; ++j1) {
          (*x_ptr)(j1, j2) = x_[new_i++];
        }
      }
    }
    return v;
  }
};

// 2D set
// function to catch types and dispatch to classes
template<class I1_, class I2_>
stoch_ind_set2_c<I1_, I2_, AnyCppAD2D<I1_, I2_>::value>
stoch_ind_set(NimArr<2, CppAD::AD<double> > &x,
                                         const I1_ &i1,
                                         const I2_ &i2) {
  return stoch_ind_set2_c<I1_, I2_, AnyCppAD2D<I1_, I2_>::value>(x, i1, i2);
}


///////////////
// 3D set

// 3D set
// base class
template<class I1_, class I2_, class I3_>
class stoch_ind_set3_base_c {
  public:
  typedef NimArr<3, CppAD::AD<double> > ArrADd;
  ArrADd *x_ptr;
  stoch_ind_set3_base_c(ArrADd &x_, const I1_ &i1_, const I2_ &i2_, const I3_ &i3_) : x_ptr(&x_), i1(i1_), i2(i2_), i3(i3_) {}
  I1_ i1; I2_ i2; I3_ i3;
};

// 3D set
// template base class for specialization
template<class I1_, class I2_, class I3_, bool UseCppAD = false>
class stoch_ind_set3_c;

// 3D set
// no indicies are CppAD
template<class I1_, class I2_, class I3_>
class stoch_ind_set3_c<I1_, I2_, I3_, false> : public stoch_ind_set3_base_c<I1_, I2_, I3_> {
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::stoch_ind_set3_base_c;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::x_ptr;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::i1;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::i2;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::i3;
  public:
  CppAD::AD<double> operator=(const CppAD::AD<double> &v) {
    //std::cout<<"In operator= for template 2D stoch_ind_set_c"<<std::endl;
    (*x_ptr)(i1, i2, i3) = v;
    return v;
  }
};

// 3D set
// at least one index is CppAD
template<class I1_, class I2_, class I3_>
class stoch_ind_set3_c<I1_, I2_, I3_, true> : public stoch_ind_set3_base_c<I1_, I2_, I3_> {
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::stoch_ind_set3_base_c;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::x_ptr;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::i1;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::i2;
  using stoch_ind_set3_base_c<I1_, I2_, I3_>::i3;
  public:
  CppAD::AD<double> operator=(const CppAD::AD<double> &v) {
    //std::cout<<"In operator= for CppAD 2D stoch_ind_set_c"<<std::endl;
    if((*x_ptr).isMap()) std::cout<<"have not implemented NimArr map case yet."<<std::endl;
    const int *strides = (*x_ptr).strides();
    CppAD::AD<double> flat_i = i1*strides[0] + i2*strides[1] + i3*strides[2]; // ignore x.offset because map case needs special handling anyway
    if(CppAD::Constant(flat_i)) {
      (*x_ptr)(CppADvalue<I1_>::v(i1), CppADvalue<I2_>::v(i2), CppADvalue<I3_>::v(i3)) = v;
    } else {
      const int *xdim = (*x_ptr).dim();
      size_t new_i = 0;
      CppAD::vector< CppAD::AD<double> > x_((*x_ptr).size());
      for(size_t j3 = 0; j3 < xdim[2]; ++j3 ) {
        for(size_t j2 = 0; j2 < xdim[1]; ++j2 ) {
          for(size_t j1 = 0; j1 < xdim[0]; ++j1) {
            x_[new_i++] = (*x_ptr)(j1, j2, j3);
          }
        }
      }
      dyn_ind_set(x_, flat_i, v);
      new_i = 0;
      for(size_t j3 = 0; j3 < xdim[2]; ++j3 ) {
        for(size_t j2 = 0; j2 < xdim[1]; ++j2 ) {
          for(size_t j1 = 0; j1 < xdim[0]; ++j1) {
            (*x_ptr)(j1, j2, j3) = x_[new_i++];
          }
        }
      }
    }
    return v;
  }
};

// 3D set
// function to catch types and dispatch to classes
template<class I1_, class I2_, class I3_>
stoch_ind_set3_c<I1_, I2_, I3_, AnyCppAD3D<I1_, I2_, I3_>::value>
stoch_ind_set(NimArr<3, CppAD::AD<double> > &x,
              const I1_ &i1,
              const I2_ &i2,
              const I3_ &i3) {
  return stoch_ind_set3_c<I1_, I2_, I3_, AnyCppAD3D<I1_, I2_, I3_>::value>(x, i1, i2, i3);
}

class atomic_dyn_ind_get_class : public CppAD::atomic_three< double >, public nimble_atomic_base {
  /*
    This is for y = x[i].
    The type_x, taylor_x, depend_x, etc will be length(x) + 1, with last element set to i.
  */
 public:
 atomic_dyn_ind_get_class(const std::string& name);
 private:
  virtual bool for_type(
                        const CppAD::vector<double>&               parameter_x ,
                        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                        CppAD::vector<CppAD::ad_type_enum>&        type_y      );

  virtual bool rev_depend(
                          const CppAD::vector<double>&          parameter_x ,
                          const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                          CppAD::vector<bool>&                depend_x    ,
                          const CppAD::vector<bool>&          depend_y
                          );

  virtual bool forward(
                       const CppAD::vector<double>&               parameter_x  ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                       size_t                              need_y       ,
                       size_t                              order_low    ,
                       size_t                              order_up     ,
                       const CppAD::vector<double>&               taylor_x     ,
                       CppAD::vector<double>&                     taylor_y     );

  virtual bool forward(
                       const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                       size_t                              need_y       ,
                       size_t                              order_low    ,
                       size_t                              order_up     ,
                       const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
                       CppAD::vector<CppAD::AD<double> >&                     taylor_y     );

  virtual bool reverse(
                       const CppAD::vector<double>&               parameter_x ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                       size_t                              order_up    ,
                       const CppAD::vector<double>&               taylor_x    ,
                       const CppAD::vector<double>&               taylor_y    ,
                       CppAD::vector<double>&                     partial_x   ,
                       const CppAD::vector<double>&               partial_y   );

  virtual bool reverse(
                       const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                       size_t                              order_up    ,
                       const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
                       const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
                       CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
                       const CppAD::vector<CppAD::AD<double> >&               partial_y   );
};

class atomic_dyn_ind_set_class : public CppAD::atomic_three< double >, public nimble_atomic_base {
  /*
    This is for y[i] = x.
    The type_x, taylor_x, depend_x, etc will have length = length(y) +  2, with elements y, then x (scalar), then index.
    This is very much like how in R, "y[i] = x" is really "y = `[<-`(y, i, x)"
    We will turn use "dyn_ind_set(y, i, x)" for "y[i] = x" and it will really do "<object of this class>(c(y, i, x), y)"
    The inputs for y will have length(y).
    Note that this operation is more costly in AD tape operations than dyn_ind_get, because y is duplicated every time.
  */
 public:
 atomic_dyn_ind_set_class(const std::string& name);
 private:
  virtual bool for_type(
                        const CppAD::vector<double>&               parameter_x ,
                        const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                        CppAD::vector<CppAD::ad_type_enum>&        type_y      );

  virtual bool rev_depend(
                          const CppAD::vector<double>&          parameter_x ,
                          const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                          CppAD::vector<bool>&                depend_x    ,
                          const CppAD::vector<bool>&          depend_y
                          );
  virtual bool forward(
                       const CppAD::vector<double>&               parameter_x  ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                       size_t                              need_y       ,
                       size_t                              order_low    ,
                       size_t                              order_up     ,
                       const CppAD::vector<double>&               taylor_x     ,
                       CppAD::vector<double>&                     taylor_y     );

  virtual bool forward(
                       const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
                       size_t                              need_y       ,
                       size_t                              order_low    ,
                       size_t                              order_up     ,
                       const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
                       CppAD::vector<CppAD::AD<double> >&                     taylor_y     );

  virtual bool reverse(
                       const CppAD::vector<double>&               parameter_x ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                       size_t                              order_up    ,
                       const CppAD::vector<double>&               taylor_x    ,
                       const CppAD::vector<double>&               taylor_y    ,
                       CppAD::vector<double>&                     partial_x   ,
                       const CppAD::vector<double>&               partial_y   );

  virtual bool reverse(
                       const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
                       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
                       size_t                              order_up    ,
                       const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
                       const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
                       CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
                       const CppAD::vector<CppAD::AD<double> >&               partial_y   );
};

#endif // _NIMDERIVS_ATOMIC_DYN_IND
