#ifndef _NIMDERIVS_ATOMIC_CLASSES
#define _NIMDERIVS_ATOMIC_CLASSES

// See end of this file for global object definitions and functions to call

template<class T>
class unary_atomic_class : public CppAD::atomic_three<T> {
  // This layer in the class hierarchy simply provides the same for_type and rev_depend
  // for all atomic classes representing unary functions below.
 public:
  unary_atomic_class() {};
 private:
   // for_type is essential.
  // This determines which elements of y are constants, dynamic parameters, or variables
  //   depending on types of x.
  // If omitted, calls to forward (and possibly reverse) will not happen.
  virtual bool for_type(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
  { 
    type_y[0] = type_x[0];
    return true;
  }
  // rev_depend is used when the optimize() method is called for a tape (an ADFun).
  virtual bool rev_depend(
			  const CppAD::vector<double>&          parameter_x ,
			  const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
			  CppAD::vector<bool>&                depend_x    ,
			  const CppAD::vector<bool>&          depend_y
			  ) {
    depend_x[0] = depend_y[0];
    return true;
  }
}

class atomic_lgamma_class : public unary_atomic_class<double>
  public:
  // From atomic_three_get_started
  atomic_lgamma_class(const std::string& name) : 
    CppAD::atomic_three<double>(name)        
    { }

private:
  virtual bool forward(
      const CppAD::vector<double>&               parameter_x  ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
      size_t                              need_y       ,
      size_t                              order_low    ,
      size_t                              order_up     ,
      const CppAD::vector<double>&               taylor_x     ,
      CppAD::vector<double>&                     taylor_y     )
  {
    if(order_low <= 0 & order_up >= 0) {
      taylor_y[0] = Rf_lgammafn(taylor_x[0]);
    }
    double fprime;
    if(order_low > 0) fprime = Rf_psigamma(taylor_x[0], 0);
    if(order_low <= 1 & order_up >= 1) {
      taylor_y[1] = fprime * taylor_x[1];
      // f'(x) (x')
    }
    if(order_low <= 2 & order_up >= 2) {
      taylor_y[2] = 0.5 * (Rf_psigamma(taylor_x[0], 1) * taylor_x[1] * taylor_x[1] + 
        fprime * 2 * taylor_x[2]);
      // 0.5 * ((f''(x)) (x')^2 + 2*f'(x) (x''))
      // Note x'' is a taylor coeff so it is 0.5*2nd deriv
    }
    return true;
  }
  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   )
  {
    partial_x[0] = 0;
    double fprime = Rf_psigamma(taylor_x[0], 0);
    if(order_up >= 1) {
      partial_x[1] += partial_y[1] * fprime;
      partial_x[0] += partial_y[1] * Rf_psigamma(taylor_x[0], 1) * taylor_x[1];
    }
    partial_x[0] += partial_y[0] * fprime;
    // dG/dx = dG/dy  dy/dx 
    return true;
  }
};


/******************************/

class atomic_pnorm1_class :public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  atomic_pnorm1_class(const std::string& name) : 
    CppAD::atomic_three<double>(name)        
    { }

private:
   virtual bool forward(
      const CppAD::vector<double>&               parameter_x  ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
      size_t                              need_y       ,
      size_t                              order_low    ,
      size_t                              order_up     ,
      const CppAD::vector<double>&               taylor_x     ,
      CppAD::vector<double>&                     taylor_y     )
  {
    // This can do more with checking sizes, need_y, orders, and
    // what is a parameter.
    // need_y determines whether all elements need calculation,
    // or only those for parameters, or only those for variables.
    if(order_low == 0) {
      taylor_y[0] = Rf_pnorm5(taylor_x[0], 0, 1, 1, 0); // (x, mu, sigma, lower_tail, log)
      // F(x)
    }
    double Fprime(0.);
    if(order_low <= 1 && order_up >= 1) {
      Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
      taylor_y[1] = Fprime * taylor_x[1];
      // F'(x) (x')
    }
    if(order_low <= 2 && order_up >= 2) {
      if(Fprime == 0) {     // only calculate Fprime if it is really needed
	if(taylor_x[2] != 0)
	  Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
      }
      taylor_y[2] = 0.5 * (- taylor_x[0]) * taylor_y[1] * taylor_x[1] +
	Fprime*taylor_x[2]; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
    }
    return true;
  }

  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   )
  {
    double Fprime = Rf_dnorm4(taylor_x[0], 0, 1, 0);
    partial_x[0] = 0;
    if(order_up >= 1) partial_x[1] = 0;
    if(order_up >= 2) return false; // not implemented
    if(order_up >= 1) { // Needed for 2nd derivs
      partial_x[0] += partial_y[1] * (- taylor_x[0]) * Fprime * taylor_x[1];
      partial_x[1] += partial_y[1] * Fprime;
    } 
    // Needed for 1st derivs
    partial_x[0] += partial_y[0] * Fprime;
    return true;
  }
};

/*******************************************/

class atomic_qnorm1_class : public unary_atomic_class<double> {
  public:
  // From atomic_three_get_started
  atomic_qnorm1_class(const std::string& name) : 
    CppAD::atomic_three<double>(name)        
    { }

private:
  virtual bool forward(
      const CppAD::vector<double>&               parameter_x  ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
      size_t                              need_y       ,
      size_t                              order_low    ,
      size_t                              order_up     ,
      const CppAD::vector<double>&               taylor_x     ,
      CppAD::vector<double>&                     taylor_y     )
  {
    // This can do more with checking sizes, need_y, orders, and
    // what is a parameter.
    // need_y determines whether all elements need calculation,
    // or only those for parameters, or only those for variables.
    if(order_low == 0) {
      taylor_y[0] = Rf_qnorm5(taylor_x[0], 0, 1, 1, 0); // (p, mu, sigma, lower_tail, log)
      // F(x)
    }
    double z = taylor_y[0];
    double invFprime(0.);
    if(order_low <= 1 && order_up >= 1) {
      invFprime = Rf_dnorm4(z, 0, 1, 0);
      taylor_y[1] = taylor_x[1] / invFprime;
      // F'(x) (x')
    }
    if(order_low <= 2 && order_up >= 2) {
      taylor_y[2] = 0.5 * z * taylor_y[1] * taylor_y[1];
      if(taylor_x[2] != 0) {
        if(invFprime == 0)
          invFprime = Rf_dnorm4(z, 0, 1, 0);
        taylor_y[2] += taylor_x[2] / invFprime; // 0.5 * (-x) * F'(x) (x') * (x') + 0.5 * 2 * F'(x) * (x'')
      }
    }
    return true;
  }

  virtual bool reverse(
      const CppAD::vector<double>&               parameter_x ,
      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
      size_t                              order_up    ,
      const CppAD::vector<double>&               taylor_x    ,
      const CppAD::vector<double>&               taylor_y    ,
      CppAD::vector<double>&                     partial_x   ,
      const CppAD::vector<double>&               partial_y   )
  {
    double z = taylor_y[0];
    double invFprime = Rf_dnorm4(z, 0, 1, 0);
    partial_x[0] = 0;
    if(order_up >= 1) partial_x[1] = 0;
    if(order_up >= 2) return false; // not implemented
    if(order_up >= 1) { // Needed for 2nd derivs
      partial_x[0] += partial_y[1] * (z/(invFprime*invFprime)) * taylor_x[1];
      partial_x[1] += partial_y[1] / invFprime;
    } 
    // Needed for 1st derivs
    partial_x[0] += partial_y[0] / invFprime;
    return true;
  }
};

#endif

template<class T>
T atomic_lfactorial(T x) {
  return atomic_lgamma(x + 1);
}

template<class T>
T atomic_lgamma(T x) {
  static atomic_lgamma_class static_atomic_lgamma("atomic_lgamma");
  CppAD::vector<T> in(1);
  CppAD::vector<T> out(1);
  in[0] = x;
  static_atomic_lgamma(in, out);
  return out[0];
}

template<class T>
T atomic_pnorm1(T x) {
  atomic_pnorm1_class static_atomic_pnorm1("atomic_pnorm1");
  CppAD::vector<T> in(1);
  CppAD::vector<T> out(1);
  in[0] = x;
  static_atomic_pnorm1(in, out);
  return out[0];
}

template<class T>
T atomic_qnorm1(T x) {
  atomic_qnorm1_class static_atomic_qnorm1("atomic_qnorm1");
  CppAD::vector<T> in(1);
  CppAD::vector<T> out(1);
  in[0] = x;
  static_atomic_qnorm1(in, out);
  return out[0];
}
