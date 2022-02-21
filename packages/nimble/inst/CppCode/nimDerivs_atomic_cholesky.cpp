#include <nimble/nimDerivs_atomic_cholesky.h>
#include <nimble/nimDerivs_atomic_cache.h>
#include <nimble/nimbleCppAD.h> // for ad_timer only

template<typename T>
void HalfDiag(T &Mat) { // There may be a pure Eigen way to do this.
  // This is the upper version of the Phi function of Iain Murray's archive paper
  size_t n = Mat.rows();
  for(size_t i = 0; i < n; ++i) {
    Mat(i, i) *= 0.5;
  }
}

//#define _TIME_AD_CHOL // also need _TIME_AD in nimbleCppAD.h
#ifdef _TIME_AD_CHOL
ad_timer derivs_chol_timer("derivs_chol");
SEXP report_AD_timers() {
  derivs_chol_timer.show_report();
  return(R_NilValue);
}
SEXP reset_AD_timers(SEXP SreportInterval) {
  derivs_chol_timer.reset();
  derivs_chol_timer.set_interval(INTEGER(SreportInterval)[0]);
  return(R_NilValue);
}
void derivs_chol_timer_start() {derivs_chol_timer.start(false);}
void derivs_chol_timer_stop() {derivs_chol_timer.stop(false);}

#endif

/*
We follow Murray (2018), but he uses lower-triangular whereas
we follow R in using upper triangular.  So we adjust Murray's results 
accordingly.

----
Value:
U = Chol(Sigma)
Sigma = L L^T = U^T U.  Here U is the "Y", defined by: U = Chol(Sigma)
-----
Forward 1:
dSigma = dU^T U + U^T dU
U^-T dSigma U^-1 = U^-T dU-T + dU U^-1
This is lower*lower + upper*upper, and the result must be symmetric, so we can eliminate one.
Define Phi (Murray's notation, changed for upper triangular)
  to be the upper_HalfDiag function above.  Phi(A) = upper triangular of A with the diagonal multiplied by 0.5.
Phi(U^-T dSigma U^-1) = dU U^-1
dU = Phi(U^-T dSigma U^-1) * U
Z = U^-T A = U^T.solve(A)
Z = A U^-1: Z^T = U^-T A^T. So Z = U^T.solve(A^T)^T
-----
Reverse 1: (Here U is "Y" and Sigma is "X")
<Yadjoint, dY> = <Yadjoint, Phi(U^-T dSigma U^-1) * U>, but we need to arrange dSigma out of this,
<Yadjoint U^T, Phi(U^-T dSigma U^-1)>
< Phi(Yadjoint U^T), U^-T dSigma U^-1 >
< U^-1 Phi(Yadjoint U^T) U^-T, dSigma >
This sort-of identifies Sigma_adjoint as U^-1 Phi(Yadjoint U^T) U^-T.  However, there is the issue of symmetry and independent elements.
This says that a "dG" is the sum of the component-wise product of (U^-1 Phi(Yadjoint U^T) U^-T) and dSigma.
Now we need that in terms of dK, the upper triangular of dSigma.
We write dSigma = dK + dK^T - diag(dK),   (following Murray here)
In other words, we view K as the "real" input, and Sigma as a function of K. Sigma = K + K^T - diag(K).
We seek
< Kadjoint, dK > = < U^-1 Phi(Yadjoint U^T) U^-T, dSigma >
   = < U^-1 Phi(Yadjoint U^T) U^-T, dK + dK^T - diag(dK) >
   = < U^-1 Phi(Yadjoint U^T) U^-T, dK > + < (U^-1 Phi(Yadjoint U^T) U^-T)^T, dK > - < diag(U^-1 Phi(Yadjoint U^T) U^-T), dK>

Giving Kadjoint = U^-1 Phi(Yadjoint U^T) U^-T + (U^-1 Phi(Yadjoint U^T) U^-T)^T - diag(U^-1 Phi(Yadjoint U^T) U^-T)
Z = A U^-T: Z^T = U^-1 A^T = U.solve(A^T).  So Z = U.solve(A^T)^T
-----
Reverse 2

<Yadjoint, dY> + <Ydot_adjoint, dYdot> = terms above + new terms

new terms: <Ydot_adjoint, dYdot>
To obtain dYdot:
Sigma_dot = Udot^T U + U^T Udot


*/

atomic_cholesky_class::atomic_cholesky_class(const std::string& name) : CppAD::atomic_three<double>(name) {};

// Output inherits type from inputs above it and/or to the left of it.
// i.e. input x(i, j) impacts output y(>=i, >=j).
// Although we are doing upper-triangular cholesky, Sigma = U^T U (X = Y^T T),
// I find it easier to think about Sigma = L L^T and then U = L^T.
bool atomic_cholesky_class::for_type(
				     const CppAD::vector<double>&               parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  int nsq = type_y.size();
  int n = static_cast<size_t>(sqrt(static_cast<double>(nsq)));
  std::vector< CppAD::ad_type_enum > row_types(n, CppAD::constant_enum);
  CppAD::ad_type_enum this_type, this_col_type, this_row_type;
  for(int j = 0; j < n; ++j) {
    this_col_type = CppAD::constant_enum;
    for(int i = 0; i <= j; ++i) {
      this_row_type = row_types[i];
      this_type = type_x[i + j*n];
      // match running column type and this type to larger
      if(this_type > this_col_type) {
	this_col_type = this_type;
      } else {
	this_type = this_col_type;
      }
      // match running row type and this type to larger
      if(this_type > this_row_type) {
	this_row_type = this_type;
      } else {
	this_type = this_row_type;
      }
      this_col_type = this_type; // pick up promotion from row in case it pulled this_type above the running col type
      row_types[i] = this_row_type;
      type_y[i + j*n] = this_type;
    }
    for(int i = j+1; i < n; ++i) { // set strict below-diagonal elements constant
      type_y[i + j*n] = CppAD::constant_enum;
    }
  }
  return true;
}

// rev_depend is used when the optimize() method is called for a tape (an ADFun).
bool atomic_cholesky_class::rev_depend(
				       const CppAD::vector<double>&          parameter_x ,
				       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				       CppAD::vector<bool>&                depend_x    ,
				       const CppAD::vector<bool>&          depend_y
				       ) {
  // We need to set any elements of depend_x to true that are inputs to elements of y with depend_y true
  // x(i,j) impacts on y(>=i, >=j).
  // Therefore y(i,j) is impacted by elements x(<=i, <=j) 
  int nsq = parameter_x.size();
  int n = static_cast<size_t>(sqrt(static_cast<double>(nsq)));
  std::vector<bool> dep_cols(n, false);
  bool this_dep, dep_this_row, dep_this_col;
  for(int i = n-1; i >= 0; --i) {
    dep_this_row = false;
    for(int j = n-1; j >= i; --j) {
      dep_this_col = dep_cols[j];
      this_dep = depend_y[i + j*n];
      if(this_dep > dep_this_row) {
	dep_this_row = this_dep;
      } else {
	this_dep = dep_this_row;
      }
      if(this_dep > dep_this_col) {
	dep_this_col = this_dep;
      } else {
	this_dep = dep_this_col;
      }
      dep_this_row = this_dep;
      dep_cols[j] = this_dep;
      depend_x[i + j * n] = this_dep;
    }
    // When treating input as upper-diagonal
    for( int j = 0; j < i; ++j) {
      depend_x[i + j*n] = false;
    }
  }

  // When treating input as symmetric
  // for(int i = 0; i < n; ++i) {
  //   for( int j = 0; j < i; ++j) {
  //     depend_x[i + j*n] = depend_x[j + i*n];
  //   }
  // }
  return true;
}

bool atomic_cholesky_class::forward(
				    const CppAD::vector<double>&               parameter_x  ,
				    const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				    size_t                              need_y       ,
				    size_t                              order_low    ,
				    size_t                              order_up     ,
				    const CppAD::vector<double>&               taylor_x     ,
				    CppAD::vector<double>&                     taylor_y     ) {
#ifdef _TIME_AD_CHOL
  derivs_chol_timer_start();
#endif
  //forward mode
  //  std::cout<<"In Cholesky forward "<<order_low<<" "<<order_up<<std::endl;
  int nrow = order_up + 1;
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
  // populate cholesky into taylor_y
  if(order_low <= 0 & order_up >= 0) {//value
    // Eigen::MatrixXd xMat(n, n);
    //     uppervec2mat(taylor_x, xMat, n, nrow); // Potential to use an Eigen map if only order is 0
    // Potential to avoid turning chol into a full matrix here, leaving it as an llt object.
    EigenMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    EigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
    Ymap = Xmap.template selfadjointView<Eigen::Upper>().llt().matrixU();
    // std::cout<<"Ymap calculated in Forward 0"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < i; ++j) std::cout<<"\t";
    //   for(int j = i; j < n; ++j) std::cout<<Ymap(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    double_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
  }
  if(order_low <= 1 & order_up >= 1) {
    //  printf("In forward >1\n");
    double_cache.check_and_set_cache(this,
				     parameter_x,
				     type_x,
				     0,
				     order_up,
				     taylor_x,
				     taylor_y.size());
    int cache_nrow = double_cache.nrow();
    EigenMap Ymap(double_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );

    EigenMap dYmap(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    EigenConstMap dXmap(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow));


    // std::cout<<"Ymap"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < i; ++j) std::cout<<"\t";
    //   for(int j = i; j < n; ++j) std::cout<<Ymap(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    
    // std::cout<<"dXmap"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<dXmap(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    // U is Y
    // Sigma is X.
    // Y = chol(X)
    // dX^T should be used as dXmap.selfAdjointView to get the symmetry, but that doesn't seem to be supported as a solve argument.
    Eigen::MatrixXd tempMat = dXmap.template selfadjointView<Eigen::Upper>(); //dSigma

    // std::cout<<"tempMat"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<tempMat(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    Eigen::MatrixXd UinvT_dSigma_Uinv = 
      Ymap.transpose().triangularView<Eigen::Lower>().solve(
							    (Ymap.transpose().triangularView<Eigen::Lower>().solve(tempMat)).eval().transpose()
							    ).eval().triangularView<Eigen::Upper>();
    // The .eval()'s are very defensive.  In some cases with maps on LHS, Eigen does not seem to respect them.
    // std::cout<<"UinvT_dSigma_Uinv"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<UinvT_dSigma_Uinv(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    HalfDiag(UinvT_dSigma_Uinv);
    // std::cout<<"UinvT_dSigma_Uinv"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<UinvT_dSigma_Uinv(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    
    dYmap = UinvT_dSigma_Uinv * Ymap.triangularView<Eigen::Upper>(); // other elements of Ymap might be nans
    // std::cout<<"dYmap"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < i; ++j) std::cout<<"\t";
    //   for(int j = i; j < n; ++j) std::cout<<dYmap(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    double_cache.set_cache( 1, 1, order_up, taylor_x, taylor_y );
  }
  //  printf("done cholesky forward\n");
#ifdef _TIME_AD_CHOL
  derivs_chol_timer_stop();
#endif
  return true;
}

bool atomic_cholesky_class::forward(
				    const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
				    const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				    size_t                              need_y       ,
				    size_t                              order_low    ,
				    size_t                              order_up     ,
				    const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
				    CppAD::vector<CppAD::AD<double> >&                     taylor_y     ) {
  //forward mode
  // printf("In cholesky forward\n");
  int nrow = order_up + 1;
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
  // populate cholesky into taylor_y
  if(order_low <= 0 & order_up >= 0) {//value
    // Eigen::MatrixXd xMat(n, n);
    //     uppervec2mat(taylor_x, xMat, n, nrow); // Potential to use an Eigen map if only order is 0
    // Potential to avoid turning chol into a full matrix here, leaving it as an llt object.
    metaEigenMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    metaEigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
    Ymap = nimDerivs_EIGEN_CHOL(Xmap);
    CppADdouble_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
  }
  if(order_low <= 1 & order_up >= 1) {
    // printf("In forward >1\n");
    CppADdouble_cache.check_and_set_cache(this,
					  parameter_x,
					  type_x,
					  0,
					  order_up,
					  taylor_x,
					  taylor_y.size());
    int cache_nrow = CppADdouble_cache.nrow();
    metaEigenMap Ymap(CppADdouble_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );

    metaEigenMap dYmap(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    metaEigenConstMap dXmap(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow));
      
    // U is Y
    // Sigma is X.
    // Y = chol(X)
    // dX^T should be used as dXmap.selfAdjointView to get the symmetry, but that doesn't seem to be supported as a solve argument.
    MatrixXd_CppAD tempMat = dXmap.template selfadjointView<Eigen::Upper>(); //dSigma
    MatrixXd_CppAD UinvT_dSigma_Uinv =
      nimDerivs_EIGEN_FS( Ymap.transpose(),
			  nimDerivs_EIGEN_FS( Ymap.transpose(), tempMat ).transpose() ).triangularView<Eigen::Upper>();
	
    // Ymap.transpose().triangularView<Eigen::Lower>().solve(
    // 						      (Ymap.transpose().triangularView<Eigen::Lower>().solve(tempMat)).transpose()
    // 						      ).triangularView<Eigen::Upper>();
    HalfDiag(UinvT_dSigma_Uinv);
    dYmap = nimDerivs_matmult(UinvT_dSigma_Uinv, Ymap.triangularView<Eigen::Upper>());
    CppADdouble_cache.set_cache( 1, 1, order_up, taylor_x, taylor_y );

  }
  //  printf("done cholesky forward\n");
  return true;
}

bool atomic_cholesky_class::reverse(
				    const CppAD::vector< double >&               parameter_x ,
				    const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				    size_t                              order_up    ,
				    const CppAD::vector< double >&               taylor_x    ,
				    const CppAD::vector< double >&               taylor_y    ,
				    CppAD::vector< double >&                     partial_x   ,
				    const CppAD::vector< double >&               partial_y   ) {
  //reverse m_ode
  // printf("In reverse\n");
  // std::cout<<"In Chol reverse "<<order_up<<std::endl;
#ifdef _TIME_AD_CHOL
  derivs_chol_timer_start();
#endif

  int nrow = order_up + 1;
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
    
  double_cache.check_and_set_cache(this,
				   parameter_x,
				   type_x,
				   order_up >= 1 ? 1 : 0, // only use cached values up to order 1
				   order_up,
				   taylor_x,
				   taylor_y.size());
  int cache_nrow = double_cache.nrow();
  EigenConstMap Ymap(double_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
  
  EigenConstMap Yadjoint_map(&partial_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
  EigenMap Xadjoint_map(&partial_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  // std::cout<<"Ymap input to reverse 1"<<std::endl;
  // for(int i = 0; i < n; ++i) {
  //   for(int j = 0; j < n; ++j) std::cout<<Ymap(i, j)<<"\t";
  //   std::cout<<std::endl;
  // }

  // std::cout<<"Yadjoint input to reverse 1"<<std::endl;
  // for(int i = 0; i < n; ++i) {
  //   for(int j = 0; j < n; ++j) std::cout<<Yadjoint_map(i, j)<<"\t";
  //   std::cout<<std::endl;
  // }
  if(order_up >= 0) {


    
    // Yadjoint_map is how final quantity changes w.r.t elements of Y.
    // By definition we treat it as upper triangular.
    // Eigen doesn't seem to support matrix multiplication of two triangular views.
    // However, when we use upper_HalfDiag, we remove lower-diagonal elements of the product,
    // which means that lower-diagonal elements of Yadjoint_map are de facto removed.
    //      Eigen::MatrixXd Yadjoint_upper = Yadjoint_map.triangularView<Eigen::Upper>();
    Eigen::MatrixXd Yadjoint_UT = (Yadjoint_map * Ymap.transpose().template triangularView<Eigen::Lower>()).template triangularView<Eigen::Upper>();
    HalfDiag(Yadjoint_UT);
    // std::cout<<"Yadjoint_UT"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<Yadjoint_UT(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    Eigen::MatrixXd Uinv_PhiStuff_UinvT = Ymap.template triangularView<Eigen::Upper>().solve(
											     Ymap.template triangularView<Eigen::Upper>().solve( Yadjoint_UT.transpose() ).eval().transpose()
											     ).eval();
    // Resulting adjoint is how final quantity changes w.r.t. elements of Sigma (i.e. X, i.e. covariance matrix)
    // This should be only upper triangular

    // std::cout<<"YUinv_PhiStuff_UinvT"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<Uinv_PhiStuff_UinvT(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    Xadjoint_map = Uinv_PhiStuff_UinvT.template triangularView<Eigen::Upper>();
    Xadjoint_map += Uinv_PhiStuff_UinvT.transpose().template triangularView<Eigen::StrictlyUpper>();      
    //Xadjoint_map = Uinv_PhiStuff_UinvT +  Uinv_PhiStuff_UinvT.transpose();
    //    Xadjoint_map.diagonal() -= Uinv_PhiStuff_UinvT.diagonal();
    // cut diagonal in half and remove lower diagonal:
    //upper_HalfDiag(Xadjoint_map);

    // std::cout<<"Xadjoint calculated in reverse 1"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<Xadjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
  }
  if(order_up >= 1) {
    //    std::cout<<"entering second order reverse"<<std::endl;
    
    //    EigenConstMap Ydot_map(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    EigenConstMap Ydot_map(double_cache.taylor_y_ptr() + 1, n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
  
    EigenConstMap Ydot_adjoint_map(&partial_y[1], n, n, EigStrDyn(nrow*n, nrow ) );

    // std::cout<<"Ydot"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<Ydot_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    // std::cout<<"Ydot_adjoint"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<Ydot_adjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    Eigen::MatrixXd P = (Ydot_adjoint_map * Ymap.transpose().template triangularView<Eigen::Lower>()).template triangularView<Eigen::Upper>();
    HalfDiag(P);

    // std::cout<<"P"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<<P(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }


    Eigen::MatrixXd Uinv_P_UinvT = Ymap.template triangularView<Eigen::Upper>().solve(
										      Ymap.template triangularView<Eigen::Upper>().solve( P.transpose() ).eval().transpose()
										      ).eval();

    // std::cout<<"Uinv_P_UinvT"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<< Uinv_P_UinvT(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    EigenMap Xdot_adjoint_map(&partial_x[1], n, n, EigStrDyn(nrow*n, nrow) );
    Xdot_adjoint_map = Uinv_P_UinvT.template triangularView<Eigen::Upper>();
    Xdot_adjoint_map += Uinv_P_UinvT.transpose().template triangularView<Eigen::StrictlyUpper>();      
      
    Eigen::MatrixXd R = (Ydot_map * Ymap.template triangularView<Eigen::Upper>().solve( P + P.transpose() ) ).eval().template triangularView<Eigen::Upper>();
    HalfDiag(R);

    // std::cout<<"R"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<< R(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    
    Eigen::MatrixXd Xadjoint_term = Ymap.template triangularView<Eigen::Upper>().solve(
										       Ymap.template triangularView<Eigen::Upper>().solve(R.transpose() ).eval().transpose()
										       ).eval();
    // std::cout<<"Xadjoint_term"<<std::endl;
    // for(int i = 0; i < n; ++i) {
    //   for(int j = 0; j < n; ++j) std::cout<< Xadjoint_term(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    Xadjoint_map -= Xadjoint_term.template triangularView<Eigen::Upper>();
    Xadjoint_map -= Xadjoint_term.transpose().template triangularView<Eigen::StrictlyUpper>(); 
  }
#ifdef _TIME_AD_CHOL
  derivs_chol_timer_start();
#endif
  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

bool atomic_cholesky_class::reverse(
				    const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				    const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				    size_t                              order_up    ,
				    const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
				    const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
				    CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
				    const CppAD::vector<CppAD::AD<double> >&               partial_y   ) {
  int nrow = order_up + 1;
  int n = static_cast<int>(sqrt(static_cast<double>(taylor_x.size()/nrow)));


  CppADdouble_cache.check_and_set_cache(this,
					parameter_x,
					type_x,
					order_up >= 1 ? 1 : 0, // only use cached values up to order 1
					order_up,
					taylor_x,
					taylor_y.size());
  //  metaEigenConstMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
  int cache_nrow = CppADdouble_cache.nrow();
  metaEigenConstMap Ymap(CppADdouble_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );

  metaEigenConstMap Yadjoint_map(&partial_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
  metaEigenMap Xadjoint_map(&partial_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  if(order_up >= 0) {
    MatrixXd_CppAD Yadjoint_UT =  Ymap.transpose().template triangularView<Eigen::Upper>();//nimDerivs_matmult(Yadjoint_map, Ymap.transpose().template triangularView<Eigen::Lower>()).template triangularView<Eigen::Upper>();
    //      MatrixXd_CppAD Yadjoint_UT = (Yadjoint_map * Ymap.transpose().template triangularView<Eigen::Lower>()).template triangularView<Eigen::Upper>();
    HalfDiag(Yadjoint_UT);
    MatrixXd_CppAD Uinv_PhiStuff_UinvT = nimDerivs_EIGEN_BS(Ymap,
							    nimDerivs_EIGEN_BS(Ymap, Yadjoint_UT.transpose() ).transpose()
							    );
    // MatrixXd_CppAD Uinv_PhiStuff_UinvT = Ymap.template triangularView<Eigen::Upper>().solve(
    // 										      Ymap.template triangularView<Eigen::Upper>().solve( Yadjoint_UT.transpose() ).transpose()
    // 										      );
    Xadjoint_map = Uinv_PhiStuff_UinvT.template triangularView<Eigen::Upper>();
    Xadjoint_map += Uinv_PhiStuff_UinvT.transpose().template triangularView<Eigen::StrictlyUpper>();      
  }
  if(order_up >= 1) {
    //    metaEigenConstMap Ydot_map(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    metaEigenConstMap Ydot_map(CppADdouble_cache.taylor_y_ptr() + 1, n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );

    metaEigenConstMap Ydot_adjoint_map(&partial_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    MatrixXd_CppAD P = nimDerivs_matmult(Ydot_adjoint_map, Ymap.transpose().template triangularView<Eigen::Lower>()).template triangularView<Eigen::Upper>();
    //MatrixXd_CppAD P = (Ydot_adjoint_map * Ymap.transpose().template triangularView<Eigen::Lower>()).template triangularView<Eigen::Upper>();
    HalfDiag(P);
    MatrixXd_CppAD Uinv_P_UinvT = nimDerivs_EIGEN_BS(Ymap,
						     nimDerivs_EIGEN_BS(Ymap, P.transpose() ).transpose()
						     );
    // MatrixXd_CppAD Uinv_P_UinvT = Ymap.template triangularView<Eigen::Upper>().solve(
    // 										       Ymap.template triangularView<Eigen::Upper>().solve( P.transpose() ).transpose()
    // 										       );
    metaEigenMap Xdot_adjoint_map(&partial_x[1], n, n, EigStrDyn(nrow*n, nrow) );
    Xdot_adjoint_map = Uinv_P_UinvT.template triangularView<Eigen::Upper>();
    Xdot_adjoint_map += Uinv_P_UinvT.transpose().template triangularView<Eigen::StrictlyUpper>();      
      
    MatrixXd_CppAD R = nimDerivs_matmult(Ydot_map,
					 nimDerivs_EIGEN_BS(Ymap, P + P.transpose() ) ).template triangularView<Eigen::Upper>();
    //      MatrixXd_CppAD R = (Ydot_map * Ymap.template triangularView<Eigen::Upper>().solve( P + P.transpose() ) ).template triangularView<Eigen::Upper>();
    HalfDiag(R);
    MatrixXd_CppAD Xadjoint_term = nimDerivs_EIGEN_BS(Ymap,
						      nimDerivs_EIGEN_BS(Ymap, R.transpose() ).transpose()
						      );
    // MatrixXd_CppAD Xadjoint_term = Ymap.template triangularView<Eigen::Upper>().solve(
    // 								Ymap.template triangularView<Eigen::Upper>().solve(R.transpose() ).transpose()
    // 								);
    Xadjoint_map -= Xadjoint_term.template triangularView<Eigen::Upper>();
    Xadjoint_map -= Xadjoint_term.transpose().template triangularView<Eigen::StrictlyUpper>(); 
  }
  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

void atomic_cholesky(const MatrixXd_CppAD &x, // This (non-template) type forces any incoming expression to be evaluated
		     MatrixXd_CppAD &y) {
  //  static atomic_cholesky_class atomic_cholesky("atomic_cholesky"); // this has no state information so the same object can be used for all cases
  atomic_cholesky_class *atomic_cholesky; // Need to do it this way for multiple compilation units
  int n = x.rows();
  CppAD::vector<CppAD::AD<double> > xVec(n*n);
  mat2vec(x, xVec); // could be mat2vec_lower_zero but it doesn't seem to matter.
  CppAD::vector<CppAD::AD<double> > yVec(n*n);
  void *tape_mgr = CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr();
  atomic_cholesky = new_atomic_cholesky(tape_mgr, "atomic_cholesky");
  (*atomic_cholesky)(xVec, yVec);
  y.resize(n, n);
  vec2mat(yVec, y);
  if(CppAD::AD<double>::get_tape_handle_nimble() == nullptr) {
    delete_atomic_cholesky(tape_mgr, atomic_cholesky);
  } else {
    track_nimble_atomic(atomic_cholesky,
			CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
			CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
}

MatrixXd_CppAD nimDerivs_EIGEN_CHOL(const MatrixXd_CppAD &x) {
  MatrixXd_CppAD ans;
  atomic_cholesky(x, ans);
  return ans;
}
