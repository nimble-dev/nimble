#include <nimble/nimDerivs_atomic_PDinverse_logdet.h>
#include <nimble/nimDerivs_atomic_cache.h>

/*
Atomic class for matrix inverse and log determinant.
See PDinverse_logdet.
This implementation follows the successful pattern of TMB.
The motivating use case is the multivariate normal distribution.

See the matinverse atomic notes as well.

Input X is a positive definite (PD) matrix.
Output Y is a vector of precision elements followed
by a scalar that is the log determinant of the covariance.
These can internally be computed with Cholesky decomposition,
but Cholesky is not placed on the AD tape. Rather, derivatives
for the matrix inverse use the rules for matrix inverse
and derivatives for the log determinant use the rules for
log determinant (Jacobi's formula). The latter uses the 
matrix inverse, giving a second reason (beyond shared use
of Cholesky) to do these steps together.

X could represent the precision or the covariance matrix.
To-Do: Figure out how to handle prec vs. cov.
----
For covariance matrix X:
Value:
Y = c(X^{-1}, log(det(X))) = c(Y1, Y2)
X is n-x-n.
Y is n*n + 1 vector. Y1 are the n*n precision elements, Y2 is log(det(X)).

----
Forward first order
dY1 = -Y1 * dX * Y1 (see atomic_matinverse notes)

// Extend to the case of symmetric matrices represented by triangular parts.

dY2 = trace(Y1 * dX) // Jacobi's formula for log determinant

----
Reverse first order
From Y1:
Xadjoint = -t(Y1) %*% Yadjoint1 %*% t(Y1);


From Y2:
<Yadjoint2, dY2> = <Yadjoint2, trace(Y1 * dX)>
                 = <Yadjoint2, sum_i inprod(Y1[i,], dX[,i])
                 = sum_i <Yadjoint2, inprod(Y1[i,], dX[,i])>
                 = sum_i <Yadjoint2, Y1[i,] * dX[,i]>
                 = sum_i <Y1[i,]^T Yadjoint2, dX[,i]>
                 = sum_i <Xadjoint[, i], dX[ ,i]>
Xadjoint[, i] += Y1[i,]^T Yadjoint2
i.e.
Xadjoint += t(Y1) * Yadjoint2

----
Reverse second order
// I am going to wait on this.
*/

/*
Storage and usage of diagonal half-matrices:
- 
*/

atomic_PDinverse_logdet_class::atomic_PDinverse_logdet_class(const std::string& name) :
  CppAD::atomic_three<double>(name) {};

bool atomic_PDinverse_logdet_class::for_type(
				       const CppAD::vector<double>&               parameter_x ,
				       const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				       CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  // printf("In PDinverse_logdet for_type\n");
  size_t n = type_x.size(); // y is one longer than x
  // All types must be the same.
  CppAD::ad_type_enum final_type(CppAD::constant_enum);
  CppAD::ad_type_enum this_type;
  for(size_t i = 0; i < n; ++i) {
    this_type = type_x[i];
    if(this_type == CppAD::dynamic_enum)
      final_type = CppAD::dynamic_enum;
    if(this_type == CppAD::variable_enum) {
      final_type = CppAD::variable_enum;
      break;
    }
  }
  for(size_t i = 0; i < n+1; ++i) type_y[i] = final_type;
  return true;
}

bool atomic_PDinverse_logdet_class::rev_depend(
					 const CppAD::vector<double>&          parameter_x ,
					 const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
					 CppAD::vector<bool>&                depend_x    ,
					 const CppAD::vector<bool>&          depend_y
					 ) {
  // printf("In PDinverse_logdet reverse_depend\n");
  bool any_depend_y(false);
  size_t ny = depend_y.size();
  for(size_t i = 0; i < ny; ++i) {
    if(depend_y[i]) {
      any_depend_y = true;
      break;
    }
  }
  if(depend_x.size() != ny-1) {
    printf("In PDinverse_logdet rev_depend, somehow size of depend_x (n*n) does not match size of depend_y (n*n+1)).  That should never happen.\n");
  }
  size_t nx = depend_x.size();
  for(size_t i = 0; i < nx; ++i) {
    depend_x[i] = any_depend_y;
  }
  return false;
}

// Forward mode
inline bool atomic_PDinverse_logdet_class::forward(
    const CppAD::vector<double>& parameter_x,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    size_t need_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& taylor_x,
    CppAD::vector<double>& taylor_y
) {
  //forward mode
  // printf("In PDinverse_logdet forward\n");
  size_t nrow = order_up + 1;
  size_t n = static_cast< size_t>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
  EigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  if((order_low <= 0) && (order_up >= 0)) { // value
    EigenMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    
    // Assisted by copilot (begin)
    // (i) Cholesky decomposition (upper-triangular, as in NIMBLE convention)
    auto chol = Xmap.template selfadjointView<Eigen::Upper>().llt();
    //Eigen::MatrixXd chol = llt.matrixU(); // Upper-triangular Cholesky factor

    // (ii) Matrix inverse from the Cholesky
    // Solve U^T U = X, so X^{-1} = U^{-1} (U^T)^{-1}
    Ymap = chol.solve(Eigen::MatrixXd::Identity(n, n)).template triangularView<Eigen::Upper>();

    // (iii) Log determinant from the Cholesky
    // log(det(X)) = 2 * sum(log(diag(U)))
    double logdet = 0.0;
    for(int i = 0; i < n; ++i)
        logdet += 2.0 * std::log(chol.matrixU()(i, i));
    // Assisted by copilot (end)

   taylor_y[n * n * nrow + 0] = logdet; // Last element is log(det(X))
    double_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
  }
 if((order_low <= 1) && (order_up >= 1)) {
    // printf("In forward >1\n");
    double_cache.check_and_set_cache(this,
				     parameter_x,
				     type_x,
				     0,
				     order_up,
				     taylor_x,
				     taylor_y.size());
    size_t cache_nrow = double_cache.nrow();
    EigenMap Ymap(double_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
    EigenMap dYmap(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    EigenConstMap dXmap(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow));
    // multiplying two selfadjointViews does not seem to work.
    // So, I will materialize at least one of them. If that works, I'll make it a member variable.
    MatrixXd Ysym = Ymap.template selfadjointView<Eigen::Upper>();
    MatrixXd dXsym = dXmap.template selfadjointView<Eigen::Upper>();
//    dYmap = -(Ymap.template selfadjointView<Eigen::Upper>() * dXmap.template selfadjointView<Eigen::Upper>() * Ymap.template selfadjointView<Eigen::Upper>()).template triangularView<Eigen::Upper>();
    dYmap = -(Ysym * dXsym * Ysym);
    dYmap = dYmap.template triangularView<Eigen::Upper>();

// std::cout<<"dXmap "<<dXmap<<std::endl;
// std::cout << "Ymap as a matrix from forward 1:\n";
// // Output dYmap as a matrix for readability
// for(int i = 0; i < dYmap.rows(); ++i) {
//     for(int j = 0; j < dYmap.cols(); ++j) {
//         std::cout << dYmap(i, j) << " ";
//     }
//     std::cout << std::endl;
// }

    // Assisted by copilot (begin)
    // Jacobi's formula for log determinant: trace(Y1 * dX)
    // Since trace(A*B) = sum_ij A_ij * B_ji, so sum of element-wise product of Ymap and dXmap.transpose()
    taylor_y[n * n * nrow + 1] = (Ysym.array() * dXsym.transpose().array()).sum();
    // Assisted by copilot (end)
 }
    return true;
}

inline bool atomic_PDinverse_logdet_class::forward(
    const CppAD::vector<CppAD::AD<double> >& parameter_x,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    size_t need_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<CppAD::AD<double> >& taylor_x,
    CppAD::vector<CppAD::AD<double> >& taylor_y
) {
  //forward mode
  // printf("In PDinverse_logdet forward\n");
  size_t nrow = order_up + 1;
  size_t n = static_cast< size_t>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
  metaEigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  if((order_low <= 0) && (order_up >= 0)) { // value
  // metaEigenMap Ymap(&taylor_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    NimArr<2, CppAD::AD<double> > NimArrX;
    NimArr<1, CppAD::AD<double> > NimArrY;
    NimArrX.setSize(n, n);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        NimArrX(i, j) = Xmap(i, j);
    NimArrY = nimDerivs_PDinverse_logdet(NimArrX);
   for (size_t i = 0; i < n*n+1; ++i)
      taylor_y[i*nrow] = NimArrY[i];
    CppADdouble_cache.set_cache( 0, 0, order_up, taylor_x, taylor_y );
  }
 if((order_low <= 1) && (order_up >= 1)) {
    // printf("In meta-forward >1\n");
    CppADdouble_cache.check_and_set_cache(this,
				     parameter_x,
				     type_x,
				     0,
				     order_up,
				     taylor_x,
				     taylor_y.size());
    size_t cache_nrow = CppADdouble_cache.nrow();
    metaEigenMap Ymap(CppADdouble_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
    metaEigenMap dYmap(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    metaEigenConstMap dXmap(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow));
    MatrixXd_CppAD Ysym = Ymap.template selfadjointView<Eigen::Upper>();
    MatrixXd_CppAD dXsym = dXmap.template selfadjointView<Eigen::Upper>();
//    dYmap = -(Ysym * dXsym * Ysym);
    dYmap = nimDerivs_matmult(-Ysym, nimDerivs_matmult( dXsym,  Ysym ) );
    dYmap = dYmap.template triangularView<Eigen::Upper>();
    // Assisted by copilot (begin)
    // Jacobi's formula for log determinant: trace(Y1 * dX)
    // Since trace(A*B) = sum_ij A_ij * B_ji, so sum of element-wise product of Ymap and dXmap.transpose()
    taylor_y[n * n * nrow + 1] = (Ysym.array() * dXsym.transpose().array()).sum();
    // Assisted by copilot (end)
 }
    return true;
}

// Reverse mode
inline bool atomic_PDinverse_logdet_class::reverse(
    const CppAD::vector<double>& parameter_x,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    size_t order_up,
    const CppAD::vector<double>& taylor_x,
    const CppAD::vector<double>& taylor_y,
    CppAD::vector<double>& partial_x,
    const CppAD::vector<double>& partial_y
) {
    size_t nrow = order_up + 1;
    size_t n = static_cast<size_t>(sqrt(static_cast<double>(taylor_x.size()/nrow)));
    
    // if(order_up >= 1) {
    //   EigenConstMap dYmap(&taylor_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    //   std::cout << "Ymap as a matrix starting reverse 2:\n";
    //   // Output dYmap as a matrix for readability
    //   for(int i = 0; i < dYmap.rows(); ++i) {
    //     for(int j = 0; j < dYmap.cols(); ++j) {
    //       std::cout << dYmap(i, j) << " ";
    //     }
    //     std::cout << std::endl;
    //   }
    // }

    double_cache.check_and_set_cache(this,
        parameter_x,
        type_x,
        0, // only use cached values up to order 0
        order_up,
        taylor_x,
        taylor_y.size());

        size_t cache_nrow = double_cache.nrow();  
        EigenConstMap Ymap(double_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
        MatrixXd Ysym = Ymap.template selfadjointView<Eigen::Upper>();
        EigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
        MatrixXd Xsym = Xmap.template selfadjointView<Eigen::Upper>();
        EigenMap Xadjoint_map(&partial_x[0], n, n, EigStrDyn(nrow*n, nrow) );
        if(order_up >= 0) {
          EigenConstMap Yadjoint_map(&partial_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
          Xadjoint_map = (-Ysym.transpose() * Yadjoint_map.template triangularView<Eigen::Upper>() * Ysym.transpose());
          for(size_t i = 0; i < n; ++i) {
            Xadjoint_map(i, i) += Ysym(i, i) * partial_y[n*n*nrow+0]; // Add the contribution from logdet
            for(size_t j = i+1; j < n; ++j) {
              Xadjoint_map(i, j) += Xadjoint_map(j, i) + // complete the contribution from inverse
              2 * Ysym(i, j) * partial_y[n*n*nrow+0]; // Add the contribution from logdet
              Xadjoint_map(j, i) = 0; 
            }
          }
          //            Xadjoint_map += partial_y[n*n*nrow+0] * Ysym.transpose(); // Add the contribution from logdet
        }
        // TO-DO: Make some of the Eigen matrix objects member variables.
        if(order_up >= 1) {
          EigenConstMap Xdot_map(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow ) ); // K-dot
          MatrixXd Xdot_sym = Xdot_map.template selfadjointView<Eigen::Upper>(); //Sigma-dot
          
          EigenConstMap Ydot_adjoint_map(&partial_y[1], n, n, EigStrDyn(nrow*n, nrow ) ); // J-dot adjoint
          MatrixXd Ydot_adjoint_sym = Ydot_adjoint_map.template selfadjointView<Eigen::Upper>(); // P-dot adjoint
          
          EigenMap Xdot_adjoint_map(&partial_x[1], n, n, EigStrDyn(nrow*n, nrow) ); // K-dot adjoint
          
          //Eigen::MatrixXd Y_Xdot_Y_transpose = (Ysym * Xdot_sym * Ysym).transpose();
          //Eigen::MatrixXd Ydot_adjoint_Ytranspose = Ydot_adjoint_sym * Ysym.transpose();
          Eigen::MatrixXd Y_Jdotadjoint_Y = Ysym * Ydot_adjoint_map.template triangularView<Eigen::Upper>() * Ysym;
          Eigen::MatrixXd Y_Xdot = Ysym * Xdot_sym;
          Eigen::MatrixXd new_Xadjoint_terms = (Y_Jdotadjoint_Y * Y_Xdot.transpose() + 
          Y_Xdot * Y_Jdotadjoint_Y);
          Eigen::MatrixXd new_Xadjoint_terms_from_logdet = -Y_Xdot * Ysym;
          for(size_t i = 0; i < n; ++i) {
            Xadjoint_map(i, i) += new_Xadjoint_terms(i, i); // from inverse to diag adj of X
            Xdot_adjoint_map(i, i) = -Y_Jdotadjoint_Y(i, i); // from inverse to diag adj of Xdot
            Xadjoint_map(i, i) += new_Xadjoint_terms_from_logdet(i, i) * partial_y[n*n*nrow+1]; // from logdet to diag adj of X
            Xdot_adjoint_map(i, i) += Ysym(i, i) * partial_y[n*n*nrow+1]; // from logdet to diag adj of Xdot

            for(size_t j = i+1; j < n; ++j) {
              Xadjoint_map(i, j) += new_Xadjoint_terms(i, j) + new_Xadjoint_terms(j, i); // from inverse to adj of X
              Xdot_adjoint_map(i, j) = -Y_Jdotadjoint_Y(i, j)-Y_Jdotadjoint_Y(j, i); // from inverse to adj of Xdot
              Xadjoint_map(i, j) += (new_Xadjoint_terms_from_logdet(i,j) +
                                     new_Xadjoint_terms_from_logdet(j,i)) * partial_y[n*n*nrow+1]; // from logdet to adj of X
              Xdot_adjoint_map(i, j) += 2* Ysym(i, j) * partial_y[n*n*nrow+1]; // from logdet to adj of Xdot
              Xdot_adjoint_map(j, i) = 0;
            }
          }
          //            Xdot_adjoint_map = -(Y_Jdotadjoint_Y.template triangularView<Eigen::Upper>());
        }
//        Xadjoint_map = Xadjoint_map.template triangularView<Eigen::Upper>();
    return true;
}

inline bool atomic_PDinverse_logdet_class::reverse(
				      const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      size_t                              order_up    ,
				      const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
				      const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
				      CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
				      const CppAD::vector<CppAD::AD<double> >&               partial_y   )
{
  size_t nrow = order_up + 1;
  size_t n = static_cast<size_t>(sqrt(static_cast<double>(taylor_x.size()/nrow)));

  CppADdouble_cache.check_and_set_cache(this,
					parameter_x,
					type_x,
					0, // only use cached values up to order 0
					order_up,
					taylor_x,
					taylor_y.size());
  size_t cache_nrow = CppADdouble_cache.nrow();
  metaEigenConstMap Ymap(CppADdouble_cache.taylor_y_ptr(), n, n, EigStrDyn(cache_nrow*n, cache_nrow ) );
  MatrixXd_CppAD Ysym = Ymap.template selfadjointView<Eigen::Upper>();
  metaEigenConstMap Xmap(&taylor_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  MatrixXd_CppAD Xsym = Xmap.template selfadjointView<Eigen::Upper>();
  metaEigenMap Xadjoint_map(&partial_x[0], n, n, EigStrDyn(nrow*n, nrow) );
  // QUESTION: Save recorded operations by hand-coding due to the triangular stuff.
  if(order_up >= 0) {
    metaEigenConstMap Yadjoint_map(&partial_y[0], n, n, EigStrDyn(nrow*n, nrow ) );
    //Xadjoint_map = (-Ysym.transpose() * Yadjoint_map.template triangularView<Eigen::Upper>() * Ysym.transpose());
    Xadjoint_map = nimDerivs_matmult(-Ysym.transpose(),
            nimDerivs_matmult(Yadjoint_map.template triangularView<Eigen::Upper>(),
            Ysym.transpose()));
    for(size_t i = 0; i < n; ++i) {
      Xadjoint_map(i, i) += Ysym(i, i) * partial_y[n*n*nrow+0]; // Add the contribution from logdet
      for(size_t j = i+1; j < n; ++j) {
        Xadjoint_map(i, j) += Xadjoint_map(j, i) + // complete the contribution from inverse
        2 * Ysym(i, j) * partial_y[n*n*nrow+0]; // Add the contribution from logdet
        Xadjoint_map(j, i) = 0; 
      }
    }
  }

  if(order_up >= 1) {
    metaEigenConstMap Xdot_map(&taylor_x[1], n, n, EigStrDyn(nrow*n, nrow ) );
    MatrixXd_CppAD Xdot_sym = Xdot_map.template selfadjointView<Eigen::Upper>(); //Sigma-dot  
    metaEigenConstMap Ydot_adjoint_map(&partial_y[1], n, n, EigStrDyn(nrow*n, nrow ) );
    MatrixXd_CppAD Ydot_adjoint_sym = Ydot_adjoint_map.template selfadjointView<Eigen::Upper>(); // P-dot adjoint
    metaEigenMap Xdot_adjoint_map(&partial_x[1], n, n, EigStrDyn(nrow*n, nrow) );
    //metaEigenMatrixXd Y_Xdot_Y_transpose = nimDerivs_matmult(Ymap, nimDerivs_matmult( Xdot_map, Ymap)).transpose();
    //metaEigenMatrixXd Ydot_adjoint_Ytranspose = nimDerivs_matmult(Ydot_adjoint_map, Ymap.transpose());
    MatrixXd_CppAD Y_Jdotadjoint_Y = nimDerivs_matmult(Ysym,
          nimDerivs_matmult(Ydot_adjoint_map.template triangularView<Eigen::Upper>(), Ysym));
    MatrixXd_CppAD Y_Xdot = nimDerivs_matmult(Ysym, Xdot_sym);
    MatrixXd_CppAD new_Xadjoint_terms = 
       nimDerivs_matmult(Y_Jdotadjoint_Y, Y_Xdot.transpose()) + 
       nimDerivs_matmult(Y_Xdot, Y_Jdotadjoint_Y);
    MatrixXd_CppAD new_Xadjoint_terms_from_logdet = nimDerivs_matmult(-Y_Xdot, Ysym);

    for(size_t i = 0; i < n; ++i) {
      Xadjoint_map(i, i) += new_Xadjoint_terms(i, i);
      Xdot_adjoint_map(i, i) = -Y_Jdotadjoint_Y(i, i);

      Xadjoint_map(i, i) += new_Xadjoint_terms_from_logdet(i, i) * partial_y[n*n*nrow+1]; // from logdet to diag adj of X

      Xdot_adjoint_map(i, i) += Ysym(i, i) * partial_y[n*n*nrow+1];
      for(size_t j = i+1; j < n; ++j) {
        Xadjoint_map(i, j) += new_Xadjoint_terms(i, j) + new_Xadjoint_terms(j, i);
        Xdot_adjoint_map(i, j) = -Y_Jdotadjoint_Y(i, j)- Y_Jdotadjoint_Y(j, i);
        Xadjoint_map(i, j) += (new_Xadjoint_terms_from_logdet(i,j) +
               new_Xadjoint_terms_from_logdet(j,i)) * partial_y[n*n*nrow+1]; // from logdet to adj of X
        Xdot_adjoint_map(i, j) += 2* Ysym(i, j) * partial_y[n*n*nrow+1];
        Xdot_adjoint_map(j, i) = 0;
      }
    }
  }
  return true;
}

/* void PDinverse_logdet(const MatrixXd_CppAD &x, // This (non-template) type forces any incoming expression to be evaluated
			 MatrixXd_CppAD &y) {
  // static PDinverse_logdet_class PDinverse_logdet("PDinverse_logdet"); // this has no state information so the same object can be used for all cases
  atomic_PDinverse_logdet_class *PDinverse_logdet;
  size_t n = x.rows();
  std::vector<CppAD::AD<double> > xVec(n*n);
  mat2vec(x, xVec);
  std::vector<CppAD::AD<double> > yVec(n*n+1);
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    PDinverse_logdet = new atomic_PDinverse_logdet_class("PDinverse_logdet");
  } else {
    void *tape_mgr = CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr();
    PDinverse_logdet = new_atomic_PDinverse_logdet(tape_mgr, "PDinverse_logdet");
  }
  (*PDinverse_logdet)(xVec, yVec);
  y.resize(n*n+1, 1);
  for(size_t i = 0; i < n*n+1; ++i) {
    y(i) = yVec[i];
  }
  if(!recording) {
    delete PDinverse_logdet;
  } else {
    track_nimble_atomic(PDinverse_logdet,
			CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
			CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
} */

/*
MatrixXd_CppAD nimDerivs_PDinverse_logdet(const MatrixXd_CppAD &x, const CppAD::AD<double> &prec_param) {
  // prec_param is ignored for now.
  MatrixXd_CppAD ans;
  PDinverse_logdet(x, ans);
  return ans;
}
*/

NimArr<1, CppAD::AD<double> > nimDerivs_PDinverse_logdet(const NimArr<2, CppAD::AD<double> > &x) {
  // prec_param was being ignored and has been removed for now.
  atomic_PDinverse_logdet_class *PDinverse_logdet;
  size_t n = x.dim()[0];
  std::vector<CppAD::AD<double> > xVec(n*n);
  for(size_t i = 0; i < n; ++i) {
    for(size_t j = 0; j < n; ++j) {
      xVec[j*n + i] = x(i, j);
    }
  }
  std::vector<CppAD::AD<double> > yVec(n*n+1);
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    PDinverse_logdet = new atomic_PDinverse_logdet_class("PDinverse_logdet");
  } else {
    void *tape_mgr = CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr();
    PDinverse_logdet = new_atomic_PDinverse_logdet(tape_mgr, "PDinverse_logdet");
  }
  (*PDinverse_logdet)(xVec, yVec);
  NimArr<1, CppAD::AD<double> > y;
  y.initialize(0, false, n*n+1);
  for(size_t i = 0; i < n*n+1; ++i) {
    y(i) = yVec[i];
  }
  if(!recording) {
    delete PDinverse_logdet;
  } else {
    track_nimble_atomic(PDinverse_logdet,
			CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
			CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
  return y;
}
