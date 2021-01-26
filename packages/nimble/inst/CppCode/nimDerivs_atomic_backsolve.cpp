#include <nimble/nimDerivs_atomic_backsolve.h>
#include <nimble/nimDerivs_atomic_forwardsolve.h>

/* 
Atomic class for backsolve
We have Y = backsolve(A, B), defined by Y = A^-1 %*% B for A upper triangular. i.e. AY = B
A is n1-x-n1
B is n1-x-n2
Y is n1-x-n2

The input X is c(A, B)
In CppAD notation:
  input length is n = n1*n1 + n1*n2
  output length is m = n1*n2

From Giles(2008)
dY = A^{-1} %*% (dB - dA %*% Y) = backsolve(A, dB - dA %*% Y) [in R-ish notation]
adjoint(B) = A^{-T} %*% adjoint(Y) = forwardsolve(t(A), adjoint(Y)) = backsolve(A, ajoint(Y), transpose = TRUE) [in R-ish notation)
adjoint(A) = -A^{-T} %*% adjoint(Y) %*% Y^T = -adjoint(B) %*% Y^T

In my notation:
-----
Value:
Y = solve(A, B) = A^{-1} * B, A upper triangular
-----
First order:
A and B are both parts of "X"
A * Y = B
dA * Y + A dY = dB
dY = A^-1 * (dB - dA * Y) = solve(A, dB - dA * Y)
------
Reverse first order:
<Yadjoint, dY> = <Aadjoint, dA> + <Badjoint, dB>
               = <Yadjoint, - A^-1 * dA * Y>  + <Yadjoint, A^-1 * dB>
               = <- A^-T * Yadjoint * Y^T, dA> + <A^-T * Yadjoint, dB > 
This gives (matching Giles)
Aadjoint = - A^-T * Yadjoint * Y^T = -solve(A^T,  Yadjoint * Y^T)
Badjoint = A^-T * Yadjoint = solve(A^T, Yadjoint)
Note that Aadjoint = -Badjoint * Y^T
And actually the strict lower-diagonal of A is 0 by definition, so really
Aadjoint = Upper(-Badjoint * Y^T)
------
Reverse second order:
Fdot = A^-1 * Bdot - A^-1 * Adot * Y
dFdot = (dA^-1 * Bdot) + (A^-1 dBdot) - (dA^-1 * Adot * Y) - (A^-1 * dAdot * Y) - (A^-1 * Adot * dY)
Is it better to do it this way:
A * Fdot = Bdot - Adot * Y
dA * Fdot + A * dFdot = dBdot - dAdot * Y - Adot * dY
dFdot = (A^-1 dBdot) - (A^-1 dAdot Y) - (A^-1 Adot dY) - (A^-1 dA Fdot)
      = (matches 2nd) - (matches 4th) - (matches 5th)  - (not clear)
I guess we use, for A = F(A) = A^-1 that dA^-1 = -A^-1 * dA * A^-1
So everything matches *if*
  (dA^-1 * Bdot) - (dA^-1 * Adot * Y) =  -(A^-1 dA Fdot)
  (dA^-1 * Bdot) - (dA^-1 * Adot * Y) =  - A^-1 dA (A^-1 * Bdot - A^-1 * Adot * Y)
                                      =  - dA^-1 *Bdot - dA^-1 * Adot * Y
  So they match.  Better to use the second version instead of forming dA^-1

  The dY term needs to be expanded
  - (A^-1 Adot dY) = - A^-1 Adot A^-1 * (dB - dA * Y) = -(A^-1 Adot A^-1 * dB) + (A^-1 Adot A^-1 dA Y)

<Yadjoint, dY> + <Ydot_adjoint, dYdot> = <reverse first order> + <Ydot_adjoint,  (A^-1 dBdot) - (A^-1 dAdot Y) -(A^-1 Adot A^-1 * dB) + (A^-1 Adot A^-1 dA Y) - (A^-1 dA Fdot) >
<Ydot_adjoint, dYdot> = <A^-T Ydot_adjoint, dBdot> +
                       < -A^-T  Ydot_adjoint Y^T, dAdot > + 
                       < -(A^-1 Adot A^-1)^T  Ydot_adjoint , dB> + 
                       < (A^-1 Adot A^-1)^T Ydot_adjoint Y^T, dA> +
                       < -A^-T Ydot_adjoint Fdot^T, dA>

Aadjoint = -(A^-T * Yadjoint * Y^T) + ((A^-1 Adot A^-1)^T Ydot_adjoint Y^T) - (A^-T Ydot_adjoint Ydot^T)
Badjoint = A^-T * Yadjoint - (A^-1 Adot A^-1)^T  Ydot_adjoint
Note
Aadjoint = -Badjoint * Y^T  - (A^-T Ydot_adjoint Ydot^T)

Adot_adjoint = -A^-T  Ydot_adjoint Y^T
Bdot_adjoint = A^-T Ydot_adjoint
Note
Adot_adjoint = -Bdot_adjoint Y^T

A^-1 Z A^-1 = solve(A, W) where
W = Z A^-1 or WA = Z or A^T W^T = Z^T, so W^T = solve(A^T, Z^T) or W = solve(A^T, Z^T)^T
So A^-1 Z A^-1 = solve(A, solve(A^T, Z^T)^T)
*/

atomic_backsolve_class::atomic_backsolve_class(const std::string& name) : CppAD::atomic_three<double>(name) {};

bool atomic_backsolve_class::for_type(
				      const CppAD::vector<double>&               parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
  //printf("In backsolve for_type\n");
  // Row i of output depends on rows i and >i in A input.
  // Element (i,j) of output depends on elements (>=i, j) of B input
  int n = type_x.size();
  int m = type_y.size();
  int n1sq = n-m;
  int n1 = sqrt( static_cast<double>(n1sq) );
  int n2 = m/n1;
  //    std::cout<<"n = "<<n<<" m = "<<m<<" n1 = "<<n1<<" n2 = "<<n2<<std::endl;

  CppAD::vector<CppAD::ad_type_enum> x1RowTypes(n1); //A
  CppAD::vector<CppAD::ad_type_enum> x2ColTypes(n2); //B

  CppAD::ad_type_enum row_type_below(CppAD::constant_enum);
  CppAD::ad_type_enum this_row_type;
  CppAD::ad_type_enum item_type;
  for(int i = n1-1; i >= 0; --i) { // go backwards because rows depend on those below them.  Must use int to get valid >= 0 comparison.
    this_row_type = row_type_below;
    if(this_row_type != CppAD::variable_enum) { // no need to check if row type is already at the "max"
      for(size_t j = i; j < n1; ++j) { // only look at upper triangle
	item_type = type_x[i + j*n1];
	if(item_type == CppAD::variable_enum) {
	  this_row_type = CppAD::variable_enum;
	  break;
	} else {
	  if(item_type == CppAD::dynamic_enum)
	    this_row_type = CppAD::dynamic_enum;
	}
      }
    }
    x1RowTypes[i] = row_type_below = this_row_type;
  }

  CppAD::ad_type_enum B_item_type;
  CppAD::ad_type_enum B_row_type_below;  // type of elements below row i, done for one column a time
  CppAD::ad_type_enum  B_this_row_type;  // type of the current row, done for one column at a time
  for(int j = 0; j < n2; ++j) {
    // for the j-th column
    B_row_type_below = CppAD::constant_enum;
    for(int i = n1-1; i >= 0; --i) { // must be int, not size_t
      B_this_row_type = B_row_type_below; // can't be lower type than what's below it in the current column
      // if B_this_row_type is not already variable (max type), update based on type of element B[i,j]
      if(B_this_row_type < CppAD::variable_enum) {
	B_item_type = type_x[n1sq + i + j*n1];
	if(B_item_type == CppAD::variable_enum) {
	  B_this_row_type = CppAD::variable_enum;
	} else {
	  if(B_item_type == CppAD::dynamic_enum)
	    B_this_row_type = CppAD::dynamic_enum;
	}
	B_row_type_below = B_this_row_type;
      }

      item_type = CppAD::constant_enum; // output type

      if(x1RowTypes[i] == CppAD::variable_enum || B_this_row_type == CppAD::variable_enum) {
	item_type = CppAD::variable_enum;
      } else {
	if(x1RowTypes[i] == CppAD::dynamic_enum || B_this_row_type == CppAD::dynamic_enum) {
	  item_type = CppAD::dynamic_enum;
	}
      }
      type_y[i + j*n1] = item_type;
    }
  }
  //  std::cout<<"forward type analysis"<<std::endl;
  // std::cout<<"type_x (A parts)"<<std::endl;
  //   for(int i = 0; i < n1; ++i) {
  //     for(int j = 0; j < n1; ++j) std::cout<< type_x[i + j*n1]<<"\t";
  //     std::cout<<std::endl;
  //   }
  //   std::cout<<"type_x (B parts)"<<std::endl;
  //   for(int i = 0; i < n1; ++i) {
  //     for(int j = 0; j < n2; ++j) std::cout<< type_x[n1sq + i + j*n1]<<"\t";
  //     std::cout<<std::endl;
  //   }
  //   std::cout<<"type_y"<<std::endl;
  //   for(int i = 0; i < n1; ++i) {
  //     for(int j = 0; j < n2; ++j) std::cout<< type_y[i + j*n1]<<"\t";
  //     std::cout<<std::endl;
  //   }

  return true;
}

bool atomic_backsolve_class::rev_depend(
					const CppAD::vector<double>&          parameter_x ,
					const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
					CppAD::vector<bool>&                depend_x    ,
					const CppAD::vector<bool>&          depend_y
					) {
  // if depend_y(i, j) is true, we need every element of depend_x that is an input to depend_y(i,j) to be set true
  // that means elements of depend_x for A(>=i, all j) (but only upper diagonal) and of B(>=i, j)
  //printf("In backsolve reverse_depend\n");
  int n = depend_x.size();
  int m = depend_y.size();
  int n1sq = n-m;
  int n1 = sqrt( static_cast<double>(n1sq) );
  int n2 = m/n1;

  int i_first_true_any_col, i_first_true_this_col;
  
  i_first_true_any_col = n1; // initialize past the end, i.e. no rows have depend_y true for this column (j)
  for(int j = 0; j < n2; ++j) {
    i_first_true_this_col = n1; // ditto
    // find the first row in the column that had depend_y true
    for(int i = 0; i < n1; ++i) {
      if(depend_y[i + j*n1]) {
	i_first_true_this_col = i;
	break;
      }
    }
    // keep track of the first row in any column that has depend_y true
    if(i_first_true_this_col < i_first_true_any_col) {
      i_first_true_any_col = i_first_true_this_col;
    }
    // set depend_x elements corresponding to B true for >= first row with depend_y true for this column
    for(int i = i_first_true_this_col; i < n1; ++i) {
      depend_x[n1sq + i + j*n1] = true;
    }
    for(int i = 0; i < i_first_true_this_col; ++i) {
      depend_x[n1sq + i + j*n1] = false; // false otherwise
    }
  }
  // set all upper_diagonal elements of A for rows >= first row with depend_y true for any column true
  for(int i = i_first_true_any_col; i < n1; ++i) {
    for(int j = i; j < n1; ++j) {
      depend_x[i + j * n1] = true;
    }
    for(int j = 0; j < i; ++j) {
      depend_x[i + j * n1] = false; // below-diagonal elements for row i
    }
  }
  for(int i = 0; i < i_first_true_any_col; ++i) { // other rows
    for(int j = 0; j < n1; ++j) {                     // all cols 
      depend_x[i + j * n1] = false;
    }    
  }
  // std::cout<<"reverse depend analysis"<<std::endl;
  //   std::cout<<"y_depend"<<std::endl;
  //   for(int i = 0; i < n1; ++i) {
  //     for(int j = 0; j < n2; ++j) std::cout<< depend_y[i + j*n1]<<"\t";
  //     std::cout<<std::endl;
  //   }
  //   std::cout<<"x_depend (A parts)"<<std::endl;
  //   for(int i = 0; i < n1; ++i) {
  //     for(int j = 0; j < n1; ++j) std::cout<< depend_x[i + j*n1]<<"\t";
  //     std::cout<<std::endl;
  //   }
  //   std::cout<<"x_depend (B parts)"<<std::endl;
  //   for(int i = 0; i < n1; ++i) {
  //     for(int j = 0; j < n2; ++j) std::cout<< depend_x[n1sq + i * j*n1]<<"\t";
  //     std::cout<<std::endl;
  //   }
  return true;
}

  
bool atomic_backsolve_class::forward(
				     const CppAD::vector<double>&               parameter_x  ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				     size_t                              need_y       ,
				     size_t                              order_low    ,
				     size_t                              order_up     ,
				     const CppAD::vector<double>&               taylor_x     ,
				     CppAD::vector<double>&                     taylor_y     ) {
  //forward mode
  //printf("In backsolve forward\n");
  int nrow = order_up + 1;
  //    std::cout<<"nrow = "<<nrow<<std::endl;
  //std::cout<<"tx size = "<<taylor_x.size()<<" ty size = "<<taylor_y.size()<<std::endl;
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq = n-m;
  int n1 = sqrt( static_cast<double>(n1sq) );
  int n2 = m/n1;
  //std::cout<<"n = "<<n<<" m = "<<m<<" n1 = "<<n1<<" n2 = "<<n2<<std::endl;

  EigenMap Ymap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
  EigenConstMap Amap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  if(order_low <= 0 & order_up >= 0) { // value
    //printf("In forward 0\n");
    // We could compile different cases depending on need for strides or not.
    EigenConstMap Bmap(&taylor_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // There is some indication in nimble EigenTypeDefs.h that we need to force a copy of B
    Ymap = Amap.template triangularView<Eigen::Upper>().solve(Bmap);      
  }
  if(order_low <= 1 & order_up >= 1) {
    //printf("In forward >1\n");
    //      solve(A, dB - dA * Y)
    EigenConstMap dA_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    EigenConstMap dB_map(&taylor_x[1 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    EigenMap dY_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    dY_map = Amap.template triangularView<Eigen::Upper>().solve(dB_map - dA_map.template triangularView<Eigen::Upper>() * Ymap).eval();// This .eval() is necessary and I don't understand why.  Normally that would be for aliasing, but there should be no over-lapping points here.  There are inter-woven maps, and that's weird but should work.

    // std::cout<<"dY"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n1; ++j) std::cout<< dY_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

  }
  return true;
}

bool atomic_backsolve_class::forward(
				     const CppAD::vector<CppAD::AD<double> >&               parameter_x  ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x       ,
				     size_t                              need_y       ,
				     size_t                              order_low    ,
				     size_t                              order_up     ,
				     const CppAD::vector<CppAD::AD<double> >&               taylor_x     ,
				     CppAD::vector<CppAD::AD<double> >&                     taylor_y     ) {
  //forward mode
  //printf("In backsolve forward\n");
  int nrow = order_up + 1;
  //    std::cout<<"nrow = "<<nrow<<std::endl;
  //std::cout<<"tx size = "<<taylor_x.size()<<" ty size = "<<taylor_y.size()<<std::endl;
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq = n-m;
  int n1 = sqrt( static_cast<double>(n1sq) );
  int n2 = m/n1;
  //std::cout<<"n = "<<n<<" m = "<<m<<" n1 = "<<n1<<" n2 = "<<n2<<std::endl;

  metaEigenMap Ymap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
  metaEigenConstMap Amap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  if(order_low <= 0 & order_up >= 0) { // value
    //printf("In forward 0\n");
    // We could compile different cases depending on need for strides or not.
    metaEigenConstMap Bmap(&taylor_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    // There is some indication in nimble EigenTypeDefs.h that we need to force a copy of B
    Ymap = nimDerivs_EIGEN_BS(Amap, Bmap);
    //    Ymap = Amap.template triangularView<Eigen::Upper>().solve(Bmap);      
  }
  if(order_low <= 1 & order_up >= 1) {
    //printf("In forward >1\n");
    //      solve(A, dB - dA * Y)
    metaEigenConstMap dA_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    metaEigenConstMap dB_map(&taylor_x[1 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    metaEigenMap dY_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );

    dY_map = nimDerivs_EIGEN_BS(Amap, dB_map - nimDerivs_matmult(dA_map.template triangularView<Eigen::Upper>(), Ymap));
    //    dY_map = Amap.template triangularView<Eigen::Upper>().solve(dB_map - dA_map * Ymap).eval();// This .eval() is necessary and I don't understand why.  Normally that would be for aliasing, but there should be no over-lapping points here.  There are inter-woven maps, and that's weird but should work.
  }
  return true;
}

bool atomic_backsolve_class::reverse(
				     const CppAD::vector<double>&               parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     size_t                              order_up    ,
				     const CppAD::vector<double>&               taylor_x    ,
				     const CppAD::vector<double>&               taylor_y    ,
				     CppAD::vector<double>&                     partial_x   ,
				     const CppAD::vector<double>&               partial_y   )
{
  //reverse mode
  int nrow = order_up + 1;
  //  std::cout<<"entering reverse----------------------------------"<<std::endl;
  //std::cout<<"nrow = "<<nrow<<std::endl;
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq = n-m;
  int n1 = sqrt( static_cast<double>(n1sq) );
  int n2 = m/n1;
    
  EigenConstMap Amap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  EigenMap Aadjoint_map(&partial_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  EigenMap Badjoint_map(&partial_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
  EigenConstMap Ymap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
  if(order_up >= 0) {
    /* EigenTemplateTypes<double>::typeEigenConstMapStrd Bmap(&taylor_x[0 + n1sq*nrow], */
    /* 							 n1, n2, EigStrDyn(nrow*n1, nrow ) ); */
    EigenConstMap Yadjoint_map(&partial_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
      
    Badjoint_map = Amap.transpose().template triangularView<Eigen::Lower>().solve(Yadjoint_map);
    if(order_up == 0) { // otherwise this gets included below
      Aadjoint_map = (-Badjoint_map * Ymap.transpose()).template triangularView<Eigen::Upper>();
    }
    // std::cout<<"Y"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Ymap(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    // std::cout<<"Yadjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Yadjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    // std::cout<<"Aadjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n1; ++j) std::cout<< Aadjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    // std::cout<<"Badjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Badjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

  }
  if(order_up >= 1) {
    //    std::cout<<"------ starting 2nd order----"<<std::endl;
    EigenConstMap Adot_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    EigenConstMap Ydot_adjoint_map(&partial_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    EigenConstMap Ydot_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    EigenMap Adot_adjoint_map(&partial_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    EigenMap Bdot_adjoint_map(&partial_x[1 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    

    // std::cout<<"Adot"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n1; ++j) std::cout<< Adot_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    // std::cout<<"Ydot_adjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Ydot_adjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    // std::cout<<"Ydot"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Ydot_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    
    // std::cout<<"Adot_adjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n1; ++j) std::cout<< Adot_adjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    // std::cout<<"Bdot_adjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Bdot_adjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

    Eigen::MatrixXd Adot_stuff = Adot_map.transpose().template triangularView<Eigen::Lower>(); // doesn't work inside solve

    /* Note: A^-1 Z A^-1 = solve(A, solve(A^T, Z^T)^T) */
      
    /* Badjoint = A^-T * Yadjoint - (A^-1 Adot A^-1)^T  Ydot_adjoint */
    Badjoint_map -= Amap.template triangularView<Eigen::Upper>().solve( Amap.transpose().template triangularView<Eigen::Lower>().solve(Adot_stuff).transpose() ).transpose() * Ydot_adjoint_map;
    // std::cout<<"Badjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Badjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    /* Bdot_adjoint = A^-T Ydot_adjoint */
    Bdot_adjoint_map = Amap.transpose().template triangularView<Eigen::Lower>().solve(Ydot_adjoint_map).eval(); // This eval is necessary.  I am not clear when evals are necessary after solves or not.  
    // std::cout<<"Bdot_adjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n2; ++j) std::cout<< Bdot_adjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    /* Aadjoint = -Badjoint * Y^T  - (A^-T Ydot_adjoint Ydot^T)= -Badjoint * Y^T  - ( Bdot_adjoint Ydot^T) , including both Badjoint terms*/
    Aadjoint_map = (-Badjoint_map * Ymap.transpose()).template triangularView<Eigen::Upper>();
    // std::cout<<"Aadjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n1; ++j) std::cout<< Aadjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    Aadjoint_map -= (Bdot_adjoint_map * Ydot_map.transpose()).template triangularView<Eigen::Upper>();  
    // std::cout<<"Aadjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n1; ++j) std::cout<< Aadjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }
    /* Adot_adjoint = -Bdot_adjoint Y^T */
    Adot_adjoint_map = (-Bdot_adjoint_map * Ymap.transpose()).template triangularView<Eigen::Upper>();
    // std::cout<<"Adot_adjoint"<<std::endl;
    // for(int i = 0; i < n1; ++i) {
    //   for(int j = 0; j < n1; ++j) std::cout<< Adot_adjoint_map(i, j)<<"\t";
    //   std::cout<<std::endl;
    // }

  }
  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

bool atomic_backsolve_class::reverse(
				     const CppAD::vector<CppAD::AD<double> >&               parameter_x ,
				     const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				     size_t                              order_up    ,
				     const CppAD::vector<CppAD::AD<double> >&               taylor_x    ,
				     const CppAD::vector<CppAD::AD<double> >&               taylor_y    ,
				     CppAD::vector<CppAD::AD<double> >&                     partial_x   ,
				     const CppAD::vector<CppAD::AD<double> >&               partial_y   )
{
  //reverse mode
  int nrow = order_up + 1;
  //std::cout<<"nrow = "<<nrow<<std::endl;
  int n = taylor_x.size()/nrow;
  int m = taylor_y.size() / nrow;
  int n1sq = n-m;
  int n1 = sqrt( static_cast<double>(n1sq) );
  int n2 = m/n1;
    
  metaEigenConstMap Amap(&taylor_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  metaEigenMap Aadjoint_map(&partial_x[0], n1, n1, EigStrDyn(nrow*n1, nrow) );
  metaEigenMap Badjoint_map(&partial_x[0 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
  metaEigenConstMap Ymap(&taylor_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );
  if(order_up >= 0) {
    /* EigenTemplateTypes<double>::typeEigenConstMapStrd Bmap(&taylor_x[0 + n1sq*nrow], */
    /* 							 n1, n2, EigStrDyn(nrow*n1, nrow ) ); */
    metaEigenConstMap Yadjoint_map(&partial_y[0], n1, n2, EigStrDyn(nrow*n1, nrow ) );

    Badjoint_map = nimDerivs_EIGEN_FS(Amap.transpose(), Yadjoint_map);
    // Badjoint_map = Amap.transpose().template triangularView<Eigen::Lower>().solve(Yadjoint_map);
    // We could potentially make a nimDerivs_matmult that takes advantage of needing only upper triangular results
    if(order_up == 0) Aadjoint_map = nimDerivs_matmult(-Badjoint_map, Ymap.transpose()).template triangularView<Eigen::Upper>();//-Badjoint_map * Ymap.transpose(); // otherwise this gets included below
  }
  if(order_up >= 1) {
    metaEigenConstMap Adot_map(&taylor_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    metaEigenConstMap Ydot_adjoint_map(&partial_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    metaEigenConstMap Ydot_map(&taylor_y[1], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    metaEigenMap Adot_adjoint_map(&partial_x[1], n1, n1, EigStrDyn(nrow*n1, nrow) );
    metaEigenMap Bdot_adjoint_map(&partial_x[1 + n1sq*nrow], n1, n2, EigStrDyn(nrow*n1, nrow ) );
    

    /* Note: A^-1 Z A^-1 = solve(A, solve(A^T, Z^T)^T) */
      
    /* Badjoint = A^-T * Yadjoint - (A^-1 Adot A^-1)^T  Ydot_adjoint */      

    MatrixXd_CppAD Adot_stuff = Adot_map.transpose().template triangularView<Eigen::Lower>(); // doesn't work inside solve

    
    Badjoint_map -= nimDerivs_matmult( nimDerivs_EIGEN_BS( Amap,
							    nimDerivs_EIGEN_FS(Amap.transpose(), Adot_stuff).transpose() ).transpose() ,
				       Ydot_adjoint_map);

    //    Badjoint_map -= Amap.template triangularView<Eigen::Upper>().solve( Amap.transpose().template triangularView<Eigen::Lower>().solve(Adot_map.transpose()).transpose() ).transpose() * Ydot_adjoint_map;
    /* Bdot_adjoint = A^-T Ydot_adjoint */
    Bdot_adjoint_map = nimDerivs_EIGEN_FS(Amap.transpose(), Ydot_adjoint_map);
    //    Bdot_adjoint_map = Amap.transpose().template triangularView<Eigen::Lower>().solve(Ydot_adjoint_map).eval(); // This eval is necessary.  I am not clear when evals are necessary after solves or not.  
    /* Aadjoint = -Badjoint * Y^T  - (A^-T Ydot_adjoint Ydot^T)= -Badjoint * Y^T  - ( Bdot_adjoint Ydot^T) , including both Badjoint terms*/
    Aadjoint_map = nimDerivs_matmult(-Badjoint_map, Ymap.transpose()).template triangularView<Eigen::Upper>();
    Aadjoint_map -= nimDerivs_matmult(Bdot_adjoint_map, Ydot_map.transpose()).template triangularView<Eigen::Upper>();
    //    Aadjoint_map = -Badjoint_map * Ymap.transpose() - Bdot_adjoint_map * Ydot_map.transpose();
    /* Adot_adjoint = -Bdot_adjoint Y^T */
    Adot_adjoint_map = nimDerivs_matmult(-Bdot_adjoint_map, Ymap.transpose()).template triangularView<Eigen::Upper>();
    //    Adot_adjoint_map = -Bdot_adjoint_map * Ymap.transpose();
  }
  if(order_up >= 2) {
    printf("Unsupported reverse order requested\n");
    return false;
  }
  return true;
}

void atomic_backsolve(const MatrixXd_CppAD &A,
		      const MatrixXd_CppAD &B,
		      MatrixXd_CppAD &Y) {
  // Y = backsolve(A, B) = A^-1 B
  // A is n1-x-n1
  // B and Y are n1-x-n2
  static atomic_backsolve_class atomic_backsolve("atomic_backsolve");
  int n1 = A.rows();
  int n2 = B.cols();
  std::vector<CppAD::AD<double> > xVec(n1*n1 + n1*n2);
  mat2vec(A, xVec);
  mat2vec(B, xVec, n1*n1);
  std::vector<CppAD::AD<double> > yVec(n1*n2);
  atomic_backsolve(xVec, yVec);
  Y.resize(n1, n2);
  vec2mat(yVec, Y);  
}

MatrixXd_CppAD nimDerivs_EIGEN_BS(const MatrixXd_CppAD &A,
				   const MatrixXd_CppAD &B) {
  MatrixXd_CppAD ans;
  atomic_backsolve(A, B, ans);
  return ans;
}
