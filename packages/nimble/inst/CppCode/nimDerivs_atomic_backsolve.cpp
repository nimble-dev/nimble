#include <nimble/nimDerivs_atomic_backsolve.h>
#include <nimble/nimDerivs_atomic_forwardsolve.h>
#include <nimble/nimDerivs_atomic_cache.h>

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

#define USE_NEW_DYNAMIC_FS

bool atomic_backsolve_class::for_type(
				      const CppAD::vector<double>&               parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
#ifdef VERBOSE_ATOMIC_FS
  printf("In backsolve for_type\n");
#endif
  // Row i of output depends on rows i and <i in A input.
  // Element (i,j) of output depends on elements (<=i, j) of B input
  int n = type_x.size();
  int m = type_y.size();
  int n1sq, n1, n2;
  if((!Aconstant()) && (!Bconstant())) n1sq = n-m;
  if(Aconstant()) n1sq = get_X_stored().size();
  if(Bconstant()) n1sq = n;
  if(Aconstant() && Bconstant())
    std::cout<<"atomic_backsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

  CppAD::vector<CppAD::ad_type_enum> x1RowTypes(n1); //A
  CppAD::vector<CppAD::ad_type_enum> x2ColTypes(n2); //B

  CppAD::ad_type_enum row_type_below(CppAD::constant_enum);
  CppAD::ad_type_enum this_row_type;
  CppAD::ad_type_enum item_type;
  for(int i = n1-1; i >= 0; --i) { // go backwards because rows depend on those below them.  Must use int to get valid >= 0 comparison.
    this_row_type = row_type_below;
    if(!Aconstant()) {
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
    }
    x1RowTypes[i] = row_type_below = this_row_type;
  }

  CppAD::ad_type_enum B_item_type;
  CppAD::ad_type_enum B_row_type_below;  // type of elements below row i, done for one column a time
  CppAD::ad_type_enum  B_this_row_type;  // type of the current row, done for one column at a time
  int Boffset = Aconstant() ? 0 : n1sq;

  for(int j = 0; j < n2; ++j) {
    // for the j-th column
    B_row_type_below = CppAD::constant_enum;
    for(int i = n1-1; i >= 0; --i) { // must be int, not size_t
      B_this_row_type = B_row_type_below; // can't be lower type than what's below it in the current column
      // if B_this_row_type is not already variable (max type), update based on type of element B[i,j]
      if(B_this_row_type < CppAD::variable_enum) {
	if(!Bconstant()) {
	  B_item_type = type_x[Boffset + i + j*n1];
	  if(B_item_type == CppAD::variable_enum) {
	    B_this_row_type = CppAD::variable_enum;
	  } else {
	    if(B_item_type == CppAD::dynamic_enum)
	      B_this_row_type = CppAD::dynamic_enum;
	  }
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
#ifdef VERBOSE_ATOMIC_FS
  printf("In backsolve rev_depend\n");
#endif
  int n = depend_x.size();
  int m = depend_y.size();
  int n1sq, n1, n2;
  if((!Aconstant()) && (!Bconstant())) n1sq = n-m;
  if(Aconstant()) n1sq = get_X_stored().size();
  if(Bconstant()) n1sq = n;
  if(Aconstant() && Bconstant())
    std::cout<<"atomic_backsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

  int Boffset = Aconstant() ? 0 : n1sq;

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
    if(!Bconstant()) {
      for(int i = i_first_true_this_col; i < n1; ++i) {
	depend_x[Boffset + i + j*n1] = true;
      }
      for(int i = 0; i < i_first_true_this_col; ++i) {
	depend_x[Boffset + i + j*n1] = false; // false otherwise
      }
    }
  }
  if(!Aconstant()) {
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
  }
  return true;
}

#define ND_SOLVE(A, B) (A).template triangularView<Eigen::Upper>().solve(B)
#define ND_MAIN_TRI(A) (A).template triangularView<Eigen::Upper>()
#define ND_TRANS_SOLVE(A, B) (A).transpose().template triangularView<Eigen::Lower>().solve(B)
#define ND_TRANS_TRI(A) (A).transpose().template triangularView<Eigen::Lower>()
#define ATOMIC_SOLVE_CLASS atomic_backsolve_class
#define ND_META_SOLVE(A, B) nimDerivs_EIGEN_BS(A, B)
#define ND_META_TRANS_SOLVE(A, B) nimDerivs_EIGEN_FS(A, B)

#include "nimDerivs_atomic_solve_generic.cpp"

#undef ND_SOLVE
#undef ND_MAIN_TRI
#undef ND_TRANS_SOLVE
#undef ND_TRANS_TRI
#undef ATOMIC_SOLVE_CLASS
#undef ND_META_SOLVE
#undef ND_META_TRANS_SOLVE

bool check_A_diagonal_upper(const MatrixXd_CppAD &A) {
  // THIS IGNORES LOWER OFF-DIAGONAL. IT IS FOR backsolve
  int n1 = A.rows();
  if(A.cols() != n1)
    cout<<"A is not square in check_A_diagonal"<<endl;
  bool diagonal(true);
  for(int iii = 0; iii < n1-1; ++iii) {
    for(int jjj = iii+1; jjj < n1; ++jjj) {
      diagonal &= CppAD::IdenticalZero(A(iii, jjj));
      if(!diagonal) break;
    }
  }
  return diagonal;
}


void atomic_backsolve(const MatrixXd_CppAD &A,
		      const MatrixXd_CppAD &B,
		      MatrixXd_CppAD &Y) {
  // Y = backsolve(A, B) = A^-1 B, with A upper-triangular
  // This will handle all of A const vs var and/or all of B const vs var.
  // Then it will cut down work based on regions of zeros.
  // For A, the case of interest is diagonal (i.e. zero on the below-diagonal).
  // For B, the case of interest is starting and ending regions of zeros
  // In the future blocks of var within const or such can be implemented.
  //
  // First set up constant cases and put all testing in place.
  int n1 = A.rows();
  int n2 = B.cols();
  //  cout<<n1<<" "<<n2<<endl;
  if(A.cols() != n1)
    cout<<"A is not square in atomic_backsolve"<<endl;
  if(B.rows() != n1)
    cout<<"incommensurate matrices in atomic_backsolve"<<endl;

  Y.resize(n1, n2);

  // Cases are handled as follows:
  // - A is 1x1: Use scalar elements and let CppAD handle all var vs const issues.
  // - A is diagonal: Again use scalar elements
  // - B has some last set of rows with all constant zeros: fill in zero results and recurse on the non-zero solve
  // - (From here down, we have recursed into a nested solve).
  // - A and B are both all constant: Solve directly (no atomic) and enter results in Y
  // - (In the future, more complicated handling of portions of A and/or B constant can be inserted here.)
  // - A or B or neither are constant: record atomic

  // A is 1x1:
  // if(n1 == 1) {
  //   for(int j = 0; j < n2; ++j)
  //     Y(0, j) = B(0, j) / A(0, 0);
  //   return;
  // }

  // A is diagonal (This subsumes the case that A is 1x1)
  if(check_A_diagonal_upper(A)) {
    for(int i = 0; i < n1; ++i) {
      for(int j = 0; j < n2; ++j) {
	Y(i, j) = B(i, j) / A(i, i);
      }
    }
    return;
  }

  // Check constant zero region for B
  int BrowStartNZ, BrowEndNZ, BcolStartNZ, BcolEndNZ;
  bool B_is_constant_zero;

  auto zero_cond = [](const CppAD::AD<double> &x)->bool {return CppAD::IdenticalZero(x);};
  B_is_constant_zero = delineate_condition_region(zero_cond, B,
						  BrowStartNZ, BrowEndNZ, BcolStartNZ, BcolEndNZ);
  //  std::cout<<"checking: "<< BrowStartNZ<<" "<< BrowEndNZ<<" "<< BcolStartNZ<<" "<< BcolEndNZ<<" "<<B_is_constant_zero<<std::endl;
  if((BrowEndNZ < n1)  || B_is_constant_zero) {
    // B has some final rows with all constant zeros: fill in zero results and recurse to solve non-zero region
    // note that if B is all constant, we'll have BrowStartNZ == BrowEndNZ == n1
    int lowest_constant_zero_row = B_is_constant_zero ? 0 : BrowEndNZ;
    for(int i = n1-1; i >= lowest_constant_zero_row; --i) {
      for(int j = 0; j < n2; ++j) {
  	Y(i, j) = 0;
      }
    }
    if(!B_is_constant_zero) {    // Recurse on a nested forward solve
      Y.block(0, 0, BrowEndNZ, n2) = nimDerivs_EIGEN_BS(A.block(0, 0, BrowEndNZ, BrowEndNZ), B.block(0, 0, BrowEndNZ, n2) );
    }
    return;
  }

  // Find constant vs non-constant regions
  int ArowStart, ArowEnd, AcolStart, AcolEnd;
  int BrowStart, BrowEnd, BcolStart, BcolEnd;
  bool A_is_constant, B_is_constant;

  auto const_cond = [](const CppAD::AD<double> &x)->bool {return CppAD::Constant(x);};
  A_is_constant = delineate_condition_region(const_cond, A,
					     ArowStart, ArowEnd, AcolStart, AcolEnd,
					     true, true, false);
  B_is_constant = delineate_condition_region(const_cond, B,
					     BrowStart, BrowEnd, BcolStart, BcolEnd);
  // std::cout<<"checking2: "<< BrowStart<<" "<< BrowEnd<<" "<< BcolStart<<" "<< BcolEnd<<" "<<B_is_constant<<std::endl;
  // std::cout<<"checking3: "<< ArowStart<<" "<< ArowEnd<<" "<< AcolStart<<" "<< AcolEnd<<" "<<A_is_constant<<std::endl;
  // - A and B are both all constant: Solve directly (no atomic) and enter results in Y
  if(A_is_constant && B_is_constant) { // The whole problem, both A and B, are constant
    Y = A.template triangularView<Eigen::Upper>().solve(B);
    return;
  }

  // - A or B or neither are constant: record atomic
  atomic_backsolve_class *atomic_backsolve;
  atomic_backsolve = new atomic_backsolve_class("atomic_backsolve");
  atomic_backsolve->Aconstant() = A_is_constant;
  atomic_backsolve->Bconstant() = B_is_constant;
  
  int xVecSize = 0;
  if(!A_is_constant) xVecSize += n1*n1;
  if(!B_is_constant) xVecSize += n1*n2;
  std::vector<CppAD::AD<double> > xVec(xVecSize);
  if(A_is_constant) {
    atomic_backsolve->set_X_stored(A);
  } else {
    mat2vec(A, xVec, 0);
  }
  if(B_is_constant) {
    atomic_backsolve->set_X_stored(B);
  } else {
    if(A_is_constant) {
      mat2vec(B, xVec, 0);
    } else {
      mat2vec(B, xVec, n1*n1);
    }
  }
  //  cout<<xVec.size()<<endl;

  std::vector<CppAD::AD<double> > yVec(n1*n2);
  (*atomic_backsolve)(xVec, yVec);
  // dummy test : for(int i = 0; i < n1*n2; ++i) yVec[i] = xVec[i];
  vec2mat(yVec, Y);
  if(CppAD::AD<double>::get_tape_handle_nimble() == nullptr) {
    delete atomic_backsolve;
  } else {
    track_nimble_atomic(atomic_backsolve,
			CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr());
  }

}



void atomic_backsolve_old(const MatrixXd_CppAD &A,
			  const MatrixXd_CppAD &B,
			  MatrixXd_CppAD &Y) {
  // Y = backsolve(A, B) = A^-1 B
  // A is n1-x-n1
  // B and Y are n1-x-n2
  //  static atomic_backsolve_class atomic_backsolve("atomic_backsolve");
  atomic_backsolve_class *atomic_backsolve;
  int n1 = A.rows();
  int n2 = B.cols();
  std::vector<CppAD::AD<double> > xVec(n1*n1 + n1*n2);
  mat2vec(A, xVec);
  mat2vec(B, xVec, n1*n1);
  std::vector<CppAD::AD<double> > yVec(n1*n2);
  atomic_backsolve = new atomic_backsolve_class("atomic_backsolve");
  (*atomic_backsolve)(xVec, yVec);
  Y.resize(n1, n2);
  vec2mat(yVec, Y);  
}

MatrixXd_CppAD nimDerivs_EIGEN_BS(const MatrixXd_CppAD &A,
				   const MatrixXd_CppAD &B) {
  MatrixXd_CppAD ans;
  atomic_backsolve(A, B, ans);
  return ans;
}
