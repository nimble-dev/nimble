#include <nimble/nimDerivs_atomic_forwardsolve.h>
#include <nimble/nimDerivs_atomic_backsolve.h>
#include <nimble/nimDerivs_atomic_cache.h>

/* 
forward solve relations are like backsolve relations, with forward and back swapped where needed.
*/

// #define VERBOSE_ATOMIC_FS

#define USE_NEW_DYNAMIC_FS


bool atomic_forwardsolve_class::for_type(
				      const CppAD::vector<double>&               parameter_x ,
				      const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
				      CppAD::vector<CppAD::ad_type_enum>&        type_y      )
{
#ifdef VERBOSE_ATOMIC_FS
  printf("In forwardsolve for_type\n");
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
    std::cout<<"atomic_forwardsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

 // HOT-WIRED FOR TESTING
  // std::cout<<"Remember to turn off hot-wired for_type and rev_depend"<<std::endl;
  // for(size_t i = 0; i < type_y.size(); ++i)
  //   type_y[i] = CppAD::variable_enum;
  // return true;
  
  CppAD::vector<CppAD::ad_type_enum> x1RowTypes(n1); //A
  CppAD::vector<CppAD::ad_type_enum> x2ColTypes(n2); //B

  CppAD::ad_type_enum row_type_above(CppAD::constant_enum);
  CppAD::ad_type_enum this_row_type;
  CppAD::ad_type_enum item_type;
  for(int i = 0; i < n1; ++i) { 
    this_row_type = row_type_above;
    if(!Aconstant()) {
      if(this_row_type != CppAD::variable_enum) { // no need to check if row type is already at the "max"
	for(size_t j = 0; j <= i; ++j) { // only look at lower triangle
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
    x1RowTypes[i] = row_type_above = this_row_type;
  }

  CppAD::ad_type_enum B_item_type;
  CppAD::ad_type_enum B_row_type_above;  // type of elements above row i, done for one column a time
  CppAD::ad_type_enum  B_this_row_type;  // type of the current row, done for one column at a time
  int Boffset = Aconstant() ? 0 : n1sq;
  for(int j = 0; j < n2; ++j) {
    // for the j-th column
    B_row_type_above = CppAD::constant_enum;
    for(int i = 0; i < n1; ++i) { 
      B_this_row_type = B_row_type_above; // can't be lower type than what's above it in the current column
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
	B_row_type_above = B_this_row_type;
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
  return true;
}

bool atomic_forwardsolve_class::rev_depend(
					const CppAD::vector<double>&          parameter_x ,
					const CppAD::vector<CppAD::ad_type_enum>&  type_x      ,
					CppAD::vector<bool>&                depend_x    ,
					const CppAD::vector<bool>&          depend_y
					) {
  // if depend_y(i, j) is true, we need every element of depend_x that is an input to depend_y(i,j) to be set true
  // that means elements of depend_x for A(<=i, all j) (but only lower diagonal) and of B(<=i, j)
  //printf("In forward reverse_depend\n");
#ifdef VERBOSE_ATOMIC_FS
  printf("In forwardsolve rev_depend\n");
#endif
  int n = depend_x.size();
  int m = depend_y.size();
  int n1sq, n1, n2;
  if((!Aconstant()) && (!Bconstant())) n1sq = n-m;
  if(Aconstant()) n1sq = get_X_stored().size();
  if(Bconstant()) n1sq = n;
  if(Aconstant() && Bconstant())
    std::cout<<"atomic_forwardsolve is being used with both A and B constant.  This should not happen."<<std::endl;
  n1 = sqrt( static_cast<double>(n1sq) );
  n2 = m/n1;

  // HOT-WIRED FOR TESTING
  // for(size_t i = 0; i < depend_x.size(); ++i)
  //   depend_x[i] = true;
  // return true;
  
  int Boffset = Aconstant() ? 0 : n1sq;

  int i_last_true_any_col, i_last_true_this_col;
  
  i_last_true_any_col = -1; // initialize past the end, i.e. no rows have depend_y true for this column (j)
  for(int j = 0; j < n2; ++j) {
    i_last_true_this_col = -1; // ditto
    // find the last row in the column that had depend_y true
    for(int i = n1-1; i >=0 ; --i) {
      if(depend_y[i + j*n1]) {
	i_last_true_this_col = i;
	break;
      }
    }
    // keep track of the last row in any column that has depend_y true
    if(i_last_true_this_col > i_last_true_any_col) {
      i_last_true_any_col = i_last_true_this_col;
    }
    // set depend_x elements corresponding to B true for <= last row with depend_y true for this column
    if(!Bconstant()) {
      for(int i = i_last_true_this_col; i >=0 ; --i) {
	depend_x[Boffset + i + j*n1] = true;
      }
      for(int i = n1-1; i > i_last_true_this_col; --i) {
	depend_x[Boffset + i + j*n1] = false;
      }
    }
  }
  if(!Aconstant()) {
    // set all lower diagonal elements of A for rows <= last row with depend_y true for any column true
    for(int i = i_last_true_any_col; i >= 0; --i) {
      for(int j = 0; j <= i; ++j) {
	depend_x[i + j * n1] = true;
      }
      for(int j = i+1; j < n1; ++j) {
	depend_x[i + j * n1] = false;
      }    
    }
    for(int i = n1-1; i > i_last_true_any_col; --i) { // all other rows
      for(int j = 0; j < n1; ++j) {                   // all cols
	depend_x[i + j * n1] = false;
      }
    }
  }
  return true;
}

#define ND_SOLVE(A, B) (A).template triangularView<Eigen::Lower>().solve(B)
#define ND_MAIN_TRI(A) (A).template triangularView<Eigen::Lower>()
#define ND_TRANS_SOLVE(A, B) (A).transpose().template triangularView<Eigen::Upper>().solve(B)
#define ND_TRANS_TRI(A) (A).transpose().template triangularView<Eigen::Upper>()
#define ATOMIC_SOLVE_CLASS atomic_forwardsolve_class
#define ND_META_SOLVE(A, B) nimDerivs_EIGEN_FS(A, B)
#define ND_META_TRANS_SOLVE(A, B) nimDerivs_EIGEN_BS(A, B)

#include "nimDerivs_atomic_solve_generic.cpp"

#undef ND_SOLVE
#undef ND_MAIN_TRI
#undef ND_TRANS_SOLVE
#undef ND_TRANS_TRI
#undef ATOMIC_SOLVE_CLASS
#undef ND_META_SOLVE
#undef ND_META_TRANS_SOLVE

bool check_A_diagonal_lower(const MatrixXd_CppAD &A) {
  // THIS IGNORES UPPER OFF-DIAGONAL. IT IS FOR FORWARD_SOLVE
  int n1 = A.rows();
  if(A.cols() != n1)
    cout<<"A is not square in check_A_diagonal_lower"<<endl;
  bool diagonal(true);
  for(int iii = 0; iii < n1; ++iii) {
    for(int jjj = 0; jjj < iii; ++jjj) {
      diagonal &= CppAD::IdenticalZero(A(iii, jjj));
      if(!diagonal) break;
    }
  }
  return diagonal;
}

void atomic_forwardsolve(const MatrixXd_CppAD &A,
			 const MatrixXd_CppAD &B,
			 MatrixXd_CppAD &Y) {
  // Y = forwardsolve(A, B) = A^-1 B, with A lower-triangular
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
    cout<<"A is not square in atomic_forwardsolve"<<endl;
  if(B.rows() != n1)
    cout<<"incommensurate matrices in atomic_forwardsolve"<<endl;

  Y.resize(n1, n2);

  // Cases are handled as follows:
  // - A is 1x1: Use scalar elements and let CppAD handle all var vs const issues.
  // - A is diagonal: Again use scalar elements
  // - B has some first set of rows with all constant zeros: fill in zero results and recurse on the non-zero solve
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
  if(check_A_diagonal_lower(A)) {
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
  // This is called to set the last four arguments by reference
  B_is_constant_zero = delineate_condition_region(zero_cond, B,
						  BrowStartNZ, BrowEndNZ, BcolStartNZ, BcolEndNZ);
  if(BrowStartNZ > 0) {
    // B has some first set of rows with all constant zeros: fill in zero results and recurse to solve non-zero region
    for(int i = 0; i < BrowStartNZ; ++i) {
      for(int j = 0; j < n2; ++j) {
	Y(i, j) = 0;
      }
    }
    if(BrowStartNZ != BrowEndNZ) {    // Recurse on a nested forward solve
      Y.block(BrowStartNZ, 0, n1 - BrowStartNZ, n2) = nimDerivs_EIGEN_FS(A.block(BrowStartNZ, BrowStartNZ, n1 - BrowStartNZ, n1 - BrowStartNZ), B.block(BrowStartNZ, 0, n1 - BrowStartNZ, n2) );
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
					     true, false, true);
  B_is_constant = delineate_condition_region(const_cond, B,
					     BrowStart, BrowEnd, BcolStart, BcolEnd);

  //  cout<<A_is_constant<<" "<<B_is_constant<<endl;
  // - A and B are both all constant: Solve directly (no atomic) and enter results in Y
  if(A_is_constant && B_is_constant) { // The whole problem, both A and B, are constant
    Y = A.template triangularView<Eigen::Lower>().solve(B);
    return;
  }

  auto dyn_cond = [](const CppAD::AD<double> &x)->bool {return CppAD::Dynamic(x);};
  bool A_is_dynamic(false);
  int a, b, c, d; // dummies
  if(!A_is_constant)
    A_is_dynamic = delineate_condition_region(dyn_cond, A, a, b, c, d);
  bool B_is_dynamic(false);
  if(!B_is_constant)
    B_is_dynamic = delineate_condition_region(dyn_cond, B, a, b, c, d);
  
  // - A or B or neither are constant: record atomic
  atomic_forwardsolve_class *atomic_forwardsolve;
  bool recording = CppAD::AD<double>::get_tape_handle_nimble() != nullptr;
  if(!recording) {
    atomic_forwardsolve = new atomic_forwardsolve_class("atomic_forwardsolve");
  } else {
    void *tape_mgr = CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr();
    atomic_forwardsolve = new_atomic_forwardsolve(tape_mgr, "atomic_forwardsolve");
  }
  atomic_forwardsolve->Aconstant() = A_is_constant;
  atomic_forwardsolve->Bconstant() = B_is_constant;

  atomic_forwardsolve->Avariable() = !(A_is_constant || A_is_dynamic);
  atomic_forwardsolve->Bvariable() = !(B_is_constant || B_is_dynamic);

  int xVecSize = 0;
  if(!A_is_constant) xVecSize += n1*n1;
  if(!B_is_constant) xVecSize += n1*n2;
  std::vector<CppAD::AD<double> > xVec(xVecSize);
  if(A_is_constant) {
    atomic_forwardsolve->set_X_stored(A);
  } else {
    mat2vec(A, xVec, 0);
  }
  if(B_is_constant) {
    atomic_forwardsolve->set_X_stored(B);
  } else {
    if(A_is_constant) {
      mat2vec(B, xVec, 0);
    } else {
      mat2vec(B, xVec, n1*n1);
    }
  }
  // cout<<xVec.size()<<endl;

  std::vector<CppAD::AD<double> > yVec(n1*n2);
  (*atomic_forwardsolve)(xVec, yVec);
  // dummy test : for(int i = 0; i < n1*n2; ++i) yVec[i] = xVec[i];
  vec2mat(yVec, Y);
  if(!recording) {
    delete atomic_forwardsolve;
  } else {
    track_nimble_atomic(atomic_forwardsolve,
			CppAD::AD<double>::get_tape_handle_nimble()->nimble_CppAD_tape_mgr_ptr(),
			CppAD::local::atomic_index_info_vec_manager_nimble<double>::manage() );
  }
}


MatrixXd_CppAD nimDerivs_EIGEN_FS(const MatrixXd_CppAD &A,
				  const MatrixXd_CppAD &B) {
#ifdef VERBOSE_ATOMIC_FS
  std::cout<<"Entering nimDerivs_EIGEN_FS"<<std::endl;
#endif
  MatrixXd_CppAD ans;
  atomic_forwardsolve(A, B, ans);
#ifdef VERBOSE_ATOMIC_FS
  std::cout<<"Leaving nimDerivs_EIGEN_FS"<<std::endl;
#endif
  return ans;
}
