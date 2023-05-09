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
void mat2vec_lower_zero(const MT &m, VT &v, size_t offset = 0) {
  size_t nrow = m.rows();
  size_t ncol = m.cols();
  for(size_t i=0; i < nrow; ++i) {
    for(size_t j=0; j < ncol; ++j) {
      if(j < i) v[i + j*nrow + offset] = 0;
      else      v[i + j*nrow + offset] = m(i,j);
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

template<typename Cond>
bool delineate_condition_region(const Cond &cond,
				const MatrixXd_CppAD &x,
				int &rowStart_arg, int &rowEnd_arg,
				int &colStart_arg, int &colEnd_arg,
				bool reset = true,
				bool ignore_lower = false,
				bool ignore_upper = false) {
  // This generically checks for a box within x such that condition cond is false.
  // In one use, the condition is that elements are CppAD constants.
  // In another use, the condition is that they have value 0.
  
  // All "end" values are C-style, that is one past the last valid index.

  // Result is returned through rowStart_arg, rowEnd_arg, colStart_arg, colEnd_arg.
  // The "arg" values are ignored upon input if reset = true.
  // If reset = false, they are used to mark the region to check.
  // In either case, they are used to return the result.
  //
  // If ignore_lower is true, then strictly lower diagonal values will not be checked.
  // This means the result behave as if the condition is true for all those values.
  // If ignore_upper is true, then strictly upper diagonal values will not be checked.

  if(ignore_lower && ignore_upper)
    std::cout<<"delineate_constant_region should not be called with both ignore_lower and ignore_upper true."<<std::endl;
  
  int nRow = x.rows();
  int nCol = x.cols();
  // These four values will store the result until the end, when the result
  // will be moved into the "[]_arg" reference variables.
  // They will be initialized as if reset = TRUE.
  // Note that they are flipped, so rowStart begins at nRow, the value it should take
  // if all rows satisfy cond(ition).
  // rowEnd starts at 0, but if all rows satisfy cond, it will be nRow upon completion.
  int rowStart(nRow);  // first row that is not all constant
  int rowEnd(0);       // one past last row that is not all constant
  int colStart(nCol);  // first col that is not all constant
  int colEnd(0);       // one past last col that is not all constant
  // If not resetting, take input regions by flipping input values.
  if(!reset) {
    rowStart = rowEnd_arg;
    rowEnd = rowStart_arg; 
    colStart = colEnd_arg;
    colEnd = colStart_arg;
  } else {
    // Use the input variables as for loop bounds.
    // If !reset, the above clause achieves this, so it is only done if reset==true.
    rowStart_arg = rowEnd; // These are deliberately flipped.  Default 0
    rowEnd_arg = rowStart; // Default nRow
    colStart_arg = colEnd; // Default 0
    colEnd_arg = colStart; // Default nCol
  }

  auto cond2 = [&](int i, int j) {
    if(ignore_lower)
      return i < j ? true : cond(x(i,j));
    if(ignore_upper)
      return i > j ? true : cond(x(i,j));
    return cond(x(i,j));
  };
  
  // For each row, counting up from start, check if all columns satisfy cond
  for(int i = rowStart_arg; i < rowEnd_arg; ++i) {
    bool this_row_true(true);
    for(int j = colStart_arg; j < colEnd_arg; ++j) {
      if(!(this_row_true &= cond2(i, j ))) break;
    }
    if(!this_row_true) { // mark the first row that has an element failing cond
      rowStart = i;
      break;
    }
  }

  // If all rows satsify cond, there is no need to iterate again for rowEnd
  if(rowStart == rowEnd_arg) {
    rowEnd = rowEnd_arg;
  } else {
    // For each row, counting down from end, check if all columns satsify cond
    for(int i = rowEnd_arg-1; i >= rowStart_arg; --i) {
      bool this_row_true(true);
      for(int j = colStart_arg; j < colEnd_arg; ++j) {
	if(!(this_row_true &= cond2(i, j ))) break;
      }
      if(!this_row_true) {
	rowEnd = i + 1; // mark one past the highest row (first from counting down) that has an element failing cond
	break;
      }
    }
  }

  // For each col, counting up from start, check if all rows satsify cond
  for(int j = colStart_arg; j < colEnd_arg; ++j) {
    bool this_col_true(true);
    for(int i = rowStart_arg; i < rowEnd_arg; ++i) {
      if(!(this_col_true &= cond2(i, j ))) break;
    }
    if(!this_col_true) {
      colStart = j; // mark first col that has an element failing cond
      break;
    }
  }

  // If all cols satisfy cond, there is no need to iterate again for colEnd
  if(colStart == colEnd_arg) {
    colEnd = colEnd_arg;
  } else {
    // For each col, counting down from end, check if all rows satisfy cond
    for(int j = colEnd_arg-1; j >= colStart_arg; --j) {
      bool this_col_true(true);
      for(int i = rowStart_arg; i < rowEnd_arg; ++i) {
	if(!(this_col_true &= cond2(i, j ))) break;
      }
      if(!this_col_true) {
	colEnd = j + 1; // mark one past last highest col (first from counting down) that has an element failing cod
	break;
      }
    }
  }
  // Does the entire block satisfy cond?
  bool all_true = (rowStart == rowEnd_arg) && (rowEnd == rowEnd_arg) && (colStart == colEnd_arg) && (colEnd == colEnd_arg);
  rowStart_arg = rowStart;
  rowEnd_arg = rowEnd;
  colStart_arg = colStart;
  colEnd_arg = colEnd;

  return all_true;
}



template<typename Cond>
bool delineate_condition_region(const Cond &cond,
				const MatrixXd_CppAD &x,
				size_t &rowStart_arg, size_t &rowEnd_arg,
				size_t &colStart_arg, size_t &colEnd_arg,
				bool reset = true,
				bool ignore_lower = false,
				bool ignore_upper = false) {
  // This generically checks for a box within x such that condition cond is false.
  // In one use, the condition is that elements are CppAD constants.
  // In another use, the condition is that they have value 0.
  
  // All "end" values are C-style, that is one past the last valid index.

  // Result is returned through rowStart_arg, rowEnd_arg, colStart_arg, colEnd_arg.
  // The "arg" values are ignored upon input if reset = true.
  // If reset = false, they are used to mark the region to check.
  // In either case, they are used to return the result.
  //
  // If ignore_lower is true, then strictly lower diagonal values will not be checked.
  // This means the result behave as if the condition is true for all those values.
  // If ignore_upper is true, then strictly upper diagonal values will not be checked.

  if(ignore_lower && ignore_upper)
    std::cout<<"delineate_constant_region should not be called with both ignore_lower and ignore_upper true."<<std::endl;
  
  size_t nRow = x.rows();
  size_t nCol = x.cols();
  // These four values will store the result until the end, when the result
  // will be moved into the "[]_arg" reference variables.
  // They will be initialized as if reset = TRUE.
  // Note that they are flipped, so rowStart begins at nRow, the value it should take
  // if all rows satisfy cond(ition).
  // rowEnd starts at 0, but if all rows satisfy cond, it will be nRow upon completion.
  size_t rowStart(nRow);  // first row that is not all constant
  size_t rowEnd(0);       // one past last row that is not all constant
  size_t colStart(nCol);  // first col that is not all constant
  size_t colEnd(0);       // one past last col that is not all constant
  // If not resetting, take input regions by flipping input values.
  if(!reset) {
    rowStart = rowEnd_arg;
    rowEnd = rowStart_arg; 
    colStart = colEnd_arg;
    colEnd = colStart_arg;
  } else {
    // Use the input variables as for loop bounds.
    // If !reset, the above clause achieves this, so it is only done if reset==true.
    rowStart_arg = rowEnd; // These are deliberately flipped.  Default 0
    rowEnd_arg = rowStart; // Default nRow
    colStart_arg = colEnd; // Default 0
    colEnd_arg = colStart; // Default nCol
  }

  auto cond2 = [&](size_t i, size_t j) {
    if(ignore_lower)
      return i < j ? true : cond(x(i,j));
    if(ignore_upper)
      return i > j ? true : cond(x(i,j));
    return cond(x(i,j));
  };
  
  // For each row, counting up from start, check if all columns satisfy cond
  for(size_t i = rowStart_arg; i < rowEnd_arg; ++i) {
    bool this_row_true(true);
    for(size_t j = colStart_arg; j < colEnd_arg; ++j) {
      if(!(this_row_true &= cond2(i, j ))) break;
    }
    if(!this_row_true) { // mark the first row that has an element failing cond
      rowStart = i;
      break;
    }
  }

  // If all rows satsify cond, there is no need to iterate again for rowEnd
  if(rowStart == rowEnd_arg) {
    rowEnd = rowEnd_arg;
  } else {
    // For each row, counting down from end, check if all columns satsify cond
    for(size_t i = rowEnd_arg; i > rowStart_arg; --i) {
      bool this_row_true(true);
      for(size_t j = colStart_arg; j < colEnd_arg; ++j) {
	if(!(this_row_true &= cond2(i-1, j ))) break;
      }
      if(!this_row_true) {
	rowEnd = i; // mark one past the highest row (first from counting down) that has an element failing cond
	break;
      }
    }
  }

  // For each col, counting up from start, check if all rows satsify cond
  for(size_t j = colStart_arg; j < colEnd_arg; ++j) {
    bool this_col_true(true);
    for(size_t i = rowStart_arg; i < rowEnd_arg; ++i) {
      if(!(this_col_true &= cond2(i, j ))) break;
    }
    if(!this_col_true) {
      colStart = j; // mark first col that has an element failing cond
      break;
    }
  }

  // If all cols satisfy cond, there is no need to iterate again for colEnd
  if(colStart == colEnd_arg) {
    colEnd = colEnd_arg;
  } else {
    // For each col, counting down from end, check if all rows satisfy cond
    for(size_t j = colEnd_arg; j > colStart_arg; --j) {
      bool this_col_true(true);
      for(size_t i = rowStart_arg; i < rowEnd_arg; ++i) {
	if(!(this_col_true &= cond2(i, j-1 ))) break;
      }
      if(!this_col_true) {
	colEnd = j ; // mark one past last highest col (first from counting down) that has an element failing cod
	break;
      }
    }
  }
  // Does the entire block satisfy cond?
  bool all_true = (rowStart == rowEnd_arg) && (rowEnd == rowEnd_arg) && (colStart == colEnd_arg) && (colEnd == colEnd_arg);
  rowStart_arg = rowStart;
  rowEnd_arg = rowEnd;
  colStart_arg = colStart;
  colEnd_arg = colEnd;

  return all_true;
}


#endif
