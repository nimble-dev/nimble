Table: (#tab:coreR) Basic R manipulation functions in NIMBLE. To find help in R for NIMBLE's version of a function, use the "nim" prefix and capitalize the next
letter. E.g. `help(nimC)` for help with `c()`.

  Function        Comments (differences from R)
  --------------- -----------------------------------
  `c()`           No `recursive` argument.
  `rep()`         No `rep.int` or `rep_len` arguments.
  `seq()` and ':' Negative integer sequences from ‘:’, e.g. , `2:1` do not work.
  `which()`       No `arr.ind` or `useNames` arguments.
  `diag()`        Works like R in three ways: `diag(vector)` returns a matrix with `vector` on the diagonal; 
                  `diag(matrix)` returns the diagonal vector of `matrix`; 
                  `diag(n)` returns an $n \times n$ identity matrix. No `nrow` or `ncol` arguments.
  `diag()<-`      Works for assigning the diagonal vector of a matrix.
  `dim()`         Works on a vector as well as higher-dimensional arguments.
  `length()`
  `is.na()`       Does not correctly handle NAs from R that are type `'logical'`, 
                  so convert these using `as.numeric()` before passing from R to NIMBLE.
  `is.nan()`
  `numeric()`     Allows additional arguments to control initialization.
  `logical()`     Allows additional arguments to control initialization.
  `integer()`     Allows additional arguments to control initialization.
  `matrix()`      Allows additional arguments to control initialization.
  `array()`       Allows additional arguments to control initialization.
  indexing        Arbitrary integer and logical indexing is supported for objects of one or two dimensions. 
                  For higher-dimensional objects, only `:` indexing works and then only to create an object
                  of at most two dimensions.


