Table: (\#tab:functionsMatrix) Functions operating on vectors and matrices. Status column indicates if the function is currently
provided in NIMBLE.

  Usage                   Description                                Comments                         Status
  ----------------------- ------------------------------------------ -------------------------------- ------  
  `inverse(x)`            matrix inverse                             $x$ symmetric, positive def.     yes
  `chol(x)`               matrix Cholesky factorization              $x$ symmetric, positive def.     yes
                                                                     returns upper triangular matrix
  `t(x)`                  matrix transpose                           $x^\top$                         yes
  `x%*%y`                 matrix multiply                            $xy$; $x$, $y$ conformant        yes
  `inprod(x, y)`          dot product                                $x^\top y$; $x$ and $y$ vectors  yes
  `solve(x,y)`            solve system of equations                  $x^{-1} y$; $y$ matrix or vector yes
  `forwardsolve(x, y)`    solve lower-triangular system of equations $x^{-1} y$; $x$ lower-triangular yes
  `backsolve(x, y)`       solve upper-triangular system of equations $x^{-1} y$; $x$ upper-triangular yes
  `logdet(x)`             log matrix determinant                     $\log|x|$                        yes
  `asRow(x)`              convert vector to 1-row matrix             sometimes automatic              yes
  `asCol(x)`              convert vector to 1-column matrix          sometimes automatic              yes
  `sum(x)`                sum of elements of `x`                                                      yes
  `mean(x)`               mean of elements of `x`                                                     yes
  `sd(x)`                 standard deviation of elements of `x`                                       yes
  `prod(x)`               product of elements of `x`                                                  yes
  `min(x), max(x)`        min. (max.) of elements of `x`                                              yes
  `pmin(x, y), pmax(x,y)` vector of mins (maxs) of elements of                                        yes
                          `x` and `y`
  `interp.lin(x, v1, v2)` linear interpolation                                                        no
