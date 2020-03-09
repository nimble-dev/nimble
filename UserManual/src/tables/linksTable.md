Table: (\#tab:links) Link functions.

  Link function     Description           Range       Inverse
  ----------------- --------------------- ----------- ------------------
  `cloglog(y) <- x` Complementary log log $0 < y < 1$ `y <- icloglog(x)`
  `log(y) <- x`     Log                   $0 < y$     `y <- exp(x)` 
  `logit(y) <- x`   Logit                 $0 < y < 1$ `y <- expit(x)` 
  `probit(y) <- x`  Probit                $0 < y < 1$ `y <- iprobit(x)`
