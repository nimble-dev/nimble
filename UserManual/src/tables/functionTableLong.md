Table: (#tab:functions) Functions operating on scalars, many of which can operate on each element (component-wise) of
vectors and matrices. Status column indicates if the function is currently provided in NIMBLE. Vector input column
indicates if the function can take a vector as an argument (i.e., if the function is vectorized).

  Usage                  Description                     Comments                     Status  Vector input
  ---------------------- ------------------------------- ---------------------------- ------  ------------
  `x | y`, `x & y`       logical OR ($|$) and AND(&)                                  yes     yes 
  `!x`                   logical not                                                  yes     yes
  `x > y, x >= y`        greater than (and or equal to)                               yes     yes 
  `x < y, x <= y`        less than (and or equal to)                                  yes     yes 
  `x != y, x == y`       (not) equals                                                 yes     yes
  `x + y, x - y, x * y`  component-wise operators        mix of scalar and vector     yes     yes
  `x / y`                component-wise division         vector $x$ and scalar $y$    yes     yes
  `x^y, pow(x, y)`       power                           $x^y$; vector $x$,scalar $y$ yes     yes
  `x %% y`               modulo (remainder)                                           yes     no
  `min(x1, x2),`         min. (max.) of two scalars                                   yes     See `pmin`, 
    `max(x1, x2)`                                                                                 `pmax`
  `exp(x)`               exponential                                                  yes     yes
  `log(x)`               natural logarithm                                            yes  	  yes
  `sqrt(x)`              square root                                                  yes     yes
  `abs(x)`               absolute value                                               yes  	  yes
  `step(x)`              step function at 0              0 if $x<0$, 1 if $x>=0$      yes     yes
  `equals(x, y)`         equality of two scalars         1 if $x==y$, 0 if $x != y$   yes	  yes
  `cube(x)`              third power                     $x^3$                        yes	  yes
  `sin(x), cos(x),`      trigonometric functions                                      yes     yes
    `tan(x)` 
  `asin(x), acos(x),`    inverse trigonometric functions                              yes     yes
    `atan(x)`
  `asinh(x), acosh(x),`  inv. hyperbolic trig. functions                              yes     yes
    `atanh(x)`
  `logit(x)`             logit                           $\log(x/(1-x))$              yes     yes
  `ilogit(x), expit(x)`  inverse logit                   $\exp(x)/(1 + \exp(x))$      yes     yes
  `probit(x)`            probit (Gaussian quantile)      $\Phi^{-1}(x)$               yes     yes
  `iprobit(x), phi(x)`   inverse probit (Gaussian CDF)   $\Phi(x)$                    yes     yes
  `cloglog(x)`           complementary log log           $\log(-\log(1-x))$           yes     yes
  `icloglog(x)`          inverse complementary log log   $1 - \exp(-\exp(x))$         yes     yes
  `ceiling(x)`           ceiling function                $\lceil(x)\rceil$            yes     yes
  `floor(x)`             floor function                  $\lfloor(x)\rfloor$          yes     yes
  `round(x)`             round to integer                                             yes     yes
  `trunc(x)`             truncation to integer                                        yes     yes
  `lgamma(x), loggam(x)` log gamma function              $\log \Gamma(x)$             yes     yes
  `besselK(k, nu,`       modified bessel function        returns vector even if       yes     yes
  `...expon.scaled)`      of the second kind             $k$ a matrix/array
  `log1p(x)`             log of 1 + x                    $\log(1+x)$                  yes     yes
  `lfactorial(x),`       log factorial                   $\log x!$                    yes     yes
    `logfact(x)`
  `qDIST(x, PARAMS)`     “q” distribution functions      canonical parameterization   yes     yes
  `pDIST(x, PARAMS)`     “p” distribution functions      canonical parameterization   yes     yes
  `rDIST(x, PARAMS)`     “r” distribution functions      canonical parameterization   yes     yes 
  `dDIST(x, PARAMS)`     “d” distribution functions      canonical parameterization   yes     yes
  `sort(x)`                                                                           no
  `rank(x, s)`                                                                        no
  `ranked(x, s)`                                                                      no
  `order(x)`                                                                          no

