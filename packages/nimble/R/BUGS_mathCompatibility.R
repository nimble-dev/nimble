## "Logical" (sensu BUGS) functions for compatibility with BUGS

linkInverses <- list('logit' = quote(expit()),
                     'log' = quote(exp()),
                     'probit' = quote(phi()),
                     'cloglog' = quote(icloglog()) )

expit <- function(x) 1 / (1 + exp(-x))
probit <- function(x) qnorm(x)
ilogit <- function(x) 1 / (1 + exp(-x)) ## same as expit
iprobit <- function(x) pnorm(x)
icloglog <- function(x) 1-exp(-exp(x))
cloglog <- function(x) log(-log(1-x))
nimEquals <- function(x1, x2) if(x1 == x2) 1 else 0 ## "equals" conflicts with a usage in testthat
logfact <- function(x) lfactorial(x)
loggam <- function(x) lgamma(x)
logit <- function(x) log(x / (1-x))
phi <- function(x) pnorm(x) ## same as iprobit
pow <- function(x1, x2) x1^x2
nimStep <- function(x) ifelse(x >= 0, 1, 0) ## We rename step to nimStep before execution to avoid masking the R step function
inverse <- function(x) solve(x)
cube <- function(x) x^3
inprod <- function(v1, v2) sum(v1 * v2)
logdet <- function(m) {
  # this ensures that R returns NaN for negative determinants, to mimic C
  out <- determinant(m, logarithm = TRUE)
  if(out$sign == 1) return(out$modulus) else return(NaN)
}
