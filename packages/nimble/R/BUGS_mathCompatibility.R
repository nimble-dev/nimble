## "Logical" (sensu BUGS) functions for compatibility with BUGS

linkInverses <- list('logit' = quote(expit()),
                     'log' = quote(exp()),
                     'probit' = quote(phi()),
                     'cloglog' = quote(icloglog()) )
#' @export
expit <- function(x) 1 / (1 + exp(-x))
#' @export
probit <- function(x) qnorm(x)
#' @export
ilogit <- function(x) 1 / (1 + exp(-x)) ## same as expit
#' @export
iprobit <- function(x) pnorm(x)
#' @export
icloglog <- function(x) 1-exp(-exp(x))
#' @export
cloglog <- function(x) log(-log(1-x))
#' @export
nimEquals <- function(x1, x2) if(x1 == x2) 1 else 0 ## "equals" conflicts with a usage in testthat
#' @export
logfact <- function(x) lfactorial(x) ## lgamma(x+1) (from R) not equivalent to lgamma1p(x) (from our C). lgamma seems to be endowed with setting numerically 0 values to zero.
#' @export
loggam <- function(x) lgamma(x)
#' @export
logit <- function(x) log(x / (1-x))
#' @export
phi <- function(x) pnorm(x) ## same as iprobit
#' @export
pow <- function(x1, x2) x1^x2
#' @export
nimStep <- function(x) ifelse(x >= 0, 1L, 0L) ## We rename step to nimStep before execution to avoid masking the R step function
#' @export
inverse <- function(x) solve(x)
#' @export
cube <- function(x) x*x*x ## Doing is this way instead of x^3 makes it numerically identical to eigen's cube
#' @export
inprod <- function(v1, v2) sum(v1 * v2)
#' @export
logdet <- function(m) {
  # this ensures that R returns NaN for negative determinants, to mimic C
  out <- determinant(m, logarithm = TRUE)
  if(out$sign == 1) return(out$modulus) else return(NaN)
}
