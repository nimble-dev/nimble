# additional distributions provided by NIMBLE
# FIXME: these should be modified to go directly to C w/o type conversion, with error-checking in C

dwish_chol <- function(x, cholesky, df, scale_param = TRUE, log = FALSE) {
  # scale_param = TRUE is the GCSR parameterization (i.e., scale matrix); scale_param = FALSE is the BUGS parameterization (i.e., rate matrix)
  .Call('C_dwish_chol', as.double(x), as.double(cholesky), as.double(df), as.double(scale_param), as.integer(log))
}

rwish_chol <- function(n = 1, cholesky, df, scale_param = TRUE) {
    if(n != 1) warning('rwish_chol only handles n = 1 at the moment')
    matrix(.Call('C_rwish_chol', as.double(cholesky), as.double(df), as.double(scale_param)), nrow = sqrt(length(cholesky)))
}

ddirch <- function(x, alpha, log = FALSE) {
    .Call('C_ddirch', as.double(x), as.double(alpha), as.logical(log))
}

rdirch <- function(n = 1, alpha) {
    if(n != 1) warning('rdirch only handles n = 1 at the moment')
  .Call('C_rdirch', as.double(alpha))
}

dmulti <- function(x, size = sum(x), prob, log = FALSE) {
  .Call('C_dmulti', as.double(x), as.double(size), as.double(prob), as.logical(log))
}

rmulti <- function(n = 1, size, prob) {
  if(n != 1) warning('rmulti only handles n = 1 at the moment')
  .Call('C_rmulti', as.double(size), as.double(prob))
}

dcat <- function(x, prob, log = FALSE) {
  .Call('C_dcat', as.double(x), as.double(prob), as.logical(log))
}

rcat <- function(n = 1, prob) {
  .Call('C_rcat', as.integer(n), as.double(prob))
}

dt_nonstandard <- function(x, df = 1, mu = 0, sigma = 1, log = FALSE) {
  .Call('C_dt_nonstandard', as.double(x), as.double(df), as.double(mu), as.double(sigma), as.logical(log))
}

rt_nonstandard <- function(n, df = 1, mu = 0, sigma = 1) {
  .Call('C_rt_nonstandard', as.integer(n), as.double(df), as.double(mu), as.double(sigma))
}

dmnorm_chol <- function(x, mean, cholesky, prec_param = TRUE, log = FALSE) {
  # cholesky should be upper triangular
  # FIXME: allow cholesky to be lower tri
  .Call('C_dmnorm_chol', as.double(x), as.double(mean), as.double(cholesky), as.logical(prec_param), as.logical(log))
}

rmnorm_chol <- function(n = 1, mean, cholesky, prec_param = TRUE) {
 ## cholesky should be upper triangular
 ## FIXME: allow cholesky to be lower tri
    if(n != 1) warning('rmnorm_chol only handles n = 1 at the moment')
    .Call('C_rmnorm_chol', as.double(mean), as.double(cholesky), as.logical(prec_param))
}

dinterval <- function(x, t, c, log = FALSE) {
    .Call('C_dinterval', as.integer(x), as.double(t), as.double(c), as.logical(log))
}

rinterval <- function(n = 1, t, c) {
    .Call('C_rinterval', as.integer(n), as.double(t), as.double(c))
}

dconstraint <- function(x, cond, log = FALSE) {
    if(is.na(x) || is.na(cond)) return(x + cond) # mimic how R's C functions handle NA and NaN inputs
    if(x == cond || x == 0) result <- 1 else result <- 0
    if(log) return(log(result)) else return(result)
}

rconstraint <- function(n = 1, cond) {
    if(n != 1) stop('rconstraint only handles n = 1 at the moment')
    if(is.na(cond)) {
        warning("NAs produced")
        return(NaN)
    }
    return(cond)
}

# exp_nimble extends R to allow rate or scale and provide common interface via 'rate' to C functions

#' The Exponential Distribution
#'
#'   Density, distribution function, quantile function and random
#'   generation for the exponential distribution with rate 'rate'
#'     (i.e., mean ‘1/rate’) or 'scale' parameterizations.
#' 
#' @aliases rexp_nimble, qexp_nimble, pexp_nimble
#' @author Christopher Paciorek
#' @details NIMBLE's exponential distribution functions use Rmath's functions
#' under the hood, but are parameterized to take both rate and scale and to
#' use 'rate' as the core parameterization in C.

dexp_nimble <- function(x, rate = 1/scale, scale = 1, log = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
    .Call('C_dexp_nimble', as.double(x), as.double(rate), as.logical(log))
}

rexp_nimble <- function(n = 1, rate = 1/scale, scale = 1) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
    .Call('C_rexp_nimble', as.integer(n), as.double(rate), as.logical(log))
}

pexp_nimble <- function(x, rate = 1/scale, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
  .Call('C_pexp_nimble', as.double(x), as.double(rate), as.logical(lower.tail), as.logical(log.p))
}

qexp_nimble <- function(p, rate = 1/scale, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
  .Call('C_qexp_nimble', as.double(p), as.double(rate), as.logical(lower.tail), as.logical(log.p))
}

