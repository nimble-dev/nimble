# additional distributions provided by NIMBLE
# in general we use doubles on C side so convert parameters to doubles here before passing to C

#' The Wishart Distribution
#'
#' Density and random generation for the Wishart distribution, using either the Cholesky factor of the scale matrix or the rate matrix.
#'
#' @aliases rwish_chol
#' 
#' @param x vector of values.
#' @param n number of observations (only \code{n=1} is handled currently).
#' @param cholesky upper-triangular Cholesky factor of either the scale matrix (when \code{scale_param} is TRUE) or rate matrix (otherwise).
#' @param df degrees of freedom.
#' @param scale_param logical; if TRUE the Cholesky factor is that of the scale matrix; otherwise, of the rate matrix.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Christopher Paciorek
#' @export
#' @details See Gelman et al., Appendix A or the BUGS manual for mathematical details. The rate matrix as used here is defined as the inverse of the scale matrix, \eqn{S^{-1}}, given in Gelman et al. 
#' @return \code{dwish_chol} gives the density and \code{rwish_chol} generates random deviates.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
#' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' df <- 40
#' ch <- chol(matrix(c(1, .7, .7, 1), 2))
#' x <- rwish_chol(1, ch, df = df)
#' dwish_chol(x, ch, df = df)
#' 
dwish_chol <- function(x, cholesky, df, scale_param = TRUE, log = FALSE) {
  # scale_param = TRUE is the GCSR parameterization (i.e., scale matrix); scale_param = FALSE is the BUGS parameterization (i.e., rate matrix)
  .Call('C_dwish_chol', as.double(x), as.double(cholesky), as.double(df), as.double(scale_param), as.logical(log))
}

rwish_chol <- function(n = 1, cholesky, df, scale_param = TRUE) {
    if(n != 1) warning('rwish_chol only handles n = 1 at the moment')
    matrix(.Call('C_rwish_chol', as.double(cholesky), as.double(df), as.double(scale_param)), nrow = sqrt(length(cholesky)))
}

#' The Dirichlet Distribution
#'
#' Density and random generation for the Dirichlet distribution
#'
#' @aliases rdirch
#' 
#' @param x vector of values.
#' @param n number of observations (only \code{n=1} is handled currently).
#' @param alpha vector of parameters of same length as \code{x}
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Christopher Paciorek
#' @export
#' @details See Gelman et al., Appendix A or the BUGS manual for mathematical details. 
#' @return \code{ddirch} gives the density and \code{rdirch} generates random deviates.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
#' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' alpha <- c(1, 10, 30)
#' x <- rdirch(1, alpha)
#' ddirch(x, alpha)
#' 
ddirch <- function(x, alpha, log = FALSE) {
    .Call('C_ddirch', as.double(x), as.double(alpha), as.logical(log))
}

rdirch <- function(n = 1, alpha) {
    if(n != 1) warning('rdirch only handles n = 1 at the moment')
  .Call('C_rdirch', as.double(alpha))
}

#' The Multinomial Distribution
#'
#' Density and random generation for the Multinomial distribution
#'
#' @aliases rdirch
#' 
#' @param x vector of values.
#' @param n number of observations (only \code{n=1} is handled currently).
#' @param prob vector of probabilities, summing to one, of same length as \code{x}
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Christopher Paciorek
#' @export
#' @details See Gelman et al., Appendix A or the BUGS manual for mathematical details. 
#' @return \code{dmulti} gives the density and \code{rmulti} generates random deviates.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
#' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' size <- 30
#' probs <- c(1/4, 1/10, 1 - 1/4 - 1/10)
#' x <- rmulti(1, size, probs)
#' dmulti(x, size, probs)
#' 
dmulti <- function(x, size = sum(x), prob, log = FALSE) {
  .Call('C_dmulti', as.double(x), as.double(size), as.double(prob), as.logical(log))
}

rmulti <- function(n = 1, size, prob) {
  if(n != 1) warning('rmulti only handles n = 1 at the moment')
  .Call('C_rmulti', as.double(size), as.double(prob))
}

#' The Categorical Distribution
#'
#' Density and random generation for the Categorical distribution
#'
#' @aliases rdirch
#' 
#' @param x non-negative integer-value numeric value
#' @param n number of observations
#' @param prob vector of probabilities, summing to one
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Christopher Paciorek
#' @export
#' @details See the BUGS manual for mathematical details. 
#' @return \code{dcat} gives the density and \code{rcat} generates random deviates.
##' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' probs <- c(1/4, 1/10, 1 - 1/4 - 1/10)
#' x <- rcat(n = 30, probs)
#' dcat(x, probs)
#' 
dcat <- function(x, prob, log = FALSE) {
  .Call('C_dcat', as.double(x), as.double(prob), as.logical(log))
}

rcat <- function(n = 1, prob) {
  .Call('C_rcat', as.integer(n), as.double(prob))
}


#' The t Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the t distribution with \code{df} degrees of freedom,
#' allowing non-zero location, \code{mu},
#' and non-unit scale, \code{sigma} 
#'
#' @aliases rt_nonstandard, qt_nonstandard, pt_nonstandard
#' 
#' @param x vector of values.
#' @param n number of observations 
#' @param df vector of degrees of freedom values
#' @param mu vector of location values
#' @param sigma vector of scale values
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @param log.p logical; if TRUE, probabilities p are given by user as log(p).
#' @param lower.tail logical; if TRUE (default) probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#' @author Christopher Paciorek
#' @export
#' @details See Gelman et al., Appendix A or the BUGS manual for mathematical details. 
#' @return \code{dt_nonstandard} gives the density, \code{pt_nonstandard} gives the distribution
#' function, \code{qt_nonstandard} gives the quantile function, and \code{rt_nonstandard}
#' generates random deviates.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
#' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' x <- rt_nonstandard(50, df = 1, mu = 5, sigma = 1)
#' dt_nonstandard(x, 3, 5, 1)
#' 
dt_nonstandard <- function(x, df = 1, mu = 0, sigma = 1, log = FALSE) {
  .Call('C_dt_nonstandard', as.double(x), as.double(df), as.double(mu), as.double(sigma), as.logical(log))
}

rt_nonstandard <- function(n, df = 1, mu = 0, sigma = 1) {
  .Call('C_rt_nonstandard', as.integer(n), as.double(df), as.double(mu), as.double(sigma))
}

pt_nonstandard <- function(x, df = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('C_pt_nonstandard', as.double(x), as.double(df), as.double(mu), as.double(sigma), as.logical(lower.tail), as.logical(log.p))
}

qt_nonstandard <- function(p, df = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('C_qt_nonstandard', as.double(p), as.double(df), as.double(mu), as.double(sigma), as.logical(lower.tail), as.logical(log.p))
}

dmnorm_chol <- function(x, mean, cholesky, prec_param = TRUE, log = FALSE) {
  # cholesky should be upper triangular
  # FIXME: allow cholesky to be lower tri
  .Call('C_dmnorm_chol', as.double(x), as.double(mean), as.double(cholesky), as.double(prec_param), as.logical(log))
}

rmnorm_chol <- function(n = 1, mean, cholesky, prec_param = TRUE) {
 ## cholesky should be upper triangular
 ## FIXME: allow cholesky to be lower tri
    if(n != 1) warning('rmnorm_chol only handles n = 1 at the moment')
    .Call('C_rmnorm_chol', as.double(mean), as.double(cholesky), as.double(prec_param))
}

dinterval <- function(x, t, c, log = FALSE) {
    .Call('C_dinterval', as.double(x), as.double(t), as.double(c), as.logical(log))
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
    return(as.integer(cond))
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
    .Call('C_rexp_nimble', as.integer(n), as.double(rate))
}

pexp_nimble <- function(q, rate = 1/scale, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
  .Call('C_pexp_nimble', as.double(q), as.double(rate), as.logical(lower.tail), as.logical(log.p))
}

qexp_nimble <- function(p, rate = 1/scale, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
  .Call('C_qexp_nimble', as.double(p), as.double(rate), as.logical(lower.tail), as.logical(log.p))
}

