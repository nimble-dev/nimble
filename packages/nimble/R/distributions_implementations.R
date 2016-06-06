# additional distributions provided by NIMBLE
# in general we use doubles on C side so convert parameters to doubles here before passing to C

#' The Wishart Distribution
#'
#' Density and random generation for the Wishart distribution, using the Cholesky factor of either the scale matrix or the rate matrix.
#'
#' @name Wishart
#' @aliases wishart
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
NULL

#' @rdname Wishart
#' @export
dwish_chol <- function(x, cholesky, df, scale_param = TRUE, log = FALSE) {
  # scale_param = TRUE is the GCSR parameterization (i.e., scale matrix); scale_param = FALSE is the BUGS parameterization (i.e., rate matrix)
  .Call('C_dwish_chol', as.double(x), as.double(cholesky), as.double(df), as.double(scale_param), as.logical(log))
}

#' @rdname Wishart
#' @export
rwish_chol <- function(n = 1, cholesky, df, scale_param = TRUE) {
    if(n != 1) warning('rwish_chol only handles n = 1 at the moment')
    matrix(.Call('C_rwish_chol', as.double(cholesky), as.double(df), as.double(scale_param)), nrow = sqrt(length(cholesky)))
}

#' The Dirichlet Distribution
#'
#' Density and random generation for the Dirichlet distribution
#'
#' @name Dirichlet
#' @aliases dirichlet
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
NULL

#' @rdname Dirichlet
#' @export
ddirch <- function(x, alpha, log = FALSE) {
    .Call('C_ddirch', as.double(x), as.double(alpha), as.logical(log))
}

#' @rdname Dirichlet
#' @export
rdirch <- function(n = 1, alpha) {
    if(n != 1) warning('rdirch only handles n = 1 at the moment')
  .Call('C_rdirch', as.double(alpha))
}

#' The Multinomial Distribution
#'
#' Density and random generation for the multinomial distribution
#'
#' @name Multinomial
#' @aliases multinomial
#' 
#' @param x vector of values.
#' @param n number of observations (only \code{n=1} is handled currently).
#' @param size number of trials.
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
NULL

#' @rdname Multinomial
#' @export
dmulti <- function(x, size = sum(x), prob, log = FALSE) {
  .Call('C_dmulti', as.double(x), as.double(size), as.double(prob), as.logical(log))
}

#' @rdname Multinomial
#' @export
rmulti <- function(n = 1, size, prob) {
  if(n != 1) warning('rmulti only handles n = 1 at the moment')
  .Call('C_rmulti', as.double(size), as.double(prob))
}

#' The Categorical Distribution
#'
#' Density and random generation for the categorical distribution
#'
#' @name Categorical
#' 
#' @param x non-negative integer-value numeric value.
#' @param n number of observations.
#' @param prob vector of probabilities, summing to one.
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
NULL

#' @rdname Categorical
#' @export
dcat <- function(x, prob, log = FALSE) {
  .Call('C_dcat', as.double(x), as.double(prob), as.logical(log))
}

#' @rdname Categorical
#' @export
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
#' @name t
#' 
#' @param x vector of values.
#' @param n number of observations.
#' @param df vector of degrees of freedom values.
#' @param p vector of probabilities.
#' @param q vector of quantiles.
#' @param mu vector of location values.
#' @param sigma vector of scale values.
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
NULL

#' @rdname t
#' @export
dt_nonstandard <- function(x, df = 1, mu = 0, sigma = 1, log = FALSE) {
  .Call('C_dt_nonstandard', as.double(x), as.double(df), as.double(mu), as.double(sigma), as.logical(log))
}

#' @rdname t
#' @export
rt_nonstandard <- function(n, df = 1, mu = 0, sigma = 1) {
  .Call('C_rt_nonstandard', as.integer(n), as.double(df), as.double(mu), as.double(sigma))
}

#' @rdname t
#' @export
pt_nonstandard <- function(q, df = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('C_pt_nonstandard', as.double(q), as.double(df), as.double(mu), as.double(sigma), as.logical(lower.tail), as.logical(log.p))
}

#' @rdname t
#' @export
qt_nonstandard <- function(p, df = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('C_qt_nonstandard', as.double(p), as.double(df), as.double(mu), as.double(sigma), as.logical(lower.tail), as.logical(log.p))
}

#' The Multivariate Normal Distribution
#'
#' Density and random generation for the multivariate normal distribution, using the Cholesky factor of either the precision matrix or the covariance matrix.
#'
#' @name MultivariateNormal
#' 
#' @param x vector of values.
#' @param n number of observations (only \code{n=1} is handled currently).
#' @param mean vector of values giving the mean of the distribution.
#' @param cholesky upper-triangular Cholesky factor of either the precision matrix (when \code{prec_param} is TRUE) or covariance matrix (otherwise).
#' @param prec_param logical; if TRUE the Cholesky factor is that of the precision matrix; otherwise, of the covariance matrix.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Christopher Paciorek
#' @export
#' @details See Gelman et al., Appendix A or the BUGS manual for mathematical details. The rate matrix as used here is defined as the inverse of the scale matrix, \eqn{S^{-1}}, given in Gelman et al. 
#' @return \code{dmnorm_chol} gives the density and \code{rmnorm_chol} generates random deviates.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
#' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' mean <- c(-10, 0, 10)
#' covmat <- matrix(c(1, .9, .3, .9, 1, -0.1, .3, -0.1, 1), 3)
#' ch <- chol(covmat)
#' x <- rmnorm_chol(1, mean, ch, prec_param = FALSE)
#' dmnorm_chol(x, mean, ch, prec_param = FALSE)
NULL

#' @rdname MultivariateNormal
#' @export
dmnorm_chol <- function(x, mean, cholesky, prec_param = TRUE, log = FALSE) {
  # cholesky should be upper triangular
  # FIXME: allow cholesky to be lower tri
  .Call('C_dmnorm_chol', as.double(x), as.double(mean), as.double(cholesky), as.double(prec_param), as.logical(log))
}

#' @rdname MultivariateNormal
#' @export
rmnorm_chol <- function(n = 1, mean, cholesky, prec_param = TRUE) {
 ## cholesky should be upper triangular
 ## FIXME: allow cholesky to be lower tri
    if(n != 1) warning('rmnorm_chol only handles n = 1 at the moment')
    .Call('C_rmnorm_chol', as.double(mean), as.double(cholesky), as.double(prec_param))
}

#' The Multivariate t Distribution
#'
#' Density and random generation for the multivariate t distribution, using the Cholesky factor of either the precision matrix (i.e., inverse scale matrix) or the scale matrix.
#'
#' @name Multivariate-t
#'
#' @aliases mvt multivariate-t 
#' 
#' @param x vector of values.
#' @param n number of observations (only \code{n=1} is handled currently).
#' @param mu vector of values giving the location of the distribution.
#' @param cholesky upper-triangular Cholesky factor of either the precision matrix (i.e., inverse scale matrix) (when \code{prec_param} is TRUE) or scale matrix (otherwise).
#' @param df degrees of freedom.
#' @param prec_param logical; if TRUE the Cholesky factor is that of the precision matrix; otherwise, of the scale matrix.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Peter Sujan
#' @export
#' @details See Gelman et al., Appendix A or the BUGS manual for mathematical details. The 'precision' matrix as used here is defined as the inverse of the scale matrix, \eqn{\Sigma^{-1}}, given in Gelman et al. 
#' @return \code{dmvt_chol} gives the density and \code{rmvt_chol} generates random deviates.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
#' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' mu <- c(-10, 0, 10)
#' scalemat <- matrix(c(1, .9, .3, .9, 1, -0.1, .3, -0.1, 1), 3)
#' ch <- chol(scalemat)
#' x <- rmvt_chol(1, mu, ch, df = 1, prec_param = FALSE)
#' dmvt_chol(x, mu, ch, df = 1, prec_param = FALSE)
#' 
NULL

#' @rdname Multivariate-t
#' @export
dmvt_chol <- function(x, mu, cholesky, df, prec_param = TRUE, log = FALSE) {
  # cholesky should be upper triangular
  # FIXME: allow cholesky to be lower tri
  .Call('C_dmvt_chol', as.double(x), as.double(mu), as.double(cholesky),
        as.double(df), as.double(prec_param), as.logical(log))
}

#' @rdname Multivariate-t
#' @export
rmvt_chol <- function(n = 1, mu, cholesky, df, prec_param = TRUE) {
  ## cholesky should be upper triangular
  ## FIXME: allow cholesky to be lower tri
  if(n != 1) warning('rmnorm_chol only handles n = 1 at the moment')
  .Call('C_rmvt_chol', as.double(mu), as.double(cholesky),
        as.double(df), as.double(prec_param))
}

#' Interval calculations 
#'
#' Calculations to handle censoring
#'
#' @name Interval
#' 
#' @param x vector of interval indices.
#' @param n number of observations.
#' @param t vector of values.
#' @param c vector of one or more values delineating the intervals.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Christopher Paciorek
#' @export
#' @details Used for working with censoring in BUGS code.
#' Taking \code{c} to define the endpoints of two or more intervals (with implicit endpoints of plus/minus infinity), \code{x} (or the return value of \code{rinterval}) gives the non-negative integer valued index of the interval in which \code{t} falls. See the NIMBLE manual for additional details. 
#' @return \code{dinterval} gives the density and \code{rinterval} generates random deviates,
#' but these are unusual as the density is 1 if \code{x} indicates the interval in which \code{t}
#' falls and 0 otherwise and the deviates are simply the interval(s) in which \code{t} falls.
##' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' endpoints <- c(-3, 0, 3)
#' vals <- c(-4, -1, 1, 5)
#' x <- rinterval(4, vals, endpoints)
#' dinterval(x, vals, endpoints)
#' dinterval(c(1, 5, 2, 3), vals, endpoints)
NULL

 
#' @rdname Interval
#' @export
dinterval <- function(x, t, c, log = FALSE) {
    .Call('C_dinterval', as.double(x), as.double(t), as.double(c), as.logical(log))
}

#' @rdname Interval
#' @export
rinterval <- function(n = 1, t, c) {
    .Call('C_rinterval', as.integer(n), as.double(t), as.double(c))
}

#' Constraint calculations in NIMBLE
#'
#' Calculations to handle censoring
#'
#' @name Constraint
#'
#' @param x value indicating whether \code{cond} is TRUE or FALSE
#' @param n number of observations (only \code{n=1} is handled currently).
#' @param cond logical value   
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Christopher Paciorek
#' @export
#' @details Used for working with constraints in BUGS code.
#' See the NIMBLE manual for additional details. 
#' @return \code{dconstraint} gives the density and \code{rconstraint} generates random deviates,
#' but these are unusual as the density is 1 if \code{x} matches \code{cond} and
#' 0 otherwise and the deviates are simply the value of \code{cond}
##' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' constr <- 3 > 2 && 4 > 0
#' x <- rconstraint(1, constr)
#' dconstraint(x, constr)
#' dconstraint(0, 3 > 4)
#' dconstraint(1, 3 > 4)
#' rconstraint(1, 3 > 4)
NULL

#' @rdname Constraint
#' @export
dconstraint <- function(x, cond, log = FALSE) {
    if(length(x) > 1 || length(cond) > 1) stop('dconstraint is not vectorized')
    if(is.na(x) || is.na(cond)) return(x + cond) # mimic how R's C functions handle NA and NaN inputs
    if(x == cond || x == 0) result <- 1 else result <- 0
    if(log) return(log(result)) else return(result)
}

#' @rdname Constraint
#' @export
rconstraint <- function(n = 1, cond) {
    if(n != 1 || length(cond) > 1) stop('rconstraint only handles n = 1 at the moment')
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
#'   generation for the exponential distribution with rate
#'     (i.e., mean of \code{1/rate}) or scale parameterizations.
#' 
#' @name Exponential
#' 
#' @param x vector of values.
#' @param n number of observations.
#' @param p vector of probabilities.
#' @param q vector of quantiles.
#' @param rate vector of rate values.
#' @param scale vector of scale values.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @param log.p logical; if TRUE, probabilities p are given by user as log(p).
#' @param lower.tail logical; if TRUE (default) probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#' @author Christopher Paciorek
#' @export
#' @details NIMBLE's exponential distribution functions use Rmath's functions
#' under the hood, but are parameterized to take both rate and scale and to
#' use 'rate' as the core parameterization in C, unlike Rmath, which uses 'scale'.
#' See Gelman et al., Appendix A or
#' the BUGS manual for mathematical details. 
#' @return \code{dexp_nimble} gives the density, \code{pexp_nimble} gives the distribution
#' function, \code{qexp_nimble} gives the quantile function, and \code{rexp_nimble}
#' generates random deviates.
#' @references Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
#' @seealso \link{Distributions} for other standard distributions
#' 
#' @examples
#' x <- rexp_nimble(50, scale = 3)
#' dexp_nimble(x, scale = 3)
NULL


#' @rdname Exponential
#' @export
dexp_nimble <- function(x, rate = 1/scale, scale = 1, log = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
    .Call('C_dexp_nimble', as.double(x), as.double(rate), as.logical(log))
}

#' @rdname Exponential
#' @export
rexp_nimble <- function(n = 1, rate = 1/scale, scale = 1) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
    .Call('C_rexp_nimble', as.integer(n), as.double(rate))
}

#' @rdname Exponential
#' @export
pexp_nimble <- function(q, rate = 1/scale, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
  .Call('C_pexp_nimble', as.double(q), as.double(rate), as.logical(lower.tail), as.logical(log.p))
}

#' @rdname Exponential
#' @export
qexp_nimble <- function(p, rate = 1/scale, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if (!missing(rate) && !missing(scale)) {
        if (abs(rate * scale - 1) < 1e-15) 
            warning("specify 'rate' or 'scale' but not both")
        else stop("specify 'rate' or 'scale' but not both")
    }
  .Call('C_qexp_nimble', as.double(p), as.double(rate), as.logical(lower.tail), as.logical(log.p))
}

