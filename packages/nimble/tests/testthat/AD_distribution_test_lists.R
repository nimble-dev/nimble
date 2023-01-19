#########################################
## distribution function parameterization
#########################################

## create a base distribution specification which will
## be extended to specify the actual distributions
distn_base <- list(
  ##
  ## name of distn_param entry should be abbreviated distribution name
  ## as listed in Rdist field of distributions_inputList.R (e.g. dnorm -> norm)
  ##
  name = '',
  ##
  ## variants should create a valid distribution function name
  ## when prepended to name
  ##
  variants = c('d'),

  args = list(
    ##
    ## the rand_variate arg must be included
    ##
    rand_variate = list( ## TODO: better naming convention
      ##
      ## input_gen_fun is a function that should generate a number on the support,
      ## or return a function that uses the other args to generate random
      ## variates from this distribution. In the latter case, the formals of the
      ## returned function should be from among the other arg names. See the
      ## definition of binom for an example of how to use the parameters of the
      ## distribution.
      input_gen_fun = NULL,
      ##
      ## The type field is required for all args.
      ##
      type = NULL
    )
    ##
    ## log is a special arg that will trigger a test to check that
    ## the _logFixed variant is not incorrectly used when included.
    ## See distributions_R.hpp. First we test with a fixed log arg
    ## in distn_tests but later we'll add a log argument to each
    ## distn_params entry in distn_with_log_tests, so don't add it here.
    ##
    ## log = list(
    ##   input_gen_fun = function(n) sample(0:1, size = n, replace = TRUE),
    ##   type = 'double(0)'
    ## )
    ##
  ),
  ##
  ## You may want to pass in more_args to distribution function expr,
  ## but not as one of the nimbleFunction args, in which case they
  ## can be included as follows:
  ##
  ## more_args = list(
  ##   d = list(log = 1),
  ##   p = list(log.p = 1)
  ## ),
  ##
  ## wrt should be a character vector and can include arg names
  ## and any distribution function first argument names among the
  ## choices 'x', 'q', and 'p'.
  ##
  wrt = character(0)
)

## all distribution specifications will be added to this list and
## converted into an AD test parameterization by make_distribution_fun_AD_test()
distn_params <- list()

########################
## Binomial distribution
########################

distn_params[['binom_base']] <- distn_base
distn_params[['binom_base']]$name <- 'binom'
distn_params[['binom_base']]$args <- list(
  rand_variate = list(
    ##
    ## `size` here refers to the size parameter of the binomial distribution.
    ## Since the support depends on a parameter, return a function which
    ## test_AD will use with the realized value of the parameter `size`.
    ##
    input_gen_fun = function(n) function(size, prob)
        if(n >= max(length(size), length(prob))) rbinom(n, size, prob) else rbinom(n, min(size), prob),
    type = c('double(0)', 'double(1)'),
    size = c(1, 5)
  ),
  size = list(
    input_gen_fun = function(n) sample(5:10, size = n, replace = TRUE),
    type = c('double(0)', 'double(1)'),
    size = list(1, 3)
  ),
  prob = list(
    input_gen_fun = function(n) runif(n, .4, .6),
    type = c('double(0)', 'double(1)'),
    size = list(1, 7)
  )
)
distn_params[['binom_base']]$wrt <- c('prob') ## Why doesn't x get baked in? Maybe the atomic?

###########################
## Categorical distribution
###########################

# dcat really won't work due to the problem in CppAD with VecAD,
# which does not successfully get handled by base2ad, which is
# used in double-taping.
## distn_params[['cat_base']] <- distn_base
## distn_params[['cat_base']]$name <- 'cat'
## distn_params[['cat_base']]$args <- list(
##   rand_variate = list(
##     ## support depends on k = length(prob)
##     input_gen_fun = function(n) function(prob) sample(1:length(prob), size = n, replace = TRUE, prob = prob),
##     type = c('double(0)', 'double(1)'),
##     size = list(1, 1) # vector is not supported
##   ),
##   prob = list(
##     input_gen_fun = function(n) {prob <- runif(n); prob/sum(prob)},
##     type = c('double(1)', 'double(1)'),
##     size = list(5, 8)
##   )
## )
## distn_params[['cat_base']]$wrt <- c('prob')

###########################
## Multinomial distribution
###########################

## no size arg
## This is supported through compilation by inserting a default
## size = sum(x) before calling dmulti.
distn_params[['multi_no_size']] <- distn_base
distn_params[['multi_no_size']]$name <- 'multi'
distn_params[['multi_no_size']]$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) sample(1:10, n, replace = TRUE),
    type = c('double(1, 5)')
  ),
  prob = list(
    input_gen_fun = function(n) {prob <- runif(n); prob/sum(prob)},
    type = c('double(1, 5)')
  )
)
distn_params[['multi_no_size']]$wrt = c('prob')

## including size
distn_params[['multi_with_size']] <- distn_params[['multi_no_size']]
distn_params[['multi_with_size']]$args[['rand_variate']]$input_gen_fun <-
  function(n) function(size, prob) nimble::rmulti(n = 1, size, prob)
distn_params[['multi_with_size']]$args[['size']]  <- list(
  input_gen_fun = function(n) sample(10:100, 1, replace = TRUE),
  type = c('double(0)')
)

#################################
## Negative Binomial distribution
#################################

distn_params[['nbinom_base']] <- distn_base
distn_params[['nbinom_base']]$name <- 'nbinom'
distn_params[['nbinom_base']]$args <- list(
  # constrain range of parameters and draws to avoid numerical problems with uncompiled derivs
  rand_variate = list(
    input_gen_fun = function(n) function(size, prob) sample(10:15, size = n, replace = TRUE), #rnbinom(n, size, prob),
    type = c('double(0)', 'double(1, 5)')
  ),
  size = list(
    input_gen_fun = function(n) runif(n, 5, 15),
    type = c('double(0)', 'double(1, 3)')
  ),
  prob = list(
    input_gen_fun = function(n) runif(n, .4, .6),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['nbinom_base']]$wrt <- c('size', 'prob')

#######################
## Poisson distribution
#######################

distn_params[['pois_base']] <- distn_base
distn_params[['pois_base']]$name <- 'pois'
distn_params[['pois_base']]$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) function(lambda) rpois(n, lambda),
    type = c('double(0)', 'double(1, 5)')
  ),
  lambda = list(
    input_gen_fun = function(n) runif(n, max = 100),
    type = c('double(0)', 'double(1, 3)')
  )
)
distn_params[['pois_base']]$wrt <- c('lambda')

####################
## Beta distribution
####################

distn_params[['beta_base']] <- distn_base
distn_params[['beta_base']]$name <- 'beta'
distn_params[['beta_base']]$args <- list(
  # It can be hard to get all reasonable values with reasonable derivatives.
  # For this reason we choose shapes from 3 to 5 and x from .2 to .8 and
  # it should all be reasonable
  rand_variate = list(
    input_gen_fun = function(n) function(shape1, shape2) runif(n, .2, .8), #rbeta(n, shape1, shape2),
    type = c('double(0)', 'double(1, 5)')
  ),
  shape1 = list(
    input_gen_fun = function(n) runif(n, 3, 5),
    type = c('double(0)', 'double(1, 3)')
  ),
  shape2 = list(
    input_gen_fun = function(n) runif(n, 3, 5),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['beta_base']]$wrt <- c('shape1', 'shape2', 'x')

###########################
## Chi-squared distribution
###########################

distn_params[['chisq_base']] <- distn_base
distn_params[['chisq_base']]$name <- 'chisq'
distn_params[['chisq_base']]$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) function(df) rchisq(n, df),
    type = c('double(0)', 'double(1, 5)')
  ),
  df = list(
    input_gen_fun = function(n) runif(n, max = 100),
    type = c('double(0)', 'double(1, 3)')
  )
)
distn_params[['chisq_base']]$wrt <- c('x')

############################################
## Double Exponential (Laplace) distribution
############################################

## Here's an example of creating a distribution base
## specification and adding alternative parameterizations.

dexp_base <- distn_base
dexp_base$name <- 'dexp'
dexp_base$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) runif(n, min = -10, max = 10),
    type = c('double(0)', 'double(1, 5)')
  ),
  location = list(
    input_gen_fun = function(n) runif(n, min = -10, max = 10),
    type = c('double(0)', 'double(1, 3)')
  )
)
dexp_base$wrt <- c('location', 'x')

## scale parameterization
distn_params[['dexp_scale']] <- dexp_base
distn_params[['dexp_scale']]$args[['scale']] <- list(
  input_gen_fun = function(n) runif(n, max = 10),
  type = c('double(0)', 'double(1, 7)')
)
distn_params[['dexp_scale']]$wrt <- c(dexp_base$wrt, 'scale')

## rate parameterization
distn_params[['dexp_rate']] <- dexp_base
distn_params[['dexp_rate']]$args[['rate']] <-
    distn_params[['dexp_scale']]$args[['scale']]
distn_params[['dexp_rate']]$wrt <- c(dexp_base$wrt, 'rate')

###########################
## Exponential distribution
###########################

exp_base <- distn_base
exp_base$args = list(
  rand_variate = list(
    type = c('double(0)')
  )
)
exp_base$wrt <- c('x')

## base R version
distn_params[['exp_R']] <- exp_base
distn_params[['exp_R']]$name <- 'exp'
distn_params[['exp_R']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(rate) rexp(n, rate)
}
distn_params[['exp_R']]$args[['rate']] = list(
  input_gen_fun = function(n) runif(n, min = 3, max = 10),
  type = c('double(0)', 'double(1, 3)')
)
distn_params[['exp_R']]$wrt <- c(exp_base$wrt, 'rate')

## NIMBLE version, rate parameterization
distn_params[['exp_nimble_rate']] <- distn_params[['exp_R']]
distn_params[['exp_nimble_rate']]$name <- 'exp_nimble'

## NIMBLE version, scale parameterization
distn_params[['exp_nimble_scale']] <- exp_base
distn_params[['exp_nimble_scale']]$name <- 'exp_nimble'
distn_params[['exp_nimble_scale']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(scale) rexp(n, 1/scale)
}
distn_params[['exp_nimble_scale']]$args[['scale']]  <- list(
  input_gen_fun = function(n) runif(n, min = 3, max = 10),
  type = c('double(0)', 'double(1, 3)')
)
distn_params[['exp_nimble_scale']]$wrt <- c(exp_base$wrt, 'scale')

#####################
## Gamma distribution
#####################

distn_params[['gamma_scale']] <- distn_params[['exp_nimble_scale']]
distn_params[['gamma_rate']] <- distn_params[['exp_nimble_rate']]

## change the name
distn_params[['gamma_scale']]$name <-
  distn_params[['gamma_rate']]$name <- 'gamma'

## add the shape parameter
distn_params[['gamma_scale']]$args[['shape']] <-
  distn_params[['gamma_rate']]$args[['shape']] <- list(
    input_gen_fun = function(n) runif(n, min = 2, max = 4),
    type = c('double(0)', 'double(1, 7)')
  )

## add shape arg to wrt
distn_params[['gamma_scale']]$wrt <- c(
  distn_params[['gamma_scale']]$wrt, 'shape'
)
distn_params[['gamma_rate']]$wrt <- c(
  distn_params[['gamma_rate']]$wrt, 'shape'
)

## change how random variates are generated
distn_params[['gamma_rate']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(shape, rate) rgamma(n, shape, rate)
}
distn_params[['gamma_scale']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(shape, scale) rgamma(n, shape, scale=scale)
}

#############################
## Inverse Gamma distribution
#############################

distn_params[['invgamma_scale']] <- distn_params[['gamma_scale']]
distn_params[['invgamma_scale']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(shape, scale) rinvgamma(n, shape, scale=scale)
}

distn_params[['invgamma_rate']] <- distn_params[['gamma_rate']]
distn_params[['invgamma_rate']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(shape, rate) rinvgamma(n, shape, rate=rate)
}

distn_params[['invgamma_scale']]$name <-
  distn_params[['invgamma_rate']]$name <- 'invgamma'

##################################
## Sqrt Inverse Gamma distribution
##################################

distn_params[['sqrtinvgamma_scale']] <- distn_params[['gamma_scale']]
distn_params[['invgamma_scale']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(shape, scale) rsqrtinvgamma(n, shape, scale=scale)
}
distn_params[['sqrtinvgamma_rate']] <- distn_params[['gamma_rate']]
distn_params[['invgamma_rate']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(shape, rate) rsqrtinvgamma(n, shape, rate=rate)
}
distn_params[['sqrtinvgamma_scale']]$name <-
  distn_params[['sqrtinvgamma_rate']]$name <- 'sqrtinvgamma'

##########################
## Log-normal distribution
##########################

distn_params[['lnorm_base']] <- distn_base
distn_params[['lnorm_base']]$name <- 'lnorm'
distn_params[['lnorm_base']]$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) function(meanlog, sdlog) rlnorm(n, meanlog, sdlog),
    type = c('double(0)', 'double(1, 5)')
  ),
  meanlog = list(
    input_gen_fun = function(n) runif(n, min = -1.5, max = 1.5),
    type = c('double(0)', 'double(1, 3)')
  ),
  sdlog = list(
    input_gen_fun = function(n) runif(n, min = 0.5, max = 2),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['lnorm_base']]$wrt <- c('meanlog', 'sd', 'x')

########################
## Logistic distribution
########################

distn_params[['logis_base']] <- distn_params[['dexp_scale']]
distn_params[['logis_base']]$name <- 'logis'
distn_params[['logis_base']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(location, scale) rlogis(n, location, scale)
}

######################
## Normal distribution
######################

distn_params[['norm_base']] <- distn_base
distn_params[['norm_base']]$name <- 'norm'
distn_params[['norm_base']]$variants = c('d')
distn_params[['norm_base']]$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) function(mean, sd) rnorm(n, mean, sd),
    type = c('double(0)', 'double(1, 5)')
  ),
  mean = list(
    input_gen_fun = function(n) runif(n, min = -100, max = 100),
    type = c('double(0)', 'double(1, 3)')
  ),
  sd = list(
    input_gen_fun = function(n) runif(n, max = 10),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['norm_base']]$wrt <- c('mean', 'sd', 'x')

#################
## t-distribution
#################

## usual R version
distn_params[['t_R']] <- distn_params[['chisq_base']]
distn_params[['t_R']]$name <- 't'
distn_params[['t_R']]$args[['rand_variate']]$input_gen_fun <-
  function(n) function(df) rt(n, df)

## non-standard version
distn_params[['t_nonstandard']] <- distn_params[['t_R']]
distn_params[['t_nonstandard']]$name <- 't_nonstandard'

distn_params[['t_nonstandard']]$args[['mu']] <-
  distn_params[['norm_base']]$args[['mean']]

distn_params[['t_nonstandard']]$args[['sigma']] <-
  distn_params[['norm_base']]$args[['sd']]

distn_params[['t_nonstandard']]$wrt <- c(
  distn_params[['t_nonstandard']]$wrt, 'mu', 'sigma'
)

distn_params[['t_nonstandard']]$args[['rand_variate']]$input_gen_fun <- function(n) {
  function(mu, sigma) rnorm(n, mu, sigma)
}

#######################
## Uniform distribution
#######################

distn_params[['unif_base']] <- distn_base
distn_params[['unif_base']]$name <- 'unif'
## For uniform, we ensure that all variates are within all ranges created by recycling-rule combinations
## NB: Uncompiled derivatives will be NaN for -Inf values, but compiled derivatives will have nonsense values.
distn_params[['unif_base']]$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) runif(n, -9.9, 9.9),
    type = c('double(0)', 'double(1, 5)')
  ),
  min = list(
    input_gen_fun = function(n) runif(n, min = -100, max = -10),
    type = c('double(0)', 'double(1, 3)')
  ),
  max = list(
    input_gen_fun = function(n) runif(n, min = 10, max = 100),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['unif_base']]$wrt <- c('min', 'max', 'x')

#######################
## Weibull distribution
#######################

distn_params[['weibull_base']] <- distn_base   ## distn_params[['gamma_scale']]
distn_params[['weibull_base']]$name <- 'weibull'
distn_params[['weibull_base']]$args <- list(
  rand_variate = list(
    input_gen_fun = function(n) {function(shape, scale) rweibull(n, shape, scale)},
    type = c('double(0)', 'double(1, 5)')
  ),
  shape = list(
    input_gen_fun = function(n) runif(n, 2, 5),
    type = c('double(0)', 'double(1, 3)')
  ),
  scale = list(
    input_gen_fun = function(n) runif(n, 2, 5),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['weibull_base']]$wrt <- c('shape', 'scale', 'x')


#########################
## Dirichlet distribution
#########################

## Dirichlet is not easy to test because the x values must sum to 1.
## The "live on the simplex."
## See test-ADfunctions for a Dirichlet-specific test.
##
## distn_params[['dirch_base']] <- distn_base
## distn_params[['dirch_base']]$name <- 'dirch'
## distn_params[['dirch_base']]$args <- list(
##   rand_variate = list(
##     input_gen_fun = function(n) {prob <- runif(n); prob/sum(prob)},
##     type = c('double(1, 5)')
##   ),
##   alpha = list(
##     input_gen_fun = function(n) runif(n, max = 3), ## large values might not work well
##     type = c('double(1, 5)')
##   )
## )
## distn_params[['dirch_base']]$wrt = c('alpha', 'x')

###################################
## Multivariate Normal distribution
###################################

chol_base <- distn_base
chol_base$args <- list(
  cholesky = list(
    input_gen_fun = gen_pos_def_matrix,
    type = c('double(2, c(5, 5))')
  )
)

distn_params[['mnorm_chol_base']] <- chol_base
distn_params[['mnorm_chol_base']]$name <- 'mnorm_chol'
distn_params[['mnorm_chol_base']]$args[['rand_variate']] <- list(
  input_gen_fun = function(n) function(mean, cholesky)
    rmnorm_chol(n = 1, mean, cholesky),
  type = c('double(1, 5)')
)
distn_params[['mnorm_chol_base']]$args[['mean']] <- list(
  input_gen_fun = function(n) runif(n, min = -10, max = 10),
  type = c('double(1, 5)')
)
distn_params[['mnorm_chol_base']]$wrt <- c('mean', 'cholesky', 'x')

##############################
## Multivariate t-distribution
##############################

distn_params[['mvt_chol_base']] <- chol_base
distn_params[['mvt_chol_base']]$name <- 'mvt_chol'
distn_params[['mvt_chol_base']]$args[['rand_variate']] <- list(
  input_gen_fun = function(n) function(mu, cholesky, df)
    rmvt_chol(n = 1, mu, cholesky, df),
  type = c('double(1, 5)')
)
distn_params[['mvt_chol_base']]$args[['mu']] <- list(
  input_gen_fun = function(n) runif(n, min = -10, max = 10),
  type = c('double(1, 5)')
)
distn_params[['mvt_chol_base']]$args[['df']] <- list(
  input_gen_fun = function(n) sample(6:9, size = n, replace = TRUE),
  type = c('double(0)')
)
distn_params[['mvt_chol_base']]$wrt <- c('mu', 'cholesky', 'x')

#######################
## Wishart distribution
#######################

## The make_AD_test2 machinery is buggy for matrix arguments, so
## this will be a one-off test in test-ADfunctions.
## distn_params[['wish_chol_base']] <- chol_base
## distn_params[['wish_chol_base']]$name <- 'wish_chol'
## distn_params[['wish_chol_base']]$args[['rand_variate']] <- list(
##   input_gen_fun = function(n) function(cholesky, df)
##     rwish_chol(n = 1, cholesky, df),
##   type = c('double(2, c(5, 5))')
## )
## distn_params[['wish_chol_base']]$args[['df']] <- list(
##   input_gen_fun = function(n) runif(n, max = 100),
##   type = c('double(0)')
## )
## distn_params[['wish_chol_base']]$wrt <- c('cholesky', 'x')

##########################################################
## add a fixed log = 1 to the distribution fun expressions
##########################################################

## add more_args to include as part of distribution function expr
distn_params_log_1 <- lapply(distn_params, function(param) {
  more_args = lapply(
    param$variants, switch,
    d = list(log = 1),
    p = list(log.p = 1),
    q = list(log.p = 1)
  )
  names(more_args) <- param$variants
  param$more_args = more_args
  param
})

#########################################################
## create distribution function tests, with fixed log = 1
#########################################################
distn_tests <- unlist(
  lapply(distn_params_log_1, make_distribution_fun_AD_test),
  recursive = FALSE
)

set.seed(456)
distn_tests2 <- unlist(
  lapply(distn_params_log_1, make_distribution_fun_AD_test, maker = make_AD_test2),
  recursive = FALSE
)

#######################################################################
## create another set of distribution functions tests, variable log arg
#######################################################################

# TO-DO: Create special cases with 0s, 1s, maxs, mins, etc.
# TO-DO: Check wrt discrete params or rv
# binom: 

## add the arg 'log' to all the distn_params
distn_params_with_log <- lapply(distn_params, function(param) {
  param$args$log <- list(
    input_gen_fun = function(n) sample(0:1, size = n, replace = TRUE),
    type = 'double(0)'
  )
  param
})

distn_with_log_tests <- unlist(
  lapply(distn_params_with_log, make_distribution_fun_AD_test),
  recursive = FALSE
)

set.seed(567)
distn_with_log_tests2 <- unlist(
  lapply(distn_params_with_log, make_distribution_fun_AD_test, maker = make_AD_test2),
  recursive = FALSE
)


# Things to watch out for:
# in meta-taping, leaving elements out of wrt means they can get baked in: good and bad
