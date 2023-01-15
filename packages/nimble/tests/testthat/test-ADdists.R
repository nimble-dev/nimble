source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

## I'll use more comprehensive tests in AD_distribution_test_lists
## and corner-case or special-case tests here.

## Status and to-do
## dcat: broken in double taping
## dbinom: looks ok.  add test case p = 0 or 1 and x = 0 or n
## dmulti with no size: add test case any p = 0 or x = 0 or n
## dmulti with size:    ditto
## dnbinom: looks ok.  add test for p = 1 and/or x = 0
## dpois: add test case of x = 0 and/or mu = 0
## dbeta: will break for x == 0 or x == 1 and the values of a and/or b where this is allowed by the dist.
## chisq: will break for x = 0
## dexp: check on breaking cases
## exp: will break for lambda = 0 (works in R)
## gamma: check for x = 0 case
## invgamma: ditto
## sqrtinvgamma: ditto
## lnorm: good
## logis: good
## norm: good
## t: good
## t_nonstandard: good
## unif: good. check on x at boundaries
## weibull: good. check on x = 0
## 
##
## dmnorm, dmvt: good but add test with prec_param flipped
## wishart: good - does it need more variants?
##
## TO-DO:
## Add warning message if prec_param or scale_param is not CppAD::Constant
## Add test for dirch
## Add test for invwishard
## Add test for lkj_chol
## Are there others that are needed now?


normDistTests <- list(
  make_AD_test2('dnorm', c(x='double()', mean='double()', sd='double()'),
                input_gen_funs = list(x=rnorm, mu=rnorm, sd=runif)),
  make_AD_test2('dnorm', c(x='double()', mean='double()', sd='double()'),
                more_args = list(log = 1),
                input_gen_funs = list(x=rnorm, mu=rnorm, sd=runif)),
  make_AD_test2('dnorm', c(x='double()', mean='double()', sd='double()', log='double()'),
                wrt = c('x', 'mean', 'sd'),
                inputs = list(record = c(x = 1.2, mean = 0.8, sd = 0.3, log = 1),
                              test   = c(x = 1.1, mean = 0.9, sd = 0.5, log = 0)))
)
result <- lapply(normDistTests, test_AD2, verbose = FALSE)

grunif <- function(min, max) function(size) runif(size, min, max)

betaDistTests <- list(
  make_AD_test2('dbeta', c(x='double()', shape1='double()', shape2='double()'),
                input_gen_funs = list(x='rbeta(size, 1, 1)', shape1='runif(size, 0.9, 1.1)', shape2='runif(size, 0.9, 1.1)')),
  make_AD_test2('dbeta', c(x='double()', shape1='double()', shape2='double()'),
                more_args = list(log = 1),
                input_gen_funs = list(x='rbeta(size, 1, 1)', shape1='runif(size, 0.9, 1.1)', shape2='runif(size, 0.9, 1.1)')),
  make_AD_test2('dbeta', c(x='double()', shape1='double()', shape2='double()', log='double()'),
                wrt = c('x', 'shape1', 'shape2'),
                inputs = list(record = c(x = 0.6, shape1 = 0.8, shape2 = 1.3, log = 1),
                              test   = c(x = 0.43, shape1 = 0.9, shape2 = 1.6, log = 0)))
)
# PROBLEM: Set tols
# WILL NOT WORK AT X = 0
result <- lapply(betaDistTests, test_AD2, verbose = FALSE)

nimbleOptions(enableDerivs= EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
