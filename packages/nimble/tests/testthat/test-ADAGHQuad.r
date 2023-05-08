# Tests of Laplace approximation
source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

library(nimble)

# check internal consistency of optim method variants
check_laplace_alternative_methods <- function(cL, # compiled laplace algorithm
                                              cm, # compiled model
                                              m,  # original model (or list with values)
                                              opt, # possibly already-run LaplaceMLE result,
                                              methods = 1:3, # methods to check
                                              summ_orig, # summarized Laplace MLE result (original)
                                              summ_trans # summarized Laplace MLE result (transformed)
) {
  vars <- cm$getVarNames()
  reset <- function() {
    for(v in vars) cm[[v]] <- m[[v]]
  }
  if(missing(opt)) {
    reset()
    opt <- cL$LaplaceMLE()
  }
  if(missing(summ_orig)){
    summ_orig <- cL$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  }
  if(missing(summ_trans)){
    summ_trans <- cL$summary(opt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  }
  ref_method <- cL$get_method()
  for(method in methods) {
    if(method != ref_method) {
      reset()
      cL$set_method(method)
      opt_alt <- cL$LaplaceMLE()
      expect_equal(opt$par, opt_alt$par, tolerance = 0.01)
      expect_equal(opt$value, opt_alt$value, tolerance = 1e-7)
      summ_orig_alt <- cL$summary(opt_alt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
      summ_trans_alt <- cL$summary(opt_alt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
      expect_equal(summ_orig$params$estimate, summ_orig_alt$params$estimate, tol = 1e-5)
      expect_equal(summ_orig$random$estimate, summ_orig_alt$random$estimate, tol = 1e-5)
      expect_equal(summ_orig$params$stdError, summ_orig_alt$params$stdError, tol = 1e-5)
      expect_equal(summ_orig$random$stdError, summ_orig_alt$random$stdError, tol = 1e-5)
      expect_equal(summ_orig$vcov, summ_orig_alt$vcov, tol = 1e-5)
      expect_equal(summ_trans$params$estimate, summ_trans_alt$params$estimate, tol = 1e-5)
      expect_equal(summ_trans$random$estimate, summ_trans_alt$random$estimate, tol = 1e-5)
      expect_equal(summ_trans$params$stdError, summ_trans_alt$params$stdError, tol = 1e-5)
      expect_equal(summ_trans$random$stdError, summ_trans_alt$random$stdError, tol = 1e-5)
      expect_equal(summ_trans$vcov, summ_trans_alt$vcov, tol = 1e-5)
    }
  }
  invisible(NULL)
}

test_that("AGH Quadrature Normal-Normal 1D works", {
set.seed(123)
n <- 50
m <- nimbleModel(nimbleCode({
  # priors 
  b0 ~ dnorm(0, 1000)
  sigma1 ~ dunif(0, 1000)
  sigma2 ~ dunif(0, 1000)
  for(i in 1:n){
	b[i] ~ dnorm(mean = 0, sd = sigma2)
	mu[i] <- b0 + b[i]
	y[i] ~ dnorm(mean = mu[i], sd = sigma1)
  }}), data = list(y=rnorm(n, 5, sqrt(1 + 0.5^2))), constants = list(n=n),
  inits = list(b0 = 3.5, sigma1 = 1, sigma2 = 1), buildDerivs = TRUE)
	  
  mQuad <- buildAGHQuad(model = m)
  mLaplace <- buildAGHQuad(model = m, nQuad = 1)
  cm <- compileNimble(m)
  cQL <- compileNimble(mQuad, mLaplace, project = m)
  cmQuad <- cQL$mQuad
  cmLaplace <- cQL$mLaplace
 
  ll.norm <- function(pars)
  {
	b0 <- pars[1]
	sigma1 <- pars[2]
	sigma2 <- pars[3]
	ll <- 0
	for( i in seq_along(m$y) ) {
		ll <- ll + dnorm( m$y[i], mean = b0, sd = sqrt(sigma1^2 + sigma2^2), log = TRUE ) 
	}
	ll
  }
  gr.ll.norm <- function(pars)
  {
	b0 <- pars[1]
	sigma1 <- pars[2]
	sigma2 <- pars[3]
	var0 <- sigma1^2+sigma2^2
	dll <- 0
	for( i in seq_along(m$y) ) {
		di <- m$y[i]-b0
		dll <- dll + -0.5*(c(0, 2*sigma1, 2*sigma2)/var0) - c(-2*di/(2*var0), -0.5*(di)^2*var0^(-2)*2*sigma1, -0.5*(di)^2*var0^(-2)*2*sigma2 )
	}
	dll
  }
  test.val1 <- c(5, 1, 0.5)
  test.val2 <- c(2, 0.5, 0.1)
  
  ## Check Marginalization
  logLik5 <- cmQuad$calcMargLogLik(test.val1)
  logLikTru <- ll.norm(test.val1)
  expect_equal(logLik5, logLikTru, tol = 1e-12)	## Should be very similar for Normal-Normal Case.

  logLik5 <- cmQuad$calcMargLogLik(test.val2)
  logLikTru <- ll.norm(test.val2)
  expect_equal(logLik5, logLikTru, tol = 1e-12)	## Should be very similar for Normal-Normal Case.

  ## Check Gradient of Marginalization
  gr_quad5 <- cmQuad$gr_margLogLik(test.val1)
  gr_laplace <- cmLaplace$gr_margLogLik(test.val1)
  expect_equal(gr_quad5, gr_laplace, tol = 1e-08)	## Should be very similar for Normal-Normal Case.

  gr_quad5 <- cmQuad$gr_margLogLik(test.val2)
  gr_laplace <- cmLaplace$gr_margLogLik(test.val2)
  expect_equal(gr_quad5, gr_laplace, tol = 1e-05)	## Approx more different for poor values.

  expect_equal(gr_quad5, gr.ll.norm(test.val2), tol = 1e-08)	## Approx more different for poor values.

  opt <- cmQuad$calculateMLE(pStart = test.val1)	## Needs decent starting values. 
  mle.tru <- optim(test.val1, ll.norm, gr.ll.norm, control = list(fnscale = -1))
  expect_equal(opt$value, mle.tru$value, tol = 1e-8) # Same log likelihood. Diff parameter values.

  ## Check covariance?
})

test_that("AGH Quadrature 1D Poisson-Gamma for checking nQuad", {
set.seed(123)
n <- 1
m <- nimbleModel(nimbleCode({
  # priors 
  a ~ dunif(0, 1000)
  b ~ dunif(0, 1000)
  for(i in 1:n){
	lambda[i] ~ dgamma(a, b)
	y[i] ~ dpois(lambda[i])
  }}), data = list(y=rpois(n, rgamma(n, 10, 2))), constants = list(n=n),
  inits = list(a = 10, b = 2), buildDerivs = TRUE)

## Marginal model is equivalent to a negative binomial with this parameterization.
m.nb <- nimbleModel(nimbleCode({
  # priors 
  a ~ dunif(0, 1000)
  b ~ dunif(0, 1000)
  for(i in 1:n){
	y[i] ~ dnbinom(size=a, prob=b/(1+b))
  }}), data = list(y=m$y), constants = list(n=n),
  inits = list(a = 10, b = 2), buildDerivs = TRUE)

  cm <- compileNimble(m)
  mQuad <- buildAGHQuad(model = m, nQuad = 20)
  cmQuad <- compileNimble(mQuad, project = m)
  test.val1 <- c(10, 2)
  test.val2 <- c(1, 5)
 
  ## Check Marginalization
  logLik5 <- cmQuad$calcMargLogLik(test.val1)
  m.nb$a <- test.val1[1]; m.nb$b <- test.val1[2]
  logLikTru <- m.nb$calculate("y")
  expect_equal(logLik5, logLikTru, tol = 1e-10)	## Should be very similar.

  ## Check Marginalization
  logLik5 <- cmQuad$calcMargLogLik(test.val2)
  m.nb$a <- test.val2[1]
  m.nb$b <- test.val2[2]
  logLikTru <- m.nb$calculate("y")
  expect_equal(logLik5, logLikTru, tol = 1e-8)	## Should be very similar.

  ## True Gradient
  tru.ll.gr <- function(pars)
  {
	a <- pars[1]
	b <- pars[2]
	da <- digamma(a+m$y) + log(b) - digamma(a) - log(b + 1)
	db <- a/b - (a + m$y)/(b+1)
	return(c(da,db))
  }
  
  tru.gr <- tru.ll.gr(c(50,2))
  m.nb$a <- 50; m.nb$b <- 2
  tru.logLik <- m.nb$calculate('y')
 
  ## Check node accuracy
  ## Tolerance should decrease in loop,
  res <- NULL; gr <- NULL
  for( i in 1:25 )
  {
	mQuad <- buildAGHQuad(model = m, nQuad = i)
    cmQuad <- compileNimble(mQuad, project = m)
	expect_equal(cmQuad$calcMargLogLik(c(50,2)), tru.logLik, tol = 0.01^(sqrt(i)))
	expect_equal(cmQuad$gr_margLogLik(c(50,2)), tru.gr, tol = 0.01^(sqrt(i)))
  }
})

test_that("AGH Quadrature 1D Binomial-Beta check 3 methods", {
  set.seed(123)
  n <- 50
  N <- 5
  m <- nimbleModel(nimbleCode({
    # priors 
    a ~ dgamma(1,1)
    b ~ dgamma(1,1)
    for(i in 1:n){
	  p[i] ~ dbeta(a, b)
	  y[i] ~ dbinom(p[i], N)
    }
  }), data = list(y = rbinom(n, N, rbeta(n, 10, 2))), 
    constants = list(N = N, n=n), inits = list(a = 10, b = 2), 
    buildDerivs = TRUE)
  
  cm <- compileNimble(m)
  mQuad <- buildAGHQuad(model = m)
  cmQuad <- compileNimble(mQuad, project = m)

  param.val <- c(7, 1)

  cmQuad$set_method(1)
  ll.11 <- cmQuad$calcMargLogLik(param.val)
  ll.12 <- cmQuad$calcMargLogLik(param.val+1)
  gr.11 <- cmQuad$gr_margLogLik(param.val)
  gr.12 <- cmQuad$gr_margLogLik(param.val+1)
  cmQuad$set_method(2)
  ll.21 <- cmQuad$calcMargLogLik(param.val)
  ll.22 <- cmQuad$calcMargLogLik(param.val+1)
  gr.21 <- cmQuad$gr_margLogLik(param.val)
  gr.22 <- cmQuad$gr_margLogLik(param.val+1)
  cmQuad$set_method(3)
  ll.31 <- cmQuad$calcMargLogLik(param.val)
  ll.32 <- cmQuad$calcMargLogLik(param.val+1)
  gr.31 <- cmQuad$gr_margLogLik(param.val)
  gr.32 <- cmQuad$gr_margLogLik(param.val+1)

  ## All the methods should return equivalent results, or at least nearly with some small
  ## numerical differences from the calls to the AD.
  expect_equal(ll.11, ll.21)
  expect_equal(ll.11, ll.31)
  expect_equal(ll.12, ll.22)
  expect_equal(ll.12, ll.32)
  expect_equal(gr.11, gr.21)
  expect_equal(gr.11, gr.31)
  expect_equal(gr.12, gr.22)  
  expect_equal(gr.12, gr.32)

  ## Check gradient and marginalization accuracy.
  ll.betabin <- function(pars){
    a <- pars[1]
    b <- pars[2]
    ll <- 0
    for( i in seq_along(m$y))
    {
      ll <- ll - lbeta(a,b) + lchoose(N, m$y[i]) + lbeta(a + m$y[i], b + N-m$y[i])
    }
    ll
  }

  dibeta <- function(a,b)
  {
    da <-  digamma(a) - digamma(a+b) 
	db <-  digamma(b) - digamma(a+b)
	c(da, db)
  }

  gr.betabin <- function(pars){
	a <- pars[1]
	b <- pars[2]
	dll <- 0
	for( i in seq_along(m$y))
	{
      dll <- dll - dibeta(a,b) + dibeta(a + m$y[i], b + N-m$y[i])
	}
	return(dll)
  }
  
  ## Check against 'exact'
  expect_equal(ll.21, ll.betabin(param.val), tol = 1e-03)	## Accuracy isn't 'great'
  expect_equal(gr.21, gr.betabin(param.val), tol = 1e-02)	

  ## Crank up the nodes to check accuracy.
  mQuad <- buildAGHQuad(model = m, nQuad = 15)
  cmQuad <- compileNimble(mQuad, project = m)
  ll.15 <- cmQuad$calcMargLogLik(param.val)
  gr.15 <- cmQuad$gr_margLogLik(param.val)

  expect_equal(ll.15, ll.betabin(param.val), tol = 1e-05)	## Accuracy is better.
  expect_equal(gr.15, gr.betabin(param.val), tol = 1e-05)	
  
  ## Lots of nodes.
  mQuad <- buildAGHQuad(model = m, nQuad = 50)
  cmQuad <- compileNimble(mQuad, project = m)
  ll.50 <- cmQuad$calcMargLogLik(param.val)
  gr.50 <- cmQuad$gr_margLogLik(param.val)

  expect_equal(ll.50, ll.betabin(param.val), tol = 1e-10)	## Accuracy is very good.
  expect_equal(gr.50, gr.betabin(param.val), tol = 1e-10)

})

test_that("AGH Quadrature 1D Check MLE.", {
  set.seed(123)
  n <- 50
  N <- 50
  m <- nimbleModel(nimbleCode({
    # priors 
    a ~ dgamma(1,1)
    b ~ dgamma(1,1)
    for(i in 1:n){
	  p[i] ~ dbeta(a, b)
	  y[i] ~ dbinom(p[i], N)
    }
  }), data = list(y = rbinom(n, N, rbeta(n, 10, 2))), 
    constants = list(N = N, n=n), inits = list(a = 10, b = 2), 
    buildDerivs = TRUE)
  
  cm <- compileNimble(m)
  mQuad <- buildAGHQuad(model = m)
  cmQuad <- compileNimble(mQuad, project = m)

  ## Check gradient and marginalization accuracy.
  ll.betabin <- function(pars){
    a <- exp(pars[1])
    b <- exp(pars[2])
    ll <- 0
    for( i in seq_along(m$y))
    {
      ll <- ll - lbeta(a,b) + lchoose(N, m$y[i]) + lbeta(a + m$y[i], b + N-m$y[i])
    }
    ll
  }

  dibeta <- function(a,b)
  {
    da <-  digamma(a) - digamma(a+b) 
	db <-  digamma(b) - digamma(a+b)
	c(da, db)
  }

  gr.betabin <- function(pars){
	a <- exp(pars[1])
	b <- exp(pars[2])
	dll <- 0
	for( i in seq_along(m$y))
	{
      dll <- dll - dibeta(a,b) + dibeta(a + m$y[i], b + N-m$y[i])
	}
	return(dll)
  }
  
  mle.tru <- optim(log(c(10,2)), ll.betabin, gr.betabin, control = list(fnscale = -1))
  mle.par <- exp(mle.tru$par)
  
  ## Check with 5 quad points.
  mle.quad <- cmQuad$calculateMLE()
  expect_equal(mle.quad$par, mle.par, tol = 1e-02)
  expect_equal(mle.quad$value, mle.tru$value, tol = 1e-03)
  
  ## Check with 51 quad points.
  mQuad <- buildAGHQuad(model = m, nQuad = 51)
  cmQuad <- compileNimble(mQuad, project = m)
  mle.quad51 <- cmQuad$calculateMLE()
  expect_equal(mle.quad51$par, mle.par, tol = 1e-04)
  expect_equal(mle.quad51$value, mle.tru$value, tol = 1e-08)

 
})