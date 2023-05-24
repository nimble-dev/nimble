# Tests of AGH Quadrature approximation
source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

library(nimble)

## Temp for Paul's unbuilt testing:
library(testthat)

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
	  
  mQuad <- buildAGHQuad(model = m, nQuad = 5)
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

  ## Test against output from buildLaplace
  expect_equal(cmLaplace$calcLogLik(test.val1), -72.5573684952759521,  tol = 1e-16 )
  expect_equal(cmLaplace$calcLogLik(test.val2), -1000.9603427653298695,  tol = 1e-16 )

  ## Values from buildLaplace with testval1 and testval2
  grTestLaplace1 <- c(1.5385734890331042, -6.3490351165007235, -3.1745175582503542)
  grTestLaplace2 <- c(584.3200662134241838, 3706.5010041520263258, 741.3001253250956779)
  
  ## Check Marginalization
  logLik5 <- cmQuad$calcLogLik(test.val1)
  logLikTru <- ll.norm(test.val1)
  expect_equal(logLik5, logLikTru, tol = 1e-14)	## Should be very similar for Normal-Normal Case.

  logLik5 <- cmQuad$calcLogLik(test.val2)
  logLikTru <- ll.norm(test.val2)
  expect_equal(logLik5, logLikTru, tol = 1e-14)	## Should be very similar for Normal-Normal Case.

  ## Check Gradient of Marginalization
  gr_quad51 <- cmQuad$gr_logLik(test.val1)
  gr_laplace1 <- cmLaplace$gr_logLik(test.val1)
  expect_equal(gr_quad51, gr_laplace1, tol = 1e-08)	## Should be very similar for Normal-Normal Case.
  expect_equal(gr_laplace1, grTestLaplace1, tol = 1e-16)	## Compare against buildLaplace

  gr_quad52 <- cmQuad$gr_logLik(test.val2)
  gr_laplace2 <- cmLaplace$gr_logLik(test.val2)
  expect_equal(gr_quad52, gr_laplace2, tol = 1e-05)	## Approx more different for poor values.
  expect_equal(gr_laplace2, grTestLaplace2, tol = 1e-16)	## Compare against buildLaplace. Should be equivalent.

  expect_equal(gr_quad52, gr.ll.norm(test.val2), tol = 1e-08)	## Approx more different for poor values.

  opt <- cmQuad$findMLE(pStart = test.val1)	## Needs decent starting values. 
  mle.tru <- optim(test.val1, ll.norm, gr.ll.norm, control = list(fnscale = -1))
  expect_equal(opt$value, mle.tru$value, tol = 1e-8) # Same log likelihood. Diff parameter values.

  ## Check covariance?

  # Values from Laplace directly.
  # mLaplace <- buildLaplace(model = m)
  # cm <- compileNimble(m)
  # cL <- compileNimble(mLaplace, project = m)
  # x1 <- cL$Laplace(test.val1)
  # x2 <- cL$Laplace(test.val2)
  # g1 <- cL$gr_Laplace(test.val1)
  # g2 <- cL$gr_Laplace(test.val2)
  # sprintf("%.16f", x1)
  # sprintf("%.16f", x2)
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
  logLik20 <- cmQuad$calcLogLik(test.val1)
  m.nb$a <- test.val1[1]; m.nb$b <- test.val1[2]
  m.nb$calculate()
  logLikTru <- m.nb$calculate("y")
  expect_equal(logLik20, logLikTru, tol = 1e-10)	## Should be very similar.

  ## Check Marginalization
  logLik20.2 <- cmQuad$calcLogLik(test.val2)
  m.nb$a <- test.val2[1]; m.nb$b <- test.val2[2]
  m.nb$calculate()
  logLikTru.2 <- m.nb$calculate("y")
  expect_equal(logLik20.2, logLikTru.2, tol = 1e-8)	## Should be very similar.

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
  m.nb$a <- 50; m.nb$b <- 2; m.nb$calculate()
  tru.logLik <- m.nb$calculate('y')
 
  ## Check node accuracy
  ## Tolerance should decrease in loop,
  res <- NULL; gr <- NULL
  for( i in 1:25 )
  {
	mQuad <- buildAGHQuad(model = m, nQuad = i)
    cmQuad <- compileNimble(mQuad, project = m)
	expect_equal(cmQuad$calcLogLik(c(50,2)), tru.logLik, tol = 0.01^(sqrt(i)))
	expect_equal(cmQuad$gr_logLik(c(50,2)), tru.gr, tol = 0.01^(sqrt(i)))
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
  mQuad <- buildAGHQuad(model = m, nQuad = 5)
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
  mle.quad <- cmQuad$findMLE(pStart = c(10,2))
  expect_equal(mle.quad$par, mle.par, tol = 1e-02)
  expect_equal(mle.quad$value, mle.tru$value, tol = 1e-03)
  
  ## Check with 51 quad points.
  mQuad <- buildAGHQuad(model = m, nQuad = 51)
  cmQuad <- compileNimble(mQuad, project = m)
  mle.quad51 <- cmQuad$findMLE(pStart = c(10,2))
  expect_equal(mle.quad51$par, mle.par, tol = 1e-04)
  expect_equal(mle.quad51$value, mle.tru$value, tol = 1e-08)
})

test_that("AGH Quadrature Comparison to LME4 1 RE", {
  set.seed(123)
  n <- 50
  J <- 10
  nobs <- n*J
  grp <- rep(1:n, each = J)
  m <- nimbleModel(nimbleCode({
    # priors 
    b0 ~ dnorm(0, 1000)
    # sigma1 ~ dunif(0, 1000)
    # sigma2 ~ dunif(0, 1000)
    sigma1 ~ dgamma(1,1)
    sigma2 ~ dgamma(1,1)
    for(i in 1:n){
	  b[i] ~ dnorm(mean = 0, sd = sigma2)
	  mu[i] <- b0 + b[i]
    }  
    for(i in 1:nobs){
      y[i] ~ dnorm(mean = mu[grp[i]], sd = sigma1)
    }}), constants = list(n=n, nobs=nobs, grp = grp),
    inits = list(b = rnorm(n, 0, 0.5), b0 = 3.5, sigma1 = 1.5, sigma2 = 0.5), buildDerivs = TRUE)
  m$simulate('y')
  m$setData('y')
  m$calculate()  

  cm <- compileNimble(m)	  
  mQuad <- buildAGHQuad(model = m, nQuad = 21, control = list(outOptimControllist = list(reltol = 1e-16)))
  mLaplace <- buildAGHQuad(model = m, nQuad = 1, control = list(outOptimControllist = list(reltol = 1e-16)))
  cQL <- compileNimble(mQuad, mLaplace, project = m)
  cmQuad <- cQL$mQuad
  cmLaplace <- cQL$mLaplace
 
    # mod.lme4 <- lme4::lmer(m$y ~ 1 + (1|grp), REML = FALSE, 
	# control = lmerControl(optimizer= "optimx", optCtrl  = list(method="L-BFGS-B")))
    # sprintf("%.16f",   summary(mod.lme4)$sigma)
    # sprintf("%.16f",   lme4::fixef(mod.lme4))
    # sprintf("%.16f",   attr(unclass(lme4::VarCorr(mod.lme4))[[1]], 'stddev'))  

    # mod.tmb <- glmmTMB::glmmTMB(m$y ~ 1 + (1|grp))
    # sprintf("%.16f",   summary(mod.tmb)$sigma)
    # sprintf("%.16f",   lme4::fixef(mod.tmb))
    # sprintf("%.16f",   attr(unclass(glmmTMB::VarCorr(mod.tmb))[[1]]$grp, 'stddev'))  

  mleLME4 <- c( 3.5679609790094040, 1.4736809813876610, 0.3925194078627622 )
  mleTMB <-  c( 3.5679629394855974, 1.4736809255475793, 0.3925215998142128 )
  mleLaplace <- cmLaplace$findMLE()$par
  mleQuad <- cmQuad$findMLE()$par

  expect_equal(mleLaplace, mleLME4, tol = 1e-7)
  expect_equal(mleQuad, mleLME4, tol = 1e-7)
  expect_equal(mleQuad, mleLaplace, tol = 1e-7)

  expect_equal(mleLaplace, mleTMB, tol = 1e-6)
  expect_equal(mleQuad, mleTMB, tol = 1e-6)

  gr_mle <- cmQuad$gr_logLik(mleLME4) ## MLE gradient check.
  expect_equal(gr_mle, c(0,0,0), tol = 1e-5)

  ## Compare MLE after running twice.
  mleLaplace2 <- cmLaplace$findMLE()$par
  mleQuad2 <- cmQuad$findMLE()$par
  expect_equal(mleLaplace, mleLaplace2, tol = 1e-8) 
  expect_equal(mleQuad, mleQuad2, tol = 1e-8)

})


## This might be better to compare for MLE as lme4 does some different
## optimization steps for LMMs.
test_that("AGH Quadrature Comparison to LME4 1 RE for Poisson-Normal", {
  set.seed(123)
  n <- 50
  J <- 10
  nobs <- n*J
  grp <- rep(1:n, each = J)
  m <- nimbleModel(nimbleCode({
    # priors 
    b0 ~ dnorm(0, 1000)
    sigma ~ dgamma(1,1)
    for(i in 1:n){
	  b[i] ~ dnorm(mean = 0, sd = sigma)
	  mu[i] <- exp(b0 + b[i])
    }  
    for(i in 1:nobs){
      y[i] ~ dpois(mu[grp[i]])
    }}), constants = list(n=n, nobs=nobs, grp = grp),
    inits = list(b = rnorm(n, 0, 0.5), b0 = 3.5, sigma = 0.5), buildDerivs = TRUE)
  m$simulate('y')
  m$setData('y')
  m$calculate()  

  cm <- compileNimble(m)	  
  mQuad <- buildAGHQuad(model = m, nQuad = 21, control = list(outOptimControllist = list(reltol = 1e-16)))
  mLaplace <- buildAGHQuad(model = m, nQuad = 1, control = list(outOptimControllist = list(reltol = 1e-16)))
  cQL <- compileNimble(mQuad, mLaplace, project = m)
  cmQuad <- cQL$mQuad
  cmLaplace <- cQL$mLaplace

    # library(lme4)
    # mod.lme4 <- lme4::glmer(m$y ~ 1 + (1|grp), family = "poisson", nAGQ = 1,	## nAGQ = 21 for nQuad = 21
	  # control = glmerControl(optimizer= "optimx", optCtrl  = list(method="L-BFGS-B")))
    # sprintf("%.16f",   lme4::fixef(mod.lme4))
    # sprintf("%.16f",   attr(unclass(lme4::VarCorr(mod.lme4))[[1]], 'stddev'))  
    # sprintf("%.16f", logLik(mod.lme4, fixed.only = TRUE) )

  lme4_laplace <- -1695.4383630192869532
	# lme4_nquad21 <- -360.6258864811468356  ## Only proportional in lme4, not the actual marginal. Can't use to compare.

  mleLME4_nquad21 <- c( 3.5136587320416126, 0.4568722479747411)
  mleLME4_laplace <- c( 3.5136586190857675, 0.4568710881066258)
  mleLaplace <- cmLaplace$findMLE()$par
  mleQuad <- cmQuad$findMLE()$par

	## Compare the marginal log likelihood for the laplace method.
  logLikLaplace <- cmLaplace$calcLogLik( mleLME4_laplace )
  expect_equal(logLikLaplace, lme4_laplace, tol = 1e-11) ## Very very similar maximization. Even if slightly different estimates.

  ## Compare mle for laplace to lme4 laplace
  ## and nQuad = 21 for both methods
  expect_equal(mleLaplace, mleLME4_laplace, tol = 1e-5)
  expect_equal(mleQuad, mleLME4_nquad21, tol = 1e-5)
  expect_equal(mleQuad, mleLaplace, tol = 1e-5)
  
  ## Compare MLE after running twice.
  mleLaplace2 <- cmLaplace$findMLE()$par
  mleQuad2 <- cmQuad$findMLE()$par
  expect_equal(mleLaplace, mleLaplace2, tol = 1e-6) 
  expect_equal(mleQuad, mleQuad2, tol = 1e-8)
})
