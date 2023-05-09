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
nimbleOptions(enableDerivs = TRUE)
source("C:/Users/Paul/Documents/Nimble/nimble/packages/nimble/R/AGHQuadrature.R")
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

  ## Test against output from buildLaplace
  expect_equal(cmLaplace$calcMargLogLik(test.val1), -72.5573684952759521,  tol = 1e-16 )
  expect_equal(cmLaplace$calcMargLogLik(test.val2), -1000.9603427653298695,  tol = 1e-16 )

  ## Values from buildLaplace with testval1 and testval2
  grTestLaplace1 <- c(1.5385734890331042, -6.3490351165007235, -3.1745175582503542)
  grTestLaplace2 <- c(584.3200662134241838, 3706.5010041520263258, 741.3001253250956779)
  
  # sprintf("%.16f", x2)

  ## Check Marginalization
  logLik5 <- cmQuad$calcMargLogLik(test.val1)
  logLikTru <- ll.norm(test.val1)
  expect_equal(logLik5, logLikTru, tol = 1e-14)	## Should be very similar for Normal-Normal Case.

  logLik5 <- cmQuad$calcMargLogLik(test.val2)
  logLikTru <- ll.norm(test.val2)
  expect_equal(logLik5, logLikTru, tol = 1e-14)	## Should be very similar for Normal-Normal Case.

  ## Check Gradient of Marginalization
  gr_quad51 <- cmQuad$gr_margLogLik(test.val1)
  gr_laplace1 <- cmLaplace$gr_margLogLik(test.val1)
  expect_equal(gr_quad51, gr_laplace1, tol = 1e-08)	## Should be very similar for Normal-Normal Case.
  expect_equal(gr_laplace1, grTestLaplace1, tol = 1e-16)	## Compare against buildLaplace


  gr_quad52 <- cmQuad$gr_margLogLik(test.val2)
  gr_laplace2 <- cmLaplace$gr_margLogLik(test.val2)
  expect_equal(gr_quad52, gr_laplace2, tol = 1e-05)	## Approx more different for poor values.
  expect_equal(gr_laplace2, grTestLaplace2, tol = 1e-16)	## Compare against buildLaplace. Should be equivalent.

  expect_equal(gr_quad52, gr.ll.norm(test.val2), tol = 1e-08)	## Approx more different for poor values.

  opt <- cmQuad$calculateMLE(pStart = test.val1)	## Needs decent starting values. 
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
    m$calculate()  

    cm <- compileNimble(m)	  
    mQuad <- buildAGHQuad(model = m, nQuad = 21, control = list(ndeps = rep(1e-8, 3)))
    mLaplace <- buildAGHQuad(model = m, nQuad = 1)
    cQL <- compileNimble(mQuad, mLaplace, project = m)
    cmQuad <- cQL$mQuad
    cmLaplace <- cQL$mLaplace
 
    mod.lme4 <- lme4::lmer(m$y ~ 1 + (1|grp), REML = FALSE, 
	control = lmerControl(optimizer= "optimx", optCtrl  = list(method="L-BFGS-B")))
    sprintf("%.16f",   summary(mod.lme4)$sigma)
    sprintf("%.16f",   lme4::fixef(mod.lme4))
    sprintf("%.16f",   attr(unclass(lme4::VarCorr(mod.lme4))[[1]], 'stddev'))  

    # mod.tmb <- glmmTMB::glmmTMB(m$y ~ 1 + (1|grp))
    # sprintf("%.16f",   summary(mod.tmb)$sigma)
    # sprintf("%.16f",   lme4::fixef(mod.tmb))
    # sprintf("%.16f",   attr(unclass(glmmTMB::VarCorr(mod.tmb))[[1]]$grp, 'stddev'))  

    mleLME4 <- c( 3.5679609790094040, 1.4736809813876610, 0.3925194078627622 )
    mleTMB <-  c( 3.5679629394855974, 1.4736809255475793, 0.3925215998142128 )
    mleLaplace <- cmLaplace$calculateMLE()$par
    mleQuad <- cmQuad$calculateMLE()$par

    expect_equal(mleLaplace, mleLME4, tol = 1e-4)
    expect_equal(mleQuad, mleLME4, tol = 1e-4)
    expect_equal(mleQuad, mleLaplace, tol = 1e-4)

    expect_equal(mleLaplace, mleTMB, tol = 1e-4)
    expect_equal(mleQuad, mleTMB, tol = 1e-3)

    cmQuad$gr_margLogLik(mleLME4) ## Pretty good MLE...
    
    print(values(cm, c('b0', 'sigma1', 'sigma2')))
    cmQuad$calcMargLogLik(round(mleLME4, 2))

    mle <- NULL
    gr <- mle
    res <- mle
    for( i in 1:25 )
    {
	  mQuad <- buildAGHQuad(model = m, nQuad = i)
      cmQuad <- compileNimble(mQuad, project = m)
      # gr <- rbind(gr, cmQuad$gr_margLogLik(round(mleLME4, 2)))
      # res <- c(res, cmQuad$calcMargLogLik(round(mleLME4, 2)))
	  mle <- rbind(mle, cmQuad$calculateMLE(method = 'Nelder-Mead')$par)
    }

    par(mfrow = c(3,1))
    plot(1:i, mle[,1], type = 'l')
    abline(h = mleLME4[1], col = 'red')
    plot(1:i, mle[,2], type = 'l')
    abline(h = mleLME4[2], col = 'red')
    plot(1:i, mle[,3], type = 'l')
    abline(h = mleLME4[3], col = 'red')

    par(mfrow = c(3,1))
    plot(1:i, gr[,1], type = 'l')
    plot(1:i, gr[,2], type = 'l')
    plot(1:i, gr[,3], type = 'l')

    par(mfrow = c(3,1))
    plot(1:i, res, type = 'l')
 
    mleQ <- grQ <- resQ <- NULL
    for(i in 1:50){
	  mleQ <- rbind(mleQ, cmQuad$calculateMLE()$par)
      resQ <- c(resQ, cmQuad$calcMargLogLik(mleQ[i,]))
      grQ <- rbind(grQ, cmQuad$gr_margLogLik(mleQ[i,]))
    }
    # values(cm, c('b0', 'sigma1', 'sigma2'))
    par(mfrow = c(3,1))
    plot(1:i, mleQ[,1], type = 'l', xlab = 'iter', ylab = 'intercept')
    abline(h = mleLME4[1], col = 'red')
    plot(1:i, mleQ[,2], type = 'l', xlab = 'iter', ylab = expression(sigma))
    abline(h = mleLME4[2], col = 'red')
    plot(1:i, mleQ[,3], type = 'l', xlab = 'iter', ylab = expression(sigma[re]))
    abline(h = mleLME4[3], col = 'red')
 
    par(mfrow = c(3,1))
    plot(1:i, grQ[,1], type = 'l')
    plot(1:i, grQ[,2], type = 'l')
    plot(1:i, grQ[,3], type = 'l')

  pstart <- values(cm, c("b0", "sigma1", "sigma2"))
  pstart[2:3] <- log(pstart[2:3])
  fit <- optim(pstart, cmQuad$p_transformed_margLogLik, cmQuad$p_transformed_gr_margLogLik, 
	method = "BFGS", control =list(fnscale = -1, trace = 10, ndeps = rep(0.000000001, 3)), hessian = TRUE)
  cmQuad$p_transformed_gr_margLogLik(fit$par)
  fit.pars <- fit$par
  fit.pars[2:3] <- exp(fit.pars[2:3])
  cmQuad$gr_margLogLik(fit.pars)

  nlminb(pstart, cmQuad$p_transformed_margLogLik, gradient = cmQuad$p_transformed_gr_margLogLik,
       lower = -Inf, upper = Inf)

  fit2 <- optim(values(cm, c("b0", "sigma1", "sigma2")), cmQuad$calcMargLogLik, cmQuad$gr_margLogLik, 
	method = "BFGS", control =list(fnscale = -1, maxit = 1000), hessian = TRUE)
  cmQuad$gr_margLogLik(fit2$par)


  mle1 <- cmQuad$calculateMLE(pStart = c(2, 1, 0.8))$par
  mle2 <- cmQuad$calculateMLE(pStart = c(3.5, 1.7, 0.5))$par
  mle3 <- cmQuad$calculateMLE(pStart = mleLME4)$par

  # Values from Laplace directly.
  mLaplace <- buildLaplace(model = m, control = list(maxit = 500))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, project = m)
  mleL1 <- cL$LaplaceMLE(pStart = c(2, 1, 0.8))$par
  mleL2 <- cL$LaplaceMLE(pStart = c(3.5, 1.7, 0.5))$par
  

  mleL <- grL <- NULL
  for(i in 1:25) {
	mleL <- rbind(mleL, cL$LaplaceMLE(method= "Nelder-Mead")$par) 
    grL <- rbind(grL, cL$gr_Laplace(mleL[i,]))
  }
  
  par(mfrow = c(3,1))
  plot(1:nrow(mleL), mleL[,1], type = 'l', xlab = 'iter', ylab = 'intercept')
  abline(h = mleLME4[1], col = 'red')
  plot(1:nrow(mleL), mleL[,2], type = 'l', xlab = 'iter', ylab = expression(sigma))
  abline(h = mleLME4[2], col = 'red')
  plot(1:nrow(mleL), mleL[,3], type = 'l', xlab = 'iter', ylab = expression(sigma[re]))
  abline(h = mleLME4[3], col = 'red')

  par(mfrow = c(3,1))
  plot(1:nrow(grL), grL[,1], type = 'l')
  plot(1:nrow(grL), grL[,2], type = 'l')
  plot(1:nrow(grL), grL[,3], type = 'l')

# lme4::lmer(m$y ~ 1 + (1|grp), REML = FALSE, control = lmerControl(optimizer= "optimx", optCtrl  = list(method="L-BFGS-B")))

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
    m$calculate()  

    cm <- compileNimble(m)	  
    mQuad <- buildAGHQuad(model = m, nQuad = 21)
    mLaplace <- buildAGHQuad(model = m, nQuad = 1)
    cQL <- compileNimble(mQuad, mLaplace, project = m)
    cmQuad <- cQL$mQuad
    cmLaplace <- cQL$mLaplace
 
    mod.lme4 <- lme4::glmer(m$y ~ 1 + (1|grp), family = "poisson", nAGQ = 21,	## Choose equal number of nodes as lme4.
	control = glmerControl(optimizer= "optimx", optCtrl  = list(method="L-BFGS-B")))
    sprintf("%.16f",   lme4::fixef(mod.lme4))
    sprintf("%.16f",   attr(unclass(lme4::VarCorr(mod.lme4))[[1]], 'stddev'))  

    mleLME4 <- c( 3.5136587320416126, 0.4568722479747411)
    mleLaplace <- cmLaplace$calculateMLE()$par
    mleQuad <- cmQuad$calculateMLE()$par

    expect_equal(mleLaplace, mleLME4, tol = 1e-5)
    expect_equal(mleQuad, mleLME4, tol = 1e-4)
    expect_equal(mleQuad, mleLaplace, tol = 1e-4)

    expect_equal(mleLaplace, mleTMB, tol = 1e-4)
    expect_equal(mleQuad, mleTMB, tol = 1e-3)

    cmQuad$gr_margLogLik(mleLME4) ## Pretty good MLE...
    
    print(values(cm, c('b0', 'sigma1', 'sigma2')))
    cmQuad$calcMargLogLik(round(mleLME4, 2))

    mle <- NULL
    gr <- mle
    res <- mle
    for( i in 1:25 )
    {
	  mQuad <- buildAGHQuad(model = m, nQuad = i)
      cmQuad <- compileNimble(mQuad, project = m)
      # gr <- rbind(gr, cmQuad$gr_margLogLik(round(mleLME4, 2)))
      # res <- c(res, cmQuad$calcMargLogLik(round(mleLME4, 2)))
	  mle <- rbind(mle, cmQuad$calculateMLE(method = 'Nelder-Mead')$par)
    }

    par(mfrow = c(3,1))
    plot(1:i, mle[,1], type = 'l')
    abline(h = mleLME4[1], col = 'red')
    plot(1:i, mle[,2], type = 'l')
    abline(h = mleLME4[2], col = 'red')
    plot(1:i, mle[,3], type = 'l')
    abline(h = mleLME4[3], col = 'red')

    par(mfrow = c(3,1))
    plot(1:i, gr[,1], type = 'l')
    plot(1:i, gr[,2], type = 'l')
    plot(1:i, gr[,3], type = 'l')

    par(mfrow = c(3,1))
    plot(1:i, res, type = 'l')
 
    mleQ <- grQ <- resQ <- NULL
    for(i in 1:50){
	  mleQ <- rbind(mleQ, cmQuad$calculateMLE()$par)
      resQ <- c(resQ, cmQuad$calcMargLogLik(mleQ[i,]))
      grQ <- rbind(grQ, cmQuad$gr_margLogLik(mleQ[i,]))
    }
    # values(cm, c('b0', 'sigma1', 'sigma2'))
    par(mfrow = c(3,1))
    plot(1:i, mleQ[,1], type = 'l', xlab = 'iter', ylab = 'intercept')
    abline(h = mleLME4[1], col = 'red')
    plot(1:i, mleQ[,2], type = 'l', xlab = 'iter', ylab = expression(sigma))
    abline(h = mleLME4[2], col = 'red')
    plot(1:i, mleQ[,3], type = 'l', xlab = 'iter', ylab = expression(sigma[re]))
    abline(h = mleLME4[3], col = 'red')
 
    par(mfrow = c(3,1))
    plot(1:i, grQ[,1], type = 'l')
    plot(1:i, grQ[,2], type = 'l')
    plot(1:i, grQ[,3], type = 'l')

  pstart <- values(cm, c("b0", "sigma1", "sigma2"))
  pstart[2:3] <- log(pstart[2:3])
  fit <- optim(pstart, cmQuad$p_transformed_margLogLik, cmQuad$p_transformed_gr_margLogLik, 
	method = "BFGS", control =list(fnscale = -1, trace = 10, ndeps = rep(0.000000001, 3)), hessian = TRUE)
  cmQuad$p_transformed_gr_margLogLik(fit$par)
  fit.pars <- fit$par
  fit.pars[2:3] <- exp(fit.pars[2:3])
  cmQuad$gr_margLogLik(fit.pars)

  nlminb(pstart, cmQuad$p_transformed_margLogLik, gradient = cmQuad$p_transformed_gr_margLogLik,
       lower = -Inf, upper = Inf)

  fit2 <- optim(values(cm, c("b0", "sigma1", "sigma2")), cmQuad$calcMargLogLik, cmQuad$gr_margLogLik, 
	method = "BFGS", control =list(fnscale = -1, maxit = 1000), hessian = TRUE)
  cmQuad$gr_margLogLik(fit2$par)


  mle1 <- cmQuad$calculateMLE(pStart = c(2, 1, 0.8))$par
  mle2 <- cmQuad$calculateMLE(pStart = c(3.5, 1.7, 0.5))$par
  mle3 <- cmQuad$calculateMLE(pStart = mleLME4)$par

  # Values from Laplace directly.
  mLaplace <- buildLaplace(model = m, control = list(maxit = 500))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, project = m)
  mleL1 <- cL$LaplaceMLE(pStart = c(2, 1, 0.8))$par
  mleL2 <- cL$LaplaceMLE(pStart = c(3.5, 1.7, 0.5))$par
  

  mleL <- grL <- NULL
  for(i in 1:25) {
	mleL <- rbind(mleL, cL$LaplaceMLE(method= "Nelder-Mead")$par) 
    grL <- rbind(grL, cL$gr_Laplace(mleL[i,]))
  }
  
  par(mfrow = c(3,1))
  plot(1:nrow(mleL), mleL[,1], type = 'l', xlab = 'iter', ylab = 'intercept')
  abline(h = mleLME4[1], col = 'red')
  plot(1:nrow(mleL), mleL[,2], type = 'l', xlab = 'iter', ylab = expression(sigma))
  abline(h = mleLME4[2], col = 'red')
  plot(1:nrow(mleL), mleL[,3], type = 'l', xlab = 'iter', ylab = expression(sigma[re]))
  abline(h = mleLME4[3], col = 'red')

  par(mfrow = c(3,1))
  plot(1:nrow(grL), grL[,1], type = 'l')
  plot(1:nrow(grL), grL[,2], type = 'l')
  plot(1:nrow(grL), grL[,3], type = 'l')

# lme4::lmer(m$y ~ 1 + (1|grp), REML = FALSE, control = lmerControl(optimizer= "optimx", optCtrl  = list(method="L-BFGS-B")))

})
