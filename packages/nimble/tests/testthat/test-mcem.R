# Set up numerically precise MLE and vcov results
# for both regular and transformed coordinates
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                         31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))


pump_llh <- function(p, trans = FALSE) {
  # for theta ~ dgamma(shape = alpha, rate = beta)
  # we get lambda = t*theta ~ dgamma(shape = alpha, rate = beta/t )
  # For Y ~ dpois(lambda), marginal negative binomial
  # has size = alpha and prob = (rate / (rate+1) = (beta/t) / (beta/t + 1) = beta / (beta+t))
  alpha <- p[1]
  beta <- p[2]
  if(trans) {
    alpha <- exp(alpha)
    beta <- exp(beta)
  }
  llh <- 0
  for(i in 1:10) {
    size <- alpha
    rate <- beta / pumpConsts$t[i]
    prob <- rate / (rate + 1)
    llh <- llh + dnbinom(pumpData$x[i], size = size, prob = prob, log=TRUE)
  }
  llh
}

pump_optim <- optim(c(.75, 1.5), pump_llh,
                  lower = c(0, 0), upper = c(Inf, Inf),
                  method = "L-BFGS-B", control = list(fnscale = -1),
                  hessian=TRUE)
pump_mle <- pump_optim$par
pump_vcov <- inverse(-pump_optim$hessian)

pump_optim_trans <- optim(log(c(.75, 1.5)), pump_llh,
                          method = "BFGS", control = list(fnscale = -1),
                          hessian=TRUE, trans = TRUE)
pump_mle_trans <- pump_optim_trans$par
pump_vcov_trans <- inverse(-pump_optim_trans$hessian)

test_that("MCEM error-trapping", {
 pumpCode <- nimbleCode({
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta);
      lambda[i] <- theta[i]*t[i];
      x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0);
  beta ~ dgamma(0.1,1.0);
 })

 pumpConsts <- list(N = 10,
               t = c(94.3, 15.7, 62.9, 126, 5.24,
                 31.4, 1.05, 1.05, 2.1, 10.5))

 pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

 pumpInits <- list(alpha = 1, beta = 1,
              theta = rep(0.1, pumpConsts$N))
 pumpModel <- nimbleModel(code = pumpCode, constants = pumpConsts,
                          data = pumpData, inits = pumpInits,
                          buildDerivs=TRUE)

 # Catch old syntax usage (for nimble < 1.2.0)
 expect_error(pumpMCEM <- buildMCEM(model = pumpModel, burnIn=200))
 # Correct message when boxConstraints will be ignored
 expect_message(pumpMCEM <- buildMCEM(model = pumpModel,
                                      control=list(boxConstraints = list(a = 1))),
                "boxConstraints will be ignored")

 # Catch when node setup is fully provided by some latent nodes don't exist.
 defaultMargNodes <- setupMargNodes(pumpModel)
 defaultMargNodes$randomEffectsNodes <- c(defaultMargNodes$randomEffectsNodes, "junk")
 expect_error(pumpMCEM <- buildMCEM(model = pumpModel, paramNodes=defaultMargNodes))
})

test_that("MCEM pump with all defaults (roxygen example)", {
 pumpCode <- nimbleCode({
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta);
      lambda[i] <- theta[i]*t[i];
      x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0);
  beta ~ dgamma(0.1,1.0);
 })

 pumpConsts <- list(N = 10,
               t = c(94.3, 15.7, 62.9, 126, 5.24,
                 31.4, 1.05, 1.05, 2.1, 10.5))

 pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

 pumpInits <- list(alpha = 1, beta = 1,
              theta = rep(0.1, pumpConsts$N))
 pumpModel <- nimbleModel(code = pumpCode, constants = pumpConsts,
                          data = pumpData, inits = pumpInits,
                          buildDerivs=TRUE)

 pumpMCEM <- buildMCEM(model = pumpModel)

 CpumpModel <- compileNimble(pumpModel)

 CpumpMCEM <- compileNimble(pumpMCEM, project=pumpModel)

 MLE <- CpumpMCEM$findMLE()
 vcov <- CpumpMCEM$vcov(MLE$par)

 expect_equal(MLE$par, pump_mle, tol = 0.03)
 expect_equal(vcov, pump_vcov, tol = 0.03)

 samples <- as.matrix(CpumpMCEM$mvSamples)
 expect_equal(ncol(samples), 10)
 expect_true(nrow(samples)>=1000)
})

test_that("MCEM with pump transform=FALSE",{
  pumpCode <- nimbleCode({
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })

  pump <- nimbleModel(code = pumpCode, name = 'pump',
                      constants = pumpConsts,
                      data = pumpData,
                      inits = pumpInits,
                      #                    check = FALSE,
                      buildDerivs = TRUE)

  cpump <- compileNimble(pump)

  ## Basic version: Use box constraints
  ## build an MCEM algorithm with Ascent-based convergence criterion
  pumpMCEM <- buildMCEM(model = pump,
                        latentNodes = 'theta',
                        control = list(burnIn = 300,
                                       mcmcControl = list(adaptInterval = 100),
                                       boxConstraints = list( list( c('alpha', 'beta'),
                                                                   limits = c(0, Inf) ) ),
                                       tol = 0.01, alpha = .01, beta = .01, gamma = .01,
                                       buffer = 1e-6,
                                       transform=FALSE))
  cpumpMCEM <- compileNimble(pumpMCEM, project = pump)

  # Run the basic MCEM that uses no transformation
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, Mfactor = 0.5)
  expect_equal(out$par, pump_mle, tol = 0.03)
  vc <- cpumpMCEM$vcov(out$par, trans = FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)
  # Although using the transformation or not is baked into the optimization
  # the vcov procedure allows to return with or without it.
  # In this case, there should be no difference since there is no
  # transormation in use.
  vct <- cpumpMCEM$vcov(out$par, trans = FALSE, returnTrans = TRUE)
  expect_identical(vct, vc)

  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 500, maxM = 500, burnin = 100, maxIter = 1)
  expect_false(abs(out$par[2] - pump_mle[2]) < .02)
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000,
                       burnin = 300, maxIter = 20, continue = TRUE)
  expect_equal(out$par, pump_mle, tol = 0.03)

  # Next check on the continue=TRUE feature and see if
  # we get to convergence going one by one from the outside
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, burnin = 300,
                       Mfactor = 0.5, maxIter = 1)
  for(i in 1:4)
    out <- cpumpMCEM$findMLE(continue = TRUE)
  expect_equal(out$par, pump_mle, tol = 0.03)
  vc <- cpumpMCEM$vcov(out$par, trans = FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)

  # Next check that it will run with ascent=FALSE
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, burnin = 300,
                       Mfactor = 0.5, ascent=FALSE)
  expect_equal(out$par, pump_mle, tol = 0.03)

  # Next check that it will run with ascent=TRUE but
  # meaningful limit set by maxM
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 2000, burnin = 300,
                       Mfactor = 0.5, ascent=TRUE)
  expect_equal(out$par, pump_mle, tol = 0.03)

  # Next check that it will run with minIter=2
  # and also with C = 3 (converge at least three times)
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 2000, burnin = 300,
                       Mfactor = 0.5, ascent=TRUE, minIter = 2,
                       C = 3)
  expect_equal(out$par, pump_mle, tol = 0.03)
})

test_that("MCEM pump with transform=TRUE", {
  # transformed version: Use parameter transformation
  ## build an MCEM algorithm with Ascent-based convergence criterion
  pumpCode <- nimbleCode({
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })

  pump <- nimbleModel(code = pumpCode, name = 'pump',
                      constants = pumpConsts,
                      data = pumpData,
                      inits = pumpInits,
                      #                    check = FALSE,
                      buildDerivs = TRUE)

  cpump <- compileNimble(pump)

  pumpMCEM <- buildMCEM(model = pump,
                        latentNodes = 'theta',
                        control = list(burnin = 300,
                                       mcmcControl = list(adaptInterval = 100),
                                       tol = 0.01, alpha = .01, beta = .01, gamma = .01,
                                       buffer = 1e-6,
                                       transform = TRUE))
  cpumpMCEM <- compileNimble(pumpMCEM, project = pump)

  ## tol changed from .001 to .01 to make the test run faster
  ## Correspondingly expect_equal tolerance
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000)
  expect_equal(out$par, pump_mle, tol = 0.01)
  vc <- cpumpMCEM$vcov(out$par, trans = FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)
  vct <- cpumpMCEM$vcov(out$par, trans = FALSE, returnTrans=TRUE)
  expect_equal(vct, pump_vcov_trans, tol=0.03)

  # Now check when we request parameters returned on the transformed
  # scale
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, returnTrans=TRUE)
  expect_equal(out$par, pump_mle_trans, tol=0.01)
  vc <- cpumpMCEM$vcov(out$par, trans = TRUE)
  expect_equal(vc, pump_vcov_trans, tol=0.02)
  vct <- cpumpMCEM$vcov(out$par, trans = TRUE, returnTrans=FALSE)
  expect_equal(vct, pump_vcov, tol=0.03)
})

test_that("MCEM pump without AD, with transform=FALSE", {
  ####
  # version without AD
  pumpCode <- nimbleCode({
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })

  pump <- nimbleModel(code = pumpCode, name = 'pump',
                      constants = pumpConsts,
                      data = pumpData,
                      inits = pumpInits,
                      #                    check = FALSE,
                      buildDerivs = TRUE)

  cpump <- compileNimble(pump)

  pumpMCEM <- buildMCEM(model = pump,
                        latentNodes = 'theta',
                        control = list(useDerivs=FALSE,
                                       burnIn = 300,
                                       mcmcControl = list(adaptInterval = 100),
                                       boxConstraints = list( list( c('alpha', 'beta'),
                                                                   limits = c(0, Inf) ) ),
                                       tol = 0.01, alpha = .01, beta = .01, gamma = .01,
                                       buffer = 1e-6,
                                       transform=FALSE))
  cpumpMCEM <- compileNimble(pumpMCEM, project = pump)
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, Mfactor = 0.5)
  expect_equal(out$par, pump_mle, tol = 0.03)

  vc <- cpumpMCEM$vcov(out$par, trans=FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)
  vc <- cpumpMCEM$vcov(out$par, trans=FALSE, resetSamples = TRUE, M = 10000)
  expect_equal(vc, pump_vcov, tol = 0.03)
  vct <- cpumpMCEM$vcov(out$par, trans=FALSE, returnTrans = TRUE)
  expect_identical(vc, vct)
})


test_that("MCEM pump without AD, with transform=TRUE", {
  ####
  # version without AD
  pumpCode <- nimbleCode({
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })

  pump <- nimbleModel(code = pumpCode, name = 'pump',
                      constants = pumpConsts,
                      data = pumpData,
                      inits = pumpInits,
                      #                    check = FALSE,
                      buildDerivs = TRUE)

  cpump <- compileNimble(pump)

  pumpMCEM <- buildMCEM(model = pump,
                        latentNodes = 'theta',
                        control = list(useDerivs=FALSE,
                                       burnIn = 300,
                                       mcmcControl = list(adaptInterval = 100),
                                       tol = 0.01, alpha = .01, beta = .01, gamma = .01,
                                       buffer = 1e-6,
                                       transform=TRUE))
  cpumpMCEM <- compileNimble(pumpMCEM, project = pump)

  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, Mfactor = 0.5)
  expect_equal(out$par, pump_mle, tol = 0.01)
  vc <- cpumpMCEM$vcov(out$par, trans = FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)
  vct <- cpumpMCEM$vcov(out$par, trans = FALSE, returnTrans=TRUE)
  expect_equal(vct, pump_vcov_trans, tol=0.03)

  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, Mfactor = 0.5, returnTrans = TRUE)
  vc <- cpumpMCEM$vcov(out$par, trans=TRUE)
  expect_equal(vc, pump_vcov_trans, tol = 0.03)
  vc <- cpumpMCEM$vcov(out$par, trans=TRUE, returnTrans=FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)
})

## Test of case where some parts of the model
## are not related to missing data or latent states

test_that("MCEM with pump transform=FALSE and some calcNodesOther",{
  pumpConsts <- list(N = 10,
                     t = c(94.3, 15.7, 62.9, 126, 5.24,
                           31.4, 1.05, 1.05, 2.1, 10.5))

  pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22),
                   theta = c(rep(NA, 8), .8, 1.2))

  pumpInits <- list(alpha = 1, beta = 1,
                    theta = c(rep(0.1, pumpConsts$N-2), NA, NA))


  pump_llh <- function(p, trans = FALSE) {
                                        # for theta ~ dgamma(shape = alpha, rate = beta)
                                        # we get lambda = t*theta ~ dgamma(shape = alpha, rate = beta/t )
                                        # For Y ~ dpois(lambda), marginal negative binomial
                                        # has size = alpha and prob = (rate / (rate+1) = (beta/t) / (beta/t + 1) = beta / (beta+t))
    alpha <- p[1]
    beta <- p[2]
    if(trans) {
      alpha <- exp(alpha)
      beta <- exp(beta)
    }
    llh <- 0
    for(i in 1:8) {
      size <- alpha
      rate <- beta / pumpConsts$t[i]
      prob <- rate / (rate + 1)
      llh <- llh + dnbinom(pumpData$x[i], size = size, prob = prob, log=TRUE)
    }
    for(i in 9:10) {
      llh <- llh + dgamma(pumpData$theta[i], shape = alpha, rate = beta, log=TRUE)
      llh <- llh + dpois(pumpData$x[i], pumpData$theta[i], log=TRUE)
    }
    llh
  }

  pump_optim <- optim(c(.75, 1.5), pump_llh,
                      lower = c(0, 0), upper = c(Inf, Inf),
                      method = "L-BFGS-B", control = list(fnscale = -1),
                      hessian=TRUE)
  pump_mle <- pump_optim$par
  pump_vcov <- inverse(-pump_optim$hessian)

  pumpCode <- nimbleCode({
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })

  pump <- nimbleModel(code = pumpCode, name = 'pump',
                      constants = pumpConsts,
                      data = pumpData,
                      inits = pumpInits,
                      #                    check = FALSE,
                      buildDerivs = TRUE)

  cpump <- compileNimble(pump)

  ## Basic version: Use box constraints
  ## build an MCEM algorithm with Ascent-based convergence criterion
  pumpMCEM <- buildMCEM(model = pump,
                        latentNodes = 'theta[1:8]',
                        control = list(burnIn = 300,
                                       mcmcControl = list(adaptInterval = 100),
                                       boxConstraints = list( list( c('alpha', 'beta'),
                                                                   limits = c(0, Inf) ) ),
                                       tol = 0.01, alpha = .01, beta = .01, gamma = .01,
                                       buffer = 1e-6,
                                       transform=FALSE))
  expect_equal(length(pumpMCEM$E$calcNodesOther), 3)
  cpumpMCEM <- compileNimble(pumpMCEM, project = pump)

  # Run the basic MCEM that uses no transformation
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, Mfactor = 0.5, C = 2)
  expect_equal(out$par, pump_mle, tol = 0.01)
  vc <- cpumpMCEM$vcov(out$par, trans = FALSE, atMLE=FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)
})


test_that("MCEM with pump transform=FALSE and some calcNodesOther, noAD",{
  pumpConsts <- list(N = 10,
                     t = c(94.3, 15.7, 62.9, 126, 5.24,
                           31.4, 1.05, 1.05, 2.1, 10.5))

  pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22),
                   theta = c(rep(NA, 8), .8, 1.2))

  pumpInits <- list(alpha = 1, beta = 1,
                    theta = c(rep(0.1, pumpConsts$N-2), NA, NA))


  pump_llh <- function(p, trans = FALSE) {
                                        # for theta ~ dgamma(shape = alpha, rate = beta)
                                        # we get lambda = t*theta ~ dgamma(shape = alpha, rate = beta/t )
                                        # For Y ~ dpois(lambda), marginal negative binomial
                                        # has size = alpha and prob = (rate / (rate+1) = (beta/t) / (beta/t + 1) = beta / (beta+t))
    alpha <- p[1]
    beta <- p[2]
    if(trans) {
      alpha <- exp(alpha)
      beta <- exp(beta)
    }
    llh <- 0
    for(i in 1:8) {
      size <- alpha
      rate <- beta / pumpConsts$t[i]
      prob <- rate / (rate + 1)
      llh <- llh + dnbinom(pumpData$x[i], size = size, prob = prob, log=TRUE)
    }
    for(i in 9:10) {
      llh <- llh + dgamma(pumpData$theta[i], shape = alpha, rate = beta, log=TRUE)
      llh <- llh + dpois(pumpData$x[i], pumpData$theta[i], log=TRUE)
    }
    llh
  }

  pump_optim <- optim(c(.75, 1.5), pump_llh,
                      lower = c(0, 0), upper = c(Inf, Inf),
                      method = "L-BFGS-B", control = list(fnscale = -1),
                      hessian=TRUE)
  pump_mle <- pump_optim$par
  pump_vcov <- inverse(-pump_optim$hessian)

  pumpCode <- nimbleCode({
    for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
  })

  pump <- nimbleModel(code = pumpCode, name = 'pump',
                      constants = pumpConsts,
                      data = pumpData,
                      inits = pumpInits,
                      #                    check = FALSE,
                      buildDerivs = FALSE)

  cpump <- compileNimble(pump)

  ## Basic version: Use box constraints
  ## build an MCEM algorithm with Ascent-based convergence criterion
  pumpMCEM <- buildMCEM(model = pump,
                        latentNodes = 'theta[1:8]',
                        control = list(burnIn = 300,
                                       mcmcControl = list(adaptInterval = 100),
                                       boxConstraints = list( list( c('alpha', 'beta'),
                                                                   limits = c(0, Inf) ) ),
                                       tol = 0.01, alpha = .01, beta = .01, gamma = .01,
                                       buffer = 1e-6,
                                       transform=FALSE,
                                       useDerivs=FALSE))
  expect_equal(length(pumpMCEM$E$calcNodesOther), 3)
  cpumpMCEM <- compileNimble(pumpMCEM, project = pump)

  # Run the basic MCEM that uses no transformation
  set.seed(0)
  cpump$alpha <- pump$alpha
  cpump$beta <- pump$beta
  out <- cpumpMCEM$findMLE(initM = 1000, maxM = 20000, Mfactor = 0.5, C = 2)
  expect_equal(out$par, pump_mle, tol = 0.01)
  vc <- cpumpMCEM$vcov(out$par, trans = FALSE, atMLE=FALSE)
  expect_equal(vc, pump_vcov, tol = 0.03)
})


#######
# Adopting some (far from all) tests from test-ADlaplace

test_that("MCEM simplest 1D works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd = 2)
      a ~ dnorm(mu, sd = 3)
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )
  mMCEM <- buildMCEM(model = m)
  cm <- compileNimble(m)
  cMCEM <- compileNimble(mMCEM, project = m)
  set.seed(0)
  opt <- cMCEM$findMLE()
  expect_equal(opt$par, 4, tol = .01)
  vc <- cMCEM$vcov(opt$par)
  # See test-ADlaplace for derivation
  expect_equal(vc[1,1], 13, tol = .02)
})

test_that("MCEM simplest 1D with a constrained parameter works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd = 2)
      a ~ dnorm(mu, sd = 3)
      mu ~ dexp(1.0)
      # I am not sure why Laplace works if init mu is 0 (log(0)=-Inf)
      # but this doesn't. Both use BFGS.
    }), data = list(y = 4), inits = list(a = -1, mu = 0.5),
    buildDerivs = TRUE
  )

  MCEM <- buildMCEM(model = m)
  cm <- compileNimble(m)
  cMCEM <- compileNimble(MCEM, project = m)

  set.seed(0)
  opt <- cMCEM$findMLE(C = 2)
  vc <- cMCEM$vcov(opt$par)
  # V[a] = 9
  # V[y] = 9 + 4 = 13
  # Cov[a, y] = V[a] = 9 (not needed)
  # y ~ N(mu, 13)
  # muhat = y = 4
  # ahat = (9*y+4*mu)/(9+4) = y = 4
  # Jacobian of ahat wrt transformed param log(mu) is 4/13*mu = 4*mu/13 = 16/13
  # Hessian of joint loglik wrt a: -(1/4 + 1/9)
  # Hessian of marginal loglik wrt transformed param log(mu) is (y*mu - 2*mu*mu)/13 = -4^2/13
  # Variance of transformed param is 13/16
  expect_equal(opt$par, 4, tol = 0.01)

  # Covariance matrix (of mu and a) on transformed scale
  vcov_transform <- matrix(c(0, 0, 0, 1/(1/4+1/9)), nrow = 2) + matrix(c(1, 16/13), ncol = 1) %*% (13/16) %*% t(matrix(c(1, 16/13), ncol = 1))
  # Covariance matrix on original scale
  vcov <- diag(c(4, 1)) %*% vcov_transform %*% diag(c(4, 1))
  expect_equal(vcov[1,1], vc[1,1], tol = 0.1)
  vct <- cMCEM$vcov(opt$par, returnTrans=TRUE)
  expect_equal(vcov_transform[1,1], vct[1,1], tol = 0.03)

  # Get results on transformed scale
  set.seed(0)
  cm$mu <- 0.5
  opt <- cMCEM$findMLE(C = 2, returnTrans = TRUE)
  expect_equal(opt$par, log(4), tol = 0.01)
  vc <- cMCEM$vcov(opt$par, trans = TRUE, returnTrans = FALSE)
  expect_equal(vcov[1,1], vc[1,1], tol = 0.1)
  vct <- cMCEM$vcov(opt$par, trans = TRUE, returnTrans = TRUE)
  expect_equal(vcov_transform[1,1], vct[1,1], tol = 0.03)
})

# Skipping " simplest 1D (constrained) with multiple data"
# The main difference there should be simply model structure
# that we don't need to test here.

test_that("MCEM simplest 1D with deterministic intermediates and multiple data works", {
  set.seed(1)
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, sd = 5)
      a ~ dexp(rate = exp(0.5 * mu))
      for (i in 1:5){
        y[i] ~ dnorm(0.2 * a, sd = 2)
      }
    }), data = list(y = rnorm(5, 1, 2)), inits = list(mu = 2, a = 1),
    buildDerivs = TRUE
  )
  MCEM <- buildMCEM(model = m)
  cm <- compileNimble(m)
  cMCEM <- compileNimble(MCEM, project = m)

  cm$mu <- 2
  opt <- cMCEM$findMLE(C = 2)
  vc <- cMCEM$vcov(opt$par)
  # Results are checked using those from TMB
  # See test-ADlaplace for checking results from TMB
  # However, for this little toy, we get a slightly different estimate
  # and so we check it with grid-based integration of the llh
  #
  #
  ##   true_llh <- function(mu) {
  ##   lambda <- exp(0.5 * mu)
  ##   a_max <- qexp(0.9999, lambda)
  ##   a_grid <- seq(0, a_max, length = 1000)
  ##   delta_a <- diff(a_grid[1:2])
  ##   a_grid <- a_grid[-1]-delta_a/2
  ##   lh <- 0
  ##   y <- cm$y
  ##   for(i in 1:999) {
  ##     llh_given_a <- dnorm(y, mean = 0.2*a_grid[i], sd = 2, log=TRUE) |>
  ##       sum()
  ##     lh <- lh + exp(llh_given_a) * dexp(a_grid[i], lambda)*delta_a
  ##   }
  ##   llh <- log(lh)
  ##   llh
  ## }
  ## true_ans <- optimize(true_llh, lower=-4, upper=4, maximum=TRUE)
  ## true_mle <- true_ans$maximum # -2.856
  ## true_hess <- nimDerivs(true_llh(true_mle))
  ## true_var <- 1/(-true_hess$hessian[1,1,1]) # 8.150
  ## true_sd <- sqrt(true_var)

  expect_equal(opt$par, -2.856, .04)
  expect_equal(vc[1,1], 8.150, 0.1)
  # Check covariance matrix for params only
})

# Case "1D with deterministic intermediates works" from Laplace seems redundant

test_that("MCEM 1D with a constrained parameter and deterministic intermediates works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(0.2 * a, sd = 2)
      a ~ dnorm(0.5 * mu, sd = 3)
      mu ~ dexp(1.0)
    }), data = list(y = 4), inits = list(a = -1, mu = 0.5),
    buildDerivs = TRUE
  )

  MCEM <- buildMCEM(model = m)
  cm <- compileNimble(m)
  cMCEM <- compileNimble(MCEM, project = m)

  set.seed(0)
  opt <- cMCEM$findMLE(tol = 0.0001)
  vc <- cMCEM$vcov(opt$par)
  # Correct answer first by direct integration and
  # next by theory
  ## true_llh <- function(mu) {
  ##   mean_a <- 0.5 * mu
  ##   a_min <- qnorm(0.00001, mean = mean_a, sd = 3)
  ##   a_max <- qnorm(0.99999, mean = mean_a, sd = 3)
  ##   a_grid <- seq(a_min, a_max, length = 1000)
  ##   delta_a <- diff(a_grid[1:2])
  ##   a_grid <- a_grid[-1]-delta_a/2
  ##   lh <- 0
  ##   y <- cm$y
  ##   for(i in 1:999) {
  ##     llh_given_a <- dnorm(y, mean = 0.2*a_grid[i], sd = 2, log=TRUE) |>
  ##       sum()
  ##     lh <- lh + exp(llh_given_a) * dnorm(a_grid[i], mean_a, sd=3)*delta_a
  ##   }
  ##   llh <- log(lh)
  ##   llh
  ## }
  ## true_ans <- optimize(true_llh, lower=20, upper=60, maximum=TRUE)
  ## true_mle <- true_ans$maximum # 40
  ## true_hess <- nimDerivs(true_llh(true_mle))
  ## true_var <- 1/(-true_hess$hessian[1,1,1]) # 436.007
  ## true_sd <- sqrt(true_var)
  # And by theory
  # V[a] = 9
  # V[y] = 0.2^2 * 9 + 4 = 4.36
  # y ~ N(0.2*0.5*mu, 4.36)
  # muhat = y/(0.2*0.5) = 40
  # ahat = (9*0.2*y + 4*0.5*mu)/(4+9*0.2^2) = 20
  # Jacobian of ahat wrt transformed param log(mu) is 4*0.5*mu/(4+9*0.2^2) = 18.34862
  # Hessian of joint loglik wrt a: -(0.2^2/4 + 1/9)
  # Hessian of marginal loglik wrt param mu is -(0.2*0.5)^2/4.36
  # Hessian of marginal loglik wrt transformed param log(mu) is (0.2*0.5*y*mu - 2*0.1^2*mu*mu)/4.36 = -3.669725
  vcov_transform <- matrix(c(0, 0, 0, 1/(0.2^2/4+1/9)), nrow = 2) + matrix(c(1, 18.34862), ncol = 1) %*% (1/3.669725) %*% t(matrix(c(1, 18.34862), ncol = 1))
  # Covariance matrix on original scale
  vcov <- diag(c(40,1)) %*% vcov_transform %*% diag(c(40,1))

  expect_equal(opt$par, 40, tol = 1)
  expect_equal(vc[1,1], vcov[1,1], tol = 50) # out of 436

  cm$mu <- 0.5
  opt <- cMCEM$findMLE(tol = 0.0001, returnTrans = TRUE)
  expect_equal(opt$par, log(40), tol = 0.02)
  vc <- cMCEM$vcov(opt$par, trans=TRUE, returnTrans = TRUE, M = 50000)
  expect_equal(vc[1,1], vcov_transform[1,1], tol = 0.1) #
})

# Skipping "with deterministic intermediates and multiple data works"

test_that("MCEM simplest 2x1D works, with multiple data for each", {
  set.seed(1)
  y <- matrix(rnorm(6, 4, 5), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) {
        mu_y[i] <- 0.8*a[i]
        for(j in 1:3)
          y[i, j] ~ dnorm(mu_y[i], sd = 2)
        a[i] ~ dnorm(mu_a, sd = 3)
      }
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = y), inits = list(a = c(-2, -1), mu = 0),
    buildDerivs = TRUE
  )

  MCEM <- buildMCEM(model = m)
  cm <- compileNimble(m)
  cMCEM <- compileNimble(MCEM, project = m)

  set.seed(0)
  opt <- cMCEM$findMLE()
  expect_equal(opt$par, mean(y)/(0.8*0.5), tol = 0.04) # optim's reltol is about 1e-8 but that is for the value, not param.
  vc <- cMCEM$vcov(opt$par)
  ## # V[a] = 9
  ## # V[y[i]] = 0.8^2 * 9 + 4 = 9.76
  ## # Cov[a, y[i]] = 0.8 * 9 = 7.2
  ## # Cov[y[i], y[j]] = 0.8^2 * 9 = 5.76, within a group
  ## Cov_A_Y <- matrix(nrow = 8, ncol = 8)
  ## Cov_A_Y[1, 1:8] <- c(  9,    0,  7.2,  7.2,  7.2,    0,    0,    0)
  ## Cov_A_Y[2, 1:8] <- c(  0,    9,    0,    0,    0,  7.2,  7.2,  7.2)
  ## Cov_A_Y[3, 1:8] <- c(7.2,    0, 9.76, 5.76, 5.76,    0,    0,    0)
  ## Cov_A_Y[4, 1:8] <- c(7.2,    0, 5.76, 9.76, 5.76,    0,    0,    0)
  ## Cov_A_Y[5, 1:8] <- c(7.2,    0, 5.76, 5.76, 9.76,    0,    0,    0)
  ## Cov_A_Y[6, 1:8] <- c(  0,  7.2,    0,    0,    0, 9.76, 5.76, 5.76)
  ## Cov_A_Y[7, 1:8] <- c(  0,  7.2,    0,    0,    0, 5.76, 9.76, 5.76)
  ## Cov_A_Y[8, 1:8] <- c(  0,  7.2,    0,    0,    0, 5.76, 5.76, 9.76)
  ## Cov_Y <- Cov_A_Y[3:8, 3:8]
  ## chol_cov <- chol(Cov_Y)
  ## res <- dmnorm_chol(as.numeric(t(y)), mean(y), cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  ## expect_equal(opt$value, res)
  ## # muhat = mean(y)/(0.8*0.5)
  ## # ahat[1] = (9*0.8*sum(y[1,]) + 4*0.5*mu)/(4+9*0.8^2*3)
  ## # ahat[2] = (9*0.8*sum(y[2,]) + 4*0.5*mu)/(4+9*0.8^2*3)
  ## # Jacobian of ahat[i] wrt mu is 4*0.5/(4+9*0.8^2*3) = 0.09398496
  ## # Hessian of joint loglik wrt a[i]a[i]: -(3*0.8^2/4 + 1/9); wrt a[i]a[j]: 0
  ## # Hessian of marginal loglik wrt mu: -0.04511278 (numerical, have not got AD work)
  muhat <- mean(y)/(0.8*0.5)
  # Covariance matrix
  vcov <- diag(c(0, rep(1/(3*0.8^2/4 + 1/9), 2))) + matrix(c(1, rep(0.09398496, 2)), ncol = 1) %*% (1/0.04511278) %*% t(matrix(c(1, rep(0.09398496, 2)), ncol = 1))
  expect_equal(vcov[1,1], vc[1,1], tol = 0.02)
  ## Check covariance matrix for params only
})

#skipping " with 2x1D random effects needing joint integration works, without intermediate nodes"
# and "Laplace with 2x1D random effects needing joint integration works, with intermediate nodes"

# This case would need more attention
# The MCEM convergence path is painfully slow.
# It might be better with a reparameterization.
# But for now it is thought that this test is not really needed
# for any specific functionality.
## test_that("MCEM with 2x2D random effects for 1D data that are separable works, with intermediate nodes", {
##   # This example gives a really bad MCEM convergence path
##   set.seed(1)
##   # y[i, j] is jth datum from ith group
##   y <- array(rnorm(8, 6, 5), dim = c(2, 2, 2))
##   cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
##   m <- nimbleModel(
##     nimbleCode({
##       for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
##       mu_a[1] <- 0.8 * mu[1]
##       mu_a[2] <- 0.2 * mu[2]
##       for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
##       for(i in 1:2) {
##         for(j in 1:2) {
##           y[1, j, i] ~ dnorm( 0.5 * a[i, 1], sd = 1.8) # this ordering makes it easier below
##           y[2, j, i] ~ dnorm( 0.1 * a[i, 2], sd = 1.2)
##         }
##       }
##     }),
##     data = list(y = y),
##     inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2), mu = c(0, .5)),
##     constants = list(cov_a = cov_a),
##     buildDerivs = TRUE
##   )

##   MCEM <- buildMCEM(model = m)
##   cm <- compileNimble(m)
##   cMCEM <- compileNimble(MCEM, project = m)

##   Laplace <- buildLaplace(m)
##   cLaplace <- compileNimble(Laplace, project = m)
##   MLE <- cLaplace$findMLE()
##   cLaplace$calcLogLik(c(12.9840, 405.9512))
##   cLaplace$calcLogLik(c(12.2434, 391.044))
##   # this takes a long time. the convergence path is poor.
##   set.seed(0)
##   cm$mu <- c(0, 0.5)
##   opt <- cMCEM$findMLE(alpha = 0.1, maxIter = 300, maxM = 50000)
## })

test_that("MCMC for simple LME case works", {
  set.seed(1)
  g <- rep(1:10, each = 5)
  n <- length(g)
  x <- runif(n)
  mc <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ dnorm((random_int[g[i]]) +
                   (random_slope[g[i]])*x[i], sd = sigma_res)
    }
    for(i in 1:ng) {
      random_int[i] ~ dnorm(fixed_int , sd = sigma_int)
      random_slope[i] ~ dnorm(fixed_slope, sd = sigma_slope)
    }
    sigma_int ~ dhalfflat() #dunif(0, 10)
    sigma_slope ~ dhalfflat() #dunif(0, 10)
    sigma_res ~ dhalfflat() #dunif(0, 10)
    fixed_int ~ dnorm(0, sd = 100)
    fixed_slope ~ dnorm(0, sd = 100)
  })
  m <- nimbleModel( mc,
    constants = list(g = g, ng = max(g), n = n, x = x),
    buildDerivs = TRUE
  )
  params <- c("fixed_int", "fixed_slope", "sigma_int", "sigma_slope", "sigma_res")
  values(m, params) <- c(10, 0.5, 3, .25, 0.2)
  m$simulate(m$getDependencies(params, downstream = TRUE, self = FALSE))
  m$setData('y')
  y <- m$y
  cm <- compileNimble(m)

  # We can do this by lme4, but instead we'll use nimble's Laplace
  #library(lme4)
  #manual_fit <- lmer(y ~ x + (1 + x || g), REML = FALSE)
  #data <- data.frame(y = y, x = x, g = g)

  # The reason to use our Laplace is to compare MLEs in terms of
  # the log likelihood itself.

  MCEM <- buildMCEM(model = m, latentNodes = c("random_int", "random_slope"),
                    control = list(
                      mcmcControl = list(adaptInterval = 20)))
  cMCEM <- compileNimble(MCEM, project = m)

  set.seed(1)
  values(cm, params) <- c(0, 0, 1, 1, 1)
  opt <- cMCEM$findMLE(maxIter = 2000, initM = 1000, burnin = 300,
                   maxM = 50000, tol = 0.001)

  m2 <- nimbleModel( mc,
                    constants = list(g = g, ng = max(g), n = n, x = x),
                    data = list(y = m$y),
    buildDerivs = TRUE
  )
  for(v in m2$getVarNames()) m2[[v]] <- m[[v]]
  m2$calculate()

  cm2 <- compileNimble(m2)
  Laplace <- buildLaplace(model=m2, randomEffectsNodes = c("random_int", "random_slope"))
  cLaplace <- compileNimble(Laplace, project = m2)
  MLE <- cLaplace$findMLE()

  expect_equal(MLE$value, cLaplace$calcLogLik(opt$par), tolerance = 0.04)
})
