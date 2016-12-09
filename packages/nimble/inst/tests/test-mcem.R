source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of MCEM")

pumpCode <- nimbleCode({
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta)
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1,1.0)
})

pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                       31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))

pump <- nimbleModel(code = pumpCode, name = 'pump',
                       constants = pumpConsts,
                       data = pumpData,
                       inits = pumpInits,
                       check = FALSE)

compileNimble(pump)

#build an MCEM algorithm with Ascent-based convergence criterion
pumpMCEM <- buildMCEM(model = pump,
                      latentNodes = 'theta', burnIn = 300,
                      mcmcControl = list(adaptInterval = 100),
                      boxConstraints = list( list( c('alpha', 'beta'),
                                                  limits = c(0, Inf) ) ),
                      C = 0.01, alpha = .01, beta = .01, gamma = .01, buffer = 1e-6)
## C changed from .001 to .01 to make the test run faster
## Correspondingly changed tolerance below from 0.01 to 0.04
set.seed(0)
out <- pumpMCEM(initM = 1000)
names(out) <- NULL

mle <- c(0.82, 1.26)
try(test_that("Test that MCEM finds the MLE: ",
              expect_equal(out, mle, tolerance = 0.04)))
