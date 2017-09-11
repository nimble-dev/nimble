source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
context("Testing of derivatives for calculate() for nimbleModels")

ADCode1 <- nimbleCode({
  x[1] ~ dnorm(0, 1)
  x[2] ~ dnorm(0, 1)
  y[1] ~ dnorm(x[1], 1)
  y[2] ~ dnorm(x[2], 1)
})
ADMod1 <- nimbleModel(code = ADCode1, data = list(y = numeric(2)), dimensions = list(y = c(2)),
                      inits = list(x = c(1,1)))
test_ADModelCalculate(ADMod1, name = "ADMod1", calcNodeNames = list(c('x', 'y'), c('y[2]'), c(ADMod1$getDependencies('x'))),
                      wrt = list(c('x', 'y'), c('x[1]', 'y[1]'), c('x[1:2]', 'y[1:2]'), c('x[1]', 'y', 'x[2]')), testR = TRUE)

test_ADModelCalculate(ADMod1, name = "ADMod1", calcNodeNames = list(c('x', 'y')),
                      wrt = list(c('x[1]', 'y[1]')), testR = TRUE)


ADCode2 <- nimbleCode({
  y[1:2] ~ dmnorm(z[1:2], diagMat[,])
  z[1:2] <- x[1:2] + c(1,1)
  x[1] ~ dnorm(1, 1)
  x[2] ~ dnorm(1, 1)
})
ADMod2 <- nimbleModel(
  code = ADCode2, dimensions = list(x = 2, y = 2, z = 2), constants = list(diagMat = diag(2)),
  inits = list(x = c(2.1, 1.2), y  = c(-.1,-.2)))
test_ADModelCalculate(ADMod2, name = "ADMod2", calcNodeNames = list(c('x', 'y'), c('y[2]'), c(ADMod2$getDependencies('x'))),
                      wrt = list(c('x[1]', 'y[1]'), c('x[1:2]', 'y[1:2]')), testR = TRUE, testCompiled = FALSE)

ADCode3 <- nimbleCode({
  for(i in 1:2){
    y[i, 1:2] ~ dmnorm(mu[1:2], sigma[1:2, 1:2])
  }
  meanVec[1:2] <- c(0,0)
  mu[1:2] ~ dmnorm(meanVec[1:2], diagMat[1:2, 1:2])
  sigma[1:2, 1:2] ~ dwish(diagMat[1:2, 1:2], 3)
})
simData <- matrix(1:4, nrow = 2, ncol = 2)
ADMod3 <- nimbleModel(
  code = ADCode3, constants = list(diagMat = diag(2)),
  data = list(y = simData), inits = list(mu = c(-1.5, 0.8), sigma = diag(2)))
test_ADModelCalculate(ADMod3, calcNodeNames = list(c('mu', 'y'), c('y[1, 2]'), c(ADMod3$getDependencies(c('mu', 'sigma'))),
                                                   c('sigma', 'y')),
                      wrt = list(c('mu', 'y'), c('sigma[1,1]', 'y[1, 2]'), c('mu[1:2]', 'sigma[1:2, 2]')), testR = TRUE,
                      testCompiled = FALSE)


### State Space Model Test

ADCode4 <- nimbleCode({
  x0 ~ dnorm(0,1)
  x[1] ~ dnorm(x0, 1)
  y[1] ~ dnorm(x[1], var = 2)
  for(i in 2:3) {
    x[i] ~ dnorm(x[i-1], 1)
    y[i] ~ dnorm(x[i], var = 2)
  }
})
testdata = list(y = c(0,1,2))
ADMod4 <- nimbleModel(
  code = ADCode4, data = list(y = 0.5+1:3), inits = list(x0 = 1.23, x = -1:3))
ADMod4$simulate(ADMod4$getDependencies('x'))
test_ADModelCalculate(ADMod4, calcNodeNames = list(c('y'), c('y[2]'), c(ADMod4$getDependencies(c('x', 'x0')))),
                      wrt = list(c('x0'), c('x[0]', 'x[1]', 'y[1]'), c('x[1:2]', 'y[1:2]')), testR = TRUE,
                      testCompiled = FALSE)


dir = nimble:::getBUGSexampleDir('equiv')
Rmodel <- readBUGSmodel('equiv', data = NULL, inits = list(tau = c(.2, .2), pi = 1, phi = 1, mu = 1), dir = dir, useInits = TRUE,
                        check = FALSE)
simulate(Rmodel, Rmodel$getDependencies('d'))

## Higher tolerance for more complex chain rule calculations in this model.
test_ADModelCalculate(Rmodel, calcNodeNames = list(Rmodel$getDependencies('tau'), Rmodel$getDependencies('sigma'),  Rmodel$getDependencies('d'),
                                                   Rmodel$getDependencies('d[1]')),
                      wrt = list(c('tau'), c('sigma'), c('d'), c('d[2]')), testR = TRUE,
                      testCompiled = FALSE, tolerance = .5)


## Ragged arrays not supported for derivs yet.
# dir = nimble:::getBUGSexampleDir('bones')
# Rmodel <- readBUGSmodel('bones', data = NULL, inits = NULL, dir = dir, useInits = TRUE,
#                         check = FALSE)
# test_ADModelCalculate(Rmodel, calcNodeNames = list(Rmodel$getDependencies('theta'), Rmodel$getDependencies('grade')),
#                       wrt = list(c('theta'), c('grade'), c('theta', 'grade'), c('grade[1, 2]')), testR = TRUE,
#                       testCompiled = FALSE)
