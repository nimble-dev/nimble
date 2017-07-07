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
                      wrt = list(c('x[1]', 'y[1]'), c('x[1:2]', 'y[1:2]')), testR = TRUE)



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

test_ADModelCalculate(ADMod3, calcNodeNames = list(c('mu', 'y'), c('y[2]'), c(ADMod3$getDependencies(c('mu', 'sigma')))),
                      wrt = list(c('mu', 'y'), c('sigma[1,1]', 'y[1]'), c('mu[1:2]', 'sigma[1:2, 2]')), testR = TRUE,
                      testCompiled = FALSE)