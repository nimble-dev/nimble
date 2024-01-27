source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

nimbleOptions(useADcholAtomic = TRUE)
nimbleOptions(useADsolveAtomic  = TRUE)                              
nimbleOptions(useADmatMultAtomic = TRUE)                             
nimbleOptions(useADmatInverseAtomic  = TRUE)                       

relTol <- eval(formals(test_ADModelCalculate)$relTol)
relTol[3] <- 1e-6
relTol[4] <- 1e-4

verbose <- TRUE

##nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for calculate() for nimbleModels")

## Start of Nick's tests ##

if(FALSE) {

test_that('Derivs of calculate function work for model ADMod1', {
  ADCode1 <- nimbleCode({
    x[1] ~ dnorm(0, 1)
    x[2] ~ dnorm(0, 1)
    y[1] ~ dnorm(x[1], 1)
    y[2] ~ dnorm(x[2], 1)
  })
  ADMod1 <- nimbleModel(code = ADCode1, data = list(y = numeric(2)), dimensions = list(y = c(2)),
                        inits = list(x = c(1,1)))
  test_ADModelCalculate_nick(ADMod1, name = 'ADMod1', calcNodeNames = list(c('x', 'y'), c('y[2]'), c(ADMod1$getDependencies('x'))),
                        wrt = list(c('x', 'y'), c('x[1]', 'y[1]'), c('x[1:2]', 'y[1:2]'), c('x[1]', 'y', 'x[2]')),
                        order = c(0, 1, 2))
})

test_that('Derivs of calculate function work for model ADMod2', {
  ADCode2 <- nimbleCode({
    y[1:2] ~ dmnorm(z[1:2], diagMat[,])
    z[1:2] <- x[1:2] + c(1,1)
    x[1] ~ dnorm(1, 1)
    x[2] ~ dnorm(1, 1)
  })
  ADMod2 <- nimbleModel(
    code = ADCode2, dimensions = list(x = 2, y = 2, z = 2), constants = list(diagMat = diag(2)),
    inits = list(x = c(2.1, 1.2), y  = c(-.1,-.2)))
  test_ADModelCalculate_nick(ADMod2, name = 'ADMod2', calcNodeNames = list(c('x', 'y'), c('y[2]'), c(ADMod2$getDependencies('x'))),
                        wrt = list(c('x[1]', 'y[1]'), c('x[1:2]', 'y[1:2]')), order = c(0, 1, 2))
})

#C++ derivs of below model won't work until we've implemented derivs of wishart
# test_that('R derivs of calculate function work for model ADMod3', {
#   ADCode3 <- nimbleCode({
#     for(i in 1:2){
#       y[i, 1:2] ~ dmnorm(mu[1:2], sigma[1:2, 1:2])
#     }
#     meanVec[1:2] <- c(0,0)
#     mu[1:2] ~ dmnorm(meanVec[1:2], diagMat[1:2, 1:2])
#     sigma[1:2, 1:2] ~ dwish(diagMat[1:2, 1:2], 3)
#   })
#   simData <- matrix(1:4, nrow = 2, ncol = 2)
#   ADMod3 <- nimbleModel(
#     code = ADCode3, constants = list(diagMat = diag(2)),
#     data = list(y = simData), inits = list(mu = c(-1.5, 0.8), sigma = diag(2)))
#   test_ADModelCalculate_nick(ADMod3, name = 'ADMod3', calcNodeNames = list(c('mu', 'y'), c('y[1, 2]'), c(ADMod3$getDependencies(c('mu', 'sigma'))),
#                                                      c('sigma', 'y')),
#                         wrt = list(c('mu', 'y'), c('sigma[1,1]', 'y[1, 2]'), c('mu[1:2]', 'sigma[1:2, 2]')),
#                         testCompiled = TRUE, order = c(0, 1))
# })


### State Space Model Test

test_that('Derivs of calculate function work for model ADMod4', {
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
    code = ADCode4, data = list(y = 0.5+1:3), inits = list(x0 = 1.23, x = 1:3))
  ADMod4$simulate(ADMod4$getDependencies('x'))
  test_ADModelCalculate_nick(ADMod4, name = 'ADMod4', calcNodeNames = list(c('y'), c('y[2]'), c(ADMod4$getDependencies(c('x', 'x0')))),
                        wrt = list(c('x0'), c('x0', 'x[1]', 'y[1]'), c('x[1:2]', 'y[1:2]')), tolerance = .1,
                        order = c(0, 1, 2))
})

test_that('Derivs of calculate function work for model ADmod5 (tricky indexing)', {
  ADCode5 <- nimbleCode({
    y[1:2] ~ dmnorm(z[1:2], diagMat[,])
    z[1] <- x[2]
    z[2:3] <- x[1:2] + c(1,1)
    x[1] ~ dnorm(.1, 10)
    x[2] ~ dnorm(1, 3)
  })
  ADMod5 <- nimbleModel(
    code = ADCode5, dimensions = list(x = 2, y = 2, z = 3), constants = list(diagMat = diag(2)),
    inits = list(x = c(1, 1.2), y  = c(-.1,-.2)))
  test_ADModelCalculate_nick(ADMod5, name = 'ADMod5', calcNodeNames = list(c('y'), c('y[2]'), c(ADMod5$getDependencies(c('x')))),
                        wrt = list(c('x'), c('x[1]', 'z[1]', 'y[1]'), c('x[1:2]', 'y[1:2]')), tolerance = .1,
                        order = c(0, 1, 2))
})


test_that('Derivs of calculate function work for model equiv', {
    ## This test gives a crash on i = 3, j = 1.  It is due to pow(x, y) with y = 0.
    ## 1. Try changing model code to use pow_int instead of pow.
    ## 2. Need to error-trap to avoid crash if pow is used.
  dir = nimble:::getBUGSexampleDir('equiv')
  Rmodel <- readBUGSmodel('equiv', data = NULL, inits = list(tau = c(.2, .2), pi = 1, phi = 1, mu = 1), dir = dir, useInits = TRUE,
                          check = FALSE)
  initModel <- initializeModel(Rmodel)
  initModel$run()
  ## Higher tolerance for more complex chain rule calculations in this model.
  test_ADModelCalculate_nick(Rmodel, name = 'equiv',
                        calcNodeNames = list(Rmodel$getDependencies('tau'),
                                             Rmodel$getDependencies('sigma'),
                                             Rmodel$getDependencies('d'),
                                             Rmodel$getDependencies('d[1]')),
                        wrt = list(c('tau'),
                                   c('sigma'),
                                   c('d'),
                                   c('d[2]')),
                        tolerance = .5,
                        order = c(0, 1, 2))
})

test_that("Derivs of calculate function work for Daniel's SSM", {
ssmCode <- nimbleCode({
  a ~ dunif(-0.9999, 0.9999)
  b ~ dnorm(0, sd = 1000)
  sigPN ~ dunif(1e-04, 1)
  sigOE ~ dunif(1e-04, 1)
  x[1] ~ dnorm(b/(1 - a), sd = sqrt(sigPN^2 + sigOE*sigOE)) #sqrt(sigPN^2 + sigOE^2))
  y[1] ~ dnorm(x[1], sd = sigOE)
  for (i in 2:t) {
    x[i] ~ dnorm(x[i - 1] * a + b, sd = sigPN)
    y[i] ~ dnorm(x[i], sd = sigOE)
  }
})
ssmConsts <- list(t = 100)
ssmData <- list(y = c(20.24405,20.57693,20.49357,20.34159,20.45759,20.43326,20.20554,20.12860,20.14756,20.20781,20.23022,20.26766,20.22984,20.37703,20.13641,20.05309,19.95709,20.19303,20.30562,20.54443,20.91010,20.70580,20.42344,20.19795,20.28816,20.31894,20.76939,20.77023,20.83486,20.29335,20.40990,20.19601,20.04083,19.76056,19.80810,19.83129,19.69174,19.90069,19.87623,19.63371,19.62360,19.72630,19.64450,19.86779,20.17104,20.34797,20.32968,20.48027,20.46694,20.47006,20.51676,20.40695,20.18715,19.97552,19.88331,19.67831,19.74702,19.47502,19.24408,19.37179,19.38277,19.15034,19.08723,19.37051,19.14274,19.46433,19.62459,19.77971,19.54194,19.39081,19.61621,19.51307,19.34745,19.17019,19.26829,19.58943,19.77143,19.83582,19.71198,19.67746,19.75053,20.40197,20.49363,20.37079,20.19005,20.55862,20.48523,20.33071,19.97069,19.79758,19.83811,19.79728,19.86277,19.86836,19.92481,19.88095,20.24899,20.55165,20.22707,20.11235))
ssmInitial <- list(a = 0.95, b=1, sigOE=0.05,sigPN = 0.2,  x= c(20.26036,20.51331,20.57057,20.35633,20.33736,20.47321,20.22002,20.14917,20.19216,20.26969,20.21135,20.22745,20.20466,20.41158,20.13408,20.08023,19.98956,20.13543,20.32709,20.55840,20.88206,20.74740,20.47671,20.14012,20.29953,20.33778,20.80916,20.75773,20.84349,20.35654,20.41045,20.20180,20.02872,19.74226,19.80483,19.81842,19.69770,19.84564,19.88211,19.70559,19.56090,19.73728,19.66545,19.88158,20.13870,20.39163,20.37372,20.47429,20.39414,20.42024,20.55560,20.40462,20.15831,19.89425,19.79939,19.72692,19.74565,19.42233,19.22730,19.36489,19.37289,19.19050,19.00823,19.35738,19.14293,19.48812,19.67329,19.82750,19.58979,19.43634,19.61278,19.56739,19.38584,19.19260,19.32732,19.65500,19.65295,19.84843,19.68285,19.69620,19.77497,20.31795,20.45797,20.32650,20.24045,20.60507,20.51597,20.30076,19.98100,19.86709,19.85965,19.74822,19.86730,19.90523,19.86970,19.87286,20.28417,20.46212,20.22618,20.13689))
Rmodel <- nimbleModel(code = ssmCode, name = 'SSMcorrelated', constants = ssmConsts, data = ssmData, inits = ssmInitial)
test_ADModelCalculate_nick(Rmodel, name = 'SSM',
                      calcNodeNames = list(Rmodel$getDependencies('a'),
                                           Rmodel$getDependencies('b'),
                                           Rmodel$getDependencies('sigPN'),
                                           Rmodel$getDependencies('x')),
                      wrt = list(c('a', 'b'),
                                 c('sigPN'),
                                 c('x[1]')),
                      tolerance = 0.0001,
                      order = c(0, 1, 2))
})

test_that("Derivs of calculate function work for rats model", {
  Rmodel <- readBUGSmodel('rats', dir = getBUGSexampleDir('rats'))
  test_ADModelCalculate_nick(Rmodel, name = 'rats', calcNodeNames = list(Rmodel$getNodeNames(),
                                                                    Rmodel$getDependencies('mu[1, 2]'),
                                                                    Rmodel$getDependencies('alpha'),
                                                                    Rmodel$getDependencies('beta')),
                        wrt = list(c('alpha[1]',
                                     'beta[1]'),
                                   c('mu[1,1]')),
                        tolerance = .0001,
                        order = c(0, 1, 2))
})

## bones uses dcat.  dcat works as of 2/25/22.
## However these tests are *very* slow for uncompiled,
## so I've only confirmed the first test passes.
## Also there are derivatives wrt variables that follow a dcat.
## In uncompiled execution, this generates endless warnings
## becuase x+delta is non-integer, but dcat returns only integers.
## In compiled execution it should would but I haven't checked.
## I am leaving these commented-out because they are so slow.
## dir = nimble:::getBUGSexampleDir('bones')
## Rmodel <- readBUGSmodel('bones', data = NULL, inits = NULL, dir = dir, useInits = TRUE,
##                         check = FALSE)
## test_ADModelCalculate_nick(Rmodel, name = 'bones', calcNodeNames = list(Rmodel$getDependencies('theta'), Rmodel$getDependencies('grade')),
##                            wrt = list(c('theta'), c('grade'), c('theta', 'grade'), c('grade[1, 2]')),
##                            order = c(0, 1, 2))

}
## end of Nick's tests ##



## basic model, with lifted nodes
set.seed(1)
inits <- list(mu0 = 1.2, tau = 1.5, tau0 = 2.2, mu = c(0.1, 1.1, 2.1))
data <- list(a = matrix(rnorm(6), 3, 2))
code <- nimbleCode({
        mu0 ~ dnorm(0.25, 1.25)
        tau0 ~ dgamma(1.5, 2.5)
        tau ~ dgamma(1.2, 2.3)
        for(j in 1:3) {
            for(i in 1:2) 
                a[j, i] ~ dnorm(mu[j], var = tau)
            mu[j] ~ dnorm(mu0, var = tau0)
        }
    })
model <- nimbleModel(code, inits = inits, data = data)
relTolTmp <- relTol
relTolTmp[3] <- 1e-4
relTolTmp[4] <- 1e-2

test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, useParamTransform = TRUE,
                      checkCompiledValuesIdentical = FALSE, name = 'basic model, lifted nodes')


## basic state space model
set.seed(1) 
code <- nimbleCode({
    x0 ~ dnorm(0,1)
    x[1] ~ dnorm(x0, 1)
    y[1] ~ dnorm(x[1], var = 2)
    for(i in 2:3) {
        x[i] ~ dnorm(x[i-1], 1)
        y[i] ~ dnorm(x[i], var = 2)
    }
})
data <- list(y = rnorm(3))
model <- nimbleModel(code, data = data)
model$simulate()
model$calculate()
relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-4
relTolTmp[4] <- 1e-2

test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, useParamTransform = TRUE,
                      checkCompiledValuesIdentical = FALSE, name = 'basic state space') 


## basic tricky indexing
set.seed(1)
code <- nimbleCode({
    y[1:2] ~ dmnorm(z[1:2], cov = covMat[,])
    z[1] <- x[2]
    z[2:3] <- x[1:2] + 1
    x[1] ~ dnorm(.1, 10)
    x[2] ~ dnorm(1, 3)
})
model <- nimbleModel(code, dimensions = list(x = 2, y = 2, z = 3), inits = list(covMat = matrix(c(1.9, .7, .7, 1.3), 2), x = c(1, 1.2)), data = list(y = c(-.1,-.2)))
relTolTmp <- relTol
relTolTmp[4] <- 1e-3

test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, name = 'basic tricky indexing', useParamTransform = TRUE,
                      newUpdateNodes = list(covMat = matrix(c(0.7, .25, .25, .7), 2)), checkCompiledValuesIdentical = FALSE)


## state space model
set.seed(1)
code <- nimbleCode({
  a ~ dunif(-0.9999, 0.9999)
  b ~ dnorm(0, sd = 1000)
  sigPN ~ dunif(1e-04, 1)
  sigOE ~ dunif(1e-04, 1)
  x[1] ~ dnorm(b/(1 - a), sd = sqrt(sigPN^2 + sigOE^2))
  y[1] ~ dnorm(x[1], sd = sigOE)
  for (i in 2:t) {
    x[i] ~ dnorm(x[i - 1] * a + b, sd = sigPN)
    y[i] ~ dnorm(x[i], sd = sigOE)
  }
})
constants <- list(t = 100)
data <- list(y = c(20.24405,20.57693,20.49357,20.34159,20.45759,20.43326,20.20554,20.12860,20.14756,20.20781,20.23022,20.26766,20.22984,20.37703,20.13641,20.05309,19.95709,20.19303,20.30562,20.54443,20.91010,20.70580,20.42344,20.19795,20.28816,20.31894,20.76939,20.77023,20.83486,20.29335,20.40990,20.19601,20.04083,19.76056,19.80810,19.83129,19.69174,19.90069,19.87623,19.63371,19.62360,19.72630,19.64450,19.86779,20.17104,20.34797,20.32968,20.48027,20.46694,20.47006,20.51676,20.40695,20.18715,19.97552,19.88331,19.67831,19.74702,19.47502,19.24408,19.37179,19.38277,19.15034,19.08723,19.37051,19.14274,19.46433,19.62459,19.77971,19.54194,19.39081,19.61621,19.51307,19.34745,19.17019,19.26829,19.58943,19.77143,19.83582,19.71198,19.67746,19.75053,20.40197,20.49363,20.37079,20.19005,20.55862,20.48523,20.33071,19.97069,19.79758,19.83811,19.79728,19.86277,19.86836,19.92481,19.88095,20.24899,20.55165,20.22707,20.11235))
inits <- list(a = 0.95, b=1, sigOE=0.05,sigPN = 0.2,  x= c(20.26036,20.51331,20.57057,20.35633,20.33736,20.47321,20.22002,20.14917,20.19216,20.26969,20.21135,20.22745,20.20466,20.41158,20.13408,20.08023,19.98956,20.13543,20.32709,20.55840,20.88206,20.74740,20.47671,20.14012,20.29953,20.33778,20.80916,20.75773,20.84349,20.35654,20.41045,20.20180,20.02872,19.74226,19.80483,19.81842,19.69770,19.84564,19.88211,19.70559,19.56090,19.73728,19.66545,19.88158,20.13870,20.39163,20.37372,20.47429,20.39414,20.42024,20.55560,20.40462,20.15831,19.89425,19.79939,19.72692,19.74565,19.42233,19.22730,19.36489,19.37289,19.19050,19.00823,19.35738,19.14293,19.48812,19.67329,19.82750,19.58979,19.43634,19.61278,19.56739,19.38584,19.19260,19.32732,19.65500,19.65295,19.84843,19.68285,19.69620,19.77497,20.31795,20.45797,20.32650,20.24045,20.60507,20.51597,20.30076,19.98100,19.86709,19.85965,19.74822,19.86730,19.90523,19.86970,19.87286,20.28417,20.46212,20.22618,20.13689))
model <- nimbleModel(code, constants = constants, data = data, inits = inits)
## large uncompiled 2d11 discrepancy occurs with new wrt for x near 0.
xNew <- list(x  = rnorm(constants$t, 20, 1))
## large uncompiled 2d11 discrepancy occurs with new Y (newConstantNodes) if not consistent with x values
## (simply seems to be that the combo of the inner deriv h and outer deriv h in numDeriv::jacobian are a bit too small)
newY <- rnorm(constants$t, 20, 1)

relTolTmp <- relTol
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-2
relTolTmp[4] <- 1e-2

test_ADModelCalculate(model, relTol = relTolTmp, xNew = xNew, newConstantNodes = list(y = newY), verbose = verbose, name = 'state space model', useFasterRderivs = TRUE)

## 2022-04-22: wrt values stored in model equal but not identical to 'x'
## 2022-04-22: various discrepancies
## Detected some values out of (relative, usually) tolerance:  rOutput12$hessian   cOutput12$hessian .
##                  [,1]                     [,2]                    [,3]
## [1,] 0.00195312500000 0.0019758109433068585595 0.011481839081673342329
## Detected some values out of (relative, usually) tolerance:  rOutput2d$value   cOutput2d$value .
##                  [,1]                     [,2]                    [,3]
## [1,] 0.00195312500000 0.0019758109433068585595 0.011481839081673342329
test_ADModelCalculate(model, relTol = relTolTmp, xNew = xNew, newConstantNodes = list(y = newY), verbose = verbose, name = 'state space model',
                    checkCompiledValuesIdentical = FALSE, useParamTransform = TRUE, useFasterRderivs = TRUE)

## 2022-04-22: both cases above
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##           [,1]       [,2]         [,3]
## [1,] -0.5344375 -0.5344375 1.329515e-14
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##            [,1]       [,2]       [,3]
## [1,] -0.1628140 -0.1545240 0.05364862




## dexp and dt, which are provided by NIMBLE to allow expanded parameterizations
set.seed(1)
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(mu[i], sd = sigma)
        mu[i] ~ dt(mu0, sigma = sd0, df = nu)
    }
    mu0 ~ dnorm(0,1)
    sd0 ~ dgamma(1.5, .7)
    nu ~ dunif(0, 50)
    sigma ~ dexp(scale = tau)
    tau ~ dexp(rate = 3.5)
})
n <- 10
## if mu is near zero, can get back uncompiled 2d11 hessian (NCT issue 350 likely)
model <- nimbleModel(code, constants = list(n = n), data = list(y = rnorm(n)),
                     inits = list(mu = rnorm(n, 1, .2), sigma = runif(1), nu = 2.5, mu0 = rnorm(1), sd0 = runif(1),
                                  tau = runif(1)))
xNew = list(mu = rnorm(n, 1, .2))
            
relTolTmp <- relTol
relTolTmp[1] <- 1e-10  
relTolTmp[2] <- 1e-7 
relTolTmp[3] <- 1e-4 
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13

test_ADModelCalculate(model, relTol = relTolTmp, xNew = xNew, verbose = verbose, useFasterRderivs = TRUE, name = 'dt and dexp model')

## 2022-04-08: wrt values stored in model equal but not identical to 'x'
test_ADModelCalculate(model, relTol = relTolTmp, xNew = xNew, verbose = verbose, useFasterRderivs = TRUE, useParamTransform = TRUE,
                      checkCompiledValuesIdentical = FALSE, name = 'dt and dexp model')

## 2022-04-21 (both cases above): 
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##             [,1]        [,2]         [,3]
## [1,] -0.12124044 -0.12124044 3.319485e-15
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##             [,1]        [,2]       [,3]
## [1,] -0.01855928 -0.01834029 0.01194040



## vectorized deterministic nodes

set.seed(1)
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dpois(mu[i])
        logmu[i] ~ dnorm(0, 1)
    }
    mu[1:n] <- exp(logmu[1:n])
})
n <- 10
model <- nimbleModel(code, constants = list(n = n), data = list(y = rpois(n, 1)),
                     inits = list(logmu = rnorm(n)))
newY <- rpois(n, 2)

relTolTmp <- relTol
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2
test_ADModelCalculate(model, newConstantNodes = list(y = newY), relTol = relTolTmp, verbose = verbose, name = 'deterministic vectorized model')
test_ADModelCalculate(model, newConstantNodes = list(y = newY), useParamTransform = TRUE, relTol = relTolTmp,
                      checkCompiledValuesIdentical = FALSE, verbose = verbose, name = 'deterministic vectorized model')
## 2022-04-12: all set


## truncation
## Note that constraints are not handled
set.seed(1)
code <- nimbleCode({
    for(i in 1:n) {
        y[i, 1:4] ~ dmnorm(mu[1:4], pr[1:4,1:4])
    }
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
    ## constraint ~ dconstraint(mu[1] + mu[2] > 0)  ## not handled (dconstraint has `if` in its density)
    mu[3] ~ T(dnorm(0,1), -3, 3)
    mu[4] ~ T(dnorm(0,1), 0, )
})
n <- 10
inits <- list(mu = c(0.35, -0.25, 1.3, 2.7), pr = diag(rep(1,4)))
y <- matrix(0, n, 4)
for(i in 1:n)
    y[i, ] <- rmnorm_chol(1, inits$mu, diag(rep(0.2, 4)), prec_param = FALSE)
model <- nimbleModel(code, constants = list(n = n), data = list(y = y), inits = inits)
newPr <- crossprod(matrix(rnorm(4*4), 4))

relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[4] <- 1e-2
## 2022-04-16: all set
test_ADModelCalculate(model, newUpdateNodes = list(pr = newPr), relTol = relTolTmp, verbose = verbose, checkCompiledValuesIdentical = FALSE, name = 'truncation model')

relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13
## 2022-04-16: various values equal but not identical
## Detected some values out of relative tolerance:  cOutput2d$value   c(cOutput012$hessian) .
## [1] 8.466706e-01 8.466706e-01 1.311281e-15
test_ADModelCalculate(model, newUpdateNodes = list(pr = newPr), useParamTransform = TRUE, relTol = relTolTmp, checkCompiledValuesIdentical = FALSE, verbose = verbose, name = 'truncation model')

if(FALSE) {   ## truncation on non-top nodes not handled for now; issue #254
    code <- nimbleCode({
        y ~ dnorm(mu, 1)
        mu ~ T(dnorm(mu0, 1), 0, 10)
        mu0 ~ dnorm(0,1)
    })
    model <- nimbleModel(code, data = list(y = 1), inits = list(mu = 0.5, mu0 = 1))
    test_ADModelCalculate(model, relTol = relTol, verbose = verbose, name = 'truncation on non-top node')
    
    
    code <- nimbleCode({
        y ~ dnorm(mu, 1)
        mu ~ T(dbeta(a,b), 0.1, 0.8)
        a ~ dunif(0, 5)
        b ~ dunif(0,5 )
    })
    model <- nimbleModel(code, data = list(y = 1), inits = list(mu = 0.5,a=1,b=1))
    test_ADModelCalculate(model, relTol = relTol, verbose = verbose, name = 'truncation with dbeta')
}

## complicated indexing 
set.seed(1)
code <- nimbleCode({
    x[2:4,3:5] <- S[1:3,1:3]
    x[1,2:5] ~ dmnorm(z[1:4], pr2[1:4,1:4])
    alphas[1:4] <- exp(x[1, 1:4]) / sum(exp(x[1,1:4]))
    x[2:5, 2] <- alphas[1:4]
    x[1,1] ~ dgamma(3,1.7)
    for(i in 2:5)
        x[i, 1] ~ dgamma(1, 1)
    for(i in 3:5)
        x[5, i] ~ dnorm(0, 1)
    w[1:3] ~ dmnorm(z[1:3], cov = x[2:4, 3:5])
    for(j in 1:5)
        y[1:5, j] ~ dmnorm(x[j, 1:5], cov = pr[1:5,1:5])
})
inits <- list(S = diag(rep(1, 3)), z = rep(0, 4), pr = diag(5), pr2 = diag(4))
model <- nimbleModel(code, inits = inits)
model$simulate()
model$calculate()
model$setData('y','w')
newPr <- crossprod(matrix(rnorm(5*5), 5))
newPr2 <- crossprod(matrix(rnorm(4*4), 4))
newS <- crossprod(matrix(rnorm(3*3), 3))

relTolTmp <- relTol
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-2
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13

test_ADModelCalculate(model, newUpdateNodes = list(S = newS, pr = newPr, pr2 = newPr2), useParamTransform = TRUE, relTol = relTolTmp, useFasterRderivs = TRUE, checkCompiledValuesIdentical = FALSE, verbose = verbose, name = 'complicated indexing')


## using different subsets of a matrix
set.seed(1)
code <- nimbleCode({
    x1[1:5] ~ dmnorm(z[1:5], pr5[1:5, 1:5])
    x2[1:4] ~ dmnorm(z[1:4], pr4[1:4, 1:4])
    for(i in 1:5)
        y1[i] ~ dnorm(x1[i], 1)
    for(i in 1:4)
        y2[i] ~ dnorm(x2[i], 1)
})
inits <- list(z = rep(0,5), pr5 = diag(5), pr4 = diag(4))
model <- nimbleModel(code, inits = inits)
model$simulate()
model$calculate()
model$setData('y1', 'y2')
newPr5 <- crossprod(matrix(rnorm(5*5), 5))
newPr4 <- crossprod(matrix(rnorm(4*4), 4))
relTolTmp <- relTol
relTolTmp[2] <- 1e-5
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2

test_ADModelCalculate(model, newUpdateNodes = list(pr5 = newPr5, pr4 = newPr4), useParamTransform = TRUE, checkCompiledValuesIdentical = TRUE, relTol = relTolTmp, useFasterRderivs = TRUE, verbose = verbose, name = 'different subsets of a matrix')

## vectorized covariance matrix
set.seed(1)
code <- nimbleCode({
    Sigma1[1:n,1:n] <- exp(-dist[1:n,1:n]/rho)
    y[1:n] ~ dmnorm(mu1[1:n], cov = Sigma1[1:n,1:n])
    mu1[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    rho ~ dgamma(2, 3)
})

n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(dist = dd, rho = rgamma(1, 1, 1), z = rep(1, n), pr = diag(n)))
model$simulate()
model$calculate()
model$setData('y')
newPr <- crossprod(matrix(rnorm(5*5), 5))
newDist <- as.matrix(dist(runif(5)))

relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-3
relTolTmp[5] <- 1e-13

test_ADModelCalculate(model, newUpdateNodes = list(pr = newPr, dist = newDist),
                      checkCompiledValuesIdentical = FALSE, relTol = relTolTmp, absTolThreshold = 1e-12, verbose = verbose, name = 'dmnorm with vectorized covariance matrix')
## 2022-04-20: various compiled values equal but not identical
test_ADModelCalculate(model, newUpdateNodes = list(pr = newPr, dist = newDist), useParamTransform = TRUE,
                      checkCompiledValuesIdentical = FALSE, relTol = relTolTmp, absTolThreshold = 1e-12, verbose = verbose, name = 'dmnorm with vectorized covariance matrix')

## 2022-04-22: (both cases above)
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##            [,1]       [,2]         [,3]
## [1,]  0.0968665  0.0968665 1.575939e-15


## vectorized covariance matrix, chol param
set.seed(1)
code <- nimbleCode({
    Sigma1[1:n,1:n] <- exp(-dist[1:n,1:n]/rho)
    y[1, 1:n] ~ dmnorm(mu1[1:n], cov = Sigma1[1:n,1:n])
    mu1[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    mu2[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    Ucov[1:n, 1:n] <- chol(Sigma1[1:n,1:n])
    y[2, 1:n] ~ dmnorm(mu2[1:n], cholesky = Ucov[1:n,1:n], prec_param = 0)
    rho ~ dgamma(2, 3)
})

n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)

## uncompiled 2d11 non-trivial discrepancies seem to correspond with near zero x values (likely NCT issue 350)
## so set z not equal to zero so don't have small mu1, mu2 values
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(dist = dd, rho = rgamma(1, 1, 1),
                                  z = rep(1, n), pr = diag(n)))
model$simulate()
model$calculate()
model$setData('y')
newPr <- crossprod(matrix(rnorm(5*5), 5))
newDist <- as.matrix(dist(runif(5)))

relTolTmp <- relTol
relTolTmp[1] <- 1e-10
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13

test_ADModelCalculate(model, newUpdateNodes = list(pr = newPr, dist = newDist),
                      relTol = relTolTmp, absTolThreshold = 1e-12, verbose = verbose, name = 'dmnorm with vectorized covariance matrix, chol param')

## 2022-04-20: various compiled values equal but not identical
test_ADModelCalculate(model, newUpdateNodes = list(pr = newPr, dist = newDist), useParamTransform = TRUE, checkCompiledValuesIdentical = FALSE,
                      relTol = relTolTmp, absTolThreshold = 1e-12, verbose = verbose, name = 'dmnorm with vectorized covariance matrix, chol param')

## 2022-04-21: (both cases above) various minor discrepancies, e.g.,
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##           [,1]      [,2]         [,3]
## [1,] 0.2314997 0.2314997 1.918315e-15
## [2,] 0.2314997 0.2314997 1.918315e-15
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##            [,1]       [,2]        [,3]
## [1,] 0.03049935 0.03059135 0.003007278



## MVN with various parameterizations and user-defined functions

## user-defined cov function with loops
## works with work-around that forces index vars to be treated as ints; see NCT issue 130
covFunLoop <- nimbleFunction(
    run = function(dist = double(2), rho = double(0)) {
        n = dim(dist)[1]
        out = nimMatrix(nrow = n, ncol = n)
        ## i <- 1L; j <- 1L  ## NCT 130 work-around
        for(i in 1:n)
            for(j in 1:n)
                out[i,j] <- exp(-dist[i,j]/rho)
        returnType(double(2))
        return(out)
    }, buildDerivs = list(run = list(ignore = c('i','j'))))
    # buildDerivs = TRUE) # if use NCT 130 work-around



## user-defined cov function vectorized
set.seed(1)
covFunVec <- nimbleFunction(
    run = function(dist = double(2), rho = double(0)) {
        out <- exp(-dist/rho)
        returnType(double(2))
        return(out)
    }, buildDerivs = TRUE)

code <- nimbleCode({
    Sigma3[1:n,1:n] <- covFunVec(dist[1:n,1:n], rho)
    y[1:n] ~ dmnorm(mu3[1:n], cov = Sigma3[1:n,1:n])
    mu3[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    rho ~ dgamma(2, 3)
})
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(dist = dd, rho = rgamma(1, 1, 1),
                                  z = rep(0, n), pr = diag(n)))
newPr <- crossprod(matrix(rnorm(n*n),n))
newDist <- as.matrix(dist(runif(n)))

model$simulate()
model$calculate()
model$setData('y')
relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2

## 2022-04-22: various compiled values equal but not identical
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTolTmp, absTolThreshold = 1e-12, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      newUpdateNodes = list(dist = newDist, pr = newPr),
                      name = 'dmnorm with user-defined vectorized fxn')

## 2022-04-23: both cases above
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##             [,1]        [,2]         [,3]
## [1,] -0.42237762 -0.42237762 2.891359e-15



## various dmvt parameterizations
set.seed(1)
code <- nimbleCode({
    y[1, 1:n] ~ dmvt(mu = mu1[1:n], prec = Q[1:n,1:n], df = nu)

    Uprec[1:n, 1:n] <- chol(Q[1:n,1:n])
    Ucov[1:n, 1:n] <- chol(Sigma[1:n,1:n])
    y[2, 1:n] ~ dmvt(mu2[1:n], cholesky = Uprec[1:n,1:n], df = nu, prec_param = 1)
    y[3, 1:n] ~ dmvt(mu3[1:n], cholesky = Ucov[1:n,1:n], df = nu, prec_param = 0)
    y[4, 1:n] ~ dmvt(mu4[1:n], scale = Sigma[1:n, 1:n], df = nu)
    mu1[1:n] ~ dmvt(z[1:n], prec = pr[1:n,1:n], df = nu)
    mu2[1:n] ~ dmvt(z[1:n], prec = pr[1:n,1:n], df = nu)
    mu3[1:n] ~ dmvt(z[1:n], prec = pr[1:n,1:n], df = nu)
    mu4[1:n] ~ dmvt(z[1:n], prec = pr[1:n,1:n], df = nu)
    nu ~ dunif(0, 50)
})

n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
Sigma <- exp(-dd/0.1)
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(z = rep(1, n), pr = diag(n),
                                  Sigma = Sigma, Q = solve(Sigma)))
newSigma <- crossprod(matrix(rnorm(n*n), n))
newQ <- crossprod(matrix(rnorm(n*n), n))

model$simulate()
model$calculate()
model$setData('y')

relTolTmp <- relTol
relTolTmp[1] <- 1e-10
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-2
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13

test_ADModelCalculate(model,
                      newUpdateNodes = list(Q = newQ, Sigma = newSigma, pr = newPr),
                      checkCompiledValuesIdentical = FALSE, relTol = relTolTmp, verbose = verbose, name = 'various dmvt parameterizations')

## 2022-04-22: various compiled values equal but not identical
test_ADModelCalculate(model, useParamTransform = TRUE,checkCompiledValuesIdentical = FALSE,
                      newUpdateNodes = list(Q = newQ, Sigma = newSigma, pr = newPr),
                      relTol = relTolTmp, verbose = verbose, name = 'various dmvt parameterizations')

## 2022-04-23: both cases above
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##             [,1]        [,2]         [,3]
## [1,]  0.01955694  0.01955694 1.930137e-13
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##           [,1]      [,2]       [,3]
## [1,] -3.348612 -3.248763 0.03073416
## Detected some values out of (relative, usually) tolerance:  rOutput012$hessian   cOutput012$hessian .
##               [,1]          [,2]       [,3]
## [1,] -4.887581e-06 -4.809407e-06 0.01625429



## dirichlet as likelihood so not differentiating wrt something with constraint.
set.seed(1)
code <- nimbleCode({
    p[1:k] ~ ddirch(alpha[1:k])
    for(i in 1:k)
        alpha[i] ~ dgamma(1.3, 1.5)
})
k <- 4
model <- nimbleModel(code, constants = list(k = k), data = list(p = c(.2, .4, .15, .25)), inits = list(alpha = runif(4)))
newP <- rdirch(1, rep(1,4))

relTolTmp <- relTol
relTolTmp[1] <- 1e-12
relTolTmp[2] <- 1e-7
relTolTmp[4] <- 1e-3

## 2022-04-23: all set
test_ADModelCalculate(model,
                      newConstantNodes = list(p = newP),
                      relTol = relTolTmp, verbose = verbose, name = 'Dirichlet likelihood')


## 2022-04-22: various compiled values equal but not identical
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
## [1] -5.8683571e-02 -5.8683571e-02  3.9020035e-15
test_ADModelCalculate(model, useParamTransform = TRUE,
                      newConstantNodes = list(p = newP),
                      relTol = relTolTmp, verbose = verbose, name = 'Dirichlet likelihood')



## dwish and dinvwish so long as not differentiating w.r.t. the random variable (since it has constraints)
## Note that nu must exceed n and can't be set to runif(0,1) via test_ADModelCalculate.
set.seed(1)
code <- nimbleCode({
    R[1:n,1:n] <- sigma2 * exp(-dist[1:n, 1:n] / rho)
    US[1:n,1:n] <- chol(inverse(R[1:n,1:n]))
    UR[1:n,1:n] <- chol(R[1:n,1:n])
    W1[1:n, 1:n] ~ dwish(R[1:n, 1:n], nu + n)
    W2[1:n, 1:n] ~ dwish(S = R[1:n, 1:n], df = nu + n)
    W3[1:n, 1:n] ~ dwish(cholesky = UR[1:n,1:n], df = nu + n, scale_param = 0)
    W4[1:n, 1:n] ~ dwish(cholesky = US[1:n,1:n], df = nu + n, scale_param = 1)
    IW1[1:n, 1:n] ~ dinvwish(R[1:n, 1:n], nu + n)
    IW2[1:n, 1:n] ~ dinvwish(R = R[1:n, 1:n], df = nu + n)
    IW3[1:n, 1:n] ~ dinvwish(cholesky = UR[1:n,1:n], df = nu + n, scale_param = 0)
    IW4[1:n, 1:n] ~ dinvwish(cholesky = US[1:n,1:n], df = nu + n, scale_param = 1)
    rho ~ dgamma(2, 3)
    sigma2 ~ dgamma(2, 3)
    nu ~ dgamma(2, 3)
})
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n), inits = list(dist = dd, nu = 5))
model$simulate()
model$calculate()
model$setData(c('W1','W2','W3','W4','IW1','IW2','IW3','IW4'))

newW1 <- crossprod(matrix(rnorm(5*5), 5))
newW2 <- crossprod(matrix(rnorm(5*5), 5))
newW3 <- crossprod(matrix(rnorm(5*5), 5))
newW4 <- crossprod(matrix(rnorm(5*5), 5))
newIW1 <- crossprod(matrix(rnorm(5*5), 5))
newIW2 <- crossprod(matrix(rnorm(5*5), 5))
newIW3 <- crossprod(matrix(rnorm(5*5), 5))
newIW4 <- crossprod(matrix(rnorm(5*5), 5))

relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2

test_ADModelCalculate(model, newUpdateNodes = list(dist = newDist),
                      newConstantNodes = list(W1 = newW1, W2 = newW2, W3 = newW3, W4 = newW4,
                                        IW1 = newIW1, IW2 = newIW2, IW3 = newIW3, IW4 = newIW4),
                      relTol = relTolTmp, verbose = verbose, name = 'dwish, dinvwish')

## 2022-04-22: various compiled values equal but not identical
test_ADModelCalculate(model, useParamTransform = TRUE, newUpdateNodes = list(dist = newDist),
                      newConstantNodes = list(W1 = newW1, W2 = newW2, W3 = newW3, W4 = newW4,
                                        IW1 = newIW1, IW2 = newIW2, IW3 = newIW3, IW4 = newIW4),
                      relTol = relTolTmp, verbose = verbose, name = 'dwish, dinvwish')

## Both of the above:
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##           [,1]      [,2]         [,3]
## [1,]  192.0026  192.0026 2.960554e-15
## Only with paramTransform:
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
## [1] -2.37781548 -2.41502054  0.01540569




## simple user-defined distribution
dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0, default = 1), 
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- log(rate) - x*rate
        if(log) return(logProb)
        else return(exp(logProb)) 
    }, buildDerivs = TRUE)

code <- nimbleCode({ 
    y ~ dmyexp(rho)
    rho ~ dgamma(2, 3)
})

set.seed(1)
model <- nimbleModel(code, data = list(y = rgamma(1,1,1)), inits = list(rho = rgamma(1, 1, 1)))
relTolTmp <- relTol
relTolTmp[1] <- 1e-14

test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, name = 'simple user-defined distribution')
## 2022-04-22: various compiled values equal but not identical
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTolTmp, verbose = verbose, name = 'simple user-defined distribution')

dtest <- nimbleFunction(
    run = function(x = double(0), mu = double(1), 
                   log = integer(0, default = 0)) {
        returnType(double(0))
        tmp <- mu^2
        out <- tmp[1]
        if(log) return(out) else return(exp(out))
    }, buildDerivs = TRUE)

rtest <- nimbleFunction(
    run = function(n = integer(0), mu = double(1)) {
        returnType(double(0))
        out <- rnorm(1, 0, 1)
        return(out)
    })

code <- nimbleCode({ 
    y ~ dtest(mu[1:2])
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
})

set.seed(1)
model <- nimbleModel(code, inits = list(mu = rnorm(2)))
model$simulate()
model$calculate()
model$setData('y')

relTolTmp <- relTol
relTolTmp[1] <- 1e-14

test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, name = 'vectorized power')
## 2022-04-22: various compiled values equal but not identical
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTolTmp, verbose = verbose, name = 'vectorized power')


dGPdist <- nimbleFunction(
    run = function(x = double(1), dist = double(2), rho = double(0),
                   log = integer(0, default = 0)) {
        returnType(double(0))
        Sigma <- exp(-dist/rho)
        L <- t(chol(Sigma))
        out <- -sum(log(diag(L)))
        qf <- forwardsolve(L, x)
        out <- out - 0.5*sum(qf^2)
        if(log) return(out) else return(exp(out))
    }, buildDerivs = TRUE)

rGPdist <- nimbleFunction(
    run = function(n = integer(0), dist = double(2), rho = double(0)) {
        returnType(double(1))
        Sigma <- exp(-dist/rho)
        U <- chol(Sigma)
        p <- dim(dist)[1]
        out <- (U %*% rnorm(p))[,1]
        return(out)
    })



## Dirichlet
code <- nimbleCode({
    y[1:k] ~ dmulti(p[1:k], n)
    p[1:k] ~ ddirch(alpha[1:k])
    for(i in 1:k)
        alpha[i] ~ dgamma(1.3, 1.5)
})
n <- 30
k <- 4
set.seed(1)
model <- nimbleModel(code, constants = list(k = k, n = n), data = list(y = rmulti(1, n, rep(1/k, k))), inits = list(p = c(.2, .4, .15, .25), alpha = runif(4)))
newP <- rdirch(1, rep(1:k))
newY <- rmulti(1, n, newP)
relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-2
relTolTmp[4] <- 1e-2

test_ADModelCalculate(model, x = 'prior', useParamTransform = TRUE, newUpdateNodes = list(p = newP), newConstantNodes = list(y = newY),
                      relTol = relTolTmp,absTolThreshold = 1e-12, verbose = verbose, name = 'Dirichlet paramTransform') 
## 2022-04-25: various cases not identical
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
## [1] -7.310951e-05 -7.310951e-05  4.145310e-12
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##               [,1]          [,2]       [,3]
## [1,] -1.800050e-05 -1.957906e-06 8.19375019
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##               [,1]          [,2]       [,3]
## [1,] -3.077516e-05 -3.311835e-05 0.07075190
## Detected some values out of (relative, usually) tolerance:  rOutput12$hessian   cOutput12$hessian .
##               [,1]          [,2]       [,3]
## [1,] -3.576279e-06 -3.778980e-06 0.05363926


code <- nimbleCode({
    Sigma1[1:n,1:n] <- exp(-dist[1:n,1:n]/rho)
    Sigma2[1:n,1:n] <- covFunLoop(dist[1:n,1:n], rho)
    Sigma3[1:n,1:n] <- covFunVec(dist[1:n,1:n], rho)
    
    y[1, 1:n] ~ dmnorm(mu1[1:n], cov = Sigma1[1:n,1:n])
    y[2, 1:n] ~ dmnorm(mu2[1:n], cov = Sigma2[1:n,1:n])
    y[3, 1:n] ~ dmnorm(mu3[1:n], cov = Sigma3[1:n,1:n])

    Q[1:n,1:n] <- inverse(Sigma1[1:n, 1:n])
    y[4, 1:n] ~ dmnorm(mu4[1:n], Q[1:n,1:n])

    Uprec[1:n, 1:n] <- chol(Q[1:n,1:n])
    Ucov[1:n, 1:n] <- chol(Sigma1[1:n,1:n])
    y[5, 1:n] ~ dmnorm(mu5[1:n], cholesky = Uprec[1:n,1:n], prec_param = 1)
    y[6, 1:n] ~ dmnorm(mu6[1:n], cholesky = Ucov[1:n,1:n], prec_param = 0)
    y[7, 1:n] ~ dGPdist(dist[1:n, 1:n], rho)
    
    W1[1:n, 1:n] ~ dinvwish(R = R[1:n,1:n], df = nu)

    UR[1:n, 1:n] <- chol(R[1:n,1:n])
    US[1:n, 1:n] <- chol(inverse(R[1:n,1:n]))
    W2[1:n, 1:n] ~ dinvwish(cholesky = UR[1:n,1:n], df = nu, scale_param = 0)
    W3[1:n, 1:n] ~ dinvwish(cholesky = US[1:n,1:n], df = nu, scale_param = 1)

    W4[1:n, 1:5] ~ dwish(R[1:n, 1:n], df = nu)
    W5[1:n, 1:5] ~ dwish(cholesky = UR[1:n, 1:n], df = nu, scale_param = 0)
    W6[1:n, 1:5] ~ dwish(cholesky = US[1:n, 1:n], df = nu, scale_param = 1)
    
    mu1[1:n] ~ dmnorm(z[1:n], cov = W1[1:n,1:n])
    mu2[1:n] ~ dmnorm(z[1:n], cov = W2[1:n,1:n])
    mu3[1:n] ~ dmnorm(z[1:n], cov = W3[1:n,1:n])
    mu4[1:n] ~ dmnorm(z[1:n], W4[1:n,1:n])
    mu5[1:n] ~ dmnorm(z[1:n], W5[1:n,1:n])
    mu6[1:n] ~ dmnorm(z[1:n], W6[1:n,1:n])
    rho ~ dgamma(2, 3)
    nu ~ dunif(0, 100)
})

set.seed(1)
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
R <- crossprod(matrix(rnorm(n^2), n, n))
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(dist = dd, R = R, nu = 8, rho = rgamma(1, 1, 1),
                                                                 z = rep(1, n)))
model$simulate()
model$calculate()
model$setData('y')

newDist <- as.matrix(dist(runif(n)))
newR <- crossprod(matrix(rnorm(n*n), n))
newW1 <- crossprod(matrix(rnorm(n*n), n))
newW2 <- crossprod(matrix(rnorm(n*n), n))
newW3 <- crossprod(matrix(rnorm(n*n), n))
newW4 <- crossprod(matrix(rnorm(n*n), n))
newW5 <- crossprod(matrix(rnorm(n*n), n))
newW6 <- crossprod(matrix(rnorm(n*n), n))
relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-2

test_ADModelCalculate(model, newUpdateNodes = list(nu = 12.1, dist = newDist, R = newR, W1 = newW1, W2 = newW2, W3 = newW3, W4 = newW4, W5 = newW5, W6 = newW6),
                      x = 'prior', absTolThreshold = 1e-12,
                      useParamTransform = TRUE, useFasterRderivs = TRUE,
                      relTol = relTolTmp, verbose = verbose,
                      name = 'various multivariate dists')

## 2022-05-05:
## I can't fully explain the magnitude of the 2d11 disparity. If I manually call jacobian, I can see
## smaller disparities than if I do my own hand-calculated 2nd derivs, but not the magnitude seen here.
## However, given the single-taped uncompiled Hessians match compiled Hessians, I don't think anything is wrong.
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##              [,1]          [,2]      [,3]
## [1,]  1.918025622  0.0000000000 1.9180256
## Detected some values out of (relative, usually) tolerance:  rOutput2d11$jacobian   cOutput2d11$jacobian .
##           [,1] [,2]      [,3]
## [1,] 12.198946    0 12.198946
## [2,] 12.198946    0 12.198946
## [3,] -6.777284    0  6.777284
## [4,] -6.777284    0  6.777284
## [5,]  2.709502    0  2.709502
## Detected some values out of (relative, usually) tolerance:  cOutput2d$value   c(cOutput012$hessian) .
##            [,1]       [,2]         [,3]
## [1,] -0.9472489 -0.9472489 1.200179e-13
## Detected some values out of  absolute  tolerance   :  rOutput012$hessian   cOutput012$hessian .
##              [,1]         [,2]      [,3]
## [1,] -0.002075195 -0.002125045 0.0234584
## Detected some values out of  absolute  tolerance   :  rOutput01$jacobian   cOutput01$jacobian .
## [1] -4.988689e-01 -4.988643e-01  9.101791e-06

## various cases equal but not identical

## 2022-07-29: have gotten seg faults in EB scenario, 'Testing new wrt values with new constantNodes' for the above when running on its own, 

## Mimic BUGS jaws model, given difficulty in getting Hessians wrt Wishart params to match in actual jaw, presumably simply because of the numerical difficulties with values being used.

code <- nimbleCode({
  for (i in 1:N) {
     Y[i,1:M] ~ dmnorm(mu[1:M], Omega[1:M,1:M]);  # The 4 measurements for each  
  }                                   # boy are multivariate normal

  for(j in 1:M) {     # location model for mean bone length at each age
     mu[j] <- beta0;                                         # constant
  }

  beta0 ~ dnorm(0.0, 0.001); 
  Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 6);    # between-child variance in length at each age  
})

set.seed(1)
M <- 4
m <- nimbleModel(code, constants = list(M=M, N = 20), inits = list(R=diag(4)))
m$beta0 <- 45
m$Omega <- rwish_chol(1, chol(m$R), 4)
m$simulate('mu')
m$simulate(m$getDependencies('Omega'), includeData = TRUE)
m$setData('Y')

newR <- crossprod(matrix(rnorm(M*M), M))
newOmega <- solve(cov(m$Y))

relTolTmp <- relTol
relTolTmp[2] <- 1e-5
relTolTmp[3] <- 1e-4
relTolTmp[4] <- 1e-3
test_ADModelCalculate(m, xNew = list(beta0 = 42, Omega = newOmega), newUpdateNodes = list(R = newR, Omega = newOmega), 
                      absTolThreshold = 1e-10,
                      useParamTransform = TRUE, useFasterRderivs = TRUE,
                      relTol = relTolTmp, verbose = verbose,
                      name = 'jaw facsimile')

## 2022-08-02: worst results are various compiled vs. uncompiled Hessians off by 0.005 - 0.009

## Mimic schools model, given difficulty in getting Hessians wrt Wishart params to match in actual jaw, presumably simply because of the numerical difficulties with values being used.

code <- nimbleCode({
   for(p in 1:N) {
      Y[p] ~ dnorm(mu[p], tau[p])
      mu[p] <- alpha[school[p],1] + alpha[school[p],2]*LRT[p] + alpha[school[p],3]*VR[p,1] + beta[1]*LRT2[p] +
          beta[2]*VR[p,2] + beta[3]*Gender[p] +
               beta[4]*School.gender[p,1] + beta[5]*School.gender[p,2]+
               beta[6]*School.denom[p,1] + beta[7]*School.denom[p,2]+
               beta[8]*School.denom[p,3]
      log(tau[p]) <- theta + phi*LRT[p]
      LRT2[p] <- LRT[p]^2
   }

   # Priors for fixed effects:
   for (k in 1:8)
       beta[k] ~ dnorm(0.0, 0.0001)   
   theta ~ dnorm(0.0, 0.0001)
   phi ~ dnorm(0.0, 0.0001)

   # Priors for random coefficients:
   for (j in 1:M) {
      alpha[j,1:D] ~ dmnorm(gamma[1:D], T[1:D,1:D]) 
   }
 
   # Hyper-priors:
   gamma[1:D] ~ dmnorm(mn[1:D], prec[1:D,1:D])

   T[1:D,1:D] ~ dwish(R[1:D,1:D],3)
})

e <- new.env()
source(system.file('classic-bugs','vol2','schools','schools-data.R', package = 'nimble'), local = e)
e$D <- 3
#e$LRT <- rnorm(length(e$LRT))  # otherwise get crazy taus after exp()
rm('Y', envir = e)

set.seed(1)
m <- nimbleModel(code, constants = as.list(e))                                        
m$theta <- .1
m$phi <- .1
m$beta <- rnorm(8)
m$beta[1] <- 0.01
m$T <- rwish_chol(1, chol(m$R), 3)
m$gamma <- rnorm(3)
m$alpha <- matrix(rnorm(e$M*e$D),e$M)
m$simulate(m$getDependencies(c('theta','phi','beta'), self = FALSE, downstream = TRUE, includeData = TRUE))
m$setData('Y')

newR <- matrix(c(.19,-.004, .004, -.004, .15, .01, .004, .01, .07), 3, 3)
newT <- solve(cov(m$alpha))
newPrec <- solve(matrix(c(9000, 100, 100, 100, 9000, 100, 100, 100, 9000), 3))
newY <- rnorm(length(m$Y), m$Y, .5) 

relTolTmp <- relTol
relTolTmp[2] <- 1e-4
relTolTmp[3] <- 1e-4
relTolTmp[3] <- 1e-4
test_ADModelCalculate(m, xNew = list(gamma = rnorm(3, m$gamma, .1), theta = .15, phi = .07, beta = m$beta + c(.005, rnorm(7, 0, .2)),
                                     alpha = rnorm(e$D*e$M, m$alpha, 0.2), T = newT), newUpdateNodes = list(T = newT, R = newR),
                      newConstantNodes = list(R = newR, prec = newPrec, Y = newY), 
                      absTolThreshold = 1e-11,
                      useParamTransform = TRUE, useFasterRderivs = TRUE,
                      relTol = relTolTmp, verbose = verbose,
                      name = 'schools facsimile')

## 2022-08-02: various bad Jacobians and Hessians involving parameter T.
## E.g., Hessian relative errors of 0.02, 0.21, 210, 52
## Verified that for three specific cases (.02,.21,210) that using numDeriv and/or modifying 'h'
## in pracma:hessian gives relative errors of ~1e-5 for first two and .00078 for third.
## So it appears that even for single-taped case, can get very bad numerical derivative estimates.

## loop through BUGS models

## example test of a BUGS model - take this out as epil is in loop below
if(FALSE) {
    model <- readBUGSmodel(model = 'epil2.bug', inits = 'epil-inits.R', data = 'epil-data.R',
                           useInits = TRUE, dir = nimble:::getBUGSexampleDir('epil'))
    set.seed(1)
    model$simulate('b1')
    model$calculate()
    newY <- matrix(rpois(length(model$y), 2), nrow = nrow(model$y))
    test_ADModelCalculate(model, newConstantNodes = list(y = newY),
                          relTol = relTol, verbose = verbose, name = 'epil', useFasterRderivs = TRUE)
}


## now loop through BUGS models.
set.seed(1)
examples <- c('blocker', 'dyes', 'epil', 'equiv', 'line', 'pump', 'rats', 'beetles', 'jaw', 'dugongs', 'seeds','oxford','schools')
bugsFile <- examples
initsFile <- dataFile <- rep(NA, length(examples))
names(bugsFile) <- names(initsFile) <- names(dataFile) <- examples

## litters uses truncation on non-top node (issue #254) so left out for now
## bones, inhaler, pigs has dcat so left out for now (but we should be able to handle this)
## biops has stoch indexing so left out for now
## lsat requires setting some indexing and doing a bunch of initialization
## kidney, mice has dinterval so left out for now
## asia has dcat and max
## stagnant has dcat and step
## eyes has dnormmix
## TODO: could add: air, alli, birats, cervix, hearts, ice, orange, but these have similar features to models already checked.

## customize file names as needed
bugsFile['beetles'] <- 'beetles-logit'
bugsFile['jaw'] <- 'jaw-constant'
bugsFile['epil'] <- 'epil2'
bugsFile <- paste0(bugsFile, '.bug')

## initsFile['epil'] <- 'epil-inits.R'
dataFile['epil'] <- 'epil-data.R'
initsFile['jaw'] <- 'jaw-inits.R'
dataFile['jaw'] <- 'jaw-data.R'
initsFile['beetles'] <- 'beetles-inits.R'
dataFile['beetles'] <- 'beetles-data.R'
initsFile['schools'] <- 'schools-inits.R'

## customize whether any nodes need to be initialized
simulateNodes <- list()
length(simulateNodes) <- length(examples)
names(simulateNodes) <- examples
simulateNodes[['epil']] <- 'b1'
simulateNodes[['dyes']] <- 'mu'
simulateNodes[['equiv']] <- 'd'
simulateNodes[['oxford']] <- 'sigma'
simulateNodes[['dugongs']] <- 'gamma'
simulateNodes[['jaw']] <- 'Omega'
## simulateNodes[['litters']] <- c('mu', 'theta')
simulateNodes[['seeds']] <- 'b'

relTol_default <- relTol
relTols <- list()
length(relTols) <- length(examples)
names(relTols) <- examples

relTols[['blocker']] <- c(1e-15, 1e-5, 1e-2, 1e-2)
relTols[['dyes']] <- c(1e-15, 1e-7, 1e-3, 1e-3)
relTols[['epil']] <- c(1e-15, 1e-6, 1e-2, 1e-2)
relTols[['equiv']] <- c(1e-15, 1e-6, 1e-4, 1e-4)
relTols[['line']] <- c(1e-15, 1e-7, 1e-5, 1e-4)
relTols[['pump']] <- c(1e-15, 1e-7, 1e-4, 1e-4)
relTols[['rats']] <- c(1e-15, 1e-5, 1e-1, 1e-1)
relTols[['beetles']] <- c(1e-14, 1e-8, 1e-4, 1e-4)
relTols[['oxford']] <- c(1e-15, 1e-4, 1e-6, 1e-4)
relTols[['schools']] <- c(1e-15, 1e-3, 1e-1, 1e-1)

newConstantNodes <- list()
newConstantNodes[['blocker']] <- list(rt = rbinom(22, 10, .5), rc = rbinom(22, 10, .5))
newConstantNodes[['oxford']] <- list(r0 = rbinom(120, 10, .5), r1 = rbinom(120, 10, .5))
newConstantNodes[['pump']] <- list(x = rpois(10, 2))
newConstantNodes[['beetles']] <- list(r = sample(0:1, 8, replace = TRUE))
newConstantNodes[['seeds']] <- list(r = rbinom(21, 4, .5))
newConstantNodes[['jaw']] <- list(R = crossprod(matrix(rnorm(4*4), 4)))
newConstantNodes[['schools']] <- list(pr = crossprod(matrix(rnorm(3*3), 3)))
newConstantNodes[['dyes']]  <- list(y = matrix(rnorm(30, 1500, 10), 6 ,5), tau.within = 1/450, tau.between = 1/420)

newUpdateNodes <- list()
newUpdateNodes[['jaw']] <- list(Omega = crossprod(matrix(rnorm(4*4), 4)))
newUpdateNodes[['schools']] <- list(T = crossprod(matrix(rnorm(3*3), 3)))
newUpdateNodes[['dyes']]  <- list(tau.within = 1/450, tau.between = 1/420)
newUpdateNodes[['oxford']] <- list(beta1 = .0035, beta2 = .000634)
newUpdateNodes[['dugongs']] <- list(gamma = .83)

wrtGeneration <- rep('given', length(examples))
names(wrtGeneration) <- examples
wrtGeneration['jaw'] <- 'prior'
# wrtGeneration['schools'] <- 'prior'

xNew <- list()
xNew[['oxford']] <- list(alpha = .11, beta1 = .12, beta2 = .13, sigma = .9, mu = rnorm(120,0,.1), b = rnorm(120, 0, .1))
xNew[['schools']] <- list(theta = .11, phi = .12, T = matrix(c(15.3, -4.1, -4.1, -4.1, 150.7, -35.23, -4.1, -35.23, 150.7), 3),
                         gamma = rnorm(3, 0, .1), beta = rnorm(8, 0, .1), alpha = matrix(rnorm(38*3, 0, .1), 38, 3))
xNew[['dyes']] <- list(mu = rnorm(6, 1500, 4), theta = rnorm(1, 1500, 5))  # otherwise uncompiled derivs are squaring large y-mu deviations
xNew[['blocker']] <- list(tau = 4.2, mu= rnorm(22, 1, 0.5), delta = rnorm(22, 1, 0.5), delta.new = .5, d = .5)
xNew[['rats']] <- list(alpha = rnorm(30, 250, 10), beta = rnorm(30, 6, 1))  # otherwise, if use much smaller values, get bad uncompiled Hessian with h=1e-4 because magnitude of function is very large O(-1e7)
xNew[['jaw']] <- list(beta0 = 35, beta1 = 0.1, beta2 = -0.1, Omega = matrix(c(1.5,1.2,-1.5,-.1,1.2,8,-3,2.5,-1.5,-3,5,-1.4,-.1,2.5,-1.4,7),4))

inits <- list()
length(inits) <- length(examples)
names(inits) <- examples

## NCT issue 350
inits[['dyes']] <- list(tau.between = 1/400, tau.within = 1/400, theta = 1500)
inits[['blocker']] <- list(tau  = 4, mu= rnorm(22, 1, 0.5), delta = rnorm(22, 1, 0.5), delta.new = 1, d = 1)
inits[['epil']] <- list(tau.b1 = 4, alpha0 = 1, alpha.Base = 0.1, alpha.Trt = -0.1, alpha.BT = 0.1, alpha.Age = -0.1, alpha.V4 = 0.1)

for(i in seq_along(examples)) {
    set.seed(1)
    cat("Testing ", examples[i], ".\n")
    if(is.na(initsFile[i])) tmpInits <- NULL else tmpInits <- initsFile[i]
    if(is.na(dataFile[i])) tmpData <- NULL else tmpData <- dataFile[i]
    if(!is.null(inits[[i]]))
        for(j in seq_along(inits[[i]])) 
            tmpInits[[names(inits[[i]])[j]]] <- inits[[i]][[j]]
    model <- readBUGSmodel(model = bugsFile[i], inits = tmpInits, data = tmpData, useInits = TRUE,
                           dir = nimble:::getBUGSexampleDir(examples[i]))
    if(!is.null(simulateNodes[[i]])) {
        model$simulate(simulateNodes[[i]])
    }
    out <- model$calculate()
    if(is.null(relTols[[examples[i]]])) {
        relTol <- relTol_default
    } else relTol = relTols[[examples[i]]]
    test_ADModelCalculate(model, newConstantNodes = newConstantNodes[[examples[i]]], newUpdateNodes = newUpdateNodes[[examples[i]]],
                          x = wrtGeneration[examples[i]], xNew = xNew[[examples[i]]], useParamTransform = TRUE, relTol = relTol, verbose = verbose, name = bugsFile[i],
                          useFasterRderivs = TRUE, absTolThreshold = 1e-12)
}

## 2022-07 results
## dyes: some Hessian errors ~0.1 for HMC partial and ML partial; probably fine given other uncompiled Hessian issues with these examples but haven't checked
## rats: some Hessian errors ~0.1; checked a case and it seems fine - just numerical errors in uncompiled Hessians
## oxford: checked one of the Hessian anomalies and messing with h for finite element brings uncompiled and compiled into alignment
## schools: various Hessian anomalies; it appears these are simply numerical errors in uncompiled Hessians.
## jaw: very large anomalies; uncompiled Hessian seems unstable but not yet able to verify that the issue is its inaccuracy.

## Running leuk and salm requires monkeying with code, which makes it a pain to deal with directories as done above, so do individually here.

## leuk
set.seed(1)
relTolTmp <- relTol
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-3
relTolTmp[4] <- 1e-3
writeLines(c("var", "Y[N,T],", "dN[N,T];"), con = file.path(tempdir(), "leuk.bug"))
system.in.dir(paste("cat leuk.bug >>", file.path(tempdir(), "leuk.bug")), dir = system.file('classic-bugs','vol1','leuk', package = 'nimble'))
system.in.dir(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
model <- readBUGSmodel(model = file.path(tempdir(), "leuk.bug"), data = system.file('classic-bugs','vol1','leuk','leuk-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','leuk','leuk-init.R', package = 'nimble'), useInits = TRUE)
out <- model$calculate()
newConstantNodes <- list(dN = matrix(rpois(17*42, 2), 17))
test_ADModelCalculate(model, newConstantNodes = newConstantNodes, useParamTransform = TRUE, relTol = relTolTmp, verbose = verbose, name = 'leuk', useFasterRderivs = TRUE, absTolThreshold = 1e-12)

## 2022-06-25: Some fairly extreme 2d11 discrepancies
## Detected some values out of tolerance   :  rOutput2d11$jacobian   cOutput2d11$jacobian .
##            [,1]     [,2]     [,3]
## [1,]  -59.30499  -15.001 2.953402
## Detected some values out of tolerance   :  rOutput12$hessian   cOutput12$hessian .
##             [,1]        [,2]         [,3]
## [1,]  0.01189137  0.01187638 1.261831e-03
## Detected some values out of  absolute  tolerance   :  rOutput2d11$jacobian   cOutput2d11$jacobian .
##           [,1]        [,2]      [,3]
## [1,] 9.9941568 -6.34527345 2.5750553
## Detected some values out of  absolute  tolerance   :  cOutput2d$value   c(cOutput012$hessian) .
## [1] 9.664767e-02 9.664767e-02 2.871831e-15
## Detected some values out of    tolerance   :  rOutput01$jacobian   cOutput01$jacobian .
## [1] -1.892888e+01 -1.892887e+01  3.018687e-07


## salm: easy to have this blow up because of exponentiation unless 'gamma' (part of wrt) is quite small
relTolTmp <- relTol
relTolTmp[2] <- 1e-5
relTolTmp[3] <- 1e-2
relTolTmp[4] <- 1e-2
set.seed(1)
writeLines(c("var","logx[doses];"), con = file.path(tempdir(), "salm.bug"))
system.in.dir(paste("cat salm.bug >>", file.path(tempdir(), "salm.bug")), dir = system.file('classic-bugs','vol1','salm', package = 'nimble'))
model <- readBUGSmodel(model = file.path(tempdir(), "salm.bug"), data = system.file('classic-bugs','vol1','salm','salm-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','salm','salm-init.R', package = 'nimble'), useInits = TRUE)
# model <- readBUGSmodel(model = file.path(tempdir(), "salm.bug"), data = system.file('classic-bugs','vol1','salm','salm-data.R', package = 'nimble'),  inits = list(tau = 0.1, alpha.star = 0.1, beta = 0.1, gamma = 0.01), useInits = TRUE)
model$simulate('lambda')
model$calculate()
newUpdateNodes <- list(gamma = 0.012)
newConstantNodes <- list(y = matrix(rpois(6*3, 2), 6))
xNew <- list(gamma = .012)
test_ADModelCalculate(model, xNew = xNew, newUpdateNodes = newUpdateNodes, newConstantNodes = newConstantNodes, useParamTransform = TRUE, relTol = relTolTmp, verbose = verbose, name = 'salm', useFasterRderivs = TRUE)
## 2022-07-02: Note that I checked the first two of these 4 cases and compiled Hessians are correct compared to hand-calculated derivs
## Detected some values out of  relative  tolerance   :  rOutput12$hessian   cOutput12$hessian .
##           [,1]       [,2]       [,3]
## [1,] 0.06298828 0.06220637 0.01256961
## Detected some values out of  relative  tolerance   :  rOutput2d11$jacobian   cOutput2d11$jacobian .
##            [,1]        [,2]      [,3]
## [1,]  31.603067 -0.05611477 564.18623
## Detected some values out of  relative  tolerance   :  rOutput01$jacobian   cOutput01$jacobian .
##             [,1]        [,2]         [,3]
## [1,]  0.04182428  0.04166845 3.739910e-03
## Detected some values out of  relative  tolerance   :  rOutput12$hessian   cOutput12$hessian .
##           [,1]        [,2]       [,3]
## [1,] -0.046875 -0.05611477 0.16465846



## 2022-03-29: some large magnitude discrepancy between R and C single-taped Hessians, though relative discrepancy not much bigger than 0.001; some cases where R 2d11 hessian values not zero when true value is zero, leading to big discrepancy; some comparisons equal but not identical, couple other minor discrepancies

## March 2022 results:

## epil2: R vs C 2d11 hessian discrepancy
## blocker: R vs C 2d11 hessian discrepancy (some large, very likely NCT issue 350), some compiled comparisons equal but not identical
## dyes: some compiled comparisons equal but not identical; HMC/MAP partial and ML partial has some non-negligible discrepancy of R vs C single-taped hessians - even a few of these after more careful setting of tau.{between,within} values; it looks like the default h for pracma::hessian may be too small in this case (interestingly numDeriv::hessian does better here)
## dyes: R 2d11 Hessian for [3,3] element (tau.between) is wrong for HMC/MAP when tau.between=1; NCT issue 350
## equiv: R vs C 2d11 hessian discrepancy (some large), some compiled comparisons equal but not identical, compiled 2d/012 hessian minor discrepancy, a few R vs C minor discrepancy
## line: R vs C 2d11 hessian discrepancy (some large),some comparisons equal but not identical
## pump: R vs C 2d11 hessian discrepancy; some comparisons equal but not identical, couple other minor discrepancies
## rats: R vs C 2d11 hessian discrepancy (some large); some comparisons equal but not identical, couple other minor discrepancies; some exceedance of thresholds for R vs C
## beetles: R vs C 2d11 hessian discrepancy; some comparisons equal but not identical, couple other minor discrepancies
## jaw: R vs C 2d11 hessian discrepancy; some comparisons equal but not identical, couple other minor discrepancies
## dugongs: R vs C 2d11 hessian discrepancy; some comparisons equal but not identical, couple other minor discrepancies
## seeds: R vs C 2d11 hessian discrepancy; some comparisons equal but not identical, couple other minor discrepancies
## schools: R vs C 2d11 hessian discrepancy; some comparisons equal but not identical
## oxford: R vs C 2d11 hessian discrepancy (some large),some comparisons equal but not identical; some other minor R vs C discrepancies;


## Feb 2021 results:

## epil2: 2d11 uncompiled hessian out of tolerance with compiled

## blocker: uncomp 2d11  hessian way off for one component in various cases - possibly issue 275
## various compiled value, 1d jacobian, 2d11 hessian ENI to compiled logProb, 

## Dec 2020 results:

## epil: some serious problems on at least one test case

## Results as of 2020-12-23
## blocker: compiled 0th deriv equal not identical to cLogProb_orig for HMC/MAP partial
## dyes: all set
## epil: compiled 0th deriv equal not identical to cLogProb_orig, cLogProb{12,01} equal not identical to cLogProb_{orig,new}
## equiv: all set
## line: all set
## oxford: compiled 0th deriv equal not identical to cLogProb_orig, cLogProb{12,01} equal not identical to cLogProb_{orig,new} 
## pump: compiled 0th deriv equal not identical to cLogProb_orig, rWrt and cWrt equal but not identical to x, cLogProb12 and cVals12 equal but not identical to cLogProb_orig and cVals_orig
## rats:
## HMC/MAP partial: compiled 0th deriv equal not identical to cLogProb_orig
## ML partial - compiled 0th deriv equal not identical to cLogProb_orig; uncompiled/compiled Jacobian values out of tolerance
## beetles: compiled 0th deriv equal not identical to cLogProb_orig (and effect on comp/uncomp 0th deriv match); cLogProb{12,01} equal not identical to cLogProb_{orig,new}, 
## jaw: cVals{12,01} equal not identical to cVals_{orig,new}; rWrt equal not identical to x.
## dugongs: issue 255
## schools: cVals{12,01} equal not identical to cVals_{orig,new}; cWrt, rWrt equal not identical to x.
## leuk: tons of NAs in compiled derivs, presumably issue #253
## seeds: compiled 0th deriv equal not identical to cLogProb_orig, cLogProb{12,01} equal not identical to cLogProb_{orig,new}


################### OLD NOTES #####################
## Issues as of spring 2020:
## oxford gives all NaN values for R derivs; C derivs have some -Inf for gradient and some NaNs for Hessian.
## Need to look at model structure to see what is going on.

## dyes has differences O(1e-3) for gradient and O(0.1) for Hessian; need to look in more detail
## dugongs has differences O(1) for gradient and hessian


## Issues as of spring 2019; keeping for now in case they help in understanding current issues listed above.
## dyes: somewhat large relative difference for Hessian - see code below
## for partial deriv wrt tau.between,tau.within, c deriv is 0 and R deriv is -0.0617
## direct numerical calculation via second difference also gives 0
## so something is weird with R deriv
## dugongs: various major mismatches, including logprob

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
