source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)
##nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for calculate() for nimbleModels")

## First set of tests are one's Nick developed

## Second set are those under development by Chris, which use the
## test_ADModelCalculate that allows for random wrt args

## Start of Nick's tests ##

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
  x[1] ~ dnorm(b/(1 - a), sd = sqrt(sigPN^2 + sigOE^2))
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
                                 c('x[1]')), tolerance = .5,
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
                        tolerance = .5,
                        order = c(0, 1, 2))
})

## Ragged arrays not supported for derivs yet.
# dir = nimble:::getBUGSexampleDir('bones')
# Rmodel <- readBUGSmodel('bones', data = NULL, inits = NULL, dir = dir, useInits = TRUE,
#                         check = FALSE)
# test_ADModelCalculate_nick(Rmodel, calcNodeNames = list(Rmodel$getDependencies('theta'), Rmodel$getDependencies('grade')),
#                       wrt = list(c('theta'), c('grade'), c('theta', 'grade'), c('grade[1, 2]')), testR = TRUE,
#                       testCompiled = FALSE)


## end of Nick's tests ##

## Start of Chris' tests ##

## the first few of these mimic and may replace Nick's tests

## basic state space model
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
test_ADModelCalculate(name = 'basic state space', model)

## basic tricky indexing
code <- nimbleCode({
    y[1:2] ~ dmnorm(z[1:2], diagMat[,])
    z[1] <- x[2]
    z[2:3] <- x[1:2] + c(1,1)
    x[1] ~ dnorm(.1, 10)
    x[2] ~ dnorm(1, 3)
})
model <- nimbleModel(code, dimensions = list(x = 2, y = 2, z = 3), constants = list(diagMat = diag(2)),
                     inits = list(x = c(1, 1.2), y  = c(-.1,-.2)))
test_ADModelCalculate(name = 'basic tricky indexing', model)

## state space model
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
test_ADModelCalculate(model, name = 'state space model')

## link functions on stochastic nodes (not covered in BUGS examples)
## plus alt params and NIMBLE-provided distributions
code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dpois(mu)
    log(mu) ~ dnorm(mu0, sd = sigma)
    sigma ~ dinvgamma(shape = a, rate = b)
    a ~ dexp(scale = 1.5)
    b ~ dexp(2.1)
    mu0 ~ dflat()
})
n <- 10
model <- nimbleModel(code, constants = list(n = n), data = list(y = rpois(n, 1)),
                     inits = c(mu0 = rnorm(1), sigma = runif(1), mu = rnorm(1), a = runif(1), b = runif(1)))
test_ADModelCalculate(model, name = 'stochastic link model')

## vectorized deterministic nodes

code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dpois(mu[i])
        logmu[i] ~ dnorm(0, 1)
    }
    mu[1:n] <- exp(logmu[1:n])
})
n <- 10
model <- nimbleModel(code, constants = list(n = n), data = list(y = rpois(n, 1)),
                     inits = c(logmu = rnorm(n)))
test_ADModelCalculate(model, name = 'deterministic vectorized model')

## truncation and constraints
set.seed(1)
code <- nimbleCode({
    for(i in 1:n) {
        y[i, 1:4] ~ dmnorm(mu[1:4], pr[1:k,1:k])
    }
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
    constraint <- dconstraint(mu[1] + mu[2] > 0)
    mu[3] ~ T(dnorm(0,1), -3, 3)
    mu[4] ~ T(dnorm(0,1), , 0)
})
n <- 10
inits <- list(mu = c(0.35, -0.25, 1.3, -2.7), pr = diag(rep(1,4)))
y <- matrix(0, n, 4)
for(i in 1:4)
    y[i, ] <- rmnorm_chol(1, inits$mu, diag(rep(0.2, 4)), prec_param = FALSE)
model <- nimbleModel(code, constants = list(n = n), data = list(y = y, constraint = 1),
                     inits = inits)
test_ADModelCalculate(model, name = 'truncation/constraints model')


## complicated indexing
## fails to compile because of dinvwish
code <- nimbleCode({
    x[2:4,3:5] ~ dinvwish(S[1:3,1:3], df = nu)
    x[1,2:5] ~ dmnorm(z[1:4], pr[1:4,1:4])
    alphas[1:4] <- exp(x[1, 1:4])
    x[2:5, 2] ~ ddirch(alphas[1:4])
    x[1,1] ~ dpois(3)
    for(i in 2:5)
        x[i, 1] ~ dgamma(1, 1)
    for(i in 3:5)
        x[5, i] ~ dnorm(0, 1)
    w[1:3] ~ dmnorm(z[1:3], cov = x[2:4, 3:5])
    for(j in 1:5)
        y[1:5, j] ~ dmnorm(x[j, 1:5], cov = pr[1:5,1:5])
})
inits <- list(nu = 7, S = diag(rep(1, 3)), z = rep(0, 4), pr = diag(rep(1,5)))
model <- nimbleModel(code, inits = inits)
model$simulate()
model$setData('y','w')
test_ADModelCalculate(model, name = 'complicated indexing')

## complicated indexing w/o Wishart
## Dirichlet causes diffs as Rderivs uses calculate and would get NaN
## but Cderivs doesn't enforce sum to 1 in actual density calc.
code <- nimbleCode({
    x[2:4,3:5] <- S[1:3,1:3]
    x[1,2:5] ~ dmnorm(z[1:4], pr[1:4,1:4])
    alphas[1:4] <- exp(x[1, 1:4])
    x[2:5, 2] ~ ddirch(alphas[1:4])
    x[1,1] ~ dpois(3)
    for(i in 2:5)
        x[i, 1] ~ dgamma(1, 1)
    for(i in 3:5)
        x[5, i] ~ dnorm(0, 1)
    w[1:3] ~ dmnorm(z[1:3], cov = x[2:4, 3:5])
    for(j in 1:5)
        y[1:5, j] ~ dmnorm(x[j, 1:5], cov = pr[1:5,1:5])
})
inits <- list(nu = 7, S = diag(rep(1, 3)), z = rep(0, 4), pr = diag(rep(1,5)))
model <- nimbleModel(code, inits = inits)
model$simulate()
model$setData('y','w')
test_ADModelCalculate(model, name = 'complicated indexing')

## complicated indexing w/o Wishart or Dirichlet
code <- nimbleCode({
    x[2:4,3:5] <- S[1:3,1:3]
    x[1,2:5] ~ dmnorm(z[1:4], pr[1:4,1:4])
    alphas[1:4] <- exp(x[1, 1:4]) / sum(exp(x[1,1:4]))
    x[2:5, 2] <- alphas[1:4]
    x[1,1] ~ dpois(3)
    for(i in 2:5)
        x[i, 1] ~ dgamma(1, 1)
    for(i in 3:5)
        x[5, i] ~ dnorm(0, 1)
    w[1:3] ~ dmnorm(z[1:3], cov = x[2:4, 3:5])
    for(j in 1:5)
        y[1:5, j] ~ dmnorm(x[j, 1:5], cov = pr[1:5,1:5])
})
inits <- list(nu = 7, S = diag(rep(1, 3)), z = rep(0, 4), pr = diag(rep(1,5)))
model <- nimbleModel(code, inits = inits)
model$simulate()
model$setData('y','w')
test_ADModelCalculate(model, name = 'complicated indexing')

## MVN with various parameterizations and user-defined functions

covFunLoop <- nimbleFunction(
    run = function(dist = double(2), rho = double(0)) {
        n = dim(dist)[1]
        out = nimMatrix(nrow = n, ncol = n)
        for(i in 1:n)
            for(j in 1:n)
                out[i,j] <- exp(-dist[i,j]/rho)
        returnType(double(2))
        return(out)
    }, enableDerivs = TRUE)

covFunVec <- nimbleFunction(
    run = function(dist = double(2), rho = double(0)) {
        out <- exp(-dist/rho)
        returnType(double(2))
        return(out)
    }, enableDerivs = TRUE)

GPdist <- nimbleFunction(
    run = function(x = double(1), dist = double(2), rho = double(0)) {
        returnType(double(0))
        Sigma <- exp(-dist/rho)
        U <- t(chol(Sigma))
        out <- -sum(log(diag(U)))
        qf <- forwardsolve(U, x)
        out <- out - 0.5*sum(qf^2)
        return(out)
    })


code <- nimbleCode({
    Sigma1[1:n,1:n] <- exp(-dist[1:n,1:n]/rho)
    Sigma2[1:n,1:n] <- covFunLoop(dist[1:n,1:n], rho)
    Sigma3[1:n,1:n] <- covFunVec(dist[1:n,1:n], rho)
    
    y[1, 1:n] ~ dmnorm(mu1[1:n], cov = Sigma1[1:n,1:n])
    y[2, 1:n] ~ dmnorm(mu2[1:n], cov = Sigma2[1:n,1:n])
    y[3, 1:n] ~ dmnorm(mu3[1:n], cov = Sigma3[1:n,1:n])

    Q[1:n,1:n] <- inverse(Sigma1[1:n, 1:n])
    y[4, 1:n] ~ dmnorm(mu4[1:n], Q[1:n,1:n])

    Uprec[1:n, 1:n] <- chol(Q)
    Ucov[1:n, 1:n] <- chol(Sigma1)
    y[5, 1:n] ~ dmnorm(mu5[1:n], cholesky = Uprec, prec_param = 1)
    y[6, 1:n] ~ dmnorm(mu6[1:n], cholesky = Ucov, prec_param = 0)
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
    nu ~ dunif(10, 20)
})


n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
R <- crossprod(matrix(rnorm(n^2), n, n))
model <- nimbleModel(code, constants = list(n = n), inits = list(dist = dd, R = R, nu = 15, rho = rgamma(1),
                                                                 z = rep(0, n)))
model$simulate()
model$setData('y')
test_ADModelCalculate(model, name = 'various multivariate dists')
                    

## loop through BUGS models

## example test of a BUGS model
model <- readBUGSmodel(model = 'epil2.bug', inits = 'epil-inits.R', data = 'epil-data.R',
                       useInits = TRUE, dir = nimble:::getBUGSexampleDir('epil'))
set.seed(1)
model$simulate('b1')
model$calculate()
test_ADModelCalculate(model, name = 'epil', order = c(0,1))

## now loop through BUGS models.
## figure out if need to simulate for any of these as we do above for 'epil'
## also will need to specify bugs,inits,data files by specific name in some cases where naming is non-standard
examples <- c('epil', 'alli', 'biops', 'blocker', 'bones', 'dyes', 'equiv', 'lsat', 'line', 'oxford',
              'pump', 'rats', 'schools', 'dugongs', 'jaw', 'beetles')
bugsFile <- examples
initsFile <- dataFile <- rep(NA, length(examples))
names(bugsFile) <- names(initsFile) <- names(dataFile) <- examples

## customize file names as needed
bugsFile['beetles'] <- 'beetles-logit'
bugsFile['jaw'] <- 'jaw-constant'
bugsFile['epil'] <- 'epil2'
bugsFile <- paste0(bugsFile, '.bug')
initsFile['epil'] <- 'epil-inits.R'
dataFile['epil'] <- 'epil-data.R'

orders <- rep(2, length(examples))

## customize whether any nodes need to be initialized
simulateNodes <- list()
length(simulateNodes) <- length(examples)
names(simulateNodes) <- examples
simulateNodes['epil'] <- 'b1'

for(i in seq_along(examples)) {
    if(is.na(initsFile[i])) tmpInits <- NULL else tmpInits <- initsFile[i]
    if(is.na(dataFile[i])) tmpData <- NULL else tmpData <- dataFile[i]
    model <- readBUGSmodel(model = bugsFile[i], inits = tmpInits, data = tmpData, useInits = TRUE,
                           dir = nimble:::getBUGSexampleDir(examples[i]))
    if(!is.null(simulateNodes[[i]])) {
        model$simulate(simulateNodes[[i]])
        model$calculate()
    }
    ## may need to just do jacobian for bigger models
    test_ADModelCalculate(model, name = bugsFile[i], order = c(0,1))
}


