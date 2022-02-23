source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)

warning("All atomics currently off as result of NCT issue 274")

##nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for calculate() for nimbleModels")


test_that('pow and pow_int work', {

  ## Following example is reduced from BUGS equiv example.
  ## It is a case from which we have had some problems and crashes in the past.
  mc <- nimbleCode({
    d ~ dnorm(0, sd = 10)
    trt <- 0
    a <- -1
    m <- pow(a, trt) + d
    y ~ dnorm(m, sd = 2)
  })
  m <- nimbleModel(mc, data = list(y = 1.5),
                   inits = list(d = .3))

  calcNodes <- c(m$getDependencies("d"))
  wrtNodes <- 'd'
  order <- 0:2
  m$calculate()
  wrapperDerivs <- nimDerivs(m$calculate(calcNodes), wrt = wrtNodes, order = order)
  testFunctionInstance <- testCompiledModelDerivsNimFxn(m, calcNodes, wrtNodes, order)
  cm <- compileNimble(m)
  ctestFunctionInstance <- compileNimble(testFunctionInstance, project =  m, resetFunctions = TRUE)
  cDerivs <- ctestFunctionInstance$run()
  expect_equal(wrapperDerivs$value, cDerivs$value)
  expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian)
  expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = 1e-4)
  m$trt <- cm$trt <- -2
  m$calculate()
  cm$calculate()
  wrapperDerivs <- nimDerivs(m$calculate(calcNodes), wrt = wrtNodes, order = order)
  cDerivs <- ctestFunctionInstance$run()
  expect_equal(wrapperDerivs$value, cDerivs$value)
  expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian)
  expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = 1e-4)

  mc <- nimbleCode({
    d ~ dnorm(0, sd = 10)
    trt <- -1
    a <- -1
    m <- pow_int(a, trt) + d
    y ~ dnorm(m, sd = 2)
  })
  m <- nimbleModel(mc, data = list(y = 1.5),
                   inits = list(d = .3))
  calcNodes <- c(m$getDependencies("d"))
  wrtNodes <- 'd'
  order <- 0:2
  m$calculate()
  wrapperDerivs <- nimDerivs(m$calculate(calcNodes), wrt = wrtNodes, order = order)
  testFunctionInstance <- testCompiledModelDerivsNimFxn(m, calcNodes, wrtNodes, order)
  cm <- compileNimble(m)
  ctestFunctionInstance <- compileNimble(testFunctionInstance, project =  m, resetFunctions = TRUE)
  cDerivs <- ctestFunctionInstance$run()
  expect_equal(wrapperDerivs$value, cDerivs$value)
  expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian)
  expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = 1e-4)
  m$trt <- cm$trt <- -2
  m$calculate()
  cm$calculate()
  wrapperDerivs <- nimDerivs(m$calculate(calcNodes), wrt = wrtNodes, order = order)
  cDerivs <- ctestFunctionInstance$run()
  expect_equal(wrapperDerivs$value, cDerivs$value)
  expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian)
  expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = 1e-4)  
})

## check logic of results
test_that('makeUpdateNodes works correctly', {
    ## updateNodes should include all nodes that affect calculate(calcNodes) that are not in wrt,
    ## including any immediate (either stoch or det) parent nodes not in either calcNodes or in wrt.
    ## Will include stoch nodes that are in calcNodes but not wrt, but not det nodes that are in calcNodes.
    ## per PdV, RHS-only nodes should be (but are not yet) in constantNodes.
    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dnorm(mu[i], var = sigma2)
            mu[i] ~ dnorm(mu0, tau)
        }
        sigma2 ~ dgamma(1.3, 0.3)
        tau ~ dgamma(0.7, 0.5)
    })
    n <- 3
    model <- nimbleModel(code, data = list(y = rnorm(n)), constants = list(n = n),
                         inits = list(mu = rnorm(n), sigma2 = 2.1, tau = 1.7, mu0 = 0.8))

    ## Various wrt and calcNodes cases; not exhaustive, but hopefully capturing most possibilities
    
    ## distinct wrt and calcNodes: updateNodes should have 'mu' and determ parents not in calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu','y'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP",
                                           model$expandNodeNames('mu')))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## distinct wrt and calcNodes, where calcNodes includes a deterministic node
    ## updateNodes should have 'mu' and determ parents not in calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu','y', "lifted_d1_over_sqrt_oPtau_cP"), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP",
                                           model$expandNodeNames('mu')))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## with option to include data nodes in updateNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu','y'), model,
                              dataAsConstantNodes = FALSE)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP",
                                       model$expandNodeNames('mu'), model$expandNodeNames('y')))
    expect_identical(result$constantNodes, character(0))

    ## scalar calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu[1]','y[1]'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP", "mu[1]"))
    expect_identical(result$constantNodes, "y[1]")

    ## subset of a variable as calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu[1:2]','y[1:2]'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP", "mu[1]", "mu[2]"))
    expect_identical(result$constantNodes, c("y[1]", "y[2]"))

    ## a wrt node included in calcNodes
    ## updateNodes should have the stoch calcNodes not in wrt, plus determ parents not in calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0'), calcNodes = c('sigma2', 'mu0', 'mu','y'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP",
                                           model$expandNodeNames('mu')))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## some overlap of wrt and calcNodes
    ## updateNodes should have determ parents not in calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'tau', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu','y'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## full overlap of wrt and calcNodes
    ## updateNodes should have determ parents not in calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu','y'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## wrt are intermediate nodes
    ## updateNodes should have determ and stoch parents not in calcNodes
    result <- makeUpdateNodes(wrt = c('mu'), calcNodes = c('mu','y'), model)
    expect_identical(result$updateNodes, c("mu0", "lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## no data nodes in calcNodes
    ## updateNodes should have determ parents not in calcNodes
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu'), model)
    expect_identical(result$updateNodes, c("lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$constantNodes, character(0))

    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dnorm(mu[i], var = sigma2)
        }
        mu[1:n] ~ dmnorm(mu0[1:n], pr[1:n, 1:n])
        sigma2 ~ dgamma(1.3, 0.3)
    })
    n <- 3
    model <- nimbleModel(code, data = list(y = rnorm(n)), constants = list(n = n),
                         inits = list(mu = rnorm(n), sigma2 = 2.1, mu0 = rnorm(n), pr = diag(n)))

    prElems <- model$expandNodeNames('pr', returnScalarComponents = TRUE)
    lftChElems <- model$expandNodeNames('lifted_chol_oPpr_oB1to3_comma_1to3_cB_cP', returnScalarComponents = TRUE)

    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0'), calcNodes = c('mu','y'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems,
                                       model$expandNodeNames('mu', returnScalarComponents = TRUE)))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0'), calcNodes = c('mu','y'), model,
                              dataAsConstantNodes = FALSE)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems,
                                           model$expandNodeNames('mu', returnScalarComponents = TRUE),
                                           model$expandNodeNames('y')))
    expect_identical(result$constantNodes, character(0))

    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0'), calcNodes = c('mu[1]','y[1]'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems, 'mu[1]'))
    expect_identical(result$constantNodes, "y[1]")

    ## NCT issue 257: should include mu[3] presumably
    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0'), calcNodes = c('mu[1:2]','y[1:2]'), model)
    expect_identical(result$updateNodes, c(prElems, "lifted_sqrt_oPsigma2_cP", lftChElems, model$expandNodeNames('mu')))
    expect_identical(result$constantNodes, c("y[1]", "y[2]"))

    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0'), calcNodes = c('sigma2', 'mu0', 'mu','y'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems,
                                       model$expandNodeNames('mu', returnScalarComponents = TRUE)))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu','y'), model)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeUpdateNodes(wrt = c('mu'), calcNodes = c('mu','y'), model)
    expect_identical(result$updateNodes, c(model$expandNodeNames('mu0', returnScalarComponents = TRUE),
                                           "lifted_sqrt_oPsigma2_cP", lftChElems))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeUpdateNodes(wrt = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu'), model)
    expect_identical(result$updateNodes, c(lftChElems))
    expect_identical(result$constantNodes, character(0))
   
}

## First set of tests are ones Nick developed

## Second set are those under development by Chris, which use the
## test_ADModelCalculate that tests various standard use cases.

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

# This causes crashes due to atomic memory mgmt.
# Removing the first calcNodeNames case fixes it, so something is being over-written there.
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

## Fails on dcat.  Potential issue of support for ragged arrays.
## dir = nimble:::getBUGSexampleDir('bones')
## Rmodel <- readBUGSmodel('bones', data = NULL, inits = NULL, dir = dir, useInits = TRUE,
##                         check = FALSE)
## test_ADModelCalculate_nick(Rmodel, name = 'bones', calcNodeNames = list(Rmodel$getDependencies('theta'), Rmodel$getDependencies('grade')),
##                            wrt = list(c('theta'), c('grade'), c('theta', 'grade'), c('grade[1, 2]')),
##                            order = c(0, 1, 2))


## end of Nick's tests ##

## Start of Chris' tests ##

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
test_ADModelCalculate(model, verbose = TRUE, name = 'basic model, lifted nodes')

## the first few of these mimic and may replace Nick's tests

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
test_ADModelCalculate(model, verbose = TRUE, name = 'basic state space') 
## with SOME random seeds, R and C 2d11 jacobian only match to only 1-2 digits with new updateNode values
## presumably just stochasticity in that the Hessian tolerance is .001 

## basic tricky indexing
code <- nimbleCode({
    y[1:2] ~ dmnorm(z[1:2], cov = covMat[,])
    z[1] <- x[2]
    z[2:3] <- x[1:2] + 1
    x[1] ~ dnorm(.1, 10)
    x[2] ~ dnorm(1, 3)
})
model <- nimbleModel(code, dimensions = list(x = 2, y = 2, z = 3), constants = list(covMat = matrix(c(1.9, .7, .7, 1.3), 2)), inits = list(x = c(1, 1.2), y  = c(-.1,-.2)))
test_ADModelCalculate(model, name = 'basic tricky indexing')

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

test_ADModelCalculate(model, name = 'state space model', useFasterRderivs = TRUE, verbose =  TRUE)
## cOutput1d$value not identical to c(cOutput012$jacobian)
## R and C 2d11 jacobian match to limited digits

## link functions on stochastic nodes (not covered in BUGS examples)
## plus alt params and NIMBLE-provided distributions
code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dpois(mu)
    log(mu) ~ dnorm(mu0, sd = sigma)
    sigma ~ dinvgamma(shape = a, rate = b)
    a ~ dexp(scale = 1.5)
    b ~ dexp(2.1)
    mu0 ~ dnorm(0, .00001)  # dflat()  # dflat not handled 
})
n <- 10
log_mu_init <- rnorm(1)
## Need to initialize lifted node.
model <- nimbleModel(code, constants = list(n = n), data = list(y = rpois(n, 1)),
                     inits = list(mu0 = rnorm(1), sigma = runif(1), mu = exp(log_mu_init),
                                log_mu = log_mu_init, a = runif(1), b = runif(1)))
test_ADModelCalculate(model, name = 'stochastic link model')

## dexp and dt, which are provided by NIMBLE to allow expanded parameterizations
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
model <- nimbleModel(code, constants = list(n = n), data = list(y = rnorm(n)),
                     inits = list(mu = rnorm(n), sigma = runif(1), nu = 2.5, mu0 = rnorm(1), sd0 = runif(1),
                                  tau = runif(1)))
## (not seeing anymore) Heisenbug: with verbose=F (the default), can get cLogProb12 equal but not identical to cLogProb_orig
test_ADModelCalculate(model, name = 'dt and dexp model')
## very slow to run rOutput2d11 (2-3 minutes) (numDeriv::jacobian slower than pracma::jacobian); might want to use fasterRderivs

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
                     inits = list(logmu = rnorm(n)))
test_ADModelCalculate(model, name = 'deterministic vectorized model')
## 2d11 hessian R and C mismatch - as bad as only 1 digit in common

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
test_ADModelCalculate(model, name = 'truncation model')

code <- nimbleCode({
    y ~ dnorm(mu, 1)
    mu ~ T(dnorm(mu0, 1), 0, 10)
    mu0 ~ dnorm(0,1)
})
model <- nimbleModel(code, data = list(y = 1), inits = list(mu = 0.5, mu0 = 1))
## compilation error: issue #254
test_ADModelCalculate(model, name = 'truncation on non-top node')


code <- nimbleCode({
    y ~ dnorm(mu, 1)
    mu ~ T(dbeta(a,b), 0.1, 0.8)
    a ~ dunif(0, 5)
    b ~ dunif(0,5 )
})
model <- nimbleModel(code, data = list(y = 1), inits = list(mu = 0.5,a=1,b=1))


## compilation error
test_ADModelCalculate(model, name = 'truncation with dbeta')


## complicated indexing 
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
test_ADModelCalculate(model, name = 'complicated indexing')
## compiled/uncompiled 01, 012, 02 values out of tolerance in various cases
## comp/unc 2d11 hessian out of tolerance
## same types of issues but not exact same numerical results when using atomics
## occasional rOutput2d11 values way off (e.g., first value in EB, idx =3 is -243 (or -251) (see issue 275 to reproduce) (not seeing as of 2021-04-29)

## if use pracma::jacobian:
## unc/compiled 01,012,02 values out of tolerance in various cases
## unc/compiled 012, 12 jacobian out of tolerance in one case

## using different subsets of a matrix
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
test_ADModelCalculate(model, name = 'different subsets of a matrix')
## comp/unc values out of tolerance for HMC/MAP partial

## vectorized covariance matrix
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
                     inits = list(dist = dd, rho = rgamma(1, 1, 1), z = rep(0, n), pr = diag(n)))
model$simulate()
model$calculate()
model$setData('y')
test_ADModelCalculate(model, excludeUpdateNodes = 'dist', name = 'dmnorm with vectorized covariance matrix')
## comp 2d/012 hessian out of tolerance, comp/unc values out of tolerance for ML partial for non-atomics
## with atomics, only one comp 2d/012 hessian out of tolerance

## vectorized covariance matrix, chol param
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
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(dist = dd, rho = rgamma(1, 1, 1),
                                  z = rep(0, n), pr = diag(n)))
model$simulate()
model$calculate()
model$setData('y')
## rOutput 01, 012, 02 values out of tolerance with compiled counterparts in one use case
test_ADModelCalculate(model, excludeUpdateNodes = 'dist',
                      name = 'dmnorm with vectorized covariance matrix, chol param')


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
    }, enableDerivs = list(run = list(noDeriv_vars = c('i','j'))))
    # enableDerivs = TRUE) # if use NCT 130 work-around

code <- nimbleCode({
    Sigma2[1:n,1:n] <- covFunLoop(dist[1:n,1:n], rho)
    y[1:n] ~ dmnorm(mu2[1:n], cov = Sigma2[1:n,1:n])
    mu2[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    rho ~ dgamma(2, 3)
})

                               
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(rho = rgamma(1, 1, 1), dist = dd, z = rep(0, n), pr = diag(n)))
model$simulate()
model$calculate()
model$setData('y')
test_ADModelCalculate(model, excludeUpdateNodes = 'dist',
                      name = 'dnorm with user-defined fxn for covariance with loops')
## 2d/012 hessian out of tolerance

## user-defined cov function vectorized

covFunVec <- nimbleFunction(
    run = function(dist = double(2), rho = double(0)) {
        out <- exp(-dist/rho)
        returnType(double(2))
        return(out)
    }, enableDerivs = TRUE)

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
model$simulate()
model$calculate()
model$setData('y')
## compiled 2d, 012 hessians out of tolerance in several cases
## compiled and uncompiled 01, 012, 02 value out of tolerance in a couple cases
test_ADModelCalculate(model, excludeUpdateNodes = 'dist',
                      name = 'dmnorm with user-defined vectorized fxn')
## 2d/012 hessian out of tolerance

## other dmnorm parameterizations
code <- nimbleCode({
    y[1, 1:n] ~ dmnorm(mu1[1:n], Q[1:n,1:n])

    Uprec[1:n, 1:n] <- chol(Q[1:n,1:n])
    Ucov[1:n, 1:n] <- chol(Sigma[1:n,1:n])
    y[2, 1:n] ~ dmnorm(mu2[1:n], cholesky = Uprec[1:n,1:n], prec_param = 1)
    y[3, 1:n] ~ dmnorm(mu3[1:n], cholesky = Ucov[1:n,1:n], prec_param = 0)
    mu1[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    mu2[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    mu3[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
})


n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
Sigma <- exp(-dd/0.1)
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(z = rep(0, n), pr = diag(n),
                                  Sigma = Sigma, Q = solve(Sigma)))
model$simulate()
model$calculate()
model$setData('y')
## comp/unc value out of tolerance for EB
test_ADModelCalculate(model, name = 'various dmnorm parameterizations')

## various dmvt parameterizations

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
                     inits = list(z = rep(0, n), pr = diag(n),
                                  Sigma = Sigma, Q = solve(Sigma)))
model$simulate()
model$calculate()
model$setData('y')
## uncompiled/compiled 2d11 hessian out of tolerance in some cases 
## unc/com 01,012,02 value out of tolerance in some cases
test_ADModelCalculate(model, name = 'various dmvt parameterizations')


## dirichlet as likelihood so not differentiating wrt something with constraint.
code <- nimbleCode({
    p[1:k] ~ ddirch(alpha[1:k])
    for(i in 1:k)
        alpha[i] ~ dgamma(1.3, 1.5)
})
k <- 4
model <- nimbleModel(code, constants = list(k = k), data = list(p = c(.2, .4, .15, .25)), inits = list(alpha = runif(4)))
test_ADModelCalculate(model, name = 'Dirichlet likelihood')

## dwish and dinvwish so long as not differentiating w.r.t. the random variable (since it has constraints)
## Note that nu must exceed n and can't be set to runif(0,1) via test_ADModelCalculate.
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

## comp/unc values out of tolerance in some cases
## 2d, 012 comp hessian out of tolerance in some cases
test_ADModelCalculate(model, excludeUpdateNodes = 'dist', verbose = TRUE, name = 'dwish, dinvwish')


## simple user-defined distribution
dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0, default = 1), 
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- log(rate) - x*rate
        if(log) return(logProb)
        else return(exp(logProb)) 
    }, enableDerivs = TRUE)

code <- nimbleCode({ 
    y ~ dmyexp(rho)
    rho ~ dgamma(2, 3)
})

model <- nimbleModel(code, data = list(y = rgamma(1,1,1)), inits = list(rho = rgamma(1, 1, 1)))
test_ADModelCalculate(model, name = 'simple user-defined distribution')


## vectorized powers cause problems, issue #253
dtest <- nimbleFunction(
    run = function(x = double(0), mu = double(1), 
                   log = integer(0, default = 0)) {
        returnType(double(0))
        tmp <- mu^2
        out <- tmp[1]
        if(log) return(out) else return(exp(out))
    }, enableDerivs = TRUE)

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

model <- nimbleModel(code, inits = list(mu = rnorm(2)))
model$simulate()
model$calculate()
model$setData('y')
## NaN values in compiled derivs
test_ADModelCalculate(model, name = 'vectorized power')


## user-defined distribution, issue #253

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
    }, enableDerivs = TRUE)

rGPdist <- nimbleFunction(
    run = function(n = integer(0), dist = double(2), rho = double(0)) {
        returnType(double(1))
        Sigma <- exp(-dist/rho)
        U <- chol(Sigma)
        p <- dim(dist)[1]
        out <- (U %*% rnorm(p))[,1]
        return(out)
    })

code <- nimbleCode({ 
    y[1:n] ~ dGPdist(dist[1:n, 1:n], rho)
    rho ~ dgamma(2, 3)
})


n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n, dist = dd), inits = list(rho = rgamma(1, 1, 1)))
model$simulate()
model$calculate()
model$setData('y')
## all compiled Jacobian/Hessian values are NaN
## issue seems to be in qf^2
test_ADModelCalculate(model, name = 'user-defined distribution')
## (formerly: model compiles but numerical differences; see later comments in NCT issue #220)

## Use additional matrix functions
## logdet() not yet allowed

code <- nimbleCode({
    y1 ~ dnorm(inprod(beta[1,1:n], x[1:n]), 1)
    # det <- logdet(Sigma[1:n,1:n])
    # y2 ~ dnorm(det, sd = 100)
    Sigma[1:n,1:n] <- sigma2 * exp(-dist[1:n,1:n] / rho)
    w[1,1:n] <- solve(Sigma[1:n,1:n], beta[2,1:n])
    w[2,1:n] <- Sigma[1:n,1:n] %*% beta[3,1:n]
    U[1:n,1:n] <- chol(Sigma[1:n,1:n])
    w[3,1:n] <- backsolve(U[1:n,1:n], beta[4,1:n])
    for(i in 1:4)
        beta[i,1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    for(i in 1:3)
        yy[i,1:n] ~ dmnorm(w[i,1:n], pr[1:n,1:n])
    rho ~ dgamma(2, 3)
    sigma2 ~ dgamma(2, 3)
})
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n, x = rnorm(n), z = rep(0, n), pr = diag(n), dist = dd),
                     inits = list(rho = rgamma(1, 1, 1), sigma2 = rgamma(1, 1, 1)),
                     data = list(y1 = rnorm(1), #  y2 = rnorm(1),
                                 yy = matrix(rnorm(n*3), 3, n)))
model$simulate()
model$calculate()
model$setData(c('y1','yy')) # 'y2'
## unc/comp value exceeds tolerance
## compiled 1d, 012 jacobian equal not identical (only for non-atomics)
## compiled 2d, 012 hessian exceeds tolerance
## compiled 2d11, 012 hessian equal not identical (only for non-atomics)
test_ADModelCalculate(model, excludeUpdateNodes = 'dist', name = 'various matrix functions')

## Various combinations of updateNodes, wrt, calcNodes

code <- nimbleCode({
    a ~ dgamma(1.1, 0.8)
    b ~ dnorm(z, var = a)
    z ~ dnorm(0, 1)
})
model <- nimbleModel(code, data = list(b = 1.2), inits = list(a = 1.3, z = 0.7))
## calcNodes excludes det intermediates
test_ADModelCalculate(model, wrt = 'a', calcNodes = c('a', 'b'), name = 'update nodes case 1a')
model <- nimbleModel(code, data = list(b = 1.2), inits = list(a = 1.3, z = 0.7))
## calcNodes includes det intermediates
test_ADModelCalculate(model, wrt = 'a', calcNodes = c('a', 'b', 'lifted_sqrt_oPa_cP'), name = 'update nodes case 1b')

code <- nimbleCode({
    a ~ dgamma(1.1, 0.8)
    for(i in 1:4) 
        b[i] ~ dnorm(z[i], var = a)
    z[1] ~ dnorm(0, 1)
    z[2] ~ dnorm(0, 1)
    z[3:4] ~ dmnorm(mu0[1:2], pr[1:2,1:2])        
})
model <- nimbleModel(code, data = list(b = rnorm(4)), inits = list(a = 1.3, z = runif(4), pr = diag(2), mu0 = rep(0, 2)))
## calcNodes excludes det intermediates
test_ADModelCalculate(model, wrt = 'a', calcNodes = c('a', 'b'), name = 'update nodes case 2a')
model <- nimbleModel(code, data = list(b = rnorm(4)), inits = list(a = 1.3, z = runif(4), pr = diag(2), mu0 = rep(0, 2)))
## calcNodes includes det intermediates
test_ADModelCalculate(model, wrt = 'a', calcNodes = c('a', 'b', "lifted_sqrt_oPa_cP"), name = 'update nodes case 2b')


## Parameter transform system and full use of ddirch, dwish, dinvwish

## basic variance component
set.seed(1)
code <- nimbleCode({
    y ~ dnorm(mu, sd = sigma)
    sigma ~ dinvgamma(1.3, 0.7)
    mu ~ dnorm(0, 1)
})
model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma = rgamma(1, 1, 1), mu = rnorm(1)))
## compiled values and logProb equal not identical (ENI) to cLogProb_new
## compiled 2d11, 012 hessian ENI
test_ADModelCalculate(model, useParamTransform = TRUE, name = 'basic param transform')


set.seed(1)
code <- nimbleCode({
    y ~ dnorm(mu, var = sigma2)
    sigma2 ~ dinvgamma(1.3, 0.7)
    mu ~ dnorm(0, 1)
})
model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = rgamma(1, 1, 1), mu = rnorm(1)))
## compiled values and logProb equal not identical (ENI) to cLogProb_new
## compiled 2d11, 012 hessian ENI
test_ADModelCalculate(model, useParamTransform = TRUE, name = 'basic param transform, with lifted')

## now check if model is out-of-state
set.seed(1)
code <- nimbleCode({
    y ~ dnorm(0, sd = sigma)
    sigma <- sqrt(sigma2)
    sigma2 ~ dinvgamma(1.3, 0.7)
})

model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = 2))
model$sigma <- 1
model$calculate('y')
model$logProb_y
## no issues
test_ADModelCalculate(model)

model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = 2))
model$sigma <- 1
model$calculate('y')
model$logProb_y
## compiled values and logProb equal not identical (ENI) to cLogProb_new
## compiled 2d11, 012 hessian ENI
test_ADModelCalculate(model, useParamTransform = TRUE)

model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = 2))
model$sigma <- 1
model$calculate('y')
model$logProb_y
test_ADModelCalculate(model, useFasterRderivs = TRUE)

model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = 2))
model$sigma <- 1
model$calculate('y')
model$logProb_y
## compiled values and logProb equal not identical (ENI) to cLogProb_new
## compiled 2d11, 012 hessian ENI
test_ADModelCalculate(model, useFasterRderivs = TRUE, useParamTransform = TRUE)

## Dirichlet
code <- nimbleCode({
    y[1:k] ~ dmulti(p[1:k], n)
    p[1:k] ~ ddirch(alpha[1:k])
    for(i in 1:k)
        alpha[i] ~ dgamma(1.3, 1.5)
})
n <- 30
k <- 4
model <- nimbleModel(code, constants = list(k = k, n = n), data = list(y = rmulti(1, n, rep(1/k, k))), inits = list(p = c(.2, .4, .15, .25), alpha = runif(4)))

## various vals, logProbs, wrt ENI
## comp 2d11/012 hessian ENI
## comp 1d/012 jacobian ENI
## comp/uncom 2d11 jacobian out of tolerance (quite a bit)

## no longer seeing this issue as of 2021-04-29
## ISSUE: ML partial: compiled wrt and x quite different
## cWrt01,012 having values be restored in ML partial wrt=alpha[1,2]
## see issue 277 - it occurs regardless of using full alpha or partial alpha
## per PdV this seems like a corner case we may not fix
test_ADModelCalculate(model, x = 'prior', useParamTransform = TRUE, excludeUpdateNodes = 'p',
                      name = 'Dirichlet paramTransform', seed = 2)
## seed = 1 (default) produces p[1]=2e-6 in new x and that causes -Inf in 0th order deriv and NaN in jac/hessian



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
    ## previous bug prevents this:
    ## y[7, 1:n] ~ dGPdist(dist[1:n, 1:n], rho)
    
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
    nu0 ~ dunif(0, 10)
    nu <- 10 + nu0  ## hack to ensure big enough to ensure p.d. when simulating new updateNodes
})

n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
R <- crossprod(matrix(rnorm(n^2), n, n))
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(dist = dd, R = R, nu0 = 5, rho = rgamma(1, 1, 1),
                                                                 z = rep(0, n)))
model$simulate()
model$calculate()
model$setData('y')

## non-atomics has seg fault in HMC/MAP partial
## *** caught segfault ***
##address 0x55da1bd80ba0, cause 'memory not mapped'
## Traceback:
##  1: .Call("CALL_nfRefClass_R_GlobalEnv46_run", x, order, .basePtr)

## atomics: 
## compvalues, logProb_new, wrt ENI
## comp/unc values out of tolerance
## comp jac/Hessians ENI
## comp/unc 2d11 jacobian quite out of tolerance
## comp/unc 2d/012 hessian out of tolerance 
test_ADModelCalculate(model, excludeUpdateNodes = 'dist', x = 'prior',
                      useParamTransform = TRUE, useFasterRderivs = TRUE,
                      verbose = TRUE,
                      name = 'various multivariate dists')


## loop through BUGS models

## example test of a BUGS model - take this out as epil is in loop below
if(FALSE) {
    model <- readBUGSmodel(model = 'epil2.bug', inits = 'epil-inits.R', data = 'epil-data.R',
                           useInits = TRUE, dir = nimble:::getBUGSexampleDir('epil'))
    set.seed(1)
    model$simulate('b1')
    model$calculate()
    test_ADModelCalculate(model, name = 'epil', useFasterRderivs = TRUE)
}


## now loop through BUGS models.
examples <- c('blocker', 'dyes', 'epil', 'equiv', 'line', 'oxford', 'pump', 'rats', 'beetles', 'jaw', 'dugongs', 'schools', 'seeds')
bugsFile <- examples
initsFile <- dataFile <- rep(NA, length(examples))
names(bugsFile) <- names(initsFile) <- names(dataFile) <- examples

## litters uses truncation on non-top node (issue #254) so left out for now
## bones, inhaler has dcat so left out for now (but we should be able to handle this)
## biops has stoch indexing so left out for now
## lsat requires setting some indexing and doing a bunch of initialization
## kidney, mice has dinterval so left out for now
## TODO: look at the various other bugs Volume 2 examples to see if there are others we should include


## customize file names as needed
bugsFile['beetles'] <- 'beetles-logit'
bugsFile['jaw'] <- 'jaw-constant'
bugsFile['epil'] <- 'epil2'
bugsFile <- paste0(bugsFile, '.bug')

initsFile['epil'] <- 'epil-inits.R'
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

relTol_default <- eval(formals(test_ADModelCalculate)$relTol)
relTols <- list()
length(relTols) <- length(examples)
names(relTols) <- examples



for(i in seq_along(examples)) {
    cat("Testing ", examples[i], ".\n")
    if(is.na(initsFile[i])) tmpInits <- NULL else tmpInits <- initsFile[i]
    if(is.na(dataFile[i])) tmpData <- NULL else tmpData <- dataFile[i]
    model <- readBUGSmodel(model = bugsFile[i], inits = tmpInits, data = tmpData, useInits = TRUE,
                           dir = nimble:::getBUGSexampleDir(examples[i]))
    if(!is.null(simulateNodes[[i]])) {
        model$simulate(simulateNodes[[i]])
    }
    out <- model$calculate()
    if(is.null(relTols[[examples[i]]])) {
        relTol <- relTol_default
    } else relTol = relTols[[examples[i]]]
    test_ADModelCalculate(model, useParamTransform = TRUE, verbose = TRUE, name = bugsFile[i], relTol = relTol, useFasterRderivs = TRUE)
}

## need to monkey with leuk and salm code, which makes it a pain to deal with directories as done above.

## leuk
writeLines(c("var", "Y[N,T],", "dN[N,T];"), con = file.path(tempdir(), "leuk.bug"))
system.in.dir(paste("cat leuk.bug >>", file.path(tempdir(), "leuk.bug")), dir = system.file('classic-bugs','vol1','leuk', package = 'nimble'))
system.in.dir(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
model <- readBUGSmodel(model = file.path(tempdir(), "leuk.bug"), data = system.file('classic-bugs','vol1','leuk','leuk-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','leuk','leuk-init.R', package = 'nimble'), useInits = TRUE)
out <- model$calculate()
## issue #253
test_ADModelCalculate(model, useParamTransform = TRUE, verbose = TRUE, name = 'leuk', relTol = relTol_default, useFasterRderivs = TRUE)

## salm
writeLines(c("var","logx[doses];"), con = file.path(tempdir(), "salm.bug"))
system.in.dir(paste("cat salm.bug >>", file.path(tempdir(), "salm.bug")), dir = system.file('classic-bugs','vol1','salm', package = 'nimble'))
model <- readBUGSmodel(model = file.path(tempdir(), "salm.bug"), data = system.file('classic-bugs','vol1','salm','salm-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','salm','salm-init.R', package = 'nimble'), useInits = TRUE)
model$simulate('lambda')
model$calculate()
## rWrt, cWrt equal not identical to x
test_ADModelCalculate(model, useParamTransform = TRUE, verbose = TRUE, name = 'salm', relTol = relTol_default, useFasterRderivs = TRUE)


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

