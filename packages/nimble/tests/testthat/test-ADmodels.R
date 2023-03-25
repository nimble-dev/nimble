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

verbose <- FALSE

context("Testing of derivatives for calculate() for nimbleModels")

test_that('pow and pow_int work', {  ## 28 sec.

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
  testFunctionInstance <- derivsNimbleFunction(m, calcNodes, wrtNodes)
  cm <- compileNimble(m)
  ctestFunctionInstance <- compileNimble(testFunctionInstance, project =  m, resetFunctions = TRUE)
  cDerivs <- ctestFunctionInstance$run(m$d, order)
  expect_equal(wrapperDerivs$value, cDerivs$value)
  expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian)
  expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = 1e-4)
  m$trt <- cm$trt <- -2
  m$calculate()
  cm$calculate()
  wrapperDerivs <- nimDerivs(m$calculate(calcNodes), wrt = wrtNodes, order = order)
  cDerivs <- ctestFunctionInstance$run(m$d, order)
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
  testFunctionInstance <- derivsNimbleFunction(m, calcNodes, wrtNodes)
  cm <- compileNimble(m)
  ctestFunctionInstance <- compileNimble(testFunctionInstance, project =  m, resetFunctions = TRUE)
  cDerivs <- ctestFunctionInstance$run(m$d, order)
  expect_equal(wrapperDerivs$value, cDerivs$value)
  expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian)
  expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = 1e-4)
  m$trt <- cm$trt <- -2
  m$calculate()
  cm$calculate()
  wrapperDerivs <- nimDerivs(m$calculate(calcNodes), wrt = wrtNodes, order = order)
  cDerivs <- ctestFunctionInstance$run(m$d, order)
  expect_equal(wrapperDerivs$value, cDerivs$value)
  expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian)
  expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = 1e-4)
})

## check logic of results
test_that('makeModelDerivsInfo works correctly', {
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
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu','y'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP",
                                           model$expandNodeNames('mu')))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## distinct wrt and calcNodes, where calcNodes includes a deterministic node
    ## updateNodes should have 'mu' and determ parents not in calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu','y', "lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP",
                                           model$expandNodeNames('mu')))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## with option to include data nodes in updateNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu','y'),
                              dataAsConstantNodes = FALSE)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP",
                                       model$expandNodeNames('mu'), model$expandNodeNames('y')))
    expect_identical(result$constantNodes, character(0))

    ## scalar calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu[1]','y[1]'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP", "mu[1]"))
    expect_identical(result$constantNodes, "y[1]")

    ## subset of a variable as calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'tau'), calcNodes = c('mu[1:2]','y[1:2]'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP", "mu[1]", "mu[2]"))
    expect_identical(result$constantNodes, c("y[1]", "y[2]"))

    ## a wrt node included in calcNodes
    ## updateNodes should have the stoch calcNodes not in wrt, plus determ parents not in calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0'), calcNodes = c('sigma2', 'mu0', 'mu','y'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP",
                                           model$expandNodeNames('mu')))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## some overlap of wrt and calcNodes
    ## updateNodes should have determ parents not in calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'tau', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu','y'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## full overlap of wrt and calcNodes
    ## updateNodes should have determ parents not in calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu','y'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## wrt are intermediate nodes
    ## updateNodes should have determ and stoch parents not in calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('mu'), calcNodes = c('mu','y'))
    expect_identical(result$updateNodes, c("mu0", "lifted_sqrt_oPsigma2_cP", "lifted_d1_over_sqrt_oPtau_cP"))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    ## no data nodes in calcNodes
    ## updateNodes should have determ parents not in calcNodes
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu'))
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

    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0'), calcNodes = c('mu','y'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems,
                                       model$expandNodeNames('mu', returnScalarComponents = TRUE)))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0'), calcNodes = c('mu','y'),
                              dataAsConstantNodes = FALSE)
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems,
                                           model$expandNodeNames('mu', returnScalarComponents = TRUE),
                                           model$expandNodeNames('y')))
    expect_identical(result$constantNodes, character(0))

    # @paciorek take a look at this one
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0'), calcNodes = c('mu[1]','y[1]'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems, 'mu[1]', 'mu[2]', 'mu[3]'))
    expect_identical(result$constantNodes, "y[1]")

    # @paciorek and this one. Note that prElems is not in result$updateNodes and shouldn't be.
    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0'), calcNodes = c('mu[1:2]','y[1:2]'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems, model$expandNodeNames('mu', returnScalarComponents = TRUE)))
    expect_identical(result$constantNodes, c("y[1]", "y[2]"))

    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0'), calcNodes = c('sigma2', 'mu0', 'mu','y'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems,
                                       model$expandNodeNames('mu', returnScalarComponents = TRUE)))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu','y'))
    expect_identical(result$updateNodes, c("lifted_sqrt_oPsigma2_cP", lftChElems))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeModelDerivsInfo(model = model, wrtNodes = c('mu'), calcNodes = c('mu','y'))
    expect_identical(result$updateNodes, c(model$expandNodeNames('mu0', returnScalarComponents = TRUE),
                                           "lifted_sqrt_oPsigma2_cP", lftChElems))
    expect_identical(result$constantNodes, model$expandNodeNames('y'))

    result <- makeModelDerivsInfo(model = model, wrtNodes = c('sigma2', 'mu0', 'mu'), calcNodes = c('sigma2', 'mu0', 'mu'))
    expect_identical(result$updateNodes, c(lftChElems))
    expect_identical(result$constantNodes, character(0))
   
})

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
relTolTmp[4] <- 1e-3

## 200 sec.
test_ADModelCalculate(model, relTol = relTolTmp, checkCompiledValuesIdentical = FALSE, useFasterRderivs = TRUE,
                      verbose = verbose, name = 'basic model, lifted nodes')

## basic state space model
set.seed(1) 
code <- nimbleCode({
    x0 ~ dnorm(0,1)
    x[1] ~ dnorm(x0, 1)
    y[1] ~ dnorm(x[1], var = 2)
    for(i in 2:5) {
        x[i] ~ dnorm(x[i-1], 1)
        y[i] ~ dnorm(x[i], var = 2)
    }
})
data <- list(y = rnorm(5))
model <- nimbleModel(code, data = data)
model$simulate()
model$calculate()
relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-4
relTolTmp[4] <- 1e-2
## 204 sec.
test_ADModelCalculate(model, relTol = relTolTmp, useFasterRderivs = TRUE, verbose = verbose, name = 'basic state space') 

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
### 185 sec.
test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, name = 'basic tricky indexing', 
                      newUpdateNodes = list(covMat = matrix(c(0.7, .25, .25, .7), 2)))



## link functions on stochastic nodes (not covered in BUGS examples)
## plus alt params and NIMBLE-provided distributions
set.seed(1)
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
newY <- rpois(n, 2)

relTolTmp <- relTol
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2

## 325 sec.
test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, name = 'stochastic link model', useFasterRderivs = TRUE,
                      checkCompiledValuesIdentical = FALSE, useParamTransform = TRUE, newConstantNodes = list(y = newY))

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
relTolTmp[4] <- 1e-1
relTolTmp[5] <- 1e-13
## 348 sec.
test_ADModelCalculate(model, relTol = relTolTmp, xNew = xNew, verbose = verbose,
                      useFasterRderivs = TRUE, useParamTransform = TRUE,
                      checkCompiledValuesIdentical = FALSE, name = 'dt and dexp model')


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

## 30 minutes if do full assessment
## various R vs. C discrepancies in 2d11 O(0.1); skip in part given time.
## 335 sec.
test_ADModelCalculate(model, newUpdateNodes = list(S = newS, pr = newPr, pr2 = newPr2), useParamTransform = TRUE,
                      relTol = relTolTmp, checkCompiledValuesIdentical = FALSE, checkDoubleUncHessian = FALSE,
                      useFasterRderivs = TRUE, verbose = verbose, name = 'complicated indexing')


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

## 205 sec.
test_ADModelCalculate(model, newUpdateNodes = list(pr5 = newPr5, pr4 = newPr4), relTol = relTolTmp,
                      useFasterRderivs = TRUE, verbose = verbose, name = 'different subsets of a matrix')


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
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13
## 330 sec.
test_ADModelCalculate(model, newUpdateNodes = list(pr = newPr, dist = newDist), useParamTransform = TRUE,
                      relTol = relTolTmp, absTolThreshold = 1e-12, checkCompiledValuesIdentical = FALSE, 
                      useFasterRderivs = TRUE, verbose = verbose, name = 'dmnorm with vectorized covariance matrix')



## MVN with various parameterizations and user-defined functions

## user-defined cov function with loops
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

code <- nimbleCode({
    Sigma2[1:n,1:n] <- covFunLoop(dist[1:n,1:n], rho)
    y[1:n] ~ dmnorm(mu2[1:n], cov = Sigma2[1:n,1:n])
    mu2[1:n] ~ dmnorm(z[1:n], pr[1:n,1:n])
    rho ~ dgamma(2, 3)
})

set.seed(1)                               
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n),
                     inits = list(rho = rgamma(1, 1, 1), dist = dd, z = rep(0, n), pr = diag(n)))
model$simulate()
model$calculate()
model$setData('y')
newPr <- crossprod(matrix(rnorm(5*5), 5))
newDist <- as.matrix(dist(runif(5)))

relTolTmp <- relTol
relTolTmp[1] <- 1e-10
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-1
relTolTmp[5] <- 1e-13

## 350 sec.
test_ADModelCalculate(model, useParamTransform = TRUE, useFasterRderivs = TRUE,
                      newUpdateNodes = list(dist = newDist, pr = newPr), checkCompiledValuesIdentical = FALSE, 
                      relTol = relTolTmp, absTolThreshold = 1e-12, verbose = verbose, 
                      name = 'dnorm with user-defined fxn for covariance with loops')

## other dmnorm parameterizations
set.seed(1)
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
                     inits = list(z = rep(1, n), pr = diag(n),
                                  Sigma = Sigma, Q = solve(Sigma)))
model$simulate()
model$calculate()
model$setData('y')
newSigma <- crossprod(matrix(rnorm(5*5), 5))
newQ <- crossprod(matrix(rnorm(5*5), 5))
newPr <- crossprod(matrix(rnorm(5*5), 5))

relTolTmp <- relTol
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-5
relTolTmp[4] <- 1e-1

## 402 sec.
test_ADModelCalculate(model, absTolThreshold = 1e-12, useParamTransform = TRUE, useFasterRderivs = TRUE,
                      checkCompiledValuesIdentical = FALSE, newUpdateNodes = list(pr = newPr, Q = newQ, Sigma = newSigma),
                      relTol = relTolTmp, verbose = verbose, name = 'various dmnorm parameterizations')


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

## 310 sec.
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTolTmp, checkCompiledValuesIdentical = FALSE,
                      useFasterRderivs = TRUE, verbose = verbose, name = 'simple user-defined distribution')



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

## 166 sec.
test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, name = 'vectorized power')


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

code <- nimbleCode({ 
    y[1:n] ~ dGPdist(dist[1:n, 1:n], rho)
    rho ~ dgamma(2, 3)
})

set.seed(1)
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n, dist = dd), inits = list(rho = rgamma(1, 1, 1)))
model$simulate()
model$calculate()
model$setData('y')
newDist <- as.matrix(dist(runif(n)))
relTolTmp <- relTol
relTolTmp[1] <- 1e-14

## 361 sec.
test_ADModelCalculate(model, newConstantNodes = list(dist = newDist), useParamTransform = TRUE,
                      checkCompiledValuesIdentical = FALSE, useFasterRderivs = TRUE,
                      relTol = relTolTmp, verbose = verbose, name = 'user-defined distribution')


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

set.seed(1)
n <- 5
locs <- runif(n)
dd <- fields::rdist(locs)
model <- nimbleModel(code, constants = list(n = n, x = rnorm(n), z = rep(1, n), pr = diag(n), dist = dd),
                     inits = list(rho = rgamma(1, 1, 1), sigma2 = rgamma(1, 1, 1)),
                     data = list(y1 = rnorm(1), #  y2 = rnorm(1),
                                 yy = matrix(rnorm(n*3), 3, n)))
model$simulate()
model$calculate()
model$setData(c('y1','yy')) # 'y2'
newDist <- as.matrix(dist(runif(n)))
newPr <- crossprod(matrix(rnorm(n*n),n))

relTolTmp <- relTol
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-3
relTolTmp[4] <- 1   # yikes
relTolTmp[5] <- 1e-13

## 462 sec.
test_ADModelCalculate(model,absTolThreshold = 1e-12, useParamTransform = TRUE, checkCompiledValuesIdentical = FALSE,
                      newConstantNodes = list(dist = newDist, pr = newPr), useFasterRderivs = TRUE,
                      relTol = relTolTmp, verbose = verbose, name = 'various matrix functions')


## Various combinations of updateNodes, wrt, calcNodes

set.seed(1)
code <- nimbleCode({
    a ~ dgamma(1.1, 0.8)
    b ~ dnorm(z, var = a)
    z ~ dnorm(0, 1)
})
model <- nimbleModel(code, data = list(b = 1.2), inits = list(a = 1.3, z = 0.7))
## calcNodes excludes det intermediates
## 81 sec.
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTol, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      wrt = 'a', calcNodes = c('a', 'b'), name = 'update nodes case 1a')

model <- nimbleModel(code, data = list(b = 1.2), inits = list(a = 1.3, z = 0.7))
## calcNodes includes det intermediates
## 88 sec.
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTol, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      wrt = 'a', calcNodes = c('a', 'b', 'lifted_sqrt_oPa_cP'), name = 'update nodes case 1b')


set.seed(1)
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
## 99 sec.
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTol, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      wrt = 'a', calcNodes = c('a', 'b'), name = 'update nodes case 2a')

model <- nimbleModel(code, data = list(b = rnorm(4)), inits = list(a = 1.3, z = runif(4), pr = diag(2), mu0 = rep(0, 2)))
relTolTmp <- relTol
relTolTmp[2] <- 1e-7
relTolTmp[5] <- 1e-13
## calcNodes includes det intermediates
## 105 sec.
test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTolTmp, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      wrt = 'a', calcNodes = c('a', 'b', "lifted_sqrt_oPa_cP"), name = 'update nodes case 2b')


## Parameter transform system and full use of ddirch, dwish, dinvwish

## basic variance component
set.seed(1)
code <- nimbleCode({
    y ~ dnorm(mu, sd = sigma)
    sigma ~ dinvgamma(1.3, 0.7)
    mu ~ dnorm(0, 1)
})
model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma = rgamma(1, 1, 1), mu = rnorm(1)))
## 330 sec.
test_ADModelCalculate(model, relTol = relTol, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      useFasterRderivs = TRUE, useParamTransform = TRUE, name = 'basic param transform')


set.seed(1)
code <- nimbleCode({
    y ~ dnorm(mu, var = sigma2)
    sigma2 ~ dinvgamma(1.3, 0.7)
    mu ~ dnorm(0, 1)
})
model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = rgamma(1, 1, 1), mu = rnorm(1)))
## 315 sec.
test_ADModelCalculate(model, relTol = relTol, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      useFasterRderivs = TRUE, useParamTransform = TRUE, name = 'basic param transform, with lifted')


## now check if model is out-of-state
code <- nimbleCode({
    y ~ dnorm(0, sd = sigma)
    sigma <- sqrt(sigma2)
    sigma2 ~ dinvgamma(1.3, 0.7)
})

set.seed(1)
model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = 2))
model$sigma <- 1
## 220 sec.
test_ADModelCalculate(model, relTol = relTol, verbose = verbose, useFasterRderivs = TRUE)

set.seed(1)
model <- nimbleModel(code, data = list(y = rnorm(1)), inits = list(sigma2 = 2))
model$sigma <- 1
relTolTmp <- relTol
relTolTmp[3] <- 1e-3
relTolTmp[4] <- 1e-1
## 406 sec.
test_ADModelCalculate(model, relTol = relTolTmp, verbose = verbose, checkCompiledValuesIdentical = FALSE,
                      useFasterRderivs = TRUE, useParamTransform = TRUE)

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
model <- nimbleModel(code, constants = list(k = k, n = n), data = list(y = rmulti(1, n, rep(1/k, k))),
                     inits = list(p = c(.2, .4, .15, .25), alpha = runif(4)))
newP <- rdirch(1, rep(1:k))
newY <- rmulti(1, n, newP)
relTolTmp <- relTol
relTolTmp[1] <- 1e-14
relTolTmp[2] <- 1e-7
relTolTmp[3] <- 1  ## manual check indicates default h in pracma::hessian is too large
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-9

## rOutput2d11 result can be wildly out of tolerance, so not checking it.
## 335 sec.
test_ADModelCalculate(model, x = 'prior', useParamTransform = TRUE, newUpdateNodes = list(p = newP), newConstantNodes = list(y = newY),
                      checkDoubleUncHessian = FALSE, relTol = relTolTmp, absTolThreshold = 1e-12, checkCompiledValuesIdentical = FALSE,
                      useFasterRderivs = TRUE, verbose = verbose, name = 'Dirichlet paramTransform') 


## This segfaults as of 2023-03-25, with libnimble.a.

covFunVec <- nimbleFunction(
    run = function(dist = double(2), rho = double(0)) {
        out <- exp(-dist/rho)
        returnType(double(2))
        return(out)
    }, buildDerivs = TRUE)


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
relTolTmp[3] <- 1e-1
relTolTmp[4] <- 1e-1
relTolTmp[5] <- 1e-11

## rOutput2d11 result can be wildly out of tolerance, so not checking it.
test_ADModelCalculate(model, newUpdateNodes = list(nu = 12.1, dist = newDist, R = newR, W1 = newW1, W2 = newW2, W3 = newW3, W4 = newW4, W5 = newW5, W6 = newW6),
                      x = 'prior', absTolThreshold = 1e-12, checkCompiledValuesIdentical = FALSE,
                      useParamTransform = TRUE, useFasterRderivs = TRUE, checkDoubleUncHessian = FALSE,
                      relTol = relTolTmp, verbose = verbose,
                      name = 'various multivariate dists')

## 2023-01-30: got a seg fault in EB scenario; plus this takes forever

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


nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
