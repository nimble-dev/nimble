source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

context("Testing of cross-validation")

###  BUGS model from Chapter 5 of Gelman and Hill
###  Below LOO-CV value from Gelman '13 "Understanding predictive 
###  information criteria for Bayesian models"

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = TRUE) ## try to prevent timeout on Travis due to lack of output

test_that("voter model cross-validation is accurate: ", {
  y <- c(44.6, 57.76, 49.91, 61.34, 49.60, 61.79, 48.95, 44.70, 59.17, 53.94,
         46.55, 54.74, 50.27, 51.24, 46.32)
  growth <- c(2.4, 2.89, .85, 4.21, 3.02, 3.62, 1.08, -.39, 3.86, 2.27, .38,
              1.04, 2.36, 1.72, .1)
  N = length(y)
  voterCode = nimbleCode({
    for(i in 1:N){
      y[i] ~ dnorm(beta_1 + growth[i]*beta_2, sd = sigma)
    }
    sigma ~ dunif(0,50)
    beta_1 ~ dnorm(0, .01)
    beta_2 ~ dnorm(0, .01)
  })
  dataList = list(y = y);
  constList = list(growth = growth, N = N)
  voterModel = nimbleModel(code = voterCode, data = dataList, 
                           constants = constList,
                           inits = list(beta_1 = 44, beta_2 = 3.75, 
                                        sigma = 4.4))
  voterConf <- configureMCMC(voterModel)
  testOut <- runCrossValidate(MCMCconfiguration = voterConf,
                              k = 15, ## do loo-cv so that we can compare against a known value
                              foldFunction = 'random',
                              lossFunction = 'predictive',
                              MCMCcontrol = list(niter = 500),
                              nBootReps = 2)
  expect_lt(abs(testOut$CVvalue*15 + 43.8), 2)
})

## Next two tests are just to ensure that the CV algo runs error-free
## with different loss functions / on different models.

test_that("voter model cross-validation runs using MSE loss: ", {
  y <- c(44.6, 57.76, 49.91, 61.34, 49.60, 61.79, 48.95, 44.70, 59.17, 53.94,
         46.55, 54.74, 50.27, 51.24, 46.32)
  growth <- c(2.4, 2.89, .85, 4.21, 3.02, 3.62, 1.08, -.39, 3.86, 2.27, .38,
              1.04, 2.36, 1.72, .1)
  N = length(y)
  voterCode = nimbleCode({
    for(i in 1:N){
      y[i] ~ dnorm(beta_1 + growth[i]*beta_2, sd = sigma)
    }
    sigma ~ dunif(0,50)
    beta_1 ~ dnorm(0, .01)
    beta_2 ~ dnorm(0, .01)
  })
  dataList = list(y = y);
  constList = list(growth = growth, N = N)
  voterModel = nimbleModel(code = voterCode, data = dataList, 
                           constants = constList,
                           inits = list(beta_1 = 44, beta_2 = 3.75, 
                                        sigma = 4.4))
  voterConf <- configureMCMC(voterModel)
  testOut <- runCrossValidate(MCMCconfiguration = voterConf,
                              k = 2, 
                              foldFunction = 'random',
                              lossFunction = 'MSE',
                              MCMCcontrol = list(niter = 100),
                              nBootReps = 2,
                              silent = TRUE)
})

test_that("cross-validation runs correctly when nBootReps = NA: ", {
  y <- c(44.6, 57.76, 49.91, 61.34, 49.60, 61.79, 48.95, 44.70, 59.17, 53.94,
         46.55, 54.74, 50.27, 51.24, 46.32)
  growth <- c(2.4, 2.89, .85, 4.21, 3.02, 3.62, 1.08, -.39, 3.86, 2.27, .38,
              1.04, 2.36, 1.72, .1)
  N = length(y)
  voterCode = nimbleCode({
    for(i in 1:N){
      y[i] ~ dnorm(beta_1 + growth[i]*beta_2, sd = sigma)
    }
    sigma ~ dunif(0,50)
    beta_1 ~ dnorm(0, .01)
    beta_2 ~ dnorm(0, .01)
  })
  dataList = list(y = y);
  constList = list(growth = growth, N = N)
  voterModel = nimbleModel(code = voterCode, data = dataList, 
                           constants = constList,
                           inits = list(beta_1 = 44, beta_2 = 3.75, 
                                        sigma = 4.4))
  voterConf <- configureMCMC(voterModel)
  testOut <- runCrossValidate(MCMCconfiguration = voterConf,
                              k = 2, 
                              foldFunction = 'random',
                              lossFunction = 'MSE',
                              MCMCcontrol = list(niter = 100),
                              nBootReps = NA)
  expect_true(is.na(testOut$CVstandardError))
})


test_that("Radon model cross-validation runs using MSE loss: ", {
  url <- "http://stat.columbia.edu/~gelman/arm/examples/arsenic/wells.dat"
  wells <- try(read.table(url), silent = TRUE)
  if (inherits(wells, 'try-error')) skip("No internet connection")
  wells$dist100 <- with(wells, dist / 100)
  X <- model.matrix(~ dist100 + arsenic, wells)
  radonCode = nimbleCode({
    for(i in 1:3){
      beta[i] ~ dnorm(0, 1)
    }
    coefs[1:N] <- X[1:N, 1:3] %*% beta[1:3]
    for(i in 1:N){
      y[i] ~ dbern(ilogit(coefs[i]))
    }
  })
  dataList = list(y = wells$switch);
  constList = list(N = nrow(X), 
                   X = X)
  radonModel = nimbleModel(code = radonCode, data = dataList, 
                           constants = constList,
                           inits = list(beta = rep(1, 3)))
  radonConf <- configureMCMC(radonModel)
  radonConf$addSampler('beta', type = 'AF_slice')
  testOut <- runCrossValidate(MCMCconfiguration = radonConf,
                              k = 2, 
                              foldFunction = 'random',
                              lossFunction = 'MSE',
                              MCMCcontrol = list(niter = 100),
                              nBootReps = 2)
})

test_that("dyes model cross-validation exact results", {
    dyesCode <- nimbleCode({
        for (i in 1:BATCHES) {
            for (j in 1:SAMPLES) {
                y[i,j] ~ dnorm(mu[i], tau.within);
            }
            mu[i] ~ dnorm(theta, tau.between);
        }
        theta ~ dnorm(0.0, 1.0E-10);
        tau.within ~ dgamma(0.001, 0.001);  sigma2.within <- 1/tau.within;
        tau.between ~ dgamma(0.001, 0.001);  sigma2.between <- 1/tau.between;
    })
    ##
    dyesData <- list(y = matrix(c(1545, 1540, 1595, 1445, 1595,
                         1520, 1440, 1555, 1550, 1440,
                         1630, 1455, 1440, 1490, 1605,
                         1595, 1515, 1450, 1520, 1560,
                         1510, 1465, 1635, 1480, 1580,
                         1495, 1560, 1545, 1625, 1445),
                         nrow = 6, ncol = 5))
    ##
    dyesConsts <- list(BATCHES = 6, SAMPLES = 5)
    ##
    dyesInits <- list(theta = 1500, tau.within = 1, tau.between =  1)
    ##
    dyesModel <- nimbleModel(dyesCode, dyesConsts, dyesData, dyesInits)
    ##
    dyesFoldFunction <- function(i){
        foldNodes_i <- paste0('y[', i, ', ]')  # will return 'y[1,]' for i = 1 e.g.
        return(foldNodes_i)
    }
    ##
    RMSElossFunction <- function(simulatedDataValues, actualDataValues){
        dataLength <- length(simulatedDataValues) # simulatedDataValues is a vector
        SSE <- 0
        for(i in 1:dataLength){
            SSE <- SSE + (simulatedDataValues[i] - actualDataValues[i])^2
        }
        MSE <- SSE / dataLength
        RMSE <- sqrt(MSE)
        return(RMSE)
    }
    ##
    dyesMCMCconfiguration <- configureMCMC(dyesModel)
    ##
    set.seed(0)
    ##
    crossValOutput <- runCrossValidate(MCMCconfiguration = dyesMCMCconfiguration, k = 6,
                                       foldFunction = dyesFoldFunction,
                                       lossFunction = RMSElossFunction,
                                       MCMCcontrol = list(niter = 5000, nburnin = 500))
    ##
    expect_equal(crossValOutput$CVvalue, 63.872534)
    expect_equal(crossValOutput$CVstandardError, 0.0023289839)
    expect_equal(unname(sapply(crossValOutput$foldCVinfo, function(x) x['foldCVvalue'])), c(56.525316, 41.326957, 73.621281, 61.466634, 109.806489, 40.488530))
    expect_equal(unname(sapply(crossValOutput$foldCVinfo, function(x) x['foldCVstandardError'])), c(0.0052012467, 0.0058034916, 0.0060221136, 0.0059338851, 0.0045952414, 0.0064763733))
})



options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
