source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

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
  expect_equal(testOut$CVvalue*15, -43.8, tolerance = 2)
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


options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
