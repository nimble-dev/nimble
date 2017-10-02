source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of cross-validation")

###  BUGS model from Chapter 5 of Gelman and Hill
###  Below LOO-CV value from Gelman '13 "Understanding predictive 
###  information criteria for Bayesian models"

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
                              MCMCcontrol = list(nMCMCiters = 500),
                              NbootReps = 2)
  expect_equal(testOut$CVvalue, -43.8, tolerance = 2)
})



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
                              k = 2, ## do loo-cv so that we can compare against a known value
                              foldFunction = 'random',
                              lossFunction = 'MSE',
                              MCMCcontrol = list(nMCMCiters = 500),
                              NbootReps = 2)
  expect_equal(testOut$CVvalue, -43.8, tolerance = 2)
})



