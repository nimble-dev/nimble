source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context("Testing of calcWAIC")

###  BUGS models from Chapter 5 of Gelman and Hill
###  Below WAIC values from Gelman '13 "Understanding predictive 
###  information criteria for Bayesian models"

test_that("school model WAIC is accurate", {
  sigma     <- c(15,10,16,11, 9,11,10,18)
  schoolobs <- c(28,8, -3, 7,-1, 1,18,12)
  schoolSATcode <- nimbleCode({
    for(i in 1:N) {
      schoolmean[i] ~ dnorm(mu,itau)
      thes[i] <- 1/(sigma[i])^2
      schoolobs[i] ~ dnorm(schoolmean[i],thes[i])
    }
    mu ~ dnorm(0,0.1) 
    itau   ~ dgamma(1e-3,0.225)
  })
  schoolSATmodel <- nimbleModel(code = schoolSATcode,
                                data=list(sigma = sigma,
                                          schoolobs =schoolobs),
                                constants = list(N = length(schoolobs)),
                                inits = list(schoolmean = rep(3, 8),
                                             mu = 3,
                                             itau = .2))
  temporarilyAssignInGlobalEnv(schoolSATmodel)
  compileNimble(schoolSATmodel)
  schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('schoolmean'))
  schoolSATmcmc <- buildMCMC(schoolSATmcmcConf, enableWAIC = TRUE)
  temporarilyAssignInGlobalEnv(schoolSATmcmc)
  CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel)
  CschoolSATmcmc$run(50000)
  expect_equal(CschoolSATmcmc$calculateWAIC(), 61.8, tolerance = 2.0)

  schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('mu', 'itau'))
  expect_error(buildMCMC(schoolSATmcmcConf, enableWAIC = TRUE), 
               "To calculate WAIC in NIMBLE, all parameters")
  ## mWAIC not enabled as of 0.10.1
  ## schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('mu'))
  ## expect_error(buildMCMC(schoolSATmcmcConf, enableWAIC = TRUE))
  ## different set of monitors than above, so different waic value expected
  ## expect_equal(nimbleMCMC(model = schoolSATmodel, WAIC = TRUE)$WAIC, 67, tolerance = 8) 

  schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('schoolmean'))
  schoolSATmcmc <- buildMCMC(schoolSATmcmcConf, enableWAIC = TRUE)
  CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel, 
                                  resetFunctions = TRUE)
  expect_equal(runMCMC(CschoolSATmcmc, WAIC = TRUE)$WAIC, 61.8, 
               tolerance = 2) 
  schoolSATmcmc <- buildMCMC(schoolSATmcmcConf, enableWAIC = FALSE)
  expect_error(runMCMC(schoolSATmcmc, WAIC = TRUE))
})

test_that("voter model WAIC is accurate", {
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
  temporarilyAssignInGlobalEnv(voterModel)
  CvoterModel <- compileNimble(voterModel)
  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_1', 'beta_2',
                                                          'sigma'))
  votermcmc <- buildMCMC(votermcmcConf, enableWAIC = TRUE)
  temporarilyAssignInGlobalEnv(votermcmc)
  Cvotermcmc <- compileNimble(votermcmc, project = voterModel)
  Cvotermcmc$run(50000)
  expect_equal(Cvotermcmc$calculateWAIC(), 87.2, tolerance = 2.0)
  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_2',
                                                          'sigma'))
  expect_error(buildMCMC(votermcmcConf, enableWAIC = TRUE))
  
  ## additional testing of validity of monitored nodes below
  voterCode = nimbleCode({
    for(i in 1:N){
      y[i] ~ dnorm(beta_12 + growth[i]*beta_2, sd = sigma_2)
    }
    beta_12 <- beta_1*2
    sigma_2 ~ dexp(sigma*3)
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

  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_1',
                                                          'beta_2',
                                                          'sigma_2'))
  expect_message(buildMCMC(votermcmcConf, enableWAIC = TRUE), 
                 "Monitored nodes are valid for WAIC",
                 all = FALSE, fixed = TRUE)

  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_12',
                                                          'beta_2',
                                                          'sigma'))
  expect_error(buildMCMC(votermcmcConf, enableWAIC = TRUE), 
                 "To calculate WAIC in NIMBLE, all parameters")
  
  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_12',
                                                          'beta_2',
                                                          'sigma_2'))
  expect_error(buildMCMC(votermcmcConf, enableWAIC = TRUE),
                 "To calculate WAIC in NIMBLE, all parameters")

  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_12',
                                                          'beta_2',
                                                          'sigma'))
  expect_error(buildMCMC(votermcmcConf, enableWAIC = TRUE), 
                 "To calculate WAIC in NIMBLE, all parameters")

  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_1',
                                                          'beta_2'))
  expect_error(buildMCMC(votermcmcConf, enableWAIC = TRUE),
               "To calculate WAIC in NIMBLE, all parameters")
})


test_that("Radon model WAIC is accurate", {
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
  temporarilyAssignInGlobalEnv(radonModel)
  CradonModel <- compileNimble(radonModel)
  radonmcmcConf <- configureMCMC(radonModel, monitors = c('beta'))
  radonmcmc <- buildMCMC(radonmcmcConf, enableWAIC = TRUE)
  temporarilyAssignInGlobalEnv(radonmcmc)
  Cradonmcmc <- compileNimble(radonmcmc, project = radonModel)
  Cradonmcmc$run(10000)
  expect_equal(Cradonmcmc$calculateWAIC(1000), 3937, tolerance = 10)
  radonmcmcConf <- configureMCMC(radonModel, monitors = c('coefs'))
  ## check to ensure monitoring deterministic nodes works
  expect_message(radonmcmc <- buildMCMC(radonmcmcConf, enableWAIC = TRUE), 
                 "Monitored nodes are valid for WAIC",
                 all = FALSE, fixed = TRUE)
  Cradonmcmc <- compileNimble(radonmcmc, project = radonModel,
                              resetFunctions = TRUE)
  Cradonmcmc$run(10000)
  ## monitoring coefs is equivalent to monitoring beta, so waic should match
  expect_equal(Cradonmcmc$calculateWAIC(1000), 3937, tolerance = 10)
})

nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
