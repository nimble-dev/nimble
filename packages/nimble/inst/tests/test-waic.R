source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of calcWAIC")

###  BUGS models from Chapter 5 of Gelman and Hill
###  Below WAIC values from Gelman '13 "Understanding predictive 
###  information criteria for Bayesian models"

test_that("school model WAIC is accurate: ", {
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
  schoolSATmcmc <- buildMCMC(schoolSATmcmcConf)
  temporarilyAssignInGlobalEnv(schoolSATmcmc)
  CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel)
  CschoolSATmcmc$run(50000)
  expect_equal(CschoolSATmcmc$calculateWAIC(), 61.8, tolerance = 2.0)
})

test_that("voter model WAIC is accurate: ", {
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
  votermcmc <- buildMCMC(votermcmcConf)
  temporarilyAssignInGlobalEnv(votermcmc)
  Cvotermcmc <- compileNimble(votermcmc, project = voterModel)
  Cvotermcmc$run(50000)
  expect_equal(Cvotermcmc$calculateWAIC(), 87.2, tolerance = 2.0)
})


test_that("Radon model WAIC is accurate: ", {
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
  radonmcmc <- buildMCMC(radonmcmcConf)
  temporarilyAssignInGlobalEnv(radonmcmc)
  Cradonmcmc <- compileNimble(radonmcmc, project = radonModel)
  Cradonmcmc$run(10000)
  expect_equal(Cradonmcmc$calculateWAIC(1000), 3937, tolerance = 10)
})
