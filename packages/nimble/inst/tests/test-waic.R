source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of calcWAIC")


###  BUGS code from Chapter 5 of Gelman and Hill
###  Below WAIC value from Gelman '13 "Understanding predictive information criteria for Bayesian models"

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
                              constants = list(N = length(schoolobs)))

compileNimble(schoolSATmodel)

schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('schoolmean'))
schoolSATmcmc <- buildMCMC(schoolSATmcmcConf)
CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel)
test_that("Test that WAIC returns -Inf if too few posterior samples are used: ",
          expect_equal(CschoolSATmcmc$calculateWAIC(), -Inf))
CschoolSATmcmc$run(50000)
test_that("Test that WAIC is accurate: ",
          expect_equal(CschoolSATmcmc$calculateWAIC(), 61.8, tolerance = 2.0))

