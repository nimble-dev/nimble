source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of calcWAIC")

sigma     <- c(15,10,16,11, 9,11,10,18)
schoolobs <- c(28,8, -3, 7,-1, 1,18,12)

schoolSATcode <- nimbleCode({
  for(i in 1:N) {
    schoolmean[i] ~ dnorm(mu,itau)
    thes[i] <- 1/(sigma[i])^2
    schoolobs[i] ~ dnorm(schoolmean[i],thes[i])
  }
  mu ~ dnorm(0,alpha)
  alpha <- .01
  itau   ~ dgamma(1e-3,0.225)
  tau <- sqrt(1/itau)
})

schoolSATmodel <- nimbleModel(code = schoolSATcode,
                              data=list(sigma = sigma,
                                        schoolobs =schoolobs),
                              constants = list(N = length(schoolobs)))

compileNimble(schoolSATmodel)

schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('schoolmean'))
schoolSATmcmc <- buildMCMC(schoolSATmcmcConf)
CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel)
CschoolSATmcmc$run(50000)
schoolWAICfunc <- calcWAIC(schoolSATmodel, schoolSATmcmc)
CschoolWAICfunc <- compileNimble(schoolWAICfunc, project = schoolSATmodel)
schoolWAIC <- CschoolWAICfunc$run(20000)

##below WAIC value from Gelman '13 "Understanding predictive information criteria for Bayesian models"
try(test_that("Test that WAIC is accurate: ",
              expect_equal(schoolWAIC, 61.8, tolerance = 2.0)))

