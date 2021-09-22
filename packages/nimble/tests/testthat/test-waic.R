source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context("Testing of WAIC")

###  BUGS models from Chapter 5 of Gelman and Hill
###  Below WAIC values from Gelman '13 "Understanding predictive 
###  information criteria for Bayesian models"

test_that("school model WAIC is accurate using original WAIC implementation", {
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
  schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('schoolmean'),
                                     enableWAIC = TRUE)
  schoolSATmcmc <- buildMCMC(schoolSATmcmcConf)
  temporarilyAssignInGlobalEnv(schoolSATmcmc)
  CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel)
  CschoolSATmcmc$run(50000)
  expect_lt(abs(CschoolSATmcmc$getWAIC()$WAIC - 61.8), 2.0)

  schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('mu', 'itau'),
                                     enableWAIC = TRUE, controlWAIC = list(online = FALSE))
  expect_error(buildMCMC(schoolSATmcmcConf), 
               "To calculate WAIC in NIMBLE, all parameters")
  ## mWAIC not enabled as of 0.10.1
  ## schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('mu'))
  ## expect_error(buildMCMC(schoolSATmcmcConf, enableWAIC = TRUE))
  ## different set of monitors than above, so different waic value expected
  ## expect_lt(abs(nimbleMCMC(model = schoolSATmodel, WAIC = TRUE)$WAIC - 67, 8) 

  schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('schoolmean'),
                                     enableWAIC = TRUE)
  schoolSATmcmc <- buildMCMC(schoolSATmcmcConf)
  CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel, 
                                  resetFunctions = TRUE)
  expect_lt(abs(runMCMC(CschoolSATmcmc, WAIC = TRUE)$WAIC$WAIC - 61.8), 2) 
  schoolSATmcmcConf <- configureMCMC(schoolSATmodel, monitors = c('schoolmean'))
  schoolSATmcmc <- buildMCMC(schoolSATmcmcConf)
  CschoolSATmcmc <- compileNimble(schoolSATmcmc, project = schoolSATmodel, 
                                  resetFunctions = TRUE)
  expect_error(runMCMC(CschoolSATmcmc, WAIC = TRUE), "mcmc argument must have been created")
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
  votermcmcConf <- configureMCMC(voterModel, monitors =
            c('beta_1', 'beta_2', 'sigma'), enableWAIC = TRUE)
  votermcmc <- buildMCMC(votermcmcConf)
  temporarilyAssignInGlobalEnv(votermcmc)
  Cvotermcmc <- compileNimble(votermcmc, project = voterModel)
  Cvotermcmc$run(50000)
  expect_lt(abs(Cvotermcmc$getWAIC()$WAIC - 87.2), 2.0)
  
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
                                                          'sigma_2'),
                                 enableWAIC = TRUE,
                                 controlWAIC = list(online = FALSE))
  expect_message(buildMCMC(votermcmcConf), 
                 "Monitored nodes are valid for WAIC",
                 all = FALSE, fixed = TRUE)

  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_12',
                                                          'beta_2',
                                                          'sigma'),
                                 enableWAIC = TRUE,
                                 controlWAIC = list(online = FALSE))
  expect_error(buildMCMC(votermcmcConf), 
                 "To calculate WAIC in NIMBLE, all parameters")
  
  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_12',
                                                          'beta_2',
                                                          'sigma_2'),
                                 enableWAIC = TRUE,
                                 controlWAIC = list(online = FALSE))
  expect_error(buildMCMC(votermcmcConf),
                 "To calculate WAIC in NIMBLE, all parameters")

  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_12',
                                                          'beta_2',
                                                          'sigma'),
                                 enableWAIC = TRUE,
                                 controlWAIC = list(online = FALSE))
  expect_error(buildMCMC(votermcmcConf), 
                 "To calculate WAIC in NIMBLE, all parameters")

  votermcmcConf <- configureMCMC(voterModel, monitors = c('beta_1',
                                                          'beta_2'),
                                 enableWAIC = TRUE,
                                 controlWAIC = list(online = FALSE))
  expect_error(buildMCMC(votermcmcConf),
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
  radonmcmcConf <- configureMCMC(radonModel, monitors = c('beta'), enableWAIC = TRUE)
  print(radonmcmcConf)
  radonmcmc <- buildMCMC(radonmcmcConf)
  temporarilyAssignInGlobalEnv(radonmcmc)
  Cradonmcmc <- compileNimble(radonmcmc, project = radonModel)
  Cradonmcmc$run(10000, nburnin = 1000)
  expect_lt(abs(Cradonmcmc$getWAIC()$WAIC - 3937), 2)
  monitors <- radonModel$getNodeNames(determOnly = TRUE)
  radonmcmcConf <- configureMCMC(radonModel, monitors = monitors, enableWAIC = TRUE,
                                 controlWAIC = list(online = FALSE))
  ## check to ensure monitoring deterministic nodes fails
  expect_error(radonmcmc <- buildMCMC(radonmcmcConf), 
                 "To calculate WAIC in NIMBLE, all parameters")
})


test_that("New WAIC implementation matches old implementation for conditional, ungrouped", {
    set.seed(1)
    J <- 5
    I <- 10
    tau <- 1
    sigma <- 1
    
    mu <- rnorm(J, 0 , tau)
    
    y <- matrix(0, J, I)
    for(j in 1:J) 
        y[j, ] <- rnorm(I, mu[j], sigma)

    code <- nimbleCode({
        tau ~ dunif(0, 10)
        sigma ~ dunif(0, 10)
        mu0 ~ dnorm(0, 10)
        for(j in 1:J) {
            mu[j] ~ dnorm(mu0, sd = tau)
            for(i in 1:I)
                y[j, i] ~ dnorm(mu[j], sd = sigma)
        }
    })

    inits <- list(mu0 = 0, tau = 0.5, sigma = 1.5,
                               mu = rnorm(J, 0, 0.5))

    m <- nimbleModel(code, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits = inits)
    cm <- compileNimble(m)

    mcmc <- buildMCMC(m, enableWAIC = TRUE, controlWAIC = list(online = FALSE),
                  monitors = c('mu0','mu','sigma','tau'))
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    out1 <- runMCMC(cmcmc, niter = 1000, WAIC = TRUE, inits = inits)
    waic1 <- cmcmc$calculateWAIC()
    expect_identical(waic1, cmcmc$getWAIC()$WAIC)

    m <- nimbleModel(code, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits <- inits)
    cm <- compileNimble(m)
    mcmc <- buildMCMC(m, enableWAIC = TRUE)
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    out2 <- runMCMC(cmcmc, niter = 1000, WAIC = TRUE, inits = inits)
    waic2 <- cmcmc$getWAIC()
    expect_equal(waic2$WAIC, waic1)
    expect_identical(c(waic2$WAIC, waic2$lppd, waic2$pWAIC), c(out2$WAIC$WAIC, out2$WAIC$lppd, out2$WAIC$pWAIC))

    waic2full <- cmcmc$getWAICdetails(returnElements = TRUE)
    expect_identical(waic2full$thin, FALSE)
    expect_identical(waic2full$online, TRUE)
    expect_identical(waic2full$marginal, FALSE)
    expect_identical(waic2full$niterMarginal, 0)
    expect_equal(sum(waic2full$WAIC_elements), waic2$WAIC)
    expect_equal(sum(waic2full$lppd_elements), waic2$lppd)
    expect_equal(sum(waic2full$pWAIC_elements), waic2$pWAIC)
    expect_identical(length(waic2full$WAIC_elements), as.integer(J*I))

    ## Use of nimbleMCMC
    set.seed(1)
    out <- nimbleMCMC(code,data = list(y = y),
                      constants = list(I = I, J = J),
                      inits = inits, niter = 1000, WAIC = TRUE)
    expect_identical(out$WAIC$WAIC, waic2$WAIC)

    ## Multiple chains; check first equivalent to first above
    m <- nimbleModel(code, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits = inits)
    cm <- compileNimble(m)

    mcmc <- buildMCMC(m, enableWAIC = TRUE, controlWAIC = list(online = FALSE),
                  monitors = c('mu0','mu','sigma','tau'))
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    out3 <- runMCMC(cmcmc, niter = 1000, nchains = 3, WAIC = TRUE, inits = inits, perChainWAIC = TRUE)
    waic3 <- cmcmc$calculateWAIC()
    expect_identical(out3$WAIC[3], waic3)
    expect_identical(out3$WAIC[1], waic1)

})


test_that("Conditional grouped WAIC matches conditional 'ungrouped' with a mv node", {
    set.seed(1)
    J <- 5
    I <- 10
    tau <- 1
    sigma <- 1
    
    mu <- rnorm(J, 0 , tau)
    
    y <- matrix(0, J, I)
    for(j in 1:J) 
        y[j, ] <- rnorm(I, mu[j], sigma)

    code <- nimbleCode({
        tau ~ dunif(0, 10)
        sigma ~ dunif(0, 10)
        mu0 ~ dnorm(0, 10)
        for(j in 1:J) {
            mu[j] ~ dnorm(mu0, sd = tau)
            for(i in 1:I)
                y[j, i] ~ dnorm(mu[j], sd = sigma)
        }
    })

    inits <- list(mu0 = 0, tau = 0.5, sigma = 1.5,
                               mu = rnorm(J, 0, 0.5))

    code2 <- nimbleCode({
        tau ~ dunif(0, 10)
        sigma ~ dunif(0, 10)
        mu0 ~ dnorm(0, 10)
        Sigma[1:I, 1:I] <- sigma^2 * Iden[1:I, 1:I]
        for(j in 1:J) {
            mu[j] ~ dnorm(mu0, sd = tau)
            mn[j, 1:I] <- ones[1:I]*mu[j]
            y[j, 1:I] ~ dmnorm(mn[j, 1:I], cov = Sigma[1:I, 1:I])
        }
    })

    inits2 <- c(inits, list(ones = rep(1, I), Iden = diag(I)))
    

    m <- nimbleModel(code2, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits = inits2)
    cm <- compileNimble(m)

    conf <- configureMCMC(m, nodes = NULL, enableWAIC = TRUE, controlWAIC = list(online = FALSE), 
                          monitors = c('mu0','mu','sigma','tau'))
    conf$addSampler('tau','RW')
    conf$addSampler('sigma','RW')
    for(j in 1:J)
        conf$addSampler(paste0('mu[', j, ']'), 'RW')
    conf$addSampler('mu0','conjugate')
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    out1 <- runMCMC(cmcmc, niter = 1000, WAIC = TRUE, inits = inits)
    waic1 <- cmcmc$calculateWAIC()

    dataGroups = list('y[1, ]', 'y[2, ]', 'y[3, ]', 'y[4, ]', 'y[5, ]') 
    m <- nimbleModel(code, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits = inits)
    cm <- compileNimble(m)
    conf <- configureMCMC(m, nodes = NULL, enableWAIC = TRUE, controlWAIC = list(dataGroups = dataGroups))
    conf$addSampler('tau','RW')
    conf$addSampler('sigma','RW')
    for(j in 1:J)
        conf$addSampler(paste0('mu[', j, ']'), 'RW')
    conf$addSampler('mu0','conjugate')
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    out2 <- runMCMC(cmcmc, niter = 1000, inits = inits)
    waic2 <- cmcmc$getWAIC()
    expect_identical(waic2$WAIC, waic1)

})

test_that("Marginal grouped WAIC implementation matches exact based on analytic integral", {
    set.seed(1)
    J <- 5
    I <- 10
    tau <- 1
    sigma <- 1
    
    mu <- rnorm(J, 0 , tau)
    
    y <- matrix(0, J, I)
    for(j in 1:J) 
        y[j, ] <- rnorm(I, mu[j], sigma)

    code <- nimbleCode({
        tau ~ dunif(0, 10)
        sigma ~ dunif(0, 10)
        mu0 ~ dnorm(0, 10)
        for(j in 1:J) {
            mu[j] ~ dnorm(mu0, sd = tau)
            for(i in 1:I)
                y[j, i] ~ dnorm(mu[j], sd = sigma)
        }
    })

    inits <- list(mu0 = 0, tau = 0.5, sigma = 1.5,
                               mu = rnorm(J, 0, 0.5))
    dataGroups = list('y[1, ]', 'y[2, ]', 'y[3, ]', 'y[4, ]', 'y[5, ]') 
    m <- nimbleModel(code, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits = inits)
    cm <- compileNimble(m)
    conf <- configureMCMC(m, enableWAIC = TRUE, 
                          controlWAIC = list(dataGroups = dataGroups,
                                             marginalizeNodes = 'mu',
                                             niterMarginal = 10000))
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    samples <- runMCMC(cmcmc, niter = 1000, inits = inits)
    waic <- cmcmc$getWAIC()

    ## Now compute WAIC using analytic marginal likelihood.
    logPredProbs <- matrix(nrow = nrow(samples), ncol = J)
    logAvgProb <- 0
    pWAIC <- 0
            
    for(i in 1:nrow(samples)) {
        cov <- matrix(samples[i, 'tau']^2, I, I)
        diag(cov) <- diag(cov) + samples[i, 'sigma']^2
        ch <- chol(cov)
        for(j in 1:J) 
            logPredProbs[i, j] <-  dmnorm_chol(y[j, 1:I], samples[i, 'mu0'],
                                cholesky = ch, prec_param = FALSE, log = TRUE)
    }
    for(j in 1:J) {
        maxLogPred <- max(logPredProbs[,j])
        thisDataLogAvgProb <- maxLogPred + log(mean(exp(logPredProbs[,j] - maxLogPred)))
        logAvgProb <- logAvgProb + thisDataLogAvgProb
        pointLogPredVar <- var(logPredProbs[,j])
        pWAIC <- pWAIC + pointLogPredVar
    }
    WAIC <- -2*(logAvgProb - pWAIC)
    expect_lt(abs(waic$lppd - logAvgProb), .002)
    expect_lt(abs(waic$pWAIC - pWAIC), .01)
})

test_that("invalid WAIC configuration is trapped", {
    set.seed(1)
    J <- 5
    I <- 10
    tau <- 1
    sigma <- 1
    
    mu <- rnorm(J, 0 , tau)
    
    y <- matrix(0, J, I)
    for(j in 1:J) 
        y[j, ] <- rnorm(I, mu[j], sigma)

    code <- nimbleCode({
        tau ~ dunif(0, 10)
        sigma ~ dunif(0, 10)
        mu0 ~ dnorm(0, 10)
        for(j in 1:J) {
            mu[j] ~ dnorm(mu0, sd = tau)
            for(i in 1:I)
                y[j, i] ~ dnorm(mu[j], sd = sigma)
        }
    })

    inits <- list(mu0 = 0, tau = 0.5, sigma = 1.5,
                               mu = rnorm(J, 0, 0.5))
    dataGroups = list('y[1, ]', 'y[2, ]', 'y[3, ]', 'y[4, ]', 'y[5, ]') 
    m <- nimbleModel(code, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits = inits)
    conf <- configureMCMC(m, enableWAIC = TRUE, 
                          controlWAIC = list(dataGroups = dataGroups[1:3],
                                             marginalizeNodes = 'mu'))
    expect_warning(mcmc <- buildMCMC(conf), "Potential problem with data grouping")

    conf <- configureMCMC(m, enableWAIC = TRUE, 
                          controlWAIC = list(dataGroups = c(dataGroups, 'tau'),
                                             marginalizeNodes = 'mu'))
    expect_warning(mcmc <- buildMCMC(conf), "Potential problem with data grouping")
    
    conf <- configureMCMC(m, enableWAIC = TRUE, 
                          controlWAIC = list(marginalizeNodes = c('mu0', 'mu')))
    expect_warning(mcmc <- buildMCMC(conf), "Potential problem with nodes to marginalize over")
})


test_that("standalone offline WAIC works", {
    set.seed(1)
    J <- 5
    I <- 10
    tau <- 1
    sigma <- 1

    mu <- rnorm(J, 0 , tau)

    y <- matrix(0, J, I)
    for(j in 1:J) 
        y[j, ] <- rnorm(I, mu[j], sigma)

    code <- nimbleCode({
        tau ~ dunif(0, 10)
        sigma ~ dunif(0, 10)
        mu0 ~ dnorm(0, 10)
        for(j in 1:J) {
            mu[j] ~ dnorm(mu0, sd = tau)
            for(i in 1:I)
                y[j, i] ~ dnorm(mu[j], sd = sigma)
        }
    })
    inits <- list(mu0 = 0, tau = 0.5, sigma = 1.5,
                  mu = rnorm(J, 0, 0.5))

    m <- nimbleModel(code, data = list(y = y),
                     constants = list(I = I, J = J),
                     inits <- inits)
    cm <- compileNimble(m)

    ## Baseline offline WAIC result.
    mcmc <- buildMCMC(m, monitors = c('mu0','mu','sigma','tau'), enableWAIC = TRUE,
            controlWAIC = list(online = FALSE))          
    cmcmc <- compileNimble(mcmc, project = m)
    set.seed(1)
    out <- runMCMC(cmcmc, niter = 1000, inits = inits, WAIC = TRUE)
    waic_gold <- out$WAIC
    

    ## Various uses of post hoc calculateWAIC for offline WAIC
    mcmc <- buildMCMC(m, monitors = c('mu0','mu','sigma','tau'))
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    set.seed(1)
    samples <- runMCMC(cmcmc, niter = 1000, inits = inits)

    out <- calculateWAIC(cmcmc)
    expect_identical(waic_gold, out$WAIC)
    out <- calculateWAIC(cmcmc, cm)
    expect_identical(waic_gold, out$WAIC)
    out <- calculateWAIC(mcmc)
    expect_identical(waic_gold, out$WAIC)

    mcmc <- buildMCMC(m, monitors = c('mu0','mu','sigma','tau'))
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    set.seed(1)
    samples <- runMCMC(cmcmc, niter = 1000, inits = inits)
    out <- calculateWAIC(samples, cm)
    expect_identical(waic_gold, out$WAIC)

    mcmc <- buildMCMC(m, monitors = c('mu0','mu','sigma','tau'))
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    set.seed(1)
    samples <- runMCMC(cmcmc, niter = 1000, inits = inits)
    out <- calculateWAIC(samples, m)
    expect_identical(waic_gold, out$WAIC)
    

    ## Use of calculateWAIC simply as a wrapper to already-enabled offline WAIC
    mcmc <- buildMCMC(m, enableWAIC = TRUE, controlWAIC = list(online = FALSE),
                      monitors = c('mu0','mu','sigma','tau'))
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    set.seed(1)
    out <- runMCMC(cmcmc, niter = 1000, inits = inits)
    waic <- calculateWAIC(cmcmc)  
    expect_identical(waic_gold, waic)
    
    ## Use of calculateWAIC simply as a wrapper to already-enabled online WAIC
    mcmc <- buildMCMC(m, enableWAIC = TRUE,
                      monitors = c('mu0','mu','sigma','tau'))
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    set.seed(1)
    out <- runMCMC(cmcmc, niter = 1000, inits = inits)
    waic <- calculateWAIC(cmcmc)  
    expect_equal(waic_gold, waic$WAIC)

    ## Missing monitors 
    mcmc <- buildMCMC(m, monitors = c('mu0'))
    cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
    out <- runMCMC(cmcmc, niter = 1000, inits = inits)
    expect_error(calculateWAIC(cmcmc), "must be monitored")
})


nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
