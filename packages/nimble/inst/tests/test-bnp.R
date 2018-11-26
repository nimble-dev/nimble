rm(list=ls())
library(nimble)
source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

library(mvtnorm)
#library(MCMCpack)

RwarnLevel <- options('warn')$warn
options(warn = 1)

#options(warn = RwarnLevel)
#nimbleOptions(verbose = nimbleVerboseSetting)
#nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context('Testing of BNP functionality')

test_that("Test that new cluster parameters are correctly updated in CRP sampler", {
  
  set.seed(1)
  
  # Data ~ Poisson(5). Starting values are extremely away fromt their true values. 
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      lambda[i] ~ dgamma(1, 0.01)
      y[i] ~ dpois(lambda[xi[i]])
    }
    alpha ~ dgamma(1, 1)
  })
  n <- 300
  Consts <- list(n = n)
  Inits <- list(xi = rep(1, n), 
                lambda = rep(200, n), 
                alpha = 200)
  Data <- list(y = rpois(n, 5))
  
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('xi', 'alpha', 'lambda'), print=FALSE)  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  
  output <- runMCMC(cMCMC, niter=1000, nburn=0, thin=1 , inits=Inits, setSeed=FALSE)
  
  xiSam <- output[500:1000, 302:601]
  lambdaSam <- output[500:1000, 2:301]
  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  lambdaUnique <- sapply(1:500, function(i) unique(lambdaSam[i, xiSam[i, ]]) )
  cond <- sapply(1:500, function(i) sum(weights[[i]]*lambdaUnique[[i]]))
  
  expect_equal(mean(cond), 5, tol=2*sd(cond), scale=1,
               info = paste0("incorrect update of cluster parameters in Poisson data"))
  
  # test for updating new cluster parameters, with conjugacy

  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], 1)
      mu[i] <- muTilde[xi[i]]
    }
    
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    
    xi[1:n] ~ dCRP(alpha, size = n)
    sd0 ~ dhalfflat()
    alpha ~ dgamma(1, 1)      
    mu0 ~ dflat()
  })
  
  n <- 30
  constants <- list(n = n)
  ## all data plausibly from first cluster except 50th data point
  data <- list(y = c(rnorm(n-1, 0, 1), 50))
  ## muTilde is good for all but last data point. muTilde[2] is bad for the last data point (so that we can see that it changes to a good value, which is what the conjugate sampler for xi should ensure)
  inits <- list(alpha = 1, mu0 = 0, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, rep(-10, n-1)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  ## now check that cmodel$muTilde[2] is near 50
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  
  xiSam <- output[, 34:63]
  muTildeSam <- output[, 3:32]
  cond <- as.numeric(muTildeSam[2])
  
  expect_equal(cond, 50, tolerance=2*1, scale=1, 
               info = paste0("incorrect update of cluster parameters in mixture of normals 1 data"))
  
  # test for updating new cluster parameters, without conjugacy
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], 1)
      mu[i] <- muTilde[xi[i]]
    }
    
    for(i in 1:n) {  
      muTilde[i] ~ dt(mu0, df = 40, sigma = sd0) # force non-conjugate
    }
    
    xi[1:n] ~ dCRP(alpha, size = n)
    sd0 ~ dhalfflat()
    alpha ~ dgamma(1, 1)      
    mu0 ~ dflat()
  })
  
  n <- 30
  constants <- list(n = n)
  ## all data plausibly from first cluster except 50th data point
  data <- list(y = c(rnorm(n-1, 0, 1), 50))
  ## muTilde is good for all but last data point. muTilde[2] is bad for the last data point (so that we can see that it changes to a good value, which is what the conjugate sampler for xi should ensure)
  inits <- list(alpha = 1, mu0 = 0, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, rep(-10, n-1)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  ## now check that 50th obs remains in initial bad cluster because muTilde[2] is even worse
  output <- runMCMC(cmcmc, niter=1, nburn=0, thin=1 , inits=inits, setSeed=FALSE)

  clusterID <- output[1, paste0('xi[', n, ']')]
  attributes(clusterID) <- NULL
  expect_identical(clusterID, 1,
                   info = 'non-conjugate incorrectly chose bad new cluster')

  ## check that 50th obs moves to better cluster now that muTilde[2] is decent
  cmodel$muTilde[2] <- 40
  cmodel$calculate()
  output <- runMCMC(cmcmc, niter=1, nburn=0, thin=1, setSeed=FALSE)
  clusterID <- output[1, paste0('xi[', n, ']')]
  attributes(clusterID) <- NULL
  expect_identical(clusterID, 2,
                   info = 'non-conjugate incorrectly did not choose good new cluster')
  value <- output[1, "muTilde[2]"]
  attributes(value) <- NULL
  expect_equal(value, 40, tolerance=3,
               info = 'non-conjugate has strange cluster parameter for new cluster')
   
  # We start with only one active component and the data is a mixture of 3 normal ditributions
  set.seed(1)
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i] ~ dnorm(0, var = s2[i]/lambda)
      s2[i] ~ dinvgamma(2, 1)
      y[i] ~ dnorm(mu[xi[i]],  var = s2[xi[i]])
    }
    lambda ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)
  })
  n <- 300
  Consts <- list(n = n)
  Inits <- list(xi = rep(1, n), 
                mu = rep(-20, n), 
                s2 = rep(0.1, n),
                alpha = 200,
                lambda = 1)
  thetas <- c(rep(-5, 100), rep(5, 100), rep(0, 100))
  Data <- list(y = rnorm(n, thetas, 1))
  
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 's2', 'alpha', 'lambda'), print=FALSE)  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  
  output <- runMCMC(cMCMC, niter=7000, nburn=2000, thin=10 , inits=Inits, setSeed=FALSE)
  expect_output(outputG <- getSamplesDPmeasure(cMCMC))
  
  Tr <- outputG$trunc
  samplesG <- outputG$samples
  grid <- seq(-15, 15, len=100)
  samF <- matrix(0, ncol=length(grid), nrow=500)
  for(i in 1:500) {
    samF[i, ] <- sapply(grid, function(x)sum(samplesG[i, 1:Tr] * dnorm(x, samplesG[i, (2*Tr+1):(3*Tr)], sqrt(samplesG[i, (Tr+1):(2*Tr)]))))
  }
  
  # distance between the estimation and truth at each point of the grid 
  f0 <- sapply(grid, function(x) 0.3*dnorm(x, -5, 1) + 0.3*dnorm(x, 0, 1) + 0.3*dnorm(x, 5, 1))
  for(i in 1:100) {
    expect_equal(mean(samF[,i]), f0[i], tol=2*sd(samF[,i]), scale=1,
                 info = paste0("incorrect update of cluster parameters in mixture of normals data. Grid_i = ", i))
  }
})


test_that("testing bivariate normal mixture models with CRP", {
  
  set.seed(1)
  
  # bivariate normal kernel with unknown mean
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i, 1:2] ~ dmnorm(mu0[1:2], cov = S0[1:2, 1:2])
      y[i, 1:2] ~ dmnorm(mu[xi[i], 1:2],  cov = Sigma0[1:2, 1:2])
    }
    alpha ~ dgamma(1, 1)
  })
  n <- 100
  Consts <- list(n = n, Sigma0 = diag(1, 2), S0 = diag(100, 2), mu0=c(0,0))
  Inits <- list(xi = sample(1:100, size=n, replace = TRUE), 
                mu = matrix(-10, ncol=2, nrow=n), 
                alpha = 1)
  Data <- list(y = rmvnorm(n, c(10, 20), diag(1, 2)))
  
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 'alpha'), print=FALSE)  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  
  output <- runMCMC(cMCMC, niter=100, nburnin=50, thin=1 , inits=Inits,setSeed=FALSE, progressBar=TRUE)
  outputG <- getSamplesDPmeasure(cMCMC)
  
  xiSam <- output[, 202:301]
  muSam1 <- output[, 2:101]
  muSam2 <- output[, 102:201]

  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  muSam1Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam1[i, xiSam[i, ]]) )
  muSam2Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam2[i, xiSam[i, ]]) )
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam1Unique[[i]]))
  expect_equal(mean(cond), 10, tol=2*sd(cond), scale=1,
               info = paste0("incorrect update of cluster parameters in bivariate normal data with unknown mean. First component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam2Unique[[i]]))
  expect_equal(mean(cond), 20, tol=2*sd(cond), scale=1,
               info = paste0("incorrect update of cluster parameters in bivariate normal data with unknown mean. Second component"))
  
  
  # bivariate normal kernel with unknown mean and unknown variance
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i, 1:2] ~ dmnorm(mu0[1:2], cov = S0[1:2, 1:2])
      Sigma[1:2, 1:2, i] ~ dinvwish(S = R0[1:2, 1:2], df = 4)
      SigmaAux[1:2, 1:2, i] <- Sigma[1:2, 1:2, xi[i]] / lambda  
      y[i, 1:2] ~ dmnorm(mu[xi[i], 1:2],  cov = SigmaAux[1:2, 1:2, i] )
    }
    alpha ~ dgamma(1, 1)
    lambda ~ dgamma(1, 1)
  })
  n <- 100
  Consts <- list(n = n, S0 = diag(1, 2), mu0=c(0,0), R0 = diag(1, 2))
  Sigma <- array(0, c(2,2,n))
  for(i in 1:n)
    Sigma[, , i] <- matrix(c(1, 0, 0, 1), 2, 2)
  Inits <- list(xi = sample(1:n, size=n, replace = FALSE), 
                mu = matrix(-10, ncol=2, nrow=n),
                Sigma=Sigma,
                alpha = 1,
                lambda = 1)
  Data <- list(y = rmvnorm(n, c(10, 20), diag(1, 2)))
  
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m, showCompilerOutput = FALSE)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 'Sigma', 'alpha', 'lambda'), print=FALSE)  #, 'lambda'
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  
  #output <- runMCMC(cMCMC, niter=5000, nburn=3000, thin=2 ,
  #                  inits=Inits, setSeed=FALSE)
  output <- runMCMC(cMCMC, niter=10, nburnin=1, thin=1 ,inits=Inits, setSeed=FALSE)
  
  outputG <- getSamplesDPmeasure(cMCMC)
  
  xiSam <- output[, 603:702]
  sigma11Sam <- output[, seq(1, 400, by=4)]
  muSam1 <- output[, 403:502]
  muSam2 <- output[, 503:602]
  
  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  muSam1Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam1[i, xiSam[i, ]]) )
  muSam2Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam2[i, xiSam[i, ]]) )
  sigma11Unique <- sapply(1:nrow(xiSam), function(i) unique(sigma11Sam[i, xiSam[i, ]]) )
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam1Unique[[i]]))
  expect_equal(mean(cond), 10, tol=2*sd(cond), scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. First component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam2Unique[[i]]))
  expect_equal(mean(cond), 20, tol=2*sd(cond), scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. Second component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*sigma11Unique[[i]]))
  expect_equal(mean(cond), 1, tol=2*sd(cond), scale=1,
               info = paste0("incorrect update of covariance matrix parameters in bivariate normal data with unknown mean and variance. [1, 1] component"))
  
  
})



test_that("testing dnorm_invgamma distribution in BNP functionality", {
  set.seed(1)
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      params[i, 1:2] ~ dnorm_invgamma(mu0, lambda, a, b)
      y[i] ~ dnorm(params[xi[i], 1],  var = params[xi[i], 2])
    }
    mu0 ~ dnorm(0, sd=20)
    lambda ~ dgamma(1, 1)
    a ~ dgamma(2, 1)
    b ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)
  })
  n <- 100
  Consts <- list(n = n)
  Inits <- list(xi = rep(1, n), 
                params = cbind(rep(0, n), rep(1, n)),
                mu0=0, a=2, b=1,
                alpha = 1,
                lambda = 1)
  Data = list( y = c(rnorm(n/2, 0,1), rnorm(n/2, 30,1)))
  
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('xi','params', 'alpha', 'lambda', 'mu0', 'a', 'b'), print=FALSE)  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  
  output <- runMCMC(cMCMC, niter=10000, nburn=5000, thin=5 , inits=Inits, setSeed=TRUE)
  
  muSam <- output[, 6:105]
  xiSam <- output[, 206:305]
  
  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  muSamUnique <- sapply(1:nrow(xiSam), function(i) unique(muSam[i, xiSam[i, ]]) )

  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSamUnique[[i]]))

  expect_equal(mean(cond), 0.5*0 + 0.5*30, tol=2*sd(cond), scale=1,
               info = paste0("incorrect results when dnorm_invgamma prior is used in normal model"))
  
  # for checking conjugacy detection:
  #class(mMCMC$samplerFunctions[[13]]$helperFunctions$contentsList[[1]])[1]
})

    
test_that("sampleDPmeasure: testing that required variables in MCMC modelValues are monitored", {
  set.seed(1)
  
  ## membership variable not being monitored
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 ~ dgamma(1, 1)
    for(i in 1:6){
      mu[i] ~ dnorm(0, 1)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6, conc0 = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The node having the dCRP distribution')
 
  mConf <- configureMCMC(m, monitors = c('xi', 'conc0', 'mu'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_output(output <- getSamplesDPmeasure(mMCMC))
  
  # cluster variable not being monitored
  set.seed(1)
  code <- nimbleCode({
    xi[1:6] ~ dCRP(1, 6)
    mu0 ~ dnorm(0, 1)
    s20 ~ dgamma(1, 1)
    for(i in 1:6){
      mu[i] ~ dnorm(mu0, s20)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6, mu0 = 0, s20 = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The node\\(s\\) representing the cluster variables') 
  
  mConf <- configureMCMC(m, monitors = c('mu', 'xi'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes')

  mConf <- configureMCMC(m, monitors = c('mu', 'xi', 'mu0', 's20'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_output(output <- getSamplesDPmeasure(mMCMC))
  
  ## concentration parameter not being monitored:
  set.seed(1)
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 ~ dgamma(a, rate=b)
    a ~ dgamma(1, rate=1)
    b ~ dgamma(1, rate=0.1)
    for(i in 1:6){
      mu[i] ~ dnorm(0, 1)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6, conc0 = 1, a = 1, b = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m, monitors = c('xi', 'mu'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'conc0'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_output(outputG <- getSamplesDPmeasure(mMCMC))

  ## concentration parameter deterministic parent not being monitored:
  set.seed(1)
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 <- a + b
    a ~ dgamma(1, rate=1)
    b <- d + 1
    d ~ dgamma(1, 1)
    for(i in 1:6){
      mu[i] ~ dnorm(0, 1)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6, conc0 = 1, a = 1, d = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m, monitors = c('xi', 'mu'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.')

  ## note that if 'b' not 'd' is initialized then this will fail because model initialization overwrites
  ## 'b' with NA based on NA in 'd' and then sampleDPmeasure setup thinks it can't calculate conc0
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'b'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_output(outputG <- getSamplesDPmeasure(mMCMC))

  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'd'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_output(outputG <- getSamplesDPmeasure(mMCMC))  
})

test_that("sampleDPmeasure: checking for uncompiled modelValues input", {
  library(nimble)
  set.seed(1)
  code <- nimbleCode({
    xi[1:6] ~ dCRP(1, 6)
    for(i in 1:6){
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = rnorm(6))
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m) 
  mConf <- configureMCMC(m, monitors = c('xi', 'mu'))
  mMCMC <- buildMCMC(mConf)
  CmMCMC <- compileNimble(mMCMC, project=m)
  output <- runMCMC(CmMCMC, 10)
  getSamplesDPmeasure(CmMCMC)
  expect_output(getSamplesDPmeasure(CmMCMC))
})




test_that("stick_breaking nimble function calculation and use is correct", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  ltruth <- log(truth)
  
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation"))
  
  expect_equal(stick_breaking(x, log=TRUE),
               ltruth,
               info = paste0("incorrect stick_breaking nimble function log calculation"))
  
  cSB <- compileNimble(stick_breaking)
  
  expect_equal(cSB(x, log=FALSE),
               truth,
               info = paste0("incorrect compiled stick_breaking nimble function calculation"))
  
  expect_equal(cSB(x, log=TRUE),
               ltruth,
               info = paste0("incorrect compiled stick_breaking nimble function log calculation"))
  
  x <- c(0.1, 0.4, -0.1, 0.3)
  expect_output(aux <- stick_breaking(x, log=FALSE), "values in 'z' have to be in", 
                info = "stick_breaking not warning of negative component")
  expect_equal(aux, rep(NaN, length(x)+1),
               info = "stick_breaking not correctly handling negative component")
  
  x <- c(0.1, 5, 0.4, 0.3)
  expect_output(aux <- stick_breaking(x, log=FALSE), "values in 'z' have to be in")
  expect_equal(aux, rep(NaN, length(x)+1),
               info = "stick_breaking incorrectly handling larger than 1 component")
  
  x <- c(0.1, 0.2, 0, 0.3, 0.8)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation with one 0 component"))
  
  x <- c(0.1, 0.2, 1, 0.3, 0.8)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation with one 1 component"))
})


test_that("Stick breaking model calculation is correct", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  
  SB_code <- nimbleCode({
    for(i in 1:5) z[i] ~ dbeta(1, 1)
    w[1:6] <- stick_breaking(z[1:5])
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  SB_model <- nimbleModel(SB_code, data=Inits)
  
  SB_model$z <- x
  SB_model$calculate()
  
  expect_equal(c(SB_model$w), truth,
               info = paste0("incorrect stick breaking weights in model"))
  
  c_SB_model <- compileNimble(SB_model)
  
  c_SB_model$z <- x
  c_SB_model$calculate()
  c_SB_model$w
  
  expect_equal(c(c_SB_model$w), truth,
               info = paste0("incorrect stick breaking weights in compiled model"))
  
})



## test simple models

model <- function() {
  for(j in 1:5) 
    z[j] ~ dbeta(1, 1)
  w[1:6] <- stick_breaking(z[1:5])
  for(i in 1:10){
    xi[i] ~ dcat(w[1:6])
  }
}

Inits <- list(z = rep(0.5,5))
Data <- list(xi = 1:10)

testBUGSmodel(example = 'test1', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)

model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    y[i] ~ dnorm( thetatilde[xi[i]], var=1)
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=1:10)
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test2', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)

model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    theta[i] <- thetatilde[xi[i]]
    y[i] ~ dnorm( theta[i], var=1)
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=1:10)
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test3', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)

model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    theta[i] <- thetatilde[xi[i]]
    y[i] ~ dnorm( theta[i], var=s2tilde[xi[i]])
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
    s2tilde[i] ~ dinvgamma(1, 1)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=1:10, s2tilde=rep(1,10))
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test4', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)



test_that("Testing conjugacy detection with bnp stick breaking models", { 
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dbeta(1, 1)
      }
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[6]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy")
  
  ## replace beta by uniform distribution: here no conjugacy is detected.
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dunif(0,1)
      }
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4)))
  conf <- configureMCMC(m)
  expect_failure(expect_match(conf$getSamplers()[[6]]$name,
                              "conjugate_dbeta_dcat",
                              info = "failed to detect categorical-beta conjugacy"))
  
  
  ## concentration parameter added in beta distribution
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dbeta(1, conc)
      }
      conc ~ dgamma(1,1)
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4), conc=1))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[7]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy")
  
})

test_that("Testing BNP model using stick breaking representation", { 
  
  Code=nimbleCode(
    {
      for(i in 1:Trunc) {
        thetatilde[i] ~ dnorm(mean=0, var=40) 
        s2tilde[i] ~ dinvgamma(shape=1, scale=0.5) 
      }
      for(i in 1:(Trunc-1)) {
        z[i] ~ dbeta(1, 1)
      }
      w[1:Trunc] <- stick_breaking(z[1:(Trunc-1)])
      
      for(i in 1:N) {
        xi[i] ~ dcat(w[1:Trunc])
        theta[i] <- thetatilde[xi[i]]
        s2[i] <- s2tilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2[i])
      }
    }
  )
  
  Consts <- list(N=50, Trunc=25)
  set.seed(1)
  Inits <- list(thetatilde = rnorm(Consts$Trunc, 0, sqrt(40)),
                s2tilde = rinvgamma(Consts$Trunc, shape=1, scale=0.5),
                z = rbeta(Consts$Trunc-1, 1, 1),
                xi = sample(1:10, size=Consts$N, replace=TRUE))
  Data = list(y = c(rnorm(Consts$N/2,5,sqrt(4)), rnorm(Consts$N/2,-5,sqrt(4))))
  
  model = nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
  cmodel = compileNimble(model)
  
  modelConf = configureMCMC(model,, thin=100)
  
  expect_match(modelConf$getSamplers()[[51]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy in BNP model")
  
  modelMCMC = buildMCMC(modelConf)
  CmodelMCMC = compileNimble(modelMCMC, project=model, resetFunctions=TRUE)
  
  CmodelMCMC$run(10000)
  
  ## results from the algorithm:
  samples = as.matrix(CmodelMCMC$mvSamples)
  s2Sam = samples[, 1:25]
  thetaSam = samples[, 26:50]
  zSam = samples[, 51:74]
  Tr = 25
  Wpost = t(apply(zSam, 1, function(x)c(x[1], x[2:(Tr-1)]*cumprod(1-x[1:(Tr-2)]), cumprod(1-x[1:(Tr-1)])[N=Tr-1])))
  
  ngrid = 302
  grid = seq(-10, 25,len=ngrid)
  
  # posterior samples of the density
  nsave = 100
  predSB = matrix(0, ncol=ngrid, nrow=nsave)
  for(i in 1:nsave) {
    predSB[i, ] = sapply(1:ngrid, function(j)sum(Wpost[i, ]*dnorm(grid[j], thetaSam[i,],sqrt(s2Sam[i,]))))
  }
  
  ## hist(Data$y, freq=FALSE, xlim=c(min(grid), max(grid)))
  ## points(grid, apply(predSB, 2, mean), col="blue", type="l", lwd=2)
  
  f0 <- function(x) 0.5*dnorm(x,5,sqrt(4)) + 0.5*dnorm(x,-5,sqrt(4))
  fhat <- apply(predSB, 2, mean)
  f0grid <- sapply(grid, f0)
  
  L1dist <- mean(abs(f0grid - fhat))
  
  expect_equal(L1dist, 0.01, tol=0.01,
               info = "wrong estimation of density in DPM of normal distrbutions")
  
})

## Testing CRP distribution

test_that("Check error given when model has no cluster variables", {
  
  set.seed(1)
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 ~ dgamma(1, 1)
    for(i in 1:6){
      y[i] ~ dnorm(xi[i], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2),  conc0 = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  
  expect_error(buildMCMC(mConf) ,
               'sampler_CRP:  The model should have at least one cluster variable.')
  
})

test_that("dCRP nimble function calculates density correctly",{
  
  x <- c(1,1,2,1,1,2)
  conc <- 1
  
  truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
    (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1))
  ltruth <- log(truth)
  
  expect_equal(dCRP(x, conc, size=length(x), log=FALSE),
               truth,
               info = paste0("incorrect dCRP nimble function calculation"))
  
  expect_equal(dCRP(x, conc, size=length(x), log=TRUE),
               ltruth,
               info = paste0("incorrect dCRP nimble function calculation in log scale"))
  
  cdCRP <- compileNimble(dCRP)
  
  expect_equal(cdCRP(x, conc, size=length(x)), (truth), 
               info = paste0("incorrect dCRP value in compiled nimble function"))
  
  expect_equal(cdCRP(x, conc, size=length(x), log=TRUE), (ltruth), 
               info = paste0("incorrect dCRP value in compiled nimble function  in log scale"))
  
  expect_equal(dCRP(x, conc=-1, size=length(x), log=FALSE),
               NaN,
               info = paste0("incorrect parameters space allowed"))
  
  expect_error(dCRP(x, conc=1, size=3, log=FALSE), "length of 'x' has to be equal to 'size'")
  
  expect_error(dCRP(x, conc=1, size=10, log=FALSE), "length of 'x' has to be equal to 'size'")
  
})

test_that("CRP model calculation and dimensions are correct:", {
  
  x <- c(1,1,2,1,1,2)
  conc <- 1
  
  truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
    (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1))
  ltruth <- log(truth)
  
  CRP_code <- nimbleCode({
    x[1:6] ~ dCRP(conc, size=6)
  })
  
  Consts <- list(conc = 1)
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model <- nimbleModel(CRP_code, data=Inits, constants=Consts)
  
  CRP_model$x <- x
  expect_equal(exp(CRP_model$calculate()), truth,
               info = paste0("incorrect likelihood value for dCRP"))
  
  c_CRP_model <- compileNimble(CRP_model)
  c_CRP_model$x
  expect_equal(exp(c_CRP_model$calculate()), truth,
               info = paste0("incorrect likelihood value for compiled dCRP"))
  
  
  # different length of x and size:
  CRP_code2 <- nimbleCode({
    x[1:6] ~ dCRP(1, size=10)
  })
  
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model2 <- nimbleModel(CRP_code2, data=Inits)
  expect_error(CRP_model2$calculate(), "length of 'x' has to be equal to 'size'")
  
  # different length of x and size:
  CRP_code3 <- nimbleCode({
    x[1:6] ~ dCRP(1, size=3)
  })
  
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model3 <- nimbleModel(CRP_code3, data=Inits)
  expect_error(CRP_model3$calculate(), "length of 'x' has to be equal to 'size'")
  
  
})


test_that("random sampling from CRP in model with additional levels", {
  
  conc <- 1
  set.seed(0)
  size <- 6
  r_samps <- t(replicate(10000, rCRP(n = 1, conc, size = size)))
  # K is the number of unique components in x of length 6
  true_EK <- sum(conc/(conc+1:size-1))
  
  expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
               tol = 0.01,
               info = "Difference in expected mean of K exceeds tolerance")
  
  # sampling from the model:
  set.seed(1)
  CRP_code <- nimbleCode({
    x[1:6] ~ dCRP(conc=1, size=6)
    for(i in 1:6){
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[x[i]], 1)
    }
  })
  Inits <- list(x = c(1,1,2,1,1,2), mu = 1:6)
  Data <- list( y =  rnorm(6))
  CRP_model <- nimbleModel(CRP_code, data=Data, inits=Inits)
  c_CRP_model <- compileNimble(CRP_model)
  
  simul_samp <- function(model) {
    model$simulate()
    return(model$x)
  }
  simul_samps <- t(replicate(10000, simul_samp(c_CRP_model)))
  
  expect_equal(mean(apply(simul_samps, 1, function(x)length(unique(x)))), true_EK, 
               tol = 0.01,
               info = "Difference in expected mean of K, from compiled model, exceeds tolerance")
  
})


test_that("Testing conjugacy detection with models using CRP", { 
  
  ## dnorm_dnorm
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm with truncation
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    for(i in 1:2)
        mu[i] ~ dnorm(0,1)
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[3]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")

  ## dnorm_dnorm one more level of hierarchy
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(beta,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
    beta ~ dnorm(0,1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4), beta =1))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[1]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  
  ## dnorm_dnorm and deterministic nodes
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dnorm(mui[i], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm and deterministic nodes and truncation
  code = nimbleCode({
    for(i in 1:4) {
      mui[i] <- mu[xi[i]]
      y[i] ~ dnorm(mui[i], sd = 1)
    }
    for(i in 1:2)
        mu[i] ~ dnorm(0,1)
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[3]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  
  ## dnorm_dpois
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dpois(10)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rpois(4, 10)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
  ## dgamma_dpois
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dpois(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rpois(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  
  ## dgamma_dexp
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dexp(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dexp")
  
  
  ## dgamma_dgamma
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dgamma(4, mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rgamma(4, 4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dgamma")
  
  ## dgamma_dnorm
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dnorm(4, tau = mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4, 4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dnorm")
  
  ## dgamma_dweib
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dweib(shape=4, lambda = mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rweibull(4, 4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dweib")
  
  ## dgamma_dinvgamma
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1, rate=1)
      y[i] ~ dinvgamma(shape=4, scale = mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rinvgamma(4, 4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dinvgamma")
  
  
  ## dbeta_dbern
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dbeta(1,1)
      y[i] ~ dbern(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rbinom(4, size=1, prob=0.5)),
                  inits = list(xi = rep(1,4), mu=rbeta(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dbern")
  
  
  ## dbeta_dbinom
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dbeta(1,1)
      y[i] ~ dbinom(size=10, prob=mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rbinom(4, size=10, prob=0.5)),
                  inits = list(xi = rep(1,4), mu=rbeta(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dbin")
  
  ## dbeta_dnegbin
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dbeta(1,1)
      y[i] ~ dnegbin(size=10, prob=mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnbinom(4, size=10, prob=0.5)),
                  inits = list(xi = rep(1,4), mu=rbeta(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dnegbin")
  
  
  ## ddirch_dmulti
  code=nimbleCode(
    {
      for(i in 1:4){
        p[i,1:3] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[xi[i],1:3], size=3)
      }
      xi[1:4] ~ dCRP(conc=1, size=4)
    }
  )
  set.seed(1)
  p0 <- matrix(0, ncol=3, nrow=4)
  y0 <- matrix(0, ncol=3, nrow=4)
  for(i in 1:4){
    p0[i,]=rdirch(1, c(1, 1, 1))
    y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
  }
  m = nimbleModel(code, 
                  data = list(y = y0),
                  inits = list(xi = rep(1,4), p=p0), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m)
  mcmc <- buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[1]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  ## dnorm_dnorm_dinvgamma
  code = nimbleCode({
    for(i in 1:4) {
      s2[i] ~ dinvgamma(1, 1)
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], var = s2[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4), s2=rinvgamma(4, 1,1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[9]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## dgamma_dexp and deterministic nodes
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dexp(mui[i])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dexp")
  
  
  ## dgamma_dexp, deterministic nodes, and conjugacy is broken
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dexp(mui[i]+3)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
  ## dgamma_dexp, deterministic nodes, and conjugacy is broken
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dexp(3*mui[i])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")

  ## non-exchangeable prior for tilde nodes
  code = nimbleCode({
    for(i in 1:4){
        mu[i] <- muTilde[xi[i]]
        y[i] ~ dnorm(mu[i], sd = 1)
        muTilde[i] ~ dnorm(mu0[i], sd = s0)
        mu0[i] ~ dnorm(0,1)
    }
    xi[1:4] ~ dCRP(1, 4)
    s0 ~ dhalfflat()
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
})

## simple tests of models

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dnorm(0,1)
    y[i] ~ dnorm(mu[xi[i]], sd = 1)
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rnorm(4))
data = list(y = rnorm(4))

testBUGSmodel(example = 'test1', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dpois(mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rgamma(4, 1, 1))
data = list(y = rpois(4, 4))

testBUGSmodel(example = 'test2', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dexp(mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rgamma(4, 1, 1))
data = list(y = rexp(4, 4))

testBUGSmodel(example = 'test3', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dgamma(4, mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rgamma(4, 1, 1))
data = list(y = rgamma(4, 4, 4))

testBUGSmodel(example = 'test4', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dbeta(1,1)
    y[i] ~ dbern(mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rbeta(4, 1, 1))
data = list(y = rbinom(4, size=1, prob=0.5))

testBUGSmodel(example = 'test5', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4){
    p[i,1:3] ~ ddirch(alpha=alpha0[1:3])
    y[i,1:3] ~ dmulti(prob=p[xi[i],1:3], size=3)
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
set.seed(1)
p0 <- matrix(0, ncol=3, nrow=4)
y0 <- matrix(0, ncol=3, nrow=4)
for(i in 1:4){
  p0[i,]=rdirch(1, c(1, 1, 1))
  y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
}
inits = list(xi = 1:4, p=p0)
data = list(y = y0)
alpha0 = c(1,1,1)

testBUGSmodel(example = 'test6', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    s2[i] ~ dinvgamma(1, 1)
    mu[i] ~ dnorm(0,1)
    y[i] ~ dnorm(mu[xi[i]], var = s2[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rnorm(4), s2=rinvgamma(4, 1,1))
data = list(y = rnorm(4))

testBUGSmodel(example = 'test7', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1, rate=1)
    y[i] ~ dinvgamma(shape=4, scale = mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
data = list(y = rinvgamma(4, 4, 4))

testBUGSmodel(example = 'test8', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dweib(shape=4, lambda = mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rgamma(4, 1, 1))
data = list(y = rweibull(4, 4, 4))

testBUGSmodel(example = 'test9', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dnorm(4, tau = mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rgamma(4, 1, 1))
data = list(y = rnorm(4, 4, 4))

testBUGSmodel(example = 'test10', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dbeta(1,1)
    y[i] ~ dbinom(size=10, prob=mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rbeta(4, 1, 1))
data = list(y = rbinom(4, size=10, prob=0.5))

testBUGSmodel(example = 'test11', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dbeta(1,1)
    y[i] ~ dnegbin(size=10, prob=mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rbeta(4, 1, 1))
data = list(y = rnbinom(4, size=10, prob=0.5))

testBUGSmodel(example = 'test11', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


test_that("sampleDPmeasure can be used for more complicated models", {
  
  ## no deterministic node, conc param is fixed
  set.seed(1)
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=a0)
      y[i] ~ dpois(lambdaTilde[xi[i]])
    }
    xi[1:10] ~ dCRP(conc = 1, size=10)
    a0 ~ dgamma(1, 1)
  })
  Inits <- list( xi = 1:10, 
                 lambdaTilde = rgamma(10, shape=1, rate=0.1), a0=1)
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  
  monitors <- c('lambdaTilde','xi', 'a0')
  mConf <- configureMCMC(m, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC, niter = 1000)
  
  # old use of getSamplesDPmeasure:
  #rdens = getSamplesDPmeasure(cMCMC)
  #cdens = compileNimble(rdens$samples, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC))
  expect_false(any(is.na(samplesG$samples)))
  
  ## no deterministic node, random conc param
  set.seed(1)
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=a0)
      y[i] ~ dpois(lambdaTilde[xi[i]])
    }
    xi[1:10] ~ dCRP(conc0, size=10)
    conc0 ~ dgamma(1,1)
    a0 ~ dgamma(1, 1)
  })
  Inits <- list( xi = 1:10, conc0=1,
                 lambdaTilde = rgamma(10, shape=1, rate=0.1), a0=1)
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  
  monitors <- c('lambdaTilde','xi', 'conc0')
  mConf <- configureMCMC(m, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
  
  ## with deterministic node, conc param is fixed
  set.seed(1)
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=a0)
      lambda[i] <- lambdaTilde[xi[i]]
      y[i] ~ dpois(lambda[i])
    }
    xi[1:10] ~ dCRP(conc = 1, size=10)
    a0 ~ dgamma(1, 1)
  })
  Inits <- list( xi = 1:10, 
                 lambdaTilde = rgamma(10, shape=1, rate=0.1), a0=1)
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  
  monitors <- c('lambdaTilde','xi', 'a0')
  mConf <- configureMCMC(m, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
  
  ## with deterministic node, random conc param
  set.seed(1)
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=a0)
      lambda[i] <- lambdaTilde[xi[i]]
      y[i] ~ dpois(lambda[i])
    }
    xi[1:10] ~ dCRP(conc0, size=10)
    conc0 ~ dgamma(1,1)
    a0 ~ dgamma(1, 1)
  })
  Inits <- list( xi = 1:10, conc0=1,
                 lambdaTilde = rgamma(10, shape=1, rate=0.1), a0=1)
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  
  monitors <- c('lambdaTilde','xi', 'conc0', 'a0')
  mConf <- configureMCMC(m, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
  
  # two cluster parameters, one deterministic parameter, fixed conc
  set.seed(1)
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=40) 
        s2tilde[i] ~ dinvgamma(1, scale=0.5)
      }
      xi[1:10] ~ dCRP( 1 , size=10)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), 
             thetatilde=rnorm(10, 0, sqrt(40)),
             s2tilde = rinvgamma(10, 1, scale=0.5))
  Data=list(y=c(rnorm(5,-5,sqrt(5)), rnorm(5,5,sqrt(4))))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  
  monitors <- c('thetatilde', 's2tilde', 'xi')
  mConf <- configureMCMC(m,, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
  
  # two cluster parameters, one deterministic parameter, random conc
  set.seed(1)
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(1, scale=0.5)
      }
      xi[1:10] ~ dCRP( conc0 , size=10)
      conc0 ~ dgamma(1, 1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    lambda ~ dgamma(1, 1)
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), conc0 =1 ,
             thetatilde=rnorm(10, 0, sqrt(40)),
             s2tilde = rinvgamma(10, 1, scale=0.5),
             lambda=1)
  Data=list(y=c(rnorm(5,-5,sqrt(5)), rnorm(5,5,sqrt(4))))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  monitors <- c('thetatilde', 's2tilde', 'xi', 'conc0', 'lambda')
  mConf <- configureMCMC(m,, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
  
  ## two cluster parameters, two deterministc parameter, random conc with random parameters
  set.seed(1)
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(2, scale=0.01)
      }
      xi[1:10] ~ dCRP( conc0 , size=10)
      conc0 ~ dgamma(a,1)
      a ~ dgamma(1, 1)
      lambda ~ dgamma(1, 1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        s2[i] <- s2tilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2[i])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), conc0 = 1 , a=1,
             thetatilde=rnorm(10, 0, sqrt(40)),
             s2tilde = rinvgamma(10, 1, scale=0.01), 
             lambda=1)
  Data=list(y=c(rnorm(5,-5,sqrt(5)), rnorm(5,5,sqrt(4))))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  monitors <- c('thetatilde', 's2tilde', 'xi', 'conc0', 'a', 'lambda')
  mConf <- configureMCMC(m,, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #if(FALSE) expect_output(cdens$run(), "Approximating the random measure")  ## very slow; need to check this
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
  
  # two cluster parameters, two deterministc parameter, 
  # random conc dfined by two random parameters
  set.seed(1)
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var = s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(1, scale=0.5)
      }
      xi[1:10] ~ dCRP( conc0 + conc1 , size=10)
      conc0 ~ dgamma(1,1)
      conc1 ~ dgamma(1,1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        s2[i] <- s2tilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2[i])
      }
      lambda ~ dgamma(1, 1)
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), conc0 =1 , conc1=1,
             thetatilde=rnorm(10, 0, sqrt(40)),
             s2tilde = rinvgamma(10, 1, scale=0.5),
             lambda = 1)
  Data=list(y=c(rnorm(5,-5,sqrt(5)), rnorm(5,5,sqrt(4))))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  monitors <- c('thetatilde', 's2tilde', 'xi', 'conc0', 'conc1', 'lambda')
  mConf <- configureMCMC(m,, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
  
  # two cluster parameters, two deterministic parameter, 
  # random conc defined by two random parameters
  set.seed(1)
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var = s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(1, scale=0.5)
      }
      xi[1:10] ~ dCRP( conc0 + conc1 , size=10)
      conc1 ~ dgamma(1,1)
      conc0 ~ dgamma(a,b)
      b ~ dgamma(1,1)
      a ~ dgamma(1,1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        s2[i] <- s2tilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2[i])
      }
      lambda ~ dgamma(1, 1)
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), conc0 =1 , conc1=1,
             a = 1, b = 1,
             thetatilde=rnorm(10, 0, sqrt(40)),
             s2tilde = rinvgamma(10, 1, scale=0.5),
             lambda = 1 )
  Data=list(y=c(rnorm(5,-5,sqrt(5)), rnorm(5,5,sqrt(4))))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  
  monitors <- c('thetatilde', 's2tilde', 'xi', 'conc0', 'conc1', 'lambda')
  mConf <- configureMCMC(m,, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, 1000)
  
  #rdens = sampleDPmeasure(m, mMCMC$mvSamples)
  #cdens = compileNimble(rdens, project = m)
  #expect_output(cdens$run(), "Approximating the random measure")
  expect_output(samplesG <- getSamplesDPmeasure(cMCMC)$samples)
  expect_false(any(is.na(samplesG)))
})


## testing sampler assigment for conc parameter

test_that("Testing sampler assignment and misspecification of priors for conc parameter", { 
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dgamma(1, 1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  conf <- configureMCMC(m)
  expect_equal(conf$getSamplers()[[5]]$name, "CRP_concentration")
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dexp(1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  conf <- configureMCMC(m)
  expect_equal(conf$getSamplers()[[5]]$name, "RW")
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dunif(0,1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  conf <- configureMCMC(m)
  expect_equal(conf$getSamplers()[[5]]$name, "RW")
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dnorm(-10,1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  ## we do not warn of negative concentration values because there could be many such
  ## warnings in certain MCMC samplers for the concentration parameter
  expect_failure(expect_output(m$simulate(), "value of concentration parameter"))
  expect_error(m$calculate())
  ## think about better way to tell the user that the prior for alpha is wrong
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(0, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4)))
  ## we do not warn of negative concentration values because there could be many such
  ## warnings in certain MCMC samplers for the concentration parameter
  expect_failure(expect_output(m$simulate(), "value of concentration parameter has to be larger than zero"))
  expect_error(m$calculate())
  
})


## testing misspecification of dimension in a model

test_that("Testing of misspecification of dimension when using CRP", { 
  
  ## more labels than observations
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:10] ~ dCRP(conc=1, size=10)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,10), mu=rnorm(4)))
  conf <- configureMCMC(m)
  expect_error(buildMCMC(conf))
  
  
  ## more observations than labels 
  code = nimbleCode({
    for(i in 1:10) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  expect_error(nimbleModel(code, data = list(y = rnorm(10)),
                           inits = list(xi = rep(1,4), mu=rnorm(10))))
  
  
  ## different obervations with same label
  code = nimbleCode({
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
    y[1] ~ dnorm(mu[xi[1]], 1)
    y[2] ~ dnorm(mu[xi[1]], 1)
    xi[1:2] ~ dCRP(conc=1, size=2)
  })
  m <- nimbleModel(code, data = list(y = rnorm(2)),
                   inits = list(xi = rep(1,2), mu=rnorm(2)))
  conf <- configureMCMC(m)
  expect_error(buildMCMC(conf))
  
  
  ## same obervation with different label
  code = nimbleCode({
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
    y[1] ~ dnorm(mu[xi[1]], 1)
    y[1] ~ dnorm(mu[xi[2]], 1)
    xi[1:2] ~ dCRP(conc=1, size=2)
  })
  expect_error(nimbleModel(code, data = list(y = rnorm(2)),
                           inits = list(xi = rep(1,2), mu=rnorm(2))))
  
  ## less tilde variables than observations
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)  
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=1)  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50)))
  conf <- configureMCMC(m)
  expect_warning(buildMCMC(conf))
  
  
  
  ## multiple tilde parameters
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)
      s2[i] ~ dinvgamma(1,1)
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=s2[xi[i]])  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50), s2=rinvgamma(50,1,1)))
  conf <- configureMCMC(m)
  expect_warning(buildMCMC(conf))
  
  
  ## multiple tilde parameters, one is common for every observation
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)
      s2[i] ~ dinvgamma(1,1)
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=s2[xi[1]])  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50), s2=rinvgamma(50,1,1)))
  conf <- configureMCMC(m)
  expect_error(buildMCMC(conf))
  
  ## more than one label used for each observation
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)  
    }
    for(i in 1:99){
      y[i] ~ dnorm(mu[xi[i]]+mu[xi[i+1]], var=1)  
    }
    y[100] ~ dnorm(mu[xi[100]], 1)
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50)))
  conf <- configureMCMC(m)
  expect_error(buildMCMC(conf))
  
  ## test that a message is sent when more tilde variables than defined are needed
  code = nimbleCode({
    for(i in 1:3){
      mu[i] ~ dnorm(0,1)  
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=1)  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = c(rnorm(20, -5) , rnorm(20, 0), rnorm(20, 5),
                                           rnorm(20, 10), rnorm(20, 20))),
                   inits = list(xi = rep(1,100), mu=rnorm(3)))
  cm <- compileNimble(m)
  conf <- configureMCMC(m)
  expect_warning(mMCMC <- buildMCMC(conf))
  cmMCMC=compileNimble(mMCMC, project=m, resetFunctions=TRUE)
  set.seed(1)
  expect_output(cmMCMC$run(1), "CRP_sampler: This MCMC is not fully nonparametric.")
  
})



## Test real BNP models:

test_that("Testing BNP model based on CRP", { 
  
  ## DPM of Poisson distribution:
  code <- nimbleCode({
    for(i in 1:n){
      lambda[i] ~ dgamma(shape=1, rate=0.01)
      y[i] ~ dpois(lambda[xi[i]])
    }
    xi[1:n] ~ dCRP(conc = conc0, size=n)
  })
  
  data0 <- list(y = c(rpois(20, 10), rpois(20, 5), rpois(60, 50)))
  Consts <- list(n = 100, conc0 = 1)
  Inits <- list( xi = 1:Consts$n, lambda = rgamma(Consts$n, shape=1, rate=0.01))
  m <- nimbleModel(code, data = data0, inits = Inits, constants = Consts,  calculate=TRUE)
  cm <- compileNimble(m) 
  
  mConf <- configureMCMC(m)
  mMCMC <- buildMCMC(mConf)
  expect_equal(class(mMCMC$samplerFunctions[[101]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  CmMCMC <- compileNimble(mMCMC, project=m, resetFunctions=TRUE)
  
  nburn <- 500
  nsave <- 100
  ntotal <- nburn + nsave
  itersave <- nburn + (1:nsave)
  set.seed(1)
  CmMCMC$run(ntotal)
  
  ## results:
  n <- Consts$n
  samples <- as.matrix(CmMCMC$mvSamples)
  lamSam <- samples[itersave, 1:n]
  xiSam <- samples[itersave, (n+1):(2*n)]
  kSam <- apply(xiSam, 1, function(x)length(unique(x)))
  
  ygrid <- 0:100
  fSam <- matrix(0, ncol = length(ygrid), nrow = nsave)
  
  lamPost <- c()
  Trunc <- 25
  vj=c(rbeta(Trunc-1, 1, Consts$conc0+Consts$n), 1)
  wj=c(vj[1], vj[2:(Trunc-1)]*cumprod(1-vj[1:(Trunc-2)]), cumprod(1-vj[1:(Trunc-1)])[Trunc-1])
  probs=c(rep(1/(Consts$conc0+Consts$n), Consts$n), Consts$conc0/(Consts$conc0+Consts$n))
  for(i in 1:nsave){
    for(j in 1:Trunc){
      index=sample(1:(Consts$n+1), size=1, prob=probs)
      if(index==(Consts$n+1)){
        lamPost[j] <- rgamma(1,shape=1, rate=0.01)
      }else{
        lamPost[j]=lamSam[i, xiSam[i,index]]
      }
    }
    fSam[i, ]=sapply(ygrid, function(x)sum(wj*dpois(x, lamPost)))
  }
  fHat=apply(fSam, 2, mean)
  
  f0 <- function(x) 0.2*dpois(x, 10) + 0.2*dpois(x, 5) + 0.6*dpois(x, 50)
  f0grid <- sapply(ygrid, f0)
  
  L1dist <- mean(abs(f0grid - fHat))
  
  expect_equal(L1dist, 0.01, tol=0.01,
               info = "wrong estimation of density in DPM of Poisson distrbutions")
  
  
  ## DPM of Poisson distribution with prior for conc parameter:
  code <- nimbleCode({
    for(i in 1:n){
      lambda[i] ~ dgamma(shape=1, rate=0.01)
      y[i] ~ dpois(lambda[xi[i]])
    }
    xi[1:n] ~ dCRP(conc = alpha, size=n)
    alpha ~ dgamma(1, 1)
  })
  
  data0 <- list(y = c(rpois(20, 10), rpois(20, 5), rpois(60, 50)))
  Consts <- list(n = 100, conc0 = 1)
  Inits <- list( xi = 1:Consts$n, lambda = rgamma(Consts$n, shape=1, rate=0.01), alpha = 1)
  m <- nimbleModel(code, data = data0, inits = Inits, constants = Consts,  calculate=TRUE)
  cm <- compileNimble(m) 
  
  mConf <- configureMCMC(m, monitors = c('xi', 'lambda', 'alpha'))
  mMCMC <- buildMCMC(mConf)
  expect_equal(mConf$getSamplers()[[101]]$name, "CRP_concentration")
  expect_equal(class(mMCMC$samplerFunctions[[102]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  CmMCMC <- compileNimble(mMCMC, project=m, resetFunctions=TRUE)
  
  nburn <- 500
  nsave <- 100
  ntotal <- nburn + nsave
  itersave <- nburn + (1:nsave)
  set.seed(1)
  CmMCMC$run(ntotal)
  
  #-- results:
  n <- Consts$n
  samples <- as.matrix(CmMCMC$mvSamples)
  concSam <- samples[itersave, 1]
  lamSam <- samples[itersave, 2:(n+1)]
  xiSam <- samples[itersave, (n+2):ncol(samples)]
  kSam <- apply(xiSam, 1, function(x)length(unique(x)))
  
  ygrid <- 0:100
  fSam <- matrix(0, ncol = length(ygrid), nrow = nsave)
  
  lamPost <- c()
  aproxError <- 10^(-10)
  Trunc <- ceiling(log(aproxError)/log(mean(concSam)/(mean(concSam)+1))+1)
  for(i in 1:nsave){
    vj=c(rbeta(Trunc-1, 1, concSam[i]+Consts$n), 1)
    wj=c(vj[1], vj[2:(Trunc-1)]*cumprod(1-vj[1:(Trunc-2)]), cumprod(1-vj[1:(Trunc-1)])[Trunc-1])
    probs=c(rep(1/(concSam[i]+Consts$n), Consts$n), concSam[i]/(concSam[i]+Consts$n))
    
    for(j in 1:Trunc){
      index=sample(1:(Consts$n+1), size=1, prob=probs)
      if(index==(Consts$n+1)){
        lamPost[j] <- rgamma(1,shape=1, rate=0.01)
      }else{
        lamPost[j]=lamSam[i, xiSam[i,index]]
      }
    }
    fSam[i, ]=sapply(ygrid, function(x)sum(wj*dpois(x, lamPost)))
  }
  fHat=apply(fSam, 2, mean)
  
  #plot(ygrid, fHat, type="l")
  #lines(ygrid, f0grid)
  
  L1dist <- mean(abs(f0grid - fHat))
  
  expect_equal(L1dist, 0.01, tol=0.01,
               info = "wrong estimation of density in DPM of Poisson distrbutions")
  
  
  ## DPM of normal densities with random means and variances
  Code=nimbleCode(
    {
      for(i in 1:50){
        thetatilde[i] ~ dnorm(mean=0, var=100) 
        s2tilde[i] ~ dinvgamma(shape=1, scale=0.01) 
      }
      for(i in 1:100){
        y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])
      }
      xi[1:100] ~ dCRP(conc=1, size=100)
    }
  )
  
  set.seed(1)
  aux <- sample(1:10, size=100, replace=TRUE)
  Inits <- list(xi = aux,
                thetatilde = rnorm(50, 0, sqrt(100)),
                s2tilde = rinvgamma(50, shape=1, scale=0.01))
  Data <- list(y = c(rnorm(50, 5,sqrt(4)), rnorm(50, -5, sqrt(2))))
  
  m <- nimbleModel(Code, data=Data, inits=Inits, calculate=TRUE)
  cmodel <- compileNimble(m)
  
  mConf <- configureMCMC(m)
  expect_warning(mMCMC <- buildMCMC(mConf))
  expect_equal(class(mMCMC$samplerFunctions[[101]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  CmMCMC=compileNimble(mMCMC, project=m, resetFunctions=TRUE)
  
  nburn <- 500
  nsave <- 100
  ntotal <- nburn + nsave
  itersave <- nburn + (1:nsave)
  set.seed(1)
  CmMCMC$run(ntotal)
  
  N <- 100
  samples <- as.matrix(CmMCMC$mvSamples)
  xisamples <- samples[, (100+1):(2*N)]
  s2samples <- samples[, 1:50]
  thetasamples <- samples[, (50+1):(100)]
  
  conc <- 1
  sec=seq(-20,20,len=100)
  Trunc=35
  thetaNew=matrix(0, ncol=2, nrow=Trunc)
  vj=c(rbeta(Trunc-1, 1, N+conc), 1)
  wj=c(vj[1], vj[2:(Trunc-1)]*cumprod(1-vj[1:(Trunc-2)]), cumprod(1-vj[1:(Trunc-1)])[Trunc-1])
  probs=c(rep(1/(N+conc), N), conc/(N+conc))
  fsam=matrix(0,ncol=length(sec), nrow=nsave)
  for(i in 1:nsave){
    for(j in 1:Trunc){
      index=sample(1:(N+1), size=1, prob=probs)
      if(index==(N+1)){
        thetaNew[j,1] <- rnorm(1, 0, sqrt(100))
        thetaNew[j,2] <- rinvgamma(1, 1, 0.01)
      }else{
        thetaNew[j,]=c(thetasamples[i, xisamples[i,index]],s2samples[i, xisamples[i,index]])
      }
    }
    fsam[i, ]=sapply(sec, function(x)sum(wj*dnorm(x,thetaNew[,1],sqrt(thetaNew[,2]))))
  }
  fHat=apply(fsam, 2, mean)
  
  f0sec <- sapply(sec, function(x) 0.5*dnorm(x, 5, sqrt(4)) + 0.5*dnorm(x, -5, sqrt(4)))
  L1dist <- mean(abs(f0sec - fHat))
  
  expect_equal(L1dist, 0.01, tol=0.01,
               info = "wrong estimation of density in DPM of Poisson distrbutions")
  
})







test_that("Testing more BNP models based on CRP", { 
  ## Avandia meta-analysis
  codeBNP <- nimbleCode({
    for(i in 1:nStudies) {
      y[i] ~ dbin(size = nStudies, prob = q[i])
      x[i] ~ dbin(size = nStudies, prob = p[i])
      q[i] <- expit(theta + gamma[i])
      p[i] <- expit(gamma[i])
      gamma[i] ~ dnorm(mu[i], var = tau[i])
      mu[i] <- muTilde[xi[i]]
      tau[i] <- tauTilde[xi[i]]
    }
    for(i in 1:nStudies) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
      tauTilde[i] ~ dinvgamma(a0, b0)
    }
    xi[1:nStudies] ~ dCRP(conc, size = nStudies)
    conc ~ dgamma(1, 1)
    mu0 ~ dflat()
    sd0 ~ dunif(0, 100)
    a0 ~ dunif(0, 100)
    b0 ~ dunif(0, 100)
    theta ~ dflat()
  })
  
  Consts=list(nStudies=10)
  set.seed(1)
  Inits=list(gamma=rep(1,10),
             muTilde=rep(1,10),
             tauTilde=rep(1,10),
             xi=rep(1,10),
             conc =1,
             mu0 = 0,
             sd0 = 1,
             a0 = 1,
             b0 = 1,
             theta = 0)
  
  Data=list(y=rbinom(10, 10, 0.5), x=rbinom(10, 10, 0.5))
  
  model<-nimbleModel(codeBNP, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
  cmodel<-compileNimble(model)
  
  mConf <- configureMCMC(model)
  mMCMC <- buildMCMC(mConf)
  
  expect_equal(mConf$getSamplers()[[1]]$name, "CRP_concentration")
  expect_equal(mConf$getSamplers()[[7]]$name, "CRP")
  expect_equal(class(mMCMC$samplerFunctions[[7]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## Using testBUGSmodel
  model <- function() {
    for(i in 1:10) {
      y[i] ~ dbin(size = 10, prob = q[i])
      x[i] ~ dbin(size = 10, prob = p[i])
      q[i] <- expit(theta + gamma[i])
      p[i] <- expit(gamma[i])
      gamma[i] ~ dnorm(mu[i], var = tau[i])
      mu[i] <- muTilde[xi[i]]
      tau[i] <- tauTilde[xi[i]]
    }
    for(i in 1:10) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
      tauTilde[i] ~ dinvgamma(a0, b0)
    }
    xi[1:10] ~ dCRP(conc, size = 10)
    conc ~ dgamma(1, 1)
    mu0 ~ dflat()
    sd0 ~ dunif(0, 100)
    a0 ~ dunif(0, 100)
    b0 ~ dunif(0, 100)
    theta ~ dflat()
  }
  testBUGSmodel(example = 'test8', dir = "",
                model = model, data = Data, inits = Inits,
                useInits = TRUE)
  
  
  ## myeloma semiparametric AFT model
  ## data from 'myeloma' in the 'emplik' package
  ## library(emplik)
  ## data(myeloma)
  time <- c(1.25,1.25,2,2,2,3,5,5,6,6,6,6,7,7,7,9,11,11,11,11,11,13,14,15,16,16,17,17,18,19,19,24,25,26,32,35,37,41,41,51,52,54,58,66,67,88,89,92,4,4,7,7,8,12,11,12,13,16,19,19,28,41,53,57,77)
  vstatus <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # 0 = alive (i.e., censored)
  logBUN <- c(2.2175,1.9395,1.5185,1.7482,1.301,1.5441,2.2355,1.6812,1.3617,2.1139,1.1139,1.415,1.9777,1.0414,1.1761,1.7243,1.1139,1.2304,1.301,1.5682,1.0792,0.7782,1.3979,1.6021,1.3424,1.3222,1.2304,1.5911,1.4472,1.0792,1.2553,1.301,1,1.2304,1.3222,1.1139,1.6021,1,1.1461,1.5682,1,1.2553,1.2041,1.4472,1.3222,1.1761,1.3222,1.4314,1.9542,1.9243,1.1139,1.5315,1.0792,1.1461,1.6128,1.3979,1.6628,1.1461,1.3222,1.3222,1.2304,1.7559,1.1139,1.2553,1.0792)
  HGB <- c(9.4,12,9.8,11.3,5.1,6.7,10.1,6.5,9,10.2,9.7,10.4,9.5,5.1,11.4,8.2,14,12,13.2,7.5,9.6,5.5,14.6,10.6,9,8.8,10,11.2,7.5,14.4,7.5,14.6,12.4,11.2,10.6,7,11,10.2,5,7.7,10.1,9,12.1,6.6,12.8,10.6,14,11,10.2,10,12.4,10.2,9.9,11.6,14,8.8,4.9,13,13,10.8,7.3,12.8,12,12.5,14)
  
  n <- length(time)
  alive <- vstatus == 0
  cens_time <- rep(NA, n)
  cens_time[alive] <- time[alive]
  cens_time[!alive] <- Inf
  time[alive] <- NA
  
  logBUN <- (logBUN - mean(logBUN)) / sd(logBUN)
  HGB <- (HGB - mean(HGB)) / sd(HGB)
  
  ## accelerated failure time model per https://www4.stat.ncsu.edu/~ghosal/papers/PMR.pdf for Bayesian semiparametric AFT models
  codeAFT <- nimbleCode({
    for(i in 1:n) {
      x[i] ~ dweib(alpha, exp(lambda[i]))   # 'data' nodes
      is_cens[i] ~ dinterval(x[i], c[i])
      lambda[i] <-  inprod(Z[i, 1:p], delta[1:p]) + eta[i] 
      eta[i] <- etaTilde[xi[i]]
    }
    xi[1:n] ~ dCRP(conc, size = n)
    conc ~ dgamma(1, 1)
    for(i in 1:n){
      etaTilde[i] ~ dunif(b0, B0)
    }
    alpha ~ dunif(a0, A0)
    for(j in 1:p){
      delta[j] ~ dflat() 
    }
  })
  
  constants = list(b0 = -10, B0 = 10, a0 = 0.1, A0 = 10, p = 2, n = n, c
                   = cens_time, Z = cbind(logBUN, HGB))
  data = list(is_cens = as.numeric(alive), x = time)
  xInit <- rep(NA, n)
  xInit[alive] <- cens_time[alive] + 10
  inits = list(alpha = 1, delta = c(0, 0), conc = 1, 
               etaTilde = runif(n,constants$b0, constants$B0),
               xi = sample(1:3, n, replace = TRUE), x = xInit)
  
  model <- nimbleModel(codeAFT, constants = constants, data = data, inits = inits)
  conf = configureMCMC(model)
  mcmc = buildMCMC(conf)
  
  expect_equal(conf$getSamplers()[[1]]$name, "CRP_concentration")
  expect_equal(conf$getSamplers()[[70]]$name, "CRP")
  expect_equal(class(mcmc$samplerFunctions[[70]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## Using testBUGSmodel
  model <- function() {
    for(i in 1:n) {
      x[i] ~ dweib(alpha, 1+exp(lambda[i]))   # 'data' nodes
      is_cens[i] ~ dinterval(x[i], c[i])
      lambda[i] <-  inprod(Z[i, 1:p], delta[1:p]) + eta[i] 
      eta[i] <- etaTilde[xi[i]]
    }
    xi[1:n] ~ dCRP(conc, size = n)
    conc ~ dgamma(1, 1)
    for(i in 1:n){
      etaTilde[i] ~ dunif(b0, B0)
    }
    alpha ~ dunif(a0, A0)
    for(j in 1:p){
      delta[j] ~ dflat() 
    }
  }
  
  Data = list(is_cens = as.numeric(alive), x = time, 
              b0 = -10, B0 = 10, a0 = 0.1, A0 = 10, p = 2, n = n,
              c = cens_time, Z = cbind(logBUN, HGB))
  xInit <- rep(NA, n)
  xInit[alive] <- cens_time[alive] + 10
  Inits = list(alpha = 1, delta = c(0, 0), conc = 1, 
               etaTilde = runif(n,Data$b0, Data$B0),
               xi = sample(1:3, n, replace = TRUE), x = xInit)
  
  testBUGSmodel(example = 'test9', dir = "",
                model = model, data = Data, inits = Inits, 
                useInits = TRUE)
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

##-- test: random sampling from a compiled model adding one more level:
test_that("random sampling from model works fine", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  
  SB_code2 <- nimbleCode({
    for(i in 1:5) 
      z[i] ~ dbeta(1, 1)
    w[1:6] <- stick_breaking(z[1:5])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  data <- list(xi = rep(1,10))
  SB_model2 <- nimbleModel(SB_code2, data=data, inits=Inits)
  
  c_SB_model2 <- compileNimble(SB_model2)
  
  c_SB_model2$z <- x
  c_SB_model2$calculate()
  
  expect_equal(c_SB_model2$w, truth,
               info = paste0("incorrect stick breaking weights in SB_model2"))
  
  #-- sampling via simulate:
  set.seed(0)
  simul_samp <- function(model) {
    model$simulate()
    return(model$w)
  }
  simul_samps <- t(replicate(10000, simul_samp(c_SB_model2)))
  
  trueE <- c(0.5^(1:5) )
  
  #-- checking the mean of the components of a vector that has a generalized dirichelt distribution
  #-- if z_i ~ beta then weights, w, defined by a SB representation have a generalized dirichlet distribution
  #-- and the expectation of w_j=0.5^j (in this case a=b=1).
  
  expect_equal(apply(simul_samps, 2, mean)[1:5], trueE, tol=0.01,
               info = paste0("incorrect weights (w) sampling  in SB_model2"))
  
  
  # wrong specification of stick variables
  SB_code3 <- nimbleCode({
    for(i in 1:5) 
      z[i] ~ dgamma(10, 10)
    w[1:6] <- stick_breaking(z[1:5])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  # wrong prior and starting values for stick variables: how to recognize this wrong specification?
  set.seed(1)
  Inits <- list(z = rgamma(5, 10, 10))
  data <- list(xi = 1:10)
  expect_output(nimbleModel(SB_code3, data=data, inits=Inits)) # message is sent because z >1.
  
  # good stating values for stick variables
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  data <- list(xi = rep(1,10))
  SB_model3 <- nimbleModel(SB_code3, data=data, inits=Inits)
  expect_output(SB_model3$simulate(), "values in 'z' have to be in") # message is sent because z >1.
  
  # wrong specification of length in stick variables, should be 5
  SB_code4 <- nimbleCode({
    for(i in 1:4) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:4])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  # wrong length for stick variables, a warning is sent in the nimbleModel function
  # how to recognize this wrong specification in the test?
  set.seed(1)
  Inits <- list(z = rbeta(4, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(nimbleModel(SB_code4, data=data, inits=Inits), message = "number of items to replace")
  
  
  # wrong specification of length in stick variables, should be 5
  # no warning in nimbleModel function
  SB_code5 <- nimbleCode({
    for(i in 1:2) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:2])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(2, 10, 10))
  data <- list(xi = rep(1,10))
  SB_model5 <- nimbleModel(SB_code5, data=data, inits=Inits)
  cSB_model5 <- compileNimble(SB_model5)
  expect_output(cSB_model5$calculate('w'), "Error in mapCopy")
  
  
  
  # longer vector of stick variables
  SB_code6 <- nimbleCode({
    for(i in 1:10) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:10])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  # wrong length for stick variables, a warning is sent in the nimbleModel function
  set.seed(1)
  Inits <- list(z = rbeta(10, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(SB_model6 <- nimbleModel(SB_code6, data=data, inits=Inits), "number of items to replace")
  cSB_model6 <- compileNimble(SB_model6)
  expect_output(cSB_model6$calculate('w'), "Error in mapCopy")
})


