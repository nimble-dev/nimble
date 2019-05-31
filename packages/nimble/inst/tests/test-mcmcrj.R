########## TO REMOVE ##########
# library(nimble)
# library(testthat)
# source("~/github/nimble/packages/nimble/R/MCMC_RJ.R")
########## TO REMOVE ##########

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
RwarnLevel <- options('warn')$warn
options(warn = 1)

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context('Testing of MCMC_RJ functionality')

test_that("Test configureRJ with no indicator variables", {

  ## Linear regression with 2 covariates, one in the model
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    beta1 ~ dnorm(0, sd = 100)
    beta2 ~ dnorm(0, sd = 100)
    sigma ~ dunif(0, 100)
    for(i in 1:50) {
      Ypred[i] <- beta0 + beta1 * x1[i] + beta2 * x2[i]
      Y[i] ~ dnorm(Ypred[i], sd = sigma)
    }
  })
  
  ## Data simulation
  set.seed(0)
  x1 <- runif(50, -1, 1)
  x2 <- runif(50, -1, 1)
  Y <- rnorm(50, 1.5 + 2 * x1, sd = 1)
  
  data   <- list(Y = Y, x1 = x1, x2 = x2)
  inits  <- list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y))
  
  m <- nimbleModel(code, data=data, inits=inits)
  mConf <- configureMCMC(m)
   
  ## One node
  nodes <-  c("beta2")
  expect_error(configureRJ(mConf, nodes), 
               'Provide indicatorNodes or priorProb vector')
  
  #####################################
  ## One node, multiple parameters 
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(fixedValue = c(0,1))), 
               'inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(mean = c(0,1))), 
               'inconsistent length of mean argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(scale = c(2,1))), 
               'inconsistent length of scale argument and specified number of RJ targetNodes')
  
  ## priorProb not probabilities
  expect_error(configureRJ(mConf, nodes, prior = -1))
  expect_error(configureRJ(mConf, nodes, prior = 2))
  
  #####################################
  ## Multiple nodes, less paramters
  nodes <-  c("beta0", "beta1", "beta2")
  expect_error(configureRJ(mConf, nodes, prior = c(0.5, 0.5)), 
               'Length of priorProb vector must match targetNodes length')
  
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(fixedValue = c(0,1))), 
               'inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(mean = c(0,1))), 
               'inconsistent length of mean argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(scale = c(2,1))), 
               'inconsistent length of scale argument and specified number of RJ targetNodes')
  
  ## priorProb not probabilities
  expect_error(configureRJ(mConf, nodes, prior = c(0.5, 2, 0.2)), 
               'Elements in priorProb must be probabilities in')
  

})
  

test_that("Test configureRJ with multivariate node - no indicator", {
  ##############################
  ## Multivariate node
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    
    mu[1:5] <- rep(0, 5)
    sigma[1:5] <- 1/rep(100, 5)
    simgma.mat[1:5, 1:5] <- diag(sigma[1:5])
    beta[1:5] ~ dmnorm(mu[1:5], sigma_mat[1:5, 1:5])
    
    for(i in 1:10) {
      Ypred[i] <- beta0 + sum(X[i,1:5]*beta[1:5])
      Y[i] ~ dnorm(Ypred[i], sd = sigma.y)
    }
    sigma.y ~ dunif(0, 100)
  })
  
  ## simulate some data
  set.seed(1)
  X <- matrix(rnorm(10*5), 10, 5) 
  betaTrue <- c(2, -2, 3, 0, 0)
  eps <- rnorm(10)
  Y <- as.vector(X%*%betaTrue + eps)
  
  data   <- list(Y = Y, X = X)
  inits <- list(beta0 = 0, beta = rep(0, 5), sigma.y = sd(Y), sigma_mat = diag(rep(1/100, 5)), mu = rep(0, 5))
  
  m <- nimbleModel(code, data=data, inits=inits)
  mConf <- configureMCMC(m)
  
  ## test multivariate node with joint sampler
  expect_error(configureRJ(mConf, "beta", prior =0.5), 
              'beta is multivariate using a joint sampler; only univariate samplers can be used')

  ## test multivariate node with univariate samplers 
  nodeAsScalar <- mConf$model$expandNodeNames("beta", returnScalarComponents = TRUE)
  ## acceptable case
  mConf$removeSamplers("beta")
  for(node in nodeAsScalar){
    mConf$addSampler(node, type = "RW")
  }
  
  targetNodes <-  c("beta")
  control <-  list(fixedValue = 0, mean = 0, scale = 2)
  ## this should work
  expect_error(configureRJ(mcmcConf = mConf, targetNodes = targetNodes, priorProb = 0.5, control = control), NA)
  
  ## test double call to configureRJ
  expect_error(configureRJ(mcmcConf = mConf, targetNodes = targetNodes, priorProb = 0.5, control = control), 
  'is already configure for reversible jump')

})


test_that("Check passing node vector - no indicator", {
  #####################################
  ## Vector node
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)

    for(i in 1:5){
      beta[i] ~ dnorm(0, sd = 100)
    }
    sigma ~ dunif(0, 100)
    for(i in 1:10) {
      Ypred[i] <- beta0 + sum(X[i,1:5]*beta[1:5])
      Y[i] ~ dnorm(Ypred[i], sd = sigma)
    }
  })

  ## simulate some data
  set.seed(1)
  X <- matrix(rnorm(10*5), 10, 5)
  betaTrue <- c(2, -2, 3, 0, 0)
  eps <- rnorm(10)
  Y <- as.vector(X%*%betaTrue + eps)

  data   <- list(Y = Y, X = X)
  inits  <- list(beta0 = 0, beta = rep(0, 10), sigma = sd(Y))

  m <- nimbleModel(code, data=data, inits=inits)
  mConf <- configureMCMC(m)

  ## no error
  expect_error(configureRJ(mConf, c("beta"), prior = 0.5), NA)
  
  mConf <- configureMCMC(m)
  expect_error(configureRJ(mConf, c("beta[1]", "beta[2:4]"), prior = 0.5), NA)
  
  mConf <- configureMCMC(m)
  expect_error(configureRJ(mConf, c("beta[1]", "beta[2:4]"), prior = c(0.5, 0.2)), NA)
})



test_that("Check sampler_RJ behaviour - no indicator", {
  
  ## Linear regression with 2 covariates, one in the model
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    beta1 ~ dnorm(0, sd = 100)
    beta2 ~ dnorm(0, sd = 100)
    sigma ~ dunif(0, 100)
    for(i in 1:50) {
      Ypred[i] <- beta0 + beta1 * x1[i] + beta2 * x2[i]
      Y[i] ~ dnorm(Ypred[i], sd = sigma)
    }
  })
  
  ## Data simulation
  set.seed(0)
  x1 <- runif(50, -1, 1)
  x2 <- runif(50, -1, 1)
  Y <- rnorm(50, 1.5 + 2 * x1, sd = 1)
  
  data   <- list(Y = Y, x1 = x1, x2 = x2)
  
  ## check sampler behaviour 
  m <- nimbleModel(code, data=data)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('beta1', 'beta2'))
  configureRJ(mConf, c('beta1', 'beta2'), prior = 0.5)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC,  niter=1000, nburnin = 900, thin=1, 
                    inits = list(beta0 = 1, beta1 = 1, beta2 = 1, sigma = sd(Y)), setSeed = 1)
  
  ## beta2 should be more likely to be 0
  expect_true(sum(output[, 'beta2'] == 0)/100 > 0.5)
  # expect_true(mean(output[which(output[, 'beta2'] != 0), 'beta2']) - coef(lm(Y ~ x1 + x2))[3] < 0.05) ## should check that beta2 is small when in the model
  
  ## beta1 should be less likely to be 0
  expect_true(sum(output[, 'beta1'] == 0)/100 < 0.5)
  ## beta1 estimate (comparison with lm estimate)
  expect_equal(mean(output[which(output[, 'beta1'] != 0), 'beta1']), as.numeric(coef(lm(Y ~ x1 + x2))[2]) , tolerance=0.1, scale = 1)

  # ## beta1 should be in the model in last 100 iterations (chain has converged)
  # expect_false(any(output[, 'beta1'] == 0))
  
  #######
  ## change proposal mean for beta1 - still reasonable even if far
  ## dnorm(1.5, 3, 1) = 0.12
  m <- nimbleModel(code, data=data)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('beta1'))
  
  configureRJ(mConf, 'beta1', prior = 0.5, control = list(mean = 3))

  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC,  niter=100, thin=1,
                    inits = list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y)), setSeed = 1)
  ## beta1 estimate (comparison with lm estimate)
  expect_equal(mean(output[which(output[, 'beta1'] != 0), 'beta1']), as.numeric(coef(lm(Y ~ x1 + x2))[2]) , tolerance=0.1, scale = 1)
  
  #######
  ## fixed value on true beta1
  m <- nimbleModel(code, data=data)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('beta1'))
  
  configureRJ(mConf, 'beta1', prior = 0.5, control = list(fixedValue = 1.5))
  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC,  niter=100, thin=1,
                    inits = list(beta0 = 0, beta1 = 0, beta2 = 0,  sigma = sd(Y)), setSeed = 1)
  expect_equal(mean(output[which(output[, 'beta1'] != 0), 'beta1']), 1.5 , tolerance=0.01, scale = 1)
  
  #######
  ## fixedValue on far value for beta2 
  m <- nimbleModel(code, data=data)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('beta2'))
  
  configureRJ(mConf, 'beta2', prior = 0.5, control = list(fixedValue = 5))
  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC,  niter=100, thin=1,
                    inits = list(beta0 = 1, beta1 = 1, beta2 = 1, sigma = sd(Y)), setSeed = 1)
  
  ## still beta2 is in the models but really small
  expect_equal(mean(output[which(output[, 'beta2'] != 0), 'beta2']), 0 , tolerance=0.1, scale = 1)
  
  
  if(.Platform$OS.type != "windows") {
    nimble:::clearCompiled(m)
  }
})


######################################
## Tests using indicator variables
######################################

test_that("Test configureRJ with indicator variables", {
  
  ## Linear regression with 2 covariates, one in the model
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    beta1 ~ dnorm(0, sd = 100)
    beta2 ~ dnorm(0, sd = 100)
    sigma ~ dunif(0, 100)
    z1 ~ dbern(psi)  ## indicator variable for including beta2
    z2 ~ dbern(psi)  ## indicator variable for including beta2
    psi ~ dbeta(1, 1)
    for(i in 1:50) {
      Ypred[i] <- beta0 + beta1 * z1 * x1[i] + beta2 * z2 * x2[i]
      Y[i] ~ dnorm(Ypred[i], sd = sigma)
    }
  })
  
  ## Data simulation
  set.seed(0)
  x1 <- runif(50, -1, 1)
  x2 <- runif(50, -1, 1)
  Y <- rnorm(50, 1.5 + 2 * x1, sd = 1)
  
  data   <- list(Y = Y, x1 = x1, x2 = x2)
  inits  <- list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y), z2 = 1, z1 = 1, psi = 0.5)
  
  m <- nimbleModel(code, data=data, inits=inits)
  mConf <- configureMCMC(m)
  
  ## One node
  nodes <-  c("beta2")
  expect_error(configureRJ(mConf, nodes), 
               'Provide indicatorNodes or priorProb vector')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2")), 
               'Length of indicatorNodes vector must match targetNodes length')
  
  ## One node, multiple parameters
  # expect_error(configureRJ(mConf, nodes, indicatorNodes = "z1", control = list(fixedValue = c(0,1))), 
  #              'inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = "z1", control = list(mean = c(0,1))), 
               'inconsistent length of mean argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = "z1", control = list(scale = c(2,1))), 
               'inconsistent length of scale argument and specified number of RJ targetNodes')
  
  ## Multiple nodes, less paramters
  nodes <-  c("beta0", "beta1", "beta2")
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2")), 
               'Length of indicatorNodes vector must match targetNodes length')

  # expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2"), control = list(fixedValue = c(0,1))), 
  #              'inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2"), control = list(mean = c(0,1))), 
               'inconsistent length of mean argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2"), control = list(scale = c(2,1))), 
               'inconsistent length of scale argument and specified number of RJ targetNodes')  


  # if(.Platform$OS.type != "windows") {
  #   nimble:::clearCompiled(m)
  # }
})


test_that("Test configureRJ with multivariate node - indicator", {
  
  ##############################
  ## Multivariate node
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    
    mu[1:5] <- rep(0, 5)
    sigma[1:5] <- 1/rep(100, 5)
    simgma.mat[1:5, 1:5] <- diag(sigma[1:5])
    beta[1:5] ~ dmnorm(mu[1:5], sigma_mat[1:5, 1:5])
    
    for(i in 1:5){
      ## indicator variables
      z[i] ~ dbern(0.5)
    }
    
    for(i in 1:10) {
      Ypred[i] <- beta0 + sum(X[i,1:5]*beta[1:5]*z[1:5])
      Y[i] ~ dnorm(Ypred[i], sd = sigma.y)
    }
    sigma.y ~ dunif(0, 100)
  })
  
  ## simulate some data
  set.seed(1)
  X <- matrix(rnorm(10*5), 10, 5) 
  betaTrue <- c(2, -2, 3, 0, 0)
  eps <- rnorm(10)
  Y <- as.vector(X%*%betaTrue + eps)
  
  data   <- list(Y = Y, X = X)
  inits <- list(beta0 = 0, beta = rep(0, 5), sigma.y = sd(Y), sigma_mat = diag(rep(1/100, 5)), mu = rep(0, 5))
  
  m <- nimbleModel(code, data=data, inits=inits)
  mConf <- configureMCMC(m)
  
  ## test multivariate node with joint sampler
  expect_error(configureRJ(mConf, "beta", indicatorNodes = "z"), 
               'beta is multivariate using a joint sampler; only univariate samplers can be used')
  
  ## test multivariate node with univariate samplers 
  nodeAsScalar <- mConf$model$expandNodeNames("beta", returnScalarComponents = TRUE)
  ## acceptable case
  mConf$removeSamplers("beta")
  for(node in nodeAsScalar){
    mConf$addSampler(node, type = "RW")
  }
  
  ## this should work
  control <-  list(fixedValue = 0, mean = 0, scale = 2)
  expect_error(configureRJ(mcmcConf = mConf, targetNodes = "beta", indicatorNodes = "z", control = control), NA)
  
  ## test double call to configureRJ
  expect_error(configureRJ(mcmcConf = mConf, targetNodes = "beta", indicatorNodes = "z", control = control), 
               'is already configure for reversible jump')
  
  
  if(.Platform$OS.type != "windows") {
    nimble:::clearCompiled(m)
  }
  
  
})


test_that("Check sampler_RJ_indicator behaviour - indicator", {

  ## Linear regression with 2 covariates, one in the model
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    beta1 ~ dnorm(0, sd = 100)
    beta2 ~ dnorm(0, sd = 100)
    sigma ~ dunif(0, 100)
    z1 ~ dbern(psi)  ## indicator variable for including beta2
    z2 ~ dbern(psi)  ## indicator variable for including beta2
    psi ~ dbeta(1, 1)
    for(i in 1:50) {
      Ypred[i] <- beta0 + beta1 * z1 * x1[i] + beta2 * z2 * x2[i]
      Y[i] ~ dnorm(Ypred[i], sd = sigma)
    }
  })
  
  ## Data simulation
  set.seed(0)
  x1 <- runif(50, -1, 1)
  x2 <- runif(50, -1, 1)
  Y <- rnorm(50, 1.5 + 2 * x1, sd = 1)
  
  data   <- list(Y = Y, x1 = x1, x2 = x2)
  inits  <- list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y), z2 = 1, z1 = 1, psi = 0.5)
  
  m <- nimbleModel(code, data=data, inits=inits)
  cm <- compileNimble(m)
  
  ## check sampler behaviour 
  m <- nimbleModel(code, data=data, inits=inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('beta1', 'beta2', 'z1', 'z2'))
  configureRJ(mConf, c('beta1', 'beta2'), indicator =c('z1', 'z2'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE, resetFunctions = TRUE)
  output <- runMCMC(cMCMC,  niter=1000, nburnin = 900, thin=1, 
                    inits = list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y)), setSeed = 1)
  
  ## beta2 should be more likely to be 0
  expect_true(mean(output[, 'z2']) < 0.5)
  ## beta2 should be 0 when z1 is 0
  expect_equal(sum(output[, 'beta2'] != 0)/100, mean(output[, 'z2']) )

  ## beta1 should be less likely to be 0
  expect_true(mean(output[, 'z1']) > 0.5)
  ## beta1 should be 0 when z1 is 0
  expect_equal(sum(output[, 'beta1'] != 0)/100, mean(output[, 'z1']) )
  
  ## check beta1 estimate
  expect_equal(mean(output[which(output[, 'z1'] != 0), 'beta1']), as.numeric(coef(lm(Y ~ x1 + x2))[2]) , tolerance=0.1, scale = 1)

  
  ## more challeging data
  set.seed(0)
  x1 <- runif(50, -1, 1)
  x2 <- runif(50, -1, 1)
  Y <- rnorm(50, 1.5 + 1 * x1 - 1 * x2, sd = 1)
  
  data   <- list(Y = Y, x1 = x1, x2 = x2)
  inits  <- list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y), z2 = 1, z1 = 1, psi = 0.5)
  
  m <- nimbleModel(code, data=data, inits=inits)
  cm <- compileNimble(m)
  
  mConf <- configureMCMC(m, monitors = c('beta1', 'beta2', 'z1', 'z2'))
  configureRJ(mConf, c('beta1', 'beta2'), indicator =c('z1', 'z2'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE, resetFunctions = TRUE)
  output <- runMCMC(cMCMC,  niter=100, nburnin = 0, thin=1, 
                    inits = list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y)), setSeed = 1)
  
  ## check toggled_sampler
  ## when indicators are zero parameters are zero 
  expect_equal(which(output[, 'beta1'] == 0), which(output[, 'z1'] == 0))
  expect_equal(which(output[, 'beta2'] == 0), which(output[, 'z2'] == 0))
  
  if(.Platform$OS.type != "windows") {
    nimble:::clearCompiled(m)
  }
  
})



test_that("Check passing node vector - indicator", {
  #####################################
  ## Vector node
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    
    for(i in 1:5){
      beta[i] ~ dnorm(0, sd = 100)
      z[i] ~ dbern(psi[i])
      psi[i] ~ dbeta(1, 1)
    }
    
    sigma ~ dunif(0, 100)
    for(i in 1:10) {
      Ypred[i] <- beta0 + sum(X[i,1:5]*beta[1:5]*z[1:5])
      Y[i] ~ dnorm(Ypred[i], sd = sigma)
    }
  })
  
  ## simulate some data
  set.seed(1)
  X <- matrix(rnorm(10*5), 10, 5)
  betaTrue <- c(2, -2, 3, 0, 0)
  eps <- rnorm(10)
  Y <- as.vector(X%*%betaTrue + eps)
  
  data   <- list(Y = Y, X = X)
  inits  <- list(beta0 = 0, beta = rep(0, 5), z = rep(0, 5), psi = rep(0.5, 5), sigma = sd(Y))
  m <- nimbleModel(code, data=data, inits=inits)
  mConf <- configureMCMC(m)
  
  ## no error
  expect_error(configureRJ(mConf, targetNodes = "beta", indicatorNodes = "z"), NA)
  
  mConf <- configureMCMC(m)
  expect_error(configureRJ(mConf, c("beta[1]", "beta[2:4]"), indicatorNodes = c("z[1]", "z[2:4]")), NA)
  
  ## throws error
  mConf <- configureMCMC(m)
  expect_error(configureRJ(mConf, c("beta[1]", "beta[2:4]"), indicatorNodes = "z"), 
               'Length of indicatorNodes vector must match targetNodes length')

    
  if(.Platform$OS.type != "windows") {
    nimble:::clearCompiled(m)
  }
  
})
