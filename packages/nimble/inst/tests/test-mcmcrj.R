source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context('Testing of MCMC_RJ functionality')

test_that("Test configureRJ with no indicator variables", {

  ## Define model 2 covariate, one in the model
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
  
  ## One node, multiple parameters
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(fixedValue = c(0,1))), 
               'inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(mean = c(0,1))), 
               'inconsistent length of mean argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(scale = c(2,1))), 
               'inconsistent length of scale argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(positive = c(FALSE, FALSE))), 
               'inconsistent length of positive argument and specified number of RJ targetNodes')
  
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
  expect_error(configureRJ(mConf, nodes, prior = 0.5, control = list(positive = c(FALSE, FALSE))), 
               'inconsistent length of positive argument and specified number of RJ targetNodes')
  

  ## SHOULD I CHECK ALSO FOR THE TYPE OF VALUE PASSED?
  # mConf$printSamplers()
  
  
  ## indicator provided instead of prior not provided
  ## THIS RUNS WITH NO ISSUES BUT OUTPUT WILL BE WRONG - SHOULD WE PREVENT IT
  # nodes <-  c("beta2")
  # configureRJ(mConf, nodes, "beta0")
  # mConf$printSamplers()
  # 
  # mMCMC <- buildMCMC(mConf)
  # cm <-  compileNimble(m)
  # cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  # output <- runMCMC(cMCMC,  niter=1000, nburnin = 900, thin=1)

  
  ## SHOULD I CHECK ALSO FOR THE TYPE OF VALUE PASSED? eg positive = number
  # mConf$printSamplers()
  
  # if(.Platform$OS.type != "windows") {
  #   nimble:::clearCompiled(m)
  # }
})
  
test_that("Test configureRJ with indicator variables", {
  
  ## Define model 2 covariate, one in the model
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
  expect_error(configureRJ(mConf, nodes, indicatorNodes = "z1", control = list(fixedValue = c(0,1))), 
               'inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = "z1", control = list(mean = c(0,1))), 
               'inconsistent length of mean argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = "z1", control = list(scale = c(2,1))), 
               'inconsistent length of scale argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = "z1", control = list(positive = c(FALSE, FALSE))), 
               'inconsistent length of positive argument and specified number of RJ targetNodes')
  
  ## Multiple nodes, less paramters
  nodes <-  c("beta0", "beta1", "beta2")
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2")), 
               'Length of indicatorNodes vector must match targetNodes length')

  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2"), control = list(fixedValue = c(0,1))), 
               'inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2"), control = list(mean = c(0,1))), 
               'inconsistent length of mean argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2"), control = list(scale = c(2,1))), 
               'inconsistent length of scale argument and specified number of RJ targetNodes')
  expect_error(configureRJ(mConf, nodes, indicatorNodes = c("z1", "z2"), control = list(positive = c(FALSE, FALSE))), 
               'inconsistent length of positive argument and specified number of RJ targetNodes')
  


  # if(.Platform$OS.type != "windows") {
  #   nimble:::clearCompiled(m)
  # }
})

test_that("Check sampler_RJ behaviour (no indicator)", {
  
  ## Define model 2 covariate, one in the model
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
  # inits  <- list(beta0 = 1, beta1 = 0, beta2 = 0, sigma = sd(Y))
  
  # m <- nimbleModel(code, data=data)
  # cm <- compileNimble(m)
  
  ## check sampler behaviour 
  m <- nimbleModel(code, data=data)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('beta1', 'beta2'))
  configureRJ(mConf, c('beta1', 'beta2'), prior = 0.5)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE, resetFunctions = TRUE)
  output <- runMCMC(cMCMC,  niter=1000, nburnin = 900, thin=1, 
                    inits = list(beta0 = 1, beta1 = 1, beta2 = 1, sigma = sd(Y)), setSeed = 1)
  
  ##-------------------------------##
  ## not sure if this is a valid test
  ##-------------------------------##
  ## beta2 should be more likely to be 0
  expect_true(sum(output[, 'beta2'] == 0)/100 > 0.5)
  # expect_true(mean(output[which(output[, 'beta2'] != 0), 'beta2']) - coef(lm(Y ~ x1 + x2))[3] < 0.05) ## should check that beta2 is small when in the model
  
  ## beta1 should be in the model
  expect_false(any(output[, 'beta1'] == 0))
  ## check beta1 estimate
  expect_equal(mean(output[which(output[, 'beta1'] != 0), 'beta1']), as.numeric(coef(lm(Y ~ x1 + x2))[2]) , tolerance=0.2, scale = 1)
  
    ## Different fixedvalue?  
#   m <- nimbleModel(code, data=data)
#   cm <- compileNimble(m)
#   mConf <- configureMCMC(m, monitors = c('beta1', 'beta2'))
#   configureRJ(mConf, c('beta1', 'beta2'), prior = 0.5, control = list(fixedValue = c(2, 0)))
})

test_that("Check sampler_RJ_indicator behaviour (with indicator)", {

  ## Define model 2 covariate, one in the model
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
                    inits = list(beta0 = 1, beta1 = 1, beta2 = 1, sigma = sd(Y)), setSeed = 1)
  
  ##-------------------------------##
  ## not sure if this is a valid test
  ##-------------------------------##
  ## beta2 should be more likely to be 0
  expect_true(sum(output[, 'beta2'] == 0)/100 > 0.5)
  
  # expect_true(mean(output[which(output[, 'beta2'] != 0), 'beta2']) - coef(lm(Y ~ x1 + x2))[3] < 0.05) ## would check that beta2 is small when in the model
  ## beta1 should be in the model
  expect_false(any(output[, 'beta1'] == 0))
  ## check beta1 estimate
  expect_equal(mean(output[which(output[, 'beta1'] != 0), 'beta1']), as.numeric(coef(lm(Y ~ x1 + x2))[2]) , tolerance=0.2, scale = 1)
  

})