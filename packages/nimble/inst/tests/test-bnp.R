source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context('Testing of BNP functionality')

test_that("Test that sampleDPmeasure can be used for more complicated models", {
  set.seed(1)
  
  ## no deterministic node, conc param is fixed
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=a0)
      y[i] ~ dpois(lambdaTilde[xi[i]])
    }
    xi[1:10] ~ dCRP(conc = 1, size=10)
    a0 ~ dgamma(1, 1)
  })
  Inits <- list( xi = sample(1:2, size=10, replace=TRUE), 
                 lambdaTilde = rgamma(10, shape=1, rate=1), a0=1)
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  monitors <- c('lambdaTilde','xi', 'a0')
  mConf <- configureMCMC(m, monitors = monitors)
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC,  niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)
  expect_false(any(is.na(samplesG$samples)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  
  ## no deterministic node, random conc param
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=a0)
      y[i] ~ dpois(lambdaTilde[xi[i]])
    }
    xi[1:10] ~ dCRP(conc0, size=10)
    conc0 ~ dgamma(1,1)
    a0 ~ dgamma(1, 1)
  })
  Inits <- list( xi = sample(1:2, size=10, replace=TRUE), conc0=1,
                 lambdaTilde = rgamma(10, shape=1, rate=1), a0=1)
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('lambdaTilde','xi', 'conc0', 'a0'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  ## with deterministic node, conc param is fixed
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=a0)
      lambda[i] <- lambdaTilde[xi[i]]
      y[i] ~ dpois(lambda[i])
    }
    xi[1:10] ~ dCRP(1, size=10)
    a0 ~ dgamma(1, 1)
  })
  Inits <- list( xi = sample(1:2, size=10, replace=TRUE), 
                 lambdaTilde = rgamma(10, shape=1, rate=1), a0=1)
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('lambdaTilde','xi', 'a0'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter = 1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  ## with deterministic node, random conc param
  code <- nimbleCode({
    for(i in 1:10){
      lambdaTilde[i] ~ dgamma(shape=1, rate=0.01)
      lambda[i] <- lambdaTilde[xi[i]]
      y[i] ~ dpois(lambda[i])
    }
    xi[1:10] ~ dCRP(conc0, size=10)
    conc0 ~ dgamma(1,1)
  })
  Inits <- list( xi = sample(1:2, size=10, replace=TRUE), conc0=1,
                 lambdaTilde = rgamma(10, shape=1, rate=0.01))
  Data <- list(y = c(rpois(10, 8)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('lambdaTilde','xi', 'conc0'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  ## two cluster parameters, one deterministic parameter, fixed conc
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=1) 
        s2tilde[i] ~ dinvgamma(2, scale=1)
      }
      xi[1:10] ~ dCRP( 1 , size=10)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), 
             thetatilde=rnorm(10, 0, 1),
             s2tilde = rinvgamma(10, 2, scale=1))
  Data=list(y=c(rnorm(5,-5, 1), rnorm(5,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  
  ## two cluster parameters, one deterministic parameter, random conc
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=1) 
        s2tilde[i] ~ dinvgamma(2, scale=1)
      }
      xi[1:10] ~ dCRP( conc0 , size=10)
      conc0 ~ dgamma(1, 1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), 
             thetatilde=rnorm(10, 0, 1),
             s2tilde = rinvgamma(10, 2, scale=1), conc0=1)
  Data=list(y=c(rnorm(5,-5, 1), rnorm(5,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi', 'conc0'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  
  ## two cluster parameters, two deterministc parameter, random conc with random parameters
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(2, scale=1)
      }
      xi[1:10] ~ dCRP( conc0 , size=10)
      conc0 ~ dgamma(a,1)
      a ~ dgamma(1, 1)
      lambda ~ dgamma(1, 1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), 
             thetatilde=rnorm(10, 0, 1),
             s2tilde = rinvgamma(10, 2, scale=1), conc0=1, a=1, lambda=1)
  Data=list(y=c(rnorm(5,-5, 1), rnorm(5, 5, 1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi', 'conc0', 'a', 'lambda'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)
  expect_false(any(is.na(samplesG$samples)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  
  
  ## two cluster parameters, one deterministc parameter, random conc defined by two random parameters
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(2, scale=1)
      }
      xi[1:10] ~ dCRP( conc0 + conc1 , size=10)
      conc0 ~ dgamma(1,1)
      conc1 ~ dgamma(1,1)
      lambda ~ dgamma(1, 1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), 
             thetatilde=rnorm(10, 0, 1),
             s2tilde = rinvgamma(10, 2, scale=1), conc0=1, conc1=1,
             lambda=1)
  Data=list(y=c(rnorm(5,-5, 1), rnorm(5,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi', 'conc0', 'conc1',
                                          'lambda'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  
  ## two cluster parameters, one deterministic parameter, random conc defined by two random parameters
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(2, scale=1)
      }
      xi[1:10] ~ dCRP( conc0 + conc1 , size=10)
      conc1 ~ dgamma(1,1)
      conc0 ~ dgamma(a,b)
      b ~ dgamma(1,1)
      a ~ dgamma(1,1)
      lambda ~ dgamma(1,1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), 
             thetatilde=rnorm(10, 0, 1),
             s2tilde = rinvgamma(10, 2, scale=1), conc0=1, conc1=1,  a=1, b=1, lambda=1)
  Data=list(y=c(rnorm(5,-5, 1), rnorm(5,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi', 'conc0', 'conc1', 'a', 'b',  'lambda'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 0, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  ## two cluster parameters, one deterministic parameter, random conc defined by two random parameters, normal inverse gamma prior
  code=nimbleCode(
    {
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=s2tilde[i]/lambda) 
        s2tilde[i] ~ dinvgamma(2, scale=1)
      }
      xi[1:10] ~ dCRP( conc0 , size=10)
      conc0 ~ dgamma(a,b)
      b ~ dgamma(1,1)
      a ~ dgamma(1,1)
      lambda ~ dgamma(1, 1)
      for(i in 1:10){
        theta[i] <- thetatilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=sample(1:10, size=10, replace=TRUE), 
             thetatilde=rnorm(10, 0, 1),
             s2tilde = rinvgamma(10, 2, scale=1), conc0=1, a=1, b=1, lambda=1)
  Data=list(y=c(rnorm(5,-5, 1), rnorm(5,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi', 'conc0', 'a', 'b',  'lambda'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  out <- runMCMC(cMCMC, niter=1000, nburnin = 900, thin=1)
  samplesG <- getSamplesDPmeasure(cMCMC)$samples
  expect_false(any(is.na(samplesG)))
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
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
  expect_silent(output <- getSamplesDPmeasure(mMCMC))
  
  ## cluster variable not being monitored
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
  expect_silent(output <- getSamplesDPmeasure(mMCMC))
  
  ## concentration parameter not being monitored:
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
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'conc0'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  outputG <- getSamplesDPmeasure(mMCMC)
  
  ## concentration parameter deterministic parent not being monitored:
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
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6,  a = 1,d=1,  conc0=1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m, monitors = c('xi', 'mu'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'conc0'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'b'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'b', 'd'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_silent(outputG <- getSamplesDPmeasure(mMCMC))
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'd'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_silent(outputG <- getSamplesDPmeasure(mMCMC))
})





test_that("check iid assumption in sampleDPmeasure", {
  set.seed(1)
  
  ## univariate cluster parameters are not iid
  code <- nimbleCode({
    for(i in 1:10){
      muTilde[i] ~ dnorm(i, 1)
      y[i] ~ dnorm(muTilde[xi[i]], 1)
    }
    xi[1:10] ~ dCRP(conc = 1, size=10)
  })
  Inits <- list( xi = sample(1:2, size=10, replace=TRUE), 
                 muTilde = rep(1, 10))
  Data <- list(y = c(rnorm(10, 0,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('muTilde','xi'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC,  niter=1, nburnin = 0, thin=1)
  expect_error(samplesG <- getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }

  ## also not IID
  code=nimbleCode({
      xi[1:3] ~ dCRP(1, size = 3)
      thetatilde[1] ~ dnorm(0, 1)
      thetatilde[2] ~ dt(0, 1, 1)
      thetatilde[3] ~ dt(0, 1, 1)
      s2tilde[1] ~ dinvgamma(2, 1)
      s2tilde[2] ~ dgamma(1, 1)
      s2tilde[3] ~ dgamma(1, 1)
      for(i in 1:3){
          y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])
      }
  }
  )
  Inits <- list(xi = rep(1, 3), thetatilde=rep(0,3), s2tilde=rep(1,3))
  Data <- list(y = rnorm(3,-5, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi'))
  expect_silent(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  output <- cMCMC$run(1)
  expect_error(getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }

  
  ## one cluster param not with same distribution
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      thetatilde[1] ~ dnorm(0, 1)
      thetatilde[2] ~ dt(0, 1, 1)
      thetatilde[3] ~ dt(0, 1, 1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), thetatilde=rep(0,3))
  Data=list(y=rnorm(10, 0,1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  cMCMC$run(1)
  expect_error(getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  
  ## bivariate cluster parameters are not iid
  code <- nimbleCode({
    for(i in 1:10){
      mutilde[i] ~ dnorm(i, s2tilde[i]/lambda)
      s2tilde[i] ~ dinvgamma(2, 1)
      y[i] ~ dnorm(mutilde[xi[i]], var=s2tilde[xi[i]])
    }
    lambda ~ dgamma(1, 1)
    xi[1:10] ~ dCRP(conc = 1, size=10)
  })
  Inits <- list( xi = sample(1:2, size=10, replace=TRUE), 
                 mutilde = rep(1, 10), s2tilde=rep(1, 10), lambda=1)
  Data <- list(y = c(rnorm(10, 0,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m, monitors = c('mutilde','s2tilde', 'lambda', 'xi'))
  expect_error(mMCMC <- buildMCMC(mConf),
               "sampler_CRP: Cluster parameters must be conditionally independent")  
  
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      thetatilde[1] ~ dnorm(0, 1)
      thetatilde[2] ~ dt(0, 1, 1)
      thetatilde[3] ~ dt(0, 1, 1)
      s2tilde[1] ~ dinvgamma(2, 1)
      s2tilde[2] ~ dgamma(1, 1)
      s2tilde[3] ~ dgamma(1, 1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=rep(1, 10), thetatilde=rep(0,3), s2tilde=rep(1,3))#
  Data=list(y=rnorm(10, 0,1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  cMCMC$run(1)
  cMCMC$run(1, reset=FALSE)  # Claudia, why do we call $run twice?
  expect_error(getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
})




test_that("check use of epsilon parameters in getSamplesDPmeasure", {
  set.seed(1)
  
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
    sd0 ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)      
    mu0 ~ dnorm(0, var=10)
  })
  
  n <- 30
  constants <- list(n = n)
  data <- list(y = rnorm(n, 0, 1))
  inits <- list(alpha = 1, mu0 = 0, sd0 = 5, xi = rep(1, n),
                muTilde = rep(0,n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  outputG <- getSamplesDPmeasure(cmcmc)
  tr1 <- outputG$trunc
  
  outputG <- getSamplesDPmeasure(cmcmc, epsilon = 0.1)
  tr2 <- outputG$trunc
  
  outputG <- getSamplesDPmeasure(cmcmc, epsilon = 0.00001)
  tr3 <- outputG$trunc
  
  expect_true(tr1 > tr2,
              info='getSamplesDPmeasure: truncation level for larger epsilon incorrectly computed')
  expect_true(tr1 < tr3,
              info='getSamplesDPmeasure: truncation level for smaller epsilon incorrectly computed')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }
})

test_that("Test that new cluster parameters are correctly updated in CRP sampler", {
    ## Note: CP thinks these checks are not all that helpful - because of the small sample
    ## size, the estimates and truth are rather different so we need a high tolerance,
    ## so it's feasible that even an incorrectly-coded algorithm could pass the tests.
    
  set.seed(1)
  
  ## Data ~ Poisson(5). Starting values are extremely away from their true values. 
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      lambda[i] ~ dgamma(1, 0.01)
      y[i] ~ dpois(lambda[xi[i]])
    }
    alpha ~ dgamma(1, 1)
  })
  n <- 30
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
  
  output <- runMCMC(cMCMC, niter=1000, nburn=900, thin=1 , inits=Inits, setSeed=FALSE)
  
  xiSam <- output[, grep('xi', colnames(output))]
  lambdaSam <- output[, grep('lambda', colnames(output))]
  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  lambdaUnique <- sapply(1:nrow(output), function(i) unique(lambdaSam[i, xiSam[i, ]]) )
  cond <- sapply(1:nrow(output), function(i) sum(weights[[i]]*lambdaUnique[[i]]))
  
  expect_equal(mean(cond), 5, tol=2*1, scale=1,
               info = paste0("incorrect update of cluster parameters in Poisson data"))

  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
 
  ## normal - independent normal - inv gamma
  ## We start with only one active component and the data is a mixture of 3 normal ditributions
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i] ~ dnorm(0, var = s20/lambda)
      s2[i] ~ dinvgamma(2, 1)
      y[i] ~ dnorm(mu[xi[i]],  var = s2[xi[i]])
    }
    lambda ~ dgamma(1, 1)
    s20 ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)
  })
  n <- 30
  Consts <- list(n = n)
  Inits <- list(xi = rep(1, n), 
                mu = rep(-20, n), 
                s2 = rep(0.1, n),
                alpha = 200,
                lambda = 1,
                s20 = 1)
  thetas <- c(rep(-5, 10), rep(5, 10), rep(0, 10))
  Data <- list(y = rnorm(n, thetas, 1))
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 's2', 'alpha', 'lambda', 's20'), print=FALSE)  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  
  output <- runMCMC(cMCMC, niter=4000, nburnin=2000, thin=10 , inits=Inits, setSeed=FALSE)
  outputG <- getSamplesDPmeasure(cMCMC)
  
  Tr <- outputG$trunc
  samplesG <- outputG$samples
  grid <- seq(-15, 15, len=100)
  samF <- matrix(0, ncol=length(grid), nrow=nrow(samplesG))
  for(i in 1:nrow(samplesG)) {
      samF[i, ] <- sapply(grid, function(x)
          sum(samplesG[i, 1:Tr] * dnorm(x, samplesG[i, (2*Tr+1):(3*Tr)],
                                        sqrt(samplesG[i, (Tr+1):(2*Tr)]))))
  }
  
  ## distance between the estimation and truth at each point of the grid 
  f0 <- sapply(grid, function(x) dnorm(x, -5, 1)/3 + dnorm(x, 0, 1)/3 + dnorm(x, 5, 1)/3)
  for(i in 1:length(grid)) {
    expect_equal(mean(samF[,i]), f0[i], tolerance=2*0.1, scale=1,
                 info = paste("incorrect update of cluster parameters in mixture of normals data. Grid_i = ", i))
  }
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  
  ## normal - normal inverse gamma
  ## We start with only one active component and the data is a mixture of 3 normal ditributions
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
  n <- 30
  Consts <- list(n = n)
  Inits <- list(xi = rep(1, n), 
                mu = rep(-20, n), 
                s2 = rep(0.1, n),
                alpha = 200,
                lambda = 1)
  thetas <- c(rep(-5, 10), rep(5, 10), rep(0, 10))
  Data <- list(y = rnorm(n, thetas, 1))
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 's2', 'alpha', 'lambda'))  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  
  output <- runMCMC(cMCMC, niter=4000, nburnin=2000, thin=10 , inits=Inits, setSeed=FALSE)
  outputG <- getSamplesDPmeasure(cMCMC)
  
  Tr <- outputG$trunc
  samplesG <- outputG$samples
  grid <- seq(-15, 15, len=100)
  samF <- matrix(0, ncol=length(grid), nrow=nrow(samplesG))
  for(i in 1:nrow(samplesG)) {
      samF[i, ] <- sapply(grid, function(x)
          sum(samplesG[i, 1:Tr] * dnorm(x, samplesG[i, (2*Tr+1):(3*Tr)],
                                        sqrt(samplesG[i, (Tr+1):(2*Tr)]))))
  }
  
  ## distance between the estimation and truth at each point of the grid 
  f0 <- sapply(grid, function(x) dnorm(x, -5, 1)/3 + dnorm(x, 0, 1)/3 + dnorm(x, 5, 1)/3)
  for(i in 1:length(grid)) {
    expect_equal(mean(samF[,i]), f0[i], tolerance=2*0.1, scale=1,
                 info = paste("incorrect update of cluster parameters in mixture of normals data and conjugate normal - N-IG. Grid_i = ", i))
  }
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  ## conjugate normal - normal-inverse gamma 
  code=nimbleCode(
    {
      xi[1:4] ~ dCRP(1 , size=4)
      for(i in 1:4){
        thetatilde[i] ~ dnorm(0 , var=s2tilde[i]/lambda)
        s2tilde[i] ~ dinvgamma(2, 1)
        y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])
      }
      lambda ~ dgamma(1, 1)
    }
  )
  Inits=list(xi=c(1,1,1,1), thetatilde=c(0, -10, -10,-10), s2tilde=c(1,10,10,10), lambda=1)
  Data=list(y=c(0, 0, 10,10))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi', 'lambda'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m) 
  out <- runMCMC(cMCMC, niter=100, nburnin = 90, thin=1)
  means <- sapply(1:10, function(i) mean(out[i, 6:9][out[i, 10:13]]) )
  expect_equal(mean(means), 5, tolerance=2*1,
               info = 'wrong results from normal - normal - invgamma conjugate CRP')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
})




test_that("Test that the nonconjugate CRP sampler works fine ", {
  ## CP:  these are short runs and small sample sizes; we should consider more robust tests.
    
  set.seed(1)
  
  ## xi=1:n, the sampler requires samples from G_0
  code=nimbleCode(
    {
      xi[1:3] ~ dCRP(1 , size=3)
      for(i in 1:3){
        thetatilde[i] ~ dt(0, 1, 1) # implies non conjugate CRP sampler 
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=1:3, thetatilde=rep(-10,3))
  Data=list(y=rnorm(3,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m) 
  out <- runMCMC(cMCMC, niter=10, nburnin = 0, thin=1)
  means <- sapply(1:10, function(i) mean(out[i, 1:3][out[i, 4:6]]) )
  expect_equal(mean(means), 0, tolerance=2*1, info = 'wrong results from non conjugate CRP 1')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
    
  ## xi=1:n, the sampler requires samples from G_0, tildeVars have dependencies
  code=nimbleCode(
    {
      xi[1:3] ~ dCRP(1 , size=3)
      for(i in 1:3){
        thetatilde[i] ~ dt(mu0, 1, 1)   # implies non conjugate CRP sampler 
        y[i] ~ dnorm(thetatilde[xi[i]] + beta, var=1)
      }
      beta ~ dnorm(0, 1)
      mu0 ~ dnorm(0,var=10)
    }
  )
  Inits=list(xi=1:3, thetatilde=rep(0,3), beta=0, mu0=0)
  Data=list(y=rnorm(3,10, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi', 'beta', 'mu0'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m) 
  out <- runMCMC(cMCMC, niter=10, nburnin = 0, thin=1)
  means <- sapply(1:10, function(i) mean(out[i, 3:5][out[i, 6:8]] + out[i, 1]) )
  expect_equal(mean(means), 10, tolerance=2*1, info = 'wrong results from non conjugate CRP 2')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
})

test_that("Test opening of new clusters in CRP sampler ", {

  ## test for updating new cluster parameters, with conjugacy
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. 
  inits <- list(alpha = 1, mu0 = 0, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  ## now check that cmodel$muTilde[2] is near 50 and first obs has moved to 2nd cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_equal(output[1, 'muTilde[2]'], 50, tolerance = 2, info = 'incorrect update of parameter for second cluster',
               check.attributes = FALSE)
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=2), 'incorrect cluster for first obs')
  expect_identical(output[1, 'muTilde[1]'], c('muTilde[1]'=0), 'incorrect update of parameter for first cluster')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }

  ## test for updating new cluster parameters, without conjugacy
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. prior generates better value for first data point so new cluster should be opened.
  inits <- list(alpha = 1, mu0 = 50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] has changed and first obs has moved to 2nd cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_true(output[1, 'muTilde[2]'] != -50, 'incorrect update of parameter for second cluster')
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=2), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }

  ## test for (not) updating new cluster parameters, without conjugacy, no movement expected
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. prior generates even worse value for new cluster in terms of first obs, so no movement expected.
  inits <- list(alpha = 1, mu0 = -50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, 50, rep(0, n-2)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] is unchanged and first obs has stayed in only cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_true(output[1, 'muTilde[2]'] == 50, 'incorrect update of parameter for second cluster')
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=1), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }
  
  ## test for updating new cluster parameters, without conjugacy, movement to existing cluster expected
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. Prior generates even worse value for new cluster in terms of first obs, so no movement expected.
  inits <- list(alpha = 1, mu0 = -50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  inits$xi[1] <- 2
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] has remained and first obs has moved to first cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_true(output[1, 'muTilde[2]'] == -50, 'incorrect update of parameter for second cluster')
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=1), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }

  ## test for updating new cluster parameters, without conjugacy, singleton with movement to new cluster, same label, expected
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. prior generates better value for new cluster in terms of first obs, so movement expected, but singleton, so new cluster should be same ID as old cluster.
  inits <- list(alpha = 1, mu0 = 50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  inits$xi[1] <- 2
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] has changed and first obs has 'stayed' in 2nd cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_equal(output[1, 'muTilde[2]'], 50, tolerance = 5,
               info = 'incorrect update of parameter for second cluster',
               check.attributes = FALSE)
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=2), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }
})

test_that("Test reset frunction in CRP sampler ", {
  set.seed(1)
  
  ## the data set has more clusters than cluster parameters are available
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:2)
        thetatilde[i] ~ dnorm(0, 1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), thetatilde=rep(0,2))#
  Data=list(y=c(rnorm(3,-5, 1), rnorm(3,5, 1), rnorm(4, 0,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  expect_output(cMCMC$run(1), info='CRP_sampler: This MCMC is for a parametric model')
  cMCMC$run(1, reset=FALSE)
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
})


test_that("Test that cluster parameters and membership variable are independent in CRP sampler ", {

  ## membership variable depends on cluster params
  code=nimbleCode({
    for(i in 1:10) {
      muTilde[i] ~ dnorm(0, 1)  
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], 1)
    }
    xi[1:10] ~ dCRP(exp(muTilde[1]) , size=10)
  })
  Inits=list(xi=rep(1, 10), muTilde=rep(0,10))
  Data=list(y=rnorm(10,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(mConf),
               'sampler_CRP: Only the variables being clustered can depend on the cluster parameters')
  
  
  ## cluster params depend on membership variable
  code=nimbleCode({
    for(i in 1:10) {
      muTilde[i] ~ dnorm(log(xi[1]), 1)  
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], 1)
    }
    xi[1:10] ~ dCRP(1 , size=10)
  })
  Inits=list(xi=rep(1, 10), muTilde=rep(0,10))
  Data=list(y=rnorm(10,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(mConf),
               'sampler_CRP: Detected that the CRP variable is used in some way not as an index')
  
  
  ## one more node depends on membership variable
  code=nimbleCode({
    mu0 ~ dnorm(xi[1], 1) 
    for(i in 1:10) {
      muTilde[i] ~ dnorm(0, 1)  
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], 1)
    }
    xi[1:10] ~ dCRP(1 , size=10)
  })
  Inits=list(xi=rep(1, 10), muTilde=rep(0,10), mu0=0)
  Data=list(y=rnorm(10,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(mConf),
               'sampler_CRP: Detected that the CRP variable is used in some way not as an index')
    
  ## non related variable depends on cluster variable and membership variable
  code=nimbleCode({
    for(i in 1:10) {
      muTilde[i] ~ dnorm(0, 1)  
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], 1)
    }
    xi[1:10] ~ dCRP(1 , size=10)
    tau ~ dnorm(muTilde[xi[1]], 1)
  })
  Inits=list(xi=rep(1, 10), muTilde=rep(0,10), tau=1)
  Data=list(y=rnorm(10,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(mConf),
               'sampler_CRP: Detected unusual indexing')
    
  ## another variable depends on variable to be clustered and membership variable
  code=nimbleCode({
    for(i in 1:10) {
      muTilde[i] ~ dnorm(0, 1)  
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], 1)
      x[i] ~ dnorm(y[i] + muTilde[xi[i]], 1)
    }
    xi[1:10] ~ dCRP(1 , size=10)
  })
  Inits=list(xi=rep(1, 10), muTilde=rep(0,10), x = rep(0,10))
  Data=list(y=rnorm(10,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(mConf),
               'sampler_CRP: Cluster membership variable used in multiple declarations')
  
})

test_that("Test only data depends on cluster variable in CRP sampler", {

  ## not only data depends on xi and mu: case is safe because length of data and xi is different
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i] ~ dnorm(0, var = s2[i]/lambda)
      s2[i] ~ dinvgamma(2, 1)
      y[i] ~ dnorm(mu[xi[i]],  var = s2[xi[i]])
      x[i] ~ dnorm(mu[xi[i]], 1)
    }
    lambda ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)
    
  })
  n <- 30
  Consts <- list(n = n, x=rep(1, n))
  Inits <- list(xi = rep(1, n), 
                mu = rep(-20, n), 
                s2 = rep(0.1, n),
                alpha = 200,
                lambda = 1)
  thetas <- c(rep(-5, 10), rep(5, 10), rep(0, 10))
  Data <- list(y = rnorm(n, thetas, 1))
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 's2', 'alpha', 'lambda'))  
  expect_error(mMCMC <- buildMCMC(mConf),
               'sampler_CRP: Cluster membership variable used in multiple declarations')
  
  ## additional node depends on cluster parameters
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i] ~ dnorm(0, var = s2[i]/lambda)
      s2[i] ~ dinvgamma(2, 1)
      y[i] ~ dnorm(mu[xi[i]],  var = s2[xi[i]])
    }
    x ~ dnorm(mu[1], 1)
    lambda ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)
    
  })
  n <- 30
  Consts <- list(n = n)
  Inits <- list(xi = rep(1, n), 
                mu = rep(-20, n), 
                s2 = rep(0.1, n),
                alpha = 200,
                lambda = 1, x=1)
  thetas <- c(rep(-5, 10), rep(5, 10), rep(0, 10))
  Data <- list(y = rnorm(n, thetas, 1))
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 's2', 'alpha', 'lambda'))  
  expect_error(mMCMC <- buildMCMC(mConf),
               'sampler_CRP: Only the variables being clustered can depend on the cluster parameters')
  
})



test_that("testing multivariate normal mixture models with CRP", {
  set.seed(1)
  
  ## bivariate normal kernel with unknown mean
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
  rho0 <- 0.8
  sigmaX <- 1; sigmaY <- 1
  muX <- 10; muY <- 20
  z1 <- rnorm(n); z2 <- rnorm(n)
  Data <- list(y = cbind( sigmaX*z1+muX, sigmaY*(rho0*z1+sqrt(1-rho0^2)*z2) +muY ))
  
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 'alpha'), print=FALSE)  
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  output <- runMCMC(cMCMC, niter=2000, nburnin=1900, thin=1 , inits=Inits, setSeed=FALSE)
  outputG <- getSamplesDPmeasure(cMCMC)
  expect_false(any(is.na(outputG$samples)))
  
  xiSam <- output[, 202:301]
  muSam1 <- output[, 2:101]
  muSam2 <- output[, 102:201]
  
  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  muSam1Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam1[i, xiSam[i, ]]) )
  muSam2Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam2[i, xiSam[i, ]]) )
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam1Unique[[i]]))
  expect_equal(mean(cond), 10, tol=2*1, scale=1,
               info = paste0("incorrect update of cluster parameters in bivariate normal data with unknown mean. First component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam2Unique[[i]]))
  expect_equal(mean(cond), 20, tol=2*1, scale=1,
               info = paste0("incorrect update of cluster parameters in bivariate normal data with unknown mean. Second component"))
  
  
  ## bivariate normal kernel with unknown mean and unknown variance
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
  Inits <- list(xi = sample(1:n, size=n, replace = TRUE), 
                mu = matrix(-10, ncol=2, nrow=n),
                Sigma=Sigma,
                alpha = 1,
                lambda = 1)
  
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m, showCompilerOutput = FALSE)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 'Sigma', 'alpha', 'lambda'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  output <- runMCMC(cMCMC, niter=2000, nburn=1900, thin=1, inits=Inits, setSeed=FALSE)
  outputG <- getSamplesDPmeasure(cMCMC)
  expect_false(any(is.na(outputG$samples)))
  
  xiSam <- output[, 603:702]
  sigma11Sam <- output[, seq(1, 400, by=4)]
  muSam1 <- output[, 403:502]
  muSam2 <- output[, 503:602]
  
  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  muSam1Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam1[i, xiSam[i, ]]) )
  muSam2Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam2[i, xiSam[i, ]]) )
  sigma11Unique <- sapply(1:nrow(xiSam), function(i) unique(sigma11Sam[i, xiSam[i, ]]) )
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam1Unique[[i]]))
  expect_equal(mean(cond), 10, tol=2*1, scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. First component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam2Unique[[i]]))
  expect_equal(mean(cond), 20, tol=2*1, scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. Second component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*sigma11Unique[[i]]))
  expect_equal(mean(cond), 1, tol=2*1, scale=1,
               info = paste0("incorrect update of covariance matrix parameters in bivariate normal data with unknown mean and variance. [1, 1] component"))
  
  
  ## 4-dimensional normal kernel with unknown mean and unknown variance
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i, 1:4] ~ dmnorm(mu0[1:4], cov = S0[1:4, 1:4])
      Sigma[1:4, 1:4, i] ~ dinvwish(S = R0[1:4, 1:4], df = 4)
      SigmaAux[1:4, 1:4, i] <- Sigma[1:4, 1:4, xi[i]] / lambda  
      y[i, 1:4] ~ dmnorm(mu[xi[i], 1:4],  cov = SigmaAux[1:4, 1:4, i] )
    }
    alpha ~ dgamma(1, 1)
    lambda ~ dgamma(1, 1)
  })
  n <- 100
  Consts <- list(n = n, S0 = diag(1, 4), mu0=c(0,0,0,0), R0 = diag(1, 4))
  Sigma <- array(0, c(4,4,n))
  for(i in 1:n)
    Sigma[, , i] <- diag(1, 4)
  Inits <- list(xi = sample(1:n, size=n, replace = TRUE), 
                mu = matrix(10, ncol=4, nrow=n),
                Sigma=Sigma,
                alpha = 1,
                lambda = 1)
  z=matrix(rnorm(4*n), ncol=4, nrow=n)
  mu = c(5, 10, 15, 20)
  Sigma=diag(2, 4)
  B <- chol(Sigma)
  x <- t(mu + B%*%t(z))
  Data <- list(y = x)
  m <- nimbleModel(code, data=Data, inits=Inits, constants = Consts)
  cm <- compileNimble(m, showCompilerOutput = FALSE)
  mConf <- configureMCMC(m, monitors = c('xi','mu', 'Sigma', 'alpha', 'lambda'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  output <- runMCMC(cMCMC, niter=50, nburn=30, thin=1, inits=Inits, setSeed=FALSE)
  outputG <- getSamplesDPmeasure(cMCMC)
  expect_false(any(is.na(outputG$samples)))
  
  xiSam <- output[, grep('xi', colnames(output))]
  sigma11Sam <- output[, seq(1, n*16, by=16)]
  muSam <- output[, grep('mu', colnames(output))]
  muSam1 <- muSam[, 1:n]
  muSam2 <- muSam[, (n+1):(2*n)]
  muSam3 <- muSam[, (2*n+1):(3*n)]
  muSam4 <- muSam[, (3*n+1):(4*n)]
  
  weights <- apply(xiSam, 1, function(x) as.numeric(table(x)/sum(table(x)) ))
  muSam1Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam1[i, xiSam[i, ]]) )
  muSam2Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam2[i, xiSam[i, ]]) )
  muSam3Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam3[i, xiSam[i, ]]) )
  muSam4Unique <- sapply(1:nrow(xiSam), function(i) unique(muSam4[i, xiSam[i, ]]) )
  sigma11Unique <- sapply(1:nrow(xiSam), function(i) unique(sigma11Sam[i, xiSam[i, ]]) )
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam1Unique[[i]]))
  expect_equal(mean(cond), 5, tol=2*sqrt(2), scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. First component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam2Unique[[i]]))
  expect_equal(mean(cond), 10, tol=2*sqrt(2), scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. Second component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam3Unique[[i]]))
  expect_equal(mean(cond), 15, tol=2*sqrt(2), scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. Third component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*muSam4Unique[[i]]))
  expect_equal(mean(cond), 20, tol=2*sqrt(2), scale=1,
               info = paste0("incorrect update of mean parameters in bivariate normal data with unknown mean and variance. Fourth component"))
  
  cond <- sapply(1:nrow(xiSam), function(i) sum(weights[[i]]*sigma11Unique[[i]]))
  expect_equal(mean(cond), 2, tol=2*1, scale=1,
               info = paste0("incorrect update of covariance matrix parameters in bivariate normal data with unknown mean and variance. [1, 1] component"))
  
})



test_that("Test that not nonparametric MCMC message in CRP sampler is printed", {
  set.seed(1)
  
  ## 2 mean cluster parameters, data is mixture of 3 normals, concentration param is fixed
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:2)
        thetatilde[i] ~ dnorm(mean=0, var=10)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), 
             thetatilde=c(0,0))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m)
  expect_output(out <- runMCMC(mcmc=cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is for a parametric model.')
  
  
  ## 2 mean cluster parameters, data is mixture of 3 normals, concentration param is random
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(conc0 , size=10)
      conc0 ~ dgamma(1, 1)
      for(i in 1:2)
        thetatilde[i] ~ dnorm(mean=0, var=10)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), 
             thetatilde=c(0,0), conc0=1)
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m)
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m)
  expect_output(out <- runMCMC(mcmc=cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is not for a proper model.')
  

  ## fewer cluster parameters than observations, concentration param is fixed, conjugate case
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:5)
        thetatilde[i] ~ dnorm(mean=0, var=10)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1:5, 2), 
             thetatilde=rep(0,5))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  expect_output(out <- runMCMC(cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is for a parametric model.')
  
  
  ## fewer cluster parameters than onservation, concentration param is fixed, non conjugate case
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:5)
        thetatilde[i] ~ dt(0,1,1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1:5, 2), 
             thetatilde=rep(0,5))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  expect_output(out <- runMCMC(cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is for a parametric model.')
  
  
   
  ## no message is sent when xi=1:n
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=10)
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=1:10, 
             thetatilde=rep(0,10))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  expect_silent(out <- runMCMC(cMCMC, niter=1, nburnin = 0, thin=1))  
  
  
  ## dirichlet-multinomial model, no message is sent when xi = 1:n
  code=nimbleCode(
    {
      for(i in 1:4){
        p[i,1:3] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[xi[i],1:3], size=3)
      }
      xi[1:4] ~ dCRP(conc=1, size=4)
    }
  )
  p0 <- matrix(0, ncol=3, nrow=4)
  y0 <- matrix(0, ncol=3, nrow=4)
  for(i in 1:4){
    p0[i,]=rdirch(1, c(1, 1, 1))
    y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
  }
  m = nimbleModel(code, 
                  data = list(y = y0),
                  inits = list(xi = 1:4, p=p0), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m, monitors=c('p', 'xi'))
  mcmc <- buildMCMC(conf)
  cm = compileNimble(m)
  cmcmc=compileNimble(mcmc,project=m)
  expect_silent(cmcmc$run(100))
  
})

test_that("Check error given when model has no cluster variables", {
    ## Originally this tested whether there are no tildeVars but with new check for 'xi' appearing
    ## in non-index role, the error is caught differently.
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
               'sampler_CRP: Detected that the CRP variable is used in some way not as an index')
  
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
  
  
  ## different length of x and size:
  CRP_code2 <- nimbleCode({
    x[1:6] ~ dCRP(1, size=10)
  })
  
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model2 <- nimbleModel(CRP_code2, data=Inits)
  expect_error(CRP_model2$calculate(), "length of 'x' has to be equal to 'size'")
  
  ## different length of x and size:
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
  ## K is the number of unique components in x of length 6
  true_EK <- sum(conc/(conc+1:size-1))
  
  expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
               tol = 0.01,
               info = "Difference in expected mean of K exceeds tolerance")
  
  ## sampling from the model:
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



test_that("Testing posterior sampling and prior predictive computation with conjugate models using CRP", { 
  
  ## dnorm_dnorm
  code = nimbleCode({
    xi[1:4] ~ dCRP(conc=1, size=4)
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    
  })
  data = list(y = rnorm(4))
  inits = list(xi = rep(1,4), mu=rnorm(4))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  dataVar <- m$getParam('y[1]', 'var') 
  priorVar <- m$getParam('mu[1]', 'var')
  priorMean <- m$getParam('mu[1]', 'mean')
  postVar <- 1 / (1 / dataVar + 1 / priorVar) # from conjugate sampler
  postMean <- postVar * (data$y[1] / dataVar + priorMean / priorVar) # from conjugate sampler
  pTgivenY <- dnorm(m$mu[1] , postMean, sqrt(postVar), log = TRUE) # from conjugate sampler
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rnorm(1 , postMean, sqrt(postVar))
  expect_identical(smp, m$mu[1])
  
  
  
  ## dnorm_invgamma_dnorm
  code = nimbleCode({
    xi[1:4] ~ dCRP(conc=1, size=4)
    for(i in 1:4) {
      mu[i] ~ dnorm(0, var = s2[i]/2)
      s2[i] ~ dinvgamma(shape=2, scale=1)
      y[i] ~ dnorm(mu[xi[i]], var=s2[xi[i]])
    }
  })
  data = list(y = rnorm(4))
  inits = list(xi = rep(1,4), mu=rnorm(4), s2=rinvgamma(4, 2, 1))
  m = nimbleModel(code, data=data, inits=inits)
  conf = configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT1 <- m$getLogProb('mu[1]')
  pT2 <- m$getLogProb('s2[1]')
  
  priorMean <- m$getParam('mu[1]', 'mean')
  kappa <- values(m, 's2[1]')[1]/m$getParam('mu[1]', 'var')
  priorShape <- m$getParam('s2[1]', 'shape')
  priorScale <- m$getParam('s2[1]',  'scale')
  pTgivenY2 <- dinvgamma(m$s2[1], shape = priorShape + 1/2,
                         scale = priorScale + kappa * (data$y[1] - priorMean)^2 / (2*(1+kappa)),
                         log=TRUE)
  pTgivenY1 <- dnorm(m$mu[1], mean = (kappa * priorMean + data$y[1])/(1 + kappa), 
                     sd = sqrt(m$s2[1] / (1+kappa)),
                     log=TRUE) 
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT1 + pT2 + pYgivenT - pTgivenY1 - pTgivenY2)
  
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- rinvgamma(1, shape = priorShape + 1/2,
                    scale = priorScale + kappa * (data$y[1] - priorMean)^2 / (2*(1+kappa)) )
  smp2 <- rnorm(1, mean = (kappa * priorMean + data$y[1])/(1 + kappa), 
                sd = sqrt(smp1 / (1+kappa))) 
  expect_identical(smp1, m$s2[1])
  expect_identical(smp2, m$mu[1])
  
  
  ## dgamma_dpois
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dpois(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rpois(4, 4))
  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape <- m$getParam('mu[1]', 'shape')
  priorRate <- m$getParam('mu[1]', 'rate')
  pTgivenY <- dgamma(m$mu[1], shape = priorShape + data$y[1], rate = priorRate + 1, log=TRUE)

  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rgamma(1 , shape = priorShape + data$y[1], rate = priorRate + 1)
  expect_identical(smp, m$mu[1])
  
  
  
  ## dbeta_dbern
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dbeta(1,1)
      y[i] ~ dbern(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rbinom(4, size=1, prob=0.5))
  inits = list(xi = rep(1,4), mu=rbeta(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape1 <- m$getParam('mu[1]', 'shape1')
  priorShape2 <- m$getParam('mu[1]', 'shape2')
  pTgivenY <- dbeta(m$mu[1], shape1=priorShape1+data$y[1], shape2=priorShape2+1-data$y[1], log=TRUE)

  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rbeta(1 , shape1=priorShape1+data$y[1], shape2=priorShape2+1-data$y[1])
  expect_identical(smp, m$mu[1])
  
  
  ## dbeta_dbin
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dbeta(1,1)
      y[i] ~ dbinom(size=10, prob=mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rbinom(4, size=10, prob=0.5))
  inits = list(xi = rep(1,4), mu=rbeta(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape1 <- m$getParam('mu[1]', 'shape1')
  priorShape2 <- m$getParam('mu[1]', 'shape2')
  dataSize <- m$getParam('y[1]', 'size')
  pTgivenY <- dbeta(m$mu[1], shape1=priorShape1+data$y[1], shape2=priorShape2+dataSize-data$y[1], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rbeta(1 , shape1=priorShape1+data$y[1], shape2=priorShape2+dataSize-data$y[1])
  expect_identical(smp, m$mu[1])
  
  
  ## dbeta_dnegbin
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dbeta(1,1)
      y[i] ~ dnegbin(size=10, prob=mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rnbinom(4, size=10, prob=0.5))
  inits = list(xi = rep(1,4), mu=rbeta(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape1 <- m$getParam('mu[1]', 'shape1')
  priorShape2 <- m$getParam('mu[1]', 'shape2')
  dataSize <- m$getParam('y[1]', 'size')
  pTgivenY <- dbeta(m$mu[1], shape1=priorShape1+dataSize, shape2=priorShape2+data$y[1], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rbeta(1 , shape1=priorShape1+dataSize, shape2=priorShape2+data$y[1])
  expect_identical(smp, m$mu[1])
  
  
  ## dgamma_dexp:
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dexp(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rexp(4, 4))
  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape <- m$getParam('mu[1]', 'shape')
  priorRate <- m$getParam('mu[1]', 'rate')
  pTgivenY <- dgamma(m$mu[1], shape=priorShape+1, rate=priorRate+data$y[1], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rgamma(1, shape=priorShape+1, rate=priorRate+data$y[1])
  expect_identical(smp, m$mu[1])
  
  
  
  ## dgamma_dgamma:
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dgamma(4, mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rgamma(4, 4, 4))
  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape <- m$getParam('mu[1]', 'shape')
  priorRate <- m$getParam('mu[1]', 'rate')
  dataShape <- m$getParam('y[1]', 'shape')
  pTgivenY <- dgamma(m$mu[1], shape=dataShape+priorShape, rate=priorRate+data$y[1], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rgamma(1, shape=dataShape+priorShape, rate=priorRate+data$y[1])
  expect_identical(smp, m$mu[1])
  
  
  ## dgamma_dweib:
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dweib(shape=4, lambda = mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rweibull(4, 4, 4))
  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape <- m$getParam('mu[1]', 'shape')
  priorRate <- m$getParam('mu[1]', 'rate')
  dataShape <- m$getParam('y[1]', 'shape')
  pTgivenY <- dgamma(m$mu[1], shape=1+priorShape, rate=priorRate+data$y[1]^dataShape, log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rgamma(1, shape=1+priorShape, rate=priorRate+data$y[1]^dataShape)
  expect_identical(smp, m$mu[1])
  
  
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
  data = list(y = y0)
  inits = list(xi = rep(1,4), p=p0)
  m = nimbleModel(code, data=data, inits=inits,
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m)
  mcmc <- buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('p[1,1:3]')
  
  priorAlpha <- m$getParam('p[1, 1:3]', 'alpha')
  pTgivenY <- ddirch(m$p[1,1:3], alpha = priorAlpha+data$y[1, 1:3], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rdirch(1, alpha = priorAlpha+data$y[1, 1:3])
  expect_identical(smp, m$p[1, 1:3])
  
  
  ## dgamma_dinvgamma:
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1, rate=1)
      y[i] ~ dinvgamma(shape=4, scale = mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rinvgamma(4, 4, 4))
  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]')
  pT <- m$getLogProb('mu[1]')
  
  priorShape <- m$getParam('mu[1]', 'shape')
  priorRate <- m$getParam('mu[1]', 'rate')
  dataShape <- m$getParam('y[1]', 'shape')
  pTgivenY <- dgamma(m$mu[1], shape=dataShape+priorShape, rate=priorRate+1/data$y[1], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rgamma(1, shape=dataShape+priorShape, rate=priorRate+1/data$y[1])
  expect_identical(smp, m$mu[1])
  
  
  ## dgamma_dnorm:
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dnorm(0, tau = mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  data = list(y = rnorm(4, 0, 4))
  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
  m = nimbleModel(code, data=data, inits= inits)
  conf = configureMCMC(m)
  mcmc = buildMCMC(conf)
  
  pYgivenT <- m$getLogProb('y[1]') ; dnorm(data$y[1], 0, sqrt(1/m$mu[1]), log=TRUE)
  pT <- m$getLogProb('mu[1]') ; dgamma(m$mu[1], 1, 1,log=TRUE)
  
  dataMean <- m$getParam('y[1]', 'mean')
  priorShape <- m$getParam('mu[1]', 'shape')
  priorRate <- m$getParam('mu[1]', 'rate')
  pTgivenY <- dgamma(m$mu[1], shape = priorShape + 0.5, rate = (priorRate + 0.5*(data$y[1]-dataMean)^2), log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rgamma(1 , shape = priorShape + 0.5, rate = priorRate + (data$y[1]-dataMean)^2/2)
  expect_identical(smp, m$mu[1])
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm with truncation
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    for(i in 1:2){
      mu[i] ~ dnorm(0,1)}
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "sampler_CRP: The number of cluster parameters is less")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm and deterministic nodes and truncation
  code = nimbleCode({
    for(i in 1:4) {
      mui[i] <- mu[xi[i]]
      y[i] ~ dnorm(mui[i], sd = 1)
    }
    for(i in 1:2){
      mu[i] ~ dnorm(0,1)}
    xi[1:4] ~ dCRP(conc=1, size=4) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm and non-standard indexing
  code = nimbleCode({
    for(i in 1:4) {
      mu[i, 2] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i], 2], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=cbind(rnorm(4),rnorm(4))))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm and non-standard indexing
  code = nimbleCode({
    for(i in 1:4) {
      mu[2, i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[2, xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=t(cbind(rnorm(4),rnorm(4)))))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dexp")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dgamma")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dnorm")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dweib")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dinvgamma")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dbern")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dbin")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dnegbin")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc <- buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  ## non-standard ordering/indexing of ddirch-dmulti
  code=nimbleCode(
    {
      for(i in 1:4){
        p[1:3, i] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[1:3, xi[i]], size=3)
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
                  inits = list(xi = rep(1,4), p=t(p0)), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc <- buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  code=nimbleCode(
    {
      for(i in 1:4){
        p[i, 2:4] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[xi[i], 2:4], size=3)
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
  p0 <- cbind(rep(0, 4), p0)
  m = nimbleModel(code, 
                  data = list(y = y0),
                  inits = list(xi = rep(1,4), p=p0), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc <- buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  ## dnorm, dinvgamma, not conjugate
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## dnorm, dinvgamma, not conjugate
  code = nimbleCode({
    for(i in 1:4) {
      sigma[i] ~ dinvgamma(1, 1)
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = sigma[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), sigma = rinvgamma(4, 1,1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## dnorm_invgamma, conjugate; we detect conjugacy
  code = nimbleCode({
    for(i in 1:4) {
      s2[i] ~ dinvgamma(a,b)
      mu[i] ~ dnorm(0, var = s2[i]/kappa)
      y[i] ~ dnorm(mu[xi[i]], var = s2[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
    kappa ~ dgamma(1, 1)
    a ~ dgamma(1, 1)
    b ~ dgamma(1, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4), s2=rinvgamma(4, 1,1), a=1, b=1, kappa=2))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(nimble:::checkCRPconjugacy(m, 'xi[1:4]'), "conjugate_dnorm_invgamma_dnorm")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  
  ## model with deterministic nodes
  code = nimbleCode({
    for(i in 1:4) {
      s2Tilde[i] ~ dinvgamma(a,b)
      s2[i] <- s2Tilde[xi[i]]
      muTilde[i] ~ dnorm(0, var = s2Tilde[i]/kappa)
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], var = s2[i])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
    kappa ~ dgamma(1, 1)
    a ~ dgamma(1, 1)
    b ~ dgamma(1, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), muTilde=rnorm(4), s2Tilde=rinvgamma(4, 1,1), a=1, b=1, kappa=2))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(nimble:::checkCRPconjugacy(m, 'xi[1:4]'), "conjugate_dnorm_invgamma_dnorm")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dexp")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
})


test_that("Testing handling (including error detection) with non-standard CRP model specification",{
 
  n <- 20
  const <- list(n = n)
  inits <- list(xi = rep(1,n), muTilde = rnorm(n), conc = 1)
  data <- list(y = rnorm(n))
  tildeNames <- paste0("muTilde[", 1:n, "]")
  target <- paste0("xi[1:", n, "]")
  
  ## basic model to check results of findClusterNodes()
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n){
      muTilde[i] ~ dnorm(0,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## more complicated indexing to check results of findClusterNodes()
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i], 2]
    }
    for(i in 1:n){muTilde[i, 2] ~ dnorm(0,1)}
    
  })
  inits2 <- inits
  inits2$muTilde <- matrix(rnorm(n*2), n)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, ", 2]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(2, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## fewer tildeNodes than observations
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:(n-2)){muTilde[i] ~ dnorm(0,1)}
    
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "less than the number of potential clusters")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n-2, clusterNodeInfo$nTilde)
  
  ## indirect indexing; we don't have a good way to handle this.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[b[i]]
    }
    for(j in 1:n)
    {b[j] <- xi[j]}
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Detected that the CRP variable is used in some way not as an index")
  
  ## cluster ID as second index
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(muTilde[2, xi[i]], var = 1)
    }
    for(i in 1:n)
    {muTilde[2, i] ~ dnorm(0,1)}
  })
  inits2 <- inits
  inits2$muTilde <- rbind(rnorm(n), rnorm(n))
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[2, ", 1:n, "]"))
  expect_equal(2, clusterNodeInfo$numIndexes)
  expect_equal(2, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  
  ## clusterID as second index, additional nodes that are not clusterNodes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(muTilde[2, xi[i]], var = 1)
    }
    for(j in 1:2)
      for(i in 1:n)
      {muTilde[j, i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[2, ", 1:n, "]"))
  expect_equal(2, clusterNodeInfo$numIndexes)
  expect_equal(2, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  
  ## cluster ID in a function
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+1]
    }
    for(i in 1:(n+1))
    {muTilde[i] ~ dnorm(0,1)}
  })
  inits2$muTilde <- rnorm(n+1)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n+1), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## cluster ID in a function, no extraneous muTilde
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+1]
    }
    for(i in 2:(n+1))
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n+1), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## reordering of muTildes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[n-xi[i]+1]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", n:1, "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## function in index (not allowed)
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[n-i+1]]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "not designed for this case")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(TRUE, clusterNodeInfo$targetIndexedByFunction)
  
  ## clusterNodes indexing doesn't begin at 1
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+2]
    }
    for(i in 3:(n+2))
    {muTilde[i] ~ dnorm(0,1)}
  })
  inits2$muTilde <- rnorm(n+2)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 3:(n+2), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## Extra nodes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:(2*n))
    {muTilde[i] ~ dnorm(0,1)}
  })
  inits2$muTilde <- rnorm(2*n)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## missing first cluster node
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 2:(n-2))
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  expect_warning(conf <- configureMCMC(m), "missing cluster parameter")
  expect_warning(mcmc <- buildMCMC(conf), "missing cluster parameter")
  expect_warning(clusterNodeInfo <- nimble:::findClusterNodes(m, target), "missing cluster parameter")
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n-2), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## cluster node indexing shifted
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+1]
    }
    for(i in 2:(n-2))
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "less than the number of potential clusters")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n-2), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n-3, clusterNodeInfo$nTilde)
  
  ## extra dependency on cluster nodes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0, 1)}
    z ~ dnorm(muTilde[1], 1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Only the variables being clustered")
  
  ## multiple observations per cluster membership; not yet handled
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(j in 1:2) {
      for(i in 1:n) {
        y[i,j] ~ dnorm(mu[i], var = 1)
      }}
    for(i in 1:n)
    {mu[i] <- muTilde[xi[i]]}
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = list(y = matrix(rnorm(2*n),n)), constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "when there is one variable being clustered")
  
  ## Extraneous node that is ok.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
    z ~ dnorm(muTilde[n+1], 1)
  })
  inits2$muTilde <- rnorm(n+1)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## Awkward trapping of observations with different distributions.
  ## We should trap this when checking conjugacy instead and simply assign non-conjugate sampler.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) 
      y[i] ~ dnorm(mu[i], var = 1)
    for(i in 1:(n-1))
    {mu[i] <- muTilde[xi[i]]}
    mu[n] <- exp(muTilde[xi[n]])
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0, 1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Detected unusual indexing in")
  
  
  ## conjugate but observations not IID  
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2[i])
      mu[i] <- muTilde[xi[i]]
      s2[i] ~ dgamma(1,1)
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## conjugacy not detected because observations have multiple declarations
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:(n/2)) {
      y[i] ~ dnorm(mu[i], var = 1)
    }
    for(i in ((n/2)+1):n)
    {y[i] ~ dnorm(mu[i], var = 1)}
    for(i in 1:n)
    {mu[i] <- muTilde[xi[i]]}
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## observations not independent
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    y[1] ~ dnorm(mu[1], var = s2[1])
    for(i in 2:n)
      y[i] ~ dnorm(mu[i]+y[i-1], var = s2[i])
    for(i in 1:n) {
      mu[i] <- muTilde[xi[i]]
      s2[i] ~ dgamma(1,1)
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Variables being clustered must be conditionally independent.")
  
  ## cluster nodes not independent
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    muTilde[1] ~ dnorm(0, 1)
    for(i in 2:n)
    {muTilde[i] ~ dnorm(muTilde[i-1],1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Cluster parameters must be conditionally independent.")
  
  ## cluster nodes not exchangeable so non-conjugate
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:(n-1) )
    {muTilde[i] ~ dnorm(mu0[i],1)}
    muTilde[n] ~ dgamma(1,1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## cluster membership variables not independent of cluster parameters
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc + muTilde[1], n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
    
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Only the variables being clustered can depend")
  
  ## cluster membership variables not independent of cluster parameters
  ## This is not detected by the check for independence but by use of 'xi' as non-index
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
      tmp[i] ~ dnorm(0,1)
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(tmp[xi[i]],1)}
    
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Only the variables being clustered can depend")
  
  inits$s2Tilde <- rep(1, n)
  
  ## Conjugate normal-invgamma
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0, var = s2Tilde[i])
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:n, "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", 1:n, "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(rep(TRUE, 2), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(rep(n,2), clusterNodeInfo$nTilde)
  
  ## nTilde < n for one of the cluster parameters. Note confusing error message.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i]]
    }
    kappa ~ dgamma(1,1)
    for(i in 1:n) 
      {muTilde[i] ~ dnorm(0, var = s2Tilde[i]/kappa)}
    for(i in 1:(n-1))
      {s2Tilde[i] ~ dinvgamma(1,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Cluster parameters must be conditionally independent")
  
  ## nTilde < n 
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i]]
    }
    kappa ~ dgamma(1,1)
    for(i in 1:(n-1)) {
      muTilde[i] ~ dnorm(0,var = s2Tilde[i]/kappa)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "less than the number of potential clusters")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:(n-1), "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", 1:(n-1), "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(rep(TRUE, 2), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(c(n-1, n-1), clusterNodeInfo$nTilde)
  
  ## CRP variable used in multiple indices; disallowing this.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i],xi[i]]
    }
    for(i in 1:n) {
      for(j in 1:n)
        {muTilde[i,j] ~ dnorm(0,var=s2Tilde[i]/3)}
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  inits2 <- inits
  inits2$muTilde <- matrix(rnorm(n^2),n)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  expect_error(conf <- configureMCMC(m), "CRP variable used multiple times")
  
  ## weird ordering of muTilde/s2Tilde but should be ok
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[n-xi[i]+1])
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,var=s2Tilde[n-i+1]/3)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:n, "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", n:1, "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(c(FALSE, TRUE), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(rep(n,2), clusterNodeInfo$nTilde)
  
  ## s2Tildes in different order than muTildes so not conjugate.
  ## CRP_sampler is INCORRECT for this because can't sample from distr of an s2Tilde given the muTilde that depends on it.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[n-xi[i]+1])
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,var=s2Tilde[i]/3)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Cluster parameters must be conditionally independent")
  
  ## Model is valid, but in trying to catch weird uses of CRP variable we don't allow this. should be ok
  inits2 <- inits
  inits2$muTilde <- cbind(rnorm(n), rgamma(n, 1, 1))
  inits2$s2Tilde <- NULL
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(muTilde[xi[i],1], var = muTilde[xi[i],2])
    }
    for(i in 1:n) {
      muTilde[i,1] ~ dnorm(0,var=s2Tilde[i]/3)
      muTilde[i,2] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Cluster membership variable used in multiple declarations")
  
  ## Non-conjugate, bivariate
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(muTilde[xi[i]], var = exp(s2Tilde[xi[i]]))
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,1)
      s2Tilde[i] ~ dgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:n, "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", 1:n, "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(rep(TRUE, 2), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(rep(n,2), clusterNodeInfo$nTilde)
  
  ## cross clustering
  
  data$y <- matrix(rnorm(n^2), n)
  
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      for(j in 1:n)
        {y[i,j] ~ dnorm(muTilde[xi[i]], var = s2Tilde[xi[j]])}
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,1)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "NIMBLE can only sample when there is one variable being clustered")
  
  inits$muTilde <- matrix(rnorm(n^2), n)
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      for(j in 1:n)
        {y[i,j] ~ dnorm(muTilde[xi[i],xi[j]], var = s2Tilde[xi[i]])}
    }
    for(i in 1:n)  {
      for(j in 1:n)
        {muTilde[i,j] ~ dnorm(0,1)}
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  expect_error(conf <- configureMCMC(m), "CRP variable used multiple times in")
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
  expect_error(buildMCMC(conf), "NIMBLE can only sample when there is one variable")
  
  
  ## more observations than labels 
  code = nimbleCode({
    for(i in 1:10) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  expect_error(nimbleModel(code, data = list(y = rnorm(10)),
                           inits = list(xi = rep(1,4), mu=rnorm(10))),
               "dimensions specified are smaller")
  
  
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
  expect_error(buildMCMC(conf), "sampler_CRP: Detected unusual indexing")
  
  
  ## same obervation with different label
  code = nimbleCode({
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
    y[1] ~ dnorm(mu[xi[1]], 1)
    y[1] ~ dnorm(mu[xi[2]], 1)
    xi[1:2] ~ dCRP(conc=1, size=2)
  })
  expect_error(nimbleModel(code, data = list(y = rnorm(2)),
                           inits = list(xi = rep(1,2), mu=rnorm(2))),
               "There are multiple definitions")
  
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
  expect_warning(buildMCMC(conf),
                 "The number of cluster parameters is less than the number of potential clusters")
  
  
  
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
  expect_warning(buildMCMC(conf),
                 "The number of cluster parameters is less than the number of potential clusters")
  
  
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
                   inits = list(xi = rep(1,100), mu=rnorm(50), s2=rinvgamma(1,1,1)))
  conf <- configureMCMC(m)
  expect_error(buildMCMC(conf), "Detected unusual indexing in")
  
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
  expect_error(buildMCMC(conf), "Detected unusual indexing in")
  
  ## test that a message is sent when truncation is hit
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
  expect_output(cmMCMC$run(1), "CRP_sampler: This MCMC is for a parametric model")
  
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
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mMCMC <- buildMCMC(mConf)
  expect_equal(class(mMCMC$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  CmMCMC <- compileNimble(mMCMC, project=m)
  samples <- runMCMC(CmMCMC, niter=600, nburnin=500)
  
  samplesG <- getSamplesDPmeasure(CmMCMC)
  trunc <- samplesG$trunc
  samplesG <- samplesG$samples
  
  ygrid <- seq(0, 100, by=1)
  fSam <- matrix(0, ncol=length(ygrid), nrow=nrow(samplesG))
  for(i in 1:nrow(samplesG)){
    fSam[i, ] <- sapply(ygrid, function(x)sum(samplesG[i, 1:trunc]*dpois(x, samplesG[i, (trunc+1):(2*trunc)])))
  }
  fHat <- apply(fSam, 2, mean)
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
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mMCMC <- buildMCMC(mConf)
  concIndex <- which(sapply(mConf$getSamplers(), function(x) x[['target']]) == 'alpha')
  expect_equal(mConf$getSamplers()[[concIndex]]$name, "CRP_concentration")
  expect_equal(class(mMCMC$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  CmMCMC <- compileNimble(mMCMC, project=m)
  samples <- runMCMC(CmMCMC, niter=600, nburnin=500)
  
  samplesG <- getSamplesDPmeasure(CmMCMC)
  trunc <- samplesG$trunc
  samplesG <- samplesG$samples
  
  ygrid <- seq(0, 100, by=1)
  fSam <- matrix(0, ncol=length(ygrid), nrow=nrow(samplesG))
  for(i in 1:nrow(samplesG)){
    fSam[i, ] <- sapply(ygrid, function(x)sum(samplesG[i, 1:trunc]*dpois(x, samplesG[i, (trunc+1):(2*trunc)])))
  }
  fHat <- apply(fSam, 2, mean)
  f0 <- function(x) 0.2*dpois(x, 10) + 0.2*dpois(x, 5) + 0.6*dpois(x, 50)
  f0grid <- sapply(ygrid, f0)
  
  L1dist <- mean(abs(f0grid - fHat))
  
  expect_equal(L1dist, 0.01, tol=0.01,
               info = "wrong estimation of density in DPM of Poisson distrbutions with random concentration param")
  
  
  ## DPM of normal densities with random means and variances
  Code=nimbleCode(
    {
      for(i in 1:50){
        thetatilde[i] ~ dnorm(mean=0, var=100) 
        s2tilde[i] ~ dinvgamma(shape=2, scale=1) 
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
                thetatilde = rnorm(50, 0, 1),
                s2tilde = rinvgamma(50, shape=2, scale=1))
  Data <- list(y = c(rnorm(50, 5,sqrt(4)), rnorm(50, -5, sqrt(2))))
  
  m <- nimbleModel(Code, data=Data, inits=Inits, calculate=TRUE)
  cmodel <- compileNimble(m)
  
  mConf <- configureMCMC(m, monitors = c('xi', 'thetatilde', 's2tilde'))
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mMCMC <- buildMCMC(mConf))
  expect_equal(class(mMCMC$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  CmMCMC=compileNimble(mMCMC, project=m)
  samples <- runMCMC(CmMCMC, niter=2000, nburnin=1900)
  
  samplesG <- getSamplesDPmeasure(CmMCMC)
  trunc <- samplesG$trunc
  samplesG <- samplesG$samples
  
  ygrid <- seq(-10, 10, len=100)
  fSam <- matrix(0, ncol=length(ygrid), nrow=nrow(samplesG))
  for(i in 1:nrow(samplesG)){
    fSam[i, ] <- sapply(ygrid, function(x)sum(samplesG[i, 1:trunc]*
                                                dnorm(x,samplesG[i, (2*trunc+1):(3*trunc)],sqrt(samplesG[i, (trunc+1):(2*trunc)]))))
  }
  fHat <- apply(fSam, 2, mean)
  f0 <- function(x) 0.5*dnorm(x, 5,sqrt(4)) + 0.5*dnorm(x, -5,sqrt(2))
  f0grid <- sapply(ygrid, f0)
  
  L1dist <- mean(abs(f0grid - fHat))
  
  expect_equal(L1dist, 0.01, tol=0.01,
               info = "wrong estimation of density in DPM of normal distrbutions")
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
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mMCMC <- buildMCMC(mConf)
  
  expect_equal(mConf$getSamplers()[[1]]$name, "CRP_concentration")
  expect_equal(class(mMCMC$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
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
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_identical(length(crpIndex), 1L)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
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
  
  modelConf = configureMCMC(model, thin=100)
  
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
  
  ## posterior samples of the density
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
  
  
  ## wrong specification of stick variables
  SB_code3 <- nimbleCode({
    for(i in 1:5) 
      z[i] ~ dgamma(10, 10)
    w[1:6] <- stick_breaking(z[1:5])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  ## wrong prior and starting values for stick variables
  set.seed(1)
  Inits <- list(z = rgamma(5, 10, 10))
  data <- list(xi = 1:10)
  expect_output(m <- nimbleModel(SB_code3, data=data, inits=Inits),
                "values in 'z' have to be in \\(0,1\\)")
  
  ## good starting values for stick variables
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  data <- list(xi = rep(1,10))
  SB_model3 <- nimbleModel(SB_code3, data=data, inits=Inits)
  expect_output(m <- SB_model3$simulate(), "values in 'z' have to be in \\(0,1\\)")
  
  ## wrong specification of length in stick variables, should be 5
  SB_code4 <- nimbleCode({
    for(i in 1:4) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:4])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  ## wrong length for stick variables
  set.seed(1)
  Inits <- list(z = rbeta(4, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(m <- nimbleModel(SB_code4, data=data, inits=Inits),
                 "number of items to replace")
  
  
  ## wrong specification of length in stick variables, should be 5
  ## no warning in nimbleModel function
  ## This is a shortcoming in NIMBLE processing of sizes in models.
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
  expect_failure(expect_error(SB_model5 <- nimbleModel(SB_code5, data=data, inits=Inits)))
  cSB_model5 <- compileNimble(SB_model5)
  expect_output(cSB_model5$calculate('w'), "Error in mapCopy")
  
  ## longer vector of stick variables
  SB_code6 <- nimbleCode({
    for(i in 1:10) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:10])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  ## wrong length for stick variables, a warning is sent in the nimbleModel function
  set.seed(1)
  Inits <- list(z = rbeta(10, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(SB_model6 <- nimbleModel(SB_code6, data=data, inits=Inits),
                 "number of items to replace")
  cSB_model6 <- compileNimble(SB_model6)
  expect_output(cSB_model6$calculate('w'), "Error in mapCopy")
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
  expect_equal(conf$getSamplers('alpha')[[1]]$name, "CRP_concentration")
  
  
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
  expect_equal(conf$getSamplers('alpha')[[1]]$name, "RW")
  
  
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
  expect_equal(conf$getSamplers('alpha')[[1]]$name, "RW")
  
  
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

test_that("Testing that cluster parameters are appropriately updated and mvSaved in good state", {
    ## Should always reject new clusters
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ T(dnorm(50, 1), -200, 200)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu = rep(0, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_identical(cm$mu[2], 0, info = 'mu[2] is changed')
    expect_identical(cm$mu[3], 0, info = 'mu[3] is changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())
    
    ## Should accept a cluster
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ T(dnorm(0, 1), -200, 200)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = c(0, rnorm(n-1, 50, 1)))
    inits <- list(xi = rep(1,n), mu = rep(50, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_equal(cm$mu[2], 0, tolerance = 3, info = 'mu[2] is not changed')
    expect_equal(cm$mu[3], 0, tolerance = 3, info = 'mu[3] is not changed')
    expect_identical(cm$mu[4], 50, info = 'mu[4] is changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())

    ## Should always reject new clusters
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(50, 1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu = rep(0, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_identical(cm$mu[2], 0, info = 'mu[2] is changed')
    expect_identical(cm$mu[3], 0, info = 'mu[3] is changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())
    
    ## Should accept a cluster
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0, 1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = c(0, rnorm(n-1, 50, 1)))
    inits <- list(xi = rep(1,n), mu = rep(50, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_equal(cm$mu[2], 0, tolerance = 3, info = 'mu[2] is not changed')
    expect_identical(cm$mu[3], 50, info = 'mu[3] is  changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())

    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            muTilde[i] ~ T(dnorm(0, 1), -200, 200)
            sigmaTilde[i] ~ dinvgamma(a,b)
            mu[i] <- muTilde[xi[i]]
            y[i] ~ dnorm(mu[i], sd = sigmaTilde[xi[i]])
        }
        a ~ dgamma(1,1)
        b ~ dgamma(1,1)
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), muTilde = rnorm(n), sigmaTilde = rinvgamma(n, 1, 1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(10)
    expect_identical(cmcmc$mvSaved[['muTilde']], cm$muTilde)
    expect_identical(cmcmc$mvSaved[['logProb_muTilde']], cm$logProb_muTilde)
    expect_identical(cmcmc$mvSaved[['sigmaTilde']], cm$sigmaTilde)
    expect_identical(cmcmc$mvSaved[['logProb_sigmaTilde']], cm$logProb_sigmaTilde)
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(cmcmc$mvSaved[['y']], cm$y)
    expect_identical(cmcmc$mvSaved[['logProb_y']], cm$logProb_y)
    expect_identical(cmcmc$mvSaved[['a']], cm$a)
    expect_identical(cmcmc$mvSaved[['logProb_a']], cm$logProb_a)
    expect_identical(cmcmc$mvSaved[['b']], cm$b)
    expect_identical(cmcmc$mvSaved[['logProb_b']], cm$logProb_b)
  
})
  
test_that("Testing wrapper sampler that avoids sampling empty clusters", {
    set.seed(1)
    code = nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rnorm(n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[2]]$control$wrapped_type, 'conjugate_dnorm_dnorm_dynamicDeps',
                     info = "cluster wrapper sampler not conjugate")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rgamma(n,1,1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[2]]$control$wrapped_type, 'RW',
                     info = "cluster wrapper sampler conjugate")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')


    ## Check that changes in cluster parameters persist even if no direct sampling of them
    
    set.seed(1)
    code = nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rnorm(n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    conf$removeSamplers('mu')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rgamma(n,1,1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    conf$removeSamplers('mu')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    ## p=2 non-conjugate case
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0,1)
            sigma[i] ~ dinvgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], var = sigma[xi[i]])
        }
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu = rnorm(n), sigma = rinvgamma(n, 1, 1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[17]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (2*n+1):(3*n)])-1  # 2nd to last so have a few transitions
    focalClusterPresent <- apply(out[ , (2*n+1):(3*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    
    focalClusterName <- paste0('mu[', focalCluster, ']')
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')
    focalClusterName <- paste0('sigma[', focalCluster, ']')
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    ## ensure resetting cluster param values works for mv cluster parameter
    set.seed(1)
    code = nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i, 1:3] ~ dmnorm(z[1:3], pr[1:3,1:3])
            tmp[i, 1:3] <- exp(mu[xi[i],1:3])
            y[i,1:3] ~ dmnorm(tmp[i,1:3], pr[1:3,1:3])
        }
        
    })
    n <- 15
    data <- list(y = matrix(rnorm(n*3, 1, 1),n))
    inits <- list(xi = rep(1,n), mu = matrix(rnorm(n*3), n), pr = diag(3), z = rep(0,3))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[2]]$control$wrapped_type, 'RW_block',
                     info = "cluster wrapper sampler conjugate")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ', 3]')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (3*n+1):(4*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ', 3]')
    focalClusterPresent <- apply(out[ , (3*n+1):(4*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    ## ensure that two different kinds of wrapped samplers compile and run
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        xi2[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu2[i] ~ dnorm(0,1)
            y2[i] ~ dnorm(mu2[xi2[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n), y2 = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rgamma(n,1,1), xi2 = rep(1,n), mu2 = rnorm(n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    ## It's hard to get testthat checking of logging to work right, so just
    ## running so that testthat catches if code errors out.
    cmcmc <- compileNimble(mcmc,project=m)
    out <- runMCMC(cmcmc, 10)
    
    
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

