# litters model based on original BUGS example (not the example in
# the JAGS version of the BUGS examples or included in the
# NIMBLE package)
# in particular no use of T(,)


modelCode <- function(){
  for (i in 1:G) {
     for (j in 1:N) {
        r[i,j] ~ dbin(p[i,j], n[i,j]);
        p[i,j] ~ dbeta(a[i], b[i]) 
     }
     mu[i] <- a[i] / (a[i] + b[i]);
     theta[i] <- 1 / (a[i] + b[i]);
     a[i] ~ dgamma(1, .001)
     b[i] ~ dgamma(1, .001)
   }
}

modelData <- list(
               G = 2,
               n = structure(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13, 7, 5, 10, 7, 6, 10, 10, 10, 7), .Dim = c(2, 16)),
               N = 16,
               r = structure(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9, 12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0), .Dim = c(2, 16))
         )

Rmodel <- readBUGSmodel(modelCode, data = modelData)

Cmodel <- compileNimble(Rmodel)
mcmcspec <- MCMCspec(Rmodel)
mcmcspec$addMonitors(c('a', 'b'))

# first MCMC: default specification

Rmcmc <- buildMCMC(mcmcspec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


niter <- 20
set.seed(0);     Rmcmc(niter)
set.seed(0);     Cmcmc(niter)

Rsamples  <- as.matrix(nfVar(Rmcmc, 'mvSamples'))
Csamples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))

if(require(testthat)) {
  context(paste0("testing litters MCMC"))
  test_that("test of equality of output from R and C versions of litters MCMC", {
    expect_that(Rsamples, equals(Csamples), info = paste("R and C posterior samples are not equal"))
  })
}

# second MCMC: use slice samplers at top-level

Rmodel2 <- readBUGSmodel(modelCode, data = modelData)

mcmcspec <- MCMCspec(Rmodel2)
for(node in c('a[1]', 'b[1]', 'a[2]', 'b[2]')) 
  mcmcspec$addSampler('slice', list(targetNode = node , adaptInterval = 100))
# remove MH for a,b

Rmcmc <- buildMCMC(mcmcspec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel, reset = TRUE)


niter <- 20
set.seed(0);     Rmcmc(niter)
set.seed(0);     Cmcmc(niter)

Rsamples  <- as.matrix(nfVar(Rmcmc, 'mvSamples'))
Csamples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))

if(require(testthat)) {
  context(paste0("testing litters MCMC"))
  test_that("test of equality of output from R and C versions of litters MCMC", {
    expect_that(Rsamples, equals(Csamples), info = paste("R and C posterior samples are not equal"))
  })
}

# third MCMC: use block samplers at top-level
mcmcspec$addSampler('RW_block', list(targetNodes = c('a[1]', 'b[1]') , adaptInterval = 100))
mcmcspec$addSampler('RW_block', list(targetNodes = c('a[2]', 'b[2]') , adaptInterval = 100))

Rmcmc <- buildMCMC(mcmcspec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 20
set.seed(0);     Rmcmc(niter)
set.seed(0);     Cmcmc(niter)

Rsamples  <- as.matrix(nfVar(Rmcmc, 'mvSamples'))
Csamples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))

if(require(testthat)) {
  context(paste0("testing litters MCMC"))
  test_that("test of equality of output from R and C versions of litters MCMC", {
    expect_that(Rsamples, equals(Csamples), info = paste("R and C posterior samples are not equal"))
  })
}

# fourth MCMC: cross-level (not debugged yet)
mcmcspec$addSampler('crossLevel', control = list(topNodes = c('a[1]', 'b[1]'), adaptInterval = 100))
mcmcspec$addSampler('crossLevel', control = list(topNodes = c('a[2]', 'b[2]'), adaptInterval = 100))

Rmcmc <- buildMCMC(mcmcspec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 10
set.seed(0);     Rmcmc(niter)
set.seed(0);     Cmcmc(niter)

Rsamples  <- as.matrix(nfVar(Rmcmc, 'mvSamples'))
Csamples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))

if(require(testthat)) {
  context(paste0("testing litters MCMC"))
  test_that("test of equality of output from R and C versions of litters MCMC", {
    expect_that(Rsamples, equals(Csamples), info = paste("R and C posterior samples are not equal"))
  })
}


### OLD MCMCski comparison with JAGS
if(FALSE) {
  jagsTime = system.time({out1 <- jags(data = modelData,
    parameters.to.save = c('p','mu','a','b','theta'), n.chains = 1,
    n.iter = 12000, n.burnin = 2000, n.thin = 1, model.file = modelCode, DIC = FALSE,
    jags.seed = 0)})
}
