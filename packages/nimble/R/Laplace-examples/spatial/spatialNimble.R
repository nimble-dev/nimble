## A spatial Poisson GLMM 
rm(list=ls())
library(nimble)
nimbleOptions(buildDerivs = TRUE)

## Model code
spatialCode <- nimbleCode({
  ## Priors
  for(i in 1:2){
    b[i] ~ dnorm(0, sd = 100)
  }
  a ~ dnorm(0, sd = 100)
  log_sigma ~ dnorm(0, sd = 100)
  
  ## Covariance matrix
  for(i in 1:n){
    cov[i, i] <- 1
    for(j in 1:(i-1)){
      cov[i, j] <- exp(-a * dd[i, j])
      cov[j, i] <- cov[i, j]
    }
  }
  ## Zero mean 
  for(i in 1:n){
    mean[i] <- 0
  }
  ## Multivariate normal
  ## prec[1:n, 1:n] <- inverse(cov[1:n, 1:n])
  ## u[1:n] ~ dmnorm(mean[1:n], prec = prec[1:n, 1:n])
  u[1:n] ~ dmnorm(mean[1:n], cov = cov[1:n, 1:n])
  
  ## Multivariate Poisson
  eta[1:n] <- X[,] %*% b[1:2] + exp(log_sigma) * u[1:n]
  for(i in 1:n){
    y[i] ~ dpois(exp(eta[i]))
  }
})
## Load data
source("spatial_data.R")
dd <- as.matrix(dist(Z))

## Build the model
spatialGLMM <- nimbleModel(spatialCode, 
                           name = "spatialGLMM",
                           constants = list(n = n, dd = dd),
                           data = list(y = y, X = X),
                           inits = list(b = c(1, 1), a = 2, log_sigma = -1, u = rep(1, n)),
                           buildDerivs = TRUE
                           )
## Compile the model
CspatialGLMM <- compileNimble(spatialGLMM)

## Setup
randomEffectsNodes <- "u"
calcNodes <- spatialGLMM$getDependencies(randomEffectsNodes)
paramNodes <- c("b", "a", "log_sigma")

## Values of parameters and random effects nodes
p  <- values(spatialGLMM, paramNodes)
re <- values(spatialGLMM, randomEffectsNodes)

## Laplace approximation
spLaplace  <- buildLaplace(spatialGLMM, paramNodes, randomEffectsNodes, calcNodes, control = list(split = FALSE))
CspLaplace <- compileNimble(spLaplace, project = spatialGLMM)

## Method 2 seems to be the best for this spatial model
method <- 2
CspLaplace$set_method(method)
## Correct answers are returned (first-time run is much slower)
system.time(nimres <- CspLaplace$Laplace(p))
system.time(tmbres <- obj$fn(p))
system.time(nimres2 <- CspLaplace$Laplace(p+1))
system.time(tmbres2 <- obj$fn(p+1))

## Calculate MLEs
nimtime1 <- system.time(nimres1 <- CspLaplace$LaplaceMLE1(p))
nimtime2 <- system.time(nimres2 <- CspLaplace$LaplaceMLE2(p))
nimtime21 <- system.time(nimres.nlminb <- nlminb(p+1, function(x) -CspLaplace$Laplace2(x),
                                          function(x) -CspLaplace$gr_Laplace2(x),
                                          lower=c(-100.0, -100.0, 0.01, -3.0),
                                          upper=c( 100.0,  100.0, 3.00,  3.0)))
nimtime3 <- system.time(nimres3 <- CspLaplace$LaplaceMLE3(p))
## Compare nimble and TMB
source("spatialTMB.R") 
c(nimtime1, nimtime2, tmbtime)
rbind(nimres1$par, nimres2$par, tmbres$par)


