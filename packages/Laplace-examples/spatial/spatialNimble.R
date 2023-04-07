rm(list=ls())
library(nimble)
source("spatialTMB.R")
## *** Note that parametrizing the model using a / log(a) causes a (big) difference to 
## the efficiency of optimization using TMB. The latter is faster.  
nimbleOptions(buildDerivs = TRUE)
## nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE) 

## Model code
code <- nimbleCode({
  ## Priors
  for(i in 1:2){
    b[i] ~ dnorm(0, sd = 100)
  }
  log_a ~ dnorm(0, sd = 100)
  a <- exp(log_a)
  ## a ~ dgamma(0.1, 1.0)
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
  u[1:n] ~ dmnorm(mean[1:n], cov = cov[1:n, 1:n])
  
  ## Multivariate Poisson
  eta[1:n] <- X[,] %*% b[1:2] + exp(log_sigma) * u[1:n]
  for(i in 1:n){
    y[i] ~ dpois(exp(eta[i]))
  }
})
## Source data
source("spatial_data.R")
dd <- as.matrix(dist(Z))
model <- nimbleModel(code,
                     constants = list(n = n, dd = dd, X = X),
                     data = list(y = y),
                     inits = list(b = c(1, 1), log_a = log(2), log_sigma = -1, u = rep(1, n)),
                     buildDerivs = TRUE)
Cmodel <- compileNimble(model)
## Cmodel$calculate() 

## Build Laplace
laplace <- buildLaplace(model)
Claplace <- compileNimble(laplace, project = model)
Claplace$Laplace(rep(1, 4))
Claplace$gr_Laplace(rep(1, 4))

## Calculate MLEs and standard errors
nimtime <- system.time(opt <- Claplace$LaplaceMLE())
nimtime2 <- system.time(nimres <- Claplace$summary(opt, calcRandomEffectsStdError = TRUE))

## Compare results:
rbind(nimtime + nimtime2, tmbtime + tmbtime2)
## Parameters
rbind(nimres$params$estimate, tmbres$par)
rbind(nimres$params$stdError, tmbsumm[1:4,"Std. Error"])
## Random effects
max(abs(nimres$random$estimate - tmbsumm[5:104, "Estimate"]))
max(abs(nimres$random$stdError - tmbsumm[5:104,"Std. Error"]))
