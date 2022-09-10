rm(list=ls())
library(nimble)
source("spatialTMB.R")
## *** Note that parametrizing the model using a / log(a) causes a (big) difference to 
## the efficiency of optimization using TMB. The latter is faster.  
tmbtime ## TMB run time: ~0.7s

nimbleOptions(buildDerivs = TRUE)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE) 

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
Cmodel$calculate() 

paramNodes <- c("b", "log_a", "log_sigma") 
model$getNodeNames(topOnly = TRUE, stochOnly = TRUE) ## Order of parameters by default

## Build Laplace: do everything using the default settings
laplace <- buildLaplace(model)
Claplace <- compileNimble(laplace, project = model, resetFunctions = TRUE)

## MLEs
p <- values(model, paramNodes)
Claplace$Laplace(p)
Claplace$gr_Laplace(p)

Claplace$set_method(1) ## Single taping for derivatives calculation with separate components
nimtime1 <- system.time(nimres1 <- Claplace$LaplaceMLE(p))
nimSumm1 <- Claplace$summary(nimres1)

Claplace$set_method(2) ## Double taping for derivatives calculation with separate components
nimtime2 <- system.time(nimres2 <- Claplace$LaplaceMLE(p))
nimSumm2 <- Claplace$summary(nimres2)

## Use nlminb
nimtime_nlmb <- system.time(
  nimres_nlmb <- nlminb(p, 
                        function(x) -Claplace$Laplace(x),
                        function(x) -Claplace$gr_Laplace(x))
)
## nimble + nlminb gives very close results to TMB for out optimization
rbind(nimres_nlmb$par, tmbres$par)
## Compare run times 
rbind(nimtime1, nimtime2, nimtime_nlmb, tmbtime)
## MLEs are very close
rbind(nimSumm1$estimate, nimSumm2$estimate, nimres_nlmb$par, tmbres$par)
## Std. errors: very close as well
rbind(nimSumm1$stdError, nimSumm2$stdError, rep[,"Std. Error"])
