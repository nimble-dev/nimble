rm(list=ls())
library(nimble)
## source("../DLaplace.R")
nimbleOptions(buildDerivs = TRUE)
source("lmmTMB.R")
## Model code
code <- nimbleCode({
  logsdu  ~ dnorm(0, sd = 10.0)
  logsd0  ~ dnorm(0, sd = 10.0)
  beta[1] ~ dnorm(0, sd = 10.0)
  beta[2] ~ dnorm(0, sd = 10.0)
  for(i in 1:nre){
    u[i] ~ dnorm(0, sd = exp(logsdu))
  }
  mean[1:nobs] <- A[1:nobs, 1:2] %*% beta[1:2] + B[1:nobs, 1:nre] %*% u[1:nre]
  for(j in 1:nobs){
    x[j] ~ dnorm(mean[j], sd = exp(logsd0))
  }
})
constants <- list(nre = length(u),
                  nobs = length(x),
                  A = A_nonsparse,
                  B = B_nonsparse)
inits <- list(logsdu = 1,
              logsd0 = 1,
              u = u*0,
              beta = beta*0)
data <- list(x = x)
model <- nimbleModel(code, constants = constants, data = data, inits = inits, buildDerivs = T)
cmodel <- compileNimble(model)
cmodel$calculate()
## Laplace
paramNodes <- c("beta", "logsdu", "logsd0")
randomEffectsNodes <- "u"
calcNodes <- model$getDependencies(randomEffectsNodes)
laplace <- buildLaplace(model, paramNodes, randomEffectsNodes, calcNodes,
                        control = list(split=TRUE))
claplace <- compileNimble(laplace, project = model)
p <- values(model, paramNodes)
claplace$set_method(2)
claplace$Laplace(p)

## Calculate MLEs
## nimtime1 <- system.time(nimres1 <- claplace$LaplaceMLE1(p))
nimtime2 <- system.time(nimres2 <- claplace$LaplaceMLE2(p))
## nimtime3 <- system.time(nimres3 <- claplace$LaplaceMLE3(p)) ## This takes long...
nimtime.nlminb <- system.time(
  nimres.nlminb <- nlminb(p, 
                          function(x) {-claplace$Laplace2(x)},
                          function(x) {-claplace$gr_Laplace2(x)})
)

## Results match well, but nimble is slower than TMB
rbind(nimtime2, nimtime.nlminb, tmbtime)
rbind(nimres2$par, nimres.nlminb$par, tmbres$par)
