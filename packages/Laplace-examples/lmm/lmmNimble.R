rm(list=ls())
library(nimble)
nimbleOptions(buildDerivs = TRUE)
## nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
source("lmmTMB.R")

## Construct a nimble model
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

## Build and compile model
model <- nimbleModel(code, constants = constants, data = data, inits = inits, buildDerivs = TRUE)
cmodel <- compileNimble(model)
cmodel$calculate()

## Build and compile Laplace 
laplace <- buildLaplace(model) 
claplace <- compileNimble(laplace, project = model)

## A single evaluation of Laplace and its gradient
claplace$calcLaplace(c(1, 1, 0, 0))
claplace$gr_Laplace(c(1, 1, 0, 0))

## Calculate MLEs and standard errors
## The warning can be safely ignored
nimtime <- system.time(opt <- claplace$findMLE())
nimtime2 <- system.time(nimres <- claplace$summary(opt, randomEffectsStdError = TRUE))

## Call nlminb for the optimisation, which is slightly faster
fn <- function(x) -claplace$calcLaplace(x)
gr <- function(x) -claplace$gr_Laplace(x)
nimtime3 <- system.time(opt2 <- nlminb(c(1, 1, 0, 0), fn, gr))

## Compare results
rbind(nimres$params$estimates, tmbres$par)
rbind(nimres$params$stdErrors, tmbsumm[1:4, "Std. Error"])
max(abs(nimres$randomEffects$estimates - tmbsumm[5:118, "Estimate"]))
max(abs(nimres$randomEffects$stdErrors - tmbsumm[5:118, "Std. Error"]))
