## A state-space model for testing nimble Laplace 
rm(list=ls())
library(nimble)
nimbleOptions(buildDerivs = TRUE)
## nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE) 

## Read in data
dat<- scan("thetalog.dat", skip=3, quiet=TRUE)
N <- length(dat) # Number of random effects (observations)

## Model code
thetalogCode <- nimbleCode({
  ## Parameters on log scale
  logr0 ~  dnorm(0, sd = 100)
  logpsi ~ dnorm(0, sd = 100)
  logK ~   dnorm(0, sd = 100)
  logQ ~   dnorm(0, sd = 100)
  logR ~   dnorm(0, sd = 100)
  
  ## Back-tranformed parameters
  r0 <-  exp(logr0)
  psi <- exp(logpsi)
  K <-   exp(logK)
  Q <-   exp(logQ)
  R <-   exp(logR)
  
  ## Initial state
  u[1] ~ dunif(-Inf, Inf)
  ## Initial observation
  y[1] ~ dnorm(u[1], sd = sqrt(R))
  
  ## States and observations for t >= 2
  for(t in 2:N){
    mean.u[t-1] <- u[t - 1] + r0 * (1 - pow(exp(u[t - 1]) / K, psi))
    u[t] ~ dnorm(mean.u[t-1], sd = sqrt(Q))
    y[t] ~ dnorm(u[t], sd = sqrt(R))
  }
})

## Code the model in nimble
thetalog <- nimbleModel(code = thetalogCode,
                        name = "thetalog",
                        constants = list(N = N), 
                        data = list(y = dat),
                        inits = list(logr0 = 0, logpsi = 0, logK = 6, logQ = 0, logR = 0,  u = rep(0, N)),
                        buildDerivs = T
                        )
## Compile the model
Cthetalog <- compileNimble(thetalog)

## Setup for Laplace
randomEffectsNodes <- "u[1:200]"
calcNodes <- thetalog$getDependencies(randomEffectsNodes)
## We exclude u[1] here
calcNodes <- setdiff(calcNodes, "u[1]") 
paramNodes <- c("logr0", "logpsi", "logK", "logQ", "logR")
## Parameters values
p <- values(thetalog, paramNodes)

## Build Laplace
## Note that currently the default setting of buildLaplace does not work properly for this model (the initial state is a headache)
thetalogLaplace <- buildLaplace(thetalog, paramNodes, randomEffectsNodes, calcNodes, 
                                control= list(split = FALSE, warn = FALSE))
CthetalogLaplace <- compileNimble(thetalogLaplace, project = thetalog)
CthetalogLaplace$Laplace(p)
CthetalogLaplace$gr_Laplace(p)

## Calculate MLEs and standard errors
nimtime <- system.time(opt <- CthetalogLaplace$LaplaceMLE(p-1))
nimtime2 <- system.time(nimres <- CthetalogLaplace$summary(opt, calcRandomEffectsStdError = TRUE))

## Run TMB code
source("thetalogTMB.R") 

## Compare runtimes, MLEs, and standard errors
rbind(nimtime + nimtime2, tmbtime + tmbtime2)
rbind(nimres$params$estimate, tmbres$par)
rbind(nimres$params$stdError, tmbsumm[1:5, "Std. Error"])
max(abs(nimres$random$estimate - tmbsumm[6:205, "Estimate"]))
max(abs(nimres$random$stdError - tmbsumm[6:205, "Std. Error"]))
