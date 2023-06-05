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
                        buildDerivs = TRUE
                        )
## Compile the model
Cthetalog <- compileNimble(thetalog)

## Setup for Laplace
randomEffectsNodes <- "u[1:200]"
calcNodes <- thetalog$getDependencies(randomEffectsNodes)
## We exclude u[1] here
calcNodes <- setdiff(calcNodes, "u[1]") 
paramNodes <- c("logr0", "logpsi", "logK", "logQ", "logR")

## Build Laplace
## Note that currently the default setting of buildLaplace does not work properly 
## for this model (the initial state is a headache)
laplace <- buildLaplace(thetalog, paramNodes, randomEffectsNodes, calcNodes, 
                        control= list(split = FALSE, check = FALSE))
cLaplace <- compileNimble(laplace, project = thetalog)

## A single evaluation of Laplace and its gradient
cLaplace$calcLaplace(c(0,0,6,0,0))
cLaplace$gr_Laplace(c(0,0,6,0,0))

## Calculate MLEs and standard errors
## Warnings here can be safely ignored
nimtime <- system.time(opt <- cLaplace$findMLE(c(0,0,6,0,0)-1))
nimtime2 <- system.time(nimres <- cLaplace$summary(opt, randomEffectsStdError = TRUE))

## Call nlminb for the optimization, which is faster
fn <- function(x) -cLaplace$calcLaplace(x)
gr <- function(x) -cLaplace$gr_Laplace(x)
nimtime3 <- system.time(opt2 <- nlminb(c(0,0,6,0,0)-1, fn, gr))

## Run TMB code
source("thetalogTMB.R") 

## Compare MLEs and standard errors
rbind(nimres$params$estimates, tmbres$par)
rbind(nimres$params$stdErrors, tmbsumm[1:5, "Std. Error"])
max(abs(nimres$randomEffects$estimates - tmbsumm[6:205, "Estimate"]))
max(abs(nimres$randomEffects$stdErrors - tmbsumm[6:205, "Std. Error"]))
