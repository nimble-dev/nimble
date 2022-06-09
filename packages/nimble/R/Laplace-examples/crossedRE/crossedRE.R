rm(list=ls())
library(nimble)
library(lme4)
library(glmmTMB)
nimbleOptions(buildDerivs = TRUE)
## Load data
data(Penicillin)
## Fit the model using glmmTMB
tmbres <- glmmTMB(diameter ~ 1 + (1|plate) + (1|sample), data=Penicillin, family=gaussian)

## Prepare data for nimble
N <- nrow(Penicillin) 
plate <- rep(1:24, 6)
np <- 24
sample <- rep(1:6, 24)
ns <- 6
## nimble model code
code <- nimbleCode({
  ## Intercept
  beta ~ dnorm(0, sd = 100)
  ## Standard deviations
  sigma  ~ dgamma(0.1, 1)
  sigmap ~ dgamma(0.1, 1)
  sigmas ~ dgamma(0.1, 1)
  ## Random effects for plate
  for(i in 1:np){
    mup[i] ~ dnorm(0, sd = sigmap)
  }
  ## Random effects for sample
  for(i in 1:ns){
    mus[i] ~ dnorm(0, sd = sigmas)
  }
  ## Observations
  for(i in 1:N){
    y[i] ~ dnorm(beta + mup[plate[i]] + mus[sample[i]], sd = sigma)
  }
})
## Build and compile the model
constants <- list(N=N, np=np, ns=ns, plate=plate, sample=sample)
data <- list(y=Penicillin$diameter)
inits <- list(beta = 1,  mup = rep(0, np), mus = rep(0,ns), 
              sigma = 1, sigmap = 1, sigmas = 1)
model <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
cmodel <- compileNimble(model)
cmodel$calculate()

## Setup for Laplace
randomEffectsNodes <- c("mup", "mus")
calcNodes <- model$getDependencies(randomEffectsNodes)
paramNodes <- c("beta", "sigmap", "sigmas", "sigma")

## Build Laplace
laplace <- buildLaplace(model, paramNodes, randomEffectsNodes, calcNodes, 
                        control = list(split = FALSE, warn = TRUE))
## Compile Laplace 
comptime <- system.time(
  Claplace <- compileNimble(laplace, project = model)
)
Claplace$set_method(2)
p <- values(model, paramNodes)
Claplace$Laplace(p)
nimtime <- system.time(nimres <- Claplace$LaplaceMLE(p))

## Fixed effect OK, but it seems we are not working correctly for the standard deviations...
allres <- rbind(nimres$par, c(tmbres$fit$par[1], exp(tmbres$fit$par[3:4]), sigma(tmbres)))
colnames(allres) <- c("beta", "sigmap", "sigmas", "sigma")
allres
