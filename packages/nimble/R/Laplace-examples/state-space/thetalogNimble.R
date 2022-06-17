## A state-space model for testing nimble Laplace 
library(nimble)
nimbleOptions(buildDerivs = TRUE)

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
    ## When pow or '^' is used here, fatal errors occur when calling Laplace.
    ## mean.u[t-1] <- u[t - 1] + r0 * (1 - pow(exp(u[t - 1]) / K, psi))
    ## mean.u[t-1] <- u[t - 1] + r0 * (1 - psi * exp(u[t - 1]) / K)
    mean.u[t-1] <- u[t - 1] + r0 * psi * (1 - exp(u[t - 1]) / K)
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
thetalogLaplace <- buildLaplace(thetalog, paramNodes, randomEffectsNodes, calcNodes, control= list(split = FALSE))
CthetalogLaplace <- compileNimble(thetalogLaplace, project = thetalog)

## Compute one Laplace approximation
CthetalogLaplace$set_method(2)
CthetalogLaplace$Laplace(p)

## Calculate MLE using different methods 
nimtime1 <- system.time(nimres1 <- CthetalogLaplace$LaplaceMLE1(p-1))
nimtime2 <- system.time(nimres2 <- CthetalogLaplace$LaplaceMLE2(p))
## nimtime3 <- system.time(nimres3 <- CthetalogLaplace$LaplaceMLE3(p))

## Use nlminb for outer optimization, which is faster, but still slower than TMB
## Convergence code 8 is annoying but results are very close to those of TMB
nimtime.nlminb <- system.time(
  nimres.nlminb <- nlminb(p, 
                          function(x) -CthetalogLaplace$Laplace2(x),
                          function(x) -CthetalogLaplace$gr_Laplace2(x))
)

source("thetalogTMB.R") ## TMB code
rbind(nimres1$par, nimres2$par, nimres.nlminb$par, tmbres$par)
c(nimres1$value, nimres2$value, tmbres$objective)