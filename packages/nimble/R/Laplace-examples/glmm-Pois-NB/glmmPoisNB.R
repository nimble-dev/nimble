## Poisson and negative binomial GLMMs
rm(list=ls())
library(nimble)
library(glmmTMB)
nimbleOptions(buildDerivs = T)

## Load data from glmmTMB
data(Salamanders)
N <- nrow(Salamanders) ## Number of obs

## Construct the design matrix (spp + mined)
X <- rep(1, N)
sppnames <- levels(Salamanders$spp)
nspp <- length(sppnames)
for(i in 2:nspp) X <- cbind(X, ifelse(Salamanders$spp==sppnames[i], 1, 0))
X <- cbind(X, ifelse(Salamanders$mined=="no", 1, 0))
n <- ncol(X)

## Use numbers instead of factors for variable site
site <- as.vector(Salamanders$site)
sitenames <- levels(Salamanders$site)
nsites <- length(sitenames)
for(i in 1:nsites) site[which(Salamanders$site == sitenames[i])] <- i
site <- as.numeric(site)

##------------------------------------------------------------------------------
## Poisson model
## Model code
pcode <- nimbleCode({
  ## Fixed effects
  for(i in 1:n){
    beta[i] ~ dnorm(0, sd = 100)
  }
  ## Random effects for sites
  sigma ~ dunif(0, 10)
  for(i in 1:nsites){
    mu[i] ~ dnorm(0, sd = sigma)
  }
  ## Linear predictor
  lp[1:N] <- X[1:N,1:n] %*% beta[1:n] 
  ## Observations
  for(i in 1:N){
    lambda[i] <- exp(lp[i] + mu[site[i]])
    y[i] ~ dpois(lambda[i])
  }
})

## Create the model in nimble
pmodel <- nimbleModel(code = pcode,
                      constants = list(N=N, n=n, site=site, nsites=nsites, X=X), 
                      data = list(y = Salamanders$count),
                      inits = list(beta = rep(0, n),  mu = rep(0, nsites), sigma = 1),
                      buildDerivs = T)
## Compile the model
Cpmodel <- compileNimble(pmodel)

## Setup for Laplace
randomEffectsNodes <- "mu"
calcNodes <- pmodel$getDependencies(randomEffectsNodes)
paramNodes <- c("beta", "sigma")

## Laplace approximation 
plaplace <- buildLaplace(pmodel, paramNodes, randomEffectsNodes, calcNodes, 
                        control = list(split = FALSE, warn = TRUE))
Cplaplace <- compileNimble(plaplace, project = pmodel)
Cplaplace$set_method(2)

## A single Laplace calculation
p <- values(pmodel, paramNodes)
Cplaplace$Laplace(p)

## MLE
nimtimep <- system.time(nimresp <- Cplaplace$LaplaceMLE(p))

## Fit the model using glmmTMB
## TMB is slightly slower here, but the time includes that for compilation etc.
tmbtimep <- system.time(tmbresp <- glmmTMB(count ~ spp + mined + (1|site), Salamanders, family=poisson))
tmbsummp <- summary(tmbresp)

rbind(nimtimep, tmbtimep)
## Very close 
rbind(nimresp$par, c(tmbsummp$coefficients$cond[,"Estimate"],
                     attributes(tmbsummp$varcor$cond$site)$stddev))

##------------------------------------------------------------------------------
## Negative binomial model
## Model code
nbcode <- nimbleCode({
  ## Fixed effects
  for(i in 1:n){
    beta[i] ~ dnorm(0, sd = 100)
  }
  ## Size param for dnbinom
  phi ~ dgamma(0.1, 1) 
  ## Random effects for sites
  sigma ~ dunif(0, 10)
  for(i in 1:nsites){
    mu[i] ~ dnorm(0, sd = sigma)
  }
  ## Linear predictor
  lp[1:N] <- X[1:N,1:n] %*% beta[1:n] 
  ## Observations
  for(i in 1:N){
    m[i] <- exp(lp[i] + mu[site[i]]) ## Mean for negative binomial
    y[i] ~ dnbinom(size = phi, prob = m[i]/(m[i]+phi))
  }
})

## Create the model in nimble
nbmodel <- nimbleModel(code = nbcode,
                       constants = list(N=N, n=n, site=site, nsites=nsites, X=X), 
                       data = list(y = Salamanders$count),
                       inits = list(beta = rep(0, n),  mu = rep(0, nsites), sigma = 1, phi = 1),
                       buildDerivs = T)
## Compile the model
Cnbmodel <- compileNimble(nbmodel)

## Setup for Laplace
randomEffectsNodes <- "mu"
calcNodes <- nbmodel$getDependencies(randomEffectsNodes)
paramNodes <- c("beta", "phi", "sigma")

## Laplace approximation 
nblaplace <- buildLaplace(nbmodel, paramNodes, randomEffectsNodes, calcNodes, 
                          control = list(split = FALSE, warn = TRUE))
Cnblaplace <- compileNimble(nblaplace, project = nbmodel, resetFunctions = T)
Cnblaplace$set_method(2)

## A single Laplace calculation
p <- values(nbmodel, paramNodes)
Cnblaplace$Laplace(p)

## MLE
nimtimenb <- system.time(nimresnb <- Cnblaplace$LaplaceMLE(p))
## TMB
tmbtimenb <- system.time(tmbresnb <- glmmTMB(count ~ spp + mined + (1|site), Salamanders, family=nbinom2))

## Results do not match... Need to investigate.
