## Poisson and negative binomial GLMMs
rm(list=ls())
library(nimble)
library(glmmTMB)
nimbleOptions(buildDerivs = TRUE)
## nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

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
  log_sigma ~ dnorm(0, sd = 100)
  sigma <- exp(log_sigma)
  ## sigma ~ dunif(0, 10)
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

## Create model
pmodel <- nimbleModel(code = pcode,
                      constants = list(N=N, n=n, site=site, nsites=nsites, X=X), 
                      data = list(y = Salamanders$count),
                      inits = list(beta = rep(1, n),  mu = rep(0, nsites), log_sigma = 0),
                      buildDerivs = TRUE)
## Compile the model
Cpmodel <- compileNimble(pmodel)
Cpmodel$calculate()

## Setup for Laplace
# randomEffectsNodes <- "mu"
# calcNodes <- pmodel$getDependencies(randomEffectsNodes)
# paramNodes <- c("beta", "log_sigma")

## Laplace approximation 
plaplace <- buildLaplace(pmodel)
Cplaplace <- compileNimble(plaplace, project = pmodel)
## A single evaluation of Laplace and its gradient
Cplaplace$calcLaplace(rep(1, 9))
Cplaplace$gr_Laplace(rep(1, 9))

## Calculate MLEs and standard errors
## Warnings here can be safely ignored
popt <- Cplaplace$findMLE()
nimresp <- Cplaplace$summary(popt, randomEffectsStdError = TRUE)

## Fit the model using glmmTMB
tmbresp <- glmmTMB(count ~ spp + mined + (1|site), Salamanders, family=poisson)
tmbsummp <- summary(tmbresp)

## MLEs are very close: 
rbind(nimresp$params$estimates, tmbresp$fit$par)

## Std. errors for regression coefficients
rbind(nimresp$params$stdErrors[1:8], tmbsummp$coefficients$cond[,"Std. Error"])

##------------------------------------------------------------------------------
## Negative binomial model
## Model code
nbcode <- nimbleCode({
  ## Fixed effects
  for(i in 1:n){
    beta[i] ~ dnorm(0, sd = 100)
  }
  ## Size param for dnbinom
  ## phi ~ dgamma(0.1, 1) 
  log_phi ~ dnorm(0, sd = 100)
  phi <- exp(log_phi)
  ## Random effects for sites
  ## sigma ~ dunif(0, 10)
  log_sigma ~ dnorm(0, sd = 100)
  sigma <- exp(log_sigma)
  for(i in 1:nsites){
    mu[i] ~ dnorm(0, sd = sigma)
  }
  ## Linear predictor
  lp[1:N] <- X[1:N,1:n] %*% beta[1:n] 
  ## Observations
  for(i in 1:N){
    m[i] <- exp(lp[i] + mu[site[i]]) ## Mean for negative binomial
    y[i] ~ dnbinom(size = phi, prob = phi/(m[i]+phi))
  }
})

## Create the model in nimble
nbmodel <- nimbleModel(code = nbcode,
                       constants = list(N=N, n=n, site=site, nsites=nsites, X=X), 
                       data = list(y = Salamanders$count),
                       inits = list(beta = rep(0, n),  mu = rep(0, nsites), log_sigma = 0, log_phi = 0),
                       buildDerivs = TRUE)
## Compile the model
Cnbmodel <- compileNimble(nbmodel)
Cnbmodel$calculate()

## Setup for Laplace
# randomEffectsNodes <- "mu"
# calcNodes <- nbmodel$getDependencies(randomEffectsNodes)
# paramNodes <- c("beta", "log_phi", "log_sigma")

## Laplace approximation 
nblaplace <- buildLaplace(nbmodel)
Cnblaplace <- compileNimble(nblaplace, project = nbmodel)
## A single evaluation of Laplace and its gradient
Cnblaplace$calcLaplace(rep(0, 10))
Cnblaplace$gr_Laplace(rep(0, 10))

## Calculate MLEs and standard errors
## Warnings here can be safely ignored
nbopt <- Cnblaplace$findMLE()
nimres <- Cnblaplace$summary(nbopt, randomEffectsStdError = TRUE)

## Use glmmTMB
tmbresnb <- glmmTMB(count ~ spp + mined + (1|site), Salamanders, family=nbinom2)
tmbsummnb <- summary(tmbresnb)

## MLEs are very close: 
rbind(nimres$params$estimates, tmbresnb$fit$par)

## Standard errors of regression coefficients: very close 
rbind(nimres$params$stdErrors[1:8], tmbsummnb$coefficients$cond[,"Std. Error"])
