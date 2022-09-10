## Poisson and negative binomial GLMMs
rm(list=ls())
library(nimble)
library(glmmTMB)
nimbleOptions(buildDerivs = TRUE)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

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

## Create the model in nimble
pmodel <- nimbleModel(code = pcode,
                      constants = list(N=N, n=n, site=site, nsites=nsites, X=X), 
                      data = list(y = Salamanders$count),
                      inits = list(beta = rep(1, n),  mu = rep(0, nsites), log_sigma = 0),
                      buildDerivs = T)
## Compile the model
Cpmodel <- compileNimble(pmodel)
Cpmodel$calculate()

## Setup for Laplace
# randomEffectsNodes <- "mu"
# calcNodes <- pmodel$getDependencies(randomEffectsNodes)
# paramNodes <- c("beta", "log_sigma")

## Laplace approximation 
plaplace <- buildLaplace(pmodel)
Cplaplace <- compileNimble(plaplace, project = pmodel, resetFunctions = TRUE)

## Initial parameter values
p <- values(pmodel, c("beta", "log_sigma"))
Cplaplace$Laplace(p)
Cplaplace$gr_Laplace(p)

## MLE: all the same
Cplaplace$set_method(1) ## Single taping for derivatives calculation with separate components
nimtimep1 <- system.time(nimresp1 <- Cplaplace$LaplaceMLE(p))
nimsummp1 <- Cplaplace$summary(nimresp1)

Cplaplace$set_method(2) ## Double taping for derivatives calculation with separate components
nimtimep2 <- system.time(nimresp2 <- Cplaplace$LaplaceMLE(p))
nimsummp2 <- Cplaplace$summary(nimresp2)

Cplaplace$set_method(3) ## Double taping with everything packed together
nimtimep3 <- system.time(nimresp3 <- Cplaplace$LaplaceMLE(p))
nimsummp3 <- Cplaplace$summary(nimresp3)

## Use nlminb
nimtimep.nlminb <- system.time(
  nimresp.nlminb <- nlminb(p,
                          function(x) -Cplaplace$Laplace(x),
                          function(x) -Cplaplace$gr_Laplace(x))
)

## Fit the model using glmmTMB
## TMB is slightly slower here, but the time includes that for compilation etc.
tmbtimep <- system.time(tmbresp <- glmmTMB(count ~ spp + mined + (1|site), Salamanders, family=poisson))
tmbsummp <- summary(tmbresp)

rbind(nimtimep1, nimtimep2, nimtimep3, nimtimep.nlminb, tmbtimep)
## MLEs are very close: 
rbind(nimsummp1$estimate, 
      nimsummp2$estimate,
      nimsummp3$estimate,
      nimresp.nlminb$par,
      tmbresp$fit$par)

## Std. errors for regression coefficients
rbind(nimsummp1$stdError[1:8], 
      nimsummp2$stdError[1:8],
      nimsummp3$stdError[1:8],
      c(tmbsummp$coefficients$cond[,"Std. Error"]))

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
                       buildDerivs = T)
## Compile the model
Cnbmodel <- compileNimble(nbmodel)
Cnbmodel$calculate()

## Setup for Laplace
# randomEffectsNodes <- "mu"
# calcNodes <- nbmodel$getDependencies(randomEffectsNodes)
# paramNodes <- c("beta", "log_phi", "log_sigma")
# nblaplace <- buildLaplace(nbmodel, paramNodes, randomEffectsNodes, calcNodes, control = list())
## Laplace approximation 
nblaplace <- buildLaplace(nbmodel)
Cnblaplace <- compileNimble(nblaplace, project = nbmodel, resetFunctions = T)

## A single Laplace calculation
paramNodes <- c("beta", "log_phi", "log_sigma")
p <- values(nbmodel, paramNodes)
Cnblaplace$Laplace(p)
Cnblaplace$gr_Laplace(p)

## MLE
Cnblaplace$set_method(1) ## Single taping for derivatives calculation with separate components
nimtimenb1 <- system.time(nimresnb1 <- Cnblaplace$LaplaceMLE(p))
nimsummnb1 <- Cnblaplace$summary(nimresnb1)

Cnblaplace$set_method(2) ## Double taping for derivatives calculation with separate components
nimtimenb2 <- system.time(nimresnb2 <- Cnblaplace$LaplaceMLE(p))
nimsummnb2 <- Cnblaplace$summary(nimresnb2)

Cnblaplace$set_method(3) ## Double taping with everything packed together
nimtimenb3 <- system.time(nimresnb3 <- Cnblaplace$LaplaceMLE(p))
nimsummnb3 <- Cnblaplace$summary(nimresnb3)

## Use nlminb
nimtimenb.nlminb <- system.time(
  nimresnb.nlminb <- nlminb(p,
                           function(x) -Cnblaplace$Laplace(x),
                           function(x) -Cnblaplace$gr_Laplace(x))
)

## Use glmmTMB
tmbtimenb <- system.time(tmbresnb <- glmmTMB(count ~ spp + mined + (1|site), Salamanders, family=nbinom2))
tmbsummnb <- summary(tmbresnb)

## Compare runtimes
rbind(nimtimenb1, nimtimenb2, nimtimenb3, nimtimenb.nlminb, tmbtimenb)
## Compare MLEs: very close
## MLEs are very close: 
rbind(nimsummnb1$estimate, 
      nimsummnb2$estimate,
      nimsummnb3$estimate,
      nimresnb.nlminb$par,
      tmbresnb$fit$par)

## Standard errors of regression coefficients: very close 
rbind(nimsummnb1$stdError[1:8], 
      nimsummnb2$stdError[1:8],
      nimsummnb3$stdError[1:8],
      tmbsummnb$coefficients$cond[,"Std. Error"])
