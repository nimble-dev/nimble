rm(list=ls())
library(nimble)
library(glmmTMB)
nimbleOptions(buildDerivs = TRUE)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

## Load data
data(cbpp, package="lme4")

## Prepare data for nimble
N <- nrow(cbpp) 
## Construct the design matrix
X <- rep(1, N)
periods <- levels(cbpp$period)
nperiod <- length(periods)
for(i in 2:nperiod) X <- cbind(X, ifelse(cbpp$period==periods[i], 1, 0))
n <- ncol(X)

## Make variable herd numeric
herd <- as.numeric(cbpp$herd)
nherd <- length(unique(herd))

## nimble model code
code <- nimbleCode({
  ## Fixed effects
  for(i in 1:n){
    beta[i] ~ dnorm(0, sd = 100)
  }
  ## sigma ~ dunif(0, 10)
  log_sigma ~ dnorm(0, sd = 100)
  sigma <- exp(log_sigma)
  ## Random effects
  for(i in 1:nherd){
    mu[i] ~ dnorm(0, sd = sigma)
  }
  ## Linear predictor
  lp[1:N] <- X[1:N,1:n] %*% beta[1:n] 
  ## Observations
  for(i in 1:N){
    logit(p[i]) <- lp[i] + mu[herd[i]]
    y[i] ~ dbinom(size = size[i], prob = p[i])
  }
})

## Build and compile the model
constants <- list(n = n, nherd = nherd, N = N, X = X, herd = herd, size = cbpp$size)
data <- list(y = cbpp$incidence)
inits <- list(beta = rep(0, n),  mu = rep(0, nherd), log_sigma = 0)
model <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
cmodel <- compileNimble(model)
cmodel$calculate()

## Setup for Laplace
# randomEffectsNodes <- "mu"
# calcNodes <- model$getDependencies(randomEffectsNodes)
# paramNodes <- c("beta", "sigma")

## Build Laplace
laplace <- buildLaplace(model, control = list())
## Compile Laplace 
comptime <- system.time(Claplace <- compileNimble(laplace, project = model, resetFunctions = TRUE))
## Check it works
p <- values(model, c("beta", "log_sigma"))
Claplace$Laplace(p)
Claplace$gr_Laplace(p)

## Calculate MLEs
claplace$set_method(1) ## Single taping for derivatives calculation with separate components
nimtime1 <- system.time(nimres1 <- Claplace$LaplaceMLE(p))
nimsumm1 <- Claplace$summary(nimres1)

claplace$set_method(2) ## Double taping for derivatives calculation with separate components
nimtime2 <- system.time(nimres2 <- Claplace$LaplaceMLE(p))
nimsumm2 <- Claplace$summary(nimres2)

claplace$set_method(3) ## Double taping with everything packed together
nimtime3 <- system.time(nimres3 <- Claplace$LaplaceMLE(p)) 
nimsumm3 <- Claplace$summary(nimres3)

## Use nlminb
nimtime.nlminb <- system.time(
  nimres.nlminb <- nlminb(p,
                          function(x) -Claplace$Laplace(x),
                          function(x) -Claplace$gr_Laplace(x))
)

## Fit the model using glmmTMB
tmbtime <- system.time(tmbres <- glmmTMB(cbind(incidence, size-incidence) ~ period + (1|herd), 
                                         family=binomial, data=cbpp))
summtmbres <- summary(tmbres)

## MLEs match well
rbind(nimtime1, nimtime2, nimtime3, nimtime.nlminb, tmbtime)
rbind(nimsumm1$estimate, 
      nimsumm2$estimate,
      nimsumm3$estimate,
      nimres.nlminb$par,
      tmbres$fit$par)

## Std. errors: very close
rbind(nimsumm1$stdError[1:4], 
      nimsumm2$stdError[1:4],
      nimsumm3$stdError[1:4],
      summtmbres$coefficients$cond[,"Std. Error"])
