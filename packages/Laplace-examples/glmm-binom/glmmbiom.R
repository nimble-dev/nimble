rm(list=ls())
library(nimble)
library(glmmTMB)
nimbleOptions(buildDerivs = TRUE)
## nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

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

## Build and compile Laplace
laplace <- buildLaplace(model)
Claplace <- compileNimble(laplace, project = model)

## Calculate MLEs
opt <- Claplace$findMLE()
## Calculate random effects estimates and standard errors
nimres <- Claplace$summary(opt, randomEffectsStdError = TRUE)

## Fit the model using glmmTMB
tmbres <- glmmTMB(cbind(incidence, size-incidence) ~ period + (1|herd), 
                  family=binomial, data=cbpp)
tmbsumm <- summary(tmbres)

## Compare results
rbind(nimres$params$estimates, tmbres$fit$par)
rbind(nimres$params$stdErrors[1:4], tmbsumm$coefficients$cond[,"Std. Error"])

