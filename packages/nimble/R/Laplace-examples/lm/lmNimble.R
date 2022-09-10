rm(list=ls())
library(nimble)
nimbleOptions(buildDerivs = TRUE)
## nimbleFunction for maximum likelihood
source("nimMaxLik.R") 
source("lmTMB.R")
## Model code
code <- nimbleCode({
  a ~ dnorm(0, sd = 10.0)
  b ~ dnorm(0, sd = 10.0)
  logSigma ~ dnorm(0, sd = 10.0)
  for(i in 1:n){
    y[i] ~ dnorm(a + b*x[i], sd = exp(logSigma))
  }
})
## Constants, initial values and data
constants <- list(n = length(data$Y))
inits <- list(a = 1,
              b = 1,
              logSigma = 1,
              x = data$x)
dat <- list(y = data$Y)
## Build and compile model
model <- nimbleModel(code, constants = constants, data = dat, inits = inits, buildDerivs = T)
cmodel <- compileNimble(model)
cmodel$calculate()
## MLE
paramNodes <- c("a", "b", "logSigma")
calcNodes <- model$getDependencies(paramNodes, self = FALSE)
mle <- nimMaxLik(model, paramNodes, calcNodes)
cmle <- compileNimble(mle, project = model, resetFunctions = T)
p <- values(model, paramNodes)
nimtime <- system.time(nimres <- cmle$maxLik(p))

## Runtimes and MLEs are very close
rbind(tmbtime, nimtime)
rbind(tmbres$par, nimres$par)


