rm(list=ls())
library(nimble)
nimbleOptions(buildDerivs = TRUE)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
source("lmmTMB.R")
## Model code
code <- nimbleCode({
  logsdu  ~ dnorm(0, sd = 10.0)
  logsd0  ~ dnorm(0, sd = 10.0)
  beta[1] ~ dnorm(0, sd = 10.0)
  beta[2] ~ dnorm(0, sd = 10.0)
  for(i in 1:nre){
    u[i] ~ dnorm(0, sd = exp(logsdu))
  }
  mean[1:nobs] <- A[1:nobs, 1:2] %*% beta[1:2] + B[1:nobs, 1:nre] %*% u[1:nre]
  for(j in 1:nobs){
    x[j] ~ dnorm(mean[j], sd = exp(logsd0))
  }
})
constants <- list(nre = length(u),
                  nobs = length(x),
                  A = A_nonsparse,
                  B = B_nonsparse)
inits <- list(logsdu = 1,
              logsd0 = 1,
              u = u*0,
              beta = beta*0)
data <- list(x = x)
model <- nimbleModel(code, constants = constants, data = data, inits = inits, buildDerivs = T)
cmodel <- compileNimble(model)
cmodel$calculate()
## Laplace
laplace <- buildLaplace(model) 
claplace <- compileNimble(laplace, project = model)
## Test it works
p <- values(model, c("logsdu", "logsd0", "beta"))
claplace$Laplace(p)
claplace$gr_Laplace(p)

## Calculate MLE
claplace$set_method(1) ## Single taping for derivatives calculation with separate components
nimtime1 <- system.time(nimres1 <- claplace$LaplaceMLE(p))
nimsumm1 <- claplace$summary(nimres1)

claplace$set_method(2) ## Double taping for derivatives calculation with separate components
nimtime2 <- system.time(nimres2 <- claplace$LaplaceMLE(p))
nimsumm2 <- claplace$summary(nimres2)

# claplace$set_method(3) ## Double taping with everything packed together
# nimtime3 <- system.time(nimres3 <- claplace$LaplaceMLE(p)) ## This takes longer: 4mins 
nimtime.nlminb <- system.time(
  nimres.nlminb <- nlminb(p,
                          function(x) -claplace$Laplace(x),
                          function(x) -claplace$gr_Laplace(x))
)

## Compare run times and MLEs
rbind(nimtime1, nimtime2, nimtime.nlminb, tmbtime)
rbind(nimsumm1$estimate, nimsumm2$estimate, nimres.nlminb$par, tmbres$par)

## Std errors
rbind(nimsumm1$stdError[1:4],
      nimsumm2$stdError[1:4],
      tmbsumm[,"Std. Error"])
