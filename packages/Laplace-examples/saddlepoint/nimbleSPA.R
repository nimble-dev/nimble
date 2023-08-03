rm(list=ls())
## Install nimble from the saddlepoint branch / source "Laplace.R" from the branch
library(nimble)
## We first define a custom normal distribution for saddlepoint approximation.
## Note that this is not a valid PDF, but we will obtain an approximate log-likelihood 
## after integrating out "s", the independent variable of the cumulant generating 
## function of the normal distribution. We can think of the custom normal distribution
## has three variables, mu, sigma, and s and we will need to "integrate out" s through
## the saddlepoint approximation (computations are very similar to those of the 
## Laplace approximation). In the saddlepoint approximation, we need to solve s from a 
## saddlepoint equation, which involves s and the model parameters mu and sigma. 
## This means that the solution of s depends on mu and sigma. Thinking about this using
## the logic of Laplace approximation, s is a latent variable to be integrated out, and 
## mu and sigma are top-level parameters. 

dsaddlepoint_normal <- nimbleFunction(
  run = function(x = double(0),
                 mu = double(0),
                 sigma = double(0),
                 s = double(0),
                 log = logical(0, default = TRUE)){
    ## Cumulant generating function of normal with mean mu and sd sigma
    cgf <- mu*s + 0.5*sigma^2*s^2 
    res <- s*x - cgf
    if(log) return(res)
    else return(exp(res))
    returnType(double())
  },
  buildDerivs = TRUE
)
registerDistributions('dsaddlepoint_normal')

## Model code using the custom normal distribution
code <- nimbleCode({
  mu ~ dnorm(0, sd = 10)
  sigma ~ dhalfflat()
  for (i in 1:10){
    ## s[i] essentially does not follow a distribution
    ## To use the Laplace functionality in nimble, we need to specify them as latent states
    ## The distribution used here does not matter
    s[i] ~ dnorm(mu, sd = sigma)
    y[i] ~ dsaddlepoint_normal(mu = mu, sigma = sigma, s = s[i])
  }
})
## Generate data
set.seed(1)
y <- rnorm(10, 0, 1)
model <- nimbleModel(code, 
                     data = list(y = y), 
                     inits = list(mu = 0, sigma = 1, s = rep(0, 10)),
                     buildDerivs = TRUE)
cmodel <- compileNimble(model)
## Note that the calcNodes needed are the data nodes "y" only; do NOT include "s".
## Because the full likelihood has all been given in the custom normal distribution
normalSpa <- buildAGHQuad(model,
                          paramNodes = c("mu", "sigma"),
                          randomEffectsNodes = "s",
                          calcNodes = "y",
                          calcNodesOther = NULL,
                          control = list(saddlepoint = TRUE, split = TRUE, check = FALSE))
cnormalSpa <- compileNimble(normalSpa, project = model)
mle <- cnormalSpa$findMLE()

## Use the true normal distribution for comparison
code2 <- nimbleCode({
  mu ~ dnorm(0, sd = 10)
  sigma ~ dhalfflat()
  for (i in 1:10){
    y[i] ~ dnorm(mu, sd = sigma)
  }
})
model2 <- nimbleModel(code2, 
                      data = list(y = y), 
                      inits = list(mu = 0, sigma = 1),
                      buildDerivs = TRUE)
                     
cmodel2 <- compileNimble(model2)
## No latent states here, so Laplace is just maximum likelihood
laplace <- buildLaplace(model2)
claplace <- compileNimble(laplace, project = model2)
mle2 <- claplace$findMLE()

## Very close results
rbind(mle$par, mle2$par)
rbind(mle$value, mle2$value)
