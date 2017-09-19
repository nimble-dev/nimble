rm(list=ls())
library(nimble)
library(parallel)
setwd('occupancy/analysis/crossValidation')
source('crossValidationFunction.R')

dyesCode <- nimbleCode({
  for (i in 1:BATCHES) {
    for (j in 1:SAMPLES) {
      y[i,j] ~ dnorm(mu[i], sd = sigma.within);
    }
    mu[i] ~ dnorm(theta, sd = sigma.between);
  }

  theta ~ dnorm(0.0, 1.0E-10);
  sigma.within ~ dunif(0, 100)
  sigma.between ~ dunif(0, 100)
})

dyesModel <- nimbleModel(dyesCode,
                         constants = list(BATCHES = 6, SAMPLES = 5))

## data <- matrix(c(1545, 1540, 1595, 1445, 1595, 1520, 1440, 1555, 1550,
##                  1440, 1630, 1455, 1440, 1490, 1605, 1595, 1515, 1450,
##                  1520, 1560, 1510, 1465, 1635, 1480, 1580, 1495, 1560,
##                  1545, 1625, 1445), nrow = 6)

data <- cbind(rnorm(6, 0, 1), rnorm(6, 6, 1), rnorm(6, 4, 1),
              rnorm(6, 7, 1),
              rnorm(6, 5, 1))

dyesModel$setData(list(y = data))

options(mc.cores=1)
niter <- 10000
output <- runCrossValidate(model=dyesModel,
                           dataNames= "y",
                           MCMCIter= niter,
                           burnInProp=0.1,
                           thin=1,
                           leaveOutIndex=2,
                           MCMCdefs=NULL)


## ****************************************************
dyesCodeSimp <- nimbleCode({
  for (i in 1:BATCHES) {
    for (j in 1:SAMPLES) {
          y[i,j] ~ dnorm(theta, sd = sigma)
    }
  }
  theta ~ dnorm(0.0, 1.0E-10);
  sigma ~ dunif(0, 100)
})

dyesModelSimp <- nimbleModel(dyesCodeSimp,
                             constants =list(BATCHES = 6, SAMPLES = 5))
dyesModelSimp$setData(list(y = data))

output.simp <- runCrossValidate(model=dyesModelSimp,
                           dataNames= "y",
                           MCMCIter= niter,
                           burnInProp=0.1,
                           thin=1,
                           leaveOutIndex=2,
                           MCMCdefs=NULL)
