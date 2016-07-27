################################################################################################################
## A simple test of convergence from very poor starting values to a MVGaussian target with strong correlated  ##
################################################################################################################
library(nimble)
library(coda)

## Just a toy example for generating a posterior with very strong correlation
## mvg = Multivariate Gaussian
mvgCode <- nimbleCode({ 
    y[1:d] ~ dmnorm(mean=mu[1:d], cov=sigma[1:d,1:d]) 
})

d         <- 2
rho       <- 0.99
mu        <- rep(0, d)
sigma     <- matrix(c(1,rho,rho,1),d,d)
yInitial  <- c(0, 0)             ## These very good starting values are replaced below
Constants <- list(d=d)
Inits     <- list(mu=mu, sigma=sigma, y=yInitial)

## Model 1, for testing Nimble's default sampler
mvg      <- nimbleModel(mvgCode, constants=Constants, inits=Inits)
mvg$calculate('y')
cmvg     <- compileNimble(mvg)
mcmcConf <- configureMCMC(mvg, print=TRUE)
mcmcConf$removeSamplers()
mcmcConf$printSamplers()
mcmcConf$addSampler('y', type="RW_block")
mcmcConf$printSamplers()
mcmcConf$resetMonitors()
mcmcConf$addMonitors(c('y', 'logProb_y'))
mcmcConf$getMonitors()
mcmc     <- buildMCMC(mcmcConf)
cmcmc    <- compileNimble(mcmc)

## Model 2, for testing a modified version
mvg2      <- nimbleModel(mvgCode, constants=Constants, inits=Inits)
mvg2$calculate('y')
cmvg2     <- compileNimble(mvg2)
mcmcConf2 <- configureMCMC(mvg2, print=TRUE)
mcmcConf2$removeSamplers()
mcmcConf2$addSampler('y', type="RW_block_ETA")
mcmcConf2$printSamplers() 
mcmcConf2$resetMonitors()
mcmcConf2$addMonitors(c('y', 'logProb_y'))
mcmcConf2$getMonitors()
mcmc2     <- buildMCMC(mcmcConf2)
cmcmc2    <- compileNimble(mcmc2)


#################################
## Set up two graphics devices ##
#################################
while(length(dev.list())>0) { dev.off() } ## Turn off all graphics devices
X11()
X11()

####################################################################################
## Scenario 1: High correlation, terrible starting values & default block sampler ##
####################################################################################
## Start adaptive MCMC far from target
set.seed(1)
yStart  <- c(-1, 1) * 1E10 ## Some really terrible starting values
cmvg$y  <- yStart 
nIter   <- 1E5
cmcmc$run(nIter, reset=TRUE) 
samples <- as.matrix(cmcmc$mvSamples); dim(samples)
samples <- tail(samples, nIter)
mc      <- as.mcmc(samples)
dev.set(dev.list()[1])
plot(mc)
dev.set(dev.list()[2])
plot(samples[,2],samples[,3], typ="l")
cmvg$y           ## -25061065 -25060434  ## These are very far from the target
cmvg$calculate() ## -3.155986e+14        ## And this is still quite terrible

## Continue without a resetting the adaptive algorithm
cmcmc$run(nIter, reset=FALSE) 
samples <- as.matrix(cmcmc$mvSamples); dim(samples)
samples <- tail(samples, nIter)
mc      <- as.mcmc(samples)
dev.set(dev.list()[1])
plot(mc) ## Don't be decived, these are terrible. The variance is microscopic, the sampler is stuck on a ridge & is trying to sample in the wrong direction
dev.set(dev.list()[2])
plot(samples[,2],samples[,3], typ="l")
as.vector(tail(samples, 1))
cmvg$y           ## -25061065 -25060434  ## Almost identical to before. The MCMC is stuck, the proposal distribution is stretched in the wrong direction & the has been down scaled to have a miniscule variance
cmvg$calculate() ## -3.155986e+14

## Continue with a reset
cmcmc$run(nIter, reset=TRUE) 
samples <- as.matrix(cmcmc$mvSamples); dim(samples)
samples <- tail(samples, nIter)
mc      <- as.mcmc(samples)
dev.set(dev.list()[1])
plot(mc)
dev.set(dev.list()[2])
plot(samples[,2],samples[,3], typ="l") ## The sampler has changed direction now
as.vector(tail(samples, 1))
cmvg$y           ## 1.116239 1.119036 ## We're close now
cmvg$calculate() ## -0.5072501
## Let's examine the hill climbing a bit closer
plot(log(samples[,1]-samples[1,1]+1), xlab="iteration (i)", ylab="log (logProb[i] - logProb[1] + 1)", typ="l") ## 
## Remove a burn-in period
dev.set(dev.list()[1])
burn  <- 1E4
mc <- as.mcmc(tail(samples, nIter-burn))
plot(mc)
dev.set(dev.list()[2])
burn  <- 2E4
mc <- as.mcmc(tail(samples, nIter-burn))
plot(mc)
## Now the sampling looks good.
## The main point here is RW_block required 1 reset before it could hit the target.
## For MCMC on a more complex model with few clues about starting values this is an undesirable characteristic.
## A method that auto-resets should be prefered.


##########################################################################################################
## Scenario 2: High correlation, terrible starting values & default block sampler + univariate samplers ##
##########################################################################################################
## Start adaptive MCMC far from target
set.seed(1)
mvg2$y  <- yStart           ## The same really terrible starting values
cmvg2$y <- yStart           ## The same really terrible starting values
nIter   <- 1E5
cmcmc2$run(nIter, reset=TRUE)  ## mcmc2$run(nIter, reset=TRUE) 
samples <- as.matrix(cmcmc2$mvSamples); dim(samples)
samples <- tail(samples, nIter)
mc2     <- as.mcmc(samples)
dev.set(dev.list()[1])
plot(mc2)
dev.set(dev.list()[2])
plot(samples[,2],samples[,3], typ="l")
cmvg2$y           ## -0.5164177 -0.6477722
cmvg2$calculate() ## -0.480978

## Let's examine the hill climbing a bit closer
plot(log(samples[,1]-samples[1,1]+1), xlab="iteration (i)", ylab="log (logProb[i] - logProb[1] + 1)", typ="l") ## 

## Remove a burn-in period
dev.set(dev.list()[1])
burn  <- 1E4
mc2 <- as.mcmc(tail(samples, nIter-burn))
plot(mc2)
dev.set(dev.list()[2])
burn  <- 2E4
mc2 <- as.mcmc(tail(samples, nIter-burn))
plot(mc2) ## Convergence already looks good. With the more flexible adaptation convergence is faster and does not require a user to reset.

par(mfrow=n2mfrow(2))
dim(as.matrix(mc[,2:3])) 
dim(as.matrix(mc2[,2:3]))
plot(as.matrix(mc[,2:3]), pch=19, cex=0.2, main="RW_block")
plot(as.matrix(mc2[,2:3]), pch=19, cex=0.2, main="RW_block_ETA") 
