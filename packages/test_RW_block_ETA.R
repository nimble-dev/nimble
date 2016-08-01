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
##
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
mcmcConf2$addSampler('y', type="RW_block_effTA", control = list(readaptability=0.1)) ## 0 gives original sampler, 1 gives greatest readaptibility but at the expense of rate, even values pretty close to zero (e.g. 0.01) work well on this example.
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
## Remove burn-in
lpy   <- samples[,"logProb_y[1]"]
ppy   <- exp(lpy - max(lpy))
ppy   <- ppy/sum(ppy)
(burn <- sum(ppy < 1/nrow(samples)))
mc    <- as.mcmc(tail(samples[,c("y[1]","y[2]")], nrow(samples)-burn))
plot(mc)
nrow(mc) ## 49819    
## 
## Now the sampling looks good.
## The main point here is RW_block required 1 reset before it could hit the target.
## For MCMC on a more complex model with few clues about starting values this is an undesirable characteristic.
## A method that auto-resets should be prefered.
##
effectiveSize(mc)
##     y[1]     y[2] 
## 6247.290 6574.205 
## 
summary(mc)
## Iterations = 1:49819
## Thinning interval = 1 
## Number of chains = 1 
## Sample size per chain = 49819 
##
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
##          Mean     SD Naive SE Time-series SE
## y[1] -0.03099 0.9980 0.004471        0.01263
## y[2] -0.02824 0.9975 0.004469        0.01230
## 
## 2. Quantiles for each variable:
##        2.5%     25%      50%    75% 97.5%
## y[1] -1.963 -0.7136 -0.02585 0.6393 1.936
## y[2] -1.972 -0.7025 -0.02659 0.6463 1.927


##########################################################################################################
## Scenario 2: High correlation, terrible starting values & default block sampler + univariate samplers ##
##########################################################################################################
## Start adaptive MCMC far from target
set.seed(1)
mvg2$y  <- yStart  ## The same really terrible starting values
cmvg2$y <- yStart  ## The same really terrible starting values
nIter   <- 2E5     ## Equivalent to two of the three runs used above
cmcmc2$run(nIter, reset=TRUE)  ## mcmc2$run(nIter, reset=TRUE) 
samples <- as.matrix(cmcmc2$mvSamples); dim(samples)
mc2     <- as.mcmc(samples)
dev.set(dev.list()[1])
plot(mc2)
dev.set(dev.list()[2])
plot(samples[,2],samples[,3], typ="l")
cmvg2$y           ##  0.6248406 0.8214114
cmvg2$calculate() ## -1.10813
## Let's examine the hill climbing a bit closer
plot(log(samples[,1]-samples[1,1]+1), xlab="iteration (i)", ylab="log (logProb[i] - logProb[1] + 1)", typ="l") ## 
## Remove burn-in
lpy   <- samples[,"logProb_y[1]"]
ppy   <- exp(lpy - max(lpy))
ppy   <- ppy/sum(ppy)
(burn <- sum(ppy < 1/nrow(samples)))
mc2   <- as.mcmc(tail(samples[,c("y[1]","y[2]")], nrow(samples)-burn))
plot(mc2)
nrow(mc2)          ## 98136
effectiveSize(mc2) ## 20719.01 20440.75 
##     y[1]     y[2] 
## 20719.01 20440.75 
## 
## summary(mc2)
## Iterations = 1:98136
## Thinning interval = 1 
## Number of chains = 1 
## Sample size per chain = 98136 
## 
## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:
##         Mean     SD Naive SE Time-series SE
## y[1] 0.01447 0.9962 0.003180       0.006921
## y[2] 0.01314 0.9979 0.003185       0.006980
## 
## 2. Quantiles for each variable:
##        2.5%     25%      50%    75% 97.5%
## y[1] -1.910 -0.6636 0.018176 0.6896 1.972
## y[2] -1.914 -0.6688 0.008595 0.6847 1.981

#############################
## Compare the two methods ##
#############################
par(mfrow=n2mfrow(2))
dim(as.matrix(mc)) 
dim(as.matrix(mc2))
plot(as.matrix(mc), pch=19, cex=0.2, main="RW_block")
plot(as.matrix(mc2), pch=19, cex=0.2, main="RW_block_ETA") 

## Heidelberger and Welch Halfwidth test suggests both chains required more samples for 2 decimal point precision in the mean.
heidel.diag(mc,  eps=0.1, pvalue=0.05)
heidel.diag(mc2, eps=0.1, pvalue=0.05) 

##################
## Further runs ##
##################
nIter       <- 1E6
Thin        <- 10
cmcmc$thin  <- Thin
cmcmc2$thin <- Thin
cmcmc$run(nIter, reset=FALSE)  
cmcmc2$run(nIter, reset=FALSE)  

samples <- tail(as.matrix(cmcmc$mvSamples)[,2:3], 1E5)
mc      <- as.mcmc(samples)
dim(mc)
##
samples <- tail(as.matrix(cmcmc2$mvSamples)[,2:3], 1E5)
mc2     <- as.mcmc(samples)
dim(mc2)
## 
dev.set(dev.list()[1])
plot(mc)
dim(mc)
dev.set(dev.list()[2])
plot(mc2)
dim(mc2)

nrow(mc)           ## 100000
nrow(mc2)          ## 100000
effectiveSize(mc)  ## 84705.17 86567.81 
effectiveSize(mc2) ## 93847.49 93054.19 

dev.set(dev.list()[1])
autocorr.plot(mc)  ## Higher because the sampler had fewer iterations at the target to adapt
dev.set(dev.list()[2])
autocorr.plot(mc2) ## Found the target earlier and so the propCov was better adapted than the other sampler

