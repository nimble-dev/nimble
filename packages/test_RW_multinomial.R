########################################################
## Missing data problem with multinomial distribution ##
########################################################

## rm(list=ls())

library(nimble)
library(coda)

## source("~/nimbleProject/multinomial/sampler.R"); source("~/nimbleProject/multinomial/test_RW_multinomial.R")

## For accessing "persistent member variables"
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE) 
nimbleOptions()

############################
## Model in BUGS language ## 
############################
codeTest <- nimbleCode ({
    X[1:nGroups] ~ dmultinom(size=N, prob=pVecX[1:nGroups])
    Y[1:nGroups] ~ dmultinom(size=N, prob=pVecY[1:nGroups])
    for (ii in 1:nGroups) {
        Z[ii] ~ dbeta(1 + X[ii], 1 + Y[ii])
    }
})

## CONSTANTS, INITS & DATA 
nGroups   <- 5
N         <- 1E6
pVecX     <- rdirch(1, rep(1, nGroups)) ## rep(1/nGroups, nGroups)
pVecY     <- rdirch(1, rep(1, nGroups)) ## rep(1/nGroups, nGroups)
X         <- rmultinom(1, N, pVecX)[,1]
Y         <- rmultinom(1, N, pVecY)[,1]
Z         <- rbeta(nGroups, 1+X, 1+Y)
                                        # Resample X and Y to make the MCMC more challenging
Xini      <- rmultinom(1, N, sample(pVecX))[,1]
Yini      <- rmultinom(1, N, sample(pVecY))[,1]
## Lists for nimbleModel
Constants <- list(nGroups=nGroups)                     ## Can't modify after compilation
Inits     <- list(X=Xini, Y=Yini, pVecX=pVecX, pVecY=pVecY, N=N) ## Can modify after compilation
Data      <- list(Z=Z)                                           ## Can modify after compilation


## ASSEMBLE NIMBLE MODEL
modelTest <- nimbleModel(codeTest, constants=Constants, inits=Inits, data=Data, check=TRUE) ## FALSE

## EXAMINE MODEL
modelTest$getNodeNames()
modelTest$X
modelTest$Y
modelTest$Z
if (FALSE) ## TRUE
    plot(modelTest$graph)

## TO COMPILE THE MCMC YOU MUST FIRST DO THIS...
cModelTest <- compileNimble(modelTest)

## CONFIGURE MCMC
mcmcTestConfig <- configureMCMC(cModelTest, print=TRUE) ## control=list(adaptInterval=1000)
mcmcTestConfig$getMonitors()
mcmcTestConfig$printSamplers()
## [1] RW_block sampler: X[1:5] ## THIS CHOICE OF SAMPLER ISN'T GOING TO WORK
## [2] RW_block sampler: Y[1:5] ## THIS CHOICE OF SAMPLER ISN'T GOING TO WORK
mcmcTestConfig$removeSamplers()
mcmcTestConfig$printSamplers()
## Add our custom-built random walk sampler on node 'x', with a fixed proposal standard deviation = 0.1
mcmcTestConfig$removeSamplers()
mcmcTestConfig$addSampler(target = 'X', type = 'RW_multinomial', control = list(adaptive=TRUE, adaptInterval=100)) 
mcmcTestConfig$addSampler(target = 'Y', type = 'RW_multinomial', control = list(adaptive=TRUE, adaptInterval=100)) 
mcmcTestConfig$printSamplers()
## BUILD & COMPILE MCMC
mcmcTest  <- buildMCMC(mcmcTestConfig) 
cMcmcTest <- compileNimble(mcmcTest, project=modelTest)

## Optionally resample data
cModelTest$N      <- N <- 1E3
(cModelTest$pVecX <- sort(rdirch(1, rep(1, nGroups)))) 
(cModelTest$pVecY <- sort(rdirch(1, rep(1, nGroups)))) 
simulate(cModelTest, "X", includeData=TRUE); (X <- cModelTest$X)
simulate(cModelTest, "Y", includeData=TRUE); (Y <- cModelTest$Y)
simulate(cModelTest, "Z", includeData=TRUE); (Z <- cModelTest$Z)

###########################
## Reset Adative Sampler ##
###########################
RESET <- TRUE

##############
## Run MCMC ##
##############
niter  <- 1E4
print(sysT <- system.time(cMcmcTest$run(niter, reset = RESET))) 
Tail   <- tail(samples <- as.matrix(cMcmcTest$mvSamples), 25)
if (nrow(samples)>niter)
    samples <- samples[-(1:niter),]
dim(as.matrix(cMcmcTest$mvSamples))
RESET  <- FALSE ## TRUE
iColsX <- 1:nGroups
iColsY <- iColsX + nGroups
##############
## ACF Plot ##
##############
if (FALSE) { ## TRUE ## 
    par(mfrow=c(1,1))
    mcSamples <- as.mcmc(samples)
    acfplot(mcSamples)
}
##################################
## Trajectory Plots & Histogram ##
##################################
plotHistograms <- N <= 1E4 ## FALSE ## TRUE
par(mfrow=c(2,1+plotHistograms))
##
for (ii in 1:2) {
    if (ii == 1) {
        yMaxX <- 0
        yMaxY <- 0
    }
    ##
    plot (samples[,1],ylim=range(samples[,iColsX]), typ="n")
    for (ii in iColsX)
        lines(samples[,ii], col=rainbow(10,alpha=0.75)[ii])
    ##
    if (plotHistograms) {
        hist(samples[,1], col=rainbow(2*nGroups, alpha=0.1)[1], breaks=min(samples):max(samples), prob=TRUE, ylim=c(0,yMaxX))
        for (ii in iColsX) {
            h <- hist(samples[,ii], prob=TRUE, 
                      col=rainbow(2*nGroups, alpha=0.1)[ii],
                      border=rainbow(2*nGroups, alpha=0.1)[ii],
                      breaks=min(samples[,iColsX]):max(samples[,iColsX]), add=TRUE)
            yMaxX <- max(yMaxX, h$density)
        }
    }
    ## 
    plot (samples[,1],ylim=range(samples[,iColsY]), typ="n")
    for (ii in iColsY)
        lines(samples[,ii], col=rainbow(10,alpha=0.75)[ii])
    ##
    if (plotHistograms) {
        hist(samples[,1+nGroups], col=rainbow(2*nGroups, alpha=0.1)[1], breaks=min(samples[,iColsY]):max(samples[,iColsY]), prob=TRUE, ylim=c(0,yMaxY))
        for (ii in iColsY) {
            h <- hist(samples[,ii], prob=TRUE, 
                      col=rainbow(2*nGroups, alpha=0.1)[ii],
                      border=rainbow(2*nGroups, alpha=0.1)[ii],
                      breaks=min(samples[,iColsY]):max(samples[,iColsY]), add=TRUE)
            yMaxY <- max(yMaxY, h$density)
        }
    }
}
print(Tail)

cMcmcTest$samplerFunctions$contentsList[[1]]$AcceptRates
cMcmcTest$samplerFunctions$contentsList[[1]]$timesAccepted /
    cMcmcTest$samplerFunctions$contentsList[[1]]$timesRan ## Some NaNs here are normal.
cMcmcTest$samplerFunctions$contentsList[[1]]$RescaleThreshold
cMcmcTest$samplerFunctions$contentsList[[1]]$ScaleShifts
cMcmcTest$samplerFunctions$contentsList[[1]]$ENSwapDeltaMatrix
cMcmcTest$samplerFunctions$contentsList[[1]]$totalAdapted
cMcmcTest$samplerFunctions$contentsList[[1]]$ENSwapDeltaMatrix /
    cMcmcTest$samplerFunctions$contentsList[[1]]$totalAdapted
cMcmcTest$samplerFunctions$contentsList[[1]]$ENSwapMatrix
cMcmcTest$samplerFunctions$contentsList[[1]]$totalAdapted
cMcmcTest$samplerFunctions$contentsList[[1]]$ScaleShifts
cMcmcTest$samplerFunctions$contentsList[[1]]$timesAccepted 
cMcmcTest$samplerFunctions$contentsList[[1]]$timesRan 

cMcmcTest$samplerFunctions$contentsList[[2]]$AcceptRates
cMcmcTest$samplerFunctions$contentsList[[2]]$timesAccepted /
    cMcmcTest$samplerFunctions$contentsList[[2]]$timesRan  ## Some NaNs here are normal.
cMcmcTest$samplerFunctions$contentsList[[2]]$RescaleThreshold
cMcmcTest$samplerFunctions$contentsList[[2]]$ScaleShifts
cMcmcTest$samplerFunctions$contentsList[[2]]$ENSwapDeltaMatrix
cMcmcTest$samplerFunctions$contentsList[[2]]$totalAdapted
cMcmcTest$samplerFunctions$contentsList[[2]]$ENSwapDeltaMatrix /
    cMcmcTest$samplerFunctions$contentsList[[2]]$totalAdapted
cMcmcTest$samplerFunctions$contentsList[[2]]$ENSwapMatrix
cMcmcTest$samplerFunctions$contentsList[[2]]$totalAdapted
cMcmcTest$samplerFunctions$contentsList[[2]]$ScaleShifts

