### Load packages
library(nimble)
library(MASS)
library(KernSmooth)
library(SamplerCompare)
library(mcmc)
library(graphics)
library(coda)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
####################################################################
### Parameter for comparisons                                                                                            ############
####################################################################
Nreplicate = 20
Scales = c(0.1,0.2,0.3,0.4)
Methods = c("coda", "MBB", "ics", "ar.act", "discrete","matrix", "d-ics", "d-ar.act")
Gridsizes = seq(90,110, by =10)
Bandwidths = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
Samplers = c("RW", "slice sampler")
RunSamplers = c("SNS", "CNS") # single no state, combine no state, "single with state", "combine with state")

####################################################################
### Create a directory to store results                                                                                  ############
####################################################################

mainDir <- "./"
subDir <- "results"

if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
} else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))

}


### Set up model
Code <- nimbleCode({
  sigma ~ dunif(0, 100)
  for (i in 1:n)
    x[i] ~ dnorm(0, sd=sigma)
  mu ~ dnorm(0, sd = 1000)
  for (i in 1:n) 
    for (j in 1:m) 
      Y[i,j] ~ dnorm(mu + x[i], sd = sigma)
})

Consts <- list(n = 4, m=4)

Data <- list(Y = matrix(c(5, 1, 5, 14,19, 1, 1, 4, 11,12,13,15, 14, 16, 22, 3), ncol=4))

Inits <- list(sigma = 1, mu = 0)

Rmodel <- nimbleModel(code = Code, name = 'RandomEffect', constants = Consts, data = Data, inits = Inits)

### Make RW_record sampler
####################################################################
### scalar RW sampler with normal proposal distribution ############
####################################################################

RW_record <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ## these lines are new:
    numSamples <- 0
    before <- c(0, 0)
    after <- c(0, 0)
    if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
    if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
    ###  control list extraction  ###
    logScale      <- control$log
    reflective    <- control$reflective
    adaptive      <- control$adaptive
    adaptInterval <- control$adaptInterval
    scale         <- control$scale
    if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    ###  node list generation  ###
    calcNodes  <- model$getDependencies(target)
    ###  numeric value generation  ###
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory          <- c(0, 0)
    acceptanceRateHistory <- c(0, 0)
    optimalAR <- 0.44
    gamma1    <- 0
    range <- getDistribution(model$getNodeDistribution(target))$range
  },
  
  run = function() {
    ## these lines are new:
    numSamples <<- numSamples + 1
    setSize(before, numSamples)
    setSize(after, numSamples)
    before[numSamples] <<- model[[target]]
    ## back to the original sampler function code:
    currentValue <- model[[target]]
    if(!logScale)    propValue <-     rnorm(1, mean =     currentValue,  sd = scale)
    else             propValue <- exp(rnorm(1, mean = log(currentValue), sd = scale))
    if(reflective)   while(propValue < range[1] | propValue > range[2]) {
      if(propValue < range[1]) propValue <- 2*range[1] - propValue
      if(propValue > range[2]) propValue <- 2*range[2] - propValue    }
    model[[target]] <<- propValue
    logMHR <- calculateDiff(model, calcNodes)
    if(logScale)     logMHR <- logMHR + log(propValue) - log(currentValue)
    jump <- decide(logMHR)
    if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    if(adaptive)     adaptiveProcedure(jump)
    ## this line new:
    after[numSamples] <<- model[[target]]
  },
  
  methods = list(
    
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        setSize(scaleHistory,          timesAdapted)
        setSize(acceptanceRateHistory, timesAdapted)
        scaleHistory[timesAdapted] <<- scale
        acceptanceRateHistory[timesAdapted] <<- acceptanceRate
        gamma1 <<- 1/((timesAdapted + 3)^0.8)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      scaleHistory          <<- scaleHistory          * 0
      acceptanceRateHistory <<- acceptanceRateHistory * 0
      gamma1 <<- 0
    }
  ), where = getLoadingNamespace()
)

####################################################################
### Perry's moving block bootstrap code `               ############
####################################################################

MBB <- function(sample, blockLength = 100, m = 100, FUN, ..., doRbind = TRUE) {
  n <- length(sample)
  numBlocks <- n / blockLength
  if(numBlocks != round(numBlocks)) stop('Must give a blockLength that divides into length(sample) with no remainder.')
  ans <- vector('list', m)
  for(iBS in 1:m) {
    blockStarts <- sample(1:(n - blockLength + 1), size = numBlocks, replace = TRUE)
    indices <- as.numeric(outer(blockStarts, 0:(blockLength-1), '+'))
    ans[[iBS]] <- FUN(sample[indices], ...)
  }
  if(doRbind) return(do.call('rbind', ans)) else return(ans)
}
### Make results to store output for comparison
results<-matrix(rep(0,160),ncol=8)
for(s4 in Bandwidths){
for(s3 in Gridsizes){
for(s2 in Scales){

### Try one sampler (S1) alone

spec <- configureMCMC(Rmodel, nodes = NULL)
## Oops, this means we have not been sampling X.  I will do it that way for now so that I can replicate the results Dao was getting.
spec$printSamplers() ## It's empty, since we set nodes=NULL
spec$addSampler(target = 'sigma', type = 'RW_record', control = list(adaptive = FALSE, scale = s2)) ## S1
spec$addSampler(target = 'mu', type = 'RW_record', control = list(adaptive = FALSE, scale = .1)) ## with X's fixed, mu is conditionally independent from sigma so this is doing nothing in relation to sigma!
spec$printSamplers() ## Now we see the samplers

spec$addMonitors(c('sigma','mu'))
Rmcmc <- buildMCMC(spec)

#### compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Nchain = 10000

for (s1 in 1:Nreplicate){
set.seed(s1)
#mcmc$run(Nchain, simulateAll=TRUE)
Cmcmc$run(Nchain)

### ESS and spec

sigmaSamples <- as.matrix(Cmcmc$mvSamples)[,'sigma']
essSigma <- effectiveSize(sigmaSamples)
essSigma

results[s1,1]=spectrum0.ar(sigmaSamples)$spec ## since order estimated at < 50, no need for "big" version


var(sigmaSamples) / Nchain  ## naive standard error for mean(sigma)
var(sigmaSamples) / essSigma  ## actual standard error for mean(sigma)
var(sigmaSamples) * (Nchain/essSigma) ## spec

### spec from moving block bootstrap
blockLength = 100
test <- MBB(sigmaSamples, blockLength = blockLength, m =200, FUN = mean) ## create 1000 MBB samples from blocks of length 200 and estimate mean of each MBB sample.  blockLength must be large enough that autocorrelation with lag=blockLength is approximately 0
effectiveSize(sigmaSamples) ## usual estimator
var(sigmaSamples) / var(test)  ## estimator of ESS based on MBB
results[s1,2]=var(sigmaSamples)^2 / var(test)  ## estimator of spec based on MBB

### ESS from Initial sequence estimators
tau <-ics.act(sigmaSamples)
results[s1,3]=var(sigmaSamples)*tau
### ESS from AR process
tau <- ar.act(sigmaSamples)$act  # from Thompson 
results[s1,4]=var(sigmaSamples)*tau  ## estimator of IACT 

### Kernel
N <- s3
beforeSamples <- Cmcmc$samplerFunctions$contentsList[[1]]$before
afterSamples <- Cmcmc$samplerFunctions$contentsList[[1]]$after

x = cbind(beforeSamples, afterSamples)
est <- bkde2D(x, bandwidth=c(s4, s4), gridsize = c(N,N), range.x = list(c(0,12), c(0,12)))
#contour(est$x1, est$x2, est$fhat)
K <- est$fhat

head(K)

#contour(est$x1, est$x2, est$fhat)
#points(beforeSamples, afterSamples, pch = '.')
K <- t(K) 
for(j in 1:N) K[,j] <- K[,j] / sum(K[,j])
colSums(K) 
## Now I think K[i,j] is probability of going to x[i] from x[j].
e = eigen(K)
firste=e$vector[,1]
firste  ##first eigen vector of K
firste = Re(firste)
if(sum(firste)<0){
  firste = -firste
}

firste[which(firste<0.0001)]=0
firste = firste/sum(firste)
firste

## Now I think K[i,j] is probability of going to x[i] from x[j].

### Generating a chain for the discrete kernel
x <- est$x1
discreteChainIndices <- integer(Nchain)
discreteChainIndices[1] <- round(length(x)/2)
for(i in 2:Nchain) {
discreteChainIndices[i] <- sample(1:N, 1,
prob = K[,discreteChainIndices[i-1]])
}
discreteChain <- x[discreteChainIndices]
#plot(discreteChain,pch='.')
#plot(beforeSamples, afterSamples, pch='.')
#points(jitter(discreteChain[1:(Nchain-1)]), jitter(discreteChain[2:Nchain]), col = 'red')

effectiveSize(discreteChain) 
results[s1,5] =spectrum0.ar(discreteChain)$spec

### ESS and spec from the discrete chain

#plot(x, firste)
m=sum(x * firste)
m
Var1 = firste*(x-m)*(x-m)
V= sum(Var1)
V
tau = 1
a <- K
a1=K
b=a1
## I note that x taken from the kernel grid above does not match d$mids
for(lag in 1:500){ ## only going up to 300 because contributions should vanish
  for (j in 1:N){
    ## in case of sum of column = 0
    if(sum(a[,j])==0){
      a[,j] = rep(1/N,N)
    }
    b[,j]=a[,j]/sum(a[,j])*firste[j] #
  }
  ## We want sum_{ij} P(x[j]) K^lag[i,j] (x[i] - m) (x[j] - m)
  ## Here a represents K^lag
  ## d$density[j] represents P(x[j])
  ## b[i,j] represents K^lag[i,j] * P(x[j])
  ## t(x-m) %*% b %*% (x-m) represents the sum above
  
  ## treating (M-j)/M as approx 1
  thisContribution <- ((t(x-m)%*% (b %*%(x-m)))[1,1]) / V
  cat("autocorrelation for lag ", lag, " = ", thisContribution, "\n")
  tau= tau+ 2*thisContribution
  a=a%*%a1
}
tau

results[s1,6] = tau * var(discreteChain)

### ICS for disrete chain:

results[s1,7] = ics.act(discreteChain)* var(discreteChain)

### AR for disrete chain:

results[s1,8] =ar.act(discreteChain)$act * var(discreteChain)

}
nameresult = paste("RW_","sc_",s2,"gs_",s3,"bw_",s4,".pdf")
pdf(nameresult)
boxplot(results, names=Methods, main="Boxplot of spectrum at 0 from different methods")
dev.off()
rm(list=c("spec","Rmcmc","Cmcmc", "Cmodel"))
}
}
}