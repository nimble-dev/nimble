dataCalcVirtual <- nimbleFunctionVirtual(
  run = function() 
    returnType(double())
)

calcDataNF <- nimbleFunction(
  contains = dataCalcVirtual,
  setup = function(model, node){},
  run = function(){
    returnType(double())
    return(model$calculate(node))
  }
)


#' Calculates WAIC from MCMC samples
#' 
#' Takes a NIMBLE model and a NIMBLE MCMC object.  Creates an algorithm that will calculate the WAIC from the posterior samples of the MCMC object.
#' 
#' @param model a NIMBLE model 
#' @param mcmc a NIMBLE MCMC object, created from a call to buildMCMC().  
#'
#'  @author Nicholas Michaud
#' @export
#' @details 
#' Calling calcWAIC will produce an uncompiled (R) function object, say \code{Rwaic}.
#' The  WAIC function has a single arugment:
#'
#' \code{burnin}: The number of iterations to subtract from the beginning of the posterior samples of the MCMC object.
#' 
#' The uncompiled (R) WAIC function may be compiled to a compiled WAIC object, taking care to compile in the same project as both the R model object and the MCMC object that were used as arguments to \code{calcWAIC()}, using: \code{Cwaic <- compileNimble(Rwaic, project=Rmodel)}

#' The compiled function will function identically to the uncompiled object, except acting on the compiled model object and compiled MCMC samples.

#' The \code{run} method of \code{calcWAIC} calculates the WAIC of a given model.  The WAIC is calculated using the posterior samples from the provided NIMBLE MCMC object.  As such, the MCMC algorithm must have been run prior to calculating the WAIC.     
#' The WAIC (Watanabe, 2010) is calculated from Equations 5, 12, and 13 in Gelman (2014).  Note that the set of all parameters monitored by \code{mcmc} will be treated as \eqn{theta} for the purposes of e.g. Equation 5 from Gelman (2014). 
#'  All parameters downstream of the monitored parameters that are necessary to calculate \eqn{p(y|theta)} will be simulated from the posterior samples of \eqn{theta}.
#'
#' @references 
#' 
#' Gelman, Andrew, Jessica Hwang, and Aki Vehtari (2014). Understanding predictive information criteria for Bayesian models. \emph{Statistics and Computing} 24(6), 997-1016.
#'  
#'  Watanabe, Sumio (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. \emph{Journal of Machine Learning Research} 11, 3571-3594.
#' @examples
#' \dontrun{
#' pumpCode <- nimbleCode({ 
#'  for (i in 1:N){
#'      theta[i] ~ dgamma(alpha,beta);
#'      lambda[i] <- theta[i]*t[i];
#'      x[i] ~ dpois(lambda[i])
#'  }
#'  alpha ~ dexp(1.0);
#'  beta ~ dgamma(0.1,1.0);
#' })
#'
#' pumpConsts <- list(N = 10,
#'               t = c(94.3, 15.7, 62.9, 126, 5.24,
#'                 31.4, 1.05, 1.05, 2.1, 10.5))
#'
#' pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
#'
#' pumpInits <- list(alpha = 1, beta = 1,
#'              theta = rep(0.1, pumpConsts$N))
#' pumpModel <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
#'                   data = pumpData, inits = pumpInits)
#'
#' Cpump <- compileNimble(pump)
#' pumpConf <- configureMCMC(pump, print = TRUE)
#' pumpMCMC <- buildMCMC(pumpConf)
#' CpumpMCMC <- compileNimble(pumpMCMC, project = pump)
#' CpumpMCMC$run(niter = 1000)
#' pumpWAIC <- calcWAIC(pumpModel, pumpMCMC)
#' CpumpWAIC <- compileNimble(pumpWAIC, project = pump)
#' WAIC <- CpumpWAIC$run(burnIn = 500)
#' }

calcWAIC <- nimbleFunction(
  setup = function(model, mcmc){
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    dataCalcNFList <- nimbleFunctionList(dataCalcVirtual)
    for(i in seq_along(dataNodes)){
      dataCalcNFList[[i]] <- calcDataNF(model, dataNodes[i])
    }
    dataNodeLength <- length(dataNodes)
    sampledNodes <- model$getVarNames(FALSE, colnames(as.matrix(mcmc$mvSamples)))
    paramDeps <- model$getDependencies(sampledNodes, self = FALSE)
    mvSamples <- mcmc$mvSamples
  },
  run = function(burnIn = integer()){
    numMCMCSamples <- getsize(mvSamples) - burnIn
    if((numMCMCSamples) < 2){
      print('Error, need more than 1 post burn-in MCMC sample.')
      return(-Inf)
    }
    logPredProbs <- numeric(numMCMCSamples)
    logAvgProb <- 0
    pWAIC <- 0
    currentVals <- values(model, sampledNodes)
    for(i in 1:dataNodeLength){
      meanPredProb <- 0
      for(j in 1:numMCMCSamples){
        copy(mvSamples, model, nodesTo = sampledNodes, row = j + burnIn)
        model$simulate(paramDeps)
        logPredProbs[j] <- dataCalcNFList[[i]]$run()
      }
      meanPredProb <- mean(exp(logPredProbs))
      pointLogPredVar <- sd(logPredProbs)^2
      pWAIC <- pWAIC + pointLogPredVar
      logAvgProb <- logAvgProb + log(meanPredProb)
    }
    WAIC <- -2*(logAvgProb- pWAIC)
    values(model, sampledNodes) <<- currentVals
    model$calculate(paramDeps)
    returnType(double())
    return(WAIC)
  }, where = getLoadingNamespace()
)
