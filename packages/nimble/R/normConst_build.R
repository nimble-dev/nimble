##  Contains code to approximate normalizing constant.
##  Algorithm found in Section 3 of Gelman and Meng, Statitsical Science Vol. 13, No. 2, 163-185, 1998
##  The getNormConst$run() funtion returns an estimate of the normalizing constant.
##  note the getValsFunc_  functions are used to extract mcmc samples from the mv object and return them.
##  The three getValsFunc_ functions correspond to scalar, vector, and matrix valued variables.

getValsVirtual <- nimbleFunctionVirtual(
  run = function(mcmcSamps = integer()){ 
    returnType(double(2))
  }
)

getValsFunc_Scalar = nimbleFunction(
  contains = getValsVirtual,
  setup = function(var, mcmc){
    mv <- mcmc$mvSamples
  },
  run = function(mcmcSamps = integer()){  
    declare(sampsOut, double(2, c(mcmcSamps,1)))
    for(i in 1:mcmcSamps){
      sampsOut[i,1] <- mv[var,i][1]
    }
    returnType(double(2))
    return(sampsOut)
  }
)
  
getValsFunc_Vector = nimbleFunction(
  contains = getValsVirtual,
  setup = function(var, varLength, mcmc){
    mv <- mcmc$mvSamples
  },
  run = function(mcmcSamps = integer()){
    declare(sampsOut, double(2, c(mcmcSamps,varLength))) 
    for(i in 1:mcmcSamps){
      for(j in 1:varLength){
        sampsOut[i,j] <- mv[var, i][j]
      }
    }
    returnType(double(2))
    return(sampsOut)
  }
)

getValsFunc_Matrix = nimbleFunction(
  contains = getValsVirtual,
  setup = function(var, dim1,dim2, mcmc){
    mv <- mcmc$mvSamples
    varLength = dim1*dim2
  },
  run = function(mcmcSamps = integer()){
    declare(sampsOut, double(2, c(mcmcSamps,varLength)))
    
    for(i in 1:mcmcSamps){
      for(row in 1:dim1){
        for(col in 1:dim2){
          sampsOut[i, (row-1)*dim2+col] <- mv[var, i][row,col]
        }
      }
    }
    returnType(double(2))
    return(sampsOut)
  }
)

impMV_simFunc = nimbleFunction(
  setup = function(model, nodes){
    paramDetermNode <- model$getDependencies(nodes, determOnly = TRUE) 
    areDetermNodes <- length(paramDetermNode) > 0
  },
  run = function(normMean = double(1), normCov = double(2)){     
    ##  for each sample, calculate likelihod of model|sample and MVN density of sample to create norm constant
    tempSamp <- rmnorm_chol(1, normMean, normCov, prec_param = 0)
    values(model, nodes) <<- tempSamp      
    if(areDetermNodes){
      calculate(model, paramDetermNode) 
    }
    numer <-exp(calculate(model))
    if(is.nan(numer)) numer <- 0
    denom <- dmnorm_chol(tempSamp, normMean, normCov, prec_param = 0)
    single_Cest <- numer/denom
    returnType(double())
    return(single_Cest)
  }
)

impScalar_simFunc = nimbleFunction(
  setup = function(model, nodes){
    paramDetermNode <- model$getDependencies(nodes, determOnly = TRUE) 
    areDetermNodes <- length(paramDetermNode) > 0
  },
  run = function(normMean = double(1), normCov = double(2)){     
    ##  for each sample, calculate likelihod of model,sample and normal density of sample to create norm constant
    tempSamp <- rnorm(1, normMean[1], normCov[1,1])
    model[[nodes]] <<- tempSamp      
    if(areDetermNodes){
      calculate(model, paramDetermNode) 
    }
    numer <-exp(calculate(model))
    if(is.nan(numer)) numer <- 0
    denom <- dnorm(tempSamp, normMean[1], normCov[1,1])
    single_Cest <- numer/denom
    returnType(double())
    return(single_Cest)
  }
)

##  A function which calculates the mean and covariance matrix of the mcmc samples for all stochastic,
##  non-data nodes in the model, and uses these to create a normal approximation to the posterior distribution of these
##  nodes.  Samples are taken from this normal approximation and importance sampling is used to estimate the normalizing
##  constant.
normImpSamp = nimbleFunction(
  setup = function(model, nodes, totalLength){
    paramDetermNode <- model$getDependencies(nodes, determOnly = TRUE) 
    areDetermNodes <- length(paramDetermNode) > 0

    ## different methods depending on whether there are one or more nodes that were sampled 
    if(totalLength>1)
      impSimFunc <- impMV_simFunc(model, nodes)
    else
      impSimFunc <- impScalar_simFunc(model, nodes)
  },
  run = function(vals = double(2), nSamps = integer(0),  nImpSamps = integer(0)){
    declare(normMean, double(1, totalLength))  # observed mean
    declare(normCov, double(2, c(totalLength, totalLength)))  # observed covariance matrix
    declare(normCEst, double(1, nImpSamps))  # vector of normalizing constant estimates
    declare(nSampsd, double())  # nSamps but in double type
    declare(out, double(1, 2))  # output
    nSampsd <- nSamps

    ## Calculate empirical mean and covariance
    normMult <- 1/(nSampsd-1)

    if(totalLength > 1){
      for(i in 1:totalLength){
        normMean[i] <- mean(vals[,i])
      }
      for(j in 1:nSamps){
        normCov <- normCov + (vals[j,]-normMean)%*%t(vals[j,]-normMean)
      }
      normCov <- normMult*normCov
      normCov <- chol(normCov)
    }
    else{
      normMean[1] <- mean(vals[,1])
      normCov[1,1] <- sd(vals[,1])^2
      normCov <- normMult*normCov
    }
    
    ##  generate nImpSamps random samples from MVN distribtuion
    ##  for each sample, calculate likelihod of model,sample and MVN density of sample to create norm constant
    for(k in 1:nImpSamps){
      normCEst[k] <- impSimFunc$run(normMean, normCov)
    }
    L <- mean(normCEst)
    varLL <- var(normCEst) / nImpSamps / L^2
    if(is.nan(varLL)) { varLL <- Inf
                        return(-Inf)  }
    return(log(L))
    out[1] <- log(L)
    out[2] <- sd(varLL)  # also return var. of log importance weights
    returnType(double(1))
    return(out)
  }
)
    
                                  
#' Calculates a normalizing constant for a NIMBLE model using a MVN approximation to the posterior.
#' 
#' Takes a nimble model and integrates over all stochastic, non-data nodes. The run() method returns an approximation of log of the normalizing constant.
#' The getVar() method returns the standard deviation of the importance weights used to approximate the normalizing constant.  
#' 
#' @param model 			A nimble model 
#' @param burnIn			burn-in used for MCMC sampler
#' @param mcmcControl		list passed to \code{MCMCSpec}, a nimble function that builds the MCMC sampler. See \code{help(MCMCSpec)} for more details
#' @author Nick Michaud
#' @export
#' @details getNormConst returns an importance sampling approximation to the log of the normalizing constant of a model (f(data|model)), using a method
#' detailed in Gelman and Meng, Statitsical Science Vol. 13, No. 2, 163-185, 1998.  The posterior distribution of all stochastic nodes is approximated 
#' by a multivariate normal distribution.  In cases where this approximation is not appropriate, the estimate of the normalizing constant will not be reliable.
#' An idea of the Monte Carlo errror in the approximation can be gotten by the getVar() function, which will return the  variance of the log of the 
#' importance sampling weights.  
#'
#' See user manual for more 
#' @section Runtime Arguments:
#'	\itemize{
#'	\item{\code{mcmcSamps}}	{
#'    number of iterations to run the MCMC algorithm for.
#'	}
#'  \item{\code{importanceSamps}}	{
#'    number of importance samples to use.
#'	}
#'
#' @examples
#' 
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
#' compileNimble(pumpModel)
#' getPumpNormConst <- getNormConst(pump, burnIn = 1000, mcmcControl = list(adaptInterval = 100))
#' cgetPumpNormConst <- compileNimble(getPumpNormConst, project = pumpModel)
#' normConstEst <- cgetPumpNormConst$run(10000,10000)
#' normConstCV <- cgetPumpNormConst$getCV()
#'
getNormConst <-nimbleFunction(
  setup = function(model, burnIn = 1000, mcmcControl = list(adaptInterval = 100)){
    
    ## Get nodes to integrate over
    intNodes <- model$getNodeNames(stochOnly = T, includeData = FALSE)
    if(length(intNodes) == 0)
      stop("no stochastic nodes")
    
    intVars <- model$getVarNames(nodes =  intNodes)  # need var names too
    
    my_initializeModel <- initializeModel(model)
    
    
    mcmc_Latent_Spec <- configureMCMC(model, nodes = intNodes, monitors = intNodes, control = mcmcControl) 
    mcmc_Latent <- buildMCMC(mcmc_Latent_Spec)
    
    numVars <- length(intVars)
    varLengths <- unname(sapply(intVars, function(x){return(length(model[[x]]))}))
    totalLength <- sum(varLengths)
    varDims <-lapply(intVars, function(x){return(dim(model[[x]]))})
    impSampFunc <- normImpSamp(model, intNodes, totalLength) 
    mv <- mcmc_Latent$mvSamples

    ## a list of functions (one for each nimble variable) which will extract and return
    ## the values from the variable.  Different functions depending on whether variable is
    ## a scalar, vector, or matrix
    getValsFuncList <- nimbleFunctionList(getValsVirtual) 
    if(numVars == 1){
      varLengths <- c(varLengths, 0)
    }

    for(var in 1:numVars){
      if(varLengths[var]==1){
        getValsFuncList[[var]] <- getValsFunc_Scalar(intVars[var], mcmc_Latent)
      }
      else{
        if(is.na(varDims[[var]][1])) varDims[[var]][1] <- 1
        if(length(varDims[[var]])==1){
          getValsFuncList[[var]] <- getValsFunc_Vector(intVars[var], varLengths[var], mcmc_Latent)
        }
        else{
          getValsFuncList[[var]] <- getValsFunc_Matrix(intVars[var], varDims[[var]][1],varDims[[var]][2], mcmc_Latent)
        } 
      }
    }
    ## initial value for variance of log importance sampling weights.  if the getVar() method
    ## is called before the normalizing constant has been estimated, a variance of -1 will be returned.
    variance <- -1
  },
  run = function(mcmcSamps = integer(0), importanceSamps = integer(0)){
    declare(vars, double(2,c(mcmcSamps, totalLength)))  # vars will hold the mcmc samples for all of our variables 
    numSamps <- mcmcSamps - burnIn  # the number of mcmc samples to actually use
    stVals <- values(model, intNodes)  # save the starting values of the model
    if(burnIn >= mcmcSamps)
      stop("Number of MCMC samples less than burnIn, MCMC quitting")
    if(mcmcSamps < 2)
      stop("Must have at least 2 MCMC samples (although more are recommended)")
    if(importanceSamps < 2)
      stop("Must have at least 2 importance samples (although more are recommended)")
    
    ##  initialize model and run mcmc
    my_initializeModel$run() 
    mcmc_Latent$run(mcmcSamps) 
  
    ##  iterate through the variables, extracting the mcmc samples for each one and storing in the vars object.
    colInd <- 1
    for(var in 1:numVars){
      vars[,colInd:(colInd+varLengths[var]-1)] <- getValsFuncList[[var]]$run(mcmcSamps)
      colInd <- colInd + varLengths[var]
    }
    
    ## use importance sampling to estimate normalizing constant and variance of importance weights
    impFuncOutput <- impSampFunc$run(vars, numSamps, importanceSamps)
    normConst <- impFuncOutput[1]
    variance <<- impFuncOutput[2]/exp(normConst)
    
    values(model, intNodes) <<- stVals 
    
    returnType(double())
    return(normConst)
  },
  
  ## A method to return the variance  of the log importance sampling weights.
  ## A high variance value may indicate that more mcmc or importance weights should be used, 
  ## or that the normal approximation to the posterior distribtuion is not appropriate.
  methods = list(
    getVar = function(){
      returnType(double())
      return(variance)
    }
  ),where = getLoadingNamespace()
)

