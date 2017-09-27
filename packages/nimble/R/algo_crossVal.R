crossValCalculate <- function(MCMCOut, realData){
  discrep <- sum(colMeans((MCMCOut - saveData)^2))
  return(discrep)
}

getSD <- function(MCMCout, sampNum, NbootReps, saveData){
  l <- ceiling(min(1000, (sampNum)/20)) ##length of each
  q <- sampNum - l + 1
  h <- ceiling(sampNum/l)
  blockcvValues <- numeric(NbootReps) ## block estimates of cv
  for(r in 1:NbootReps){
    randNum <- runif(h, 0, 1)
    randIndexStart <- ceiling(randNum*q)
    indexes <- cbind(randIndexStart, randIndexStart + l - 1)
    chunks <- c(apply(indexes, 1, function(x){
      seq(x[1], x[2], 1)}
    ))
    blockcvValues[r] <- crossValCalculate(MCMCout[chunks,], saveData)
  }
  return(sd(blockcvValues))
}

calcCrossVal <- function(i,
                         prepModel,
                         leaveOutNames,
                         currentDataNames,
                         mcmcIter,
                         burnInProp,
                         returnSamples,
                         NbootReps){
  saveData <- values(prepModel, leaveOutNames[[i]])
  newModel <- prepModel$newModel()
  newModel$resetData()
  browser()
  values(newModel, leaveOutNames[[i]]) <- NA
  newModel$setData(currentDataNames[[i]])
  compileNimble(newModel)
  modelMCMCConf <- configureMCMC(newModel,
                                 monitors = leaveOutNames[[i]])
  modelMCMC <- buildMCMC(modelMCMCConf)
  C.modelMCMC <- compileNimble(modelMCMC,
                               project = newModel,
                               resetFunctions = (i != 1))
  C.modelMCMC$run(MCMCIter)
  MCMCout <- as.matrix(C.modelMCMC$mvSamples)[,leaveOutNames[[i]]]
  sampNum <- dim(MCMCout)[1]
  startIndex <- max(c(ceiling(burnInProp*MCMCIter), 1))
  message(paste("dropping data block", i))
  crossValAveSD <- getSD(MCMCout=MCMCout, sampNum=sampNum,
                         NbootReps=NbootReps,
                         saveData=saveData)
  crossVal <- crossValCalculate(MCMCout[startIndex:sampNum,], saveData)
  nimble:::clearCompiled(C.modelMCMC)
  return(list(crossValAverage = crossVal,
              crossValAveSD = crossValAveSD,
              samples= if(returnSamples) MCMCout else NA))
}

runCrossValidate <- function(model,
                             dataNames,
                             k,
                             MCMCIter,
                             burnInProp, 
                             MCMCdefs=NULL,
                             returnSamples=FALSE,
                             NbootReps=200){
  

#' Performs k-fold cross-validation for a given NIMBLE model.  
#' 
#' EXPERIMENTAL.  Takes a NIMBLE model and conducts k-fold cross validation on it, returning a measure of the model's predictive performance. 
#' 
#' @param model an uncompiled nimble model 
#' @param k an integer argument determining the number of folds that should be used for cross-validation.
#' @param MCMCIter the number of MCMC iterations to run for each fold of the cross-validation algorithm. 
#' @param burnInProp a number between 0 and 1 specifying the proportion of MCMC iterations that should be discarded from the beginning of each run.
#' @param discrepancyFunctionArgs an R \code{list} containing informaton to be used in the \code{discrepancyFunction}.  See details.
#' @param cpppControl	(optional) an R \code{list} with parameters governing the CPPP calculation.  See details for specific parameters.
#' @param mcmcControl (optional) an R \code{list} with parameters governing the MCMC  algorithm,  See details for specific parameters.
#' @param origMCMCOutput (optional) a matrix with samples from the posterior distribution of the model for all model parameters.
#' Each matrix column should contain samples from a different model parameter, and should be named as that parameter.  This argument should be supplied if
#' a user has already run an MCMC and obtained posterior samples from the model, and wishes to use these samples in the CPPP algorithm. 
#' @param returnSamples a logical argument indicating whether to return all posterior samples from all simulated data sets.  This can result in a very large 
#' returned object,  Defaults to \code{FALSE}.
#' @param nCores an integer argument indicating the number of cpu cores to use for CPPP calculation.  In order to use >1 core, the \code{parallel} R package must be 
#' installed.  Only Mac and Linux operating systems support multiple cores at this time.  
#' 
#' @author Nicholas Michaud and Lauren Ponisio
#' @export
#' @details The CPPP algorithm used is based on that introduced by Hjoort et al. (2006), with some modifications.  Specfically, for each simulated data set and posterior predictive p-value
#' calculated from that data set, a bootstrap estimate of the Monte-Carlo error in the p-value is calculated as well.  These errors are used to produce an error-adjusted calibrated posterior predictive p-value. 
#'
#'  @section The \code{discrepancyFunction} Argument:
#' The \code{discrepancyFunction} argument is a \code{nimbleFunction} that will calculate the discrepancy between the model and the observations.  The choice of discrepancy function is highly model-specific, 
#' and as such no default is provided for this argument.  The \code{discrepancyFunction} provided must follow these rules:
#' \enumerate{
#'   \item The \code{discrepancyFunction} must be a \code{nimbleFunction} with \code{setup} code and {run} code.  See the User's Manual for more information on writing \code{nimbleFunctions}.
#'   \item The \code{setup} code to \code{discrepancyFunction} must accept two arguments: \code{model}, 
#'   which is the model the CPPP is being calculated for, and \code{discrepancyFunctionArgs}, 
#'   which is an R list with any needed information for calculating the discrepancy (e.g. node names, tuning parameters, etc.).
#'   \item The \code{run} code of the discrepancy function must take no arguments and return a scalar double, and must \code{contain} the \code{discrepancyFunction_BASE} virtual \code{nimbleFunction}.
#' }
#' The \code{discrepancyFunctionArgs} argument to \code{runCPPP} is passed to the \code{discrepancyFunction}.  See the examples section below for an example of a discrepancy function.
#' 
#' 
#' @section The \code{cpppControl} Argument:
#' The \code{cpppControl} argument is a list with the following elements:
#'	\itemize{
#'	  \item{\code{nSimVals}}	{
#'      an integer argument determining how many simulated posterior predictive p-values should be calculated to estimate the distribution of p-values.  Defaults to 100.
#'	  }
#'	  \item{\code{nCalcSamples}}{
#'	    an integer argument determining how many posterior samples should be used to calculate the posterior predictive p-value for each simulated data set.  Defaults to 1000.
#'	  }
#'	  \item{\code{nBootstrapSDReps}} {
#'	    an integer argument determining how many bootstrap repititions should be used to estimate the standard error of each posterior predictive p-value.  Defaults to 200.
#'	  }
#'  }
#' @section The \code{mcmcControl} Argument:
#' The \code{mcmcControl} argument is a list with the following elements:
#'	\itemize{
#'	  \item{\code{nMCMCIters}}	{
#'      an integer argument determining how many MCMC iterations should be run for each posterior predictive p-value calculation.  Defaults to 10000, but should probably be manually set.
#'	  }
#'	\item{\code{burnInProp}}{
#'	   the proportion of the MCMC chain to discard as burn-in for each posterior predictive p-value calculation.  Must be between 0 and 1.  Defaults to 0.
#'	}
#'  }
#' 
#' @return
#' an R list with four elements:
#' \itemize{
#' \item{\code{CPPP}} {The CPPP value, measuring the amount of "surprise" the data has given the model.  Lower values indicate greater "surprise".  Will be between 0 and 1.}
#' \item{\code{observed.PPP}} {The posterior predictive p-value for the observed data, along with an estimate of the Monte Carlo error for this p-value.}
#' \item{\code{simulated.PPPs}} {Posterior predictive p-values for each of the simulated data sets, along with error estimates for each p-value.  The set of simulated posterior predictive p-values provide an empirical distribution of p-values, which is used to calculate the CPPP.}
#' \item{\code{samples}} {An R list, only returned when \code{returnSamples = TRUE}.  Each element of this list will be a matrix of posterior samples from the model given a set of simulated data.  There will be \code{nSimVals} sets of samples.}
#' }
#'
#' @references 
#' Hjort, N. L., Dahl, F. A., & Steinbakk, G. H. (2006). Post-processing posterior predictive p values. \emph{Journal of the American Statistical Association}, 101(475), 1157-1174.
#' @examples
#' \dontrun{
#'
#'  
#' pumpCode <- nimbleCode({
#'   for (i in 1:N){
#'     theta[i] ~ dgamma(alpha, beta)
#'     lambda[i] <- theta[i]*t[i]
#'     x[i] ~ dpois(lambda[i])
#'   }
#'   alpha ~ dexp(1.0)
#'   beta ~ dgamma(0.1,1.0)
#' })
#' 
#' pumpConsts <- list(N = 10,
#'                    t = c(94.3, 15.7, 62.9, 126, 5.24,
#'                          31.4, 1.05, 1.05, 2.1, 10.5))
#' 
#' pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
#' 
#' pumpModel <- nimbleModel(code=pumpCode,
#'                          constants=pumpConsts,
#'                          data=pumpData)
#' 
#' # Define the discrepancy function.
#' # Note that the discrepancy function must be a nimbleFunction that takes
#' # model and discFunctionArgs as its only two arguments.
#' 
#' pumpDiscFunction <- nimbleFunction(
#'   setup = function(model, discFunctionArgs){ 
#'     dataNames <- discFunctionArgs[['dataNames']]
#'   },
#'   run = function(){
#'     dataVals <- values(model, dataNames)
#'     lambdaVals <- values(model, 'lambda')
#'     disc <- 0
#'     for(i in 1:dim(dataVals)[1]){
#'       disc <- disc + ((dataVals[i] - lambdaVals[i])/sqrt(lambdaVals[i])) # subtract the expected value, divide by the sd.
#'                                                                          # Note here the expected value and sd are the same,
#'                                                                          # as the data comes from a poisson distribution.
#'     }
#'     returnType(double(0))  # must return a double(0)
#'     return(disc)
#'   },
#'   contains = discrepancyFunction_BASE
#' )
#' 
#' pumpCPPP <- runCPPP(pumpModel,
#'                     dataNames = 'x',
#'                     discrepancyFunction = pumpDiscFunction,
#'                     discrepancyFunctionArgs = list(dataNames = 'x'),
#'                     mcmcControl = list(nMCMCIters = 20000))
#' }


runCrossValidate <- function(model,
                             dataNames,
                             k,
                             MCMCIter,
                             burnInProp, 
                             MCMCdefs=NULL,
                             returnSamples=FALSE,
                             NbootReps=200){
  # simpSetMCMCDefs <- function(Rmodel, MCMCdefs, MCMCname) {
  #   eval(MCMCdefs[[MCMCname]])
  # }
  # if(!is.null(MCMCdefs)){
    ## eval(MCMCdefs.opt2[['nimbleOpt2']])
    ## customSpec <- simpSetMCMCDefs(occ.R.model, MCMCdefs.opt2, 'nimbleOpt2')
  # }
  if(!inherits(model, "RmodelBaseClass")){
    stop("model is not an Rmodel")
  }
  # if(thin == 0 | thin > MCMCIter){
  #   stop("thin needs to be an integer > 0 and < MCMCIter")
  # }
  if(burnInProp >= 1 | burnInProp < 0){
    stop("burnInProp needs to be between 0 and 1")
  }
  testDataNames <- try(model[[dataNames]], silent=TRUE)
  if(inherits(testDataNames, "try-error")){
    stop(paste("dataNames", dataNames,
               "is not the name of the data in model"))
  } else{
    test2DataNames <- all(model$expandNodeNames(dataNames) %in%
                            model$getNodeNames(dataOnly=TRUE))
    if(test2DataNames == FALSE){
      stop(paste("dataNames", dataNames,
                 "is not the name of the data in model"))
    }
  }
  allDataNodeNames <- unique(sapply(dataNames, function(x){return(model$expandNodeNames(nodes = x))}))
  approxNumPerFold <- ceiling(length(allDataNodeNames)/k)
  leaveOutNames <- list()
  currentDataNames <- list()
  reducedDataNodeNames <- allDataNodeNames
  for(i in 1:(k-1)){
    leaveOutNames[[i]] <- sample(reducedDataNodeNames, approxNumPerFold)
    reducedDataNodeNames <- reducedDataNodeNames[-which(reducedDataNodeNames %in% leaveOutNames[[i]])]
    currentDataNames[[i]] <- allDataNodeNames[-which(allDataNodeNames %in% leaveOutNames[[i]])]
  }
  leaveOutNames[[k]] <- reducedDataNodeNames
  currentDataNames[[k]] <- allDataNodeNames[-which(allDataNodeNames %in% leaveOutNames[[k]])]
  
  prepModel <- model$newModel()
  crossValOut <- mclapply(1:k, calcCrossVal,
                          prepModel,
                          leaveOutNames,
                          currentDataNames,
                          mcmcIter,
                          burnInProp,
                          returnSamples,
                          NbootReps)
  crossVal <- sum(sapply(crossValOut, function(x) x$crossValAverage),
                  na.rm=TRUE)
  crossValSD <- sqrt(sum(sapply(crossValOut,
                                function(x) x$crossValAveSD^2),
                         na.rm=TRUE))
  out <- list(crossVal=crossVal,
              crossValSD=crossValSD)
  if(!returnSamples){
    out$samples <- NULL
  }
  else{
    out$samples <- lapply(crossValOut, function(x) x$samples)
  }
  
  return(out)
}
