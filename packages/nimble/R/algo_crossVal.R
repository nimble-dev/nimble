calcCrossValSD <- function(MCMCout, sampNum, NbootReps, saveData, lossFunction, predLoss){
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
    blockcvValues[r] <- mean(apply(MCMCout[chunks,], 1, lossFunction, saveData))
    if(predLoss){
      blockcvValues[r] <- log(blockcvValues[r])
    }
  }
  return(sd(blockcvValues))
}


calcCrossVal <- function(i,
                         conf,
                         foldFunction,
                         lossFunction,
                         MCMCiter,
                         burnInProp,
                         returnSamples,
                         NbootReps){
  model <- conf$model
  leaveOutNames <- model$expandNodeNames(foldFunction(i))
  currentDataNames <- model$getNodeNames(dataOnly = T)
  currentDataNames <- currentDataNames[!(currentDataNames %in% leaveOutNames)]
  saveData <- values(model, leaveOutNames)
  newModel <- model$newModel(check = FALSE)
  newModel$resetData()
  values(newModel, leaveOutNames) <- NA
  newModel$setData(model$getVarNames(nodes = currentDataNames))
  compileNimble(newModel)
  predLoss <- FALSE
  if(lossFunction == 'predictive'){
    paramNames <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    dependencies <- model$getDependencies(paramNames)
    missingDataNames <- leaveOutNames
    lossFunction <- function(posteriorSamples, observedData){
      values(model, paramNames) <- posteriorSamples
      model$calculate(dependencies)
      return(exp(model$calculate(missingDataNames)))
    }
    leaveOutNames <- paramNames
    predLoss <- TRUE
  }
  modelMCMCConf <- configureMCMC(newModel, nodes = NULL, monitors = leaveOutNames)
  modelMCMCConf$samplerConfs <- conf$samplerConfs
  modelMCMC <- buildMCMC(modelMCMCConf)
  C.modelMCMC <- compileNimble(modelMCMC,
                               project = newModel)
  C.modelMCMC$run(MCMCiter)
  MCMCout <- as.matrix(C.modelMCMC$mvSamples)[,leaveOutNames]
  sampNum <- dim(MCMCout)[1]
  startIndex <- max(c(ceiling(burnInProp*MCMCiter), 1))
  message(paste("dropping data fold", i))
  if(predLoss){
    values(model, missingDataNames) <- saveData
  }
  crossValAveSD <- calcCrossValSD(MCMCout=MCMCout, sampNum=sampNum,
                                  NbootReps=NbootReps,
                                  saveData=saveData,
                                  lossFunction = lossFunction,
                                  predLoss = predLoss)
  crossVal <- mean(apply(MCMCout[startIndex:sampNum,], 1, lossFunction, saveData))
  if(predLoss){
    crossVal <- log(crossVal)
  }
  nimble:::clearCompiled(C.modelMCMC)
  return(list(crossValAverage = crossVal,
              crossValAveSD = crossValAveSD,
              samples= if(returnSamples) MCMCout else NA))
}

MSElossFunction <- function(simulatedDataValues, actualDataValues){
  MSE <- mean((simulatedDataValues - actualDataValues)^2)
  return(MSE)
}

generateRandomFoldFunction <- function(model, k){
  dataNodes <- model$expandNodeNames(model$getNodeNames(dataOnly = TRUE))
  nDataNodes <- length(dataNodes)
  if(k > nDataNodes){
    stop('Cannot have more folds than there are data nodes in the model.')
  }
  approxNumPerFold <- ceiling(nDataNodes/k)
  leaveOutNames <- list()
  reducedDataNodeNames <- dataNodes
  for(i in 1:(k-1)){
    leaveOutNames[[i]] <- sample(reducedDataNodeNames, approxNumPerFold)
    reducedDataNodeNames <- reducedDataNodeNames[-which(reducedDataNodeNames %in% leaveOutNames[[i]])]
  }
  leaveOutNames[[k]] <- reducedDataNodeNames
  randomFoldFunction <- function(i){
    return(leaveOutNames[[i]])
  }
}



#' Perform k-fold cross-validation on a NIMBLE model  
#' 
#' EXPERIMENTAL.  Takes a NIMBLE model and conducts k-fold cross validation on 
#' it, returning a measure of the model's predictive performance. 
#' 
#'  @param MCMCconfiguration a NIMBLE MCMC configuration object, returned by a 
#'  call to configureMCMC(). 
#'  @param k an integer argument determining the number of folds that should be
#'  used for cross-validation.
#'  @param foldFunction one of (1) an R function taking a single integer argument
#'  \code{i}, and returning a character vector with the names of the data nodes 
#'  to leave out of the model for fold \code{i}, or (2) the character string 
#'  \code{"random"}, indicating that data nodes will be randomly partitioned 
#'  into \code{k} folds.  Note that choosing "random" and setting \code{k} equal to the total number of data nodes in the model will perform leave-one-out cross-validation.  Defaults to \code{"random"}.  See 'Details'.
#'  @param lossFunction one of (1) an R function taking a set of simulated data
#'  and a set of observed data, and calculating the loss from those, or (2) a
#'  character string naming one of NIMBLE's built-in loss functions.  If a 
#'  character string, must be one of \code{"preictive"} to use the log 
#'  predictive density as a loss function or \code{"MSE"} to use the mean squared
#'  error as a loss function.  Defaults to \code{"MSE"}.  See 'Details'
#'  for information on creating a user-defined loss function.
#'  @param MCMCcontrol (optional) an R \code{list} with parameters governing the 
#'  MCMC algorithm,  See 'Details' for specific parameters.
#'  @param returnSamples a logical argument indicating whether to return all 
#'  posterior samples from all MCMC runs (note there will be \code{k} sets of
#'  posterior samples returned).  This can result in a very large returned
#'  object.  Defaults to \code{FALSE}.
#'  @param nCores an integer argument indicating the number of cpu cores to use
#'  for CV calculation.  In order to use >1 core, the \code{parallel} R package
#'  must be installed.  Only Mac and Linux operating systems support multiple
#'  cores at this time.  Defaults to 1.  
#'  @param nBootReps an integer argument indicating the number of bootstrap samples
#'  to use when estimating the Monte Carlo error of the cross validation metric.
#' 
#' @author Nicholas Michaud and Lauren Ponisio
#' @export
#' @details k-fold CV in NIMBLE proceeds by separating the data in a \code{nimbleModel} into k folds, as determined by the
#' \code{foldFunction} argument.  For each fold, the corresponding data is held out of the model, 
#' and MCMC is run to simulate data values for the left out data (these values are drawn from
#' the data's posterior predictive distribution).  Then, the simlated data values are compared to the 
#' known, left-out data values via the specified \code{lossFunction}.  The loss values are averaged over the
#' posterior samples for reach fold, and these averaged values for each fold are then averaged over all folds to produce a single 
#' loss estimate.  Additionally, estimates of the Monte Carlo error for each fold are returned.
#' 
#'  @section The \code{foldFunction} Argument:
#'  If the default 'random' method is not used, the \code{foldFunction} argument is an R function that takes a single integer argument \code{i}.  \code{i} is guaranteed to be within the range \code{[1, k]}.
#'  For each integer value \code{i}, the function should return a character vector of node names corresponding to the data nodes that will be left out of the model for that fold.
#'  The returned node names can be expanded, but don't need to be.   For example, if fold \code{i} is inteded to leave out the model nodes \code{x[1]}, \code{x[2]} and \code{x[3]} then the function could return either \code{c('x[1]', 'x[2]', 'x[3]')} or \code{c( 'x[1:3]')}.
#'   
#'  @section The \code{lossFunction} Argument:
#'  If you don't wish to use \ NIMBLE's built in "MSE" loss functions, you may provide your own R function
#'  as the \code{lossFunction} argument to \code{runCrossValidate}.  A user-supplied  \code{lossFunction} must be an R function
#'  that takes two arguments: the first, named \code{simulatedDataValues}, will be a vector of simulated data values.  The second,
#'  named \code{actualDataValues}, will be a vector of observed data values corresponding to the simulated data values in \code{simulatedDataValues}.  The loss function should return a single scalar number.
#'  See 'Examples' for an example of a user-defined loss function.
#'  
#' @section The \code{MCMCcontrol} Argument:
#' The \code{MCMCcontrol} argument is a list with the following elements:
#'	\itemize{
#'	  \item{\code{nMCMCiters}}	{
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
#' \item{\code{CVvalue}} {The CV value, measuring the model's ability to predict new data.  Smaller relative values indicate better model performance.}
#' \item{\code{CVstandardError}} {The standard error of the CV value, giving an indication of the total Monte Carlo error in the CV estimate.}
#' \item{\code{foldCVInfo}} {An list of fold CV values and standard errors for each fold.}
#' \item{\code{samples}} {An R list, only returned when \code{returnSamples = TRUE}.  The i'th element of this list will be a matrix of posterior samples from the model with the i'th fold of data left out.  There will be \code{k} sets of samples.}
#' }
#'
#' @examples
#' \dontrun{
#'
#' ## We conduct CV on the classic "dyes" BUGS model.
#' 
#' dyesCode <- nimbleCode({
#'   for (i in 1:BATCHES) {
#'     for (j in 1:SAMPLES) {
#'       y[i,j] ~ dnorm(mu[i], tau.within);
#'     }
#'     mu[i] ~ dnorm(theta, tau.between);
#'   }
#' 
#'   theta ~ dnorm(0.0, 1.0E-10);
#'   tau.within ~ dgamma(0.001, 0.001);  sigma2.within <- 1/tau.within;
#'   tau.between ~ dgamma(0.001, 0.001);  sigma2.between <- 1/tau.between;
#' })
#' 
#' dyesData <- list(y = matrix(c(1545, 1540, 1595, 1445, 1595,
#'                               1520, 1440, 1555, 1550, 1440,
#'                               1630, 1455, 1440, 1490, 1605,
#'                               1595, 1515, 1450, 1520, 1560, 
#'                               1510, 1465, 1635, 1480, 1580,
#'                               1495, 1560, 1545, 1625, 1445), 
#'                               nrow = 6, ncol = 5))
#'                               
#' dyesConsts <- list(BATCHES = 6,
#'                    SAMPLES = 5)
#'                    
#' dyesInits <- list(theta = 1500, tau.within = 1, tau.between =  1)
#' 
#' dyesModel <- nimbleModel(code = dyesCode,
#'                          constants = dyesConsts,
#'                          data = dyesData,
#'                          inits = dyesInits)
#' 
#' # Define the fold function.
#' # This function defines the data to leave out for the i'th fold
#' # as the i'th row of the data matrix y.  This implies we will have
#' # 6 folds.
#' 
#' dyesFoldFunction <- function(i){
#'   foldNodes_i <- paste0('y[', i, ', ]')  # will return 'y[1,]' for i = 1 e.g.
#'   return(foldNodes_i)
#' }
#' 
#' # We define our own loss function as well.
#' # The function below will compute the root mean squared error.
#' 
#' RMSElossFunction <- function(simulatedDataValues, actualDataValues){
#'   dataLength <- length(simulatedDataValues) # simulatedDataValues is a vector
#'   SSE <- 0
#'   for(i in 1:dataLength){
#'     SSE <- SSE + (simulatedDataValues[i] - actualDataValues[i])^2
#'   }
#'   MSE <- SSE / dataLength
#'   RMSE <- sqrt(MSE)
#'   return(RMSE)
#' }
#' 
#' dyesMCMCconfiguration <- configureMCMC(dyesModel)
#'   
#' crossValOutput <- runCrossValidate(MCMCconfiguration = dyesMCMCconfiguration,
#'                                    k = 6,
#'                                    foldFunction = dyesFoldFunction,
#'                                    lossFunction = RMSElossFunction,
#'                                    MCMCcontrol = list(nMCMCiters = 5000,
#'                                                       burnInProp = .1))
#'   
#' 

runCrossValidate <- function(MCMCconfiguration,
                             k,
                             foldFunction = 'random',
                             lossFunction = 'MSE',
                             MCMCcontrol = list(), 
                             returnSamples = FALSE,
                             nCores = 1,
                             NbootReps = 200){
  
  model <- MCMCconfiguration$model
  nMCMCiters <- MCMCcontrol[['nMCMCiters']]  
  if(is.null(nMCMCiters)){
    nMCMCiters <- 10000
    warning("Defaulting to 10,000 MCMC iterations for each MCMC run.")
  }
  burnInProp <-  MCMCcontrol[['burnInProp']]  
  if(is.null(burnInProp)){
    burnInProp <- 0
  }
  else if(burnInProp >= 1 | burnInProp < 0){
    stop("burnInProp needs to be between 0 and 1")
  }
  if(is.character(foldFunction) && foldFunction == 'random'){
    foldFunction <- generateRandomFoldFunction(model, k)
  }
  if(is.character(lossFunction) && lossFunction == 'MSE'){
    lossFunction <- MSElossFunction
  }
  allLeaveoutDataNames <- lapply(1:k, function(x){
    thisDataNames <- try(model$expandNodeNames(foldFunction(x)), silent = TRUE)
    if(inherits(thisDataNames, "try-error")){
      stop(paste("foldFunction is returning incorrect node names for fold", x))
    }
    else{
      if(!all(model$expandNodeNames(thisDataNames) %in% model$expandNodeNames(model$getNodeNames(dataOnly = TRUE)))){
        stop(paste("foldFunction is returning names of non-data nodes for fold", x))
      }
    }
  })
  
  libError <-  try(library('parallel'), silent = TRUE)
  if(inherits(libError, 'try-error') && nCores > 1){
    warning("The 'parallel' package must be installed to use multiple cores for CPPP calculation.  Defaulting to one core.")
    nCores <- 1
  }
  
  if(nCores > 1){
    ## simulate the ppp values
    crossValOut <- mclapply(1:k, calcCrossVal,
                            MCMCconfiguration,
                            foldFunction,
                            lossFunction,
                            nMCMCiters,
                            burnInProp,
                            returnSamples,
                            NbootReps,
                            mc.cores = ncores)
  }
  else{
    crossValOut <- lapply(1:k, calcCrossVal,
                          MCMCconfiguration,
                            foldFunction,
                            lossFunction,
                            nMCMCiters,
                            burnInProp,
                            returnSamples,
                            NbootReps)
  }
  CVvalue <- sum(sapply(crossValOut, function(x) x$crossValAverage),
                  na.rm=TRUE)
  CVstandardError <- sqrt(sum(sapply(crossValOut,
                                function(x) x$crossValAveSD^2),
                         na.rm=TRUE))
  foldCVinfo <- lapply(crossValOut, function(x){return(c(foldCVvalue = x$crossValAverage,
                                                    foldCVstandardError = x$crossValAveSD))})
  out <- list(CVvalue=CVvalue,
              CVstandardError=CVstandardError,
              foldCVinfo = foldCVinfo)
  if(!returnSamples){
    out$samples <- NULL
  }
  else{
    out$samples <- lapply(crossValOut, function(x) x$samples)
  }
  return(out)
}
