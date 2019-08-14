##  Contains code to run a Bayes map iterated filter.  The filter has a build
##  function (buildIteratedFilter2) and a step function (IF2stepNS). Also contains 
##  a function for calculating useful quantities. The doPars2() function
##  is specialized to different variables in the mv objects, and allows
##  these variables to be set / retrieved within the mv objects without
##  explicitly calling the variable name.
##  
##
##  Implementing follows Ionides et. all 2015 and adapted Liu and West filter of Nicholas Michaud

IF2StepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer(), n = integer(), coolingRate = double(),
                 useStoredSamples = integer())
    returnType(double())
)

IF2Step <- nimbleFunction(
  contains = IF2StepVirtual,
  setup = function(model, mvWSamples, mvEWSamples, latentNodes, latentVar,
                   iNode, paramNodes, numParams, sigma, silent = FALSE) {
    notFirst <- iNode != 1
    isSecond <- iNode == 2
    prevNode <- latentNodes[if(notFirst) iNode-1 else iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- latentNodes[iNode]
    parDeterm <- model$getDependencies(paramNodes, determOnly=TRUE)
    parAndPrevDeterm <- c(parDeterm, prevDeterm)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    t <- iNode  # current time point
    totalTime <- length(latentNodes)
    isLast <- (t == totalTime)
  },
  run = function(m = integer(), j = integer(), coolingRate = double(), useStoredSamples = integer()) {
    returnType(double())
    l <- numeric(m, init=FALSE)
    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    coolParam <- (coolingRate)^((t-1+(j - 1)*totalTime)/(50*totalTime))
    coolSigma <- coolParam*sigma
    
    for(i in 1:m) {
      ids[i] <- i ## for initial weights copying
    }
    for(i in 1:m) {
        if(useStoredSamples == 1) { 
            nimCopy(mvEWSamples, model, nodes = paramNodes, row = i)
            if(notFirst)
                copy(mvEWSamples, model, nodes = latentVar, nodesTo = prevNode, row = i)
        }
        currentValues <- values(model, paramNodes)
        for(j in 1:numParams)
            currentValues[j] <- rnorm(1, currentValues[j], coolSigma[j])
        values(model, paramNodes) <<- currentValues
        calculate(model, parAndPrevDeterm)
        simulate(model, thisNode)
        calculate(model, thisDeterm)
        wts[i]  <- exp(calculate(model, thisData))
        if(is.nan(wts[i])) wts[i] <- 0
        logProb <- calculate(model, paramNodes)
        if(is.na(logProb) | logProb == -Inf)
            wts[i] <- 0
        nimCopy(model, mvWSamples, nodes = thisNode, nodesTo = latentVar, rowTo = i)
        nimCopy(model, mvWSamples, nodes = paramNodes, rowTo = i)
    }
    wts <- wts/sum(wts)
    rankSample(wts, m, ids, silent)
    for(i in 1:m){
      copy(mvWSamples, mvEWSamples, nodes = latentVar, row = ids[i], rowTo = i)
      copy(mvWSamples, mvEWSamples, nodes = paramNodes, row = ids[i], rowTo = i)
      mvWSamples['wts',i][1] <<- log(wts[i])
    }
  return(0)
  },  where = getLoadingNamespace()
)

#' Create an IF2 algorithm.  
#' 
#' @description Create an IF2 algorithm for a given NIMBLE state space model.  
#'
#' @param model A NIMBLE model object, typically representing a state 
#'  space model or a hidden Markov model.
#' @param nodes A character vector specifying the latent model nodes 
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.
#' @param params A character vector specifying the top-level parameters to obtain maximum likelihood estimates of. 
#'   If unspecified, parameter nodes are specified as all stochastic top level nodes which
#'  are not in the set of latent nodes specified in \code{nodes}.
#' @param control  A list specifying different control options for the IF2 algorithm.  Options are described in the \sQuote{details} section below.

#' @author Nicholas Michaud, Dao Nguyen, and Christopher Paciorek
#' @family particle filtering methods
#' @details 
#' 
#' Each of the \code{control()} list options are described in detail below:
#' \describe{
#'  \item{sigma}{A vector specifying a non-negative perturbation magnitude for each element of the \code{params} argument.  Defaults to a vector of 1's.}
#'  \item{inits}{A vector specifying an initial value for each element of the \code{params} argument.  Defaults to the parameter values in the model at the time the model is built.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.  
#'  Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to \code{TRUE}.}
#' }
#' 
#'  The IF2 agorithm uses iterated filtering to estimate maximum likelihood values for top-level parameters for a state space model.  
#'        
#'  The resulting specialized IF2 algorithm will accept the following arguments:
#'  
#'  \describe{
#'    \item{m}{A single integer specifying the number of particles to use for each run of the filter.  Defaults to 10,000.}
#'    \item{n}{A single integer specifying the number of overall filter iterations to run.  Defaults to 5.}
#'    \item{coolingRate}{A double specifying the cooling rate to use for the IF2 algorithm.  Defaults to 0.2.}
#'
#'  The \code{run} fuction will return a vector with MLE estimates.  Additionally, once the specialized algorithm has been run, it can be continued for additional iterations by calling the \code{continueRun} method.
#'  
#' @export
#' 
#' @references Ionides, E.L., D. Nguyen, Y. Atchad{\'e}, S. Stoev, and A.A. King (2-15). Inference for dynamic and latent variable models via iterated, perturbed Bayes maps. \emph{Proceedings of the National Academy of Sciences}, 112(3), 719-724.
#' 
#' @examples
#' \dontrun{
#' model <- nimbleModel(code = ...)
#' my_IF2 <- buildIteratedFilter2(model, 'x[1:100]', params = 'sigma_x')
#' Cmodel <- compileNimble(model)
#' Cmy_IF2 <- compileNimble(my_IF2, project = model)
#' # MLE estimate of a top level parameter named sigma_x:
#' sigma_x_MLE <- Cmy_IF2$run(m = 10000, n = 10)
#' # Continue running algorithm for more precise estimate:
#' sigma_x_MLE <- Cmy_IF2$continueRun(n = 10) 
#' }
buildIteratedFilter2 <- nimbleFunction(
  setup = function(model, nodes, params = NULL, control = list()){
    
    #control list extraction
    silent <- control[['silent']]
    inits <- control[['inits']]
    sigma <- control[['sigma']]
    timeIndex <- control[['timeIndex']]
    initModel <- control[['initModel']]
    if(is.null(silent)) silent <- TRUE
    if(is.null(initModel)) initModel <- TRUE

    # if unspecified, parameter nodes are specified as all stochastic top level nodes which
    # are not in the set of latent nodes above
    if(all(is.null(params))){
      params <-  model$getNodeNames(stochOnly = TRUE, includeData = FALSE,
                                    topOnly = TRUE)
      params <- params[!params %in% nodes]
    }
    params <- model$expandNodeNames(params)
    numParams <- length(params)
    parDeterm <- model$getDependencies(params, determOnly=TRUE)

    if(identical(params, character(0)))
      stop('buildIteratedFilter2: There must be at least one parameter for IF2 to optimize with respect to.')
    if(any(params %in% nodes))
      stop('buildIteratedFilter2: Parameters cannot be latent states.')
    if(!all(params %in% model$getNodeNames(stochOnly = TRUE)))
        stop('buildIteratedFilter2: Parameters must be stochastic nodes.')
      
    paramVars <-  model$getVarNames(nodes = params)

    if(is.null(sigma)){
      sigma <- rep(1, numParams)
    }
    if(length(sigma) != numParams)
       stop("buildIteratedFilter2: The 'sigma' control list argument must be a vector specifying a
            non-negative perturbation magnitude for each element of 'params'. The length of 'sigma'
            does not match the length of 'params'.")
    
    if(any(sigma < 0))
      stop("buildIteratedFilter2: All values of 'sigma' should be non-negative.")
    
    ## Check for inits values
    if(is.null(inits)){
      inits <- values(model, params)
    }
    if(length(inits) < numParams)
      stop("buildIteratedFilter2: The 'inits' control list argument must be a vector specifying initial values for each element of 'params'. The length of 'inits' does not match the number of parameters.")
    
    #get latent state info
    latentVar <- model$getVarNames(nodes = nodes) 
    if(length(unique(latentVar)) > 1){
      stop("buildIteratedFilter2: All latent nodes must be in the same variable.")
    }

    info <- model$getVarInfo(latentVar)
    latentDims <- info$nDim
    if(is.null(timeIndex)){
      timeIndex <- which.max(info$maxs)
      timeLength <- max(info$maxs)
      if(sum(info$maxs == timeLength) > 1) # check if multiple dimensions share the max index size
         stop("buildIteratedFilter2: Unable to determine which dimension indexes time. Specify manually using the 'timeIndex' control list argument.")
    } else{
      timeLength <- info$maxs[timeIndex]
    }
     ## CJP note: this assumes nodes in variable are in time order. Should we tell the user this is our assumption?  
    nodes <- paste(info$varName, "[", rep(",", timeIndex-1), 1:timeLength,
                   rep(",", info$nDim - timeIndex), "]", sep="")
             
    dims <- lapply(nodes, function(node) nimDim(model[[node]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent 
                                      states varies')

    my_initializeModel <- initializeModel(model, silent = silent)

    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[c(latentVar, paramVars)]
    names <- sapply(modelSymbolObjects, function(x) return(x$name))
    type <- sapply(modelSymbolObjects, function(x) return(x$type))
    size <- lapply(modelSymbolObjects, function(x) {
      if(identical(x$size, numeric(0))) return(1)
      return(x$size)})

    ## check why size for latent is not already ok
    size[[latentVar]] <- as.numeric(dims[[1]])      
    mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))
    names <- c(names, "wts")
    type <- c(type, "double")
    size$wts <- 1
    mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))
    
    IF2StepFunctions <- nimbleFunctionList(IF2StepVirtual)
    for(iNode in seq_along(nodes))
      IF2StepFunctions[[iNode]] <- IF2Step(model,  mvWSamples, mvEWSamples, nodes, latentVar,
                                           iNode, params, numParams, sigma, silent)

    estimate <- nimNumeric(numParams)  
  
    oldJ <- 0
    oldM <- 0
  },
  run = function(m = integer(default = 10000), n = integer(default = 5), 
                 coolingRate = double(default = 0.2)) {
    
    values(model, params) <<- inits
    my_initializeModel$run()
    resize(mvWSamples, m)
    resize(mvEWSamples, m)

    for(i in 1:m)
      mvWSamples['wts',i][1] <<- log(1/m)
   
    for(j in 1:n){
      for(iNode in seq_along(IF2StepFunctions)) {
       useStoredSamples <- ((iNode > 1) | j > 1)
       IF2StepFunctions[[iNode]]$run(m, j, coolingRate, useStoredSamples)
      }
    }

      ## Calculate mean as final estimate.
    estimate <<- rep(0, numParams)
      for(i in 1:m) {
          nimCopy(mvWSamples, model, nodes = params, row = i)
          estimate <<- estimate + values(model, params)
      }
      estimate <<- estimate / m

    ## Leave model with parameter estimates.  
      values(model, params) <<- estimate
      model$calculate(parDeterm)
    oldM <<- m
    oldJ <<- j 
    
    returnType(double(1))
    return(estimate)
  },
  methods = list(
    continueRun = function(n = integer(default = 5), coolingRate = double(default = 0.2)){
      useStoredSamples <- 1
      newN <- oldJ + n
      for(j in (oldJ+1):newN){
        for(iNode in seq_along(IF2StepFunctions)) {
          IF2StepFunctions[[iNode]]$run(oldM, j, coolingRate, 
                                        useStoredSamples)
        }
      }

      ## Calculate mean as final estimate.
    estimate <<- rep(0, numParams)
      for(i in 1:oldM) {
          nimCopy(mvWSamples, model, nodes = params, row = i)
          estimate <<- estimate + values(model, params)
      }
      estimate <<- estimate / oldM

    ## Leave model with parameter estimates.  
      values(model, params) <<- estimate
      model$calculate(parDeterm)

      oldJ <<- newN

      returnType(double(1))
      return(estimate)
    }
  ),
  where = getLoadingNamespace()
)
