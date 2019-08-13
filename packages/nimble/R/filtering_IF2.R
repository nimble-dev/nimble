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

IF2SetParVirtual <- nimbleFunctionVirtual(
  methods = list(
    scalarSet = function(scalars = double(1), mvWset = integer(), m = integer(),
                         ids = integer(1)){},
    vectorSet = function(vectors = double(2), mvWset = integer(), m = integer(),
                         ids = integer(1)){},
    scalarGet = function(mvWset = integer(), m = integer()){
      returnType(double(1))},
    vectorGet = function(mvWset = integer(), m = integer(), length = integer()){
      returnType(double(2))},
    calcScalarMean = function(m = integer()){returnType(double())},
    calcVectorMean = function(m = integer(), length = integer()){
      returnType(double(1))}
  )
)

doPars1 <- nimbleFunction(
  contains = IF2SetParVirtual,
  setup = function(parName, mvWSamples, mvEWSamples) {
  },
  methods = list(
    scalarSet = function(scalars = double(1), mvWset = integer(), m = integer(), ids = integer(1)){
      if(mvWset == 1){
        for(i in 1:m){
          mvWSamples[parName, i][1] <<- scalars[ids[i]]
        }
      }
      else{
        for(i in 1:m){
          mvEWSamples[parName, i][1] <<- scalars[ids[i]]
        }
      }
    },
    vectorSet = function(vectors = double(2), mvWset = integer(), m = integer(), ids = integer(1)){
      if(mvWset == 1){
        for(i in 1:m){
          mvWSamples[parName, i] <<- vectors[,ids[i]]
        }
      }
      else{
        for(i in 1:m){
          mvEWSamples[parName, i] <<- vectors[,ids[i]]
        }
      }
    },
    scalarGet = function(mvWset = integer(), m = integer()){
        vecOut <- numeric(m, init=FALSE)
        if(mvWset == 1){
            for(i in 1:m){
          vecOut[i] <- mvWSamples[parName, i][1]
        }
      }
      else{
        for(i in 1:m){
          vecOut[i] <- mvEWSamples[parName, i][1]
        }
      }
      returnType(double(1))
      return(vecOut)
    },
    vectorGet = function(mvWset = integer(), m = integer(), length = integer()){
        matOut <- matrix(nrow = length, ncol = m, init=FALSE)
        if(mvWset == 1){
        for(i in 1:m){
          matOut[,i] <- mvWSamples[parName, i]
        }
      }
      else{
        for(i in 1:m){
          matOut[,i] <- mvEWSamples[parName, i]
        }
      }
      returnType(double(2))
      return(matOut)
    },
    calcScalarMean = function(m = integer()){
      returnType(double())
      sumMean = 0
      for(i in 1:m){
        sumMean <- sumMean + mvWSamples[parName, i][1]*1/m

      }
      return(sumMean)
    },
    calcVectorMean = function(m = integer(), length = integer()){
      returnType(double(1))
      sumMean <- numeric(length)
      for(i in 1:m){
        sumMean[1:length] <- sumMean[1:length] + 
          mvWSamples[parName, i][1:length]*1/m 
      }
      return(sumMean)
    }
  ), where = getLoadingNamespace())


IF2Step <- nimbleFunction(
  contains = IF2StepVirtual,
  setup = function(model, mvWSamples, mvEWSamples, nodes, paramVarDims,
                   iNode, paramNodes, paramVars, names,
                   sigma, silent = FALSE) {
    notFirst <- iNode != 1
    isSecond <- iNode == 2
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    parDeterm <- model$getDependencies(paramNodes, determOnly=TRUE)
    parAndPrevDeterm <- c(parDeterm, prevDeterm)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    t <- iNode  # current time point
    totalTime <- length(nodes)     
    isLast <- (t == totalTime)
    
    # Get names of xs node for current and previous time point (used in copy)

    prevXName <- names    
    thisXName <- names
    currInd <- 1
    prevInd <- 1 

    numParams <- sum(paramVarDims)
    numParamVars <- length(paramVars)
    paramInds <- c(0,cumsum(paramVarDims)) 
    singleParam <- (numParams == 1)
    
    ##  varSize keeps track of the size of each parameter we are estimating
    varSize <- rep(0, length(paramInds)-1)
    doVarList <- nimbleFunctionList(IF2SetParVirtual)
    for(i in 1:numParamVars){
      doVarList[[i]] <- doPars1(paramVars[i], mvWSamples, mvEWSamples)
      varSize[i] <- paramInds[i+1]-paramInds[i]
    }
    if(singleParam)
      varSize = c(varSize, 0)  # ensure that varSize is treated as a vector even with only one parameter
  },
  run = function(m = integer(), j = integer(), coolingRate = double(), useStoredSamples = integer()) {
    returnType(double())
    l <- numeric(m, init=FALSE)
    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    tmpPars <- matrix(nrow = numParams, ncol = m, init=FALSE)
    coolParam <- (coolingRate)^((t-1+(j - 1)*totalTime)/(50*totalTime))
    coolSigma <- coolParam*sigma
    
    for(i in 1:m) {
      ids[i] <- i ## for initial weights copying
    }
    if(useStoredSamples == 1){
     for(j in 1:numParamVars){
       if(varSize[j] == 1){
         tmpPars[paramInds[j+1], ] <- doVarList[[j]]$scalarGet(0,m)
       }
       else{
         tmpPars[(paramInds[j]+1):paramInds[j+1], ] <- doVarList[[j]]$vectorGet(0, m, varSize[j])
       }
     }
    }
    for(i in 1:m) {
       if(useStoredSamples == 1) {
         if(singleParam == 1){
           tmpPars[1,i] <- rnorm(1, tmpPars[1,i],  coolSigma[1])
           values(model, paramNodes) <<- tmpPars[,i]
         }
         else{
           for(j in 1:numParamVars){
              tmpPars[j, i] <- rnorm(1, tmpPars[j, i],
                                                coolSigma[j])
          
          }
          
           
           values(model, paramNodes) <<- tmpPars[,i] 
         }
         if(notFirst){
           copy(mvEWSamples, model, nodes = prevXName, nodesTo = prevNode, row = i)
         }
       }
       else{
         for(j in 1:numParamVars){
           if(varSize[j] == 1){
             tmpPars[paramInds[j+1], i] <- rnorm(1, values(model, paramVars[j])[1],
                                                 coolSigma[paramInds[j+1]])
           }
           else{
             tmpPars[paramInds[j]+j, i] <- rnorm(1, values(model, paramVars[j])[j],
                                                coolSigma[j])
             
           }
         }
         values(model, paramNodes) <<- tmpPars[,i]
      }
      calculate(model, parAndPrevDeterm)
      simulate(model, thisNode)
      copy(model, mvWSamples, nodes = thisNode, nodesTo = thisXName, rowTo = i)
      calculate(model, thisDeterm)
      wts[i]  <- exp(calculate(model, thisData))
      if(is.nan(wts[i])) wts[i] <- 0
      if(calculate(model, paramNodes) == -Inf) 
        wts[i] <- 0
    }
    for(j in 1:numParamVars){
      if(varSize[j] == 1){
        doVarList[[j]]$scalarSet(tmpPars[paramInds[j+1], ], 1, m, ids)
      }
      else{
        doVarList[[j]]$vectorSet(tmpPars[(paramInds[j]+1):paramInds[j+1], ], 1, m, ids)
      }
    }
    wts <- wts/sum(wts)
    rankSample(wts, m, ids, silent)
    for(i in 1:m){
      copy(mvWSamples, mvEWSamples, nodes = thisXName, nodesTo = thisXName, row = ids[i], rowTo = i)
      mvWSamples['wts',i][currInd] <<- log(wts[i])
    }
    for(j in 1:numParamVars){
      if(varSize[j] == 1){
         doVarList[[j]]$scalarSet(tmpPars[paramInds[j+1], ], 0, m, ids)
       }
       else{
         doVarList[[j]]$vectorSet(tmpPars[(paramInds[j]+1):paramInds[j+1], ], 0, m, ids)
       }
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

    paramsLength <- length(model$expandNodeNames(params, 
                                                 returnScalarComponents = TRUE))
    if(is.null(sigma)){
      sigma <- rep(1, paramsLength)
    }
    if(length(sigma) < paramsLength)
       stop("buildIteratedFilter2: The 'sigma' control list argument must be a vector specifying a
            non-negative perturbation magnitude for each element of 'params'. The length of 'sigma'
            does not match the length of 'params'.")
    
    if(any(sigma < 0))
      stop("buildIteratedFilter2: All values of 'sigma' should be non-negative.")
    
    ## Check for inits values
    if(is.null(inits)){
      inits <- values(model, params)
    }
    if(length(inits) < paramsLength)
      stop("buildIteratedFilter2: The 'inits' control list argument must be a vector specifying initial values for each element of 'params'. The length of 'inits' does not match the number of parameters.")
    
    #get latent state info
    varName <- sapply(nodes, function(x) { return(model$getVarNames(nodes = x)) } )
    if(length(unique(varName)) > 1){
      stop("buildIteratedFilter2: All latent nodes must be in the same variable.")
    }
    varName <- varName[1]
    info <- model$getVarInfo(varName)
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
    latentVars <- model$getVarNames(nodes = nodes)
      
    # if unspecified, parameter nodes are specified as all stochastic top level nodes which
    # are not in the set of latent nodes above
    if(all(is.null(params))){
      params <-  model$getNodeNames(stochOnly = TRUE, includeData = FALSE,
                                           topOnly = TRUE)
      parLatents <- sapply(params, function(x){ return(model$getVarNames(nodes = x) %in% latentVars)})
      params <- params[!parLatents]
    }
    unsortParams <- model$expandNodeNames(params, sort = FALSE)
    params <- model$expandNodeNames(params, sort = TRUE)
    sortingOrder <- sapply(unsortParams, function(x) return(which(x == params)))
    if(identical(params, character(0)))
      stop('buildIteratedFilter2: There must be at least one higher level parameter for IF2 to work.')
    if(any(params %in% nodes))
      stop('buildIteratedFilter2: Parameters cannot be latent states.')
    if(!all(params%in%model$getNodeNames(stochOnly = TRUE)))
      stop('buildIteratedFilter2: Parameters must be stochastic nodes.')
    paramVars <-  model$getVarNames(nodes =  params)  # need var names too
    pardimcheck <- sapply(paramVars, function(var){
      if(length(nimDim(model[[var]]))>1)
        stop("buildIteratedFilter2: IF2 doesn't work for matrix-valued top-level parameters.")
    })
    dims <- lapply(nodes, function(var) nimDim(model[[var]]))
    if(length(unique(dims)) > 1) stop("buildIteratedFilter2: sizes or dimension of latent states cannot vary.")
    paramDims <-   sapply(params, function(n) nimDim(model[[n]]))
    
    my_initializeModel <- initializeModel(model, silent = silent)
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[c(latentVars, paramVars)]
    names <- sapply(modelSymbolObjects, function(x) return(x$name))
    type <- sapply(modelSymbolObjects, function(x) return(x$type))
    size <- lapply(modelSymbolObjects, function(x) {
      if(identical(x$size, numeric(0))) return(1)
      return(x$size)})
    
    size[[latentVars]] <- as.numeric(dims[[1]])      
    mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))
    names <- c(names, "wts")
    type <- c(type, "double")
    size$wts <- 1
    mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))
    
    names <- names[1]
    varSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[paramVars]
    paramVarDims <- unlist(lapply(varSymbolObjects, function(x){
      if(identical(x$size, numeric(0))) return(1)
      return(x$size)}))

    IF2StepFunctions <- nimbleFunctionList(IF2StepVirtual)
    for(iNode in seq_along(nodes))
      IF2StepFunctions[[iNode]] <- IF2Step(model,  mvWSamples, mvEWSamples, nodes, paramVarDims, 
                                         iNode, params, paramVars, names, sigma, silent)
    
    numParams <- sum(paramVarDims)
    numParamVars <- length(paramVars)
    paramInds <- c(0, cumsum(paramVarDims)) 
    singleParam <- (numParams == 1)
    
    ##  varSize keeps track of the size of each parameter we are estimating
    varSize <- rep(0, length(paramInds)-1)
    doVarList <- nimbleFunctionList(IF2SetParVirtual)
    for(i in 1:numParamVars){
      doVarList[[i]] <- doPars1(paramVars[i], mvWSamples, mvEWSamples)
      varSize[i] <- paramInds[i+1]-paramInds[i]
    }
    if(singleParam)
      varSize = c(varSize, 0)  # ensure that varSize is treated as a vector even with only one parameter
    
    oldJ <- 0
    oldM <- 0
  },
  run = function(m = integer(default = 10000), n = integer(default = 5), 
                 coolingRate = double(default = 0.2)) {
    
    values(model, params) <<- inits
    my_initializeModel$run()
    resize(mvWSamples, m)
    resize(mvEWSamples, m)
    outPars <- numeric(numParams)
    sortedOutPars <- numeric(numParams)
    for(i in 1:m)
      mvWSamples['wts',i][1] <<- log(1/m)
   
    for(j in 1:n){
      for(iNode in seq_along(IF2StepFunctions)) {
       useStoredSamples <- ((iNode > 1) | j > 1)
       IF2StepFunctions[[iNode]]$run(m, j, coolingRate, useStoredSamples)
      }
    }
    for(k in 1:numParamVars){
      if(varSize[k] == 1){
        outPars[paramInds[k+1]] <- doVarList[[k]]$calcScalarMean(m)
      }
      else{
        outPars[(paramInds[k]+1):paramInds[k+1]] <- 
          doVarList[[k]]$calcVectorMean(m, varSize[k])
      }
    }
    oldM <<- m
    oldJ <<- j 
    
    for(i in 1:paramsLength)
      sortedOutPars[i] <- outPars[i] 
    
    returnType(double(1))
    return(sortedOutPars)
  },
  methods = list(
    continueRun = function(n = integer(default = 5), coolingRate = double(default = 0.2)){
      outPars <- numeric(numParams)
      sortedOutPars <- numeric(numParams)
      useStoredSamples <- 1
      newN <- oldJ + n
      for(j in (oldJ+1):newN){
        for(iNode in seq_along(IF2StepFunctions)) {
          IF2StepFunctions[[iNode]]$run(oldM, j, coolingRate, 
                                        useStoredSamples)
        }
      }
      for(k in 1:numParamVars){
        if(varSize[k] == 1){
          outPars[paramInds[k+1]] <- doVarList[[k]]$calcScalarMean(oldM)
        }
        else{
          outPars[(paramInds[k]+1):paramInds[k+1]] <- 
            doVarList[[k]]$calcVectorMean(oldM, varSize[k])
        }
      }
      oldJ <<- newN
      for(i in 1:paramsLength)
        sortedOutPars[i] <- outPars[i]
      returnType(double(1))
      return(sortedOutPars)
    }
  ),
  where = getLoadingNamespace()
)
