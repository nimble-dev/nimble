##  Contains code to run bootstrap particle filters.
##  We have a build function (buildBootstrapFilter),
##  and step function.

bootStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer(), threshNum=double(), prevSamp = logical()) {
    returnType(double(1))
  },
  methods = list(
    returnESS = function() {
      returnType(double())
    }
  )
)

# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.

bootFStep <- nimbleFunction(
  name = 'bootFStep',
  contains = bootStepVirtual,
  setup = function(model,
                   mvEWSamples,
                   mvWSamples,
                   nodes,
                   iNode,
                   names,
                   saveAll,
                   smoothing,
                   resamplingMethod,
                   silent = FALSE) {
    notFirst <- iNode != 1
    modelSteps <- particleFilter_splitModelSteps(model, nodes, iNode, notFirst)
    prevDeterm <- modelSteps$prevDeterm
    calc_thisNode_self <- modelSteps$calc_thisNode_self
    calc_thisNode_deps <- modelSteps$calc_thisNode_deps
    
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    thisNode <- nodes[iNode]
    ## prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    ## thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    ## thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    ## t is the current time point.
    t <- iNode
    ## Get names of xs node for current and previous time point (used in copy)
    if(saveAll == 1){
      allPrevNodes <- model$expandNodeNames(nodes[1:(iNode-1)])
      prevXName <- prevNode    
      thisXName <- thisNode
      currInd <- t
      prevInd <- t-1
      if(isTRUE(smoothing)){
        currInd <- 1
        prevInd <- 1
      }
    }
    else{
      allPrevNodes <- names
      prevXName <- names    
      thisXName <- names
      currInd <- 1
      prevInd <- 1 
    }
    isLast <- (iNode == length(nodes))
    ess <- 0
    resamplerFunctionList <- nimbleFunctionList(resamplerVirtual)
    if(resamplingMethod == 'default'){
      resamplerFunctionList[[1]] <- residualResampleFunction()
      defaultResamplerFlag <- TRUE
    }
    if(resamplingMethod == 'residual')
      resamplerFunctionList[[1]] <- residualResampleFunction()
    if(resamplingMethod == 'multinomial')
      resamplerFunctionList[[1]] <- multinomialResampleFunction()
    if(resamplingMethod == 'stratified')
      resamplerFunctionList[[1]] <- stratifiedResampleFunction()
    if(resamplingMethod == 'systematic')
      resamplerFunctionList[[1]] <- systematicResampleFunction()
  },
  run = function(m = integer(),
                 threshNum = double(),
                 prevSamp = logical()) {
    returnType(double(1))
    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    llEst <- numeric(m, init=FALSE)
    out <- numeric(2, init=FALSE)
    
    for(i in 1:m) {
      if(notFirst) {  
        if(smoothing == 1){
          copy(mvEWSamples, mvWSamples, nodes = allPrevNodes, 
               nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        copy(mvEWSamples, model, nodes = prevXName, nodesTo = prevNode, row = i)
        model$calculate(prevDeterm) 
      }
      model$simulate(calc_thisNode_self)
      ## The logProbs of calc_thisNode_self are, correctly, not calculated.
      copy(model, mvWSamples, nodes = thisNode, nodesTo = thisXName, row = i)
      wts[i]  <- model$calculate(calc_thisNode_deps)
      if(is.nan(wts[i])){
        out[1] <- -Inf
        out[2] <- 0
        return(out)
      }
      if(prevSamp == 0){ ## defaults to 1 for first step and then comes from the previous step
        wts[i] <- wts[i] + mvWSamples['wts',i][prevInd]
        llEst[i] <- wts[i]
      }
      else{
          llEst[i] <- wts[i] - log(m)
      }
    }
    
    stepllEst <- log(sum(exp(llEst)))
    if(is.nan(stepllEst)){
      out[1] <- -Inf
      out[2] <- 0
      return(out)
    }
    if(stepllEst == Inf | stepllEst == -Inf){
      out[1] <- -Inf
      out[2] <- 0
      return(out)
    }
    
    out[1] <- stepllEst
    
    ## Normalize weights and calculate effective sample size .
    wts <- exp(wts)/sum(exp(wts))
    ess <<- 1/sum(wts^2) 
    
    ## Determine whether to resample by weights or not.
    if(ess < threshNum){
      if(defaultResamplerFlag){
        rankSample(wts, m, ids, silent)	 
      }
      else{
        ids <- resamplerFunctionList[[1]]$run(wts)  
      }
      ## out[2] is an indicator of whether resampling takes place.
      ## Resampling affects how ll estimate is calculated at next time point.
      out[2] <- 1
      for(i in 1:m){
        copy(mvWSamples, mvEWSamples, nodes = thisXName, nodesTo = thisXName,
             row = ids[i], rowTo = i)
        if(smoothing == 1){
          copy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = ids[i], rowTo=i)
        }
        mvWSamples['wts',i][currInd] <<- log(wts[i])
      }
    }
    else{
      out[2] <- 0
      for(i in 1:m){
        copy(mvWSamples, mvEWSamples, nodes = thisXName, nodesTo = thisXName,
             row = i, rowTo = i)
        if(smoothing == 1){
          copy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        mvWSamples['wts',i][currInd] <<- log(wts[i])
      }
    }
    return(out)
  },
  methods = list(
    returnESS = function(){
      returnType(double(0))
      return(ess)
    }
  )
)

#' Create a bootstrap particle filter algorithm to estimate log-likelihood.
#'
#'@description Create a bootstrap particle filter algorithm for a given NIMBLE state space model.  
#'
#' @param model A nimble model object, typically representing a state 
#'  space model or a hidden Markov model.
#' @param nodes  A character vector specifying the latent model nodes 
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.  All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author Daniel Turek and Nicholas Michaud
#' @details 
#' 
#' Each of the \code{control()} list options are described in detail here:
#' \describe{
#'  \item{thresh}{ A number between 0 and 1 specifying when to resample: the resampling step will occur when the
#'   effective sample size is less than \code{thresh} times the number of particles. Defaults to 0.8. Note that at the last time step, resampling will always occur so that the \code{mvEWsamples} \code{modelValues} contains equally-weighted samples.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter. Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function),  \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}.  Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (TRUE), or only for the most recent time point (FALSE)}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}.  \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.  
#'  Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#' }
#' 
#'  The bootstrap filter starts by generating a sample of estimates from the 
#'  prior distribution of the latent states of a state space model.  At each time point, these particles are propagated forward 
#'  by the model's transition equation.  Each particle is then given a weight 
#'  proportional to the value of the observation equation given that particle. 
#'  The weights are used to draw an equally-weighted sample of the particles at this time point.
#'  The algorithm then proceeds
#'  to the next time point.  Neither the transition nor the observation equations are required to 
#'  be normal for the bootstrap filter to work.   
#'  
#'  The resulting specialized particle filter algorthm will accept a
#'  single integer argument (\code{m}, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm 
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the \code{mvWSamples} modelValues object, with corresponding logged weights in \code{mvWSamples['wts',]}.
#'  An equally weighted sample from the posterior can be found in the \code{mvEWSamples} \code{modelValues} object.
#'  
#'  Note that if the \code{thresh} argument is set to a value less than 1, resampling may not take place at every time point.  
#'  At time points for which resampling did not take place, \code{mvEWSamples} will not contain equally weighted samples.
#'  To ensure equally weighted samples in the case that \code{thresh < 1}, we recommend resampling from \code{mvWSamples} at each time point 
#'  after the filter has been run, rather than using \code{mvEWSamples}.
#'
#' @section \code{returnESS()} Method:
#'  Calling the \code{returnESS()} method of a bootstrap filter after that filter has been \code{run()} for a given model will return a vector of ESS (effective
#'  sample size) values, one value for each time point.
#'  
#' @export
#' 
#' @family particle filtering methods
#' @references Gordon, N.J., D.J. Salmond, and A.F.M. Smith. (1993). Novel approach to nonlinear/non-Gaussian Bayesian state estimation. \emph{IEEE Proceedings F (Radar and Signal Processing)}. Vol. 140. No. 2. IET Digital Library, 1993.
#' @examples
#' \dontrun{
#' model <- nimbleModel(code = ...)
#' my_BootF <- buildBootstrapFilter(model, 'x[1:100]', control = list(thresh  = 1))
#' Cmodel <- compileNimble(model)
#' Cmy_BootF <- compileNimble(my_BootF, project = model)
#' logLike <- Cmy_BootF$run(m = 100000)
#' ESS <- Cmy_BootF$returnESS()
#' boot_X <- as.matrix(Cmy_BootF$mvEWSamples)
#' }
buildBootstrapFilter <- nimbleFunction(
    name = 'buildBootstrapFilter',
  setup = function(model, nodes, control = list()) {
    
    #control list extraction
    thresh <- control[['thresh']]
    saveAll <- control[['saveAll']]
    smoothing <- control[['smoothing']]
    silent <- control[['silent']]
    timeIndex <- control[['timeIndex']]
    initModel <- control[['initModel']]
    resamplingMethod <- control[['resamplingMethod']]
    if(is.null(thresh)) thresh <- .8
    if(is.null(silent)) silent <- TRUE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(initModel)) initModel <- TRUE
    if(is.null(resamplingMethod)) resamplingMethod <- 'default'
    if(!(resamplingMethod %in% c('default', 'multinomial', 'systematic', 'stratified',
                                 'residual')))
      stop('buildBootstrapFilter: "resamplingMethod" must be one of: "default", "multinomial", "systematic", "stratified", or "residual". ')
    ## latent state info
    nodes <- findLatentNodes(model, nodes, timeIndex)  
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1)
      stop('buildBootstrapFilter: sizes or dimensions of latent states varies.')
    vars <- model$getVarNames(nodes =  nodes)  

    my_initializeModel <- initializeModel(model, silent = silent)
    
    if(0 > thresh || 1 < thresh || !is.numeric(thresh))
      stop('buildBootstrapFilter: "thresh" must be between 0 and 1.')
    if(!saveAll & smoothing)
      stop("buildBootstrapFilter: must have 'saveAll = TRUE' for smoothing to work.")
    ## Create mv variables for x state and sampled x states.  If saveAll=TRUE, 
    ## the sampled x states will be recorded at each time point.
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]
    if(saveAll){
      
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      
      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                              types = type,
                                              sizes = size))
      
      names <- c(names, "wts")
      type <- c(type, "double")
      size$wts <- length(dims)
      ##  Only need one weight per particle (at time T) if smoothing == TRUE.
      if(smoothing)
        size$wts <- 1 
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                              types = type,
                                              sizes = size))
      
    }
    else{
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size[[1]] <- as.numeric(dims[[1]])
      
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
    }
    
    bootStepFunctions <- nimbleFunctionList(bootStepVirtual)
    for(iNode in seq_along(nodes)){
      bootStepFunctions[[iNode]] <- bootFStep(model, mvEWSamples, mvWSamples,
                                              nodes, iNode, names, saveAll,
                                              smoothing, resamplingMethod,
                                              silent) 
    }
    essVals <- rep(0, length(nodes))
    lastLogLik <- -Inf
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    if(initModel) my_initializeModel$run()
    resize(mvWSamples, m)
    resize(mvEWSamples, m)
    threshNum <- ceiling(thresh*m)
    logL <- 0   
    ## prevSamp indicates whether resampling took place at the previous time point.
    prevSamp <- 1
    for(iNode in seq_along(bootStepFunctions)) {
      if(iNode == length(bootStepFunctions))
        threshNum <- m  ## always resample at last time step so mvEWsamples is equally-weighted
      out <- bootStepFunctions[[iNode]]$run(m, threshNum, prevSamp)
      logL <- logL + out[1]
      prevSamp <- out[2]
      essVals[iNode] <<- bootStepFunctions[[iNode]]$returnESS()
      if(logL == -Inf) {lastLogLik <<- logL; return(logL)} 
      if(is.nan(logL)) {lastLogLik <<- -Inf; return(-Inf)}
      if(logL == Inf)  {lastLogLik <<- -Inf; return(-Inf)} 
    }
    lastLogLik <<- logL
    return(logL)
  },
  methods = list(
    getLastLogLik = function() {
      return(lastLogLik)
      returnType(double())
    },
    setLastLogLik = function(lll = double()) {
      lastLogLik <<- lll
    },
    returnESS = function(){
      returnType(double(1))
      return(essVals)
    }
  )
)


