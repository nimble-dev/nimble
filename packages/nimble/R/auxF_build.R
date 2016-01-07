##  Contains code to run auxiliary particle filters.
##  We have a build function (buildAuxF),
##  and step function (auxFStep)
##
##  This version of the APF is based on 
##  Pitt et al., 2012

auxStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer()) 
    returnType(double())
)

auxFuncVirtual <- nimbleFunctionVirtual(
  methods = list(
    lookahead = function(){}
  )
)

auxMLookFunc = nimbleFunction(
  contains = auxFuncVirtual,
  setup = function(model, node){
    nodeFuncList <- nimbleFunctionList(node_stoch_dmnorm)
    nodeFuncList[[1]] <- model$nodeFunctions[[node]]
  },
  methods = list(
    lookahead = function(){
      model[[node]] <<- nodeFuncList[[1]]$get_mean()
    }
  ), where = getLoadingNamespace()
)

auxLookFunc = nimbleFunction(
  contains = auxFuncVirtual,
  setup = function(model, node){
    nodeFuncList <- nimbleFunctionList(node_stoch_dnorm) 
    nodeFuncList[[1]] <- model$nodeFunctions[[node]]
  },
  methods = list(
    lookahead = function(){
      model[[node]] <<- nodeFuncList[[1]]$get_mean()
    }
  ), where = getLoadingNamespace()
)

auxSimFunc = nimbleFunction(
  contains = auxFuncVirtual,
  setup = function(model, node){},
  methods = list(
    lookahead = function(){
      simulate(model, node)
    }), where = getLoadingNamespace()
)

auxFStep <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mvEWSamp, mvWSamp, nodes, iNode, names, saveAll, smoothing,
                   lookahead, silent = TRUE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    allPrevNodes <- nodes[1:(iNode-1)]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE) 
    #current time point
    t <- iNode
    # Get names of x and xs node for current and previous time point,
    # will be different depending on whether we are saving all time points
    # or only the most recent
    if(saveAll == 1){
      prevXName <- prevNode    
      thisXName <- thisNode
      currInd <- t
      prevInd <- t-1
      if(smoothing == T){
        currInd <- 1
        prevInd <- 1
      }
    }
    else{
      prevXName <- names    
      thisXName <- names
      currInd <- 1
      prevInd <- 1 
    }
    isLast <- (iNode == length(nodes))
    xDim <- nimDim(model[[thisXName]])
    
    auxFuncList <- nimbleFunctionList(auxFuncVirtual) 
    if(lookahead == "mean"){
    if(xDim > 1){
      auxFuncList[[1]] <- auxMLookFunc(model, thisNode)
    }
    else{
      auxFuncList[[1]] <- auxLookFunc(model, thisNode)
    } 
    }
    else{
      auxFuncList[[1]] <- auxSimFunc(model, thisNode)
    }
  },
  run = function(m = integer()) {
    returnType(double())
    declare(auxll, double(1,m))
    declare(auxWts, double(1,m))
    declare(normAuxWts, double(1,m))
    declare(wts, double(1,m))
    declare(normWts, double(1,m))
    declare(ids, integer(1, m))
    declare(ll, double(1,m))
    declare(LL, double(1,m))

    ## This is the look-ahead step, not conducted for first time-point
    if(notFirst){ 
      for(i in 1:m) {
        if(smoothing == 1){
          copy(mvEWSamp, mvWSamp, nodes = allPrevNodes, nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        copy(mvWSamp, model, prevXName, prevNode, row=i)        
        calculate(model, prevDeterm)
        auxFuncList[[1]]$lookahead()
        calculate(model, thisDeterm)
        auxll[i] <- calculate(model, thisData)  # get p(y_t+1 | x_t+1)
        if(is.nan(auxll[i])){
          return(-Inf)
        }
        auxll[i] <- auxll[i]+calculate(model, thisNode) # multiply by p(x_t+1 | x_t)
        auxWts[i] <- auxll[i] + mvWSamp['wts',i][prevInd] # multiply by weight from time t

      }
      normAuxWts <- exp(auxWts)/sum(exp(auxWts))  # normalize weights and resample
      rankSample(normAuxWts, m, ids, silent)
    }   
    for(i in 1:m) {
      if(notFirst) {
        copy(mvWSamp, model, nodes = prevXName, nodesTo = prevNode, row = ids[i])
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)  # simulate from x_t+1 | x_t
      copy(model, mvWSamp, nodes = thisNode, nodesTo = thisXName, row=i)
      calculate(model, thisDeterm)
      ll[i]  <- calculate(model, thisData)  # get p(y_t+1 | x_t+1)
      if(is.nan(ll[i])){
        return(-Inf)
      }
      if(notFirst){
        wts[i] <- ll[i]-auxll[ids[i]]  # construct weight following step 4 of paper    
      }
      else{
        wts[i] <- ll[i]  # First step has no auxiliary weights
      }
    }
    
    normWts <- exp(wts)/sum(exp(wts))
    
    for(i in 1:m){
      ##  save weights for use in next timepoint's look-ahead step
      mvWSamp['wts', i][currInd] <<- log(normWts[i])   
    }
    rankSample(normWts, m, ids, silent)
    for(i in 1:m){
      if(smoothing == 1){
        copy(mvWSamp, mvEWSamp, nodes = allPrevNodes, nodesTo = allPrevNodes, row = ids[i], rowTo=i)
      }
      copy(mvWSamp, mvEWSamp, thisXName, thisXName, ids[i], i)
    }
    
    
    ##  Calculate likelihood p(y_t+1 | y_1:t) as in equation (3) of paper
    if(notFirst){
      outLL <- sum(exp(wts))/m
      outLL <- outLL*sum(exp(auxWts))
    }
    else{
      outLL <- sum(exp(wts))/m
    }
    return(log(outLL))
  }, where = getLoadingNamespace()
)


#' Creates an auxiliary particle filter algorithm to estimate log-likelihood.
#'
#' @param model A nimble model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes over which the particle filter will stochastically integrate over to estimate the log-likelihood function
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author  Nicholas Michaud
#' @family particle filtering methods
#' @details 
#' 
#' \describe{
#'  \item{"saveAll"}{Indicates whether to save state samples for all time points (T), or only for the most recent time point (F)}
#' }
#' 
#'      The auxiliary particle filter modifies the bootstrap filter (\code{\link{buildBootF}})
#'      by adding a lookahead step to 
#'      the algorithm: before propagating particles from one time point to the next via the transition equation, 
#'      the auxiliary filter calculates a weight for each pre-propogated particle by predicting how well
#'      the particle will agree with the next data point.  These pre-weights are used to conduct an initial 
#'      resampling step before propagation. 
#' 
#'  The resulting specialized particle filter algorthm will accept a
#'  single integer argument (m, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm 
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the mvWSamp model values object, with corresponding logged weights in mvWSamp['wts',].
#'  An equally weighted sample from the posterior can be found in mvEWsamp.  
#'  
#'   The auxiliary particle filter uses a lookahead function to select promising particles before propogation.  Currently, the lookahead
#'   funciton uses the expected value of the latent state at the next time point given the current particle, e E(x[t+1]|x[t]).
#'   The auxiliary particle filter currently only works for models with univariate normal transition densities. 
#'   @references Pitt, Michael K., and Neil Shephard. "Filtering via simulation: Auxiliary particle filters."
#'    Journal of the American statistical association 94.446 (1999): 590-599.
#'   @references Pitt, Michael K., et al. "On some properties of Markov chain Monte Carlo simulation methods based on the particle filter." 
#'   Journal of Econometrics 171.2 (2012): 134-151.
#' @examples
#' model <- nimbleModel(code = ...)
#' my_AuxF <- buildAuxF(model, 'x[1:100]', control = list(saveAll = T))
#' Cmodel <- compileNimble(model)
#' Cmy_AuxF <- compileNimble(my_AuxF, project = model)
#' logLike <- Cmy_AuxF(m = 100000)
#' hist(as.matrix(Cmy_Auxf$mvEWSamp, 'x'))
#' @export
buildAuxF <- nimbleFunction(
  setup = function(model, nodes, control = list()) {
    
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    vars <- model$getVarNames(nodes =  nodes)  # need var names too
    
    if(length(unique(dims)) > 1) 
      stop('sizes or dimension of latent states varies')
    
    saveAll <- control[['saveAll']]
    silent <- control[['silent']]
    smoothing <- control[['smoothing']]
    lookahead <- control[['lookahead']]
    
    if(is.null(silent)) silent <- TRUE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(lookahead)) lookahead <- "simulate"
    
    my_initializeModel <- initializeModel(model, silent = silent)
    
    
    if(lookahead == "mean"){
      errors <- sapply(model$expandNodeNames(nodes), function(node){tryCatch(model$nodes[[node]]$get_mean(), error=function(a){return("error")})})
      if(any(errors == "error", na.rm=T)) stop("transition equation must be normal to use mean as lookahead function")
    } 
    else if(lookahead != "simulate"){
      stop("lookahead must be either simulate or mean")
    }
    
    # Create mv variables for x state and sampled x states.  If saveAll=T, 
    # the sampled x states will be recorded at each time point. 
    modelSymbolObjects = model$getSymbolTable()$getSymbolObjects()[vars]
    if(saveAll){
      
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      mvEWSamp <- modelValues(modelValuesSpec(vars = names,
                                              type = type,
                                              size = size))
      
      names <- c(names, "wts")
      type <- c(type, "double")
      size$wts <- length(dims)
      if(smoothing == T){
        size$wts <- 1 ##  only need one weight per particle (at time T) if smoothing == T
      }
      mvWSamp  <- modelValues(modelValuesSpec(vars = names,
                                              type = type,
                                              size = size))
      
    }
    else{
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size[[1]] <- as.numeric(dims[[1]])
      
      
      mvEWSamp <- modelValues(modelValuesSpec(vars = names,
                                              type = type,
                                              size = size))
      
      names <- c(names, "wts")
      type <- c(type, "double")
      size$wts <- 1
      mvWSamp  <- modelValues(modelValuesSpec(vars = names,
                                              type = type,
                                              size = size))
    }
    
    names <- names[1]
    auxStepFunctions <- nimbleFunctionList(auxStepVirtual)
    for(iNode in seq_along(nodes))
      auxStepFunctions[[iNode]] <- auxFStep(model, mvEWSamp, mvWSamp, nodes,
                                            iNode, names, saveAll, smoothing, lookahead, silent)
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    my_initializeModel$run()
    resize(mvEWSamp, m) 
    resize(mvWSamp, m)  
    logL <- 0
    for(iNode in seq_along(auxStepFunctions)) {
      logL <- logL + auxStepFunctions[[iNode]]$run(m)

      # when all particles have 0 weight, likelihood becomes NAN
      # this happens if top-level params have bad values - possible
      # during pmcmc for example
      if(is.nan(logL)) return(-Inf)
      if(logL == -Inf) return(logL) 
      if(logL == Inf) return(-Inf) 
    }
  
    return(logL)
  },  where = getLoadingNamespace()
)

