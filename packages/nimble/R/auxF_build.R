##  Contains code to run auxiliary particle filters.
##  We have a build function (buildAuxF),
##  and step function (auxFStep)
##  The details of the specific APF algorithm are in
##  "A Tutorial on Particle Filtering and Smoothing:
##   Fifteen years later" by Doucet and Johansen

auxStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer(), thresh_num=double()) 
    returnType(double())
)




auxFStep <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mvEWSamp, mvWSamp, nodes, iNode, names, saveAll, silent) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
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
    }
    else{
      prevXName <- names    
      thisXName <- names
      currInd <- 1
      prevInd <- 1 
    }
    isLast <- (iNode == length(nodes))
    meanFuncList <- nimbleFunctionList(node_stoch_dnorm)
    meanFuncList[[1]] <- model$nodeFunctions[[thisNode]]
  },
  run = function(m = integer(), thresh_num = double()) {
    returnType(double())
    declare(auxl, double(1,m))
    declare(auxWts, double(1,m))
    declare(wts, double(1,m))
    declare(ids, integer(1, m))
    declare(l, double(1,m))
    declare(LL, double(1,m))
    ess <- 0
    resamp <- 0
    if(notFirst){ 
      for(i in 1:m) {
        copy(mvWSamp, model, prevXName, prevNode, row=i)        
        calculate(model, prevDeterm) 
        model[[thisNode]] <<-meanFuncList[[1]]$get_mean() # returns E(x_t+1 | x_t)
        calculate(model, thisDeterm)
        auxl[i] <- calculate(model, thisData)
        auxWts[i] <- auxl[i] + mvWSamp['wts',i][prevInd]        
      }
      auxWts <- exp(auxWts)/sum(exp(auxWts))
      ess <- 1/sum(auxWts^2)
      
      if(ess<thresh_num){
        resamp <- 1
        rankSample(auxWts, m, ids, silent)
        for(i in 1:m){
          if(saveAll == 1)
            copy(mvWSamp, mvEWSamp, prevXName, prevXName, ids[i],i)
          mvWSamp['wts',i][prevInd]   <<- log(1/m)
        }
      }
      else{
        for(i in 1:m){
          copy(mvWSamp, mvEWSamp, prevXName, prevXName, i,i)
          mvWSamp['wts',i][prevInd]   <<- log(auxWts[i] )     
        }
      }
    }   
    for(i in 1:m) {
      if(notFirst) {
        copy(mvEWSamp, model, nodes = prevXName, nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mvWSamp, nodes = thisNode, nodesTo = thisXName, row=i)
      calculate(model, thisDeterm)
      l[i]  <- calculate(model, thisData)
      if(notFirst){
        if(resamp == 1){
          mvWSamp['wts',i][currInd] <<- l[i]-auxl[ids[i]]
          LL[i] <- l[i]-auxl[ids[i]]
        }
        else{
          mvWSamp['wts',i][currInd] <<- l[i] - auxl[i]
          LL[i] <- l[i] - auxl[i]
        }
      }
      
      else{
        mvWSamp['wts',i][currInd] <<- l[i]
        LL[i] <- l[i]
      }
    }
    
    if(isLast){
      for(i in 1:m){
        wts[i] <-  exp(mvWSamp['wts',i][currInd] )
      }
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mvWSamp, mvEWSamp, thisXName, thisXName, ids[i], i)
      }
    }
    
    return(log(mean(exp(LL))))
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
#'  \item{"thresh"}{ A number between 0 and 1 specifying when to resample: the resampling step will occur when the
#'   effective sample size is less than thresh*(number of particles).  Defaults to 0.5.}
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
#'   The auxiliary particle filter uses a lookeahead function to select promising particles before propogation.  Currently, the lookahead
#'   funciton uses the expected valu of the latent state at the next time point given the current particle, e E(x[t+1]|x[t]).
#'   The auxiliary particle filter currently only works for models with univariate normal transition densities. 
#'   @references Pitt, Michael K., and Neil Shephard. "Filtering via simulation: Auxiliary particle filters."
#'    Journal of the American statistical association 94.446 (1999): 590-599.
#' @examples
#' model <- nimbleModel(code = ...)
#' my_AuxF <- buildAuxF(model, 'x[1:100]', control = list(thresh = 0.9))
#' Cmodel <- compileNimble(model)
#' Cmy_AuxF <- compileNimble(my_AuxF, project = model)
#' logLike <- Cmy_AuxF(m = 100000)
#' hist(as.matrix(Cmy_Auxf$mvEWSamp, 'x'))
#' @export
buildAuxF <- nimbleFunction(
  setup = function(model, nodes, control = list()) {
    
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    vars <- model$getVarNames(nodes =  nodes)  # need var names too
    
    if(length(unique(dims)) > 1) 
      stop('sizes or dimension of latent states varies')
    
    
    thresh <- control[['thresh']]
    saveAll <- control[['saveAll']]
    silent <- control[['silent']]
    if(is.null(silent)) silent <- FALSE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(thresh)) thresh <- .5
    
    if(thresh<0 || thresh>1 || !is.numeric(thresh)) 
      stop('thresh must be between 0 and 1')
    
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
                                            iNode, names, saveAll, silent)
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    my_initializeModel$run()
    resize(mvEWSamp, m) 
    resize(mvWSamp, m)  
    thresh_num <- ceiling(thresh*m)
    logL <- 0
    for(iNode in seq_along(auxStepFunctions)) {
      logL <- logL + auxStepFunctions[[iNode]]$run(m, thresh_num)      
      if(logL == -Inf) return(logL) }
    return(logL)
  },  where = getLoadingNamespace()
)

