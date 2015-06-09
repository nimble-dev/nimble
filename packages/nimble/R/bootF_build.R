##  Contains code to run bootstrap particle filters.
##  We have a build function (buildBootF),
##  and step functions.  There are two step functions: 
##  one which ends with NS (only saves latent state samples from
##  most recent time point), one which ends with S (saves latent state 
##  samples from all time points).



auxStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer(), thresh_num=double()) 
    returnType(double())
)


# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.

bootFStepNS <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mv, nodes, iNode, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    thisNode <- nodes[iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    isLast <- (iNode == length(nodes))
  },
  run = function(m = integer(), thresh_num=double()) {
    returnType(double())
    declare(wts, double(1, m))
    declare(ids, integer(1, m))
    declare(ess, double())
    declare(llEst, double(1,m))
    for(i in 1:m) {
      if(notFirst) {  
        copy(mv, model, nodes = 'xs', nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = 'x', row = i)
      calculate(model, thisDeterm)
      wts[i]  <- exp(calculate(model, thisData))
      if(notFirst){
        llEst[i] <- wts[i]*mv['wts',i][1]
      }
      else{
        llEst[i] <- wts[i]/m
      }
    }
    # Normalize weights and calculate effective sample size 
    wts <- wts/sum(wts)
    ess <- 1/sum(wts^2) 
    
    # Determine whether to resample by weights or not
    if(ess < thresh_num){
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mv, mv, nodes = 'x', nodesTo = 'xs', row = ids[i], rowTo = i)
        mv['wts',i][1] <<- 1/m
      }
    }
    else{
      for(i in 1:m){
        copy(mv, mv, "x", "xs", i, i)
        mv['wts',i][1] <<- wts[i]
      }
    }
    return(log(sum(llEst)))
  },  where = getLoadingNamespace()
)


bootFStepS <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mv, nodes, iNode, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    thisNode <- nodes[iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    t <- iNode  # current time point
    # Get names of xs node for current and previous time point (used in copy)
    thisXSName <- paste("xs[,", t, "]", sep = "") 
    prevXSName <- paste("xs[,", t-1, "]", sep = "")
    thisXName <- paste("x[,", t, "]", sep = "")
    isLast <- (iNode == length(nodes))
  },
  run = function(m = integer(), thresh_num=double()) {
    returnType(double())
    declare(wts, double(1, m))
    declare(ids, integer(1, m))
    declare(ess, double())
    declare(llEst, double(1,m))
    for(i in 1:m) {
      if(notFirst) {  
        copy(mv, model, nodes = prevXSName, nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = thisXName, row = i)
      calculate(model, thisDeterm)
      wts[i]  <- exp(calculate(model, thisData))
      if(notFirst){
        llEst[i] <- wts[i]*mv['wts',i][1,(t-1)]
      }
      else{
        llEst[i] <- wts[i]/m
      }
    }
    
    # Normalize weights and calculate effective sample size 
    wts <- wts/sum(wts)
    ess <- 1/sum(wts^2) 
    
    # Determine whether to resample by weights or not
    if(ess < thresh_num){
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mv, mv, nodes = thisXName, nodesTo = thisXSName, row = ids[i], rowTo = i)
        mv['wts',i][1,t] <<- 1/m
      }
    }
    else{
      for(i in 1:m){
        copy(mv, mv, nodes = thisXName, nodesTo = thisXSName, row = i, rowTo = i)
        mv['wts',i][1,t] <<- wts[i]
      }
    }
    return(log(sum(llEst)))
  },  where = getLoadingNamespace()
)

#' Creates a bootstrap particle filter algorithm to estimate log-likelihood.
#'
#' @param model A nimble model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes over which the particle filter will stochastically integrate over to estimate the log-likelihood function
#' @param thresh A number between 0 and 1 specifying when to resample: the resampling step will occur when the effective sample size is less than thresh*(number of particles)
#' @param saveAll  Whether to save state samples for all time points (T), or only for the most recent time points (F)
#' @author Daniel Turek
#' @details The resulting specialized particle filter algorthm will accept a
#'  single integer argument (m, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm 
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in mv['x',], with corresponding unlogged weights in mv['wts',].
#'  An equally weighted sample from the posterior can be found in mv['xs',]. 
#' @family filtering methods
#' @examples
#' model <- nimbleModel(code = ...)
#' my_BootF <- buildBootF(model, 'x[1:100]')
#' Cmodel <- compileNimble(model)
#' Cmy_BootF <- compileNimble(my_BootF, project = model)
#' logLike <- Cmy_BootF(m = 100000)
#' boot_X <- Cmy_BootF$mv['xs',]
#' @export
buildBootF <- nimbleFunction(
  setup = function(model, nodes, thresh = 0.5, silent = FALSE, saveAll = FALSE) {
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimbleDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent states varies')
    if(0>thresh || 1<thresh || !is.numeric(thresh)) stop('thresh must be between 0 and 1')
    
    # Create mv variables for x state and sampled x states. 
    # Size depends on whethersaveAll=T or not 
    if(!saveAll){
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs','wts'),  
                                        type = c('double', 'double','double'),
                                        size = list(x = dims[[1]],
                                                    xs = dims[[1]],
                                                    wts = dims[[1]])))
      bootStepFunctions <- nimbleFunctionList(auxStepVirtual)
      for(iNode in seq_along(nodes)){
        bootStepFunctions[[iNode]] <- bootFStepNS(model, mv, nodes, 
                                                  iNode, silent)  
      }
    }
    
    else{
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs', 'wts'),
                                        type = c('double', 'double', 'double'),
                                        size = list(x = c(dims[[1]],
                                                          length(dims)),
                                                    xs = c(dims[[1]],
                                                           length(dims)),
                                                    wts = c(dims[[1]],
                                                            length(dims))
                                        )))
      bootStepFunctions <- nimbleFunctionList(auxStepVirtual)
      for(iNode in seq_along(nodes)){
        bootStepFunctions[[iNode]] <- bootFStepS(model, mv, nodes, 
                                                 iNode, silent)   
      }
    }
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    my_initializeModel$run()
    resize(mv, m)
    thresh_num <- ceiling(thresh*m)
    logL <- 0    
    for(iNode in seq_along(bootStepFunctions)) { 
      logL <- logL + bootStepFunctions[[iNode]]$run(m,thresh_num)
      if(logL == -Inf) return(logL)
    }
    return(logL)
  },  where = getLoadingNamespace()
)
