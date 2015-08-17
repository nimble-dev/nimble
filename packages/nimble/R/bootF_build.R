##  Contains code to run bootstrap particle filters.
##  We have a build function (buildBootF),
##  and step function.


bootStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer(), threshNum=double(), prevSamp = logical()) 
    returnType(double(1))
)


# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.

bootFStep <- nimbleFunction(
  contains = bootStepVirtual,
  setup = function(model, mv, nodes, iNode, saveAll, smoothing, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    thisNode <- nodes[iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    t <- iNode  # current time point
    # Get names of xs node for current and previous time point (used in copy)
    if(saveAll == 1){
      prevXSName <- paste("xs[,",t-1,"]",sep="")    
      thisXSName <- paste("xs[,",t,"]",sep="")
      prevXName <- paste("x[,",t-1,"]",sep="")
      thisXName <- paste("x[,",t,"]",sep="")
      currInd <- t
      prevInd <- t-1
    }
    else{
      prevXSName <- "xs"
      thisXSName <- "xs"
      prevXName <- "x"
      thisXName <- "x"
      currInd <- 1
      prevInd <- 1 
    }
    isLast <- (iNode == length(nodes))
  },
  run = function(m = integer(), threshNum = double(), prevSamp = logical()) {
    returnType(double(1))
    declare(wts, double(1, m))
    declare(ids, integer(1, m))
    declare(ess, double())
    declare(llEst, double(1,m))
    declare(out, double(1,2))
    
    for(i in 1:m) {
      if(notFirst) {  
        if(smoothing == 1){
          copy(mv, mv, nodes = 'xs', nodesTo = 'x', row = i)
        }
        copy(mv, model, nodes = prevXSName, nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      ## below line is used to get around the R version of this funciton shrinking mv[wts, i] to a one dimensional object when saveAll ==F,
      ## when it must remain two dimensional.  We add an unnecessary second column, which we fill with 0's
      if(!notFirst & !saveAll) mv['wts', i][1,2]  <<- 0
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = thisXName, row = i)
      calculate(model, thisDeterm)
      wts[i]  <- exp(calculate(model, thisData))
      if(prevSamp == 0){
        llEst[i] <- wts[i]*mv['wts',i][1,prevInd]
      }
      else{
        llEst[i] <- wts[i]/m
      }
    }
    
    out[1] <- log(sum(llEst))
    # Normalize weights and calculate effective sample size 
    wts <- wts/sum(wts)
    ess <- 1/sum(wts^2) 
    
    # Determine whether to resample by weights or not
    if(ess < threshNum){
      rankSample(wts, m, ids, silent)
      # out[2] is an indicator of whether resampling takes place
      # affects how ll estimate is calculated at next time point.
      out[2] <- 1
      for(i in 1:m){
        if(smoothing == 0)
          copy(mv, mv, nodes = thisXName, nodesTo = thisXSName, row = ids[i], rowTo = i)
        else
          copy(mv, mv, nodes = 'x', nodesTo = 'xs', row = ids[i], rowTo = i)
        mv['wts',i][1,currInd] <<- wts[i]
      }
    }
    else{
      out[2] <- 0
      for(i in 1:m){
        if(smoothing == 0)
          copy(mv, mv, nodes = thisXName, nodesTo = thisXSName, row = i, rowTo = i)
        else
          copy(mv, mv, nodes = 'x', nodesTo = 'xs', row = i, rowTo = i)
        mv['wts',i][1,currInd] <<- wts[i]
      }
    }
    return(out)
  }, where = getLoadingNamespace()
)

#' Creates a bootstrap particle filter algorithm to estimate log-likelihood.
#'
#' @param model A nimble model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes over which the particle filter will stochastically integrate over to estimate the log-likelihood function
#' @param thresh A number between 0 and 1 specifying when to resample: the resampling step will occur when the effective sample size is less than thresh*(number of particles)
#' @param saveAll  Whether to save state samples for all time points (T), or only for the most recent time points (F)
#' @param filtering  If saveAll = T, should we save filtered estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]), or save samples from f(x[t]|y[t]) at each time point.
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
  setup = function(model, nodes, filterControl = list(thresh = 0.5,  silent = FALSE, saveAll = FALSE, smoothing = FALSE)) {
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent states varies')
    
    thresh <- filterControl[['thresh']]
    saveAll <- filterControl[['saveAll']]
    smoothing <- filterControl[['smoothing']]
    silent <- filterControl[['silent']]
    if(is.null(silent)) silent <- FALSE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(thresh)) thresh <- .5
    
    if(0>thresh || 1<thresh || !is.numeric(thresh)) stop('thresh must be between 0 and 1')
    if(!saveAll) smoothing <- FALSE
    # Create mv variables for x state and sampled x states.  If saveAll=T, 
    # the sampled x states will be recorded at each time point. 
    # Note: the algorithm will no longer run in R, since R automatically
    # reduces the dimension of the "wts"  variable in the mv, which
    # breaks the step function.  C keeps 2 dimensions and works fine.
    if(!saveAll){
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs','wts'),  
                                        type = c('double', 'double','double'),
                                        size = list(x = dims[[1]],
                                                    xs = dims[[1]],
                                                    wts =  c(1,2))))
    }
    
    else{
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs', 'wts'),
                                        type = c('double', 'double', 'double'),
                                        size = list(x = c(dims[[1]],
                                                          length(dims)),
                                                    xs = c(dims[[1]],
                                                           length(dims)),
                                                    wts = c(1,
                                                            length(dims)))))
     
    }
    
    bootStepFunctions <- nimbleFunctionList(bootStepVirtual)
    for(iNode in seq_along(nodes)){
      bootStepFunctions[[iNode]] <- bootFStep(model, mv, nodes,
                                               iNode, saveAll, smoothing, silent) 
    }
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(out, double(1,2))
    my_initializeModel$run()
    resize(mv, m)
    threshNum <- ceiling(thresh*m)
    logL <- 0   
    # prevSamp indicates whether resampling took place at the
    # previous time point.
    prevSamp <- 1
    for(iNode in seq_along(bootStepFunctions)) { 
      out <- bootStepFunctions[[iNode]]$run(m,threshNum, prevSamp)
      logL <- logL + out[1]
      prevSamp <- out[2]
      if(logL == -Inf) return(logL)
    }
    return(logL)
  }, where = getLoadingNamespace()
)
