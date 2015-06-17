##  Contains code to run auxiliary particle filters.
##  We have a build function (buildAuxF),
##  and step function (auxFStep.  



auxStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer(), thresh_num=double()) 
    returnType(double())
)




auxFStep <- nimbleFunction(
  contains = auxStepVirtual,
  setup = function(model, mv, nodes, iNode, saveAll, silent = FALSE) {
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
    if(saveAll==1){
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
    ess <- 0
    if(notFirst){
      for(i in 1:m) {
        copy(mv, model, prevXName, prevNode, row=i)
        calculate(model, prevDeterm) 
        model[[thisNode]] <<-meanFuncList[[1]]$get_mean() # returns E(x_t+1 | x_t)
        calculate(model, thisDeterm)
        auxl[i] <- exp(calculate(model, thisData))
        auxWts[i] <- auxl[i]*mv['wts',i][1,prevInd]
      }
      auxWts <- auxWts/sum(auxWts)
      ess <- 1/sum(auxWts^2)
      if(ess<thresh_num){
        resamp <- 1
        rankSample(auxWts, m, ids, silent)
        for(i in 1:m){
          copy(mv, mv, prevXName, prevXSName, ids[i],i)
        }
      }
      else{
        resamp <- 0
        for(i in 1:m){
          copy(mv, mv, prevXName, prevXSName, i,i)
        }
      }
    }   
    for(i in 1:m) {
      if(notFirst) {
        copy(mv, model, nodes = prevXSName, nodesTo = prevNode, row = i)
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = thisXName, row = i)
      calculate(model, thisDeterm)
      l[i]  <- exp(calculate(model, thisData))
      if(notFirst){
        if(resamp == 1){
          mv['wts',i][1,currInd] <<- l[i]/auxl[ids[i]]
        }
        else{
          mv['wts',i][1,currInd] <<- l[i]/auxl[i]
        }
      }
      
      else{
        mv['wts',i][1,currInd] <<- l[i]
      }
    }
    
    if(isLast){
      for(i in 1:m){
        wts[i] <-  mv['wts',i][1,currInd] 
      }
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mv, mv, thisXName, thisXSName, ids[i], i)
      }
    }
    return(log(mean(l)))
  }, where = getLoadingNamespace()
)


#' Creates an auxiliary particle filter 
#' 
#' @param model A nimble model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes over which the Auxiliary particle filter will stochastically integrate over to estimate the log-likelihood function
#' @param thresh A number between 0 and 1 specifying when to resample: the resampling step will occur when the effective sample size is less than thresh*(number of particles)
#' @param saveAll  Whether to save state samples for all time points (T), or only for the most recent time points (F)
#' @author Nick Michaud
#' @family filtering methods
#' @details The resulting specialized particle filter algorthm will accept a
#'  single integer argument (m, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm 
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in mv['x',], with corresponding unlogged weights in mv['wts',].
#'  An equally weighted sample from the posterior can be found in mv['xs',]. 
#' @examples
#' model <- nimbleModel(code = ...)
#' my_AuxF <- buildAuxF(model, 'x[1:100]', .9)
#' Cmodel <- compileNimble(model)
#' Cmy_AuxF <- compileNimble(my_AuxF, project = model)
#' logLike <- Cmy_AuxF(m = 100000)
#' hist(as.matrix(Cmy_Auxf$mv, 'xs'))
#' @export
buildAuxF <- nimbleFunction(
  setup = function(model, nodes, thresh=.5, silent = FALSE, saveAll = FALSE) {
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimbleDim(model[[n]]))
    if(length(unique(dims)) > 1) 
      stop('sizes or dimension of latent states varies')
    if(thresh<0 || thresh>1 || !is.numeric(thresh)) 
      stop('thresh must be between 0 and 1')
    
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
                                                    wts = c(dims[[1]],1))))
    }
    else{
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs','wts'),  
                                        type = c('double', 'double','double'),
                                        size = list(x = c(dims[[1]],
                                                          length(dims)),
                                                    xs = c(dims[[1]],
                                                           length(dims)),
                                                    wts = c(dims[[1]],
                                                               length(dims))))) 
    }
    auxStepFunctions <- nimbleFunctionList(auxStepVirtual)
    for(iNode in seq_along(nodes))
      auxStepFunctions[[iNode]] <- auxFStep(model, mv, nodes,
                                             iNode, saveAll, silent)
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    my_initializeModel$run()
    resize(mv, m)  
    thresh_num <- ceiling(thresh*m)
    logL <- 0
    for(iNode in seq_along(auxStepFunctions)) {
      logL <- logL + auxStepFunctions[[iNode]]$run(m, thresh_num)      
      if(logL == -Inf) return(logL) }
    return(logL)
  },  where = getLoadingNamespace()
)

