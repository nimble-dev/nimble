##  Contains code to run a Liu and West filter.  The filter has a build
##  function (buildLWF) and a step function (LWstepNS). Also contains 
##  a function for caluclating useful quantities. The doPars() function
##  is specialized to different variables in the mv objects, and allows
##  these variables to be set / retrieved within the mv objects without
##  explicitly calling the variable name.
##  
##
##  Still to implement:
##  option to skip auxiliary step in filter (Ionides '11), which will also allow for threshold resampling




LWStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer())
    returnType(double())
)



# Has a return_mean method which returns the mean of a normally distributed nimble node.
normMean <- nimbleFunction(
  setup = function(model, node) {
    nfList <- nimbleFunctionList(node_stoch_dnorm)
    nfList[[1]] <- model$nodeFunctions[[node]]
  },  
  methods = list(                        
    return_mean = function() {         
      returnType(double())           
      return(nfList[[1]]$get_mean()) 
    }                                  
  ), where = getLoadingNamespace()                                    
)

LWSetParVirtual <- nimbleFunctionVirtual(
  methods = list(
    scalarSet = function(scalars = double(1), mvWset = integer(), m = integer()){},
    vectorSet = function(vectors = double(2), mvWset = integer(), m = integer()){},
    scalarGet = function(mvWset = integer(), m = integer()){returnType(double(1))},
    vectorGet = function(mvWset = integer(), m = integer(), length = integer()){returnType(double(2))}
  )
)

doPars <- nimbleFunction(
  contains = LWSetParVirtual,
  setup = function(parName, mvWSamp, mvEWSamp) {
  },
  methods = list(
    scalarSet = function(scalars = double(1), mvWset = integer(), m = integer()){
      if(mvWset == 1){
        for(i in 1:m){
          mvWSamp[parName, i][1] <<- scalars[i]
        }
      }
      else{
        for(i in 1:m){
          mvEWSamp[parName, i][1] <<- scalars[i]
        }
      }
    },
    vectorSet = function(vectors = double(2), mvWset = integer(), m = integer()){
      if(mvWset == 1){
        for(i in 1:m){
          mvWSamp[parName, i] <<- vectors[,i]
        }
      }
      else{
        for(i in 1:m){
          mvEWSamp[parName, i] <<- vectors[,i]
        }
      }
    },
    scalarGet = function(mvWset = integer(), m = integer()){
      declare(vecOut, double(1,m))
      if(mvWset == 1){
        for(i in 1:m){
          vecOut[i] <- mvWSamp[parName, i][1]
        }
      }
      else{
        for(i in 1:m){
          vecOut[i] <- mvEWSamp[parName, i][1]
        }
      }
      returnType(double(1))
      return(vecOut)
    },
    vectorGet = function(mvWset = integer(), m = integer(), length = integer()){
      declare(matOut, double(2,c(length,m)))
      if(mvWset == 1){
        for(i in 1:m){
          matOut[,i] <- mvWSamp[parName, i]
        }
      }
      else{
        for(i in 1:m){
          matOut[,i] <- mvEWSamp[parName, i]
        }
      }
      returnType(double(2))
      return(matOut)
    } 
  ), where = getLoadingNamespace())


LWStep <- nimbleFunction(
  contains = LWStepVirtual,
  setup = function(model, mvWSamp, mvEWSamp, nodes, paramVarDims, iNode, paramNodes, paramVars, names, saveAll, d, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    parDeterm <- model$getDependencies(paramNodes, determOnly=TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    getmean <- normMeanA(model, thisNode)
    isLast <- (iNode == length(nodes))
    
    t <- iNode  # current time point
    # Get names of xs node for current and previous time point (used in copy)
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
    numParams <- sum(paramVarDims)
    numParamVars <- length(paramVars)
    paramInds <- c(0,cumsum(paramVarDims)) 
    parCalc <- LWparFunc(d, numParams, prevInd) #d is discount factor, default to 0.99  
    
    doVarList <- nimbleFunctionList(LWSetParVirtual)
    for(i in 1:numParamVars)
      doVarList[[i]] <- doPars(paramVars[i], mvWSamp, mvEWSamp)
  },
  run = function(m = integer()) {
    returnType(double())
    declare(auxWts, double(1,m))
    declare(l, double(1,m))
    declare(wts, double(1, m))
    declare(ids, integer(1, m))
    declare(preWts, double(1,m))
    declare(tmpPars, double(2, c(numParams, m)))
    ess <- 0
    
    if(notFirst){
      for(j in 1:numParamVars){
        tmpSize <- paramInds[j+1]-paramInds[j]
        if(tmpSize == 1){
          tmpPars[paramInds[j+1], ] <- doVarList[[j]]$scalarGet(1,m)
        }
        else{
          tmpPars[(paramInds[j]+1):paramInds[j+1], ] <- doVarList[[j]]$vectorGet(1, m, tmpSize)
        }
      }
      for(i in 1:m){
        wts[i] <- mvWSamp['wts',i][prevInd]
      }
      meanVec <- parCalc$shrinkMean(m, wts, tmpPars )  # shrink parameter particles towards mean  
      for(i in 1:m) {
        copy(mvWSamp, model, prevXName, prevNode, row=i)
        values(model, paramNodes) <<- meanVec[,i]
        calculate(model, parDeterm) 
        calculate(model, prevDeterm) 
        model[[thisNode]] <<- getmean$return_mean()
        calculate(model, thisDeterm)
        auxWts[i] <- exp(calculate(model, thisData))
        if(is.nan(auxWts[i])) auxWts[i] <- 0 #check for ok param values
        preWts[i] <- exp(log(auxWts[i])+wts[i]) #first resample weights
      }
      rankSample(preWts, m, ids, silent)
      # Reassign weights and pars so that cov matrix is calculated correctly
      for(i in 1:m){
        copy(mvWSamp, mvEWSamp, nodes = prevXName, nodesTo = prevXName, row = ids[i], rowTo = i)
        wts[i] <- mvWSamp['wts',ids[i]][prevInd] 
      }
      for(j in 1:numParamVars){
        tmpSize <- paramInds[j+1]-paramInds[j]
        if(tmpSize == 1){
          tmpPars[paramInds[j+1], ] <- doVarList[[j]]$scalarGet(1,m)
        }
        else{
          tmpPars[(paramInds[j]+1):paramInds[j+1], ] <- doVarList[[j]]$vectorGet(1, m, tmpSize)
        }
      }
      parChol <- parCalc$cholesVar(m,wts,tmpPars) #  calculate MC cov matrix
    }
    for(i in 1:m) {
      if(notFirst) {  
        tmpPars[,i] <- rmnorm_chol(1, meanVec[,ids[i]],parChol, prec_param=0) 
        values(model, paramNodes) <<- tmpPars[,i]
        copy(mvWSamp, model, nodes = prevXName, nodesTo = prevNode, row = ids[i]) 
      }
      else{
        simulate(model, paramNodes)
        tmpPars[,i] <- values(model,paramNodes)
      }
      calculate(model, parDeterm) 
      simulate(model, thisNode)
      copy(model, mvWSamp, nodes = thisNode, nodesTo = thisXName, rowTo = i)
      calculate(model, thisDeterm)
      l[i]  <- exp(calculate(model, thisData))
      if(is.nan(l[i])) l[i] <- 0
      ess <- ess+(l[i])^2
      #rescale weights by pre-sampling weight
      if(notFirst)
        wts[i] <- log(l[i])-log(auxWts[ids[i]])
      else 
        wts[i] <- log(l[i])
      mvWSamp['wts',i][currInd] <<- wts[i] 
    }
    #     #normalizing weights and calculating effective sample size 
    ess <- 1/(ess/(sum(l)^2))
    
    for(j in 1:numParamVars){
      tmpSize <- paramInds[j+1]-paramInds[j]
      if(tmpSize == 1){
        doVarList[[j]]$scalarSet(tmpPars[paramInds[j+1], ], 1 ,m)
      }
      else{
        doVarList[[j]]$vectorSet(tmpPars[(paramInds[j]+1):paramInds[j+1], ], 1 ,m)
      }
    }   
    
    if(isLast){
      for(i in 1:m){
        wts[i] <-  exp(mvWSamp['wts',i][currInd])
      }
      rankSample(wts, m, ids, silent)
      for(i in 1:m){
        copy(mvWSamp, mvEWSamp, thisXName, thisXName, ids[i], i)
      }
      for(j in 1:numParamVars){
        tmpSize <- paramInds[j+1]-paramInds[j]
        if(tmpSize == 1){
          doVarList[[j]]$scalarSet(tmpPars[paramInds[j+1], ids[1:m]], 0 ,m)
        }
        else{
          doVarList[[j]]$vectorSet(tmpPars[(paramInds[j]+1):paramInds[j+1], ids[1:m]], 0 ,m)
        }
      } 
    }
    return(log(mean(l)))
  },  where = getLoadingNamespace()
)

# Has two methods: shrinkMean, which shrinks each parameter particle towards
# the mean of all particles, and cholesVar, which returns the cholesky 
# decomposition of the weighted MC covariance matrix
LWparFunc <- nimbleFunction(
  setup = function(d, parDim, prevInd){
    # Calculate h^2 and a using specified discount factor d
    hsq <- 1-((3*d-1)/(2*d))^2 
    a <- sqrt(1-hsq)
  },  
  methods = list(                        
    shrinkMean = function(m = integer(), wts = double(1), pars = double(2)) {
      declare(parMean, double(2))
      declare(wtMean,double(2, c(parDim, m)))
      declare(oneMat, double(1, m))
      for(i in 1:m){
        oneMat[i] <- 1
      }
      wts <- wts/sum(wts)
      parMean <- (pars%*%wts)%*%oneMat
      wtMean <- a*pars + (1-a)*parMean  
      returnType(double(2))
      return(wtMean)
    },
    cholesVar = function(m = integer(), wts = double(1), pars = double(2)){
      declare(varMat, double(2, c(parDim, parDim)))
      returnType(double(2)) 
      wts <- wts/sum(wts)
      wMean <- pars%*%wts 
      for(i in 1:m)
        varMat <- varMat + wts[i]*((pars[,i]-wMean)%*%t(pars[,i]-wMean))
      varMat <- hsq*varMat
      return(chol(varMat))
    }
  ), where = getLoadingNamespace() 
)

#' Creates a Liu and West filter.  
#'
#' @param model A nimble model object, typically representing a state 
#'  space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes 
#'  over which the particle filter will stochastically integrate over to
#'  estimate the log-likelihood function
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.

#' @author Nicholas Michaud
#' @family particle filtering methods
#' @details 
#' 
#' Each of the control() list options are described in detail below:
#' \describe{
#'  \item{"params"}{A character vector sepcifying the parameters you would 
#'  like to estimate the posterior distribution of.  If unspecified, parameter nodes are specified as all stochastic top level nodes which
#'  are not in the set of latent nodes specified in 'nodes'.}
#'  \item{"d"}{A discount factor for the Liu-West filter.  Should be close to,
#'  but not above, 1.}
#'  \item{"saveAll"}{Indicates whether to save state samples for all time points (T), or only for the most recent time point (F)}
#' }
#' 
#'  The Liu and West filter samples from the posterior 
#'  distribution of both the latent states and top-level parameters for a state space model.  
#'  Each particle in the Liu and West filter contains values not only for latent states, 
#'  but also for top level parameters.  Latent states are propogated via an auxiliary step, 
#'  as in the auxiliary particle filter (\code{\link{buildAuxF}}).
#'  Top-level parameters are propagated from one 
#'  time point to the next through a smoothed kernel density based on previous particle values.  
#'        
#'  The resulting specialized particle filter algorthm will accept a
#'  single integer argument (m, default 10,000), which specifies the number
#'  of random \'particles\' to use for sampling from the posterior distributions.  The algorithm  saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states and top-level parameters in mvWSamp, with corresponding logged weights in mvWSamp['wts',].
#'  An equally weighted sample from the posterior can be found in mvEWSamp. 
#'  
#' @references Liu, Jane, and Mike West. "Combined parameter and state estimation in simulation-based filtering." 
#' Sequential Monte Carlo methods in practice. Springer New York, 2001. 197-223.
#' 
#' @examples
#' model <- nimbleModel(code = ...)
#' my_LWF <- buildPF(model, 'x[1:100]')
#' Cmodel <- compileNimble(model)
#' Cmy_LWF <- compileNimble(my_LWF, project = model)
#' logLike <- Cmy_LWF(m = 100000)
#' lw_X <- as.matrix(Cmy_LWF$mvEWSamp, 'x')
#' 
#' #samples from posterior of a top level parameter named sigma:
#' lw_sigma <- as.matrix(Cmy_LWF$mvEWSamp, 'sigma')
#' 
#' 
#' @export
buildLWF <- nimbleFunction(
  setup = function(model, nodes, control = list()){
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    
    
    
    
    saveAll <- control[['saveAll']]
    silent <- control[['silent']]
    params <- control[['params']]
    d <- control[['d']]
    if(is.null(silent)) silent <- FALSE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(d)) d <- .99
    # if unspecified, parameter nodes are specified as all stochastic top level nodes which
    # are not in the set of latent nodes above
    if(is.null(params)|is.na(params)){
      params <- setdiff(model$getNodeNames(stochOnly=T, includeData=F,
                                           topOnly=T),nodes)
    }
    params <- model$expandNodeNames(params, sort = TRUE)
    
    if(identical(params, character(0)))
      stop('must be at least one higher level parameter for Liu and West filter to work')
    if(any(params %in% nodes))
      stop('parameters cannot be latent states')
    if(!all(params%in%model$getNodeNames(stochOnly=T)))
      stop('parameters must be stochastic nodes')
    
    latentVars <- model$getVarNames(nodes = nodes)
    paramVars <-  model$getVarNames(nodes =  params)  # need var names too
    
    
    if(!saveAll) smoothing <- FALSE
    
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent 
                                      states varies')
    paramDims <-   sapply(params, function(n) nimDim(model[[n]]))
    
    
    modelSymbolObjects = model$getSymbolTable()$getSymbolObjects()[c(latentVars, paramVars)]
    if(saveAll){
      
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x){
        if(identical(x$size, numeric(0))) return(1)
        return(x$size)})      
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
      size <- lapply(modelSymbolObjects, function(x){
        if(identical(x$size, numeric(0))) return(1)
        return(x$size)})
      
      
      size[[latentVars]] <- as.numeric(dims[[1]])
      
      
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
    varSymbolObjects = model$getSymbolTable()$getSymbolObjects()[paramVars]
    paramVarDims <- unlist(lapply(varSymbolObjects, function(x){
      if(identical(x$size, numeric(0))) return(1)
      return(x$size)}))
    #names <- model$getSymbolTable()$getSymbolObjects()[[latentVars]]$name
    LWStepFunctions <- nimbleFunctionList(LWStepVirtual)
    for(iNode in seq_along(nodes))
      LWStepFunctions[[iNode]] <- LWStep(model,  mvWSamp, mvEWSamp, nodes, paramVarDims, 
                                         iNode, params, paramVars, names, saveAll, d, silent)
    
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    my_initializeModel$run()
    resize(mvWSamp, m)
    resize(mvEWSamp, m)
    
    for(i in 1:m)
      mvWSamp['wts',i][1] <<- log(1/m)
    logL <- 0
    for(iNode in seq_along(LWStepFunctions)) {
      logL <- logL + LWStepFunctions[[iNode]]$run(m)
      if(logL == -Inf) return(logL) }
    return(logL)
  },  where = getLoadingNamespace()
)


