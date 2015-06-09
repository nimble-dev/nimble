##  Contains code to run a Liu and West filter.  The filter has a build
##  function (buildLWF) and a step function (LWstepNS). Also contains 
##  a function for caluclating useful quantities.  
##
##  Still to implement: LWstepS, better parameter name checking, 
##  option to skip auxiliary step in filter (Ionides '11)
##  getpars() R function which returns nicely formatted matrix with 
##  parameter names (rows) and posteriors.



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

LWStepNS <- nimbleFunction(
  contains = LWStepVirtual,
  setup = function(model, mv, nodes, parDim, iNode, paramNodes, d, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    parDeterm <- model$getDependencies(paramNodes, determOnly=TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    parCalc <- LWparFunc(mv, d, parDim) #d is discount factor, default to 0.99
    getmean <- normMean(model, thisNode)
  },
  run = function(m = integer()) {
    returnType(double())
    declare(auxWts, double(1,m))
    declare(l, double(1,m))
    declare(wts, double(1, m))
    declare(ids, integer(1, m))
    declare(preWts, double(1,m))
    declare(tmpPars, double(2, c(parDim, m)))
    ess <- 0
    if(notFirst){
      meanVec <- parCalc$shrinkMean(m)  # shrink parameter particles towards mean  
      for(i in 1:m) {
        copy(mv, model, 'xs', prevNode, row=i)
        values(model, paramNodes) <<- meanVec[,i]
        calculate(model, parDeterm) 
        calculate(model, prevDeterm) 
        model[[thisNode]] <<- getmean$return_mean()
        calculate(model, thisDeterm)
        auxWts[i] <- exp(calculate(model, thisData))
        if(is.nan(auxWts[i])) auxWts[i] <- 0 #check for ok param values
        preWts[i] <- auxWts[i]*mv['wts',i][1] #first resample weights
      }
      rankSample(preWts, m, ids, silent)
      # Reassign weights and pars so that cov matrix is calculated correctly
      for(i in 1:m){
        wts[i] <- mv['wts',ids[i]][1] 
        tmpPars[,i] <- mv['pars', ids[i]] 
      }
      parChol <- parCalc$cholesVar(m,wts,tmpPars) #  calculate MC cov matrix
    }

    for(i in 1:m) {
      if(notFirst) {  
        tmpPars[,i] <- rmnorm_chol(1, meanVec[,ids[i]],parChol, prec_param=0) 
        values(model, paramNodes) <<- tmpPars[,i]
        copy(mv, model, nodes = 'xs', nodesTo = prevNode, row = ids[i]) 
      }
      else{
        simulate(model, paramNodes)
        mv['pars',i] <<- values(model,paramNodes)
      }
      calculate(model, parDeterm) 
      simulate(model, thisNode)
      copy(model, mv, nodes = thisNode, nodesTo = 'x', row = i)
      calculate(model, thisDeterm)
      l[i]  <- exp(calculate(model, thisData))
      if(is.nan(l[i])) l[i] <- 0
      ess <- ess+(l[i])^2
      #rescale weights by pre-sampling weight
      if(notFirst)
        wts[i] <- l[i]/auxWts[ids[i]]
      else 
        wts[i] <- l[i]
      if(is.nan(wts[i])) wts[i] <- 0
      
    }
    #normalizing weights and calculating effective sample size 
    ess <- 1/(ess/(sum(l)^2))
  
    for(i in 1:m){
      copy(mv, mv, "x", "xs", i, i)  
      if(notFirst) mv['pars',i] <<- tmpPars[,i]
      mv['wts',i][1] <<- wts[i] 
    }
    return(log(mean(l)))
  },  where = getLoadingNamespace()
)

# Has two methods: shrinkMean, which shrinks each parameter particle towards
# the mean of all particles, and cholesVar, which returns the cholesky 
# decomposition of the weighted MC covariance matrix
LWparFunc <- nimbleFunction(
  setup = function(mv,d, parDim){
    # Calculate h^2 and a using specified discount factor d
    hsq <- 1-((3*d-1)/(2*d))^2 
    a <- sqrt(1-hsq)
  },  
  methods = list(                        
    shrinkMean = function(m = integer()) {
      declare(parMean, double(2))
      declare(wtMean,double(2, c(parDim, m)))
      declare(pars, double(2, c(parDim, m)))
      declare(wts, double(1, m))
      declare(oneMat, double(1, m))
      for(i in 1:m){
        pars[,i] <- mv['pars',i]  
        wts[i] <- mv['wts',i][1]
        oneMat[i] <- 1}
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
#' @param params A character vector sepcifying the parameters you would 
#'  like to estimate.
#' @param d  A discount factor for the LW filter.  Should be close to,
#'  but not above, 1.
#' @author Nick Michaud
#' @family filtering methods
#' @details The resulting specialized particle filter algorthm will accept a
#'  single integer argument (m, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm 
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in mv['x',], with corresponding unlogged weights in mv['wts',].
#'  An equally weighted sample from the posterior can be found in mv['xs',]. 
#'  Samples from the posterior distribution of the specified parameters are in
#'  mv['params',].
#' @examples
#' model <- nimbleModel(code = ...)
#' my_LWF <- buildPF(model, 'x[1:100]')
#' Cmodel <- compileNimble(model)
#' Cmy_LWF <- compileNimble(my_LWF, project = model)
#' logLike <- Cmy_LWF(m = 100000)
#' lw_X <- Cmy_LWF$mv['xs',]
#' lw_pars <- Cmy_LWF$mv['pars',]
#' @export
buildLWF <- nimbleFunction(
  setup = function(model, nodes, params=NA, d = .99, silent = FALSE) {
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    # if unspecified, parameter nodes are specified as all stochastic top level nodes which
    # are not in the set of latent nodes above
    
    if(is.na(params)){
      params <- setdiff(model$getNodeNames(stochOnly=T, includeData=F,
                                           topOnly=T),nodes)
    }
    else{
      if(any(params %in% nodes))
        stop('parameters cannot be latent states')
      if(!(params%in%model$getNodeNames(stochOnly=T)))
        stop('parameters must be stochastic nodes')
    }
    dims <- lapply(nodes, function(n) nimbleDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent 
                                      states varies')
    parDim <- length(params)
    
    #  Create mv variables for x state, sampled x states, params, and weights.
      mv <- modelValues(modelValuesSpec(vars = c('x', 'xs', 'pars', 'wts'),  
                                        type = c('double', 'double', 'double',
                                                 'double'),
                                        size = list(x = dims[[1]],
                                                    xs = dims[[1]],
                                                    pars = parDim,
                                                    wts = 1)))
    LWStepFunctions <- nimbleFunctionList(LWStepVirtual)
    for(iNode in seq_along(nodes))
       LWStepFunctions[[iNode]] <- LWStepNS(model, mv, nodes, parDim, 
                                            iNode, params, d, silent)
    
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    my_initializeModel$run()
    resize(mv, m)
    for(i in 1:m)
      mv['wts',i][1] <<- 1/m
    logL <- 0
    for(iNode in seq_along(LWStepFunctions)) {
      logL <- logL + LWStepFunctions[[iNode]]$run(m)
      if(logL == -Inf) return(logL) }
    return(logL)
  },  where = getLoadingNamespace()
)
