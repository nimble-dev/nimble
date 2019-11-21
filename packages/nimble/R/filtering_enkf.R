##  Contains code to run Ensemble Kalman Filter.
##  Algorithm found in Evenson '03.
##  The buildEnsembleKF() funtion builds and runs the ENKF.
##  The ENKFStep() function gets samples from latent states
##  for one timepoint. The ENKFStep function uses either the 
##  ENKFMultFunc or ENFKScalFunc depending on whether the dependent data
##  is a scalar or vector at each time point.


ENKFStepVirtual <- nimbleFunctionVirtual(
  run = function(m = integer()){ 
    returnType()
  }
)

ENKFFuncVirtual <- nimbleFunctionVirtual(
  methods = list(
    getMean = function(){returnType(double(1))},
    getVar = function(){returnType(double(2))}
  )
)



#  returns mean and cov matrix for MVN data node
enkfMultFunc = nimbleFunction(
    name = 'enkfMultFunc',
  contains = ENKFFuncVirtual,
  setup = function(model, thisData){},
  methods = list(
    getMean = function(){
      returnType(double(1))
      return(model$getParam(thisData, 'mean'))
    },
    getVar = function(){
      returnType(double(2))
      return(model$getParam(thisData, 'cov'))
    }
  ), where = getLoadingNamespace()
)

#  returns mean (as vector) and var (as matrix) for normal data node  
enkfScalFunc = nimbleFunction(
    name = 'enkfScalFunc',
  contains = ENKFFuncVirtual,
  setup = function(model, thisData){},
  methods = list(
    getMean = function(){
      returnType(double(1))
      #  output is length 2 to ensure it is not recast as a scalar, only first element is used
      ##declare(output, double(1, 2))  
      output <- numeric(2)
      output[1] <- model$getParam(thisData, 'mean')
      return(output)
    },
    getVar = function(){
      returnType(double(2))
      ##declare(outMat, double(2, c(1,1)))
      outMat <- matrix(nrow = 1, ncol = 1, init=FALSE)
      outMat[1,1] <- model$getParam(thisData, 'var')
      return(outMat)
    }
  ), where = getLoadingNamespace()
)

# Ensemble Kalman filter step
# Currently assumes that model is of form specified in Gillijns et al. '06
# Does not check to verify this.

ENKFStep <- nimbleFunction(
    name = 'ENKFStep',
  contains = ENKFStepVirtual,
  setup = function(model, mvSamples, nodes, iNode, xDim, yDim, saveAll, names, silent = FALSE) {
    notFirst <- iNode != 1
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]
    thisNode <- nodes[iNode]
    prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    t <- iNode  # current time point
    yObs <- c(sapply(thisData, function(x) model[[x]])) #all dependent data for this time point
    # Get names of xs node for current and previous time point (used in copy)
    if(saveAll == 1){
      prevXSName <- prevNode    
      thisXSName <- thisNode
      currInd <- t
      prevInd <- t-1
    }
    else{
      prevXSName <- names    
      thisXSName <- names
      currInd <- 1
      prevInd <- 1
    }
    yLength <- sum(yDim)  #total number of dependent data points for this time point
    meanVec <- rep(0, yLength)
    
    if(yLength==1){
      yObs<-c(yObs,0)
      meanVec <- c(0,0)
    }

    ## indices for locations of data from different nodes within our yf vector
    ## e.g. if there were two y nodes depending on x[t], with the first being a vector of length 2, and the second being 
    ## a scalar , yInds would be (0, 2, 3)  - the first node would go in yf[1:2] and the second node in yf[3] 
    yInds <- c(0,cumsum(yDim)) 
    # 0 vector
    nNodes <- length(yDim)
    #determine whether each data node is multivariate or scalar and assign correct function
    ENKFFuncList <- nimbleFunctionList(ENKFFuncVirtual) 
    for(yNode in 1:length(yDim)){
      if(yDim[yNode] > 1){
        ENKFFuncList[[yNode]] <- enkfMultFunc(model, thisData[yNode])
      }
      else{
        ENKFFuncList[[yNode]] <- enkfScalFunc(model, thisData[yNode])
      } 
    }
  },
  run = function(m = integer()) {
    #declare(xf, double(2, c(xDim, m))) 
      xf <- matrix(nrow = xDim, ncol = m, init=FALSE)
      ##declare(yf, double(2, c(yLength, m)))
      yf <- matrix(nrow = yLength, ncol = m, init=FALSE)
      ##declare(varMat, double(2, c(yLength, yLength))) # combined covariance matrix for all depndent nodes
      varMat <- matrix(nrow = yLength, ncol = yLength, init=FALSE)
      ##declare(perturb, double(2, c(yLength, m)))  # matrix for perturbed observations  
      perturb <- matrix(nrow = yLength, ncol = m, init=FALSE)
      declare(md, double())
    md <- m  # make m a double for calculations
    if(yLength == 1){
      setSize(yObs, 1)
      setSize(meanVec, 1)
    }  
    
    
    # girst get var, doesn't depend on x value
    for(j in 1:nNodes){
      yVar <- ENKFFuncList[[j]]$getVar() 
      ySize <- yInds[j+1]-yInds[j]
      if(ySize == 1){
        varMat[yInds[j+1], yInds[j+1]] <-yVar[1,1]
      }
      else{
        varMat[(yInds[j]+1):yInds[j+1], (yInds[j]+1):yInds[j+1]] <- yVar
      }
    }
    
    
    #  Forecast Step
    #  cycle through data nodes and particles,
    #  forecast x values, get mean of y, add to mean vector
    for(i in 1:m) {
      if(notFirst) {
        copy(from = mvSamples, to = model, nodes = prevXSName, nodesTo = prevNode, row = i)          
        calculate(model, prevDeterm) 
      }
      simulate(model, thisNode)
      calculate(model, thisDeterm) 
      xf[,i] <- values(model,thisNode)
      for(j in 1:nNodes){
        yfVal <- ENKFFuncList[[j]]$getMean()
        ySize <- yInds[j+1]-yInds[j]
        if(ySize == 1){
          yf[yInds[j+1], i] <- yfVal[1]
        }
        else{
          yf[(yInds[j]+1):yInds[j+1], i] <- yfVal
        }
      }
    }
    
    #  Analysis step
    #  first calculate approx Kalman gain matrix (pg. 350)
    oneVec <- numeric(m,1)
    efx <- xf - (1/md)*(xf%*%oneVec%*%t(oneVec))
    efy <- yf -  (1/md)*(yf%*%oneVec%*%t(oneVec))
    kMat <- (1/(md-1))*efx%*%t(efy)%*%inverse((1/(md-1))*(efy%*%t(efy))+varMat)
    
    #  next, cycle through particles, create perturbed observations,
    #  and store in model values object
    if(yLength == 1){
      for(i in 1:m){
        perturb[1,i] <-yObs[1] + rnorm(1, 0, sqrt(varMat[1,1]))
        xFilterSamp <- xf[,i] + kMat%*%(perturb[,i] -yf[,i]) 
        values(model, thisNode) <<- xFilterSamp
        copy(model, mvSamples, thisNode, thisXSName, rowTo = i)
      }
    }
    else{
      cholesky <- chol(varMat)    
      for(i in 1:m){
        perturb[,i] <- yObs + rmnorm_chol(1, meanVec, cholesky, 0) 
        xFilterSamp <- xf[,i] + kMat%*%(perturb[,i] -yf[,i]) 
        values(model, thisNode) <<- xFilterSamp
        copy(model, mvSamples, thisNode, thisXSName, rowTo = i)      
      }
    }
  }, where = getLoadingNamespace()
)

#' Create an Ensemble Kalman filter algorithm to sample from latent states.
#' 
#' @description Create an Ensemble Kalman filter algorithm for a given NIMBLE state space model.  
#'
#' @param model A NIMBLE model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes 
#'  the Ensemble Kalman Filter will estimate. All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author Nicholas Michaud
#' @family particle filtering methods

#' @details 
#' The \code{control()} list option is described in detail below:
#' \describe{
#'  \item{saveAll}{Indicates whether to save state samples for all time points (TRUE), or only for the most recent time point (FALSE)}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.  
#'  Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#' }
#' 
#' Runs an Ensemble Kalman filter to estimate a latent state given observations at each time point.  The ensemble Kalman filter
#' is a Monte Carlo approximation to a Kalman filter that can be used when the model's transition euqations do not follow a normal distribution.  
#' Latent states (x[t]) and observations (y[t]) can be scalars or vectors at each time point, 
#' and sizes of observations can vary from time point to time point.
#' In the BUGS model, the observations (y[t]) must be equal to some (possibly nonlinear) deterministic function
#' of the latent state (x[t]) plus an additive error term.  Currently only normal and multivariate normal
#' error terms are supported.
#' The transition from x[t] to x[t+1] does not have to be normal or linear.  Output from the posterior distribution of the latent
#' states is stored in \code{mvSamples}.
#' @family filtering methods
#' @references Houtekamer, P.L., and H.L. Mitchell. (1998). Data assimilation using an ensemble Kalman filter technique. \emph{Monthly Weather Review}, 126(3), 796-811.
#' @examples
#' \dontrun{
#' model <- nimbleModel(code = ...)
#' my_ENKFF <- buildEnsembleKF(model, 'x')
#' Cmodel <- compileNimble(model)
#' Cmy_ENKF <- compileNimble(my_ENKF, project = model)
#' Cmy_ENKF$run(m = 100000)
#' ENKF_X <- as.matrix(Cmy_ENKF$mvSamples, 'x')
#' hist(ENKF_X)
#' }
#' @export
buildEnsembleKF <- nimbleFunction(
    name = 'buildEnsembleKF',
  setup = function(model, nodes, control = list()) {

    #control list extraction
    saveAll <- control[['saveAll']]
    silent <- control[['silent']]
    timeIndex <- control[['timeIndex']]
    initModel <- control[['initModel']]
    if(is.null(silent)) silent <- FALSE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(initModel)) initModel <- TRUE 
    ## get latent state info
    nodes <- findLatentNodes(model, nodes, timeIndex)
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent states varies')
    xDim <- dims[[1]]
     
    #  get list of y nodes depening on each x node
    #  necessary if the bugs model specifies something like:
    #  y[1,1] ~ dnorm(x[1], 1)
    #  y[1,2] ~ dnorm(5*pow(x[1],2),1)
    #  where multiple separately specified y nodes depend on the same x node(s)
    yNodes <- lapply(nodes, function(n) model$getDependencies(n, dataOnly = TRUE))  
    yDim <- unlist(sapply(yNodes, function(x){unname(sapply(x, function(z) nimDim(model[[z]])))})) #dimensions of each dependent y node
    
    errors <- sapply(model$expandNodeNames(yNodes), function(node){tryCatch(getParam(model, node, 'mean'), error=function(a){return("error")})})
    if(any(errors == "error", na.rm=TRUE)) stop("cannot use ENKF for this model, cannot obtain mean of data nodes ")
    
    
    # Create mv variables for x state.  If saveAll=TRUE, 
    # the  x states will be recorded at each time point.
    vars <- model$getVarNames(nodes =  nodes)  # need var names too
    modelSymbolObjects = model$getSymbolTable()$getSymbolObjects()[vars]
    if(saveAll){
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x){
        if(identical(x$size, numeric(0))) return(1)
        return(x$size)})
      
      mvSamples <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))
    }
    else{
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x){
        if(identical(x$size, numeric(0))) return(1)
        return(x$size)})
      
      size[[1]] <- as.numeric(dims[[1]])
      
      mvSamples <- modelValues(modelValuesConf(vars = names,
                                            types = type,
                                            sizes = size))
    }
    
    my_initializeModel <- initializeModel(model)
    names <- names[1]
    ENKFStepFunctions <- nimbleFunctionList(ENKFStepVirtual)
    for(iNode in seq_along(nodes)){
      ENKFStepFunctions[[iNode]] <- ENKFStep(model, mvSamples, nodes, 
                                             iNode, xDim, yDim[[iNode]], saveAll, names, silent) 
    }
  },
  run = function(m = integer(default = 100)) {
    if(initModel == TRUE) my_initializeModel$run()
    resize(mvSamples, m) 
    for(iNode in seq_along(ENKFStepFunctions)) { 
      ENKFStepFunctions[[iNode]]$run(m)
    }
  },  where = getLoadingNamespace()
)
