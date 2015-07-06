##  Contains code to run Ensemble Kalman Filter.
##  Algorithm found in Evenson '03.
##  The buildENKF() funtion builds and runs the ENKF.
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
  contains = ENKFFuncVirtual,
  setup = function(model, thisData, mv, thisXName){
    data <- model[[thisData]]
    dataFuncList <- nimbleFunctionList(node_stoch_dmnorm)
    dataFuncList[[1]] <- model$nodeFunctions[[thisData]]
  },
  methods = list(
    getMean = function(){
      returnType(double(1))
      return(dataFuncList[[1]]$get_mean())
    },
    getVar = function(){
      returnType(double(2))
      return(dataFuncList[[1]]$get_cov())
    }
  ), where = getLoadingNamespace()
)

#  returns mean (as vector) and var (as matrix) for normal data node  
enkfScalFunc = nimbleFunction(
  contains = ENKFFuncVirtual,
  setup = function(model, thisData, mv, thisXName){
    data <- model[[thisData]]
    output <- c(0,0)
    dataFuncList <- nimbleFunctionList(node_stoch_dnorm) 
    dataFuncList[[1]] <- model$nodeFunctions[[thisData]]
  },
  methods = list(
    getMean = function(){
      returnType(double(1))
      setSize(output, 1)
      output[1] <<- dataFuncList[[1]]$get_mean()
      return(output)
    },
    getVar = function(){
      returnType(double(2))
      declare(outMat, double(2, c(1,1)))
      outMat[1,1] <- dataFuncList[[1]]$get_var()
      return(outMat)
    }
  ), where = getLoadingNamespace()
)


# Ensemble Kalman filter step
# Currently assumes that model is of form specified in Gillijns et al. '06
# Does not check to verify this.

ENKFStep <- nimbleFunction(
  contains = ENKFStepVirtual,
  setup = function(model, mv, nodes, iNode, xDim, yDim, saveAll, silent = FALSE) {
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
      prevXName <- paste("x[,",t-1,"]",sep="")
      thisXName <- paste("x[,",t,"]",sep="")
    }
    else{
      prevXName <- "x"
      thisXName <- "x"
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
        ENKFFuncList[[yNode]] <- enkfMultFunc(model, thisData[yNode],  mv, thisXName)
      }
      else{
        ENKFFuncList[[yNode]] <- enkfScalFunc(model, thisData[yNode],  mv, thisXName)
      } 
    }
  },
  run = function(m = integer()) {
    declare(xf, double(2, c(xDim, m)))  # currently only works for x of length 1!
    declare(yf, double(2, c(yLength, m)))
    declare(varMat, double(2, c(yLength, yLength))) # combined covariance matrix for all depndent nodes
    declare(preturb, double(2, c(yLength, m)))  # matrix for preturbed observations  
    declare(md, double())
    md <- m  # make m a double for calculations
    if(yLength == 1){
      setSize(yObs, 1)
      setSize(meanVec, 1)
    }
    #  Forecast Step
    #  cycle through data nodes and particles,
    #  forecast x values, get mean and variance of y, add to mean vector and cov matrix
    for(j in 1:nNodes){
      yVar <- ENKFFuncList[[j]]$getVar() # var doesn't depend on x value
      for(i in 1:m) {
        if(notFirst) {
          copy(mv, model, nodes = prevXName, nodesTo = prevNode, row = i)
          calculate(model, prevDeterm) 
        }
        simulate(model, thisNode)
        calculate(model, thisDeterm) 
        xf[,i] <- values(model,thisNode)
        yfVal <- ENKFFuncList[[j]]$getMean()
        ySize <- yInds[j+1]-yInds[j]
        if(ySize == 1){
          yf[yInds[j+1], i] <-yfVal[1]
          varMat[yInds[j+1], yInds[j+1]] <-yVar[1,1]
        }
        else{
          yf[(yInds[j]+1):yInds[j+1], i] <- yfVal
          varMat[(yInds[j]+1):yInds[j+1], (yInds[j]+1):yInds[j+1]] <- yVar
        }
      }
    }
    

    #  Analysis step
    #  first calculate approx Kalman gain matrix (pg. 350)
    oneVec <- nimVector(1,m)
    efx <- xf - (1/md)*(xf%*%oneVec%*%t(oneVec))
    efy <- yf -  (1/md)*(yf%*%oneVec%*%t(oneVec))
    kMat <- (1/(md-1))*efx%*%t(efy)%*%inverse((1/(md-1))*(efy%*%t(efy)+varMat))

    
    #  next, cycle through particles, create preturbed observations,
    #  and store in model values object
    if(yLength == 1){
      for(i in 1:m){
        preturb[1,i] <-yObs[1] + rnorm(1, 0, sqrt(varMat[1,1]))
        mv[thisXName, i] <<-xf[,i] + kMat%*%(preturb[,i] -yf[,i]) 
      }
    }
    else{
      cholesky <- chol(varMat)    
      for(i in 1:m){
        preturb[,i] <- yObs + rmnorm_chol(1, meanVec, cholesky, 0) 
        mv[thisXName, i] <<-xf[,i] + kMat%*%(preturb[,i] -yf[,i]) 
      }
    }
  }, where = getLoadingNamespace()
)

#' Creates an Ensemble Kalman filter algorithm to sample from latent states.
#'
#' @param model A nimble model object, typically representing a state space model or a hidden Markov model
#' @param nodes A character vector specifying the latent model nodes which the ENKF will estimate.
#' @param saveAll  Whether to save state samples for all time points (T), or only for the most recent time points (F)
#' @author Nick Michaud
#' @details Runs an Ensemble Kalman filter to estimate a latent state given observations at each time point.  
#' Latent state (x[t]) and Observations (y[t]) can be scalars or vectors at each time point, 
#' and sizes of observations can vary from time point to time point.
#' In the BUGS model, the observations (y[t]) must be equal to some (possibly nonlinear) deterministic function
#' of the latent state (x[t]) plus an additive error term.  Currently only normal/mvn error terms are supported.
#' The transition from x(t) to x(t+1) does not have to be normal or linear.
#' @family filtering methods
#' @examples
#' model <- nimbleModel(code = ...)
#' my_ENKFF <- buildENKF(model, 'x')
#' Cmodel <- compileNimble(model)
#' Cmy_ENKF <- compileNimble(my_ENKF, project = model)
#' Cmy_ENKF$run(m = 100000)
#' ENKF_X <- as.matrix(Cmy_ENKF$mv, 'x')
#' hist(ENKF_X)
#' @export
buildENKF <- nimbleFunction(
  setup = function(model, nodes,  silent = FALSE, saveAll = FALSE) {
    my_initializeModel <- initializeModel(model)
    nodes <- model$expandNodeNames(nodes, sort = TRUE)
    dims <- lapply(nodes, function(n) nimbleDim(model[[n]]))
    if(length(unique(dims)) > 1) stop('sizes or dimension of latent states varies')
    xDim <- dims[[1]]
    #  get list of y nodes depening on each x node
    #  necessary if the bugs model specifies something like:
    #  y[1,1] ~ dnorm(x[1], 1)
    #  y[1,2] ~ dnorm(5*pow(x[1],2),1)
    #  where multiple separately specified y nodes depend on the same x node(s)
    yNodes <- lapply(nodes, function(n) model$getDependencies(n, dataOnly = TRUE))  
    yDim <- lapply(yNodes, function(x){unname(sapply(x, function(z) nimbleDim(model[[z]])))}) #dimensions of each dependent y node
    
    # Create mv variables for x state.  If saveAll=T, 
    # the  x states will be recorded at each time point. 
    if(!saveAll){
      mv <- modelValues(modelValuesSpec(vars = c('x'),  
                                        type = c('double'),
                                        size = list(x = xDim)))
    }
    
    else{
      mv <- modelValues(modelValuesSpec(vars = c('x'),
                                        type = c('double'),
                                        size = list(x = c(xDim,
                                                          length(dims)))))
    }
    ENKFStepFunctions <- nimbleFunctionList(ENKFStepVirtual)
    for(iNode in seq_along(nodes)){
     ENKFStepFunctions[[iNode]] <- ENKFStep(model, mv, nodes, 
                                               iNode, xDim, yDim[[iNode]], saveAll,  silent) 
    }
  },
  run = function(m = integer(default = 100)) {
    my_initializeModel$run()
    resize(mv, m) 
    for(iNode in seq_along(ENKFStepFunctions)) { 
      ENKFStepFunctions[[iNode]]$run(m)
    }
  },  where = getLoadingNamespace()
)
