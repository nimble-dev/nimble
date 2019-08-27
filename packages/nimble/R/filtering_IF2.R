##  Contains code to run a Bayes map iterated filter.  The filter has a build
##  function (buildIteratedFilter2) and a step function (IF2stepNS). Also contains 
##  a function for calculating useful quantities. The doPars2() function
##  is specialized to different variables in the mv objects, and allows
##  these variables to be set / retrieved within the mv objects without
##  explicitly calling the variable name.
##  
##
##  Implementing follows Ionides et. all 2015 and adapted Liu and West filter of Nicholas Michaud

initializeParamSwarm <- nimbleFunction(
    setup = function(model, mvEWSamples, paramNodes, numParams, sigma) {
    },
    run = function(m = integer()) {
        returnType(double())
        initialValues <- values(model, paramNodes)
        currentValues <- initialValues
        for(i in 1:m) {
            for(j in 1:numParams)
                currentValues[j] <- rnorm(1, initialValues[j], sigma[j])
            values(model, paramNodes) <<- currentValues
            nimCopy(model, mvEWSamples, nodes = paramNodes, row = i)
        }
        return(0)
    }
)


IF2Step0 <- nimbleFunction(
    setup = function(model, mvEWSamples, baselineNode, latentVar,
                     paramNodes, numParams, sigma, timeLength, silent = FALSE) {

        thisNodeExpanded <- model$expandNodeNames(baselineNode, sort = TRUE)
        ## Code simplified from particleFilter_splitModelSteps.
        thisDeterm <- model$getDependencies(baselineNode, determOnly = TRUE)
        if(length(thisDeterm) > 0) {
            thisDeterm_is_intermediate <- logical(length(thisDeterm))
            for(i in seq_along(thisDeterm)) {
                theseDeps <- model$getDependencies(thisDeterm[i], stochOnly = TRUE)
                thisDeterm_is_intermediate[i] <- any(theseDeps %in% thisNodeExpanded)
            }
            thisDeterm_self <- thisDeterm[ thisDeterm_is_intermediate ]
            thisDeterm <- thisDeterm[ !thisDeterm_is_intermediate ]
            calc_thisNode_self <-  model$expandNodeNames(c(thisNodeExpanded, thisDeterm_self), ## only for the sort
                                                         sort = TRUE)
        } else {
            calc_thisNode_self <- thisNodeExpanded
        }
        ## Probably not necessary, as thisDeterm should already be expanded and sorted.
        calc_thisNode_deps <- model$expandNodeNames(thisDeterm, sort = TRUE) 
        
        parDeterm <- model$getDependencies(paramNodes, determOnly=TRUE)

    },
    run = function(m = integer(), j = integer(), coolingRate = double()) {
        returnType(double())
        l <- numeric(m, init=FALSE)
        ## use same sigma as for t=1
        coolParam <- (coolingRate)^(((j - 1)*timeLength)/(50*timeLength))
        coolSigma <- coolParam*sigma
        
        for(i in 1:m) {
            nimCopy(mvEWSamples, model, nodes = paramNodes, row = i)
            currentValues <- values(model, paramNodes)
            for(j in 1:numParams)
                currentValues[j] <- rnorm(1, currentValues[j], coolSigma[j])
            values(model, paramNodes) <<- currentValues
            calculate(model, parDeterm)
            simulate(model, calc_thisNode_self)
            calculate(model, calc_thisNode_deps)
            nimCopy(model, mvEWSamples, nodes = baselineNode, nodesTo = latentVar, rowTo = i)
            nimCopy(model, mvEWSamples, nodes = paramNodes, rowTo = i)
        }
        return(0)
    },  where = getLoadingNamespace()
)

IF2StepVirtual <- nimbleFunctionVirtual(
    run = function(m = integer(), n = integer(), coolingRate = double())
        returnType(double())
)

IF2Step <- nimbleFunction(
    contains = IF2StepVirtual,
    setup = function(model, mvWSamples, mvEWSamples, latentNodes, latentVar, baselineNode, baseline,
                     iNode, paramNodes, numParams, sigma, timeLength, silent = FALSE) {
        notFirst <- iNode != 1
        isSecond <- iNode == 2
        prevNode <- latentNodes[if(notFirst) iNode-1 else iNode]

        modelSteps <- particleFilter_splitModelSteps(model, latentNodes, iNode, notFirst)
        prevDeterm <- modelSteps$prevDeterm
        calc_thisNode_self <- modelSteps$calc_thisNode_self
        calc_thisNode_deps <- modelSteps$calc_thisNode_deps
        
        thisNode <- latentNodes[iNode]
        parDeterm <- model$getDependencies(paramNodes, determOnly=TRUE)
        parAndPrevDeterm <- c(parDeterm, prevDeterm)

        isLast <- (iNode == timeLength)
        coolSigma <- numeric(numParams)
        coolParam <- 0
    },
    run = function(m = integer(), j = integer(), coolingRate = double()) {
        returnType(double())
        l <- numeric(m, init=FALSE)
        wts <- numeric(m, init=FALSE)
        ids <- integer(m, 0)
        coolParam <<- (coolingRate)^((iNode-1+(j - 1)*timeLength)/(50*timeLength))
        coolSigma <<- coolParam*sigma
        
        for(i in 1:m) {
            ids[i] <- i  ## use of 1:m causes type mismatch
            nimCopy(mvEWSamples, model, nodes = paramNodes, row = i)
            if(notFirst) {
                copy(mvEWSamples, model, nodes = latentVar, nodesTo = prevNode, row = i)
            } else {
                if(baseline)
                    copy(mvEWSamples, model, nodes = latentVar, nodesTo = baselineNode, row = i)
            }
            currentValues <- values(model, paramNodes)
            for(j in 1:numParams)
                currentValues[j] <- rnorm(1, currentValues[j], coolSigma[j])
            values(model, paramNodes) <<- currentValues
            calculate(model, parAndPrevDeterm)
            simulate(model, calc_thisNode_self)
            logProb <- calculate(model, calc_thisNode_deps)
            wts[i]  <- exp(logProb)
            if(is.nan(wts[i])) wts[i] <- 0
            logProb <- calculate(model, paramNodes)
            if(is.na(logProb) | logProb == -Inf)
                wts[i] <- 0
            nimCopy(model, mvWSamples, nodes = thisNode, nodesTo = latentVar, rowTo = i)
            nimCopy(model, mvWSamples, nodes = paramNodes, rowTo = i)
            mvWSamples['wts', i][1] <<- wts[i]
        }
        lik <- mean(wts)
        wts <- wts/(m*lik) ## = wts / sum(wts)
        rankSample(wts, m, ids, silent)
        for(i in 1:m){
            copy(mvWSamples, mvEWSamples, nodes = latentVar, row = ids[i], rowTo = i)
            copy(mvWSamples, mvEWSamples, nodes = paramNodes, row = ids[i], rowTo = i)
        }
        return(lik)
    },  where = getLoadingNamespace()
)

#' Create an IF2 algorithm.  
#' 
#' @description Create an IF2 algorithm for a given NIMBLE state space model.  
#'
#' @param model A NIMBLE model object, typically representing a state 
#'  space model or a hidden Markov model.
#' @param nodes A character vector specifying the latent model nodes 
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.
#' @param params A character vector specifying the top-level parameters to obtain maximum likelihood estimates of. 
#'   If unspecified, parameter nodes are specified as all stochastic top level nodes which
#'  are not in the set of latent nodes specified in \code{nodes}.
#' @param baselineNode A character vector specifying the node that is the latent node at the "0th" time step. The first node in \code{nodes} should depend on this baseline, but it should have no data depending on it. If \code{NULL} (the default), any initial state is taken to be fixed at the values present in the model at the time the algorithm is run.
#' @param control  A list specifying different control options for the IF2 algorithm.  Options are described in the \sQuote{details} section below.

#' @author Nicholas Michaud, Dao Nguyen, and Christopher Paciorek
#' @family particle filtering methods
#' @details 
#' 
#' Each of the \code{control()} list options are described in detail below:
#' \describe{
#'  \item{sigma}{A vector specifying a non-negative perturbation magnitude for each element of the \code{params} argument.  Defaults to a vector of 1's.}
#'  \item{initParamSigma}{A vector specifying a vector of standard deviations to use when simulating an initial particle swarm centered on the initial value of the parameters.}
#'  \item{inits}{A vector specifying an initial value for each element of the \code{params} argument.  Defaults to the parameter values in the model at the time the model is built.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.  
#'  Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to \code{TRUE}.}
#' }
#' 
#'  The IF2 agorithm uses iterated filtering to estimate maximum likelihood values for top-level parameters for a state space model.  
#'        
#'  The resulting specialized IF2 algorithm will accept the following arguments:
#'  
#'  \describe{
#'    \item{m}{A single integer specifying the number of particles to use for each run of the filter.  Defaults to 10,000.}
#'    \item{n}{A single integer specifying the number of overall filter iterations to run.  Defaults to 5.}
#'    \item{coolingRate}{A double specifying the cooling rate to use for the IF2 algorithm.  Defaults to 0.2.}
#'
#'  The \code{run} fuction will return a vector with MLE estimates.  Additionally, once the specialized algorithm has been run, it can be continued for additional iterations by calling the \code{continueRun} method.
#'  
#' @export
#' 
#' @references Ionides, E.L., D. Nguyen, Y. Atchad{\'e}, S. Stoev, and A.A. King (2-15). Inference for dynamic and latent variable models via iterated, perturbed Bayes maps. \emph{Proceedings of the National Academy of Sciences}, 112(3), 719-724.
#' 
#' @examples
#' \dontrun{
#' model <- nimbleModel(code = ...)
#' my_IF2 <- buildIteratedFilter2(model, 'x[1:100]', params = 'sigma_x')
#' Cmodel <- compileNimble(model)
#' Cmy_IF2 <- compileNimble(my_IF2, project = model)
#' # MLE estimate of a top level parameter named sigma_x:
#' sigma_x_MLE <- Cmy_IF2$run(m = 10000, n = 10)
#' # Continue running algorithm for more precise estimate:
#' sigma_x_MLE <- Cmy_IF2$continueRun(n = 10) 
#' }
buildIteratedFilter2 <- nimbleFunction(
    setup = function(model, nodes, params = NULL, baselineNode = NULL, control = list()){
        
        ## control list extraction
        silent <- control[['silent']]
        inits <- control[['inits']]
        sigma <- control[['sigma']]
        ## Still need to allow user to provide initial param particles
        initParamSigma <- control[['initParamSigma']]  
        timeIndex <- control[['timeIndex']]
        initModel <- control[['initModel']]
        if(is.null(silent)) silent <- TRUE
        if(is.null(initModel)) initModel <- TRUE

        if(!is.null(baselineNode))
            if(length(model$getDependencies(baselineNode, dataOnly = TRUE)))
                stop("buildIteratedFilter2: 'baselineNode' should not have any data nodes as dependents.")

        ## if unspecified, parameter nodes are specified as all stochastic top level nodes which
        ## are not in the set of latent nodes above
        if(all(is.null(params))){
            params <-  model$getNodeNames(stochOnly = TRUE, includeData = FALSE,
                                          topOnly = TRUE)
            params <- params[!params %in% nodes]
        }
        params <- model$expandNodeNames(params)
        numParams <- length(params)
        parDeterm <- model$getDependencies(params, determOnly=TRUE)

        if(identical(params, character(0)))
            stop('buildIteratedFilter2: There must be at least one parameter for IF2 to optimize with respect to.')
        if(any(params %in% nodes))
            stop('buildIteratedFilter2: Parameters cannot be latent states.')
        if(!all(params %in% model$getNodeNames(stochOnly = TRUE)))
            stop('buildIteratedFilter2: Parameters must be stochastic nodes.')
        
        paramVars <-  model$getVarNames(nodes = params)

        if(is.null(sigma)){
            sigma <- rep(1, numParams)
        }
        if(length(sigma) != numParams)
            stop("buildIteratedFilter2: The 'sigma' control list argument must be a vector specifying a
            non-negative perturbation magnitude for each element of 'params'. The length of 'sigma'
            does not match the length of 'params'.")
        
        if(any(sigma < 0))
            stop("buildIteratedFilter2: All values of 'sigma' should be non-negative.")
        
        ## Check for inits values
        if(is.null(inits)){
            inits <- values(model, params)
        }
        if(length(inits) < numParams)
            stop("buildIteratedFilter2: The 'inits' control list argument must be a vector specifying initial values for each element of 'params'. The length of 'inits' does not match the number of parameters.")
        
        latentVar <- model$getVarNames(nodes = nodes) 
        if(length(unique(latentVar)) > 1){
            stop("buildIteratedFilter2: All latent nodes must be in the same variable.")
        }

        info <- model$getVarInfo(latentVar)
        latentDims <- info$nDim
        if(is.null(timeIndex)){
            timeIndex <- which.max(info$maxs)
            timeLength <- max(info$maxs)
            if(sum(info$maxs == timeLength) > 1) # check if multiple dimensions share the max index size
                stop("buildIteratedFilter2: Unable to determine which dimension indexes time. Specify manually using the 'timeIndex' control list argument.")
        } else{
            timeLength <- info$maxs[timeIndex]
        }
        ## CJP note: this assumes nodes in variable are in time order.
        ## Should we tell the user this is our assumption?  
        nodes <- paste(info$varName, "[", rep(",", timeIndex-1), 1:timeLength,
                       rep(",", info$nDim - timeIndex), "]", sep="")
        
        dims <- lapply(nodes, function(node) nimDim(model[[node]]))
        if(length(unique(dims)) > 1) stop('sizes or dimension of latent 
                                      states varies')

        my_initializeModel <- initializeModel(model, silent = silent)

        modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[c(latentVar, paramVars)]
        names <- sapply(modelSymbolObjects, function(x) return(x$name))
        types <- sapply(modelSymbolObjects, function(x) return(x$type))
        sizes <- lapply(modelSymbolObjects, function(x) {
            if(identical(x$size, numeric(0))) return(1)
            return(x$size)})
        sizes[[latentVar]] <- as.numeric(dims[[1]])      

        mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                   types = types,
                                                   sizes = sizes))
        names <- c(names, "wts")
        types <- c(types, "double")
        sizes$wts <- 1
        mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                   types = types,
                                                   sizes = sizes))

        initializeParamSwarmFunction <- initializeParamSwarm(model, mvEWSamples, params,
                                                       numParams, initParamSigma)
        if(is.null(baselineNode)) {
            ## use params[1] simply as a dummy as IF2Step0Function needs to exist
            baselineNode <- params[1]
            IF2Step0Function <- IF2Step0(model,  mvEWSamples, params[1], latentVar,
                                          params, numParams, sigma, timeLength, silent)
            baseline <- FALSE
        } else {
            IF2Step0Function <- IF2Step0(model,  mvEWSamples, baselineNode, latentVar,
                                         params, numParams, sigma, timeLength, silent)
            baseline <- TRUE
        }
        
        IF2StepFunctions <- nimbleFunctionList(IF2StepVirtual)
        for(iNode in seq_along(nodes))
            IF2StepFunctions[[iNode]] <- IF2Step(model,  mvWSamples, mvEWSamples, nodes, latentVar, baselineNode,
                                                 baseline, iNode, params, numParams, sigma, timeLength, silent)

        estimate <- nimNumeric(numParams)  
        
        oldJ <- 0
        oldM <- 0
        logLik <- nimNumeric(2)
        estimates <- nimMatrix(0, 1, numParams)
        estSD <- nimMatrix(0, 1, numParams)
    },
    run = function(m = integer(default = 10000), n = integer(default = 5), 
                   coolingRate = double(default = 0.2)) {
        
        values(model, params) <<- inits
        my_initializeModel$run()
        resize(mvWSamples, m)
        resize(mvEWSamples, m)
        estSD <<- nimMatrix(0, n, numParams)
        estimates <<- nimMatrix(0, n, numParams)
        logLik <<- nimNumeric(n, value = 0)

        initializeParamSwarmFunction$run(m)
        for(j in 1:n){
            ## Initialize latent process and params at time 0.
            if(baseline)
                IF2Step0Function$run(m, j, coolingRate)
            for(iNode in seq_along(IF2StepFunctions)) {
                logLik[j] <<- logLik[j] + log(IF2StepFunctions[[iNode]]$run(m, j, coolingRate))
            }
            
            ## Compute estimate and sd of particles at each iteration for diagnostics.
            ## Revisit possibility of weighted sample (when changed to that, got strange results.)
            for(i in 1:m) {
                nimCopy(mvEWSamples, model, nodes = params, row = i)
                estimates[j, ] <<- estimates[j, ] + values(model, params) 
            }
            estimates[j, ] <<- estimates[j, ] / m
            for(i in 1:m) {
                nimCopy(mvEWSamples, model, nodes = params, row = i)
                estSD[j, ] <<- estSD[j, ] + (values(model, params) - estimates[j,])^2
            }
            estSD[j, ] <<- sqrt(estSD[j, ] / m)
        }

        estimate <<- estimates[n, ]

        ## Put final parameter estimates into model.
        values(model, params) <<- estimate
        model$calculate(parDeterm)
        oldM <<- m
        oldJ <<- j 
        
        returnType(double(1))
        return(estimate)
    },
    methods = list(
        continueRun = function(n = integer(default = 5), coolingRate = double(default = 0.2)){
            ## add storage of logLik, estimates and particle SDs as in main $run
            useStoredSamples <- 1
            newN <- oldJ + n
            for(j in (oldJ+1):newN){
                if(baseline) 
                    IF2Step0Function$run(oldM, j, coolingRate)
                for(iNode in seq_along(IF2StepFunctions)) {
                    IF2StepFunctions[[iNode]]$run(oldM, j, coolingRate)
                }
            }

            ## Calculate mean as final estimate.
            estimate <<- rep(0, numParams)
            for(i in 1:oldM) {
                nimCopy(mvWSamples, model, nodes = params, row = i)
                estimate <<- estimate + values(model, params)
            }
            estimate <<- estimate / oldM

            ## Leave model with parameter estimates.  
            values(model, params) <<- estimate
            model$calculate(parDeterm)

            oldJ <<- newN

            returnType(double(1))
            return(estimate)
        }
    ),
    where = getLoadingNamespace()
)
