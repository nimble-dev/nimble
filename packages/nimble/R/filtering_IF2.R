##  Contains code to run a Bayes map iterated filter.  The filter has a build
##  function (buildIteratedFilter2) and a step function (IF2stepNS). Also contains 
##  a function for calculating useful quantities. The doPars2() function
##  is specialized to different variables in the mv objects, and allows
##  these variables to be set / retrieved within the mv objects without
##  explicitly calling the variable name.
##  
##
##  Implementation follows Ionides et. all 2015 and adapted Liu and West filter of Nicholas Michaud

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
    run = function(m = integer(), j = integer(), alpha = double()) {
        returnType(double())
        l <- numeric(m, init=FALSE)
        ## use same sigma as for t=1
        coolParam <- (alpha)^(((j - 1)*timeLength)/(50*timeLength))
        coolSigma <- coolParam*sigma
        
        for(i in 1:m) {
            nimCopy(mvEWSamples, model, nodes = paramNodes, row = i)
            currentValues <- values(model, paramNodes)
            for(j in 1:numParams)
                currentValues[j] <- rnorm(1, currentValues[j], coolSigma[j])
            values(model, paramNodes) <<- currentValues
            model$calculate(parDeterm)
            model$simulate(calc_thisNode_self)
            model$calculate(calc_thisNode_deps)
            nimCopy(model, mvEWSamples, nodes = baselineNode, nodesTo = latentVar, rowTo = i)
            nimCopy(model, mvEWSamples, nodes = paramNodes, rowTo = i)
        }
        return(0)
    },  where = getLoadingNamespace()
)

IF2StepVirtual <- nimbleFunctionVirtual(
    run = function(m = integer(), n = integer(), alpha = double())
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
        ## Remove any element of parDeterm that is not needed by prevDeterm, calc_thisNode_self or calc_thisNode_deps
        nodes_that_matter <- c(prevDeterm, calc_thisNode_self, calc_thisNode_deps)
        if(length(parDeterm) > 0) {
            thisDeterm_is_intermediate <- logical(length(parDeterm))
            for(i in seq_along(parDeterm)) {
                theseDeps <- model$getDependencies(parDeterm[i])
                thisDeterm_is_intermediate[i] <- !(parDeterm[i] %in% nodes_that_matter) & any(theseDeps %in% nodes_that_matter)
            }
            parDeterm <- parDeterm[thisDeterm_is_intermediate]
        }
        parAndPrevDeterm <- c(parDeterm, prevDeterm)
        parAndPrevDeterm <- model$expandNodeNames(parAndPrevDeterm, sort = TRUE) ## Ensure sorting, because a node in parDeterm that is also in prevDeterm will have been removed from parDeterm, but that is where it potentially comes first.
        
        isLast <- (iNode == timeLength)
        coolSigma <- numeric(numParams)
        if(numParams == 1)
            coolSigma <- c(coolSigma, 0)
        coolParam <- 0
    },
    run = function(m = integer(), j = integer(), alpha = double()) {
        returnType(double())
        l <- numeric(m, init=FALSE)
        wts <- numeric(m, init=FALSE)
        ids <- integer(m, 0)
        coolParam <<- (alpha)^((iNode-1+(j - 1)*timeLength)/(50*timeLength))
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
            model$calculate(parAndPrevDeterm)
            model$simulate(calc_thisNode_self)
            logProb <- model$calculate(calc_thisNode_deps)
            wts[i]  <- exp(logProb)
            if(is.nan(wts[i])) wts[i] <- 0
            logProb <- model$calculate(paramNodes)
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
#' @param nodes  A character vector specifying the latent model nodes 
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.  All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param params A character vector specifying the top-level parameters to obtain maximum likelihood estimates of. 
#'   If unspecified, parameter nodes are specified as all stochastic top level nodes which
#'  are not in the set of latent nodes specified in \code{nodes}.
#' @param baselineNode A character vector specifying the node that is the latent node at the "0th" time step. The first node in \code{nodes} should depend on this baseline, but \code{baselineNode} should have no data depending on it. If \code{NULL} (the default), any initial state is taken to be fixed at the values present in the model at the time the algorithm is run. 
#' @param control  A list specifying different control options for the IF2 algorithm.  Options are described in the \sQuote{details} section below.

#' @author Nicholas Michaud, Dao Nguyen, and Christopher Paciorek
#' @family particle filtering methods
#' @details 
#' 
#' Each of the \code{control()} list options are described in detail below:
#' \describe{
#'  \item{sigma}{A vector specifying a non-negative perturbation magnitude for each element of the \code{params} argument.  Defaults to a vector of 1's.}
#'  \item{initParamSigma}{An optional vector specifying a vector of standard deviations to use when simulating an initial particle swarm centered on the initial value of the parameters. Defaults to \code{sigma}.}
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
#'    \item{m}{A single integer specifying the number of particles to use for each run of the filter. }
#'    \item{n}{A single integer specifying the number of overall filter iterations to run. }
#'    \item{alpha}{A double specifying the cooling factor to use for the IF2 algorithm. }
#'
#'  The \code{run} fuction will return a vector with the estimated MLE.  Additionally, once the specialized algorithm has been run, it can be continued for additional iterations by calling the \code{continueRun} method.
#'
#' @section Reparameterization:
#'
#' The IF2 algorithm perturbs the parameters using a normal distribution, which may not be optimal for parameters whose support is not the whole real line, such as variance parameters, which are restricted to be positive. We recommend that users reparameterize the model in advance, e.g., writing variances and standard deviations on the log scale and probabilities on the logit scale. This requires specifying priors directly on the transformed parameters. 
#'  
#' @section Parameter prior distributions:
#'
#' While NIMBLE's IF2 algorithm requires prior distributions on the parameters, the IF2 algorithm produces maximum likelihood estimates and does not directly use those prior distributions. We require the prior distributions to be stated only so that we can automatically determine which model nodes are the parameters. The IF2 algorithm also makes use of any bounds on the parameters.
#'
#' @section Diagnostics and information stored in the algorithm object:
#'
#' The IF2 algorithm stores the estimated MLEs, one from each iteration, in \code{estimates}. It also stores standard deviation of the particles from each iteration, one per parameter, in \code{estSD}. Finally it stores the estimated log-likelihood at the estimated MLE from each iteration in \code{logLik}.
#'  
#' @export
#' 
#' @references Ionides, E.L., D. Nguyen, Y. Atchad{\'e}, S. Stoev, and A.A. King (2015). Inference for dynamic and latent variable models via iterated, perturbed Bayes maps. \emph{Proceedings of the National Academy of Sciences}, 112(3), 719-724.
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
#' # visualize progression of the estimated log-likelihood
#' ts.plot(CmyIF2$logLik)
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

        if(!is.null(baselineNode)) {
            baselineNode <- model$expandNodeNames(baselineNode)
            if(length(model$getDependencies(baselineNode, dataOnly = TRUE)))
                stop("buildIteratedFilter2: 'baselineNode' should not have any data nodes as dependents.")
        }
        
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
        if(!all(params %in% model$getNodeNames(stochOnly = TRUE)))
            stop('buildIteratedFilter2: Parameters must be stochastic nodes.')
        
        paramVars <-  model$getVarNames(nodes = params)
        
        ## Be sure baselineNode is not in paramVars
        if(!is.null(baselineNode)) {
            if(any(baselineNode %in% params)) {
                params <- params[!params %in% baselineNode]
                cat("buildIteratedFilter2: removing baselineNode from parameters.\n")
            }
        }
        cat("buildIteratedFilter2: IF2 algorithm built to maximize parameters: ", paste(params, collapse = ', '), ".\n", sep = '')
        
        if(is.null(sigma)){
            sigma <- rep(1, numParams)
        }
        if(length(sigma) != numParams)
            stop("buildIteratedFilter2: The 'sigma' control list argument must be a vector specifying a
            non-negative perturbation magnitude for each element of 'params'. The length of 'sigma'
            does not match the length of 'params'.")
        
        if(any(sigma < 0))
            stop("buildIteratedFilter2: All values of 'sigma' should be non-negative.")

        if(is.null(initParamSigma)) initParamSigma <- sigma
        if(length(initParamSigma) != numParams)
            stop("buildIteratedFilter2: If non-NULL, the 'initParamSigma' control list argument must be a vector specifying a non-negative perturbation magnitude for each element of 'params'. The length of 'initParamSigma' does not match the length of 'params'.")


        ## Check for inits values
        if(is.null(inits)){
            inits <- values(model, params)
        }
        if(length(inits) < numParams)
            stop("buildIteratedFilter2: The 'inits' control list argument must be a vector specifying initial values for each element of 'params'. The length of 'inits' does not match the number of parameters.")
        if(length(inits) == 1)
            inits <- c(inits, 0)
        
        latentVar <- model$getVarNames(nodes = nodes)
        ## This check may not strictly be needed.
        if(length(unique(latentVar)) > 1){
            stop("buildIteratedFilter2: All latent nodes must be in the same variable.")
        }
        nodes <- findLatentNodes(model, nodes, timeIndex)
        timeLength <- length(nodes)
        
        expandedNodes <- model$expandNodeNames(nodes)
        if(any(baselineNode %in% expandedNodes))
            stop("buildIteratedFilter2: 'baselineNode' cannot be part of latent states.")
        if(any(params %in% expandedNodes))
            stop('buildIteratedFilter2: Parameters cannot be latent states.')

        dims <- lapply(nodes, function(node) nimDim(model[[node]]))
        if(length(unique(dims)) > 1) stop('buildIteratedFilter2: sizes or dimension of latent states varies; this is not allowed.')

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

        ## need vector for compilation
        if(length(initParamSigma) == 1)
            initParamSigma <- c(initParamSigma, 0)
        if(length(sigma) == 1)
            sigma <- c(sigma, 0)
        
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

        if(numParams == 1) {
            estimate <- nimNumeric(2)  # need vector for compilation
        } else estimate <- nimNumeric(numParams)  
        
        oldJ <- 0
        oldM <- 0
        logLik <- nimNumeric(2)
        estimates <- nimMatrix(0, 1, numParams)
        estSD <- nimMatrix(0, 1, numParams)
    },
    run = function(m = integer(),
                   niter = integer(), 
                   alpha = double()) {
        
        values(model, params) <<- inits[1:numParams]
        my_initializeModel$run()
        resize(mvWSamples, m)
        resize(mvEWSamples, m)
        estSD <<- nimMatrix(0, niter, numParams)
        estimates <<- nimMatrix(0, niter, numParams)
        logLik <<- nimNumeric(niter, value = 0)

        initializeParamSwarmFunction$run(m)
        for(j in 1:niter){
            checkInterrupt()
            ## Initialize latent process and params at time 0.
            if(baseline)
                IF2Step0Function$run(m, j, alpha)
            for(iNode in seq_along(IF2StepFunctions)) {
                logLik[j] <<- logLik[j] + log(IF2StepFunctions[[iNode]]$run(m, j, alpha))
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

        estimate <<- estimates[niter, ]

        ## Put final parameter estimates into model.
        values(model, params) <<- estimate  
        model$calculate(parDeterm)
        oldM <<- m
        oldJ <<- j 
        
        returnType(double(1))
        return(estimate)
    },
    methods = list(
        continueRun = function(niter = integer(), alpha = double()) {
            tmpEstimates <- estimates
            tmpEstSD <- estSD
            tmpLogLik <- logLik
            niter_old <- dim(estimates)[1]
            niter_total <- niter_old + niter
            setSize(estimates, niter_total, numParams, fillZeros = TRUE)
            setSize(estSD, niter_total, numParams, fillZeros = TRUE)
            setSize(logLik, niter_total, fillZeros = TRUE)
            estimates[1:niter_old, ] <<- tmpEstimates
            estSD[1:niter_old, ] <<- tmpEstSD
            logLik[1:niter_old] <<- tmpLogLik
            
            useStoredSamples <- 1
            newN <- oldJ + niter
            for(j in (oldJ+1):newN){
                checkInterrupt()
                if(baseline) 
                    IF2Step0Function$run(oldM, j, alpha)
                for(iNode in seq_along(IF2StepFunctions)) {
                    logLik[j] <<- logLik[j] + log(IF2StepFunctions[[iNode]]$run(oldM, j, alpha))
                }
                ## Compute estimate and sd of particles at each iteration for diagnostics.
                for(i in 1:oldM) {
                    nimCopy(mvEWSamples, model, nodes = params, row = i)
                    estimates[j, ] <<- estimates[j, ] + values(model, params) 
                }
                estimates[j, ] <<- estimates[j, ] / oldM
                for(i in 1:oldM) {
                    nimCopy(mvEWSamples, model, nodes = params, row = i)
                    estSD[j, ] <<- estSD[j, ] + (values(model, params) - estimates[j,])^2
                }
                estSD[j, ] <<- sqrt(estSD[j, ] / oldM)
            }

            estimate <<- estimates[niter_total, ]

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
