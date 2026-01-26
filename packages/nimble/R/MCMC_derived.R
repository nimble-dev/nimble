
####################################################################
### virtual nimbleFunction template for all derived quantities #####
####################################################################

#' @rdname derived
#' @export
derived_BASE <- nimbleFunctionVirtual(
    name = 'derived_BASE',
    run = function(timesRan = double()) { },
    methods = list(
        set_interval = function(newInterval = double()) { },
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), chain = double()) { },
        after_chain  = function() { },
        get_results  = function() { returnType(double(2))    },
        get_names    = function() { returnType(character(1)) },
        reset        = function() { }
    ),
    methodControl = list(
        before_chain = list(required = FALSE),
        after_chain  = list(required = FALSE),
        reset        = list(required = FALSE)
    )
)



####################################################################
### derived quantity: mean #########################################
####################################################################

#' @rdname derived
#' @export
derived_mean <- nimbleFunction(
    name = 'derived_mean',
    contains = derived_BASE,
    setup = function(model, mcmc, interval, control) {
        ## control list extraction
        nodes              <- extractControlElement(control, 'nodes',              defaultValue = character())
        recordingFrequency <- extractControlElement(control, 'recordingFrequency', defaultValue = 0)
        ## node list generation
        if(is.list(nodes))   nodes <- unlist(nodes)
        nodes <- model$expandNodeNames(nodes)
        ## names generation
        names <- if(length(nodes) < 2) c(nodes,'','') else nodes     ## vector
        ## numeric value generation
        nSamples <- 0
        nResults <- length(nodes)
        vals       <- numeric(max(nResults, 2))    ## vector
        onlineMean <- numeric(max(nResults, 2))    ## vector
        results <- array(0, c(1, nResults))
    },
    run = function(timesRan = double()) {
        if(nResults == 0)   return()
        nSamples <<- nSamples + 1
        vals <<- values(model, nodes)
        if(nSamples == 1) {
            onlineMean <<- vals
        } else {
            onlineMean <<- onlineMean + (vals - onlineMean) / nSamples
        }
        if(recordingFrequency != 0 & timesRan %% recordingFrequency == 0) {
           results[timesRan/recordingFrequency,] <<- onlineMean
        }
    },
    methods = list(
        set_interval = function(newInterval = double()) {
            interval <<- newInterval
        },
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), chain = double()) {
            if(recordingFrequency == 0) {
                nKeep <- 1
            } else {
                nKeep <- floor(niter / (interval * recordingFrequency))
            }
            results <<- nimArray(NA, c(nKeep, nResults))
        },
        after_chain = function() {
            if(recordingFrequency == 0) {
                results[1,] <<- onlineMean
            }
        },
        get_results = function() {
            returnType(double(2))
            return(results)
        },
        get_names = function() {
            returnType(character(1))
            return(names)
        },
        reset = function() {
            nSamples <<- 0
            vals       <<- nimNumeric(nResults, value = NA)
            onlineMean <<- nimNumeric(nResults, value = NA)
        }
    )
)



####################################################################
### derived quantity: variance #####################################
####################################################################

#' @rdname derived
#' @export
derived_variance <- nimbleFunction(
    name = 'derived_variance',
    contains = derived_BASE,
    setup = function(model, mcmc, interval, control) {
        ## control list extraction
        nodes              <- extractControlElement(control, 'nodes',              defaultValue = character())
        recordingFrequency <- extractControlElement(control, 'recordingFrequency', defaultValue = 0)
        ## node list generation
        if(is.list(nodes))   nodes <- unlist(nodes)
        nodes <- model$expandNodeNames(nodes)
        ## names generation
        names <- if(length(nodes) < 2) c(nodes,'','') else nodes  ## vector
        ## numeric value generation
        nSamples <- 0
        nResults <- length(nodes)
        vals    <- numeric(max(nResults, 2))          ## vector
        prvMean <- numeric(max(nResults, 2))          ## vector
        newMean <- numeric(max(nResults, 2))          ## vector
        sumSqur <- numeric(max(nResults, 2))          ## vector
        results <- array(0, c(1, nResults))
    },
    run = function(timesRan = double()) {
        if(nResults == 0)   return()
        nSamples <<- nSamples + 1
        vals <<- values(model, nodes)
        ## Welford's algorithm for stable online variance
        if(nSamples == 1) {
            newMean <<- vals
            sumSqur <<- nimNumeric(nResults)
        } else {
            prvMean <<- newMean
            newMean <<- prvMean + (vals - prvMean) / nSamples
            sumSqur <<- sumSqur + (vals - prvMean) * (vals - newMean)
        }
        if(recordingFrequency != 0 & timesRan %% recordingFrequency == 0) {
            if(nSamples == 1) {
                results[timesRan/recordingFrequency,] <<- rep(NA, nResults)
            } else {
                results[timesRan/recordingFrequency,] <<- sumSqur / (nSamples-1)
            }
        }
    },
    methods = list(
        set_interval = function(newInterval = double()) {
            interval <<- newInterval
        },
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), chain = double()) {
            if(recordingFrequency == 0) {
                nKeep <- 1
            } else {
                nKeep <- floor(niter / (interval * recordingFrequency))
            }
            results <<- nimArray(NA, c(nKeep, nResults))
        },
        after_chain = function() {
            if(recordingFrequency == 0) {
                if(nSamples == 1) {
                    results[1,] <<- rep(NA, nResults)
                } else {
                    results[1,] <<- sumSqur / (nSamples-1)
                }
            }
        },
        get_results = function() {
            returnType(double(2))
            return(results)
        },
        get_names = function() {
            returnType(character(1))
            return(names)
        },
        reset = function() {
            nSamples <<- 0
            vals     <<- nimNumeric(nResults, value = NA)
            prvMean  <<- nimNumeric(nResults, value = NA)
            newMean  <<- nimNumeric(nResults, value = NA)
            sumSqur  <<- nimNumeric(nResults, value = NA)
        }
    )
)



####################################################################
### derived quantity: logProb ######################################
####################################################################

getLogProb_virtual <- nimbleFunctionVirtual(
    run = function() { returnType(double()) }
)

getLogProbNF <- nimbleFunction(
    contains = getLogProb_virtual,
    setup = function(model, nodes) {},
    run = function() {
        returnType(double())
        return(model$getLogProb(nodes))
    }
)

#' @rdname derived
#' @export
derived_logProb <- nimbleFunction(
    name = 'derived_logProb',
    contains = derived_BASE,
    setup = function(model, mcmc, interval, control) {
        ## control list extraction
        nodes  <- extractControlElement(control, 'nodes',  defaultValue = '.all')
        silent <- extractControlElement(control, 'silent', defaultValue = FALSE)
        ## node list generation
        nodeList <- if(is.character(nodes)) {
                        as.list(unlist(lapply(nodes, function(x) if(identical(x,'.all')) '.all' else model$expandNodeNames(x))))
                    } else nodes
        allBool <- sapply(nodeList, function(x) identical(x, '.all'))
        nodeList <- lapply(nodeList, function(x) if(identical(x,'.all')) model$getNodeNames(stochOnly=TRUE) else x)
        ## names generation
        if(is.list(nodes) && !is.null(names(nodes))) {
            names <- names(nodes)
        } else {
            names <- character(length(nodeList))
            sumIndex <- 0
            for(i in seq_along(nodeList)) {
                if(allBool[i]) {
                    names[i] <- '_all_nodes_'
                    next }
                if(length(nodeList[[i]]) > 1) {
                    sumIndex <- sumIndex + 1
                    names[i] <- paste0('sum', sumIndex)
                    next }
                names[i] <- nodeList[[i]]
            }
        }
        if(length(names) < 2)   names <- c(names, '', '')     ## vector
        ## numeric value generation
        nResults <- length(nodeList)
        results <- array(0, c(1, nResults))
        ## nested function and function list definitions
        getLogProbNFL <- nimbleFunctionList(getLogProb_virtual)
        for(i in seq_along(nodeList))   getLogProbNFL[[i]] <- getLogProbNF(model, nodeList[[i]])
        ## checks
        uniqueNodes <- unique(unlist(nodeList))
        missingInd <- sapply(uniqueNodes, function(x) length(model$expandNodeNames(x)) == 0)
        missingNodes <- uniqueNodes[missingInd]
        if(length(missingNodes) && !silent)
            warning('logProb derived quantity function is using node names which are not in the model: ',
                    paste0(missingNodes, collapse=', '), call. = FALSE)
    },
    run = function(timesRan = double()) {
        if(nResults == 0)   return()
        for(i in 1:nResults) {
            results[timesRan, i] <<- getLogProbNFL[[i]]$run()
        }
    },
    methods = list(
        set_interval = function(newInterval = double()) {
            interval <<- newInterval
        },
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), chain = double()) {
            nKeep <- floor(niter / interval)
            results <<- nimArray(NA, c(nKeep, nResults))
        },
        get_results = function() {
            returnType(double(2))
            return(results)
        },
        get_names = function() {
            returnType(character(1))
            return(names)
        },
        reset = function() {
            results <<- nimArray(NA, c(1, nResults))
        }
    )
)



####################################################################
### derived quantity: predictive ###################################
####################################################################

#' @rdname derived
#' @export
derived_predictive <- nimbleFunction(
    name = 'derived_predictive',
    contains = derived_BASE,
    setup = function(model, mcmc, interval, control) {
        ## control list extraction
        nodes    <- extractControlElement(control, 'nodes',    defaultValue = '.missing')
        simNodes <- extractControlElement(control, 'simNodes', defaultValue = '.missing')
        sort     <- extractControlElement(control, 'sort',     defaultValue = TRUE)
        silent   <- extractControlElement(control, 'silent',   defaultValue = FALSE)
        ## node list generation
        ## all predictive stochastic nodes, and their deterministic dependencies
        ppNodesAndDeps <- model$getDependencies(model$getNodeNames(predictiveOnly = TRUE), downstream = TRUE)
        if(length(nodes) == 1 && nodes == '.missing') {
            saveNodes <- model$getNodeNames(endOnly = TRUE, includeData = FALSE)
            isEnd <- sapply(saveNodes, function(x) length(model$getDependencies(x, self = FALSE)) == 0)
            saveNodes <- saveNodes[isEnd]
        } else if(length(nodes) == 1 && nodes == '.all') {
            ## this approach missed predicted deterministic nodes, which also have a stochastic dependency:
            ##saveNodes <- unique(c(
            ##    ppNodesAndDeps,                                            ## predictive stochastic nodes and deterministic dependencies
            ##    model$getNodeNames(endOnly = TRUE, includeData = FALSE)    ## deterministic derived quantities
            ##))
            allModelNodes <- model$getNodeNames()
            allNodeDownstreamDeps <- lapply(allModelNodes, function(x) model$getDependencies(x, downstream = TRUE))
            haveDataDepsBool <- sapply(allNodeDownstreamDeps, function(x) any(model$isData(x)))
            saveNodes <- allModelNodes[!haveDataDepsBool]
            saveNodes <- model$topologicallySortNodes(saveNodes)
        } else {
            saveNodes <- model$expandNodeNames(nodes)
        }
        sampledNodesList <- lapply(
            mcmc$samplerFunctions$contentsList,
            function(x) {
                nodes <- character()
                if('target'         %in% ls(x))   nodes <- c(nodes, x$target)
                if('targetAsScalar' %in% ls(x))   nodes <- c(nodes, x$targetAsScalar)
                if('simNodes'       %in% ls(x))   nodes <- c(nodes, x$simNodes)
                return(nodes)
            })
        sampledNodes <- model$expandNodeNames(unlist(sampledNodesList))
        if(simNodes == '.missing') {
            ## figuring this case out was not easy -DT May 2025
            upToDateNodes <- setdiff(     ## nodes which should be kept up-to-date by the MCMC
                model$getNodeNames(includePredictive = FALSE),   ## all model nodes, excluding predictive stochastic nodes
                ppNodesAndDeps                                   ## predictive stochastic nodes and deterministic dependencies
            )
            ## now we add any predictive nodes (and their deterministic dependencies), which might have samplers assigned
            sampledNodeDeps <- model$getDependencies(sampledNodes)
            sampledNodeDetermDeps <- sampledNodeDeps[model$isDeterm(sampledNodeDeps)]
            upToDateNodes <- unique(c(
                upToDateNodes,
                sampledNodes,             ## sampled stochastic nodes
                sampledNodeDetermDeps     ## deterministic dependencies of sampled nodes
            ))
            saveNodesParents <- model$getParents(saveNodes, self = TRUE, upstream = TRUE)   ## everything upstream from (and including) saveNodes
            simNodes <- setdiff(saveNodesParents, upToDateNodes)
        }
        if(sort)   simNodes <- model$topologicallySortNodes(simNodes)
        calcNodes <- model$getDependencies(simNodes)
        mvSaved <- mcmc$mvSaved
        ## names generation
        names <- if(length(saveNodes) < 2) c(saveNodes,'','') else saveNodes     ## vector
        ## numeric value generation
        nResults <- length(saveNodes)
        results <- array(0, c(1, nResults))
        ## checks
        otherwiseSampledSimNodes <- intersect(simNodes, sampledNodes)
        if(length(otherwiseSampledSimNodes) & !silent) {
            message('  [Warning] predictive derived quantity function is simulating nodes being updated by other MCMC samplers: ',
                    paste0(otherwiseSampledSimNodes, collapse=', '))
        }
        stochSaveNodes <- saveNodes[model$isStoch(saveNodes)]
        nonPPsaveNodes <- setdiff(stochSaveNodes, ppNodesAndDeps)
        if(length(nonPPsaveNodes) & !silent) {
            message('  [Warning] predictive derived quantity function is operating on non-posterior-predictive nodes: ',
                    paste0(nonPPsaveNodes, collapse=', '))
        }
    },
    run = function(timesRan = double()) {
        if(nResults == 0)   return()
        model$simulate(simNodes)
        model$calculate(calcNodes)
        results[timesRan,] <<- values(model, saveNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        set_interval = function(newInterval = double()) {
            interval <<- newInterval
        },
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), chain = double()) {
            nKeep <- floor(niter / interval)
            results <<- nimArray(NA, c(nKeep, nResults))
        },
        get_results = function() {
            returnType(double(2))
            return(results)
        },
        get_names = function() {
            returnType(character(1))
            return(names)
        },
        reset = function() {
            results <<- nimArray(0, c(1, nResults))
        }
    )
)



#' MCMC Derived Quantities
#'
#' Details of the NIMBLE MCMC engine handling of derived quantities, which are deterministic functions that can be calaculated and recorded after each MCMC sampling iteration.
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' @param mcmc (uncompiled) MCMC object
#' @param interval interval (of MCMC iterations) at which the derived quantity is calculated
#' @param control named list that controls the precise behavior of the derived quantity calculation, with elements specific to  type of derived quantity.
#' 
#' @section Mean and variance:
#'
#' The \code{mean} and \code{variance} derived quantity functions calculate the running mean and variance, respectively, for each node specified in the \code{nodes} argument.  If added to an MCMC configuration object using the \code{addDerivedQuantity} method, then a value of the \code{interval} argument may also be provided to \code{addDerivedQuantity}. In that case, the value of \code{interval} specifies the number of MCMC iterations between calculations of the statistic.  When the statistic is calculated, only the current value of each node is used to update the statistic.  For example, if \code{interval} is 2, then every other MCMC iteration is used to calculate an updated value of the statistic.  If no value of \code{interval} is provided as an argument to \code{addDerivedQuantity}, then the default value is the thinning interval \code{thin} of the MCMC.
#'
#' The \code{mean} and \code{variance} derived quantity functions both accept the following control list elements:
#' \itemize{
#' \item nodes. The set of model nodes used for tracking the statistic.
#' \item recordingFrequency. The frequency (number of calculations of the statistic) afer which the value of the statistic is saved.  For example, if \code{recordingFrequency} is 1, then the value of the statistic is saved after every update of its value.  But if \code{recordingFrequency} is 10, then the value of the statistic is only saved after every tenth update of its value.  The dafault value of \code{recordingFrequency} is 0, which corresponds to a special case: the value of the statistic is only recorded a single time, which is on the final iteration of the MCMC chain.
#' }
#'
#' @section Model log-densities:
#'
#' The \code{logProb} derived quantity function calculates and records values of the log-density of individual nodes or (summed) groups of nodes.   If added to an MCMC configuration object using the \code{addDerivedQuantity} method, then a value of the \code{interval} argument may also be provided to \code{addDerivedQuantity}. In that case, the value of \code{interval} specifies the number of MCMC iterations between recordings of the log-density values.  For example, if \code{interval} is 2, then log-density values will be recorded upon every other MCMC iteration.  If no value of \code{interval} is provided as an argument to \code{addDerivedQuantity}, then the default value is the thinning interval \code{thin} of the MCMC.
#'
#' The \code{logProb} derived quantity function accepts the following control list elements:
#' \itemize{
#' \item nodes. The \code{nodes} argument determines the individual nodes, or (summed) groups of nodes, for recording log-density values.  When provided as a character vector, the individual log-density of each node in this vector will be recorded.  When provided as a list, each list element may contain one or mode node names, and separately for the node(s) in each element of the list, the summed log-density list will be calculated.  In addition, the keyword \code{".all"} may also be provided in either the vector or list argument, which corresponds to the set of all stochastic model nodes (including data).
#' \item silent.  By default, the \code{logProb} derived quantity function will issue a warning when the \code{nodes} argument includes node names which are not present in the model.  This warning may be suppressed by setting \code{silent} to \code{TRUE}.
#' }
#'
#' @section Posterior predictive nodes and derived quantities:
#'
#' The \code{predictive} derived quantity function simulates the values of posterior predictive nodes in the model and stores these simulated values.  This may be useful when a model structure includes posterior predictive nodes (or deterministically defined posterior derived quantities), but for reasons of efficiency, these nodes may not undergo MCMC sampling.  In such cases, the \code{predictive} derived quantity function may be assigned to these nodes, and when executed it will simulate new values for these nodes and record the simulated values.
#'
#' When added to an MCMC configuration object using the \code{addDerivedQuantity} method, a value of the \code{interval} argument may also be provided to \code{addDerivedQuantity}. In that case, the value of \code{interval} specifies the number of MCMC iterations between operations of the \code{predictive} function.  For example, if \code{interval} is 2, then prediction and storing values takes place every other MCMC iteration.  If no value of \code{interval} is provided as an argument to \code{addDerivedQuantity}, then the default value is the thinning interval \code{thin} of the MCMC.
#'
#' The \code{predictive} derived quantity function accepts the following control list elements:
#' \itemize{
#' \item nodes. The \code{nodes} argument defines the predictive nodes which will have their values recorded and returned.  The \code{predictive} function will automatically determine which nodes upstream of \code{nodes} need to be simulated, prior to storing the values of \code{nodes}.  If omitted,  the default set of \code{nodes} will be all (non-data) terminal model nodes.  In addition, specifying the value \code{nodes = ".all"} will take \code{nodes} to be all predictive nodes in the model (all deterministic and stochastic nodes which have no downstream data dependencies).
#' \item simNodes. The \code{simNodes} specifies the set of model nodes which will be simulated, prior to recording the values of \code{nodes}.  If ommited, which would be the usual case, \code{simNodes} will automatically represent the set of all model nodes which must be simulated in order to reach the predictive \code{nodes}.  By providing \code{simNodes}, more fine-grained control of the simulation process is possible, including omitting simulation of some nodes.
#' \item sort. The \code{sort} argument determines whether the simulation of \code{simNodes} takes place in topological order.  This argument has a default value of \code{TRUE}.  When specified as \code{FALSE} the simulation of \code{simNodes} will take place in the order in which they were specified in the \code{simNodes} argument, which may not be in their natural order of dependency.
#' \item silent.  By default, the \code{predictive} derived quantity function will issue a warning when the \code{simNodes} argument includes node names which are being updated by some sampler function in the MCMC.  This warning may be suppressed by setting \code{silent} to \code{TRUE}.
#' }
#'
#' @name derived
#'
#' @aliases derived_mean derived_variance derived_logProb
#'
#' @examples
#' \dontrun{
#' conf <- configureMCMC(model)
#'
#' conf$addDerivedQuantity("mean", nodes = c("a", "b"))
#' 
#' conf$addDerivedQuantity("mean", nodes = "theta", interval = 5)
#' 
#' conf$addDerivedQuantity("variance", nodes = "x[1:4]", control = list(recordingFrequency = 10))
#'
#' conf$addDerivedQuantity("logProb", nodes = c('alpha', 'beta'))
#'
#' conf <- configureMCMC(model, mean = 'a', variance = 'b', logProb = TRUE)
#' }
#' 
#' @seealso \code{\link{configureMCMC}} \code{\link{addDerivedQuantity}} \code{\link{buildMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
#'
#' @author Daniel Turek
#'
NULL

