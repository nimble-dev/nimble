
####################################################################
### virtual nimbleFunction template for all derived quantities #####
####################################################################

#' @rdname derived
#' @export
derived_BASE <- nimbleFunctionVirtual(
    name = 'derived_BASE',
    run = function(iter = double()) { },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) { },
        after_chain  = function() { },
        getResults   = function() { returnType(double(2))    },
        getNames     = function() { returnType(character(1)) },
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
#' @rdname derived
#' @export
derived_mean <- nimbleFunction(
    name = 'derived_mean',
    contains = derived_BASE,
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        nodes              <- extractControlElement(control, 'nodes',              defaultValue = character())
        recordingFrequency <- extractControlElement(control, 'recordingFrequency', defaultValue = 1)
        ## node list generation
        nodes <- model$expandNodeNames(nodes)
        ## names generation
        names <- if(length(nodes) == 1) c(nodes,'') else nodes    ## vector
        ## numeric value generation
        nSamples <- 0
        saveFrequency <- 0
        nextResultsRow <- 1
        nResults <- length(nodes)
        vals       <- numeric(max(nResults, 2))    ## vector
        onlineMean <- numeric(max(nResults, 2))    ## vector
        results <- array(0, c(1, nResults))
    },
    run = function(iter = double()) {
        if(nResults == 0)   return()
        nSamples <<- nSamples + 1
        vals <<- values(model, nodes)
        onlineMean <<- onlineMean + (vals - onlineMean) / nSamples
        if(iter %% saveFrequency == 0) {
            results[nextResultsRow,] <<- onlineMean
            nextResultsRow <<- nextResultsRow + 1
        }
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            saveFrequency <<- interval * recordingFrequency
            nKeep <- floor(niter / saveFrequency)
            setSize(results, nKeep, nResults)
        },
        getResults = function() {
            returnType(double(2))
            return(results)
        },
        getNames = function() {
            returnType(character(1))
            return(names)
        },
        reset = function() {
            nSamples <<- 0
            nextResultsRow <<- 1
            vals       <<- numeric(nResults)
            onlineMean <<- numeric(nResults)
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
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        nodes              <- extractControlElement(control, 'nodes',              defaultValue = character())
        recordingFrequency <- extractControlElement(control, 'recordingFrequency', defaultValue = 1)
        ## node list generation
        nodes <- model$expandNodeNames(nodes)
        ## names generation
        names <- if(length(nodes) == 1) c(nodes,'') else nodes    ## vector
        ## numeric value generation
        nSamples <- 0
        saveFrequency <- 0
        nextResultsRow <- 1
        nResults <- length(nodes)
        vals    <- numeric(max(nResults, 2))                      ## vector
        prvMean <- numeric(max(nResults, 2))                      ## vector
        newMean <- numeric(max(nResults, 2))                      ## vector
        sumSqur <- numeric(max(nResults, 2))                      ## vector
        results <- array(0, c(1, nResults))
    },
    run = function(iter = double()) {
        if(nResults == 0)   return()
        nSamples <<- nSamples + 1
        vals <<- values(model, nodes)
        ## Welford's algorithm for stable online variance
        if(nSamples == 1) {
            newMean <<- vals
            sumSqur <<- numeric(nResults)
        } else {
            prvMean <<- newMean
            newMean <<- prvMean + (vals - prvMean) / nSamples
            sumSqur <<- sumSqur + (vals - prvMean) * (vals - newMean)
        }
        if(iter %% saveFrequency == 0) {
            if(nSamples == 1) {
                results[nextResultsRow,] <<- rep(NA, nResults)
            } else {
                results[nextResultsRow,] <<- sumSqur / (nSamples-1)
            }
            nextResultsRow <<- nextResultsRow + 1
        }
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            saveFrequency <<- interval * recordingFrequency
            nKeep <- floor(niter / saveFrequency)
            setSize(results, nKeep, nResults)
        },
        getResults = function() {
            returnType(double(2))
            return(results)
        },
        getNames = function() {
            returnType(character(1))
            return(names)
        },
        reset = function() {
            nSamples       <<- 0
            nextResultsRow <<- 1
            setSize(vals,    nResults)
            setSize(prvMean, nResults)
            setSize(newMean, nResults)
            setSize(sumSqur, nResults)
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
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        nodes  <- extractControlElement(control, 'nodes',    defaultValue = '.all')
        silent <- extractControlElement(control, 'silent',   defaultValue = FALSE)
        ## node list generation
        nodeList <- if(is.character(nodes)) {
                        as.list(unlist(lapply(nodes, function(n) if(n=='.all') '.all' else Rmodel$expandNodeNames(n))))
                    } else nodes
        allBool <- sapply(nodeList, function(n) identical(n, '.all'))
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
        if(length(names) == 1)   names <- c(names, '')    ## vector
        ## numeric value generation
        nextResultsRow <- 1
        nResults <- length(nodeList)
        results <- array(0, c(1, nResults))
        ## nested function and function list definitions
        getLogProbNFL <- nimbleFunctionList(getLogProb_virtual)
        for(i in seq_along(nodeList))   getLogProbNFL[[i]] <- getLogProbNF(model, nodeList[[i]])
        ## checks
        uniqueNodes <- unique(unlist(nodeList))
        missingInd <- sapply(uniqueNodes, function(n) length(model$expandNodeNames(n)) == 0)
        missingNodes <- uniqueNodes[missingInd]
        if(length(missingNodes) && !silent)
            warning('logProb derived quantity function is using node names which are not in the model: ',
                    paste0(missingNodes, collapse=', '), call. = FALSE)
    },
    run = function(iter = double()) {
        if(nResults > 0) {
            for(i in 1:nResults)   results[nextResultsRow, i] <<- getLogProbNFL[[i]]$run()
            nextResultsRow <<- nextResultsRow + 1
        }
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            nKeep <- floor(niter / interval)
            setSize(results, nKeep, nResults)
        },
        getResults = function() {
            returnType(double(2))
            return(results)
        },
        getNames = function() {
            returnType(character(1))
            return(names)
        },
        reset = function() {
            nextResultsRow <<- 1
            results <<- array(0, c(1, nResults))
        }
    )
)



#' MCMC Derived Quantities
#'
#' Details of the NIMBLE MCMC engine handles derived quantities, which are deterministic functions that can be calaculated and recorded after each MCMC sampling iteration.
#'
#' @section Running Mean and Variance
#'
#' The \code{mean} and \code{variance} derived quantity functions calculate the running mean and variance, respectively, for each node specified in the \code{nodes} argument.  If added to an MCMC configuration object using the \code{addDerivedQuantity} method, then a value of the \code{interval} argument may also be provided to \code{addDerivedQuantity}. In that case, the value of \code{interval} specifies the number of MCMC iterations between calculations of the running statistic.  When the statistic is calculated, only the current value of each node is used to update the statistic.  For example, if \code{interval} is 2, then every other MCMC iteration is used to calculate an updated value of the running statistic.
#'
#' The \code{mean} and \code{variance} derived quantity functions both accept the following control list elements:
#' \itemize{
#' \item nodes. The set of model nodes used for tracking the running statistic.
#' \item recordingFrequency. The frequency (number of calculations of the running statistic) with which the value of the statistic is saved.  For example, itf \code{recordingFrequency} is 10, then the value of the running statistic is only saved after every tenth update of the running value.
#' }
#'
#' @section Model Log-Densities
#'
#' The \code{logProb} derived quantity function calculates and records values of the log-density of individual nodes or (summed) groups of nodes.   If added to an MCMC configuration object using the \code{addDerivedQuantity} method, then a value of the \code{interval} argument may also be provided to \code{addDerivedQuantity}. In that case, the value of \code{interval} specifies the number of MCMC iterations between recordings of the log-density values.  For example, if \code{interval} is 2, then log-density values will be recorded upon every other MCMC iteration.
#'
#' The \code{logProb} derived quantity function accepts the following control list elements:
#' \itemize{
#' \item nodes. The \code{nodes} argument determines the individual nodes, or (summed) groups of nodes, for recording log-density values.  When provided as a character vector, the individual log-density of each node in this vector will be recorded.  When provided as a list, each list element may contain one or mode node names, and separately for the node(s) in each element of the list, the summed log-density list will be calculated.  In addition, the keyword \code{".all"} may also be provided in either the vector or list argument, which corresponds to the set of all stochastic model nodes (including data).
#' }
#'
#' @name derived
#'
#' @aliases derived_mean derived_variance derived_logProb
#'
#' @examples
#' conf$addDerivedQuantity("mean", nodes = c("a", "b"))
#' 
#' conf$addDerivedQuantity("mean", nodes = "theta", interval = 5)
#' 
#' conf$addDerivedQuantity("variance", nodes = "x[1:4]", control = list(recordingFrequency = 10))
#'
#' conf$addDerivedQuantity("logProb", nodes = c('alpha', 'beta'))
#'
#' conf <- configureMCMC(model, mean = 'a', variance = 'b', logProb = TRUE)
#' 
#' @seealso \code{\link{configureMCMC}} \code{\link{addDerivedQuantity}} \code{\link{buildMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
#'
#' @author Daniel Turek
#'
NULL

