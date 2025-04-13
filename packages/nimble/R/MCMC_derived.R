
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
### derived quantity: runningMean ##################################
####################################################################

#' @rdname derived
#' @export
derived_runningMean <- nimbleFunction(
    name = 'derived_runningMean',
    contains = derived_BASE,
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        nodes <- extractControlElement(control, 'nodes', defaultValue = character())
        ## node list generation
        nodes <- model$expandNodeNames(nodes)
        ## numeric value generation
        count <- 1
        nResults <- length(nodes)
        vals <- numeric(nResults)
        results <- array(0, c(1, nResults))
    },
    run = function(iter = double()) {
        if(nResults > 0) {
            vals <<- values(model, nodes)
            if(count == 1) {
                results[count,] <<- vals
            } else {
                results[count,] <<- results[count-1,] + (1/count) * (vals - results[count-1,])
            }
            count <<- count + 1
        }
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            nKeep <- floor(niter / interval)
            setSize(results, nKeep, nResults) },
        getResults = function() {
            returnType(double(2))
            return(results) },
        getNames = function() {
            returnType(character(1))
            return(nodes) },
        reset = function() {
            count <<- 1
            results <<- array(0, c(1, nResults)) }
    )
)



####################################################################
### derived quantity: runningVariance ##############################
####################################################################

#' @rdname derived
#' @export
derived_runningVariance <- nimbleFunction(
    name = 'derived_runningVariance',
    contains = derived_BASE,
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        nodes <- extractControlElement(control, 'nodes', defaultValue = character())
        ## node list generation
        nodes <- model$expandNodeNames(nodes)
        ## numeric value generation
        count <- 1
        nResults <- length(nodes)
        vals <- numeric(nResults)
        runMean <- array(0, c(1, nResults))
        sumSqur <- array(0, c(1, nResults))
        results <- array(0, c(1, nResults))
    },
    run = function(iter = double()) {
        if(nResults > 0) {
            vals <<- values(model, nodes)
            if(count == 1) {
                runMean[count,] <<- vals
                sumSqur[count,] <<- rep(0,  nResults)
                results[count,] <<- rep(NA, nResults)
            } else {
                ## Welford's algorithm for stable online variance
                runMean[count,] <<- runMean[count-1,] + (1/count) * (vals - runMean[count-1,])
                sumSqur[count,] <<- sumSqur[count-1,] + (vals - runMean[count-1,]) * (vals - runMean[count,])
                results[count,] <<- sumSqur[count,] / (count-1)
            }
            count <<- count + 1
        }
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            nKeep <- floor(niter / interval)
            setSize(runMean, nKeep, nResults)
            setSize(sumSqur, nKeep, nResults)
            setSize(results, nKeep, nResults) },
        getResults = function() {
            returnType(double(2))
            return(results) },
        getNames = function() {
            returnType(character(1))
            return(nodes) },
        reset = function() {
            count <<- 1
            runMean <<- array(0, c(1, nResults))
            sumSqur <<- array(0, c(1, nResults))
            results <<- array(0, c(1, nResults)) }
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
        nodes  <- extractControlElement(control, 'nodes', defaultValue = '.all')
        silent <- extractControlElement(control, 'silent',   defaultValue = FALSE)
        ## node list generation
        nodeList <- if(is.character(nodes)) list(nodes) else nodes
        names <- character(length(nodeList))
        sumIndex <- 0
        for(i in seq_along(nodeList)) {
            if(identical(nodeList[[i]], '.all')) {
                names[i] <- '_all_nodes_'
                nodeList[[i]] <- model$getNodeNames(stochOnly = TRUE)
                next }
            if(length(nodeList[[i]]) > 1) {
                sumIndex <- sumIndex + 1
                names[i] <- paste0('sum', sumIndex)
                next }
            names[i] <- nodeList[[i]]
        }
        ## numeric value generation
        count <- 1
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
            for(i in 1:nResults)   results[count, i] <<- getLogProbNFL[[i]]$run()
            count <<- count + 1
        }
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            nKeep <- floor(niter / interval)
            setSize(results, nKeep, nResults) },
        getResults = function() {
            returnType(double(2))
            return(results) },
        getNames = function() {
            returnType(character(1))
            return(names) },
        reset = function() {
            count <<- 1
            results <<- array(0, c(1, nResults)) }
    )
)

