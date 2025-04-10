
####################################################################
### virtual nimbleFunction template for all derived functions ######
####################################################################

#' @rdname derived
#' @export
derived_BASE <- nimbleFunctionVirtual(
    name = 'derived_BASE',
    run = function(iter = double()) { },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) { },
        after_chain  = function() { },
        getResults   = function() { returnType(double(2)) },
        reset        = function() { }
    ),
    methodControl = list(
        before_chain = list(required = FALSE),
        after_chain  = list(required = FALSE),
        reset        = list(required = FALSE)
    )
)



####################################################################
### derived: test ##################################################
####################################################################

#' @rdname derived
#' @export
derived_test <- nimbleFunction(
    name = 'derived_test',
    contains = derived_BASE,
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        count <- 1
        results <- array(0, c(1, 1))
    },
    run = function(iter = double()) {
        results[count, 1] <<- iter
        count <<- count + 1
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            nKeep <- floor(niter/interval)
            setSize(results, nKeep, 1)
        },
        after_chain = function() { },
        getResults = function() {
            returnType(double(2))
            return(results)
        },
        reset = function() {
            results <<- results * 0
            count <<- 1
        }
    )
)



####################################################################
### derived: logProb ###############################################
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
        nodeList <- extractControlElement(control, 'nodeList', defaultValue = list('.all'))
        silent   <- extractControlElement(control, 'silent',   defaultValue = FALSE)
        ## node list generation
        allBool <- sapply(nodeList, function(x) identical('.all', x))
        allInd <- if(length(allBool)) which(allBool) else numeric()
        for(ind in allInd)   nodeList[[ind]] <- model$getNodeNames(stochOnly = TRUE)
        ## numeric value generation
        count <- 1
        nResults <- length(nodeList)
        results <- array(0, c(1, nResults))
        ## nested function and function list definitions
        getLogProbNFL <- nimbleFunctionList(getLogProb_virtual)
        for(i in seq_along(nodeList))   getLogProbNFL[[i]] <- getLogProbNF(model, nodeList[[i]])
        ## checks
        nodes <- unique(unlist(nodeList))
        missingInd <- sapply(nodes, function(n) length(model$expandNodeNames(n)) == 0)
        missingNodes <- nodes[missingInd]
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
            setSize(results, nKeep, nResults)
        },
        getResults = function() {
            returnType(double(2))
            return(results)
        },
        reset = function() {
            count <<- 1
            results <<- array(0, c(1, nResults))
        }
    )
)



####################################################################
### derived: runningMean ###########################################
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
        results <- array(0, c(1, nResults))
    },
    run = function(iter = double()) {
        if(nResults > 0) {
            if(count == 1) {
                results[count, 1:nResults] <<- values(model, nodes)
            } else {
                results[count, 1:nResults] <<- ((count-1)/count) * results[count-1, 1:nResults] + (1/count) * values(model, nodes)
            }
            count <<- count + 1
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
        reset = function() {
            count <<- 1
            results <<- array(0, c(1, nResults))
        }
    )
)
