
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
derived_mean <- nimbleFunction(
    name = 'derived_mean',
    contains = derived_BASE,
    setup = function(model, mvSaved, mvSamples, mvSamples2, interval, control) {
        ## control list extraction
        nodes <- extractControlElement(control, 'nodes', defaultValue = character())
        ## node list generation
        nodes <- model$expandNodeNames(nodes)
        names <- if(length(nodes) == 1) c(nodes,'') else nodes    ## vector
        ## numeric value generation
        burnin <- 0
        thin1 <- 0
        nextRow <- 1
        nSamplesUsed <- 0
        nResults <- length(nodes)
        thisMean <- numeric(max(nResults, 2))    ## vector
        results     <- array(0, c(1, nResults))
        tempStorage <- array(0, c(1, nResults))
    },
    run = function(iter = double()) {
        if(iter < burnin) return()
        if(nResults == 0) return()
        nSamples <- floor((iter-burnin) / thin1)
        nNewSamples <- nSamples - nSamplesUsed
        if(nNewSamples == 0) return()
        if(nextRow == 1)   setSize(tempStorage, nNewSamples, nResults)
        for(i in 1:nNewSamples) {
            nimCopy(from = mvSamples, to = model, row = nSamplesUsed+i, nodes = nodes)
            tempStorage[i,] <<- values(model, nodes)
        }
        nimCopy(from = mvSaved, to = model, row = 1, nodes = nodes)
        for(i in 1:nResults)   thisMean[i] <<- mean(tempStorage[,i])
        if(nextRow == 1) {
            results[nextRow,] <<- thisMean
        } else {
            results[nextRow,] <<- results[nextRow-1,] + (1/nextRow) * (thisMean - results[nextRow-1,])
        }
        nextRow <<- nextRow + 1
        nSamplesUsed <<- nSamples
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            burnin <<- nburnin
            thin1 <<- thin[1]
            setSize(thisMean, nResults)
            nKeep <- floor((niter-nburnin) / interval)
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
            nextRow <<- 1
            nSamplesUsed <<- 0
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
        nodes <- extractControlElement(control, 'nodes', defaultValue = character())
        ## node list generation
        nodes <- model$expandNodeNames(nodes)
        names <- if(length(nodes) == 1) c(nodes,'') else nodes    ## vector
        ## numeric value generation
        burnin <- 0
        thin1 <- 0
        nextRow <- 1
        nSamplesUsed <- 0
        nResults <- length(nodes)
        vals    <- numeric(max(nResults, 2))    ## vector
        prvMean <- numeric(max(nResults, 2))    ## vector
        newMean <- numeric(max(nResults, 2))    ## vector
        sumSqur <- numeric(max(nResults, 2))    ## vector
        results <- array(0, c(1, nResults))
    },
    run = function(iter = double()) {
        if(iter < burnin) return()
        if(nResults == 0) return()
        nSamples <- floor((iter-burnin) / thin1)
        nNewSamples <- nSamples - nSamplesUsed
        if(nNewSamples == 0) return()
        ## Welford's algorithm for stable online variance
        for(i in 1:nNewSamples) {
            nimCopy(from = mvSamples, to = model, row = nSamplesUsed+i, nodes = nodes)
            vals <<- values(model, nodes)
            if(i==1 & nSamplesUsed==0) {
                newMean <<- vals
                sumSqur <<- numeric(nResults)
                next
            }
            prvMean <<- newMean
            newMean <<- prvMean + (vals - prvMean) / (nSamplesUsed+i)
            sumSqur <<- sumSqur + (vals - prvMean) * (vals - newMean)
        }
        if(nSamples == 1) {
            results[nextRow,] <<- rep(NA, nResults)
        } else {
            results[nextRow,] <<- sumSqur / (nSamples-1)
        }
        nextRow <<- nextRow + 1
        nSamplesUsed <<- nSamples
        nimCopy(from = mvSaved, to = model, row = 1, nodes = nodes)
    },
    methods = list(
        before_chain = function(niter = double(), nburnin = double(), thin = double(1), nchains = double()) {
            burnin <<- nburnin
            thin1 <<- thin[1]
            setSize(vals,    nResults)
            setSize(prvMean, nResults)
            setSize(newMean, nResults)
            setSize(sumSqur, nResults)
            nKeep <- floor((niter-nburnin) / interval)
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
            nextRow <<- 1
            nSamplesUsed <<- 0
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
        nodeList <- if(is.character(nodes)) as.list(nodes) else nodes
        names <- character(max(length(nodeList),2))      ## vector
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
        nextRow <- 1
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
            for(i in 1:nResults)   results[nextRow, i] <<- getLogProbNFL[[i]]$run()
            nextRow <<- nextRow + 1
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
            nextRow <<- 1
            results <<- array(0, c(1, nResults))
        }
    )
)

