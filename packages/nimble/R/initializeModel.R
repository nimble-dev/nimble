

#' Performs initialization of nimble model node values and log probabilities
#'
#' @param model A setup argument, which specializes an instance of this nimble function to a particular model.
#' @param silent logical indicating whether to suppress logging information
#' @author Daniel Turek
#' @details This nimbleFunction may be used at the beginning of nimble algorithms to perform model initialization.
#' The intended usage is to specialize an instance of this nimbleFunction in the setup function of an algorithm,
#' then execute that specialied function at the beginning of the algorithm run function.
#' The specialized function takes no arguments.
#'
#' Executing this function ensures that all right-hand-side only nodes have been assigned real values,
#' that all stochastic nodes have a real value, or otherwise have their simulate() method called,
#' that all deterministic nodes have their simulate() method called,
#' and that all log-probabilities have been calculated with the current model values.
#' An error results if model initialization encounters a problem, for example a missing right-hand-side only
#' node value.
#' @examples
#' myNewAlgorithm <- nimbleFunction(
#'    setup = function(model, ...) {
#'       my_initializeModel <- initializeModel(model)
#'       ....
#'    },
#'    run = function(...) {
#'       my_initializeModel()
#'       ....
#'    }
#' )
#' @export
initializeModel <- nimbleFunction(
    name = 'initializeModel',
    setup = function(model, silent = FALSE) {
        
        ##determDepsOfRHSonly <- setdiff(
        ##    model$getDependencies(model$getMaps('nodeNamesRHSonly'), determOnly = TRUE),
        ##    model$getDependencies(model$getNodeNames(stochOnly = TRUE), determOnly = TRUE))
        
        initFunctionList <- nimbleFunctionList(nodeInit_virtual)
        startInd <- 1

        RHSonlyNodes <- model$getMaps('nodeNamesRHSonly')
        if(length(RHSonlyNodes) > 0) {
            initFunctionList[[startInd]] <- checkRHSonlyInit(model = model, nodes = RHSonlyNodes)
            startInd <- startInd + 1
        }

        topDetermNodes <- model$getNodeNames(topOnly = TRUE, determOnly = TRUE)
        for(i in seq_along(topDetermNodes)) {
            initFunctionList[[startInd]] <- topDetermNodeInit(model = model, node = topDetermNodes[i])
            startInd <- startInd + 1
        }

        stochNonDataNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        for(i in seq_along(stochNonDataNodes)) {
            initFunctionList[[startInd]] <- stochNodeInit(model = model, node = stochNonDataNodes[i], silent = silent)
            startInd <- startInd + 1
        }

    },
    
    run = function() {
        ##model$calculate(determDepsOfRHSonly)
        for(i in seq_along(initFunctionList)) {
            initFunctionList[[i]]$run()
        }
        calculate(model)
    },  where = getLoadingNamespace()
)


nodeInit_virtual <- nimbleFunctionVirtual()

checkRHSonlyInit <- nimbleFunction(
    name = 'checkRHSonlyInit',
    contains = nodeInit_virtual,
    setup = function(model, nodes) {},
    run = function() {
        vals <- values(model, nodes)
        if(is.na.vec(vals) | is.nan.vec(vals)) print('warning: value of right hand side only node not initialized')
    },    where = getLoadingNamespace()
)

topDetermNodeInit <- nimbleFunction(
    name = 'topDetermNodeInit',
    contains = nodeInit_virtual,
    setup = function(model, node) {},
    run = function() {
        theseVals <- values(model, node)
        if(is.na.vec(theseVals) | is.nan.vec(theseVals)) calculate(model, node)
        theseVals <- values(model, node)
        if(is.na.vec(theseVals) | is.nan.vec(theseVals)) print('warning: value of top-level deterministic node ',node,': value is NA or NaN even after trying to calculate.')
    },    where = getLoadingNamespace()
)

stochNodeInit <- nimbleFunction(
    name = 'stochNodeInit',
    contains = nodeInit_virtual,
    setup = function(model, node, silent) {
        thisDetermNodes <- model$getDependencies(node, determOnly=TRUE)
    },
    run = function() {
        theseVals <- values(model, node)
        if(is.na.vec(theseVals) | is.nan.vec(theseVals)) simulate(model, node)
        theseVals <- values(model, node)
        if(is.na.vec(theseVals) | is.nan.vec(theseVals)) print('warning: value of stochastic node ',node,': value is NA or NaN even after trying to simulate.')
        lp <- calculate(model, node)
        if(is.na(lp) | is.nan(lp)) print('warning: problem initializing stochastic node ', node, ': logProb is NA or NaN.')
        if(!(is.na(lp) | is.nan(lp))) {
            if(lp < -1e12) {
                if(!silent) print('warning: problem initializing stochastic node ', node, ': logProb less than -1e12.')
            }
        }
        model$calculate(thisDetermNodes)
    },    where = getLoadingNamespace()
)
