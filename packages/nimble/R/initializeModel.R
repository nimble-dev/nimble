

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
    setup = function(model, silent = FALSE) {
        initFunctionList <- nimbleFunctionList(nodeInit_virtual)
        iter <- 1

        RHSonlyNodes <- model$getMaps('nodeNamesRHSonly')
        if(length(RHSonlyNodes) > 0) {
            initFunctionList[[iter]] <- checkRHSonlyInit(model = model, nodes = RHSonlyNodes)
            iter <- iter + 1
        }
        
        stochNonDataNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        for(i in seq_along(stochNonDataNodes))
            initFunctionList[[iter + i - 1]] <- stochNodeInit(model, stochNonDataNodes[i], silent)

        allDetermNodes <- model$getNodeNames(determOnly = TRUE)
        determNodesNodeFxnVector <- nodeFunctionVector(model = model, nodeNames = allDetermNodes)
    },
    
    run = function() {
        for(i in seq_along(initFunctionList))  {   
            calculate(nodeFxnVector = determNodesNodeFxnVector)
            initFunctionList[[i]]$run()
        }
        calculate(model)
    },  where = getLoadingNamespace()
)


nodeInit_virtual <- nimbleFunctionVirtual()

checkRHSonlyInit <- nimbleFunction(
    contains = nodeInit_virtual,
    setup = function(model, nodes) {},
    run = function() {
        vals <- values(model, nodes)
        if(is.na.vec(vals)) print('warning: value of right hand side only node not initialized')
    },    where = getLoadingNamespace()
)

stochNodeInit <- nimbleFunction(
    contains = nodeInit_virtual,
    setup = function(model, node, silent) {},
    run = function() {
        theseVals <- values(model, node)
        if(is.na.vec(theseVals)) simulate(model, node)
        theseVals <- values(model, node)
        if(is.na.vec(theseVals)) print('warning: value of stochastic node is NA')
        lp <- calculate(model, node)
        if(is.na(lp)) print('warning: problem initializing stochastic node ', node, ', logProb is NA')
        if(!is.na(lp)) {
            if(lp < -1e12) {
                if(!silent) print('warning: problem initializing stochastic node, logProb less than -1e12')
            }
        }
    },    where = getLoadingNamespace()
)
