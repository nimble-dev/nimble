

#' Performs initialization nimble model node values and log probabilities
#'
#' @param model A setup argument, which specializes an instance of this nimble function to a particular model.
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
initializeModel <- nimbleFunction(
    setup = function(model) {
        RHSonlyNodes <- model$getMaps('nodeNamesRHSonly')
        hasRHSonlyNodes <- length(RHSonlyNodes) > 0

        allDetermNodes <- model$getNodeNames(determOnly = TRUE)
        determNodesNodeFxnVector <- nodeFunctionVector(model = model, nodeNames = allDetermNodes)
        
        stochNonDataNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        
        initFunctionList <- nimbleFunctionList(mcmcNodeInit_virtual)

        iter <- 1
        if(hasRHSonlyNodes) {
            initFunctionList[[iter]] <- mcmcCheckRHS_Init(model = model, node = RHSonlyNodes)
            iter <- iter + 1
        }
        
        for(i in seq_along(stochNonDataNodes))
            initFunctionList[[iter + i - 1]] <- mcmcStochNode_Init(model, stochNonDataNodes[i])
    },
    
    run = function() {
        
        for(i in seq_along(initFunctionList))  {   
            calculate(nodeFxnVector = determNodesNodeFxnVector)
            initFunctionList[[i]]$run()
        }
        calculate(model)

    },  where = getLoadingNamespace()
)


