

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
        initFunctionList <- nimbleFunctionList(nodeInit_virtual)
        iter <- 1

        RHSonlyNodes <- model$getMaps('nodeNamesRHSonly')
        if(length(RHSonlyNodes) > 0) {
            initFunctionList[[iter]] <- checkRSHonlyInit(model = model, nodes = RHSonlyNodes)
            iter <- iter + 1
        }
        
        stochNonDataNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        for(i in seq_along(stochNonDataNodes))
            initFunctionList[[iter + i - 1]] <- stochNodeInit(model, stochNonDataNodes[i])

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

checkRSHonlyInit <- nimbleFunction(
    contains = nodeInit_virtual,
    setup = function(model, nodes) {},
    run = function() {
        vals <- values(model, nodes)
        if(is.na.vec(vals)) print('Value of right hand side only node not initialized')
    },    where = getLoadingNamespace()
)

stochNodeInit <- nimbleFunction(
    contains = nodeInit_virtual,
    setup = function(model, node) {},
    run = function() {
        theseVals <- values(model, node)
        if(is.na.vec(theseVals)) simulate(model, node)
        lp <- calculate(model, node)
        if(is.na(lp)) print('Problem initializing stochastic node, logProb is NA')
        if(lp < -1e12) print('Problem initializing stochastic node, logProb less than -1e12')
    },    where = getLoadingNamespace()
)







## Commenting out, these seem to not be used any longer.
## -DT May 2015
##
## RHSonlyInit_virtual <- nimbleFunctionVirtual()
##
## RHSonlyInit <- nimbleFunction(
##     contains = RHSonlyInit_virtual,
##     setup = function(model, node) {},
##     run = function() {
##         nv <- values(model, node)
##         if(is.na.vec(nv) | is.nan.vec(nv))     print('missing value in right-hand-side only node; cannot initialize model')
##     }, where = getLoadingNamespace()
## )
##
## nodeInit <- nimbleFunction(
##     contains = nodeInit_virtual,
##     setup = function(model, node) {
##         gID <- model$modelDef$nodeName2GraphIDs(node)
##         type <- model$modelDef$maps$types[gID]
##         isDeterm = FALSE
##         isStoch = FALSE
##         if(type == 'stoch')
##             isStoch = TRUE
##         else if(type == 'determ')
##             isDeterm = TRUE
##     },
##     run = function() {
##         if(isDeterm) {
##             calculate(model, node)
##             nv <- values(model, node)
##             if(is.na.vec(nv) | is.nan.vec(nv))     print('deterministic model node is NA or NaN in model initialization')
##         }
##         if(isStoch) {
##             nv <- values(model, node)
##             if(is.na.vec(nv)) {
##                 simulate(model, node)
##                 nv <- values(model, node)
##             }
##             if(is.na.vec(nv) | is.nan.vec(nv))     print('stochastic model node is NA or NaN in model initialization')
##             lp <- calculate(model, node)
##             if(is.na(lp) | is.nan(lp) | lp < -1e12)              print('stochastic model value is NA, NaN or too small in model initialization')
##         }
##     }, where = getLoadingNamespace()
## )
##
## fillDetermTop_Init <- nimbleFunction(
##     contains = nodeInit_virtual,
##     setup = function(model, node){},
##     run = function(){
##         nil <- calculate(model, node)
##     }, where = getLoadingNamespace()
## )




