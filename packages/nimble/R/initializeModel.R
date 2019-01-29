

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
##' @export
initializeModel <- nimbleFunction(
    name = 'initializeModel',
    setup = function(model, silent = FALSE) {
        initFunctionList <- nimbleFunctionList(nodeInit_virtual)
        startInd <- 1

        RHSonlyNodes <- model$getMaps('nodeNamesRHSonly')
        RHSonlyVarNames <- removeIndexing(RHSonlyNodes)
        RHSonlyVarNamesUnique <- unique(RHSonlyVarNames)
        RHSonlyNodesListByVariable <- lapply(RHSonlyVarNamesUnique, function(var) RHSonlyNodes[RHSonlyVarNames==var])
        if(length(RHSonlyNodesListByVariable) > 0) {
            lengths <- sapply(RHSonlyNodesListByVariable, length)
            if(any(lengths == 0)) stop('something went wrong in RHS node model initialization')
            if(sum(lengths) != length(RHSonlyNodes)) stop('something went wrong in RHS node model initialization')
            for(i in seq_along(RHSonlyNodesListByVariable)) {
                initFunctionList[[startInd]] <- checkRHSonlyInit(model = model, nodes = RHSonlyNodesListByVariable[[i]], variable = RHSonlyVarNamesUnique[i])
                startInd <- startInd + 1
            }
        }
        
        LHSnodes <- model$getNodeNames()
        nodeTypes <- model$getNodeType(LHSnodes)
        typeCode <- ifelse(nodeTypes == "determ", 1,
                    ifelse(nodeTypes == "stoch", 2,
                           3)) ## Neither "stoch" nor "determ"
        if(any(typeCode == 3)) message(paste0("found typeCode == 3 for "), paste(LHSnodes[typeCode==3], sep=", "))
        isData <- model$isData(LHSnodes)
    },
    
    run = function() {
        for(i in seq_along(initFunctionList)) {
            initFunctionList[[i]]$run()
        }
        for(i in seq_along(typeCode)) {
            if(typeCode[i] == 1) { ## determ
                initialize_deterministic(i)
            } else {
                if(typeCode[i] == 2) {
                    if(isData[i]) {
                        initialize_stoch_data_node(i)
                    } else {
                        initialize_stoch_non_data_node(i)
                    }                    
                }
            }
        }
        model$calculate(LHSnodes)
    },
    methods = list(
        initialize_deterministic = function(i = integer()) {
            nodeValue <- values(model, LHSnodes[i])
            if(is.na.vec(nodeValue) | is.nan.vec(nodeValue)) {
                model$calculate(LHSnodes[i])
                nodeValue <- values(model, LHSnodes[i])
                if(is.na.vec(nodeValue) | is.nan.vec(nodeValue))
                    print('warning: value of deterministic node ',LHSnodes[i],': value is NA or NaN even after trying to calculate.')
            }
        },
        initialize_stoch_data_node = function(i = integer()) {
            nodeValue <- values(model, LHSnodes[i])
            if(is.na.vec(nodeValue) | is.nan.vec(nodeValue))
                print('warning: value of data node ', LHSnodes[i],': value is NA or NaN.')
            lp <- model$calculate(LHSnodes[i])
            if(is.na(lp) | is.nan(lp)) {
                print('warning: logProb of data node ', LHSnodes[i], ': logProb is NA or NaN.')
            } else {
                if(lp == -Inf) {
                    if(!silent) print('warning: logProb of data node ', LHSnodes[i], ': logProb is -Inf.')
                } else if(lp < -1e12) {
                    if(!silent) print('warning: logProb of data node ', LHSnodes[i], ': logProb less than -1e12.')
                }
            }
        },
        initialize_stoch_non_data_node = function(i = integer()) {
            nodeValue <- values(model, LHSnodes[i])
            if(is.na.vec(nodeValue) | is.nan.vec(nodeValue)) {
                model$simulate(LHSnodes[i], includeData = TRUE) ## includeData = TRUE suppresses a warning
                nodeValue <- values(model, LHSnodes[i])
                if(is.na.vec(nodeValue) | is.nan.vec(nodeValue))
                    print('warning: value of stochastic node ', LHSnodes[i],': value is NA or NaN even after trying to simulate.')
            }
            lp <- model$calculate(LHSnodes[i])
            if(is.na(lp) | is.nan(lp)) {
                print('warning: problem initializing stochastic node ', LHSnodes[i], ': logProb is NA or NaN.')
            } else {
                if(lp == -Inf) {
                    if(!silent) print('warning: problem initializing stochastic node ', LHSnodes[i], ': logProb is -Inf.')
                } else if(lp < -1e12) {
                    if(!silent) print('warning: problem initializing stochastic node ', LHSnodes[i], ': logProb less than -1e12.')
                }
            }
        }
    ),
    where = getLoadingNamespace()
)

nodeInit_virtual <- nimbleFunctionVirtual()

checkRHSonlyInit <- nimbleFunction(
    name = 'checkRHSonlyInit',
    contains = nodeInit_virtual,
    setup = function(model, nodes, variable) {},
    run = function() {
        vals <- values(model, nodes)
        if(is.na.vec(vals) | is.nan.vec(vals)) print('warning: value in right-hand-side-only variable is NA or NaN, in variable: ', variable)
    },    where = getLoadingNamespace()
)

## determNodeInit <- nimbleFunction(
##     name = 'determNodeInit',
##     contains = nodeInit_virtual,
##     setup = function(model, node, silent) {},
##     run = function() {
##         nodeValue <- values(model, node)
##         if(is.na.vec(nodeValue) | is.nan.vec(nodeValue)) calculate(model, node)
##         nodeValue <- values(model, node)
##         if(is.na.vec(nodeValue) | is.nan.vec(nodeValue)) print('warning: value of deterministic node ',node,': value is NA or NaN even after trying to calculate.')
##     },    where = getLoadingNamespace()
## )

## stochDataNodeInit <- nimbleFunction(
##     name = 'stochDataNodeInit',
##     contains = nodeInit_virtual,
##     setup = function(model, node, silent) {},
##     run = function() {
##         nodeValue <- values(model, node)
##         if(is.na.vec(nodeValue) | is.nan.vec(nodeValue)) print('warning: value of data node ',node,': value is NA or NaN.')
##         lp <- calculate(model, node)
##         if(is.na(lp) | is.nan(lp)) print('warning: logProb of data node ', node, ': logProb is NA or NaN.')
##         if(!(is.na(lp) | is.nan(lp))) {
##             if(lp == -Inf) {
##                 if(!silent) print('warning: logProb of data node ', node, ': logProb is -Inf.')
##             } else if(lp < -1e12) {
##                 if(!silent) print('warning: logProb of data node ', node, ': logProb less than -1e12.')
##             }
##         }
##     },    where = getLoadingNamespace()
## )

## stochNonDataNodeInit <- nimbleFunction(
##     name = 'stochNonDataNodeInit',
##     contains = nodeInit_virtual,
##     setup = function(model, node, silent) {},
##     run = function() {
##         nodeValue <- values(model, node)
##         if(is.na.vec(nodeValue) | is.nan.vec(nodeValue)) simulate(model, node)
##         nodeValue <- values(model, node)
##         if(is.na.vec(nodeValue) | is.nan.vec(nodeValue)) print('warning: value of stochastic node ',node,': value is NA or NaN even after trying to simulate.')
##         lp <- calculate(model, node)
##         if(is.na(lp) | is.nan(lp)) print('warning: problem initializing stochastic node ', node, ': logProb is NA or NaN.')
##         if(!(is.na(lp) | is.nan(lp))) {
##             if(lp == -Inf) {
##                 if(!silent) print('warning: problem initializing stochastic node ', node, ': logProb is -Inf.')
##             } else if(lp < -1e12) {
##                 if(!silent) print('warning: problem initializing stochastic node ', node, ': logProb less than -1e12.')
##             }
##         }
##     },    where = getLoadingNamespace()
## )

