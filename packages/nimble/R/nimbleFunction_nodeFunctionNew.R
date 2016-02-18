nodeFunctionNew <- function(LHS, RHS, name = NA, altParams, logProbNodeExpr, type, setupOutputExprs, evaluate = TRUE, where = globalenv()) {
    if(!(type %in% c('stoch', 'determ')))       stop(paste0('invalid argument to nodeFunction(): type = ', type))
    browser()
    LHSrep <- replaceSetupOutputsWithIndexedNodeInfo(LHS, setupOutputExprs)
    RHSrep <- replaceSetupOutputsWithIndexedNodeInfo(RHS, setupOutputExprs)
    altParamsRep <- lapply(altParams, replaceSetupOutputsWithIndexedNodeInfo, RHS, setupOutputExprs)
    logProbNodeExprRep <- replaceSetupOutputsWithIndexedNodeInfo(logProbNodeExpr, setupOutputExprs)
    nodeFunctionTemplate <-
        substitute(
            nimbleFunction(##contains      = CONTAINS,
                           setup         = SETUPFUNCTION,
                           methods       = METHODS,
                           name          = name,
                           where = where),
            list(##CONTAINS      = nndf_createContains(RHS, type), ## this was used for intermediate classes for get_scale style parameter access, prior to getParam
                 SETUPFUNCTION = nndf_createSetupFunction(),  ##nndf = new node function
                 METHODS       = nndf_createMethodList(LHSrep, RHSrep, altParamsRep, logProbNodeExprRep, type),
                 where         = where)
        )
    if(evaluate)    return(eval(nodeFunctionTemplate))     else       return(nodeFunctionTemplate)
}

replaceSetupOutputsWithIndexedNodeInfo <- function(code, setupOutputExprs) {
    browser()
    ## use name INDEXEDNODEINFO_
    invisible(NULL)
}

## creates a function object for use as setup argument to nimbleFunction()
nndf_createSetupFunction <- function() {
    setup <- function(indexedNodeInfoEnv) {
        indexedNodeInfoTable <- indexedNodeInfoTableClass(indexedNodeInfoEnv, setupOutputExprs)
        invisible(NULL)
    }
    return(setup)
}

## creates a list of the methods calculate, simulate, and getLogProb, corresponding to LHS, RHS, and type arguments
nndf_createMethodList <- function(LHS, RHS, altParams, logProbNodeExpr, type) {
    if(type == 'determ') {
        methodList <- eval(substitute(
            list(
                simulate   = function(indexedNodeInfo = indexedNodeInfoClass()) { LHS <<- RHS                                                 },
                calculate  = function(indexedNodeInfo = indexedNodeInfoClass()) { simulate(indexedNodeInfo);    returnType(double());   return(invisible(0)) },
                calculateDiff = function(indexedNodeInfo = indexedNodeInfoClass()) {simulate(indexedNodeInfo);  returnType(double());   return(invisible(0)) },
                getLogProb = function(indexedNodeInfo = indexedNodeInfoClass()) {                returnType(double());   return(0)            }
            ),
            list(LHS=LHS, 
                 RHS=RHS)))
    }
    if(type == 'stoch') {
        methodList <- eval(substitute(
            list(
                simulate   = function(indexedNodeInfo = indexedNodeInfoClass()) { LHS <<- STOCHSIM                                                         },
                calculate  = function(indexedNodeInfo = indexedNodeInfoClass()) { STOCHCALC_FULLEXPR;   returnType(double());   return(invisible(LOGPROB)) },
                calculateDiff = function(indexedNodeInfo = indexedNodeInfoClass()) {STOCHCALC_FULLEXPR_DIFF; LocalAns <- LocalNewLogProb - LOGPROB;  LOGPROB <<- LocalNewLogProb;
                                            returnType(double());   return(invisible(LocalAns))},
                getLogProb = function(indexedNodeInfo = indexedNodeInfoClass()) {                       returnType(double());   return(LOGPROB)            }
            ),
            list(LHS       = LHS,
                 LOGPROB   = logProbNodeExpr,
                 STOCHSIM  = ndf_createStochSimulate(RHS),
                 STOCHCALC_FULLEXPR = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS),
                 STOCHCALC_FULLEXPR_DIFF = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, diff = TRUE))))
        if(FALSE) {
        if(nimbleOptions()$compileAltParamFunctions) {
            distName <- as.character(RHS[[1]])
            ## add accessor function for node value; used in multivariate conjugate sampler functions
            typeList <- getDistribution(distName)$types[['value']]
            methodList[['get_value']] <- ndf_generateGetParamFunction(LHS, typeList$type, typeList$nDim)
            ## add accessor functions for stochastic node distribution parameters
            for(param in names(RHS[-1])) {
                if(!param %in% c("lower", "upper")) {
                    typeList <- getDistribution(distName)$types[[param]]
                    methodList[[paste0('get_',param)]] <- ndf_generateGetParamFunction(RHS[[param]], typeList$type, typeList$nDim)
                }
            }
            for(i in seq_along(altParams)) {
                altParamName <- names(altParams)[i]
                typeList <- getDistribution(distName)$types[[altParamName]]
                methodList[[paste0('get_',altParamName)]] <- ndf_generateGetParamFunction(altParams[[altParamName]], typeList$type, typeList$nDim)
            }
        }
        } ## if(FALSE) to cut out old get_XXX param system
            ## new for getParam, eventually to replace get_XXX where XXX is each param name
            ## TO-DO: unfold types and nDims more thoroughly (but all types are implemented as doubles anyway)
            ## understand use of altParams vs. all entries in typesListAllParams
            ## need a value Entry
            allParams <- c(list(value = LHS), as.list(RHS[-1]), altParams)
            typesListAllParams <- getDistribution(distName)$types
            ##numParams <- length(typesListAllParams)
            typesNDims <- unlist(lapply(typesListAllParams, `[[`, 'nDim'))
            typesTypes <- unlist(lapply(typesListAllParams, `[[`, 'type'))
            paramIDs <- getDistribution(distName)$paramIDs
            ## rely on only double for now
            for(nDimSupported in c(0, 1, 2)) {
                boolThisCase <- typesNDims == nDimSupported & typesTypes == 'double'
                paramNamesToUse <- names(typesListAllParams)[boolThisCase]
                caseName <- paste0("getParam_",nDimSupported,"D_double")
                if(length(paramNamesToUse) > 0) 
                    methodList[[caseName]] <- ndf_generateGetParamSwitchFunction(allParams[paramNamesToUse], paramIDs[paramNamesToUse], type = 'double', nDim = nDimSupported) 
            }
        
    }
    ## add model$ in front of all names, except the setupOutputs
    methodList <- nndf_addModelDollarSignsToMethods(methodList, exceptionNames = c("LocalAns", "LocalNewLogProb","PARAMID_","PARAMANSWER_", "INDEXEDNODEINFO_"))
    return(methodList)
}

nndf_addModelDollarSignsToMethods <- function(methodList, exceptionNames = character()) {
    for(i in seq_along(methodList)) {
        body(methodList[[i]]) <-addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames))
    }
    return(methodList)
}
