nodeFunctionNew <- function(LHS, RHS, name = NA, altParams, logProbNodeExpr, type, setupOutputExprs, evaluate = TRUE, where = globalenv()) {
    if(!(type %in% c('stoch', 'determ')))       stop(paste0('invalid argument to nodeFunction(): type = ', type))
    setupOutputLabels <- nndf_makeNodeFunctionIndexLabels(setupOutputExprs) ## should perhaps move to the declInfo for preservation
    LHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(LHS, setupOutputLabels)
    RHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(RHS, setupOutputLabels)
    altParamsRep <- lapply(altParams, nndf_replaceSetupOutputsWithIndexedNodeInfo, setupOutputLabels)
    logProbNodeExprRep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(logProbNodeExpr, setupOutputLabels)
    nodeFunctionTemplate <-
        substitute(
            nimbleFunction(##contains      = CONTAINS,
                           setup         = SETUPFUNCTION,
                           methods       = METHODS,
                           name          = name,
                           check         = FALSE,
                           where = where),
            list(##CONTAINS      = nndf_createContains(RHS, type), ## this was used for intermediate classes for get_scale style parameter access, prior to getParam
                 SETUPFUNCTION = nndf_createSetupFunction(),  ##nndf = new node function
                 METHODS       = nndf_createMethodList(LHSrep, RHSrep, altParamsRep, logProbNodeExprRep, type),
                 where         = where)
        )
    if(evaluate)    return(eval(nodeFunctionTemplate))     else       return(nodeFunctionTemplate)
}

nndf_makeNodeFunctionIndexLabels <- function(setupOutputExprs) {
    setupOutputLabels <- 1:length(setupOutputExprs)
    names(setupOutputLabels) <- names(setupOutputExprs)
    setupOutputLabels
}

nndf_makeNodeFunctionIndexAccessCall <- function(index) {
    ## use name INDEXEDNODEINFO_
    substitute(getNodeFunctionIndexedInfo(INDEXEDNODEINFO_, INDEX), list(INDEX = index)) ## still needs unity decrement for c++
}

nndf_replaceSetupOutputsWithIndexedNodeInfo <- function(code, setupOutputLabels) {
    cLength <- length(code)
    if(cLength == 1) {
        if(is.name(code)) {
            varName <- as.character(code)
            if(varName %in% names(setupOutputLabels))
                return(nndf_makeNodeFunctionIndexAccessCall(setupOutputLabels[[varName]]))
        }
        return(code)
    }
    if(is.call(code)) {
        for(i in 2:cLength)
            code[[i]] <- nndf_replaceSetupOutputsWithIndexedNodeInfo(code[[i]], setupOutputLabels)
        return(code)
    }
    return(code)
}

# need function to be defined to pass CRAN but setupOutputs is
# never called - it is processed out of nimbleFunction setup code
setupOutputs <- function(...) NULL

## creates a function object for use as setup argument to nimbleFunction()
nndf_createSetupFunction <- function() {
    setup <- function(model, BUGSdecl) {
        indexedNodeInfoTable <- indexedNodeInfoTableClass(BUGSdecl)
        setupOutputs(indexedNodeInfoTable)
        invisible(NULL)
    }
    return(setup)
}

indexedNodeInfoTableClass <- function(BUGSdecl) {
    structure(
        list(unrolledIndicesMatrix = BUGSdecl$unrolledIndicesMatrix),
             class = 'indexedNodeInfoTableClass')
}

## creates a list of the methods calculate, simulate, and getLogProb, corresponding to LHS, RHS, and type arguments
nndf_createMethodList <- function(LHS, RHS, altParams, logProbNodeExpr, type) {
    if(type == 'determ') {
        methodList <- eval(substitute(
            list(
                simulate   = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { LHS <<- RHS                                                 },
                calculate  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { simulate(INDEXEDNODEINFO_ = INDEXEDNODEINFO_);    returnType(double());   return(invisible(0)) },
                calculateDiff = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {simulate(INDEXEDNODEINFO_ = INDEXEDNODEINFO_);  returnType(double());   return(invisible(0)) },
                getLogProb = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {                returnType(double());   return(0)            }
            ),
            list(LHS=LHS,
                 RHS=RHS)))
    }
    if(type == 'stoch') {
        methodList <- eval(substitute(
            list(
                simulate   = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { LHS <<- STOCHSIM                                                         },
                calculate  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { STOCHCALC_FULLEXPR;   returnType(double());   return(invisible(LOGPROB)) },
                calculateDiff = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {STOCHCALC_FULLEXPR_DIFF; LocalAns <- LocalNewLogProb - LOGPROB;  LOGPROB <<- LocalNewLogProb;
                                            returnType(double());   return(invisible(LocalAns))},
                getLogProb = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {                       returnType(double());   return(LOGPROB)            }
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
        distName <- as.character(RHS[[1]])

        allParams <- c(list(value = LHS), as.list(RHS[-1]), altParams)
        typesListAllParams <- getDistribution(distName)$types
        ##numParams <- length(typesListAllParams)
        typesNDims <- unlist(lapply(typesListAllParams, `[[`, 'nDim'))
        typesTypes <- unlist(lapply(typesListAllParams, `[[`, 'type'))
        paramIDs <- getDistribution(distName)$paramIDs
        ## rely on only double for now
        for(nDimSupported in c(0, 1, 2)) {
            boolThisCase <- typesNDims == nDimSupported ## & typesTypes == 'double' ## until (if ever) we have separate handling of integer params, these should be folded in with doubles.  We don't normally have any integer params, because we handle integers as doubles
            paramNamesToUse <- names(typesListAllParams)[boolThisCase]
            caseName <- paste0("getParam_",nDimSupported,"D_double")
            if(length(paramNamesToUse) > 0)
                methodList[[caseName]] <- nndf_generateGetParamSwitchFunction(allParams[paramNamesToUse], paramIDs[paramNamesToUse], type = 'double', nDim = nDimSupported)
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

nndf_generateGetParamSwitchFunction <- function(typesListAll, paramIDs, type, nDim) {
    if(any(unlist(lapply(typesListAll, is.null)))) stop(paste('problem creating switch function for getParam from ', paste(paste(names(typesListAll), as.character(typesListAll), sep='='), collapse=',')))
    paramIDs <- as.integer(paramIDs)
    answerAssignmentExpressions <- lapply(typesListAll, function(x) substitute(PARAMANSWER_ <- ANSEXPR, list(ANSEXPR = x)))
    switchCode <- as.call(c(list(quote(nimSwitch), quote(PARAMID_), paramIDs), answerAssignmentExpressions))
    if(nDim == 0) {
        answerInitCode <- quote(PARAMANSWER_ <- 0)  ## this avoids a Windows compiler warning about a possibly unassigned return variable

        ans <- try(eval(substitute(
            function(PARAMID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                ANSWERINITCODE
                SWITCHCODE
                return(PARAMANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, ANSWERINITCODE = answerInitCode, SWITCHCODE = switchCode)
        )))
    } else {
        ans <- try(eval(substitute(
            function(PARAMID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                SWITCHCODE
                return(PARAMANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, SWITCHCODE = switchCode)
        )))
    }
    if(inherits(ans, 'try-error')) browser()
    attr(ans, 'srcref') <- NULL
    ans
}

