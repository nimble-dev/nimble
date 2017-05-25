nodeFunctionNew <- function(LHS, RHS, name = NA, altParams, bounds, parents, 
                            logProbNodeExpr, type, setupOutputExprs, evaluate = TRUE, where = globalenv()) {    
    if(!(type %in% c('stoch', 'determ')))       stop(paste0('invalid argument to nodeFunction(): type = ', type))
    setupOutputLabels <- nndf_makeNodeFunctionIndexLabels(setupOutputExprs) ## should perhaps move to the declInfo for preservation
    LHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(LHS, setupOutputLabels)
    RHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(RHS, setupOutputLabels)
    altParamsRep <- lapply(altParams, nndf_replaceSetupOutputsWithIndexedNodeInfo, setupOutputLabels)
    boundsRep <- lapply(bounds, nndf_replaceSetupOutputsWithIndexedNodeInfo, setupOutputLabels)
    logProbNodeExprRep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(logProbNodeExpr, setupOutputLabels)
    nodeFunctionTemplate <-
        substitute(
            nimbleFunction(##contains      = CONTAINS,
                           setup         = SETUPFUNCTION,
                           methods       = METHODS,
                           name          = name,
                           check         = FALSE,
                           enableDerivs  = CALCAD_LIST,
                           where = where),
            list(##CONTAINS      = nndf_createContains(RHS, type), ## this was used for intermediate classes for get_scale style parameter access, prior to getParam
                 SETUPFUNCTION = nndf_createSetupFunction(),  ##nndf = new node function
                 METHODS       = nndf_createMethodList(LHSrep, RHSrep, parents, altParamsRep, boundsRep, logProbNodeExprRep, type),
                 CALCAD_LIST   = (if(type == 'stoch') list(getCalcADFunName()) else list()),
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

## creates a list of the methods calculate, simulate, getParam, getBound, and getLogProb, corresponding to LHS, RHS, and type arguments
nndf_createMethodList <- function(LHS, RHS, parents, altParams, bounds, logProbNodeExpr, type) {
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
                CALCADFUNNAME  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { STOCHCALC_FULLEXPR_AD;   returnType(double());   return(invisible(LOGPROB_AD)) },
                calculateDiff = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {STOCHCALC_FULLEXPR_DIFF; LocalAns <- LocalNewLogProb - LOGPROB;  LOGPROB <<- LocalNewLogProb;
                                            returnType(double());   return(invisible(LocalAns))},
                getLogProb = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {                       returnType(double());   return(LOGPROB)            }
            ),
            list(LHS       = LHS,
                 LOGPROB   = logProbNodeExpr,
                 LOGPROB_AD = as.name('logProb'),
                 STOCHSIM  = ndf_createStochSimulate(RHS),
                 STOCHCALC_FULLEXPR = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS),
                 STOCHCALC_FULLEXPR_AD = ndf_createStochCalculate(as.name('logProb'), LHS, RHS, ADFunc = TRUE),
                 STOCHCALC_FULLEXPR_DIFF = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, diff = TRUE))))
        
        names(methodList)[3] <-  getCalcADFunName() ## replase CALCADFUNNAME with real name
        parentsArgs <- if(length(parents) > 0) list() else NULL
        for(i in seq_along(parents)){
          parentsArgs[[parents[i]]] <- quote(double(1))  
        }
        # parentsArgs <- rep(quote(double), length(parents))
        # if(!is.null(parents)){
        #   names(parentsArgs) <- parents
        # }
        selfWithNoInds <- c(quote(double(1)))
        names(selfWithNoInds) <-  strsplit(deparse(LHS), '[', fixed = TRUE)[[1]][1]
        formals(methodList[[getCalcADFunName()]]) <- c( formals(methodList[[getCalcADFunName()]]), 
                                                        selfWithNoInds, parentsArgs)
        if(FALSE) {
        if(nimbleOptions()$compileAltParamFunctions) {
            distName <- as.character(RHS[[1]])
            ## add accessor function for node value; used in multivariate conjugate sampler functions
            type = getType(distName)
            nDim <- getDimension(distName)
            methodList[['get_value']] <- ndf_generateGetParamFunction(LHS, type, nDim)
            ## add accessor functions for stochastic node distribution parameters
            for(param in names(RHS[-1])) {
                if(!param %in% c("lower", "upper")) {
                    type = getType(distName, param)
                    nDim <- getDimension(distName, param)
                    methodList[[paste0('get_',param)]] <- ndf_generateGetParamFunction(RHS[[param]], type, nDim)
                }
            }
            for(i in seq_along(altParams)) {
                altParamName <- names(altParams)[i]
                type = getType(distName, altParamName)
                nDim <- getDimension(distName, altParamName)
                methodList[[paste0('get_',altParamName)]] <- ndf_generateGetParamFunction(altParams[[altParamName]], type, nDim)
            }
        }
        } ## if(FALSE) to cut out old get_XXX param system
            ## new for getParam, eventually to replace get_XXX where XXX is each param name
            ## TO-DO: unfold types and nDims more thoroughly (but all types are implemented as doubles anyway)
            ## understand use of altParams vs. all entries in typesListAllParams (i.e., getDistributionInfo(distName)$types
        ## need a value Entry
        distName <- as.character(RHS[[1]])

        allParams <- c(list(value = LHS), as.list(RHS[-1]), altParams)

        typesNDims <- getDimension(distName, includeParams = TRUE)
        typesTypes <- getType(distName, includeParams = TRUE)
        paramIDs <- getParamID(distName, includeParams = TRUE)
        ## rely on only double for now
        for(nDimSupported in c(0, 1, 2)) {
            boolThisCase <- typesNDims == nDimSupported ## & typesTypes == 'double' ## until (if ever) we have separate handling of integer params, these should be folded in with doubles.  We don't normally have any integer params, because we handle integers as doubles
            paramNamesToUse <- getParamNames(distName)[boolThisCase]
            caseName <- paste0("getParam_",nDimSupported,"D_double")
            if(length(paramNamesToUse) > 0)
                methodList[[caseName]] <- nndf_generateGetParamSwitchFunction(allParams[paramNamesToUse], paramIDs[paramNamesToUse], type = 'double', nDim = nDimSupported)
        }
        nDimSupported <- 0  # even multivar nodes have single lower and single upper since truncation not supported for mv distributions, so lower and upper come from distribution 'range'
        caseName <- paste0("getBound_",nDimSupported,"D_double")
        methodList[[caseName]] <- nndf_generateGetBoundSwitchFunction(bounds, seq_along(bounds), type = 'double', nDim = nDimSupported)
    }
    ## add model$ in front of all names, except the setupOutputs
    methodList <- nndf_addModelDollarSignsToMethods(methodList, exceptionNames = c("LocalAns", "LocalNewLogProb","PARAMID_","PARAMANSWER_", "BOUNDID_", "BOUNDANSWER_", "INDEXEDNODEINFO_"), 
                                                    ADexceptionNames = c(names(parentsArgs), names(selfWithNoInds), 'logProb'))
    return(methodList)
}

nndf_addModelDollarSignsToMethods <- function(methodList, exceptionNames = character(), ADexceptionNames = character()) {
    for(i in seq_along(methodList)) {
      if(names(methodList)[i] == getCalcADFunName()){
        # body(methodList[[i]]) <- removeIndices(body(methodList[[i]]))
        body(methodList[[i]]) <- addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames, ADexceptionNames))
      }
      else  body(methodList[[i]]) <- addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames))
    }
    return(methodList)
}

nndf_generateGetParamSwitchFunction <- function(typesListAll, paramIDs, type, nDim) {
    if(any(unlist(lapply(typesListAll, is.null)))) stop(paste('problem creating switch function for getParam from ', paste(paste(names(typesListAll), as.character(typesListAll), sep='='), collapse=',')))
    paramIDs <- as.integer(paramIDs)
    answerAssignmentExpressions <- lapply(typesListAll, function(x) substitute(PARAMANSWER_ <- ANSEXPR, list(ANSEXPR = x)))
    switchCode <- as.call(c(list(quote(nimSwitch), quote(PARAMID_), paramIDs), answerAssignmentExpressions))
    # avoid arg name mismatch based on R partial arg name matching
    names(answerAssignmentExpressions) <- NULL
    names(switchCode)[2:3] <- c('paramID', 'IDoptions')
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

nndf_generateGetBoundSwitchFunction <- function(typesListAll, boundIDs, type, nDim) {
    if(any(unlist(lapply(typesListAll, is.null)))) stop(paste('problem creating switch function for getBound from ', paste(paste(names(typesListAll), as.character(typesListAll), sep='='), collapse=',')))
    boundIDs <- as.integer(boundIDs)
    answerAssignmentExpressions <- lapply(typesListAll, function(x) substitute(BOUNDANSWER_ <- ANSEXPR, list(ANSEXPR = x)))
    switchCode <- as.call(c(list(quote(nimSwitch), quote(BOUNDID_), boundIDs), answerAssignmentExpressions))
    if(nDim == 0) {
        answerInitCode <- quote(BOUNDANSWER_ <- 0)  ## this avoids a Windows compiler warning about a possibly unassigned return variable

        ans <- try(eval(substitute(
            function(BOUNDID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                ANSWERINITCODE
                SWITCHCODE
                return(BOUNDANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, ANSWERINITCODE = answerInitCode, SWITCHCODE = switchCode)
        )))
    } else {
        ans <- try(eval(substitute(
            function(BOUNDID_ = integer(), INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {
                returnType(TYPE(NDIM))
                SWITCHCODE
                return(BOUNDANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, SWITCHCODE = switchCode)
        )))
    }
    if(inherits(ans, 'try-error')) browser()
    attr(ans, 'srcref') <- NULL
    ans
}
