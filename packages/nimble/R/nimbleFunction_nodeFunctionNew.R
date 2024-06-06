nodeFunctionNew <- function(LHS,
                            RHS,
                            name = NA,
                            altParams,
                            bounds,
                            parentsSizeAndDims,
                            logProbNodeExpr,
                            type,
                            setupOutputExprs,
                            dynamicIndexInfo = NULL,
                            unrolledIndicesMatrix = NULL,
                            evaluate = TRUE,
                            nodeDim = NULL,
                            buildDerivs = FALSE,
                            where = globalenv()) {
    if(!(type %in% c('stoch', 'determ')))       stop(paste0('invalid argument to nodeFunction(): type = ', type))
    setupOutputLabels <- nndf_makeNodeFunctionIndexLabels(setupOutputExprs) ## should perhaps move to the declInfo for preservation
    LHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(LHS, setupOutputLabels)
    RHSrep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(RHS, setupOutputLabels)
    altParamsRep <- lapply(altParams, nndf_replaceSetupOutputsWithIndexedNodeInfo, setupOutputLabels)
    boundsRep <- lapply(bounds, nndf_replaceSetupOutputsWithIndexedNodeInfo, setupOutputLabels)
    logProbNodeExprRep <- nndf_replaceSetupOutputsWithIndexedNodeInfo(logProbNodeExpr, setupOutputLabels)
    
    if(getNimbleOption('allowDynamicIndexing')) {
        if(length(dynamicIndexInfo)) {
            for(i in seq_along(dynamicIndexInfo))
                dynamicIndexInfo[[i]]$indexCode <- nndf_replaceSetupOutputsWithIndexedNodeInfo(dynamicIndexInfo[[i]]$indexCode, setupOutputLabels)
            dynamicIndexLimitsExpr <- nndf_generateDynamicIndexLimitsExpr(dynamicIndexInfo)
        } else dynamicIndexLimitsExpr <- NULL
    } else dynamicIndexLimitsExpr <- NULL
    if(isTRUE(buildDerivs)) {
        parents <- names(parentsSizeAndDims)
        ADconstantsInfo <- makeSizeAndDimListForIndexedInfo(LHSrep, parents, unrolledIndicesMatrix)
        ADconstantsInfo <- makeSizeAndDimListForIndexedInfo(RHSrep, parents, unrolledIndicesMatrix,
                                                            allSizeAndDimList = ADconstantsInfo)
      parentIndexInfoList <- nndf_extractNodeIndices(LHSrep, parents)
      parentIndexInfoList <- nndf_extractNodeIndices(RHSrep, parents, indexExprList = parentIndexInfoList)
      for(i in seq_along(parentIndexInfoList)){
        for(j in seq_along(parentIndexInfoList[[i]])){
          parentsSizeAndDims[[names(parentIndexInfoList)[i]]][[j]]$indexExpr <- parentIndexInfoList[[i]][[j]]$indexExpr
        }
      }
    }
    else{
      ADconstantsInfo <- list()
    }

    nodeFunctionTemplate <-
        substitute(
            nimbleFunction(##contains      = CONTAINS,
                           setup         = SETUPFUNCTION,
                           methods       = METHODS,
                           name          = name,
                           check         = FALSE,
                           buildDerivs  = CALCAD_LIST,
                           where = where)
          ,
            list(##CONTAINS      = nndf_createContains(RHS, type), ## this was used for intermediate classes for get_scale style parameter access, prior to getParam
                 SETUPFUNCTION = nndf_createSetupFunction(buildDerivs, RHS),  ##nndf = new node function
                METHODS       = nndf_createMethodList(LHSrep,
                                                      RHSrep,
                                                      parentsSizeAndDims,
                                                      altParamsRep,
                                                      boundsRep,
                                                      logProbNodeExprRep,
                                                      type,
                                                      dynamicIndexLimitsExpr,
                                                      RHS,
                                                      nodeDim,
                                                      buildDerivs),
                CALCAD_LIST   = if(isTRUE(buildDerivs)) list(calculate_ADproxyModel = list(isNode = TRUE)) else list(),
                where         = where)
        )
    if(evaluate){
      returnFunc <- eval(nodeFunctionTemplate)
      assign('parentsSizeAndDims', parentsSizeAndDims, envir = environment(returnFunc))    
      assign('ADconstantsInfo', ADconstantsInfo, envir = environment(returnFunc))      
      return(returnFunc)
    }
    else       return(nodeFunctionTemplate)
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
nndf_createSetupFunction <- function(buildDerivs = FALSE, RHS) {
    if(isTRUE(buildDerivs)) {
        setup <- function(model, BUGSdecl) {
            indexedNodeInfoTable <- indexedNodeInfoTableClass(BUGSdecl)
            ADproxyModel <- model$ADproxyModel
            setupOutputs(indexedNodeInfoTable)
            invisible(NULL)
        }
    } else {
        setup <- function(model, BUGSdecl) {
            indexedNodeInfoTable <- indexedNodeInfoTableClass(BUGSdecl)
            setupOutputs(indexedNodeInfoTable)
            invisible(NULL)
        }
    }

    if(isTRUE(getNimbleOption('allowNFobjInModel'))) {
      # We will find things like a$foo and insert a line of setup code
      # so that a is created in setup code and thus will be identified by the
      # nimble compiler as a setup output, which later makes it a class member variable in C++.
      # This line of code is "a <- eval(a)". (The eval environment should be correct
      # because the nodeFunction [nimbleFunction] sets that up.)
      find_NFs_recurse <- function(code, result = list()) {
        if(length(code)==1 && !is.call(code)) return(result)
        if(code[[1]]=="$") return(c(result, list(code[[2]])))
        for(i in seq_along(code)) {
          result <- c(result, find_NFs_recurse(code[[i]]))
        }
        result
      }
      new_NFs <- find_NFs_recurse(RHS)
      new_NFs <- unique(new_NFs)
      if(length(new_NFs)) {
        any_nested_dollar_signs <- any(unlist(lapply(new_NFs, function(x) length(x) > 1 && x[[1]]=='$')))
        if(any_nested_dollar_signs) {
          warning("  [warning] when using a nimbleFunction within model code, there can only be one '$' (e.g. a$b, not a$b$c).")
        }
        # One can get a $ in a line of model code from something like
        # a[] <- eigen(B[,])$values
        # so we must ignore anything like a$foo where the a part is itself a call.
        ignore <- unlist(lapply(new_NFs, is.call))
        new_NFs <- new_NFs[!ignore]
        new_lines <- lapply(new_NFs, function(x) substitute(X <- eval(X),
                                                       list(X = x)))
        new_lines <- c(new_lines, quote(invisible(NULL)))
        new_body <- body(setup)
        new_body <- as.call(c(as.list(new_body), new_lines))
        body(setup) <- new_body
      }
    }

    return(setup)
}

indexedNodeInfoTableClass <- function(BUGSdecl) {
    structure(
        list(unrolledIndicesMatrix = BUGSdecl$unrolledIndicesMatrix),
             class = 'indexedNodeInfoTableClass')
}

## creates a list of the methods calculate, simulate, getParam, getBound, and getLogProb, corresponding to LHS, RHS, and type arguments

nndf_createMethodList <- function(LHS,
                                  RHS,
                                  parentsSizeAndDims,
                                  altParams,
                                  bounds,
                                  logProbNodeExpr,
                                  type,
                                  dynamicIndexLimitsExpr,
                                  RHSnonReplaced,
                                  nodeDim,
                                  buildDerivs = FALSE) {
    if(type == 'determ') {
        methodList <- eval(substitute(
            list(
                simulate   = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {  DETERMSIM                                                },
                calculate  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { simulate(INDEXEDNODEINFO_ = INDEXEDNODEINFO_);    returnType(double());   return(invisible(0)) },
                calculateDiff = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {simulate(INDEXEDNODEINFO_ = INDEXEDNODEINFO_);  returnType(double());   return(invisible(0)) },
                getLogProb = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {                returnType(double());   return(0)            }
            ),
            list(LHS=LHS, # no longer used but kept for reference
                 RHS=RHS, # no longer used but kept for reference
                 DETERMSIM = ndf_createDetermSimulate(LHS, RHS, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced, nodeDim = nodeDim)
                 )))
    }
    if(type == 'stoch') {
        methodList <- eval(substitute(
            list(
                simulate   = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { STOCHSIM                                                         },
                calculate  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { STOCHCALC_FULLEXPR;   returnType(double());   return(invisible(LOGPROB)) },
                calculateDiff = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {STOCHCALC_FULLEXPR_DIFF; LocalAns <- LocalNewLogProb - LOGPROB;  LOGPROB <<- LocalNewLogProb;
                                            returnType(double());   return(invisible(LocalAns))},
                getLogProb = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) {                       returnType(double());   return(LOGPROB)            }
            ),
            list(LHS       = LHS, # no longer used but kept for reference
                 LOGPROB   = logProbNodeExpr,
                 STOCHSIM  = ndf_createStochSimulate(LHS, RHS, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced, nodeDim = nodeDim),
                 STOCHCALC_FULLEXPR = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced),
                 STOCHCALC_FULLEXPR_DIFF = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, diff = TRUE, dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced))))
        if(FALSE) {

        if(getNimbleOption('compileAltParamFunctions')) {
            distName <- safeDeparse(RHS[[1]])

            ## add accessor function for node value; used in multivariate conjugate sampler functions
            type = getType(distName)
            nDim <- getDimension(distName)
            methodList[['get_value']] <- ndf_generateGetParamFunction(LHS, type, nDim)
            ## add accessor functions for stochastic node distribution parameters
            for(param in names(RHS[-1])) {
                if(!param %in% c("lower_", "upper_")) {
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
        distName <- safeDeparse(RHS[[1]])

        allParams <- c(list(value = LHS), as.list(RHS[-1]), altParams)

        typesNDims <- getDimension(distName, includeParams = TRUE)
        typesTypes <- getType(distName, includeParams = TRUE)
        paramIDs <- getParamID(distName, includeParams = TRUE)
        ## rely on only double for now
        for(nDimSupported in c(0, 1, 2)) {
            boolThisCase <- typesNDims == nDimSupported ## & typesTypes == 'double' ## until (if ever) we have separate handling of integer params, these should be folded in with doubles.  We don't normally have any integer params, because we handle integers as doubles
            paramNamesToUse <- getParamNames(distName)[boolThisCase]
            caseName <- paste0("getParam_",nDimSupported,"D_double")

            ## Special handling needed for dinterval
            ## Second parameter is a vector but is expected to handle a scalar
            ## It goes in nDimSupported == 1, but we need it's return type cast to vector
            ## if it is in fact only a scalar. We do so by wrapping in c() if
            ## it doesn't any `:` in it
            if(nDimSupported == 1) {
                allParams[paramNamesToUse] <- lapply(allParams[paramNamesToUse],
                                                     function(x) {
                                                         if(':' %in% all.names(x))
                                                             x
                                                         else
                                                             substitute(c(X), list(X = x))
                                                     })
            }
            
            if(length(paramNamesToUse) > 0)
                methodList[[caseName]] <- nndf_generateGetParamSwitchFunction(allParams[paramNamesToUse], paramIDs[paramNamesToUse], type = 'double', nDim = nDimSupported)
        }
        nDimSupported <- 0  # even multivar nodes have single lower and single upper since truncation not supported for mv distributions, so lower and upper come from distribution 'range'
        caseName <- paste0("getBound_",nDimSupported,"D_double")
        methodList[[caseName]] <- nndf_generateGetBoundSwitchFunction(bounds, seq_along(bounds), type = 'double', nDim = nDimSupported)
    }
    ## The next three lines are left over from old approach to AD.  It is not clear if these are still needed.
    parentsArgs <-c()
    ADexceptionNames <- c(names(parentsArgs), deparse(logProbNodeExpr[[2]]))
    exceptionNames <- c("LocalAns", "LocalNewLogProb","PARAMID_","PARAMANSWER_", "BOUNDID_", "BOUNDANSWER_", "INDEXEDNODEINFO_")
    methodList <- nndf_addModelDollarSignsToMethods(methodList, exceptionNames = exceptionNames,
                                                    ADexceptionNames = ADexceptionNames)
    if(isTRUE(buildDerivs)) {
        methodList[['calculate_ADproxyModel']] <- methodList[['calculate']]
        if(type=="stoch") {
          if(!is.null(dynamicIndexLimitsExpr)) {# If we have stochastic indexing, make the AD version not include the error-trapping of lower and upper bounds, because it won't work.
            ## Twist here is to introduce the name AD_return_value_MPQVRBHW_ (which should be unique enough)
            ## so that it can be returned without a second indexing call just to get the returned value.
            tempMethodList <- eval(substitute(
              list(
                calculate  = function(INDEXEDNODEINFO_ = internalType(indexedNodeInfoClass)) { AD_return_value_MPQVRBHW_ <- STOCHCALC_FULLEXPR;   returnType(double());   return(invisible(AD_return_value_MPQVRBHW_)) }
              ),
              list(STOCHCALC_FULLEXPR = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, dynamicIndexLimitsExpr = NULL, RHSnonReplaced = RHSnonReplaced),
                   LOGPROB = logProbNodeExpr)))
            tempMethodList <- nndf_addModelDollarSignsToMethods(tempMethodList, exceptionNames = c(exceptionNames, "AD_return_value_MPQVRBHW_"),
                                                                ADexceptionNames = ADexceptionNames)

            newBody <- body(tempMethodList[['calculate']])
          } else {
            newBody <- body(methodList[['calculate']])
          }
        } else {
            newBody <- body(methodList[['simulate']])
            newBody[[length(newBody) + 1]] <- quote(return(0))
            newBody[[length(newBody) + 1]] <- quote(returnType(double()))
        }
            
        body( methodList[['calculate_ADproxyModel']] ) <-
            eval(substitute(
                substitute(BODY,
                           list(model = as.name("ADproxyModel"))),
                list(BODY = newBody)))
        ## we could also do calculateDiff as follows,
        ##   although derivatives would only be to newly-calculated logProbs, not cached ones.
        ## methodList[['calculateDiff_ADproxyModel']] <- methodList[['calculateDiff']]
        ## body( methodList[['calculateDiff_ADproxyModel']] ) <-
        ##     eval(substitute(
        ##         substitute(BODY,
        ##                    list(model = as.name("ADproxyModel"))),
        ##         list(BODY = body( methodList[['calculateDiff_ADproxyModel']]))))
    }
    return(methodList)
}

nndf_extractNodeIndices <- function(code, nodesToExtract, indexExprList = list()){
  if(is.call(code)){
    if(safeDeparse(code[[1]]) == '[') {
      if(safeDeparse(code[[2]]) %in% nodesToExtract){
        thisIndexExpr <- list()
        for(i in 1:(length(code)-2)){
          if(is.call(code[[i + 2]]) && safeDeparse(code[[i+2]][[1]]) == ':'){
            thisIndexExpr[[i]]   <- list(code[[i+2]][[2]], code[[i+2]][[3]])
          }
          else{
            thisIndexExpr[[i]] <- code[[i+2]]
          }
        }
        if(is.null(indexExprList[[safeDeparse(code[[2]])]])) indexExprList[[safeDeparse(code[[2]])]][[1]]$indexExpr <- thisIndexExpr
        else{
          indexExprList[[safeDeparse(code[[2]])]][[length(indexExprList[[safeDeparse(code[[2]])]]) + 1]] <- list()
          indexExprList[[safeDeparse(code[[2]])]][[length(indexExprList[[safeDeparse(code[[2]])]])]]$indexExpr <- thisIndexExpr
        }
        return(indexExprList)
      }
    }
    if(length(code) > 1){
      for(i in 2:length(code)){
        indexExprList <- nndf_extractNodeIndices(code[[i]], nodesToExtract, indexExprList)
      }
    }
  }
  return(indexExprList)
}

nndf_makeParentSizeExpr <- function(sizeInfoList){
  if(sizeInfoList$nDim == 0) return(parse(text = 'c(1)')[[1]])
  else{
    dimInds <- sizeInfoList$lengths[sizeInfoList$lengths > 1]
    return(parse(text = paste0('c(', paste(dimInds, collapse = ', '), ')'))[[1]])
  }
}

nndf_addModelDollarSignsToMethods <- function(methodList, exceptionNames = character(), ADexceptionNames = character()) {
    for(i in seq_along(methodList)) {
      ## if(names(methodList)[i] == getCalcADFunName()){
      ##   body(methodList[[i]]) <- removeIndices(body(methodList[[i]]))
      ##   body(methodList[[i]]) <- addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames, ADexceptionNames))
      ## }
      ##  else
        body(methodList[[i]]) <- addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames))
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
    if(inherits(ans, 'try-error')) stop('In nndf_generateGetParamSwitchFunction.', call. = FALSE)
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
    if(inherits(ans, 'try-error')) stop('In nndf_generateGetBoundSwitchFunction.', call. = FALSE)
    attr(ans, 'srcref') <- NULL
    ans
}

nndf_generateDynamicIndexLimitsExpr <- function(dynamicIndexInfo) {
    indivLimits <- lapply(dynamicIndexInfo, function(x) substitute(VAREXPR >= LOWER & VAREXPR <= UPPER,
                              list(VAREXPR =  x$indexCode,
                             LOWER = x$lower,
                             UPPER = x$upper)))
    dynamicIndexLimitsExpr <- indivLimits[[1]]
    if(length(indivLimits) > 1)
        for(i in 2:length(indivLimits))
            dynamicIndexLimitsExpr <- substitute(FIRST & SECOND,
                                                 list(FIRST = dynamicIndexLimitsExpr,
                                                      SECOND = indivLimits[[i]]))
    ## FIXME: that puts extra () in the expression; could potentially also construct as
    ## tmp <- quote(1 & 1)
    ## tmp[[2]] <- dynamicIndexLimitsExpr
    ## tmp[[3]] <- indivLimits[[i]]
    ## check if () are stripped out in C++; if so then it doesn't matter anyway
    return(dynamicIndexLimitsExpr)
}

