
## helper function that adds an argument to a call
## used to add needed arguments for C versions of {d,p,q}${dist} functions
addArg <- function(code, value, name) {
    newArgIndex <- length(code) + 1
    code[newArgIndex] <- value
    names(code)[newArgIndex] <-  name
    return(code)
}


ndf_createDetermSimulate <- function(LHS, RHS, dynamicIndexLimitsExpr, RHSnonReplaced, nodeDim) {
    if(nimbleOptions()$allowDynamicIndexing && !is.null(dynamicIndexLimitsExpr)) {
        ## error gracefully if dynamic index too small or large; we don't catch non-integers within the bounds though
        if(is.null(nodeDim)) {
            nanExpr <- NaN
        } else nanExpr <- substitute(nimRep(NaN, LENGTH),
                                     list(LENGTH = prod(nodeDim)))
        code <- substitute(if(isTRUE(CONDITION)) LHS <<- RHS
                           else {
                               LHS <<- NANEXPR
                               print(TEXT) },  
                           list(LHS = LHS,
                                RHS = RHS,
                                NANEXPR = nanExpr,
                                CONDITION = dynamicIndexLimitsExpr,
                                TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))))
    } else code <- substitute(LHS <<- RHS,
                              list(LHS = LHS,
                                   RHS = RHS))
    return(code)
}

## changes 'dnorm(mean=1, sd=2)' into 'rnorm(1, mean=1, sd=2)'
ndf_createStochSimulate <- function(LHS, RHS, dynamicIndexLimitsExpr, RHSnonReplaced, nodeDim) {
    BUGSdistName <- as.character(RHS[[1]])
    RHS[[1]] <- as.name(getDistributionInfo(BUGSdistName)$simulateName)   # does the appropriate substituion of the distribution name
    if(length(RHS) > 1) {    for(i in (length(RHS)+1):3)   { RHS[i] <- RHS[i-1];     names(RHS)[i] <- names(RHS)[i-1] } }    # scoots all named arguments right 1 position
    RHS[[2]] <- 1;     names(RHS)[2] <- ''    # adds the first (unnamed) argument '1'    
    if("lower_" %in% names(RHS) || "upper_" %in% names(RHS)) {
        RHS <- ndf_createStochSimulateTrunc(RHS, discrete = getAllDistributionsInfo('discrete')[BUGSdistName])
    } 
    if(nimbleOptions()$allowDynamicIndexing && !is.null(dynamicIndexLimitsExpr)) {
        if(is.null(nodeDim)) {
            nanExpr <- NaN
        } else nanExpr <- substitute(rep(NaN, LENGTH),
                                     list(LENGTH = prod(nodeDim)))
        code <- substitute(if(isTRUE(CONDITION)) LHS <<- RHS
                           else {
                               LHS <<- NANEXPR
                               print(TEXT) },  
                           list(LHS = LHS,
                                RHS = RHS,
                                NANEXPR = nanExpr,
                                CONDITION = dynamicIndexLimitsExpr,
                                TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))))
    } else {
        code <- substitute(LHS <<- RHS,
                           list(LHS = LHS,
                                RHS = RHS))
    }
    return(code)
}


## changes 'rnorm(mean=1, sd=2, lower_=0, upper_=3)' into correct truncated simulation
##   using inverse CDF
ndf_createStochSimulateTrunc <- function(RHS, discrete = FALSE) {
    lowerPosn <- which("lower_" == names(RHS))
    upperPosn <- which("upper_" == names(RHS))
    lower <- RHS[[lowerPosn]]
    upper <- RHS[[upperPosn]]
    RHS <- RHS[-c(lowerPosn, upperPosn)]
    dist <- substring(as.character(RHS[[1]]), 2, 1000)
    
    lowerTailName <- 'lower.tail' 
    logpName <- 'log.p' 
    logName <- 'log' 
    # setup for runif(1, pdist(lower,...), pdist(upper,...))
    # pdist() expression template for inputs to runif()
    pdistTemplate <- RHS
    pdistTemplate[[1]] <- as.name(paste0("p", dist))
    pdistTemplate <- addArg(pdistTemplate, 1, lowerTailName)
    pdistTemplate <- addArg(pdistTemplate, 0, logpName)

    if(discrete) {
        ddistTemplate <- RHS
        ddistTemplate[[1]] <- as.name(paste0("d", dist))
        ddistTemplate <- addArg(ddistTemplate, 0, logName)
        ceilTemplate <- quote(ceiling(x))
    } else ddistTemplate <- NULL
    
    # create bounds for runif() using pdist expressions
    MIN_EXPR <- 0
    MAX_EXPR <- 1
    VALUE_EXPR <- 0
    if(lower != -Inf) {
        pdistTemplate[[2]] <- lower
        if(discrete) {
            ceilTemplate[[2]] <- lower
            ddistTemplate[[2]] <- ceilTemplate
            pdistTemplate[[2]] <- ceilTemplate
        }
        VALUE_EXPR <- ddistTemplate
        MIN_EXPR <- pdistTemplate
    } 
    if(upper != Inf) {
        pdistTemplate[[2]] <- upper
        MAX_EXPR <- pdistTemplate
    }
    
    # now create full runif() expression
    if(discrete && lower != -Inf) {
        substCode <- quote(runif(1, MIN - VALUE, MAX))
    } else {
          substCode <- quote(runif(1, MIN, MAX))
      }
    RUNIF_EXPR <- eval(substitute(substitute(e, list(
        MIN = MIN_EXPR,
        MAX = MAX_EXPR,
        VALUE = VALUE_EXPR)), list(e = substCode)))

    # create full qdist(runif(...),...) expression
    RHS[[1]] <- as.name(paste0("q", dist))
    RHS[[2]] <- RUNIF_EXPR
    RHS <- addArg(RHS, 1, lowerTailName)
    RHS <- addArg(RHS, 0, logpName)

    return(RHS)
}

## changes 'dnorm(mean=1, sd=2)' into 'dnorm(LHS, mean=1, sd=2, log=TRUE)'
ndf_createStochCalculate <- function(logProbNodeExpr, LHS, RHS, diff = FALSE, ADFunc = FALSE, dynamicIndexLimitsExpr, RHSnonReplaced) {
    BUGSdistName <- as.character(RHS[[1]])
    RHS[[1]] <- as.name(getDistributionInfo(BUGSdistName)$densityName)   # does the appropriate substituion of the distribution name
    if(length(RHS) > 1) {    for(i in (length(RHS)+1):3)   { RHS[i] <- RHS[i-1];     names(RHS)[i] <- names(RHS)[i-1] } }    # scoots all named arguments right 1 position
    RHS[[2]] <- LHS;     names(RHS)[2] <- ''    # adds the first (unnamed) argument LHS
    if("lower_" %in% names(RHS) || "upper_" %in% names(RHS)) {
        return(ndf_createStochCalculateTrunc(logProbNodeExpr, LHS, RHS, diff = diff, discrete = getAllDistributionsInfo('discrete')[BUGSdistName], dynamicIndexLimitsExpr = dynamicIndexLimitsExpr, RHSnonReplaced = RHSnonReplaced))
    } else {
        userDist <- BUGSdistName %in% getAllDistributionsInfo('namesVector', userOnly = TRUE)
        RHS <- addArg(RHS, 1, 'log')  # adds the last argument log=TRUE (log_value for user-defined) # This was changed to 1 from TRUE for easier C++ generation
        if(nimbleOptions()$allowDynamicIndexing && !is.null(dynamicIndexLimitsExpr)) {
            if(diff) {
                code <- substitute(if(isTRUE(CONDITION)) LocalNewLogProb <- STOCHCALC
                                   else {LocalNewLogProb <- NaN
                                       print(TEXT)}, # LocalNewLogProb <- -Inf,
                                   list(STOCHCALC = RHS,
                                        CONDITION = dynamicIndexLimitsExpr,
                                        TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))))
            } else if(ADFunc){  ## don't want global assignment for _AD_ functions.
                code <- substitute(if(isTRUE(CONDITION)) LOGPROB <- STOCHCALC
                                   else {LOGPROB <- NaN
                                       print(TEXT)},
                                   list(LOGPROB = logProbNodeExpr,
                                        STOCHCALC = RHS,
                                        CONDITION = dynamicIndexLimitsExpr,
                                        TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))))
            } else {

                code <- substitute(if(isTRUE(CONDITION)) LOGPROB <<- STOCHCALC
                                   else {LOGPROB <<- NaN
                                       print(TEXT)}, # LOGPROB <<- -Inf,
                                   list(LOGPROB = logProbNodeExpr,
                                        STOCHCALC = RHS,
                                        CONDITION = dynamicIndexLimitsExpr,
                                        TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))))
            }
        } else {
            if(diff) {
                code <- substitute(LocalNewLogProb <- STOCHCALC,
                                   list(
                                       STOCHCALC = RHS))
            } else if(ADFunc){  ## don't want global assignment for _AD_ functions.
                code <- substitute(LOGPROB <- STOCHCALC,
                                   list(LOGPROB = logProbNodeExpr,
                                        STOCHCALC = RHS))
            } else{
                code <- substitute(LOGPROB <<- STOCHCALC,
                                   list(LOGPROB = logProbNodeExpr,
                                        STOCHCALC = RHS))
            }
        }
        return(code)
    }
}

## changes 'dnorm(mean=1, sd=2, lower=0, upper=3)' into correct truncated calculation
ndf_createStochCalculateTrunc <- function(logProbNodeExpr, LHS, RHS, diff = FALSE, discrete = FALSE, dynamicIndexLimitsExpr, RHSnonReplaced) {
    lowerPosn <- which("lower_" == names(RHS))
    upperPosn <- which("upper_" == names(RHS))
    lower <- RHS[[lowerPosn]]
    upper <- RHS[[upperPosn]]
    RHS <- RHS[-c(lowerPosn, upperPosn)]
    dist <- substring(as.character(RHS[[1]]), 2, 1000)
    lowerTailName <- 'lower.tail' 
    logpName <- 'log.p' 
    logName <- 'log' 

    pdistTemplate <- RHS
    pdistTemplate[[1]] <- as.name(paste0("p", dist))
    pdistTemplate <- addArg(pdistTemplate, 1, lowerTailName)
    pdistTemplate <- addArg(pdistTemplate, 0, logpName)

    if(discrete) {
        ddistTemplate <- RHS
        ddistTemplate[[1]] <- as.name(paste0("d", dist))
        ddistTemplate <- addArg(ddistTemplate, 0, logName)
        ceilTemplate <- quote(ceiling(x))
    } else ddistTemplate <- NULL
    
    PDIST_LOWER <- 0
    PDIST_UPPER <- 1
    DDIST_LOWER <- 0
    if(lower != -Inf) {
        pdistTemplate[[2]] <- lower
        if(discrete) {
            ceilTemplate[[2]] <- lower
            ddistTemplate[[2]] <- ceilTemplate
            pdistTemplate[[2]] <- ceilTemplate
        }
        PDIST_LOWER <- pdistTemplate
        DDIST_LOWER <- ddistTemplate
    } 
    if(upper != Inf) {
        pdistTemplate[[2]] <- upper
        PDIST_UPPER <- pdistTemplate
    }

    RHS <- addArg(RHS, 1, logName)  # add log=1 now that pdist() created without 'log'

                code <- substitute(if(isTRUE(CONDITION)) LocalNewLogProb <- STOCHCALC
                                   else { LocalNewLogProb <- NaN
                                       print(TEXT)}, # LocalNewLogProb <- -Inf,
                                   list(STOCHCALC = RHS,
                                        CONDITION = dynamicIndexLimitsExpr,
                                        TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))))

    
    if(nimbleOptions()$allowDynamicIndexing && !is.null(dynamicIndexLimitsExpr)) {
        if(diff) {
            if(discrete && lower != -Inf) {
                substCode <- quote(if(isTRUE(CONDITION)) {
                                       if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                           LocalNewLogProb <- DENSITY - log(PDIST_UPPER - PDIST_LOWER + DDIST_LOWER)
                                       else LocalNewLogProb <- -Inf
                                   } else {
                                       LocalNewLogProb <- NaN
                                       cat(TEXT) }
                                   )
            } else {
                substCode <- quote(if(isTRUE(CONDITION)) {
                                       if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                           LocalNewLogProb <- DENSITY - log(PDIST_UPPER - PDIST_LOWER)
                                       else LocalNewLogProb <- -Inf
                                   } else {
                                       LocalNewLogProb <- NaN
                                       print(TEXT) }
                                           )
            }
            code <- eval(substitute(substitute(e, 
                                               list(
                                                   LOWER = lower,
                                                   UPPER = upper,
                                                   VALUE = LHS,
                                                   DENSITY = RHS,
                                                   PDIST_LOWER = PDIST_LOWER,
                                                   PDIST_UPPER = PDIST_UPPER,
                                                   DDIST_LOWER = DDIST_LOWER,
                                                   CONDITION = dynamicIndexLimitsExpr,
                                                   TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))
                                               )), list( e = substCode)))
        } else {
            if(discrete && lower != -Inf) {
                substCode <- quote(if(isTRUE(CONDITION)) {
                                       if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                           LOGPROB <<- DENSITY - log(PDIST_UPPER - PDIST_LOWER + DDIST_LOWER)
                                       else LOGPROB <<- -Inf
                                   } else {
                                       LOGPROB <<- NaN
                                       cat(TEXT) }
                                   )
            } else {
                substCode <- quote(if(isTRUE(CONDITION)) {
                                       if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                         LOGPROB <<- DENSITY - log(PDIST_UPPER - PDIST_LOWER)
                                   else LOGPROB <<- -Inf
                                   } else {
                                       LOGPROB <<- NaN
                                       print(TEXT) }
                                   )
            }
            code <- eval(substitute(substitute(e, 
                                               list(
                                                   LOWER = lower,
                                                   UPPER = upper,
                                                   VALUE = LHS,
                                                   LOGPROB = logProbNodeExpr,
                                                   DENSITY = RHS,
                                                   PDIST_LOWER = PDIST_LOWER,
                                                   PDIST_UPPER = PDIST_UPPER,
                                                   DDIST_LOWER = DDIST_LOWER,
                                                   CONDITION = dynamicIndexLimitsExpr,
                                                   TEXT = paste0("  [Warning] Dynamic index out of bounds: ", safeDeparse(RHSnonReplaced))
                                               )), list(e = substCode)))
        }
    } else {
        if(diff) {
            if(discrete && lower != -Inf) {
                substCode <- quote(if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                       LocalNewLogProb <- DENSITY - log(PDIST_UPPER - PDIST_LOWER + DDIST_LOWER)
                                   else LocalNewLogProb <- -Inf)
            } else {
                substCode <- quote(if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                       LocalNewLogProb <- DENSITY - log(PDIST_UPPER - PDIST_LOWER)
                                   else LocalNewLogProb <- -Inf)
            }
            code <- eval(substitute(substitute(e, 
                                               list(
                                                   LOWER = lower,
                                                   UPPER = upper,
                                                   VALUE = LHS,
                                                   DENSITY = RHS,
                                                   PDIST_LOWER = PDIST_LOWER,
                                                   PDIST_UPPER = PDIST_UPPER,
                                                   DDIST_LOWER = DDIST_LOWER
                                               )), list( e = substCode)))
        } else {
            if(discrete && lower != -Inf) {
                substCode <- quote(if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                       LOGPROB <<- DENSITY - log(PDIST_UPPER - PDIST_LOWER + DDIST_LOWER)
                                   else LOGPROB <<- -Inf)
            } else {
                substCode <- quote(if(isTRUE(LOWER <= VALUE & VALUE <= UPPER))
                                       LOGPROB <<- DENSITY - log(PDIST_UPPER - PDIST_LOWER)
                                   else LOGPROB <<- -Inf)
            }
            code <- eval(substitute(substitute(e, 
                                               list(
                                                   LOWER = lower,
                                                   UPPER = upper,
                                                   VALUE = LHS,
                                                   LOGPROB = logProbNodeExpr,
                                                   DENSITY = RHS,
                                                   PDIST_LOWER = PDIST_LOWER,
                                                   PDIST_UPPER = PDIST_UPPER,
                                                   DDIST_LOWER = DDIST_LOWER
                                               )), list(e = substCode)))
        }
    }
    return(code)
}

## creates the accessor method to return value 'expr'
ndf_generateGetParamFunction <- function(expr, type, nDim) {
    type <- 'double'  ## (NOTE) node values and paramters are always implemented as doubles
    ans <- try(eval(substitute(
        function() {
            returnType(TYPE(NDIM))
            return(PARAMEXPR)
        },
        list(TYPE=as.name(type), NDIM=nDim, PARAMEXPR=expr)
    )))
    if(inherits(ans, 'try-error')) stop('In ndf_generateGetParamFunction.', .call = FALSE)
    ans
}

ndf_createSingleMethod <- function(type, nDim) {
    type <- 'double'  ## (NOTE) node values and paramters are always implemented as doubles
    methodDef <- substitute(function() { returnType(TYPE(NDIM)) }, list(TYPE=as.name(type), NDIM=nDim))
    methodDef[[4]] <- NULL
    eval(methodDef)
}

ndf_createVirtualNodeFunctionDefinition <- function(types = list()) {
    methodsList <- lapply(types, function(singleType) ndf_createSingleMethod(type=singleType$type, nDim=singleType$nDim))
    if(length(methodsList) > 0)     names(methodsList) <- paste0('get_', names(methodsList))
    virtualFunctionDef <- substitute(
        nimbleFunctionVirtual(
            contains = 'nodeFun',
            methods = METHODS
        ),
        list(METHODS = methodsList)
    )
    return(virtualFunctionDef)
}

ndf_createVirtualNodeFunctionDefinitionsList <- function(userAdded = FALSE) {
    defsList <- list()
    if(!userAdded) {
        defsList$node_determ <- ndf_createVirtualNodeFunctionDefinition()
        for(distName in getAllDistributionsInfo('namesVector', nimbleOnly = TRUE)) {
            defsList[[paste0('node_stoch_', distName)]] <- ndf_createVirtualNodeFunctionDefinition(getDistributionInfo(distName)$types)
        }
    } else {
        # this deals with user-provided distributions
        if(exists('distributions', nimbleUserNamespace, inherits = FALSE)) {
            for(distName in getAllDistributionsInfo('namesVector', userOnly = TRUE))
                defsList[[paste0('node_stoch_', distName)]] <- ndf_createVirtualNodeFunctionDefinition(getDistributionInfo(distName)$types)
        } else stop("ndf_createVirtualNodeFunctionDefinitionsList: no 'distributions' list in nimbleUserNamespace.")
    }
    return(defsList)
}

virtualNodeFunctionDefinitions <- ndf_createVirtualNodeFunctionDefinitionsList()
createNamedObjectsFromList(virtualNodeFunctionDefinitions)

