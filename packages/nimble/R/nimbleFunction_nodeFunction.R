

nodeFunction <- function(LHS, RHS, name = NA, altParams, logProbNodeExpr, type, setupOutputExprs, evaluate = TRUE, where = globalenv()) {
    if(!(type %in% c('stoch', 'determ')))       stop(paste0('invalid argument to nodeFunction(): type = ', type))
    nodeFunctionTemplate <-
        substitute(
            nimbleFunction(contains      = CONTAINS,
                           setup         = SETUPFUNCTION,
                           methods       = METHODS,
                           name          = name,
                           where = where),
            list(CONTAINS      = ndf_createContains(RHS, type),
                 SETUPFUNCTION = ndf_createSetupFunction(setupOutputExprs),
                 METHODS       = ndf_createMethodList(LHS, RHS, altParams, logProbNodeExpr, type, setupOutputExprs),
                 where         = where)
        )
    if(evaluate)    return(eval(nodeFunctionTemplate))     else       return(nodeFunctionTemplate)
}


## creates the name of the node class inheritance (nimbleFunction(contains = ....)
ndf_createContains <- function(RHS, type) {
    if(nimbleOptions()$compileAltParamFunctions) {
        if(type == 'determ')   tag <- 'determ'
        if(type == 'stoch')    tag <- paste0('stoch_', RHS[[1]])
        containsText <- paste0('node_', tag)
        return(as.name(containsText))
    }
    return(NULL)
}

## creates a function object for use as setup argument to nimbleFunction()
ndf_createSetupFunction <- function(setupOutputExprs) {
    setup <- function() {}
    allSetupOutputExprs <- c(quote(model), setupOutputExprs)
    formals(setup) <- nf_createAList(allSetupOutputExprs)
    return(setup)
}


## creates a list of the methods calculate, simulate, and getLogProb, corresponding to LHS, RHS, and type arguments
ndf_createMethodList <- function(LHS, RHS, altParams, logProbNodeExpr, type, setupOutputExprs) {
    if(type == 'determ') {
        methodList <- eval(substitute(
            list(
                simulate   = function() { LHS <<- RHS                                                 },
                calculate  = function() { simulate();    returnType(double());   return(invisible(0)) },
                calculateDiff = function() {simulate();  returnType(double());   return(invisible(0)) },
                getLogProb = function() {                returnType(double());   return(0)            }
            ),
            list(LHS=LHS, 
                 RHS=RHS)))
    }
    if(type == 'stoch') {
        methodList <- eval(substitute(
            list(
                simulate   = function() { LHS <<- STOCHSIM                                                         },
                calculate  = function() { STOCHCALC_FULLEXPR;   returnType(double());   return(invisible(LOGPROB)) },
                calculateDiff = function() {STOCHCALC_FULLEXPR_DIFF; LocalAns <- LocalNewLogProb - LOGPROB;  LOGPROB <<- LocalNewLogProb;
                                            returnType(double());   return(invisible(LocalAns))},
                getLogProb = function() {                       returnType(double());   return(LOGPROB)            }
            ),
            list(LHS       = LHS,
                 LOGPROB   = logProbNodeExpr,
                 STOCHSIM  = ndf_createStochSimulate(RHS),
                 STOCHCALC_FULLEXPR = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS),
                 STOCHCALC_FULLEXPR_DIFF = ndf_createStochCalculate(logProbNodeExpr, LHS, RHS, diff = TRUE))))
        if(nimbleOptions()$compileAltParamFunctions) {
            distName <- as.character(RHS[[1]])
            ## add accessor function for node value; used in multivariate conjugate sampler functions
            type <- getType(distName)
            nDim <- getDimension(distName)
            methodList[['get_value']] <- ndf_generateGetParamFunction(LHS, type, nDim)
            ## add accessor functions for stochastic node distribution parameters
            for(param in names(RHS[-1])) {
                if(!param %in% c("lower", "upper")) {
                    type <- getType(distName, param)
                    nDim <- getDimension(distName, param)
                    methodList[[paste0('get_',param)]] <- ndf_generateGetParamFunction(RHS[[param]], type, nDim)
                }
            }
            for(i in seq_along(altParams)) {
                altParamName <- names(altParams)[i]
                type <- getType(distName, altParamName)
                nDim <- getDimension(distName, altParamName)
                methodList[[paste0('get_',altParamName)]] <- ndf_generateGetParamFunction(altParams[[altParamName]], type, nDim)
            }
            ## new for getParam, eventually to replace get_XXX where XXX is each param name
            ## TO-DO: unfold types and nDims more thoroughly (but all types are implemented as doubles anyway)
            ## understand use of altParams vs. all entries in getDistributionInfo(distName)$types
            ## need a value Entry
            allParams <- c(list(value = LHS), as.list(RHS[-1]), altParams)
            typesNDims <- getDimension(distName, includeParams = TRUE)
            typesTypes <- getType(distName, includeParams = TRUE)
            paramIDs <- getParamID(distName, includeParams = TRUE)

            ## rely on only double for now
            for(nDimSupported in c(0, 1, 2)) {
                boolThisCase <- typesNDims == nDimSupported & typesTypes == 'double'
                paramNamesToUse <- getParamNames(distName)[boolThisCase]
                caseName <- paste0("getParam_",nDimSupported,"D_double")
                if(length(paramNamesToUse) > 0) 
                    methodList[[caseName]] <- ndf_generateGetParamSwitchFunction(allParams[paramNamesToUse], paramIDs[paramNamesToUse], type = 'double', nDim = nDimSupported) 
            }
        }
    }
    ## add model$ in front of all names, except the setupOutputs
    methodList <- ndf_addModelDollarSignsToMethods(methodList, setupOutputExprs, exceptionNames = c("LocalAns", "LocalNewLogProb","PARAMID_","PARAMANSWER_"))
    return(methodList)
}

replaceDistributionAliasesNameOnly <- function(dist) {
    if (as.character(dist) %in% names(distributionAliases)) {
        dist <- as.name(distributionAliases[dist])
    }
    return(dist)
}


## helper function that adds an argument to a call
## used to add needed arguments for C versions of {d,p,q}${dist} functions
addArg <- function(code, value, name) {
    newArgIndex <- length(code) + 1
    code[newArgIndex] <- value
    names(code)[newArgIndex] <-  name
    return(code)
}



## changes 'dnorm(mean=1, sd=2)' into 'rnorm(1, mean=1, sd=2)'
ndf_createStochSimulate <- function(RHS) {
    BUGSdistName <- as.character(RHS[[1]])
    RHS[[1]] <- as.name(getDistributionInfo(BUGSdistName)$simulateName)   # does the appropriate substituion of the distribution name
    if(length(RHS) > 1) {    for(i in (length(RHS)+1):3)   { RHS[i] <- RHS[i-1];     names(RHS)[i] <- names(RHS)[i-1] } }    # scoots all named arguments right 1 position
    RHS[[2]] <- 1;     names(RHS)[2] <- ''    # adds the first (unnamed) argument '1'
    if("lower" %in% names(RHS) || "upper" %in% names(RHS))
        RHS <- ndf_createStochSimulateTrunc(RHS, discrete = getAllDistributionsInfo('discrete')[BUGSdistName])
    return(RHS)
}


## changes 'rnorm(mean=1, sd=2, lower=0, upper=3)' into correct truncated simulation
##   using inverse CDF
ndf_createStochSimulateTrunc <- function(RHS, discrete = FALSE) {
    lowerPosn <- which("lower" == names(RHS))
    upperPosn <- which("upper" == names(RHS))
    lower <- RHS[[lowerPosn]]
    upper <- RHS[[upperPosn]]
    RHS <- RHS[-c(lowerPosn, upperPosn)]
    dist <- substring(as.character(RHS[[1]]), 2, 1000)
    # userDist <- sum(BUGSdist %in% getAllDistributionsInfo('namesVector', userOnly = TRUE))
    # back to using periods in name because we now mangle the nf arg names
    lowerTailName <- 'lower.tail' # ifelse(userDist, 'lower_tail', 'lower.tail')
    logpName <- 'log.p'  # ifelse(userDist, 'log_p', 'log.p')
    logName <- 'log' # ifelse(userDist, 'log_value', 'log')
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
ndf_createStochCalculate <- function(logProbNodeExpr, LHS, RHS, diff = FALSE) {
    BUGSdistName <- as.character(RHS[[1]])
    RHS[[1]] <- as.name(getDistributionInfo(BUGSdistName)$densityName)   # does the appropriate substituion of the distribution name
    if(length(RHS) > 1) {    for(i in (length(RHS)+1):3)   { RHS[i] <- RHS[i-1];     names(RHS)[i] <- names(RHS)[i-1] } }    # scoots all named arguments right 1 position
    RHS[[2]] <- LHS;     names(RHS)[2] <- ''    # adds the first (unnamed) argument LHS
    if("lower" %in% names(RHS) || "upper" %in% names(RHS)) {
        return(ndf_createStochCalculateTrunc(logProbNodeExpr, LHS, RHS, diff = diff, discrete = getAllDistributionsInfo('discrete')[BUGSdistName]))
    } else {
          userDist <- BUGSdistName %in% getAllDistributionsInfo('namesVector', userOnly = TRUE)
          RHS <- addArg(RHS, 1, 'log')  # adds the last argument log=TRUE (log_value for user-defined) # This was changed to 1 from TRUE for easier C++ generation
          if(diff) {
              code <- substitute(LocalNewLogProb <- STOCHCALC,
                                 list(##LOGPROB = logProbNodeExpr,
                                      STOCHCALC = RHS))
          } 
          else{
            code <- substitute(LOGPROB <<- STOCHCALC,
                               list(LOGPROB = logProbNodeExpr,
                                    STOCHCALC = RHS))
          }
          return(code)
    }
}

## changes 'dnorm(mean=1, sd=2, lower=0, upper=3)' into correct truncated calculation
ndf_createStochCalculateTrunc <- function(logProbNodeExpr, LHS, RHS, diff = FALSE, discrete = FALSE) {
    lowerPosn <- which("lower" == names(RHS))
    upperPosn <- which("upper" == names(RHS))
    lower <- RHS[[lowerPosn]]
    upper <- RHS[[upperPosn]]
    RHS <- RHS[-c(lowerPosn, upperPosn)]
    dist <- substring(as.character(RHS[[1]]), 2, 1000)
    # userDist <- sum(as.character(RHS[[1]]) %in% getAllDistributionsInfo('namesVector', userOnly = TRUE))
    # back to using periods in name because we now mangle the nf arg names
    lowerTailName <- 'lower.tail' # ifelse(userDist, 'lower_tail', 'lower.tail')
    logpName <- 'log.p' # ifelse(userDist, 'log_p', 'log.p')
    logName <- 'log' # ifelse(userDist, 'log_value', 'log')

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

    if(diff) {
        if(discrete && lower != -Inf) {
            substCode <- quote(if(LOWER <= VALUE & VALUE <= UPPER)
                                   LocalNewLogProb <- DENSITY - log(PDIST_UPPER - PDIST_LOWER + DDIST_LOWER)
                               else LocalNewLogProb <- -Inf)
        } else {
              substCode <- quote(if(LOWER <= VALUE & VALUE <= UPPER)
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
              substCode <- quote(if(LOWER <= VALUE & VALUE <= UPPER)
                                     LOGPROB <<- DENSITY - log(PDIST_UPPER - PDIST_LOWER + DDIST_LOWER)
                                 else LOGPROB <<- -Inf)
          } else {
                substCode <- quote(if(LOWER <= VALUE & VALUE <= UPPER)
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
    return(code)
}


ndf_generateGetParamSwitchFunction <- function(typesListAll, paramIDs, type, nDim) {
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
            function(PARAMID_ = integer()) {
                returnType(TYPE(NDIM))
                ANSWERINITCODE
                SWITCHCODE
                return(PARAMANSWER_)
            },
            list(TYPE = as.name(type), NDIM=nDim, ANSWERINITCODE = answerInitCode, SWITCHCODE = switchCode)
        )))
    } else {
        ans <- try(eval(substitute(
            function(PARAMID_ = integer()) {
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
    if(inherits(ans, 'try-error')) browser()
    ans
}

## adds model$ on front of all node names, in the bodys of methods in methodList
ndf_addModelDollarSignsToMethods <- function(methodList, setupOutputExprs, exceptionNames = character()) {
    for(i in seq_along(methodList)) {
        body(methodList[[i]]) <-addModelDollarSign(body(methodList[[i]]), exceptionNames = c(exceptionNames, as.character(setupOutputExprs)))
    }
    return(methodList)
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
        if(exists('distributions', nimbleUserNamespace)) {
            for(distName in getAllDistributionsInfo('namesVector', userOnly = TRUE))
                defsList[[paste0('node_stoch_', distName)]] <- ndf_createVirtualNodeFunctionDefinition(getDistributionInfo(distName)$types)
        } else stop("ndf_createVirtualNodeFunctionDefinitionsList: no 'distributions' list in nimbleUserNamespace.")
    }
    return(defsList)
}



virtualNodeFunctionDefinitions <- ndf_createVirtualNodeFunctionDefinitionsList()
createNamedObjectsFromList(virtualNodeFunctionDefinitions)

