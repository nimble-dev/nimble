

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
    }
    ## add model$ in front of all names, except the setupOutputs
    methodList <- ndf_addModelDollarSignsToMethods(methodList, setupOutputExprs, exceptionNames = c("LocalAns", "LocalNewLogProb"))
    return(methodList)
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
    RHS[[1]] <- as.name(getDistribution(as.character(RHS[[1]]))$simulateName)   # does the appropriate substituion of the distribution name
    if(length(RHS) > 1) {    for(i in (length(RHS)+1):3)   { RHS[i] <- RHS[i-1];     names(RHS)[i] <- names(RHS)[i-1] } }    # scoots all named arguments right 1 position
    RHS[[2]] <- 1;     names(RHS)[2] <- ''    # adds the first (unnamed) argument '1'
    if("lower" %in% names(RHS) || "upper" %in% names(RHS))
        RHS <- ndf_createStochSimulateTrunc(RHS)
    return(RHS)
}


## changes 'rnorm(mean=1, sd=2, lower=0, upper=3)' into correct truncated simulation
##   using inverse CDF
ndf_createStochSimulateTrunc <- function(RHS) {
    lowerPosn <- which("lower" == names(RHS))
    upperPosn <- which("upper" == names(RHS))
    lower <- RHS[[lowerPosn]]
    upper <- RHS[[upperPosn]]
    RHS <- RHS[-c(lowerPosn, upperPosn)]
    dist <- substring(as.character(RHS[[1]]), 2, 1000)
    # userDist <- sum(BUGSdist %in% getDistributionsInfo('namesVector', userOnly = TRUE))
    # back to using periods in name because we now mangle the nf arg names
    lowerTailName <- 'lower.tail' # ifelse(userDist, 'lower_tail', 'lower.tail')
    logpName <- 'log.p'  # ifelse(userDist, 'log_p', 'log.p')
    # setup for runif(1, pdist(lower,...), pdist(upper,...))
    # pdist() expression template for inputs to runif()
    pdistTemplate <- RHS
    pdistTemplate[[1]] <- as.name(paste0("p", dist))
    pdistTemplate <- addArg(pdistTemplate, 1, lowerTailName)
    pdistTemplate <- addArg(pdistTemplate, 0, logpName)
    # create bounds for runif() using pdist expressions
    MIN_EXPR <- 0
    MAX_EXPR <- 1
    if(lower != -Inf) {
        pdistTemplate[[2]] <- lower
        MIN_EXPR <- pdistTemplate
    } 
    if(upper != Inf) {
        pdistTemplate[[2]] <- upper
        MAX_EXPR <- pdistTemplate
    }
    
    # now create full runif() expression
    RUNIF_EXPR <- substitute(runif(1, MIN, MAX), list(
        MIN = MIN_EXPR,
        MAX = MAX_EXPR))

    # create full qdist(runif(...),...) expression
    RHS[[1]] <- as.name(paste0("q", dist))
    RHS[[2]] <- RUNIF_EXPR
    RHS <- addArg(RHS, 1, lowerTailName)
    RHS <- addArg(RHS, 0, logpName)

    return(RHS)
}

## changes 'dnorm(mean=1, sd=2)' into 'dnorm(LHS, mean=1, sd=2, log=TRUE)'
ndf_createStochCalculate <- function(logProbNodeExpr, LHS, RHS, diff = FALSE) {
    RHS[[1]] <- as.name(getDistribution(as.character(RHS[[1]]))$densityName)   # does the appropriate substituion of the distribution name
    if(length(RHS) > 1) {    for(i in (length(RHS)+1):3)   { RHS[i] <- RHS[i-1];     names(RHS)[i] <- names(RHS)[i-1] } }    # scoots all named arguments right 1 position
    RHS[[2]] <- LHS;     names(RHS)[2] <- ''    # adds the first (unnamed) argument LHS
    if("lower" %in% names(RHS) || "upper" %in% names(RHS)) {
        return(ndf_createStochCalculateTrunc(logProbNodeExpr, LHS, RHS))
    } else {
          userDist <- as.character(RHS[[1]]) %in% getDistributionsInfo('namesVector', userOnly = TRUE)
          RHS <- addArg(RHS, 1, 'log')  # adds the last argument log=TRUE (log_value for user-defined) # This was changed to 1 from TRUE for easier C++ generation
          if(diff) {
              code <- substitute(LocalNewLogProb <- STOCHCALC,
                                 list(LOGPROB = logProbNodeExpr,
                                      STOCHCALC = RHS))
          } else {
              code <- substitute(LOGPROB <<- STOCHCALC,
                                 list(LOGPROB = logProbNodeExpr,
                                      STOCHCALC = RHS))
          }
          return(code)
    }
}

## changes 'dnorm(mean=1, sd=2, lower=0, upper=3)' into correct truncated calculation
ndf_createStochCalculateTrunc <- function(logProbNodeExpr, LHS, RHS) {
    lowerPosn <- which("lower" == names(RHS))
    upperPosn <- which("upper" == names(RHS))
    lower <- RHS[[lowerPosn]]
    upper <- RHS[[upperPosn]]
    RHS <- RHS[-c(lowerPosn, upperPosn)]
    dist <- substring(as.character(RHS[[1]]), 2, 1000)
    # userDist <- sum(as.character(RHS[[1]]) %in% getDistributionsInfo('namesVector', userOnly = TRUE))
    # back to using periods in name because we now mangle the nf arg names
    lowerTailName <- 'lower.tail' # ifelse(userDist, 'lower_tail', 'lower.tail')
    logpName <- 'log.p' # ifelse(userDist, 'log_p', 'log.p')
    logName <- 'log' # ifelse(userDist, 'log_value', 'log')

    pdistTemplate <- RHS
    pdistTemplate[[1]] <- as.name(paste0("p", dist))
    pdistTemplate <- addArg(pdistTemplate, 1, lowerTailName)
    pdistTemplate <- addArg(pdistTemplate, 0, logpName)

    PDIST_LOWER <- 0
    PDIST_UPPER <- 1
    if(lower != -Inf) {
        pdistTemplate[[2]] <- lower
        PDIST_LOWER <- pdistTemplate
    } 
    if(upper != Inf) {
        pdistTemplate[[2]] <- upper
        PDIST_UPPER <- pdistTemplate
    }

    RHS <- addArg(RHS, 1, logName)  # add log=1 now that pdist() created without 'log'

    # unlike JAGS we have (L < X <= U), i.e., (L,U], as otherwise we would need
    # machinations to deal with the X=L case (pdist functions provide only P(X<=L),P(X>L))
    # this only matters for discrete distributions
    # one option if necessary would be to check for discrete and then use ddist too
    code <- substitute(if(LOWER <= VALUE & VALUE <= UPPER)
                           LOGPROB <<- DENSITY - log(PDIST_UPPER - PDIST_LOWER)
                       else LOGPROB <<- -Inf,
                       list(
                           LOWER = lower,
                           UPPER = upper,
                           VALUE = LHS,
                           LOGPROB = logProbNodeExpr,
                           DENSITY = RHS,
                           PDIST_LOWER = PDIST_LOWER,
                           PDIST_UPPER = PDIST_UPPER
                       ))
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
    virtualFuncionDef <- substitute(
        nimbleFunctionVirtual(
            contains = 'nodeFun',
            methods = METHODS
        ),
        list(METHODS = methodsList)
    )
    return(virtualFuncionDef)
}

ndf_createVirtualNodeFunctionDefinitionsList <- function(userAdded = FALSE) {
    defsList <- list()
    if(!userAdded) {
        defsList$node_determ <- ndf_createVirtualNodeFunctionDefinition()
        for(distName in getDistributionsInfo('namesVector', nimbleOnly = TRUE)) {
            defsList[[paste0('node_stoch_', distName)]] <- ndf_createVirtualNodeFunctionDefinition(getDistribution(distName)$types)
        }
    } else {
        # this deals with user-provided distributions
        if(exists('distributions', nimbleUserNamespace)) {
            for(distName in getDistributionsInfo('namesVector', userOnly = TRUE))
                defsList[[paste0('node_stoch_', distName)]] <- ndf_createVirtualNodeFunctionDefinition(getDistribution(distName)$types)
        } else stop("ndf_createVirtualNodeFunctionDefinitionsList: no 'distributions' list in nimbleUserNamespace.")
    }
    return(defsList)
}



virtualNodeFunctionDefinitions <- ndf_createVirtualNodeFunctionDefinitionsList()
createNamedObjectsFromList(virtualNodeFunctionDefinitions)

