

distributionsClass <- setRefClass(
    Class = 'distributionsClass',
    
    fields = list(
        distObjects   = 'ANY',		#'list',      ## a list of distClass objects, names of each element are the BUGS distribution name
        namesVector   = 'ANY',		#'character',      ## a character vector of the (BUGS) names of all distributions
        namesExprList = 'ANY',		#'list',   ## a list of the expressions of the (BUGS) names of all distributions
        matchCallEnv  = 'ANY',		#'environment',   ## an environment containing distribution functions which run match.call()
        translations  = 'ANY'		#'list'   ## a list of the (R) d-dist and r-dist function names. element names are BUGS distributions
    ),
    
    methods = list(
        initialize = function(dil) {
            distObjects <<- list()
            namesExprList <<- list()
            translations <<- list()
            for(i in seq_along(dil))     distObjects[[i]] <<- distClass(dil[[i]], names(dil)[i])
            names(distObjects) <<- names(dil)
            namesVector <<- names(dil)
            namesExprList <<- lapply(namesVector, as.name)
            matchCallEnv <<- new.env()
            for(distName in namesVector)     assign(distName, distObjects[[distName]]$makeMatchCallFunction(), matchCallEnv)
            translations <<- lapply(distObjects, function(d) c(d$densityName, d$simulateName))
        },

        add = function(dil) {
              distObjectsNew <- list()
              nms <- names(dil)
              dupl <- which(nms %in% getAllDistributionsInfo('namesVector', userOnly = TRUE))
              if(length(dupl)) {
                  for(i in seq_along(dupl)) {
                      remove(nms[dupl[i]])
                  }
                  ## distObjects[dupl] <<- NULL
                  ## namesVector <<- namesVector[-dupl]
                  ## namesExprList[dupl] <<- NULL
                  ## translations[dupl] <<- NULL
                  nmsDuplicated <- paste0(nms[dupl], collapse = ', ')
                  cat(paste0("Overwriting the following user-supplied distributions: ", nmsDuplicated, "\n"))
              }
              for(i in seq_along(dil))     distObjectsNew[[i]] <- distClass(dil[[i]], nms[i])
              names(distObjectsNew) <- nms
              translations <<- c(translations, lapply(distObjectsNew, function(d) c(d$densityName, d$simulateName)))

              distObjects <<- c(distObjects, distObjectsNew)
              namesVector <<- c(namesVector, nms)
              namesExprList <<- c(namesExprList, lapply(nms, as.name))
              for(distName in nms) assign(distName, distObjects[[distName]]$makeMatchCallFunction(), matchCallEnv)
          },

        remove = function(dn) {
            namesVector <<- namesVector[!namesVector %in% dn]
            namesExprList[namesExprList == as.name(dn)] <<- NULL
            eval(substitute(rm(x, envir = matchCallEnv), list(x = dn)))
            translations[dn] <<- NULL
            distObjects[dn] <<- NULL
        }
        )
    )
              

setMethod('[[',   'distributionsClass',
          function(x, i) {
              return(x$distObjects[[i]])
          }
)

setMethod('[',   'distributionsClass',
          function(x, i) {
              return(x$distObjects[i])
          }
)


distClass <- setRefClass(
    Class = 'distClass',
    
    fields = list(
        BUGSdistName = 'ANY',	#'character',   ## the (BUGS) name of the distribution
        BUGSdistExpr = 'ANY',   # the BUGS distribution expression, as provided in the original inputs list, with all possible parameter names
        RdistExprList = 'ANY',	#'list',  ## a list of the R distribution expressions, along with their parameters and re-parametrizations
        numAlts = 'ANY',	#'numeric',   ## the number of alternate reparametrizations provided
        alts = 'ANY',		#'list',
        exprs = 'ANY',		#'list',
        reqdArgs = 'ANY',	#'character',   ## chracter vector of the required arguments in our R implementation of each distribution; we always reparametrize to this
        densityName = 'ANY',	#'character',   ## the (R) name of the d-dist function, e.g. 'dnorm'
        simulateName = 'ANY',	#'character',   ## the (R) name of the r-dist function, e.g. 'rnorm'
        altParams = 'ANY',	#'list',    ## the (named) list of alternate parameters we'll have available, list elements are the expressions for each parameter 
        discrete = 'ANY',	#'logical',   ## logical, if the distribution is discrete
        pqAvail = 'ANY',        #'logical', ## if the p (CDF) and q (inverse CDF/quantile) functions are available
        mixedSizes = 'ANY',     ##   if TRUE, then parameters of this distribution could have varied sizes, and is exempted from this check in model$checkBasics()
        range = 'ANY',          #'numeric',  ## lower and upper limits of distribution domain
        types = 'ANY',		#'list',     ## named list (names are 'node', ALL reqdArgs, and ALL altParams), each element is a named list: list(type = 'double', nDim = 0) <- default values
        paramIDs = 'ANY'        #'integer'   ## named vector of unique integer ID for each parameter
### typesForVirtualNodeFunction = 'ANY'		#'list'  ## version of 'types' for making the virtualNodeFunction definiton.  same as above, except without 'value'
    ),
    
    methods = list(
        initialize = function(distInputList, BUGSdistName) {
            RdistExprList <<- list()
            altParams <<- list()
            types <<- list()
            BUGSdistName <<- BUGSdistName
            BUGSdistExpr <<- parse(text=distInputList$BUGSdist)[[1]]
            if(BUGSdistExpr[[1]] != BUGSdistName)   stop(paste0('inconsistent BUGS distribution names for distribution: ', BUGSdistName))
            RdistTextVector <- if(is.null(distInputList$Rdist)) character() else distInputList$Rdist
            RdistExprList <<- lapply(RdistTextVector, function(t) parse(text=t)[[1]])
            init_altsExprsReqdArgs()
            numAlts <<- length(alts)
            simulateName <<- sub('^d', 'r', densityName)
            init_altParams(distInputList)
            discrete <<- if(is.null(distInputList$discrete))    FALSE    else    distInputList$discrete
            pqAvail <<- if(is.null(distInputList$pqAvail))    FALSE    else    distInputList$pqAvail
            mixedSizes <<- if(is.null(distInputList$mixedSizes))    FALSE    else    distInputList$mixedSizes
            init_range(distInputList)
            init_types(distInputList)
            init_paramIDs()
        },
        
        init_altsExprsReqdArgs = function() {
            alts <<- list()
            exprs <<- list()
            if(length(RdistExprList) == 0) {
                params <- as.list(BUGSdistExpr[-1])   # removes the distribution name
                paramsText <- lapply(params, deparse)
                reqdArgs <<- sapply(paramsText, function(pt) init_getReqdArgs(pt))
                densityName <<- as.character(BUGSdistExpr[[1]])
            } else {
                params <- lapply(RdistExprList, `[`, -1)        # removes the distribution names
                paramsText <- lapply(params, function(x) lapply(x, deparse))
                reqdArgsList <- lapply(paramsText, function(pt) init_getReqdArgs(pt))
                densityNamesList <- lapply(RdistExprList, function(expr) as.character(expr[[1]]))
                if(length(unique(lapply(reqdArgsList, sort))) > 1)
                    stop('R/NIMBLE parameter names and order not consistent across alternative parameterizations')
                if(length(unique(densityNamesList)) > 1)
                    stop('R/NIMBLE density names not consistent across alternative parameterizations')
                reqdArgs <<- reqdArgsList[[1]]
                densityName <<- densityNamesList[[1]]
                for(i in seq_along(params)) {
                    boolNoDefault <- if (is.null(names(paramsText[[i]]))) rep(TRUE, length(paramsText[[i]])) else names(paramsText[[i]]) == ''
                    if(sum(!boolNoDefault)) {
                        exprs[[i]] <<- lapply(params[[i]][!boolNoDefault], function(x) {names(x) <- NULL; x})
                        BUGSargs <- unique(unlist(c(lapply(exprs[[i]], all.vars), paramsText[[i]][boolNoDefault])))
                        names(BUGSargs) <- NULL
                        if(!identical(sort(BUGSargs), sort(reqdArgs))) alts[[i]] <<- BUGSargs
                    } else {
                        if(!identical(sort(unlist(paramsText[[i]])), sort(reqdArgs)))   stop(paste0('reparametization number ', i, ' for ', BUGSdistName, ' with no default argument values must exactly match arguments of canonical parameterization'))
                    }
                }
            }
        },
        
        init_getReqdArgs = function(x) {
            args <- if(is.null(names(x))) rep('', length(x)) else names(x)
            args[args == ''] <- unlist(x[args == ''])
            return(args)
        },
        
        init_range = function(distInputList) {
            if(!is.null(distInputList$range)) {
                if(length(distInputList$range) != 2)
                    stop("'Range' element of ", BUGSdistExpr[[1]], " must be a vector of length two.")
                if(is.numeric(distInputList$range)) {
                    range <<-list(lower = distInputList$range[1], upper = distInputList$range[2])
                } else {  
                    parsedRangeArg <- lapply(distInputList$range, function(x) parse(text=x)[[1]])
                    range <<- lapply(parsedRangeArg, function(x) x[[3]])
                    names(range) <<- unlist(lapply(parsedRangeArg, function(x) x[[2]]))
                    if(!identical(names(range), c('lower', 'upper')))
                        stop("'Range' element of ", BUGSdistExpr[[1]], " expected to contain 'lower' and 'upper'.")
                }
            } else range <<- list(lower = -Inf, upper = Inf)
        },

        init_altParams = function(distInputList) {
            altParams <<- list()
            if(!is.null(distInputList$altParams)) {
                parsedAltParamArg <- lapply(distInputList$altParams, function(x) parse(text=x)[[1]])
                altParams <<- lapply(parsedAltParamArg, function(x) x[[3]])
                names(altParams) <<- unlist(lapply(parsedAltParamArg, function(x) x[[2]]))
            }
        },
        
        init_types = function(distInputList) {
            typeArgCharVector <- if(!is.null(distInputList$types)) distInputList$types else character(0)
            typeArgList <- init_types_makeArgList(typeArgCharVector)
            if('value' %in% c(reqdArgs, names(altParams)))    stop('going to have a name conflict with \'value\' in distribution declaration')
            allTypeNames <- c('value', reqdArgs, names(altParams))
            for(typeName in allTypeNames) {
                typeList <- if(typeName %in% names(typeArgList))     typeArgList[[typeName]]     else     list(type='double', nDim=0)   # default type
                if(!(typeList$type %in% c('double', 'integer', 'logical')))     stop(paste0('unknown type specified in distribution: ', typeList$type))
                if(typeList$nDim > 0 && typeList$type != 'double') 
                    stop("Non-scalar integer or logical found in distribution function.\nPlease use type 'double' for all non-scalars in distribution functions.")
                if(!(typeList$nDim %in% 0:1000))     stop(paste0('unknown nDim specified in distribution: ', typeList$nDim))  ## yes, specificying maximum dimension of 1000
                types[[typeName]] <<- typeList
            }            
        },

        init_paramIDs = function() {
            paramIDs <<- seq_along(types)
            names(paramIDs) <<- names(types)
        },
        
        init_types_makeArgList = function(typeArgCharVector) {
            parsedArgList <- try(lapply(typeArgCharVector, function(x) parse(text=x, keep.source = FALSE)[[1]]))
            if(is(parsedArgList, 'try-error'))
                stop("init_types_makeArgList: problem with arguments ", paste(typeArgCharVector, collapse = ","), ". Perhaps you didn't define types for your user-defined distribution nimbleFunctions?")
            allNames <- unlist(lapply(parsedArgList, function(pa) as.character(pa[[2]])))
            if('x' %in% allNames) {
                warning("init_types_makeArgList: Found 'x' in 'types', changing to 'value'.")
                allNames[which('x' == allNames)] <- 'value'
            }
            declExprs <- lapply(parsedArgList, function(pa) pa[[3]])
            allTypes <- unlist(lapply(parsedArgList, function(pa) as.character(pa[[3]][[1]])))
            allDims <- unlist(lapply(parsedArgList, function(pa) if(length(pa[[3]]) == 1) 0 else as.numeric(pa[[3]][[2]])))
            argList <- list()
            for(i in seq_along(allNames)) {
                argList[[allNames[i]]] <- list(type = allTypes[i], nDim = allDims[i])
            }
            return(argList)
        },
        
        makeMatchCallFunction = function() {
            vars <- BUGSdistExpr[-1]
            functionText <- paste0('function(', paste0(vars, collapse=', '), ') { match.call() }')
            functionDef <- parse(text = functionText)[[1]]
            eval(functionDef)
        }
    )
)

#####################################################################################################
#####################################################################################################
#####  process user-supplied distributions ##########################################################
#####################################################################################################
#####################################################################################################

checkDistributionInput <- function(distributionInput) {
    allowedFields <- unique(unlist(sapply(distributionsInputList, names)))
    if(sum(!names(distributionInput) %in% allowedFields)) 
        stop(paste0(names(distributionInput), " has unknown field."))
    if(!sum(is.character(distributionInput$BUGSdist))) stop(paste0(distributionInput$BUGSdist, ": field 'BUGSdist' is not of type character."))
    if(exists("Rdist", distributionInput, inherits = FALSE) && !sum(is.character(distributionInput$Rdist))) stop(paste0(distributionInput$BUGSdist, ": field 'Rdist' is not type of character."))
    if(exists("discrete", distributionInput, inherits = FALSE) && !sum(is.logical(distributionInput$discrete))) stop(paste0(distributionInput$BUGSdist, ": field 'discrete' is not type logical."))
    if(exists("pqAvail", distributionInput, inherits = FALSE) && !sum(is.logical(distributionInput$pqAvail))) stop(paste0(distributionInput$BUGSdist, ": field 'pqAvail' is not of type logical."))
    if(exists("range", distributionInput, inherits = FALSE) && (!is.numeric(distributionInput$range) || length(distributionInput$range) != 2)) stop(paste0(distributionInput$BUGSdist, ": field 'range' is not a vector of two numeric values."))
    if(exists("types", distributionInput, inherits = FALSE) && !sum(is.character(distributionInput$types))) stop(paste0(distributionInput$BUGSdist, ": field 'types' is not of type character."))
    if(exists("altParams", distributionInput, inherits = FALSE) && !sum(is.character(distributionInput$altParams))) stop(paste0(distributionInput$BUGSdist, ": field 'altParams' is not of type character."))
    if(length(distributionInput$BUGSdist) > 1 || (exists('discrete', distributionInput, inherits = FALSE) && length(distributionInput$discrete) > 1) || (exists('pqAvail', distributionInput, inherits = FALSE) && length(distributionInput$pqAvail) > 1))
        stop(paste0(names(distributionInput), " field 'BUGSdist', 'discrete', 'altParams', or 'pqAvail' is not of length one."))
    invisible(NULL)
}


# check for log last in d, n as first in r, lower.tail, log.p in p,q
checkDistributionFunctions <- function(distributionInput, userEnv) {
    if(is.list(distributionInput)) {
        if(exists('Rdist', distributionInput, inherits = FALSE))
            inputString <- distributionInput$Rdist else inputString <- distributionInput$BUGSdist
        densityName <- as.character(parse(text = inputString)[[1]][[1]])
    } else densityName <- distributionInput
    simulateName <- sub('^d', 'r', densityName)
    nofun <- FALSE
    if(exists(densityName, where = userEnv)) {
        rcf <- get(densityName, pos = userEnv)
        if(!is.rcf(rcf)) {
            nofun <- TRUE
        } else if(environment(rcf)$nfMethodRCobject$returnType != quote(double()) &&
                  environment(rcf)$nfMethodRCobject$returnType != quote(double(0)))
              stop(paste0("checkDistributionFunctions: density function for ", densityName,
                          " has invalid or missing returnType, which must be 'double(0)' (or equivalently 'double()')."))
    } else nofun <- TRUE
    if(nofun) {
        if(densityName %in% c('+','-','*','/','%%','%*%','[','[[','$','^','|','||','&','&&',':','<','<=','>','>=','!=','==')) {
            stop(paste0("checkDistributionFunctions: expression '", densityName,
                        "' found where a density function is expected. Did you mistakenly use `~` instead of `<-`?"))
        } else stop(paste0("checkDistributionFunctions: density function for ", densityName,
                    " is not available.  It must be a nimbleFunction (with no setup code)."))
    }
    if(!exists(simulateName, where = userEnv) || !is.rcf(get(simulateName, pos = userEnv))) {
        cat(paste0("Warning: random generation function for ", densityName,
                    " is not available. NIMBLE is generating a placeholder function, ", simulateName, ", that will invoke an error if an algorithm needs to simulate from this distribution. Some algorithms (such as random-walk Metropolis MCMC sampling) will work without the ability to simulate from the distribution.  If simulation is needed, provide a nimbleFunction (with no setup code) to do it.\n"))
        rargInfo <- environment(get(densityName, pos = userEnv))$nfMethodRCobject$argInfo
        returnType <- deparse(unlist(rargInfo[[1]]))
        returnDim <- 0
        if(length(rargInfo[[1]]) > 1)
            returnDim <- rargInfo[[1]][[2]]
        rargInfo <- rargInfo[-length(rargInfo)]  # remove 'log' argument
        rargInfo[[1]] <- quote(integer(0))
        names(rargInfo)[1] <- 'n'
        args <- paste(names(rargInfo), as.character(rargInfo), sep = "=", collapse = ', ')
        if(returnDim == 0)
            returnCreation <- "x <- 0" else if(returnDim == 1) returnCreation <- "x <- nimNumeric()" else
                                                          returnCreation <- "x <- nimMatrix()"
        # build nf from text as unclear how to pairlist info in rargInfo with substitute
        nfCode <- paste0("nimbleFunction(run = function(", args, ") { stop('user-defined distribution ", densityName, " provided without random generation function.')\nreturnType(", returnType, ")\n", returnCreation, "\nreturn(x)})")
	## Want to assign to same environmet as the 'd' function.
        ## If user does use `assign` to put in GlobalEnv, we don't
        ## do that here automatically as CRAN policy says packages should not modify the GlobalEnv.
        ## Should be ok in terms of running the model inside a function, so long as
        ## the simulate function is not called, which it shouldn't be as it is a dummy.
        assign(simulateName, eval(parse(text = nfCode)), userEnv)

    }

    dargs <- args <- formals(get(densityName, pos = userEnv))
    nArgs <- length(args)
    if(nArgs < 2) stop(paste0("checkDistributionFunctions: expecting at least two arguments ('x', 'log') as arguments for the density function for ", densityName, "."))
    if(names(args)[1] != "x") stop(paste0("checkDistributionFunctions: expecting 'x' as the first argument for the density function for ", densityName, "."))
    if(names(args)[nArgs] != "log") stop(paste0("checkDistributionFunctions: expecting 'log' as the last argument for the density function for ", densityName, "."))
    dargs <- dargs[-c(1,nArgs)]
    
    rargs <- args <- formals(get(simulateName, pos = userEnv))
    nArgs <- length(args)
    if(nArgs < 1) stop(paste0("checkDistributionFunctions: expecting at least one argument ('n') as arguments for the simulation function for ", densityName, "."))
    if(names(args)[1] != "n") stop(paste0("checkDistributionFunctions: expecting 'n' as the first argument for the simulation function for ", densityName, "."))
    rargs <- rargs[-1]
    
    if(!identical(dargs, rargs))
        warning(paste0("checkDistributionFunctions: parameter arguments not the same amongst density and simulation functions for ", densityName, ". Continuing anyway based on arguments to the density function; algorithms using the simulation function are unlikely to function properly."))
    if(!is.null(distributionInput) && is.list(distributionInput) && exists("pqAvail", distributionInput, inherits = FALSE) &&
       distributionInput$pqAvail) {
        cdfName <- sub('^d', 'p', densityName)
        quantileName <- sub('^d', 'q', densityName)
        if(!is.rcf(get(cdfName, pos = userEnv)) || !is.rcf(get(quantileName, pos = userEnv)))
            stop(paste0("checkDistributionFunctions: Either distribution (CDF) or quantile (inverse CDF) functions for ", densityName,
                        " are not available.  If needed, they must be nimbleFunction (with no setup code)."))

        pargs <- args <- formals(get(cdfName, pos = userEnv))
        nArgs <- length(args)
        if(nArgs < 3) stop(paste0("checkDistributionFunctions: expecting at least three arguments ('q', 'lower.tail', and 'log.p') as arguments for the distribution function for ", densityName, "."))
        if(names(args)[1] != "q") stop(paste0("checkDistributionFunctions: expecting 'q' as the first argument for the distribution function for ", densityName, "."))
        if(names(args)[nArgs] != "log.p") stop(paste0("checkDistributionFunctions: expecting 'log.p' as the last argument for the distribution function for ", densityName, "."))
        if(names(args)[nArgs-1] != "lower.tail") stop(paste0("checkDistributionFunctions: expecting 'lower.tail' as the last argument for the distribution function for ", densityName, "."))
        pargs <- pargs[-c(1,nArgs-1,nArgs)]
        
        qargs <- args <- formals(get(quantileName, pos = userEnv))
        nArgs <- length(args)
        if(nArgs < 3) stop(paste0("checkDistributionFunctions: expecting at least three arguments ('p', 'lower.tail', and 'log.p') as arguments for the quantile function for ", densityName, "."))
        if(names(args)[1] != "p") stop(paste0("checkDistributionFunctions: expecting 'p' as the first argument for the quantile function for ", densityName, "."))
        if(names(args)[nArgs] != "log.p") stop(paste0("checkDistributionFunctions: expecting 'log.p' as the last argument for the quantile function for ", densityName, "."))
        if(names(args)[nArgs-1] != "lower.tail") stop(paste0("checkDistributionFunctions: expecting 'lower.tail' as the last argument for the quantile function for ", densityName, "."))
        qargs <- qargs[-c(1,nArgs-1,nArgs)]

        if(!identical(dargs, pargs) || !identical(dargs, qargs))
            stop(paste0("checkDistributionFunctions: parameter arguments not the same amongst density, distribution, and quantile functions for ", densityName, "."))
    }
    invisible(NULL)
}


getMaxDim <- function(typeList) 
    max(sapply(typeList, '[[', 'nDim'))

getValueDim <- function(distObject) 
    distObject$types$value$nDim

prepareDistributionInput <- function(densityName, userEnv) {
    out <- list()
    dist <- get(densityName, pos = userEnv)
    args <- formals(dist)
    args <- args[!names(args) %in% c('x', 'log')]
    if(!length(args))
        argInfo <- NULL
    if(length(args) == 1)
        argInfo <- names(args)
    if(length(args) > 1)
        argInfo <- paste0(paste0(names(args)[-length(args)], sep = ',', collapse = ''), names(args)[length(args)],
                          collapse = '')
    out$BUGSdist <- paste0(densityName, "(", argInfo, ")", collapse = '')
    typeInfo <- get('nfMethodRCobject', environment(eval(as.name(densityName), envir = userEnv)))$argInfo
    out$types <- paste0('value = ', deparse(typeInfo$x))
    typeInfo <- typeInfo[!names(typeInfo) %in% c('x', 'log')]
    if(length(typeInfo))
        out$types <- c(out$types, paste0(names(typeInfo), ' = ', sapply(typeInfo, deparse)))

    # check consistent types
    simulateName <- sub('^d', 'r', densityName)
    typeInfo <- get('nfMethodRCobject', environment(eval(as.name(simulateName), envir = userEnv)))$argInfo
    typeInfo <- typeInfo[names(typeInfo) != "n"]
    rtypes <- character(0)
    if(length(typeInfo))
        rtypes <- c(rtypes, paste0(names(typeInfo), ' = ', sapply(typeInfo, deparse)))
    if(!identical(out$types[-1], rtypes))
        stop(paste0("prepareDistributionInfo: types/dimensions of parameters are not the same in density and simulation functions: '", densityName, "' and '", simulateName, "'."))

    # check for p and q functions
    cdfName <- sub('^d', 'p', densityName)
    quantileName <- sub('^d', 'q', densityName)
    out$pqAvail <- exists(cdfName, where = userEnv) && exists(quantileName, where = userEnv) && is.rcf(get(cdfName, pos = userEnv)) &&
        is.rcf(get(quantileName, pos = userEnv))

    # check consistent types
    if(out$pqAvail) {
        typeInfo <- get('nfMethodRCobject', environment(eval(as.name(cdfName), envir = userEnv)))$argInfo
        typeInfo <- typeInfo[!names(typeInfo) %in% c('q', 'log.p', 'lower.tail')]
        ptypes <- numeric(0)
        if(length(typeInfo))
            ptypes <- c(ptypes, paste0(names(typeInfo), ' = ', sapply(typeInfo, deparse)))
        if(!identical(out$types[-1], ptypes))
            stop(paste0("prepareDistributionInfo: types/dimensions of parameters are not the same in the density and distribution functions: '", densityName, "' and '", cdfName, "'."))
  
        typeInfo <- get('nfMethodRCobject', environment(eval(as.name(quantileName), envir = userEnv)))$argInfo
        typeInfo <- typeInfo[!names(typeInfo) %in% c('p', 'log.p', 'lower.tail')]
        qtypes <- numeric(0)
        if(length(typeInfo))
            qtypes <- c(qtypes, paste0(names(typeInfo), ' = ', sapply(typeInfo, deparse)))
        if(!identical(out$types[-1], qtypes))
            stop(paste0("prepareDistributionInfo: types/dimensions of parameters are not the same in the density and quantile functions: '", densityName, "' and '", quantileName, "'."))
    }
    return(out)
}
    
#' Add user-supplied distributions for use in NIMBLE BUGS models
#'
#' Register distributional information so that NIMBLE can process
#' user-supplied distributions in BUGS model code
#'
#' @param distributionsInput either a list or character vector specifying the user-supplied distributions. If a list, it should be a named list of lists in the form of that shown in \code{nimble:::distributionsInputList} with each list having required field \code{BUGSdist} and optional fields \code{Rdist}, \code{altParams}, \code{discrete}, \code{pqAvail}, \code{types}, and with the name of the list the same as that of the density function. Alternatively, simply a character vector providing the names of the density functions for the user-supplied distributions.
#' @param userEnv environment in which to look for the nimbleFunctions that provide the distribution; this will generally not need to be set by the user as it will default to the environment from which this function was called.
#' @param verbose logical indicating whether to print additional logging information
#' 
#' @author Christopher Paciorek
#' @export
#' @details
#' When \code{distributionsInput} is a list of lists, see below for more information on the structure of the list. When \code{distributionsInput} is a character vector, the distribution is assumed to be of standard form, with parameters assumed to be the arguments provided in the density nimbleFunction, no alternative parameterizations, and the distribution assumed to be continuous with range from minus infinity to infinity. The availability of distribution and quantile functions is inferred from whether appropriately-named functions exist in the global environment.
#'
#' One usually does not need to explicitly call \code{registerDistributions} as it will be called automatically when the user-supplied distribution is used for the first time in BUGS code. However, if one wishes to provide alternative parameterizations, to provide a range, or to indicate a distribution is discrete, then one still must explicitly register the distribution using \code{registerDistributions} with the argument in the list format.
#'
#' Format of the component lists when \code{distributionsInput} is a list of lists:
#' \itemize{
#' \item{\code{BUGSdist}} {
#' a character string in the form of the density name (starting with 'd') followed by the names of the parameters in parentheses. When alternative parameterizations are given in \code{Rdist}, this should be an exhaustive list of the unique parameter names from all possible parameterizations, with the default parameters specified first.
#' }
#' \item{\code{Rdist}} {
#' an optional character vector with one or more alternative specifications of the density; each alternative specification can be an alternative name for the density, a different ordering of the parameters, different parameter name(s), or an alternative parameterization. In the latter case, the character string in parentheses should provide a given reparameterization as comma-separated name = value pairs, one for each default parameter, where name is the name of the default parameter and value is a mathematical expression relating the default parameter to the alternative parameters or other default parameters. The default parameters should correspond to the input arguments of the nimbleFunctions provided as the density and random generation functions. The mathematical expression can use any of the math functions allowed in NIMBLE (see the \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual}) as well as user-supplied nimbleFunctions (which must have no setup code). The names of your nimbleFunctions for the distribution functions must match the function name in the \code{Rdist} entry (or if missing, the function name in the \code{BUGSdist} entry
#' }
#' \item{\code{discrete}} {
#' a optional logical indicating if the distribution is that of a discrete random variable. If not supplied, distribution is assumed to be for a continuous random variable.
#' }
#' \item{\code{pqAvail}} {
#' an optional logical indicating if distribution (CDF) and quantile (inverse CDF) functions are provided as nimbleFunctions. These are required for one to be able to use truncated versions of the distribution. Only applicable for univariate distributions. If not supplied, assumed to be FALSE.
#' }
#' \item{\code{altParams}} {
#' a character vector of comma-separated 'name = value' pairs that provide the mathematical expressions relating non-canonical parameters to canonical parameters (canonical parameters are those passed as arguments to your distribution functions). These inverse functions are used for MCMC conjugacy calculations when a conjugate relationship is expressed in terms of non-default parameters (such as the precision for normal-normal conjugacy). If not supplied, the system will still function but with a possible loss of efficiency in certain algorithms.
#' }
#' \item{\code{types}} {
#' a character vector of comma-separated 'name = input' pairs indicating the type and dimension of the random variable and parameters (including default and alternative parameters). 'input' should take the form 'double(d)' or 'integer(d)', where 'd' is 0 for scalars, 1 for vectors, 2 for matrices. Note that since NIMBLE uses doubles for numerical calculations and the default type  is \code{double(0)}, one should generally use 'double' and one need only specify the type for non-scalars. 'name' should be either 'value' to indicate the random variable itself or the parameter name to indicate a given parameter.  
#' }
#' \item{\code{range}} {
#' a vector of two values giving the range of the distribution for possible use in future algorithms (not used currently). When the lower or upper limit involves a strict inequality (e.g., $x>0$), you should simply treat it as a non-strict inequality ($x>=0$, and set the lower value to 0). Also we do not handle ranges that are functions of parameters, so simply use the smallest/largest possible values given the possible parameter values. If not supplied this is taken to be \code{(-Inf, Inf)}.
#' }
#' }
#' @examples
#' dmyexp <- nimbleFunction(
#'    run = function(x = double(0), rate = double(0), log = integer(0)) {
#'        returnType(double(0))
#'        logProb <- log(rate) - x*rate
#'        if(log) {
#'            return(logProb)
#'        } else {
#'            return(exp(logProb))
#'        }
#'    })
#' rmyexp <- nimbleFunction(
#'    run = function(n = integer(0), rate = double(0)) {
#'        returnType(double(0))
#'        if(n != 1) nimPrint("rmyexp only allows n = 1; using n = 1.")
#'        dev <- runif(1, 0, 1)
#'        return(-log(1-dev) / rate)
#'    }
#'    )
#' registerDistributions(list(
#'     dmyexp = list(
#'               BUGSdist = "dmyexp(rate, scale)",
#'               Rdist = "dmyexp(rate = 1/scale)",
#'               altParams = "scale = 1/rate",
#'               pqAvail = FALSE)))
#' code <- nimbleCode({
#'     y ~ dmyexp(rate = r)
#'     r ~ dunif(0, 100)
#' })
#' m <- nimbleModel(code, inits = list(r = 1), data = list(y = 2))
#' m$calculate('y')
#' m$r <- 2
#' m$calculate('y')
#' m$resetData()
#' m$simulate('y')
#' m$y
#'
#' # alternatively, simply specify a character vector with the
#' # name of one or more 'd' functions
#' deregisterDistributions('dmyexp')
#' registerDistributions('dmyexp')
#'
#' # or simply use in BUGS code without registration
#' deregisterDistributions('dmyexp')
#' m <- nimbleModel(code, inits = list(r = 1), data = list(y = 2))
#'
#' # example of Dirichlet-multinomial registration to illustrate
#' # use of 'types' (note that registration is not actually needed
#' # in this case)
#' ddirchmulti <- nimbleFunction(
#'     run = function(x = double(1), alpha = double(1), size = double(0), 
#'                    log = integer(0, default = 0)) {
#'         returnType(double(0))
#'         logProb <- lgamma(size) - sum(lgamma(x)) + lgamma(sum(alpha)) - 
#'             sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + 
#'                                                                  size)
#'         if(log) return(logProb)
#'         else return(exp(logProb))
#'     })
#'
#' rdirchmulti <- nimbleFunction(
#'     run = function(n = integer(0), alpha = double(1), size = double(0)) {
#'         returnType(double(1))
#'         if(n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
#'         p <- rdirch(1, alpha)
#'         return(rmulti(1, size = size, prob = p))
#'     })
#'
#' registerDistributions(list(
#'     ddirchmulti = list(
#'         BUGSdist = "ddirchmulti(alpha, size)",
#'         types = c('value = double(1)', 'alpha = double(1)')
#'         )
#'     ))
registerDistributions <- function(distributionsInput, userEnv = parent.frame(), verbose = nimbleOptions('verbose')) {
    if(missing(distributionsInput)) {
        cat("No distribution information supplied\n")
    } else {
        if(!(is.character(distributionsInput) || (is.list(distributionsInput) &&
                                                  (length(distributionsInput) == 1 || is.list(distributionsInput[[1]])))))
                                                   stop("'distributionsInput' should be a named list of lists or a character vector")
        if(is.character(distributionsInput)) {
           nms <- distributionsInput
         } else {
            nms <- names(distributionsInput)
          }
        if(verbose) {
            nmsTogether <- paste0(nms, collapse = ', ')
            cat(paste0("Registering the following user-provided distributions: ", nmsTogether, "\n"))
        }
        dupl <- nms[nms %in% getAllDistributionsInfo('namesVector', nimbleOnly = TRUE)]
        if(length(dupl)) {
            distributionsInput[dupl] <- NULL
            duplTogether <- paste0(dupl, collapse = ', ')
            cat(paste0("Ignoring the following user-supplied distributions as they have the same names as default NIMBLE distributions: ", duplTogether, ". Please rename to avoid the conflict.\n"))
          }

        if(is.list(distributionsInput)) 
           sapply(distributionsInput, checkDistributionInput)
        sapply(distributionsInput, checkDistributionFunctions, userEnv = userEnv)
        if(is.character(distributionsInput)) {
          distributionsInput <- lapply(distributionsInput, prepareDistributionInput, userEnv = userEnv)
          names(distributionsInput) <- nms
        }
           
        if(exists('distributions', nimbleUserNamespace, inherits = FALSE)) {
            nimbleUserNamespace$distributions$add(distributionsInput)
        } else 
            nimbleUserNamespace$distributions <- distributionsClass(distributionsInput)
        virtualNodeFunctionDefinitions <- ndf_createVirtualNodeFunctionDefinitionsList(userAdded = TRUE)
        createNamedObjectsFromList(virtualNodeFunctionDefinitions, envir = .GlobalEnv)

    # note don't use rFunHandler as rUserDist nimbleFunction needs n as first arg so it works on R side, therefore we have n in the C version of the nimbleFunction and don't want to strip it out in Cpp generation
      }
    invisible(NULL)
        
}


#' Remove user-supplied distributions from use in NIMBLE BUGS models
#'
#' Deregister distributional information originally supplied by the user
#' for use in BUGS model code
#'
#' @param distributionsNames a character vector giving the names of the distributions to be dergistered
#' @author Christopher Paciorek
#' @export
deregisterDistributions <- function(distributionsNames) {
    if(!exists('distributions', nimbleUserNamespace, inherits = FALSE)) 
        cat("No user-supplied distributions are registered\n")
    matched <- distributionsNames %in% getAllDistributionsInfo('namesVector', userOnly = TRUE)
    if(sum(matched)) {
        distsMatched <- paste0(distributionsNames[matched], collapse = ', ')
        cat(paste0("Deregistering ", distsMatched, " from user-registered distributions\n"))
    }
    if(sum(!matched))
        for(nm in distributionsNames[!matched]) {
            cat(paste0("Cannot deregister ", nm, " as it is not registered as a user-defined distribution\n"))
        }
    
    distributionsNames <- distributionsNames[matched]
    if(length(distributionsNames)) {
        if(sum(!nimbleUserNamespace$distributions$namesVector %in% distributionsNames)) {
            sapply(distributionsNames, function(x) nimbleUserNamespace$distributions$remove(x))
        } else {  # all distributions to be removed
              rm(distributions, envir = nimbleUserNamespace)
        }
    }
    invisible(NULL)
}
    
#####################################################################################################
#####################################################################################################
#####  API for accessing info about distributions ###################################################
#####################################################################################################
#####################################################################################################


getDistributionList <- function(dists) {
    boolNative <- dists %in% distributions$namesVector
    if(all(boolNative)) return(distributions[dists])
    missingDists <- dists[!boolNative]
    allFound <- FALSE
    if(exists('distributions', nimbleUserNamespace, inherits = FALSE)) {
        if(all(missingDists %in% nimbleUserNamespace$distributions$namesVector))
            allFound <- TRUE
    }
    if(allFound) {
        ans <- vector('list', length(dists))
        ans[boolNative] <- distributions[dists[boolNative]]
        ans[!boolNative] <- nimbleUserNamespace$distributions[missingDists]
        return(ans)
    }
    notFound <- missingDists[ !(missingDists %in% nimbleUserNamespace$distributions$namesVector) ]
    stop(paste0('In getDistributions, distributions named ', paste(notFound, sep = ',', collapse = ","), ' could not be found.')) 
}

# note that getDimension and isDiscrete are not included as aliases below because they have the same name as modelBaseClass methods so we are having help for them direct to help(modelBaseClass) as we expect more usage of the modelBaseClass methods by users

#' Get information about a distribution
#'
#' Give information about each BUGS distribution
#'
#' @name distributionInfo
#' @aliases isUserDefined pqDefined getType getParamNames getDistributionInfo
#' 
#' @param dist a character vector of length one, giving the name of the distribution (as used in BUGS code), e.g. \code{'dnorm'}
#'
#' @param params an optional character vector of names of parameters for which dimensions are desired (possibly including \'value\' and alternate parameters)
#'
#' @param valueOnly a logical indicating whether to only return the dimension of the value of the node
#'
#' @param includeParams a logical indicating whether to return dimensions of parameters. If TRUE and \'params\' is NULL then dimensions of all parameters, including the dimension of the value of the node, are returned
#'
#' @param includeValue a logical indicating whether to return the string 'value', which is the name of the node value
#'
#' @author Christopher Paciorek
#' @details
#' NIMBLE provides various functions to give information about a BUGS distribution. In some cases, functions of the same name and similar functionality operate on the node(s) of a model as well (see \code{help(modelBaseClass)}).
#' 
#' \code{getDistributionInfo} returns an internal data structure (a reference class object) providing various information about the distribution. The output is not very user-friendly, but does contain all of the information that NIMBLE has about the distribution.
#'
#' \code{isDiscrete} tests if a BUGS distribution is a discrete distribution.
#'
#' \code{isUserDefined} tests if a BUGS distribution is a user-defined distribution.
#'
#' \code{pqAvail} tests if a BUGS distribution provides distribution ('p') and quantile ('q') functions.
#' 
#' \code{getDimension} provides the dimension of the value and/or parameters of a BUGS distribution. The return value is a numeric vector with an element for each parameter/value requested.
#'
#' \code{getType} provides the type (numeric, logical, integer) of the value and/or parameters of a BUGS distribution. The return value is a character vector with an element for each parameter/value requested. At present, all quantities are stored as numeric (double) values, so this function is of little practical use but could be exploited in the future.
#'
#' \code{getParamNames} provides the value and/or parameter names of a BUGS distribution.
#' 
#' @examples
#' distInfo <- getDistributionInfo('dnorm')
#' distInfo
#' distInfo$range
#'
#' isDiscrete('dbin')
#' 
#' isUserDefined('dbin')
#' 
#' pqDefined('dgamma')
#' pqDefined('dmnorm')
#' 
#' getDimension('dnorm')
#' getDimension('dnorm', includeParams = TRUE)
#' getDimension('dnorm', c('var', 'sd'))
#' getDimension('dcat', includeParams = TRUE)
#' getDimension('dwish', includeParams = TRUE)
#' 
#' getType('dnorm')
#' getType('dnorm', includeParams = TRUE)
#' getType('dnorm', c('var', 'sd'))
#' getType('dcat', includeParams = TRUE)
#' getType('dwish', includeParams = TRUE)
#'
#' getParamNames('dnorm', includeValue = FALSE)
#' getParamNames('dmnorm')
#'
NULL

#' @rdname distributionInfo
#' @export
getDistributionInfo <- function(dist) {
    if(is.na(dist)) return(NA)
    ans <- distributions[[dist]]
    if(!is.null(ans)) return(ans)
    ##    if(dist %in% distributions$namesVector) return(distributions[[dist]])
    ans <- nimbleUserNamespace$distributions[[dist]]
    if(!is.null(ans)) return(ans)
    ##if(exists('distributions', nimbleUserNamespace, inherits = FALSE) && dist %in% nimbleUserNamespace$distributions$namesVector)
    ##    return(nimbleUserNamespace$distributions[[dist]])
    stop(paste0("getDistributionInfo: ", dist, " is not a distribution provided by NIMBLE or supplied by the user."))
}

getAllDistributionsInfo <- function(kind, nimbleOnly = FALSE, userOnly = FALSE) {
    if(kind %in% c('namesVector', 'namesExprList', 'translations')) {
        if(userOnly) out <- NULL else out <- get(kind, distributions)
        if(!nimbleOnly && exists('distributions', nimbleUserNamespace, inherits = FALSE))
            out <- c(out, get(kind, nimbleUserNamespace$distributions))
        return(out)
    }

if(kind %in% c('pqAvail', 'discrete')) {
        if(userOnly) out <- NULL else out <- sapply(distributions$distObjects, '[[', kind)
        if(!nimbleOnly && exists('distributions', nimbleUserNamespace, inherits = FALSE))
            out <- c(out, sapply(nimbleUserNamespace$distributions$distObjects, '[[', kind))
        return(out)
    }
    stop(paste0("getAllDistributionInfo: ", kind, " is not available from the distributions information."))
}

evalInDistsMatchCallEnv <- function(expr) {
    dist <- as.character(expr[[1]])
    if(dist %in% distributions$namesVector)
        return(eval(expr, distributions$matchCallEnv))
    if(exists('distributions', nimbleUserNamespace, inherits = FALSE) &&
       dist %in% nimbleUserNamespace$distributions$namesVector)
        return(eval(expr, nimbleUserNamespace$distributions$matchCallEnv))
    stop(paste0("evalInDistsMatchCallEnv: ", dist, " is not a distribution provided by NIMBLE or supplied by the user."))
}

stripPrefix <- function(vec, prefix = "d")
    return(gsub(paste0("^", prefix), "", vec))

BUGSdistToRdist <- function(BUGSdists, dIncluded = FALSE) {
    Rdists <- lapply(getAllDistributionsInfo('translations'), `[[`, 1)
    if(!dIncluded) names(Rdists) <- stripPrefix(names(Rdists))
    results <- unlist(Rdists[BUGSdists])
    names(results) <- NULL
    if(!dIncluded) return(stripPrefix(results)) else return(results)
}

#' @export
isDiscrete <- function(dist) {
    if(is.na(dist)) return(NA)
    if(length(dist) > 1 || !inherits(dist, 'character'))
        stop("isDiscrete: 'dist' should be a character vector of length 1")
    return(getDistributionInfo(dist)$discrete)
}

#' @rdname distributionInfo
#' @export 
isUserDefined <- function(dist) {
    if(is.na(dist)) return(dist)
    if(length(dist) > 1 || !inherits(dist, 'character'))
        stop("isUserDistribution: 'dist' should be a character vector of length 1")
    if(exists('distributions', nimbleUserNamespace, inherits = FALSE) &&
       dist %in% getAllDistributionsInfo('namesVector', userOnly = TRUE))
      return(TRUE) else return(FALSE)
}

#' @rdname distributionInfo
#' @export
pqDefined <- function(dist) {
    if(is.na(dist)) return(NA)
    if(length(dist) > 1 || !inherits(dist, 'character'))
        stop("pqDefined: 'dist' should be a character vector of length 1")
   return(getDistributionInfo(dist)$pqAvail)
} 

## not user-facing. only for use in model$checkBasics(),
## to avoid "same size check" for distribution parameters
isMixedSizes <- function(dist) {
    if(is.na(dist)) return(NA)
    if(length(dist) > 1 || !inherits(dist, 'character'))
        stop("isMixedSizes: 'dist' should be a character vector of length 1")
   return(getDistributionInfo(dist)$mixedSizes)
}

#' @export
getDimension <- function(dist, params = NULL, valueOnly = is.null(params) &&
                         !includeParams, includeParams = !is.null(params)) {
    if(length(dist) == 1 && is.na(dist)) return(NA)  # in case of passing a determ node
    if(length(dist) > 1 || !inherits(dist, 'character'))
      stop("getDimension: 'dist' should be a character vector of length 1")
  distInfo <- getDistributionInfo(dist)
  
  if(!includeParams && !valueOnly)
    stop("getDimension: no parameters or value requested")
  if(valueOnly && (!is.null(params) || includeParams))
    stop("getDimension: 'valueOnly' cannot be TRUE if parameters also requested")
  if(!includeParams && !is.null(params))
    stop("getDimension: 'params' is not NULL but 'includeParams' is FALSE")
  if(valueOnly) {
    params <- 'value'
  } else {  
    if(includeParams && is.null(params)) 
      params <- getParamNames(dist, includeValue = TRUE)
  }
  notFound <- which(! params %in% getParamNames(dist))
  if(length(notFound)) {
    if('x' %in% params[notFound]) warning("getDimension: use 'value' instead of 'x'.")
    stop("getDimension: these parameter names not found: ", params[notFound])
  }
  out <- sapply(params, function(p) distInfo$types[[p]]$nDim)
  return(out)
}

getParamID <- function(dist, params = NULL, valueOnly = is.null(params) &&
                       !includeParams, includeParams = !is.null(params)) {
    if(length(dist) == 1 && is.na(dist)) return(NA)
    if(length(dist) > 1 || !inherits(dist, 'character'))
    stop("getType: 'dist' should be a character vector of length 1")
  distInfo <- getDistributionInfo(dist)
  
  if(!includeParams && !valueOnly)
    stop("getDimension: no parameters or value requested")
  if(valueOnly && (!is.null(params) || includeParams))
    stop("getDimension: 'valueOnly' cannot be TRUE if parameters also requested")
  if(!includeParams && !is.null(params))
    stop("getDimension: 'params' is not NULL but 'includeParams' is FALSE")
  if(valueOnly) {
    params <- 'value'
  } else {  
    if(includeParams && is.null(params)) 
      params <- getParamNames(dist, includeValue = TRUE)
  }
  notFound <- which(! params %in% getParamNames(dist))
  if(length(notFound)) {
    if('x' %in% params[notFound]) warning("getParamID: use 'value' instead of 'x'.")
    stop("getParamID: these parameter names not found: ", params[notFound])
  }
  out <- distInfo$paramIDs[params]
  return(out)
}

#' @rdname distributionInfo
#' @export
getType <- function(dist, params = NULL, valueOnly = is.null(params) &&
                       !includeParams, includeParams = !is.null(params)) {
    if(length(dist) == 1 && is.na(dist)) return(NA)
    if(length(dist) > 1 || !inherits(dist, 'character'))
        stop("getType: 'dist' should be a character vector of length 1")
    distInfo <- getDistributionInfo(dist)
    
  if(!includeParams && !valueOnly)
    stop("getType: no parameters or value requested")
  if(valueOnly && (!is.null(params) || includeParams))
    stop("getType: 'valueOnly' cannot be TRUE if parameters also requested")
  if(!includeParams && !is.null(params))
    stop("getType: 'params' is not NULL but 'includeParams' is FALSE")
  if(valueOnly) {
    params <- 'value'
  } else {  
    if(includeParams && is.null(params)) 
      params <- getParamNames(dist, includeValue = TRUE)
  }
    notFound <- which(! params %in% getParamNames(dist))
    if(length(notFound)) {
        if('x' %in% params[notFound]) warning("getParamID: use 'value' instead of 'x'.")
        stop("getType: these parameter names not found: ", params[notFound])
    }
    out <- sapply(params, function(p) distInfo$types[[p]]$type)
    return(out)
}

# perhaps have args to allow only reqdArgs or only altParams?

#' @rdname distributionInfo
#' @export
getParamNames <- function(dist, includeValue = TRUE) {
    if(length(dist) == 1 && is.na(dist)) return(NA)
    if(length(dist) > 1 || !inherits(dist, 'character'))
        stop("getParamNames: 'dist' should be a character vector of length 1")
    distInfo <- getDistributionInfo(dist)
    names <- names(distInfo$paramIDs)
    if(!includeValue)
        names <- names[!names == 'value']
    return(names)
}

#####################################################################################################
#####################################################################################################
#####  executable code, creates global system variable 'distributions' and 'distribution_aliases  ###
#####################################################################################################
#####################################################################################################

distributions <- distributionsClass(distributionsInputList)

# removed by CJP as getDistribution() and getDistributionsInfo() make it unneeded
# getDistributionsObject <- function() {
#   distributions
# }

processDistributionAliases <- function(distributionsInputList) {
    tmp <- sapply(distributionsInputList, function(x) if(length(x$alias)) x$alias else NULL)
    # next two lines avoid need for regex processing if we used unlist() when a dist has multiple aliases
    aliases <- rep(names(tmp), sapply(tmp, length))
    names(aliases) <- unlist(tmp, use.names = FALSE) 
   
    return(aliases)
}
                  
distributionAliases <- processDistributionAliases(distributionsInputList) 


distribution_dFuns <- BUGSdistToRdist(getAllDistributionsInfo('namesVector'), dIncluded = TRUE)
distribution_rFuns <- gsub("^d", "r", distribution_dFuns)

pqAvail <- names(which(getAllDistributionsInfo('pqAvail')))
pqDists <- BUGSdistToRdist(pqAvail, dIncluded = TRUE)

distribution_pFuns <- gsub("^d", "p", pqDists)
distribution_qFuns <- gsub("^d", "q", pqDists)

distributionFuns <- c(distribution_dFuns, distribution_rFuns, distribution_pFuns, distribution_qFuns)

## following sections are added for use in genCpp_operatorLists and other places.  Slightly different need is to have separate list of scalar distributions and to use Rdist names

## dCRP is causing warnings in 'make man' though doesn't seem to cause errors in build package, but filter out dCRP to avoid the warnings.
nms <- getAllDistributionsInfo('namesVector')
nms <- nms[nms != 'dCRP']

scalar_distribution_bool <- unlist(lapply(nms, function(x) all(unlist(lapply(getDistributionInfo(x)$types, function(y) y$nDim == 0 )))))
scalar_distribution_dFuns <- BUGSdistToRdist(nms[scalar_distribution_bool], dIncluded = TRUE)
scalar_distribution_rFuns <- gsub("^d", "r", scalar_distribution_dFuns)

scalar_pqAvail_bool <- nimble:::getAllDistributionsInfo('pqAvail')[nms] & scalar_distribution_bool
scalar_pqAvail_dFuns <- BUGSdistToRdist(nms[scalar_pqAvail_bool], dIncluded = TRUE)
scalar_distribution_pFuns <- gsub("^d", "p", scalar_pqAvail_dFuns)
scalar_distribution_qFuns <- gsub("^d", "q", scalar_pqAvail_dFuns)

rm(nms, scalar_distribution_bool, scalar_pqAvail_bool, scalar_pqAvail_dFuns)


