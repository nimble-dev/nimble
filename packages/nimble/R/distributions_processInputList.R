

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
              dupl <- which(nms %in% getDistributionsInfo('namesVector', userOnly = TRUE))
              if(length(dupl)) {
                  distObjects[dupl] <<- NULL
                  namesVector <<- namesVector[-dupl]
                  namesExprList[dupl] <<- NULL
                  translations[dupl] <<- NULL
                  cat("Overwriting the following user-supplied distributions:", nms[dupl], ".\n", sep = " ")
              }
              for(i in seq_along(dil))     distObjectsNew[[i]] <- distClass(dil[[i]], nms[i])
              names(distObjectsNew) <- nms
              translations <<- c(translations, lapply(distObjectsNew, function(d) c(d$densityName, d$simulateName)))

              distObjects <<- c(distObjects, distObjectsNew)
              namesVector <<- c(namesVector, nms)
              namesExprList <<- c(namesExprList, lapply(namesVector, as.name))
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
        types = 'ANY'		#'list',     ## named list (names are 'node', ALL reqdArgs, and ALL altParams), each element is a named list: list(type = 'double', nDim = 0) <- default values
        ### typesForVirtualNodeFunction = 'ANY'		#'list'  ## version of 'types' for making the virtualNodeFunction definiton.  same as above, except without 'value'
    ),
    
    methods = list(
        initialize = function(distInputList, BUGSdistName) {
            RdistExprList <<- list()
            alts <<- list()
            exprs <<- list()
            altParams <<- list()
            types <<- list()
            BUGSdistName <<- BUGSdistName
            BUGSdistExpr <<- parse(text=distInputList$BUGSdist)[[1]]
            if(BUGSdistExpr[[1]] != BUGSdistName)   stop(paste0('inconsistent BUGS distribution names for distribution: ', BUGSdistName))
            RdistTextVector <- if(is.null(distInputList$Rdist)) character() else distInputList$Rdist
            RdistExprList <<- lapply(RdistTextVector, function(t) parse(text=t)[[1]])
            numAlts <<- length(RdistExprList)
            init_altsExprsReqdArgs()
            simulateName <<- sub('^d', 'r', densityName)
            init_altParams(distInputList)
            discrete <<- if(is.null(distInputList$discrete))    FALSE    else    distInputList$discrete
            pqAvail <<- if(is.null(distInputList$pqAvail))    FALSE    else    distInputList$pqAvail
            init_types(distInputList)
        },
        
        init_altsExprsReqdArgs = function() {
            alts <<- list()
            exprs <<- list()
            if(numAlts == 0) {
                params <- as.list(BUGSdistExpr[-1])   # removes the distribution name
                paramsText <- lapply(params, deparse)
                reqdArgs <<- sapply(paramsText, function(pt) init_getReqdArgs(pt))
                densityName <<- as.character(BUGSdistExpr[[1]])
            } else {
                params <- lapply(RdistExprList, `[`, -1)        # removes the distribution names
                paramsText <- lapply(params, function(x) lapply(x, deparse))
                reqdArgsList <- lapply(paramsText, function(pt) init_getReqdArgs(pt))
                densityNamesList <- lapply(RdistExprList, function(expr) as.character(expr[[1]]))
                if(length(unique(reqdArgsList)) > 1)
                    stop('R/NIMBLE parameter names and order not consistent across alternative parameterizations')
                if(length(unique(densityNamesList)) > 1)
                    stop('R/NIMBLE density names not consistent across alternative parameterizations')
                reqdArgs <<- reqdArgsList[[1]]
                densityName <<- densityNamesList[[1]]
                for(i in seq_along(params)) {
                    boolNoDefault <- if (is.null(names(paramsText[[i]]))) seq_along(paramsText[[i]]) else names(paramsText[[i]]) == ''
                    if(sum(!boolNoDefault)) {
                        exprs[[i]] <<- lapply(params[[i]][!boolNoDefault], function(x) {names(x) <- NULL; x})
                        BUGSargs <- unique(unlist(c(lapply(exprs[[i]], all.vars), paramsText[[i]][boolNoDefault])))
                        names(BUGSargs) <- NULL
                        if(!identical(sort(BUGSargs), sort(reqdArgs))) alts[[i]] <<- BUGSargs
                    }
                }
            }
        },
        
        init_getReqdArgs = function(x) {
            args <- if(is.null(names(x))) rep('', length(x)) else names(x)
            args[args == ''] <- unlist(x[args == ''])
            return(args)
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
                if(!(typeList$nDim %in% 0:1000))     stop(paste0('unknown nDim specified in distribution: ', typeList$nDim))  ## yes, specificying maximum dimension of 1000
                types[[typeName]] <<- typeList
            }
        },
        
        init_types_makeArgList = function(typeArgCharVector) {
            parsedArgList <- lapply(typeArgCharVector, function(x) parse(text=x)[[1]])
            allNames <- unlist(lapply(parsedArgList, function(pa) as.character(pa[[2]])))
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

checkDistributionsInput <- function(distributionsInput) {
    allowedFields <- unique(unlist(sapply(nimble:::distributionsInputList, names)))
    if(sum(!names(distributionsInput) %in% allowedFields)) 
        stop(paste0(names(distributionsInput), " has unknown field."))
    if(!sum(is.character(distributionsInput$BUGSdist))) stop(paste0(distributionsInput$BUGSdist, ": field 'BUGSdist' is not of type character."))
    if(exists("Rdist", distributionsInput) && !sum(is.character(distributionsInput$Rdist))) stop(paste0(distributionsInput$BUGSdist, ": field 'Rdist' is not type of character."))
    if(exists("discrete", distributionsInput) && !sum(is.logical(distributionsInput$discrete))) stop(paste0(distributionsInput$BUGSdist, ": field 'discrete' is not type logical."))
    if(exists("pqAvail", distributionsInput) && !sum(is.logical(distributionsInput$pqAvail))) stop(paste0(distributionsInput$BUGSdist, ": field 'pqAvail' is not of type logical."))
    if(exists("types", distributionsInput) && !sum(is.character(distributionsInput$types))) stop(paste0(distributionsInput$BUGSdist, ": field 'types' is not of type character."))
    if(exists("altParams", distributionsInput) && !sum(is.character(distributionsInput$altParams))) stop(paste0(distributionsInput$BUGSdist, ": field 'altParams' is not of type character."))
    if(length(distributionsInput$BUGSdist) > 1 || (exists('discrete', distributionsInputList) && length(distributionsInputList$discrete) > 1) || (exists('pqAvail', distributionsInputList) && length(distributionsInputList$pqAvail) > 1))
        stop(paste0(names(distributionsInput), " field 'BUGSdist', 'discrete', 'altParams', or 'pqAvail' is not of length one."))
    invisible(NULL)
}

checkDistributionsFunctions <- function(distributionsInput) {
    if(exists('Rdist', distributionsInput))
        inputString <- distributionsInput$Rdist else inputString <- distributionsInput$BUGSdist
    densityName <- as.character(parse(text = inputString)[[1]][[1]])
    simulateName <- sub('^d', 'r', densityName)
    if(!is.rcf(get(densityName)) || !is.rcf(get(simulateName)))
        stop(paste0("Either density or random generation functions for ", densityName,
                    " are not available as nimbleFunctions without setup code."))
    # FIXME: deal with finding by R's scoping rules here and in genCpp_sizeProcessing (currently line 139)
    if(!is.null(distributionsInput) && exists("pqAvail", distributionsInput) && distributionsInput$pqAvail) {
        cdfName <- sub('^d', 'p', densityName)
        quantileName <- sub('^d', 'q', densityName)
        if(!is.rcf(get(cdfName)) || !is.rcf(get(quantileName)))
            stop(paste0("Either distribution (CDF) or quantile (inverse CDF) functions for ", densityName,
                        " are not available as nimbleFunctions without setup code."))
    }
    invisible(NULL)
}


getMaxDim <- function(typeList) 
    max(sapply(typeList, '[[', 'nDim'))

getValueDim <- function(distObject) 
    distObject$types$value$nDim



#' Add user-supplied distributions for use in NIMBLE BUGS models
#'
#' Register distributional information so that NIMBLE can process
#' user-supplied distributions in BUGS model code
#'
#' @param distributionsInputList a list in the form of that shown in \code{nimble:::distributionsInputList} with required field \code{BUGSdist} and optional fields \code{Rdist}, \code{altParams}, \code{discrete}, \code{pqAvail}, \code{types}. See Details for more information. If only one distribution is supplied it may be a list rather than a list containing a list.
#' @author Christopher Paciorek
#' @export
#' @details
#' \itemize{
#' \item{\code{BUGSdist}} {
#' a character string in the form of the density name (starting with 'd') followed by the names of the parameters in parentheses. When alternative parameterizations are given in \code{Rdist}, this should be an exhaustive list of the unique parameter names from all possible parameterizations.
#' }
#' \item{\code{Rdist}} {
#' an optional character vector with one or more alternative specifications of the density; each alternative specification can be an alternative name for the density, a different ordering of the parameters, different parameter name(s), or an alternative parameterization. In the latter case, the character string in parentheses should provide a given reparameterization as comma-separated name = value pairs, one for each default parameter, where name is the name of the default parameter and value is a mathematical expression relating the default parameter to the alternative parameters or other default parameters. The default parameters should correspond to the input arguments of the nimbleFunctions provided as the density and random generation functions. The maethematical expression can use any of the math functions allowed in NIMBLE (see the User Manual) as well as user-supplied nimbleFunctions without setup code. 
#' }
#' \item{\code{discrete}} {
#' a optional logical indicating if the distribution is that of a discrete random variable. If not supplied, distribution is assumed to be for a continuous random variable.
#' }
#' \item{\code{pqAvail}} {
#' a optional logical indicating if distribution (CDF) and quantile (inverse CDF) functions are provided as nimbleFunctions. These are required for one to be able to use truncated versions of the distribution. Only applicable for univariate distributions. If not supplied, assumed to be FALSE.
#' }
#' \item{\code{altParams}} {
#' a character vector of comma-separated name = value pairs that provide the mathematical expressions relating non-default parameters to default parameters. These inverse functions are used for MCMC conjugacy calculations when a conjugate relationship is expressed in terms of non-default parameters (such as the precision for normal-normal conjugacy). If not supplied, the system will still functional but at a possible loss of efficiency in certain algorithms.
#' }
#' \item{\code{types}} {
#' a character vector of comma-separated 'name = input' pairs indicating the type and dimension of the random variable and parameters (including default and alternative parameters). 'input' should take the form 'integer(d)' or 'double(d)', where 'd' is 0 for scalars, 1 for vectors, 2 for matrices. 'name' should be either 'value' to indicate the random variable itself or the parameter name to indicate a given parameter.  Required for multivariate distributions or when parameters are not scalars. By default all parameters are assumed double(0) and values (the random variable) are assumed either double(0) or integer(0) depending on whether the distribution is for a discrete random variable.
#' }
registerDistributions <- function(distributionsInputList) {
    if(missing(distributionsInputList)) {
        cat("No distribution information supplied.\n")
    } else {
        if(!is.list(distributionsInputList[[1]])) stop("'distributionsInputList' should be a named list of lists")

        nms <- names(distributionsInputList)
        dupl <- nms[nms %in% getDistributionsInfo('namesVector', nimbleOnly = TRUE)]
        if(length(dupl)) {
            distributionsInputList[dupl] <- NULL
            cat("Ignoring the following user-supplied distributions as they have the same names as default NIMBLE distributions:", nms, ". Please rename to avoid the conflict.\n", sep = "")
        }
        sapply(distributionsInputList, checkDistributionsInput)
        sapply(distributionsInputList, checkDistributionsFunctions)

        if(exists('distributions', nimbleUserNamespace)) {
            nimbleUserNamespace$distributions$add(distributionsInputList)
        } else 
            nimbleUserNamespace$distributions <- distributionsClass(distributionsInputList)
    }
    virtualNodeFunctionDefinitions <- ndf_createVirtualNodeFunctionDefinitionsList(userAdded = TRUE)
    createNamedObjectsFromList(virtualNodeFunctionDefinitions, envir = .GlobalEnv)

    # note don't use rFunHandler as rUserDist nimbleFunction needs n as first arg so it works on R side, therefore we have n in the C version of the nimbleFunction and don't want to strip it out in Cpp generation

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
    if(!exists('distributions', nimbleUserNamespace)) 
        cat("No user-supplied distributions are registered.\n")
    matched <- distributionsNames %in% getDistributionsInfo('namesVector', userOnly = TRUE)
    if(sum(matched)) 
        cat(paste("Deregistering ", distributionsNames[matched], " from user-registered distributions.\n"))
    if(sum(!matched))
        cat(paste(distributionsNames[!matched], " are not user-registered distributions; ignoring.\n"))
    
    distributionsNames <- distributionsNames[matched]
    if(length(distributionsNames)) {
        if(sum(!nimbleUserNamespace$distributions$namesVector %in% distributionsNames)) {
            sapply(distributionsNames, function(x) nimbleUserNamespace$distributions$remove(x))
        } else {  # all distributions to be removed
              nimbleUserNamespace$distributions <- NULL
          }
    }
    return(NULL)
}
    
#####################################################################################################
#####################################################################################################
#####  API for accessing info about distributions ###################################################
#####################################################################################################
#####################################################################################################

# this could be more elegant, e.g., providing
# - isDiscrete(distName)
# - pqAvail(distName)
# - getDistributionNames()
# at the moment it is still somewhat tied to the internal structure of our distributionsClass
# - Chris

# this is a hack because having trouble calling getDistribution() from within nodeInfoClass$isDiscrete
getDistribution2 <- function(distName) {
    getDistribution(distName)
}

getDistribution <- function(distName) {
    if(distName %in% distributions$namesVector) return(distributions[[distName]])
    if(exists('distributions', nimbleUserNamespace) && distName %in% nimbleUserNamespace$distributions$namesVector)
        return(nimbleUserNamespace$distributions[[distName]])
    stop(paste0("getDistribution: ", distName, " is not a distribution provided by NIMBLE or supplied by the user."))
}

getDistributionsInfo <- function(kind, nimbleOnly = FALSE, userOnly = FALSE) {
    if(kind %in% c('namesVector', 'namesExprList', 'translations')) {
        if(userOnly) out <- NULL else out <- get(kind, distributions)
        if(!nimbleOnly && exists('distributions', nimbleUserNamespace))
            out <- c(out, get(kind, nimbleUserNamespace$distributions))
        return(out)
    }
    if(kind == 'pqAvail') {
        if(userOnly) out <- NULL else out <- sapply(distributions$distObjects, '[[', 'pqAvail')
        if(!nimbleOnly && exists('distributions', nimbleUserNamespace))
            out <- c(out, sapply(nimbleUserNamespace$distributions$distObjects, '[[', 'pqAvail'))
        return(out)
    }
    stop(paste0("getDistributionInfo: ", kind, " is not available from the distributions information."))
}

evalInDistsMatchCallEnv <- function(expr) {
    distName <- as.character(expr[[1]])
    if(distName %in% distributions$namesVector)
        return(eval(expr, distributions$matchCallEnv))
    if(exists('distributions', nimbleUserNamespace) &&
       distName %in% nimbleUserNamespace$distributions$namesVector)
        return(eval(expr, nimbleUserNamespace$distributions$matchCallEnv))
    stop(paste0("evalInDistsMatchCallEnv: ", distName, " is not a distribution provided by NIMBLE or supplied by the user."))
}

#####################################################################################################
#####################################################################################################
#####  executable code, creates global system variable 'distributions' ##############################
#####################################################################################################
#####################################################################################################

distributions <- distributionsClass(distributionsInputList)

getDistributionsObject <- function() {
    distributions
}






