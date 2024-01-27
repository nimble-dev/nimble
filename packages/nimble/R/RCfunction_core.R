# for use in DSL code check:
otherDSLcalls <- c("{",
                   "[[",
                   "$",
                   "resize",
                   "declare",
                   "returnType",
                   "seq_along",
                   "double",
                   "rankSample",
                   "new",
                   "nimEigen",
                   "nimSvd",
                   "nimOptim",
                   "nimIntegrate",
                   "nimOptimDefaultControl",
                   "nimDerivs",
                   "any_na",
                   "any_nan",
                   "void")

nimKeyWords <- list(copy = 'nimCopy',
                    print = 'nimPrint',
                    cat = 'nimCat',
                    step = 'nimStep',
                    equals = 'nimEquals',
                    dim = 'nimDim',
                    stop = 'nimStop',
                    numeric = 'nimNumeric',
                    logical = 'nimLogical',
                    integer = 'nimInteger',
                    matrix = 'nimMatrix',
                    array = 'nimArray',
                    round = 'nimRound',
                    c = 'nimC',
                    rep = 'nimRep',
                    seq = 'nimSeq',
                    eigen = 'nimEigen',
                    svd = 'nimSvd',
                    optim = 'nimOptim',
                    integrate = 'nimIntegrate',
                    optimDefaultControl = 'nimOptimDefaultControl',
                    min.bound = 'carMinBound',
                    max.bound = 'carMaxBound',
                    derivs = 'nimDerivs')

distsNotAllowedInAD <- c(
  paste0('d', c('interval', 'constraint'))
)

fxnsNotAllowedInAD <- c(
  '%%',
  'nimMod',
  'trace',
#  'nimStep',
  'nimEigen',
  'nimSvd',
  'nimOptim',
  'getParam',
  'getBound',
  'bessel_k',
  paste0('p', c('binom', 'nbinom', 'pois','beta','chisq',
                'dexp','exp_nimble','gamma','invgamma','lnorm',
                'logis','norm','t_nonstandard','unif','weibull', 't','exp')),
  paste0('q', c('binom', 'nbinom', 'pois','beta','chisq',
                'dexp','exp_nimble','gamma','invgamma','lnorm',
                'logis','norm','t_nonstandard','unif','weibull', 't','exp')),
  paste0('r', c('binom', 'nbinom', 'pois','beta','chisq',
                'dexp','exp_nimble','flat','halfflat','gamma','invgamma','sqrtinvgamma','lnorm',
                'logis','norm','t_nonstandard','unif','weibull', 't','exp')),
  paste0('r', c('cat', 'interval', 'car_normal', 'car_proper',
                'dirch','mnorm_chol','multi','mvt_chol','lkj_corr_cholesky','wish_chol',
                'invwish_chol')),
  paste0('nimArr_d', c('interval','constraint')),
  paste0('nimArr_r', c('mnorm_chol','mvt_chol', 'lkj_corr_cholesky','wish_chol',
                       'invwish_chol', 'car_normal','car_proper','multi','dirch') ),
  'getLogProb',
  'decide',
  'rankSample',
  'any_na', 
  'any_nan',
  'is.na',
  'is.nan',
  'nimCopy','carMinBound','carMaxBound'
)

callsNotAllowedInAD <- c(distsNotAllowedInAD, fxnsNotAllowedInAD)

nfMethodRCinterface <- setRefClass(
    Class = 'nfMethodRCinterface',
    fields = list(
        argInfo    = 'ANY',
        arguments  = 'ANY',
        returnType = 'ANY',
        uniqueName = 'character'))

nfMethodRC <- setRefClass(
    Class   = 'nfMethodRC',
    contains = 'nfMethodRCinterface',
    fields  = list(
        template   = 'ANY',
        code       = 'ANY',
        neededRCfuns = 'ANY',		#list
        externalHincludes = 'ANY',
        externalCPPincludes = 'ANY',
        buildDerivs = 'ANY'
    ),
    methods = list(
        initialize = function(method,
                              name,
                              check = FALSE,
                              methodNames = NULL,
                              setupVarNames = NULL,
                              buildDerivs = FALSE,
                              where = NULL) {
            ## uniqueName is only needed for a pure RC function.
            ## It is not needed for a nimbleFunction method.
            if(!missing(name)) uniqueName <<- name 
            neededRCfuns <<- list()	
            argInfo <<- formals(method)
             ## nf_changeNimKeywords changes all nimble keywords,
             ## e.g. 'print' to 'nimPrint'; see 'nimKeyWords' list at
             ## bottom
            code <<- nf_changeNimKeywords(body(method)) 
            if(code[[1]] != '{')
                code <<- substitute({CODE}, list(CODE=code))
            ## check all code except nimble package nimbleFunctions
            generateArgs()

            if(check && "package:nimble" %in% search()) {
                nf_checkDSLcode(code, methodNames, setupVarNames, names(arguments), where)
                if(isTRUE(nimbleOptions("doADerrorTraps"))) {
                   if(!isFALSE(buildDerivs) && !is.null(buildDerivs))
                       nf_checkDSLcode_derivs(code, names(arguments), callsNotAllowedInAD)
                }
            }
            generateTemplate() ## used for argument matching
            removeAndSetReturnType(check = check)
            ## Includes for .h and .cpp files when making external calls:
            ## If needed, these will be populated by nimbleExternalCall
            externalHincludes <<- externalCPPincludes <<- list()
            buildDerivs <<- buildDerivs
        },
        generateArgs = function() {
            argsList <- nf_createAList(names(argInfo))
            for(i in seq_along(argsList)) {
                if('default' %in% names(argInfo[[i]]))
                    argsList[[i]] <- argInfo[[i]]$default }
            arguments <<- as.pairlist(argsList)
        },
        generateTemplate = function() {
            functionAsList <- list(as.name('function'))
            functionAsList[2] <- list(NULL)
            if(!is.null(args)) functionAsList[[2]] <- arguments
            functionAsList[[3]] <- quote({})
            template <<- eval(as.call(functionAsList))
        },
        removeAndSetReturnType = function(check = TRUE) {
            ## check will be FALSE in the case of a virtualNimbleFunction method
            returnLineNum <- 0
            for(i in seq_along(code)) {
                if(length(code[[i]]) > 1) {
                    if(is.name(code[[i]][[1]])) {
                        if(code[[i]][[1]] == 'returnType') {
                            returnLineNum <- i
                            break;
                        }
                    }
                }
            }
            if(sum(all.names(code) == 'returnType') > 1)
                stop('multiple returnType() declarations in nimbleFunction method; only one allowed', call. = FALSE)
            if(returnLineNum == 0) {
                ## no returnType() declaration was found;
                ## create default behavior.
                returnTypeDeclaration <- quote(void())
            } else {
                ## returnType() declaration was found
                returnTypeDeclaration <- code[[returnLineNum]][[2]]
                ## Check that there is at least one return() statement
                if(check) {
                    if(!("return" %in% all.names(code))) {
                        writeLines("Problem in this code:")
                        writeLines(deparse(code))
                        stop('returnType() declaration was provided but there is no return() statement.  At least one return() statement must be provided', call. = FALSE)
                    }
                }
                code[returnLineNum] <<- NULL
            }
            ## a very patchy solution: switch nimInteger back to integer
            if(as.character(returnTypeDeclaration[[1]]) == 'nimInteger')
                returnTypeDeclaration[[1]] <- as.name('integer')
            if(as.character(returnTypeDeclaration[[1]]) == 'nimLogical')
                returnTypeDeclaration[[1]] <- as.name('logical')
            returnType <<- returnTypeDeclaration
        },
        generateFunctionObject = function(keep.nfMethodRC = FALSE, where = NULL) {
            functionAsList <- list(as.name('function'))
            functionAsList[2] <- list(NULL)
            if(!is.null(args)) functionAsList[[2]] <- arguments
            functionAsList[[3]] <- code
            ans <- eval(
                parse(text=deparse(as.call(functionAsList)),
                      keep.source = FALSE)[[1]])
            if(!is.null(where)) {
                environment(ans) <- new.env(parent = where)
            } else environment(ans) <- new.env()
            if(keep.nfMethodRC) {
                environment(ans)$nfMethodRCobject <- .self
            }
            ans
        },
        getArgInfo       = function() { return(argInfo)    },
        getReturnType    = function() { return(returnType) })
)

## Helper function that walks an exprClass tree and finds methods based on positioning wrt dollar signs
findMethodsInExprClass <- function(expr) {
    if(is(expr, 'exprClass')) {
        if(expr$isAssign) {
            return(findMethodsInExprClass(expr$args[[2]]))
        } else {
            if(expr$isCall) {
                if(expr$name == '$') {
                    ## Next check ensures that RHS of $ is a call and not a field
                    ## Check for whether it is first arg of chainedCall prevents finding
                    ## 'foo' in m$cc(m$foo) as a method.
                    if(!is.null(expr$caller) && expr$caller$name == 'chainedCall' &&
                       identical(expr, expr$caller$args[[1]]))
                        tmp <- expr$args[[2]]$name else tmp <- NULL
                    return(c(tmp, findMethodsInExprClass(expr$args[[1]])))
                } else {
                    return(unlist(lapply(expr$args, findMethodsInExprClass)))
                }
            } 
        }
    } 
    return(NULL)
}

nf_checkDSLcode_derivs <- function(code, args, calls_not_allowed) {
    calls <- setdiff(all.names(code),
                     c(all.vars(code), args))
    problem_calls <- calls[ calls %in% calls_not_allowed ]
    if(length(problem_calls)) {
        message("  [Note] Detected use of function(s) that are not supported for derivative tracking in a function or method for which `buildDerivs` has been requested: ", paste(unique(problem_calls), collapse = ", "), ".")
    }
    NULL
}

nf_checkDSLcode_buildDerivs <- function(code, buildDerivs) {
    code <- body(code)
    codeNames <- all.names(code)
    derivsLocn <- which(codeNames %in% c('derivs', 'nimDerivs'))
    if(length(derivsLocn)) {
        for(i in seq_along(derivsLocn)) {
            if(!(length(codeNames) >= derivsLocn[i]+3 && codeNames[derivsLocn[i]+1] == '$'
                && codeNames[derivsLocn[i]+3] == 'calculate')) {
                methodName <- codeNames[derivsLocn[i]+1]
                if(isFALSE(buildDerivs) || !length(buildDerivs) || is.null(buildDerivs) ||
                    (is.character(buildDerivs) && !methodName %in% buildDerivs) ||
                   (is.list(buildDerivs) && !methodName %in% names(buildDerivs)))
                    message("  [Note] Detected use of `nimDerivs` with a function or method, `", methodName, "`, for which `buildDerivs` has not been set. This nimbleFunction cannot be compiled.") 
            }

        }
    }
    invisible(NULL)
}

nf_checkDSLcode <- function(code, methodNames, setupVarNames, args, where = NULL) {
    validCalls <- c(names(sizeCalls),
                    otherDSLcalls,
                    names(specificCallReplacements),
                    nimKeyWords,
                    methodNames,
                    setupVarNames)
    calls <- setdiff(all.names(code),
                     c(all.vars(code), args))
    
    ## Find the 'y' in cases of x$y() and x[]$y() and x[[]]$y().

    nfMethods <- findMethodsInExprClass(RparseTree2ExprClasses(code))
    
    ## don't check RHS of $ to ensure it is a valid nf method because no current way to easily find the methods of nf's defined in setup code
    nonDSLcalls <- calls[!(calls %in% c(validCalls, nfMethods))]
    if(length(nonDSLcalls)) {
        objInR <- sapply(nonDSLcalls, exists, where = where)
        nonDSLnonR <- nonDSLcalls[!objInR]
        nonDSLinR <- nonDSLcalls[objInR]
        if(length(nonDSLinR))
            ## nf and nimbleFunctionList cases probably will never
            ## occur as these need to be passed as setup args or
            ## created in setup
            ##
            ## problem with passing inputIsName when run through roxygen...
            nonDSLinR <-
                nonDSLinR[!(sapply(nonDSLinR,
                                   function(x)
                                       is.nf(x, inputIsName = TRUE, where = where)) |
                            sapply(nonDSLinR,
                                   function(x)
                                       is(get(x, envir = where), 'nimbleFunctionList')) |
                            sapply(nonDSLinR,
                                   function(x)
                                       is.rcf(x, inputIsName = TRUE, where = where)) |
                            sapply(nonDSLinR,
                                   function(x)
                                       is.nlGenerator(x, inputIsName = TRUE, where = where)))]
        if(length(nonDSLinR))
            message("  [Note] Detected possible use of R functions in nimbleFunction run code.\n         For this nimbleFunction to compile, '", paste(nonDSLinR, collapse = ', '), "' must be defined as a nimbleFunction, nimbleFunctionList, or nimbleList.")
        if(length(nonDSLnonR))
            message("  [Note] For this nimbleFunction to compile, '", paste(nonDSLnonR, collapse = ', '), "' must be defined as a nimbleFunction, nimbleFunctionList, or nimbleList before compilation.")
    }
    return(0)
}

nf_changeNimKeywords <- function(code){
    if(length(code) > 0){
        for(i in seq_along(code) ) {
            if(is.call(code) ) {
                if(!is.null(code[[i]]) ) {
                    code[[i]] <- nf_changeNimKeywordsOne(code[[i]])
                }
            }
        }
    }
    return(code)
}


nf_changeNimKeywordsOne <- function(code, first = FALSE){
    if(length(code) == 1){
        if(as.character(code) %in% names(nimKeyWords)) {
            if(is.call(code)) {
                code[[1]] <- as.name( nimKeyWords[[as.character(code)]] )
            } else {
                if(!is.character(code) & first)
                    code <- as.name( nimKeyWords[[as.character(code)]] )
            }
        }
    }
    else if(length(code) > 1){ 
        for(i in seq_along(code) ) {
            if(!is.null(code[[i]]) )
                code[[i]] <- nf_changeNimKeywordsOne(code[[i]], first = i == 1)
        }
    }
    return(code)
}
