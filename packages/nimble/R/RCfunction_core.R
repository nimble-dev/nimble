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
                    optimDefaultControl = 'nimOptimDefaultControl',
                    min.bound = 'carMinBound',
                    max.bound = 'carMaxBound',
                    derivs = 'nimDerivs')

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
        externalCPPincludes = 'ANY'
    ),
    methods = list(
        initialize = function(method,
                              name,
                              check = FALSE,
                              methodNames = NULL,
                              setupVarNames = NULL) {
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

            if(check && "package:nimble" %in% search()) 
                nf_checkDSLcode(code, methodNames, setupVarNames, names(arguments))

            generateTemplate() ## used for argument matching
            removeAndSetReturnType(check = check)
            ## Includes for .h and .cpp files when making external calls:
            ## If needed, these will be populated by nimbleExternalCall
            externalHincludes <<- externalCPPincludes <<- list() 
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

nf_checkDSLcode <- function(code, methodNames, setupVarNames, args) {
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
        objInR <- sapply(nonDSLcalls, exists)
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
                                       is.nf(x, inputIsName = TRUE)) |
                            sapply(nonDSLinR,
                                   function(x)
                                       is(get(x), 'nimbleFunctionList')) |
                            sapply(nonDSLinR,
                                   function(x)
                                       is.rcf(x, inputIsName = TRUE)) |
                            sapply(nonDSLinR,
                                   function(x)
                                       is.nlGenerator(x, inputIsName = TRUE)))]
        if(length(nonDSLinR))
            warning(paste0("Detected possible use of R functions in nimbleFunction run code. For this nimbleFunction to compile, these objects must be defined as nimbleFunctions, nimbleFunctionLists, or nimbleLists: ", paste(nonDSLinR, collapse = ', '), "."), call. = FALSE)
        if(length(nonDSLnonR))
            warning(paste0("For this nimbleFunction to compile, these objects must be defined as nimbleFunctions, nimbleFunctionLists, or nimbleLists before compilation: ", paste(nonDSLnonR, collapse = ', '), "."), call. = FALSE)
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
