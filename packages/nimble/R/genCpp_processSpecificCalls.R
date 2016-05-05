## Many cases handled here are simple replacements, e.g. ^ to pow.
## These could easily be done elsewhere, but it works cleanly to do them here.
##
## A more substantial task is to handle declare() statements. see declareHandler below.

####################################
## System to processSpecificCalls ##
####################################

specificCallReplacements <- list(
    '^' = 'pow',
    '%%' = 'nimMod',
    length = 'size',
    run = 'operator()',
    is.nan = 'ISNAN',
    is.nan.vec = 'ISNAN',
    is.na = 'ISNA',
    is.na.vec = 'ISNA',
    lgamma = 'lgammafn',
    logfact = 'lfactorial',
    loggam = 'lgammafn',
    gamma = 'gammafn',
    expit = 'ilogit',
    phi = 'iprobit',
    round = 'nimRound',
    ceiling = 'ceil',
    trunc = 'ftrunc',
    nimDim = 'dim',
    checkInterrupt = 'R_CheckUserInterrupt')
    

specificCallHandlers = c(
    list(invisible = 'removeLayerHandler',
         seq_along = 'seqAlongHandler',
         mvAccess = 'mvAccessHandler',
         numListAccess = 'mvAccessHandler',
         declare = 'declareHandler',
         nimArrayGeneral = 'nimArrayGeneralHandler',
         nfMethod = 'nfMethodErrorHandler',
         min = 'minMaxHandler',
         max = 'minMaxHandler'),
    makeCallList(names(specificCallReplacements), 'replacementHandler'),
    makeCallList(c(distribution_rFuns, 'rt', 'rexp'), 'rFunHandler'),  # exp and t allowed in DSL because in R and Rmath, but t_nonstandard and exp_nimble are the Nimble distributions for nodeFunctions
    makeCallList(c('dmnorm_chol', 'dwish_chol', 'dmulti', 'dcat', 'dinterval', 'ddirch'), 'dmFunHandler')
         )
specificCallHandlers[['rmnorm_chol']] <- 'rmFunHandler'
specificCallHandlers[['rwish_chol']] <- 'rmFunHandler'
specificCallHandlers[['rmulti']] <- 'rmFunHandler'
specificCallHandlers[['rcat']] <- 'rmFunHandler' ## not really multivar, but same processing
specificCallHandlers[['rdirch']] <- 'rmFunHandler'
specificCallHandlers[['rinterval']] <- 'rmFunHandler' ## not really multivar, but same processing

exprClasses_processSpecificCalls <- function(code, symTab) {
    if(code$isName) return(invisible())
    if(code$isCall) {
        for(i in seq_along(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
                exprClasses_processSpecificCalls(code$args[[i]], symTab)
            }
        }
        handler <- specificCallHandlers[[code$name]]
        if(!is.null(handler)) eval(call(handler, code, symTab))
    }
}

nfMethodErrorHandler <- function(code, symTab) {
    if(code$caller$name != 'chainedCall')
        stop(paste0('Error with ', nimDeparse(code), ': an nfMethod call always needs arguments. e.g. nfMethod(mf, \'foo\')().'), call. = FALSE)
    NULL
}

seqAlongHandler <- function(code, symTab) {
    code$name <- ':'
    insertExprClassLayer(code, 1, 'size')
    oldArg <- code$args[[1]]
    code$args <- vector('list', 2)
    code$args[[1]] <- 1
    setArg(code, 2, oldArg)
    NULL
}

## processes something like declare(Z, double(1, c(3, 4))) where the first argument to double is the number of dimensions and next (optional)
## give size expressions
## A difference between this and the syntax for run-time arguments is that here the sizes can be expressions.  For run-time arguments, they must be integers.  The expressions must refer to sizes of OTHER variables.  E.g. declare(Z, double(1, c(3, dim(A)[2])))
## What is not currently supported is a mix of known and unknown sizes, e.g. declare(Z, double(3, NA))
declareHandler <- function(code, symTab) {
    if(code$args[[1]]$isName) {
        newNames <- code$args[[1]]$name
    } else { ## case of declaring multiple names at once.  This is not documented in User Manual and not supported by R declare() 
        if(code$args[[1]]$name == 'c') {
            newNames <- unlist(lapply(code$args[[1]]$args, `[[`, 'name'))
            if(length(newNames) == 0) stop('Error: no names provided to declare')
        } else stop('Error: first arg to declare should be a variable name or c(name1, name2, etc)')
    }
    for(i in seq_along(newNames)) {
        if(symTab$symbolExists(newNames[i])) stop(paste0('Error: declaring ',newNames[i],' which already exists.'))
    }
    typeDeclExpr <- code$args[[2]]
    type <- typeDeclExpr$name
    if(length(typeDeclExpr$args) == 0) {
        nDim <- 0
        newSizes <- numeric(0)
    } else {
        nDim <- typeDeclExpr$args[[1]]
        newSizes <- rep(as.numeric(NA), nDim)
        if(length(typeDeclExpr$args) > 1) { ## sizes were provided
            typeSizeExpr <- typeDeclExpr$args[[2]]
            if(inherits(typeSizeExpr, 'exprClass')) { ## first size are is an expression
                if(typeSizeExpr$name == 'c')  ## it's a concatenation
                    sizeExprs <- typeSizeExpr$args ## record the args
                else { ## it's not a concatenation
                    if(nDim != 1) stop('confused in declareHandler')
                    sizeExprs <- list(typeSizeExpr) ## take single expression
                }
            } else { ## it better be numeric
                if(!is.numeric(typeSizeExpr)) stop('confused 2 in declareHandler')
                if(length(typeDeclExpr$args) != 1 + nDim) stop('confused 3 in declareHandler')
                sizeExprs <- list()
                for(i in 1:nDim)
                    sizeExprs[[i]] <- typeDeclExpr$args[[i+1]]
            }
            
            ## if(is.numeric(typeSizeExpr) & nDim == 1) sizeExprs <- list(typeSizeExpr)
            ## else if(nDim ==1 & typeSizeExpr$name != 'c') sizeExprs <- list(typeSizeExpr)
            ## else sizeExprs <- typeSizeExpr$args
            if(length(sizeExprs) != nDim) stop(paste('Error in declare for', paste(newNames, collapse = ','), ': wrong number of dimensions provided'))
            for(i in 1:nDim) {
                if(is.numeric(sizeExprs[[i]])) newSizes[i] <- sizeExprs[[i]]
            }
        }
    }
    ## The symbolTable object will be used in exprClasses_initSizes to create an entry in typeEnv when the symbol is first encountered
    for(i in seq_along(newNames)) {
        newSym <- symbolBasic(name = newNames[i], type = type, nDim = nDim, size = newSizes)
        symTab$addSymbol(newSym)
    }
    ## If sizeExprs have been provided, then, we generate a setSize call.  This will then over-ride the typeEnv entry
    if(length(typeDeclExpr$args) > 1) {
        newExprs <- vector('list', length(newNames))
        for(i in seq_along(newNames)) {
            nameExpr <- exprClass(name = newNames[i], isCall = FALSE, isName = TRUE, isAssign = FALSE, args = list())
            newExprs[[i]] <- exprClass(name = 'setSize', isCall = TRUE, isName = FALSE,
                                       isAssign = FALSE, args = c(list(nameExpr), sizeExprs),
                                       caller = code$caller, callerArgID = code$callerArgID)
            for(j in seq_along(newExprs[[i]]$args)) {
                if(inherits(newExprs[[i]]$args[[j]], 'exprClass')) {
                    newExprs[[i]]$args[[j]]$caller <- newExprs[[i]]
                    newExprs[[i]]$args[[j]]$callerArgID <- j
                }
            }
        }
        newBrackExpr <- newBracketExpr(newExprs)
        setArg(code$caller, code$callerArgID, newBrackExpr)
    } else {
        code$caller$args[[code$callerArgID]] <- exprClass(name = 'blank', isCall = TRUE, isName = FALSE, isAssign = FALSE, args = list(),
                                                          caller = code$caller, callerArgID = code$callerArgID)
    }
}


## handler for nimArrayGeneral()
## the prototype for nimArrayGeneral is *always*
## newName <- nimArrayGeneral(typeCharString, nDim,      c(sizeExpr1, ...), initialValue, initializeLogical)
## newName <- nimArrayGeneral(args[[1]],      args[[2]], args[[3]],         args[[4]],    args[[5]])
nimArrayGeneralHandler <- function(code, symTab) {
    assignmentCall <- code$caller
    newName <- assignmentCall$args[[1]]
    if(!newName$isName) stop('numeric, vector, matrix, array can only create complete variables (not part of an indexed variable)')
    newName <- newName$name
    if(symTab$symbolExists(newName)) stop(paste0('Error: declaring ',newName,' which already exists.'))
    type <- code$args[[1]]
    nDim <- code$args[[2]]
    sizeExprs <- as.list(code$args[[3]]$args)  ## strip the leading c() call, and get the rest as an actual R list
    if(length(sizeExprs) != nDim) stop('something we wrong processing in nimArrayGeneralHandler')
    newSizes <- rep(as.numeric(NA), nDim)
    for(i in 1:nDim)
        if(is.numeric(sizeExprs[[i]])) newSizes[i] <- sizeExprs[[i]]
    newSym <- symbolBasic(name = newName, type = type, nDim = nDim, size = newSizes)
    symTab$addSymbol(newSym)
    ## generate setSize() call
    nameExpr <- exprClass(name = newName, isCall = FALSE, isName = TRUE, isAssign = FALSE, args = list())
    setSizeExpr <- exprClass(name = 'setSize', isCall = TRUE, isName = FALSE,isAssign = FALSE, args = c(list(nameExpr), sizeExprs), caller = code$caller, callerArgID = code$callerArgID)
    for(j in seq_along(setSizeExpr$args)) {
        if(inherits(setSizeExpr$args[[j]], 'exprClass')) {
            setSizeExpr$args[[j]]$caller <- setSizeExpr
            setSizeExpr$args[[j]]$callerArgID <- j
        }
    }
    ## generate initialize() call
    nameExpr2 <- exprClass(name = newName, isCall = FALSE, isName = TRUE, isAssign = FALSE, args = list())
    initialValue <- code$args[[4]]
    initializeLogical <- code$args[[5]]
    ## NOTE: prototype for R call to initialize() is:
    ## initialize(nimArrName, initialValue, initializeLogical)
    initializeExpr <- exprClass(name = 'initialize', isCall = TRUE, isName = FALSE, isAssign = FALSE, args = list(nameExpr2, initialValue, initializeLogical), caller = code$caller, callerArgID = code$callerArgID)
    for(j in seq_along(initializeExpr$args)) {
        if(inherits(initializeExpr$args[[j]], 'exprClass')) {
            initializeExpr$args[[j]]$caller <- initializeExpr
            initializeExpr$args[[j]]$callerArgID <- j
        }
    }
    ## make setSize() and initialize() into a single bracketExpr, and set the caller / callerIDs
    newBrackExpr <- newBracketExpr(list(setSizeExpr, initializeExpr))
    ## taking your exprClass to new levels.....
    ## yes, this is a chained code$caller$caller lookup, to over-write the original <- assignment
    setArg(code$caller$caller, code$caller$callerArgID, newBrackExpr)
}

## for min(V), no change.  for min(v1, v2), change to pairmin(v1, v2)
minMaxHandler <- function(code, symTab) {
    if(length(code$args) == 2) code$name <- paste0('pair',code$name)
}

replacementHandler <- function(code, symTab) {
    repl <- specificCallReplacements[[code$name]]
    if(is.null(repl)) stop(paste0("No valid replacement for ", code$name), call. = FALSE)
    code$name <- repl
}
    
dmFunHandler <- function(code, symTab) {
    code$name <- paste0('nimArr_', code$name)
}

rmFunHandler <- function(code, symTab) {
    dmFunHandler(code, symTab)
    rFunHandler(code, symTab)
}

rFunHandler <- function(code, symTab) {
    ## strip the 1 from the first argument.  In the future we will need to condition on whether this is 1 or >1
    ## This is a funny step because building the R function *inserted* the 1 from the BUGS code
    notOK <- if(!is.numeric(code$args[[1]])) TRUE else code$args[[1]] != 1
    if(notOK) writeLines(paste('Warning: we currently expect to see a 1 and only a 1 as the first argument to an r[dist] call'))
    code$args[[1]] <- NULL
    for(i in seq_along(code$args)) {
        if(inherits(code$args[[i]], 'exprClass')) {
            code$args[[i]]$callerArgID <- code$args[[i]]$callerArgID - 1
        }
    }
}

removeLayerHandler <- function(code, symTab) {
    removeExprClassLayer(code)
}

mvAccessHandler <- function(code, symTab) {
	
    if(code$caller$caller$caller$name != '[') {
        insertExprClassLayer(code$caller$caller$caller, code$callerArgID, '[')
        code$caller$caller$caller$args[[2]] <- 1
    }
}
