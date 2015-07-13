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
    round = 'nimbleRound',
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
    } else {
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
        if(length(typeDeclExpr$args) > 1) {
            typeSizeExpr <- typeDeclExpr$args[[2]]
            if(is.numeric(typeSizeExpr) & nDim == 1) sizeExprs <- list(typeSizeExpr)
            else if(nDim ==1 & typeSizeExpr$name != 'c') sizeExprs <- list(typeSizeExpr)
            else sizeExprs <- typeSizeExpr$args
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
    ## If sizeExprs have been provided, then, we generate a sizeSize call.  This will then over-ride the typeEnv entry
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
