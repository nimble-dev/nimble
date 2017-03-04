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
         nfMethod = 'nfMethodErrorHandler',
         min = 'minMaxHandler',
         max = 'minMaxHandler'),
    makeCallList(names(specificCallReplacements), 'replacementHandler'),
    makeCallList(c('nimNumeric', 'nimLogical', 'nimInteger', 'nimMatrix', 'nimArray'), 'nimArrayGeneralHandler' ),
    ##makeCallList(c(distribution_rFuns, 'rt', 'rexp'), 'rFunHandler'),  # exp and t allowed in DSL because in R and Rmath, but t_nonstandard and exp_nimble are the Nimble distributions for nodeFunctions
    makeCallList(c('dmnorm_chol', 'dmvt_chol', 'dwish_chol', 'dmulti', 'dcat', 'dinterval', 'ddirch'), 'dmFunHandler')
         )
specificCallHandlers[['rmnorm_chol']] <- 'rmFunHandler'
specificCallHandlers[['rmvt_chol']] <- 'rmFunHandler'
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
        if(code$args[[1]]$name == 'nimC') {
            newNames <- unlist(lapply(code$args[[1]]$args, `[[`, 'name'))
            if(length(newNames) == 0) stop('Error: no names provided to declare')
        } else stop('Error: first arg to declare should be a variable name or c(name1, name2, etc)')
    }
    for(i in seq_along(newNames)) {
        if(symTab$symbolExists(newNames[i])) stop(paste0('Error: declaring ',newNames[i],' which already exists.'))
    }
    typeDeclExpr <- code$args[[2]]
    ## a very patchy solution: go inside declare's processSpecificCalls handler and switch nimInteger back to integer
    if(typeDeclExpr$name == 'nimInteger') typeDeclExpr$name <- 'integer'
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
                if(typeSizeExpr$name == 'nimC')  ## it's a concatenation
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

            ## We used to put specific sizes into the symbol entry, which puts them into the typeEnv
            ## But now we keep the typeEnv entries generic to avoid unknown run-time evaluation path problems
            ## Therefore, no longer copy sizeExprs into newSizes.  the newSizes defaults set above will work in the symbol
            ## The sizeExprs are still needed below for the resize.
            
            ## for(i in 1:nDim) {
            ##     if(is.numeric(sizeExprs[[i]])) newSizes[i] <- sizeExprs[[i]]
            ## }
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


## process specific calls handler for nimArrayGeneral()
## the resulting function nimArrayGeneral() looks like:
## nimArrayGeneral(typeCharString, nDim, c(sizeExpr1, ...), initializeValue, initializeLogical)
## nimArrayGeneral(     arg1,      arg2,       arg3,              arg4,            arg5       )
nimArrayGeneralHandler <- function(code, symTab) {
    ## check to see if we're inside a declare statement()
    if(code$caller$name == 'declare') {
        if(code$name != 'nimInteger') stop('something going wrong with backwards compatibility fix for declare(..., integer())')
        code$name <- 'integer'
        code$args <- if(code$args[[2]] == 0) code$args[1] else code$args[1:2]
        return()
    }
    ## collectSizes is just a parse tree annotation, not a real function
    switch(code$name,
           ##nimNumeric(length = 0, value = 0, init = TRUE)
           nimNumeric = {
##               if(inherits(code$args[[1]], 'exprClass') && code$args[[1]]$isCall && code$args[[1]]$name == 'c') stop('numeric doesnt handle c() in length')
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[1], caller=code, callerArgID=3)
               newArgs <- list('double', 1, sizeExprs, code$args[[2]], code$args[[3]])
           },
           ##nimInteger(length = 0, value = 0, init = TRUE)
           nimInteger = {
##               if(inherits(code$args[[1]], 'exprClass') && code$args[[1]]$isCall && code$args[[1]]$name == 'c') stop('integer doesnt handle c() in length')
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[1], caller=code, callerArgID=3)
               newArgs <- list('integer', 1, sizeExprs, code$args[[2]], code$args[[3]])
           },
           ##nimLogical(length = 0, value = 0, init = TRUE)
           nimInteger = {
##               if(inherits(code$args[[1]], 'exprClass') && code$args[[1]]$isCall && code$args[[1]]$name == 'c') stop('integer doesnt handle c() in length')
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[1], caller=code, callerArgID=3)
               newArgs <- list('logical', 1, sizeExprs, code$args[[2]], code$args[[3]])
           },
           ##nimVector(type = 'double', length = 0, value = 0, init = TRUE)
           ##nimVector = {},
           ##nimMatrix(value = 0, nrow = 1, ncol = 1, init = TRUE, type = 'double')
           nimMatrix = { 
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[2:3], caller=code, callerArgID=3)
               newArgs <- list(code$args[[5]], 2, sizeExprs, code$args[[1]], code$args[[4]])
           },
           ##nimArray(value = 0, dim = c(1, 1), init = TRUE, type = 'double')
           nimArray = {
               ## nimArray will handle dim=c(...), and also dim=EXPR
               if(inherits(code$args[[2]], 'exprClass') && code$args[[2]]$isCall && code$args[[2]]$name == 'nimC') {
                   ## dim argument is c(...)
                   code$args[[2]]$name <- 'collectSizes'
                   newArgs <- list(code$args[[4]], length(code$args[[2]]$args), code$args[[2]], code$args[[1]], code$args[[3]])
               } else {
                   ## dim argument is a single number or expression
                   sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[2], caller=code, callerArgID=3)
                   newArgs <- list(code$args[[4]], 1, sizeExprs, code$args[[1]], code$args[[3]])
               }
           },
           stop('should never get here')
           )
    code$name <- 'nimArrayGeneral'
    code$args <- newArgs
    ## set caller / callerArgID for the new nimArraryGeneral() call arguments
    for(i in seq_along(code$args)) {
        if(inherits(code$args[[i]], 'exprClass')) {
            code$args[[i]]$callerArgID <- i
            code$args[[i]]$caller <- code
        }
    }
    ## set caller / callerArgID for the c(size expressions...) call
    for(i in seq_along(code$args[[3]]$args)) {
        if(inherits(code$args[[3]]$args[[i]], 'exprClass')) {
            code$args[[3]]$args[[i]]$callerArgID <- i
            code$args[[3]]$args[[i]]$caller <- code$args[[3]]
        }
    }
    if(!(code$args[[1]] %in% c('double', 'integer'))) stop('unknown type in nimArrayGeneral')
    if(code$args[[2]] != length(code$args[[3]]$args)) stop('mismatch between nDim and number of size expressions in nimArrayGeneral')
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
