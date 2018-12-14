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
    is.nan = 'ISNAN',
    is.nan.vec = 'ISNAN',
    is.na = 'ISNA',
    is.na.vec = 'ISNA',
    lgamma = 'lgammafn',
    logfact = 'lfactorial',
    loggam = 'lgammafn',
    besselK = 'bessel_k',
    gamma = 'gammafn',
    expit = 'ilogit',
    phi = 'iprobit',
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
         min = 'minMaxHandler',
         max = 'minMaxHandler',
         nimSvd = 'svdHandler'),
    makeCallList(names(specificCallReplacements), 'replacementHandler'),
    makeCallList(c('nimNumeric', 'nimLogical', 'nimInteger', 'nimMatrix', 'nimArray'), 'nimArrayGeneralHandler' ),
    makeCallList(c('dmnorm_chol', 'dmvt_chol', 'dwish_chol', 'dinvwish_chol', 'dcar_normal', 'dcar_proper', 'dmulti', 'dcat', 'dinterval', 'ddirch'), 'dmFunHandler')
         )
specificCallHandlers[['rmnorm_chol']] <- 'rmFunHandler'
specificCallHandlers[['rmvt_chol']] <- 'rmFunHandler'
specificCallHandlers[['rwish_chol']] <- 'rmFunHandler'
specificCallHandlers[['rinvwish_chol']] <- 'rmFunHandler'
specificCallHandlers[['rcar_normal']] <- 'rmFunHandler'
specificCallHandlers[['rcar_proper']] <- 'rmFunHandler'
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

seqAlongHandler <- function(code, symTab) {
    code$name <- ':'
    insertExprClassLayer(code, 1, 'size')
    oldArg <- code$args[[1]]
    code$args <- vector('list', 2)
    code$args[[1]] <- 1
    setArg(code, 2, oldArg)
    NULL
}

svdHandler <- function(code, symTab){
  code$args[[2]] <- switch(tolower(code$args[[2]]),
                         none = 0,
                         thin = 1,
                         full = 2)
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
            if(length(sizeExprs) != nDim) stop(paste('Error in declare for', paste(newNames, collapse = ','), ': wrong number of dimensions provided'))
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
            nameExpr <- exprClass$new(name = newNames[i], isCall = FALSE, isName = TRUE, isAssign = FALSE, args = list())
            newExprs[[i]] <- exprClass$new(name = 'setSize', isCall = TRUE, isName = FALSE,
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
        code$caller$args[[code$callerArgID]] <- exprClass$new(name = 'blank', isCall = TRUE, isName = FALSE, isAssign = FALSE, args = list(),
                                                          caller = code$caller, callerArgID = code$callerArgID)
    }
}


## process specific calls handler for nimArrayGeneral()
## the resulting function nimArrayGeneral() looks like:
## nimArrayGeneral(typeCharString, nDim, c(sizeExpr1, ...), initializeValue, initializeLogical, fillZeros, recycle)
## nimArrayGeneral(     arg1,      arg2,       arg3,              arg4,            arg5       ,   arg6   ,   arg7 )
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
           ##case: nimNumeric(length = 0, value = 0, init = TRUE)
           nimNumeric = {
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[1], caller=code, callerArgID=3)
               newArgs <- list(type = 'double', nDim = 1, dim = sizeExprs, value = code$args[['value']], init = code$args[['init']],  fillZeros = code$args[['fillZeros']], recycle = code$args[['recycle']])
           },
           ##case: nimInteger(length = 0, value = 0, init = TRUE)
           nimInteger = {
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[1], caller=code, callerArgID=3)
               newArgs <- list(type = 'integer', nDim = 1, dim = sizeExprs, value = code$args[['value']], init = code$args[['init']],  fillZeros = code$args[['fillZeros']], recycle = code$args[['recycle']])
           },
           nimLogical = {
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[1], caller=code, callerArgID=3)
               newArgs <- list(type = 'logical', nDim = 1, dim = sizeExprs, value = code$args[['value']], init = code$args[['init']],  fillZeros = code$args[['fillZeros']], recycle = code$args[['recycle']])
           },
           ##cases: nimVector(type = 'double', length = 0, value = 0, init = TRUE)
           ##       nimVector = {},
           ##       nimMatrix(value = 0, nrow = 1, ncol = 1, init = TRUE, type = 'double')
           nimMatrix = { 
               sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args[c('nrow','ncol')], caller=code, callerArgID=3)
               newArgs <- list(type = code$args[['type']], nDim = 2, dim = sizeExprs, value = code$args[['value']], init = code$args[['init']], fillZeros = code$args[['fillZeros']], recycle = code$args[['recycle']])
           },
           ##case: nimArray(value = 0, dim = c(1, 1), init = TRUE, type = 'double')
           nimArray = {
               ## nimArray will handle dim=c(...), and also dim=EXPR
               if(inherits(code$args[['dim']], 'exprClass') && code$args[['dim']]$isCall && code$args[['dim']]$name == 'nimC') {
                   ## dim argument is c(...)
                   code$args[['dim']]$name <- 'collectSizes'
                   newArgs <- list(type = code$args[['type']], nDim = length(code$args[['dim']]$args), dim = code$args[['dim']], value = code$args[['value']], init = code$args[['init']], fillZeros = code$args[['fillZeros']], recycle = code$args[['recycle']])
               } else {
                   ## dim argument is a single number or expression
                   sizeExprs <- exprClass$new(isName=FALSE, isCall=TRUE, isAssign=FALSE, name='collectSizes', args=code$args['dim'], caller=code, callerArgID=3)
                   if(!inherits(code$args[['dim']], 'exprClass')) {
                       ## a constant was given
                       newArgs <- list(type = code$args[['type']], nDim = 1, dim = sizeExprs, value = code$args[['value']], init = code$args[['init']], fillZeros = code$args[['fillZeros']], recycle = code$args[['recycle']])
                   } else {
                       ## an expression was given: use nDim = -1 as first step to flagging it via unpackNDim for resolution during size processing
                       newArgs <- list(type = code$args[['type']], nDim = -1, dim = sizeExprs, value = code$args[['value']], init = code$args[['init']], fillZeros = code$args[['fillZeros']], recycle = code$args[['recycle']])
                   }
               }
               ## Check if an nDim argument was provided
               if(!is.null(code$args[['nDim']])) { ## Typically nDim would be provided for array if length of dim argument is not known at compile time.  e.g. nimArray(value = v, dim = myDim, type = 'double', nDim = 3).  nDim = 3 is necessary because we often won't know at compile time what the length of myDim is.
                   nDim <- code$args[['nDim']]
                   ## If it was provided in a situation where we haven't flagged it as useful, emit a warning
                   if(newArgs[['nDim']] != -1) if(nDim != newArgs[['nDim']]) warning("Possible nDim mismatch")
                   ## Enter the nDim in the right place and also flag unpackNDim=TRUE to handling during size processing
                   newArgs[['nDim']] <- nDim
                   newArgs[['unpackNDim']] <- TRUE
               }
           },
           stop('There is some problem processing a call to numeric, integer, logical, matrix or array.')
           )
    code$name <- 'nimArrayGeneral'
    code$args <- newArgs
    ## fix AST caller/argument relationships:
    ## set caller / callerArgID for the new nimArraryGeneral() call arguments
    for(i in seq_along(code$args)) {
        if(inherits(code$args[[i]], 'exprClass')) {
            code$args[[i]]$callerArgID <- i
            code$args[[i]]$caller <- code
        }
    }
    ## fix AST caller/argument relationships for the size expressions:
    ## set caller / callerArgID for the c(size expressions...) call
    for(i in seq_along(code$args[['dim']]$args)) {
        if(inherits(code$args[['dim']]$args[[i]], 'exprClass')) {
            code$args[['dim']]$args[[i]]$callerArgID <- i
            code$args[['dim']]$args[[i]]$caller <- code$args[[3]]
        }
    }
    if(!(code$args[['type']] %in% c('double', 'integer', 'logical'))) stop('unknown type in nimArrayGeneral')
    if(code$args[['nDim']] != -1)
        if(is.null(newArgs[['unpackNDim']]))
            if(code$args[['nDim']] != length(code$args[['dim']]$args)) stop('mismatch between nDim and number of size expressions in nimArrayGeneral')
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
    if(code$name %in% c('dwish_chol', 'dinvwish_chol', 'dmnorm_chol', 'dmvt_chol'))
        code$args$overwrite_inputs <- 0
    code$name <- paste0('nimArr_', code$name)
}

rmFunHandler <- function(code, symTab) {
    if(code$name %in% c('rwish_chol', 'rinvwish_chol'))
        code$args$overwrite_inputs <- 0
    dmFunHandler(code, symTab)
    rFunHandler(code, symTab)
}

rFunHandler <- function(code, symTab) {
    ## strip the 1 from the first argument.  In the future we will need to condition on whether this is 1 or >1
    ## This is a funny step because building the R function *inserted* the 1 from the BUGS code
    notOK <- if(!is.numeric(code$args[[1]])) TRUE else code$args[[1]] != 1
    if(notOK) writeLines(paste('Warning: we currently expect to see a 1 and only a 1 as the first argument to an rcat or rinterval call'))
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
