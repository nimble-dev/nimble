assignmentAsFirstArgFuns <- c('nimArr_rmnorm_chol', 'nimArr_rwish_chol', 'nimArr_rmulti', 'nimArr_rdirch', 'getValues')


getAssignmentAsFirstArgFuns <- function() {
    if(!exists('assignmentAsFirstArgFuns'), nimbleUserObjects) {
        return(assignmentAsFirstArgFuns)
    } else {
        return(c(assignmentAsFirstArgFuns, nimbleUserObjects$assignmentAsFirstArgFuns))
    }
}              

sizeCalls <- c(makeCallList(binaryOperators, 'sizeBinaryCwise'),
                makeCallList(binaryMidLogicalOperators, 'sizeBinaryCwiseLogical'),
                makeCallList(binaryOrUnaryOperators, 'sizeBinaryUnaryCwise'),
               makeCallList(unaryOperators, 'sizeUnaryCwise'), 
               makeCallList(unaryOrNonaryOperators, 'sizeUnaryNonaryCwise'),
               makeCallList(assignmentOperators, 'sizeAssign'), 
               makeCallList(reductionUnaryOperators, 'sizeUnaryReduction'), 
               makeCallList(matrixSquareReductionOperators, 'sizeMatrixSquareReduction'),
               makeCallList(reductionBinaryOperators, 'sizeBinaryReduction'),
               makeCallList(matrixMultOperators, 'sizeMatrixMult'), 
               makeCallList(matrixFlipOperators, 'sizeTranspose'),
               makeCallList(matrixSolveOperators, 'sizeSolveOp'), ## TO DO
               makeCallList(matrixEigenOperators, 'sizeEigenOp'), ## TO DO
               makeCallList(matrixSquareOperators, 'sizeUnaryCwiseSquare'), 
               list('return' = 'sizeReturn',
                    'asRow' = 'sizeAsRowOrCol',
                    'asCol' = 'sizeAsRowOrCol',
                    asDoublePtr = 'sizeasDoublePtr',
                   '[' = 'sizeIndexingBracket',
                 ## '[[' for nimbleFunctionList goes through chainedCall
                    chainedCall = 'sizeChainedCall',
                    nfVar = 'sizeNFvar',
                    map = 'sizemap', 
                    ':' = 'sizeColonOperator',
                    dim = 'sizeDimOperator',
                    'if' = 'recurseSetSizes', ##OK
                    'while' = 'recurseSetSizes',
                    callC = 'sizecallC', 
                    'for' = 'sizeFor', 
##                    mvAccess = 'sizeUnaryCwise',
 ##                   numListAccess = 'sizeNumListAccess',                   
#                    numListAccess = 'sizemvAccessBracket',
                    
                    values = 'sizeValues',
                    '(' = 'sizeUnaryCwise',
                    setSize = 'sizeSetSize', ## OK but not done for numericLists
                    resizeNoPtr = 'sizeResizeNoPtr', ## 
                    nimArr_rcat = 'sizeScalarRecurse',
                    nimArr_rinterval = 'sizeScalarRecurse',
                    nimPrint = 'sizeforceEigenize',
                    as.integer = 'sizeUnaryCwise', ## Note as.integer and as.numeric will not work on a non-scalar yet
                    as.numeric = 'sizeUnaryCwise',
                    setAll = 'sizeOneEigenCommand',
                    voidPtr = 'sizeVoidPtr'),
               makeCallList(distributionFuns, 'sizeScalarRecurse'),
               makeCallList(c('isnan','ISNAN','!','ISNA'), 'sizeScalarRecurse'),
               makeCallList(c('nimArr_dmnorm_chol', 'nimArr_dwish_chol', 'nimArr_dmulti', 'nimArr_dcat', 'nimArr_dinterval', 'nimArr_ddirch'), 'sizeScalarRecurse'),
               makeCallList(c('nimArr_rmnorm_chol', 'nimArr_rwish_chol', 'nimArr_rmulti', 'nimArr_rdirch'), 'sizeRmultivarFirstArg'),
               makeCallList(c('calculate', 'getLogProb', 'decide', 'size', 'getsize'), 'sizeScalar'),
               makeCallList(c('simulate', 'blank', 'nfMethod', 'nimFunListAccess', 'getPtr'), 'sizeUndefined')
               )


getSizeCall <- function(codeName) {
    result <- sizeCalls[[codeName]]
    if(is.null(result) && exists('sizeCalls', nimbleUserObjects))
        result <- nimbleUserObjects$sizeCalls[[codeName]]
    return(result)
}


scalarOutputTypes <- list(decide = 'logical', size = 'integer', isnan = 'logical', ISNA = 'logical', '!' = 'logical', nimArr_rcat = 'integer', nimArr_rinterval = 'integer')

## exprClasses_setSizes fills in the type information of exprClass code
## code is an exprClas object
## typeEnv is an environment returned by exprClasses_initSizes
## allowUnknown says whether it is ok to have unknown type. This will be true for LHS of assignments
##
## This returns a set of type assertions collected for each line of code
## This function operates recursively, so the type assertions returned from recursive calls are
## put into the exprClass object for which they were recursed.
##
## For example, if we have A2 <- mean(B + C)
## then typeEnv must have size expressions for B and C to get started.
## If these are matrices, the generic size expressions (for B) will be dim(B)[1] and dim(B)[2]
## Then the exprClass object for `+`(B, C) will generate assertions that dim(B)[1] == dim(C)[1] and dim(B[2]) == dim(C)[2]
## and it will copy the size expressions for B as its own size expressions
## Then the exprClass object for mean(`+`(B, C)) will create a size expression of 1 (with the same dimensions as B+C)
## Then the exprClass object for `<-`(A, mean(`+`(B, C))) will generate assertions that the size of A must be 1
## and it will set the size expressions for A and for itself to 1.
exprClasses_setSizes <- function(code, symTab, typeEnv) { ## input code is exprClass
    ## name:
    if(code$isName) {
        ## If it doesn't exist and must exist, stop
        if(code$name != "") { ## e.g. In A[i,], second index gives name==""
            if(!exists(code$name, envir = typeEnv, inherits = FALSE)) {
                if(symTab$symbolExists(code$name, TRUE)) {
                    code$type <- class(symTab$getSymbolObject(code$name, TRUE))[1]
                } else {
                    code$type <- 'unknown'
                }
            } else {
                ## otherwise fill in type fields from typeEnv object
                info <- get(code$name, envir = typeEnv)
                if(inherits(info, 'exprTypeInfoClass')) {
                    code$type <- info$type
                    code$sizeExprs <- info$sizeExprs
                    code$nDim <- info$nDim
                    code$toEigenize <- 'maybe'

                    ## If it is a vector of known length 1 and it has no indexing, insert a [1] after it
                    if(code$nDim == 1) {
                        if(code$sizeExprs[[1]] == 1) {
                            if(code$caller$name != '[') {
                                insertIndexingBracket(code$caller, code$callerArgID, 1)
                                code$caller$nDim <- 0
                                code$caller$sizeExprs <- list()
                                code$caller$toEigenize <- 'maybe'
                                code$caller$type <- code$type
                            }
                        }
                    }
                }
            }
            ## Note that generation of a symbol for LHS of an assignment is done in the sizeAssign function, which is the handler for assignments
            return(NULL)
        }
    }

    if(code$isCall) {
        if(code$name == '{') {
            ## recurse over lines
            for(i in seq_along(code$args)) {
                if(inherits(code$args[[i]], 'exprClass')) {
                    newAsserts <- exprClasses_setSizes(code$args[[i]], symTab, typeEnv)
                    code$args[[i]]$assertions <- if(is.null(newAsserts)) list() else newAsserts
                }
            }
            return(invisible(NULL))
        }
        
        sizeCall <- getSizeCalls(code$name)

        if(!is.null(sizeCall)) {
            return(eval(call(sizeCall, code, symTab, typeEnv)))
        }

        ## if(!is.null(tempSizeHandlers[[code$name]])) {
        ##     return(tempSizeHandlers[[code$name]](code, symTab, typeEnv))
        ## }

        if(symTab$symbolExists(code$name, TRUE)) { ## could be a nimbleFunction object
            return(sizeNimbleFunction(code, symTab, typeEnv) )
        }
        ## Finally, it could be an RCfunction (a nimbleFunction with no setup == a simple function) {

        if(exists(code$name)) {
            obj <- get(code$name)
            if(is.rcf(obj)) { ## it is an RC function
                nfmObj <- environment(obj)$nfMethodRCobject
                uniqueName <- nfmObj$uniqueName
                if(length(uniqueName)==0) stop("Error, a no-setup nimbleFunction with no internal name is being called:", nimDeparse(code), ". Something wasn't set up right.", call. = FALSE)
                if(is.null(typeEnv$neededRCfuns[[uniqueName]])) {
                    typeEnv$neededRCfuns[[uniqueName]] <- nfmObj
                }
                return(sizeRCfunction(code, symTab, typeEnv, nfmObj))
            }
        }
        
    }
    invisible(NULL)
}

sizeRMVNorm <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(code$args[[1]]$nDim != 1) {
        writeLines(paste0('Uh oh, first arg to rmvnorm in ', nimDeparse(code), ' is not 1D.'))
        browser()
    }
    if(code$args[[2]]$nDim != 1) {
        writeLines(paste0('Uh oh, second arg to rmvnorm in ', nimDeparse(code), ' is not 1D.'))
        browser()
    }
    if(!(code$caller %in% assignmentOperators)) {
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }
    code$type <- 'double'
    code$nDim <- 1
    code$sizeExprs <- code$args[[1]]$sizeExprs ## already checked that it is 1D
    code$toEigenize <- 'no'
}

sizemap <- function(code, symTab, typeEnv) {
    ## This will only be called on a map generated from setup
    ## Maps created from indexing in nimble code don't go through this function
    sym <- symTab$getSymbolObject(code$args[[1]], TRUE)
    code$type <- sym$type
    code$nDim <- code$args[[2]]
    code$sizeExprs <- code$args[[4]]
    code$toEigenize <- 'maybe'
    invisible(NULL)
}

sizeAsRowOrCol <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    if(!inherits(a1, 'exprClass')) stop(paste0('Error in sizeAsRowOrCol for ', nimDeparse(code), '. Argument must be an expression.'), call. = FALSE)
    if(a1$nDim == 0) stop(paste0('Error in sizeAsRowOrCol for ', nimDeparse(code), '. Argument cannot be scalar (could be fixed).'), call. = FALSE)
    code$type <- a1$type
    code$toEigenize <- 'yes'
    if(!code$name %in% c('asRow', 'asCol')) stop('Weird, how did we get to sizeAsRowOrCol without a call to asRow or asCol?', call. = FALSE)

    if(a1$nDim == 1) {
        if(code$name == 'asRow') {
            ##setAsRowOrCol(code, 1, 'asRow', type = a1$type)
            ##code <- removeExprClassLayer(code, 1)
            code$nDim <- 2
            code$sizeExprs <- c(list(1), a1$sizeExprs[[1]])
        } else {
            ##setAsRowOrCol(code, 1, 'asCol', type = a1$type)
            ##code <- removeExprClassLayer(code, 1)
            code$nDim <- 2
            code$sizeExprs <- c(a1$sizeExprs[[1]], list(1))
        }
        return(asserts)
    }

    warning(paste0(' asRow or asCol used on something with more than 1 dimension in ', nimDeparse(code)), call. = FALSE)

}

sizeNFvar <- function(code, symTab, typeEnv) {
    nfName <- code$args[[1]]$name
    nfSym <- symTab$getSymbolObject(nfName, inherits = TRUE)
    if(!inherits(nfSym, 'symbolNimbleFunction')) stop(paste0('Error in ', nimDeparse(code), '. First argument is not a nimbleFunction.'), call. = FALSE)
    nfProc <- nfSym$nfProc
    if(is.null(nfProc)) stop(paste0('Error in ', nimDeparse(code), '. Symbols in this nimbleFunction generation function not set up.'), call. = FALSE)
    objName <- code$args[[2]]
    if(!is.character(objName)) stop(paste0('Error in ', nimDeparse(code), '. Second argument must be a character string.'), call. = FALSE)
    objSym <- nfProc$setupSymTab$getSymbolObject(objName)
    if(is.null(objSym)) stop(paste0('Error in ', nimDeparse(code), '. Symbol not found in the nimbleFunction.'), call. = FALSE)
    code$nDim <- objSym$nDim
    code$type <- objSym$type
    if(code$nDim > 0) {
        code$sizeExprs <- makeSizeExpressions(objSym$size,
                                              parse(text = nimDeparse(code))[[1]])
    } else {
        code$sizeExprs <- list()
    }
    code$toEigenize <- 'maybe'
    NULL
}

sizeChainedCall <- function(code, symTab, typeEnv) { ## at the moment we have only nimFunList[[i]](a), nfMethod(nf, 'foo')(a), or nfMethod(nf[[i]], 'foo')(a)
    ## where actually the [[ would have already been replaced with nimFunListAccess
    ## In other places we generate chainedCalls for static_cast<int>(a), but those shouldn't be seen here
    a1 <- code$args[[1]] 
    if(!inherits(a1, 'exprClass')) stop(paste0('Error in sizeChainedCall for ', nimDeparse(code), '. First arg is not an expression.'), call. = FALSE)
    nfMethodRCobj <- NULL
    if(a1$name == '[[') {
        ## nimFunList[[i]](a)
        recurseSetSizes(a1, symTab, typeEnv, c(FALSE, rep(TRUE, length(a1$args)-1))) ## recursion on this is not done in generalFunSizeHandler because it skips arg1 for chainedCall = TRUE
        if(a1$args[[2]]$nDim != 0) stop(paste0('Error in sizeChainedCall for ', nimDeparse(code), '. Index for nimbleFunction list is not scalar.'), call. = FALSE)
        
        sym <- symTab$getSymbolObject(a1$args[[1]]$name, TRUE)
        if(!inherits(sym, 'symbolNimbleFunctionList')) {
            stop(paste0('Error, in ', nimDeparse(code), '. Expecting a nimbleFunction list.'), call. = FALSE)
        }
        nfMethodRCobj <- getFunctionEnvVar(nf_getGeneratorFunction(sym$baseClass), 'methodList')$run
    }
    else if(a1$name == 'nfMethod') {
        a11 <- a1$args[[1]]
        methodName <- a1$args[[2]]
        if(a11$isName) { ## e.g. in nfMethod(nf, 'foo'), a11 is nf
            sym <- symTab$getSymbolObject(a11$name, TRUE)
         
            nfMethodRCobj <- sym$nfProc$origMethods[[methodName]]
        } else {
            if(a11$name != '[[') stop(paste0('Error, in ', nimDeparse(code), '. Expecting a nimbleFunction list or a nimFun as first arg of nfMethod.'), call. = FALSE)
            ## should look like nfMethod(nflist[[i]], 'foo')
            a111 <- a11$args[[1]]
            sym <- symTab$getSymbolObject(a111$name, TRUE)
            if(!inherits(sym, 'symbolNimbleFunctionList')) {
                stop(paste0('Error, in ', nimDeparse(code), '. Expecting a nimbleFunction list.'), call. = FALSE)
            }
            nfMethodRCobj <- getFunctionEnvVar(nf_getGeneratorFunction(sym$baseClass), 'methodList')[[methodName]]
        }
    }
    else warning(paste0('Warning that we did not know what to do in sizeChainedCall for ', nimDeparse(code)))
    if(!is.null(nfMethodRCobj)) {
        returnType <- nfMethodRCobj$returnType
        argInfo <- nfMethodRCobj$argInfo
        asserts <- generalFunSizeHandler(code, symTab, typeEnv, returnType, argInfo, chainedCall = TRUE)
        return(asserts)
    }
    invisible(NULL)    
    writeLines('Warning')
}

sizeValues <- function(code, symTab, typeEnv) {
    code$nDim <- 1
    code$type <- 'double'
    code$toEigenize <- 'no'
    sym <- symTab$getSymbolObject(code$args[[1]]$name, TRUE)
    if(length(sym$lengthName)==0) stop(paste0("Error the size information for ", nimDeparse(code), " is missing."), call. = FALSE) 
    code$sizeExprs <- list(as.name(sym$lengthName))
    asserts <- list()
   
    if(code$caller$name %in% assignmentOperators) {
        if(code$callerArgID == 2) { ## ans <- values(...)
            code$name <- 'getValues'
            LHS <- code$caller$args[[1]]
            if(LHS$isName) { ## It is a little awkward to insert setSize here, but this is different from other cases in sizeAssignAfterRecursing
                assertSS <- list(substitute(setSize(LHS), list(LHS = as.name(LHS$name))))
                sym <- symTab$getSymbolObject(code$args[[1]]$name, TRUE)
                assertSS[[1]][[3]] <- as.name(sym$lengthName)
                asserts <- c(assertSS, asserts)
            }
        } ## values(...) <- P, don't change it
    } else { ## values(...) embedded in a RHS expression
        code$name <- 'getValues'
        code$toEigenize <- 'yes' ## This tricks sizeAssignAfterRecursing to generate the setSize in asserts
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        code$toEigenize <- 'no'
    }
    
    if(length(asserts)==0) NULL else asserts
}

sizeRCfunction <- function(code, symTab, typeEnv, nfmObj) {
    returnType <- nfmObj$returnType
    argInfo <- nfmObj$argInfo
    code$name <- nfmObj$uniqueName
    asserts <- generalFunSizeHandler(code, symTab, typeEnv, returnType, argInfo)
    return(asserts)
}

sizeNimbleFunction <- function(code, symTab, typeEnv) { ## This will handle other nimbleFunction run (operator()) calls or other methods of this nimbleFunction
    sym <- symTab$getSymbolObject(code$name, TRUE)
    ok <- FALSE
    if(inherits(sym, 'symbolNimbleFunction')) {
        nfMethodRCobj <- environment(sym$nfProc$nfGenerator)$methodList$run
        returnType <- nfMethodRCobj$returnType
        argInfo <- nfMethodRCobj$argInfo
        ok <- TRUE
    }
    if(inherits(sym, 'symbolMemberFunction')) {
        returnType <- sym$nfMethodRCobj$returnType
        argInfo <- sym$nfMethodRCobj$argInfo
        ok <- TRUE
    }
    if(ok) {
        asserts <- generalFunSizeHandler(code, symTab, typeEnv, returnType, argInfo)
        return(asserts)
    }
    stop(paste0('Error, in ', nimDeparse(code), ', the function name is not known and is not a nimbleFunction or a member function.'), call. = FALSE)
}

recurseSetSizes <- function(code, symTab, typeEnv, useArgs = rep(TRUE, length(code$args))) {
    ## won't be here unless code is a call.  It will not be a {
    asserts <- list()
    for(i in seq_along(code$args)) {
        if(useArgs[i]) {
            if(inherits(code$args[[i]], 'exprClass')) {
                asserts <- c(asserts, exprClasses_setSizes(code$args[[i]], symTab, typeEnv))
            }
        }
    }
    if(length(asserts)==0) NULL else asserts
}

## promote numeric output to most information-rich type, double > integer > logical
## Note this will not be correct for logical operators, where output type should be logical
arithmeticOutputType <- function(t1, t2) {
    if(t1 == 'double') return('double')
    if(t2 == 'double') return('double')
    if(t1 == 'integer') return('integer')
    if(t2 == 'integer') return('integer')
    return('logical')
}

## Generate R code for an equality assertion
identityAssert <- function(lhs, rhs, msg = "") {
    if(identical(lhs, rhs)) return(NULL)
    msg <- gsub("\"", "\\\\\"", msg)
    substitute(if(lhs != rhs) nimPrint(msg), list(lhs = lhs, rhs = rhs, msg = msg))
}
## Determine if LHS is less information-rich that RHS and issue a warning.
## e.g. if LHS is int but RHS is double.
assignmentTypeWarn <- function(LHS, RHS) {
    if(LHS == 'int' & RHS == 'double') return(TRUE)
    if(LHS == 'logical' & RHS != 'logical') return(TRUE)
    return(FALSE)
}

## used for setAll
sizeOneEigenCommand <- function(code, symTab, typeEnv) {
    if(!code$args[[1]]$isName) stop(paste0('Error for in sizeOneEigenCommand for ', nimDeparse(code), '. First arg should be a name.'), call. = FALSE)
    recurseSetSizes(code, symTab, typeEnv)
    if(code$args[[1]]$nDim != 2) stop(paste0('Error for in sizeOneEigenCommand for ', nimDeparse(code), '. At the moment only works for 2D objects.'), call. = FALSE)
    code$toEigenize <- 'yes'
    invisible(NULL)
}

## This is used for nimPrint
## If anything has toEigenize == "maybe", the whole expression gets "yes"
## That way cout<<X;  will use an eigen map for X
sizeforceEigenize <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    toEigs <- lapply(code$args, function(x) {
        if(inherits(x, 'exprClass')) x$toEigenize else 'unknown'
    })
    code$toEigenize <- if(any( unlist(toEigs) %in% c('maybe', 'yes'))) 'yes' else 'no'
    code$type <- 'unknown'
    if(length(asserts) == 0) NULL else asserts
}

sizecallC <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code$args[[1]], symTab, typeEnv)
    asserts
}

## This is for when the programmer has directly written "resize(Z, 3, dim(A)[1])".
## When the resize is automatically generated, it skips size inference
sizeSetSize <- function(code, symTab, typeEnv) {
    sym <- symTab$getSymbolObject(code$args[[1]]$name, inherits = TRUE)
    if(!inherits(sym, 'symbolNumericList')) {
        if(sym$nDim == 0) stop(paste0('Error, resizing a scalar does not make sense in: ', nimDeparse(code$args[[1]])), call. = FALSE)
        if(length(code$args) != 1 + sym$nDim) stop(paste0('Error, number of dimensions provided in resize for ', nimDeparse(code$args[[1]]), ' is wrong'), call. = FALSE)
        asserts <- recurseSetSizes(code, symTab, typeEnv, c(FALSE, rep(TRUE, sym$nDim) ) )
        ## May need intermediates if any size provided requires eigenization, although that's hard to imagine
        assign(code$name, exprTypeInfoClass$new(nDim = sym$nDim, sizeExprs = lapply(code$args[-1], nimDeparse), type = sym$type), envir = typeEnv)
        if(length(asserts)==0) NULL else asserts
    }
    if(inherits(sym, 'symbolNumericList') ) {
    	if(length(code$args) != 2 + sym$nDim) stop(paste0('Error, number of dimensions provided in resize for ', nimDeparse(code$args[[1]]), ' is wrong'), call. = FALSE )
    	assign(code$name, exprTypeInfoClass$new(nDim = sym$nDim, sizeExprs = lapply(code$args[-1], nimDeparse), type = sym$type), envir = typeEnv)
    	invisible(NULL)
    }
}

## This was redundant and we should eventually be able to remove it
sizeResizeNoPtr <- function(code, symTab, typeEnv){
    sym <- symTab$getSymbolObject(code$args[[1]]$name, inherits = TRUE)
    if(length(code$args[[2]]) != 1)  stop( paste0('Error, number of dimensions provided in resize for ', nimDeparse(code$args[[1]]), ' is wrong' ), call. = FALSE )
    assign(code$name, exprTypeInfoClass$new(nDim = 1, sizeExprs = lapply(code$args[-1], nimDeparse), type = sym$type), envir = typeEnv)
    invisible(NULL)
}


## Handler for for-loops: a fairly special case
## e.g. for(i in 1:10) {do(i)}
sizeFor <- function(code, symTab, typeEnv) {
    if(length(code$args) != 3) stop('Error in sizeFor: expected 3 arguments to a for-loop', call. = FALSE)
    ## first handle type of the indexing variable
    if(!inherits(code$args[[2]], 'exprClass')) stop('Error in sizeFor: expecting the index range to be an expression (exprClass)', call. = FALSE)
    asserts <- exprClasses_setSizes(code$args[[2]], symTab, typeEnv)
    code$args[[1]]$nDim <- 0
    code$args[[1]]$sizeExprs <- list()
    code$args[[1]]$type <- code$args[[2]]$type
    code$args[[1]]$toEigenize <- 'no'
    ## if index is unknown, create it in typeEnv and in the symTab
    if(!exists(code$args[[1]]$name, envir = typeEnv, inherits = FALSE)) {
        assign(code$args[[1]]$name, exprTypeInfoClass$new(nDim = 0, type = code$args[[1]]$type), envir = typeEnv)
        symTab$addSymbol(symbolBasic(name = code$args[[1]]$name, nDim = 0, type = code$args[[1]]$type))
    }
    typeEnv[[code$args[[1]]$name]]$sizeExprs <- list()

    ## Now the 3rd arg, the body of the loop, can be processed
    asserts <- c(asserts, exprClasses_setSizes(code$args[[3]], symTab, typeEnv))
    ## I think there shouldn't be any asserts returned since the body should be a bracket expression.
    return(invisible(NULL))
}

sizeInsertIntermediate <- function(code, argID, symTab, typeEnv, forceAssign = FALSE) {
    newName <- IntermLabelMaker()

    ## I think it is valid and general to catch maps here.
    ## For most variables, creating an intermediate involves interN <- expression being lifted
    ## But for map, which will be using a NimArr if it is lifted here, what we need to generate is setMap call
    mapcase <- if(is.numeric(code$args[[argID]])) FALSE else (code$args[[argID]]$name == 'map' & !forceAssign) 
    if(mapcase) {
        ans <- nimArrMapExpr(code$args[[argID]], symTab, typeEnv, newName)
        ## That should create the symTab entry
        ans <- RparseTree2ExprClasses(ans)
        newArgExpr <- RparseTree2ExprClasses(as.name(newName))
        newArgExpr$type <- code$args[[argID]]$type
        newArgExpr$sizeExprs <- code$args[[argID]]$sizeExprs
        newArgExpr$toEigenize <- 'maybe'
        newArgExpr$nDim <- code$args[[argID]]$nDim
    } else {

        ## One may wonder where the new variable is added to the
        ## symbolTable.  That happens when we do
        ## sizeAssignAfterRecursing, which identifies unknown LHS and
        ## creates the symTab entry.
        
        newExpr <- newAssignmentExpression()
        setArg(newExpr, 1, RparseTree2ExprClasses(as.name(newName))) 
        setArg(newExpr, 2, code$args[[argID]]) ## The setArg function should set code$caller (to newExpr) and code$callerArgID (to 3)
        ans <- c(sizeAssignAfterRecursing(newExpr, symTab, typeEnv, NoEigenizeMap = TRUE), list(newExpr))

        newArgExpr <- RparseTree2ExprClasses(as.name(newName))
        newArgExpr$type <- newExpr$args[[1]]$type
        newArgExpr$sizeExprs <- newExpr$args[[1]]$sizeExprs
        newArgExpr$toEigenize <- 'maybe'
        newArgExpr$nDim <- newExpr$args[[1]]$nDim
    }
    setArg(code, argID, newArgExpr)
    return(ans) ## This is to be inserted in a list of asserts, even though it is really core code, not just an a test or assertion
}

sizeAssign <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    asserts <- c(asserts, sizeAssignAfterRecursing(code, symTab, typeEnv))
    if(length(asserts) == 0) NULL else asserts
}

## Handler for assignment
sizeAssignAfterRecursing <- function(code, symTab, typeEnv, NoEigenizeMap = FALSE) {
    LHS <- code$args[[1]]
    RHS <- code$args[[2]]
    if(inherits(RHS, 'exprClass')) {
        RHSname <- RHS$name
        RHSnDim <- RHS$nDim
        RHStype <- RHS$type
        RHSsizeExprs <- RHS$sizeExprs
    } else {
        if(is.numeric(RHS)) {
            RHSname = ''
            RHSnDim <- 0
            RHStype <- storage.mode(RHS)
            RHSsizeExprs <- list() 
        } else {
            stop(paste("Error in sizeAssign: don't know what to do with", nimDeparse(RHS), "in", nimDeparse(code)), call. = FALSE) 
        }
    }

    test <- try(if(inherits(RHStype, 'uninitializedField') | length(RHStype)==0) {
        stop('Error with assignment in ', nimDeparse(code), '. Type of RHS is unknown.', call. = FALSE)
    })
    if(inherits(test, 'try-error')) browser()
    
    if(LHS$isName) {
        if(!exists(LHS$name, envir = typeEnv, inherits = FALSE)) { ## not in typeEnv
            ## If LHS unknown, create it in typeEnv
            if(!symTab$symbolExists(LHS$name, TRUE)) { ## not in symTab
                if(RHStype %in% c('double','integer', 'logical')) {  ## valid type to create here
                    assign(LHS$name, exprTypeInfoClass$new(nDim = RHSnDim, type = RHStype), envir = typeEnv)
                    symTab$addSymbol(symbolBasic(name = LHS$name, nDim = RHSnDim, type = RHStype))
                } else { ## not valid type to create here
                    if(RHStype == 'voidPtr') {
                        assign(LHS$name, exprTypeInfoClass$new(nDim = RHSnDim, type = RHStype), envir = typeEnv)
                        symTab$addSymbol(symbolVoidPtr(name = LHS$name, type = RHStype))
                    } 
                    else stop(paste('Error, LHS of ', nimDeparse(code),' is not in typeEnv or symTab but it cannot be added now.'), call. = FALSE)
                }
            } else { ## yes in symTab
                ## This case is ok.  It is in the symbol table but not the typeEnv.  So it is something like ptr <- getPtr(A)
                code$toEigenize <- 'no'
                code$nDim <- 0
                code$type <- 'unknown'
                code$sizeExprs <- list()
                return(invisible(NULL))
            }
        } else { ## yes in typeEnv.  must be symTab too.
            ## If LHS known, check if nDim matches RHS
            if(length(LHS$nDim) == 0) stop(paste0('Error for ', nimDeparse(code), '. nDim for LHS not set.'), call. = FALSE)
            if(length(RHSnDim) == 0) stop(paste0('Error for ', nimDeparse(code), '. nDim for RHS not set.'), call. = FALSE)
            if(LHS$nDim != RHSnDim) {
                warning(paste0('Warning, mismatched dimensions in assignment: ', nimDeparse(code), '. Going to browser(). Press Q to exit'), call. = FALSE )
                browser()
            }
            ## and warn if type issue e.g. int <- double
            if(assignmentTypeWarn(LHS$type, RHStype)) {
                writeLines(paste0('Warning, RHS numeric type is losing information in assignment to LHS.', nimDeparse(code)))
            }
        }
    }
    ## update size info in typeEnv

    assert <- NULL
    
    if(LHS$name == 'values' && length(LHS$args) == 1) { ## It is values(model_values_accessor) <- STUFF
        if(is.numeric(RHS)) stop("Cannot assign into values() from numeric", call. = FALSE)
        code$name <- 'setValues'
        
        code$args <- list(2)
        setArg(code, 1, RHS)
        setArg(code, 2, LHS$args[[1]])
        if(!(RHS$isName)) assert <- c(assert, sizeInsertIntermediate(code, 1, symTab, typeEnv) )
        return( if(length(assert) == 0) NULL else assert )
    }
    
    if(any(unlist(lapply(RHSsizeExprs, is.null)))) RHSsizeExprs <- makeSizeExpressions(rep(NA, RHSnDim), LHS$name) ## reset sizeExprs for the LHS var. re-using RHSsizeExprs for LHS.  This would only be valid if it is a nimbleFunction returning something on the RHS.  For assignment to be executed in Eigen, the RHS sizes MUST be known
    typeEnv[[LHS$name]]$sizeExprs <- RHSsizeExprs

    if(LHS$toEigenize == 'yes') writeLines('Warning from sizeAssign: not expecting LHS to have toEigenize == yes')
    code$toEigenize <-if(inherits(RHS, 'exprClass')) {
        if(RHS$toEigenize == 'no') 'no' else {
            if(RHS$toEigenize == 'unknown') 'no' else {
                if(RHS$toEigenize != 'yes' & (RHS$nDim == 0 | RHS$isName | (RHS$name == 'map' & NoEigenizeMap))) 'no' ## if it is scalar or is just a name or a map, we will do it via NimArr operator= .  Used to have "| RHS$name == 'map'", but this allowed X[1:3] <- X[2:4], which requires eigen, with eval triggered, to get right
                else 'yes' ## if it is 'maybe' and non-scalar and not just a name, default to 'yes'
            }
        }
    } else 'no'
    

    if(code$toEigenize == 'yes') { ## this would make more sense in eigenize_assign
    ## generate setSize(LHS, ...) where ... are dimension expressions
        if(length(RHSnDim) == 0) browser()
        if(RHSnDim > 0) {
            if(TRUE) { ## !identical(LHSdrop$sizeExprs, RHSdrop$sizeExprs)) {## This was too clever: it was to prevent redundant calls to setSize, but the problem is the previous call could have been generated inside an if-then-else, so we can't rely on it
                if(LHS$isName) {
                    assert <- list(substitute(setSize(LHS), list(LHS = as.name(LHS$name))))
                    for(i in seq_along(RHSsizeExprs)) {
                        test <- try(assert[[1]][[i + 2]] <- RHS$sizeExprs[[i]])
                        if(inherits(test, 'try-error')) browser()
                    }
                } else {
##                    cat('Warning: we do not yet generate size matching assertions for indexed LHS expressions')
                }
            }
        }
    } else {
        if(inherits(RHS, 'exprClass')) {
            ## If we have A <- map(B, ...), we need to generate a setMap for the RHS, which will be done by sizeInsertIntermediate 
            if(RHS$name == 'map') assert <- c(assert, sizeInsertIntermediate(code, 2, symTab, typeEnv) )
        }
        if(inherits(LHS, 'exprClass')) {
            # ditto
            if(LHS$name == 'map') assert <- c(assert, sizeInsertIntermediate(code, 1, symTab, typeEnv) )
        }
    }
    
    code$nDim <- code$args[[1]]$nDim <- RHSnDim
    code$type <- code$args[[1]]$type <- RHStype
    code$sizeExprs <- code$args[[1]]$sizeExprs <- RHSsizeExprs

    if(RHSname %in% getAssignmentAsFirstArgFuns()) {
        code$name <- RHS$name
        oldArgs <- RHS$args
        code$args <- list(length(oldArgs) + 1)
        for(i in seq_along(oldArgs)) {
            setArg(code, i+1, oldArgs[[i]])
        }
        setArg(code, 1, LHS)
    }
    return(assert)
}



sizeasDoublePtr <- function(code, symTab, typeEnv) {
    ## This could also handle copies from ints to doubles, which would ALWAYS require a copy
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    nDim <- code$args[[1]]$nDim
    
    targetString <- nimDeparse(code$args[[1]])
    targetName <- Rname2CppName(targetString)
    targetExpr <- parse(text = targetString, keep.source = FALSE)[[1]]
    ptrName <- paste0(targetName, '_DoublePtr')
    copyName <- paste0(targetName, '_contigCopy')
    subList <- list(var = targetExpr, copy = as.name(copyName), ptr = as.name(ptrName))
    codeBefore <- substitute( if(isMap(var) ) { copy <- var; ptr <- getPtr(copy)} else {ptr <- getPtr(var)},
                             subList )
    codeAfter <- substitute( after( if(isMap(var)) { mapCopy(var, copy) } ),  ## after() tags the assertion to go after the code line
                            subList )

    if(!symTab$symbolExists( ptrName ))
        symTab$addSymbol( symbolPtr(name =  ptrName, type = 'double') )
    if(!symTab$symbolExists( copyName )) {
        symTab$addSymbol( symbolBasic(name = copyName, type = 'double', nDim = nDim) )
    }
    
    
    codeBefore <- RparseTree2ExprClasses(codeBefore)
    exprClasses_initSizes(codeBefore, symTab, NULL, typeEnv)
    asserts <- c(asserts, exprClasses_setSizes(codeBefore, symTab, typeEnv))

    codeAfter <- RparseTree2ExprClasses(codeAfter)
    asserts <- c(asserts, exprClasses_setSizes(codeAfter, symTab, typeEnv))
    
    newArgExpr <- RparseTree2ExprClasses( substitute( ptr, subList) )
    setArg(code$caller, code$callerArgID, newArgExpr)
    
    c(asserts, list(codeBefore, codeAfter))
}

sizeScalar <- function(code, symTab, typeEnv) {
    ## use something different for distributionFuns
    if(code$args[[1]]$toEigenize == 'yes') {
        asserts <- sizeInsertIntermediate(code, 1, symTab, typeEnv)
    } else {
        asserts <- NULL
    }
    code$nDim <- 0
    outputType <- scalarOutputTypes[[code$name]]
    if(is.null(outputType)) code$type <- 'double'
    else code$type <- outputType
    code$sizeExprs <- list()
    code$toEigenize <- 'maybe' ## a scalar can be eigenized or not
    invisible(NULL)
}

sizeScalarRecurse <- function(code, symTab, typeEnv) {
    ## use something different for distributionFuns
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    ## This just forces any argument expression to be lifted.  Can we life only things to be eigenized?
    for(i in seq_along(code$args)) {
        if(inherits(code$args[[i]], 'exprClass')) {
            if(!code$args[[i]]$isName) {
                asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv) )
            }
        }
    }
    
    code$nDim <- 0
    outputType <- scalarOutputTypes[[code$name]]
    if(is.null(outputType)) code$type <- 'double'
    else code$type <- outputType
    code$sizeExprs <- list()
    code$toEigenize <- 'maybe' ## a scalar can be eigenized or not
    asserts
}

#sizeNumListAccess <- function(code, symTab, typeEnv){
#	sym = symTab$getSymbolObject(code$args[[1]]$name, TRUE)
#	code$type <- sym$type
#	
#}

sizeUndefined <- function(code, symTab, typeEnv) {
    code$nDim <- 0
    code$type <- as.character(NA)
    code$sizeExprs <- list()
    code$toEigenize <- 'maybe'
    invisible(NULL)
}

sizeBinaryUnaryCwise <- function(code, symTab, typeEnv) {
    if(length(code$args) == 1) return(sizeUnaryCwise(code, symTab, typeEnv))
    if(length(code$args) == 2) return(sizeBinaryCwise(code, symTab, typeEnv))
    
    stop('Error with sizeBinaryUnarycWise: length of arguments is not 1 or 2', call. = FALSE)
}

sizemvAccessBracket <- function(code, symTab, typeEnv) {
    ## this gets called from sizeIndexingBracket, so recurse has already been done
    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, TRUE))
    if(length(code$args) != 2) {
        stop(paste0('Error in sizemvAccessBracket for ', nimDeparse(code),'. Wrong number of indices provided.'), call. = FALSE)
    }
    if(inherits(code$args[[2]], 'exprClass')) {
        if(code$args[[2]]$nDim != 0) stop(paste0('Error in sizemvAccessBracket for ', nimDeparse(code),'. Index is not a scalar.'), call. = FALSE)
    }
    sym <- symTab$getSymbolObject(code$args[[1]]$name, TRUE) ## This is the symbolVecNimArrPtr
    code$type = sym$type
    code$nDim = sym$nDim
    code$sizeExprs <- as.list(sym$size)
    code$toEigenize <- 'maybe'
    code$name <- 'mvAccessRow'
    if(length(asserts)==0) NULL else asserts
}

sizeIndexingBracket <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    if(code$args[[1]]$type == 'symbolVecNimArrPtr') return(c(asserts, sizemvAccessBracket(code, symTab, typeEnv)))
    if(code$args[[1]]$type == 'symbolNumericList') return(c(asserts, sizemvAccessBracket(code, symTab, typeEnv)))
  
    nDimVar <- code$args[[1]]$nDim
    
    if(nDimVar != length(code$args) - 1) {
        stop(paste('Error, wrong number of indices provided for ', nimDeparse(code)), call. = FALSE)
    }
    code$nDim <- nDimVar
    code$type <- code$args[[1]]$type
    code$sizeExprs <- vector('list', length = nDimVar)
    ## We could generate asserts here to ensure sub-indexing is within bounds
    needMap <- FALSE
    iSizes <- 1
    for(i in 1:nDimVar) {
        dropThisDim <- FALSE
        if(is.numeric(code$args[[i+1]])) dropThisDim <- TRUE
        else if((code$args[[i+1]]$name != "") & (length(dropSingleSizes(code$args[[i+1]]$sizeExprs)$sizeExprs) == 0)) dropThisDim <- TRUE
        if(dropThisDim) {
            if(nimbleOptions()$indexDrop) {
                code$sizeExprs[[iSizes]] <- NULL
                code$nDim <- code$nDim - 1
            } else { 
                code$sizeExprs[[iSizes]] <- 1; iSizes <- iSizes + 1
            }
            next
        }
        needMap <- TRUE
        if(inherits(code$args[[i+1]], 'exprClass')) {
            if(code$args[[i+1]]$name != "") {
                ## An entry that is a variable possible with a length
                code$sizeExprs[[iSizes]] <- code$args[[i+1]]$sizeExprs[[1]]
            } else {
                ## blank entry (e.g. A[,i]) is an exprClass with isName = TRUE and name = ""
                code$sizeExprs[[iSizes]] <- code$args[[1]]$sizeExprs[[i]]
            }
            iSizes <- iSizes + 1
            next
        }
    }
    code$toEigenize <- 'maybe'
    if(needMap) {
        ## If this is a map on an *expression* that is not a map, lift it
        ## e.g. (A + B)[1:4] must become (Interm <- A + B; Interm[1:4])
        if(!code$args[[1]]$isName) if(code$args[[1]]$name != 'map') asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
        ## Replace with a map expression if needed
        newExpr <- makeMapExprFromBrackets(code)
        newExpr$sizeExprs <- code$sizeExprs
        newExpr$type <- code$type
        newExpr$nDim <- code$nDim
        newExpr$toEigenize <- code$toEigenize
        setArg(code$caller, code$callerArgID, newExpr)
        
    }
    if(length(asserts)==0) NULL else asserts

}

sizeColonOperator <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(length(code$args) != 2) stop(paste('Error determining size for : without two argument, in expression', nimDeparse(code)), call. = FALSE)
    if(inherits(code$args[[1]], 'exprClass')) {
        if(code$args[[1]]$nDim != 0) writeLines('WARNING, first arg to : is not scalar in expression', nimDeparse(code))
    }
    if(inherits(code$args[[2]], 'exprClass')) {
        if(code$args[[2]]$nDim != 0) writeLines('WARNING, second arg to : is not scalar in expression', nimDeparse(code))
    }

    code$type <- 'integer'
    code$nDim <- 1
    code$toEigenize <- 'maybe' ## doesn't really matter since are used in maps
    
    ## could generate an assertiong that second arg is >= first arg
    if(is.numeric(code$args[[1]]) & is.numeric(code$args[[2]])) {
        code$sizeExprs <- list(code$args[[2]] - code$args[[1]] + 1)
    } else { ## at least one part is an expression
        ## This is an awkward case:
        ## sizeExprs are R parse trees, not exprClasses
        ## But in this case, we want the expression from an exprClass.
        ## so we need to nimDeparse and then parse them
        code$sizeExprs <- list(substitute( A - B + 1, list(A = parse(text = nimDeparse(code$args[[2]]), keep.source = FALSE)[[1]],
                                                           B = parse(text = nimDeparse(code$args[[1]]), keep.source = FALSE)[[1]] ) ) )
    }
    invisible(asserts)
}

sizeDimOperator <- function(code, symTab, typeEnv) {
    ## a special case since the resulte is stored as a vector<int>
    ## recurse and set Intermediate is if it not a name or map
    ## we'll want a NimArr<1, int> constructor vector<int> 
    recurseSetSizes(code, symTab, typeEnv)
    code$type <- 'integer'
    code$nDim <- 1
    code$sizeExprs <- list( code$args[[1]]$nDim )
    code$toEigenize <- 'no'
    invisible(NULL)
}

sizeTranspose <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1) warning(paste0('More than one argument to transpose in ', nimDeparse(code), '.'), call. = FALSE)
    ans <- sizeUnaryCwise(code, symTab, typeEnv)
    if(is.numeric(code$args[[1]])) {
        warning(paste0('Confused by transpose of a numeric scalar in ', nimDeparse(code), '.  Will remove transpose.'), call. = FALSE)
        removeExprClassLayer(code$caller, 1)
        return(ans)
    }
    code$toEigenize <- 'yes'
    code$type <- code$args[[1]]$type
    if(length(code$sizeExprs) == 2) {
        if(code$nDim != 2) warning(paste0('In sizeTranspose, there are 2 sizeExprs but nDim != 2'), call. = FALSE)
        code$sizeExprs <- c(code$sizeExprs[2], code$sizeExprs[1])
    } else if(length(code$sizeExprs) == 1) {
        if(code$nDim != 1) warning(paste0('In sizeTranspose, there is 1 sizeExpr but nDim != 1'), call. = FALSE)
        ##setAsRowOrCol(code, 1, 'asRow', code$args[[1]]$type)
        ##code <- removeExprClassLayer(code, 1)
        ##code$toEigenize <- 'yes'
        code$name <- 'asRow'
        code$sizeExprs <- c(list(1), code$sizeExprs[[1]])
        code$nDim <- 2
    }
    return(ans)
}

## Handler for unary functions that operate component-wise
sizeUnaryCwise <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1){
    	stop('Error, sizeUnaryCwise called with argument length != 1', call. = FALSE)
    }

    asserts <- recurseSetSizes(code, symTab, typeEnv)
    ## lift intermediates
    a1 <- code$args[[1]]
    if(inherits(a1, 'exprClass')) {
        if(a1$toEigenize == 'no') {
            asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            a1 <- code$args[[1]]
        }
        code$nDim <- a1$nDim
        code$sizeExprs <- a1$sizeExprs
        code$type <- a1$type
    } else {
        code$nDim <- 0
        code$sizeExprs <- list()
        code$type <- 'double'
    }
    if(length(code$nDim) != 1) stop(paste0('Error for ', nimDeparse(code), '. nDim is not set.'), call. = FALSE)
    code$toEigenize <- if(code$nDim > 0) 'yes' else 'maybe'
    return(asserts)
}

## currently only inprod(v1, v2)
sizeBinaryReduction <- function(code, symTab, typeEnv) {
    if(length(code$args) != 2) stop('Error, sizeBinaryReduction called with argument length != 2', call. = FALSE);
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]] 
    a2 <- code$args[[2]]

    ok <- TRUE
    if(inherits(a1, 'exprClass')) {
        if(a1$toEigenize == 'no') {
            asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            a1 <- code$args[[1]]
        }
    } else {
        ok <- FALSE
    }
    if(inherits(a2, 'exprClass')) {
        if(a2$toEigenize == 'no') {
            asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
            a2 <- code$args[[2]]
        }
    } else {
        ok <- FALSE
    }
    if(!ok) stop(paste0("Cannot call inprod or other binary reduction operator with constant argument.  Problem is with: ", nimDeparse(code)), call. = FALSE)
    
    code$nDim <- 0
    code$sizeExprs <- list()
    code$type <- 'double'
    code$toEigenize <- 'yes'

    if(!(code$caller$name %in% c('{','<-','<<-','='))) {
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }
    if(length(asserts) == 0) NULL else asserts
}


## things like trace, det, logdet
sizeMatrixSquareReduction <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1){
    	stop('Error, sizeMatrixSquareReduction called with argument length != 1', call. = FALSE)
    }
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    if(!inherits(a1, 'exprClass')) stop('Error, sizeUnaryCwiseSquare called with an argument that is not an expression', call. = FALSE)
    if(a1$toEigenize == 'no') {
        asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
        a1 <- code$args[[1]]
    }
    if(a1$nDim != 2) stop('Error, sizeUnaryCwiseSquare was called with an argument that is not a matrix', call. = FALSE)
    code$nDim <- 0
    code$sizeExprs <- list()
    code$type <- if(code$name == 'trace') code$args[[1]]$type else 'double'
    code$toEigenize <- 'yes'

    if(!(code$caller$name %in% c('{','<-','<<-','='))) {
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }
    if(length(asserts) == 0) NULL else asserts
}

sizeUnaryCwiseSquare <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1){
    	stop('Error, sizeUnaryCwiseSquare called with argument length != 1', call. = FALSE)
    }
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    if(!inherits(a1, 'exprClass')) stop('Error, sizeUnaryCwiseSquare called with an argument that is not an expression', call. = FALSE)
    if(a1$toEigenize == 'no') {
        asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
        a1 <- code$args[[1]]
    }
    if(a1$nDim != 2) stop('Error, sizeUnaryCwiseSquare was called with an argument that is not a matrix', call. = FALSE)
    if(!identical(a1$sizeExprs[[1]], a1$sizeExprs[[2]])) {
        asserts <- identityAssert(a1$sizeExprs[[1]], a1$sizeExprs[[2]], paste0("Run-time size error: expected ", nimDeparse(a1), " to be square.") )
        if(is.integer(a1$sizeExprs[[1]])) {
            newSize <- a1$sizeExprs[[1]]
        } else {
            if(is.integer(a1$sizeExprs[[2]])) {
                newSize <- a1$sizeExprs[[2]]
            } else {
                newSize <- a1$sizeExprs[[1]]
            }
        }
    } else {
        asserts <- NULL
        newSize <- a1$sizeExprs[[1]]
    }
    code$nDim <- 2
    code$sizeExprs <- list(newSize, newSize)
    code$type <- a1$type
    code$toEigenize <- if(code$nDim > 0) 'yes' else 'maybe'
    invisible(asserts)
}

sizeUnaryNonaryCwise <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) stop('Error, sizeUnaryNonarCwise called with argument length > 1', call. = FALSE);
    if(length(code$args) == 1) return(sizeUnaryCwise(code, symTab, typeEnv))
    ## default behavior for a nonary (no-argument) function:
    code$type <- 'double'
    code$nDim <- 0
    code$sizeExprs <- list()
    code$toEigenize <- 'maybe'
    invisible(NULL)
}

## things like min, max, mean, sum
sizeUnaryReduction <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1) stop('Error, sizeUnaryReduction called with argument length != 1', call. = FALSE)
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(inherits(code$args[[1]], 'exprClass')) {
        if(!code$args[[1]]$isName) {
            if(code$args[[1]]$toEigenize == 'no') {
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            }
        }
    }

    code$nDim <- 0
    code$sizeExprs <- list()
    code$type <- code$args[[1]]$type
    code$toEigenize <- 'yes'

    if(!(code$caller$name %in% c('{','<-','<<-','='))) {
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }

    if(length(asserts) == 0) NULL else asserts
}

## There's no real point in annotating return.  Just need to recurse and lift 
sizeReturn <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) stop('Error, return has argument length > 1', call. = FALSE)
    code$toEigenize <- 'no'
    if(length(code$args) == 0) return(invisible(NULL))
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(inherits(code$args[[1]], 'exprClass')) {
        if(!code$args[[1]]$isName) {
##            if(code$args[[1]]$toEigenize == 'yes') {
            if(anyNonScalar(code$args[[1]])) {
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            }
            ##            }
        }
    }
    invisible(asserts)
}

sizeMatrixMult <- function(code, symTab, typeEnv) {
    if(length(code$args) != 2) stop('Error, sizeMatrixMult called with argument length != 2', call. = FALSE);
    a1 <- code$args[[1]]
    a2 <- code$args[[2]]

    if(!(inherits(a1, 'exprClass') & inherits(a2, 'exprClass'))) stop('Error in sizeMatrixMult: expecting both arguments to be expressions.', call. = FALSE)

    ## need promotion from vectors to matrices with asRow or asCol
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    a1 <- code$args[[1]]
    a2 <- code$args[[2]]
    
    if(a1$nDim == 0 | a2$nDim == 0) stop(paste0('Error for ', nimDeparse(code),'. Cannot do matrix multiplication with a scalar.'), call. = FALSE)

    if(a1$toEigenize == 'no') {
        asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
        a1 <- code$args[[1]]
    }
    if(a2$toEigenize == 'no') {
        asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
        a2 <- code$args[[2]]
    }
    
    ## Note that we could insert RUN-TIME adaptation of mat %*% vec and vec %*% mat
    ## but to do so we would need to generate trickier sizeExprs
    ## For now, a vector on the right will be turned into a column
    ## and a vector on the left will be turned into a row
    ## The programmer can always use asRow or asCol to control it explicitly

    if(a1$nDim == 1 & a2$nDim == 1) {
        ##a1 <- setAsRowOrCol(code, 1, 'asRow', type = a1$type)
        origSizeExprs <- a1$sizeExprs[[1]]
        a1 <- insertExprClassLayer(code, 1, 'asRow', type = a1$type)
        a1$sizeExprs <- c(list(1), origSizeExprs)
        ##a2 <- setAsRowOrCol(code, 2, 'asCol', type = a2$type)
        origSizeExprs <- a2$sizeExprs[[1]]
        a2 <- insertExprClassLayer(code, 2, 'asCol', type = a2$type)
        a2$sizeExprs <- c(origSizeExprs, list(1))
    } else {
        if(a1$nDim == 1) {
            if(a2$nDim != 2) stop(paste0('Problem in ', nimDeparse(code), '.  First arg has nDim = 1 and 2nd arg has nDim = ', a2$nDim, '.'), call. = FALSE) 
            origSizeExprs <- a1$sizeExprs[[1]]
            ## For first argument, default to asRow unless second argument has only one row, in which case make first asCol
            if(identical(a2$sizeExprs[[1]], 1)) {
                a1 <- insertExprClassLayer(code, 1, 'asCol', type = a1$type)
                a1$sizeExprs <- c(origSizeExprs, list(1))
            }
            else {
                a1 <- insertExprClassLayer(code, 1, 'asRow', type = a1$type)
                a1$sizeExprs <- c(list(1), origSizeExprs)
            }
        } else if(a2$nDim == 1) {
            origSizeExprs <- a2$sizeExprs[[1]]
            if(a1$nDim != 2) stop(paste0('Problem in ', nimDeparse(code), '.  Second arg has nDim = 1 and 1st arg has nDim = ', a1$nDim, '.'), call. = FALSE) 
            if(identical(a1$sizeExprs[[2]], 1)) {
                ##a2 <- setAsRowOrCol(code, 2, 'asRow', type = a2$type)
                a2 <- insertExprClassLayer(code, 2, 'asRow', type = a2$type)
                a2$sizeExprs <- c(list(1), origSizeExprs)
           }
            else { 
##               a2 <- setAsRowOrCol(code, 2, 'asCol', type = a2$type)
                a2 <- insertExprClassLayer(code, 2, 'asCol', type = a2$type)
                a2$sizeExprs <- c(origSizeExprs, list(1))
            }
        }
    }
    code$nDim <- 2
    code$sizeExprs <- list(a1$sizeExprs[[1]], a2$sizeExprs[[2]])
    code$type <- arithmeticOutputType(a1$type, a2$type)
    code$toEigenize <- 'yes'
    assertMessage <- paste0("Run-time size error: expected ", deparse(a1$sizeExprs[[2]]), " == ", deparse(a2$sizeExprs[[1]]))
    newAssert <- identityAssert(a1$sizeExprs[[2]], a2$sizeExprs[[1]], assertMessage)
    if(is.null(newAssert))
        return(asserts)
    else
        return(c(list(identityAssert(a1$sizeExprs[[2]], a2$sizeExprs[[1]], assertMessage)), asserts))
}

sizeSolveOp <- function(code, symTab, typeEnv) { ## this is for solve(A, b) or forwardsolve(A, b). For inverse, use inverse(A), not solve(A)
    if(length(code$args) != 2) stop('Error, sizeMatrixMult called with argument length != 2', call. = FALSE)
    ## need promotion from vectors to matrices with asRow or asCol
    a1 <- code$args[[1]]
    a2 <- code$args[[2]]
    if(!(inherits(a1, 'exprClass') & inherits(a2, 'exprClass'))) stop('Error in sizeSolveOp: expecting both arguments to be exprClasses', call. = FALSE)
    code$type <- 'double'
    if(a1$nDim != 2) stop('Error, first argument to a matrix solver must be a matrix', call. = FALSE)
    if(a2$nDim == 1) {
        a2 <- insertExprClassLayer(code, 2, 'asCol', type = as$type)
        a2$sizeExprs <- list(code$args[[2]]$sizeExprs[[1]], 1)
        a1$nDim <- 2
    }
    code$nDim <- 2
    code$sizeExprs <- c(a1$sizeExprs[[2]], a2$sizeExprs[[2]])
    assertMessage <- paste0("Run-time size error: expected ", deparse(a1$sizeExprs[[1]]), " == ", deparse(a2$sizeExprs[[1]]))
    return(identityAssert(a1$sizeExprs[[1]], a2$sizeExprs[[1]], assertMessage))
}

## deprecated and will be removed
setAsRowOrCol <- function(code, argID, rowOrCol, type ) {
    recurse <- TRUE
    if(is.numeric(code$args[[argID]])) return(NULL)
    if(code$args[[argID]]$isName) {
        recurse <- FALSE
    } else {
        if(code$args[[argID]]$name == 'map') recurse <- FALSE
    }
    if(!recurse) {
        if(code$args[[argID]]$nDim == 2) {
##            stop('Error in setAsRowOrCol: should only be used one 1-dimensional things', call. = FALSE)
            if(rowOrCol == 'asRow') {
                if(is.numeric(code$args[[argID]]$sizeExprs[[1]])) {
                    if(code$args[[argID]]$sizeExprs[[1]] == 1) {
                        return(invisible(NULL)) ## it is already a row
                    }
                }
                rowOK <- if(is.numeric(code$args[[argID]]$sizeExprs[[2]])) {  ## only ok if a1 2nd size is 1
                    if(code$sizeExprs[[2]] == 1) TRUE else FALSE
                } else FALSE
                if(!rowOK) stop(paste0('Error in setAsRowOrCol for ', nimDeparse(code), '. Cannot convert to row.'), call. = FALSE)
                lengthExpr <- code$args[[argID]]$sizeExprs[[1]]
                insertExprClassLayer(code$caller, code$callerArgID, 'transpose', type = type)
                code$nDim <- 2
                code$sizeExprs <- c(list(1), lengthExpr)
                return(code$args[[argID]])
            } else {
                if(is.numeric(code$args[[argID]]$sizeExprs[[1]])) {
                    if(code$args[[argID]]$sizeExprs[[2]] == 1) {
                        return(invisible(NULL)) ## it is already a col
                    }
                }
                colOK <- if(is.numeric(code$args[[argID]]$sizeExprs[[1]])) {  ## only ok if a1 1st size is 1
                    if(code$sizeExprs[[1]] == 1) TRUE else FALSE
                } else FALSE
                if(!colOK) stop(paste0('Error in setAsRowOrCol for ', nimDeparse(code), '. Cannot convert to col.'), call. = FALSE)
                lengthExpr <- code$args[[argID]]$sizeExprs[[1]]
                insertExprClassLayer(code$caller, code$callerArgID, 'transpose', type = type)
                code$nDim <- 2
                code$sizeExprs <- c(lengthExpr, list(1))
                return(code$args[[argID]])
            }
        } else if(code$args[[argID]]$nDim == 1) {        
            oldSizeExprs <- code$args[[argID]]$sizeExprs
            insertExprClassLayer(code, argID, rowOrCol, type = type)
            if(rowOrCol == 'asRow') {
                code$sizeExprs <- c(list(1), oldSizeExprs)
            } else {
                code$sizeExprs <- c(oldSizeExprs, list(1))
            }
            code$nDim <- 2
            code$type <- type
            ans <- code$args[[argID]]
        }
    } else {
        for(i in seq_along(code$args[[argID]]$args)) {
            setAsRowOrCol(code$args[[argID]], i, rowOrCol, type)
        }
        ans <- code$args[[argID]]
    }
    ans
}

sizeBinaryCwiseLogical <- function(code, symTab, typeEnv) {
    ans <- sizeBinaryCwise(code, symTab, typeEnv)
    code$type <- 'logical'
    return(ans)
}

## Handler for binary component-wise operators
sizeBinaryCwise <- function(code, symTab, typeEnv) {
    if(length(code$args) != 2) stop('Error, sizeBinaryCwise called with argument length != 2', call. = FALSE);

    asserts <- recurseSetSizes(code, symTab, typeEnv)
    
    ## sizes of arguments must have already been set

    ## pull out the two arguments
    a1 <- code$args[[1]] 
    a2 <- code$args[[2]]

    ## pull out aXDropNdim, aXnDim, aXsizeExprs, and aXtype (X = 1 or 2)
    if(inherits(a1, 'exprClass')) {
        if(a1$toEigenize == 'no') {
            asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            a1 <- code$args[[1]]
        }
        a1Drop <- dropSingleSizes(a1$sizeExprs)
        a1DropNdim <- length(a1Drop$sizeExprs)
        a1nDim <- a1$nDim
        a1sizeExprs <- a1$sizeExprs
        a1type <- a1$type
        a1toEigenize <- a1$toEigenize
    } else {
        a1DropNdim <- 0
        a1nDim <- 0
        a1sizeExprs <- list()
        a1type <- 'double'
    }
    if(inherits(a2, 'exprClass')) {
        if(a2$toEigenize == 'no') {
            asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
            a2 <- code$args[[2]]
        }
        a2Drop <- dropSingleSizes(a2$sizeExprs)
        a2DropNdim <- length(a2Drop$sizeExprs)
        a2nDim <- a2$nDim
        a2sizeExprs <- a2$sizeExprs
        a2type <- a2$type
    } else {
        a2DropNdim <- 0
        a2nDim <- 0
        a2sizeExprs <- list()
        a2type <- 'double'
    }
    
    ## Choose the output type by type promotion
    if(length(a1type) == 0) {warning('Problem with type of arg1 in sizeBinaryCwise', call. = FALSE); browser()}
    if(length(a2type) == 0) {warning('Problem with type of arg2 in sizeBinaryCwise', call. = FALSE); browser()}
    code$type <- arithmeticOutputType(a1type, a2type)

    code$toEigenize <- if(a1DropNdim == 0 & a2DropNdim == 0) 'maybe' else 'yes'
    
    ## Catch the case that there is at least one scalar-equivalent (all lengths == 1)
    if(a1DropNdim == 0 | a2DropNdim == 0) { 
        ## Here we will process effective scalar additions
        ## and not do any other type of size promotion/dropping
        if(a1DropNdim == 0) { ## a1 is scalar-equiv
            if(a2DropNdim == 0) { ##both are scalar-equiv
                code$nDim <- max(a1nDim, a2nDim) ## use the larger nDims
                code$sizeExprs <- rep(list(1), code$nDim) ## set sizeExprs to all 1
                code$toEigenize <- 'maybe'
            } else {
                ## a2 is not scalar equiv, so take nDim and sizeExprs from it
                code$nDim <- a2nDim
                code$sizeExprs <- a2sizeExprs
                code$toEigenize <- 'yes'
            }
        } else { ## a2 is scalar-equiv, and a1 is not
            code$nDim <- a1nDim
            code$sizeExprs <- a1sizeExprs
            code$toEigenize <- 'yes'
        }
        return(if(length(asserts) == 0) NULL else asserts)
    }
    
    if(is.null(asserts)) asserts <- list()
    ## Catch the case that the number of dimensions is not equal.
    ## This case doesn't arise as much as it used to because [ (sizeIndexingBracket) now drops single dimensions
    if(a1nDim != a2nDim) {
        ## Catch the case that one is 2D and the other is 1D-equivalent.
        ## This allows e.g. X[1,1:5] + Y[1,1,1:5].  First arg is 2D. 2nd arg is 1D-equivalent. An assertion will check that dim(X)[1] == 1  
        ## If so, wrap the 1D is asRow or asCol to orient it later for Eigen
        if(a1DropNdim == 1 & a2DropNdim == 1) {

            ## Hey, I think this is wrong: I think we should check the aXDropNdims
            if(a1nDim > 2 | a2nDim > 2) stop(paste0('Error for ', nimDeparse(code), '. Dimensions do not match and we do not match Array + Vector for dim(Array) > 2'), call. = FALSE)

            ## a1 is 2D and a2 is 1D
            if(a1nDim == 2 & a2nDim == 1) {
                a1IsCol <- identical(a1sizeExprs[[2]], 1)
                asFun <- if(a1IsCol) 'asCol' else 'asRow'
                ##a2 <- setAsRowOrCol(code, 2, asFun, type = a2type)
                a2 <- insertExprClassLayer(code, 2, asFun, type = a2type)
                a2$sizeExprs <- a1sizeExprs
                a2$nDim <- a1nDim
                a1ind <- if(a1IsCol) 1 else 2
                if(!is.numeric(a1sizeExprs[[a1ind]]) | !is.numeric(a2sizeExprs[[1]])) { ## Really do want original a2sizeExprs
                    assertMessage <- paste0("Run-time size error: expected ", deparse(a1sizeExprs[[a1ind]]), " == ", deparse(a2sizeExprs[[1]]))
                    thisAssert <- identityAssert(a1sizeExprs[[a1ind]], a2sizeExprs[[1]], assertMessage)
                    if(!is.null(thisAssert)) asserts[[length(asserts) + 1]] <- thisAssert
                } else {
                    if(a1sizeExprs[[a1ind]] != a2sizeExprs[[1]]) stop('Fixed size mismatch', call. = FALSE)
                }
                code$nDim <- a1nDim
                code$sizeExprs <- a1sizeExprs
            } else {
                a2IsCol <- identical(a2sizeExprs[[2]], 1)
                asFun <- if(a2IsCol) 'asCol' else 'asRow'
                ##a1 <- setAsRowOrCol(code, 1, asFun, type = a1type)
                a1 <- insertExprClassLayer(code, 1, asFun, type = a1type)
                a1$sizeExprs <- a2sizeExprs
                a1$type <- a1type
                a1$nDim <- a1nDim <- a2nDim
                a2ind <- if(a2IsCol) 1 else 2
                if(!is.numeric(a1sizeExprs[[1]]) | !is.numeric(a2sizeExprs[[a2ind]])) { ## Really do want the original a1sizeExprs[[1]], not the modified one.
                    assertMessage <- paste0("Run-time size error: expected ", deparse(a1sizeExprs[[1]]), " == ", deparse(a2sizeExprs[[a2ind]]))
                    thisAssert <- identityAssert(a1sizeExprs[[1]], a2sizeExprs[[a2ind]], assertMessage)
                    if(!is.null(thisAssert)) asserts[[length(asserts) + 1]] <- thisAssert 
                } else {
                    if(a1sizeExprs[[1]] != a2sizeExprs[[a2ind]]) stop('Fixed size mismatch', call. = FALSE)
                }
                code$nDim <- a2nDim
                code$sizeExprs <- a2sizeExprs
            }
        } else
            stop(paste0('Error for ', nimDeparse(code), '. Dimensions do not match and neither is scalar-equivalent.'), call. = FALSE)
    } else {
        ## dimensions match at the outset
        nDim <- a1nDim
        if(nDim > 0) {
            for(i in 1:nDim) {
                if(!is.numeric(a1sizeExprs[[i]]) | !is.numeric(a2sizeExprs[[i]])) {
                    assertMessage <- paste0("Run-time size error: expected ", deparse(a1sizeExprs[[i]]), " == ", deparse(a2sizeExprs[[i]]))
                    thisAssert <- identityAssert(a1sizeExprs[[i]], a2sizeExprs[[i]], assertMessage)
                    if(!is.null(thisAssert)) asserts[[length(asserts) + 1]] <- thisAssert 
                } else {
                    if(a1sizeExprs[[i]] != a2sizeExprs[[i]]) stop('Fixed size mismatch', call. = FALSE)
                }
            }
        }
        code$nDim <- a1$nDim
        code$sizeExprs <- vector('list', code$nDim)
        for(i in seq_along(code$sizeExprs)) code$sizeExprs[[i]] <- if(is.numeric(a1sizeExprs[[i]])) a1sizeExprs[[i]] else a2sizeExprs[[i]]
    }
    if(length(asserts) == 0) NULL else asserts
}


mvFirstArgCheckLists <- list(nimArr_rmnorm_chol = list(c(1, 2, 0), ## dimensionality of ordered arguments AFTER the first, which is for the return value.  e.g. mean (1D), chol(2D), prec_param(scalar)
                                 1, 'double'), ## 1 = argument from which to take answer size, double = answer type
                             nimArr_rwish_chol = list(c(2, 0, 0), ## chol, df, prec_param
                                 1, 'double'),
                             nimArr_rmulti = list(c(0, 1), ## size, probs
                                 2, 'double'), ## We treat integer rv's as doubles
                             nimArr_rdirch = list(c(1), 1, 'double')) ## alpha

sizeRmultivarFirstArg <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    notOK <- FALSE
    checkList <- mvFirstArgCheckLists[[code$name]]
    if(!is.null(checkList)) {
        if(length(code$args) < length(checkList[[1]])) stop(paste0("Error, not enough arguments provided for ", nimDeparse(code),"."), call. = FALSE)
        for(i in seq_along(checkList[[1]])) {
            notOK <- if(inherits(code$args[[i]], 'exprClass')) code$args[[i]]$nDim != checkList[[1]][i] else notOK            
        }
        returnSizeArgID <- checkList[[2]]
        returnType <- checkList[[3]]
    } else {
        returnSizeArgID <- 1
        returnType <- 'double'
    }

    if(notOK) {
        stop(paste0('Error for ', nimDeparse(code), '. Some argument(s) have the wrong dimension.'), call. = FALSE)
    }

    code$type <- returnType
    code$nDim <- code$args[[returnSizeArgID]]$nDim
    code$toEigenize <- 'maybe'
    code$sizeExprs <- code$args[[returnSizeArgID]]$sizeExprs
    
    for(i in seq_along(code$args)) {
        if(inherits(code$args[[i]], 'exprClass')) {
            if(!code$args[[i]]$isName) {
                asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv) )
            }
        }
    }
    if(code$nDim > 0) {
        if(!(code$caller$name %in% c('{','<-','<<-','='))) {
            asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        }
    }

    return(asserts)
}

sizeVoidPtr <- function(code, symTab, typeEnv) {
	
	
    ## lift any argument that is an expression or scalar.  
    ## We expect only one argument
    ## Lift it if it is an expression, a numeric, or a scalar
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    lift <- TRUE
    if(inherits(code$args[[1]], 'exprClass')) {
    	if(code$args[[1]]$type == 'symbolNimbleFunction') lift <- FALSE
        else if(code$args[[1]]$isName & code$args[[1]]$nDim > 0) lift <- FALSE ## will already be a pointer
        
    }
    if(lift) {
         asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv) )
     }
    code$type <- 'voidPtr'
    code$nDim <- 0
    code$toEigenize <- 'no'
    return(asserts)
}

###
## This function would be called with arguments from an RCfunction or nimbleFunction
## the functions dim and length would be taken over to work on the sizeExprs.
## but for now it can just return NAs for size expressions, and then the new returned value will have default size expressions (dim(name)[1], etc)
##
generalFunSizeHandler <- function(code, symTab, typeEnv, returnType, args, chainedCall = FALSE) {
    useArgs <- unlist(lapply(args, function(x) as.character(x[[1]]) %in% c('double', 'integer', 'logical')))

    if(chainedCall) useArgs <- c(FALSE, useArgs)
    if(length(code$args) != length(useArgs)) {
        stop(paste0('Error, wrong number of args for ',code$name), call. = FALSE)
    }
    ## Note this is NOT checking the dimensions of each arg. useArgs just means it will recurse on that and lift or do as needed

    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs)

    ## lift any argument that is an expression
    for(i in seq_along(code$args)) {
        if(useArgs[i]) {
            if(inherits(code$args[[i]], 'exprClass')) {
                if(!code$args[[i]]$isName) {
                    asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv) )
                }
            }
        }
    }

    returnTypeLabel <- as.character(returnType[[1]])
    if(returnTypeLabel == 'void') {
        code$type <- returnTypeLabel
        code$toEigenize <- 'unknown'
        return(asserts)
    }
    returnNDim <- if(length(returnType) > 1) as.numeric(returnType[[2]]) else 0
    returnSizeExprs <- vector('list', returnNDim) ## This stays blank (NULLs), so if assigned as a RHS, the LHS will get default sizes

    
    code$type <- returnTypeLabel
    code$nDim <- returnNDim
    code$sizeExprs <- returnSizeExprs
    code$toEigenize <- if(code$nDim == 0) 'maybe' else 'no'

    if(code$nDim > 0) {
        if(!(code$caller$name %in% c('{','<-','<<-','='))) {
            asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        }
    }
    
    return(asserts)
}
