assignmentAsFirstArgFuns <- c('nimArr_rmnorm_chol', 'nimArr_rmvt_chol', 'nimArr_rwish_chol', 'nimArr_rmulti', 'nimArr_rdirch', 'getValues', 'getValuesIndexRange', 'initialize', 'setWhich', 'setRepVectorTimes', 'assignVectorToNimArr', 'dimNimArr', 'assignNimArrToNimArr')
setSizeNotNeededOperators <- c('setWhich', 'setRepVectorTimes')
operatorsAllowedBeforeIndexBracketsWithoutLifting <- c('map','dim','mvAccessRow','nfVar')

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
               makeCallList(matrixSolveOperators, 'sizeSolveOp'), 
               makeCallList(matrixSquareOperators, 'sizeUnaryCwiseSquare'),
               makeCallList(matrixEigenListOperators, 'sizeMatrixEigenList'),
               list('debugSizeProcessing' = 'sizeProxyForDebugging',
                    diag = 'sizeDiagonal',
                    dim = 'sizeDim',
                    RRtest_add = 'sizeRecyclingRule',
                    which = 'sizeWhich',
                    nimC = 'sizeConcatenate',
                    nimRep = 'sizeRep',
                    nimSeqBy = 'sizeSeq',
                    nimSeqLen = 'sizeSeq',
                    nimSeqByLen = 'sizeSeq',
                    'return' = 'sizeReturn',
                    'asRow' = 'sizeAsRowOrCol',
                    'asCol' = 'sizeAsRowOrCol',
                    makeNewNimbleListObject = 'sizeNewNimbleList',
                    getParam = 'sizeGetParam',
                    getBound = 'sizeGetBound',
                    nimSwitch = 'sizeSwitch',
                    asDoublePtr = 'sizeasDoublePtr',
                   '[' = 'sizeIndexingBracket',
                   '[[' = 'sizeDoubleBracket', ## for nimbleFunctionList, this will always  go through chainedCall(nfList[[i]], 'foo')(arg1, arg2)
                    chainedCall = 'sizeChainedCall',
                    nfVar = 'sizeNFvar',
                    map = 'sizemap', 
                    ':' = 'sizeColonOperator',
                    ##dim = 'sizeDimOperator',
                    'if' = 'recurseSetSizes', ##OK
                    'while' = 'recurseSetSizes',
                    callC = 'sizecallC', 
                    'for' = 'sizeFor', 
                    
                    values = 'sizeValues',
                    '(' = 'sizeUnaryCwise',
                    setSize = 'sizeSetSize', ## OK but not done for numericLists
                    resizeNoPtr = 'sizeResizeNoPtr', ## may not be used any more 
                    nimArr_rcat = 'sizeScalarRecurse',
                    nimArr_rinterval = 'sizeScalarRecurse',
                    nimPrint = 'sizeforceEigenize',
                    ##nimCat = 'sizeforceEigenize',
                    as.integer = 'sizeUnaryCwise', ## Note as.integer and as.numeric will not work on a non-scalar yet
                    as.numeric = 'sizeUnaryCwise',
                    nimArrayGeneral = 'sizeNimArrayGeneral',
                    setAll = 'sizeOneEigenCommand',
                    voidPtr = 'sizeVoidPtr',
                    run.time = 'sizeRunTime'),
               makeCallList(scalar_distribution_dFuns, 'sizeRecyclingRule'),
               makeCallList(scalar_distribution_pFuns, 'sizeRecyclingRule'),
               makeCallList(scalar_distribution_qFuns, 'sizeRecyclingRule'),
               makeCallList(scalar_distribution_rFuns, 'sizeRecyclingRuleRfunction'),
               makeCallList(distributionFuns[!(distributionFuns %in% c(scalar_distribution_dFuns, scalar_distribution_pFuns, scalar_distribution_qFuns, scalar_distribution_rFuns))], 'sizeScalarRecurse'),
               # R dist functions that are not used by NIMBLE but we allow in DSL
               makeCallList(paste0(c('d','q','p'), 't'), 'sizeRecyclingRule'),
               rt = 'sizeRecyclingRuleRfunction',
               makeCallList(paste0(c('d','q','p'), 'exp'), 'sizeRecyclingRule'),
               rexp = 'sizeRecyclingRuleRfunction',
               makeCallList(c('isnan','ISNAN','ISNA'), 'sizeScalarRecurse'),
               makeCallList(c('nimArr_dmnorm_chol', 'nimArr_dmvt_chol', 'nimArr_dwish_chol', 'nimArr_dmulti', 'nimArr_dcat', 'nimArr_dinterval', 'nimArr_ddirch'), 'sizeScalarRecurse'),
               makeCallList(c('nimArr_rmnorm_chol', 'nimArr_rmvt_chol', 'nimArr_rwish_chol', 'nimArr_rmulti', 'nimArr_rdirch'), 'sizeRmultivarFirstArg'),
               makeCallList(c('decide', 'size', 'getsize','getNodeFunctionIndexedInfo', 'endNimbleTimer'), 'sizeScalar'),
               makeCallList(c('calculate','calculateDiff', 'getLogProb'), 'sizeScalarModelOp'),
               simulate = 'sizeSimulate',
               makeCallList(c('blank', 'nfMethod', 'getPtr', 'startNimbleTimer'), 'sizeUndefined')
               )

scalarOutputTypes <- list(decide = 'logical', size = 'integer', isnan = 'logical', ISNA = 'logical', '!' = 'logical', getNodeFunctionIndexedInfo = 'double', endNimbleTimer = 'double') # , nimArr_rcat = 'double', nimArr_rinterval = 'double')

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
                    ##if(exists('.AllowUnknowns', envir = typeEnv)) 
                        if(!typeEnv$.AllowUnknowns)
                            warning(paste0("variable '",code$name,"' has not been created yet."), call.=FALSE) 
                }
            } else {
                ## otherwise fill in type fields from typeEnv object
                info <- get(code$name, envir = typeEnv)
                if(inherits(info, 'exprTypeInfoClass')) {
                    code$type <- info$type
                    code$sizeExprs <- info$sizeExprs
                    code$nDim <- info$nDim
                    code$toEigenize <- 'maybe'
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
        sizeCall <- sizeCalls[[code$name]]
        if(!is.null(sizeCall)) {
            if(.nimbleOptions$debugSizeProcessing) {
                browser()
                eval(substitute(debugonce(XYZ), list(XYZ = as.name(sizeCall))))
            }
          test0 <- eval(call(sizeCall, code, symTab, typeEnv))
            # if(sizeCall == 'sizeAssign'){
            #   print(test0)
            #   if(length(test0)>0) browser()
            # }
            return(test0)
        }
        if(symTab$symbolExists(code$name, TRUE)) { ## could be a nimbleFunction object
            return(sizeNimbleFunction(code, symTab, typeEnv) )
        }
        ## Finally, it could be an RCfunction (a nimbleFunction with no setup == a simple function) {

        if(exists(code$name)) {
            obj <- get(code$name)
            if(is.rcf(obj)) { ## it is an RC function
                nfmObj <- environment(obj)$nfMethodRCobject
                uniqueName <- nfmObj$uniqueName
                if(length(uniqueName)==0) stop(exprClassProcessingErrorMsg(code, 'In size processing: A no-setup nimbleFunction with no internal name is being called.'), call. = FALSE)
                if(is.null(typeEnv$neededRCfuns[[uniqueName]])) {
                    typeEnv$neededRCfuns[[uniqueName]] <- nfmObj
                }
                return(sizeRCfunction(code, symTab, typeEnv, nfmObj))
            }
        }
    }
    invisible(NULL)
}

sizeProxyForDebugging <- function(code, symTab, typeEnv) {
    browser()
    origValue <- .nimbleOptions$debugSizeProcessing
    message('Entering into size processing debugging. You may need to do nimbleOptions(debugSizeProcessing = FALSE) if this exits in any non-standard way.')
    setNimbleOption('debugSizeProcessing', TRUE)
    ans <- recurseSetSizes(code, symTab, typeEnv)
    removeExprClassLayer(code$caller, 1)
    setNimbleOption('debugSizeProcessing', origValue)
    return(ans)
}

productSizeExprs <- function(sizeExprs) {
    if(length(sizeExprs)==0) return(1)
    if(length(sizeExprs)==1) return(sizeExprs[[1]])
    ans <- substitute( (A), list(A = sizeExprs[[1]]))
    for(i in 2:length(sizeExprs)) {
        ans <- substitute(A * (B), list(A = ans, B = sizeExprs[[i]]))
    }
    ans
}

multiMaxSizeExprs <- function(code, useArgs = rep(TRUE, length(code$args))) {
    if(length(code$args)==0) return(list()) ## probably something wrong
    codeArgsUsed <- code$args[useArgs]
    totalLengthExprs <- lapply(codeArgsUsed, function(x) if(inherits(x, 'exprClass')) productSizeExprs(x$sizeExprs) else 1)
    if(length(codeArgsUsed)==1) return(totalLengthExprs) ## a list of length 1
    numericTotalLengths <- unlist(lapply(totalLengthExprs, is.numeric))
    if(sum(numericTotalLengths) > 0) {
        maxKnownSize <- max(unlist(totalLengthExprs[numericTotalLengths]))
        if(sum(numericTotalLengths)==length(totalLengthExprs)) return(list(maxKnownSize))
        totalLengthExprs <- c(list(maxKnownSize), totalLengthExprs[-which(numericTotalLengths)])
    }
    numArgs <- length(totalLengthExprs) ## must be > 1 or it would have returned two lines above
    if(numArgs == 1) return(totalLengthExprs[[1]]) ## but check anyway
    lastMax <- substitute(max(A, B), list(A = totalLengthExprs[[numArgs]], B = totalLengthExprs[[numArgs-1]]))
    if(numArgs > 2) {
        for(i in (numArgs-2):1) {
            lastMax <- substitute(max(A, B), list(A = totalLengthExprs[[i]], B = lastMax))
        }
    }
    return(list(lastMax))
}

addDIB <- function(name, type) {
    paste0(name, switch(type, double = 'D', integer = 'I', logical = 'B'))
}

sizeDim <- function(code, symTab, typeEnv) {
 ##   if(code$caller$name != '[') return(list()) ## This gets specially handled in sizeIndexingBracket
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(!inherits(code$args[[1]], 'exprClass')) {
        stop(exprClassProcessingErrorMsg(code, paste0('Argument of dim is not valid')), call. = FALSE)
    }
    if(code$args[[1]]$nDim == 0) {
        stop(exprClassProcessingErrorMsg(code, paste0('dim() cannot take a scalar as its argument.')), call. = FALSE)
    }
    if(code$caller$name != '[') code$name <- 'dimNimArr'
    code$nDim <- 1
    code$type <- 'integer'
    code$toEigenize <- 'no'
    code$sizeExprs <- list( code$args[[1]]$nDim )
    return(if(is.null(asserts)) list() else asserts)
}

sizeDiagonal <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    argIsExprClass <- inherits(code$args[[1]], 'exprClass')
    nDimArg <- if(argIsExprClass) code$args[[1]]$nDim else 0
    if(nDimArg == 0) {
        code$nDim <- 2
        code$type <- 'double'
        newSizeExpr <- parse(text = nimDeparse(code$args[[1]]), keep.source = FALSE)[[1]]
        code$sizeExprs <- list(newSizeExpr, newSizeExpr)
        code$toEigenize <- 'yes'
        code$name <- addDIB('nimDiagonal', code$type) ## These all go to double anyway
        return( if(length(asserts) == 0) NULL else asserts )
    }
    if(nDimArg == 1) {
        code$nDim <- 2
        code$type <- 'double'
        newSizeExpr <- code$args[[1]]$sizeExprs[[1]]
        code$sizeExprs <- list(newSizeExpr, newSizeExpr)
        code$toEigenize <- 'yes'
        code$name <- addDIB('nimDiagonal', code$args[[1]]$type) ## double anyway
        return( if(length(asserts) == 0) NULL else asserts )
    }

    if(nDimArg == 2) {
        code$nDim <- 1
        code$type <- code$args[[1]]$type
        code$sizeExprs <- list(substitute(min(X, Y), list(X = code$args[[1]]$sizeExprs[[1]], Y = code$args[[1]]$sizeExprs[[2]])))
        code$toEigenize <- 'yes'
        code$name <- 'diagonal'
        return( if(length(asserts) == 0) NULL else asserts )
    } 
    stop(exprClassProcessingErrorMsg(code, paste0('Something is wrong with this usage of diag()')), call. = FALSE)
}

sizeWhich <- function(code, symTab, typeEnv) {
    ## which is a somewhat unique construction.
    ## It should only appear as
    ## answer <- which(boolExpr)
    ## and should be lifted to an intermediate if necessary
    ## The sizeExprs on "which" in the syntax tree will be NULL
    ## which will trigger sizeAssignAfterRecursing to make default size expressions on "answer"
    ## and then it will transform to
    ## setWhich(answer, boolExpr) for C++ output
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    code$type = 'integer'
    code$sizeExprs <- list(NULL)
    code$nDim <- 1
    code$toEigenize <- 'yes'
    code$name <- 'setWhich'

    if(!(code$caller$name %in% assignmentOperators)) {
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }
    if(length(asserts) == 0) NULL else asserts
}

sizeRecyclingRule <- function(code, symTab, typeEnv) { ## also need an entry in eigenization.
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    ## for now this is entirely for d, p and q distribution functions, so we'll look up number of arguments for recycling
    numArgs <- length(code$args)
    if(numArgs == 0) return(asserts)

    recycleArgs <- rep(TRUE, numArgs)

    dFunName <- code$name
    substr(dFunName, 1, 1) <- 'd' 
    
    thisDist <- distributions$distObjects[[dFunName]]
    if(!is.null(thisDist)) {
        numReqdArgs <- length(thisDist$reqdArgs)
        recycleArgs[-(1:(numReqdArgs+1))] <- FALSE
    }
    newSizeExprs <- multiMaxSizeExprs(code, recycleArgs)
    if(length(newSizeExprs)==1)
        if(is.numeric(newSizeExprs[[1]]))
            if(newSizeExprs[[1]] == 1)
                return(c(asserts, sizeScalarRecurse(code, symTab, typeEnv, recurse = FALSE))) ## ALSO NEED ALL ARGS TO HAVE nDim 0
    code$sizeExprs <- newSizeExprs
    code$type <- 'double' ## will need to look up from a list
    code$nDim <- 1
    code$toEigenize <- TRUE
    return(asserts)
}

sizeRecyclingRuleRfunction <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    ## for now this is entirely for r distribution functions, so we'll look up number of arguments for recycling
    numArgs <- length(code$args)
    if(numArgs == 0) return(asserts)

    ## Size determined by first arg
    ## If scalar, that gives size
    ## If vector, size is length of first argument.
    ## Problem is vector of length 1, where size should be value of first element, not length of 1.
    if(inherits(code$args[[1]], 'exprClass')) {
        if(!code$args[[1]]$isName) {
            asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
        }
        newSizeExprs <- list(substitute(rFunLength(X), list(X = as.name(code$args[[1]]$name))))
    } else {
        newSizeExprs <- list(code$args[[1]])
    }
    
    if(length(newSizeExprs)==1)
        if(is.numeric(newSizeExprs[[1]]))
            if(newSizeExprs[[1]] == 1) {
                code$args[[1]] <- NULL ## strip the first argument, which should be a 1 if we are here
                for(i in seq_along(code$args)) {
                    if(inherits(code$args[[i]], 'exprClass')) {
                        code$args[[i]]$callerArgID <- code$args[[i]]$callerArgID - 1
                    }
                }
                return(c(asserts, sizeScalarRecurse(code, symTab, typeEnv, recurse = FALSE)))
            }
    code$sizeExprs <- newSizeExprs
    code$type <- 'double' ## will need to look up from a list
    code$nDim <- 1
    code$toEigenize <- TRUE
    return(asserts)
}


sizeConcatenate <- function(code, symTab, typeEnv) { ## This is two argument version
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    ## must recurse to get nDims set
    isScalar <- unlist(lapply(code$args, function(x) if(inherits(x, 'exprClass')) x$nDim == 0 else TRUE))
    ##isExprClass <- unlist(lapply(code$args, function(x) inherits(x, 'exprClass')))
    argRLE <- rle(isScalar)
    newNumArgs <- sum(argRLE$values) + sum(argRLE$lengths[!argRLE$values]) ## number of scalar runs + sum of non-scalar runs * run-lengths
    newArgs <- vector(length(newNumArgs), mode = 'list')
    iInput <- 1
    iOutput <- 1
    concatenateIntermLabelMaker <- labelFunctionCreator("ConcatenateInterm")
    asserts <- NULL
    for(i in seq_along(argRLE$values)) {
        thisLength <- argRLE$lengths[i]
        if(!(argRLE$values[i])) {
            newArgs[(iOutput-1) + (1:thisLength)] <- code$args[(iInput-1) + (1:thisLength)]
            iInput <- iInput + thisLength
            iOutput <- iOutput + thisLength
        } else {
            newTempFixedName <- concatenateIntermLabelMaker()
            newTempVecName <- concatenateIntermLabelMaker()
            newExpr <- exprClass(isName = FALSE, isCall = TRUE, isAssign = FALSE, name = "concatenateTemp", nDim = 1, sizeExprs = list(thisLength), type = 'double')
            setArg(newExpr, 1, exprClass(isName = TRUE, isCall = FALSE, isAssign = FALSE, name = newTempVecName, nDim = 1, sizeExprs = list(thisLength), type = 'double'))
            valuesExpr <- quote(hardCodedVectorInitializer())
            thisType <- 'logical'
            for(j in 1:thisLength) {
                valuesExpr[[j+1]] <- parse(text = nimDeparse(code$args[[iInput - 1 + j]]), keep.source = FALSE)[[1]]
                if(inherits(code$args[[iInput - 1 + j]], 'exprClass'))
                    thisType <- arithmeticOutputType(thisType, code$args[[iInput - 1 + j]]$type)
                else
                    thisType <- 'double'
            }
            newExpr$type <- thisType
            newExpr$args[[1]]$type <- thisType
            iInput <- iInput + thisLength
            if(thisType == 'integer') thisType <- 'int'
            if(thisType == 'logical') thisType <- 'bool'
            newAssert <- substitute(MAKE_FIXED_VECTOR(newTempVecName, newTempFixedName, thisLength, valuesExpr, thisType),
                                    list(newTempVecName = newTempVecName, newTempFixedName = newTempFixedName,
                                         thisLength = as.numeric(thisLength), valuesExpr = valuesExpr, thisType = thisType))
            newAssert <- as.call(newAssert)
            asserts <- c(asserts, list(newAssert))
            newArgs[[iOutput]] <- newExpr
            iOutput <- iOutput + 1
        }
    }

    maxArgsOneCall <- 4
    numArgGroups <- ceiling(newNumArgs / (maxArgsOneCall-1))
    splitArgIDs <- split(1:newNumArgs, rep(1:numArgGroups, each = maxArgsOneCall-1, length.out = newNumArgs))
    ## if last is a singleton it can be put with previous group
    if(length(splitArgIDs[[numArgGroups]]) == 1) {
        if(numArgGroups > 1) {
            splitArgIDs[[numArgGroups-1]] <- c(splitArgIDs[[numArgGroups-1]], splitArgIDs[[numArgGroups]])
            splitArgIDs[[numArgGroups]] <- NULL
            numArgGroups <- numArgGroups-1
        }
    }

    newExprList <- vector(numArgGroups, mode = 'list')
    for(i in seq_along(splitArgIDs)) {
        newExprList[[i]] <- exprClass(isName = FALSE, isCall = TRUE, isAssign = FALSE, name = 'nimC', nDim = 1, toEigenize = 'yes', type = 'double')
        for(j in seq_along(splitArgIDs[[i]])) setArg(newExprList[[i]], j, newArgs[[splitArgIDs[[i]][j]]])
    }
    ## Last step is to set up nesting and make sizeExprs
    for(i in seq_along(splitArgIDs)) {
        if(i != length(splitArgIDs)) {
            setArg(newExprList[[i]], maxArgsOneCall, newExprList[[i+1]])
        }
    }
    for(i in rev(seq_along(splitArgIDs))) {
        if(inherits(newExprList[[i]]$args[[1]], 'exprClass')) {
            thisSizeExpr <-productSizeExprs(newExprList[[i]]$args[[1]]$sizeExprs) 
            thisType <- newExprList[[i]]$args[[1]]$type
        } else {
            thisSizeExpr <- 1
            thisType <- 'double'
        }
        for(j in seq_along(newExprList[[i]]$args)) {
            if(j == 1) next
            if(inherits(newExprList[[i]]$args[[j]], 'exprClass')) {
                thisSizeExpr <- substitute( (A) + (B),
                                           list(A = thisSizeExpr,
                                                B = productSizeExprs(newExprList[[i]]$args[[j]]$sizeExprs)
                                                ))
                thisType <- arithmeticOutputType(thisType, newExprList[[i]]$args[[j]]$type)
            } else {
                thisSizeExpr <- substitute( (A) + 1,
                                           list(A = thisSizeExpr
                                                ))
                thisType <- 'double'
            }
        }
        if(thisType == 'double') newExprList[[i]]$name <- 'nimCd' ## this change could get moved to genCpp_generateCpp 
        if(thisType == 'integer') newExprList[[i]]$name <- 'nimCi'
        if(thisType == 'logical') newExprList[[i]]$name <- 'nimCb'
        newExprList[[i]]$type <- thisType
        newExprList[[i]]$sizeExprs <- list(thisSizeExpr)
    }
    setArg(code$caller, code$callerArgID, newExprList[[1]])
    return(asserts)
    
    ## code$type <- arithmeticOutputType(code$args[[1]]$type, code$args[[2]]$type)
    ## if(code$type == 'double') code$name <- 'nimCd' ## this change could get moved to genCpp_generateCpp 
    ## if(code$type == 'integer') code$name <- 'nimCi'
    ## if(code$type == 'logical') code$name <- 'nimCb'
    
    ## if(code$args[[1]]$nDim > 2 | code$args[[2]]$nDim > 2) stop(exprClassProcessingErrorMsg(code, paste0('Arguments to c() must have dimension <= 2 for now.')), call. = FALSE)
    ## ## need to deal with lifting size expressions to get them fully compiled.
    ## ## otherwise we could end up with an assertion output somewhat incorrectly
    ## ## insertAssertions does convert to exprClasses but does not eigenize etc
    ## thisSizeExpr <- substitute( AAA_ + BBB_,
    ##                            list(AAA_ = productSizeExprs(code$args[[1]]$sizeExprs),   ## need to create products of sizeExprs in general
    ##                                 BBB_ = productSizeExprs(code$args[[2]]$sizeExprs)) )
    ## code$sizeExprs <- list(thisSizeExpr)
    ## code$nDim <- 1
    ## code$toEigenize <- 'yes'
    ## return(asserts)
}

sizeRep <- function(code, symTab, typeEnv) {
    ## if times is a vector: If length.out is provided, times is always ignored
    ## otherwise lift and use assignEigenToNIMBLE
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    xIsExpr <- inherits(code$args[[1]], 'exprClass')
    code$type <- if(xIsExpr) code$args[[1]]$type else 'double'

    ##    if((!inherits(code$args[[1]], 'exprClass')) || code$args[[1]]$nDim != 1) stop(exprClassProcessingErrorMsg(code, paste0('First argument to rep() must be a vector for now.')), call. = FALSE)
    includesLengthOut <- length(code$args) > 3
    if(inherits(code$args[[2]], 'exprClass')) if(code$args[[2]]$nDim != 0 && !includesLengthOut) { ## times is a vector and length.out not provided
        if(!(code$caller$name %in% assignmentOperators)) {
            asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        }
        if(length(code$args) > 2) code$args[[3]] <- NULL
        code$name <- 'setRepVectorTimes'
        code$sizeExprs <- list(NULL)
        code$nDim <- 1
        code$toEigenize <- 'yes'
        return(asserts)
    }
    
    if(code$type == 'double') code$name <- 'nimRepd' ## this change could get moved to genCpp_generateCpp 
    if(code$type == 'integer') code$name <- 'nimRepi'
    if(code$type == 'logical') code$name <- 'nimRepb'

    ## requiring for now that times and each arguments are given as integers, not expressions
    ## Since these will go into sizeExprs, which are then processed as R expressions, then as exprClasses but not fully size processed,
    ## any expression should be lifted
    if(includesLengthOut) { ## there is a "length.out" argument        ## need to lift length.out if it is more than a name or constant
        if(inherits(code$args[[2]], 'exprClass')) if(code$args[[2]]$nDim > 0) stop(exprClassProcessingErrorMsg(code, paste0('times argument to rep() must be scalar is length.out is also provided.')), call. = FALSE)
        if(inherits(code$args[[3]], 'exprClass')) { ## if length.out is present, it is argument 3
            if(!is.name(code$args[[3]])) 
                asserts <- c(asserts, sizeInsertIntermediate(code, 3, symTab, typeEnv))
            if(code$args[[3]]$nDim > 0)
                code$sizeExprs <- list( parse(text = paste0(nimDeparse(code$args[[3]]),'[1]'), keep.source = FALSE)[[1]])
            else
                code$sizeExprs <- list( parse(text = nimDeparse(code$args[[3]]), keep.source = FALSE)[[1]])
        } else {
            code$sizeExprs <- list(code$args[[3]])
        }
    } else { ## length.out absent, so times is second and each is third
        for(i in 2:3) {
            if(inherits(code$args[[i]], 'exprClass'))
                if(!is.name(code$args[[i]])) 
                    asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv))
        }
        part2 <- nimDeparse(code$args[[2]])
        if(inherits(code$args[[2]], 'exprClass')) if(code$args[[2]]$nDim > 0) part2 <- paste0(part2, '[1]')
        part3 <- nimDeparse(code$args[[3]])
        if(inherits(code$args[[3]], 'exprClass')) if(code$args[[3]]$nDim > 0) part3 <- paste0(part3, '[1]')

        thisSizeExpr <- substitute( (AAA_) * (BBB_) * (CCC_),
                                   list(AAA_ = if(xIsExpr) productSizeExprs(code$args[[1]]$sizeExprs) else 1,
                                        BBB_ = parse(text = part2, keep.source = FALSE)[[1]],
                                        CCC_ = parse(text = part3, keep.source = FALSE)[[1]] ))
        code$sizeExprs <- list(thisSizeExpr)
    }
    code$nDim <- 1
    code$toEigenize <- 'yes'
    return(asserts)
}

sizeNewNimbleList <- function(code, symTab, typeEnv){
  ## code looks like: nimListDef$new(a = '', b = 12)
  ## want to change code$caller to :
  ## { nimList <- nimListDef$new()
  ## nimList$a <- 10
  ## nimList$b <- 12 }
  ## accomplish this by copying code, getting arguments (e.g. a = 10, b = 12) from copied code and turning them into assignment 
  ## exprs in code$caller, and setting first argument of code$caller to be nimList <- nimListDef$new()
  listDefName <- code$args[[1]]$name
  if(symTab$parentST$symbolExists(listDefName)){
    listST <- symTab$getSymbolObject(listDefName, inherits = TRUE)
    code$type <- "symbolNimbleList"
    code$sizeExprs <- listST
    code$toEigenize <- "maybe"
    code$nDim <- 0
  }
  else stop('Error in sizeNewNimbleList: listGenerator not found in parentST', call. = FALSE)
  
  asserts <- list()
  if(!(code$caller$name %in% assignmentOperators)){
    asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
  } else {
      asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, useArgs = c(TRUE, rep(FALSE, length(code$args)-1))), asserts)
  }
  if(length(code$args)>1){   
    RnewExprs <- list()
    newExprs <- list()
    RnfVarExprs <- list()
    nfVarExprs <- list()
    exprCounter <- 1
    originalCode <- code 
    listElements <- listST$nlProc$symTab$getSymbolNames() ##getSymbolObjects()
    RlistNameExpr <- nimbleGeneralParseDeparse(originalCode$caller$args[[1]])    
    for(i in seq_along(listElements)) {
        thisVarName <- listElements[i]
        ##      if(!inherits(originalCode$args[[i+1]], 'exprClass') ||  (originalCode$args[[i+1]]$name != "")){  ## skip first arg, which will be name of nlDef, then check if value is ""
        if(!is.null(originalCode$args[[thisVarName]])) {
            if(!inherits(originalCode$args[[thisVarName]], 'exprClass') ||  (originalCode$args[[thisVarName]]$name != "")){  ## skip first arg, which will be name of nlDef, then check if value is ""
                ## nfVar(A, 'x') for whichever element name it's on ('x')
                ##    RnfVarExprs[[exprCounter]] <- substitute(nfVar(A, X), list(A = RlistNameExpr, X = listElements[[i]]$name))
                RnfVarExprs[[exprCounter]] <- substitute(nfVar(A, X), list(A = RlistNameExpr, X = thisVarName))
                ## nfVar(A, 'x') <- y or whatever code was provided (already recursed for size processing)
                RnewExprs[[exprCounter]] <- substitute(A <- B, list(A = RnfVarExprs[[exprCounter]],
                                                                    B = nimbleGeneralParseDeparse(originalCode$args[[thisVarName]])))
                exprCounter <- exprCounter + 1
            }
        }
    }
    
    ## embed RnewExprs in a '{' expression
    if(length(RnewExprs) != 0) {
        RbracketNewExprs <- quote(after({}))
        RbracketNewExprs[[2]][2:(length(RnewExprs) + 1)] <- RnewExprs
        bracketNewExprs <- RparseTree2ExprClasses(RbracketNewExprs)
        ## Need to install assignment target in symTab if necessary so that it
        ## will be there for recursion in the following step
        assignmentTarget <- code$caller$args[[1]]
        if(assignmentTarget$isName) {
            if(!symTab$symbolExists(assignmentTarget$name, TRUE)) {
                symTab$addSymbol(symbolNimbleList(name = assignmentTarget$name, type = code$type, nlProc = code$sizeExprs$nlProc))
            }
        }
        ## recurse into element assignments
        exprClasses_setSizes(bracketNewExprs$args[[1]], symTab, typeEnv)
        asserts <- c(asserts, list(bracketNewExprs))
        if(length(code$args) > 1) ## always if we make it this far
            code$args <- code$args[1]
    }
  }
  return(asserts)
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

## size handler for nimArrayGeneral()
## nimArrayGeneral(typeCharString, nDim, c(sizeExpr1, ...), initializeValue, initializeLogical, unpackNDim(optional))
## nimArrayGeneral(     arg1,      arg2,       arg3,              arg4,            arg5       ,       arg6     )
sizeNimArrayGeneral <- function(code, symTab, typeEnv) {
    useArgs <- c(FALSE, FALSE, FALSE, TRUE, TRUE)
    if(length(code$args) > 5) useArgs <- c(useArgs, TRUE)
    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = useArgs)  ## recurse on initialValue and initialLogical only

    ## some checking
    if(inherits(code$args[[5]], 'exprClass'))
        if(!(code$args[[5]]$nDim == 0)) stop(exprClassProcessingErrorMsg(code, paste0('init argument to numeric, logical, integer, matrix or array must be scalar')), call. = FALSE)
    
    type <- code$args[[1]]    ## args[[1]]: 'type' argument
    nDim <- code$args[[2]]    ## args[[2]]: 'nDim' argument
    unpackNDim <- if(length(code$args) > 5) code$args[[6]] else FALSE

    cSizeExprs <- code$args[[3]] ## these are the size expressions encompassed by collectSizes(), needed for purposes of the C++ line to be generated
    if(!inherits(cSizeExprs, 'exprClass'))        stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (i) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)
    if(cSizeExprs$name != 'collectSizes')         stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (ii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)

    if(unpackNDim) { 
        asserts <- c(asserts, recurseSetSizes(cSizeExprs, symTab, typeEnv))
        if(!cSizeExprs$args[[1]]$isName)
            asserts <- c(asserts, sizeInsertIntermediate(cSizeExprs, 1, symTab, typeEnv)) ## this intermediate goes a layer down the AST, but works
        if(length(cSizeExprs$args[[1]]$sizeExprs) == 0) { ## The argument expression evaluates to scalar
            if(nDim == -1) {
                nDim <- 1
                code$args[[2]] <- 1
            }
            if(nDim == 1) unpackNDim <- FALSE             ## and that's ok because nDim given as 1
        }
        if(unpackNDim) {
            if(length(cSizeExprs$args[[1]]$sizeExprs) != 1) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (ii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)
            if(nDim == -1) {## code for nDim not given but dim given as expression
                if(!is.numeric(cSizeExprs$args[[1]]$sizeExprs[[1]])) stop()
                nDim <- cSizeExprs$args[[1]]$sizeExprs[1]
                if(nDim < 1 | nDim > 4) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (iii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)
                code$args[[2]] <- nDim
            }
            varName <- as.name(cSizeExprs$args[[1]]$name)
            for(i in 1:nDim) {               
                newIndexedSizeExpr <- RparseTree2ExprClasses( substitute(X[I], list(X = varName, I = i) ) )
                setArg(cSizeExprs, i, newIndexedSizeExpr)
            }
        }
    } else {
        asserts <- c(asserts, recurseSetSizes(cSizeExprs, symTab, typeEnv))
        nonScalarWhereNeeded <- FALSE
            if(inherits(cSizeExprs$args[[1]], 'exprClass'))
                if(cSizeExprs$args[[1]]$nDim != 0)
                    nonScalarWhereNeeded <- TRUE
        if(nDim == -1) {  ## nDim wasn't provided (to nimArray) and dim was an expression, so it out to be a scalar
            if(nonScalarWhereNeeded) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (iv) with sizes or dim to numeric, logical, integer, matrix or array.  It looks like dim argument was non-scalar but nDim was not provided.')), call. = FALSE)
            nDim <- code$args[[2]] <- 1 
        } else { ## call was from numeric, integer or logical
            if(nonScalarWhereNeeded) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (v) with sizes or dim to numeric, logical, integer, matrix or array.  It looks like length argument was non-scalar.')), call. = FALSE)
        }
        
    }
        
    ## if it is a call to matrix() and the value argument is non-scalar,
    ## we will generate it in C++ as nimNewMatrix
    useNewMatrix <- FALSE
    if(nDim == 2) {
        if(inherits(code$args[[4]], 'exprClass'))
            if(code$args[[4]]$nDim > 0)
                useNewMatrix <- TRUE
    }
                
     if(code$args[[2]] != length(cSizeExprs$args)) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (iii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)

    annotationSizeExprs <- lapply(cSizeExprs$args, nimbleGeneralParseDeparse) ## and this is for purposes of the sizeExprs in the AST exprClass object
    
    missingSizes <- unlist(lapply(cSizeExprs$args, identical, as.numeric(NA)))
    ## only case where we do something useful with missingSizes is matrix(value = non-scalar, ...)
    if(any(missingSizes)) {
        ## modify sizes in generated C++ line
        if(useNewMatrix) cSizeExprs$args[missingSizes] <- -1
        else cSizeExprs$args[missingSizes] <- 1

        ## modify annotation sizeExprs
        totalInputLengthExpr <- if(inherits(code$args[[4]], 'exprClass')) productSizeExprs(code$args[[4]]$sizeExprs) else 1 ## should always be exprClass in here anyway
        ## see newMatrixClass in nimbleEigen.h
        if(missingSizes[1]) { ## missing nrow
            if(missingSizes[2]) { ## missing both
                annotationSizeExprs[[1]] <- totalInputLengthExpr
                annotationSizeExprs[[2]] <- 1
            } else { ## ncol provided
                annotationSizeExprs[[1]] <- substitute(floor( ((A)-1) / (B)) + 1, ## avoids errors from integer arithmetic
                                              list(A = totalInputLengthExpr,
                                                   B = annotationSizeExprs[[2]]))
            }
        } else { ## nrow provided, ncol missing (is both provided, we wouldn't be in this code
                annotationSizeExprs[[2]] <- substitute(floor( ((A)-1) / (B)) + 1,
                                              list(A = totalInputLengthExpr,
                                                   B = annotationSizeExprs[[1]]))
        }
    }
    
    asserts <- c(asserts, recurseSetSizes(cSizeExprs, symTab, typeEnv))
    if(!(type %in% c('double', 'integer', 'logical')))       stop('unknown type in nimArrayGeneral')
    code$name <- 'initialize' ## may be replaced below if useNewMatrix
    if(inherits(code$args[[4]], 'exprClass'))
        if(code$args[[4]]$nDim > 0)
            code$name <- 'assignNimArrToNimArr'
    code$args <- c(code$args[4:5], cSizeExprs$args)  ##  args: initialize(initializeValue, initializeLogical, sizeExpr1, sizeExpr2, etc...)
    for(i in seq_along(code$args)) {
        if(inherits(code$args[[i]], 'exprClass')) {
            code$args[[i]]$callerArgID <- i
            code$args[[i]]$caller <- code
        }
    }
    code$type <- type
    code$nDim <- nDim
    code$toEigenize <- 'no'
        ## insert intermediate unless it will be newNimMatrix
    code$sizeExprs <- annotationSizeExprs

    if(!useNewMatrix) 
        if(inherits(code$caller, 'exprClass'))
            if(!(code$caller$name %in% assignmentOperators)) {
                if(!is.null(code$caller$name)) {
                    asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
                }
            } else
                typeEnv$.ensureNimbleBlocks <- TRUE
    
    ## check for nimNewMatrix case
    if(useNewMatrix) {
        suffix <- 'D'
        if(code$type == 'integer') suffix <- 'I'
        if(code$type == 'logical') suffix <- 'B'
        code$name <- paste0("nimNewMatrix", suffix)
        code$toEigenize <- "yes"
    } else {
        ## otherwise, lift values arg if necessary
        if(inherits(code$args[[1]], 'exprClass')) ## was re-ordered here
            if(!(code$args[[1]]$isName))
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
    }
    
    return(asserts)
}

sizeRunTime <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1) stop(exprClassProcessingErrorMsg(code, paste0('run.time must take exactly 1 argument')), call. = FALSE)
    origCaller <- code$caller
    origCallerArgID <- code$callerArgID

    if(!code$caller$isAssign) { ## e.g. a + run.time({foo(y)}), should have already been lifted by buildIntermediates
        message('Problem in sizeRunTime: run.time is not in a simple assignment at this stage of processing.')
    }

    ## this is the case ans <- run.time({foo(y)})
    lhsName <- code$caller$args[[1]]$name
    timerName <- IntermLabelMaker() ##paste0(lhsName,'_TIMER_')
    newSym <- symbolNimbleTimer(name = timerName, type = 'symbolNimbleTimer')
    symTab$addSymbol(newSym)
    startTimerAssert <- RparseTree2ExprClasses(substitute(startNimbleTimer(TIMERNAME), list(TIMERNAME = as.name(timerName))))
    recurseAsserts <- recurseSetSizes(code, symTab, typeEnv) ## arg to run.time should be in {} so any nested asserts should be done by the time this finishes and this should return NULL
    if(!is.null(recurseAsserts)) {
        message('issue in sizeRunTime: recurseAsserts is not NULL')
    }
    asserts <- list(startTimerAssert, code$args[[1]])
    newCode <- RparseTree2ExprClasses(substitute(endNimbleTimer(TIMERNAME), list(TIMERNAME = as.name(timerName))))
    newCode$type <- 'double'
    newCode$nDim <- 0
    newCode$sizeExprs <- list()
    newCode$toEigenize <- 'no'
    setArg(origCaller, origCallerArgID, newCode)
    return(asserts)
}

sizeGetParam <- function(code, symTab, typeEnv) {
    if(length(code$args) > 3) {
        asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, FALSE, FALSE, rep(TRUE, length(code$args)-3)))
        for(i in 4:length(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
                if(code$args[[i]]$toEigenize=='yes') stop(exprClassProcessingErrorMsg(code, 'In sizeGetParam: There is an expression beyond the third argument that cannot be handled.  If it involve vectorized math, you need to do it separately, not in this expression.'), call. = FALSE)
            }
        }
    } else {
        asserts <- list()
    }
 
    
    paramInfoSym <- symTab$getSymbolObject(code$args[[3]]$name, inherits = TRUE)
    code$type <- paramInfoSym$paramInfo$type
    code$nDim <- paramInfoSym$paramInfo$nDim
    code$sizeExprs <- vector(mode = 'list', length = code$nDim)
    code$toEigenize <- 'no'

    if(!(code$caller$name %in% assignmentOperators)) {
        if(!is.null(code$caller$name))
            if(!(code$caller$name == '{')) ## could be on its own line -- useless but possible
                asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        code$toEigenize <- 'maybe'
    } else
        typeEnv$.ensureNimbleBlocks <- TRUE
    return(asserts)
}

sizeGetBound <- function(code, symTab, typeEnv) {
    if(length(code$args) > 3) {
        asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, FALSE, FALSE, rep(TRUE, length(code$args)-3)))
        for(i in 4:length(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
                if(code$args[[i]]$toEigenize=='yes') stop(exprClassProcessingErrorMsg(code, 'In sizeGetParam: There is an expression beyond the third argument that cannot be handled.  If it involve vectorized math, you need to do it separately, not in this expression.'), call. = FALSE)
            }
        }
    } else {
        asserts <- list()
    }
 
    
    boundInfoSym <- symTab$getSymbolObject(code$args[[3]]$name, inherits = TRUE)
    code$type <- boundInfoSym$boundInfo$type
    code$nDim <- boundInfoSym$boundInfo$nDim
    code$sizeExprs <- vector(mode = 'list', length = code$nDim)
    code$toEigenize <- 'no'

    if(!(code$caller$name %in% assignmentOperators)) {
        if(!is.null(code$caller$name))
            if(!(code$caller$name == '{')) ## could be on its own line -- useless but possible
                asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        code$toEigenize <- 'maybe'
    }
    return(asserts)
}


sizeSwitch <- function(code, symTab, typeEnv) {
    if(length(code$args) <= 2) return(invisible(NULL))
    for(i in 3:length(code$args)) { ## just like the '{' clause of exprClasses_setSizes.  This treats each of the outcomes as if it was a new line or block of code
        if(inherits(code$args[[i]], 'exprClass')) {
            newAsserts <- exprClasses_setSizes(code$args[[i]], symTab, typeEnv)
            code$args[[i]]$assertions <- if(is.null(newAsserts)) list() else newAsserts
        }
    }
    return(invisible(NULL))
}

sizeAsRowOrCol <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    if(!inherits(a1, 'exprClass')) stop(exprClassProcessingErrorMsg(code, 'In sizeAsRowOrCol: Argument must be an expression.'), call. = FALSE)
    if(a1$nDim == 0) stop(exprClassProcessingErrorMsg(code, 'In sizeAsRowOrCol:  Argument cannot be scalar (could be fixed).'), call. = FALSE)
    code$type <- a1$type
    code$toEigenize <- 'yes'
    if(!code$name %in% c('asRow', 'asCol')) stop(exprClassProcessingErrorMsg(code, 'Somehow the system got to sizeAsRowOrCol without a call to asRow or asCol.'), call. = FALSE)

    if(a1$nDim == 1) {
        if(code$name == 'asRow') {
            code$nDim <- 2
            code$sizeExprs <- c(list(1), a1$sizeExprs[[1]])
        } else {
            code$nDim <- 2
            code$sizeExprs <- c(a1$sizeExprs[[1]], list(1))
        }
        return(asserts)
    }

    warning(paste0(' asRow or asCol used on something with more than 1 dimension in ', nimDeparse(code)), call. = FALSE)

}


recurseExtractNimListArg <- function(code, symTab){
  if(length(code$args) == 0){  ## reached end level list
    listSym <- symTab$getSymbolObject(code$name, inherits = TRUE)
  }
  else{ ## progress though lower lists, passing appropriate symbol down each time
    listSym <- recurseExtractNimListArg(code$args[[1]], symTab)
    nestedObjInd <- which(sapply(listSym$nlProc$neededTypes, function(x){return(x$name == code$args[[2]])}) == TRUE)
    if(length(nestedObjInd) == 0)
      listSym <- listSym$nlProc$symTab$getSymbolObject(code$args[[2]])
    else
      listSym <- listSym$nlProc$neededTypes[[nestedObjInd]]
  }
    return(listSym)
}

## a$b becomes nfVar(a, 'b')
sizeNFvar <- function(code, symTab, typeEnv) {
    asserts <- list()
    if(!inherits(code$args[[1]], 'exprClass'))
        stop(exprClassProcessingErrorMsg(code, 'Problem using $: no name on the right?'), call. = FALSE)
    if(length(code$args) != 2)
        stop(exprClassProcessingErrorMsg(code, 'Problem using $: wrong number of arguments?'), call. = FALSE)
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    if(code$args[[1]]$isName) {
        objectName <- code$args[[1]]$name
        symbolObject <- symTab$getSymbolObject(objectName, inherits = TRUE)
        objectType <- symbolObject$type
    } else { ## if there is nesting, A$B$C, figure out what to do
        objectType <- code$args[[1]]$type
        symbolObject <- code$args[[1]]$sizeExprs ## repurposed for this role
    }
    
    isSymFunc <- objectType == 'nimbleFunction'   ## minor inconsistency in naming style here
    isSymList <- objectType == 'symbolNimbleList'
    
    ## Cases to handle (nl for nimbleList, nf for nimbleFunction):
    ## nl$a <- x     ## NimArr assignment (not setSize needed)
    ## nl$a <- x + 1 ## eigen assignment  (setSize needed)
    ## nl1$nl2$ <- x or x + 1
    ## x <- foo(nl$a)
    ## x <- foo(nl1$nl2$b)
    ## same with nf instead of any nl, in any order
    ## nl$new()$a , which becomes makeNewNimbleListObject(nl1)$a
    ## nl in nlEigenReferenceList

    if(!(isSymFunc || isSymList))
        stop(exprClassProcessingErrorMsg(code, 'In sizeNFvar: First argument is not a nimbleFunction or a nimbleList'), call. = FALSE)

    nfProc <- if(isSymFunc) symbolObject$nfProc else symbolObject$nlProc
    
    if(is.null(nfProc)) {
         browser()
        stop(exprClassProcessingErrorMsg(code, 'In handling X$Y: Symbols for X have not been set up.'), call. = FALSE)
    }
    
    memberName <- code$args[[2]]
    if(!is.character(memberName)) stop(exprClassProcessingErrorMsg(code, 'In handling X$Y: Something is wrong with Y.'), call. = FALSE)

    memberSymbolObject <- nfProc$getSymbolTable()$getSymbolObject(memberName)
    if(!is.null(memberSymbolObject)) code$type <- memberSymbolObject$type
    
    if(isSymList | isSymFunc) {
        ## symbolNimbleList
        ## We need (*nl) in C++, represented by cppPointerDereference(nl)
        if(code$args[[1]]$name != 'cppPointerDereference') {
            ## a1 <- nimble:::insertExprClassLayer(code, 1, 'cppPointerDereference')
            ## a1$type <- a1$args[[1]]$type
            ## a1$nDim <- a1$args[[1]]$nDim ## may be (always?) uninitialized
            ## a1$sizeExprs <- a1$args[[1]]$sizeExprs ## may be list()
            a1 <- nimble:::insertExprClassLayer(code, 1, 'cppPointerDereference',
                                                type = code$args[[1]]$type,
                                                nDim = code$args[[1]]$nDim,
                                                sizeExprs = code$args[[1]]$sizeExprs)
        }
    }
     
    ## following checks are on type of A$B (isSymList and isSymFunc refer to type of A)
    
    if(code$type == 'symbolNimbleList') {
        ## for a nimbleList, sizeExprs slot is used for symbol object
        ## of nlGenerator of member object
        code$sizeExprs <- memberSymbolObject
    } else if(code$type == 'nimbleFunction') {
        ## nimbleFunction
        code$sizeExprs <- memberSymbolObject
    } else if(code$type == 'nimbleFunctionList') {
        code$sizeExprs <- memberSymbolObject
    } else {
        ## a numeric etc. type
        code$nDim <- memberSymbolObject$nDim
        code$sizeExprs <- if(code$nDim > 0)
                              makeSizeExpressions(memberSymbolObject$size,
                                                  parse(text = nimDeparse(code))[[1]])
                          else
                              list()
    }
    return(asserts)
}


sizeNFvar_old <- function(code, symTab, typeEnv) {
  topLevel <- code$caller$name != 'nfVar'
  nfName <- code$args[[1]]$name
  if(nfName == 'cppPointerDereference'){
    tmpArg <- code$args[[1]]$args[[1]]
    nfName <- tmpArg$name
    while(nfName == 'cppPointerDereference'){
      tmpArg <- tmpArg$args[[1]]
      nfName <- tmpArg$name
    }
  }
  else if(nfName == 'makeNewNimbleListObject'){
    nfName <- code$args[[1]]$args[[1]]$name
  }
  asserts <- NULL
  code$toEigenize <- 'maybe'
  if(nfName == 'nfVar'){ ## accessing nested nimbleList or nested nimbleList element
    isSymList <- TRUE
    objSym <- recurseExtractNimListArg(code, symTab)
    if(is.null(objSym)) stop(exprClassProcessingErrorMsg(code, 'In sizeNFvar: Symbol not found in the nimbleFunction.'), call. = FALSE)
    asserts <- recurseSetSizes(code, symTab, typeEnv)
  }
  else{
    nfSym <- symTab$getSymbolObject(nfName, inherits = TRUE)
    isSymFunc <- inherits(nfSym, 'symbolNimbleFunction')
    isSymList <- (inherits(nfSym, 'symbolNimbleList') || inherits(nfSym, 'symbolNimbleListGenerator'))
    eigListFuncNames <- sapply(nlEigenReferenceList, function(x){return(x$nimFuncName)})
    if(nfName %in% eigListFuncNames || code$args[[1]]$name == 'makeNewNimbleListObject'){
      asserts <- recurseSetSizes(code, symTab, typeEnv)
      asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
    }
    if(!(isSymFunc || isSymList))
      stop(exprClassProcessingErrorMsg(code, 'In sizeNFvar: First argument is not a nimbleFunction or a nimbleList'), call. = FALSE)
    if(isSymFunc) nfProc <- nfSym$nfProc ## Now more generally this should be an interface
    if(isSymList) nfProc <- nfSym$nlProc
    if(is.null(nfProc)) stop(exprClassProcessingErrorMsg(code, 'In sizeNFvar: Symbols in this nimbleFunction generation function not set up.'), call. = FALSE)
    objName <- code$args[[2]]
    if(!is.character(objName)) stop(exprClassProcessingErrorMsg(code, 'In sizeNFvar: Second argument must be a character string.'), call. = FALSE)
    objSym <- nfProc$getSymbolTable()$getSymbolObject(objName)  ##nfProc$setupSymTab$getSymbolObject(objName)
    if(is.null(objSym)) stop(exprClassProcessingErrorMsg(code, 'In sizeNFvar: Symbol not found in the nimbleFunction.'), call. = FALSE)
    if(inherits(objSym, 'symbolNimbleList')) code$toEigenize <- 'no'
  }
  if(!is.null(objSym)) code$type <- objSym$type
  if(code$type != 'symbolNimbleList') code$nDim <- objSym$nDim
  if(isSymList){
    if(code$args[[1]]$name != 'cppPointerDereference'){
      a1 <- nimble:::insertExprClassLayer(code, 1, 'cppPointerDereference')
      a1$type <- a1$args[[1]]$type
      a1$nDim <- a1$args[[1]]$nDim
      a1$sizeExprs <- a1$args[[1]]$sizeExprs
      code$args[[1]] <- a1
    }
  }
  else{
    code$nDim <- objSym$nDim
    code$type <- objSym$type
  }
  if((code$type != 'symbolNimbleList') && code$nDim > 0) {
    code$sizeExprs <- makeSizeExpressions(objSym$size,
                                          parse(text = nimDeparse(code))[[1]])
  } 
  else if(code$type == 'symbolNimbleList'){
    code$sizeExprs$nlProc <-objSym$nlProc
  }
  else{
    code$sizeExprs <- list()
  }
  return(asserts)
}

sizeDoubleBracket <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(code$args[[1]]$isName) {
        objectName <- code$args[[1]]$name
        symbolObject <- symTab$getSymbolObject(objectName, inherits = TRUE)
        objectType <- symbolObject$type
    } else { ## if there is nesting, A$B$C, figure out what to do
        objectType <- code$args[[1]]$type
        symbolObject <- code$args[[1]]$sizeExprs ## repurposed for this role
    }
    isSymFuncList <- objectType == 'nimbleFunctionList'
    if(!isSymFuncList) stop('nfList[[i]] must use a nimbleFunctionList')
    code$sizeExprs <- symbolObject
    code$type <- objectType
    return(if(is.null(asserts)) list() else asserts)
}

sizeChainedCall <- function(code, symTab, typeEnv) { ## options include nfMethod(nf, 'foo')(a), or nfMethod(nf[[i]], 'foo')(a) [which arises from nf[[i]]$foo(a), where nf is a local nflist, where nf could need recursion, in which case it will be wrapped in nfVar
    ## In other places we generate chainedCalls for static_cast<int>(a), but those shouldn't be seen here
    a1 <- code$args[[1]] 
    if(!inherits(a1, 'exprClass')) stop(exprClassProcessingErrorMsg(code, 'In sizeChainedCall.  First arg is not an expression.'), call. = FALSE)
    nfMethodRCobj <- NULL

    if(a1$name != 'nfMethod') stop(exprClassProcessingErrorMsg(code, 'Some problem processing a chained call.'), call. = FALSE)

    asserts <- recurseSetSizes(a1, symTab, typeEnv, useArgs = c(TRUE, rep(FALSE, length(a1$args)-1)))

    a11 <- a1$args[[1]]
    methodName <- a1$args[[2]]

    if(a1$args[[1]]$isName) {
        objectName <- a1$args[[1]]$name
        symbolObject <- symTab$getSymbolObject(objectName, inherits = TRUE)
        objectType <- symbolObject$type
    } else { ## if there is nesting, A$B$C, figure out what to do
        objectType <- a1$args[[1]]$type
        symbolObject <- a1$args[[1]]$sizeExprs ## repurposed for this role
    }

    isSymFun <- objectType == 'nimbleFunction'
    isSymFunList <- objectType == 'nimbleFunctionList'
    
    if(! (isSymFun | isSymFunList)) stop('Problem processing what looks like a member function call.')


    if(!is.character(methodName)) stop(exprClassProcessingErrorMsg(code, 'In handling X$Y: Something is wrong with Y.'), call. = FALSE)

    if(isSymFun) {
        nfProc <- symbolObject$nfProc 
        if(is.null(nfProc)) {
            stop(exprClassProcessingErrorMsg(code, 'In handling X$Y(): Symbols for X have not been set up.'), call. = FALSE)
        }
        if(a1$args[[1]]$name != 'cppPointerDereference') {
            nimble:::insertExprClassLayer(a1, 1, 'cppPointerDereference') ## not annotated, but not needed
        }

    }

    ##sym <- symTab$getSymbolObject(a11$name, TRUE)
    
    if(isSymFun) {
        nfMethodRCobj <- nfProc$getMethodInterfaces()[[methodName]] ##sym$nfProc$origMethods[[methodName]]
        returnType <- nfProc$compileInfos[[methodName]]$returnSymbol
    } 
    if(isSymFunList) {
        nfMethodRCobj <- getFunctionEnvVar(nf_getGeneratorFunction(symbolObject$baseClass), 'methodList')[[methodName]]
        returnType <- nfMethodRCobj$returnType
    }
    ## } else {
    ##     if(a11$name != '[[') stop(exprClassProcessingErrorMsg(code, 'In sizeChainedCall. Expecting a nimbleFunction list or a nimFun as first arg of nfMethod.'), call. = FALSE)
    ##     ## should look like nfMethod(nflist[[i]], 'foo')
    ##     a111 <- a11$args[[1]]
    ##     sym <- symTab$getSymbolObject(a111$name, TRUE)
    ##     if(!inherits(sym, 'symbolNimbleFunctionList')) {
    ##         stop(exprClassProcessingErrorMsg(code, 'In sizeChainedCall. Expecting a nimbleFunction list.'), call. = FALSE)
    ##     }
    ##     nfMethodRCobj <- getFunctionEnvVar(nf_getGeneratorFunction(sym$baseClass), 'methodList')[[methodName]]
    ## }
    
 ##   }
 ##   else warning(paste0('Warning that we did not know what to do in sizeChainedCall for ', nimDeparse(code)))
    if(!is.null(nfMethodRCobj)) {
     ##   returnType <- nfMethodRCobj$returnType
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
##    if(length(sym$lengthName)==0) stop(paste0("Error the size information for ", nimDeparse(code), " is missing."), call. = FALSE) 
    ##    code$sizeExprs <- list(as.name(sym$lengthName))
    indexRangeCase <- FALSE
    if(length(code$args) == 1) {  # full vector of nodes
        code$sizeExprs <- list(substitute(cppMemberFunction(getTotalLength(ACCESSNAME)), list(ACCESSNAME = as.name(code$args[[1]]$name))))
        asserts <- list()
    } else {  # there must be index on the node
        asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, rep(TRUE, length(code$args)-1)))
        if(is.numeric(code$args[[2]])) {
            code$sizeExprs <- list(substitute(cppMemberFunction(getNodeLength(ACCESSNAME, ACCESSINDEX)), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = code$args[[2]])))
        } else {
            if(code$args[[2]]$nDim > 0)
                if(!(code$args[[2]]$isName))
                    asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
            code$sizeExprs <- list(substitute(getNodesLength_Indices(ACCESSNAME, ACCESSINDEX), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = as.name(code$args[[2]]$name))))
            indexRangeCase <- TRUE
        }
    }
   
    if(code$caller$name == "[" & code$caller$callerArgID == 1) # values(...)[.] <-
        if(typeEnv$.AllowUnknowns) ## a surrogate for being on LHS of an assignment. values(...)[] should work on RHS
            stop(exprClassProcessingErrorMsg(code, 'In sizeValues: indexing of values() on left-hand size of an assignment is not allowed.'), call. = FALSE)
     if(code$caller$name %in% assignmentOperators) {
        if(code$callerArgID == 2) { ## ans <- values(...)
            code$name <- if(!indexRangeCase) 'getValues' else 'getValuesIndexRange'
            LHS <- code$caller$args[[1]]
            if(LHS$isName) { ## It is a little awkward to insert setSize here, but this is different from other cases in sizeAssignAfterRecursing
                assertSS <- list(substitute(setSize(LHS), list(LHS = as.name(LHS$name))))
                if(length(code$args) == 1) {  # full vector of nodes
                    assertSS[[1]][[3]] <- substitute(cppMemberFunction(getTotalLength(ACCESSNAME)), list(ACCESSNAME = as.name(code$args[[1]]$name)))
                } else {  # there must be index on the node
                    if(is.numeric(code$args[[2]])) {
                        assertSS[[1]][[3]] <- substitute(cppMemberFunction(getNodeLength(ACCESSNAME, ACCESSINDEX)), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = code$args[[2]]))
                    } else {
                        assertSS[[1]][[3]] <- substitute(getNodesLength_Indices(ACCESSNAME, ACCESSINDEX), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = as.name(code$args[[2]]$name))) ## intermediate has already been inserted above, if needed
                    }
                }

                # assertSS[[1]][[3]] <- substitute(cppMemberFunction(getTotalLength(ACCESSNAME)), list(ACCESSNAME = as.name(code$args[[1]]$name)))
                asserts <- c(asserts, assertSS)
            } else
                typeEnv$.ensureNimbleBlocks <- TRUE
        }   # values(...) <- P, don't change it
    } else { ## values(...) embedded in a RHS expression
        code$name <- if(!indexRangeCase) 'getValues' else 'getValuesIndexRange'
        code$toEigenize <- 'yes' ## This tricks sizeAssignAfterRecursing to generate the setSize in asserts, in getValues case (getValuesIndexRange is in set of names to skip for that)
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
        stop(exprClassProcessingErrorMsg(code, 'In sizeNimbleFunction: A nimbleFunction method should not be processed here.'), call. = FALSE)
        ## HANDLING OF myNF$run() HERE IS DEFUNCT.  ALL SHOULD GO THROUGH sizeChainedCall now (chainedCall(nfMethod(myNF,'run'), arg1, arg2).
        ## nfMethodRCobj <- sym$nfProc$getMethodInterfaces()$run ##environment(sym$nfProc$nfGenerator)$methodList$run
        ## returnType <- nfMethodRCobj$returnType
        ## eigListClasses <- sapply(nlEigenReferenceList, function(x){return(x$className)})  
        ## if(!(as.character(returnType[1]) %in% c('double', 'integer', 'character', 'logical', 'void',
        ##                                         eigListClasses))){  
        ##   ## if we have a nl return type, find class name and match with nlGenerator in symTab ## wrong 'return' ???
        ##   outClassName <- get('return', envir = typeEnv)$sizeExprs$name
        ##   parentNLGenName <- lapply(symTab$parentST$symbols, function(x){
        ##     symType <- x$type
        ##     if(symType == 'Ronly'){
        ##       symClassName <- x[['nlProc']][['name']]
        ##       if(!is.null(symClassName) && symClassName == outClassName){
        ##         return(x$name)
        ##       }
        ##     }
        ##     return(NULL)
        ##   })
        ##   returnType[[1]] <- as.name(unlist(parentNLGenName))
        ## }
        ## argInfo <- nfMethodRCobj$argInfo
        ## ok <- TRUE
    }
    if(inherits(sym, 'symbolMemberFunction')) {
        eigListClasses <- sapply(nlEigenReferenceList, function(x){return(x$className)})  
        returnType <- sym$nfMethodRCobj$returnType ## now nfMethodRCobj could be an interface
        if(!(as.character(returnType[1]) %in% c('double', 'integer', 'character', 'logical', 'void',
                                                eigListClasses))){  
          ## if we have a nl return type, find class name and match with nlGenerator in symTab
          outClassName <- get('return', envir = typeEnv)$sizeExprs$name
          parentNLGenName <- lapply(symTab$parentST$symbols, function(x){
            symType <- x$type
            if(symType == 'Ronly'){
              symClassName <- x[['nlProc']][['name']]
              if(!is.null(symClassName) && symClassName == outClassName){
                return(x$name)
              }
            }
            return(NULL)
          })
          returnType[[1]] <- as.name(unlist(parentNLGenName))
        }
        argInfo <- sym$nfMethodRCobj$argInfo
        ok <- TRUE
    }
    if(ok) {
        asserts <- generalFunSizeHandler(code, symTab, typeEnv, returnType, argInfo)
        return(asserts)
    }
    stop(exprClassProcessingErrorMsg(code, 'In sizeNimbleFunction: The function name is not known and is not a nimbleFunction or a member function.'), call. = FALSE)
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
    if(!code$args[[1]]$isName) stop(exprClassProcessingErrorMsg(code, 'In sizeOneEigenCommand:  First arg should be a name.'), call. = FALSE)
    recurseSetSizes(code, symTab, typeEnv)
    if(code$args[[1]]$nDim != 2) stop(exprClassProcessingErrorMsg(code, 'In sizeOneEigenCommand:  At the moment only works for 2D objects.'), call. = FALSE)
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

nimbleGeneralParseDeparse <- function(code) {
    if(inherits(code,'exprClass'))
        parse(text = nimDeparse(code), keep.source = FALSE)[[1]]
    else
        code
}


sizeSetSize <- function(code, symTab, typeEnv) {
    #go inside nfVar call if resizing nimbleList element
    if(code$args[[1]]$name == 'nfVar'){
      useArg1 <- TRUE
      sym <- symTab$getSymbolObject(code$args[[1]]$args[[1]]$name)
      if(sym$type == 'symbolNimbleList'){
        sym <- sym$nlProc$symTab$getSymbolObject(code$args[[1]]$args[[2]])
      }
    } else {
      sym <- symTab$getSymbolObject(code$args[[1]]$name, inherits = TRUE)
      useArg1 <- FALSE
    } 

    if(!inherits(sym, 'symbolNumericList')) {
        if(sym$nDim == 0) stop(exprClassProcessingErrorMsg(code, 'In sizeSetSize: Resizing a scalar does not make sense.'), call. = FALSE)
        firstSizeExpr <- code$args[[2]]
        if(inherits(firstSizeExpr, 'exprClass')) {
            if(firstSizeExpr$name == 'nimC') { ## handle syntax of resize(Z, c(3, dim(A)[1]))
                if(length(firstSizeExpr$args) != sym$nDim) stop(exprClassProcessingErrorMsg(code, 'In sizeSetSize: Problem with number of dimensions provided in resize.'), call. = FALSE)
                asserts <- recurseSetSizes(firstSizeExpr, symTab, typeEnv) ## may set intermediates if needed
                ## see comment below
                ##                assign(code$args[[1]]$name, exprTypeInfoClass$new(nDim = sym$nDim, sizeExprs = lapply(firstSizeExpr$args, nimbleGeneralParseDeparse), type = sym$type), envir = typeEnv)
                for(i in 1:length(firstSizeExpr$args)) {
                    code$args[[i+1]] <- firstSizeExpr$args[[i]]
                    if(inherits(firstSizeExpr$args[[i]], 'exprClass')) {
                        firstSizeExpr$args[[i]]$caller <- code
                        firstSizeExpr$args[[i]]$callerArgID <- i+1
                    }
                }
                return(if(length(asserts)==0) NULL else asserts)
            }
        }

        asserts <- recurseSetSizes(code, symTab, typeEnv, c(useArg1, rep(TRUE, sym$nDim) ) )

        if(inherits(code$args[[2]], 'exprClass')) {
            if(code$args[[2]]$nDim > 0) {
                if(length(code$args) > 2) stop(exprClassProcessingErrorMsg(code, 'In sizeSetSize: Non-scalar argument for sizes provided, but more than one size argument also provided.  This does not look valid.'), call. = FALSE)
        
                if(!(code$args[[2]]$isName)) asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
                code$name <- 'setSizeNimArrToNimArr'
            } else {
                if(length(code$args) != 1 + sym$nDim) stop(exprClassProcessingErrorMsg(code, 'In sizeSetSize: Problem with number of dimensions provided in setSize.'), call. = FALSE)
            }
        }

        
        ## We used to update typeEnv here with the new sizes, but it is not safe to do so because the setSize might appear inside a conditional (if-then)
        ## and hence one can't know until run-time if the size will actually be changed as given.  Thus typeEnv sizeExprs are set when a variable first appears
        ## and should be either constants (and not ever setSized again, which we should check for but don't) or remain generic (dim(x)[1], etc)
        ## assign(code$args[[1]]$name, exprTypeInfoClass$new(nDim = sym$nDim, sizeExprs = lapply(code$args[-1], nimbleGeneralParseDeparse), type = sym$type), envir = typeEnv)
        return(if(length(asserts)==0) NULL else asserts)
    }
    if(inherits(sym, 'symbolNumericList') ) { ## these are deprecated
    	if(length(code$args) != 2 + sym$nDim) stop(exprClassProcessingErrorMsg(code, 'In sizeSetSize: Problem with number of dimensions provided in resize.'), call. = FALSE)
        ## no longer modify typeEnv
        ##    	assign(code$name, exprTypeInfoClass$new(nDim = sym$nDim, sizeExprs = lapply(code$args[-1], nimbleGeneralParseDeparse), type = sym$type), envir = typeEnv)
    	invisible(NULL)
    }
}


## This was redundant and we should eventually be able to remove it
sizeResizeNoPtr <- function(code, symTab, typeEnv){
    sym <- symTab$getSymbolObject(code$args[[1]]$name, inherits = TRUE)
    if(length(code$args[[2]]) != 1)  stop(exprClassProcessingErrorMsg(code, 'In sizeResizeNoPtr: Problem with number of dimensions provided in resize.'), call. = FALSE)
    ## no longer modify typeEnv
    ## assign(code$name, exprTypeInfoClass$new(nDim = 1, sizeExprs = lapply(code$args[-1], nimDeparse), type = sym$type), envir = typeEnv)
    invisible(NULL)
}


## Handler for for-loops: a fairly special case
## e.g. for(i in 1:10) {do(i)}
sizeFor <- function(code, symTab, typeEnv) {
    if(length(code$args) != 3) stop('Error in sizeFor: expected 3 arguments to a for-loop', call. = FALSE)
    ## first handle type of the indexing variable
    if(!inherits(code$args[[2]], 'exprClass')) stop(exprClassProcessingErrorMsg(code, 'In sizeFor: expected the index range to be an expression (exprClass).'), call. = FALSE)
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
    return(if(length(asserts) == 0) invisible(NULL) else asserts)
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
    ##asserts <- recurseSetSizes(code, symTab, typeEnv)
    typeEnv$.AllowUnknowns <- FALSE
    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, TRUE))
    typeEnv$.AllowUnknowns <- TRUE
    if(length(code$args) > 2){
      asserts <- c(asserts, exprClasses_setSizes(code, symTab, typeEnv))
    }
    else{
      asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, useArgs = c(TRUE, FALSE)))
      typeEnv[['.ensureNimbleBlocks']] <- FALSE ## may have been true from RHS of rmnorm etc.
      asserts <- c(asserts, sizeAssignAfterRecursing(code, symTab, typeEnv))
    }
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
        }
        else if(is.character(RHS)){
          RHSname = ''
          RHSnDim <- 0
          RHStype <- 'character'
          RHSsizeExprs <- list() 
        }
        else {
            stop(exprClassProcessingErrorMsg(code, "In sizeAssignAfterRecursing: don't know what to do with a provided expression."), call. = FALSE)
        }
    }
    test <- try(if(inherits(RHStype, 'uninitializedField') | length(RHStype)==0) {
        stop(exprClassProcessingErrorMsg(code, paste0("In sizeAssignAfterRecursing: '",RHSname, "' is not available or its output type is unknown.")), call. = FALSE)
    })
    if(inherits(test, 'try-error')) browser()
    if(LHS$isName) {
        if(!exists(LHS$name, envir = typeEnv, inherits = FALSE)) { ## not in typeEnv
            ## If LHS unknown, create it in typeEnv
            if(!symTab$symbolExists(LHS$name, TRUE)) { ## not in symTab
                if(RHStype %in% c('double','integer', 'logical')) {  ## valid type to create here
                    ## We used to delay creating sizeExprs until below, but now it always generic
                    ## assign(LHS$name, exprTypeInfoClass$new(nDim = RHSnDim, type = RHStype), envir = typeEnv)
                    assign(LHS$name, exprTypeInfoClass$new(nDim = RHSnDim, type = RHStype, sizeExprs = makeSizeExpressions(rep(NA, RHSnDim), LHS$name)), envir = typeEnv)
                    symTab$addSymbol(symbolBasic(name = LHS$name, nDim = RHSnDim, type = RHStype))
                } else { ## not valid type to create here
                    if(RHStype == 'voidPtr') {
                        ## This should be ok without sizeExprs content
                        assign(LHS$name, exprTypeInfoClass$new(nDim = RHSnDim, type = RHStype), envir = typeEnv)
                        symTab$addSymbol(symbolVoidPtr(name = LHS$name, type = RHStype))
                    } 
                    else if(RHStype == "symbolNimbleList") {
                        ## I think we have the nlProc in the RHS sizeExprs in some cases?
                      LHSnlProc <- symTab$getSymbolObject(RHS$name)$nlProc
                      if(is.null(LHSnlProc)) LHSnlProc <- RHS$sizeExprs$nlProc
                      if(is.null(LHSnlProc)) LHSnlProc <- symTab$getSymbolObject(RHS$name, inherits = TRUE)$nlProc
                      symTab$addSymbol(symbolNimbleList(name = LHS$name, type = RHStype, nlProc = LHSnlProc))
                    }
                    else if(symTab$symbolExists(RHStype, TRUE)){  ## this is TRUE if a nested nimbleFunction returns a nimbleList - the type of
                                                                  ## the returned nimbleList will be a symbolNimbleListGenerator that exists
                                                                  ## in the parent ST.
                      LHSnlProc <- symTab$getSymbolObject(RHStype, TRUE)$nlProc
                      symTab$addSymbol(symbolNimbleList(name = LHS$name, type = 'symbolNimbleList', nlProc = LHSnlProc))
                    }
                    else
                        stop(exprClassProcessingErrorMsg(code, paste0('In sizeAssignAfterRecursing: LHS is not in typeEnv or symTab and cannot be added now.')), call. = FALSE)
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
            if(length(LHS$nDim) == 0) stop(exprClassProcessingErrorMsg(code, paste0('In sizeAssignAfterRecursing: nDim for LHS not set.')), call. = FALSE)
            if(length(RHSnDim) == 0) stop(exprClassProcessingErrorMsg(code, paste0('In sizeAssignAfterRecursing: nDim for RHS not set.')), call. = FALSE)
            if(LHS$nDim != RHSnDim) {
                message(paste0('Warning, mismatched dimensions in assignment: ', nimDeparse(code), '. Going to browser(). Press Q to exit'), call. = FALSE )
                browser()
            }
            ## and warn if type issue e.g. int <- double
            if(assignmentTypeWarn(LHS$type, RHStype)) {
                message(paste0('Warning, RHS numeric type is losing information in assignment to LHS.', nimDeparse(code)))
            }
        }
    }
    ## update size info in typeEnv
    assert <- NULL

    if(LHS$name == 'values' && length(LHS$args) %in% c(1,2)) { ## It is values(model_values_accessor[, index]) <- STUFF
        # triggered when we have simple assignment into values() without indexing of values()
        if(is.numeric(RHS)) stop(exprClassProcessingErrorMsg(code, paste0('In sizeAssignAfterRecursing: Cannot assign into values() from numeric.')), call. = FALSE)
        code$name <- 'setValues' 
        code$args <- list(1 + length(LHS$args))
        setArg(code, 1, RHS)
        setArg(code, 2, LHS$args[[1]])
        if(length(LHS$args) == 2) setArg(code, 3, LHS$args[[2]])  # for indexed of form values(model, nodes[i])
        if(!(RHS$isName)) assert <- c(assert, sizeInsertIntermediate(code, 1, symTab, typeEnv) )
        return( if(length(assert) == 0) NULL else assert )
    }

    ## Note this can use LHS$name for RHSsizeExprs when returning from a nimbleFunction on RHS.  But this is probably not needed any more.
    if(any(unlist(lapply(RHSsizeExprs, is.null)))) RHSsizeExprs <- makeSizeExpressions(rep(NA, RHSnDim), LHS$name) ## reset sizeExprs for the LHS var. re-using RHSsizeExprs for LHS.  This would only be valid if it is a nimbleFunction returning something on the RHS.  For assignment to be executed in Eigen, the RHS sizes MUST be known

    ## We used to update typeEnv sizeExprs, but in some cases it is not safe to do so
    ## Hence they are created generically above if the LHS$name is new
    ## typeEnv[[LHS$name]]$sizeExprs <- RHSsizeExprs


    if(LHS$toEigenize == 'yes') {
        code$toEigenize <- 'yes'
##        message('Warning from sizeAssign: not expecting LHS to have toEigenize == yes')
    } else {
        code$toEigenize <-if(inherits(RHS, 'exprClass')) {
            if(RHS$toEigenize == 'no') 'no' else {
                if(RHS$toEigenize == 'unknown') 'no' else {
                  # if(RHS$toEigenize != 'yes' & (RHS$nDim == 0 | (RHS$isName & LHS$name != "nfVar") | (RHS$name == 'map' & NoEigenizeMap))) 'no' ## if it is scalar or is just a name or a map, we will do it via NimArr operator= .  Used to have "| RHS$name == 'map'", but this allowed X[1:3] <- X[2:4], which requires eigen, with eval triggered, to get right
                    if(RHS$toEigenize != 'yes' & (!(LHS$name %in% c('eigenBlock', 'diagonal', 'coeffSetter'))) & (RHS$nDim == 0 | RHS$isName | (RHS$name == 'map' & NoEigenizeMap))) 'no' ## if it is scalar or is just a name or a map, we will do it via NimArr operator= .  Used to have "| RHS$name == 'map'", but this allowed X[1:3] <- X[2:4], which requires eigen, with eval triggered, to get right
                    else 'yes' ## if it is 'maybe' and non-scalar and not just a name, default to 'yes'
                }
            }
        } else {
            if(is.numeric(LHS$nDim))
                if(LHS$nDim > 0) 'yes' ## This is for cases like out[1:4] <- scalar
                else 'no'
            else 'no'
        }
    }

    if(code$toEigenize == 'yes') { ## this would make more sense in eigenize_assign
    ## generate setSize(LHS, ...) where ... are dimension expressions
        if(length(RHSnDim) == 0) {
            message("confused about trying to eigenize something with nDim = 0")
            browser()
        }
        if(RHSnDim > 0) {
            # if(TRUE) { ## !identical(LHSdrop$sizeExprs, RHSdrop$sizeExprs)) {## This was too clever: it was to prevent redundant calls to setSize, but the problem is the previous call could have been generated inside an if-then-else, so we can't rely on it
            #     if(LHS$isName | LHS$name == "nfVar") {
            #         assert <- list(substitute(setSize(LHS), list(LHS = parse(text = nimDeparse(LHS), keep.source = FALSE)[[1]])))
            if(!(RHS$name %in% setSizeNotNeededOperators)) {
                if(LHS$isName | LHS$name == "nfVar") {
                    assert <- substitute(setSize(LHS), list(LHS = nimbleGeneralParseDeparse(LHS)))
                ##if(LHS$isName) {
                ##    assert <- list(substitute(setSize(LHS), list(LHS = as.name(LHS$name))))
                    for(i in seq_along(RHSsizeExprs)) {
                        test <- try(assert[[i + 2]] <- RHS$sizeExprs[[i]])
                        if(inherits(test, 'try-error')) browser()
                    }
                    assert[[ length(assert) + 1]] <- 0 ## copyValues = false
                    assert[[ length(assert) + 1]] <- 0 ## fillZeros  = false
                    assert <- list(assert)
                } else { ## We have an indexed LHS of an eigenizable expression
                    ## need special handling if it is a row assignment like x[i,] <- ...
                    ## also need to generate size assertions                    
                    if(LHS$nDim == 1) {
                        if(RHS$nDim == 2) {
                            if(is.numeric(RHS$sizeExprs[[1]])) {
                                if(RHS$sizeExprs[[1]] == 1) {
                                    newExpr <- insertExprClassLayer(code, 1, 'asRow', type = LHS$type)
                                    newExpr$sizeExprs <- RHS$sizeExprs 
                                    newExpr$type <- LHS$type
                                    newExpr$nDim <- RHS$nDim
                                    if(!is.numeric(LHS$sizeExprs[[1]]) | !is.numeric(RHS$sizeExprs[[2]])) {
                                        assertMessage <- paste0("Run-time size error: expected ", deparse(LHS$sizeExprs[[1]]), " == ", deparse(RHS$sizeExprs[[2]]))
                                        thisAssert <- identityAssert(LHS$sizeExprs[[1]], RHS$sizeExprs[[2]], assertMessage)
                                        if(!is.null(thisAssert)) assert[[length(assert) + 1]] <- thisAssert 
                                    } else {
                                        if(LHS$sizeExprs[[1]] != RHS$sizeExprs[[2]]) stop(exprClassProcessingErrorMsg(code, paste0('In sizeAssignAfterRecursing: Fixed size mismatch.')), call. = FALSE)
                                    }
                                }
                            }
                        }
                    }
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

    if(!(LHS$name %in% c('eigenBlock', 'diagonal', 'coeffSetter', 'nimNonseqIndexedd', 'nimNonseqIndexedi','nimNonseqIndexedb'))) {
        ## should already be annotated if it is an indexed assignment.
        ## It should be harmless to re-annotated EXCEPT in case like out[1:5] <- scalar
        code$nDim <- code$args[[1]]$nDim <- RHSnDim
        code$type <- code$args[[1]]$type <- RHStype
        code$sizeExprs <- code$args[[1]]$sizeExprs <- RHSsizeExprs
    }
    if(RHSname %in% assignmentAsFirstArgFuns) {
        code$name <- RHS$name
        oldArgs <- RHS$args
        LHS <- code$args[[1]] ## could have been reset by LHS$name == 'map' situation above
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

sizeScalarModelOp <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) {
        asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, rep(TRUE, length(code$args)-1)))
        ## This used to error-trap attempts at index expressions.
        ## Now they are allowed
        ## I think length(code$args) should only ever be 1 or 2 but will write more generally
        for(i in 2:length(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
##                if(code$args[[i]]$toEigenize=='yes') stop(exprClassProcessingErrorMsg(code, 'In sizeScalarModelOp: There is an expression beyond the first argument that cannot be handled.  If it involve vectorized math, you need to do it separately, not in this expression.'), call. = FALSE)
                ## Now instead we see if there is any eigenization we will lift it
                ## The C++ code is now flexible enough that we shouldn't have to lift eigen expressions
                ## The limitation at this moment is our eigenization processing, which isn't set up to eigenize just part of a line.
                if(code$args[[i]]$toEigenize=='yes')
                    asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv))
            }
        }
        if(inherits(code$args[[2]], 'exprClass')) { ## There is an index expression that may be non-scalar
            if(code$args[[2]]$nDim > 0) { ## It is non-scalar so we need to set a logical argument about whether is it a logical or numeric vector
                code$args[[ length(code$args)+1 ]] <- as.integer(code$args[[2]]$type == 'logical')
            }
        }
    } else {
        asserts <- list()
    }
    if(code$args[[1]]$toEigenize == 'yes') { ## not sure when this would be TRUE
        asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
    }
    code$nDim <- 0
    outputType <- scalarOutputTypes[[code$name]]
    if(is.null(outputType)) code$type <- 'double'
    else code$type <- outputType
    code$sizeExprs <- list()
    code$toEigenize <- 'maybe'
    asserts
}

sizeSimulate <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) {
        asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, rep(TRUE, length(code$args)-1)))
        for(i in 2:length(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
                ## changes similar to sizeScalarModelOp to allow general vectorized indexing
                ##       if(code$args[[i]]$toEigenize=='yes') stop(exprClassProcessingErrorMsg(code, 'In sizeSimulate: There is an expression beyond the first argument that cannot be handled.  If it involve vectorized math, you need to do it separately, not in this expression.'), call. = FALSE)
                if(code$args[[i]]$toEigenize=='yes')
                    asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv))##toEigenize <- 'yes'

            }
        }
        if(inherits(code$args[[2]], 'exprClass')) { ## There is an index expression that may be non-scalar
            if(code$args[[2]]$nDim > 0) { ## It is non-scalar so we need to set a logical argument about whether is it a logical or numeric vector
                code$args[[ length(code$args)+1 ]] <- as.integer(code$args[[2]]$type == 'logical')
            }
        }
    } else {
        asserts <- list()
    }

    code$nDim <- 0
    code$type <- as.character(NA)
    code$sizeExprs <- list()
    code$toEigenize <- 'maybe'
    return(asserts)
}

sizeScalarRecurse <- function(code, symTab, typeEnv, recurse = TRUE) {
    ## use something different for distributionFuns
    asserts <- if(recurse) recurseSetSizes(code, symTab, typeEnv) else list()
    ## This just forces any argument expression to be lifted.  Can we lift only things to be eigenized?
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
    if(length(asserts)==0) NULL else asserts
}

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
    stop(exprClassProcessingErrorMsg(code, paste0('In sizeBinaryUnarycWise: Length of arguments is not 1 or 2.')), call. = FALSE)
}

sizemvAccessBracket <- function(code, symTab, typeEnv) {
    ## this gets called from sizeIndexingBracket, so recurse has already been done
    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, TRUE))
    if(length(code$args) != 2) {
        stop(exprClassProcessingErrorMsg(code, paste0('In sizemvAccessBracket:  Wrong number of indices provided.')), call. = FALSE)
    }
    if(inherits(code$args[[2]], 'exprClass')) {
        if(code$args[[2]]$nDim != 0) stop(exprClassProcessingErrorMsg(code, paste0('In sizemvAccessBracket:  Index is not a scalar.')), call. = FALSE)
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

    for(i in seq_along(code$args)) {
        if(i == 1) next
        if(inherits(code$args[[i]], 'exprClass')) {
            if(code$args[[i]]$name != "")
                if(code$args[[i]]$type == 'logical') {
                    ## first insert which, then lift to intermediate
                    newExpr <- insertExprClassLayer(code, i, 'which')
                    useBool <- rep(FALSE, length(code$args))
                    useBool[i] <- TRUE
                    asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, useBool))
                    ## sizeWhich will lift it to an intermediate and annotate it 
##                    asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv))
                    
                }
        }
    }
    
    nDimVar <- code$args[[1]]$nDim

    dropBool <- TRUE
    dropArgProvided <- FALSE
    if(!is.null(names(code$args)))
        if('drop' %in% names(code$args)) {
            dropArgProvided <- TRUE
            iDropArg <- which(names(code$args) == 'drop')
        }
    if(nDimVar != length(code$args) - 1 - dropArgProvided) {
        ## only valid case with fewer index arguments than source dimensions is matrix[indices]
        if(!( (nDimVar == 2) & (length(code$args) - dropArgProvided) == 1)) {
            msg <- paste0('Error, wrong number of indices provided for ', nimDeparse(code),'.')
            stop(exprClassProcessingErrorMsg(code, msg), call. = FALSE)
        }
    }

    if(dropArgProvided) {
        dropBool <- code$args[[iDropArg]]
        if(!is.logical(dropBool)) {
            msg <- paste0(msg, "(A drop argument must be hard-coded as TRUE or FALSE, not given as a variable.)")
            stop(exprClassProcessingErrorMsg(code, msg), call. = FALSE)
        }
    }
    code$nDim <- nDimVar
    code$type <- code$args[[1]]$type
    code$sizeExprs <- vector('list', length = nDimVar)
    ## We could generate asserts here to ensure sub-indexing is within bounds
    needMap <- FALSE
    ## Track whether if all index ranges are defined by `:` or by scalar
    simpleBlockOK <- TRUE
    iSizes <- 1
    for(i in 1:nDimVar) {
        dropThisDim <- FALSE

        if(is.numeric(code$args[[i+1]])) dropThisDim <- TRUE
        else if((code$args[[i+1]]$name != "") & (length(dropSingleSizes(code$args[[i+1]]$sizeExprs)$sizeExprs) == 0)) dropThisDim <- TRUE

        isExprClass <- inherits(code$args[[i+1]], 'exprClass') ## code$args[[1]] ???
        
        if(dropThisDim) { ## The index is a scalar
            if(nimbleOptions()$indexDrop & dropBool) {
                code$sizeExprs[[iSizes]] <- NULL
                code$nDim <- code$nDim - 1
            } else { 
                code$sizeExprs[[iSizes]] <- 1; iSizes <- iSizes + 1
            }
            next
        } else {
            if(isExprClass) ## Need to think through cases here more
                if((code$args[[i+1]]$name != ':') && (code$args[[i+1]]$name != "")) simpleBlockOK <- FALSE
        }
        needMap <- TRUE ## If the "next" in if(dropThisDim) {} is always hit, then needMap will never be set to TRUE
        if(isExprClass) {
            if(code$args[[i+1]]$name != "") {
                ## An entry that is a variable possibly with a length
                code$sizeExprs[[iSizes]] <- code$args[[i+1]]$sizeExprs[[1]]
            } else {
                ## blank entry (e.g. A[,i]) is an exprClass with isName = TRUE and name = ""
                code$sizeExprs[[iSizes]] <- code$args[[1]]$sizeExprs[[i]]
                ## also at this point we will fill in a `:` expression for the indices
                newIndexExpr <- RparseTree2ExprClasses(substitute(1:N, list(N = code$args[[1]]$sizeExprs[[i]])))
                setArg(code, i+1, newIndexExpr)
                useArgs <- rep(FALSE, length(code$args))
                useArgs[i+1] <- TRUE
                asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, useArgs))
            }
            iSizes <- iSizes + 1
            next
        }
    }
    ## did all dims get dropped?
    if(length(code$sizeExprs)==0) {
        code$sizeExprs <- list() ## it was a named, list.  this creates consistency. maybe unnecessary
        ##needMap will be FALSE if we are in this clause
        if(!code$args[[1]]$isName)
            if(!(code$args[[1]]$name %in% operatorsAllowedBeforeIndexBracketsWithoutLifting)) {## 'mvAccessRow'){
                ## At this point we have decided to lift, and the next two if()s determine if that is weird due to being on LHS of assignment
                if(code$caller$name %in% assignmentOperators)
                    if(code$callerArgID == 1)
                        stop(exprClassProcessingErrorMsg(code, 'There is a problem on the left-hand side of an assignment'), call. = FALSE)
                
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            }
    }
    
    code$toEigenize <- 'maybe'
    if(needMap) {
        ## If this is a map on an *expression* that is not a map, we used to always lift it
        ## e.g. (A + B)[1:4] must become (Interm <- A + B; Interm[1:4])
        ## Now we only need to lift it if the map will not be impemented via eigenBlock

        ## Or maybe we never need to
        ## liftArg <- FALSE
        ## mappingViaEigenBlock <- simpleBlockOK && code$args[[1]]$nDim <= 2
        ## if(!code$args[[1]]$isName)                   ## In X[I], X is an expression, not just a name
        ##     if(!mappingViaEigenBlock) {              ## Handling of `[` will NOT be via .block() in Eigen
        ##         if(code$args[[1]]$name != 'map')     ## Lift unless it is already a map, which can be compounded
        ##             liftArg <- TRUE
        ##     }
        ## if(liftArg) asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))

        ## for nested blocking, we have (nonseq | eigenBlock | map) x (nonseq | eigenBlock | map)
        ##      [ coeffSetter is a version of nonseq ]
        ## (map) x (eigenBlock | map) is already handled
        ## (map) x (nonseq)           is already handled
        ## (eigenBlock) x (eigenBlock) is already handled
        ## 
        ## check whether to nest the indexing directly


        nestIndexing <- FALSE
        if(!code$args[[1]]$isName) { ## In X[i], X is an expression
            ## X is an indexing expression of some kind (other than a map, which is already a new object)
            ## It can't be coeffSetter at this point in processing flow, because the nestedness implies its caller was not <-
            if(code$args[[1]]$name %in% c('eigenBlock', 'nimNonseqIndexedd' ,'nimNonseqIndexedi' ,'nimNonseqIndexedb' )) {
                ## if it is not (eigenBlock) x (eigenBlock)
                if(!( (code$args[[1]]$name == 'eigenBlock') & (simpleBlockOK)))
                    nestIndexing <- TRUE
            }
        }
        
        if(nestIndexing) {
            if(code$args[[1]]$name == 'eigenBlock') {
                ## reach down to X and rename it
                ## put `:`(start, finish) back together.
                if(code$caller$name %in% assignmentOperators & code$callerArgID == 1) {
                    code$args[[1]]$name <- 'coeffSetter'
                } else {
                    if(code$type == 'double') code$args[[1]]$name <- 'nimNonseqIndexedd'
                    if(code$type == 'integer') code$args[[1]]$name <- 'nimNonseqIndexedi'
                    if(code$type == 'logical') code$args[[1]]$name <- 'nimNonseqIndexedb'
                }
            } else {
                ## it was already nonseq, but it might need to become coeffSetter
                if(code$caller$name %in% assignmentOperators & code$callerArgID == 1) {
                    code$args[[1]]$name <- 'coeffSetter'
                }
            }
            nestedNinds <- length(code$args[[1]]$args)-1
            nestedNdim <- code$args[[1]]$nDim
            nestedDropBool <- TRUE
            nestedDropArgProvided <- FALSE
            if(!is.null(names(code$args[[1]]$args)))
                if("drop" %in% names(code$args[[1]]$args)) {
                    nestedDropArgProvided <- TRUE
                    nestedDropBool <- code$args[[1]]$args[[ which(names(code$args[[1]]$args) == 'drop') ]]
                    nestedNinds <- nestedNinds - 1
                }
            nestedBlockBool <- rep(TRUE, nestedNinds)    ## is it preserved as a block (can still be scalar if nestedDropBool is FALSE)
            nestedScalarIndex <- rep(FALSE, nestedNinds) 
            
            for(iInd in 1:nestedNinds) {
                if(is(code$args[[1]]$args[[iInd+1]], 'exprClass'))  {
                    if(code$args[[1]]$args[[iInd+1]]$nDim == 0) {
                        nestedScalarIndex[iInd] <- TRUE
                        if(nestedDropBool) nestedBlockBool[iInd] <- FALSE
                    }
                } else {
                    nestedScalarIndex[iInd] <- TRUE
                    if(nestedDropBool) nestedBlockBool[iInd] <- FALSE
                }
            }
            
            code$args[[1]]$sizeExprs <- code$sizeExprs
            code$args[[1]]$nDim <- code$nDim
            code$args[[1]]$type <- code$type
            numIndices <- length(code$args) - 1 - dropArgProvided
            
            ## NEED TO SKIP SCALARS IF dropBool = TRUE for nested case.
            nestedInds <- which(nestedBlockBool)
            if(length(nestedInds) != numIndices) stop(exprClassProcessingErrorMsg(code, 'Wrong number of nested indices.'), call.=FALSE)
            for(iInd in 1:numIndices) {
                nestedIind <- nestedInds[iInd]
                nestedIndexIsScalar <- if(inherits(code$args[[1]]$args[[nestedIind + 1]], 'exprClass')) code$args[[1]]$args[[nestedIind + 1]]$nDim == 0 else TRUE
                if(nestedIndexIsScalar) {
                    indexIsScalar <- if(inherits(code$args[[iInd+1]], 'exprClass')) code$args[[iInd+1]]$nDim == 0 else TRUE
                    if(!indexIsScalar) warning("There is nested indexing with drop=FALSE where an index must be scalar but isn't")
                } else {
                    newExpr <- nimble:::exprClass(name = 'nimNonseqIndexedi', isName = FALSE, isCall = TRUE, isAssign = FALSE)
                    newExpr$type <- 'integer'
                    indexIsScalar <- if(inherits(code$args[[iInd+1]], 'exprClass')) code$args[[iInd+1]]$nDim == 0 else TRUE
                    newExpr$sizeExprs <- if(!indexIsScalar) c(code$args[[iInd + 1]]$sizeExprs) else list(1)
                    newExpr$nDim <- 1
                    newExpr$toEigenize <- 'yes'
                    ## sizeExprs, nDim, toEigenize?
                    setArg(newExpr, 1, code$args[[1]]$args[[nestedIind + 1]])
                    setArg(newExpr, 2, code$args[[iInd + 1]])
                    setArg(newExpr, 3, 1)
                    setArg(code$args[[1]], nestedIind + 1, newExpr)
                }
            }
            code$args[1+(1:numIndices)] <- NULL
            codeCaller <- code$caller
            codeCallerArgID <- code$callerArgID
            removeExprClassLayer(code)
            code <- codeCaller$args[[codeCallerArgID]]
            return(if(length(asserts)==0) NULL else asserts)
        }
        
        ## Replace with a map expression if needed
        if(!simpleBlockOK) {
            if(typeEnv$.ensureNimbleBlocks) {
                stop(exprClassProcessingErrorMsg(code, "LHS indexing for a multivariate random draw can only use sequential blocks (via ':')."), call. = FALSE)
            }
            ##   if(nDimVar != length(code$args) - 1) code$args[[length(code$args)]] <- NULL 
            if(code$caller$name %in% assignmentOperators & code$callerArgID == 1) {
                code$name <- 'coeffSetter'
            } else {
                if(code$type == 'double') code$name <- 'nimNonseqIndexedd' ## this change could get moved to genCpp_generateCpp 
                if(code$type == 'integer') code$name <- 'nimNonseqIndexedi'
                if(code$type == 'logical') code$name <- 'nimNonseqIndexedb'
            }
            if(length(code$args) - 1 - dropArgProvided == 1) ## only 1 index
                code$args[[3]] <- 1 ## fill in extra 1 for a second dimension
        }
        else {
            if(code$args[[1]]$nDim > 2 | typeEnv$.ensureNimbleBlocks) { ## old-style blocking from >2D down to 2D or 1D, or this is LHS for something like rmnorm, requiring a non-eigen map on LHS.
                if(dropArgProvided) code$args[[iDropArg]] <- NULL
                newExpr <- makeMapExprFromBrackets(code, dropBool)
                newExpr$sizeExprs <- code$sizeExprs
                newExpr$type <- code$type
                newExpr$nDim <- code$nDim
                newExpr$toEigenize <- code$toEigenize
                setArg(code$caller, code$callerArgID, newExpr)
            }
            else { ## blocking via Eigen
                ## newExpr <- makeEigenBlockExprFromBrackets(code, dropBool) ## at this point it is ok that code exprClass is messed up (first arg re-used in newExpr)
                ## newExpr$sizeExprs <- code$sizeExprs
                ## newExpr$type <- code$type
                ## newExpr$nDim <- code$nDim
                ## newExpr$toEigenize <- 'yes'
                ## ## note that any expressions like sum(A) in 1:sum(A) should have already been lifted
                code$name <- 'eigenBlock'
                code$toEigenize <- 'yes'
            }
            ##setArg(code$caller, code$callerArgID, newExpr)
        }
    }
    if(length(asserts)==0) NULL else asserts
}

isIntegerEquivalent <- function(code) {
    if(inherits(code, 'exprClass')) {
        if(code$type == 'integer') return(TRUE)
        return(FALSE)
    }
    if(is.logical(code)) return(FALSE)
    if(storage.mode(code) == 'integer') return(TRUE)
    code == floor(code) ## storage.mode must be 'double' so check if it's equivalent to an integer 
}

sizeSeq <- function(code, symTab, typeEnv, recurse = TRUE) {
    message('still need to handle -1L sequences')
    message('check that arguments are scalars')
    asserts <- if(recurse) recurseSetSizes(code, symTab, typeEnv) else list()
    byProvided <- code$name == 'nimSeqBy' | code$name == 'nimSeqByLen'
    lengthProvided <- code$name == 'nimSeqLen' | code$name == 'nimSeqByLen'
    integerFrom <- isIntegerEquivalent(code$args[[1]])
    integerTo <- isIntegerEquivalent(code$args[[2]])
    liftExprRanges <- TRUE
    if(integerFrom && integerTo) {
        if((!byProvided && !lengthProvided) || (byProvided && !lengthProvided && is.numeric(code$args[[3]]) && code$args[[3]] == 1)) {
            code$name = ':'
            asserts <- c(asserts, sizeColonOperator(code, symTab, typeEnv, recurse = FALSE))
            return(if(length(asserts)==0) NULL else asserts)
        }
    } else {
        if(!byProvided && !lengthProvided) {
            code$args[[3]] <- 1
            byProvided <- TRUE
        }
        if(byProvided) {
            code$name <- 'nimSeqByD'
            ## lift any expression arguments
            for(i in 1:2) {
                if(inherits(code$args[[i]], 'exprClass')) {
                    if(!code$args[[i]]$isName) {
                        asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv) )
                    }
                }
            }
            if(lengthProvided) {
                code$name <- 'nimSeqByLenD'
                thisSizeExpr <- parse(text = nimDeparse(code$args[[4]]), keep.source = FALSE)[[1]]
            } else {
                thisSizeExpr <- substitute(1 + floor((TO_ - FROM_) / BY_),
                                           list(FROM_ = parse(text = nimDeparse(code$args[[1]]), keep.source = FALSE)[[1]],
                                                TO_ = parse(text = nimDeparse(code$args[[2]]), keep.source = FALSE)[[1]],
                                                BY_ = parse(text = nimDeparse(code$args[[3]]), keep.source = FALSE)[[1]]))
            }
        } else { ## must be lengthProvided
            code$name <- 'nimSeqLenD'
            thisSizeExpr <- parse(text = nimDeparse(code$args[[4]]), keep.source = FALSE)[[1]]
        }
    }
    code$type <- 'double' ## only remaining case to catch here is -1 integer sequences, which we don't move to `:`
    code$sizeExprs <- list(thisSizeExpr)
    code$toEigenize <- 'yes'
    code$nDim <- 1
    return(if(length(asserts)==0) NULL else asserts)
}


sizeColonOperator <- function(code, symTab, typeEnv, recurse = TRUE) {
    asserts <- if(recurse) recurseSetSizes(code, symTab, typeEnv) else list()
    if(length(code$args) != 2) stop(exprClassProcessingErrorMsg(code, 'In sizeColonOperator: Problem determining size for : without two arguments.'), call. = FALSE)

    ## Had the idea that some cases be shifted to seq(), but I don't think that makes sense
    ## moveToSeqBy1 <- FALSE
    ## if(inherits(code$args[[1]], 'exprClass')) {
    ##     if(code$args[[1]]$nDim != 0) message('WARNING, first arg to : is not scalar in expression', nimDeparse(code))
    ##     if(code$args[[1]]$type == 'double') moveToSeqBy1 <- TRUE
    ## } else {
    ##     if(storage.mode(code$args[[1]]) == 'double') if(code$args[[1]] != floor(code$args[[1]])) moveToSeqBy1 <- TRUE
    ## }
    ## if(inherits(code$args[[2]], 'exprClass')) {
    ##     if(code$args[[2]]$nDim != 0) message('WARNING, second arg to : is not scalar in expression', nimDeparse(code))
    ##     if(code$args[[2]]$type == 'double') moveToSeqBy1 <- TRUE
    ## } else {
    ##     if(storage.mode(code$args[[2]]) == 'double') if(code$args[[2]] != floor(code$args[[2]])) moveToSeqBy1 <- TRUE
    ## }

    ## if(moveToSeqBy1) {
    ##     code$name <- 'nimSeq'
    ##     asserts <- c(asserts, sizeSeq(code, symTab, typeEnv, recurse = FALSE))
    ##     return(if(length(asserts)==0) NULL else asserts)
    ## }
    for(i in 1:2) {
        if(inherits(code$args[[i]], 'exprClass')) {
            if(!code$args[[i]]$isName) {
              if(! (code$args[[i]]$name == '[' && (code$args[[i]]$args[[1]]$name == 'dim' && code$args[[i]]$args[[1]]$args[[1]]$name == 'nfVar'))){
                asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv) )
              }
            }
        }
    }
    
    code$type <- 'double'
    code$nDim <- 1
    code$toEigenize <- 'maybe' 
    
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

## sizeDimOperator <- function(code, symTab, typeEnv) {
##     ## a special case since the resulte is stored as a vector<int>
##     ## recurse and set Intermediate is if it not a name or map
##     ## we'll want a NimArr<1, int> constructor vector<int> 
##     recurseSetSizes(code, symTab, typeEnv)
##     code$type <- 'integer'
##     code$nDim <- 1
##     code$sizeExprs <- list( code$args[[1]]$nDim )
##     code$toEigenize <- 'no'
##     invisible(NULL)
## }

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
        code$name <- 'asRow'
        code$sizeExprs <- c(list(1), code$sizeExprs[[1]])
        code$nDim <- 2
    }
    return(ans)
}

getArgumentType <- function(expr) {
    if(inherits(expr, 'exprClass')) {
        expr$type
    } else
        storage.mode(expr)
}

setReturnType <- function(keyword, argType) {
    handling <- returnTypeHandling[[keyword]]
    if(is.null(handling)) return('double')
    switch(handling,
           'double', ##1
           'integer', ##2
           'logical', ##3
           argType, ##4
           if(argType == 'logical') 'integer' else argType ##5
           )
}

## Handler for unary functions that operate component-wise
sizeUnaryCwise <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1){
        stop(exprClassProcessingErrorMsg(code, 'sizeUnaryCwise called with argument length != 1.'), call. = FALSE)
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
  ##      code$type <- a1$type
    } else {
        code$nDim <- 0
        code$sizeExprs <- list()
##        code$type <- 'double'
    }
    code$type <- setReturnType(code$name, getArgumentType(a1))
    if(length(code$nDim) != 1) stop(exprClassProcessingErrorMsg(code, 'In sizeUnaryCwise: nDim is not set.'), call. = FALSE)
    code$toEigenize <- if(code$nDim > 0) 'yes' else 'maybe'
    return(asserts)
}

## currently only inprod(v1, v2)
sizeBinaryReduction <- function(code, symTab, typeEnv) {
    if(length(code$args) != 2) stop(exprClassProcessingErrorMsg(code, 'In sizeBinaryReduction: argument length != 2'), call. = FALSE)
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
    if(!ok) stop(exprClassProcessingErrorMsg(code, 'Cannot call inprod or other binary reduction operator with constant argument.'), call. = FALSE)
    
    code$nDim <- 0
    code$sizeExprs <- list()
    code$type <- 'double'
    code$toEigenize <- 'yes'
    if(length(asserts) == 0) NULL else asserts
}


## things like trace, det, logdet
sizeMatrixSquareReduction <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1){
        stop(exprClassProcessingErrorMsg(code, 'sizeMatrixSquareReduction called with argument length != 1.'), call. = FALSE)
    }
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    if(!inherits(a1, 'exprClass')) stop(exprClassProcessingErrorMsg(code, 'sizeMatrixSquareReduction called with argument that is not an expression.'), call. = FALSE)
    if(a1$toEigenize == 'no') {
        asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
        a1 <- code$args[[1]]
    }
    if(a1$nDim != 2) stop(exprClassProcessingErrorMsg(code, 'sizeMatrixSquareReduction called with argument that is not a matrix.'), call. = FALSE)
    code$nDim <- 0
    code$sizeExprs <- list()
    code$type <- if(code$name == 'trace') code$args[[1]]$type else 'double'
    code$toEigenize <- 'yes'

    if(!(code$caller$name %in% c('{','<-','<<-','='))) {
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }
    if(length(asserts) == 0) NULL else asserts
}


#for eigen() and svd() function
sizeMatrixEigenList <- function(code, symTab, typeEnv){
  if(code$name == 'EIGEN_EIGEN'){
    if(length(code$args) != 2){
      stop(exprClassProcessingErrorMsg(code, 'eigen() called with inappropriate argument length.'), call. = FALSE)
    }
  }
  if(code$name == 'EIGEN_SVD'){
    if(length(code$args) != 2){
      stop(exprClassProcessingErrorMsg(code, 'svd() called with inappropriate argument length.'), call. = FALSE)
    }
  }
  asserts <- recurseSetSizes(code, symTab, typeEnv)
  a1 <- code$args[[1]]
  
  if(!inherits(a1, 'exprClass')) stop(exprClassProcessingErrorMsg(code, 'sizeMatrixEigenList called with argument that is not an expression.'), call. = FALSE)
  if(a1$nDim != 2) stop(exprClassProcessingErrorMsg(code, 'sizeMatrixEigenList called with argument that is not a matrix.'), call. = FALSE)

  code$type <- 'symbolNimbleList'
  listST <- symTab$getParentST()$getSymbolObject(paste0(code$name, 'CLASS'))
  code$sizeExprs <- listST
  code$toEigenize <- "yes"
  code$nDim <- 0
  # if(!(code$caller$name %in% c('{','<-','<<-','='))) {
  #   asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
  # }
  if(length(asserts) == 0) NULL else asserts
}

sizeUnaryCwiseSquare <- function(code, symTab, typeEnv) {
    if(length(code$args) != 1){
    	stop(exprClassProcessingErrorMsg(code, 'sizeUnaryCwiseSquare called with argument length != 1.'), call. = FALSE)
    }
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    if(!inherits(a1, 'exprClass')) stop(exprClassProcessingErrorMsg(code, 'sizeUnaryCwiseSquare called with argument that is not an expression.'), call. = FALSE)
    if(a1$toEigenize == 'no') {
        asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
        a1 <- code$args[[1]]
    }
    if(a1$nDim != 2) stop(exprClassProcessingErrorMsg(code, 'sizeUnaryCwiseSquare called with argument that is not a matrix.'), call. = FALSE)
    if(!identical(a1$sizeExprs[[1]], a1$sizeExprs[[2]])) {
        asserts <- c(asserts, identityAssert(a1$sizeExprs[[1]], a1$sizeExprs[[2]], paste0("Run-time size error: expected ", nimDeparse(a1), " to be square.") ))
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
        newSize <- a1$sizeExprs[[1]]
    }
    code$nDim <- 2
    code$sizeExprs <- list(newSize, newSize)
    code$type <- setReturnType(code$name, a1$type)
    code$toEigenize <- if(code$nDim > 0) 'yes' else 'maybe'
    invisible(asserts)
}

sizeUnaryNonaryCwise <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) stop(exprClassProcessingErrorMsg(code, 'sizeUnaryNonaryCwise called with argument length > 1'), call. = FALSE) 
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
    if(length(code$args) != 1) stop(exprClassProcessingErrorMsg(code, 'sizeUnaryReduction called with argument length != 1.'), call. = FALSE)
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(inherits(code$args[[1]], 'exprClass')) {
        ## Kludgy catch of var case here.  Can't do var(matrix) because in R that is interpreted as cov(data.frame)
        if(code$args[[1]]$nDim >= 2) {
            if(code$name == 'var') {
                stop(exprClassProcessingErrorMsg(code, 'NIMBLE compiler does not support var with a matrix (or higher dimensional) argument.'), call. = FALSE) 
            }
        }
        if(!code$args[[1]]$isName) {
            if(code$args[[1]]$toEigenize == 'no') {
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            }
        }
    }

    code$nDim <- 0
    code$sizeExprs <- list()
    code$type <- setReturnType(code$name, code$args[[1]]$type)
    code$toEigenize <- 'yes'

    if(!(code$caller$name %in% c('{','<-','<<-','='))) {
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }

    if(length(asserts) == 0) NULL else asserts
}

## There's no real point in annotating return.  Just need to recurse and lift 
sizeReturn <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) stop(exprClassProcessingErrorMsg(code, 'return has argument length > 1.'), call. = FALSE)
    code$toEigenize <- 'no'
    if(length(code$args) == 0) return(invisible(NULL))
    
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(inherits(code$args[[1]], 'exprClass')) {
        if(!code$args[[1]]$isName) {
##            if(code$args[[1]]$toEigenize == 'yes') {
            if(anyNonScalar(code$args[[1]])) {
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv, forceAssign = TRUE))
            }
            ##            }
        }
    }
    invisible(asserts)
}

sizeMatrixMult <- function(code, symTab, typeEnv) {
    if(length(code$args) != 2) stop(exprClassProcessingErrorMsg(code, 'sizeMatrixMult called with argument length != 2.'), call. = FALSE)
    a1 <- code$args[[1]]
    a2 <- code$args[[2]]

    if(!(inherits(a1, 'exprClass') & inherits(a2, 'exprClass'))) stop(exprClassProcessingErrorMsg(code, 'In sizeMatrixMult: expecting both arguments to be expressions.'), call. = FALSE)

    ## need promotion from vectors to matrices with asRow or asCol
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    a1 <- code$args[[1]]
    a2 <- code$args[[2]]
    
    if(a1$nDim == 0 | a2$nDim == 0) stop(exprClassProcessingErrorMsg(code, 'In sizeMatrixMult: Cannot do matrix multiplication with a scalar.'), call. = FALSE) 

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
        origSizeExprs <- a1$sizeExprs[[1]]
        a1 <- insertExprClassLayer(code, 1, 'asRow', type = a1$type)
        a1$sizeExprs <- c(list(1), origSizeExprs)
        origSizeExprs <- a2$sizeExprs[[1]]
        a2 <- insertExprClassLayer(code, 2, 'asCol', type = a2$type)
        a2$sizeExprs <- c(origSizeExprs, list(1))
    } else {
        if(a1$nDim == 1) {
            if(a2$nDim != 2) stop(exprClassProcessingErrorMsg(code, paste0('In sizeMatrixMult: First arg has nDim = 1 and 2nd arg has nDim = ', a2$nDim, '.')), call. = FALSE)
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
            if(a1$nDim != 2) stop(exprClassProcessingErrorMsg(code, paste0('In sizeMatrixMult: Second arg has nDim = 1 and 1st arg has nDim = ', a1$nDim, '.')), call. = FALSE)
            if(identical(a1$sizeExprs[[2]], 1)) {
                a2 <- insertExprClassLayer(code, 2, 'asRow', type = a2$type)
                a2$sizeExprs <- c(list(1), origSizeExprs)
           }
            else { 
                a2 <- insertExprClassLayer(code, 2, 'asCol', type = a2$type)
                a2$sizeExprs <- c(origSizeExprs, list(1))
            }
        }
    }
    code$nDim <- 2
    code$sizeExprs <- list(a1$sizeExprs[[1]], a2$sizeExprs[[2]])
    code$type <- setReturnType(code$name, arithmeticOutputType(a1$type, a2$type))
    code$toEigenize <- 'yes'
    assertMessage <- paste0("Run-time size error: expected ", deparse(a1$sizeExprs[[2]]), " == ", deparse(a2$sizeExprs[[1]]))
    newAssert <- identityAssert(a1$sizeExprs[[2]], a2$sizeExprs[[1]], assertMessage)
    if(is.null(newAssert))
        return(asserts)
    else
        return(c(asserts, list(newAssert)))
     ##   return(c(list(identityAssert(a1$sizeExprs[[2]], a2$sizeExprs[[1]], assertMessage)), asserts))
}

sizeSolveOp <- function(code, symTab, typeEnv) { ## this is for solve(A, b) or forwardsolve(A, b). For inverse, use inverse(A), not solve(A)
    if(length(code$args) != 2) stop(exprClassProcessingErrorMsg(code, 'sizeSolveOp called with argument length != 2.'), call. = FALSE)
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    a2 <- code$args[[2]]
    if(!(inherits(a1, 'exprClass') & inherits(a2, 'exprClass'))) stop(exprClassProcessingErrorMsg(code, 'In sizeSolveOp: expecting both arguments to be exprClasses.'), call. = FALSE)
    if(a1$nDim != 2)  stop(exprClassProcessingErrorMsg(code, 'In sizeSolveOp: first argument to a matrix solver must be a matrix.'), call. = FALSE)
    if(!any(a2$nDim == 1:2)) stop(exprClassProcessingErrorMsg(code, 'In sizeSolveOp: second argument to a matrix solver must be a vector or matrix.'), call. = FALSE)
    code$type <- setReturnType(code$name, 'double')
    code$nDim <- a2$nDim  ## keep the same dimension as the 2nd argument
    if(code$nDim == 1) { code$sizeExprs <- c(a1$sizeExprs[[1]])
                     } else { code$sizeExprs <- c(a1$sizeExprs[[1]], a2$sizeExprs[[2]]) }
    code$toEigenize <- 'yes'
    assertMessage <- paste0("Run-time size error: expected ", deparse(a1$sizeExprs[[1]]), " == ", deparse(a1$sizeExprs[[2]]))
    assert1 <- identityAssert(a1$sizeExprs[[1]], a1$sizeExprs[[2]], assertMessage)
    assertMessage <- paste0("Run-time size error: expected ", deparse(a1$sizeExprs[[1]]), " == ", deparse(a2$sizeExprs[[1]]))
    assert2 <- identityAssert(a1$sizeExprs[[1]], a2$sizeExprs[[1]], assertMessage)
    asserts <- c(asserts, assert1, assert2)
    return(asserts)
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
                if(!rowOK) stop(exprClassProcessingErrorMsg(code, 'In setAsRowOrCol: Cannot convert to row.'), call. = FALSE)
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
                if(!colOK) stop(exprClassProcessingErrorMsg(code, 'In setAsRowOrCol: Cannot convert to col.'), call. = FALSE)
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
    if(length(code$args) != 2) stop(exprClassProcessingErrorMsg(code, 'sizeBinaryCwise called with argument length != 2.'), call. = FALSE)

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
        a1type <- storage.mode(a1)
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
        a2type <- storage.mode(a2)
    }
    
    ## Choose the output type by type promotion
    if(length(a1type) == 0) {warning('Problem with type of arg1 in sizeBinaryCwise', call. = FALSE); browser()}
    if(length(a2type) == 0) {warning('Problem with type of arg2 in sizeBinaryCwise', call. = FALSE); browser()}
    code$type <- setReturnType(code$name, arithmeticOutputType(a1type, a2type))

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
        ## If so, wrap the 1D in asRow or asCol to orient it later for Eigen
        if(a1DropNdim == 1 & a2DropNdim == 1) {

            ## Hey, I think this is wrong: I think we should check the aXDropNdims
            if(a1nDim > 2 | a2nDim > 2) stop(exprClassProcessingErrorMsg(code, 'In sizeBinaryCwise: Dimensions do not match and/or NIMBLE will not handle Array + Vector for dim(Array) > 2.'), call. = FALSE) 

            ## a1 is 2D and a2 is 1D
            if(a1nDim == 2 & a2nDim == 1) {
                a1IsCol <- identical(a1sizeExprs[[2]], 1)
                asFun <- if(a1IsCol) 'asCol' else 'asRow'
                a2 <- insertExprClassLayer(code, 2, asFun, type = a2type)
                a2$sizeExprs <- a1sizeExprs
                a2$nDim <- a1nDim
                a1ind <- if(a1IsCol) 1 else 2
                if(!is.numeric(a1sizeExprs[[a1ind]]) | !is.numeric(a2sizeExprs[[1]])) { ## Really do want original a2sizeExprs
                    assertMessage <- paste0("Run-time size error: expected ", deparse(a1sizeExprs[[a1ind]]), " == ", deparse(a2sizeExprs[[1]]))
                    thisAssert <- identityAssert(a1sizeExprs[[a1ind]], a2sizeExprs[[1]], assertMessage)
                    if(!is.null(thisAssert)) asserts[[length(asserts) + 1]] <- thisAssert
                } else {
                    if(a1sizeExprs[[a1ind]] != a2sizeExprs[[1]]) stop(exprClassProcessingErrorMsg(code, 'In sizeBinaryCwise: Fixed size mismatch.'), call. = FALSE)
                }
                code$nDim <- a1nDim
                code$sizeExprs <- a1sizeExprs
            } else {
                a2IsCol <- identical(a2sizeExprs[[2]], 1)
                asFun <- if(a2IsCol) 'asCol' else 'asRow'
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
                    if(a1sizeExprs[[1]] != a2sizeExprs[[a2ind]]) stop(exprClassProcessingErrorMsg(code, 'In sizeBinaryCwise: Fixed size mismatch.'), call. = FALSE)
                }
                code$nDim <- a2nDim
                code$sizeExprs <- a2sizeExprs
            }
        } else {
            ## If at least one arg is a known scalar-equivalent, that case was handled above
            ## (But it's still not complete)
            ## Here is the case that nDims aren't equal and dropNdims aren't equal
            ## either.  We used to rely on typeEnv to keep track of when a size resulting from an operation is known to be 1 but realized that isn't safe if that operation is only conditionally executed at run time.
            ## Hence what will do now is assume the user has written valid code
            ## but add run-time size checks of which dimension must match
            ## This is currently limited in what it will handle
            ## Specifically, it assumes things should be columns
            assertMessage <- paste0("Run-time size error: expected ", deparse(a1sizeExprs[[1]]), " == ", deparse(a2sizeExprs[[1]]))
                thisAssert <- identityAssert(a1sizeExprs[[1]], a2sizeExprs[[1]], assertMessage)
                if(!is.null(thisAssert)) asserts[[length(asserts) + 1]] <- thisAssert

            if(a1nDim == 1 & a2nDim == 2) {
                assertMessage <- paste0("Run-time size error: expected ", deparse(a2sizeExprs[[2]]), " == ", 1)
                thisAssert <- identityAssert(a2sizeExprs[[2]], 1, assertMessage)
                if(!is.null(thisAssert)) asserts[[length(asserts) + 1]] <- thisAssert                
                code$sizeExprs <- a2sizeExprs
            } else {
                if(a1nDim == 2 & a2nDim == 1) {
                    assertMessage <- paste0("Run-time size error: expected ", deparse(a1sizeExprs[[2]]), " == ", 1)
                    thisAssert <- identityAssert(a1sizeExprs[[2]], 1, assertMessage)
                    if(!is.null(thisAssert)) asserts[[length(asserts) + 1]] <- thisAssert
                    code$sizeExprs <- a1sizeExprs
                } else {
                    stop(exprClassProcessingErrorMsg(code, 'In sizeBinaryCwise: Dimensions do not matchin a way that can be handled.'), call. = FALSE)
                }
            }
            code$nDim <- 2
        }
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
                    if(a1sizeExprs[[i]] != a2sizeExprs[[i]]) stop(exprClassProcessingErrorMsg(code, 'In sizeBinaryCwise: Fixed size mismatch.'), call. = FALSE)
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
                             nimArr_rmvt_chol = list(c(1, 2, 0, 0), ## dimensionality of ordered arguments AFTER the first, which is for the return value.  e.g. mean (1D), chol(2D), df(scalar), prec_param(scalar)
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
        if(length(code$args) < length(checkList[[1]])) stop(exprClassProcessingErrorMsg(code, 'In sizeRmultivarFirstArg: Not enough arguments provided.'), call. = FALSE)
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
        stop(exprClassProcessingErrorMsg(code, 'In sizeRmultivarFirstArg: Some argument(s) have the wrong dimension.'), call. = FALSE) 
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
        } else
            typeEnv$.ensureNimbleBlocks <- TRUE ## for purposes of sizeAssign, which recurses on assignment target after RHS
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
        stop(exprClassProcessingErrorMsg(code, 'In generalFunSizeHandler: Wrong number of arguments.'), call. = FALSE)
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
    if(inherits(returnType, 'symbolNimbleList')) {
        code$type <- 'symbolNimbleList'
        code$sizeExprs <- returnType
        code$toEigenize <- 'maybe'
        code$nDim <- 0
        liftIfAmidExpression <- TRUE
    } else {
        returnSymbolBasic <- inherits(returnType, 'symbolBasic')
        returnTypeLabel <- if(returnSymbolBasic) returnType$type else as.character(returnType[[1]])
        
        if(returnTypeLabel == 'void') {
            code$type <- returnTypeLabel
            code$toEigenize <- 'unknown'
            return(asserts)
        }
        returnNDim <- if(returnSymbolBasic) returnType$nDim
                      else if(length(returnType) > 1) as.numeric(returnType[[2]]) else 0
                                                                
                                        # returnSizeExprs <- if(returnTypeLabel == 'symbolNimbleList') symTab$getSymbolObject('return') else vector('list', returnNDim) ## This stays blank (NULLs), so if assigned as a RHS, the LHS will get default sizes
        returnSizeExprs <- vector('list', returnNDim) ## This stays blank (NULLs), so if assigned as a RHS, the LHS will get default sizes
        code$type <- returnTypeLabel
        code$nDim <- returnNDim
        code$sizeExprs <- returnSizeExprs
        code$toEigenize <- if(code$nDim == 0) 'maybe' else 'no'
        liftIfAmidExpression <- code$nDim > 0
    }
    
    if(liftIfAmidExpression) {
        if(!(code$caller$name %in% c('{','<-','<<-','='))) {
            asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        } else
            typeEnv$.ensureNimbleBlocks <- TRUE
    }
    return(asserts)
}
