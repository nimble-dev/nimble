sizeProc_storage_mode <- function(x) {
    if(is.na(x))
        if(storage.mode(x) == 'logical')
            return('double') ## promote logical NA to double for compiled code.
    storage.mode(x)
}
    
assignmentAsFirstArgFuns <- c('nimArr_rmnorm_chol',
                              'nimArr_rmvt_chol',
                              'nimArr_rlkj_corr_cholesky',
                              'nimArr_rwish_chol',
                              'nimArr_rinvwish_chol',
                              'nimArr_rcar_normal',
                              'nimArr_rcar_proper',
                              'nimArr_rmulti',
                              'nimArr_rdirch',
                              'getValues',
                              'getValuesIndexRange',
                              'initialize',
                              'setWhich',
                              'setRepVectorTimes',
                              'assignVectorToNimArr',
                              'dimNimArr',
                              'assignNimArrToNimArr')

setSizeNotNeededOperators <- c('setWhich',
                               'setRepVectorTimes',
                               'SEXP_2_NimArr',
                               'nimVerbatim')

operatorsAllowedBeforeIndexBracketsWithoutLifting <- c('map',
                                                       'dim',
                                                       'mvAccessRow',
                                                       'nfVar')

sizeCalls <- c(
    makeCallList(binaryOperators, 'sizeBinaryCwise'),
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
    makeCallList(nimbleListReturningOperators, 'sizeNimbleListReturningFunction'),
    nimOptim = 'sizeOptim',
    nimOptimDefaultControl = 'sizeOptimDefaultControl',
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
         '[' = 'sizeIndexingBracket',
         '[[' = 'sizeDoubleBracket', ## for nimbleFunctionList, this will always  go through chainedCall(nfList[[i]], 'foo')(arg1, arg2)
         chainedCall = 'sizeChainedCall',
         nfVar = 'sizeNFvar',
         map = 'sizemap', 
         ':' = 'sizeColonOperator',
         'if' = 'recurseSetSizes', ##OK
         'while' = 'recurseSetSizes',
         'for' = 'sizeFor', 
         cppPointerDereference = 'sizeCppPointerDereference',
         values = 'sizeValues',
         '(' = 'sizeUnaryCwise',
         setSize = 'sizeSetSize', 
         resizeNoPtr = 'sizeResizeNoPtr', ## may not be used any more 
         nimArr_rcat = 'sizeScalarRecurse',
         nimArr_rinterval = 'sizeScalarRecurse',
         nimPrint = 'sizeforceEigenize',
         nimDerivs = 'sizeNimDerivs',
         nimDerivs_calculate = 'sizeNimDerivsCalculate',
         as.integer = 'sizeUnaryCwise', 
         as.numeric = 'sizeUnaryCwise',
         nimArrayGeneral = 'sizeNimArrayGeneral',
         setAll = 'sizeOneEigenCommand',
         voidPtr = 'sizeVoidPtr',
         run.time = 'sizeRunTime',
         bessel_k = 'sizeRecyclingRuleBesselK',
         PROTECT = 'sizePROTECT',
         NimArr_2_SEXP = 'sizePROTECT', 
         Reval = 'sizeReval',
         Rf_eval = 'sizeReval',
         nimbleConvert = 'sizeNimbleConvert',
         nimbleUnconvert = 'sizeNimbleUnconvert',
         asReturnSymbol = 'sizeAsReturnSymbol'),
    makeCallList(scalar_distribution_dFuns, 'sizeRecyclingRule'),
    makeCallList(scalar_distribution_pFuns, 'sizeRecyclingRule'),
    makeCallList(scalar_distribution_qFuns, 'sizeRecyclingRule'),
    makeCallList(scalar_distribution_rFuns, 'sizeRecyclingRuleRfunction'),
    makeCallList(distributionFuns[
        !(distributionFuns %in% c(scalar_distribution_dFuns,
                                  scalar_distribution_pFuns,
                                  scalar_distribution_qFuns,
                                  scalar_distribution_rFuns))
    ], 'sizeScalarRecurseAllowMaps'),
    ## R dist functions that are not used by NIMBLE but we allow in DSL
    makeCallList(paste0(c('d','q','p'), 't'), 'sizeRecyclingRule'),
    rt = 'sizeRecyclingRuleRfunction',
    makeCallList(paste0(c('d','q','p'), 'exp'), 'sizeRecyclingRule'),
    rexp = 'sizeRecyclingRuleRfunction',
    makeCallList(c('nimAnyNA','nimAnyNaN'), 'sizeScalarRecurse'),
    makeCallList(c('nimArr_dmnorm_chol',
                   'nimArr_dmvt_chol',
                   'nimArr_dlkj_corr_cholesky',
                   'nimArr_dwish_chol',
                   'nimArr_dinvwish_chol',
                   'nimArr_dcar_normal',
                   'nimArr_dcar_proper',
                   'nimArr_dmulti',
                   'nimArr_dcat',
                   'nimArr_dinterval',
                   'nimArr_ddirch'), 'sizeScalarRecurseAllowMaps'),
    makeCallList(c('nimArr_rmnorm_chol',
                   'nimArr_rmvt_chol',
                   'nimArr_rlkj_corr_cholesky',
                   'nimArr_rwish_chol',
                   'nimArr_rinvwish_chol',
                   'nimArr_rcar_normal',
                   'nimArr_rcar_proper',
                   'nimArr_rmulti',
                   'nimArr_rdirch'), 'sizeRmultivarFirstArg'),
    makeCallList(c('decide',
                   'size',
                   'getsize',
                   'getNodeFunctionIndexedInfo',
                   'endNimbleTimer'), 'sizeScalar'),
    makeCallList(c('calculate',
                   'calculateDiff',
                   'getLogProb'), 'sizeScalarModelOp'),
    simulate = 'sizeSimulate',
    makeCallList(c('blank',
                   'nfMethod',
                   'getPtr',
                   'startNimbleTimer'), 'sizeUndefined'), ##'nimFunListAccess'
    passByMap = 'sizePassByMap',
    ADbreak = 'sizeADbreak')

scalarOutputTypes <- list(decide = 'logical',
                          size = 'integer',
                          nimAnyNA = 'logical',
                          nimAnyNaN = 'logical',
                          '!' = 'logical',
                          getNodeFunctionIndexedInfo = 'double',
                          endNimbleTimer = 'double')

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
expressionSymbolTypeReplacements <- c('symbolNimbleListGenerator', 'symbolNimbleList', 'symbolNimbleFunction', 'symbolMemberFunction')

exprClasses_setSizes <- function(code, symTab, typeEnv) { ## input code is exprClass
    ## name:
   if(code$isName) {
        ## If it doesn't exist and must exist, stop
        if(code$name != "") { ## e.g. In A[i,], second index gives name==""
            if(!exists(code$name, envir = typeEnv, inherits = FALSE)) {
                if(symTab$symbolExists(code$name, TRUE)) {
                    thisSymbolObject <- symTab$getSymbolObject(code$name, TRUE)
                    code$type <- class(thisSymbolObject)[1]
                    if(code$type %in% expressionSymbolTypeReplacements){
                      code$type <- thisSymbolObject$type
                      code$sizeExprs <- thisSymbolObject
                    }
                } else {
                    code$type <- 'unknown'
                    if(!typeEnv$.AllowUnknowns)
                        if(identical(code$name, 'pi')) { ## unique because it may be encountered anew on on RHS and be valid
                            assign('pi',
                                   exprTypeInfoClass$new(nDim = 0,
                                                         type = 'double',
                                                         sizeExprs = list()),
                                   envir = typeEnv)
                            symTab$addSymbol(
                                symbolBasic(name = 'pi',
                                            type = 'double',
                                            nDim = 0))
                            code$nDim <- 0
                            code$type <- 'double'
                            code$sizeExprs <- list()
                            code$toEigenize <- 'maybe'
                        } else {
                            if(nimbleOptions('errorIfMissingNFVariable')) {
                                stop("variable `",
                                     code$name,
                                     "` is not available.",
                                     call.=FALSE)
                            } else
                                messageIfVerbose("  [Warning] Variable `", code$name, "` is not available.")
                        }
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
            ## Add RCfunctions to neededRCfuns.
            if(typeEnv[['.allowFunctionAsArgument']]) { 
                if(exists(code$name) && is.rcf(get(code$name))) {
                    nfmObj <- environment(get(code$name))$nfMethodRCobject
                    uniqueName <- nfmObj$uniqueName
                    if (is.null(typeEnv$neededRCfuns[[uniqueName]])) {
                        typeEnv$neededRCfuns[[uniqueName]] <- nfmObj
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
                    newAsserts <-
                        exprClasses_setSizes(code$args[[i]], symTab, typeEnv)
                    code$args[[i]]$assertions <-
                        if(is.null(newAsserts)) list() else newAsserts
                }
            }
            return(invisible(NULL))
        }
        sizeCall <- sizeCalls[[code$name]]
        if(!is.null(sizeCall)) {
            if(.nimbleOptions$debugSizeProcessing) {
                browser()
                eval(
                    substitute(
                        debugonce(XYZ),
                        list(XYZ = as.name(sizeCall))
                    )
                )
            }
            test0 <- eval(call(sizeCall, code, symTab, typeEnv))
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
                if(length(uniqueName)==0)
                    stop(
                        exprClassProcessingErrorMsg(
                            code,
                            'In size processing: A no-setup nimbleFunction with no internal name is being called.'),
                        call. = FALSE)
                ## new with nimbleLists: we need to initiate compilation here so we can get full returnType information, including of nimbleLists
                RCfunProc <-
                    typeEnv$.nimbleProject$compileRCfun(obj,
                                                        initialTypeInference = TRUE)
                
                if(is.null(typeEnv$neededRCfuns[[uniqueName]])) {
                    if(!identical(RCfunProc$RCfun$uniqueName, typeEnv$.myUniqueName))
                        typeEnv$neededRCfuns[[uniqueName]] <- nfmObj
                }
                
                return(sizeRCfunction(code, symTab, typeEnv, nfmObj, RCfunProc))
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

sizeADbreak <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(length(code$args) != 1)
        stop(exprClassProcessingErrorMsg(code, paste0('ADbreak must have exactly one argument.')), call. = FALSE)
    if(code$args[[1]]$nDim != 0)
        stop(exprClassProcessingErrorMsg(code, paste0('The argument to ADbreak must be scalar.')), call. = FALSE)
    code$nDim <- 0
    code$type <- code$args[[1]]$type
    code$toEigenize <- 'no'
    code$sizeExprs <- list()
    return(if(is.null(asserts)) list() else asserts)
}

## This is used by nimbleExternalCall.
## When the external call is provided as foo, returning e.g. double(0),
## we end up needing a line of code RETURNVALUE <- foo(args).
## To get the type of RETURNVALUE, we wrap that as RETURNVALUE <- asReturnSymbol(foo(args), type, nDim)
sizeAsReturnSymbol <- function(code, symTab, typeEnv) {
    returnType <- code$args[[2]]
    returnNDim <- code$args[[3]]
    code$args <- list(code$args[[1]])
    code$args[[1]]$type <- returnType
    code$args[[1]]$nDim <- returnNDim
    code$args[[1]]$toEigenize <- 'no'
    code$args[[1]]$sizeExprs <- NULL
    removeExprClassLayer(code, 1)
    list()
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
    ## experimentalNewSizeProcessing: code$name change step stays here
    ## experimentalNewSizeProcessing: because the 3 cases are not implementation-specific
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

    if(!nimbleOptions('experimentalSelfLiftStage')) {
        if(!(code$caller$name %in% assignmentOperators)) {
            asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        }
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
    code$toEigenize <- 'yes' ## toEigen: N.B. This had TRUE
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
    ## toEigen: keep this lift here for now, since it sets up sizes.
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
    code$toEigenize <- 'yes'
    return(asserts)
}

sizeRecyclingRuleBesselK <- function(code, symTab, typeEnv) { ## also need an entry in eigenization.
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    numArgs <- length(code$args)

    # this is easily relaxed but not clear the functionality would ever be needed...
    if(!is.numeric(code$args[[3]]) && !identical(code$args[[3]]$nDim, 0))
        stop("In besselK, 'expon.scaled' must be a single value.")

    if(numArgs != 3) stop("Expecting two or three arguments for besselK function.")
    recycleArgs <- c(TRUE, TRUE, FALSE)

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

concatenateIntermLabelMaker <- labelFunctionCreator("ConcatenateInterm")

sizeConcatenate <- function(code, symTab, typeEnv) { ## This is two argument version
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    ## overall strategy is to separate runs of scaalrs and non-scalars
    ## also in C++ we don't take arbitrary arguments.  Instead we chain together calls in groups of 4
    ##     e.g. c(a1, a2, a3, a4, a5) will become c( c(a1, a2, a3, a4), a5)
    
    ## first puzzle is with nimC(scalar1, scalar2, vector1, scalar3)
    ## we need to extract the runs of scalars like (scalar1, scalar2), so they can be packed up in an object together.
    isScalar <- unlist(lapply(code$args, function(x) if(inherits(x, 'exprClass')) x$nDim == 0 else TRUE))
    ## run length encoding: This native R function returns information about repeats, so we can figure out how long each run of scalars is
    argRLE <- rle(isScalar)
    ## How many arguments will we have after packing scalars together into single objects:
    newNumArgs <- sum(argRLE$values) + sum(argRLE$lengths[!argRLE$values]) ## number of scalar runs + sum of non-scalar runs * run-lengths
    newArgs <- vector(length(newNumArgs), mode = 'list')
    iInput <- 1
    iOutput <- 1
    for(i in seq_along(argRLE$values)) {
        thisLength <- argRLE$lengths[i]
        if(!(argRLE$values[i])) { ## it is a run of non-scalars, so pack them into the new argument list, newArgs
            newArgs[(iOutput-1) + (1:thisLength)] <- code$args[(iInput-1) + (1:thisLength)]
            iInput <- iInput + thisLength
            iOutput <- iOutput + thisLength
        } else { ## it is a run of scalars, so construct an object for them
            newTempFixedName <- concatenateIntermLabelMaker()
            newTempVecName <- concatenateIntermLabelMaker()
            ## Construct:
            ## concatenateTemp(ConcatenateInterm_1),
            ##   concatenateTemp is not output to C++. It is a placeholder
            newExpr <- exprClass$new(isName = FALSE, isCall = TRUE, isAssign = FALSE, name = "concatenateTemp", nDim = 1, sizeExprs = list(thisLength), type = 'double')
            setArg(newExpr, 1, exprClass$new(isName = TRUE, isCall = FALSE, isAssign = FALSE, name = newTempVecName, nDim = 1, sizeExprs = list(thisLength), type = 'double'))

            ## hardCodedVectorInitializer is a wrapper for the "contents1, contents2, ..." below 
            valuesExpr <- quote(hardCodedVectorInitializer())
            thisType <- 'logical'
            for(j in 1:thisLength) {
                thisArgIndex <- iInput - 1 + j
                if(inherits(code$args[[thisArgIndex]], 'exprClass')) {
                    if(!code$args[[thisArgIndex]]$isName) ## a little heavy-handed: lift any expression of any kind
                        ## to avoid dealing with eigen or other handling inside initialization values
                        ## This is necessary for cases like nimC(model[[node]][2], 1.2)
                        ## because model[[node]] is a map
                        asserts <- c(asserts, sizeInsertIntermediate(code, thisArgIndex, symTab, typeEnv))
                    thisType <- arithmeticOutputType(thisType, code$args[[thisArgIndex]]$type)
                } else {
                    thisType <- sizeProc_storage_mode(code$args[[thisArgIndex]]) ##'double'
                }
                ## Putting a map, or a values access, through parse(nimDeparse) won't work
                ## So we lift any expression element above.
                ## This could be done more cleanly with more coding work.
                  valuesExpr[[j+1]] <- parse(text = nimDeparse(code$args[[thisArgIndex]]), keep.source = FALSE)[[1]]
            }
            newExpr$type <- thisType
            newExpr$args[[1]]$type <- thisType
            iInput <- iInput + thisLength
            if(thisType == 'integer') thisType <- 'int'
            if(thisType == 'logical') thisType <- 'bool'
            ## MAKE_FIXED_VECTOR("ConcatenateInterm_2", "ConcatenateInterm_1", numArgs, values, type, allowAD) goes through a customized output generator
            ##  to create something like
            ##    double ConcatenateIterm_1[] = {contents1, contents2}
            ##    std::vector<double> ConcatenateInterm_2(ConcatenateInterm_1, ConcatenateInterm_1 + length)
            ##  so there is one intermediate whose only purpose is to achieve initialization by value and a second intermediate copied from the first.
            ##     The second intermediate can later be used in the templated nimCd/nimCi/nimCb
            ##
            ## allowAD flags whether double should be changed to CppAD::AD<double> for _AD_ versions
            newAssert <- substitute(MAKE_FIXED_VECTOR(newTempVecName, newTempFixedName, thisLength, valuesExpr, thisType, allowAD),
                                    list(newTempVecName = newTempVecName, newTempFixedName = newTempFixedName,
                                         thisLength = as.numeric(thisLength),
                                         valuesExpr = valuesExpr,
                                         thisType = thisType,
                                         allowAD = !isTRUE(typeEnv$.avoidAD)))
            newAssert <- as.call(newAssert)
            asserts <- c(asserts, list(newAssert))
            newArgs[[iOutput]] <- newExpr
            iOutput <- iOutput + 1
        }
    }

    ## Next step: chain together multiple calls:
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
        newExprList[[i]] <- exprClass$new(isName = FALSE, isCall = TRUE, isAssign = FALSE, name = 'nimC', nDim = 1, toEigenize = 'yes', type = 'double')
        for(j in seq_along(splitArgIDs[[i]])) setArg(newExprList[[i]], j, newArgs[[splitArgIDs[[i]][j]]])
    }

    ## Last step is to set up nesting and make sizeExprs for each constructed argument
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
}

sizeRep <- function(code, symTab, typeEnv) {
    ## if times is a vector: If length.out is provided, times is always ignored
    ## otherwise lift and use assignEigenToNIMBLE
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    xIsExpr <- inherits(code$args[[1]], 'exprClass')
    code$type <- if(xIsExpr) code$args[[1]]$type else 'double'

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
  old_avoidAD <- typeEnv$.avoidAD
  typeEnv$.avoidAD <- TRUE  # The makes the lifted nodes for times and/or length.out be put on the ignore list for AD
  on.exit(typeEnv$.avoidAD <- old_avoidAD)
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
    ## The code looks like: nimListDef$new(a = 10, b = 12).
    ## We want to change code$caller to :
    ## { nimList <- nimListDef$new()
    ## nimList$a <- 10
    ## nimList$b <- 12 }.
    ## We accomplish this by copying code, getting arguments (e.g. a = 10, b = 12) from copied code and turning them into assignment 
    ## exprs in code$caller, and setting first argument of code$caller to be nimList <- nimListDef$new().
    
    listDefName <- code$args[[1]]$name
    if(symTab$symbolExists(listDefName, inherits = TRUE)){
        listST <- symTab$getSymbolObject(listDefName, inherits = TRUE)
    } else {
        ## We need to establish the symbol and needed type.        
        nlDef <- get(listDefName)
        ## Need the nimbleProject!
        nlp <- typeEnv$.nimbleProject$compileNimbleList(nlDef, initialTypeInference = TRUE)
        className <- nl.getListDef(nlDef)$className
        if(is.null(typeEnv$neededRCfuns[[className]])) {
            newSym <- symbolNimbleList(name = listDefName, nlProc = nlp)
            typeEnv$neededRCfuns[[className]] <- newSym
        }
        newDefSym <- symbolNimbleListGenerator(name = listDefName, nlProc = nlp)
        symTab$addSymbol(newDefSym)
        listST <- newDefSym
    }
    code$type <- "nimbleList"
    code$sizeExprs <- listST
    code$toEigenize <- "maybe"
    code$nDim <- 0
    
    asserts <- list()
    asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, useArgs = c(TRUE, rep(FALSE, length(code$args)-1))))
    if(!(code$caller$name %in% assignmentOperators)){
        intermediateAsserts <- sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv)
        ## intermediateAsserts can potentially have size setting stuff from sizeAssignAfterRecursing.
        ## Not sure if that would ever happen in this context, but to be safe we'll use last element as the actual intermediate assignment.
        ## Embed the intermediate assignment in a '{' (so insertAssertions will recurse on it) and recurse on it.
        numIntermAsserts <- length(intermediateAsserts)
        bracketedIntermAssert <- newBracketExpr(intermediateAsserts[numIntermAsserts])
        exprClasses_setSizes(bracketedIntermAssert, symTab, typeEnv)
        intermediateAsserts[[numIntermAsserts]] <- bracketedIntermAssert
        asserts <- c(asserts, intermediateAsserts)
        return(asserts)
    }
    if(length(code$args) <= 1) return(asserts)  ## There are no args to process.

    RnewExprs <- list()
    newExprs <- list()
    RnfVarExprs <- list()
    nfVarExprs <- list()
    exprCounter <- 1
    originalCode <- code 
    listElements <- listST$nlProc$symTab$getSymbolNames()
    RlistNameExpr <- nimbleGeneralParseDeparse(originalCode$caller$args[[1]])    
    for(i in seq_along(listElements)) {
        thisVarName <- listElements[i]
        if(!is.null(originalCode$args[[thisVarName]])) {
            ## Skip first arg, which will be name of nlDef, then check if value is "".
            if(!inherits(originalCode$args[[thisVarName]], 'exprClass') || (originalCode$args[[thisVarName]]$name != "")) {
                ## nfVar(A, 'x') for whichever element name it's on ('x')
                RnfVarExprs[[exprCounter]] <- substitute(nfVar(A, X), list(A = RlistNameExpr, X = thisVarName))
                ## nfVar(A, 'x') <- y or whatever code was provided (already recursed for size processing)
                RnewExprs[[exprCounter]] <- substitute(A <- B, list(A = RnfVarExprs[[exprCounter]],
                                                                    B = nimbleGeneralParseDeparse(originalCode$args[[thisVarName]])))
                exprCounter <- exprCounter + 1
            }
        }
    }
    if(length(RnewExprs) == 0) return(asserts)  ## All args have already been specified.
    
    ## Embed RnewExprs in a '{' expression.
    RbracketNewExprs <- quote(after({}))
    RbracketNewExprs[[2]][2:(length(RnewExprs) + 1)] <- RnewExprs
    bracketNewExprs <- RparseTree2ExprClasses(RbracketNewExprs)
    ## Need to install assignment target in symTab if necessary so that it
    ## will be there for recursion in the following step.
    assignmentTarget <- code$caller$args[[1]]
    if(assignmentTarget$isName) {
        if(!symTab$symbolExists(assignmentTarget$name, TRUE)) {
            symTab$addSymbol(symbolNimbleList(name = assignmentTarget$name, type = code$type, nlProc = code$sizeExprs$nlProc))
        }
    }
    ## Recurse into element assignments.
    exprClasses_setSizes(bracketNewExprs$args[[1]], symTab, typeEnv)
    asserts <- c(asserts, list(bracketNewExprs))
    if(length(code$args) > 1) ## TODO Remove this conditional, since this should always be true if we make it this far.
        code$args <- code$args[1]
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
## nimArrayGeneral(type(character), nDim, dim (c(sizeExpr1, ...)), value, init (logical), fillZeros, recycle, unpackNDim(optional))
## nimArrayGeneral(     arg1,       arg2,       arg3,              arg4,     arg5       ,    arg6  ,  arg7   ,    arg8            )
sizeNimArrayGeneral <- function(code, symTab, typeEnv) {
    useArgs <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
    if(!is.null(code$args[['unpackNDim']])) useArgs <- c(useArgs, TRUE)
    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = useArgs)  ## recurse on initialValue and initialLogical only

    ## some checking
    if(inherits(code$args[['init']], 'exprClass'))
        if(!(code$args[['init']]$nDim == 0)) stop(exprClassProcessingErrorMsg(code, paste0('init argument to numeric, logical, integer, matrix or array must be scalar')), call. = FALSE)
    
    type <- code$args[['type']] 
    nDim <- code$args[['nDim']] 
    unpackNDim <- if(!is.null(code$args[['unpackNDim']])) code$args[['unpackNDim']] else FALSE ##if(length(code$args) > 5) code$args[[6]] else FALSE

    cSizeExprs <- code$args[['dim']] ## these are the size expressions encompassed by collectSizes(), needed for purposes of the C++ line to be generated
    if(!inherits(cSizeExprs, 'exprClass'))        stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (i) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)
    if(cSizeExprs$name != 'collectSizes')         stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (ii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)

    if(unpackNDim) { ## This means length of dim unknown at compile time but nDim explicitly provided, so we construct c(dim[1], dim[2]), etc.
        asserts <- c(asserts, recurseSetSizes(cSizeExprs, symTab, typeEnv))
        if(!cSizeExprs$args[[1]]$isName)
            asserts <- c(asserts, sizeInsertIntermediate(cSizeExprs, 1, symTab, typeEnv)) ## this intermediate goes a layer down the AST, but works
        if(length(cSizeExprs$args[[1]]$sizeExprs) == 0) { ## The argument expression evaluates to scalar
            if(nDim == -1) {
                nDim <- 1
                code$args[['nDim']] <- 1
            }
            if(nDim == 1) unpackNDim <- FALSE             ## and that's ok because nDim given as 1
        }
        if(unpackNDim) {
            if(length(cSizeExprs$args[[1]]$sizeExprs) != 1) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (ii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)
            if(nDim == -1) {## code for nDim not given but dim given as expression
                if(!is.numeric(cSizeExprs$args[[1]]$sizeExprs[[1]])) stop()
                nDim <- cSizeExprs$args[[1]]$sizeExprs[1]
                if(nDim < 1 | nDim > 4) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (iii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)
                code$args[['nDim']] <- nDim
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
        if(nDim == -1) {  ## nDim wasn't provided (to nimArray) and dim was an expression, so it ought to be a scalar
            if(nonScalarWhereNeeded) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (iv) with sizes or dim to numeric, logical, integer, matrix or array.  It looks like dim argument was non-scalar but nDim was not provided.  If the dim argument to array (or nimArray) is a vector, you must also provide nDim argument to say how many dimensions will be used.')), call. = FALSE)
            nDim <- code$args[['nDim']] <- 1 
        } else { ## call was from numeric, integer or logical
            if(nonScalarWhereNeeded) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (v) with sizes or dim to numeric, logical, integer, matrix or array.  It looks like length argument was non-scalar.')), call. = FALSE)
        }       
    }
        
    ## if it is a call to matrix() and the value argument is non-scalar,
    ## we will generate it in C++ as nimNewMatrix
    useNewMatrix <- FALSE
    if(nDim == 2) {
        if(inherits(code$args[['value']], 'exprClass'))
            if(code$args[['value']]$nDim > 0)
                useNewMatrix <- TRUE  ## use eigen-compatible C++
    }
                
     if(code$args[['nDim']] != length(cSizeExprs$args)) stop(exprClassProcessingErrorMsg(code, paste0('Something wrong (iii) with sizes or dim to numeric, logical, integer, matrix or array')), call. = FALSE)

    annotationSizeExprs <- lapply(cSizeExprs$args, nimbleGeneralParseDeparse) ## and this is for purposes of the sizeExprs in the AST exprClass object
    
    missingSizes <- unlist(lapply(cSizeExprs$args, function(x) identical(x, as.numeric(NA)) | identical(x, NA))) ## note is.na doesn't work b/c the argument can be an expression and is.na warns on that ##old: identical, as.numeric(NA)))
    ## only case where we do something useful with missingSizes is matrix(value = non-scalar, ...)
    if(any(missingSizes)) {
        ## modify sizes in generated C++ line
        if(useNewMatrix) cSizeExprs$args[missingSizes] <- -1
        else cSizeExprs$args[missingSizes] <- 1

        ## modify annotation sizeExprs
        totalInputLengthExpr <- if(inherits(code$args[['value']], 'exprClass')) productSizeExprs(code$args[['value']]$sizeExprs) else 1 ## should always be exprClass in here anyway
        ## see newMatrixClass in nimbleEigen.h
        if(missingSizes[1]) { ## missing nrow
            if(missingSizes[2]) { ## missing both
                annotationSizeExprs[[1]] <- totalInputLengthExpr
                annotationSizeExprs[[2]] <- 1
            } else { ## ncol provided
                annotationSizeExprs[[1]] <- substitute(calcMissingMatrixSize(A, B), 
                                              list(A = totalInputLengthExpr,
                                                   B = annotationSizeExprs[[2]]))
            }
        } else { ## nrow provided, ncol missing (is both provided, we wouldn't be in this code
                annotationSizeExprs[[2]] <- substitute(calcMissingMatrixSize(A, B), 
                                              list(A = totalInputLengthExpr,
                                                   B = annotationSizeExprs[[1]]))
        }
    }
    
    asserts <- c(asserts, recurseSetSizes(cSizeExprs, symTab, typeEnv))
    if(!(type %in% c('double', 'integer', 'logical')))       stop('unknown type in nimArrayGeneral')
    ## Three possible calls can be emitted by choice of code$name: initialize (this becomes a NimArr member function call.  It is used if initialization is scalar, to be repeated); assignNimArrToNimArr (this becomes a call to assignNimArrToNimArr.  It is used if initialization is non-scalar and the object being created is not a matrix; nimNewMatrix[D|I|B] (this has the same name in C++.  It is used if initialization is non-scalar and the object being created is a matrix.  It creates an eigen-compatible object within an expression).
    code$name <- 'initialize' ## may be replaced below if useNewMatrix
    if(inherits(code$args[['value']], 'exprClass'))
        if(code$args[['value']]$nDim > 0)
            code$name <- 'assignNimArrToNimArr' ## could be replaced by nimNewMatrix[D|B|I] below

    ## rearrange arguments
    if(code$name == 'assignNimArrToNimArr')
        if(!useNewMatrix) 
            code$args <- c(code$args[4:7], cSizeExprs$args)  ##  args: initialize(value, init, fillZeros, recycle, sizeExpr1, sizeExpr2, etc...)
        else
            code$args <- c(code$args[c(4,5,7)], cSizeExprs$args)  ##  fillZeros has no role in this case.  nimNewMatrix creates an eigen object that has to return something for each element, so it will use a zero anyway.
    else
        code$args <- c(code$args[4:7], cSizeExprs$args) ## actually this turned out the same as for assignNimArrToNimArr.  

    ## fix code/caller relationships in AST
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

    ## check for nimNewMatrix case
    if(useNewMatrix) {
        suffix <- 'D'
        if(code$type == 'integer') suffix <- 'I'
        if(code$type == 'logical') suffix <- 'B'
        code$name <- paste0("nimNewMatrix", suffix)
        code$toEigenize <- "yes"
    } else {
        ## otherwise, lift values arg if necessary
        if(inherits(code$args[['value']], 'exprClass')) ## was re-ordered here
            if(!(code$args[['value']]$isName))
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
    }

    ## if(!useNewMatrix)
    ##     if(inherits(code$caller, 'exprClass'))
    ##         if(!(code$caller$name %in% assignmentOperators)) {
    ##             if(!is.null(code$caller$name)) {
    ##                 asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    ##             }
    ##         } else
    ##             typeEnv$.ensureNimbleBlocks <- TRUE

    if(!useNewMatrix)
      if(inherits(code$caller, 'exprClass')) {
        liftNewArr <- FALSE
        if(!is.null(code$caller$name)) {
          if(!(code$caller$name %in% assignmentOperators)) {
            liftNewArr <- TRUE
          } else {
            if(!isTRUE(code$caller$args[[1]]$isName))
              liftNewArr <- TRUE
          }
        }
        if(liftNewArr)
          asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        else
          typeEnv$.ensureNimbleBlocks <- TRUE
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
    timerName <- IntermLabelMaker()
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



## a$b becomes nfVar(a, 'b')
sizeNFvar <- function(code, symTab, typeEnv) {
    ## toEigen: Is it correct that this does not mark toEigen?
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
    isSymList <- objectType == 'nimbleList'
    
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
        stop(exprClassProcessingErrorMsg(code, 'In sizeNFvar: First argument is not a nimbleFunction or a nimbleList.\nList-like syntax with `[[` and `$` in nimbleFunction run code can only be done with a nimbleFunction or nimbleList.'), call. = FALSE)
    nfProc <- if(isSymFunc) symbolObject$nfProc else symbolObject$nlProc
    
    if(is.null(nfProc)) {
        stop(exprClassProcessingErrorMsg(code, 'In handling X$Y: Symbols for X have not been set up.'), call. = FALSE)
    }
    memberName <- code$args[[2]]
    if(!is.character(memberName)) stop(exprClassProcessingErrorMsg(code, 'In handling X$Y: Something is wrong with Y.'), call. = FALSE)

    memberSymbolObject <- nfProc$getSymbolTable()$getSymbolObject(memberName)
    if(!is.null(memberSymbolObject)) code$type <- memberSymbolObject$type
    
    if(isSymList | isSymFunc) {
        ## nimbleList
        ## We need (*nl) in C++, represented by cppPointerDereference(nl)
        if(code$args[[1]]$name != 'cppPointerDereference') {
            a1 <- insertExprClassLayer(code, 1, 'cppPointerDereference',
                                                type = code$args[[1]]$type,
                                                nDim = code$args[[1]]$nDim,
                                                sizeExprs = code$args[[1]]$sizeExprs)
        }
    }
     
    ## following checks are on type of A$B (isSymList and isSymFunc refer to type of A)
    
    if(code$type == 'nimbleList') {
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

sizeNimDerivs <- function(code, symTab, typeEnv){
  code$name <- "nimDerivs_dummy"
  ## static <- code$args[['static']]
  ## code$args[['calcNodes']] <- NULL
  ## code$args[['static']] <- NULL ## Ok since these two are last arguments.  Otherwise we need to shift args.
  updateNodesName <- code$args[['updateNodesName']]
  code$args[['updateNodesName']] <- NULL

  ## next section is adapted from sizeNimbleListReturningFunction
  ## skip first argument, which is the call for which derivatives are requested.
  typeEnv$.avoidAD <- TRUE # tells sizeConcatenate to protect lifted vectors from AD. Also tells sizeInsertIntermediate to tag intermediates as non-AD
  on.exit({typeEnv$.avoidAD <- NULL})
  asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, rep(TRUE, length(code$args)-1)))
  code$type <- 'nimbleList'
  nlGen <- nimbleListReturningFunctionList[[code$name]]$nlGen
  nlDef <- nl.getListDef(nlGen)
  className <- nlDef$className
  symbolObject <- symTab$getSymbolObject(className, inherits = TRUE)
  if(is.null(symbolObject)) {
    nlp <- typeEnv$.nimbleProject$compileNimbleList(nlGen, initialTypeInference = TRUE)
    symbolObject <- symbolNimbleListGenerator(name = className, nlProc = nlp)
    symTab$addSymbol(symbolObject)
  }
  code$sizeExprs <- symbolObject
  code$toEigenize <- "no"  # This is specialized for nimSvd and nimEigen.
  code$nDim <- 0
  
  ## asserts <- sizeNimbleListReturningFunction(code, symTab, typeEnv)
  ## lift wrt if needed.  I'm not sure why sizeNimbleListReturningFunction doesn't handle lifting
  if(inherits(code$args[['wrt']], 'exprClass')) {
    if(!code$args[['wrt']]$isName) {
      iWrt <- which(names(code$args) == 'wrt')
      if(length(iWrt) != 1) stop("problem working on wrt argument to nimDerivs")
      asserts <- c(asserts, sizeInsertIntermediate(code, iWrt, symTab, typeEnv) )
    }
  }
  if(inherits(code$args[['order']], 'exprClass')) {
    if(!code$args[['order']]$isName) {
      iOrder <- which(names(code$args) == 'order')
      if(length(iOrder) != 1) stop("problem working on order argument to nimDerivs")
      newAssert <- sizeInsertIntermediate(code, iOrder, symTab, typeEnv)
      for(iNewA in seq_along(newAssert)) {
          if(inherits(newAssert[[iNewA]], "exprClass")) {
              if(is.null(newAssert[[iNewA]]$aux))
                  newAssert[[iNewA]]$aux <- list()
              newAssert[[iNewA]]$aux$.avoidAD <- TRUE
          }
      }
      asserts <- c(asserts, newAssert)
      ## asserts <- c(asserts, sizeInsertIntermediate(code, iOrder, symTab, typeEnv) )
    }
  }
  
  insertExprClassLayer(code, which(names(code$args)=='wrt'), 'make_vector_if_necessary',
                       type = 'double',
                       nDim = 1,
                       sizeExprs = list())
  
  a1 <- insertExprClassLayer(code, which(names(code$args)=='order'), 'make_vector_if_necessary',
                             type = 'double',
                             nDim = 1,
                             sizeExprs = list())
  newADinfoName <- ADinfoLabel()
##  symTab$addSymbol(symbolADinfo$new(name = newADinfoName))
  if(!is.list(code$aux))
    code$aux <- list()

  code$aux[['ADinfoName']] <- newADinfoName
  if(!is.null(updateNodesName))
    code$aux[['updateNodesName']] <- updateNodesName
  
  ADinfoNames <- 'ADinfoNames'
  ##  ADinfoNames <- if(isTRUE(static)) 'ADstaticInfoNames' else 'ADinfoNames'
  if(is.null(typeEnv[[ADinfoNames]])) {
    typeEnv[[ADinfoNames]] <- newADinfoName
  } else {
    typeEnv[[ADinfoNames]] <- c(typeEnv[[ADinfoNames]],
                                newADinfoName)
  }
  if(!nimbleOptions('experimentalSelfLiftStage')) {
    if(!(code$caller$name %in% assignmentOperators))
      asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
  }
  if(length(asserts) == 0) NULL else asserts
}

sizeNimDerivsCalculate <- function(code, symTab, typeEnv){
  ## this contains some logic from sizeNimbleListReturningFunction,
  ## and some from sizeScalarModelOp.  For now, simplest to make a new size processor.
  nlGen <- nimbleListReturningFunctionList[[code$name]]$nlGen
  nlDef <- nl.getListDef(nlGen)
  className <- nlDef$className
  symbolObject <- symTab$getSymbolObject(className, inherits = TRUE)
  if(is.null(symbolObject)) {
    nlp <- typeEnv$.nimbleProject$compileNimbleList(nlGen, initialTypeInference = TRUE)
    symbolObject <- symbolNimbleListGenerator(name = className, nlProc = nlp)
    symTab$addSymbol(symbolObject)
  }
    
  newADinfoName <- ADinfoLabel()
  if(!is.list(code$aux))
    code$aux <- list()

  code$aux[['ADinfoName']] <- newADinfoName
  
  ADinfoNames <- 'ADinfoNames_calculate'
  if(is.null(typeEnv[[ADinfoNames]])) {
    typeEnv[[ADinfoNames]] <- newADinfoName
  } else {
    typeEnv[[ADinfoNames]] <- c(typeEnv[[ADinfoNames]],
                                newADinfoName)
  }

  ## if(is.null(typeEnv[['numNimDerivsCalculate']])) { ## running count of nimDerivs_calculate calls in the nimbleFunction
  ##   typeEnv[['numNimDerivsCalculate']] <- 1
  ## } else {
  ##   typeEnv[['numNimDerivsCalculate']] <- typeEnv[['numNimDerivsCalculate']] + 1
  ## }
  code$sizeExprs <- symbolObject
  code$type <- 'nimbleList'
  code$nDim <- 0
  code$toEigenize <- 'maybe'
  asserts <- recurseSetSizes(code, symTab, typeEnv) # useArgs = c(FALSE, rep(TRUE, code$args)))
  for(i in 2:length(code$args)) {
    if(inherits(code$args[[i]], 'exprClass')) {
      if(code$args[[i]]$toEigenize=='yes')
        asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv))
    }
  }
  if(inherits(code$args[[2]], 'exprClass') && names(code$args)[2] != 'orderVector') { ## There is an index expression that may be non-scalar
    if(code$args[[2]]$nDim > 0) { ## It is non-scalar so we need to set a logical argument about whether is it a logical or numeric vector
      code$args[[ length(code$args)+1 ]] <- as.integer(code$args[[2]]$type == 'logical')
    }
  }
  insertExprClassLayer(code, which(names(code$args)=='orderVector'), 'make_vector_if_necessary',
                       type = 'double',
                       nDim = 1,
                       sizeExprs = list())  
  if(code$args[[1]]$toEigenize == 'yes') { ## not sure when this would be TRUE
    asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
  }
  if(!nimbleOptions('experimentalSelfLiftStage')) {
    if(!(code$caller$name %in% assignmentOperators))
      asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
  }
  if(length(asserts) == 0) NULL else asserts
}

sizeNimbleListReturningFunction <- function(code, symTab, typeEnv) {
  asserts <- recurseSetSizes(code, symTab, typeEnv)
  code$type <- 'nimbleList'
  nlGen <- nimbleListReturningFunctionList[[code$name]]$nlGen
  nlDef <- nl.getListDef(nlGen)
  className <- nlDef$className
  symbolObject <- symTab$getSymbolObject(className, inherits = TRUE)
  if(is.null(symbolObject)) {
      nlp <- typeEnv$.nimbleProject$compileNimbleList(nlGen, initialTypeInference = TRUE)
      symbolObject <- symbolNimbleListGenerator(name = className, nlProc = nlp)
      symTab$addSymbol(symbolObject)
  }
  code$sizeExprs <- symbolObject
  code$toEigenize <- "yes"  # This is specialized for nimSvd and nimEigen.
  if(code$name == 'getDerivs_wrapper'){
      code$toEigenize <- 'no'  ## Temp. solution to ensure that derivsOrders argument is a nimArray and not an eigen type.
  }
  code$nDim <- 0
  if(!nimbleOptions('experimentalSelfLiftStage')) {
      if(!(code$caller$name %in% assignmentOperators))
          asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
  }
  if(length(asserts) == 0) NULL else asserts
}

## sizeOptimModel <- function(code, symTab, typeEnv) {
##     ## No call to recurseSetSizes.
##     code$type <- 'nimbleList'
##     nlGen <- nimbleListReturningFunctionList[[code$name]]$nlGen
##     nlDef <- nl.getListDef(nlGen)
##     className <- nlDef$className
##     symbolObject <- symTab$getSymbolObject(className, inherits = TRUE)
##     if(is.null(symbolObject)) {
##         nlp <- typeEnv$.nimbleProject$compileNimbleList(nlGen, initialTypeInference = TRUE)
##         symbolObject <- symbolNimbleListGenerator(name = className, nlProc = nlp)
##         symTab$addSymbol(symbolObject)
##     }
##     code$sizeExprs <- symbolObject
##     code$toEigenize <- "no"
##     code$nDim <- 0
##     list()
## }

sizeOptim <- function(code, symTab, typeEnv) {
    typeEnv$.allowFunctionAsArgument <- TRUE
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    typeEnv$.allowFunctionAsArgument <- FALSE
    code$type <- 'nimbleList'
    nlGen <- nimbleListReturningFunctionList[[code$name]]$nlGen
    nlDef <- nl.getListDef(nlGen)
    className <- nlDef$className
    symbolObject <- symTab$getSymbolObject(className, inherits = TRUE)
    if(is.null(symbolObject)) {
        nlp <- typeEnv$.nimbleProject$compileNimbleList(nlGen, initialTypeInference = TRUE)
        symbolObject <- symbolNimbleListGenerator(name = className, nlProc = nlp)
        symTab$addSymbol(symbolObject)
    }
    code$sizeExprs <- symbolObject
    code$toEigenize <- "no"
    code$nDim <- 0

    if(inherits(code$args[[1]], 'exprClass')) {
        if(!(code$args[[1]]$isName))
            asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
    }
    
    fnCode <- code$args$fn
    if(!inherits(fnCode, 'exprClass')) {
        stop(exprClassProcessingErrorMsg(code, '`fn` argument to `optim` is not valid.'), call. = FALSE)
    }
    if (fnCode$name == 'nfMethod') {
        # This is handled in cppOutputNFmethod.
    } else if(identical(fnCode$type, 'Ronly') & identical(class(fnCode$sizeExprs)[1], 'symbolMemberFunction')) {
        fnCode$name <- fnCode$sizeExprs$RCfunProc$name
        newCode <- substitute(nfMethod(this, FUN), list(FUN = fnCode$name))
        newExpr <- RparseTree2ExprClasses(newCode)
        newExpr$args[[1]]$type <- symTab$getSymbolObject(".self", TRUE)$baseType
        setArg(code, 2, newExpr)
    } else if(exists(fnCode$name) && is.rcf(get(fnCode$name))) {
        # Handle fn arguments that are RCfunctions.
        fnCode$name <- environment(get(fnCode$name))$nfMethodRCobject$uniqueName
    } else {
        stop('in `optim`, the `fn` argument, `', fnCode$name, '`, is not available or is not a nimbleFunction or nimbleFunction method.')
    }

    grCode <- code$args$gr
    if (identical(grCode, "NULL")) {
        # We simply emit "NULL".
    } else if (grCode$name == 'nfMethod') {
        # This is handled in cppOutputNFmethod.
    } else if(identical(grCode$type, 'Ronly') & identical(class(grCode$sizeExprs)[1], 'symbolMemberFunction')) {
        grCode$name <- grCode$sizeExprs$RCfunProc$name
        newCode <- substitute(nfMethod(this, FUN), list(FUN = grCode$name))
        newExpr <- RparseTree2ExprClasses(newCode)
        newExpr$args[[1]]$type <- symTab$getSymbolObject(".self", TRUE)$baseType
        setArg(code, 3, newExpr)
    } else if(exists(grCode$name) && is.rcf(get(grCode$name))) {
        # Handle gr arguments that are RCfunctions.
        grCode$name <- environment(get(grCode$name))$nfMethodRCobject$uniqueName
    } else {
        stop(paste0('unsupported gr argument in optim(par, gr = ', grCode$name, '); try an RCfunction or nfMethod instead'))
    }

    for(arg in c(code$args$lower, code$args$upper)) {
        if(inherits(arg, 'exprClass') && arg$toEigenize=='yes') {
            asserts <- c(asserts, sizeInsertIntermediate(code, arg$callerArgID, symTab, typeEnv))
        }
    }

    if(length(asserts) == 0) NULL else asserts
}

sizeOptimDefaultControl <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    code$type <- 'nimbleList'
    nlGen <- nimbleListReturningFunctionList[[code$name]]$nlGen
    nlDef <- nl.getListDef(nlGen)
    className <- nlDef$className
    symbolObject <- symTab$getSymbolObject(className, inherits = TRUE)
    if(is.null(symbolObject)) {
        nlp <- typeEnv$.nimbleProject$compileNimbleList(nlGen, initialTypeInference = TRUE)
        symbolObject <- symbolNimbleListGenerator(name = className, nlProc = nlp)
        symTab$addSymbol(symbolObject)
    }
    code$sizeExprs <- symbolObject
    code$toEigenize <- "no"
    code$nDim <- 0

    if(length(asserts) == 0) NULL else asserts
}

sizeCppPointerDereference <- function(code, symTab, typeEnv) {
  asserts <- recurseSetSizes(code, symTab, typeEnv)
  code$type <- code$args[[1]]$type
  code$sizeExprs <- code$args[[1]]$sizeExprs
  code$toEigenize <- code$args[[1]]$toEigenize
  code$nDim <- code$args[[1]]$nDim
  if(length(asserts) == 0) NULL else asserts
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

    nfProc <- symbolObject$nfProc
    if(is.null(nfProc)) {
        stop(exprClassProcessingErrorMsg(code, 'In handling X$Y(): Symbols for X have not been set up.'), call. = FALSE)
    }
    if(isSymFun) {
        
        if(a1$args[[1]]$name != 'cppPointerDereference') {
            insertExprClassLayer(a1, 1, 'cppPointerDereference') ## not annotated, but not needed
        }

    }

    if(isSymFun) {
        returnSymbol <- nfProc$compileInfos[[methodName]]$returnSymbol
        argSymTab <- nfProc$compileInfos[[methodName]]$origLocalSymTab
    } 
    if(isSymFunList) {
        returnSymbol <- nfProc$compileInfos[[methodName]]$returnSymbol
        argSymTab <- nfProc$compileInfos[[methodName]]$origLocalSymTab
    }
    if(!is.null(returnSymbol)) {
        asserts <- generalFunSizeHandlerFromSymbols(code, symTab, typeEnv, returnSymbol, argSymTab, chainedCall = TRUE)
        return(asserts)
    }
    writeLines('Warning: no return type determined from a chained call such as nimbleFunction call.')
    invisible(NULL)
}

sizeValues <- function(code, symTab, typeEnv) {
    code$nDim <- 1
    code$type <- 'double'
    code$toEigenize <- 'no'
    sym <- symTab$getSymbolObject(code$args[[1]]$name, TRUE)
    indexRangeCase <- FALSE
    if(length(code$args) == 1) {  # full vector of nodes
        code$sizeExprs <- list(substitute(cppMemberFunction(getTotalLength(ACCESSNAME)), list(ACCESSNAME = as.name(code$args[[1]]$name))))
        asserts <- list()
    } else {  # there must be index on the node
        asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, rep(TRUE, length(code$args)-1)))
        if(is.numeric(code$args[[2]])) {
            code$sizeExprs <- list(substitute(cppMemberFunction(getNodeLength(ACCESSNAME, ACCESSINDEX)), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = code$args[[2]])))
        } else {
            if(!(code$args[[2]]$isName))
                asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))

            if(code$args[[2]]$nDim > 0) {
                code$sizeExprs <- list(substitute(getNodesLength_Indices(ACCESSNAME, ACCESSINDEX), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = as.name(code$args[[2]]$name))))
                indexRangeCase <- TRUE
            } else {
                code$sizeExprs <- list(substitute(cppMemberFunction(getNodeLength(ACCESSNAME, ACCESSINDEX)), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = as.name(code$args[[2]]$name))))
            }
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
                        if(code$args[[2]]$nDim > 0) {
                            assertSS[[1]][[3]] <- substitute(getNodesLength_Indices(ACCESSNAME, ACCESSINDEX), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = as.name(code$args[[2]]$name))) ## intermediate has already been inserted above, if needed
                        } else {
                            assertSS[[1]][[3]] <- substitute(cppMemberFunction(getNodeLength(ACCESSNAME, ACCESSINDEX)), list(ACCESSNAME = as.name(code$args[[1]]$name), ACCESSINDEX = as.name(code$args[[2]]$name)))
                        }
                    }
                }

                asserts <- c(asserts, assertSS)
            } else
                typeEnv$.ensureNimbleBlocks <- TRUE
         } else {   # values(...) <- P, don't change it
             if(indexRangeCase) code$name <- 'valuesIndexRange'
         }
    } else { ## values(...) embedded in a RHS expression
        code$name <- if(!indexRangeCase) 'getValues' else 'getValuesIndexRange'
        code$toEigenize <- 'yes' ## This tricks sizeAssignAfterRecursing to generate the setSize in asserts, in getValues case (getValuesIndexRange is in set of names to skip for that)
        asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        code$toEigenize <- 'no'
    }
    if(length(asserts)==0) NULL else asserts
}

sizeRCfunction <- function(code, symTab, typeEnv, nfmObj, RCfunProc) {
    returnType <- nfmObj$returnType
    ## argInfo <- nfmObj$argInfo
    ## Insert buildDerivs label into code$aux
    if(is.list(nfmObj$buildDerivs)) {
        if(is.null(code$aux))
            code$aux <- list(buildDerivs = TRUE)
        else
            code$aux[['buildDerivs']] <- TRUE
    }
    code$name <- nfmObj$uniqueName
    returnSymbol <- RCfunProc$compileInfo$returnSymbol
    argSymTab <- RCfunProc$compileInfo$origLocalSymTab
    asserts <- generalFunSizeHandlerFromSymbols(code, symTab, typeEnv, returnSymbol, argSymTab)
    return(asserts)
}

sizeNimbleFunction <- function(code, symTab, typeEnv) { ## This will handle other nimbleFunction run calls or other methods of this nimbleFunction
    sym <- symTab$getSymbolObject(code$name, TRUE)
    ok <- FALSE
    if(inherits(sym, 'symbolNimbleFunction')) {
        stop(exprClassProcessingErrorMsg(code, 'In sizeNimbleFunction: A nimbleFunction method should not be processed here.'), call. = FALSE)
        ## HANDLING OF myNF$run() HERE IS DEFUNCT.  ALL SHOULD GO THROUGH sizeChainedCall now (chainedCall(nfMethod(myNF,'run'), arg1, arg2).
    }
    if(inherits(sym, 'symbolMemberFunction')) {
        memberRCfunProc <- sym$RCfunProc
        returnSymbol <- memberRCfunProc$compileInfo$returnSymbol
        argSymTab <- memberRCfunProc$compileInfo$origLocalSymTab
        ok <- TRUE
    }
    if(ok) {
        asserts <- generalFunSizeHandlerFromSymbols(code, symTab, typeEnv, returnSymbol, argSymTab)
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
## toEigen: N.B. This may be deprecated.
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
    toLift <- lapply(code$args, function(x) {
        if(inherits(x, 'exprClass')) (identical(x$type, 'logical') & !x$isName) else FALSE
    })
    for(i in seq_along(toLift)) {
        if(toLift[[i]])
            asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv)) 
    }
    code$toEigenize <- if(any( unlist(toEigs) %in% c('maybe', 'yes'))) 'yes' else 'no'
    code$type <- 'unknown'
    if(length(asserts) == 0) NULL else asserts
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
      if(sym$type == 'nimbleList'){
        sym <- sym$nlProc$symTab$getSymbolObject(code$args[[1]]$args[[2]])
      }
    } else {
      sym <- symTab$getSymbolObject(code$args[[1]]$name, inherits = TRUE)
      useArg1 <- FALSE
    } 
    asserts <- list()
    
    if(!inherits(sym, 'symbolNumericList')) {
        if(sym$nDim == 0) stop(exprClassProcessingErrorMsg(code, 'In sizeSetSize: Resizing a scalar does not make sense.'), call. = FALSE)
        firstSizeExpr <- code$args[[2]]

        ## first two arguments are variable to be resized and new sizes
        ## extra arguments would be fillZeros and recycle
        ## need to determine if any extra arguments were provided in order to repack arguments correctly below
        if(length(code$args) > 2)
            nExtraArgs <- length(code$args)-2
        else
            nExtraArgs <- 0

        if(nExtraArgs > 0)
            asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, c(rep(FALSE, 2), rep(TRUE, nExtraArgs))))

        if(inherits(firstSizeExpr, 'exprClass')) {
            if(firstSizeExpr$name == 'nimC') { ## handle syntax of resize(Z, c(3, dim(A)[1]))
                if(length(firstSizeExpr$args) != sym$nDim) stop(exprClassProcessingErrorMsg(code, 'In sizeSetSize: Problem with number of dimensions provided in resize.'), call. = FALSE)
                asserts <- c(asserts, recurseSetSizes(firstSizeExpr, symTab, typeEnv)) ## may set intermediates if needed
                if(nExtraArgs > 0) {
                    origExtraArgs <- code$args[3:length(code$args)] ## preserve extra arguments
                    code$args <- code$args[1:2]
                }
                for(i in 1:length(firstSizeExpr$args)) {
                    code$args[[i+1]] <- firstSizeExpr$args[[i]]
                    if(inherits(firstSizeExpr$args[[i]], 'exprClass')) {
                        firstSizeExpr$args[[i]]$caller <- code
                        firstSizeExpr$args[[i]]$callerArgID <- i+1
                    }
                }
                if(nExtraArgs > 0) { ## reinsert extra arguments on end.
                    for(i in 1:nExtraArgs) {
                        setArg(code, length(code$args) + 1, origExtraArgs[[i]])
                    }
                }
                return(if(length(asserts)==0) NULL else asserts)
            }
        }

        useArgs <- c(useArg1, TRUE )
        if(nExtraArgs > 0) useArgs <- c(useArgs, rep(FALSE, nExtraArgs))
        asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, useArgs) )

        if(inherits(code$args[[2]], 'exprClass')) {
            if(code$args[[2]]$nDim > 0) {
                 if(!(code$args[[2]]$isName)) asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
                code$name <- 'setSizeNimArrToNimArr'
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
        invisible(NULL)
    }
}


## This was redundant and we should eventually be able to remove it
## toEigen: N.B. omitting this
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

sizeInsertIntermediate <- function(code, argID, symTab, typeEnv, forceAssign = FALSE,
                                   forceType = NULL) { ## only used by sizeColonOperator to force lifted range ends to be integer
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
        if(!nimbleOptions('experimentalNewSizeProcessing')) {
            newArgExpr$toEigenize <- 'maybe'
        }
        newArgExpr$nDim <- code$args[[argID]]$nDim
    } else {

        ## One may wonder where the new variable is added to the
        ## symbolTable.  That happens when we do
        ## sizeAssignAfterRecursing, which identifies unknown LHS and
        ## creates the symTab entry.
        
        newExpr <- newAssignmentExpression()
        setArg(newExpr, 1, RparseTree2ExprClasses(as.name(newName))) 
        setArg(newExpr, 2, code$args[[argID]]) ## The setArg function should set code$caller (to newExpr) and code$callerArgID (to 3)
        if(!is.null(forceType))
            newExpr$args[[2]]$type <- forceType
        ans <- c(sizeAssignAfterRecursing(newExpr, symTab, typeEnv, NoEigenizeMap = TRUE), list(newExpr))

        newArgExpr <- RparseTree2ExprClasses(as.name(newName))
        newArgExpr$type <- newExpr$args[[1]]$type
        newArgExpr$sizeExprs <- newExpr$args[[1]]$sizeExprs
        if(!nimbleOptions('experimentalNewSizeProcessing')) {
            newArgExpr$toEigenize <- 'maybe'
        }
        newArgExpr$nDim <- newExpr$args[[1]]$nDim
    }
    setArg(code, argID, newArgExpr)
    if( isTRUE(typeEnv$.avoidAD) ) {
        typeEnv$.new_ignore <- c(typeEnv$.new_ignore, newName)
    }
    return(ans) ## This is to be inserted in a list of asserts, even though it is really core code, not just an a test or assertion
}

nimbleAliasRiskFxns <- c("t", "[", "eigenBlock", ## all "[" will be replaced by eigenBlock anyway
                         names(nimble:::sizeCalls)[ grepl("RecyclingRule", unlist(nimble:::sizeCalls) ) ],
                         "nimRep", "nimRepd", "nimRepi", "nimRepb", "nimC", "nimCd", "nimCi", "nimCb")

detectNimbleAliasRisk <- function(code, LHSname, insideRiskFxn = FALSE) {
    if(!inherits(code, "exprClass")) return(FALSE)
    if(!insideRiskFxn) {
        if(code$name %in% nimbleAliasRiskFxns)
            insideRiskFxn <- TRUE
    }
    if(insideRiskFxn)
        if(code$name == LHSname)
            return(TRUE)
    if(length(code$args) > 0) {
        for(i in seq_along(code$args)) {
            if(detectNimbleAliasRisk(code$args[[i]], LHSname, insideRiskFxn))
                return(TRUE)
        }
    }
    FALSE
}

# Inspect an expression to extract any nfVar(A, x), meaning A$x
# This passes over '[' at the top level, so that A$x[i] goes to A$x
find_nfVar_info <- function(code) {
  if(!inherits(code, 'exprClass')) return(NULL)
  if(code$name == "[") return(find_nfVar_info(code$args[[1]]))
  if(code$name == "nfVar") { # using A$x which is nfVar(A, x)
    A <- code$args[[1]]$name
    x <- code$args[[2]]
    if(is.character(x)) {
      inner <- x
    } else {
      if(!(x$name == "nfVar"))
        return(list(var = NULL, objs = NULL))
      x <- find_nfVar_info(x)
      A <- c(A, x$objs)
      inner <- x$var
    }
    result <- list(var = inner,
                   objs = A )
    return(result)
  }
  NULL
}

detect_nfVar_info_aliasRisk <- function(code, LHS_nfVar_info) {
  if(!inherits(code, "exprClass")) return(FALSE)
  if(code$name == "nfVar") {
    this_nfVar_info <- find_nfVar_info(code)
    if(identical(this_nfVar_info$var, LHS_nfVar_info$var))
      return(TRUE)
    else
      return(FALSE)
  }
  if(length(code$args) > 0) {
    for(i in seq_along(code$args)) {
      if(detect_nfVar_info_aliasRisk(code$args[[i]], LHS_nfVar_info))
        return(TRUE)
    }
  }
  FALSE
}

sizeAssign <- function(code, symTab, typeEnv) {
    typeEnv$.AllowUnknowns <- FALSE
    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, TRUE))
    typeEnv$.AllowUnknowns <- TRUE
    if(length(code$args) > 2){
      asserts <- c(asserts, exprClasses_setSizes(code, symTab, typeEnv))
    }
    else{
      asserts <- c(asserts, recurseSetSizes(code, symTab, typeEnv, useArgs = c(TRUE, FALSE)))
      typeEnv[['.ensureNimbleBlocks']] <- FALSE ## may have been true from RHS of rmnorm etc.
      LHS <- code$args[[1]]
      RHS <- code$args[[2]]
      if(inherits(LHS, 'exprClass')) {
          if(LHS$isName) {
              if(detectNimbleAliasRisk(RHS, LHS$name)) {
                  asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
              }
          } else {
            LHS_nfVar_info <- find_nfVar_info(LHS)
            if(!is.null(LHS_nfVar_info)) {
              if(detect_nfVar_info_aliasRisk(RHS, LHS_nfVar_info)) {
                asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
              }
            }
          }
      }
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
        if(is.numeric(RHS) || is.logical(RHS)) {
            RHSname = ''
            RHSnDim <- 0
            RHStype <- sizeProc_storage_mode(RHS)
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
    if(is.null(RHStype) || length(RHStype)==0) {
        if(startsWith(RHSname, "r") && gsub( "^r", "d", RHSname) %in% nimbleUserNamespace$distributions$namesVector)  # Fix issue 1355.
            stop("Missing simulation function '", RHSname, "', perhaps because it was deleted. Please use `deregisterDistributions` to deregister the distribution.")
        stop(exprClassProcessingErrorMsg(code, paste0("In sizeAssignAfterRecursing: '", RHSname, "' is not available or its output type is unknown.")), call. = FALSE)
    }
    if(LHS$isName) {
        if(!exists(LHS$name, envir = typeEnv, inherits = FALSE)) { ## not in typeEnv
            ## If LHS unknown, create it in typeEnv
            if(!symTab$symbolExists(LHS$name, TRUE)) { ## not in symTab
                if(RHStype %in% c('double','integer', 'logical')) {  ## valid type to create here
                    ## We used to delay creating sizeExprs until below, but now it always generic
                    assign(LHS$name, exprTypeInfoClass$new(nDim = RHSnDim, type = RHStype, sizeExprs = makeSizeExpressions(rep(NA, RHSnDim), LHS$name)), envir = typeEnv)
                    symTab$addSymbol(symbolBasic(name = LHS$name, nDim = RHSnDim, type = RHStype))
                } else { ## not valid type to create here
                    if(RHStype == 'voidPtr') {
                        ## This should be ok without sizeExprs content
                        assign(LHS$name, exprTypeInfoClass$new(nDim = RHSnDim, type = RHStype), envir = typeEnv)
                        symTab$addSymbol(symbolVoidPtr(name = LHS$name, type = RHStype))
                    }
                    ## a path for arbitrary symbols
                    else if(RHStype == "custom") {                        
                        ConlySym <- RHS$sizeExprs$copy() ## trick to put a symbol object here. use a copy in case this expr is from simple assignment, not creation
                        ConlySym$name <- LHS$name
                        symTab$addSymbol(ConlySym)

                        code$type <- "custom"
                        code$sizeExprs <- ConlySym ## in case there is chained assignment

                        return(invisible(NULL))
                    }
                    else if(RHStype == "nimbleList") {
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
                      symTab$addSymbol(symbolNimbleList(name = LHS$name, nlProc = LHSnlProc))
                    }
                    else
                        stop(exprClassProcessingErrorMsg(code, paste0('In sizeAssignAfterRecursing: LHS is not in typeEnv or symTab and cannot be added now.')), call. = FALSE)
                }
            } else { ## yes in symTab
                ## this is another path for arbitrary symbols, but not sure it's used.
                ## This case is ok.  It is in the symbol table but not the typeEnv.  So it is something like ptr <- getPtr(A)
                 if(!nimbleOptions('experimentalNewSizeProcessing')) {
                code$toEigenize <- 'no'
                 } ##experimentalNewSizeProcessing
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
                stop('Mismatched dimensions in assignment: ', nimDeparse(code), '.', call. = FALSE )
            }
            ## and warn if type issue e.g. int <- double
            if(assignmentTypeWarn(LHS$type, RHStype)) {
                message(paste0('Warning, RHS numeric type is losing information in assignment to LHS in line:\n', nimDeparse(code)))
            }
        }
    }
    ## update size info in typeEnv
    assert <- NULL

    if((LHS$name == 'values' || LHS$name == 'valuesIndexRange') && length(LHS$args) %in% c(1,2)) { ## It is values(model_values_accessor[, index]) <- STUFF
        # triggered when we have simple assignment into values() without indexing of values()
        if(is.numeric(RHS)) stop(exprClassProcessingErrorMsg(code, paste0('In sizeAssignAfterRecursing: Cannot assign into values() from numeric.')), call. = FALSE)
        code$name <- if(LHS$name == 'values') 'setValues' else 'setValuesIndexRange' 
        code$args <- list(1 + length(LHS$args))
        setArg(code, 1, RHS)
        setArg(code, 2, LHS$args[[1]])
        if(length(LHS$args) == 2) setArg(code, 3, LHS$args[[2]])  # for indexed of form values(model, nodes[i])
        if(!(RHS$isName)) assert <- c(assert, sizeInsertIntermediate(code, 1, symTab, typeEnv) )
        return( if(length(assert) == 0) NULL else assert )
    }

    ## Note this can use LHS$name for RHSsizeExprs when returning from a nimbleFunction on RHS.  But this is probably not needed any more.
    if(any(unlist(lapply(RHSsizeExprs, is.null)))) RHSsizeExprs <- makeSizeExpressions(rep(NA, RHSnDim), LHS$name) ## reset sizeExprs for the LHS var. re-using RHSsizeExprs for LHS.  This would only be valid if it is a nimbleFunction returning something on the RHS.  For assignment to be executed in Eigen, the RHS sizes MUST be known

     if(!nimbleOptions('experimentalNewSizeProcessing')) {
                
    if(LHS$toEigenize == 'yes') {
        code$toEigenize <- 'yes'
##        message('Warning from sizeAssign: not expecting LHS to have toEigenize == yes')
    } else {
        code$toEigenize <-if(inherits(RHS, 'exprClass')) {
            if(RHS$toEigenize == 'no') 'no' else {
                if(RHS$toEigenize == 'unknown') 'no' else {
                    if(RHS$toEigenize != 'yes' && (!(LHS$name %in% c('eigenBlock', 'diagonal', 'coeffSetter'))) && (RHS$isName || RHS$nDim == 0 || (RHS$name == 'map' && NoEigenizeMap))) 'no' ## if it is scalar or is just a name or a map, we will do it via NimArr operator= .  Used to have "| RHS$name == 'map'", but this allowed X[1:3] <- X[2:4], which requires eigen, with eval triggered, to get right
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
            if(!(RHS$name %in% setSizeNotNeededOperators)) {
                if(LHS$isName || LHS$name == "nfVar") {
                    assert <- substitute(setSize(LHS), list(LHS = nimbleGeneralParseDeparse(LHS)))
                    for(i in seq_along(RHSsizeExprs)) {
                        test <- try(assert[[i + 2]] <- RHS$sizeExprs[[i]])
                        if(inherits(test, 'try-error')) stop(paste0('In sizeAssignAfterRecursing: Error in assert[[i + 2]] <- RHS$sizeExprs[[i]] for i = ', i), call. = FALSE)
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
                                    if(!is.numeric(LHS$sizeExprs[[1]]) || !is.numeric(RHS$sizeExprs[[2]])) {
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
     } ##experimentalNewSizeProcessing
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

sizePROTECT <- function(code, symTab, typeEnv) {
    ## Do not recurse.
    code$type <- "custom"
    code$sizeExprs <- symbolSEXP(type = 'custom') ## trick to put a symbol object into sizeExprs for later use
    return(invisible(NULL))
}

sizeReval <- function(code, symTab, typeEnv) {
    code$name <- 'Rf_eval'
    return(sizePROTECT(code, symTab, typeEnv))
}

sizeNimbleConvert <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv) ## should not normally have an expression other than variable name as the argument, but do this for safety
    nDim <- code$args[[1]]$nDim
    type <- code$args[[1]]$type
    if(!code$caller$name %in% assignmentOperators) stop(exprClassProcessingErrorMsg(code, 'nimbleConvert can only be used in simple assignment.'), call. = FALSE)

    targetString <- nimDeparse(code$args[[1]])
    targetName <- Rname2CppName(targetString)
    targetExpr <- parse(text = targetString, keep.source = FALSE)[[1]]
    copyName <- paste0(targetName, '_nimbleContigCopy')
    subList <- list(var = targetExpr, copy = as.name(copyName))
    newCode <- substitute( nimArrPtr_copyIfNeeded(var, copy),
                             subList )
    ## only necessary if the result is needed
    if(!symTab$symbolExists( copyName )) {
        symTab$addSymbol(  symbolBasic(name = copyName, type = type, nDim = nDim) )
        assign(copyName, exprTypeInfoClass$new(nDim = nDim, type = type), envir = typeEnv)
    }
    newCode <- RparseTree2ExprClasses(newCode)
    newCode$type <- "custom"
    newCode$sizeExprs <- symbolPtr(type = type) ## trick to put a symbol object into sizeExprs for later use
    setArg(code$caller, code$callerArgID, newCode)
    
    asserts
}

sizeNimbleUnconvert <- function(code, symTab, typeEnv) {
    ptrString <- nimDeparse(code$args[[1]])
    ptrName <- Rname2CppName(ptrString)
    ptrExpr <- parse(text = ptrString, keep.source = FALSE)[[1]]

    targetString <- nimDeparse(code$args[[2]])
    targetName <- Rname2CppName(targetString)
    targetExpr <- parse(text = targetString, keep.source = FALSE)[[1]]

    copyName <- paste0(targetName, '_nimbleContigCopy')
    subList <- list(ptr = ptrExpr, var = targetExpr, copy = as.name(copyName))
    newCode <- substitute( nimArrPtr_copyBackIfNeeded(ptr, var, copy),
                             subList )

    newCode <- RparseTree2ExprClasses(newCode)
    setArg(code$caller, code$callerArgID, newCode)
    NULL
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

    ## length(model[[node]]) wasn't working because we were not doing recurseSetSize here
    ## However I am not sure if that is because there are cases where size expects a special argument we don't want to process (a modelValues?)
    ## So I'm going to wrap it in a try() and suppress messages
    asserts <- try(recurseSetSizes(code, symTab, typeEnv), silent = TRUE)
    if(inherits(asserts, 'try-error')) asserts <- list()
    if(code$args[[1]]$toEigenize == 'yes') {
        asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
    }
    code$nDim <- 0
    outputType <- scalarOutputTypes[[code$name]]
    if(is.null(outputType)) code$type <- 'double'
    else code$type <- outputType
    code$sizeExprs <- list()
    code$toEigenize <- 'maybe' ## a scalar can be eigenized or not
    asserts
}

sizeScalarModelOp <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) {
        asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs = c(FALSE, rep(TRUE, length(code$args)-1)))
        for(i in 2:length(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
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

sizeScalarRecurseAllowMaps <- function(code, symTab, typeEnv, recurse = TRUE) {
    ## use something different for distributionFuns

    ## Ensure that simple maps being passed will be passed without extra
    ## copy that would occur from lifting an Eigen expression.
    for(i in seq_along(code$args)) {
        ## check if each argument is purely of the form x[...]
        ## (note that x[...][...] might also be valid for passByMap
        ## but it is not handled that way currently.
        if(inherits(code$args[[i]], 'exprClass')) {
            if(code$args[[i]]$name == "[") {
                if(inherits(code$args[[i]]$args[[1]],
                            'exprClass')) { ## must be true, but I'm being defensive
                    if(code$args[[i]]$args[[1]]$isName) {
                        insertExprClassLayer(code, i, 'passByMap')
                    }
                }
            }
        }
    }

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
    ## This is for X[i, j], viewed as `[`(X, i, j), where there may be different numbers of indices, and they may be scalars, sequences defined by `:`, or arbitrary (nonSequence) vectors of integers or logicals.
    ## X itself could be Y[k, l] (or the result of processing it) or map(Y, k, l), which is created if Y is a model variable and we know we need a map into but at the point it is created there is no processing of how it should be represented, so it is just represented as an abstract map.

    ## recurse into arguments
    asserts <- recurseSetSizes(code, symTab, typeEnv)

    ## Check two special cases
    ## This is from modelValues:
    if(code$args[[1]]$type == 'symbolVecNimArrPtr') return(c(asserts, sizemvAccessBracket(code, symTab, typeEnv)))
    ## This is deprecated:
    if(code$args[[1]]$type == 'symbolNumericList') return(c(asserts, sizemvAccessBracket(code, symTab, typeEnv)))

    ## Iterate over arguments,  lifting any logical indices into which()
    ## e.g. X[i, bool] becomes X[i, Interm1], with Interm1 <- which(bool) as an assert.
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
                }
        }
    }

    ## Collect information about the number of dimensions and value of a drop argument if provided

    ## nDimVar is nDim of X
    nDimVar <- code$args[[1]]$nDim
    
    dropBool <- TRUE
    dropArgProvided <- FALSE
    if(!is.null(names(code$args)))
        if('drop' %in% names(code$args)) {
            dropArgProvided <- TRUE
            iDropArg <- which(names(code$args) == 'drop')
        }

    if(is.null(nDimVar)) 
        stop(paste0("Error in '", nimDeparse(code), "'; has '", nimDeparse(code$args[[1]]), "' been created?"))
    if(nDimVar != length(code$args) - 1 - dropArgProvided) { ## check if number of indices is correct
        ## only valid case with fewer index arguments than source dimensions is matrix[indices], where matrix can be treated as a vector
        if(!( (nDimVar == 2) & (length(code$args) - dropArgProvided) == 1)) {
            msg <- paste0('Error, wrong number of indices provided for ', nimDeparse(code),'.')
            stop(exprClassProcessingErrorMsg(code, msg), call. = FALSE)
        }
    }

    ## pick out the drop argument and check if it is logical
    if(dropArgProvided) {
        dropBool <- code$args[[iDropArg]]
        if(!is.logical(dropBool)) {
            msg <- paste0(msg, "(A drop argument must be hard-coded as TRUE or FALSE, not given as a variable.)")
            stop(exprClassProcessingErrorMsg(code, msg), call. = FALSE)
        }
    }
    ## These initial annotations may change later
    code$nDim <- nDimVar
    code$type <- code$args[[1]]$type
    ## Initialize sizeExprs
    code$sizeExprs <- vector('list', length = nDimVar)
    ## (We could generate asserts here to ensure sub-indexing is within bounds)

    ## needMap will become TRUE below unless all indices are scalars
    needMap <- FALSE

    ## Track whether if all index ranges are defined by `:` or by scalar
    ## simpleBlockOK will be TRUE if all index vectors and sequential, defined by `:`
    simpleBlockOK <- TRUE
    iSizes <- 1
    ## Iternate over dimensions of X and see which dimensions will be dropped from X[i,j,k] due to scalar indices, if drop = TRUE
    for(i in 1:nDimVar) {
        dropThisDim <- FALSE

        ## If the index is numeric, drop this dimension
        if(is.numeric(code$args[[i+1]])) dropThisDim <- TRUE
        ## If the index is not numeric but it is not a blank and its sizeExprs reveal it is a scalar-equivalent, drop this dimension
        else if((code$args[[i+1]]$name != "") & (length(dropSingleSizes(code$args[[i+1]]$sizeExprs)$sizeExprs) == 0)) dropThisDim <- TRUE

        ## Is this indices an expression?
        isExprClass <- inherits(code$args[[i+1]], 'exprClass') ## 

        if(dropThisDim) { ## The index is a scalar
            if(nimbleOptions()$indexDrop & dropBool) {  ## And flags allow dropping
                code$sizeExprs[[iSizes]] <- NULL        ## Remove that sizeExpr element
                code$nDim <- code$nDim - 1              ## reduce dimensions of result by 1
            } else { 
                code$sizeExprs[[iSizes]] <- 1; iSizes <- iSizes + 1  ## If we are not droping dimensions, set sizeExpr to 1
            }
            next
        } else {        ## not dropping a dimension, so the index is non-scalar
            if(isExprClass) ## If it is an expression that is not `:` or blank, then a simple block is not allowed
                if((code$args[[i+1]]$name != ':') && (code$args[[i+1]]$name != "")) simpleBlockOK <- FALSE
        }
        needMap <- TRUE ## If the "next" in if(dropThisDim) {} is always hit, then needMap will never be set to TRUE

        ## Update sizeExprs
        if(isExprClass) {
            if(code$args[[i+1]]$name != "") {
                ## An entry that is a variable possibly with a length
                code$sizeExprs[[iSizes]] <- code$args[[i+1]]$sizeExprs[[1]]
            } else {
                ## blank entry (e.g. A[,i]) is an exprClass with isName = TRUE and name = ""
                code$sizeExprs[[iSizes]] <- code$args[[1]]$sizeExprs[[i]]
                ## also at this point we will fill in a `:` expression for the indices,
                ## so now we have e.g. A[ 1:dim(A)[1], i ]
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

        ## We need to check whether X is an expression that needs to be lifted, say (A + B)[2, 3]
        ## We could do better for these cases 
        if(!code$args[[1]]$isName) ## It's not a name
            if(!(code$args[[1]]$name %in% operatorsAllowedBeforeIndexBracketsWithoutLifting)) {## e.g. 'mvAccessRow'
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

        ## for nested blocking, we have (nonseq | eigenBlock | map) x (nonseq | eigenBlock | map)
        ##      [ coeffSetter is a version of nonseq ]
        ## where nonseq means non-sequential indices, eigenBlock means sequential indices, and map means a model or modelValues variable marked abstractly for a map
        
        ## (map) x (eigenBlock | map) is already handled
        ## (map) x (nonseq)           is already handled
        ## (eigenBlock) x (eigenBlock) is already handled
        ## 
        ## check whether to nest the indexing directly

        ## nestIndexing TRUE means we will convert X[i, j][k, l] into X[ i[k], j[l] ] (while we are working on `[`(X[i, j], k, l)
        ## We do this for nested indexing except (eigenBlock) x (eigenBlock), which means all indices are sequential
        ## Then we just generate .block(..).block(..)
        nestIndexing <- FALSE

        ## code$args[[1]] is the X[i, j]
        if(!code$args[[1]]$isName) { ## In X[i], X is an expression
            ## X is an indexing expression of some kind (other than a map, which is already a new object)
            ## It can't be coeffSetter at this point in processing flow, because the nestedness implies its caller was not <-
            if(code$args[[1]]$name %in% c('eigenBlock', 'nimNonseqIndexedd' ,'nimNonseqIndexedi' ,'nimNonseqIndexedb' )) {
                ## if it is not (eigenBlock) x (eigenBlock)
                if(!( (code$args[[1]]$name == 'eigenBlock') & (simpleBlockOK)))
                    nestIndexing <- TRUE
            }
        }

        ## implement nestIndexing
        if(nestIndexing) {
            ## We have something like `[`( eigenBlock(X, i, j), k, l) or `[`( nimNonseqIndexedd(X, i, j), k, l)
            ## We will gradually take over the first argument to construct something that will end up like nimNonseqIndexedd(X, nimNonseqIndexedi(i, k), nimNonseqIndexedi(j, l))
            
            ## The first one was an eigenBlock (all sequential integer indices defined by `:` or blank imputed with `:`)
            if(code$args[[1]]$name == 'eigenBlock') {
                ## We have `[`( eigenBlock(X, i, j), k, l)
                ##
                ## reach down to X and rename it
                ## put `:`(start, finish) back together.
                ##
                ## If we are in  `[`( eigenBlock(X, i, j), k, l) <- Z,
                ## convert to  `[`( coeffSetter(X, i, j), k, l) <- Z,
                if(code$caller$name %in% assignmentOperators & code$callerArgID == 1) {
                    code$args[[1]]$name <- 'coeffSetter'
                } else {
                    ## otherwise, convert  `[`( eigenBlock(X, i, j), k, l) to `[`( nimNonseqIndexedd(X, i, j), k, l), e.g.
                    if(code$type == 'double') code$args[[1]]$name <- 'nimNonseqIndexedd'
                    if(code$type == 'integer') code$args[[1]]$name <- 'nimNonseqIndexedi'
                    if(code$type == 'logical') code$args[[1]]$name <- 'nimNonseqIndexedb'
                }
            } else { ## The first one was a nonSeq
                ## it was already nonseq, but it might need to become coeffSetter
                ## If we are in `[`( nimNonseqIndexedd(X, i, j), k, l) <- Z,
                ## convert to `[`( coeffSetter(X, i, j), k, l) <- Z,
                if(code$caller$name %in% assignmentOperators & code$callerArgID == 1) {
                    code$args[[1]]$name <- 'coeffSetter'
                }
            }
            ## Now construct the nesting i[k], j[l], etc.
            nestedNinds <- length(code$args[[1]]$args)-1
            nestedNdim <- code$args[[1]]$nDim
            nestedDropBool <- TRUE
            nestedDropArgProvided <- FALSE
            if(!is.null(names(code$args[[1]]$args))) ## does nimNonseqIndexedd(X, i, j) or coeffSetter(X, i, j) have named arguments?
                if("drop" %in% names(code$args[[1]]$args)) { ## is drop among the names?
                    nestedDropArgProvided <- TRUE
                    nestedDropBool <- code$args[[1]]$args[[ which(names(code$args[[1]]$args) == 'drop') ]]
                    nestedNinds <- nestedNinds - 1
                }
            nestedBlockBool <- rep(TRUE, nestedNinds)    ## is it preserved as a block (can still be scalar if nestedDropBool is FALSE)
            nestedScalarIndex <- rep(FALSE, nestedNinds) 

            ## Of the indices of nimNonseqIndexedd(X, i, j) or coeffSetter(X, i, j)
            ## which are scalars, and which are blocks
            ## If we have  nimNonseqIndexedd(X, i, j, drop = FALSE) or coeffSetter(X, i, j, drop = FALSE),
            ##    then we treat all dimensions as blocks, even if scalar indices
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

            ## Re-annotate first arg
            code$args[[1]]$sizeExprs <- code$sizeExprs
            code$args[[1]]$nDim <- code$nDim
            code$args[[1]]$type <- code$type
            numIndices <- length(code$args) - 1 - dropArgProvided

            ## Do we need to set drop carefully?
            
            ## NEED TO SKIP SCALARS IF dropBool = TRUE for nested case.
            nestedInds <- which(nestedBlockBool)
            if(length(nestedInds) != numIndices) stop(exprClassProcessingErrorMsg(code, 'Wrong number of nested indices.'), call.=FALSE)
            ## iterate over indices, constructing i[j] if necessary
            for(iInd in 1:numIndices) {
                nestedIind <- nestedInds[iInd]
                nestedIndexIsScalar <- if(inherits(code$args[[1]]$args[[nestedIind + 1]], 'exprClass')) code$args[[1]]$args[[nestedIind + 1]]$nDim == 0 else TRUE
                if(nestedIndexIsScalar) {
                    ## check:
                    ## In X[i, j][k, l], if i is scalar, k should also be scalar (can't check its value now, but should be 1 at run-time)
                    indexIsScalar <- if(inherits(code$args[[iInd+1]], 'exprClass')) code$args[[iInd+1]]$nDim == 0 else TRUE
                    if(!indexIsScalar) warning("There is nested indexing with drop=FALSE where an index must be scalar but isn't")
                } else {
                    ## construct i[k], which is really nimNonseqIndexedi(i, k)
                    newExpr <- exprClass$new(name = 'nimNonseqIndexedi', isName = FALSE, isCall = TRUE, isAssign = FALSE)
                    newExpr$type <- 'integer'
                    indexIsScalar <- if(inherits(code$args[[iInd+1]], 'exprClass')) code$args[[iInd+1]]$nDim == 0 else TRUE
                    newExpr$sizeExprs <- if(!indexIsScalar) c(code$args[[iInd + 1]]$sizeExprs) else list(1)
                    newExpr$nDim <- 1
                    newExpr$toEigenize <- 'yes'
                
                    setArg(newExpr, 1, code$args[[1]]$args[[nestedIind + 1]])
                    setArg(newExpr, 2, code$args[[iInd + 1]])
                    setArg(newExpr, 3, 1)
                    setArg(code$args[[1]], nestedIind + 1, newExpr)
                }
            }
            ## The only remaining use of a drop argument is during eigenization to determine if 1xn needs a transpose to become nx1
            ## For that purpose, the drop arg of X[ i[k], j[l] ] should be from the outer part of `[`(X[i, j, drop = TRUE|FALSE], k, l, drop = [TRUE|FALSE]), not from the X[i,j]
            code$args[[1]]$args[['drop']] <- if(dropArgProvided) dropBool else TRUE
            
            ## clear remaining indices
            ## i.e. turn `[`( nimNonseqIndexedd(X, nimNonseqIndexedi(i, k), nimNonseqIndexedi(j, l)), k, l)
            ## into `[`( nimNonseqIndexedd(X, nimNonseqIndexedi(i, k), nimNonseqIndexedi(j, l)))
            code$args[1+(1:numIndices)] <- NULL
            codeCaller <- code$caller
            codeCallerArgID <- code$callerArgID
            ## remove the `[` layer of the current processing
            ## i.e. turn `[`( nimNonseqIndexedd(X, nimNonseqIndexedi(i, k), nimNonseqIndexedi(j, l))) into
            ## imNonseqIndexedd(X, nimNonseqIndexedi(i, k), nimNonseqIndexedi(j, l))
            removeExprClassLayer(code)
            code <- codeCaller$args[[codeCallerArgID]]
            return(if(length(asserts)==0) NULL else asserts)
        }

        ## Now we are in the case where there is no nested indexing, or if there is X[i, j][k, l], it can be chained eigen blocks
        ## Replace with a map expression if needed
        if(!simpleBlockOK) {
            if(typeEnv$.ensureNimbleBlocks) {
                stop(exprClassProcessingErrorMsg(code, "LHS indexing for a multivariate random draw can only use sequential blocks (via ':')."), call. = FALSE)
            }
            ## If this is part of X[i, j] <- Z, convert to coeffSetter(X, i, j) <- Z
            if(code$caller$name %in% assignmentOperators & code$callerArgID == 1) {
                code$name <- 'coeffSetter'
            } else {
                ## otherwise convert `[`(X, i, j) to e.g. nimNonseqIndexedd(X, i, j)
                if(code$type == 'double') code$name <- 'nimNonseqIndexedd' ## this change could get moved to genCpp_generateCpp 
                if(code$type == 'integer') code$name <- 'nimNonseqIndexedi'
                if(code$type == 'logical') code$name <- 'nimNonseqIndexedb'
            }
            ## If we have nimNonseqIndexedd(X, i), make it nimNonseqIndexedd(X, i, 1) for Eigen
            if(length(code$args) - 1 - dropArgProvided == 1) ## only 1 index
                code$args[[3]] <- 1 ## fill in extra 1 for a second dimension.  ## should the index depend on dropArgProvided?
        }
        else { ## a simpleBlock is ok
            if(code$args[[1]]$nDim > 2 | typeEnv$.ensureNimbleBlocks) { ## old-style blocking from >2D down to 2D or 1D, or this is LHS for something like rmnorm, requiring a non-eigen map on LHS.
                ## We have X[i, j, k] where X has dimension > 2
                if(dropArgProvided) code$args[[iDropArg]] <- NULL
                newExpr <- makeMapExprFromBrackets(code, dropBool)
                newExpr$sizeExprs <- code$sizeExprs
                newExpr$type <- code$type
                newExpr$nDim <- code$nDim
                newExpr$toEigenize <- code$toEigenize
                setArg(code$caller, code$callerArgID, newExpr)
            }
            else { ## blocking via Eigen
                ## ## note that any expressions like sum(A) in 1:sum(A) should have already been lifted
                code$name <- 'eigenBlock'
                code$toEigenize <- 'yes'
            }
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
    asserts <- if(recurse) recurseSetSizes(code, symTab, typeEnv) else list()
    byProvided <- code$name == 'nimSeqBy' | code$name == 'nimSeqByLen'
    lengthProvided <- code$name == 'nimSeqLen' | code$name == 'nimSeqByLen'
    integerFrom <- isIntegerEquivalent(code$args[[1]])
    integerTo <- isIntegerEquivalent(code$args[[2]])
    liftExprRanges <- TRUE
    if(integerFrom && integerTo) {
        if((!byProvided && !lengthProvided) ||
           (byProvided && !lengthProvided && is.numeric(code$args[[3]]) && code$args[[3]] == 1)) {
            code$name = ':'
            if(length(code$args) > 2) {
                for(i in length(code$args):3)
                    setArg(code, i, NULL)
            }
            asserts <- c(asserts, sizeColonOperator(code, symTab, typeEnv, recurse = FALSE))
            return(if(length(asserts)==0) NULL else asserts)
        }
    }
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
        thisSizeExpr <- substitute(calcSeqLength(FROM_, TO_, BY_),##1 + floor((TO_ - FROM_) / BY_),
                                   list(FROM_ = parse(text = nimDeparse(code$args[[1]]), keep.source = FALSE)[[1]],
                                        TO_ = parse(text = nimDeparse(code$args[[2]]), keep.source = FALSE)[[1]],
                                        BY_ = parse(text = nimDeparse(code$args[[3]]), keep.source = FALSE)[[1]]))
      }
    } else { ## must be lengthProvided
      code$name <- 'nimSeqLenD'
      thisSizeExpr <- parse(text = nimDeparse(code$args[[4]]), keep.source = FALSE)[[1]]
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

    for(i in 1:2) {
        if(inherits(code$args[[i]], 'exprClass')) {
            if(!code$args[[i]]$isName) {
              if(! (code$args[[i]]$name == '[' && (code$args[[i]]$args[[1]]$name == 'dim' && code$args[[i]]$args[[1]]$args[[1]]$name == 'nfVar'))){
                asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv, forceType = "integer") )
              }
            }
        }
    }
    
    code$type <- 'double'
    code$nDim <- 1
    code$toEigenize <- 'maybe' 
    
    ## could generate an assertion that second arg is >= first arg
    if(is.numeric(code$args[[1]]) & is.numeric(code$args[[2]])) {
        if(code$args[[1]] > code$args[[2]])
            stop("sizeColonOperator: negative indexing cannot be compiled.")
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
        sizeProc_storage_mode(expr)
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
        if(!nimbleOptions('experimentalNewSizeProcessing') ) {
            if(a1$nDim == 0) {
                ## Argument is scalar.
                ## If it results from vector operation (e.g. inprod)
                ## lift that to an intermediate
                if(a1$toEigenize == 'yes') {
                    asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
                    a1 <- code$args[[1]]
                }
            } else {
                ## Argument is non-scalar.  In this case, the
                ## expression will be eigenized, so we must lift the
                ## argument to an intermediate if it *can't* be
                ## eigenized.
                if(a1$toEigenize == 'no') {
                    asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
                    a1 <- code$args[[1]]
                }
            }
        }
        code$nDim <- a1$nDim
        code$sizeExprs <- a1$sizeExprs
    } else {
        code$nDim <- 0
        code$sizeExprs <- list()
    }
    code$type <- setReturnType(code$name, getArgumentType(a1))
    if(length(code$nDim) != 1) stop(exprClassProcessingErrorMsg(code, 'In sizeUnaryCwise: nDim is not set.'), call. = FALSE)
    if(!nimbleOptions('experimentalNewSizeProcessing') ) code$toEigenize <- if(code$nDim > 0) 'yes' else 'maybe'
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
        if(code$args[[1]]$nDim == 0) 
            stop(exprClassProcessingErrorMsg(code, 'NIMBLE compiler does not support reduction operations on scalar arguments.'), call. = FALSE)
        if(!nimbleOptions('experimentalNewSizeProcessing') ) {
            if(!code$args[[1]]$isName) {
                if(code$args[[1]]$toEigenize == 'no') {
                    asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
                }
            }
        }
    }

    code$nDim <- 0
    code$sizeExprs <- list()
    code$type <- setReturnType(code$name, code$args[[1]]$type)
    if(!nimbleOptions('experimentalNewSizeProcessing') ) code$toEigenize <- 'yes'

    if(!nimbleOptions('experimentalNewSizeProcessing') ) {
        if(!(code$caller$name %in% c('{','<-','<<-','='))) {
            asserts <- c(asserts, sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
        }
    }

    if(length(asserts) == 0) NULL else asserts
}

## There's no real point in annotating return.  Just need to recurse and lift 
sizeReturn <- function(code, symTab, typeEnv) {
    if(length(code$args) > 1) stop(exprClassProcessingErrorMsg(code, 'return has argument length > 1.'), call. = FALSE)
    code$toEigenize <- 'no'
    if(!exists('return', envir = typeEnv)) stop(exprClassProcessingErrorMsg(code, 'There was no returnType declaration and the default is missing.'), call. = FALSE)

    if(length(code$args) == 0) {
        if(!identical(typeEnv$return$type, 'void'))
            stop(exprClassProcessingErrorMsg(code, 'return() with no argument can only be used with returnType(void()), which is the default if there is no returnType() statement.'), call. = FALSE)
        return(invisible(NULL))
    }
    if(identical(typeEnv$return$type, 'void'))
        stop(exprClassProcessingErrorMsg(code, 'returnType was declared void() (default) (or something invalid), which is not consistent with the object you are trying to return.'), call. = FALSE)
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    if(inherits(code$args[[1]], 'exprClass')) {
        if(typeEnv$return$type == 'nimbleList' || isTRUE(code$args[[1]]$type == 'nimbleList')) {
            if(typeEnv$return$type != 'nimbleList') stop(exprClassProcessingErrorMsg(code, paste0('return() argument is a nimbleList but returnType() statement gives a different type')), call. = FALSE)
            if(code$args[[1]]$type != 'nimbleList') stop(exprClassProcessingErrorMsg(code, paste0('returnType statement gives a nimbleList type but return() argument is not the right type')), call. = FALSE)
            ## equivalent to symTab$getSymbolObject(code$args[[1]]$name)$nlProc, if it is a name
            if(!identical(code$args[[1]]$sizeExprs$nlProc, typeEnv$return$sizeExprs$nlProc)) stop(exprClassProcessingErrorMsg(code, paste0('nimbleList given in return() argument does not match nimbleList type declared in returnType()')), call. = FALSE)
        } else { ## check numeric types and nDim
            fail <- FALSE
            if(is.null(code$args[[1]]$type)) {  # Issue 1364
                failMsg <- paste0(code$args[[1]]$name, " is not available or its output type is unknown.")
                fail <- TRUE
            } else {
                if(!identical(code$args[[1]]$type, typeEnv$return$type)) {
                    if(typeEnv$return$nDim > 0) { ## allow scalar casting of returns without error
                        failMsg <- paste0('Type ', code$args[[1]]$type, ' of the return() argument does not match type ',  typeEnv$return$type, ' given in the returnType() statement (void is default).')
                        fail <- TRUE
                    }
                }
                if(!isTRUE(all.equal(code$args[[1]]$nDim, typeEnv$return$nDim))) {
                    failMsg <- paste0( if(exists("failMsg", inherits = FALSE)) paste0(failMsg,' ') else character(),
                                      paste0('Number of dimensions ', code$args[[1]]$nDim, ' of the return() argument does not match number ',  typeEnv$return$nDim, ' given in the returnType() statement.'))
                    fail <- TRUE
                }
            }
            if(fail)
                stop(exprClassProcessingErrorMsg(code, failMsg), call. = FALSE)
        }
        if(!code$args[[1]]$isName) {
            liftArg <- FALSE
            if(code$args[[1]]$toEigenize == 'yes')
                liftArg <- TRUE
            else if(anyNonScalar(code$args[[1]]))
                liftArg <- TRUE
            if(liftArg)
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv, forceAssign = TRUE))
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
    if(!nimbleOptions('experimentalNewSizeProcessing') ) {
        if(a1$toEigenize == 'no') {
            asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
            a1 <- code$args[[1]]
        }
        if(a2$toEigenize == 'no') {
            asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
            a2 <- code$args[[2]]
        }
    }
    
    ## Note that we could insert RUN-TIME adaptation of mat %*% vec and vec %*% mat
    ## but to do so we would need to generate trickier sizeExprs
    ## For now, a vector on the right will be turned into a column
    ## and a vector on the left will be turned into a row
    ## The programmer can always use asRow or asCol to control it explicitly

    if(a1$nDim == 1 & a2$nDim == 1) {
        origSizeExprs <- a1$sizeExprs[[1]]
        a1 <- insertExprClassLayer(code, 1, 'asRow', type = a1$type, nDim = 2)
        a1$sizeExprs <- c(list(1), origSizeExprs)
        origSizeExprs <- a2$sizeExprs[[1]]
        a2 <- insertExprClassLayer(code, 2, 'asCol', type = a2$type, nDim = 2)
        a2$sizeExprs <- c(origSizeExprs, list(1))
    } else {
        if(a1$nDim == 1) {
            if(a2$nDim != 2) stop(exprClassProcessingErrorMsg(code, paste0('In sizeMatrixMult: First arg has nDim = 1 and 2nd arg has nDim = ', a2$nDim, '.')), call. = FALSE)
            origSizeExprs <- a1$sizeExprs[[1]]
            ## For first argument, default to asRow unless second argument has only one row, in which case make first asCol
            if(identical(a2$sizeExprs[[1]], 1)) {
                a1 <- insertExprClassLayer(code, 1, 'asCol', type = a1$type, nDim = 2)
                a1$sizeExprs <- c(origSizeExprs, list(1))
            }
            else {
                a1 <- insertExprClassLayer(code, 1, 'asRow', type = a1$type, nDim = 2)
                a1$sizeExprs <- c(list(1), origSizeExprs)
            }
        } else if(a2$nDim == 1) {
            origSizeExprs <- a2$sizeExprs[[1]]
            if(a1$nDim != 2) stop(exprClassProcessingErrorMsg(code, paste0('In sizeMatrixMult: Second arg has nDim = 1 and 1st arg has nDim = ', a1$nDim, '.')), call. = FALSE)
            if(identical(a1$sizeExprs[[2]], 1)) {
                a2 <- insertExprClassLayer(code, 2, 'asRow', type = a2$type, nDim = 2)
                a2$sizeExprs <- c(list(1), origSizeExprs)
           }
            else { 
                a2 <- insertExprClassLayer(code, 2, 'asCol', type = a2$type, nDim = 2)
                a2$sizeExprs <- c(origSizeExprs, list(1))
            }
        }
    }
    code$nDim <- 2
    code$sizeExprs <- list(a1$sizeExprs[[1]], a2$sizeExprs[[2]])
    code$type <- setReturnType(code$name, arithmeticOutputType(a1$type, a2$type))
    if(!nimbleOptions('experimentalNewSizeProcessing') ) code$toEigenize <- 'yes'
    assertMessage <- paste0("Run-time size error: expected ", deparse(a1$sizeExprs[[2]]), " == ", deparse(a2$sizeExprs[[1]]))
    newAssert <- identityAssert(a1$sizeExprs[[2]], a2$sizeExprs[[1]], assertMessage)
    if(is.null(newAssert))
        return(asserts)
    else
        return(c(asserts, list(newAssert)))
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
        if(!nimbleOptions('experimentalNewSizeProcessing') ) {
            if(a1$toEigenize == 'no') {
                asserts <- c(asserts, sizeInsertIntermediate(code, 1, symTab, typeEnv))
                a1 <- code$args[[1]]
            }
        }
        a1Drop <- dropSingleSizes(a1$sizeExprs)
        a1DropNdim <- length(a1Drop$sizeExprs)
        a1nDim <- a1$nDim
        a1sizeExprs <- a1$sizeExprs
        a1type <- a1$type
        if(!nimbleOptions('experimentalNewSizeProcessing') ) a1toEigenize <- a1$toEigenize
    } else {
        a1DropNdim <- 0
        a1nDim <- 0
        a1sizeExprs <- list()
        a1type <- sizeProc_storage_mode(a1)
        if(!nimbleOptions('experimentalNewSizeProcessing') ) a1toEigenize <- 'maybe'
    }
    if(inherits(a2, 'exprClass')) {
        if(!nimbleOptions('experimentalNewSizeProcessing') ) {
            if(a2$toEigenize == 'no') {
                asserts <- c(asserts, sizeInsertIntermediate(code, 2, symTab, typeEnv))
                a2 <- code$args[[2]]
            }
        }
        a2Drop <- dropSingleSizes(a2$sizeExprs)
        a2DropNdim <- length(a2Drop$sizeExprs)
        a2nDim <- a2$nDim
        a2sizeExprs <- a2$sizeExprs
        a2type <- a2$type
        if(!nimbleOptions('experimentalNewSizeProcessing') )  a2toEigenize <- a2$toEigenize
    } else {
        a2DropNdim <- 0
        a2nDim <- 0
        a2sizeExprs <- list()
        a2type <- sizeProc_storage_mode(a2)
        if(!nimbleOptions('experimentalNewSizeProcessing') )  a2toEigenize <- 'maybe'
    }
    
    ## Choose the output type by type promotion
    if(length(a1type) == 0) {stop('Problem with type of arg1 in sizeBinaryCwise', call. = FALSE)}
    if(length(a2type) == 0) {stop('Problem with type of arg2 in sizeBinaryCwise', call. = FALSE)}
    code$type <- setReturnType(code$name, arithmeticOutputType(a1type, a2type))

    if(!nimbleOptions('experimentalNewSizeProcessing') ) {
        forceYesEigenize <- identical(a1toEigenize, 'yes') | identical(a2toEigenize, 'yes')
        code$toEigenize <- if(a1DropNdim == 0 & a2DropNdim == 0)
                               if(forceYesEigenize)
                                   'yes'
                               else
                                   'maybe'
                           else 'yes'
    }
    
    ## Catch the case that there is at least one scalar-equivalent (all lengths == 1)
    ## experimentalNewSizeProcessing: The 3 'code$toEigenize <- ' should be redundant with above and could be removed during refactor
    if(a1DropNdim == 0 | a2DropNdim == 0) { 
        ## Here we will process effective scalar additions
        ## and not do any other type of size promotion/dropping
        if(a1DropNdim == 0) { ## a1 is scalar-equiv
            if(a2DropNdim == 0) { ##both are scalar-equiv
                code$nDim <- max(a1nDim, a2nDim) ## use the larger nDims
                code$sizeExprs <- rep(list(1), code$nDim) ## set sizeExprs to all 1
                if(!nimbleOptions('experimentalNewSizeProcessing') ) code$toEigenize <- if(forceYesEigenize) 'yes' else 'maybe'
            } else {
                ## a2 is not scalar equiv, so take nDim and sizeExprs from it
                code$nDim <- a2nDim
                code$sizeExprs <- a2sizeExprs
                if(!nimbleOptions('experimentalNewSizeProcessing') ) code$toEigenize <- 'yes'
            }
        } else { ## a2 is scalar-equiv, and a1 is not
            code$nDim <- a1nDim
            code$sizeExprs <- a1sizeExprs
            if(!nimbleOptions('experimentalNewSizeProcessing') ) code$toEigenize <- 'yes'
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
                             nimArr_rlkj_corr_cholesky = list(c(0, 0), ## eta, p
                                                              function(code) {
                                                                  ## example
                                                                  ## Rcode <- quote(nimArr_rlkj_corr_cholesky(eta = k * rho, p = k + a) )
                                                                  ## code <- nimble:::RparseTree2ExprClasses(Rcode)
                                                                  dimCode <- parse(text = nimDeparse(code$args[[2]]), keep.source = FALSE)[[1]]
                                                                  sizeExprs <- list(dimCode, dimCode)
                                                                  ## dimCode should be quote( k + a )
                                                                  list(nDim = 2,
                                                                       sizeExprs = sizeExprs)
                                                              },
                                                              'double'),  # '1' won't work here; problem is that no matrices are input but matrix is output
                             nimArr_rwish_chol = list(c(2, 0, 0, 0), ## chol, df, prec_param, overwrite_inputs
                                 1, 'double'),
                             nimArr_rinvwish_chol = list(c(2, 0, 0), ## chol, df, prec_param
                                 1, 'double'),
			     nimArr_rcar_normal = list(c(1, 1, 1, 0, 0, 0), 3, 'double'), ## adj, wgts, num, tau, c, zero_mean, answer size comes from num
			     nimArr_rcar_proper = list(c(1, 1, 1, 1, 1, 0, 0, 1), 1, 'double'), ## mu, C, adj, num, M, tau, gamma, evs, answer size comes from mu
                             nimArr_rmulti = list(c(0, 1), ## size, probs
                                 2, 'double'), ## We treat integer rv's as doubles
                             nimArr_rdirch = list(c(1), 1, 'double')) 

sizeRmultivarFirstArg <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    notOK <- FALSE
    checkList <- mvFirstArgCheckLists[[code$name]]
    if(!is.null(checkList)) {
        if(length(code$args) < length(checkList[[1]])) stop(exprClassProcessingErrorMsg(code, 'Not enough arguments provided.'), call. = FALSE)
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
        stop(exprClassProcessingErrorMsg(code, 'Some argument(s) have the wrong dimension.'), call. = FALSE) 
    }

    if(!(is.function(returnSizeArgID)))
        if(!inherits(code$args[[returnSizeArgID]], 'exprClass'))
            stop(exprClassProcessingErrorMsg(code, paste0('Expected ', nimDeparse(code$args[[returnSizeArgID]]) ,' to be an expression or function.')), call. = FALSE) 

    code$type <- returnType
    code$toEigenize <- 'maybe'
    if(!is.function(returnSizeArgID)) {
        code$nDim <- code$args[[returnSizeArgID]]$nDim
        code$sizeExprs <- code$args[[returnSizeArgID]]$sizeExprs
    } else {
        sizeInfo <- returnSizeArgID(code)
        code$nDim <- sizeInfo$nDim
        code$sizeExprs <- sizeInfo$sizeExprs
    }
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
    	if(code$args[[1]]$type == 'nimbleFunction') lift <- FALSE
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

sizePassByMap <- function(code, symTab, typeEnv) {
    ensureNimbleBlocks <- typeEnv$.ensureNimbleBlocks
    typeEnv$.ensureNimbleBlocks <- TRUE
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    code <- removeExprClassLayer(code, 1)
    typeEnv$.ensureNimbleBlocks <- ensureNimbleBlocks
    asserts
}

generalFunSizeHandlerFromSymbols <- function(code, symTab, typeEnv, returnSymbol, argSymTab, chainedCall = FALSE) {
    ## symbols should be in order
    useArgs <- unlist(lapply(argSymTab$symbols, function(x) {
        if(!is.null(x[['type']]))
            as.character(x$type) %in% c('double', 'integer', 'logical')
        else
            FALSE
        }))
    
    if(chainedCall) useArgs <- c(FALSE, useArgs)
    if(length(code$args) != length(useArgs)) {
        stop(exprClassProcessingErrorMsg(code, 'In generalFunSizeHandlerFromSymbols: Wrong number of arguments.'), call. = FALSE)
    }
    ## Note this is NOT checking the dimensions of each arg. useArgs just means it will recurse on that and lift or do as needed

    ## Ensure that simple maps being passed will be passed without extra
    ## copy that would occur from lifting an Eigen expression.
    for(i in seq_along(code$args)) {
        ## check if each argument is purely of the form x[...]
        ## (note that x[...][...] might also be valid for passByMap
        ## but it is not handled that way currently.
        if(inherits(code$args[[i]], 'exprClass')) {
            if(code$args[[i]]$name == "[") {
                if(inherits(code$args[[i]]$args[[1]],
                            'exprClass')) { ## must be true, but I'm being defensive
                    if(code$args[[i]]$args[[1]]$isName) {
                        insertExprClassLayer(code, i, 'passByMap')
                    }
                }
            }
        }
    }

    asserts <- recurseSetSizes(code, symTab, typeEnv, useArgs)

    ## lift any argument that is an expression
    for(i in seq_along(code$args)) {
        if(useArgs[i]) {
            if(inherits(code$args[[i]], 'exprClass')) {
                if(!code$args[[i]]$isName) {
                    forceType <- NULL
                    iSym <- i - chainedCall
                    if(argSymTab$symbols[[iSym]]$nDim == 0) ## We're only here if useArgs[i] is TRUE, which means nDim and type should be set
                        forceType <- argSymTab$symbols[[iSym]]$type
                    asserts <- c(asserts, sizeInsertIntermediate(code, i, symTab, typeEnv, forceType = forceType) )
                }
            }
        }
    }
    if(inherits(returnSymbol, 'symbolNimbleList')) {
        code$type <- 'nimbleList'
        code$sizeExprs <- returnSymbol
        code$toEigenize <- 'maybe'
        code$nDim <- 0
        liftIfAmidExpression <- TRUE
    } else {
        returnSymbolBasic <- inherits(returnSymbol, 'symbolBasic')
        returnTypeLabel <- if(returnSymbolBasic)
                               returnSymbol$type
                           else {
                               stop(exprClassProcessingErrorMsg(code, 'In generalFunSizeHandlerFromSymbols: Problem with return type.'), call. = FALSE)
                           }
        if(returnTypeLabel == 'void') {
            code$type <- returnTypeLabel
            code$toEigenize <- 'unknown'
            return(asserts)
        }
        returnNDim <- if(returnSymbolBasic) returnSymbol$nDim
                      else if(length(returnType) > 1) as.numeric(returnType[[2]]) else 0
                                                                
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
