## It is hoped substantial sets of these calls can be combined
## or implemented by a rule system.

toEigenizeNoCalls <- c('dim', 'run.time','nimOptimDefaultControl')

toEigenizeYesCalls <- c(paste0('nimDiagonal', c('D','I','B')),
                        'diagonal',
                        'which',
                        scalar_distribution_dFuns,
                        scalar_distribution_pFuns,
                        scalar_distribution_qFuns,
                        paste0(c('d','q','p'), 't'),
                        paste0(c('d','q','p'), 'exp'),
                        scalar_distribution_rFuns,
                        'rt',
                        'rexp',
                        'asRow','asCol',
                        nimbleListReturningOperators, ## should these be in 'maybe'?
                        'setAll' ## deprecate?
                        )

toEigenizeMaybeCalls <- c('map',
                          c('decide', 'size', 'getsize','getNodeFunctionIndexedInfo', 'endNimbleTimer'),
                          c('blank', 'nfMethod', 'getPtr', 'startNimbleTimer'))

toEigenizeUseRuleCalls <- c('nimPrint')

toEigenCalls <- c(makeCallList(binaryOperators, 'toEigenBinaryCwise'),
                  makeCallList(binaryMidLogicalOperators, 'sizeBinaryCwiseLogical'),
                  makeCallList(binaryOrUnaryOperators, 'toEigenBinaryUnaryCwise'),
                  makeCallList(unaryOperators, 'toEigenUnaryCwise'), 
                  makeCallList(unaryOrNonaryOperators, 'sizeUnaryNonaryCwise'),
                  makeCallList(assignmentOperators, 'toEigenAssign'), 
                  makeCallList(reductionUnaryOperators, 'sizeUnaryReduction'), 
                  makeCallList(matrixSquareReductionOperators, 'sizeMatrixSquareReduction'),
                  makeCallList(reductionBinaryOperators, 'sizeBinaryReduction'),
                  makeCallList(matrixMultOperators, 'sizeMatrixMult'), 
                  makeCallList(matrixFlipOperators, 'sizeTranspose'),
                  makeCallList(matrixSolveOperators, 'sizeSolveOp'), 
                  makeCallList(matrixSquareOperators, 'sizeUnaryCwiseSquare'),
                  ##               makeCallList(nimbleListReturningOperators, 'sizeNimbleListReturningFunction'),
                  nimOptim = 'toEigenOptim',
                  ##nimOptimDefaultControl = 'sizeOptimDefaultControl',
                  
                  makeCallList(paste0('nimC',c('d','i','b')), 'toEigenConcatenate'),
                  makeCallList(paste0('nimRep',c('d','i','b')), 'toEigenRep'),
                  list('debugSizeProcessing' = 'sizeProxyForDebugging',
                       ##      diag = 'toEigenYes',
                       ##      dim = 'toEigenNo',
                       ##       RRtest_add = 'toEigenYes', ## testing!
                       ##      which = 'toEigenYes',
                       
                       nimRep = 'sizeRep',
                       nimSeqBy = 'sizeSeq',
                       nimSeqLen = 'sizeSeq',
                       nimSeqByLen = 'sizeSeq',
                       'return' = 'sizeReturn',
                       ##  'asRow' = 'sizeAsRowOrCol',
                       ##  'asCol' = 'sizeAsRowOrCol',
                       makeNewNimbleListObject = 'toEigenNewNimbleList',
                       getParam = 'toEigenGetParam',
                       getBound = 'toEigenGetBound',
                       ## nimSwitch = 'sizeSwitch', ## NEVER GETS MARKED WITH TOEIGENIZE
                       asDoublePtr = 'toEigenIgnoreForNow',
                       '[' = 'sizeIndexingBracket',
                       ##'[[' = 'sizeDoubleBracket', ## for nimbleFunctionList, this will always  go through chainedCall(nfList[[i]], 'foo')(arg1, arg2) ## NEVER GETS MARKED WITH TOEIGENIZE
                       chainedCall = 'toEigenChainedCall',
                       ## nfVar = 'sizeNFvar', ## NEVER GETS MARKED WITH TOEIGENIZE
                       ## map = 'sizemap',
                       ':' = 'sizeColonOperator',
                       ##dim = 'sizeDimOperator',
                       'if' = 'recurseSetSizes', ##OK
                       'while' = 'recurseSetSizes',
                       callC = 'sizecallC', 
                       'for' = 'toEigenFor', 
                       cppPointerDereference = 'toEigenCppPointerDereference',
                       values = 'toEigenValues',
                       '(' = 'sizeUnaryCwise',
                       setSize = 'toEigenSetSize', ## OK but not done for numericLists
                       ## resizeNoPtr = 'sizeResizeNoPtr', ## may not be used any more 
                       nimArr_rcat = 'toEigenScalarRecurse',
                       nimArr_rinterval = 'toEigenScalarRecurse',
                       ##nimPrint = 'sizeforceEigenize',
                       ##nimCat = 'sizeforceEigenize',
                       as.integer = 'sizeUnaryCwise', ## Note as.integer and as.numeric will not work on a non-scalar yet
                       as.numeric = 'sizeUnaryCwise',
                       nimArrayGeneral = 'toEigenNimArrayGeneral',
                       ##setAll = 'sizeOneEigenCommand',
                       voidPtr = 'sizeVoidPtr'
                       ),
                    ##run.time = 'sizeRunTime'),
                  ##               makeCallList(scalar_distribution_dFuns, 'toEigenYes'),
                  ##               makeCallList(scalar_distribution_pFuns, 'toEigenYes'),
                  ##               makeCallList(scalar_distribution_qFuns, 'sizeRecyclingRule'),
                  ##               makeCallList(scalar_distribution_rFuns, 'sizeRecyclingRuleRfunction'),
                  makeCallList(distributionFuns[!(distributionFuns %in% c(scalar_distribution_dFuns, scalar_distribution_pFuns, scalar_distribution_qFuns, scalar_distribution_rFuns))], 'toEigenScalarRecurse'),
                                        # R dist functions that are not used by NIMBLE but we allow in DSL
                  ## makeCallList(paste0(c('d','q','p'), 't'), 'sizeRecyclingRule'),
                  ## rt = 'sizeRecyclingRuleRfunction',
                  ## makeCallList(paste0(c('d','q','p'), 'exp'), 'sizeRecyclingRule'),
                  ## rexp = 'sizeRecyclingRuleRfunction',
                  makeCallList(c('isnan','ISNAN','ISNA'), 'toEigenScalarRecurse'),
                  makeCallList(c('nimArr_dmnorm_chol', 'nimArr_dmvt_chol', 'nimArr_dwish_chol', 'nimArr_dinvwish_chol', 'nimArr_dmulti', 'nimArr_dcat', 'nimArr_dinterval', 'nimArr_ddirch'), 'toEigenScalarRecurse'),
                  makeCallList(c('nimArr_rmnorm_chol', 'nimArr_rmvt_chol', 'nimArr_rwish_chol', 'nimArr_rinvwish_chol', 'nimArr_rmulti', 'nimArr_rdirch'), 'sizeRmultivarFirstArg'),
                  ##      makeCallList(c('decide', 'size', 'getsize','getNodeFunctionIndexedInfo', 'endNimbleTimer'), 'sizeScalar'),
                  ## makeCallList(c('calculate','calculateDiff', 'getLogProb'), 'toEigenScalarModelOp'),
                  simulate = 'toEigenIgnoreForNow'
                  ##               makeCallList(c('blank', 'nfMethod', 'getPtr', 'startNimbleTimer'), 'sizeUndefined')
                  )


exprClasses_setToEigenize <- function(code, symTab, typeEnv) { ## input code is exprClass
    ## name:
    if(code$isName) {
        return(list()) ## names don't have toEigen annotation
    }
    if(code$isCall) {
        if(code$name == '{') {
            ## recurse over lines
            for(i in seq_along(code$args)) {
                if(inherits(code$args[[i]], 'exprClass')) {
                    newAsserts <- exprClasses_setToEigenize(code$args[[i]], symTab, typeEnv)
                    code$args[[i]]$assertions <- if(is.null(newAsserts)) list() else newAsserts
                }
            }
            return(invisible(NULL))
        }
        thisCall <- toEigenCalls[[code$name]]
        if(!is.null(thisCall)) {
            test0 <- eval(call(thisCall, code, symTab, typeEnv))
            return(test0)
        }
        if(symTab$symbolExists(code$name, TRUE)) { ## could be a nimbleFunction object
            return(toEigenNimbleFunction(code, symTab, typeEnv) )
        }
        ## Finally, it could be an RCfunction (a nimbleFunction with no setup == a simple function) {

        if(exists(code$name)) {
            obj <- get(code$name)
            if(is.rcf(obj)) { ## it is an RC function
                message('figure out how to extract the RCfunProc that should already exist')
                browser()
                RCfunProc <- typeEnv$.nimbleProject$compileRCfun(obj, initialTypeInference = TRUE)
                
                return(toEigenRCfunction(code, symTab, typeEnv, nfmObj, RCfunProc))
            }
        }
    }
    invisible(NULL)
}

toEigenNo <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'no'
    return(list())
}

toEigenYes <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'yes'
    return(list())
}

toEigenNo <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'no'
    return(list())
}

toEigenMaybe <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'maybe'
    return(list())
}

toEigenUseRule <- function(code, symTab, typeEnv) {
    message('need to implement rules for toEigenUseRule')
    return(list())
}

toEigenConcatenate <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'yes'
    return(list())
}

toEigenRep <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'yes'
    return(list())
}

toEigenNewNimbleList <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'maybe'
    return(list())
}

toEigenNimArrayGeneral <- function(code, symTab, typeEnv) {
    ## won't even have the same name at this point
    stop('setting toEigenize for nimArrayGeneral must be handled')
}

toEigenGetParam <- function(code, symTab, typeEnv) {
    stop('setting toEigenize for getParam must be handled')
    return(list())
}

toEigenGetBound <- function(code, symTab, typeEnv) {
    stop('setting toEigenize for getBound must be handled')
    return(list())
}

toEigenOptim <- function(code, symTab, typeEnv) {
    code$toEigenize <- 'no'
    return(list())
}

toEigenCppPointerDereference <- function(code, symTab, typeEnv) {
    code$toEigenize <- code$args[[1]]$toEigenize
    return(list())
}

toEigenChainedCall <- function(code, symTab, typeEnv) {
    message('setting toEigenize for chained call must be handled')
    return(list())
}

toEigenValues <- function(code, symTab, typeEnv) {
    message('setting toEigenize for values must be handled')
    return(list())
}

toEigenRCfunction <- function(code, symTab, typeEnv) {
    message('setting toEigenize for RCfunction must be handled')
    return(list())
}

toEigenNimbleFunction <- function(code, symTab, typeEnv) {
    message('setting toEigenize for RCfunction must be handled')
    return(list())
}

toEigenSetSize <- function(code, symTab, typeEnv) {
    message('setting toEigenize for toEigenSetSize must be handled')
    return(list())
}

toEigenFor <- function(code, symTab, typeEnv) {
    message('setting toEigenize for toEigenFor must be handled')
    return(list())
}

toEigenAssign <- function(code, symTab, typeEnv) {
    message('setting toEigenize for toEigenAssign must be handled')
    return(list())
}

toEigenIgnoreForNow <- function(code, symTab, typeEnv) {
    message('Can not set toEigenize: to ignore for now')
    return(list())
}

toEigenScalarModelOp <- function(code, symTab, typeEnv) {
    message('setting toEigenize for toEigenScalarModelOp must be handled')
    return(list())
}

toEigenScalarRecurse <- function(code, symTab, typeEnv) {
    message('setting toEigenize for toEigenScalarRecurse must be handled')
    return(list())
}

toEigenBinaryUnaryCwise <- function(code, symTab, typeEnv) {
    if(length(code$args) == 1) return(toEigenUnaryCwise(code, symTab, typeEnv))
    if(length(code$args) == 2) return(toEigenBinaryCwise(code, symTab, typeEnv))
    stop(exprClassProcessingErrorMsg(code, paste0('In toEigenBinaryUnarycWise: Length of arguments is not 1 or 2.')), call. = FALSE)
}

toEigenBinaryCwise <- function(code, symTab, typeEnv) {
    message('setting toEigenize for toEigenBinaryCwise must be handled')
    return(list())
}

