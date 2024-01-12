## It is hoped substantial sets of these calls can be combined
## or implemented by a rule system.

toEigenizeNoCalls <- c('dim',
                       'run.time',
                       'nimOptimDefaultControl')

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
                          c('decide',
                            'size',
                            'getsize',
                            'getNodeFunctionIndexedInfo',
                            'endNimbleTimer'),
                          c('blank',
                            'nfMethod',
                            'getPtr',
                            'startNimbleTimer'))

toEigenizeUseRuleCalls <- c('nimPrint')

toEigenCalls <- c(
    makeCallList(binaryOperators, 'toEigenBinaryCwise'),             
    makeCallList(binaryMidLogicalOperators, 'sizeBinaryCwiseLogical'),
    makeCallList(binaryOrUnaryOperators, 'toEigenBinaryUnaryCwise'), 
    makeCallList(unaryOperators, 'toEigenUnaryCwise'),               
    makeCallList(unaryOrNonaryOperators, 'sizeUnaryNonaryCwise'),
    makeCallList(assignmentOperators, 'toEigenAssign'),              
    makeCallList(reductionUnaryOperators, 'sizeUnaryReduction'),     
    makeCallList(matrixSquareReductionOperators, 'sizeMatrixSquareReduction'),
    makeCallList(reductionBinaryOperators, 'sizeBinaryReduction'),
    makeCallList(matrixMultOperators, 'toEigenMatrixMult'),          
    makeCallList(matrixFlipOperators, 'sizeTranspose'),
    makeCallList(matrixSolveOperators, 'sizeSolveOp'), 
    makeCallList(matrixSquareOperators, 'sizeUnaryCwiseSquare'),
    nimOptim = 'toEigenOptim',
    
    makeCallList(paste0('nimC',c('d','i','b')), 'toEigenConcatenate'),
    makeCallList(paste0('nimRep',c('d','i','b')), 'toEigenRep'),
    list('debugSizeProcessing' = 'sizeProxyForDebugging',
         nimRep = 'sizeRep',
         nimSeqBy = 'sizeSeq',
         nimSeqLen = 'sizeSeq',
         nimSeqByLen = 'sizeSeq',
         'return' = 'sizeReturn',
         makeNewNimbleListObject = 'toEigenNewNimbleList',
         getParam = 'toEigenGetParam',
         getBound = 'toEigenGetBound',
         asDoublePtr = 'toEigenIgnoreForNow',
         '[' = 'sizeIndexingBracket',
         chainedCall = 'toEigenChainedCall',
         ':' = 'sizeColonOperator',
         'if' = 'recurseSetSizes', 
         'while' = 'recurseSetSizes',
         callC = 'sizecallC', 
         'for' = 'toEigenFor', 
         cppPointerDereference = 'toEigenCppPointerDereference',
         values = 'toEigenValues',
         '(' = 'sizeUnaryCwise',
         setSize = 'toEigenSetSize', 
         nimArr_rcat = 'toEigenScalarRecurse',
         nimArr_rinterval = 'toEigenScalarRecurse',
         as.integer = 'sizeUnaryCwise', 
         as.numeric = 'sizeUnaryCwise',
         nimArrayGeneral = 'toEigenNimArrayGeneral',
         voidPtr = 'sizeVoidPtr'
         ),
    makeCallList(distributionFuns[
        !(distributionFuns %in% c(scalar_distribution_dFuns,
                                  scalar_distribution_pFuns,
                                  scalar_distribution_qFuns,
                                  scalar_distribution_rFuns))
    ], 'toEigenScalarRecurse'),
    makeCallList(c('isnan',  # Needed?
                   'nim_IsNA',
                   'nim_isnancpp'), 'toEigenScalarRecurse'),
    makeCallList(c('nimArr_dmnorm_chol',
                   'nimArr_dmvt_chol',
                   'nimArr_dlkj_corr_cholesky',
                   'nimArr_dwish_chol',
                   'nimArr_dinvwish_chol',
                   'nimArr_dmulti',
                   'nimArr_dcat',
                   'nimArr_dinterval',
                   'nimArr_ddirch'), 'toEigenScalarRecurse'),
    makeCallList(c('nimArr_rmnorm_chol',
                   'nimArr_rmvt_chol',
                   'nimArr_rlkj_corr_cholesky',
                   'nimArr_rwish_chol',
                   'nimArr_rinvwish_chol',
                   'nimArr_rmulti',
                   'nimArr_rdirch'), 'sizeRmultivarFirstArg'),
    simulate = 'toEigenIgnoreForNow'
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
                    newAsserts <-
                        exprClasses_setToEigenize(code$args[[i]],
                                                  symTab,
                                                  typeEnv)
                    code$args[[i]]$assertions <-
                        c(code$args[[i]]$assertions,
                          if(is.null(newAsserts)) list() else newAsserts
                          )
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
                RCfunProc <-
                    typeEnv$.nimbleProject$compileRCfun(obj,
                                                        initialTypeInference = TRUE)
                
                return(toEigenRCfunction(code, symTab, typeEnv)) ## These additional arguments are not part of the function signature: ', nfmObj, RCfunProc))'
            }
        }
    }
    invisible(NULL)
}


recurseSetToEigenize <- function(code, symTab, typeEnv, useArgs = rep(TRUE, length(code$args))) {
    asserts <- list()
    for(i in seq_along(code$args)) {
        if(useArgs[i]) {
            if(inherits(code$args[[i]], 'exprClass')) {
                asserts <- c(asserts,
                             exprClasses_setToEigenize(code$args[[i]],
                                                       symTab,
                                                       typeEnv)
                             )
            }
        }
    }
    if(length(asserts)==0) NULL else asserts
}

toEigenInsertIntermediate <- function(code, argID, symTab, typeEnv) {
    newLine <- sizeInsertIntermediate(code, argID, symTab, typeEnv)
    exprClasses_setToEigenize(newLine, symTab, typeEnv)
    newLine
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

toEigenAssign <- function(code, symTab, typeEnv) {
    asserts <- recurseSetToEigenize(code, symTab, typeEnv, useArgs = c(FALSE, TRUE))
    asserts <- c(asserts, toEigenAssignLHS(code, symTab, typeEnv))
    asserts
}

toEigenAssignLHS <- function(code, symTab, typeEnv, NoEigenizeMap = FALSE) {
    ## In original size processing, NoEigenizeMap exists to be set as TRUE from sizeInsertIntermediate.
    if(identical(code$type, 'unknown')) {
        code$toEigenize <- 'no'
        return(list())
    }
    assert <- list()
    LHS <- code$args[[1]]
    RHS <- code$args[[2]]
    
    if(LHS$toEigenize == 'yes') { ## unexpected?
        code$toEigenize <- 'yes'
    } else {
        code$toEigenize <-
            if(inherits(RHS, 'exprClass')) {
                if(RHS$toEigenize == 'no') 'no'
                else {
                    if(RHS$toEigenize == 'unknown') 'no'
                    else {
                        if(RHS$toEigenize != 'yes' &
                           (!(LHS$name %in% c('eigenBlock',
                                              'diagonal',
                                              'coeffSetter'))) &
                           (RHS$nDim == 0 |
                            RHS$isName |
                            (RHS$name == 'map' & NoEigenizeMap)))
                            'no' ## if it is scalar or is just a name or a map, we will do it via NimArr operator= .  Used to have "| RHS$name == 'map'", but this allowed X[1:3] <- X[2:4], which requires eigen, with eval triggered, to get right
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
    ## if LHS is set to eigenize 'yes set code to Eigenize 'yes'
    
    if(code$toEigenize == 'yes') { ## this would make more sense in eigenize_assign
        ## generate setSize(LHS, ...) where ... are dimension expressions
        if(length(RHS$nDim) == 0) {
            message("confused about trying to eigenize something with nDim = 0")
            browser()
        }
        if(RHS$nDim > 0) {
            if(!(RHS$name %in% setSizeNotNeededOperators)) {
                if(LHS$isName | LHS$name == "nfVar") {
                    assert <-
                        substitute(setSize(LHS),
                                   list(LHS = nimbleGeneralParseDeparse(LHS)))
                    for(i in seq_along(RHS$sizeExprs)) { ## N.B. RHS$sizeExprs may be modified by sizeAssignAfterRecursing
                        test <- try(
                            assert[[i + 2]] <- RHS$sizeExprs[[i]]
                        )
                        if(inherits(test, 'try-error')) stop(paste0('In toEigenAssignLHS: Error in assert[[i + 2]] <- RHS$sizeExprs[[i]] for i = ', i), call. = FALSE)
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
                                    newExpr <-
                                        insertExprClassLayer(code,
                                                             1,
                                                             'asRow',
                                                             type = LHS$type)
                                    newExpr$sizeExprs <- RHS$sizeExprs 
                                    newExpr$type <- LHS$type
                                    newExpr$nDim <- RHS$nDim
                                    if(!is.numeric(LHS$sizeExprs[[1]]) |
                                       !is.numeric(RHS$sizeExprs[[2]])) {
                                        assertMessage <-
                                            paste0("Run-time size error: expected ",
                                                   deparse(LHS$sizeExprs[[1]]), " == ",
                                                   deparse(RHS$sizeExprs[[2]]))
                                        thisAssert <-
                                            identityAssert(LHS$sizeExprs[[1]],
                                                           RHS$sizeExprs[[2]],
                                                           assertMessage)
                                        if(!is.null(thisAssert))
                                            assert[[length(assert) + 1]] <-
                                                thisAssert 
                                    } else {
                                        if(LHS$sizeExprs[[1]] != RHS$sizeExprs[[2]])
                                            stop(exprClassProcessingErrorMsg(code,
                                                                             paste0('In sizeAssignAfterRecursing: Fixed size mismatch.')),
                                                 call. = FALSE)
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
            if(RHS$name == 'map')
                assert <- c(assert,
                            toEigenInsertIntermediate(code, 2, symTab, typeEnv) )
        }
        if(inherits(LHS, 'exprClass')) {
                                        # ditto
            if(LHS$name == 'map')
                assert <- c(assert,
                            toEigenInsertIntermediate(code, 1, symTab, typeEnv) )
        }
    }
    assert
    
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
    asserts <- recurseSetToEigenize(code, symTab, typeEnv)
    a1 <- code$args[[1]] 
    a2 <- code$args[[2]]
    if(inherits(a1, 'exprClass')) {
        if(a1$toEigenize == 'no') {
            asserts <- c(asserts,
                         toEigenInsertIntermediate(code, 1, symTab, typeEnv))
            a1 <- code$args[[1]]
        }
        a1DropNdim <- length(dropSingleSizes(a1$sizeExprs)$sizeExprs)
        a1toEigenize <- a1$toEigenize
    } else {
        a1DropNdim <- 0
        a1toEigenize <- 'maybe'
    }
    if(inherits(a2, 'exprClass')) {
        if(a2$toEigenize == 'no') {
            asserts <- c(asserts,
                         toEigenInsertIntermediate(code, 2, symTab, typeEnv))
            a2 <- code$args[[2]]
        }
        a2DropNdim <- length(dropSingleSizes(a2$sizeExprs)$sizeExprs)
        a2toEigenize <- a2$toEigenize
    } else {
        a2DropNdim <- 0
        a2toEigenize <- 'maybe'
    }
    forceYesEigenize <-
        identical(a1toEigenize, 'yes') |
        identical(a2toEigenize, 'yes')
    code$toEigenize <- if(a1DropNdim == 0 & a2DropNdim == 0)
                           if(forceYesEigenize)
                               'yes'
                           else
                               'maybe'
                       else 'yes'
    asserts
}

toEigenUnaryCwise <- function(code, symTab, typeEnv) {
    asserts <- recurseSetSizes(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    ## lifting rule: lift if not eigenizable
    if(inherits(a1, 'exprClass')) {
        if(a1$toEigenize == 'no') {
            asserts <- c(asserts,
                         toEigenInsertIntermediate(code, 1, symTab, typeEnv))
            a1 <- code$args[[1]]
        }
    }
    ## propagation rule: yes or maybe
    code$toEigenize <- if(code$nDim > 0) 'yes' else 'maybe'
    asserts
}

toEigenMatrixMult <- function(code, symTab, typeEnv) {
    asserts <- recurseSetToEigenize(code, symTab, typeEnv)
    a1 <- code$args[[1]]
    a2 <- code$args[[2]]
    ## POSSIBLE DIFFERENT BEHAVIOR:
    ## sizeMatrixMult may insert an asRow() or asCol()
    ## It would do so *after* sizeInsertIntermediate
    ## But now that would happen *before* sizeInsertIntermediate
    ## Options: (1) The new code made be just as valid.  (2) We could make toEigenInsertIntermediate be smart about asRow or asCol
    if(a1$toEigenize == 'no') {
        asserts <- c(asserts, toEigenInsertIntermediate(code, 1, symTab, typeEnv))
    }
    if(a2$toEigenize == 'no') {
        asserts <- c(asserts, toEigenInsertIntermediate(code, 2, symTab, typeEnv))
    }
    code$toEigenize <- 'yes'
    asserts
}

toEigenUnaryReduction <- function(code, symTab, typeEnv) {
    asserts <- recurseSetToEigenize(code, symTab, typeEnv)
    if(inherits(code$args[[1]], 'exprClass')) {
        if(!code$args[[1]]$isName) {
            if(code$args[[1]]$toEigenize == 'no') {
                asserts <- c(asserts, toEigenInsertIntermediate(code, 1, symTab, typeEnv))
            }
        }
    }

    code$toEigenize <- 'yes'

    if(!(code$caller$name %in% c('{','<-','<<-','='))) {
        asserts <- c(asserts, toEigenInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv))
    }
    
    if(length(asserts) == 0) NULL else asserts
}
