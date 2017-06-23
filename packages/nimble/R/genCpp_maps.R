## Code for maps in the compilation system

## New map format will be
## map( target, map nDim, offsetExpr, sizeExpr1, ..., sizeExprNDim, strideExpr1, ..., strideExprNDim )
## e.g. a simple map for X[1:3, 2, 4:8] would give
## map (X, 2, 1 * strides(X)[1] + 2 * strides(X)[2] + 4 * strides(X)[3], 3, 5, strides(X)[1], strides(X)[3] )

## we compound maps "on the fly", meaning there is no reason to keep them separated.


## test <- RparseTree2ExprClasses(quote(x[i:4, 2, 2:j]))
## test$args[[1]]$sizeExprs <- list(quote( dim(x)[1] ), quote(dim(x)[2]), quote(dim(x)[3]) )
## test$sizeExprs <- list(quote(4-i +1), quote(j-2+1))
## test$nDim <- 2

## ans <- makeMapExprFromBrackets(test)
## test2 <- RparseTree2ExprClasses(quote(y[k:3, l]))
## setArg(test2, 1, ans)
## debug(makeMapExprFromBrackets)
## ans2 <- makeMapExprFromBrackets(test2)

makeStrideRexprs <- function(varRexpr, nDim) {
    ans <- vector('list', nDim)
    for(i in 1:nDim) {
        ans[[i]] <- substitute(strides( VE )[ ND ], list(VE = varRexpr, ND = i) )
    }
    ans
}

makeOffsetRexpr <- function(firstIndexRexprs, sourceStrideRexprs) {
    nDim <- length(firstIndexRexprs)
    if(nDim == 0) return(quote(0))
    if(nDim != length(sourceStrideRexprs)) stop('Error in makeOffsetRexpr: nDims do not match.', call. = FALSE)
    sumExprs <- vector('list', nDim)
    for(i in 1:nDim) {
        if(is.numeric(firstIndexRexprs[[i]])) firstIndexRexprs[[i]] <- firstIndexRexprs[[i]]-1
        else firstIndexRexprs[[i]] <- substitute(fIR - 1, list(fIR = firstIndexRexprs[[i]]))
    }
    iFirst <- 0
    for(i in 1:nDim) {
        if(is.numeric(firstIndexRexprs[[i]])) {
            if(firstIndexRexprs[[i]] == 0) {
                sumExprs[i] <- list(NULL)
                next
            }
        }
        sumExprs[[i]] <- substitute(A * B, list(A = firstIndexRexprs[[i]], B = sourceStrideRexprs[[i]]))
        if(iFirst == 0) iFirst <- i
    }
    if(iFirst == 0) return(quote(0))
    
    allSums <- sumExprs[[iFirst]]
    if(nDim > iFirst) {
        for(i in (iFirst+1):nDim) {
            if(!is.null(sumExprs[[i]])) {
                newSums <- substitute(A + B, list(A = allSums, B = sumExprs[[i]]))
                allSums <- newSums
            }
        }
    }
    allSums
}


blockIndexInfo <- function(code) {
    ## targetIndexExprs begin at mapExpr arg 2
    nArgs <- length(code$args)
    nDim <- nArgs-1
    
    ## iterate and set up
    blockBool <- rep(FALSE, nDim)
    firstIndexRexprs <- vector('list', nDim)
    for(i in 2:nArgs) { ## 1 is the var
        im1 <- i-1
        if(inherits(code$args[[i]], 'exprClass')) {
            if(code$args[[i]]$name != ':') {
                if(code$args[[i]]$name != "") {
                    firstIndexRexprs[[im1]] <- parse(text = nimDeparse(code$args[[i]]), keep.source = FALSE)[[1]]
                } else {
                    ## It is a blank
                    blockBool[im1] <- TRUE
                    firstIndexRexprs[[im1]] <- 1
                }
            } else {
                ## It is a ":"
                blockBool[im1] <- TRUE
                if(!is.null(code$args[[i]]$sizeExprs))
                    if(length(code$args[[i]]$sizeExprs[[1]])==1)
                        if(is.numeric(code$args[[i]]$sizeExprs[[1]]))
                            if(code$args[[i]]$sizeExprs[[1]]==1)
                                blockBool[im1] <- FALSE
                firstIndexRexprs[[im1]] <- parse(text = nimDeparse(code$args[[i]]$args[[1]]), keep.source = FALSE)[[1]]
            }
        } else {
            firstIndexRexprs[[im1]] <- code$args[[i]]
        }
    }
    list(firstIndexRexprs = firstIndexRexprs, blockBool = blockBool)
}

addTransposeIfNeededForNonSeqBlock <- function(code, drop = TRUE) {
    ## somewhat like makeEigenBlockExprFromBrackets, but the only
    ## purpose is to add a transpose step for the case X[ scalar, vector ]
    if(!drop) return(code);
    nArgs <- length(code$args)
    nDim <- nArgs-1
    blockBool <- rep(FALSE, nDim)
    for(i in 1:nDim) {
        if(inherits(code$args[[i+1]], 'exprClass')) {
            if(code$args[[i+1]]$nDim > 0) blockBool[i] <- TRUE
        }
    }
    if(identical(blockBool, c(FALSE, TRUE))) {
        code <- insertEigenTranspose(code)
    }
    code
}

insertEigenTranspose <- function(code) {
    newExpr2 <- exprClass$new(isName = FALSE, isCall = TRUE, isAssign = FALSE, name = 'eigTranspose', args = list(1))
    newExpr2$sizeExprs <- c(code$sizeExprs[2], code$sizeExprs[1])
    newExpr2$toEigenize <- 'yes' 
    newExpr2$nDim <-  code$nDim
    newExpr2$type <- code$type
    setArg(newExpr2, 1, code)
    newExpr2
}

makeEigenBlockExprFromBrackets <- function(code, drop = TRUE) {
    thisBlockIndexInfo <- blockIndexInfo(code)
    blockBool <- thisBlockIndexInfo$blockBool
    firstIndexRexprs <- thisBlockIndexInfo$firstIndexRexprs

    nArgs <- length(code$args)
    nDim <- nArgs-1
    for(i in 1:nDim) {
        if(is.numeric(firstIndexRexprs[[i]])) firstIndexRexprs[[i]] <- firstIndexRexprs[[i]]-1
        else firstIndexRexprs[[i]] <- substitute(fIR - 1, list(fIR = firstIndexRexprs[[i]]))
    }

    needTranspose <- FALSE
    if(length(blockBool)==1) {
        blockBool <- c(blockBool, FALSE)
        firstIndexRexprs[[2]] <- 0
    }
    if(identical(blockBool, c(TRUE, TRUE))) {
        P <- code$sizeExprs[[1]]
        Q <- code$sizeExprs[[2]]
    } else 
        if(identical(blockBool, c(TRUE, FALSE))) {
            P <- code$sizeExprs[[1]]
            Q <- 1
        } else
            if(identical(blockBool, c(FALSE, TRUE))) {
                P <- 1
                Q <- code$sizeExprs[[ if(drop) 1 else 2 ]]
                if(drop) needTranspose <- TRUE ## This is a case like A[3, 2:4] where we need to make Eigen handle it like a 3x1 matrix implementing a nimble vector internally, not 1x3
            } else
                if(identical(blockBool, c(FALSE, FALSE))) {
                    P <- 1
                    Q <- 1
                    if(drop) warning("Possible error handling a 1x1 block with drop = TRUE")
                }


    ## this will make the I J P Q come out as properly arranged exprClasses
    newExpr <- RparseTree2ExprClasses( substitute(eigenBlock(1, I, J, P, Q), list(I = firstIndexRexprs[[1]], J = firstIndexRexprs[[2]], P = P, Q = Q)))
    ## but the target variable already is an exprClass and so needs to be inserted properly separately:
    setArg(newExpr, 1, code$args[[1]])
    for(i in 2:5) {
        if(inherits(newExpr$args[[i]], 'exprClass')) {
            newExpr$args[[i]]$nDim <- 0
            newExpr$args[[i]]$type <- 'integer'
            newExpr$args[[i]]$sizeExprs <- list(1)
            newExpr$args[[i]]$toEigenize <- 'maybe'
        }
    }
    if(needTranspose) { 
        newExpr$nDim <- code$nDim
        newExpr$type <- code$type
        newExpr$sizeExprs <- code$sizeExprs
        newExpr$toEigenize <- 'yes' ## not really needed since will be called from eigenization
        newExpr <- insertEigenTranspose(newExpr)
    }
    newExpr
}

## this is used by sizeIndexingBracket when it hits a need for a map
makeMapExprFromBrackets <- function(code, drop = TRUE) {
    ## code nDim, type and sizeExprs have already been set, and toEigenize will be set to 'maybe'
    if(code$args[[1]]$name == 'map') {
        sourceVarName <- code$args[[1]]$args[[1]]
        sourceVarExpr <- as.name(sourceVarName)
        sourceOffsetRexpr <- code$args[[1]]$args[[3]]
        sourceSizeExprs <- code$args[[1]]$args[[4]]
        sourceNdim <- length(sourceSizeExprs)
        sourceStrideRexprs <- code$args[[1]]$args[[5]]
    } else {
        sourceVarName <- code$args[[1]]$name
        if(!code$args[[1]]$isName) writeLines(paste0('Watch out, in makeMapExprFromBrackets for ', nimDeparse(code), ', there is an expression instead of a name.'))
        sourceVarExpr <- as.name(sourceVarName)
        sourceOffsetRexpr <- quote(0)
        sourceSizeExprs <- code$args[[1]]$sizeExprs
        sourceNdim <- length(sourceSizeExprs)
        sourceStrideRexprs <- makeStrideRexprs(sourceVarExpr, sourceNdim)
    }

    thisBlockIndexInfo <- blockIndexInfo(code)
    blockBool <- thisBlockIndexInfo$blockBool
    firstIndexRexprs <- thisBlockIndexInfo$firstIndexRexprs
    if(identical(sourceOffsetRexpr, 0)) {
        targetOffsetRexpr <- makeOffsetRexpr(firstIndexRexprs, sourceStrideRexprs)
    } else {
        targetOffsetRexpr <- substitute(A + B, list(A = sourceOffsetRexpr, B = makeOffsetRexpr(firstIndexRexprs, sourceStrideRexprs)))
    }
    targetSizeExprs <- code$sizeExprs
    targetStrideRexprs <- if(drop) sourceStrideRexprs[blockBool] else sourceStrideRexprs
    targetNdim <- length(targetSizeExprs)

    ## this is an unusual exprClass object because args is just a regular list.  Its elements are not exprClass objects
    newExpr <- exprClass$new(isName = FALSE, isCall = TRUE, isAssign = FALSE, name = 'map', args = list(sourceVarName, targetNdim, targetOffsetRexpr, targetSizeExprs, targetStrideRexprs))

    newExpr

}

## This is used to build the setMap expression for a NimArr
nimArrMapExpr <- function(code, symTab, typeEnv, newName) {
    mapName <- newName
    varName <- code$args[[1]]
    needStartOffset <- !is.null(typeEnv$passedArgumentNames[[varName]])
    targetSym <- symTab$getSymbolObject(varName, TRUE)
    if(targetSym$nDim == 0) {
        writeLines("Strange, in nimArrMap, there is a case of nDim == 0")
        browser()
    }
    if(!symTab$symbolExists(mapName, TRUE)) {
        newSym <- symbolBasic(name = mapName,
                              nDim = code$nDim,
                              type = targetSym$type)
        symTab$addSymbol(newSym)
    }
    if(length(code$sizeExprs) != length(code$args[[5]])) {
        stop('Error, length(code$sizeExprs) != length(code$args[[5]]), which should be the strideRexprs')
    }
    sizeExprs <- code$args[[4]]
    strides <- code$args[[5]]

    if(needStartOffset) offsetRexpr <- substitute(getOffset(A) + chainedCall(template(static_cast, int), B), list(A =  as.name(varName), B = code$args[[3]]))
    else offsetRexpr <- substitute( chainedCall(template(static_cast, int), OE), list(OE = code$args[[3]]) )
    
    ## Need to build up this expression
    ans <- list(as.name('setMap'), as.name(newName), as.name(varName), offsetRexpr) ## varName used to be targetSym$name. Should be identical
    ans <- c(ans, strides, sizeExprs)
    ans <- as.call(ans)
    return(ans)
}

## This is used to build the AssignEigenMap expression
eigenizeNameStrided <- function(code, symTab, typeEnv, workEnv) {
    varName <- code$args[[1]]
    
    ## this is a map on a passed argument.  It may be itself be a map, so the offset from it will be needed
    needStartOffset <- !is.null(typeEnv$passedArgumentNames[[varName]])

    EigenName <- paste0(Rname2CppName(makeEigenName(varName)), IntermLabelMaker())
    targetSym <- symTab$getSymbolObject(varName, TRUE)
    if(targetSym$nDim == 0) {
        writeLines("Strange, in eigenizeNameStrided, there is a case of nDim == 0")
        browser()
    }
    
    mapSizeExprs <- code$args[[4]]
    mapStrideExprs <- code$args[[5]]
    if(length(code$args[[4]]) != length(code$args[[5]])) {
        stop('Error, length(code$args[[4]]) != length(code$args[[5]]), which give the stride exprs')
    }
    if(length(mapSizeExprs) > 2) {stop('Error, cannot eigenize a map of dimensions > 2')}
    if(length(mapSizeExprs) == 2) {
        nrowExpr <- mapSizeExprs[[1]]
        ncolExpr <- mapSizeExprs[[2]]
        strides <- mapStrideExprs
    }
    if(length(mapSizeExprs) == 1) {
        if(code$caller$name == 'asRow') {
            nrowExpr <- 1
            ncolExpr <- mapSizeExprs[[1]]
            strides <- c(list(0), mapStrideExprs)
            EigenName <- paste0(EigenName,'_asRow')
        } else {
            ## default to column
            nrowExpr <- mapSizeExprs[[1]]
            ncolExpr <- 1
            strides <- c(mapStrideExprs, list(0))
        }
    }
    if(length(mapSizeExprs) == 0) {
        nrowExpr <- 1
        nColExpr <- 1
        strides <- list(0,0)
    }

    thisMapAlreadySet <- FALSE
    if(!is.null(workEnv[['OnLHSnow']])) { ## this is the var on the LHS
        if(!is.null(workEnv[['LHSeigenName']])) stop(paste0('Error for map of ', varName, '. LHSeigenName already exists'), call. = FALSE)
        workEnv$LHSeigenName <- list(EigenName = EigenName, targetVar = varName)
        workEnv[[EigenName]] <- TRUE
    } else { ## This is on the RHS
        if(EigenName %in% ls(workEnv) ) {
            thisMapAlreadySet <- TRUE
        } else {
            workEnv[[EigenName]] <- TRUE
        }
        if(!is.null(workEnv[['LHSeigenName']])) { ## There was a LHS (there may not be for nimPrint(x), for example)
            alreadyAliased <- !is.null(workEnv[['mustAddEigenEval']])
            if(!alreadyAliased) {
                aliasRisk <- !is.null(workEnv[['aliasRisk']])
                if(varName == workEnv[['LHSeigenName']]$targetVar) { ## this uses the same targetVar as the LHS
                    if(aliasRisk || EigenName != workEnv[['LHSeigenName']]$EigenName) {
                        workEnv[['mustAddEigenEval']] <- TRUE
                    }
                }
            }
        }
    }

    
    if(!symTab$symbolExists(EigenName, TRUE)) {
        if(thisMapAlreadySet) warning(paste0('Weird, it looks like a strided Eigen map for ', varName, 'was already set but the symbol did not exist.'), call. = FALSE)
        newSym <- symbolEigenMap(name = EigenName,
                                 eigMatrix = TRUE, ## default to matrix
                                 type = targetSym$type,
                                 strides = as.numeric(c(NA, NA))) ## Not sure strides are really used from this object
        symTab$addSymbol(newSym)
    }
    
    newExprClass <- RparseTree2ExprClasses(as.name(EigenName) )
    newExprClass$eigMatrix <- TRUE
    newExprClass$sizeExprs <- code$sizeExprs
    newExprClass$nDim <- code$nDim
    newExprClass$type <- code$type
    setArg(code$caller, code$callerArgID, newExprClass)

    if(!thisMapAlreadySet) {
        if(needStartOffset) offsetRexpr <- substitute(getOffset(A) + chainedCall(template(static_cast, int), B), list(A =  as.name(varName), B = code$args[[3]]))
        else offsetRexpr <- substitute(chainedCall(template(static_cast, int), B), list(B = code$args[[3]]))
        return(RparseTree2ExprClasses(
            EigenNewExpr(EigenName, varName, offsetRexpr, makeEigenTypeLabel(TRUE, targetSym$type),
                         nrowExpr, ncolExpr, strides = rev(strides))) ## Eigen takes strides as (outer, inner). 
               ) ## varName was targetSym$name.  They should be identical
    } else {
        return(NULL)
    }
}

