
## CASE 3: modelValuesAccessorVector
## copy(modelValues, xxx, from.row=i, nodes) becomes:
## modelValues_nodes_accessors <- modelValuesAccessorVector(modelValues, nodes)
## copy(modelValues_nodes_accessors$row(i), xxx)
## modelValuesAccessorVector$getAccessors() returns a list of the modelValuesAccessor objects



copierVector <- function(accessFrom_name, accessTo_name, isFromMV, isToMV) {
    ans <- list(deparse(substitute(accessFrom_name)), deparse(substitute(accessTo_name)), isFromMV, isToMV)
    class(ans) <- 'copierVectorClass'
    ans
}

modelValuesAccessorVector <- function(mv, nodeNames, logProb = FALSE) {
    ans <- list(mv, substitute(nodeNames), logProb, parent.frame())
    class(ans) <- c('modelValuesAccessorVector', 'valuesAccessorVector')
    ans
}

makeSetCodeFromAccessorVector <- function(accessorVector) {
    nodeNames <- eval(accessorVector[[2]], envir = accessorVector[[4]])
    modelOrModelValues <- accessorVector[[1]]
    if(accessorVector[[3]]) {## logProb == TRUE
        isLogProbName <- grepl('logProb_', nodeNames)
        nodeNames <- c(nodeNames, modelOrModelValues$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
    }
    if(inherits(accessorVector, 'modelValuesAccessorVector'))
        setCode <- lapply(nodeNames, function(nn) {
            temp <- parse(text = nn, keep.source = FALSE)[[1]]
            if(is.name(temp)) return(substitute(sourceToObject$B[[rowTo]] <- oneValue, list(B = temp)))
            temp[[2]] <- substitute(sourceToObject$B[[rowTo]], list(B = temp[[2]]))
            substitute(A <- oneValue, list(A = temp))
        })
    else
        setCode <- lapply(nodeNames, function(nn) {
            temp <- parse(text = nn, keep.source = FALSE)[[1]]
            if(is.name(temp)) return(substitute(sourceToObject$B <- oneValue, list(B = temp)))
            temp[[2]] <- substitute(sourceToObject$B, list(B = temp[[2]]))
            substitute(A <- oneValue, list(A = temp))
        })
    setCode
}

makeGetCodeFromAccessorVector <- function(accessorVector) {
    nodeNames <- eval(accessorVector[[2]], envir = accessorVector[[4]])
    modelOrModelValues <- accessorVector[[1]]
    if(accessorVector[[3]]) {## logProb == TRUE
        isLogProbName <- grepl('logProb_', nodeNames)
        nodeNames <- c(nodeNames, modelOrModelValues$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
    }
    if(inherits(accessorVector, 'modelValuesAccessorVector'))
        getCode <- lapply(nodeNames, function(nn) {
            temp <- parse(text = nn, keep.source = FALSE)[[1]]
            if(is.name(temp)) return(substitute(sourceFromObject$B[[row]], list(B = temp)))
            temp[[2]] <- substitute(sourceFromObject$B[[row]], list(B = temp[[2]]))
            temp
        })
    else
        getCode <- lapply(nodeNames, function(nn) {
            temp <- parse(text = nn, keep.source = FALSE)[[1]]
            if(is.name(temp)) return(substitute(sourceFromObject$B, list(B = temp)))
            temp[[2]] <- substitute(sourceFromObject$B, list(B = temp[[2]]))
            temp
        })
    getCode
}

## getValues.modelValuesAccessorVector <- function(accessorVector, i, row = NA) {
##     if(is.na(row)) stop('Can no longer use row = NA in getValues.modelValuesAccessorVector')
##     nodeName <- eval(accessorVector[[2]], envir = accessorVector[[4]])[i]
##     nodeCode <- parse(text = nodeName, keep.source = FALSE)[[1]]
##     sourceObject <- accessorVector[[1]]
##     if(is.name(nodeCode)) return(sourceObject[[nodeName]][[row]])
##     nodeCode[[2]] <- substitute(sourceObject[[VN]][[row]], list(VN = as.character(nodeCode[[2]]) ) )
##     eval(nodeCode)
## }
## setValues.modelValuesAccessorVector <- function(accessorVector, i, vals, row = NA) {
##     if(is.na(row)) stop('Can no longer use row = NA in setValues.modelValuesAccessorVector')
##     nodeName <- eval(accessorVector[[2]], envir = accessorVector[[4]])[i]
##     nodeCode <- parse(text = nodeName, keep.source = FALSE)[[1]]
##     sourceObject <- accessorVector[[1]]
##     if(is.name(nodeCode)) nodeCode <- quote(sourceObject[[nodeName]][[row]])
##     nodeCode[[2]] <- substitute(sourceObject[[VN]][[row]], list(VN = as.character(nodeCode[[2]]) ) )
##     eval(substitute(A <- vals), list(A = nodeCode))
## }

## A wrapper on one lapply to get clearer Rprof results
## SINGLE_LAPPLY_FOR_PROFILING <- function( varNames, symTab ) {
##     ans <- lapply(varNames, function(x) {symObj <- symTab$getSymbolObject(x); list(symObj$size, symObj$nDim)})
##     ans
## }

## A wrapper on one member function to get clearer Rprof results
## nodeName2LogProbName_FOR_PROFILING <- function(md, arg1) {
##     ans <- md$nodeName2LogProbName(arg1)
##     ans
## }


## 5 and 2 are the big ones

timingFunction1 <- function(a, envir) {
    eval(a, envir = envir)
}

timingFunction2 <- function(a, b) {
    a$modelDef$nodeName2LogProbName(b)
}

timingFunction3 <- function(a) .Call('parseVar', a)

timingFunction4 <- function(a) a$getSymbolTable()

timingFunction6 <- function(symTab, x) symTab$getSymbolObject(x)

timingFunction5 <- function(varNames, symTab) lapply(varNames, function(x) {symObj <- timingFunction6(symTab, x); list(symObj$size, symObj$nDim)})

makeMapInfoFromAccessorVectorFaster <- function(accessorVector ) {
##    length <- 0
    nodeNames <- timingFunction1(accessorVector[[2]], envir = accessorVector[[4]]) ## eval(accessorVector[[2]], envir = accessorVector[[4]])
    sourceObject <- accessorVector[[1]] ## a model or modelValues
    
    if(accessorVector[[3]]) {## logProb == TRUE
        isLogProbName <- grepl('logProb_', nodeNames)
        nodeNames <- c(nodeNames, timingFunction2(sourceObject, nodeNames[!isLogProbName])) ##sourceObject$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
##        nodeNames <- c(nodeNames, nodeName2LogProbName_FOR_PROFILING(sourceObject$modelDef, nodeNames[!isLogProbName]))
    }

## time these also
    
    varNames <- timingFunction3(nodeNames) ##.Call('parseVar', nodeNames)
    symTab <- timingFunction4(sourceObject) ##sourceObject$getSymbolTable()
    varSizesAndNDims <- timingFunction5(varNames, symTab) ## lapply(varNames, function(x) {symObj <- symTab$getSymbolObject(x); list(symObj$size, symObj$nDim)})
    varSizesAndNDims2 <- symTab$makeDimAndSizeList(varNames)
    test <- varSizesAndNDims
    if(length(varNames)>0) names(test) <- varNames
    else {
        names(test) <- NULL
        names(varSizesAndNDims2) <- NULL
    }
    if(!identical(test, varSizesAndNDims2)) browser()

    varSizesAndNDims <- varSizesAndNDims2
##    varSizesAndNDims <- SINGLE_LAPPLY_FOR_PROFILING(varNames, symTab)
    
    list(nodeNames, varSizesAndNDims)
    ## mapInfo <- lapply(nodeNames, function(z) {
    ##     x <- parse(text = z, keep.source = FALSE)[[1]]
    ##     varAndIndices <- getVarAndIndices(x)
    ##     varName <- as.character(varAndIndices$varName)
    ##     varSym <- sourceObject$getSymbolTable()$getSymbolObject(varName) ## previously from model$getVarInfo(varName)
    ##     ans <- varAndIndices2mapParts(varAndIndices, varSym$size, varSym$nDim)
    ##     ans$varName <- varName
    ##     anslength <- prod(ans$sizes)
    ##     length <<- length + anslength
    ##     ans$singleton <- anslength == 1
    ##     ans$length <- anslength ## putting it last so it doesn't mess up current C++ code
    ##     ans
    ## }) ## list elements will be offset, sizes, strides, varName, and singleton in that order.  Any changes must be propogated to C++

    ## if(length(accessorVector) > 4) { ## set the length variable in the calling (setup) environment if needed
    ##     assign(accessorVector[[5]], length, envir = accessorVector[[4]]) 
    ## }
    ## mapInfo
}

makeMapInfoFromAccessorVector <- function(accessorVector ) {
    length <- 0
    nodeNames <- eval(accessorVector[[2]], envir = accessorVector[[4]])
    sourceObject <- accessorVector[[1]] ## a model or modelValues
    
    if(accessorVector[[3]]) {## logProb == TRUE
        isLogProbName <- grepl('logProb_', nodeNames)
        nodeNames <- c(nodeNames, sourceObject$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
    }
    
    mapInfo <- lapply(nodeNames, function(z) {
        x <- parse(text = z, keep.source = FALSE)[[1]]
        varAndIndices <- getVarAndIndices(x)
        varName <- as.character(varAndIndices$varName)
        varSym <- sourceObject$getSymbolTable()$getSymbolObject(varName) ## previously from model$getVarInfo(varName)
        ans <- varAndIndices2mapParts(varAndIndices, varSym$size, varSym$nDim)
        ans$varName <- varName
        anslength <- prod(ans$sizes)
        length <<- length + anslength
        ans$singleton <- anslength == 1
        ans$length <- anslength ## putting it last so it doesn't mess up current C++ code
        ans
    }) ## list elements will be offset, sizes, strides, varName, and singleton in that order.  Any changes must be propogated to C++

    if(length(accessorVector) > 4) { ## set the length variable in the calling (setup) environment if needed
        assign(accessorVector[[5]], length, envir = accessorVector[[4]]) 
    }
    mapInfo
}

valuesAccessorVector <- setRefClass( ## new implementation
    Class = 'valuesAccessorVector',
    fields = list(
        sourceObject = 'ANY',
        nodeNames = 'ANY',
        logProb = 'ANY', 
        length = 'ANY',
        code = 'ANY',
        accessCode = 'ANY',
        setCode = 'ANY',
        numAccessors = 'ANY',
        mapInfo = 'ANY'
    ),
    methods = list(
        initialize = function(modelOrModelValues, nodeNames, logProb = FALSE) { ## will need to work with IDs too
            sourceObject <<- modelOrModelValues
            inputNodeNames <- nodeNames
             if(!logProb)
                nodeNames <<- inputNodeNames
             else {
                 if(!inherits(modelOrModelValues$modelDef,'modelDefClass'))
                     stop('creating valuesAccessorVector with logProb = TRUE for an object that was not built from a model')
                 isLogProbName <- grepl('logProb_', inputNodeNames)
                 nodeNames <<- c(inputNodeNames, modelOrModelValues$modelDef$nodeName2LogProbName(inputNodeNames[!isLogProbName]))
            }
            logProb <<- logProb
            code <<- lapply(.self$nodeNames, function(x) parse(text = x, keep.source = FALSE)[[1]])
            numAccessors <<- length(code)
            length <<- NA
            mapInfo <<- NULL
        },
        makeMapInfo = function() {
            length <<- 0
            mapInfo <<- lapply(code, function(x) {
                varAndIndices <- getVarAndIndices(x)
                varName <- as.character(varAndIndices$varName)
                varSym <- sourceObject$getSymbolTable()$getSymbolObject(varName) ## previously from model$getVarInfo(varName)
                ans <- varAndIndices2mapParts(varAndIndices, varSym$size, varSym$nDim)
                ans$varName <- varName
                anslength <- prod(ans$sizes)
                length <<- length + anslength
                ans$singleton <- anslength == 1
                ans$length <- anslength ## putting it last so it doesn't mess up current C++ code
                ans
            }) ## list elements will be offset, sizes, strides, varName, and singleton in that order.  Any changes must be propogated to C++
            invisible(NULL)
        },
        getMapInfo = function() {
            if(is.null(mapInfo)) makeMapInfo()
            mapInfo
        },
        getLength = function(i = NA) {
            if(is.null(mapInfo)) makeMapInfo()
            if(is.na(i)) return(length)
            mapInfo[[i]]$length
        }
    )
    )


## I don't think this is used any more
## makeSingleModelValuesAccessor <- function(ids, modelValues){
## 	return(modelValuesAccessor(id = id, modelValues = modelValues))
##     }


## modelValuesAccessorVector <- setRefClass( ## new implementation
##    Class = 'modelValuesAccessorVector',
##     contains = 'valuesAccessorVector',
##     fields = list(row = 'ANY'),
##     methods = list(
##         initialize = function(...) {
##             callSuper(...)
##             row <<- as.integer(NA)
##             makeAccessAndSetCode()
##         },
##         makeAccessAndSetCode = function() {
##             accessCode <<- lapply(code, function(temp) {
##                 if(is.name(temp)) return(substitute(sourceObject$B[[localrow]], list(B = temp)))
##                 temp[[2]] <- substitute(sourceObject$B[[localrow]], list(B = temp[[2]]))
##                 temp
##             })
##             setCode <<- lapply(accessCode, function(x) substitute(A <- vals, list(A = x)))
##         },
##         setRow = function(row) {
##             row <<- row
##         },
##         getValues = function(i, row = NA) {
##             ## model version
##             localrow <- if(is.na(row)) .self$row else row
##             eval(accessCode[[i]])
##         },
##         setValues = function(i, vals, row = NA) {
##             ## model version
##             localrow <- if(is.na(row)) .self$row else row
##             eval(setCode[[i]])
##         }
##     ))

## length.modelValuesAccessorVector <- function(access)
## 	return(access$length)

## I don't think this is used any more
## modelValuesAccessor <- setRefClass(
##     Class = 'modelValuesAccessor',
##     fields = list(modelValues = 'ANY',
##     			   id = 'ANY'
##    #               var         = 'ANY', 		#'character',
##    #               first       = 'ANY', 		#'numeric',
##    #               last        = 'ANY', 		#'numeric',
##    #               length	  = 'ANY' 		#'numeric'
##     ),
##     methods = list(toStr = function() paste0(var, '[', first, ':', last, ']'),
##                    show  = function() cat(paste0(toStr(), '\n'))
##     )
## )
