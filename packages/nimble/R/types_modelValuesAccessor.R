
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

modelValuesAccessorVector <- function(mv, nodeNames, logProb = FALSE, logProbOnly = FALSE) {
    ans <- list(mv, substitute(nodeNames), logProb, parent.frame(), logProbOnly)
    class(ans) <- c('modelValuesAccessorVector', 'valuesAccessorVector')
    ans
}

makeSetCodeFromAccessorVector <- function(accessorVector) {
    nodeNames <- eval(accessorVector[[2]], envir = accessorVector[[4]])
    modelOrModelValues <- accessorVector[[1]]
    if(accessorVector[[3]]) {## logProb == TRUE
        isLogProbName <- grepl('logProb_', nodeNames)
        if(accessorVector[[5]]) ## logProbOnly == TRUE
            nodeNames <- c(modelOrModelValues$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
        else
            nodeNames <- c(nodeNames, modelOrModelValues$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
    }
    if(inherits(accessorVector, 'modelValuesAccessorVector'))
        setCode <- lapply(nodeNames, function(nn) {
            temp <- parse(text = nn, keep.source = FALSE)[[1]]
            if(is.name(temp)) return(substitute(sourceToObject[[B]][[rowTo]] <- oneValue, list(B = as.character(temp))))
            temp[[2]] <- substitute(sourceToObject[[B]][[rowTo]], list(B = as.character(temp[[2]])))
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
        if(accessorVector[[5]]) ## logProbOnly == TRUE
            nodeNames <- c(modelOrModelValues$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
        else            
            nodeNames <- c(nodeNames, modelOrModelValues$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
    }
    if(inherits(accessorVector, 'modelValuesAccessorVector'))
        getCode <- lapply(nodeNames, function(nn) {
            temp <- parse(text = nn, keep.source = FALSE)[[1]]
            if(is.name(temp)) return(substitute(sourceFromObject[[B]][[row]], list(B = as.character(temp))))
            temp[[2]] <- substitute(sourceFromObject$B[[row]], list(B = as.character(temp[[2]])))
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

makeMapInfoFromAccessorVectorFaster <- function(accessorVector ) {
    nodeNames <- eval(accessorVector[[2]], envir = accessorVector[[4]])
    sourceObject <- accessorVector[[1]] ## a model or modelValues
    
    if(accessorVector[[3]]) {## logProb == TRUE
        ## efficiency note for posterity:
        ## grepl is inefficient when used extensively.
        ## Using fixed = TRUE is much more efficient.
        ## In this case, substr comparison is even more efficient
        ## isLogProbName <- grepl('logProb_', nodeNames)
        isLogProbName <- substr(nodeNames, 1, 8) == 'logProb_'
        if(accessorVector[[5]]) ## logProbOnly == TRUE
            nodeNames <- c(sourceObject$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
        else
            nodeNames <- c(nodeNames, sourceObject$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
    }
    varNames <- .Call(parseVar, nodeNames)
    symTab <- sourceObject$getSymbolTable()
    varSizesAndNDims2 <- symTab$makeDimAndSizeList(varNames)    
    varSizesAndNDims <- varSizesAndNDims2
    list(nodeNames, varSizesAndNDims)
}

makeMapInfoFromAccessorVector <- function(accessorVector ) {
    length <- 0
    nodeNames <- eval(accessorVector[[2]], envir = accessorVector[[4]])
    sourceObject <- accessorVector[[1]] ## a model or modelValues
    
    if(accessorVector[[3]]) {## logProb == TRUE
        isLogProbName <- grepl('logProb_', nodeNames)
        if(accessorVector[[5]]) ## logProbOnly == TRUE
            nodeNames <- c(sourceObject$modelDef$nodeName2LogProbName(nodeNames[!isLogProbName]))
        else
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

    ## Might need to go back and check if this was ever called.
    ## if(length(accessorVector) > 4) { ## set the length variable in the calling (setup) environment if needed
    ##     assign(accessorVector[[5]], length, envir = accessorVector[[4]]) 
    ## }
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

