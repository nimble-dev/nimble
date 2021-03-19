

## adds explicit indexing to variables (using symtab), then expands indexing present on any variables
nl_expandNodeNames <- function(nodeNames, symtab, env) {
    nodeNames <- nl_addIndicesToVariables(nodeNames, symtab)
    nodeNames <- unlist(lapply(nodeNames, function(node) if(is.indexed(node)) nl_expandNodeIndex(node, env) else node))
    if(is.null(nodeNames))   return(character(0))
    return(nodeNames)
}


## Expands variables into their fully indexed form, e.g., 'y' is expanded to 'y[1:10]', using information in the symbolTable
nl_addIndicesToVariables <- function(nodeNames, symtab) {
    for(i in seq_along(nodeNames)) {
        nodeName <- nodeNames[i]
        varName <- nl_getVarNameFromNodeName(nodeName)
        if(!(varName %in% symtab$getSymbolNames()))   stop('variable not in symbol table')
        if(!is.indexed(nodeName) && (symtab$getSymbolField(varName, 'nDim') > 0)) {    ## nodeName has no indexing, and has dimension > 0
            maxs <- symtab$getSymbolField(varName, 'size')
            mins <- rep(1, length(maxs))
            indexStuff <- paste(mins, maxs, sep=':', collapse = ', ')
            nodeNames[i] <- paste0(varName, '[', indexStuff, ']')
        }
    }
    return(nodeNames)
}


## This is the same as nl_ExpandNodeIndex, except it takes a nodeExpr instead of a node char string
nl_expandNodeIndexExpr <- function(nodeExpr, env = parent.frame()) {
    if(length(nodeExpr)==1)  if(is.name(nodeExpr)) return(as.character(nodeExpr)) else stop('node expression with only one element, but not a variable name')
    indexExprs <- nodeExpr[-c(1,2)]
    indexStrs <- lapply(indexExprs, function(ind) as.character(eval(ind, envir=env)))
    numInd <- length(indexStrs)
    indexStrsCombined <- indexStrs[[numInd]]
    if(numInd > 1) for(i in seq(numInd-1, 1, by = -1)) {
        indexStrsCombined <- outer(indexStrs[[i]], indexStrsCombined, paste, sep=', ') }
    nodeStr <- as.character(nodeExpr[[2]])
    expandedNodes <- paste0(nodeStr, paste0('[',indexStrsCombined,']'))
    return(expandedNodes)
}

## same as nl_vectorizedExpandNodeIndex but takes nodeExprs instead of nodes
nl_vectorizedExpandNodeIndexExprs <- function(nodeExprs, env = parent.frame()) {
    return(unlist(lapply(nodeExprs, function(n) nl_expandNodeIndexExpr(n, env))))
}

## Expands the indexing of a single node name string, e.g., 'x[1:3]' is expanded to c('x[1]', 'x[2]', 'x[3]')
nl_expandNodeIndex <- function(node, env = parent.frame()) {
    nodeExpr <- parse(text=node, keep.source = FALSE)[[1]]
    if(length(nodeExpr)==1)  if(is.name(nodeExpr)) return(as.character(nodeExpr)) else stop('node expression with only one element, but not a variable name')
    indexExprs <- nodeExpr[-c(1,2)]
    indexStrs <- lapply(indexExprs, function(ind) as.character(eval(ind, envir=env)))
    numInd <- length(indexStrs)
    indexStrsCombined <- indexStrs[[numInd]]
    if(numInd > 1) for(i in seq(numInd-1, 1, by = -1)) {
        indexStrsCombined <- outer(indexStrs[[i]], indexStrsCombined, paste, sep=', ') }
    nodeStr <- as.character(nodeExpr[[2]])
    expandedNodes <- paste0(nodeStr, paste0('[',indexStrsCombined,']'))
    return(expandedNodes)
}


nl_vectorizedExpandNodeIndex <- function(nodes, env = parent.frame()) {
    return(unlist(lapply(nodes, function(n) nl_expandNodeIndex(n, env))))
}


## checks that all nodeNames are present in model
nl_checkNodeNamesInModel <- function(model, nodeNames, determOnly = FALSE, stochOnly = FALSE) {
    # this function not used in package, but if it were, would be good to have error messages indicate which nodes are the issue; see below for analogous situation for nl_checkVarNamesInModel
    if(!all(nodeNames %in% model$getNodeNames()))      stop('all node names not in model')
    if(determOnly) if(!all(nodeNames %in% model$getNodeNames(determOnly = TRUE)))      stop('all node names are not deterministic')
    if(stochOnly)  if(!all(nodeNames %in% model$getNodeNames(stochOnly = TRUE)))       stop('all node names are not stochastic')
}

## checks that all varNames are present in model
nl_checkVarNamesInModel <- function(model, varNames) {
    found <- varNames %in% model$getVarNames(includeLogProb = TRUE)
    if(!all(found)) 
        stop('These variables are not in model: ', paste(varNames[!found], collapse = ','), '.')
}

nl_nodeVectorReadyNodes <- function(model, nodeNames, includeData = TRUE){
	GIDs <- model$modelDef$nodeName2GraphIDs(nodeNames)
	sortedIDs <- sort(GIDs)
	if(includeData == FALSE)
		sortedIDs = sortedIDs[!model$isDataFromGraphID(GIDs)]
	return(model$modelDef$maps$graphID_2_nodeFunctionName[sortedIDs])
}




## returns a list of variable names, and flat index ranges
nl_createVarsAndFlatIndexRanges <- function(nodeNames, symtab) {
    varAndIndexes <- lapply(nodeNames, function(nn) list(var = nl_getVarNameFromNodeName(nn), ind = if(is.indexed(nn)) unlist(as.list(parse(text=nn, keep.source = FALSE)[[1]][-c(1,2)])) else numeric(0)))
    varAndFlatIndexes <- lapply(varAndIndexes, function(vai) list(var = vai$var, flatIndex = nl_determineFlatIndex(ind=vai$ind, maxs=symtab$getSymbolField(vai$var, 'size'))))
    listsOfIndexes <- list()
    for(vafi in varAndFlatIndexes) listsOfIndexes[[vafi$var]] <- list(var = vafi$var, inds = c(listsOfIndexes[[vafi$var]]$inds, vafi$flatIndex))
    listsOfAggIndexes <- lapply(listsOfIndexes, function(loi) list(var=loi$var, agInd=nl_aggregateConsecutiveBlocks(loi$inds)))
    separatedListsOfAggIndexes <- list()
    for(loai in listsOfAggIndexes) for(ind in loai$agInd) separatedListsOfAggIndexes <- c(separatedListsOfAggIndexes, list(list(var=loai$var, ind=ind)))
    return(separatedListsOfAggIndexes)
}


## determines the (scalar) flatIndex, from numeric vectors of indicies, and max_indicies
nl_determineFlatIndex = function(ind, maxs) {
		
    if(length(ind) != length(maxs))   stop('number of indices and dimensions are different')
    if(any(ind > maxs))               stop('some indices exceed dimensions')
    if(any(ind < 1))                  stop('non-positive index')
    nDim <- length(ind)
    if(nDim==0) return(1)
    if(nDim==1) return(ind)
    max_prods <- 1
    for(iDim in 2:nDim) max_prods[iDim] <- prod(maxs[1:(iDim-1)])
    indShifted <- c(ind[1], ind[-1]-1)
    flatIndex <- sum(indShifted * max_prods)
    return(flatIndex)
}


## takes a vector of flat indices, e.g. c(1,2,3,6,7,100),
## returns a list of pairs: (index_start, index_finish), e.g. [[1]] c(1,3)  [[2]] c(6,7)  [[3]] c(100,100)
## If anyone has an implementation which beats O(n) linear time, then let's use it.
nl_aggregateConsecutiveBlocks <- function(ind) {
    if(length(ind) == 0)     return(list())
    indDiffLogical <- c(TRUE, diff(ind) != 1)
    aggregated <- list()
    for(i in seq_along(indDiffLogical)) {
        if(indDiffLogical[i])     aggregated[[length(aggregated)+1]] <- c(ind[i], ind[i])
        if(!indDiffLogical[i])    aggregated[[length(aggregated)]][2] <- ind[i]
    }
    return(aggregated)
}

nl_removeNodeNamesNotInSymbolTable <- function(nodeNames, st) {
    return(nodeNames[nl_getVarNameFromNodeName(nodeNames) %in% st$getSymbolNames()])
}


nl_getVarNameFromNodeName <- function(nodeName)    gsub('\\[.*', '', nodeName)

expandMVNames <- function(mv, varNames){
	sizeList = mv$sizes
	nodeNames = NA
	nodeIndex = 0 
	for(i in seq_along(varNames) ){
		baseName = varNames[i]
		dims = mv$symTab$getSymbolObject(baseName)$size
		if(length(dims) == 0){
			nodeNames[nodeIndex+1] = baseName
			nodeIndex <- nodeIndex + 1
			}
		else{
            mins <- rep(1, length(dims))
            indexStuff <- paste(mins, dims, sep=':', collapse = ', ')
            compactNodeNames <- paste0(baseName, '[', indexStuff, ']')
            expandedNodeNames <- nl_expandNodeIndex(compactNodeNames)
            nodeNames[nodeIndex + 1:length(expandedNodeNames)] <- expandedNodeNames
            nodeIndex = nodeIndex + length(expandedNodeNames)
		}
	}
	return(nodeNames)
}

as.matrix.modelValuesBaseClass <- function(x, varNames, ...){
	if(missing(varNames))
			varNames <- x$varNames
	nrows = getsize(x)
	flatNames = expandMVNames(x, varNames)
	ans <- matrix(0.1, nrow = nrows, ncol = length(flatNames))
	colIndex = 1
	for(i in seq_along(varNames)){
		.Call(fastMatrixInsert, ans, modelValuesElement2Matrix(x, varNames[i]) , as.integer(1), as.integer(colIndex) ) 		
		colIndex = colIndex + prod(x$sizes[[varNames[i] ]])
		}
	colnames(ans) <- flatNames
	return(ans)

}

as.matrix.CmodelValues <- function(x, varNames, ...){
	if(missing(varNames))
			varNames <- x$varNames
	nrows = getsize(x)
	flatNames = expandMVNames(x, varNames)
	ans <- matrix(0.1, nrow = nrows, ncol = length(flatNames))
	colIndex = 1
	for(i in seq_along(varNames)){
		.Call(fastMatrixInsert, ans, modelValuesElement2Matrix(x, varNames[i]) , as.integer(1), as.integer(colIndex) ) 		
		colIndex = colIndex + prod(x$sizes[[varNames[i] ]])
		}
	colnames(ans) <- flatNames
	return(ans)
}

as.list.modelValuesBaseClass <- function(x, varNames, iterationAsLastIndex = FALSE, ...) {
  if(missing(varNames))
    varNames <- x$varNames
  nrows <- getsize(x)
  results <- list()
  for(v in varNames) {
    samples <- x[[v]]
    dims <- dimOrLength(samples[[1]])
    matrixVersion <- do.call("c", lapply(samples, as.numeric))
    ansDims <- c(dims, nrows)
    results[[v]] <- array(matrixVersion, dim = ansDims)
    if(!iterationAsLastIndex) {
      nDim <- length(ansDims)
      results[[v]] <- aperm(results[[v]], c(nDim, 1:(nDim-1)))
    }
  }
  results
}

as.list.CmodelValues <- function(x, varNames, iterationAsLastIndex = FALSE, ...) {
  if(missing(varNames))
    varNames <- x$varNames
  nrows <- getsize(x)
  results <- list()
  for(v in varNames) {
    samples <- x[[v]]
    dims <- dimOrLength(samples[[1]])
    matrixVersion <- do.call("c", lapply(samples, as.numeric))
    ansDims <- c(dims, nrows)
    results[[v]] <- array(matrixVersion, dim = ansDims)
    if(!iterationAsLastIndex) {
      nDim <- length(ansDims)
      results[[v]] <- aperm(results[[v]], c(nDim, 1:(nDim-1)))
    }
  }
  results
}

modelValuesElement2Matrix <- function(mv, varName){
	if(length(varName) != 1)
		stop('modelValuesElement2Matrix is a call for a single variable. For multiple variables, use modelValues2Matrix')
	matrix( as.numeric(unlist(mv[[varName]]) ), ncol = prod(mv$sizes[[varName]]), byrow = TRUE )
}


Cmatrix2mvOneVar <- function(mat, mv, varName, k){
	ptr <- mv$componentExtptrs[[varName]]
	if(inherits(ptr, 'externalptr')) {
            eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$matrix2VecNimArr, ptr, mat, rowStart = as.integer(1), rowEnd = k ))
        } else
            stop('varName not found in modelValues')
}

#' Set values of one variable of a modelValues object from an R matrix
#'
#' Normally a modelValues object is accessed one "row" at a time.  This function allows all rows for one variable to set from a matrix with one dimension more than the variable to be set.
#'
#' @param mat Input matrix
#'
#' @param mv  modelValues object to be modified.
#'
#' @param varName Character string giving the name of the variable on \code{mv} to be set
#'
#' @param k Number of rows to use
#'
#' @export
#' @details
#' This function may be deprecated in the future when a more natural system for interacting with modelValues objects is developed.
Rmatrix2mvOneVar <- function(mat, mv, varName, k){
	if( mv$symTab$symbols[[varName]][['type']] == 'double'){
		storage.mode(mat) <- 'double'
		len <- ncol(mat)
		.Call(matrix2ListDouble, mat, mv[[varName]], listStartIndex = as.integer(1), RnRows = k, Rlength = as.integer(mv$sizes[[varName]]) )
	}
	if( mv$symTab$symbols[[varName]][['type']] == 'int'){
		storage.mode(mat) <- 'integer'
		len <- ncol(mat)
		.Call(matrix2ListInt, mat, mv[[varName]], listStartIndex = as.integer(1), RnRows = k, Rlength = as.integer(mv$sizes[[varName]]) )
	}
}

matrix2mv <- function(mat, mv){
	k <- nrow(mat)
	if(mv$getSize() < k)
		mv$resize(k)
	mvVarNames <- mv$varNames
	mvSizes <- mv$sizes
	colNames <- colnames(mat)
	colNames <- removeIndexing(colNames)
	uniqueColNames <- unique(colNames)
	if(inherits(mv, 'CmodelValues')	){
		for(vN in uniqueColNames){
			totVals <- prod(mvSizes[[vN]])
			varInds <- colNames == vN
			if(totVals != sum(varInds) )
				stop('matrix2mv halted because dimensions of variables do not match')
			subMatrix <- as.matrix(mat[, varInds])
			Cmatrix2mvOneVar(subMatrix, mv, vN, k)
		}
	}
	else if(inherits(mv, 'modelValuesBaseClass') ){
		for(vN in uniqueColNames){
			totVals <- prod(mvSizes[[vN]])
			varInds <- colNames == vN
			if(totVals != sum(varInds) )
				stop('matrix2mv halted because dimensions of variables do not match')
			subMatrix <- as.matrix(mat[, varInds])
			Rmatrix2mvOneVar(subMatrix, mv, vN, k)
		}
	}
	else
		stop('argument mv is neither a CmodelValues or RmodelValues object')
}

