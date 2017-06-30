makeSingleArgWrapper <- function(nf, wrt, fxnEnv) {
  formalNames <- formalArgs(eval(nf[[1]], envir = fxnEnv)@.Data)
  wrtNames <- strsplit(wrt, '\\[')
  wrtArgIndices <- c(sapply(wrtNames, function(x){return(which(x[1] == formalNames))}))
  flatteningInfo <- list()
  for(i in seq_along(wrt)){
    arg <- nf[[wrtArgIndices[i] + 1]]
    indText <- ''
    if(length(wrtNames[[i]]) > 1){
      indText <-  paste0('[', wrtNames[[i]][[2]])
    }
    argSym <- parse(text = paste0(deparse(arg), indText))[[1]]
    dimInfo <- dim(eval(argSym, envir = fxnEnv))
    if(is.null(dimInfo)) dimInfo <- length(eval(argSym, envir = fxnEnv))
    flatteningInfo[[i]] <- list()
    flatteningInfo[[i]][[1]] <- dimInfo
    flatteningInfo[[i]][[2]] <- indText
  }
  wrappedFun <- function(x) {
    args <- list()
    nextInd <- 0
    for(i in 1:(length(nf)-1)){
      if(i %in% wrtArgIndices){
        thisWrtIndex <- which(i == wrtArgIndices)
        thisSize <- prod(flatteningInfo[[thisWrtIndex]][[1]]) ## total length
        if(length(wrtNames[[thisWrtIndex]]) > 1){
          args[[i]] <-  eval(nf[[i+1]], envir = fxnEnv)
          argBracketExpr <- parse(text = paste0('args[[i]]', flatteningInfo[[thisWrtIndex]][[2]], ' <- x[nextInd + 1:thisSize]'))[[1]]
          eval(argBracketExpr)
        }
        else{
          args[[i]] <-  x[nextInd + 1:thisSize]
        }
        if(length( flatteningInfo[[thisWrtIndex]][[1]]) > 1 ) dim(  args[[i]]) <-  flatteningInfo[[thisWrtIndex]][[1]]
        nextInd <- nextInd + thisSize
      }
      else{
       args[[i]] <- nf[[i+1]] 
      }
    }
    c(do.call(paste(nf[[1]]), args, envir = fxnEnv))
  }
  makeSingleArg <- function(){
    singleArg <- c()
    for(i in seq_along(wrtNames)){
      if(length(wrtNames[[i]]) > 1){
        arg <- eval(nf[[wrtArgIndices[i]+1]], envir = fxnEnv)
        singleArg <- c(singleArg, c(eval(parse(text =  paste0('arg', flatteningInfo[[i]][[2]]))[1])))
      }
      else{
        singleArg <- c(singleArg, c(eval(nf[[wrtArgIndices[i]+1]], envir = fxnEnv)))
      }
    }
    return(singleArg)
  }
  return(list(wrappedFun, makeSingleArg))
}


#' Nimble Derivatives
#' 
#' EXPERIMENTAL Computes the value, gradient, and Hessian of a given  \code{nimbleFunction} method.  
#' 
#' @param nimFxn a call to a \code{nimbleFunction} method with arguments included.
#' @param order an integer vector with values within the set {0, 1, 2}, corresponding to whether the function value, gradient, and Hessian should be returned respectively.
#' 
#' @export
nimDerivs <- function(nimFxn = NA, order = nimC(0,1,2), dropArgs = NULL, wrt = NULL){
  fxnEnv <- parent.frame()
  fxnCall <- match.call(function(nimFxn, order, dropArgs, wrt){})
  if(is.null(fxnCall[['order']])) fxnCall[['order']] <- order
  derivFxnCall <- fxnCall[['nimFxn']]
  if(deparse(derivFxnCall[[1]]) == 'calculate'){
    derivFxnCall <- match.call(calculate, derivFxnCall)
    return(nimDerivs_calculate(model = eval(derivFxnCall[['model']], envir = fxnEnv),
                               nodes = eval(derivFxnCall[['nodes']], envir = fxnEnv),
                               nodeFxnVector = eval(derivFxnCall[['nodeFxnVector']], envir = fxnEnv),
                               nodeFunctionIndex = eval(derivFxnCall[['nodeFunctionIndex']], envir = fxnEnv),
                               order, wrt))
  }
  if(deparse(derivFxnCall[[1]]) == 'model$calculate'){
    return(nimDerivs_calculate(model = eval(derivFxnCall[[1]][[2]], envir = fxnEnv),
                               nodes = eval(derivFxnCall[[2]], envir = fxnEnv),
                               order = order, wrt = wrt))
  }
  libError <- try(library('numDeriv'), silent = TRUE)
  if(inherits(libError, 'try-error')){
    stop("The 'numDeriv' package must be installed to use derivatives in uncompiled nimbleFunctions.")
  }
  if(is.null(wrt)){
    wrt <- formalArgs(eval(derivFxnCall[[1]], envir = fxnEnv)@.Data)
    wrt <- wrt[- which(wrt == dropArgs)]
  }
  derivFxnList <- makeSingleArgWrapper(derivFxnCall, wrt, fxnEnv)
  singleArg <- derivFxnList[[2]]()
  derivList <- genD(derivFxnList[[1]], singleArg)
  outGrad <- derivList$D[,1:derivList$p, drop = FALSE]
  outHessVals <- derivList$D[,(derivList$p + 1):dim(derivList$D)[2], drop = FALSE]
  outHess <- array(NA, dim = c(derivList$p, derivList$p, length(derivList$f0)))
  singleDimMat <- matrix(NA, nrow = derivList$p, ncol = derivList$p)
  singleDimMatUpperTriDiag <- upper.tri(singleDimMat, diag = TRUE)
  for(outDim in seq_along(derivList$f0)){
    singleDimMat[singleDimMatUpperTriDiag] <- outHessVals[outDim,]
    outHess[,,outDim] <- singleDimMat
  }
  outHess[lower.tri(outHess[,,1])] <-   outHess[upper.tri(outHess[,,1])]
  return(nimble:::ADNimbleList$new(value = derivList$f0,
                                   gradient = outGrad,
                                   hessian = outHess))
}

rDeriv_CalcNodes <- function(model, nfv, derivInfo, calcNodesLineNums, wrtLineInfo){
  model <- nfv$model
  indexingInfo <- nfv$indexingInfo
  declIDs <- indexingInfo$declIDs
  numNodes <- length(declIDs)
  unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
  chainRuleDerivList <- list()
  chainRuleHessianList <- list()
  ## totalWRTSize is the sum of the lengths of all wrt parameters.
  totalWRTSize <- sum(sapply(wrtLineInfo, function(x){return(x$lineSize)}))
  ## outDerivList will be returned from this function.
  outDerivList <- nimble:::ADNimbleList$new(value = 0,
                                            gradient = matrix(0, 
                                                              ncol = totalWRTSize,
                                                              nrow = 1),
                                            hessian = array(0, dim = c(totalWRTSize,
                                                                       totalWRTSize,
                                                                       1)))
  ## Below we get the start and end indices of all wrt parameters.
  ## For example, if we are taking derivatives with respect to 'a' and 'b', both of 
  ## length 2, outIndexStartPoints would be c(1, 3) and outIndexEndPoints would be c(2,4)
  outIndexEndPoints <- cumsum(sapply(wrtLineInfo, function(x){return(x$lineSize)}))
  outIndexStartPoints <- c(1, (outIndexEndPoints[-length(outIndexEndPoints)]+1))
  for(i in seq_along(wrtLineInfo)){
    wrtLineInfo[[i]]$lineIndices <- outIndexStartPoints[i]:outIndexEndPoints[i]
  }
  for(i in 1:numNodes) {
    if(length(calcNodesLineNums) > 0){
      declID <- declIDs[i]
      isDeterm <- model$modelDef$declInfo[[declID]]$type == 'determ'
      thisWrtLine <- which(sapply(wrtLineInfo, function(x){return(x$lineNum == i)}))
      isWrtLine <- length(thisWrtLine) > 0
      isCalcNodeLine <- i %in% calcNodesLineNums
      sizeAndDimInfo <- environment(model$nodeFunctions[[declID]]$.generatorFunction)[['parentsSizeAndDims']]
      formalArgNames <- formals(model$nodeFunctions[[ declID ]]$calculateWithArgs)
      unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ unrolledIndicesMatrixRows[i], ]
      calcWithArgs <- model$nodeFunctions[[ declID ]]$calculateWithArgs
      calcWithArgsCall <- as.call(c(list(as.name('calcWithArgs'), unrolledIndicesMatrixRow), lapply(names(formalArgNames)[-1],
                                                                                                    function(x){parse(text = convertCalcArgNameToModelNodeName(x, sizeAndDimInfo))[[1]]})))
      thisNodeSize <- length(eval(calcWithArgsCall[[3]]))   
      ## If this node is either a calulate node or a deterministic dependency of a wrt node,
      ## we need to take derivatives of its calculateWithArgs function.  The derivative function
      ## call is evaluated below.
      if(isCalcNodeLine || isDeterm){
        derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS),
                                     list(CALCCALL = calcWithArgsCall,
                                          DERIVORDERS = c(0, 1, 2),
                                          DROPARGS = 'INDEXEDNODEINFO_')))
        argSizeInfo <- sapply(3:length(calcWithArgsCall), function(x){length(eval(calcWithArgsCall[[x]]))})
        if(isDeterm){
          derivList$value <- 0
          model$nodeFunctions[[declID]]$calculate(unrolledIndicesMatrixRow)
        }
        
        ## The derivOutputFlag determines whether the derivatives of this node (node i): 
        ## should be calculated for inclusion in the chain rule output (TRUE),
        ## should be calculated for later use in the chain rule (FALSE), 
        derivOutputFlag <- if(isDeterm) FALSE else TRUE
        
        ## Below we construct two lists:
        ## parentGradients, a list of all the gradients of the parent nodes of each argument of this node.
        ## parentHessians, a list of all the hessians of the parent nodes of this node.  
        parentGradients <- vector('list', length = length(derivInfo[[2]][[i]]))
        parentHessians <- vector('list', length = length(derivInfo[[2]][[i]]))
        for(k in seq_along(derivInfo[[2]][[i]])){
          if(k == 1 && isWrtLine){
            ## The first argument (k = 1) of a node's calculateWithArgs function will always be the node itself.
            ## If this node is a wrt node, we want to set the parent gradient of the first arg (the derivative of this node wrt itself)
            ## to the identity.
            parentGradients[[k]] <- matrix(0, nrow = length(wrtLineInfo[[thisWrtLine]]$lineIndices),
                                           ncol = totalWRTSize)
            parentGradients[[k]][,wrtLineInfo[[thisWrtLine]]$lineIndices] <- diag(length(wrtLineInfo[[thisWrtLine]]$lineIndices))
          }
          else if(derivInfo[[2]][[i]][[k]][1] > 0){
            ## Otherwise, if argument k has parents that depend on a wrt node, grab the parent gradients (which will have already been calculated)
            ## and combine them into a single matrix.
            parentGradientsList <- chainRuleDerivList[derivInfo[[2]][[i]][[k]]]
            parentGradients[[k]] <- parentGradientsList[[1]]
            if(length(parentGradientsList) > 1){
              for(i_listEntry in 2:length(parentGradientsList)){
                parentGradients[[k]] <- rbind(parentGradients[[k]], parentGradientsList[[i_listEntry]])
              }
            }
            
            ## Similarly, grab the parent Hessians (which will have already been calculated)
            ## and combine them into a single array.
            parentHessiansList <- chainRuleHessianList[derivInfo[[2]][[i]][[k]]]
            parentHessiansDim <- sum(sapply(parentHessiansList, function(x){return(dim(x)[3])}))
            parentHessians[[k]] <- array(NA, dim = c(dim(parentHessiansList[[1]])[1], dim(parentHessiansList[[1]])[1],
                                                     parentHessiansDim))
            thisDim <- 1
            for(i_listEntry in 1:length(parentHessiansList)){
              parentHessians[[k]][, , thisDim:(thisDim + dim(parentHessiansList[[i_listEntry]])[3] - 1)] <- parentHessiansList[[i_listEntry]]
              thisDim <- thisDim +  dim(parentHessiansList[[i_listEntry]])[3]
            }
          }
        }
        
        if(derivOutputFlag == TRUE){
          ## If these derivatives will be included in output, we are taking derivative of a log prob. calculation,
          ## so ouput will be length 1, and input args will be all wrt args.
          chainRuleDerivList[[i]] <- matrix(0, 
                                            nrow = 1, 
                                            ncol = totalWRTSize)
          chainRuleHessianList[[i]] <- array(0, dim = c(totalWRTSize,
                                                        totalWRTSize,
                                                        1))
        }
        else{
          ## otherwise, we are taking derivative of a node value calculation (not log prob. calculation).
          ## so ouput will be the length of this node, and input args will be all wrt args.
          chainRuleDerivList[[i]] <- matrix(0,
                                            nrow = thisNodeSize,  
                                            ncol = totalWRTSize)
          chainRuleHessianList[[i]] <- array(0, dim = c(totalWRTSize,
                                                        totalWRTSize,
                                                        thisNodeSize))
        }
        
        ## Iterate over all wrt params.
        for(j in seq_along(wrtLineInfo)){
          thisArgIndex <- 0
          ## Iterate over this line's parent nodes.
          for(k in seq_along(derivInfo[[2]][[i]])){
              if(!is.null(parentGradients[[k]])){
                ## Calculate derivs of this node (i) wrt this parameter (j) for this parent node (k) via chain rule. 
                chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices] <- chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices] +
                  derivList$gradient[,(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]), drop = FALSE]%*%parentGradients[[k]][,wrtLineInfo[[j]]$lineIndices, drop = FALSE]
              }
            thisArgIndex <- thisArgIndex + argSizeInfo[k]
          }
          if(derivOutputFlag == TRUE){
            ## If this line is included in output, add the derivative of this line (i) wrt this param (j).
            outDerivList$gradient[, wrtLineInfo[[j]]$lineIndices] <- outDerivList$gradient[, wrtLineInfo[[j]]$lineIndices]  +  
              chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices]
          }
          
          ## The Hessian is calculated below using Faà di Bruno's formula.
          ## Second iteration over wrt parameters 
          for(j_2 in j:length(wrtLineInfo)){
            thisArgIndex <- 0
            ## Iterate over this line's parent nodes.
            for(k in seq_along(derivInfo[[2]][[i]])){
              if(!is.null(parentHessians[[k]])){
                for(dim1 in wrtLineInfo[[j]]$lineIndices){
                  for(dim2 in wrtLineInfo[[j_2]]$lineIndices){
                    chainRuleHessianList[[i]][dim1, dim2, ] <- chainRuleHessianList[[i]][dim1, dim2, ] +
                      c(derivList$gradient[ ,(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]), drop = FALSE]%*%parentHessians[[k]][dim1, dim2, , drop = FALSE])
                  }
                }
              }
              thisArgIndex_2 <- 0
              for(k_2 in seq_along(derivInfo[[2]][[i]])){
                    if(!is.null(parentGradients[[k]]) && !is.null(parentGradients[[k_2]])){
                      for(dim3 in 1:dim(derivList$hessian)[3]){
                        chainRuleHessianList[[i]][wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j_2]]$lineIndices, dim3] <- chainRuleHessianList[[i]][wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j_2]]$lineIndices, dim3] +
                          t(parentGradients[[k]][, wrtLineInfo[[j]]$lineIndices, drop = FALSE])%*%derivList$hessian[(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]),(thisArgIndex_2 + 1):(thisArgIndex_2 + argSizeInfo[k_2]), dim3]%*%
                          parentGradients[[k_2]][, wrtLineInfo[[j_2]]$lineIndices, drop = FALSE]
                      }
                    }
                thisArgIndex_2 <- thisArgIndex_2 + argSizeInfo[k_2]
              }
              thisArgIndex <- thisArgIndex + argSizeInfo[k]
            }
            if(derivOutputFlag == TRUE){
              ## If this line is included in output, add the Hessian of this line (i) wrt this param #1 (j) and this param #2 (j_2).
              outDerivList$hessian[wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j_2]]$lineIndices, ] <-  outDerivList$hessian[wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j_2]]$lineIndices, ]   +
                chainRuleHessianList[[i]][wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j_2]]$lineIndices,]
            }
          }
        }
      }
      if(isWrtLine){
        ## If this is a wrt node, we need to set the chainRule lists appropriately so that the chain rule
        ## will work for dependent nodes of this node.  That means taking the first and second derivs of the
        ## function f(x) = x, which will be the identity matrix and 0 respectively.
        chainRuleDerivList[[i]] <- matrix(0,
                                          nrow = thisNodeSize,  
                                          ncol = totalWRTSize)
        chainRuleHessianList[[i]] <- array(0, dim = c(totalWRTSize,
                                                      totalWRTSize,
                                                      thisNodeSize))
        chainRuleDerivList[[i]][,wrtLineInfo[[thisWrtLine]]$lineIndices] <- diag(thisNodeSize)
      }
      if(isCalcNodeLine){
        outDerivList$value <- outDerivList$value + derivList$value
      }
    }
  }
  wrtIndex <- 1
  ## Finally, we check to see which parts of the wrt parameters are returned.
  ## E.g. if x[1:2] is a node in the model, but derivs were taken wrt x[1], we drop the gradient
  ## and hessian information corresponding to x[2]
  for(i in seq_along(wrtLineInfo)){
    if(!is.na(wrtLineInfo[[i]]$correctIndices)){
      theseIndices <- (wrtIndex:(wrtIndex + wrtLineInfo[[i]]$lineSize -1))[wrtLineInfo[[i]]$correctIndices]
      allIndices <- c(1:(wrtIndex - 1), theseIndices)
      if(!i == length(wrtLineInfo)){
        allIndices <- c(allIndices, (wrtIndex + wrtLineInfo[[i]]$lineSize):length(outDerivList$gradient[1,]))
      }
      allIndices <- unique(allIndices)
      outDerivList$gradient <- outDerivList$gradient[,allIndices, drop = FALSE]
      outDerivList$hessian <- outDerivList$hessian[allIndices, allIndices, , drop = FALSE]
      wrtIndex <- wrtIndex + length(theseIndices)
    }
    else{
      wrtIndex <- wrtIndex + wrtLineInfo[[i]]$lineSize
    }
  }
  ## Reflect hessian across the diagonal
  outDerivList$hessian[lower.tri(outDerivList$hessian[,,1])] <-  outDerivList$hessian[upper.tri(outDerivList$hessian[,,1])]
  return(outDerivList)
}

convertCalcArgNameToModelNodeName <- function(calcArgName, sizeAndDimInfo){
  thisModelElementNum <- as.numeric(gsub(".*([0-9]+)$", "\\1", calcArgName)) ## Extract 1, 2, etc. from end of arg name.
  thisName <- sub("_[0-9]+$","",calcArgName) ## Extract node name from beginning of arg name.
  indexBracketInfo <- paste0('[',
                             paste0(sapply(sizeAndDimInfo[[thisName]][[thisModelElementNum]]$indexExpr, function(x){
                               if(length(x) == 1) return(deparse(x[[1]]))
                               else if( deparse(x[[1]]) == 'getNodeFunctionIndexedInfo'){
                                 x[[2]] <- parse(text = 'unrolledIndicesMatrixRow')[[1]]
                               } 
                               else{
                                 return(paste0(deparse(x[[1]]), ':', deparse(x[[2]])))
                               }}), collapse = ', '),
                             ']')
  return(paste0('model$', thisName, indexBracketInfo))
}

nimDerivs_calculate <- function(model, nodes = NA, nodeFxnVector = NULL, nodeFunctionIndex = NULL, order, wrtPars){
  if(!is.null(nodeFxnVector)){
    stop('nfv case of nimDerivs_calculate not implemented yet.')
  }
  wrtParsDeps <- model$getDependencies(wrtPars)
  nodes <- model$expandNodeNames(nodes)
  if(!all(nodes %in% wrtParsDeps)){
    warning('not all calculate nodes depend on a wrtNode')
    wrtParsDeps <- model$topologicallySortNodes(c(wrtParsDeps, nodes[which(!(nodes%in%wrtParsDeps))]))
  }
  derivInfo <- nimble:::enhanceDepsForDerivs(model$expandNodeNames(wrtPars), wrtParsDeps, model)
  stochNodes <- nodes[model$getNodeType(nodes) == 'stoch']
  calcNodesLineNums <- sapply(stochNodes, function(x){which(x == derivInfo[[1]])})
  wrtLineInfo <- list()
  for(i in seq_along(model$expandNodeNames(wrtPars))){
    wrtLineInfo[[i]] <- list()
    wrtLineInfo[[i]]$lineNum <- which(model$expandNodeNames(wrtPars)[i] == derivInfo[[1]])
    wrtLineInfo[[i]]$lineSize <- length(model[[model$expandNodeNames(wrtPars[i])]])
    if(!identical(model$expandNodeNames(wrtPars[i], returnScalarComponents = TRUE),
                  model$expandNodeNames(model$expandNodeNames(wrtPars[i]), returnScalarComponents = TRUE))){
      expandWrt <- model$expandNodeNames(wrtPars[i], returnScalarComponents = TRUE)
      correctIndices <- numeric(length(expandWrt))
      for(j in seq_along(expandWrt)){
        correctIndices[j] <- which(model$expandNodeNames(expandWrt[j], returnScalarComponents = TRUE) == 
                                     model$expandNodeNames(model$expandNodeNames(wrtPars[i]), returnScalarComponents = TRUE))
      }
    }
    else{
      correctIndices <- NA
    }
    wrtLineInfo[[i]]$correctIndices <- correctIndices
  }
  if(inherits(model, 'modelBaseClass') ){
    if(missing(nodes) ) 
      nodes <- model$getMaps('nodeNamesLHSall')
    nfv <- nodeFunctionVector(model, wrtParsDeps, sortUnique = FALSE)
    return(rDeriv_CalcNodes(model, nfv, derivInfo, calcNodesLineNums, wrtLineInfo))
  }	
}

