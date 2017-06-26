
#' Nimble Derivatives
#' 
#' EXPERIMENTAL Computes the value, gradient, and Hessian of a given  \code{nimbleFunction} method.  The R version is currently unimplemented.
#' 
#' @param nimFxn a call to a \code{nimbleFunction} method with arguments included.
#' @param order an integer vector with values within the set {0, 1, 2}, corresponding to whether the function value, gradient, and Hessian should be returned respectively.
#' 
#' @export
nimDerivs <- function(nimFxn = NA, order = nimC(0,1,2), dropArgs = NULL, wrtPars = NULL){
  fxnEnv <- parent.frame()
  fxnCall <- match.call(function(nimFxn, order, dropArgs, wrtPars){})
  if(is.null(fxnCall[['order']])) fxnCall[['order']] <- order
  derivFxnCall <- fxnCall[['nimFxn']]
  if(deparse(derivFxnCall[[1]]) == 'calculate'){
    derivFxnCall <- match.call(calculate, derivFxnCall)
    return(nimDerivs_calculate(model = eval(derivFxnCall[['model']], envir = fxnEnv),
                               nodes = eval(derivFxnCall[['nodes']], envir = fxnEnv),
                               nodeFxnVector = eval(derivFxnCall[['nodeFxnVector']], envir = fxnEnv),
                               nodeFunctionIndex = eval(derivFxnCall[['nodeFunctionIndex']], envir = fxnEnv),
                               order, wrtPars))
  }
  useArgs <- 1:(length(derivFxnCall) - 1)
  if(!is.null(dropArgs)) useArgs <- useArgs[-dropArgs]
  fxnArgLengths <- sapply(derivFxnCall, function(x){return(length(eval(x, envir = fxnEnv)))})[-1]
  totalFxnArgLength <- sum(fxnArgLengths[useArgs])
  delta <- .0001
  origValue <- eval(derivFxnCall, envir = fxnEnv)
  outLength <- sum(nimDim(origValue))
  fxph <- matrix(nrow = outLength, ncol = totalFxnArgLength)
  fxmh <- matrix(nrow = outLength, ncol = totalFxnArgLength)
  grad <- matrix(nrow = outLength, ncol = totalFxnArgLength)
  derivxy <- array(0, dim = c(totalFxnArgLength, totalFxnArgLength, outLength))
  #if(length(origValue) > 1) stop('Currently only have R derivs for functions that return scalars.')
  for(i in useArgs){
    origDFxnCall <-  derivFxnCall[[i + 1]] 
    deltaVec <- rep(0, fxnArgLengths[i])
    for(j in 1:fxnArgLengths[i]){
      thisArgNum <- if(i > 1) sum(fxnArgLengths[intersect(useArgs, c(1:(i-1)))]) + j else j
      deltaVec[j] <- delta
      derivFxnCall[[i + 1]] <- substitute(DFXNCALL + DELTAVEC, 
                                          list(DFXNCALL = derivFxnCall[[i + 1]],
                                               DELTAVEC = deltaVec))
      fxph[,thisArgNum] <- c(eval(derivFxnCall, envir = fxnEnv))
      derivFxnCall[[i + 1]] <- substitute(DFXNCALL - 2*DELTAVEC, 
                                          list(DFXNCALL = derivFxnCall[[i + 1]],
                                               DELTAVEC = deltaVec))    
      fxmh[,thisArgNum] <- c(eval(derivFxnCall, envir = fxnEnv))
      grad[,thisArgNum] <- (fxph[,thisArgNum] - fxmh[,thisArgNum])/(2*delta)
      derivxy[thisArgNum, thisArgNum, ] <- (fxph[,thisArgNum] -2*origValue + fxmh[,thisArgNum])/(delta^2)
      derivFxnCall[[i + 1]] <- origDFxnCall
      deltaVec[j] <- 0
    }
  }
  for(i in useArgs){
    origDFxnCall <-  derivFxnCall[[i + 1]] 
    deltaVec <- rep(0, fxnArgLengths[i])
    for(j in 1:fxnArgLengths[i]){
      deltaVec[j] <- delta
      thisArgNum <- if(i > 1) sum(fxnArgLengths[intersect(useArgs, c(1:(i-1)))]) + j else j
      if(j != fxnArgLengths[i]){
        for(j_2 in (j+1):fxnArgLengths[i]){
          thisArgNum_2 <- if(i > 1) sum(fxnArgLengths[intersect(useArgs, c(1:(i-1)))]) + j_2 else j_2
          deltaVec[j_2] <- delta
          derivFxnCall[[i + 1]] <- substitute(DFXNCALL + DELTAVEC, 
                                              list(DFXNCALL = derivFxnCall[[i + 1]],
                                                   DELTAVEC = deltaVec))
          fxyph <-  c(eval(derivFxnCall, envir = fxnEnv))
          derivFxnCall[[i + 1]] <- substitute(DFXNCALL - 2*DELTAVEC, 
                                              list(DFXNCALL = derivFxnCall[[i + 1]],
                                                   DELTAVEC = deltaVec))
          fxymh <-  c(eval(derivFxnCall, envir = fxnEnv))
          derivxy[thisArgNum, thisArgNum_2, ] <- 
            (fxyph - fxph[,thisArgNum] - fxph[,thisArgNum_2] + 2*origValue - fxmh[,thisArgNum] - fxmh[,thisArgNum_2] + fxymh)/(2*delta^2)
          derivFxnCall[[i + 1]] <- origDFxnCall
          deltaVec[j_2] <- 0
        }
      }
      if(i != useArgs[length(useArgs)]){
        for(i_2 in useArgs[-c(1:which(useArgs == i))]){
          origDFxnCall_2 <-  derivFxnCall[[i_2 + 1]] 
          deltaVec_2 <- rep(0, fxnArgLengths[i_2])
          for(j_2 in 1:fxnArgLengths[i_2]){
            thisArgNum_2 <- sum(fxnArgLengths[intersect(useArgs, c(1:(i_2-1)))]) + j_2
            deltaVec_2[j_2] <- delta
            derivFxnCall[[i + 1]] <- substitute(DFXNCALL + DELTAVEC, 
                                                list(DFXNCALL = derivFxnCall[[i + 1]],
                                                     DELTAVEC = deltaVec))
            derivFxnCall[[i_2 + 1]] <- substitute(DFXNCALL + DELTAVEC, 
                                                  list(DFXNCALL = derivFxnCall[[i_2 + 1]],
                                                       DELTAVEC = deltaVec_2))
            fxyph <-  c(eval(derivFxnCall, envir = fxnEnv))
            derivFxnCall[[i + 1]] <- substitute(DFXNCALL - 2*DELTAVEC, 
                                                list(DFXNCALL = derivFxnCall[[i + 1]],
                                                     DELTAVEC = deltaVec))
            derivFxnCall[[i_2 + 1]] <- substitute(DFXNCALL - 2*DELTAVEC, 
                                                  list(DFXNCALL = derivFxnCall[[i_2 + 1]],
                                                       DELTAVEC = deltaVec_2))
            fxymh <-  c(eval(derivFxnCall, envir = fxnEnv))
            derivxy[thisArgNum, thisArgNum_2, ] <- 
              (fxyph - fxph[,thisArgNum] - fxph[,thisArgNum_2] + 2*origValue - fxmh[,thisArgNum] - fxmh[,thisArgNum_2] + fxymh)/(2*delta^2)
            derivFxnCall[[i + 1]] <- origDFxnCall
            derivFxnCall[[i_2 + 1]] <- origDFxnCall_2
            deltaVec_2[j_2] <- 0
          }
        }
      }
      deltaVec[j] <- 0
    }
  }
  derivxy[lower.tri(derivxy[,,1])] <-   derivxy[upper.tri(derivxy[,,1])]
  return(nimble:::ADNimbleList$new(value = origValue,
                                   gradient = grad,
                                   hessian = derivxy))
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
                                          DROPARGS = 1)))
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

nimDerivs_calculate <- function(model, nodes, nodeFxnVector, nodeFunctionIndex, order, wrtPars){
  if(!is.null(nodeFxnVector)){
    stop('nfv case of nimDerivs_calculate not implemented yet.')
  }
  wrtParsDeps <- model$getDependencies(wrtPars)
  nodes <- model$expandNodeNames(nodes)
  if(!all(nodes %in% wrtParsDeps)){
    print('Warning: not all calculate nodes depend on a wrtNode')
  }
  derivInfo <- nimble:::enhanceDepsForDerivs(model$expandNodeNames(wrtPars), wrtParsDeps, model)
  stochNodes <- nodes[model$getNodeType(nodes) == 'stoch']
  calcNodesLineNums <- sapply(stochNodes, function(x){which(x == derivInfo[[1]])})
  wrtLineInfo <- list()
  for(i in seq_along(model$expandNodeNames(wrtPars))){
    wrtLineInfo[[i]] <- list()
    wrtLineInfo[[i]]$lineNum <- which(model$expandNodeNames(wrtPars)[i] == derivInfo[[1]])
    wrtLineInfo[[i]]$lineSize <- length(values(model, model$expandNodeNames(wrtPars)[i]))
  }
  if(inherits(model, 'modelBaseClass') ){
    if(missing(nodes) ) 
      nodes <- model$getMaps('nodeNamesLHSall')
    nfv <- nodeFunctionVector(model, wrtParsDeps)
    return(rDeriv_CalcNodes(model, nfv, derivInfo, calcNodesLineNums, wrtLineInfo))
  }	
}

