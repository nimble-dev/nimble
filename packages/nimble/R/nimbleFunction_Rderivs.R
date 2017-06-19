
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
      derivxy[thisArgNum, thisArgNum,] <- (fxph[,thisArgNum] -2*origValue + fxmh[,thisArgNum])/(delta^2)
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
          derivxy[thisArgNum, thisArgNum_2,] <- 
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
            derivxy[thisArgNum, thisArgNum_2,] <- 
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
  return(nimble:::ADNimbleList$new(value = origValue,
                                   gradient = grad,
                                   hessian = derivxy))
}


rDeriv_CalcNodes <- function(model, nfv, derivInfo, nodesLineNums, wrtLineInfo){
  model <- nfv$model
  # useCompiledNonNestedInterface <- inherits(model, 'CmodelBaseClass') & !getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')
  indexingInfo <- nfv$indexingInfo
  declIDs <- indexingInfo$declIDs
  numNodes <- length(declIDs)
  # if(numNodes < 1) return(l_Prob)
  unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
  chainRuleDerivList <- list()
  chainRuleHessianList <- list()
  outValue <- 0
  outGradient <- matrix(0, ncol = sum(sapply(wrtLineInfo, function(x){return(x$lineSize)})),
                        nrow = 1)
  outHessian <- matrix(0, nrow = sum(sapply(wrtLineInfo, function(x){return(x$lineSize)})),
                       ncol = sum(sapply(wrtLineInfo, function(x){return(x$lineSize)})))
  outIndexEndPoints <- cumsum(sapply(wrtLineInfo, function(x){return(x$lineSize)}))
  outIndexStartPoints <- c(1, (outIndexEndPoints[-length(outIndexEndPoints)]+1))
  for(i in seq_along(wrtLineInfo)){
    wrtLineInfo[[i]]$lineIndices <- outIndexStartPoints[i]:outIndexEndPoints[i]
  }
  for(i in 1:numNodes) {
    declID <- declIDs[i]
    sizeAndDimInfo <- environment(model$nodeFunctions[[declID]]$.generatorFunction)[['parentsSizeAndDims']]
    formalArgNames <- formals(model$nodeFunctions[[ declID ]]$calculateWithArgs)
    unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ unrolledIndicesMatrixRows[i], ]
    calcWithArgs <- model$nodeFunctions[[ declID ]]$calculateWithArgs
    calcWithArgsCall <- as.call(c(list(as.name('calcWithArgs'), unrolledIndicesMatrixRow), lapply(names(formalArgNames)[-1],
                                                                                                  function(x){parse(text = convertCalcArgNameToModelNodeName(x, sizeAndDimInfo))[[1]]})))
    derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS),
                                 list(CALCCALL = calcWithArgsCall,
                                      DERIVORDERS = c(0, 1, 2),
                                      DROPARGS = 1)))
    argSizeInfo <- sapply(3:length(calcWithArgsCall), function(x){length(eval(calcWithArgsCall[[x]]))})
    names(argSizeInfo) <- names(formals(calcWithArgs))[2:length(formals(calcWithArgs))]  ## used only for debugging
    if(model$modelDef$declInfo[[declID]]$type == 'determ'){
      derivList$value <- 0
      model$nodeFunctions[[declID]]$calculate(unrolledIndicesMatrixRow)
    }
    isWrtLine <- any(sapply(wrtLineInfo, function(x){return(x$lineNum == i)}))
    isNodeLine <- i %in% nodesLineNums
    
    ## The derivCalcFlags vector constructed below determines whether this node (i) should be treated as an output node (FALSE),
    ## a node that is necessary for later use in the chain rule (TRUE), or both (FALSE, TRUE).
    if(model$modelDef$declInfo[[declID]]$type == 'determ'){
      derivCalcFlags <- TRUE
    }
    else if(isNodeLine){
      derivCalcFlags <- FALSE
    }
    if(isWrtLine){
      derivCalcFlags <- c(derivCalcFlags, TRUE)
    }
    
    chainRuleHessianList[[i]] <- array(0, dim = c(sum(sapply(wrtLineInfo, function(x){return(x$lineSize)})),
                                                  sum(sapply(wrtLineInfo, function(x){return(x$lineSize)})),
                                                  length(derivList$value)))
    
    for(derivCalcFlag in derivCalcFlags){
      if(derivCalcFlag == FALSE){
        chainRuleDerivList[[i]] <- matrix(0, 
                                          nrow = length(derivList$value),
                                          ncol = sum(sapply(wrtLineInfo, function(x){return(x$lineSize)})))
      }
      else{
        chainRuleDerivList[[i]] <- matrix(0,
                                          nrow = length(eval(calcWithArgsCall[[3]])),
                                          ncol = sum(sapply(wrtLineInfo, function(x){return(x$lineSize)})))
      }
      ## Iterate over all wrt params.
      for(j in seq_along(wrtLineInfo)){
        if(wrtLineInfo[[j]]$lineNum == i){
          if(derivCalcFlag == FALSE){
            ## If this line is included in output, add the derivative of this node's calc. function wrt itself.
            outGradient[,wrtLineInfo[[j]]$lineIndices] <- outGradient[,wrtLineInfo[[j]]$lineIndices] + derivList$gradient[1:wrtLineInfo[[j]]$lineSize] 
          }
          else{
            ## If this line (i) represents this wrt param (j), set chain rule derivs to 1 (or identity).
            chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices] <- diag(length(eval(calcWithArgsCall[[3]])))
          }
        }
        ## If this line does not represent this wrt param, calc derivatives wrt this param via chain rule.
        else{
          thisArgIndex <- 0
          ## Iterate over this line's parent nodes.
          for(k in seq_along(derivInfo[[2]][[i]])){
            if(derivInfo[[2]][[i]][[k]][1] > 0){
              ## Get derivs of this parent node.
              parentGradientsList <- chainRuleDerivList[derivInfo[[2]][[i]][[k]]]
              parentGradients <- parentGradientsList[[1]]
              if(length(parentGradientsList) > 1){
                for(i_listEntry in 2:length(parentGradientsList)){
                  parentGradients <- rbind(parentGradients, parentGradientsList[[i_listEntry]])
                }
              }
              ## Calculate derivs of this node (i) wrt this parameter (j) for this parent node (k) via chain rule. 
              chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices] <- chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices] +
                derivList$gradient[,(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]), drop = FALSE]%*%parentGradients[,wrtLineInfo[[j]]$lineIndices, drop = FALSE]
            }
            thisArgIndex <- thisArgIndex + argSizeInfo[k]
          }
          if(derivCalcFlag == FALSE){
            ## If this line is included in output, add the derivative of this line (i) wrt this param (j).
            outGradient[, wrtLineInfo[[j]]$lineIndices] <- outGradient[, wrtLineInfo[[j]]$lineIndices]  +  
              chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices]
          }
        }
      }
    }
    if(isNodeLine){
      outValue <- outValue + derivList$value
    }
    
    ##HESSIAN BELOW
    ## For each node i, we want array of d^2i/dWrtnode_1 dWrtnode_2
    ## in above iteration:
    ## in the j loop, but oustide of the k loop.
    ## iterate over j_2 as well.  for each j, j_2, first iterate over k.  Will need this gradient and parent hessians.
    ##  then iterate over k and k_2.  will need this hessian and parent gradients
    # for(j in seq_along(wrtLineInfo)){
    #   for(j_2 in j:length(wrtLineInfo)){
    #       if(wrtLineInfo[[j]]$lineNum == i && wrtLineInfo[[j_2]]$lineNum){
    #         chainRuleHessianList[[i]][wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j]]$lineIndices, ] <- 1   
    #         ## If this line is also the line of a calculate node, add the derivative of this node's calc. function wrt itself.
    #         if(i %in% nodesLineNums){
    #           outHessian[wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j]]$lineIndices] <- outHessian[wrtLineInfo[[j]]$lineIndices, wrtLineInfo[[j]]$lineIndices] + 
    #             derivList$hessian[1:wrtLineInfo[[j]]$lineSize, 1:wrtLineInfo[[j]]$lineSize,] 
    #         }
    #       }
    #       else{
    #       thisArgIndex <- 0
    #     ## Iterate over this line's parent nodes.
    #       for(k in seq_along(derivInfo[[2]][[i]])){
    #         if(derivInfo[[2]][[i]][[k]][1] > 0){
    #           parentGradients <- chainRuleDerivList[derivInfo[[2]][[i]][[k]]]
    #           if(length(parentGradients) > 1){
    #             parentGradients <- t(sapply(parentGradients, as.numeric))
    #           }
    #           else{
    #             parentGradients <- parentGradients[[1]]
    #           }
    #           parentHessian <- chainRuleHessianList[derivInfo[[2]][[i]][[k]]]
    #           if(length(parentHessian) > 1){
    #             parentGradients <- t(sapply(parentGradients, as.numeric))
    #           }
    #           else{
    #             parentGradients <- parentGradients[[1]]
    #           }
    #           
    #         }
    #       }
    #       }
    #   ## If this line does not represent this wrt param, calc derivatives wrt this param via chain rule.
    #   else{
    #     thisArgIndex <- 0
    #     ## Iterate over this line's parent nodes.
    #     for(k in seq_along(derivInfo[[2]][[i]])){
    #       if(derivInfo[[2]][[i]][[k]][1] > 0){
    #         ## Get derivs of this parent node.
    #         parentGradients <- chainRuleDerivList[derivInfo[[2]][[i]][[k]]]
    #         if(length(parentGradients) > 1){
    #           parentGradients <- t(sapply(parentGradients, as.numeric))
    #         }
    #         else{
    #           parentGradients <- parentGradients[[1]]
    #         }
    #         ## Calculate derivs of this node (i) wrt this parameter (j) for this parent node (k) via chain rule. 
    #         chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices] <- chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices] +
    #           derivList$gradient[,(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]), drop = FALSE]%*%parentGradients[,wrtLineInfo[[j]]$lineIndices, drop = FALSE]
    #       }
    #       thisArgIndex <- thisArgIndex + argSizeInfo[k]
    #       if(i %in% nodesLineNums){
    #         ## If this line is also the line of a calculate node, add the derivative of this line (i) wrt this param (j).
    #         outGradient[wrtLineInfo[[j]]$lineIndices] <- outGradient[wrtLineInfo[[j]]$lineIndices]  +  
    #           chainRuleDerivList[[i]][,wrtLineInfo[[j]]$lineIndices]
    #       }
    #     }
    #   }
    # }
  }
  return(ADNimbleList$new(value = outValue,
                          gradient = outGradient))
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
  nodesLineNums <- sapply(stochNodes, function(x){which(x == derivInfo[[1]])})
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
    return(rDeriv_CalcNodes(model, nfv, derivInfo, nodesLineNums, wrtLineInfo))
  }	
}

