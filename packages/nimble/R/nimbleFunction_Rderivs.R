

makeSingleArgWrapper <- function(nf, wrt, fxnEnv) {
  ## The code below matches the wrt arguemnts provided to the call to nimDerivs
  ## to the formal arguments taken by the nimbleFunction, and gets dimension and
  ## indexing information about each of these wrt arguments. 
  formalNames <- formalArgs(eval(nf[[1]], envir = fxnEnv)@.Data)
  wrtNames <- strsplit(wrt, '\\[')
  wrtArgIndices <- c(sapply(wrtNames, function(x){return(which(x[1] == formalNames))}))
  flatteningInfo <- list()
  thisIndex <- 1
  for(i in seq_along(wrt)){
    arg <- nf[[wrtArgIndices[i] + 1]]
    indText <- ''
    if(length(wrtNames[[i]]) > 1){
      indText <-  paste0('[', wrtNames[[i]][[2]])
    }
    argSym <- parse(text = paste0(deparse(arg), indText))[[1]]
    dimInfo <- dim(eval(argSym, envir = fxnEnv))
    if(is.null(dimInfo)) dimInfo <- length(eval(argSym, envir = fxnEnv))
    lineInds <- thisIndex:(thisIndex + prod(dimInfo) - 1)
    argDimInfo <- dim(eval(nf[[wrtArgIndices[i] + 1]],envir = fxnEnv))
    flatteningInfo[[i]] <- list()
    flatteningInfo[[i]][[1]] <- argDimInfo
    flatteningInfo[[i]][[2]] <- indText
    flatteningInfo[[i]][[3]] <- lineInds
    thisIndex <- thisIndex + prod(dimInfo)
  }
  ## The wrappedFun created below is a wrapper for a call to the nimFxn argument
  ## to nimDerivs.  The wrapped version takes a single vector argument x, as 
  ## required by the numDeriv package.  It then unpacks the elements of the 
  ## vector x, using them as arguments to the nimFxn as appropriate.
  wrappedFun <- function(x) {
    args <- list()
    for(i in 1:(length(nf)-1)){
      if(i %in% wrtArgIndices){
        thisWrtIndex <- which(i == wrtArgIndices)
        thisSize <- sum(sapply(flatteningInfo[thisWrtIndex], function(x){return(prod(x[[1]]))})) ## total length
        if(length(wrtNames[[thisWrtIndex[1]]]) > 1){
            args[[i]] <-  eval(nf[[i+1]], envir = fxnEnv)
            argBracketExpr <- parse(text = paste0('args[[i]]', flatteningInfo[[thisWrtIndex[1]]][[2]], ' <- x[flatteningInfo[[thisWrtIndex[1]]][[3]]]'))[[1]]
            eval(argBracketExpr)
            if(length(thisWrtIndex) > 1){
              for(j in 2:length(wrtNames[thisWrtIndex])){
                argBracketExpr <- parse(text = paste0('args[[i]]', flatteningInfo[[thisWrtIndex[j]]][[2]], ' <- x[flatteningInfo[[thisWrtIndex[j]]][[3]]]'))[[1]]
                eval(argBracketExpr)
              }
            }
          }
          else{
            args[[i]] <-  x[flatteningInfo[[thisWrtIndex[1]]][[3]]]
          }
          if(length( flatteningInfo[[thisWrtIndex[1]]][[1]]) > 1 ) dim(  args[[i]]) <-  flatteningInfo[[thisWrtIndex[1]]][[1]]
      }
      else{
       args[[i]] <- nf[[i+1]] 
      }
    }
    c(do.call(paste(nf[[1]]), args, envir = fxnEnv))
  }
  ## The makeSingleArg function created below extracts the values of the 
  ## variables specified in the wrt arguemnt to nimDerivs and concatenates them
  ## into a single vector.
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
#' @param order an integer vector with values within the set {0, 1, 2}, corresponding to whether the function value, gradient, and Hessian should be returned respectively.  Defaults to \code{c(0, 1, 2)}.
#' @param dropArgs a vector of integers specifying any arguments to \code{nimFxn} that derivatives should not be taken with respect to.  For example, \code{dropArgs = 2} means that the second argument to \code{nimFxn} will not have derivatives taken with respect to it.  Defaults to an empty vector. 
#' @param wrt a character vector of node names to take derivatives with respect to.  If left empty, derivatives will be taken with respect to all arguments to \code{nimFxn}.
#' 
#' @details Derivatives for uncompiled nimbleFunctions are calculated using the
#' \code{numDeriv} package.  If this package is not installed, an error will
#' be issued.  Derivatives for matrix valued arguments will be returned in column-major order.
#' 
#' @return a \code{nimbleList} with elements \code{value}, \code{gradient}, and \code{hessian}.
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
  }
  if(!is.null(dropArgs)){
    removeArgs <- which(wrt == dropArgs)
    if(length(removeArgs) > 0)
      wrt <- wrt[-removeArgs]
  }
  derivFxnList <- makeSingleArgWrapper(derivFxnCall, wrt, fxnEnv)
  singleArg <- derivFxnList[[2]]()
  
  hessianFlag <- 2 %in% order
  gradientFlag <- 1 %in% order
  valueFlag <- 0 %in% order
  if(hessianFlag){
    derivList <- genD(derivFxnList[[1]], singleArg)
    if(valueFlag) outVal <- derivList$f0 
    if(gradientFlag) outGrad <- derivList$D[,1:derivList$p, drop = FALSE]
    outHessVals <- derivList$D[,(derivList$p + 1):dim(derivList$D)[2], drop = FALSE]
    outHess <- array(NA, dim = c(derivList$p, derivList$p, length(derivList$f0)))
    singleDimMat <- matrix(NA, nrow = derivList$p, ncol = derivList$p)
    singleDimMatUpperTriDiag <- upper.tri(singleDimMat, diag = TRUE)
    for(outDim in seq_along(derivList$f0)){
      singleDimMat[singleDimMatUpperTriDiag] <- outHessVals[outDim,]
      singleDimMat[lower.tri(singleDimMat)] <-   t(singleDimMat)[lower.tri(singleDimMat)]
      outHess[,,outDim] <- singleDimMat
    }
  }
  else if(gradientFlag){
    if(valueFlag) outVal <- nimFxn 
    outGrad <- jacobian(derivFxnList[[1]], singleArg)
  }
  else if(valueFlag){
    outVal <- nimFxn
  }
  outList <- nimble:::ADNimbleList$new()
  if(valueFlag) outList$value = outVal
  if(gradientFlag) outList$gradient = outGrad
  if(hessianFlag) outList$hessian = outHess
  return(outList)
}
  

nimDerivs_calculate <- function(model, nodes = NA, nodeFxnVector = NULL, nodeFunctionIndex = NULL, order, wrtPars){
  if(all(order == 0)) return(calculate(model, nodes))
  if(!inherits(model, 'modelBaseClass') ){
    stop('model argument must be a nimbleModel.')
  }
  if(missing(nodes) ) 
    nodes <- model$getMaps('nodeNamesLHSall')
  if(!is.null(nodeFxnVector)){
    stop('nfv case of nimDerivs_calculate not implemented yet.')
  }
  wrtParsDeps <- model$getDependencies(wrtPars)
  nodes <- model$expandNodeNames(nodes, sort = TRUE)
  if(!all(nodes %in% wrtParsDeps)){
    warning('not all calculate nodes depend on a wrtNode')
  }
  # nodesAndWrt <- model$expandNodeNames(c(nodes, model$expandNodeNames(wrtPars)), sort = TRUE)
  # nfv <- nodeFunctionVector(model, nodesAndWrt, sortUnique = FALSE)
  derivInfo <- nimDerivsInfoClass$new(wrtPars, nodes, model)
  hessianFlag <- 2 %in% order
  gradientFlag <- 1 %in% order || hessianFlag ## note that to calculate Hessians using chain rule, must have gradients. 
  if(gradientFlag && !(1 %in% order)) order <- c(order, 1)
  valueFlag <- 0 %in% order
  indexingInfo <- nodeFunctionVector(model, model$expandNodeNames(c(nodes, 
                                     model$expandNodeNames(wrtPars)), 
                                     sort = TRUE), 
                                     sortUnique = FALSE)$indexingInfo
  declIDs <- indexingInfo$declIDs
  unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
  chainRuleDerivList <- list()
  if(hessianFlag) chainRuleHessianList <- list()
  ## totalOutWRTSize is the sum of the lengths of all ouput wrt parameters.
  totalOutWRTSize <- sum(sapply(derivInfo$wrtToIndices, function(x){return(length(x))}))
  ## outDerivList will be returned from this function.
  outDerivList <- nimble:::ADNimbleList$new()
  if(valueFlag) outDerivList$value = 0
  outDerivList$gradient = matrix(0, ncol = totalOutWRTSize, nrow = 1)
  if(hessianFlag) outDerivList$hessian = array(0, dim = c(totalOutWRTSize, totalOutWRTSize, 1))
  totalWRTSize <-  sum(sapply(derivInfo$wrtLineSize, function(x){return(x)}))
  allWRTLines <-derivInfo$wrtLineNums
  ## Below we get the start and end indices of all wrt parameters.
  ## For example, if we are taking derivatives with respect to 'a' and 'b', both of 
  ## length 2, outIndexStartPoints would be c(1, 3) and outIndexEndPoints would be c(2,4)
  calcNodesLineNums <- which(derivInfo$calcNodeIndicators == 1)
  for(i in seq_along(declIDs)){
    if(length(calcNodesLineNums) > 0){
      declID <- declIDs[i]
      isDeterm <- derivInfo$stochasticIndicators[i] == 0
      thisWrtLine <- which(allWRTLines == i)
      isWrtLine <-  derivInfo$wrtNodeIndicators[i] == 1
      isCalcNodeLine <- i %in% derivInfo$calcNodeIndicators[i]
      ## Below shouldn't be necessary in c++ chain rule code. 
      thisNodeSize <- length(eval(calcWithArgsCall[[3]]))   
      ## If this node is either a calulate node or a deterministic dependency of a wrt node,
      ## we need to take derivatives of its calculateWithArgs function.  The derivative function
      ## call is evaluated below.
      if(isCalcNodeLine){
        
        ## Below we construct two lists:
        ## parentGradients, a list of all the gradients of the parent nodes of each argument of this node.
        ## parentHessians, a list of all the hessians of the parent nodes of this node.  
        parentGradients <- vector('list', length = length(derivInfo[[2]][[i]]))
        if(hessianFlag) parentHessians <- vector('list', length = length(derivInfo[[2]][[i]]))
        thisWRTArgs <- c()
        argSizeInfo <- numeric(length(derivInfo[[2]][[i]]))
        for(k in seq_along(derivInfo[[2]][[i]])){
          if(k == 1 && isWrtLine){
            ## The first argument (k = 1) of a node's calculateWithArgs function will always be the node itself.
            ## If this node is a wrt node, we want to set the parent gradient of the first arg (the derivative of this node wrt itself)
            ## to the identity.
            parentGradients[[k]] <- matrix(0, nrow = derivInfo$wrtLineSize[[thisWrtLine]],
                                           ncol = totalWRTSize)
            parentGradients[[k]][, derivInfo$wrtLineIndices[[thisWrtLine]]] <- diag(derivInfo$wrtLineSize[[thisWrtLine]])
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
            if(hessianFlag){
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
        }
        if(!is.null(thisWRTArgs)){
          derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS, WRT),
                                       list(CALCCALL = calcWithArgsCall,
                                            DERIVORDERS = order,
                                            DROPARGS = 'INDEXEDNODEINFO_',
                                            WRT = derivInfo$lineWRTArgsAsCharacters[[i]])))
          if(isDeterm){
            derivList$value <- 0
            model$nodeFunctions[[declID]]$calculate(unrolledIndicesMatrixRow)
          }
          
          ## The derivOutputFlag determines whether the derivatives of this node (node i): 
          ## should be calculated for inclusion in the chain rule output (TRUE),
          ## should be calculated for later use in the chain rule (FALSE), 
          derivOutputFlag <- if(isDeterm) FALSE else TRUE
          
          if(derivOutputFlag == TRUE){
            ## If these derivatives will be included in output, we are taking derivative of a log prob. calculation,
            ## so ouput will be length 1, and input args will be all wrt args.
            chainRuleDerivList[[i]] <- matrix(0, nrow = 1, ncol = totalWRTSize)
            if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWRTSize, totalWRTSize, 1))
          }
          else{
            ## otherwise, we are taking derivative of a node value calculation (not log prob. calculation).
            ## so ouput will be the length of this node, and input args will be all wrt args.
            chainRuleDerivList[[i]] <- matrix(0, nrow = thisNodeSize, ncol = totalWRTSize)
            if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWRTSize, totalWRTSize, thisNodeSize))
          }
          ## Iterate over all wrt params.
          for(j in seq_along(allWRTLines)){
            thisArgIndex <- 0
            ## Iterate over this line's parent nodes.
            for(k in seq_along(derivInfo[[2]][[i]])){
              if(!is.null(parentGradients[[k]])){
                ## Calculate derivs of this node (i) wrt this parameter (j) for this parent node (k) via chain rule. 
                chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[j]]] <- chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[j]]] +
                  derivList$gradient[,(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]), drop = FALSE]%*%parentGradients[[k]][,derivInfo$wrtLineIndices[[j]], drop = FALSE]
                
              }
              thisArgIndex <- thisArgIndex + argSizeInfo[k]
            }
            if(derivOutputFlag == TRUE){
              ## If this line is included in output, add the derivative of this line (i) wrt this param (j).
              outDerivList$gradient[, derivInfo$wrtToIndices[[j]]] <- outDerivList$gradient[, derivInfo$wrtToIndices[[j]]]  +  
                chainRuleDerivList[[i]][,derivInfo$wrtFromIndices[[j]]]
            }
            if(hessianFlag){
              ## The Hessian is calculated below using Faà di Bruno's formula.
              ## Second iteration over wrt parameters 
              for(j_2 in j:length(allWRTLines)){
                thisArgIndex <- 0
                ## Iterate over this line's parent nodes.
                for(k in seq_along(derivInfo[[2]][[i]])){
                  if(!is.null(parentHessians[[k]])){
                    for(dim1 in derivInfo$wrtLineIndices[[j]]){
                      for(dim2 in derivInfo$wrtLineIndices[[j_2]]){
                        chainRuleHessianList[[i]][dim1, dim2, ] <- chainRuleHessianList[[i]][dim1, dim2, ] +
                          c(derivList$gradient[ ,(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]), drop = FALSE]%*%parentHessians[[k]][dim1, dim2, , drop = FALSE])
                      }
                    }
                  }
                  thisArgIndex_2 <- 0
                  for(k_2 in seq_along(derivInfo[[2]][[i]])){
                    if(!is.null(parentGradients[[k]])){
                      if(!is.null(parentGradients[[k_2]])){
                        for(dim3 in 1:dim(derivList$hessian)[3]){
                          chainRuleHessianList[[i]][derivInfo$wrtLineIndices[[j]], derivInfo$wrtLineIndices[[j_2]], dim3] <- chainRuleHessianList[[i]][derivInfo$wrtLineIndices[[j]], derivInfo$wrtLineIndices[[j_2]], dim3] +
                            t(parentGradients[[k]][, derivInfo$wrtLineIndices[[j]], drop = FALSE])%*%derivList$hessian[(thisArgIndex + 1):(thisArgIndex + argSizeInfo[k]),(thisArgIndex_2 + 1):(thisArgIndex_2 + argSizeInfo[k_2]), dim3]%*%
                            parentGradients[[k_2]][, derivInfo$wrtLineIndices[[j_2]], drop = FALSE]
                        }
                      }
                      
                    }
                    thisArgIndex_2 <- thisArgIndex_2 + argSizeInfo[k_2]
                  }
                  thisArgIndex <- thisArgIndex + argSizeInfo[k]
                }
                if(derivOutputFlag == TRUE){
                  ## If this line is included in output, add the Hessian of this line (i) wrt this param #1 (j) and this param #2 (j_2).
                  outDerivList$hessian[derivInfo$wrtToIndices[[j]], derivInfo$wrtToIndices[[j_2]], ] <-  outDerivList$hessian[derivInfo$wrtToIndices[[j]], derivInfo$wrtToIndices[[j_2]], ]   +
                    chainRuleHessianList[[i]][derivInfo$wrtFromIndices[[j]], derivInfo$wrtFromIndices[[j_2]],]
                }
              }
            }
          }
        }
        else{ 
          if(valueFlag) derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS),
                                                     list(CALCCALL = calcWithArgsCall,
                                                          DERIVORDERS = c(0),
                                                          DROPARGS = 'INDEXEDNODEINFO_')))
          chainRuleDerivList[[i]] <- matrix(0, nrow = thisNodeSize, ncol = totalWRTSize)
          if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWRTSize, totalWRTSize, thisNodeSize))
        }
      }
      
      if(!isDeterm){
        ## If this is a wrt node, we need to set the chainRule lists appropriately so that the chain rule
        ## will work for dependent nodes of this node.  That means taking the first and second derivs of the
        ## function f(x) = x, which will be the identity matrix and 0 respectively.
        chainRuleDerivList[[i]] <- matrix(0,nrow = thisNodeSize, ncol = totalWRTSize)
        if(isWrtLine)   chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[thisWrtLine]]] <- diag(thisNodeSize)
        if(isWrtLine && hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWRTSize, totalWRTSize, thisNodeSize))
        if(isCalcNodeLine){
          if(valueFlag) outDerivList$value <- outDerivList$value + derivList$value
        }
      }
    }
  }
  
  ## Reflect hessian across the diagonal
  if(hessianFlag){
    upperTriHess <- outDerivList$hessian[,,1]
    upperTriHess[lower.tri(upperTriHess)] <-   t(upperTriHess)[lower.tri(upperTriHess)]
    outDerivList$hessian[,,1] <-  upperTriHess
  }
  return(outDerivList)
}

