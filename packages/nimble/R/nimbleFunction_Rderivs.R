## Creates a wrapper function for a call to calculate(model, nodes).
## Arguments:
##  calcNodeName: A character vector of node names that will be used as the nodes argument in a call to calculate(model, nodes).
##  wrtName:      A character vector of node names that will be arguments to the wrapper function.
##
## The returned wrapper function takes as arguments values for all parameters specified by the 'wrtName' argument,
## and returns the value of a call to calculate(model, nodes) evaluated at those parameter values.  The
## 'nodes' argument to the calculate function is specified by the 'calcNodeName' argument to makeADCalcWrapperFunction.
makeDerivCalcWrapperFunction <- function(calcNodeName, wrtName){
  calcFunctionArgNames <- list('model' = NA)
  for(i in seq_along(wrtName)){
    calcFunctionArgNames[[ wrtName[[i]][1]]] <- NA
  }
  argSaveLines <- list()
  argSaveLines[[1]] <- substitute(origVals <- list())
  for(i in seq_along(calcFunctionArgNames[-1])){
    argSaveLines[[i+1]] <- substitute({origVals[[I]] <- model[[ARGNAME]];
    model[[ARGNAME]] <- ARGNAMEEXPR},
    list(I = i,
         ARGNAME = names(calcFunctionArgNames)[i+1],
         ARGNAMEEXPR = as.name(names(calcFunctionArgNames)[i+1])))
  }
  callLine <-  list(substitute(outVal <- calculate(model, CALCNODENAME),
                               list(CALCNODENAME = calcNodeName)))
  argReplaceLines <- list()
  for(i in seq_along(calcFunctionArgNames[-1])){
    argReplaceLines[[i]] <- substitute(model[[ARGNAME]] <- origVals[[I]],
                                       list(I = i, ARGNAME = names(calcFunctionArgNames)[i+1] ))
  }
  returnLine <- substitute(return(outVal))
  calcWrapperFunction <- function(){}
  body(calcWrapperFunction) <- putCodeLinesInBrackets(c(argSaveLines, callLine, argReplaceLines, returnLine))
  formals(calcWrapperFunction) <- calcFunctionArgNames
  return(calcWrapperFunction)
}





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
    argSym <- parse(text = paste0(paste0(deparse(arg), collapse = ''), indText))[[1]]
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
#' EXPERIMENTAL Computes the value, Jacobian, and Hessian of a given  \code{nimbleFunction} method.  
#' 
#' @param nimFxn a call to a \code{nimbleFunction} method with arguments included.  Can also be a call to  \code{model$calculate(nodes)}, or to \code{calculate(model, nodes)}.
#' @param order an integer vector with values within the set {0, 1, 2}, corresponding to whether the function value, Jacobian, and Hessian should be returned respectively.  Defaults to \code{c(0, 1, 2)}.
#' @param dropArgs a vector of integers specifying any arguments to \code{nimFxn} that derivatives should not be taken with respect to.  For example, \code{dropArgs = 2} means that the second argument to \code{nimFxn} will not have derivatives taken with respect to it.  Defaults to an empty vector. 
#' @param wrt a character vector of either: names of function arguments (if taking derivatives of a \code{nimbleFunction} method), or node names (if taking derivatives of \code{model$calculate(nodes)}) to take derivatives with respect to.  If left empty, derivatives will be taken with respect to all arguments to \code{nimFxn}.
#' @param chainRuleDerivs a logical argument indicating whether to take derivatives of calls to \code{model$calculate(nodes)} using NIMBLE's implementation of the chain rule for testing purposes.  Defaults to \code{FALSE}.    
#' @param silent a logical argument that determines whether warnings will be displayed.
#' @details Derivatives for uncompiled nimbleFunctions are calculated using the
#' \code{numDeriv} package.  If this package is not installed, an error will
#' be issued.  Derivatives for matrix valued arguments will be returned in column-major order.
#' 
#' @return a \code{nimbleList} with elements \code{value}, \code{jacobian}, and \code{hessian}.
#' 
#' @examples 
#' 
#' \dontrun{
#' model <- nimbleModel(code = ...)
#' calcDerivs <- nimDerivs(model$calculate(model$getDependencies('x')), wrt = 'x')
#' }
#' 
#' @export
nimDerivs <- function(nimFxn = NA, order = nimC(0,1,2), dropArgs = NA, wrt = NULL, chainRuleDerivs = FALSE, silent = TRUE){
  fxnEnv <- parent.frame()
  fxnCall <- match.call(function(nimFxn, order, dropArgs, wrt, chainRuleDerivs){})
  if(is.null(fxnCall[['order']])) fxnCall[['order']] <- order
  derivFxnCall <- fxnCall[['nimFxn']]
  wrtNames <- sapply(wrt, function(x){strsplit(x, '\\[')[[1]][1]})
  ## Below are two checks to see if the nimFxn argument is a model calculate call.
  ## If so, one of two paths will be taken:
  ##  1) If chainRuleDerivs = FALSE, a wrapper function will be made and nimDerivs
  ##     will be called again on that wrapper function.
  ##  2) If chainRuleDerivs = TRUE, derivatives will be calculated using the chain rule,
  ##     via the nimDerivs_calculate function.
  if(deparse(derivFxnCall[[1]]) == 'calculate'){
    if(chainRuleDerivs == TRUE){
      derivFxnCall <- match.call(calculate, derivFxnCall)
      return(nimDerivs_calculate(model = eval(derivFxnCall[['model']], envir = fxnEnv),
                                 nodes = eval(derivFxnCall[['nodes']], envir = fxnEnv),
                                 order, wrt, silent))
    }
    else{
      model <- eval(derivFxnCall[['model']], envir = fxnEnv)
      RCalcDerivsFunction <- makeDerivCalcWrapperFunction(eval(derivFxnCall[['nodes']], envir = fxnEnv), wrtNames)
      argList <- list('model' = quote(model))
      for(k in seq_along(wrtNames)){
        argList[[wrtNames[[k]][1]]] <- model[[wrtNames[[k]][1]]]
      }
      argList <- c(list('RCalcDerivsFunction'), argList)
      fxnCall <- as.call(argList)
      fxnCall[[1]] <- quote(RCalcDerivsFunction)
      return(eval(substitute(nimDerivs(FXNCALL, wrt = WRT, order = order),
                             list(FXNCALL = fxnCall,
                                  WRT = wrt))))
    }
  }
  if(length(derivFxnCall[[1]]) == 3 && deparse(derivFxnCall[[1]][[1]]) == '$'){
    if(deparse(derivFxnCall[[1]][[3]]) == 'calculate'){
      modelError <- tryCatch(get('modelDef', pos = eval(derivFxnCall[[1]][[2]], envir = fxnEnv)), error = function(e) e)
      if(!inherits(modelError, 'error')){
        model <-  eval(derivFxnCall[[1]][[2]], envir = fxnEnv)
        if(length(derivFxnCall) < 2){
          nodes <- model$getMaps('nodeNamesLHSall')
        }
        else{
          nodes <- eval(derivFxnCall[[2]], envir = fxnEnv)
        }
        if(chainRuleDerivs == TRUE){
          return(nimDerivs_calculate(model = model,
                                     nodes = nodes,
                                     order = order, wrt = wrt))
        }
        else{
          RCalcDerivsFunction <- makeDerivCalcWrapperFunction(nodes, wrtNames)
          argList <- list('model' = quote(model))
          for(k in seq_along(wrtNames)){
            argList[[wrtNames[[k]][1]]] <- model[[wrtNames[[k]][1]]]
          }
          argList <- c(list('RCalcDerivsFunction'), argList)
          fxnCall <- as.call(argList)
          fxnCall[[1]] <- quote(RCalcDerivsFunction)
          return(eval(substitute(nimDerivs(FXNCALL, wrt = WRT, order = order),
                                 list(FXNCALL = fxnCall,
                                      WRT = wrt))))
        }
      }
    }
  }
  if(!all(wrtNames %in% formalArgs(eval(derivFxnCall[[1]], envir = fxnEnv)@.Data))){
    stop('Error:  the wrt argument to nimDerivs() contains names that are not arguments to the nimFxn argument.')
  }
  libError <- try(library('numDeriv'), silent = TRUE)
  if(inherits(libError, 'try-error')){
    stop("The 'numDeriv' package must be installed to use derivatives in uncompiled nimbleFunctions.")
  }
  if(is.null(wrt)){
    wrt <- formalArgs(eval(derivFxnCall[[1]], envir = fxnEnv)@.Data)
  }
  if(!is.na(dropArgs)){
    removeArgs <- which(wrt == dropArgs)
    if(length(removeArgs) > 0)
      wrt <- wrt[-removeArgs]
  }
  ## derivFxnList will be a list with two elements:
  ## The first element is a wrapper function to the nimFxn that the numDeriv 
  ## package will actually take derivatives of.
  ## The second element is a function that will take all arguments to the nimFxn
  ## and wrap them into a single vector argument.
  derivFxnList <- makeSingleArgWrapper(derivFxnCall, wrt, fxnEnv)
  singleArg <- derivFxnList[[2]]()
  hessianFlag <- 2 %in% order
  jacobianFlag <- 1 %in% order
  valueFlag <- 0 %in% order
  if(hessianFlag){
    ## If hessians are requested, derivatives taken using numDeriv's genD() 
    ## function.  After that, we extract the various derivative elements and 
    ## arrange them properly.
    derivList <- genD(derivFxnList[[1]], singleArg)
    if(valueFlag) outVal <- derivList$f0 
    if(jacobianFlag) outGrad <- derivList$D[,1:derivList$p, drop = FALSE]
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
  else if(jacobianFlag){
    ## If jacobians are requested, derivatives taken using numDeriv's jacobian() 
    ## function.  After that, we extract the various derivative elements and 
    ## arrange them properly.
    if(valueFlag) outVal <- nimFxn 
    outGrad <- jacobian(derivFxnList[[1]], singleArg)
  }
  else if(valueFlag){
    ## Otherwise just evaluate the function.
    outVal <- nimFxn
  }
  outList <- ADNimbleList$new()
  if(valueFlag) outList$value = outVal
  if(jacobianFlag) outList$jacobian = outGrad
  if(hessianFlag) outList$hessian = outHess
  return(outList)
}
  

nimDerivs_calculate <- function(model, nodes = NA, order, wrtPars, silent = TRUE){
  outDerivList <- ADNimbleList$new()
  if(all(order == 0)){
    outDerivList$value <- calculate(model, nodes)
    return(outDerivList)
  }
  if(!inherits(model, 'modelBaseClass') ){
    stop('model argument must be a nimbleModel.')
  }
  if(missing(nodes) ) 
    nodes <- model$getMaps('nodeNamesLHSall')
  wrtParsDeps <- model$getDependencies(wrtPars)
  nodes <- model$expandNodeNames(nodes, sort = TRUE)
  if(!all(nodes %in% wrtParsDeps)){
    if(!silent){
      warning('not all calculate nodes depend on a wrtNode')
    }
  }
  derivInfo <- nimDerivsInfoClass(wrtPars, nodes, model)
  hessianFlag <- 2 %in% order
  jacobianFlag <- 1 %in% order || hessianFlag ## note that to calculate Hessians using chain rule, must have jacobians. 
  if(jacobianFlag && !(1 %in% order)) order <- c(order, 1)
  valueFlag <- 0 %in% order
  nfv <-  nodeFunctionVector(model, model$expandNodeNames(c(nodes, 
                                                            model$expandNodeNames(wrtPars)), 
                                                          sort = TRUE), 
                             sortUnique = FALSE)
  declIDs <- nfv$indexingInfo$declIDs
  chainRuleDerivList <- list()
  if(hessianFlag) chainRuleHessianList <- list()
  ## totalOutWrtSize is the sum of the lengths of all ouput wrt parameters.
  totalOutWrtSize <- sum(sapply(derivInfo$wrtToIndices, function(x){return(length(x))}))
  ## outDerivList will be returned from this function.
  if(valueFlag) outDerivList$value = 0
  outDerivList$jacobian = matrix(0, ncol = totalOutWrtSize, nrow = 1)
  if(hessianFlag) outDerivList$hessian = array(0, dim = c(totalOutWrtSize, totalOutWrtSize, 1))
  totalWrtSize <-  sum(sapply(derivInfo$wrtLineSize, function(x){return(x)}))
  allWrtLines <- derivInfo$wrtLineNums
  ## Below we get the start and end indices of all wrt parameters.
  ## For example, if we are taking derivatives with respect to 'a' and 'b', both of 
  ## length 2, outIndexStartPoints would be c(1, 3) and outIndexEndPoints would be c(2,4)
  calcNodesLineNums <- which(derivInfo$calcNodeIndicators == 1)
  for(i in seq_along(derivInfo$allWrtAndCalcNodeNames)){
    if(length(calcNodesLineNums) > 0){
      declID <- declIDs[i]
      isDeterm <- derivInfo$stochNodeIndicators[i] == 0
      thisWrtLine <- which(allWrtLines == i)
      isWrtLine <-  derivInfo$wrtNodeIndicators[i] == 1
      isCalcNodeLine <- i %in% calcNodesLineNums
      ## Below shouldn't be necessary in c++ chain rule code.
      ## If this node is a calulate node,
      ## we need to take derivatives of its calculateWithArgs function.  The derivative function
      ## call is evaluated below.
      thisNodeSize <- length(values(model, derivInfo$allWrtAndCalcNodeNames[i]))  
      unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[nfv$indexingInfo$unrolledIndicesMatrixRows[i], ]
      if(isCalcNodeLine){
        calcWithArgs <- model$nodeFunctions[[declID]]$calculateWithArgs
        
        ## Below we construct two lists:
        ## parentjacobians, a list of all the jacobians of the parent nodes of each argument of this node.
        ## parentHessians, a list of all the hessians of the parent nodes of this node.  
        parentjacobians <- vector('list', length = length(derivInfo$parentIndicesList[[i]]))
        if(hessianFlag) parentHessians <- vector('list', length = length(derivInfo$parentIndicesList[[i]]))
        for(k in seq_along(derivInfo$parentIndicesList[[i]])){
          if(k == 1 && isWrtLine){
            ## The first argument (k = 1) of a node's calculateWithArgs function will always be the node itself.
            ## If this node is a wrt node, we want to set the parent jacobian of the first arg (the derivative of this node wrt itself)
            ## to the identity.
            parentjacobians[[k]] <- matrix(0, nrow = derivInfo$wrtLineSize[[thisWrtLine]],
                                           ncol = totalWrtSize)
            parentjacobians[[k]][, derivInfo$wrtLineIndices[[thisWrtLine]]] <- diag(derivInfo$wrtLineSize[[thisWrtLine]])
          }
          else if(derivInfo$parentIndicesList[[i]][[k]][1] > 0){
            ## Otherwise, if argument k has parents that depend on a wrt node, grab the parent jacobians (which will have already been calculated)
            ## and combine them into a single matrix.
            parentjacobiansList <- chainRuleDerivList[derivInfo$parentIndicesList[[i]][[k]]]
            parentjacobians[[k]] <- parentjacobiansList[[1]]
            if(length(parentjacobiansList) > 1){
              for(i_listEntry in 2:length(parentjacobiansList)){
                parentjacobians[[k]] <- rbind(parentjacobians[[k]], parentjacobiansList[[i_listEntry]])
              }
            }
            ## Similarly, grab the parent Hessians (which will have already been calculated)
            ## and combine them into a single array.
            if(hessianFlag){
              parentHessiansList <- chainRuleHessianList[derivInfo$parentIndicesList[[i]][[k]]]
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
        if(!is.na(derivInfo$lineWrtArgsAsCharacters[[i]][1])){
          derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS, WRT),
                                       list(CALCCALL = derivInfo$calcWithArgsCalls[[i]],
                                            DERIVORDERS = order,
                                            DROPARGS = 'INDEXEDNODEINFO_',
                                            WRT = derivInfo$lineWrtArgsAsCharacters[[i]])))
          if(isDeterm){
            derivList$value <- 0
            model$nodeFunctions[[declID]]$calculate(unrolledIndicesMatrixRow)
          }
          print(paste(i, "hessian:"))
          print(derivList$jacobian)
          ## The derivOutputFlag determines whether the derivatives of this node (node i): 
          ## should be calculated for inclusion in the chain rule output (TRUE),
          ## should be calculated for later use in the chain rule (FALSE), 
          derivOutputFlag <- if(isDeterm) FALSE else TRUE
          if(derivOutputFlag == TRUE){
            ## If these derivatives will be included in output, we are taking derivative of a log prob. calculation,
            ## so ouput will be length 1, and input args will be all wrt args.
            chainRuleDerivList[[i]] <- matrix(0, nrow = 1, ncol = totalWrtSize)
            if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, 1))
          }
          else{
            ## otherwise, we are taking derivative of a node value calculation (not log prob. calculation).
            ## so ouput will be the length of this node, and input args will be all wrt args.
            chainRuleDerivList[[i]] <- matrix(0, nrow = thisNodeSize, ncol = totalWrtSize)
            if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, thisNodeSize))
          }
          ## Iterate over all wrt params.
          for(j in seq_along(allWrtLines)){
            thisArgIndex <- 0
            ## Iterate over this line's parent nodes.
            for(k in seq_along(derivInfo$parentIndicesList[[i]])){
              if(!is.null(parentjacobians[[k]])){
                ## Calculate derivs of this node (i) wrt this parameter (j) for this parent node (k) via chain rule.
                chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[j]]] <- chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[j]]] +
                  derivList$jacobian[,(thisArgIndex + 1):(thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]), drop = FALSE]%*%parentjacobians[[k]][,derivInfo$wrtLineIndices[[j]], drop = FALSE]
                
              }
              thisArgIndex <- thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]
            }
            if(derivOutputFlag == TRUE){
              ## If this line is included in output, add the derivative of this line (i) wrt this param (j).
              outDerivList$jacobian[, derivInfo$wrtToIndices[[j]]] <- outDerivList$jacobian[, derivInfo$wrtToIndices[[j]]]  +  
                chainRuleDerivList[[i]][,derivInfo$wrtFromIndices[[j]]]
            }
            if(hessianFlag){
              ## The Hessian is calculated below using Faà di Bruno's formula.
              ## Second iteration over wrt parameters 
              for(j_2 in j:length(allWrtLines)){
                thisArgIndex <- 0
                ## Iterate over this line's parent nodes.
                for(k in seq_along(derivInfo$parentIndicesList[[i]])){
                  if(!is.null(parentHessians[[k]])){
                    for(dim1 in derivInfo$wrtLineIndices[[j]]){
                      for(dim2 in derivInfo$wrtLineIndices[[j_2]]){
                        chainRuleHessianList[[i]][dim1, dim2, ] <- chainRuleHessianList[[i]][dim1, dim2, ] +
                          c(derivList$jacobian[ ,(thisArgIndex + 1):(thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]), drop = FALSE]%*%parentHessians[[k]][dim1, dim2, , drop = FALSE])
                      }
                    }
                  }
                  thisArgIndex_2 <- 0
                  for(k_2 in seq_along(derivInfo$parentIndicesList[[i]])){
                    if(!is.null(parentjacobians[[k]])){
                      if(!is.null(parentjacobians[[k_2]])){
                        for(dim3 in 1:dim(derivList$hessian)[3]){
                          chainRuleHessianList[[i]][derivInfo$wrtLineIndices[[j]], derivInfo$wrtLineIndices[[j_2]], dim3] <- chainRuleHessianList[[i]][derivInfo$wrtLineIndices[[j]], derivInfo$wrtLineIndices[[j_2]], dim3] +
                            t(parentjacobians[[k]][, derivInfo$wrtLineIndices[[j]], drop = FALSE])%*%derivList$hessian[(thisArgIndex + 1):(thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]),(thisArgIndex_2 + 1):(thisArgIndex_2 + derivInfo$lineWrtArgSizeInfo[[i]][k_2]), dim3]%*%
                            parentjacobians[[k_2]][, derivInfo$wrtLineIndices[[j_2]], drop = FALSE]
                        }
                      }
                      
                    }
                    thisArgIndex_2 <- thisArgIndex_2 + derivInfo$lineWrtArgSizeInfo[[i]][k_2]
                  }
                  thisArgIndex <- thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]
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
                                                     list(CALCCALL = derivInfo$calcWithArgsCalls[[i]],
                                                          DERIVORDERS = c(0),
                                                          DROPARGS = 'INDEXEDNODEINFO_')))
          chainRuleDerivList[[i]] <- matrix(0, nrow = thisNodeSize, ncol = totalWrtSize)
          if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, thisNodeSize))
        }
      }
      
      if(!isDeterm){
        ## If this is a wrt node, we need to set the chainRule lists appropriately so that the chain rule
        ## will work for dependent nodes of this node.  That means taking the first and second derivs of the
        ## function f(x) = x, which will be the identity matrix and 0 respectively.
        chainRuleDerivList[[i]] <- matrix(0,nrow = thisNodeSize, ncol = totalWrtSize)
        if(isWrtLine)   chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[thisWrtLine]]] <- diag(thisNodeSize)
        if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, thisNodeSize))
        if(isCalcNodeLine){
          if(valueFlag){
            outDerivList$value <- outDerivList$value + derivList$value
          } 
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


convertWrtArgToIndices <- function(wrtArgs, nimFxnArgs, fxnName){
  ## If a vector of wrt args, get individual args.
  if(deparse(wrtArgs[[1]]) == 'nimC'){ 
    wrtArgs <- sapply(wrtArgs[-1], function(x){as.character(x)})
  }
  if(all(is.na(wrtArgs))){
    return(-1)
  }
  else if(is.null(wrtArgs[[1]])){
    wrtArgs <- names(nimFxnArgs)
  }
  else if(!is.character(wrtArgs)){
    wrtArgs <- deparse(wrtArgs)
  }
  ## Get names of wrt args with indexing removed.
  wrtArgNames <- sapply(wrtArgs, function(x){strsplit(x, '\\[')[[1]][1]})
  nimFxnArgNames <- names(nimFxnArgs)
  wrtMatchArgs <- which(nimFxnArgNames %in% wrtArgNames)
  argNameCheck <- wrtArgNames %in% nimFxnArgNames
  ## Compare wrt args to actual function args and make sure no erroneous args are present.
  if(any(argNameCheck != TRUE)) stop('Incorrect names passed to wrt argument of nimDerivs: ', fxnName, 
                                     ' does not have arguments named: ', paste(wrtArgNames[!argNameCheck], collapse = ', ' ), '.')
  ## Make sure all wrt args have type double.
  nameCheck <- sapply(wrtMatchArgs, function(x){return(class(nimFxnArgs[[x]]))})
  if(any(nameCheck == 'name')) stop('Derivatives of ', fxnName, ' being taken WRT an argument that does not have type double().')
  doubleCheck <- sapply(wrtMatchArgs, function(x){return(deparse(nimFxnArgs[[x]][[1]]) == 'double')})
  if(any(doubleCheck != TRUE)) stop('Derivatives of ', fxnName, ' being taken WRT an argument that does not have type double().')
  ## Make sure that all wrt arg dims are < 2.
  fxnArgsDims <- sapply(wrtMatchArgs, function(x){
    if(length(nimFxnArgs[[x]]) == 1) outDim <- 0
    else outDim <- nimFxnArgs[[x]][[2]]
    names(outDim) <- names(nimFxnArgs)[x]
    return(outDim)})
  
  if(any(fxnArgsDims > 2)) stop('Derivatives cannot be taken WRT an argument with dimension > 2')
  ## Determine sizes of each function arg.
  fxnArgsDimSizes <- lapply(nimFxnArgs, function(x){
    if(length(x) == 1) return(1)
    if(x[[2]] == 0) return(1)
    if(length(x) < 3) stop('Sizes of arguments to nimbleFunctions must be explictly specified (e.g. x = double(1, 4)) in order to take derivatives.')
    if(x[[2]] == 1){
      if(length(x[[3]]) == 1) return(x[[3]])
      else return(x[[3]][[2]])
    }
    else if(x[[2]] == 2) return(c(x[[3]][[2]], x[[3]][[3]]))
  })
  ## Same as above sizes, except that matrix sizes are reported as nrow*ncol instead of c(nrow, ncol).
  fxnArgsTotalSizes <- sapply(fxnArgsDimSizes, function(x){
    return(prod(x))
  }
  )
  fxnArgsTotalSizes <- c(0,fxnArgsTotalSizes)
  ## fxnArgsIndexVector is a named vector with the starting index of each argument to the nimFxn,
  ## if all arguments were flattened and put into a single vector.
  ## E.g. if nimFxn has arguments x = double(2, c(2, 2)), y = double(1, 3), and z = double(0),
  ## then fxnArgsIndexVector will be:    
  ##   x  y  z
  ##   1  5  8
  fxnArgsIndexVector <- cumsum(fxnArgsTotalSizes) + 1
  names(fxnArgsIndexVector) <- c(names(fxnArgsIndexVector)[-1], '')
  fxnArgsIndexVector <- fxnArgsIndexVector[-length(fxnArgsIndexVector)]
  wrtArgsIndexVector <- c()
  for(i in seq_along(wrtArgNames)){
    if(fxnArgsDims[wrtArgNames[i]] == 0){
      wrtArgsIndexVector <- c(wrtArgsIndexVector, fxnArgsIndexVector[wrtArgNames[i]])
    }
    else if(fxnArgsDims[wrtArgNames[i]] == 1){
      hasIndex <- grepl('^[a-z]*\\[', wrtArgs[i])
      if(hasIndex){
        argIndices <- eval(parse(text = sub('\\]', '', sub('^[a-z]*\\[', '', wrtArgs[i]))))
        if(is.null(argIndices)) argIndices <- 1:fxnArgsTotalSizes[wrtArgNames[i]]
        wrtArgsIndexVector <- c(wrtArgsIndexVector, fxnArgsIndexVector[wrtArgNames[i]] + argIndices - 1)
      }
      else{
        wrtArgsIndexVector <- c(wrtArgsIndexVector, fxnArgsIndexVector[wrtArgNames[i]] + 0:(fxnArgsTotalSizes[wrtArgNames[i]] - 1))
      }
    }
    else if(fxnArgsDims[wrtArgNames[i]] == 2){
      hasIndex <- grepl('^[a-z]*\\[', wrtArgs[i])
      if(hasIndex){
        argIndicesText <- strsplit(sub('\\]', '', sub('^[a-z]*\\[', '', wrtArgs[i])), ',', fixed = TRUE)
        if(length(argIndicesText[[1]]) != 2)  stop(paste0('Incorrect indexing provided for wrt argument to nimDerivs(): ', wrtArgs[i]))
        argIndicesRows <- eval(parse(text = argIndicesText[[1]][1]))
        argIndicesCols <- eval(parse(text = argIndicesText[[1]][2]))
        if(is.null(argIndicesRows)) argIndicesRows <- 1:fxnArgsDimSizes[[wrtArgNames[i]]][1]
        if(is.null(argIndicesCols)) argIndicesCols <- 1:fxnArgsDimSizes[[wrtArgNames[i]]][2]
        ## Column major ordering
        for(col in argIndicesCols){
          wrtArgsIndexVector <- c(wrtArgsIndexVector, fxnArgsIndexVector[wrtArgNames[i]] + (col - 1)*fxnArgsDimSizes[[wrtArgNames[i]]][2] + argIndicesRows - 1)
        }
      }
      else{
        wrtArgsIndexVector <- c(wrtArgsIndexVector,  fxnArgsIndexVector[wrtArgNames[i]] + 0:(fxnArgsTotalSizes[wrtArgNames[i]] - 1))
      }
    }
  }
  return(unname(wrtArgsIndexVector))
}

