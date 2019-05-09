
## Creates a wrapper function for a call to calculate(model, nodes).
## Arguments:
##  calcNodeName: A character vector of node names that will be used as the 
##                nodes argument in a call to calculate(model, nodes).
##  wrtName:      A character vector of node names that will be arguments to the
##                wrapper function.
##
## The returned wrapper function takes as arguments values for all parameters
## specified by the 'wrtName' argument, and returns the value of a call to
## calculate(model, nodes) evaluated at those parameter values.  The 'nodes' 
## argument to the calculate function is specified by the 'calcNodeName'
## argument to makeADCalcWrapperFunction.
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
                                       list(I = i, ARGNAME = names(
                                         calcFunctionArgNames)[i+1] ))
  }
  returnLine <- substitute(return(outVal))
  calcWrapperFunction <- function(){}
  body(calcWrapperFunction) <- putCodeLinesInBrackets(c(argSaveLines, callLine, 
                                                        argReplaceLines, 
                                                        returnLine))
  formals(calcWrapperFunction) <- calcFunctionArgNames
  return(calcWrapperFunction)
}

makeSingleArgWrapper <- function(nf, wrt, fxnEnv) {
  ## The code below matches the wrt arguemnts provided to the call to nimDerivs
  ## to the formal arguments taken by the nimbleFunction, and gets dimension and
    ## indexing information about each of these wrt arguments.
    
    fA <- formalArgs(eval(nf[[1]], envir = fxnEnv))
    wrtNames <- strsplit(wrt, '\\[')
    wrtArgIndices <- c(sapply(wrtNames, function(x){
        return(which(x[1] == fA))}))
    flatteningInfo <- list()
    thisIndex <- 1
    for(i in seq_along(wrt)){
        arg <- nf[[wrtArgIndices[i] + 1]]
        indText <- ''
        if(length(wrtNames[[i]]) > 1){
            indText <-  paste0('[', wrtNames[[i]][[2]])
        }
        argSym <- parse(text = paste0(paste0(deparse(arg), collapse = ''),
                                      indText))[[1]]
        dimInfo <- dimOrLength(eval(argSym, envir = fxnEnv))
        ##if(is.null(dimInfo)) dimInfo <- length(eval(argSym, envir = fxnEnv))
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
                thisSize <- sum(sapply(flatteningInfo[thisWrtIndex],
                                       function(x){return(prod(x[[1]]))})) #total length
                if(length(wrtNames[[thisWrtIndex[1]]]) > 1){
                    args[[i]] <-  eval(nf[[i+1]], envir = fxnEnv)
                    argBracketExpr <- parse(
                        text = paste0('args[[i]]', flatteningInfo[[thisWrtIndex[1]]][[2]],
                                      ' <- x[flatteningInfo[[thisWrtIndex[1]]][[3]]]')
                    )[[1]]
                    eval(argBracketExpr)
                    if(length(thisWrtIndex) > 1){
                        for(j in 2:length(wrtNames[thisWrtIndex])){
                            argBracketExpr <- parse(
                                text = paste0('args[[i]]', 
                                              flatteningInfo[[thisWrtIndex[j]]][[2]],
                                              ' <- x[flatteningInfo[[thisWrtIndex[j]]][[3]]]')
                            )[[1]]
                            eval(argBracketExpr)
                        }
                    }
                }
                else{
                    args[[i]] <-  x[flatteningInfo[[thisWrtIndex[1]]][[3]]]
                }
                if(length( flatteningInfo[[thisWrtIndex[1]]][[1]]) > 1 )
                    dim(  args[[i]]) <-  flatteningInfo[[thisWrtIndex[1]]][[1]]
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
                singleArg <- c(singleArg, c(
                                              eval(parse(text =  paste0('arg', flatteningInfo[[i]][[2]]))[1])))
            }
            else{
                singleArg <- c(singleArg, c(eval(nf[[wrtArgIndices[i]+1]],
                                                 envir = fxnEnv)))
            }
        }
        return(singleArg)
    }
    return(list(wrappedFun, makeSingleArg))
}


#' Nimble Derivatives
#' 
#' EXPERIMENTAL Computes the value, Jacobian, and Hessian of a given  
#' \code{nimbleFunction} method.  
#' 
#' @param nimFxn a call to a \code{nimbleFunction} method with arguments 
#' included.  Can also be a call to  \code{model$calculate(nodes)}, or to 
#' \code{calculate(model, nodes)}.
#' @param order an integer vector with values within the set {0, 1, 2}, 
#' corresponding to whether the function value, Jacobian, and Hessian should be
#'  returned respectively.  Defaults to \code{c(0, 1, 2)}.
#' @param dropArgs a vector of integers specifying any arguments to 
#' \code{nimFxn} that derivatives should not be taken with respect to.  For 
#' example, \code{dropArgs = 2} means that the second argument to \code{nimFxn}
#' will not have derivatives taken with respect to it.  Defaults to an empty
#' vector. 
#' @param wrt a character vector of either: names of function arguments 
#' (if taking derivatives of a \code{nimbleFunction} method), or node names 
#' (if taking derivatives of \code{model$calculate(nodes)}) to take derivatives 
#' with respect to.  If left empty, derivatives will be taken with respect to 
#' all arguments to \code{nimFxn}.
#' @param silent a logical argument that determines whether warnings will be
#' displayed.
#' @details Derivatives for uncompiled nimbleFunctions are calculated using the
#' \code{numDeriv} package.  If this package is not installed, an error will
#' be issued.  Derivatives for matrix valued arguments will be returned in 
#' column-major order.
#' 
#' @return a \code{nimbleList} with elements \code{value}, \code{jacobian},
#' and \code{hessian}.
#' 
#' @examples 
#' 
#' \dontrun{
#' model <- nimbleModel(code = ...)
#' calcDerivs <- nimDerivs(model$calculate(model$getDependencies('x')),
#'  wrt = 'x')
#' }
#' 
#' @export
nimDerivs <- function(nimFxn = NA, order = nimC(0,1,2), dropArgs = NA,
                      wrt = NULL, silent = TRUE){
  fxnEnv <- parent.frame()
  fxnCall <- match.call(function(nimFxn, order, dropArgs, wrt
                                 ){})
                                        #,chainRuleDerivs){})
  if(is.null(fxnCall[['order']])) fxnCall[['order']] <- order
  derivFxnCall <- fxnCall[['nimFxn']]
  wrtNames <- sapply(wrt, function(x){strsplit(x, '\\[')[[1]][1]})
  ## Below are two checks to see if the nimFxn argument is a model calculate 
  ## call. If so,  a wrapper function will be made and nimDerivs will be called 
  ## again on that wrapper function.
  if(deparse(derivFxnCall[[1]]) == 'calculate'){
    model <- eval(derivFxnCall[['model']], envir = fxnEnv)
    RCalcDerivsFunction <- makeDerivCalcWrapperFunction(eval(derivFxnCall[[
      'nodes']], envir = fxnEnv), wrtNames)
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
  if(length(derivFxnCall[[1]]) == 3 && deparse(derivFxnCall[[1]][[1]]) == '$'){
    if(deparse(derivFxnCall[[1]][[3]]) == 'calculate'){
      modelError <- tryCatch(get('modelDef', pos = eval(derivFxnCall[[1]][[2]],
                                                        envir = fxnEnv)),
                             error = function(e) e)
      if(!inherits(modelError, 'error')){
        model <-  eval(derivFxnCall[[1]][[2]], envir = fxnEnv)
        if(length(derivFxnCall) < 2){
          nodes <- model$getMaps('nodeNamesLHSall')
        }
        else{
          nodes <- eval(derivFxnCall[[2]], envir = fxnEnv)
        }
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
  fA <- formalArgs(eval(derivFxnCall[[1]], envir = fxnEnv))
  if(!all(wrtNames %in% fA)){
    stop('Error:  the wrt argument to nimDerivs() contains names that are not
         arguments to the nimFxn argument.')
  }
  ## libError <- try(library('numDeriv'), silent = TRUE)
  ## if(inherits(libError, 'try-error')){
  ##   stop("The 'numDeriv' package must be installed to use derivatives in
  ##        uncompiled nimbleFunctions.")
  ## }
  if(!require('numDeriv'))
      stop("The 'numDeriv' package must be installed to use derivatives in
         uncompiled nimbleFunctions.")

  if(is.null(wrt)){
    wrt <- fA
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
    outHessVals <- derivList$D[,(derivList$p + 1):dim(derivList$D)[2],
                               drop = FALSE]
    outHess <- array(NA, dim = c(derivList$p, derivList$p, length(derivList$f0)))
    singleDimMat <- matrix(NA, nrow = derivList$p, ncol = derivList$p)
    singleDimMatUpperTriDiag <- upper.tri(singleDimMat, diag = TRUE)
    for(outDim in seq_along(derivList$f0)){
      singleDimMat[singleDimMatUpperTriDiag] <- outHessVals[outDim,]
      singleDimMat[lower.tri(singleDimMat)] <-   t(singleDimMat)[
        lower.tri(singleDimMat)]
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
    ## wrtArgNames <- sapply(wrtArgs, function(x){strsplit(x, '\\[')[[1]][1]})
    wrtArgParsed <- lapply(wrtArgs, function(x) {
        pvl <- .Call(nimble:::makeParsedVarList, x)[[2]][[2]]
        if(is.name(pvl))
            list(argName = deparse(pvl),
                 argIndices = NULL)
        else
            list(argName = deparse(pvl[[2]]),
                 argIndices = as.list(pvl)[-c(1:2)])
    }
    )
    wrtArgNames <- unlist(lapply(wrtArgParsed, `[[`, 'argName'))
  nimFxnArgNames <- names(nimFxnArgs)
  wrtMatchArgs <- which(nimFxnArgNames %in% wrtArgNames)
  argNameCheck <- wrtArgNames %in% nimFxnArgNames
  ## Compare wrt args to actual function args and make sure no erroneous args
  ## are present.
  if(!all(argNameCheck)) stop('Incorrect names passed to wrt argument of
                                     nimDerivs: ', fxnName, 
                              ' does not have arguments named: ',
                              paste(wrtArgNames[!argNameCheck], 
                                    collapse = ', ' ), '.')
  ## Make sure all wrt args have type declarations (i.e., are not just names without types)
    nameCheck <- sapply(wrtMatchArgs, function(x){return(class(nimFxnArgs[[x]]))})
    if(any(nameCheck == 'name')) stop('Derivatives of ', fxnName, ' being taken 
                                    WRT an argument that does not have a valid type.')
    ## Make sure all wrt args have type double.
    doubleCheck <- sapply(wrtMatchArgs, function(x){
    return(deparse(nimFxnArgs[[x]][[1]]) == 'double')})
    if(!all(doubleCheck)) stop('Derivatives of ', fxnName, 
                               ' being taken WRT an argument that does not
                                    have type double().')
  ## Make sure that all wrt arg dims are < 2.
  fxnArgsDims <- sapply(wrtMatchArgs, function(x){
    if(length(nimFxnArgs[[x]]) == 1) outDim <- 0
    else outDim <- nimFxnArgs[[x]][[2]]
    names(outDim) <- names(nimFxnArgs)[x]
    return(outDim)})
  
  if(any(fxnArgsDims > 2)) stop('Derivatives cannot be taken WRT an argument
                                with dimension > 2')
  ## Determine sizes of each function arg.
  fxnArgsDimSizes <- lapply(nimFxnArgs, function(x){
    if(length(x) == 1) return(1)
    if(x[[2]] == 0) return(1)
    if(length(x) < 3) stop('Sizes of arguments to nimbleFunctions must be
                           explictly specified (e.g. x = double(1, 4)) in order
                           to take derivatives.')
    if(x[[2]] == 1){
      if(length(x[[3]]) == 1) return(x[[3]])
      else return(x[[3]][[2]])
    }
    else if(x[[2]] == 2) return(c(x[[3]][[2]], x[[3]][[3]]))
  })
  ## Same as above sizes, except that matrix sizes are reported as nrow*ncol
  ## instead of c(nrow, ncol).
  fxnArgsTotalSizes <- sapply(fxnArgsDimSizes, function(x){
    return(prod(x))
  }
  )
  fxnArgsTotalSizes <- c(0,fxnArgsTotalSizes)
  ## fxnArgsIndexVector is a named vector with the starting index of each
  ## argument to the nimFxn,
  ## if all arguments were flattened and put into a single vector.
  ## E.g. if nimFxn has arguments x = double(2, c(2, 2)), y = double(1, 3), 
  ## and z = double(0),
  ## then fxnArgsIndexVector will be:    
  ##   x  y  z
  ##   1  5  8
  fxnArgsIndexVector <- cumsum(fxnArgsTotalSizes) + 1
  names(fxnArgsIndexVector) <- c(names(fxnArgsIndexVector)[-1], '')
  fxnArgsIndexVector <- fxnArgsIndexVector[-length(fxnArgsIndexVector)]
  wrtArgsIndexVector <- c()
  for(i in seq_along(wrtArgNames)){
    if(fxnArgsDims[wrtArgNames[i]] == 0){
      wrtArgsIndexVector <- c(wrtArgsIndexVector, 
                              fxnArgsIndexVector[wrtArgNames[i]])
    }
    else if(fxnArgsDims[wrtArgNames[i]] == 1){
      hasIndex <- !is.null(wrtArgParsed[[i]]$argIndices) #grepl('^[a-z]*\\[', wrtArgs[i])
      if(hasIndex){
          if(length(wrtArgParsed[[i]]$argIndices) != 1)  stop(paste0('Incorrect indexing ',
                                                          'provided for wrt ',
                                                          'argument to ',
                                                          'nimDerivs(): ',
                                                          wrtArgs[i]))
           ##argIndices <- eval(parse(text = sub('\\]', '', sub('^[a-z]*\\[',
          ##'', wrtArgs[i]))))
          if(nimble:::is.blank(wrtArgParsed[[i]]$argIndices[[1]]))
              argIndices <- 1:fxnArgsTotalSizes[wrtArgNames[i]]
          else {
              argIndices <- eval(wrtArgParsed[[i]]$argIndices[[1]])
              if(any(argIndices < 1 | argIndices > fxnArgsTotalSizes[wrtArgNames[i]]))
                  stop(paste0('Index too large (or < 0) provided in wrt argument: ',
                              wrtArgs[i], ' for derivatives of ', fxnName))
          }
        wrtArgsIndexVector <- c(wrtArgsIndexVector, 
                                fxnArgsIndexVector[wrtArgNames[i]] +
                                  argIndices - 1)
      }
      else{
        wrtArgsIndexVector <- c(wrtArgsIndexVector, 
                                fxnArgsIndexVector[wrtArgNames[i]] +
                                  0:(fxnArgsTotalSizes[wrtArgNames[i]] - 1))
      }
    }
    else if(fxnArgsDims[wrtArgNames[i]] == 2){
      hasIndex <- !is.null(wrtArgParsed[[i]]$argIndices) #grepl('^[a-z]*\\[', wrtArgs[i])
      if(hasIndex){
##        argIndicesText <- strsplit(sub('\\]', '', sub('^[a-z]*\\[', '',
##                                                      wrtArgs[i])), ',',
##                                   fixed = TRUE)
        if(length(wrtArgParsed[[i]]$argIndices) != 2)  stop(paste0('Incorrect indexing ',
                                                          'provided for wrt ',
                                                          'argument to ',
                                                          'nimDerivs(): ',
                                                          wrtArgs[i]))
          if(nimble:::is.blank(wrtArgParsed[[i]]$argIndices[[1]]))
              argIndicesRows <- 1:fxnArgsDimSizes[[wrtArgNames[i]]][1]
          else {
              argIndicesRows <- eval(wrtArgParsed[[i]]$argIndices[[1]])
              if(any(argIndicesRows < 1 | argIndicesRows > fxnArgsDimSizes[[wrtArgNames[i]]][1]))
                  stop(paste0('A row index that is too large (or < 0) ',
                              ' was provided in wrt argument: ',
                              wrtArgs[i], ' for derivatives of ', fxnName))
          }
          if(nimble:::is.blank(wrtArgParsed[[i]]$argIndices[[2]]))
              argIndicesCols <- 1:fxnArgsDimSizes[[wrtArgNames[i]]][2]
          else { 
              argIndicesCols <- eval(wrtArgParsed[[i]]$argIndices[[2]])
              if(any(argIndicesCols < 1 | argIndicesCols > fxnArgsDimSizes[[wrtArgNames[i]]][2]))
                  stop(paste0('A column index that is too large (or < 0) ',
                              ' was provided in wrt argument: ',
                              wrtArgs[i], ' for derivatives of ', fxnName))
          
          }
        ## Column major ordering
        for(col in argIndicesCols){
          wrtArgsIndexVector <- c(wrtArgsIndexVector, 
                                  fxnArgsIndexVector[wrtArgNames[i]] + 
                                    (col - 1)*
                                    fxnArgsDimSizes[[wrtArgNames[i]]][2] +
                                    argIndicesRows - 1)
        }
      }
      else{
        wrtArgsIndexVector <- c(wrtArgsIndexVector,  
                                fxnArgsIndexVector[wrtArgNames[i]] +
                                  0:(fxnArgsTotalSizes[wrtArgNames[i]] - 1))
      }
    }
  }
  return(unname(wrtArgsIndexVector))
}

