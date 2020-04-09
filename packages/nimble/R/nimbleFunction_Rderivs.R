derivInfo <- function(x, ...) x

## dummy for recursing in sizeNimDerivs
nimDerivs_dummy <- nimbleFunction(
    run = function(call = void(),
                   order = integer(1),
                   dropArgs = void(),
                   wrt = integer(1)) {
    }
)

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
nimDerivs <- function(nimFxn = NA,
                      order = nimC(0,1,2),
                      dropArgs = NA,
                      wrt = NULL, calcNodes = NULL, static = FALSE){
  fxnEnv <- parent.frame()
  fxnCall <- match.call()
  if(is.null(fxnCall[['order']])) fxnCall[['order']] <- order
  derivFxnCall <- fxnCall[['nimFxn']]
  ## Below are two checks to see if the nimFxn argument is a model calculate 
  ## call. If so,  a wrapper function will be made and nimDerivs will be called 
  ## again on that wrapper function.
  model_calculate_case <- FALSE
  if(deparse(derivFxnCall[[1]]) == 'calculate'){
      model_calculate_case <- TRUE
  }
  if(length(derivFxnCall[[1]]) == 3 && deparse(derivFxnCall[[1]][[1]]) == '$'){
    if(deparse(derivFxnCall[[1]][[3]]) == 'calculate'){
      modelError <- tryCatch(get('modelDef', pos = eval(derivFxnCall[[1]][[2]],
                                                        envir = fxnEnv)),
                             error = function(e) e)
      if(!inherits(modelError, 'error')){
          model_calculate_case <- TRUE
      }
    }
  }
  ## Handle case of nimFxn argument being a calcNodes call for faster
  ## AD testing of model calculate calls. 
  model_calcNodes_case <- FALSE
  if(length(derivFxnCall[[1]]) == 3 && deparse(derivFxnCall[[1]][[1]]) == '$'){
      if(deparse(derivFxnCall[[1]][[3]]) == 'run'){
          nf <- eval(derivFxnCall[[1]][[2]], envir = fxnEnv)
          if(is(nf, "calcNodes_refClass")) {
              model <- eval(nf$Robject$model, envir = fxnEnv)
              if(!(exists('CobjectInterface', model) && !is(model$CobjectInterface, 'uninitializedField'))) { ## but, model should always be compiled because existence of calcNodes_refClass implies that.
                  cModel <- compileNimble(model)
              } else cModel <- eval(quote(model$CobjectInterface))
              model_calcNodes_case <- TRUE
          }
      }
  }

  if(!is.na(dropArgs)){
    removeArgs <- which(wrt == dropArgs)
    if(length(removeArgs) > 0)
      wrt <- wrt[-removeArgs]
  }
  if(model_calculate_case) {
      ans <- nimDerivs_model( order = order, wrt = wrt, derivFxnCall = derivFxnCall, fxnEnv = fxnEnv )
  } else {
      if(model_calcNodes_case) {
          ans <- nimDerivs_calcNodes( order = order, wrt = wrt, derivFxnCall = derivFxnCall, fxnEnv = fxnEnv, model = cModel, nodes = calcNodes)
      } else ans <- nimDerivs_nf( order = order, wrt = wrt, derivFxnCall = derivFxnCall, fxnEnv = fxnEnv )
  }
  ans
}

calcDerivs_hessian <- function(func, X, n, deriv_package = "pracma") {
    use_pracma <- deriv_package == "pracma"
    if(missing(n)) {
        check <- func(X)
        n <- length(check)
    }
    p <- length(X)
    outHess <- array(NA, dim = c(p, p, n))
    for(i in 1:n) {
        wrapped_func <- function(X) as.numeric(func(X))[i]
        if(use_pracma)
            outHess[, , i] <- pracma::hessian(wrapped_func, X)
        else
            outHess[, , i] <- numDeriv::hessian(wrapped_func, X)
    }
    outHess
}

calcDerivs_internal <- function(func, X, order, resultIndices ) {
    if(!require('numDeriv'))
        stop("The 'numDeriv' package must be installed to use derivatives in
         uncompiled nimbleFunctions.")

    hessianFlag <- 2 %in% order
    jacobianFlag <- 1 %in% order
    ## When called for a model$calculate, valueFlag will always be 0 here
    ## because value will be obtained later (if requested) after restoring
    ## model variables.
    valueFlag <- 0 %in% order
    outVal <- NULL

    if(hessianFlag) {
        ## determine output size
        returnType_scalar <- FALSE
        funcEnv <- environment(func)
        if (deparse(funcEnv$fxnCall[[1]]) == 'nimDerivs_nf') {
            fxnEnv <- funcEnv[['fxnEnv']]
            derivFxnCall <- funcEnv[['derivFxnCall']]
            genFunEnv <- environment(
                environment(
                    eval(derivFxnCall[[1]], envir = fxnEnv)
                )[['.generatorFunction']]
            )
            returnType <- argType2symbol(
                genFunEnv$methodList[[funcEnv$derivFxnName]]$returnType
            )
            returnType_scalar <- all(returnType$size == 1)
        } else {
            outVal <- func(X)
            returnType_scalar <- length(outVal) == 1
        }

        use_pracma <- requireNamespace("pracma")
        if(!use_pracma) {
            message("Install package pracma for more accurate numerical Hessians in uncompiled (R) execution (experimental).")
        }
        
        if (returnType_scalar) {
            ## numDeriv's hessian() only works for scalar-valued functions
            if(use_pracma) {
                hess <- pracma::hessian(func, X)
            } else {
                hess <- numDeriv::hessian(func, X)
            }
            outHess <- array(NA, dim = c(dim(hess)[1], dim(hess)[2], 1))
            outHess[,,1]  <- hess
            if(jacobianFlag) outGrad <- jacobian(func, X)
            if(valueFlag && is.null(outVal)) outVal <- func(X)
        } else {
            ## If hessians are requested for a non-scalar valued function,
            ## derivatives taken using numDeriv's genD() function.  After that,
            ## we extract the various derivative elements and arrange them properly.
            if(use_pracma) {
                value <- func(X)
                if(valueFlag && is.null(outVal)) outVal <- value
                if(jacobianFlag) outGrad <- jacobian(func, X)
                outHess <- calcDerivs_hessian(func, X, n = length(value), deriv_package = "pracma") 
            } else {
                derivList <- genD(func, X)
                if(valueFlag && is.null(outVal)) outVal <- derivList$f0
                if(jacobianFlag) outGrad <- derivList$D[,1:derivList$p, drop = FALSE]
                outHessVals <- derivList$D[,(derivList$p + 1):dim(derivList$D)[2],
                                           drop = FALSE]
                outHess <- array(NA, dim = c(derivList$p, derivList$p, length(derivList$f0)))
                singleDimMat <- matrix(NA, nrow = derivList$p, ncol = derivList$p)
                singleDimMatUpperTriDiag <- upper.tri(singleDimMat, diag = TRUE)
                singleDimMatLowerTriDiag <- lower.tri(singleDimMat)
                for(outDim in seq_along(derivList$f0)){
                    singleDimMat[singleDimMatUpperTriDiag] <- outHessVals[outDim,]
                    singleDimMat[singleDimMatLowerTriDiag] <- t(singleDimMat)[singleDimMatLowerTriDiag]
                    outHess[,,outDim] <- singleDimMat
                }
            }
        }
    } else if(jacobianFlag){
        ## If jacobians are requested, derivatives taken using numDeriv's jacobian()
        ## function.
        if(valueFlag) outVal <- func(X)
        outGrad <- jacobian(func, X)
    } else if(valueFlag)
        outVal <- func(X)
    
    outList <- ADNimbleList$new()
    if(!missing(resultIndices)) {
        if(jacobianFlag) outGrad <- outGrad[, resultIndices, drop=FALSE]
        if(hessianFlag) outHess <- outHess[resultIndices, resultIndices, , drop=FALSE]
    }
    if(valueFlag) outList$value <- outVal
    if(jacobianFlag) outList$jacobian <- outGrad
    if(hessianFlag) outList$hessian <- outHess
    return(outList)
}

nimDerivs_model <- function( nimFxn, order = c(0, 1, 2), wrt = NULL, derivFxnCall = NULL, fxnEnv = parent.frame()) {
    fxnCall <- match.call()
    ## This may be called directly (for testing) or from nimDerivs (typically).
    ## In the former case, we get derivFxnCall from the nimFxn argument.
    ## In the latter case, it is already match.call()ed and is passed here via derivFxnCall
    if(is.null(derivFxnCall))
        derivFxnCall <- fxnCall[['nimFxn']] ## either calculate(model, nodes) or model$calculate(nodes)
    else
        if(is.null(fxnCall[['order']]))
            fxnCall[['order']] <- order

    valueFlag <- 0 %in% order
    
    if(deparse(derivFxnCall[[1]]) == 'calculate'){
        ## calculate(model, nodes) format
        derivFxnCall <- match.call(nimble:::calculate, derivFxnCall)
        model <- eval(derivFxnCall[['model']], envir = fxnEnv)
        nodes <- eval(derivFxnCall[['nodes']], envir = fxnEnv)
        if(is.null(nodes))
            nodes <- model$getMaps('nodeNamesLHSall')
    } else {
        ## model$calculate(nodes) format, where nodes could be missing
        model <-  eval(derivFxnCall[[1]][[2]], envir = fxnEnv)
        if(length(derivFxnCall) < 2){ ## nodes is missing: default to all nodes
          nodes <- model$getMaps('nodeNamesLHSall')
        } else {
          nodes <- eval(derivFxnCall[[2]], envir = fxnEnv)
        }
    }
    ## We rely on fact that model$expandNodeNames will return column-major order.
    ## That makes it suitable for the canonical ordering of any blocks in wrt.

    ## We need to expandNodeNames (among other reasons) to determine
    ## any duplicates in wrt, such as "x" and "x[2]".
    
    ## It might also be possible to make use of values() and values()<- here.
    ## However, we instead are using a single model utility function (expandNodeNames)
    ## and basing every step on the result of that.  The hope is this will
    ## minimize risk of errors in case values() would somehow handle wrt in
    ## a different order than expandNodeNames.  (It shouldn't, but values() and values()<-
    ## have implementations that are hard to follow.)

    wrt_unique_names <- unique(.Call(nimble:::parseVar, wrt))
    if(!all(wrt_unique_names %in% model$getVarNames())){
        stop('Error:  the wrt argument to nimDerivs() contains variable names that are not
         in the model.')
    }
    
    ## scalar node elements, in order for results
    wrt_expanded <- model$expandNodeNames(wrt, returnScalarComponents = TRUE, unique = FALSE)
    ## unique scalar node elements
    wrt_expanded_unique <- unique(wrt_expanded)
    ## indices of original node elements in unique elements, for construction of results after getting derivatives wrt unique elements
    wrt_orig_indices <- match(wrt_expanded, wrt_expanded_unique)
    
    ## List of lines like "model$beta[2] <- x[i]", where "beta[2]" is an expanded wrt node and x[i] is a variable inside func.
    wrt_expanded_unique_assign_code <- mapply(
        function(node, I)
            parse(text = paste0("model$", node, " <- x[", I, "]"), keep.source = FALSE)[[1]],
        wrt_expanded_unique,
        seq_along(wrt_expanded_unique))
    ## We could consider doing the above step by substitute, but it would take some care to handle nodes with indices.
    
    length_wrt <- length(wrt_expanded_unique_assign_code)
    
    ## func could be made more efficient if we re-block wrt nodes that are in the same variable.
    ## makeVertexNamesFromIndexArray2 might be helpful for that purpose.
    func <- function(x) {
        do.call("{", wrt_expanded_unique_assign_code)
        ## equivalent to:
        ##        for(i in 1:length_wrt) {
        ##            eval(wrt_expanded_unique_assign_code[[i]])
        ##        }
        model$calculate(nodes)
    }
    
    ## make current values in single vector.
    ## This list will have code like `currentX[3] <- model$beta[2]`
    wrt_expanded_unique_get_code <- mapply(
        function(node, I)
            parse(text = paste0("currentX[",I,"] <- model$", node), keep.source = FALSE)[[1]],
        wrt_expanded_unique,
        seq_along(wrt_expanded_unique))
    
    currentX <- numeric(length_wrt)
    do.call("{", wrt_expanded_unique_get_code)
    
    ## Determine what nodes to restore.  We will do this by variable because in many cases
    ## it will be faster to copy and restore a whole variable in one line that to do elements
    ## with one line for each element.
    ##
    ## If value is not being returned, we restore wrt vars, nodes vars, and any logProb vars for nodes
    ## If value is being returned, we recalculate at the end, using func.  This restores any wrt nodes that
    ## may have been perturbed during finite-element derivatives as well as updating all nodes
    if(!valueFlag) {
        wrtVars <- unique(.Call(nimble:::parseVar, wrt))
        calcNodeVars <- unique(.Call(nimble:::parseVar, nodes))
                                        # determine any stochastic nodes and include their logProb nodes
        anyStoch_calcNodeVars <- unlist(lapply(calcNodeVars, function(x) model$getVarInfo(x)$anyStoch))
        logProbVars <- nimble:::makeLogProbName( calcNodeVars[anyStoch_calcNodeVars] )
        varsToRestore <- c(union(wrtVars, calcNodeVars), logProbVars)
        savedVars <- new.env()
    } else {
        varsToRestore <- character()
    }
    for(i in varsToRestore)
        savedVars[[i]] <- model[[i]]
    
    ## outGrad <- jacobian(func, currentX)
    ## We do not want calcDerivs_internal to calculate the value,
    ## because we want to restore variables first.
    if(valueFlag) order <- order[ order != 0 ]
    if(length(order) > 0) {
        ans <- calcDerivs_internal(func, currentX, order, wrt_orig_indices)
    } else ans <- NULL
    for(i in varsToRestore)
        model[[i]] <- savedVars[[i]]
    
    if(valueFlag) {
        ans$value <- func(currentX)
    }
    ans
}

nimDerivs_calcNodes <- function( nimFxn, order = c(0, 1, 2), wrt = NULL, derivFxnCall = NULL, fxnEnv = parent.frame(), model, nodes) {
    ## Handles nimDerivs(calcNodes$run()) case, used in AD testing. This function could be merged
    ## into nimDerivs_model as much of the code is duplicated but given nimDerivs for calculate
    ## is user facing and this is not, it's standalone code for now.
    
    fxnCall <- match.call()
    ## This may be called directly (for testing) or from nimDerivs (typically).
    ## In the former case, we get derivFxnCall from the nimFxn argument.
    ## In the latter case, it is already match.call()ed and is passed here via derivFxnCall
    if(is.null(derivFxnCall))
        derivFxnCall <- fxnCall[['nimFxn']] ## either calculate(model, nodes) or model$calculate(nodes)
    else
        if(is.null(fxnCall[['order']]))
            fxnCall[['order']] <- order

    valueFlag <- 0 %in% order
    
    ## We rely on fact that model$expandNodeNames will return column-major order.
    ## That makes it suitable for the canonical ordering of any blocks in wrt.

    ## We need to expandNodeNames (among other reasons) to determine
    ## any duplicates in wrt, such as "x" and "x[2]".
    
    ## It might also be possible to make use of values() and values()<- here.
    ## However, we instead are using a single model utility function (expandNodeNames)
    ## and basing every step on the result of that.  The hope is this will
    ## minimize risk of errors in case values() would somehow handle wrt in
    ## a different order than expandNodeNames.  (It shouldn't, but values() and values()<-
    ## have implementations that are hard to follow.)

    wrt_unique_names <- unique(.Call(nimble:::parseVar, wrt))
    if(!all(wrt_unique_names %in% model$getVarNames())){
        stop('Error:  the wrt argument to nimDerivs() contains variable names that are not
         in the model.')
    }
    
    ## scalar node elements, in order for results
    wrt_expanded <- model$expandNodeNames(wrt, returnScalarComponents = TRUE, unique = FALSE)
    ## unique scalar node elements
    wrt_expanded_unique <- unique(wrt_expanded)
    ## indices of original node elements in unique elements, for construction of results after getting derivatives wrt unique elements
    wrt_orig_indices <- match(wrt_expanded, wrt_expanded_unique)
    
    ## List of lines like "model$beta[2] <- x[i]", where "beta[2]" is an expanded wrt node and x[i] is a variable inside func.
    wrt_expanded_unique_assign_code <- mapply(
        function(node, I)
            parse(text = paste0("model$", node, " <- x[", I, "]"), keep.source = FALSE)[[1]],
        wrt_expanded_unique,
        seq_along(wrt_expanded_unique))
    ## We could consider doing the above step by substitute, but it would take some care to handle nodes with indices.
    
    length_wrt <- length(wrt_expanded_unique_assign_code)
    
    ## func could be made more efficient if we re-block wrt nodes that are in the same variable.
    ## makeVertexNamesFromIndexArray2 might be helpful for that purpose.
    func <- eval(substitute(function(x) {
        do.call("{", wrt_expanded_unique_assign_code)
        ## equivalent to:
        ##        for(i in 1:length_wrt) {
        ##            eval(wrt_expanded_unique_assign_code[[i]])
        ##        }
        CALCCALL
    }, list(CALCCALL = derivFxnCall)))
    
    ## make current values in single vector.
    ## This list will have code like `currentX[3] <- model$beta[2]`
    wrt_expanded_unique_get_code <- mapply(
        function(node, I)
            parse(text = paste0("currentX[",I,"] <- model$", node), keep.source = FALSE)[[1]],
        wrt_expanded_unique,
        seq_along(wrt_expanded_unique))
    
    currentX <- numeric(length_wrt)
    do.call("{", wrt_expanded_unique_get_code)
    
    ## Determine what nodes to restore.  We will do this by variable because in many cases
    ## it will be faster to copy and restore a whole variable in one line that to do elements
    ## with one line for each element.
    ##
    ## If value is not being returned, we restore wrt vars, nodes vars, and any logProb vars for nodes
    ## If value is being returned, we recalculate at the end, using func.  This restores any wrt nodes that
    ## may have been perturbed during finite-element derivatives as well as updating all nodes
    if(!valueFlag) {
        wrtVars <- unique(.Call(nimble:::parseVar, wrt))
        calcNodeVars <- unique(.Call(nimble:::parseVar, nodes))
                                        # determine any stochastic nodes and include their logProb nodes
        anyStoch_calcNodeVars <- unlist(lapply(calcNodeVars, function(x) model$getVarInfo(x)$anyStoch))
        logProbVars <- nimble:::makeLogProbName( calcNodeVars[anyStoch_calcNodeVars] )
        varsToRestore <- c(union(wrtVars, calcNodeVars), logProbVars)
        savedVars <- new.env()
    } else {
        varsToRestore <- character()
    }
    for(i in varsToRestore)
        savedVars[[i]] <- model[[i]]
    
    ## outGrad <- jacobian(func, currentX)
    ## We do not want calcDerivs_internal to calculate the value,
    ## because we want to restore variables first.
    if(valueFlag) order <- order[ order != 0 ]
    if(length(order) > 0) {
        ans <- calcDerivs_internal(func, currentX, order, wrt_orig_indices)
    } else ans <- NULL
    for(i in varsToRestore)
        model[[i]] <- savedVars[[i]]
    
    if(valueFlag) {
        ans$value <- func(currentX)
    }
    ans
}

nimDerivs_nf_charwrt <- function(derivFxn,
                                 derivFxnCall,
                                 order,
                                 wrt,
                                 fxnEnv) {
  fA <- formalArgs(derivFxn)

  ## convert 'x[2]' to quote(x[2]).
  wrt_code <- lapply(wrt,
                     function(x) parse(text = x, keep.source = FALSE)[[1]])
  
  ## convert quote(x[2]) to quote(x)
  wrt_names <- lapply(wrt_code,
                      function(x) if(is.name(x)) x else x[[2]])
  
  wrt_name_strings <- as.character(wrt_names)
  
  ## Get unique names and track indices from wrt to unique names
  wrt_unique_names <- unique(wrt_names)
  
  if(!all(wrt_unique_names %in% fA)){
    stop('Error:  the wrt argument to nimDerivs() contains names that are not
         arguments to the nimFxn argument.')
  }
  
  wrt_names_orig_indices <- match(wrt_names, wrt_unique_names)
  wrt_unique_name_strings <- as.character(wrt_unique_names)

  ## get the user-supplied arguments to the nFunction
  derivFxnCall_args <- as.list(derivFxnCall)[-1]

  ## Get the value of the user-supplied args
  ## These will be modified in func.
  ## This is ordered by fA, the arguments of the target function
  fxnArgs <- lapply(derivFxnCall_args, function(x) eval(x, envir = fxnEnv))
  names(fxnArgs) <- as.character(fA)
  
  ## Get length or dimensions
  ## ordered by fxnArgs
  unique_dims <- sapply(fxnArgs,
                 function(x) 
                     nimble:::dimOrLength(x))

  ## Get product of dimensions, which we call size
  ## ordered by fxnArgs
  unique_sizes <- sapply(unique_dims, prod)

  ## 
  scalar_index_objects <- lapply(unique_dims,
                                 function(x)
                                     array(1:prod(x), dim = x))

  eval_env <- new.env()
  ## iterate over names in any wrt.

  current_x_index <- 1
  fxnArgs_assign_code <- list()
  get_init_values_code <- list()
  result_x_indices_by_wrt <- list()
  for(i in seq_along(wrt_unique_names)) {
      ## Below, x represents the input argument to func
      this_unique_wrt_string <- wrt_unique_name_strings[i]
      ##
      dims <- nimble:::dimOrLength(fxnArgs[[this_unique_wrt_string]])
      ##
      assign(this_unique_wrt_string, array(1:prod(dims), dim = dims), envir = eval_env)
      
      ## which wrt arguments use this wrt variable
      i_wrt_orig <- which(this_unique_wrt_string == wrt_name_strings)     

      flat_indices <- lapply(wrt_code[i_wrt_orig], ## quote(x[2]) and so on for any wrt using x
                             function(wrt_code_) {
                                 eval(wrt_code_, envir = eval_env)
                             })
      unique_flat_indices <- unique(unlist(flat_indices))
      ## I is this_unique_wrt_string
      ## J is an element of unique_flat_indices
      ## K is a running element of x

      x_indices <- current_x_index - 1 + 1:length(unique_flat_indices)
      
      fxnArgs_assign_code[[i]] <- mapply(
          function(jval, kval)
              substitute(fxnArgs[[ I ]][J] <<- x[K], list(I = this_unique_wrt_string,
                                                         J = jval,
                                                         K = kval)),
          unique_flat_indices,
          x_indices)
      get_init_values_code[[i]] <- mapply(
          function(jval, kval)
              substitute(currentX[K] <- fxnArgs[[ I ]][J], list(I = this_unique_wrt_string,
                                                             J = jval,
                                                             K = kval)),
          unique_flat_indices,
          x_indices)

      result_x_indices <- lapply(flat_indices,
                                      function(fi) x_indices[match(fi, unique_flat_indices)]
                                      )
      result_x_indices_by_wrt[i_wrt_orig] <- result_x_indices
      
     current_x_index <- current_x_index + length(unique_flat_indices)
  }
  fxnArgs_assign_code <- unlist(fxnArgs_assign_code, recursive = FALSE)
  get_init_values_code <- unlist(get_init_values_code, recursive = FALSE)
  result_x_indices_all <- unlist(result_x_indices_by_wrt)

  length_x <- current_x_index - 1
  currentX <- numeric(length_x)
  do.call("{", get_init_values_code)
  ## equivalent to:
  ##  for(i in 1:length_x) {
  ##      eval(get_init_values_code[[i]])
  ##  }

  derivFxnName <- as.character(derivFxnCall[[1]])
  func <- function(x) {
      do.call("{", fxnArgs_assign_code)
      ## for(i in 1:length_x) {
      ##     ## Each line is like fxnArgs[[ 2 ]][23] <- x[5]
      ##     eval(fxnArgs_assign_code[[i]])
      ## }
      c(do.call(derivFxnName, fxnArgs, envir = fxnEnv)) #c() unrolls any answer to a vector
  }

  ans <- calcDerivs_internal(func, currentX, order, result_x_indices_all)
##  jacobian(func, currentX)
  ans 
}


nimDerivs_nf_numericwrt <- function(derivFxn,
                                    derivFxnCall,
                                    order,
                                    wrt,
                                    fxnEnv) {
  fA <- formalArgs(derivFxn)
  derivFxnCall_args <- as.list(derivFxnCall)[-1]

  fxnArgs <- lapply(derivFxnCall_args, function(x) eval(x, envir = fxnEnv))
  names(fxnArgs) <- as.character(fA)

  argSizes <- unlist(lapply(fxnArgs, dimOrLength))

  flatStartIndices <- cumsum(c(1, argSizes[-length(argSizes)]))
  names(flatStartIndices) <- names(fxnArgs)
  wrtIndexToArgName <- character()
  for(i in seq_along(fxnArgs)) {
    wrtIndexToArgName <- c(wrtIndexToArgName, rep(names(fxnArgs)[i], dimOrLength(fxnArgs[[i]])))
  }

  ## We need lines like
  ## fxnArgs[[ argName ]][j] <- x[k]
  ## j is the flat index of the variable
  ## k is the index of the concatenated x
  wrtUnique <- unique(wrt)
  fxnArgs_assign_code <- list()
  get_init_values_code <- list()
  length_result <- length(wrt)
  result_x_indices <- integer(length_result)
  for(iWrt in seq_along(wrtUnique)) {
    wrtIndex <- wrtUnique[iWrt]
    i_wrt_orig <- which(wrtIndex == wrt) ## indices of original, allowing repeats
    argName <- wrtIndexToArgName[wrtIndex]
    j <- wrtIndex - flatStartIndices[argName] + 1
    fxnArgs_assign_code[[ length(fxnArgs_assign_code) + 1]] <-
      substitute(fxnArgs[[ ARGNAME ]][J] <- x[K],
                 list(ARGNAME = argName,
                      J = as.integer(j),
                      K = iWrt))
    get_init_values_code[[ length(fxnArgs_assign_code) + 1]] <-
      substitute(currentX[K] <- fxnArgs[[ ARGNAME ]][J],
                 list(ARGNAME = argName,
                      J = as.integer(j),
                      K = iWrt))
    result_x_indices[i_wrt_orig] <- iWrt
    ## still need to add result_x_indices to deal with redundant or disordered requests
  }
  length_x <- length(wrtUnique)
  currentX <- numeric(length_x)
  do.call("{", get_init_values_code)
  
  derivFxnName <- as.character(derivFxnCall[[1]])
  func <- function(x) {
    do.call("{", fxnArgs_assign_code)
    c(do.call(derivFxnName, fxnArgs, envir = fxnEnv)) #c() unrolls any answer to a vector
  }
  ans <- calcDerivs_internal(func, currentX, order, result_x_indices)
  ans 
}

nimDerivs_nf <- function(nimFxn = NA, order = nimC(0,1,2),
                         wrt = NULL, derivFxnCall = NULL, fxnEnv = parent.frame()){
  fxnCall <- match.call()
  ## This may be called directly (for testing) or from nimDerivs (typically).
  ## In the former case, we get derivFxnCall from the nimFxn argument.
  ## In the latter case, it is already match.call()ed and is passed here via derivFxnCall
  if(is.null(derivFxnCall))
    derivFxnCall <- fxnCall[['nimFxn']] ## either calculate(model, nodes) or model$calculate(nodes)
  else
    if(is.null(fxnCall[['order']]))
      fxnCall[['order']] <- order ## I don't think this is meaningful any more. fxnCall is not used further

  ## get the actual function
  derivFxn <- eval(derivFxnCall[[1]], envir = fxnEnv)

  ## standardize the derivFxnCall arguments
  derivFxnCall <- match.call(derivFxn, derivFxnCall)

  fA <- formalArgs(derivFxn)
  if(is.null(wrt)) {
    wrt <- fA
  }
  if(is.character(wrt)) {
    ans <- nimDerivs_nf_charwrt(derivFxn,
                                derivFxnCall,
                                order,
                                wrt,
                                fxnEnv)
  } else if(is.numeric(wrt)) {
    ans <- nimDerivs_nf_numericwrt(derivFxn,
                                   derivFxnCall,
                                   order,
                                   wrt,
                                   fxnEnv)
  } else {
    stop("Invalid wrt argument in nimDerivs_nf")
  }
  ans
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
    if(length(x) < 3) stop(paste0('Problem with derivative-enabled nimbleFunction method ', fxnName,
                                 ': If the mode of enableDerivs is static=TRUE, sizes of arguments ',
                                 'must be explicitly specified (e.g. x = double(1, 4)).  If the mode of ',
                                 'enableDerivs is static=FALSE, then wrt cannot use argument ',
                                 'names.'))
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
                                    fxnArgsDimSizes[[wrtArgNames[i]]][1] +
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

