###		These functions are used for calculate/sim/getLP for the nodeFunctionVectors
###		Can either enter model, nodes or model_nodes

#' NIMBLE language function to break tracking of derivatives
#'
#' This function is used in a method of a nimbleFunction that has derivatives enabled.  It returns its value but breaks tracking of derivatives.
#'
#' @param x scalar value
#'
#' @details
#' This funcion only works with scalars.
#'
#' @export
ADbreak <- function(x) x

#' Power function for integer-valued exponent
#'
#' pow function with exponent required to be integer
#'
#' @param a Base
#' @param b Exponent
#'
#' @return a^b
#'
#' @details This is required in nimble models and nimbleFunctions if derivatives will be tracked
#' but tracked only with respect to a, such that b might be any (positive, 0, or negative) integer.
#' This contrasts with pow(a, b) (equivalent to a^b), which requires b > 0 if derivatives will be
#' tracked, even if they will only be requested with respect to a.
#' 
#' @export
pow_int <- function(a, b) a^round(b) # b should be integer.  rounding it makes finite-element derivatives 0, as needed.

#' NIMBLE language functions for R-like vector construction
#'
#' The functions \code{c}, \code{rep}, \code{seq}, \code{which}, \code{diag}, \code{length}, \code{seq_along}, \code{is.na}, \code{is.nan}, \code{any}, and \code{all} can be used in nimbleFunctions and compiled using \code{compileNimble}.
#'
#' @name nimble-R-functions
#'
#' @param ... values to be concatenated.
#' @param x vector of values to be replicated (\code{rep}), or logical array or vector (\code{which}), or object whose length is wanted (\code{length}), or input value (\code{diag}), or vector of values to be tested/checked (\code{is.na}, \code{is.nan}, \code{any}, \code{all}).
#' @param from starting value of sequence.
#' @param to end value of sequence.
#' @param by increment of the sequence.
#' @param length.out desired length of the sequence.
#'
#' @aliases nimC nimRep nimSeq c rep seq which diag length seq_along is.na is.nan any all
#'
#' @details
#' For \code{c}, \code{rep}, \code{seq}, these functions are NIMBLE's version of similar R functions, e.g., \code{nimRep} for \code{rep}.   In a \code{nimbleFunction}, either the R name (e.g., \code{rep}) or the NIMBLE name (e.g., \code{nimRep}) can be used.  If the R name is used, it will be converted to the NIMBLE name. For \code{which}, \code{length}, \code{diag}, \code{seq_along}, \code{is.na}, \code{is.nan}, \code{any}, \code{all} simply use the standard name without \code{"nim"}. These functions largely mimic (see exceptions below) the behavior of their R counterparts, but they can be compiled in a \code{nimbleFunction} using \code{compileNimble}.
#' 
#' \code{nimC} is NIMBLE's version of \code{c} and behaves identically.
#'
#' \code{nimRep} is NIMBLE's version of \code{rep}.  It should behave identically to \code{rep}.  There are no NIMBLE versions of \code{rep.int} or \code{rep_len}.
#'
#' \code{nimSeq} is NIMBLE's version of \code{seq}.  It behaves like \code{seq} with support for \code{from}, \code{to}, \code{by} and \code{length.out} arguments.  The \code{along.with} argument is not supported.  There are no NIMBLE versions of \code{seq.int}, \code{seq_along} or \code{seq_len}, with the exception that \code{seq_along} can take a nimbleFunctionList as an argument to provide the index range of a for-loop (\href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} Ch. 13). 
#'
#' \code{which} behaves like the R version but without support for \code{arr.ind} or \code{useNames} arguments.
#'
#' \code{diag} behaves like the R version but without support for the \code{nrow} and \code{ncol} arguments.
#'
#' \code{length} behaves like the R version.
#' 
#' \code{seq_along} behaves like the R version.
#'
#' \code{is.na} behaves like the R version but does not correctly handle \code{NA} values from R that are type 'logical', so convert these using \code{as.numeric()} before passing from R to NIMBLE.
#' 
#' \code{is.nan} behaves like the R version, but treats \code{NA} of type 'double' as being \code{NaN} and \code{NA} of type 'logical' as not being \code{NaN}. 
#' 
#' \code{any} behaves like the R version but takes only one argument and treats NAs as \code{FALSE}.
#'
#' \code{all} behaves like the R version but takes only one argument and treats NAs as \code{FALSE}.
#'
NULL

#' @rdname nimble-R-functions
#' @export
nimC <- function(...) {
    c(...)
}

#' @rdname nimble-R-functions
#' @export
nimRep <- function(x, ...) {
    rep(x, ...)
}

#' @rdname nimble-R-functions
#' @export
nimSeq <- function(from, to, by, length.out) { ## this creates default arguments filled in at keyword matching that are then useful in cppOutput step to determine if C++ should use nimSeqBy or nimSeqLen
    haveBy <- !missing(by)
    haveLen <- !missing(length.out)
    if(haveBy)
        if(haveLen)
            seq(from, to, by, length.out)
        else
            seq(from, to, by)
    else
        if(haveLen)
            seq(from, to, length.out = length.out)
        else
            seq(from, to)
}

#' Explicitly declare objects created in setup code to be preserved and compiled as member data
#'
#' Normally a nimbleFunction determines what objects from setup code need to be preserved for run code or other member functions.  \code{setupOutputs} allows explicit declaration for cases when an object created in setup code is not used in member functions.
#' 
#' @name setupOutputs
#'
#' @param ... An arbitrary set of names
#'
#' @details
#' Normally any object created in \code{setup} whose name appears in \code{run} or another member function is included in the saved results of setup code.  When the nimbleFunction is compiled, such objects will become member data of the resulting C++ class.  If it is desired to force an object to become member data even if it does not appear in a member function, declare it using \code{setupOutputs}.  E.g., \code{setupOutputs(a, b)} declares that \code{a} and \code{b} should be preserved.
#'
#' The \code{setupOutputs} line will be removed from the setup code.  It is really a marker during nimbleFunction creation of what should be preserved.
#' 
NULL

#' Halt execution of a nimbleFunction function method.  Part of the NIMBLE language
#'
#'
#' @param msg Character object to be output as an error message
#'
#' @aliases stop
#'
#' @author Perry de Valpine
#' @export
#' @details
#' The NIMBLE \code{stop} is similar to the native R \code{stop}, but it takes only one argument, the error message to be output.  During uncompiled NIMBLE execution, \code{nimStop} simply calls R's stop funtion. During compiled execution it calls the error function from the R headers.  \code{stop} is an alias for \code{nimStop} in the NIMBLE language  
nimStop <- function(msg) stop(msg, call. = FALSE)
# we use call.=FALSE because otherwise the error msg indicates the
# error itself occurs in nimStop() and not in the calling frame

#' Time execution of NIMBLE code
#'
#' @param code code to be timed
#'
#' @author NIMBLE Development Team
#' @details
#' Function for use in nimbleFunction run code; when nimbleFunctions are run in R, this simply wraps \code{system.time}.
#' @export
run.time <- function(code) {
    as.numeric(system.time(code)[3])
}

#' Check for interrupt (e.g. Ctrl-C) during nimbleFunction execution. Part of the NIMBLE language.
#'
#' @author Perry de Valpine
#' @export
#' @details
#' During execution of nimbleFunctions that take a long time, it is nice to occassionally check if the user has entered an interrupt and bail out of execution if so.  This function does that.  During uncompiled nimbleFunction execution, it does nothing.  During compiled execution, it calls R_checkUserInterrupt() of the R headers. 
checkInterrupt <- function() {}

#' Turn a numeric vector into a single-row or single-column matrix
#'
#' Turns a numeric vector into a matrix that has 1 row or 1 column.  Part of NIMBLE language.
#'
#' @aliases asCol
#'
#' @param x       Numeric to be turned into a single row or column matrix
#'
#' @author Perry de Valpine
#' @export
#' @details
#' In the NIMBLE language, some automatic determination of how to turn vectors into single-row or single-column matrices is done.
#' For example, in \code{A \%*\% x}, where A is a matrix and x a vector, x will be turned into a single-column matrix unless
#' it is known at compile time that A is a single column, in which case x will be turned into a single-row matrix.
#' However, if it is desired that x be turned into a single row but A cannot be determined at compile time to be a single column,
#' then one can use \code{A \%*\% asRow(x)} to force this conversion.
asRow <- function(x) {
    matrix(x, nrow = 1)
}

## Aliased in asRow
#' @rdname asRow
#' @export
asCol <- function(x) {
    matrix(x, ncol = 1)
}

#' Make an object of information about a model-parameter pairing for getParam.  Used internally
#'
#' Creates a simple getParam_info object, which has a list with a paramID and a type
#'
#' @param model A model such as returned by \code{\link{nimbleModel}}.
#'
#' @param nodes A character string naming one one or more stochastic nodes, such as "mu", "c('mu', 'beta[2]')", or "eta[1:3, 2]".  getParam only works for one node at a time, but if it is indexed (nodes[i]), then makeParamInfo sets up the information for the entire vector nodes.  The processing pathway is used by the NIMBLE compiler.
#'
#' @param param A character string naming a parameter of the distribution followed by node, such as "mean", "rate", "lambda", or whatever parameter names are relevant for the distribution of the node.
#'
#' @param vector A logical indicating whether nodes should definitely be treated as a vector in compiled code, even if it has length = 1.  For type consistency, the compiler needs this option.  If nodes has length > 1, this argument is ignored.
#' 
#' @export
#' @details This is used internally by \code{\link{getParam}}.  It is not intended for direct use by a user or even a nimbleFunction programmer.
makeParamInfo <- function(model, nodes, param, vector = FALSE) {
    ## updating to allow nodes to be a vector. getParam only works for a scalar but in a case like nodes[i] the param info is set up for the entire vector.

    ## this allows for(i in seq_along(nodes)) a <- a + model$getParam(nodes[i], 'mean') through compilation even if some instances have nodes empty and so won't be called.
  if(length(nodes) == 0) return(structure(c(list(paramID = integer()), type = NA, nDim = NA), class = 'getParam_info'))
  
  if(length(param) != 1) stop(paste0(paste0('Problem with param(s) ', paste0(param, collapse = ','), ' while setting up getParam for node ', nodes,
                                            '\nOnly one parameter is allowed.')))
  
  nodeIDs <- model$expandNodeNames(nodes, returnType = 'ids')
  nodeDeclIDs <- model$modelDef$maps$graphID_2_declID[nodeIDs]
  declID2nodeIDs <- split(nodeIDs, nodeDeclIDs)
  numDeclIDs <- length(declID2nodeIDs)
  paramIDs <- integer(numDeclIDs)
  types <- character(numDeclIDs)
  nDims <- integer(numDeclIDs)
  for(i in seq_along(declID2nodeIDs)) {
    nodeIDsFromOneDecl <- declID2nodeIDs[[i]]
    firstNodeName <- model$modelDef$maps$graphID_2_nodeName[nodeIDsFromOneDecl[1]]
    dist <- model$getDistribution(firstNodeName)
    distInfo <- getDistributionInfo(dist)
    paramIDs[i] <- distInfo$paramIDs[param]
    if(is.na(paramIDs[i]))
      stop(paste0("getParam: parameter '", param, "' not found in distribution ",
                  dist, "."))
    types[i] <- distInfo$types[[param]]$type
    nDims[i] <- distInfo$types[[param]]$nDim
  }
  if(length(unique(types)) != 1 || length(unique(nDims)) != 1)
    stop('cannot have an indexed vector of nodes used in getParam if they have different types or dimensions for the same parameter.')
    
  ## on C++ side, we always work with double
  if(types[1] %in% c('integer', 'logical')) types[1] <- 'double'

  if(length(paramIDs) == 1) { ## We could shortcut on this case earlier
    if(vector)
      paramIDvec <- c(-1L, paramIDs[1]) # paramIDs is length 1 anyway
    else
      paramIDvec <- paramIDs[1]
  } else {
    if(length(unique(paramIDs)) == 1) {
      # They are all the same.
      # We encode this in a sparse way
      paramIDvec <- c(-1L, paramIDs[1])
    } else {
      ## Otherwise, create a full vector of paramIDs
      paramIDvec <- rep(1L, length(nodeIDs))
      sourceIDs <- split(seq_along(nodeIDs), nodeDeclIDs)
      for(i in seq_along(declID2nodeIDs)) { #unsplit() would be another approach to this step.
        paramIDvec[sourceIDs[[i]] ] <- paramIDs[i]
      }
    }
  }
    
  ans <- c(list(paramID = paramIDvec), type = types[1], nDim = nDims[1])
  class(ans) <- 'getParam_info'
  ans
}

defaultParamInfo <- function() {
    ans <- list(paramID = integer(), type = 'double', nDim = 0)
    class(ans) <- 'getParam_info'
    ans
}


#' Get value of a parameter of a stochastic node in a model
#'
#' Get the value of a parameter for any single stochastic node in a model.
#'
#' @param model A NIMBLE model object
#'
#' @param node  The name of a stochastic node in the model
#'
#' @param param The name of a parameter for the node
#'
#' @param nodeFunctionIndex For internal NIMBLE use only
#'
#' @export
#' @details
#' Standard usage is as a method of a model, in the form \code{model$getParam(node, param)}, but the usage as a simple function with the model as the first argument as above is also allowed.
#'
#' For example, suppose node 'x[1:5]' follows a multivariate
#' normal distribution (dmnorm) in a model declared by BUGS code.
#' model$getParam('x[1:5]', 'mean') would return the current value of
#' the mean parameter (which may be determined from other nodes).  The
#' parameter requested does not have to be part of the
#' parameterization used to declare the node.  Rather, it can be any
#' parameter known to the distribution.  For example, one can request
#' the scale or rate parameter of a gamma distribution, regardless of
#' which one was used to declare the node.
getParam <- function(model, node, param, nodeFunctionIndex) {
    if(missing(param)) { ## already converted by keyword conversion
        stop('Missing either the node or the parameter. Please specify both the node and the parameter of interest. This error may also be triggered if this case of getParam (after keyword replacement) has not been updated for R execution with newNodeFunction system')
        ## nodeFunctionIndex would only be used here, when we make this part work
        nodeFunction <- model
        paramInfo <- node
    } else {
        ## not already converted; this is regular execution
        if(length(node) != 1) stop(paste0("getParam only works for one node at a time, but ", length(node), " were provided."))
        tmp <- model$expandNodeNames(node)
        if(length(tmp) != 1) stop(paste0("getParam only works for one node at a time, but ", node, " includes multiple nodes."))
        ## makeParamInfo, called by nodeFunctionVector, will check on length of param
        ## nodeFunctionIndex should never be used.
        nfv <- nodeFunctionVector(model, node)
        indexingInfo <- nfv$indexingInfo
        declID <- indexingInfo$declIDs[1] ## should only be one
        nodeFunction <- model$nodeFunctions[[ declID ]] 
        paramInfo <- makeParamInfo(model, node, param)
        if(is.na(paramInfo$type)) stop(paste('getParam called with empty or invalid node:', as.character(node)))
    }
    paramID <- paramInfo$paramID
    if(paramID[1]==-1)
        paramID <- paramID[2]
    nDim <- paramInfo$nDim
    type <- paramInfo$type
    unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ indexingInfo$unrolledIndicesMatrixRows[1], ]
    if(nDim > 2) {
        warning("Note that getParam is not available for parameters of dimension greater than two, but other calculations with models are unaffected")
        return(NULL)
    }
    funName <- paste0('getParam_',nDim,'D_',type)

    useCompiledNonNestedInterface <- inherits(model, 'CmodelBaseClass') & !getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')

    if(useCompiledNonNestedInterface) {
        ans <- model$nodeFunctions[[ declID ]][[1]]$callMemberFunction(model$nodeFunctions[[ declID ]][[2]], funName, paramID, unrolledIndicesMatrixRow)
    } else 
        ans <- eval(substitute(nodeFunction$FUNNAME(paramID, unrolledIndicesMatrixRow), list(FUNNAME = as.name(funName))))
    return(ans)
}

#' Make an object of information about a model-bound pairing for getBound.  Used internally
#'
#' Creates a simple getBound_info object, which has a list with a boundID and a type.
#' Unlike \code{makeParamInfo} this is more bare-bones, but keeping it for parallelism
#' with \code{getParam}.
#'
#' @param model A model such as returned by \code{\link{nimbleModel}}.
#'
#' @param nodes A character string naming a stochastic nodes, such as \code{'mu'}.
#'
#' @param bound A character string naming a bound of the distribution, either \code{'lower'} or \code{'upper'}.
#'
#' @export
#' @details This is used internally by \code{\link{getBound}}.  It is not intended for direct use by a user or even a nimbleFunction programmer.
makeBoundInfo <- function(model, nodes, bound) {
    # note: would be good to work through this and makeParamInfo to ensure that both give good error messages and handle situation where user tries to pass in multiple nodes or multiple params/bounds
    # note: also, it should be fine to combine symbolTable and sizeProcessing stuff for param and bound to avoid repetitive code
    if(length(nodes) > 1) stop("'getBound' only set up to work with a single node'") # I think this is consistent with current getParam behavior
    distInfo <- getDistributionList(model$getDistribution(nodes))
    boundIDvec <- c('lower'=1,'upper'=2)[bound]
    typeVec <- unlist(lapply(distInfo, function(x) x$types[['value']]$type)) 
    ## on C++ side, we always work with double
    if(typeVec[1] %in% c('integer', 'logical')) typeVec[1] <- 'double'
    
    # at moment have bound be scalar regardless of dimension of parameter, as assumed to be same value for all elements of a multivariate distribution; we are not fully handling truncation/bounds for multivariate nodes anyway:
    nDimVec <- 0  # unlist(lapply(distInfo, function(x) x$types[['value']]$nDim))
    ans <- c(list(boundID = boundIDvec, type = typeVec, nDim = nDimVec))
    class(ans) <- 'getBound_info'
    ans
}

#' Get value of bound of a stochastic node in a model
#'
#' Get the value of the lower or upper bound for a single stochastic node in a model.
#'
#' @param model A NIMBLE model object
#'
#' @param node  The name of a stochastic node in the model
#'
#' @param bound Either \code{'lower'} or \code{'upper'} indicating the desired bound for the node
#'
#' @param nodeFunctionIndex For internal NIMBLE use only
#'
#' @export
#' @details
#' Standard usage is as a method of a model, in the form \code{model$getBound(node, bound)}, but the usage as a simple function with the model as the first argument as above is also allowed.
#'
#' For nodes that do not involve truncation of the distribution
#' this will return the lower or upper bound of the distribution, which
#' may be a constant or for a limited number of distributions a parameter
#' or functional of a parameter (at the moment in NIMBLE, the only case
#' where a bound is a parameter is for the uniform distribution. For nodes
#' that are truncated, this will return the desired bound, which may be
#' a functional of other quantities in the model or may be a constant. 
getBound <- function(model, node, bound, nodeFunctionIndex) {
    if(!bound %in% c('lower', 'upper'))
        stop("getBound: 'bound' must be either 'lower' or 'upper'")
    nfv <- nodeFunctionVector(model, node)
    indexingInfo <- nfv$indexingInfo
    declID <- indexingInfo$declIDs[1] ## should only be one
    nodeFunction <- model$nodeFunctions[[ declID ]] 
    boundInfo <- makeBoundInfo(model, node, bound)
    boundID <- boundInfo$boundID
    nDim <- boundInfo$nDim
    type <- boundInfo$type
    unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ indexingInfo$unrolledIndicesMatrixRows[1], ]
    funName <- paste0('getBound_',nDim,'D_',type)

    useCompiledNonNestedInterface <- inherits(model, 'CmodelBaseClass') & !getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')

    if(useCompiledNonNestedInterface) {
        ans <- model$nodeFunctions[[ declID ]][[1]]$callMemberFunction(model$nodeFunctions[[ declID ]][[2]], funName, boundID, unrolledIndicesMatrixRow)
    } else 
        ans <- eval(substitute(nodeFunction$FUNNAME(boundID, unrolledIndicesMatrixRow), list(FUNNAME = as.name(funName))))
    return(ans)
}


#' @export
nimSwitch <- function(paramID, IDoptions, ...) {
    dotsList <- eval(substitute(alist(...)))
    iUse <- which(IDoptions == paramID)
    eval(dotsList[[iUse]], envir = parent.frame())
    invisible(NULL)
}

rCalcNodes <- function(model, nfv){ ##nodeFunctionVector
    l_Prob = 0
    model <- nfv$model
    useCompiledNonNestedInterface <- inherits(model, 'CmodelBaseClass') & !getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')
    indexingInfo <- nfv$indexingInfo
    declIDs <- indexingInfo$declIDs
    numNodes <- length(declIDs)
    if(numNodes < 1) return(l_Prob)
    unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
    for(i in 1:numNodes) {
        declID <- declIDs[i]
        unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ unrolledIndicesMatrixRows[i], ]
        if(useCompiledNonNestedInterface) {
            l_Prob = l_Prob + model$nodeFunctions[[ declID ]][[1]]$callMemberFunction(model$nodeFunctions[[ declID ]][[2]], 'calculate', unrolledIndicesMatrixRow)
        } else
            l_Prob = l_Prob + model$nodeFunctions[[ declID ]]$calculate(unrolledIndicesMatrixRow) ## must use nodeFunctions to have declID ordering
    }
    return(l_Prob)
}

getNodeFunctionIndexedInfo <- function(indexedNodeInfo, iCol) indexedNodeInfo[iCol]

rCalcDiffNodes <- function(model, nfv){
    l_Prob <- 0
    model <- nfv$model
    useCompiledNonNestedInterface <- inherits(model, 'CmodelBaseClass') & !getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')
    indexingInfo <- nfv$indexingInfo
    declIDs <- indexingInfo$declIDs
    numNodes <- length(declIDs)
    if(numNodes < 1) return(l_Prob)
    unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
    for(i in 1:numNodes) {
        declID <- declIDs[i]
        unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ unrolledIndicesMatrixRows[i], ]
        if(useCompiledNonNestedInterface) {
            l_Prob = l_Prob + model$nodeFunctions[[ declID ]][[1]]$callMemberFunction(model$nodeFunctions[[ declID ]][[2]], 'calculateDiff', unrolledIndicesMatrixRow)
        } else
            l_Prob = l_Prob + model$nodeFunctions[[ declID ]]$calculateDiff(unrolledIndicesMatrixRow) ## must use nodeFunctions to have declID ordering
    }
    return(l_Prob)
}

#' calculate, calculateDiff, simulate, or get the current log probabilities (densities) a set of nodes in a NIMBLE model
#'
#' calculate, calculateDiff, simulate, or get the current log probabilities (densities) of one or more nodes of a NIMBLE model and (for calculate and getLogProb) return the sum of their log probabilities (or densities).  Part of R and NIMBLE.
#' @name nodeFunctions
#' 
#' @param model        A NIMBLE model, either the compiled or uncompiled version
#' @param nodes        A character vector of node names, with index blocks allowed, such as 'x', 'y[2]', or 'z[1:3, 2:4]'
#' @param nodeFxnVector An optional vector of nodeFunctions on which to operate, in lieu of \code{model} and \code{nodes}
#' @param nodeFunctionIndex For internal NIMBLE use only
#' @param includeData  A logical argument specifying whether \code{data} nodes should be simulated into (only relevant for \code{\link{simulate}})
#' @author NIMBLE development team
#' @details
#' Standard usage is as a method of a model, in the form \code{model$calculate(nodes)}, but the usage as a simple function with the model as the first argument as above is also allowed.
#' 
#' These functions expands the nodes and then process them in the model in the order provided.  Expanding nodes means turning 'y[1:2]' into c('y[1]','y[2]') if y is a vector of scalar nodes.
#' Calculation is defined for a stochastic node as executing the log probability (density) calculation and for a deterministic node as calculating whatever function was provided on the right-hand side of the model declaration.
#'
#' Difference calculation (calculateDiff) executes the operation(s) on the model as calculate, but it returns the sum of the difference between the new log probabilities and the previous ones.
#' 
#' Simulation is defined for a stochastic node as drawing a random value from its distribution, and for deterministic node as equivalent to calculate.
#'
#' getLogProb collects and returns the sum of the log probabilities of nodes, using the log probability values currently stored in the model (as generated from the most recent call to calculate on each node)
#'
#' These functions can be used from R or in NIMBLE run-time functions that will be compiled.  When executed in R (including when an uncompiled nimbleFunction is executed), they can be slow because the nodes are expanded each time.  When compiled in NIMBLE, the nodes are expanded only once during compilation, so execution will be much faster.
#'
#' It is common to want the nodes to be provided in topologically sorted order, so that they will be calculated or simulated following the order of the model graph.  Functions such as model$getDependencies(nodes, ...) return nodes in topologically sorted order.  They can be directly sorted by model$topologicallySortNodes(nodes), but if so it is a good idea to expand names first by model$topologicallySortNodes(model$expandNodeNames(nodes))
#'
#' @return calculate and getLogProb return the sum of the log probabilities (densities) of the calculated nodes, with a contribution of 0 from any deterministic nodes
#'
#' @return calculateDiff returns the sum of the difference between the new and old log probabilities (densities) of the calculated nodes, with a contribution of 0 from any deterministic nodes.
#'
#' simulate returns NULL.
#' 
#' @rdname nodeFunctions
#' @export
calculate <- function(model, nodes, nodeFxnVector, nodeFunctionIndex)	
{
    if(!missing(nodeFxnVector)){
        return(rCalcNodes(model, nodeFxnVector))
    }
    if(inherits(model, 'modelBaseClass') ){
        if(missing(nodes) ) 
            nodes <- model$getMaps('nodeNamesLHSall')
        nfv <- nodeFunctionVector(model, nodes)
        return(rCalcNodes(model, nfv))
    }	
}

#' @rdname nodeFunctions
#' @export
calculateDiff <- function(model, nodes, nodeFxnVector, nodeFunctionIndex)		
{
    if(!missing(nodeFxnVector)){
        return(rCalcDiffNodes(model, nodeFxnVector))
    }
    if(inherits(model, 'modelBaseClass') ){
        if(missing(nodes) ) 
            nodes <- model$getMaps('nodeNamesLHSall')
        nfv <- nodeFunctionVector(model, nodes)
        return(rCalcDiffNodes(model, nfv))
    }	
}

rGetLogProbsNodes <- function(model, nfv){
    l_Prob = 0
    model <- nfv$model
    useCompiledNonNestedInterface <- inherits(model, 'CmodelBaseClass') & !getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')
    indexingInfo <- nfv$indexingInfo
    declIDs <- indexingInfo$declIDs
    numNodes <- length(declIDs)
    if(numNodes < 1) return(l_Prob)
    unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
    for(i in 1:numNodes) {
        declID <- declIDs[i]
        unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ unrolledIndicesMatrixRows[i], ]
        if(useCompiledNonNestedInterface) {
            l_Prob = l_Prob + model$nodeFunctions[[ declID ]][[1]]$callMemberFunction(model$nodeFunctions[[ declID ]][[2]], 'getLogProb', unrolledIndicesMatrixRow)
        } else
            l_Prob = l_Prob + model$nodeFunctions[[ declID ]]$getLogProb(unrolledIndicesMatrixRow) ## must use nodeFunctions to have declID ordering
    }
    return(l_Prob)
}

#' @rdname nodeFunctions
#' @export
getLogProb <- function(model, nodes, nodeFxnVector, nodeFunctionIndex)		
{
    if(!missing(nodeFxnVector)){
        return(rGetLogProbsNodes(model, nodeFxnVector))
    }
    if( inherits(model, "modelBaseClass") ){		
        if(missing(nodes) ) 
            nodes <- model$getMaps('nodeNamesLHSall')
        nfv <- nodeFunctionVector(model, nodes)
    	return(rGetLogProbsNodes(model, nfv))
    }        
}


rSimNodes <- function(model, nfv){
    model <- nfv$model
    useCompiledNonNestedInterface <- inherits(model, 'CmodelBaseClass') & !getNimbleOption('buildInterfacesForCompiledNestedNimbleFunctions')
    indexingInfo <- nfv$indexingInfo
    declIDs <- indexingInfo$declIDs
    numNodes <- length(declIDs)
    if(numNodes < 1) return()
    unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
    for(i in 1:numNodes) {
        declID <- declIDs[i]
        unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[ unrolledIndicesMatrixRows[i], ]
        if(useCompiledNonNestedInterface) {
            model$nodeFunctions[[ declID ]][[1]]$callMemberFunction(model$nodeFunctions[[ declID ]][[2]], 'simulate', unrolledIndicesMatrixRow)
            
        } else
            model$nodeFunctions[[ declID ]]$simulate(unrolledIndicesMatrixRow) ## must use nodeFunctions to have declID ordering
    }
}

#' @rdname nodeFunctions
#' @export
simulate <- function(model, nodes, includeData = FALSE, nodeFxnVector, nodeFunctionIndex)		
{
    if(!missing(nodeFxnVector)){
        rSimNodes(model, nodeFxnVector)
    }
    else if( inherits(model, "modelBaseClass") ) {
        if(missing(nodes) ) 
            nodes <- model$getMaps('nodeNamesLHSall')
        nfv <- nodeFunctionVector(model, nodes, excludeData = !includeData)
        rSimNodes(nfv$model, nfv)
    }
}


# Fill a vector with flattened values from a set of model nodes
#
# Take a vector of node names for a model and fill a vector by concatenating their values.  Works in R and NIMBLE.
#
# @param vals        the variable in the calling function to receive the values
# @param model       a NIMBLE model object, either compiled or uncompiled
# @param nodes       a vector of node names, allowing index blocks that will be expanded
#
# @author NIMBLE development team
# @export
# @details
# Calling \code{getValues(P, model, nodes)} will modify P in the calling function.  This is being deprecated by
# \code{P <- values(model, nodes)}, but the development syntax of \code{getValues(P, model, nodes)} is still supported.
#
# The result P will be the concatenation of values from the nodes requested.  When requested nodes are from matrices or arrays, the values will be flattened into a vector following column-wise order.
#
# The reverse of \code{getValues(P, model, nodes)} is \code{setValues(P, model, nodes)}, which is being deprecated by \code{values(model, nodes) <- P}.
#
# These functions work in R and in NIMBLE run-time code that can be compiled.
#
# @return NULL, but this function works by the side-effect of modifying P in the calling environment.
getValues <- function(vals, model, nodes, envir = parent.frame()) {
    valsExp = substitute(vals)
    access = modelVariableAccessorVector(model, nodes, logProb = FALSE)
    output = getValuesAccess(access)
    if(is.name(valsExp))
        assign(as.character(valsExp), output, envir = parent.frame(n = 1) ) ## again this assumes valsExp has no indexing -- need to fix
    else {
        varName <- valsExp[[2]]
        orig <- get(varName, envir = envir)
        assignTarget <- valsExp
        assignTarget[[2]] <- quote(orig)
        eval(substitute(AT <- output, list(AT = assignTarget)))
        assign(varName, orig, envir = envir)
    }
}

getValuesAccess <- function(access) {
    fromCode <- makeGetCodeFromAccessorVector(access)
    if(length(fromCode)==0) return(numeric())
    sourceFromObject <- access[[1]]
    unlist(lapply(fromCode, function(i) eval(i)))

}


setValuesAccess <- function(input, access) {
    toCode <- makeSetCodeFromAccessorVector(access)
    mapInfo <- makeMapInfoFromAccessorVector(access) ## inefficient to do both of these, but we need the lengths!
    if(length(toCode)==0) return(NULL)
    sourceToObject <- access[[1]]
    nextIndex <- 0
    ## not easy to check lengths any more
    for(i in 1:length(toCode)) {
        nextLength <- mapInfo[[i]]$length
        oneValue <- input[nextIndex + (1:nextLength)]
        eval(toCode[[i]])
        nextIndex <- nextIndex + nextLength
    }
    invisible(NULL)
}


# Fill a set of nodes in a model from a vector of values
#
# Take a vector of node names for a model and fill them sequentially from the values in a vector.  Works in R and NIMBLE.
#
# @param input        a vector of values to be put in the model's nodes
# @param model       a NIMBLE model object, either compiled or uncompiled
# @param nodes       a vector of node names, allowing index blocks that will be expanded
#
# @author NIMBLE development team
# @export
# @details
# Calling \code{setValues(P, model, nodes)} will place values from P, in order, into the nodes provided for the model.  This is deprecated by
# \code{values(model, nodes) <- P}, but the development syntax of \code{setValues(P, model, nodes)} is still supported.
#
# When provided nodes are from matrices or arrays, the values will filled following column-wise order.
#
# The reverse of \code{setValues(P, model, nodes)} is \code{getValues(P, model, nodes}, which is deprecated by \code{P <- values(model, nodes)}.
#
# These functions work in R and in NIMBLE run-time code that can be compiled.
#
# @return NULL, but this function works by the side-effect of modifying the model.
setValues <- function(input, model, nodes){
	access = modelVariableAccessorVector(model, nodes, logProb = FALSE)
	setValuesAccess(input, access)
}


#' Access or set values for a set of nodes in a model
#'
#' Access or set values for a set of nodes in a NIMBLE model. 
#'
#' @aliases values values<-
#' 
#' @description
#' Get or set values for a set of nodes in a model
#'
#' @param model       a NIMBLE model object, either compiled or uncompiled
#' @param nodes       a vector of node names, allowing index blocks that will be expanded
#' @param value       value to set the node(s) to
#' @param accessorIndex For internal NIMBLE use only
#' 
#' @author NIMBLE development team
#' @export
#' @details
#' Calling \code{values(model, nodes)} returns a vector of the concatenation of values from the nodes requested
#' \code{P <- values(model, nodes)} is a newer syntax for \code{getValues(P, model, values)}, which still works and modifies P in the calling environment.
#'
#' Calling \code{values(model, nodes) <- P} sets the value of the nodes in the model, in sequential order from the vector P.
#'
#' In both uses, when requested nodes are from matrices or arrays, the values will be handled following column-wise order.
#'
#' The older function \code{getValues(P, model, nodes)} is equivalent to \code{P <- values(model, nodes)}, and the older function \code{setValues(P, model, nodes)} is equivalent to \code{values(model, nodes) <- P}
#'
#'
#' These functions work in R and in NIMBLE run-time code that can be compiled.
#'
#' @return A vector of values concatenated from the provided nodes in the model
values <- function(model, nodes, accessorIndex){
	ans <- NA
	getValues(ans, model, nodes, parent.frame())
	ans
}

## # @rdname values

#' @rdname values 
#' @export
`values<-` <- function(model, nodes, accessorIndex, value){
	setValues(value, model, nodes)
	return(model)
}

#' Copying function for NIMBLE
#' 
#' Copies values from a NIMBLE model or modelValues object to another NIMBLE model or modelValues. Work in R and NIMBLE.  The NIMBLE keyword \code{copy} is identical to \code{nimCopy}
#' 
#' @param from		Either a NIMBLE model or modelValues object
#' @param to		Either a NIMBLE model or modelValues object
#' @param nodes		Vector of one or more node names of object \code{from} that will be copied from
#' @param nodesTo	Vector of one or more node names of object \code{to} that will be copied to. If \code{nodesTo} is \code{NULL}, will automatically be set to \code{nodes}
#' @param row		If \code{from} is a modelValues, the row that will be copied from
#' @param rowTo		If \code{to} is a modelValues, the row which will be copied to. If \code{rowTo} is \code{NA}, will automatically be set to \code{row}
#' @param logProb	A logical value indicating whether the log probabilities of the given nodes should also be copied (i.e. if \code{nodes = 'x'}
#' and \code{logProb = TRUE}, then both \code{'x'} and \code{'logProb_x'} will be copied)
#' @param logProbOnly   A logical value indicating whether only the log probabilities of the given nodes should be copied (i.e. if \code{nodes = 'x'}
#' and \code{logProbOnly = TRUE}, then only \code{'logProb_x'} will be copied)
#'
#' @aliases copy
#'
#' @author Clifford Anderson-Bergman
#' @export
#' @details
#'
#' This function copies values from one or more nodes (possibly including log probabilities for nodes) between models and modelValues objects. For modelValues objects, the row must be specified. This function allows one to conveniently copy multiple nodes, avoiding having to write a loop. 
#'
#' @examples
#'	# Building model and modelValues object
#' simpleModelCode <- nimbleCode({
#'	for(i in 1:100)
#'		x[i] ~ dnorm(0,1)
#'})
#' rModel <- nimbleModel(simpleModelCode)
#' rModelValues <- modelValues(rModel)
#'
#' #Setting model nodes
#' rModel$x <- rnorm(100)
#' #Using nimCopy in R.
#' nimCopy(from = rModel, to = rModelValues, nodes = 'x', rowTo = 1)
#' 
#' #Use of nimCopy in a simple nimbleFunction
#' cCopyGen <- nimbleFunction(
#' 	setup = function(model, modelValues, nodeNames){},
#' 	run = function(){
#' 		nimCopy(from = model, to = modelValues, nodes = nodeNames, rowTo = 1)
#' 	}
#' )
#' 
#' rCopy <- cCopyGen(rModel, rModelValues, 'x')
#' \dontrun{
#' cModel <- compileNimble(rModel)
#' cCopy <- compileNimble(rCopy, project = rModel)
#' cModel[['x']] <- rnorm(100)
#' 
#' cCopy$run() ## execute the copy with the compiled function
#' }
nimCopy <- function(from, to, nodes = NULL, nodesTo = NULL, row = NA, rowTo = NA, logProb = FALSE, logProbOnly = FALSE){
    if(is.null(nodes) )
        nodes = from$getVarNames(includeLogProb = logProb) ## allNodeNames(from)
    if( inherits(from, "modelBaseClass") ){
        accessFrom = modelVariableAccessorVector(from, nodes, logProb = logProb, logProbOnly = logProbOnly)
    } else
        if(inherits(from, "modelValuesBaseClass") || inherits(from, "CmodelValues")) {
            accessFrom = modelValuesAccessorVector(from, nodes, logProb = logProb, logProbOnly = logProbOnly)
            if(is.na(row))
                stop("Error: need to supply 'row' for a modelValues copy")
            ##accessFrom$setRow(row) ## NEW ACCESSORS
        }
        else stop('argument "from" in nimCopy is neither a model nor modelValues')

    if( inherits(to, "modelBaseClass") ){
        if(is.null(nodesTo) ) 
            accessTo = modelVariableAccessorVector(to, nodes, logProb = logProb, logProbOnly = logProbOnly)
        else
            accessTo = modelVariableAccessorVector(to, nodesTo, logProb = logProb, logProbOnly = logProbOnly)
    } else
        if(inherits(to, "modelValuesBaseClass") || inherits(to, "CmodelValues")) {
            if(is.null(nodesTo) ) 
                accessTo = modelValuesAccessorVector(to, nodes, logProb = logProb, logProbOnly = logProbOnly)
            else
                accessTo = modelValuesAccessorVector(to, nodesTo, logProb = logProb, logProbOnly = logProbOnly)
            if(is.na(rowTo))
                rowTo = row
            ##accessTo$setRow(rowTo) ## NEW ACCESSORS
        } else stop('argument "to" in nimCopy is neither a model nor modelValues')

    sourceToObject <- accessTo[[1]]
    sourceFromObject <- accessFrom[[1]]
    setCode <- makeSetCodeFromAccessorVector(accessTo)
    getCode <- makeGetCodeFromAccessorVector(accessFrom)
    lengthTo <- length(setCode)
    if(length(getCode) != lengthTo)
    stop('unequal number of entries in nimCopy') 
    if(lengthTo > 0){
        for(i in 1:lengthTo){
            oneValue <- eval(getCode[[i]]) ## may have row hardwired in
            eval(setCode[[i]]) ## oneValue is hard-wired in. may have rowTo hardwired in
        }
    }
}

#' Access or set a member variable of a nimbleFunction
#' 
#' Internal way to access or set a member variable of a nimbleFunction created during \code{setup}.  Normally in NIMBLE code you would use \code{nf$var} instead of \code{nfVar(nf, var)}. 
#'
#' @aliases nfVar nfVar<-
#'
#' @description
#' Access or set a member variable of a specialized nimbleFunction, i.e. a variable passed to or created during the \code{setup} function that is used in run code or preserved by \code{setupOutputs}.  Works in R for any variable and in NIMBLE for numeric variables.
#'
#' @param nf      a specialized nimbleFunction, i.e. a function returned by executing a function returned from \code{nimbleFunction} with \code{setup} arguments
#' @param varName a character string naming a variable in the \code{setup} function.
#' @param value value to set the variable to.
#' @author NIMBLE development team
#' @export
#' @details
#' When \code{nimbleFunction} is called and a \code{setup} function is provided, then \code{nimbleFunction} returns a function.  That function is a generator that should be called with arguments to the \code{setup} function and returns another function with \code{run} and possibly other member functions.  The member functions can use objects created or passed to \code{setup}.  During internal processing, the NIMBLE compiler turns some cases of \code{nf$var} into \code{nfVar(nf, var)}. These provide direct access to setup variables (member data).  \code{nfVar} is not typically called by a NIMBLE user or programmer.
#'
#'
#' For internal access to methods of \code{nf}, see \code{\link{nfMethod}}.
#' 
#' For more information, see \code{?nimbleFunction} and the NIMBLE \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual}.
#'
#' @return whatever varName is in the nimbleFunction nf.
#' @examples
#' nfGen1 <- nimbleFunction(
#'     setup = function(A) {
#'     B <- matrix(rnorm(4), nrow = 2)
#'     setupOutputs(B) ## preserves B even though it is not used in run-code
#'    },
#'    run = function() {
#'       print('This is A', A, '\n')
#' })
#'
#' nfGen2 <- nimbleFunction(
#'   setup = function() {
#'     nf1 <- nfGen1(1000)
#'   },
#'   run = function() {
#'       print('accessing A:', nfVar(nf1, 'A'))
#'       nfVar(nf1, 'B')[2,2] <<- -1000
#'       print('accessing B:', nfVar(nf1, 'B'))
#'    })
#'        
#' nf2 <- nfGen2()
#' nf2$run()
nfVar <- function(nf, varName) {
    refClassObj <- nf_getRefClassObject(nf)
    v <- refClassObj[[varName]]
    return(v)
}

#' @rdname nfVar
#' @export
`nfVar<-` <- function(nf, varName, value) {
    refClassObj <- nf_getRefClassObject(nf)
    refClassObj[[varName]] <- value
    return(nf)
}

#' access (call) a member function of a nimbleFunction
#' 
#' Internal function for accessing a member function (method) of a nimbleFunction.  Normally a user will write \code{nf$method(x)} instead of \code{nfMethod(nf, method)(x)}.
#'
#' @param nf          a specialized nimbleFunction, i.e. one that has already had setup parameters processed
#' @param methodName  a character string giving the name of the member function to call
#'
#' @author NIMBLE development team
#' @export
#' @details
#' nimbleFunctions have a default member function called \code{run}, and may have other member functions provided via the \code{methods} argument to \code{nimbleFunction}.
#' As an internal step, the NIMBLE compiler turns \code{nf$method(x)} into \code{nfMethod(nf, method)(x)}, but a NIMBLE user or programmer would not normally need to use \code{nfMethod} directly.
#'
#' @return a function that can be called.
nfMethod <- function(nf, methodName) {
    # refClassObj <- nf_getRefClassObject(nf)
    # methodCall <- substitute(refClassObj$METHODNAME, list(METHODNAME = methodName))
    # eval(methodCall)
    eval(substitute(nf_getRefClassObject(nf)$METHODNAME, list(METHODNAME = methodName)))
}

#' Generates a weighted sample (with replacement) of ranks
#'
#' Takes a set of non-negative \code{weights} (do not need to sum to 1) and 
#' returns a sample with \code{size} elements of the integers \code{1:length(weights)}, where the probability of being sampled is proportional
#' to the value of \code{weights}. An important note is that the output vector
#' will be sorted in ascending order. Also, right now it works slightly odd syntax (see example below). Later releases of NIMBLE will contain more natural syntax. 
#'
#' If invalid weights provided (i.e. negative weights or weights sum to 1), sets output = rep(1, size) and prints warning. 
#' \code{rankSample} can be used inside nimble functions.
#'
#' @param weights		A vector of numeric weights. Does not need to sum to 1, but must be non-negative
#' @param size			Size of sample
#' @param output		An R object into which the values will be placed. See example below for proper use
#' @param silent Logical indicating whether to suppress logging information
#' @author	Clifford Anderson-Bergman
#' @export
#' @details		
#' \code{rankSample} first samples from the joint distribution \code{size} uniform(0,1) distributions by conditionally sampling from the rank statistics. This leads to 
#' a sorted sample of uniform(0,1)'s. Then, a cdf vector is constructed from weights. Because the sample of uniforms is sorted, \code{rankSample} walks
#' down the cdf in linear time and fills out the sample.
#'  
#' @examples
#' set.seed(1)
#' sampInts = NA	#sampled integers will be placed in sampInts
#' rankSample(weights = c(1, 1, 2), size = 10, sampInts)
#' sampInts
#'# [1] 1 1 2 2 2 2 2 3 3 3
#' rankSample(weights = c(1, 1, 2), size = 10000, sampInts)
#' table(sampInts)
#' #sampInts
#' #   1    2    3 
#' #2434 2492 5074 
#'
#' #Used in a nimbleFunction
#' sampGen <- nimbleFunction(setup = function(){
#' 	x = 1:2
#' },
#' run = function(weights = double(1), k = integer() ){
#' 	rankSample(weights, k, x)
#' 	returnType(integer(1))
#' 	return(x)
#' })
#' rSamp <- sampGen()
#' rSamp$run(1:4, 5)
#' #[1] 3 3 4 4 4
rankSample <- function(weights, size, output, silent = FALSE) {
    ##cat('in R version rankSample\n')
    assign(as.character(substitute(output)), .Call(C_rankSample, as.numeric(weights), as.integer(size), as.integer(output), as.logical(silent)), envir = parent.frame())
}

#' print function for use in nimbleFunctions
#'
#' @param ...  an abitrary set of arguments that will be printed in sequence.
#'
#' @details The keyword \code{print} in nimbleFunction run-time code will be automatically turned into \code{nimPrint}.  This is a function that prints its arguments in order using \code{cat} in R, or using \code{std::cout} in C++ code generated by compiling nimbleFunctions.
#' Non-scalar numeric objects can be included, although their output will be formatted slightly different in uncompiled and compiled nimbleFunctions.
#'
#' @aliases print
#'
#' @examples
#' ans <- matrix(1:4, nrow = 2) ## R code, not NIMBLE code
#' nimPrint('Answer is ', ans) ## would work in R or NIMBLE
#'
#' @seealso \code{\link{cat}}
#' @export
nimPrint <- function(...) {
    items <- list(...)
    for(i in seq_along(items)) {if(is.array(items[[i]])) print(items[[i]]) else cat(items[[i]])}
    cat('\n')
}

#' cat function for use in nimbleFunctions
#'
#' @param ...  an arbitrary set of arguments that will be printed in sequence.
#'
#' @aliases cat
#'
#' @details
#'
#' \code{cat} in nimbleFunction run-code imitates the R function \code{\link{cat}}.  It prints its arguments in order.  No newline is inserted, so include \code{"\n"} if one is desired.
#'
#' When an uncompiled nimbleFunction is executed, R's \code{cat} is used.  In a compiled nimbleFunction, a C++ output stream is used that will generally format output similarly to R's \code{cat}. Non-scalar numeric objects can be included, although their output will be formatted slightly different in uncompiled and compiled nimbleFunctions.
#'
#' In nimbleFunction run-time code, \code{cat} is identical to \code{print} except the latter appends a newline at the end.
#'
#' \code{nimCat} is the same as \code{cat}, and the latter is converted to the former when a nimbleFunction is defined.
#'
#' @examples
#' ans <- matrix(1:4, nrow = 2) ## R code, not NIMBLE code
#' nimCat('Answer is ', ans) ## would work in R or NIMBLE
#'
#' @seealso \code{\link{print}}
#' @export
nimCat <- function(...) {
    items <- list(...)
    for(i in seq_along(items)) {if(is.array(items[[i]])) print(items[[i]]) else cat(items[[i]])}
}


#' Creates numeric, integer or logical vectors for use in nimbleFunctions
#'
#' In a \code{nimbleFunction}, \code{numeric}, \code{integer} and \code{logical} are identical to \code{nimNumeric}, \code{nimInteger} and \code{nimLogical}, respectively.
#'
#' @aliases nimInteger nimLogical numeric integer logical
#'
#' @param length the length of the vector (default = 0)
#' @param value value(s) for initializing the vector (default = 0).  This may be a vector, matrix or array but will be used as a vector.
#' @param init logical, whether to initialize elements of the vector (default = TRUE)
#' @param fillZeros logical, whether to initialize any elements not filled by (possibly recycled) \code{value} with 0 (or FALSE for \code{nimLogical}) (default = TRUE)
#' @param recycle logical, whether \code{value} should be recycled to fill the entire \code{length} of the new vector (default = TRUE)
#'
#' @details
#' These functions are similar to R's \code{\link{numeric}}, \code{\link{integer}}, \code{\link{logical}} functions, but they can be used in a nimbleFunction and then compiled using \code{compileNimble}.  Largely for compilation purposes, finer control is provided over initialization behavior.  If \code{init = FALSE}, no initialization will be done, and \code{value}, \code{fillZeros} and \code{recycle} will be ignored.  If \code{init=TRUE} and \code{recycle=TRUE}, then \code{fillZeros} will be ignored, and \code{value} will be repeated (according to R's recycling rule) as much as necessary to fill a vector of length \code{length}.  If \code{init=TRUE} and \code{recycle=FALSE}, then if \code{fillZeros=TRUE}, values of 0 (or FALSE for \code{nimLogical}) will be filled in after \code{value} up to length \code{length}.  Compiled code will be more efficient if unnecessary initialization is not done, but this may or may not be noticeable depending on the situation.
#' 
#' When used in a \code{nimbleFunction} (in \code{run} or other member function), \code{numeric}, \code{integer} and \code{logical} are immediately converted to \code{nimNumeric}, \code{nimInteger} and \code{nimLogical}, respectively.  
#' 
#' @author Daniel Turek, Chris Paciorek, Perry de Valpine
#' @aliases numeric
#' @seealso \code{\link{nimMatrix}}, \code{\link{nimArray}}
#' @export
nimNumeric <- function(length = 0, value = 0, init = TRUE, fillZeros = TRUE, recycle = TRUE) {
    fillValue <- makeFillValue(value, 'double', init)
    makeReturnVector(fillValue, length, recycle)
}

#' @rdname nimNumeric
#' @export
nimInteger <- function(length = 0, value = 0, init = TRUE, fillZeros = TRUE, recycle = TRUE) {
    fillValue <- makeFillValue(value, 'integer', init)
    makeReturnVector(fillValue, length, recycle)
}

#' @rdname nimNumeric
#' @export
nimLogical <- function(length = 0, value = 0, init = TRUE, fillZeros = TRUE, recycle = TRUE) {
    fillValue <- makeFillValue(value, 'logical', init)
    makeReturnVector(fillValue, length, recycle)
}

#' Creates matrix or array objects for use in nimbleFunctions
#' 
#' In a \code{nimbleFunction}, \code{matrix} and \code{array} are identical to \code{nimMatrix} and \code{nimArray}, respectively
#'
#' @aliases nimArray matrix array
#'
#' @param value value(s) for initialization (default = 0).  This can be a vector, matrix or array, but it will be used as a vector.
#' @param nrow the number of rows in a matrix (default = 1)
#' @param ncol the number of columns in a matrix (default = 1)
#' @param dim vector of dimension sizes in an array (default = \code{c(1, 1)})
#' @param init logical, whether to initialize values (default = \code{TRUE})
#' @param fillZeros logical, whether to initialize any elements not filled by (possibly recycled) \code{value} with 0 (or \code{FALSE} for \code{nimLogical}) (default = \code{TRUE})
#' @param recycle logical, whether \code{value} should be recycled to fill the entire contents of the new object (default = \code{TRUE})
#' @param type character representing the data type, i.e. \code{'double'}, \code{'integer'}, or \code{'logical'} (default = \code{'double'})
#' @param nDim number of dimensions in an array.  This is only necessary for \code{compileNimble} if the length of \code{dim} cannot be determined during compilation.
#'
#' @details
#' These functions are similar to R's \code{\link{matrix}} and \code{\link{array}} functions, but they can be used in a nimbleFunction and compiled using \code{compileNimble}.  Largely for compilation purposes, finer control is provided over initialization behavior, similarly to \code{\link{nimNumeric}}, \code{\link{nimInteger}}, and \code{\link{nimLogical}}. If \code{init = FALSE}, no initialization will be done, and \code{value}, \code{fillZeros} and \code{recycle} will be ignored.  If \code{init=TRUE} and \code{recycle=TRUE}, then \code{fillZeros} will be ignored, and \code{value} will be repeated (according to R's recycling rule) as much as necessary to fill the object.  If \code{init=TRUE} and \code{recycle=FALSE}, then if \code{fillZeros=TRUE}, values of 0 (or FALSE for \code{nimLogical}) will be filled in after \code{value}.  Compiled code will be more efficient if unnecessary initialization is not done, but this may or may not be noticeable depending on the situation.
#'
#' When used in a \code{nimbleFunction} (in \code{run} or other member function), \code{matrix} and \code{array} are immediately converted to \code{nimMatrix} and \code{nimArray}, respectively.
#'
#' The \code{nDim} argument is only necessary for a use like \code{dim <- c(2, 3, 4); A <- nimArray(0, dim = dim, nDim = 3)}.  It is necessary because the NIMBLE compiler must determine during compilation that \code{A} will be a 3-dimensional numeric array.  However, the compiler doesn't know for sure what the length of \code{dim} will be at run time, only that it is a vector.  On the other hand,   \code{A <- nimArray(0, dim = c(2, 3, 4))} is allowed because the compiler can directly determine that a vector of length three is constructed inline for the \code{dim} argument.
#' 
#' @author Daniel Turek and Perry de Valpine
#' @seealso \code{\link{nimNumeric}} \code{\link{nimInteger}} \code{\link{nimLogical}}
#' @export
nimMatrix <- function(value = 0, nrow = NA, ncol = NA, init = TRUE, fillZeros = TRUE, recycle = TRUE, type = 'double') {
    ## the -1's are used because nimble does not allow both missingness and default value
    ## but R's matrix function relies on both possibilities
    fillValue <- makeFillValue(value, type, init)
    mnrow <- missing(nrow) || is.na(nrow)
    mncol <- missing(ncol) || is.na(ncol)
    if(mnrow)
        if(mncol) {
            base::matrix(fillValue)
        } else {
            nrow <- ceiling( length(fillValue) / ncol )
            fillValue <- makeReturnVector(fillValue, nrow * ncol, recycle)
            base::matrix(fillValue, ncol = ncol, nrow = nrow)
        }
    else
        if(mncol) {
            ncol <- ceiling( length(fillValue) / nrow )
            fillValue <- makeReturnVector(fillValue, nrow * ncol, recycle)
            base::matrix(fillValue, nrow = nrow, ncol = ncol)
        } else {
            fillValue <- makeReturnVector(fillValue, ncol*nrow, recycle)
            base::matrix(fillValue, nrow = nrow, ncol = ncol)
        }
}


#' @rdname nimMatrix
#' @export
nimArray <- function(value = 0, dim = c(1, 1), init = TRUE, fillZeros = TRUE, recycle = TRUE, nDim, type = 'double') {
    if(!missing(nDim)) dim <- dim[1:nDim]
    fillValue <- makeFillValue(value, type, init)
    fillValue <- makeReturnVector(fillValue, prod(dim), recycle)
    if(length(dim) == 1) fillValue
    else base::array(fillValue, dim)
}

makeFillValue <- function(value, type, init) {
    fillValue <- if(init) value else 0
    fillValueTyped <- switch(type,
                             double = as.numeric(fillValue),
                             integer = as.integer(fillValue),
                             logical = as.logical(fillValue),
                             stop('unknown type argument'))
    return(fillValueTyped)
}

makeReturnVector <- function(fillValue, length, recycle) {
    if(length(fillValue) == 1) {
        if(recycle)
            rep(fillValue, length)
        else
            c(fillValue, as(rep(0, max(length-1, 0)), class(fillValue)))
    }
    else {
        if(length(fillValue) != length) {
            if(length(fillValue) < length) {
                ##warning(paste0("Not enough values provided for vector of length ",length, ".")) 
                if(recycle)
                    rep(fillValue, length.out = length)
                else
                    c(fillValue, as(rep(0, length-length(fillValue)), class(fillValue)))
            } else {
                ##warning(paste0("Too many values provided for vector of length ",length, ".")) 
                fillValue[1:length]
            }
        } else {
            fillValue
        }
    }
}


#' Explicitly declare a variable in run-time code of a nimbleFunction
#'
#' Explicitly declare a variable in run-time code of a nimbleFunction, for cases when its dimensions cannot be inferred before it is used.  Works in R and NIMBLE.
#'
#' @param name     Name of a variable to declare, without quotes
#' @param def      NIMBLE type declaration, of the form \code{TYPE(nDim, sizes)}, where \code{TYPE} is \code{integer}, \code{double}, or \code{logical}, \code{nDim} is the number of dimensions, and \code{sizes} is an optional vector of sizes concatenated with \code{c}.  If \code{nDim} is omitted, it defaults to 0, indicating a scalar.  If sizes are provided, they should not be changed subsequently in the function, including by assignment.  Omitting \code{nDim} results in a scalar.  For \code{logical}, only scalar is currently supported.
#'
#' @author NIMBLE development team
#' @export
#' @details
#' In a run-time function of a nimbleFunction (either the \code{run} function or a function provided in \code{methods} when calling \code{nimbleFunction}), the dimensionality and numeric type of a variable is inferred when possible from the statement first assigning into it.  E.g. \code{A <- B + C} infers that \code{A} has numeric types, dimensions and sizes taken from \code{B + C}.  However, if the first appearance of \code{A} is e.g. \code{A[i] <- 5}, \code{A} must have been explicitly declared.  In this case, \code{declare(A, double(1))} would make \code{A} a 1-dimensional (i.e. vector) double.
#'
#' When sizes are not set, they can be set by a call to \code{setSize} or by assignment to the whole object.  Sizes are not automatically extended if assignment is made to elements beyond the current sizes.  In compiled nimbleFunctions doing so can cause a segfault and crash the R session.
#'
#'
#' This part of the NIMBLE language is needed for compilation, but it also runs in R.  When run in R, is works by the side effect of creating or modifying \code{name} in the calling environment.
#'
#' @examples
#' declare(A, logical())             ## scalar logical, the only kind allowed
#' declare(B, integer(2, c(10, 10))) ## 10 x 10 integer matrix
#' declare(C, double(3))             ## 3-dimensional double array with no sizes set.
declare <- function(name, def){
    defCode <- substitute(def)
    name <- substitute(name)
    if(exists(as.character(name), parent.frame(), inherits = FALSE)) return(invisible(NULL))
    value <- if(defCode[[1]] == 'logical') FALSE else 0
    if(length(defCode) == 1){  ## no arg, like double()
        assign(as.character(name), value, envir = parent.frame() )
        return()
    }
    nDim = eval(defCode[[2]], envir = parent.frame() )
    if(nDim == 0 ){ ## like double(0)
        assign(as.character(name), value, envir = parent.frame() )
        return()
    }
    dims = rep(1, nDim)
    if(length(defCode) == 3) ## notation like double(2, c(3, 5))
        dims = eval(defCode[[3]], envir = parent.frame() )
    else {
        if(length(defCode) == 2 + nDim) {
            dims <- numeric(nDim)
            for(i in 1:nDim)
                dims[i] <- eval(defCode[[2 + i]], envir = parent.frame())
        }
    }
    if(length(dims) != nDim)
        stop('in declare, dimensions are not declared properly')
    if(prod(dims)==1) {
        if(length(dims) == 1)
            return(assign(as.character(name), value, envir = parent.frame()))
    }
    assign(as.character(name), array(value, dim = dims), envir = parent.frame() )
}

#' Determine if any values in a vector are NA or NaN
#'
#' NIMBLE language functions that can be used in either compiled or uncompiled
#' nimbleFunctions to detect if there are any NA or NaN values in a vector.
#'
#' @param x vector of values
#'
#' @aliases any_nan
#' @author NIMBLE Development Team
#' @export
any_na <- function(x) any(is.na(x))

## Retained for backwards compatibility
is.na.vec <- function(x) any(is.na(x))

#' @rdname any_na
#' @export
any_nan <- function(x) any(is.nan(x))

## Retained for backwards compatibility
is.nan.vec <- function(x) any(is.nan(x))

#' @export
nimRound <- round

#' Nimble wrapper around R's builtin \code{\link{optim}}.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param fn  A function to be minimized (or maximized), with first argument the
#'            vector of parameters over which minimization is to take place. It
#'            should return a scalar result.
#' @param gr  A function to return the gradient for the "BFGS", "CG" and "L-BFGS-B" methods.
#' @param ... IGNORED
#' @param method The method to be used. See `Details` section of \code{\link{optim}}. One of:
#'               "Nelder-Mead", "BFGS", "CG", "L-BFGS-B".
#'               Note that the R methods "SANN", "Brent" are not supported.
#' @param lower Vector or scalar of lower bounds for parameters.
#' @param upper Vector or scalar of upper bounds for parameters.
#' @param control A list of control parameters. See \code{Details} section of \code{\link{optim}}.
#' @param hessian Logical. Should a Hessian matrix be returned?
#'
#' @return \code{\link{optimResultNimbleList}}
#' @seealso \code{\link{optim}}
#' @export
#' @examples
#' \dontrun{
#' objectiveFunction <- nimbleFunction(
#'     run = function(par = double(1)) {
#'         return(sum(par) * exp(-sum(par ^ 2) / 2))
#'         returnType(double(0))
#'     }
#' )
#' optimizer <- nimbleFunction(
#'     run = function(method = character(0), fnscale = double(0)) {
#'         control <- optimDefaultControl()
#'         control$fnscale <- fnscale
#'         par <- c(0.1, -0.1)
#'         return(optim(par, objectiveFunction, method = method, control = control))
#'         returnType(optimResultNimbleList())
#'     }
#' )
#' cOptimizer <- compileNimble(optimizer)
#' cOptimizer(method = 'BFGS', fnscale = -1)
#' }
nimOptim <- function(par, fn, gr = "NULL", ..., method = "Nelder-Mead", lower = -Inf, upper = Inf,
                     control = nimOptimDefaultControl(), hessian = FALSE) {
    # Tweak parameters.
    if (identical(gr, "NULL")) gr <- NULL
    defaultControl <- nimOptimDefaultControl()
    Rcontrol <- list()
    ## Only enter non-default values into Rcontrol.
    ## R's optim will fill in the rest.
    ## This should match compiled handling of control list
    ## parameters in NimOptimProblem::solve in nimOptim.cpp.
    for (name in defaultControl$nimbleListDef$types$vars) {
        if (!identical(control[[name]], defaultControl[[name]])) {
            Rcontrol[[name]] <- control[[name]]
        }
    }
    result <- optim(par, fn, gr = gr, ..., method = method,
                    lower = lower, upper = upper, control = Rcontrol, hessian = hessian)
    nimResult <- do.call(optimResultNimbleList$new, result)
    # Tweak result value to exactly match C++ behavior.
    nimResult$counts <- unname(nimResult$counts)
    if (is.null(nimResult$message)) {
        nimResult$message <- ''
    }
    if (is.null(nimResult$hessian)) {
        nimResult$hessian <- matrix(nrow = 0, ncol = 0)
    }
    return(nimResult)
}


## Interface for optimizing model logProb with respect to parameters
nimOptim_model <- function(model, wrt, nodes, use.gr = TRUE, method = "BFGS",
                           lower = -Inf, upper = Inf, control = nimOptimDefaultControl(), hessian = FALSE) {
    par <- values(model, wrt)
    fn <- function(par) {
        values(model, wrt) <- par
        model$calculate(nodes)
    }
    if(use.gr) {
        gr <-  function(par) {
            numDeriv::grad(fn, par)
        }
    } else gr <- "NULL"
    control$fnscale <- -1
    nimOptim(par = par, fn = fn, gr = gr,
             method = method, lower = lower, upper = upper,
             control = control, hessian = hessian)
}


#' Creates a deafult \code{control} argument for \code{\link{optim}} (just an empty list).
#'
#' @return an empty \code{list}.
#' @seealso \code{\link{nimOptim}}, \code{\link{optim}}
#' @export
optimDefaultControl <- function() {
    return(list())
}

#' Creates a deafult \code{control} argument for \code{\link{nimOptim}}.
#'
#' @return \code{\link{optimControlNimbleList}}
#' @seealso \code{\link{nimOptim}}, \code{\link{optim}}
#' @export
nimOptimDefaultControl <- function() {
  control <- optimControlNimbleList$new()
  control$trace <- 0
  control$fnscale <- 1
  control$parscale <- NA # Must be filled in to length of par
  control$ndeps <- NA    # Ditto
  control$maxit <- NA  ## The default value depends on method.
  control$abstol = -Inf
  control$reltol <- sqrt(.Machine$double.eps)
  control$alpha <- 1.0
  control$beta <- 0.5
  control$gamma <- 2.0
  control$REPORT <- NA # Method dependent and not used in compiled version
  control$type <- 1
  control$lmm <- 5
  control$factr <- 1e7
  control$pgtol <- 0
  control$tmax <- 10
  control$temp <- 10
  return(control)
}
