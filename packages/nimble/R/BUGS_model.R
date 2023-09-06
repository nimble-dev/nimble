#' Class \code{modelBaseClass}
#' @aliases modelBaseClass getVarNames getNodeNames topologicallySortNodes resetData setData isData isEndNode getDistribution isDiscrete isBinary isStoch isDeterm isTruncated isUnivariate isMultivariate getDimension getDependencies getDependenciesList getDownstream expandNodeNames setInits checkConjugacy getCode getMacroParameters getMacroInits getConstants newModel getParents [[,modelBaseClass-method [[<-,modelBaseClass-method initializeInfo
#' @export
#' @description
#' This class underlies all NIMBLE model objects: both R model objects created from the return value of nimbleModel(), and compiled model objects.
#' The model object contains a variety of member functions, for providing information about the model structure, setting or querying properties of the model, or accessing various internal components of the model.
#' These member functions of the modelBaseClass are commonly used in the body of the \code{setup} function argument to nimbleFunction(), to aid in preparation of node vectors, nimbleFunctionLists, and other runtime inputs.
#' See documentation for \code{\link{nimbleModel}} for details of creating an R model object.
#' @author Daniel Turek
#' @examples
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x[1] ~ dnorm(mu, 1)
#'     x[2] ~ dnorm(mu, 1)
#' })
#' Rmodel <- nimbleModel(code)
#' modelVars <- Rmodel$getVarNames()   ## returns 'mu' and 'x'
#' modelNodes <- Rmodel$getNodeNames()   ## returns 'mu', 'x[1]' and 'x[2]'
#' Rmodel$resetData()
#' Rmodel$setData(list(x = c(1.2, NA)))   ## flags only 'x[1]' node as data
#' Rmodel$isData(c('mu', 'x[1]', 'x[2]'))   ## returns c(FALSE, TRUE, FALSE)
#' @seealso \code{\link{initializeModel}}
modelBaseClass <- setRefClass('modelBaseClass',
                              fields = list(
                                  modelDef = 'ANY',
                                  nodes = 'ANY',       #list
                                  vars = 'ANY',
                                  graph = 'ANY',
                                  defaultModelValues = 'ANY',
                                  name = 'ANY', 		#character
                                  isDataVars = 'ANY', #list           ## list with the dimensions of isData_vars
                                  isDataEnv = 'ANY',	#environment      ## environment holding 'logical' objects, with isData flags
                                  predictiveNodeIDs = 'ANY',          ## integer vector giving the graph IDs of predictive nodes
                                  predictiveRootNodeIDs = 'ANY',      ## integer vector giving the graph IDs of predictive root nodes
                                  classEnvironment = 'ANY', # environment in which the reference classes will be defined
                                  origData = 'ANY',
                                  origInits = 'ANY',
                                  nimbleProject = 'ANY'
                                  ),
                              methods = list(
                                  calculate = function(nodes) {
'
See `help(calculate)`
'
                                      nimble::calculate(.self, nodes)
                                  },
                                  calculateDiff = function(nodes) {
'
See `help(calculateDiff)`
'
                                      nimble::calculateDiff(.self, nodes)
                                  },
                                  getLogProb = function(nodes) {
'
See `help(getLogProb)`
'
                                      nimble::getLogProb(.self, nodes)
                                  },
                                  simulate = function(nodes, includeData = FALSE) {
'
See `help(simulate)`                                                                                                                                                                                                                                
'
                                      nimble::simulate(.self, nodes, includeData)
                                  },
                                  getParam = function(node, param) {                                                                                                                                                                                 '
See `help(getParam)`                                                                                                                                                                                                                               
'
                                      nimble::getParam(.self, node, param)
                                  },
                                  getBound = function(node, bound) {
'                                                                                                                                                                                                                                                    
See `help(getBound)`                                                                                                                                                                                                                               
'
                                      nimble::getBound(.self, node, bound)
                                  },
                                  getGraph = function() graph,
                                  setGraph = function(value) graph <<- value,
                                  plotGraph = function() igraph::plot.igraph(graph),
                                  plot      = function() plotGraph(),
                                  getModelDef = function() modelDef,
                                  setModelDef = function(value) modelDef <<- value,
                                  getMaps = function(mapName, all = FALSE){
                                      if(all == TRUE)		return(modelDef$maps)
                                      return(modelDef$maps[[mapName]])
                                  },
                                  getCode = function() {
                                      '
Return the code for a model after processing if-then-else statements, expanding macros, and replacing some keywords (e.g. nimStep for step) to avoid R ambiguity.
'
                                      modelDef$BUGScode
                                  },
                                  getMacroParameters = function(includeLHS = TRUE, includeRHS = TRUE, includeDeterm = TRUE, 
                                                                includeStoch = TRUE, includeIndices = FALSE) {
                                      '
Return list of new parameters generated by macros.
'
                                  dc <- modelDef$declInfo
                                  out <- modelDef$macroParameters

                                  #RHS/LHS
                                  if(!includeRHS){
                                    out <- lapply(out, function(x){
                                      lapply(x, function(z){
                                        z$RHS <- NULL
                                        z
                                      })
                                    })
                                  }
                                  if(!includeLHS){
                                    out <- lapply(out, function(x){
                                      lapply(x, function(z){
                                        z$LHS <- NULL
                                        z
                                      })
                                    })
                                  }

                                  # Indices
                                  if(!includeIndices){
                                    idx <- sapply(unlist(lapply(dc, function(x) x$indexExpr)), deparse)
                                    out <- lapply(out, function(x){
                                      lapply(x, function(y){
                                        yout <- lapply(y, function(z){
                                          zout <- z[! z %in% idx]
                                          if(length(zout) == 0) zout <- NULL
                                          zout
                                        })
                                      })
                                    })
                                  }

                                  # Type
                                  varnames <- sapply(unlist(lapply(dc, function(x) x$targetVarExpr)), deparse)
                                  type <- unlist(lapply(dc, function(x) x$type))
                                  names(type) <- varnames
                                  determ_pars <- names(type)[type == "determ"]
                                  stoch_pars <- names(type)[type == "stoch"]
                                  if(!includeDeterm){
                                    out <- lapply(out, function(x){
                                      lapply(x, function(y){
                                        yout <- lapply(y, function(z){
                                          zout <- z[! z %in% determ_pars]
                                          if(length(zout) == 0) zout <- NULL
                                          zout
                                        })
                                      })
                                    })
                                  }
                                  if(!includeStoch){
                                    out <- lapply(out, function(x){
                                      lapply(x, function(y){
                                        yout <- lapply(y, function(z){
                                          zout <- z[! z %in% stoch_pars]
                                          if(length(zout) == 0) zout <- NULL
                                          zout
                                        })
                                      })
                                    })
                                  }

                                  # Cleanup
                                  lapply(out, function(x){
                                    lapply(x, function(y){
                                      y <- unlist(y)
                                      if(length(y) == 0) x <- NULL
                                      names(y) <- NULL
                                      unique(y)
                                    })
                                  })
                                  },
                                  getMacroInits = function(){
                                    modelDef$macroInits
                                  },
                                  getConstants = function(){
                                    modelDef$constantsList
                                  },
                                  isEndNode = function(nodes){  #Note: it says nodes, but graphIDs are fine too. Actually they are better
                                                                          '
Determines whether one or more nodes are end nodes (nodes with no stochastic dependences)

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is logical vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.
'

                                      nodeNames <- nodes  # needed so don't have local assignment into 'nodes'
                                      nms <- nodeNames
				      if(is.character(nodeNames)) {
                                          nms <- expandNodeNames(nodeNames, unique = FALSE)
                                          nodeNames = expandNodeNames(nodeNames, returnType = 'ids', unique = FALSE)
                                      }
                                      out <- modelDef$maps$isEndNode_byGID[nodeNames]
                                      names(out) <- nms
                                      return(out)
                                  },

                                  ## returns the type of one or more node names, e.g., 'stoch' or 'determ'
                                  getNodeType = function(nodes) {
                                      graphIDs <- modelDef$nodeName2GraphIDs(nodes, unique = FALSE)
                                      types <- getMaps('types')[graphIDs]
                                      return(types)
                                  },

                                  ## returns the declaration ID corresponding to nodes
                                  getDeclID = function(nodes) {
                                      graphIDs <- if(is.character(nodes))
                                                      modelDef$nodeName2GraphIDs(nodes, unique = FALSE)
                                                  else
                                                      nodes
                                      declIDs <- getMaps('graphID_2_declID')[graphIDs]
                                      return(declIDs)
                                  },

                                  ## returns a list of the declInfo objects corresponding to nodes
                                  getDeclInfo = function(nodes) {
                                      declIDs <- getDeclID(nodes)
                                      declInfos <- modelDef$declInfo[declIDs]
                                      return(declInfos)
                                  },

                                  getUnrolledIndicesList = function(node) {
                                      di <- getDeclInfo(node)[[1]]
                                      if(length(which(di$nodeFunctionNames == node)) != 1)
                                          stop('something went wrong with Daniel\'s understanding of newNimbleModel')
                                      unrolledRowNumber <- which(di$nodeFunctionNames == node)
                                      indicesMatrix <- di$unrolledIndicesMatrix
                                      if(nrow(indicesMatrix) == 0) {
                                          if(unrolledRowNumber > 1) stop('something went wrong with Daniel\'s understanding of newNimbleModel')
                                          return(list())
                                      }
                                      unrolledIndices <- as.list(indicesMatrix[unrolledRowNumber, ])
                                      return(unrolledIndices)
                                  },

                                  ## returns the text for the distribution of a stochastic node, e.g., 'dnorm'
                                  getDistribution = function(nodes) {
                                                                          '
Returns the names of the distributions for the requested node or nodes

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.
'
                                      nodeNames <- expandNodeNames(nodes, unique = FALSE, returnType = "ids")
                                      out <- sapply(nodeNames, function(x)
				      	getDeclInfo(x)[[1]]$getDistributionName())
                                      names(out) <- modelDef$maps$graphID_2_nodeName[nodeNames]
                                      return(out)
                                  },

                                  ## returns the expr corresponding to 'param' in the distribution of 'node'
                                  getParamExpr = function(node, param) {
                                      di <- getDeclInfo(node)[[1]]
                                      if(di$type != 'stoch')  stop('getting parameter expression for a non-stochastic node')
                                      if(param %in% names(di$valueExprReplaced)) {
                                          expr <- di$valueExprReplaced[[param]]
                                      } else if(param %in% names(di$altParamExprs)) {
                                          expr <- di$altParamExprs[[param]]
                                      } else stop('getting a parameter not present in stochastic node')
                                      unrolledIndices <- getUnrolledIndicesList(node)
                                      subExpr <- eval(substitute(substitute(EXPR, unrolledIndices), list(EXPR = expr)))
                                      return(subExpr)
                                  },

                                  ##  returns the entire RHS valueExpr for 'node'
                                  getValueExpr = function(node) {
                                      expr <- getDeclInfo(node)[[1]]$valueExprReplaced
                                      unrolledIndices <- getUnrolledIndicesList(node)
                                      subExpr <- eval(substitute(substitute(EXPR, unrolledIndices), list(EXPR = expr)))
                                      return(subExpr)
                                  },

                                  isMultivariate = function(nodes) {
                                      '
Determines whether one or more nodes represent multivariate nodes

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a logical vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.
'
                                      nodeNames <- expandNodeNames(nodes, unique = FALSE)
                                      multi <- sapply(nodeNames, function(node) getDistributionInfo(getDistribution(node))$types$value$nDim > 0)
                                      ##multi <- rep(FALSE, length(nodeNames))
                                      ##for(i in seq_along(nodeNames)) {
                                      ##    nodeExpanded <- expandNodeNames(nodeNames[i], returnScalarComponents = TRUE)
                                      ##    if(length(nodeExpanded) > 1) multi[i] <- TRUE
                                      ##}
                                      return(multi)
                                  },

                                  isDiscrete = function(nodes) {
                                                                          '
Determines whether one or more nodes represent discrete random variables

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.
'
                                      dist <- getDistribution(nodes)
                                      # explicit reference to namespace needed as class definition objects inheriting from modelBaseClass not in namespace
                                      discrete <- sapply(dist, nimble::isDiscrete)
                                      #discrete <- nimble:::getDistributionInfo(dist)$discrete
                                      return(discrete)
                                  },

                                  isBinary = function(nodes) {
                                    '
Determines whether one or more nodes represent binary random variables

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.
'
                                      nodeNames <- expandNodeNames(nodes, unique = FALSE)  # needed below but duplicates what happens in getDistribution
                                      dist <- getDistribution(nodeNames)

                                      binary <- rep(FALSE, length(dist))
                                      names(binary) <- names(dist)
                                      binary[dist == 'dbern'] <- TRUE
                                      binomInds <- which(dist == 'dbin')
                                      if(length(binomInds)) {
                                          tmp <- sapply(binomInds, function(ind) getParamExpr(nodeNames[ind], 'size') == 1)
                                          binary[binomInds[tmp]] <- TRUE
                                      }
                                      binary[is.na(dist)] <- NA
                                      return(binary)
                                  },

                                # user-facing, in contrast to getNodeTypes
                                isStoch = function(nodes, nodesAlreadyExpanded = FALSE) {
                                  '
Determines whether one or more nodes are stochastic

Arguments:

nodes: A character vector specifying one or more node or variable names.
nodesAlreadyExpanded: Boolean argument indicating whether `nodes` should be expanded. Generally intended for internal use. Default is `FALSE`.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.
'
                                  ## This check handles strange case of overlapping RHSonly nodes
                                  ## when called from `setData` and checking deterministic data nodes;
                                  ## e.g., see test-mcmc.R 'conjugate MVN with ragged dependencies'
                                  if(!nodesAlreadyExpanded)
                                      nodeNames <- expandNodeNames(nodes, unique = FALSE) else nodeNames <- nodes
                                  type <- getNodeType(nodeNames)
                                  out <- type == "stoch"
                                  names(out) <- nodeNames
                                  return(out)
                                },

                                isDeterm = function(nodes, nodesAlreadyExpanded = FALSE) {
                                  '
Determines whether one or more nodes are deterministic

Arguments:

nodes: A character vector specifying one or more node or variable names.
nodesAlreadyExpanded: Boolean argument indicating whether `nodes` should be expanded. Generally intended for internal use. Default is `FALSE`.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.
'
                                  ## See comment in `isStoch`.
                                  if(!nodesAlreadyExpanded)
                                      nodeNames <- expandNodeNames(nodes, unique = FALSE) else nodeNames <- nodes
                                  type <- getNodeType(nodeNames)
                                  out <- type == "determ"
                                  names(out) <- nodeNames
                                  return(out)
                                },

                                isTruncated = function(nodes) {
                                                                      '
Determines whether one or more nodes are truncated

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent nodes names, so the length of the output may be longer than that of the input
'
                                      nodeNames <- expandNodeNames(nodes, unique = FALSE)
                                      out <- sapply(nodeNames, function(x)
	    					getDeclInfo(x)[[1]]$isTruncated())
                                      names(out) <- nodeNames
                                      return(out)
                                  },

                                  isUnivariate = function(nodes) {
                                                                                                       '
Determines whether one or more nodes represent univariate random variables

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent nodes names, so the length of the output may be longer than that of the input
'

                                      nodeNames <- expandNodeNames(nodes, unique = FALSE)
                                      dists <- getDistribution(nodeNames)
			   	      dims <- sapply(dists, getDimension)
                                      out <- dims == 1
                                      names(out) <- nodeNames
                                      return(out)
                                  },

                                  getDimension = function(node, params = NULL, valueOnly = is.null(params)
                                    && !includeParams, includeParams = !is.null(params)) {
                                                                                                                                           '
Determines the dimension of the value and/or parameters of the node

Arguments:

node: A character vector specifying a single node

params: an optional character vector of names of parameters for which dimensions are desired (possibly including \'value\' and alternate parameters)

valueOnly: a logical indicating whether to only return the dimension of the value of the node

includeParams: a logical indicating whether to return dimensions of parameters. If TRUE and \'params\' is NULL then dimensions of all parameters, including the dimension of the value of the node, are returned

Details: The return value is a numeric vector with an element for each parameter/value requested.
'

                                      dist <- getDistribution(node)
                                      if(length(dist) > 1)
                                          stop("getDimension: 'node' should be a single node in the model")
                                      dim <- nimble::getDimension(dist, params, valueOnly, includeParams)
                                      return(dim)
                                  },

                                  getVarNames = function(includeLogProb = FALSE, nodes) {
                                      '
Returns the names of all variables in a model, optionally including the logProb variables

Arguments:

logProb: Logical argument specifying whether or not to include the logProb variables.  Default is FALSE.

nodes: An optional character vector supplying a subset of nodes for which to extract the variable names and return the unique set of variable names
'
                                      if(missing(nodes)){
                                          if(includeLogProb) ans <- modelDef$varNames
                                          else ans <- names(modelDef$varInfo)
    	                              } else {
                                          ans <- unique(nimble:::removeIndexing(nodes))
                                          if(!all(ans %in% modelDef$varNames))
                                              stop(c('invalid node names provided to model$getVarNames') )
                                      }
                                      ## "includeData" argument to getVarNames (with default = TRUE)
                                      ## was removed by consensus, March 2017.
    	                              return(ans)
                                  },

                                  getNodeFunctions = function(nodes) {
                                      gids <- modelDef$nodeName2GraphIDs(nodes, unique = FALSE)
                                      dclids <- modelDef$graphIDs2indexedNodeInfo(gids)$declIDs
                                      if(length(dclids) == 1)
                                          return(nodeFunctions[[dclids]])
                                      else
                                          return(nodeFunctions[dclids])
                                  },

                                  getNodeNames = function(determOnly = FALSE, stochOnly = FALSE,
                                                          includeData = TRUE, dataOnly = FALSE, includeRHSonly = FALSE,
                                                          topOnly = FALSE, latentOnly = FALSE, endOnly = FALSE,
                                                          includePredictive = TRUE, predictiveOnly = FALSE,
                                                          returnType = 'names', returnScalarComponents = FALSE) {
                                      '
Returns a character vector of all node names in the model, in topologically sorted order.  A variety of logical arguments allow for flexible subsetting of all model nodes.

Arguments:

determOnly: Logical argument specifying whether to return only deterministic nodes.  Default is FALSE.

stochOnly: Logical argument specifying whether to return only stochastic nodes.  Default is FALSE.

includeData: Logical argument specifying whether to include \'data\' nodes (set via the member method setData).  Default is TRUE.

dataOnly: Logical argument specifying whether to return only \'data\' nodes.  Default is FALSE.

includeRHSonly: Logical argument specifying whether to include right-hand-side-only nodes (model nodes which never appear on the left-hand-side of ~ or <- in the model code).  Default is FALSE.

topOnly: Logical argument specifying whether to return only top-level nodes from the hierarchical model structure.

latentOnly: Logical argument specifying whether to return only latent (mid-level) nodes from the hierarchical model structure.

endOnly: Logical argument specifying whether to return only end nodes from the hierarchical model structure.

includePredictive: Logical argument specifying whether to include predictive nodes (stochastic nodes, which themselves are not data and have no downstream stochastic dependents which are data) from the hierarchical model structure.

predictiveOnly: Logical argument specifying whether to return only predictive nodes (stochastic nodes, which themselves are not data and have no downstream stochastic dependents which are data) from the hierarchical model structure.

returnType: Character argument specific type object returned. Options are \'names\' (returns character vector) and \'ids\' (returns numeric graph IDs for model)

returnScalar Componenets: Logical argument specifying whether multivariate nodes should return full node name (i.e. \'x[1:2]\') or should break down into scalar componenets (i.e. \'x[1]\' and \'x[2]\')

Details: Multiple logical input arguments may be used simultaneously.  For example, `model$getNodeNames(endOnly = TRUE, dataOnly = TRUE)` will return all end-level nodes from the model which are designated as \'data\'.
'
                                      ## Part of fix for Issue #340:
                                      ## In previous versions, we started with validValues all TRUE
                                      ## and LHSinferred nodes (which should never be returned) were filtered out in the final expandNodeNames call
                                      ## Now a LHSinferred node can have a name reflecting that it is a split vertex, like mu[ 1 %.s$ 5]
                                      ##   which means its elements begin at 1 and end at 5 but are not contiguous (some elements are not included)
                                      ## Such a name can be parsed but evaluating it in one of the "vars2..." environments will fail
                                      ## so it needs to be omitted from the call to expandNodeNames.
                                      ## It turns out ot be easy to do that by filtering out LHSinferred nodes at initialization of validValues
                                      ## validValues is a boolean vector aligned with modelDef$maps$nodesNames and allied vectors
                                      validValues <- modelDef$maps$types != "LHSinferred"
                                      ## Apply a series of filters of what can be included
                                      if(!includeRHSonly)               validValues[modelDef$maps$types == 'RHSonly'] <- FALSE
                                      if(determOnly)                    validValues[modelDef$maps$types != 'determ']  <- FALSE
                                      if(stochOnly)                     validValues[modelDef$maps$types != 'stoch']   <- FALSE

                                      if(!includeData | dataOnly) {
                                          boolIsData <- rep(FALSE, length(modelDef$maps$graphIDs))
                                          possibleDataIDs <- modelDef$maps$graphIDs[modelDef$maps$types == 'RHSonly' | modelDef$maps$types == 'stoch']
                                          boolIsData[possibleDataIDs] <- isDataFromGraphID(possibleDataIDs)
                                          if(!includeData)              validValues[boolIsData]  <- FALSE
                                          if(dataOnly)                  validValues[!boolIsData] <- FALSE
                                      }
                                      
                                      if(!includePredictive)         validValues <- safeUpdateValidValues(validValues, idsVec_exclude = predictiveNodeIDs)
                                      if(predictiveOnly)             validValues <- safeUpdateValidValues(validValues, idsVec_only    = predictiveNodeIDs)
                                      if(topOnly)                       validValues <- safeUpdateValidValues(validValues, idsVec_only = modelDef$maps$top_IDs)
                                      if(latentOnly)                    validValues <- safeUpdateValidValues(validValues, idsVec_only = modelDef$maps$latent_IDs)
                                      if(endOnly)                       validValues <- safeUpdateValidValues(validValues, idsVec_only = modelDef$maps$end_IDs)

                                      ## Part of fix for Issue #340.
                                      ## In general the flow of model/node querying functions sometimes flips between IDs and names multiple times
                                      ## which is inefficienty/redudant.  I am adding some logic here to avoid such a flip when it would be
                                      ## redundant.  In the future it may make sense to push this logic to be internal to expandNodesNames and/or
                                      ## nodeName2GraphIDs, but I am leaving that for a future step.

                                      ## New logic, part of fix for Issue #340:
                                      ## If returnScalarComponents is FALSE, we should not need to call expandNodeNames (stoch and determ) and RHSonly,
                                      ## which flips to names and back to IDs.  Instead we can work directly with the IDs
                                      if(!returnScalarComponents) {
                                          ans <- expandNodeNamesFromGraphIDs(which(validValues),
                                                                             returnScalarComponents = returnScalarComponents,
                                                                             returnType = returnType)
                                      } else {
                                          ## nodeNames2graphID is called inside expandNodeNames
                                          ans <- expandNodeNames(modelDef$maps$graphID_2_nodeName[validValues],
                                                                 returnScalarComponents = returnScalarComponents,
                                                                 returnType = returnType)
                                      }
                                      return(ans)
                                  },
                                  safeUpdateValidValues = function(validValues, idsVec_only, idsVec_exclude) {
                                      if(!missing(idsVec_only) && !missing(idsVec_exclude)) stop()
                                      if(!missing(idsVec_only)) {
                                          if(length(idsVec_only) == 0) {
                                              validValues <- rep(FALSE, length(validValues))
                                          } else   validValues[-idsVec_only] <- FALSE
                                      }
                                      if(!missing(idsVec_exclude)) {
                                          if(length(idsVec_exclude) > 0)   validValues[idsVec_exclude] <- FALSE
                                      }
                                      return(validValues)
                                  },
expandNodeNamesFromGraphIDs = function(graphID, returnScalarComponents = FALSE, returnType = 'names', sort = FALSE) {
    if(length(graphID)==0) return(if(returnType=='names') character() else numeric())
    if(sort)
        graphID <- sort(graphID)
    if(returnType == 'names'){
        if(returnScalarComponents) nodeNames <- modelDef$maps$elementNames[graphID] ## these are really elementIDs
        else nodeNames <- modelDef$maps$graphID_2_nodeName[graphID]
        return(nodeNames)
    }
    if(returnType == 'ids'){
        if(returnScalarComponents) warning("NIMBLE development warning: returning IDs of scalar components may not be meaningful.  Checking to see if we ever see this message.")
        return(graphID)
    }
    if(!(returnType %in% c('ids','names')))
        stop('instead expandNodeNames, imporper returnType chosen')
},
expandNodeNames = function(nodes, env = parent.frame(), returnScalarComponents = FALSE, returnType = 'names', sort = FALSE, unique = TRUE){
                                      '
Takes a vector of names of nodes or variables and returns the unique and expanded names in the model, i.e. \'x\' expands to \'x[1]\', \'x[2]\', ...

Arguments:

nodes: a vector of names of nodes (or variables) to be expanded. Alternatively, can be a vector of integer graph IDs, but this use is intended only for advanced users

returnScalarComponents: should multivariate nodes (i.e. dmnorm or dmulti) be broken up into scalar components?

returnType: return type. Options are \'names\' (character vector) or \'ids\' (graph IDs)

sort: should names be topologically sorted before being returned?

unique: should names be the unique names or should original ordering of nodes (after expansion of any variable names into node names) be preserved
'

                                      if(length(nodes) == 0) return(if(returnType=='names') character() else numeric())
                                      graphID <- modelDef$nodeName2GraphIDs(nodes, !returnScalarComponents, unique = unique, ignoreNotFound = TRUE)
                                      expandNodeNamesFromGraphIDs(graphID, returnScalarComponents, returnType, sort)
                                  },

                                  topologicallySortNodes = function(nodes, returnType = 'names') {
                                      '
Sorts the input list of node names according to the topological dependence ordering of the model structure.

Arguments:

nodes: A character vector of node or variable names, which is to be topologically sorted. Alternatively can be a numeric vector of graphIDs

returnType: character vector indicating return type. Choices are "names" or "ids"

Details: This function merely reorders its input argument.  This may be important prior to calls such as model$simulate(nodes) or model$calculate(nodes), to enforce that the operation is performed in topological order.
'
                                      nodeIDs <- expandNodeNames(nodes, returnType = 'ids')
                                      nodeIDs <- sort(nodeIDs)
                                      nodeNames <- expandNodeNames(nodeIDs, returnType = returnType)
                                      return(nodeNames)
                                  },

                                  getVarInfo = function(name, includeLogProb = TRUE) {
                                      if(missing(name)) return(modelDef$varInfo)
                                      ans <- modelDef$varInfo[[name]]
                                      if(is.null(ans) & includeLogProb) ans <- modelDef$logProbVarInfo[[name]]
                                      return(ans)
                                  },
                                  getSymbolTable = function() modelDef$symTab,

                                  init_isDataEnv = function() {
                                      ## initializes the 'isDataEnv' to logical arrays of 'FALSE', based on dimensions in 'isDataVars' object
                                      list2env(lapply(isDataVars, nimble:::createDefault_isDataObj), isDataEnv)
                                      ## at the time of initializing the isDataEnv, also set predictiveNodeIDs and predictiveRootNodeIDs
                                      ##setPredictiveNodeIDs()
                                      if(!getNimbleOption('determinePredictiveNodesInModel')) {
                                          predictiveNodeIDs     <<- integer()
                                          predictiveRootNodeIDs <<- integer()
                                      } else {
                                          predictiveNodeIDs     <<- as.integer(getNodeNames(stochOnly = TRUE,                 returnType = 'ids'))
                                          predictiveRootNodeIDs <<- as.integer(getNodeNames(stochOnly = TRUE, topOnly = TRUE, returnType = 'ids'))
                                      }
                                  },

                                  resetData = function() {
'
Resets the \'data\' property of ALL model nodes to FALSE.  Subsequent to this call, the model will have no nodes flagged as \'data\'.
'
                                      ## re-initializes the 'isDataEnv', setting everything to 'FALSE'
                                      init_isDataEnv()
                                      return(invisible(NULL))
                                  },

                                  setData = function(..., warnAboutMissingNames = TRUE) {
'
Sets the \'data\' flag for specified stochastic nodes to TRUE, and also sets the value of these nodes to the value provided.  This is the exclusive method for specifying \'data\' nodes in a model object.  When a \'data\' argument is provided to \'nimbleModel()\', it uses this method to set the data nodes.

Arguments:

...:  Arguments may be provided as named elements with numeric values or as character names of model variables.  These may be provided in a single list, a single character vector, or as multiple arguments.  When a named element with a numeric value is provided, the size and dimension must match the corresponding model variable.  This value will be copied to the model variable and any non-NA elements will be marked as data.  When a character name is provided, the value of that variable in the model is not changed but any currently non-NA values are marked as data.  Examples: setData(\'x\', y = 1:10) will mark both x and y as data and will set the value of y to 1:10.  setData(list(\'x\', y = 1:10)) is equivalent.  setData(c(\'x\',\'y\')) or setData(\'x\',\'y\') will mark both x and y as data.

Details: If a provided value (or the current value in the model when only a name is specified) contains some NA values, then the model nodes corresponding to these NAs will not have their value set, and will not be designated as \'data\'.  Only model nodes corresponding to numeric values in the argument list elements will be designated as data.  Designating a deterministic model node as \'data\' will result in an error.  Designating part of a multivariate node as \'data\' and part as non-data (NA) is allowed, but \'isData()\' will report such a node as being \'data\', calculations with the node will generally return NA, and MCMC samplers will not be assigned to such nodes.
'
                                          ## new functionality for setData():
                                          ## ... can be a list, a character vector of variable names, or a mix of both
                                          ## intention is to flag these variables as 'data', and not change any model values.
                                          ## some inefficiency here (accesses model values, then re-sets the same model values),
                                          ## but this simplifies the addition without changing exisiting code.
                                           data = list(...)
                                           ## Check if a single list or character vector was provided
                                           if(length(data) == 0 || (length(data) == 1 && is.null(data[[1]]))) ## NULL catches case that occurs with test_size
                                               return()
                                           if(is.null(names(data))) {  ## if no names, should be either characters or single list
                                               if(length(data) == 1) {
                                                   if(is.character(data[[1]])) { # e.g., c('a','b','c')
                                                       data <- as.list(data[[1]])
                                                   } else {
                                                       if(is.list(data[[1]])) {
                                                           data <- data[[1]]
                                                       } else stop("setData: single input should be a list or character strings")
                                                   }
                                               } else if(!all(sapply(data, function(x) is.character(x) && length(x) == 1)))
                                                   stop("setData: multiple inputs must be named or be individual variable names")
                                           } ##  otherwise, treat each element of ... as variable
                                           dataNames <- names(data)
                                           if(is.null(dataNames)) dataNames <- rep("", length(data))
                                           for(i in seq_along(data)) {
                                               if(dataNames[i] == "" && is.character(data[[i]])) {
                                                   dataNames[i] <- data[[i]]
                                                   if(exists(dataNames[i]))
                                                       data[[i]] <- get(dataNames[i])
                                               }
                                           }
                                           names(data) <- dataNames
                                           if(any(names(data) == ""))
                                               warning("setData: one or more elements of 'data' is unnamed.")
                                      origData <<- data
                                      ## argument is a named list of data values.
                                      ## all nodes specified (except with NA) are set to that value, and have isDataEnv$VAR set to TRUE
                                      for(iData in seq_along(data)) {
                                          varName <- dataNames[iData]
                                          varValue <- data[[varName]]
                                          if(is.data.frame(varValue)) {
                                              if(!all(sapply(varValue, is.numeric)))
                                                  stop("setData: '", varName, "' must be numeric")
                                              varValue <- as.matrix(varValue)
                                              dimnames(varValue) <- NULL
                                          }
                                          if(!(varName %in% names(isDataVars))) {
                                              ## when data is from modelDef$constantsList,
                                              ## it is possible that the constants don't exist on LHS of BUGS decls
                                              ## and so are not variables in the model.  In that case we don't want to issue the warning.
                                              if(warnAboutMissingNames
                                                 && nimble::nimbleOptions("verbose")) {
                                                  if(varName == '') {
                                                      warning('setData: unnamed element provided to setData.')
                                                  } else 
                                                      messageIfVerbose("  [Note] '", varName, "' is provided in 'data' but is not a variable in the model and is being ignored.")
                                              }
                                              ## Removing unnecessary
                                              ## elements does not
                                              ## seem to be necessary
                                              ## for later processing
                                              ## steps, but we do it
                                              ## as a cleanup step.
                                              data[[varName]] <- NULL
                                              next
                                          }
                                          if(length(isDataVars[[varName]]))
                                              scalarize <- FALSE else scalarize <- TRUE  ## if non-scalar, check actual dimensionality of input
                                          if(length(nimble::nimbleInternalFunctions$dimOrLength(varValue, scalarize = scalarize)) != length(isDataVars[[varName]]))   stop(paste0('incorrect size or dim in data: ', varName))
                                          if(!(all(nimble::nimbleInternalFunctions$dimOrLength(varValue, scalarize = scalarize) == isDataVars[[varName]])))   stop(paste0('incorrect size or dim in data: ', varName))

                                          expandedNodeNames <- expandNodeNames(varName, returnScalarComponents = TRUE)
                                          determElements <- .self$isDeterm(expandedNodeNames, nodesAlreadyExpanded = TRUE)
                                          if(any(determElements))
                                              if(any(!is.na(varValue[which(determElements)])))
                                                  stop("setData: '", varName, "' contains deterministic nodes. Deterministic nodes cannot be specified as 'data' or 'constants'.")
                                          
                                          .self[[varName]] <- varValue
                                          ## Values set as NA are not flagged as data nor are RHSonly elements.
                                          isDataVarValue <- !is.na(varValue)
                                          isDataVarValue[!.self$isStoch(expandedNodeNames, nodesAlreadyExpanded = TRUE)] <- FALSE
                                          names(isDataVarValue) <- NULL
                                          assign(varName, isDataVarValue, envir = isDataEnv)
                                      }
                                   ##   testDataFlags()  ## this is slow for large models.  it could be re-written if we want to use it routinely
                                      setPredictiveNodeIDs()
                                      return(invisible(NULL))
                                  },

                                  setPredictiveNodeIDs = function() {
                                      if(!getNimbleOption('determinePredictiveNodesInModel')) {
                                          predictiveNodeIDs     <<- integer()
                                          predictiveRootNodeIDs <<- integer()
                                          return()
                                      }
                                      ## determine all current predictive and predictive root nodes in the model
                                      predNodeIDs <- integer()
                                      predRootNodeIDs <- integer()
                                      stochNonDataIDs <- getNodeNames(stochOnly = TRUE, includeData = FALSE, returnType = 'ids')
                                      anyPredictiveNodes <- any(isEndNode(stochNonDataIDs))
                                      if(anyPredictiveNodes) {
                                          ## for starters, all stochNonData nodes, which are end nodes, are predictive nodes
                                          ## (these don't get included, otherwise):
                                          predNodeIDs <- stochNonDataIDs[isEndNode(stochNonDataIDs)]
                                          ## now, find all potential (candidate) non-end predictive nodes:
                                          candidateNonEndPredNodeIDs <- stochNonDataIDs[!isEndNode(stochNonDataIDs)]
                                          dataNodeIDs <- getNodeNames(dataOnly = TRUE, returnType = 'ids')
                                          dataNodeParentIDs <- expandNodeNames(getParents(dataNodeIDs, stochOnly = TRUE), returnType = 'ids')
                                          ## remove from candidate non-end predictive nodes all direct parents of data nodes:
                                          candidateNonEndPredNodeIDs <- setdiff(candidateNonEndPredNodeIDs, dataNodeParentIDs)
                                          nCandidate <- length(candidateNonEndPredNodeIDs)
                                          nextCandInd <- 1
                                          while(nextCandInd <= nCandidate) {
                                              thisCandNodeID <- as.numeric(candidateNonEndPredNodeIDs[nextCandInd])
                                              stochDownstreamNoSelfIDs <- getDependencies(thisCandNodeID, self = FALSE, stochOnly = TRUE, downstream = TRUE, returnType = 'ids')
                                              ## skip candidate nodes that have any downstream data nodes:
                                              if(length(intersect(stochDownstreamNoSelfIDs, dataNodeIDs)) > 0)   { nextCandInd <- nextCandInd + 1;   next }
                                              ## found a predictive root node:
                                              predRootNodeIDs <- c(predRootNodeIDs, thisCandNodeID)
                                              ## everything downstream from (and including) this root node are predictive nodes:
                                              predNodeIDs <- c(predNodeIDs, thisCandNodeID, stochDownstreamNoSelfIDs)
                                              ## remove all downstream dependencies of this node from the candidate set,
                                              ## which are already marked as being predictive nodes:
                                              candidateNonEndPredNodeIDs <- candidateNonEndPredNodeIDs[-(1:nextCandInd)]
                                              candidateNonEndPredNodeIDs <- setdiff(candidateNonEndPredNodeIDs, stochDownstreamNoSelfIDs)
                                              nCandidate <- length(candidateNonEndPredNodeIDs)
                                              nextCandInd <- 1
                                          }
                                      }
                                      ## set into the model object's fields:
                                      predictiveNodeIDs     <<- as.integer(sort(unique(predNodeIDs)))       ## can contain duplicates
                                      predictiveRootNodeIDs <<- as.integer(sort(unique(predRootNodeIDs)))
                                  },

                                  getPredictiveNodeIDs     = function() return(predictiveNodeIDs),
                                  getPredictiveRootNodeIDs = function() return(predictiveRootNodeIDs),

                                  testDataFlags = function() {
                                      ## this function tests for *mixed* T/F flags in the isData flag of all nodes.
                                      ## Only really tests when there's a vectorized node declaration, which gives rise to >1 isData flag.
                                      lapply(as.list(getNodeNames()),
                                             function(nn) { ## check if there are *both* TRUE and FALSE flags
                                                 isDataVals <- as.vector(eval(parse(text=nn, keep.source = FALSE)[[1]], envir=isDataEnv))
                                                 if(!(all(isDataVals) || all(!isDataVals))) stop(paste0('it seems we have mixed isData flags for the vectorized node: ', nn)) })
                                      lapply(as.list(getNodeNames(determOnly = TRUE)),
                                             function(nn) { ## check if any determ nodes have isData flag TRUE
                                                 isDataVals <- as.vector(eval(parse(text=nn, keep.source = FALSE)[[1]], envir=isDataEnv))
                                                 if(isDataVals[1]) stop(paste0('it seems we have isData=TRUE for a deterministic node: ', nn)) })
                                      return(invisible(NULL))
                                  },

                                  isData = function(nodes) {
'
Returns a vector of logical TRUE / FALSE values, corresponding to the \'data\' flags of the input node names.

Arguments:

nodes: A character vector of node or variable names.

Details: The variable or node names specified is expanded into a vector of model node names. A logical vector is returned, indicating whether each model node has been flagged as containing \'data\'. Multivariate nodes for which any elements are flagged as containing \'data\' will be assigned a value of TRUE.
'
                                                g_id = modelDef$nodeName2GraphIDs(nodes, unique = FALSE)
                                  		return(isDataFromGraphID(g_id))
                                  },

                                  isDataFromGraphID = function(g_id){
                                      ## returns TRUE if any elements are flagged as data
                                      nodeNames <- modelDef$maps$graphID_2_nodeName[g_id]
                                  	ret <- unlist(lapply(as.list(nodeNames),
                                                  function(nn)
                                                    return(any(eval(parse(text=nn, keep.source = FALSE)[[1]],
                                                                                envir=isDataEnv)))))
                                    if(is.null(ret))   ret <- logical(0)
                                    return(ret)
                                  },

                                  getDependenciesList = function(returnNames = TRUE, sort = TRUE) {
'
Returns a list of all dependent neighbor relationships.  Each list element gives the one-step dependencies of one vertex, and the element name is the vertex label (integer ID or character node name)

Arguments:

returnNames: If TRUE (default), list names and element contents are returns as character node names, e.g. \'x[1]\'.  If FALSE, everything is returned using graph IDs, which are unique integer labels for each node.

sort: If TRUE (default), each list element is returned in topologically sorted order.  If FALSE, they are returned in arbitrary order.

Details: This provides a fairly raw representation of the graph (model) structure that may be useful for inspecting what NIMBLE has created from model code.
'
                                      if(!returnNames)
                                          if(!sort) return(modelDef$maps$edgesFrom2To)
                                          else return(lapply(modelDef$maps$edgesFrom2To, sort))
                                      else {
                                          if(!sort) ans <- lapply(modelDef$maps$edgesFrom2To, function(x) modelDef$maps$graphID_2_nodeName[x])
                                          else ans <- lapply(modelDef$maps$edgesFrom2To, function(x) modelDef$maps$graphID_2_nodeName[sort(x)])
                                          names(ans) <- modelDef$maps$graphID_2_nodeName[as.numeric(names(ans))]
                                          return(ans)
                                      }
                                  },

                                  getDependencyPathCountOneNode = function(node) {
                                      if(length(node) > 1)
                                          stop("getDependencyPathCountOneNode: argument 'node' should provide a single node.")
                                      if(is.character(node)) {
                                          node <- modelDef$nodeName2GraphIDs(node)
                                      }
                                      if(!is.numeric(node))
                                          stop("getDependencyPathCountOneNode: argument 'node' should be a character node name or a numeric node ID.")
                                      modelDef$maps$nimbleGraph$getDependencyPathCountOneNode(node = node)
                                  },
getDependencyPaths = function(node) {
    if(length(node) > 1)
        stop("getDependencyPaths: argument 'node' should provide a single node.")
    if(is.character(node)) {
        node <- modelDef$nodeName2GraphIDs(node)
    }
    if(!is.numeric(node))
        stop("getDependencyPaths: argument 'node' should be a character node name or a numeric node ID.")
    modelDef$maps$nimbleGraph$getDependencyPaths(node = node)
},
getParentsList = function(returnNames = TRUE, sort = TRUE) {
'
Returns a list of all parent neighbor relationships.  Each list element gives the one-step parents of one vertex, and the element name is the vertex label (integer ID or character node name)

Arguments:

returnNames: If TRUE (default), list names and element contents are returns as character node names, e.g. \'x[1]\'.  If FALSE, everything is returned using graph IDs, which are unique integer labels for each node.

sort: If TRUE (default), each list element is returned in topologically sorted order.  If FALSE, they are returned in arbitrary order.

Details: This provides a fairly raw representation of the graph (model) structure that may be useful for inspecting what NIMBLE has created from model code.
'
  maps <- modelDef$maps
  # Future option: We could put edgesTo2From in the maps.
  # For now we create it on the fly here
  maxNodeID <- length(maps$vertexID_2_nodeID) ## should be same as length(maps$nodeNames)
  edgesLevels <- if(maxNodeID > 0) 1:maxNodeID else numeric(0)
  fedgesTo <- factor(maps$edgesTo, levels = edgesLevels) ## setting levels ensures blanks inserted into the splits correctly
  edgesTo2From <- split(maps$edgesFrom, fedgesTo)

  if(!returnNames)
    if(!sort) return(edgesTo2From)
  else return(lapply(edgesTo2From, sort))
  else {
    if(!sort) ans <- lapply(edgesTo2From, function(x) modelDef$maps$graphID_2_nodeName[x])
    else ans <- lapply(edgesTo2From, function(x) modelDef$maps$graphID_2_nodeName[sort(x)])
    names(ans) <- modelDef$maps$graphID_2_nodeName[as.numeric(names(ans))]
    return(ans)
  }
},

getParents = function(nodes, omit = character(), self = FALSE,
                      determOnly = FALSE, stochOnly = FALSE,
                      includeData = TRUE, dataOnly = FALSE,
                      includeRHSonly = FALSE, upstream = FALSE,
                      immediateOnly = FALSE,
                      returnType = 'names', returnScalarComponents = FALSE) {
  '
 Returns a character vector of the nodes on which the input nodes depend, sorted topologically according to the model graph, by default recursing and stopping at stochastic parent nodes.  In the genealogical metaphor for a graphical model, this function returns the "parents" of the input nodes. In the river network metaphor, it returns upstream nodes.  By default, the returned nodes omit the input nodes. Additional input arguments provide flexibility in the values returned.

Arguments:

nodes: Character vector of node names, with index blocks allowed, and/or variable names, the parents of which will be returned.

omit: Character vector of node names, which will be omitted from the nodes returned.  In addition, parent nodes beyond these omitted nodes will not be returned.  The omitted nodes argument serves to stop the upward search through the hierarchical model structure, and excludes the specified node.

self: Logical argument specifying whether to include the input argument nodes in the return vector of dependent nodes.  Default is FALSE.

determOnly: Logical argument specifying whether to return only deterministic nodes.  Default is FALSE.

stochOnly: Logical argument specifying whether to return only stochastic nodes.  Default is FALSE.  If both determOnly and stochOnly are TRUE, no nodes will be returned.

includeData: Logical argument specifying whether to include \'data\' nodes (set via nimbleModel or the setData method).  Default is TRUE.

dataOnly: Logical argument specifying whether to return only \'data\' nodes.  Default is FALSE.

includeRHSonly: Logical argument specifying whether to include right-hand-side-only nodes (model nodes which never appear on the left-hand-side of ~ or <- in the model code).  These nodes are neither stochastic nor deterministic, but instead function as variable inputs to the model.  Default is FALSE.

upstream: Logical argument specifying whether the upward search through the hierarchical model structure should continue beyond the first and subsequent stochastic nodes encountered, hence returning all nodes upstream of the input nodes.  Default is FALSE.

immediateOnly: Logical argument specifying whether only the immediate parent nodes should be returned, even if they are deterministic.  If FALSE, getParents recurses and stops at stochastic nodes.  Default is FALSE.

returnType: Character argument specifying type of object returned. Options are \'names\' (returns character vector) and \'ids\' (returns numeric graph IDs for model).

returnScalarComponenets: Logical argument specifying whether multivariate nodes should be returned as full node names (i.e. \'x[1:2]\') or as scalar componenets (i.e. \'x[1]\' and \'x[2]\').

Details: The upward search for dependent nodes propagates through deterministic nodes, but by default will halt at the first level of stochastic nodes encountered.  Use getParentsList for a list of one-step parent nodes of each node in the model.
'
  if(is.character(nodes)) {
    elementIDs <- modelDef$nodeName2GraphIDs(nodes, FALSE)
    nodeIDs <- unique(modelDef$maps$elementID_2_vertexID[elementIDs],     ## turn into IDs in the graph
                      FALSE,
                      FALSE,
                      NA)
  }
  else if(is.numeric(nodes))
    nodeIDs <- nodes
  if(is.character(omit)) {
    elementIDs <- modelDef$nodeName2GraphIDs(omit, FALSE)
    omitIDs <- unique(modelDef$maps$elementID_2_vertexID[elementIDs],     ## turn into IDs in the graph
                      FALSE,
                      FALSE,
                      NA)
  }
  else if(is.numeric(omit))
    omitIDs <- omit

  parentIDs <- modelDef$maps$nimbleGraph$getParents(nodes = nodeIDs,
                                                    omit = if(is.null(omitIDs)) integer() else omitIDs,
                                                    upstream = upstream,
                                                    immediateOnly = immediateOnly)
  if(self)	{ # The C++ call does *not* return self nodes
    nodeFunIDs <- unique(modelDef$maps$vertexID_2_nodeID[ nodeIDs ])
    parentIDs <- sort(c(parentIDs, nodeFunIDs))
  }
  if(!includeRHSonly) parentIDs <- parentIDs[modelDef$maps$types[parentIDs] != 'RHSonly']
  if(determOnly)      parentIDs <- parentIDs[modelDef$maps$types[parentIDs] == 'determ']
  if(stochOnly)	      parentIDs <- parentIDs[modelDef$maps$types[parentIDs] == 'stoch']
  if(!includeData)	parentIDs <- parentIDs[!isDataFromGraphID(parentIDs)]
  if(dataOnly)		parentIDs <- parentIDs[isDataFromGraphID(parentIDs)]
  
  parentIDs <- modelDef$nodeName2GraphIDs(modelDef$maps$graphID_2_nodeName[parentIDs], !returnScalarComponents)
  if(returnScalarComponents)
    parentIDs = unique(parentIDs, FALSE, FALSE, NA)
  if(returnType == 'ids') {
    if(returnScalarComponents) warning("NIMBLE development warning: calling getParents with returnType = ids and returnScalarComponents may not be meaningful.")
    return(depIDs)
  }
  if(returnType == 'names') {
    if(returnScalarComponents)
      return(modelDef$maps$elementNames[parentIDs])
    retVal <- modelDef$maps$nodeNames[parentIDs]
    return(retVal)
  }
  if(!(returnType %in% c('ids', 'names')))
    stop('instead getDependencies, improper returnType chosen')
},
getConditionallyIndependentSets = function(nodes,
                                           givenNodes,
                                           omit = integer(),
                                           inputType = c("latent", "param", "data"),
                                           stochOnly = TRUE,
                                           returnType = 'names',
                                           returnScalarComponents = FALSE) {
  '
Get a list of conditionally independent sets of nodes in a nimble model.

Conditionally independent sets of nodes are typically groups of latent states whose joint conditional probability (density) will not change even if any other non-fixed node is changed.  Default fixed nodes are data nodes and parameter nodes (with no parent nodes), but this can be controlled.

model: A nimble model object (uncompiled or compiled).

nodes: A vector of node names or their graph IDs that are the starting nodes from which conditionally independent sets of nodes should be found.  If omitted, the default will be all latent nodes, defined as stochastic nodes that are not data and have at least one stochastic parent node (possible with deterministic nodes in between).  Note that this will omit latent states that have no hyperparameters.  An example is the first latent state in some state-space (time-series) models, which is sometimes declared with known prior.  See type because it relates to the interpretation of nodes.

givenNodes: A vector of node names or their graph IDs that should be considered as fixed and hence can be conditioned on.  If omitted, the default will be all data nodes and all parameter nodes, the latter defined as nodes with no stochastic parent nodes (skipping over deterministic parent nodes).

omit: A vector of node names or their graph IDs that should be omitted and should block further graph exploration. 

intputType: Type of input nodes provided in nodes argument.  For \'latent\', the input nodes are interpreted as latent states, from which both downstream and upstream exploration should be done to find nodes in the same set (nodes that are not conditionally independent from each other).  For \'param\', the input nodes are interpreted as parameters, so graph exploration begins from the  top (input) and proceeds downstream.  For \'data\', the input nodes are interpreted and data nodes, so graph exploration begins from the bottom (input) and proceeds upstream.

stochOnly: Logical for whether only stochastic nodes should be returned (default = TRUE).  If FALSE, both deterministic and stochastic nodes are returned.

returnType: Either \'names\' for returned nodes to be node names or \'ids\' for returned nodes to be graph IDs.

returnScalarComponents: If FALSE (default), multivariate nodes are returned as full names (e.g. \'x[1:3]\').  If TRUE, they are returned as scalar elements (e.g. \'x[1]\',  \'x[2]\',  \'x[3]\').

Details: This function returns sets of conditionally independent nodes.  Multiple input nodes might be in the same set or different sets, and other nodes (not in codes) will be included.

By default, deterministic dependencies of givenNodes are also counted as given nodes.  This is relevant only for parent nodes. This allows the givenNodes to include only stochastic nodes.  Say we have A -> B -> C -> D.  A and D are givenNodes.  C is a latent node.  B is a deterministic node.  By default, B is considered given.  Otherwise, other dependent networks of nodes that depend on B would be grouped in the same output set as C, but this is usually not what is wanted.  Any use of the resulting output must ensure that B is calculated when necessary, as usual with nimble\'s model-generic programming.  To turn off this feature, set nimbleOptions(groupDetermWithGivenInCondIndSets = FALSE).

There is a non-exported function `nimble:::testConditionallyIndependentSets(model, sets, initialize = TRUE)` that tests whether the conditional independence of sets is valid.  It should be the case that `nimble:::testConditionallyIndependentSets(model, getConditionallyIndependentSets(model), initialize = TRUE)` returns `TRUE`.

Return value: List of nodes that are in conditionally independent sets.  Within each set, nodes are returned in topologically sorted order.  The sets themselves are returned in topologically sorted order of their first nodes.
'
  inputType <- match.arg(inputType)
  nimble:::getConditionallyIndependentSets(.self, nodes, givenNodes, omit, inputType,
                                  stochOnly, returnType, returnScalarComponents)
},
                                  getDependencies = function(nodes, omit = character(), self = TRUE,
                                      determOnly = FALSE, stochOnly = FALSE,
                                      includeData = TRUE, dataOnly = FALSE,
                                      includePredictive = getNimbleOption('getDependenciesIncludesPredictiveNodes'), predictiveOnly = FALSE,
                                      includeRHSonly = FALSE, downstream = FALSE,
                                      returnType = 'names', returnScalarComponents = FALSE) {
'
Returns a character vector of the nodes dependent upon the input argument nodes, sorted topologically according to the model graph. In the genealogical metaphor for a graphical model, this function returns the "children" of the input nodes.  In the river network metaphor, it returns downstream nodes. By default, the returned nodes include the input nodes, include both deterministic and stochastic nodes, and stop at stochastic nodes. Additional input arguments provide flexibility in the values returned.

Arguments:

nodes: Character vector of node names, with index blocks allowed, and/or variable names, the dependents of which will be returned.

omit: Character vector of node names, which will be omitted from the nodes returned.  In addition, dependent nodes subsequent to these omitted nodes will not be returned.  The omitted nodes argument serves to stop the downward search within the hierarchical model structure, and excludes the specified node.

self: Logical argument specifying whether to include the input argument nodes in the return vector of dependent nodes.  Default is TRUE.

determOnly: Logical argument specifying whether to return only deterministic nodes.  Default is FALSE.

stochOnly: Logical argument specifying whether to return only stochastic nodes.  Default is FALSE.  If both determOnly and stochOnly are TRUE, no nodes will be returned.

includeData: Logical argument specifying whether to include \'data\' nodes (set via nimbleModel or the setData method).  Default is TRUE.

dataOnly: Logical argument specifying whether to return only \'data\' nodes.  Default is FALSE.

includePredictive: Logical argument specifying whether to include predictive nodes (stochastic nodes, which themselves are not data and have no downstream stochastic dependents which are data). Used primarily to exclude predictive node calculations when setting up MCMC samplers on model parameters. Default value is controlled by the NIMBLE system option `getDependenciesIncludesPredictiveNodes`, which itself has a default value of `TRUE`.

predictiveOnly: Logical argument specifying whether to return only predictive nodes (stochastic nodes, which themselves are not data and have no downstream stochastic dependents which are data).  Default is FALSE.

includeRHSonly: Logical argument specifying whether to include right-hand-side-only nodes (model nodes which never appear on the left-hand-side of ~ or <- in the model code).  These nodes are neither stochastic nor deterministic, but instead function as variable inputs to the model.  Default is FALSE.

downstream: Logical argument specifying whether the downward search through the hierarchical model structure should continue beyond the first and subsequent stochastic nodes encountered, hence returning all nodes downstream of the input nodes.  Default is FALSE.

returnType: Character argument specifying type of object returned. Options are \'names\' (returns character vector) and \'ids\' (returns numeric graph IDs for model).

returnScalarComponenets: Logical argument specifying whether multivariate nodes should be returned as full node names (i.e. \'x[1:2]\') or as scalar componenets (i.e. \'x[1]\' and \'x[2]\').

Details: The downward search for dependent nodes propagates through deterministic nodes, but by default will halt at the first level of stochastic nodes encountered.  Use getDependenciesList for a list of one-step dependent nodes of each node in the model.
'
                                      if(is.character(nodes)) {
                                          ## always start from scalar components
                                          elementIDs <- modelDef$nodeName2GraphIDs(nodes, FALSE)
                                          nodeIDs <- unique(modelDef$maps$elementID_2_vertexID[elementIDs],     ## turn into IDs in the graph
                                                            FALSE,
                                                            FALSE,
                                                            NA)
                                      }
                                      else if(is.numeric(nodes))
                                          nodeIDs <- nodes

                                      if(is.character(omit)) { ## mimic above
                                          elementIDs <- modelDef$nodeName2GraphIDs(omit, FALSE)
                                          omitIDs <- unique(modelDef$maps$elementID_2_vertexID[elementIDs],     ## turn into IDs in the graph
                                                            FALSE,
                                                            FALSE,
                                                            NA)
                                      }
                                      else if(is.numeric(omit))
                                          omitIDs <- omit
## Go into C++
depIDs <- modelDef$maps$nimbleGraph$getDependencies(nodes = nodeIDs, omit = if(is.null(omitIDs)) integer() else omitIDs, downstream = downstream)
## ## Uncomment these lines to catch discrepancies between the C++ and R systems.
## depIDs <- nimble:::gd_getDependencies_IDs(graph = getGraph(), maps = getMaps(all = TRUE), nodes = nodeIDs, omit = omitIDs, downstream = downstream)
## if(!identical(as.numeric(depIDsOld), as.numeric(depIDs))) {
##     cat('caught a discrepancy for depIDs')
##     browser()
## }
                                      if(!includeRHSonly) depIDs <- depIDs[modelDef$maps$types[depIDs] != 'RHSonly']
                                      if(determOnly)	depIDs <- depIDs[modelDef$maps$types[depIDs] == 'determ']
                                      if(stochOnly)	depIDs <- depIDs[modelDef$maps$types[depIDs] == 'stoch']
                                      if(!self) {
                                          nodeFunIDs <- unique(modelDef$maps$vertexID_2_nodeID[ nodeIDs ])
                                          depIDs <- setdiff(depIDs, nodeFunIDs)
                                      }
                                      if(!includeData)       depIDs <- depIDs[!isDataFromGraphID(depIDs)]
                                      if(dataOnly)           depIDs <- depIDs[isDataFromGraphID(depIDs)]
                                      if(!includePredictive)    depIDs <- setdiff(  depIDs, predictiveNodeIDs)
                                      if(predictiveOnly)        depIDs <- intersect(depIDs, predictiveNodeIDs)

                                      depIDs <- modelDef$nodeName2GraphIDs(modelDef$maps$graphID_2_nodeName[depIDs], !returnScalarComponents)
                                      if(returnScalarComponents)
                                          depIDs = unique(depIDs, FALSE, FALSE, NA)
                                      if(returnType == 'ids'){
                                          if(returnScalarComponents) warning("NIMBLE development warning: calling getDependencies with returnType = ids and returnScalarComponents may not be meaningful.")
                                          return(depIDs)
                                      }
                                      if(returnType == 'names') {
                                          if(returnScalarComponents)
                                              return(modelDef$maps$elementNames[depIDs])
                                          retVal <- modelDef$maps$nodeNames[depIDs]
                                          return(retVal)
                                      }
                                      if(!(returnType %in% c('ids', 'names')))
                                          stop('instead getDependencies, improper returnType chosen')
                                  },

                                  getDownstream = function(...) {
'
Identical to getDependencies(..., downstream = TRUE)

Details: See documentation for member method getDependencies.
'
                                      getDependencies(..., downstream = TRUE)
                                  },

                                  setInits = function(inits) {
'
Sets initial values (or more generally, any named list of value elements) into the model

Arguments:

inits: A named list.  The names of list elements must correspond to model variable names.  The elements of the list must be of class numeric, with size and dimension each matching the corresponding model variable.
'
                                      if(!is.list(inits) || (is.list(inits) && length(inits) > 0 && is.null(names(inits)))) {
                                          warning('inits argument should be a named list; not adding initial values to model.')
                                          return(invisible(NULL))
                                      }

                                      origInits <<- inits
                                      for(i in seq_along(inits)) {
                                          if(names(inits)[i] == '') {
                                              warning(paste0('setInits: element ', i, ' of inits list is not named; skipping this element'))
                                              next
                                          }
                                          if(!(names(inits)[i] %in% .self$getVarNames())) {
                                              messageIfVerbose("  [Note] '", names(inits)[i], "' has initial values but is not a variable in the model and is being ignored.")
                                              next
                                          }
                                          dataVals <- .self$isDataEnv[[names(inits)[[i]] ]]
                                          if(any(dataVals)) {
                                              .self[[names(inits)[i]]][!dataVals] <- inits[[i]][!dataVals]
                                              if(any(!is.na(inits[[i]][dataVals])))
                                                  messageIfVerbose("  [Note] Ignoring non-NA values in inits for data nodes: ", names(inits)[[i]], ".")
                                          } else {
                                              inputDim <- nimbleInternalFunctions$dimOrLength(inits[[i]])
                                              varInfo <- .self$modelDef$varInfo[[names(inits)[i]]]
                                              mismatch <- FALSE
                                              if(length(inputDim) == 1 && inputDim == 1) {  # scalar could be scalar or vector of length 1
                                                  if(!(varInfo$nDim == 0 || (varInfo$nDim > 0 && identical(varInfo$maxs, rep(1, varInfo$nDim)))))
                                                      mismatch <- TRUE
                                              } else {
                                                  if(length(inputDim) != varInfo$nDim || any(inputDim != varInfo$maxs))
                                                      mismatch <- TRUE
                                              }
                                              if(mismatch)
                                                  message("  [Warning] Incorrect size or dimension of initial value for '", names(inits)[i], "'.\n         Initial value will not be used in compiled model.")
                                              .self[[names(inits)[i]]] <- inits[[i]]
                                          }
                                      }
                                  },
                                  checkConjugacy = function(nodeVector, restrictLink = NULL) {
                                      '
Determines whether or not the input nodes appear in conjugate relationships

Arguments:

nodeVector: A character vector specifying one or more node or variable names.  If omitted, all stochastic non-data nodes are checked for conjugacy.

Details: The return value is a named list, with an element corresponding to each conjugate node.  The list names are the conjugate node names, and list elements are the control list arguments required by the corresponding MCMC conjugate sampler functions.  If no model nodes are conjugate, an empty list is returned.
'
                                      if(missing(nodeVector)) nodeVector <- getNodeNames(stochOnly=TRUE, includeData=FALSE)
                                      nodeIDs <- expandNodeNames(nodeVector, returnType = 'ids')
                                      nimble:::conjugacyRelationshipsObject$checkConjugacy(.self, nodeIDs, restrictLink = restrictLink)
                                  },
                                  checkBasics = function() {
                                      '
Checks for size/dimension mismatches and for presence of NAs in model variables (the latter is not an error but a note of this is given to the user)
'
                                      # first do size checking; note that LHS of deterministic expressions are not necessarily filled in

                                      for(j in seq_along(.self$modelDef$declInfo)) {
                                              declInfo <- .self$modelDef$declInfo[[j]]
                                              nn <- length(declInfo$nodeFunctionNames)
                                              nfn <- declInfo$nodeFunctionNames[nn]
                                              nf <- .self$nodeFunctions[[j]]

                                              if(declInfo$type == 'determ') {
                                                  # check LHS and RHS are same size/dim
                                                  # need to eval within nf; constants not present otherwise
                                                  RHSsize <- try(nimble::nimbleInternalFunctions$dimOrLength(eval(codeSubstitute(declInfo$valueExprReplaced, as.list(nf)))), silent = TRUE)

                                                  LHSsize <- try(nimble::nimbleInternalFunctions$dimOrLength(eval(codeSubstitute(declInfo$targetExprReplaced, as.list(nf)))), silent = TRUE)
                                                  # apparently implicit dropping of size 1 dimensions is ok in determ node calcs
                                                  if(!is(RHSsize, 'try-error') && !is(LHSsize, 'try-error')) {
                                                      if(length(RHSsize) > 1 && any(RHSsize == 1))
                                                          RHSsize <- RHSsize[RHSsize != 1]
                                                      if(length(LHSsize) > 1 && any(LHSsize == 1))
                                                          LHSsize <- LHSsize[LHSsize != 1]

                                                      if(!identical(LHSsize, RHSsize))
                                                          stop("Size/dimension mismatch between left-hand side and right-hand size of BUGS expression: ", nimble:::safeDeparse(declInfo$code))
                                                  }
                                                  ## these lines caused a problem for functions such as chol() in BUGS code
                                                  ## removed by DT April 2016
                                                  ##if(is(RHSsize, 'try-error'))
                                                  ##    stop("Problem evaluating: ", deparse(declInfo$valueExprReplaced))
                                                  ##if(is(LHSsize, 'try-error'))
                                                  ##    stop("Problem evaluating: ", deparse(declInfo$targetExprReplaced))
                                              } else {
                                                  # check:
                                                  #   1) dims of param args match those in distInputList based on calculation
                                                  #   2) dims of param args match those in distInputList based on varInfo
                                                  #   3) sizes of vecs and row/column sizes all match for non-scalar quantities (only for Nimble-provided distributions)
                                                  dist <- nimble:::safeDeparse(declInfo$valueExprReplaced[[1]], warn = TRUE)

							# nimble:::getDimension so uses function not model method
                                                  distDims <- nimble::getDimension(dist, includeParams = TRUE)
                                                  nms <- names(distDims)
                                                  distDims <- as.integer(distDims); names(distDims) <- nms

                                                  sizes <- list(); length(sizes) <- length(nms); names(sizes) <- nms

                                                  for(k in seq_along(nms)) {
                                        # sometimes get_foo not found in env of nf (and doesn't appear in ls(nf) )
                                                      ##fun <- as.call(parse(text = paste0("nf$get_", nms[k])))
                                                      ##e = try(eval(fun))
                                                      ## NEWNODEFXN
                                                      e <- try(.self$getParam(nfn, nms[k]))

                                                      if(!is(e, "try-error")) {
                                                          if(!is.null(e)) {
                                                              sizes[[nms[k]]] <- nimble::nimbleInternalFunctions$dimOrLength(e)
                                                              if(prod(sizes[[nms[[k]]]]) == 1) sizes[[nms[[k]]]] <- numeric()
                                                          } else sizes[[nms[[k]]]] <- NA # when have param with dim > 2
                                                      } else messageIfVerbose("  [Warning] Unable to calculate parameter '", nms[k], "' for distribution ", dist, " for node '", nfn, "'; this may simply reflect that there are missing values in model variables.")
                                                  }
                                        # check dimensions based on varInfo
                                                  if(length(declInfo$targetExprReplaced) > 1) {
                                                      LHSvar <- nimble:::safeDeparse(declInfo$targetExprReplaced[[2]], warn = TRUE)
                                                  } else LHSvar <- nimble:::safeDeparse(declInfo$targetExprReplaced, warn = TRUE)
                                                  if(.self$modelDef$varInfo[[LHSvar]]$nDim < distDims['value'])
                                                      stop("Dimension of '", LHSvar, "' does not match required dimension for the distribution '", dist, "'. Necessary dimension is ", distDims['value'], ".", ifelse(distDims['value'] > 0, paste0(" You may need to include explicit indexing information, e.g., variable_name", ifelse(distDims['value'] < 2, "[1:2].", "[1:2,1:2].")), ''))
                                                  nms2 <- nms[nms%in%names(declInfo$valueExprReplaced)]
                                                  for(k in seq_along(nms2)) {
                                                      if(!is.numeric(declInfo$valueExprReplaced[[nms2[k]]]) &&
                                                         !(dist == 'dinterval' && nms2[k] == 'c') &&
                                                         !(dist == 'dcat' && nms2[k] == 'prob') &&
                                                         ( length(declInfo$valueExprReplaced[[nms2[k]]]) ==1
                                                             || nimble:::safeDeparse(declInfo$valueExprReplaced[[nms2[k]]][[1]], warn = TRUE) == '[' )) {  # can only check variables not expressions or constants
                                                          ## also dinterval can take 'c' param as scalar or vector
                                                          ## and dcat same for 'prob', so don't check
                                                          if(length(declInfo$valueExprReplaced[[nms2[k]]]) > 1) {
                                                              var <- nimble:::safeDeparse(declInfo$valueExprReplaced[[nms2[k]]][[2]], warn = TRUE)
                                                          } else var <- nimble:::safeDeparse(declInfo$valueExprReplaced[[nms2[k]]], warn = TRUE)
                                                          if(var %in% names(.self$modelDef$varInfo) && .self$modelDef$varInfo[[var]]$nDim < distDims[nms2[k]])
                                                              # check less than because variable dim can be bigger than node dim
                                                              stop("Dimension of '", nms2[k], "' does not match required dimension for the distribution '", dist, "'. Necessary dimension is ", distDims[nms2[k]], ".", ifelse(distDims[nms2[k]] > 0, paste0(" You may need to include explicit indexing information, e.g., variable_name", ifelse(distDims[nms2[k]] < 2, "[1:2].", "[1:2,1:2].")), ''))
                                                      }
                                                  }

                                        # actual dimensions
                                                  dims <- sapply(sizes, length)
                                                  toCheck <- names(dims[!is.na(sizes) & sapply(sizes, function(x) !is.null(x))])
                                                  if(dist == 'dinterval') toCheck <- toCheck[toCheck != 'c']
                                                  if(dist == 'dcat') {  # either scalar or vector is allowed (NCT issue 1251)
                                                      wh <- which(toCheck == 'prob')
                                                      if(dims[wh] > 1) 
                                                          stop("Dimension of distribution argument 'prob' does not match required dimension for the distribution 'dcat'. Necessary dimension is one (zero is also allowed).")
                                                      toCheck <- toCheck[toCheck != 'prob']
                                                  }
                                        # check dimensions based on empirical size of variables
                                                  if(!identical(dims[toCheck], distDims[toCheck])) {
                                                      mismatches <- which(dims[toCheck] != distDims[toCheck])
                                                      valueMismatch <- which('value' %in% names(mismatches))
                                                      if(length(valueMismatch))  ## Catch this first and refer to node name rather than 'value' (NCT issue 397)
                                                          stop("Dimension of '", nimble:::safeDeparse(declInfo$targetExpr),
                                                               "' does not match required dimension for the distribution '", dist,
                                                               "'. Necessary dimension is ", paste(distDims[toCheck][valueMismatch], collapse = ","), ".")
                                                      stop("Dimension of distribution argument(s) '", paste(names(mismatches), collapse = ","),
                                                           "' does not match required dimension(s) for the distribution '", dist, "'. Necessary dimension(s) are ",
                                                           paste(distDims[toCheck][mismatches], collapse = ","), ".",
                                                           ifelse(any(distDims[toCheck][mismatches] == 1),
                                                                  " You may need to ensure that you have explicit vectors and not one-row or one-column matrices.", ""))
                                                  }

                                        # check sizes
                                                  ## can exempt distributions from check that all non-scalar parameters have
                                                  ## the same length or size in each dimension (e.g., for dcar_normal)
                                                  if(!nimble:::isMixedSizes(dist)) {
                                                      mats <- dims == 2
                                                      vecs <- dims == 1
                                                      matRows <- unlist(sapply(sizes[mats], `[`, 1))
                                                      matCols <- unlist(sapply(sizes[mats], `[`, 2))
                                                      if(!length(unique(c(matRows, matCols, unlist(sizes[vecs])))) <= 1)
                                                          if(dist %in% names(nimble:::distributionsInputList)) {
                                                              stop("Size/dimension mismatch amongst vectors and matrices in BUGS expression: ", nimble:::safeDeparse(declInfo$code))
                                                          } else {
                                                              messageIfVerbose("  [Warning] Possible size/dimension mismatch amongst vectors and matrices in BUGS expression: '", nimble:::safeDeparse(declInfo$code), "'. Ignore this warning if the user-provided distribution has multivariate parameters with distinct sizes or if size of variable differs from sizes of parameters.")                                                                                                                                   }
                                                  }

                                              }
                                      }
                                    
                                      if(isTRUE(nimble::nimbleOptions('verbose'))){
                                        varsWithNAs <- NULL
                                        for(v in .self$getVarNames()){
                                          if(!nimble:::isValid(.self[[v]])){
                                            messageIfVerbose('  [Note] This model is not fully initialized. This is not an error.\n         To see which variables are not initialized, use model$initializeInfo().\n         For more information on model initialization, see help(modelInitialization).')
                                            break
                                          }
                                        }
                                      }

                                  },
                                  initializeInfo = function(stochasticLogProbs = FALSE) {
                                    '
Provides more detailed information on which model nodes are not initialized.

Arguments:

stochasticLogProbs: Boolean argument. If TRUE, the log-density value associated with each stochastic model variable is calculated and printed.
'
                                    if(stochasticLogProbs) {
                                        stochVars <- unique(nimble:::removeIndexing(.self$getNodeNames(stochOnly = TRUE)))
                                        for(v in stochVars)   cat(paste0(v, ': ', .self$calculate(v), '\n'))
                                        return(invisible(NULL))
                                    }
                                    varsWithNAs <- NULL
                                    for(v in .self$getVarNames()){
                                      if(!nimble:::isValid(.self[[v]])){
                                        varsWithNAs <- c(varsWithNAs, v)
                                      }
                                    }
                                    if(!is.null(varsWithNAs)){
                                      message('  [Note] Missing values (NAs) or non-finite values were found in model variables: ', paste(varsWithNAs, collapse = ', '), 
                                              '.\n  [Note] This is not an error, but some or all variables may need to be initialized for certain algorithms to operate properly.\n  [Note] For more information on model initialization, see help(modelInitialization).')
                                    }
                                    else{
                                      message('  [Note] All model variables are initialized.')
                                    }
                                  },
                                  check = function() {
                                      '
Checks for errors in model specification and for missing values that prevent use of calculate/simulate on any nodes
'
                                      # check for missing values and inability to calculate/simulate
                                      lp <- try(calculate())
                                      if(!nimble:::isValid(lp)) {
                                          varsToCheck <- character()
                                          for(v in .self$getVarNames())
                                              if(!nimble:::isValid(.self[[v]]) || !nimble:::isValid(getLogProb(setdiff(expandNodeNames(v), modelDef$maps$nodeNamesRHSonly))))
                                                  varsToCheck <- c(varsToCheck, v)
                                          badVars <- list(na=character(), nan=character(), inf=character())
                                          nns <- expandNodeNames(varsToCheck)
                                          nns <- topologicallySortNodes(nns)   ## should be unnecessary; just in case
                                          for(nn in nns) {
                                              val <- .self[[nn]]
                                              type <- getNodeType(nn)
                                              if(length(type) > 1) stop('check: Unexpectedly found "type" to have multiple values.')
                                              if(type == 'RHSonly') {
                                                  if(!nimble:::isValid(val)) badVars[[nimble:::whyInvalid(val)]] <- c(badVars[[nimble:::whyInvalid(val)]], nn)
                                              } else if(type == 'determ') {
                                                  test <- try(calculate(nn))
                                                  if(inherits(test, 'try-error'))
                                                      messageIfVerbose("  [Note] Cannot calculate logProb for node ", nn, ".")
                                                  val <- .self[[nn]]
                                                  if(!nimble:::isValid(val)) badVars[[nimble:::whyInvalid(val)]] <- c(badVars[[nimble:::whyInvalid(val)]], nn)
                                              } else if(type == 'stoch') {
                                                  if(!nimble:::isValid(val)) badVars[[nimble:::whyInvalid(val)]] <- c(badVars[[nimble:::whyInvalid(val)]], nn)
                                                  test <- try(val <- calculate(nn))
                                                  if(inherits(test, 'try-error'))
                                                      cat("  [Note] Cannot calculate logProb for node ", nn, ".")

                                                  if(!nimble:::isValid(val)) badVars[[nimble:::whyInvalid(val)]] <- c(badVars[[nimble:::whyInvalid(val)]], paste0('logProb_', nn))
                                              } else stop('Unknown node type: ', type)
                                          }
                                          badVars <- lapply(badVars, nimble:::removeIndexing)
                                          badVars <- lapply(badVars, unique)
                                          badVars <- lapply(badVars, function(nns) if(length(nns>0)) paste0(nns, collapse=', '))
                                          conds <- list(c('na','NAs'), c('nan','NaNs'), c('inf','Infinite values'))
                                          for(i in seq_along(conds)) {
                                              v <- badVars[[conds[[i]][1]]]
                                              m <- conds[[i]][2]
                                              if(!is.null(v)) messageIfVerbose(  "[Note] ", m, ' were detected in model variable', if(grepl(',',v)) 's' else '', ': ', v, ".")
                                          }
                                      }

                                  },

newModel = function(data = NULL, inits = NULL, modelName = character(), replicate = FALSE, check = getNimbleOption('checkModel'), calculate = TRUE) {
                                      '
Returns a new R model object, with the same model definiton (as defined from the original model code) as the existing model object.

Arguments:

data: A named list specifying data nodes and values, for use in the newly returned model.  If not provided, the data argument from the creation of the original R model object will be used.

inits: A named list specifying initial valuse, for use in the newly returned model.  If not provided, the inits argument from the creation of the original R model object will be used.

modelName: An optional character string, used to set the internal name of the model object.  If provided, this name will propagate throughout the generated C++ code, serving to improve readability.

replicate: Logical specifying whether to replicate all current values and data flags from the current model in the new model.  If TRUE, then the data and inits arguments are not used.  Default value is FALSE.

check: A logical indicating whether to check the model object for missing or invalid values.  Default is given by the NIMBLE option \'checkModel\', see help on \'nimbleOptions\' for details. 

calculate: A logical indicating whether to run \'calculate\' on the model; this will calculate all deterministic nodes and logProbability values given the current state of all nodes. Default is TRUE. For large models, one might want to disable this, but note that deterministic nodes, including nodes introduced into the model by NIMBLE, may be NA. 

Details: The newly created model object will be identical to the original model in terms of structure and functionality, but entirely distinct in terms of the internal values.
'
                                      if(replicate) {
                                          newlyCreatedModel <- modelDef$newModel(check = FALSE, calculate = FALSE)
                                          nimCopy(from = .self, to = newlyCreatedModel, logProb = TRUE)
                                          for(var in ls(isDataEnv)) newlyCreatedModel$isDataEnv[[var]] <- isDataEnv[[var]]
                                          if(check) newlyCreatedModel$check()
                                          if(calculate) newlyCreatedModel$calculate()
                                          return(newlyCreatedModel)
                                      }
                                      if(is.null(data)) data <- if( inherits(origData, 'uninitializedField') ) list() else origData
                                      if(is.null(inits)) inits <- if( inherits(origInits, 'uninitializedField') ) list() else origInits
                                      modelDef$newModel(data = data, inits = inits, modelName = modelName, check = check, calculate = calculate)
                                  }
                              )
                              )



setMethod('[[', 'modelBaseClass',
          function(x, i) {
              if(length(i) != 1) stop(paste0("Only one node can be accessed from a model using '[['."), call. = FALSE)
              if(!is.indexed(i)) {
                  eval(substitute(x$VAR, list(VAR=i)))
              } else {
                  parsedNode <- parse(text=i, keep.source = FALSE)[[1]]
                  parsedNode[[2]] <- substitute(x$VAR, list(VAR=parsedNode[[2]]))
                  eval(parsedNode)
              }
          }
)

setMethod('[[<-', 'modelBaseClass',
          function(x, i, value) {
              if(length(i) != 1) stop(paste0("Only one node can be accessed from a model using '[['."), call. = FALSE)
              if(!is.indexed(i)) {
                  eval(substitute(x$VAR <- value, list(VAR=i)))
              } else {
                  parsedNode <- parse(text=i, keep.source = FALSE)[[1]]
                  parsedNode[[2]] <- substitute(x$VAR, list(VAR=parsedNode[[2]]))
                  assignmentExpr <- substitute(TARGET <- value, list(TARGET=parsedNode))
                  eval(assignmentExpr)
              }
              return(x)
          }
)

insertSingleIndexBrackets <- function(code, varInfo) {
    cLength <- length(code)
    if(cLength == 1) {
        if(is.name(code)) {
            varName <- as.character(code)
            thisVarInfo <- varInfo[[varName]]
            if(!is.null(thisVarInfo)) {
                if(thisVarInfo$nDim == 0)
                    return(substitute(VAR[1], list(VAR = code)))
            }
        }
        return(code)
    }
    if(is.call(code)) {
        recurseIndices <- 2:cLength
        if(code[[1]] == '[') {
            if(length(code[[2]]) == 1) {
                if(is.name(code[[2]])) recurseIndices <- 3:cLength
            }
        }
        for(i in recurseIndices) {
            code[[i]] <- insertSingleIndexBrackets(code[[i]], varInfo)
        }
        return(code)
    }
    if(!is.null(code)) warning('Unexpectedly reached end of insertSingleBrackets with ', safeDeparse(code))
    return(code)
}

# for now export this as R<3.1.2 give warnings if don't

#' Class \code{RmodelBaseClass}
#' @aliases RmodelBaseClass
#' @export
#' @description
#' Classes used internally in NIMBLE and not expected to be called directly by users.
RmodelBaseClass <- setRefClass("RmodelBaseClass",
                               contains = "modelBaseClass",
                               fields = list(
                                   nodeFunctions = 'ANY',	#list
                                   nodeFunctionGeneratorNames = 'ANY', #character, for efficiency in nimbleProject$addNimbleFunctionMulti
                                   nodeGenerators = 'ANY',	#list
                                   Cname = 'ANY',		#character
                                   CobjectInterface = 'ANY'
                                   ),
                               methods = list(
                                   initialize = function( ...) {
                                       callSuper(...)
                                   },
                                   setupDefaultMV = function(where = NULL) {
                                       defaultModelValues <<- modelDef$modelValuesClass(1)
                                       nimble:::pointAt(.self, defaultModelValues, index = 1)
                                   },

                                   buildNodeFunctions = function(where = globalenv(), debug = FALSE) {
                                       if(debug) browser()
                                       iNextNodeFunction <- 1
                                       numDecls <- length(modelDef$declInfo)
                                       nodeFunctions <<- vector('list', length = numDecls)  ## for the specialized instances
                                       nodeFunctionGeneratorNames <<- character(numDecls)
                                       nodeGenerators <<- vector('list', length = numDecls) ## for the nimbleFunctions
                                       for(i in seq_along(modelDef$declInfo)) {
                                           BUGSdecl <- modelDef$declInfo[[i]]
                                           if(BUGSdecl$numUnrolledNodes == 0) next
                                           ## extract needed pieces
                                           type <- BUGSdecl$type
                                           code <- BUGSdecl$codeReplaced
                                           code <- nimble:::insertSingleIndexBrackets(code, modelDef$varInfo)
                                           LHS <- code[[2]]
                                           RHS <- code[[3]]
                                           if(nimble::nimbleOptions('experimentalEnableDerivs')){
                                             parents <- BUGSdecl$allParentVarNames()
                                             selfWithNoInds <-  strsplit(nimble:::safeDeparse(LHS, warn = TRUE), '[', fixed = TRUE)[[1]][1]
                                             parents <- c(selfWithNoInds, parents)
                                             parentsSizeAndDims <- nimble:::makeSizeAndDimList(LHS, parents, BUGSdecl$unrolledIndicesMatrix, checkRagged = TRUE)
                                             parentsSizeAndDims <- nimble:::makeSizeAndDimList(RHS, parents, BUGSdecl$unrolledIndicesMatrix,
                                                                                               allSizeAndDimList = parentsSizeAndDims, checkRagged = TRUE)
                                           } else parentsSizeAndDims <- list()

                                           if(nimble::nimbleOptions()$allowDynamicIndexing && length(BUGSdecl$dynamicIndexInfo)) {  ## need dim for node for generating NaN with invalid dynamic indexes
                                               nodeSizeAndDims <- nimble:::makeSizeAndDimList(LHS, nimble:::safeDeparse(BUGSdecl$targetVarExpr, warn = TRUE),
                                                                                              BUGSdecl$unrolledIndicesMatrix,
                                                                                              checkRagged = FALSE)
                                               nodeDim <- nodeSizeAndDims[[nimble:::safeDeparse(BUGSdecl$targetVarExpr, warn = TRUE)]][[1]]$lengths
                                               nodeDim <- nodeDim[nodeDim != 1] ## will be NULL if scalar
                                               if(!length(nodeDim)) nodeDim <- NULL
                                           } else nodeDim <- NULL

                                           altParams <- BUGSdecl$altParamExprs
                                           altParams <- lapply(altParams, nimble:::insertSingleIndexBrackets, modelDef$varInfo)
                                           bounds <- BUGSdecl$boundExprs
                                           bounds <- lapply(bounds, nimble:::insertSingleIndexBrackets, modelDef$varInfo)
                                           logProbNodeExpr <- BUGSdecl$logProbNodeExpr
                                           logProbNodeExpr <- nimble:::insertSingleIndexBrackets(logProbNodeExpr, modelDef$logProbVarInfo)
                                           setupOutputExprs <- BUGSdecl$replacementNameExprs
                                           ## ensure they are in the same order as the columns of the unrolledIndicesMatrix, because that is assumed in nodeFunctionNew
                                           ## This can be necessary in a case like for(j in ...) for(i in ...) x[i,j] ~ ...; because x uses inner index first

                                           dynamicIndexInfo <- BUGSdecl$dynamicIndexInfo
                                           if(nimble::nimbleOptions()$allowDynamicIndexing) {
                                               for(iI in seq_along(dynamicIndexInfo))
                                                   dynamicIndexInfo[[iI]]$indexCode <- nimble:::insertSingleIndexBrackets(dynamicIndexInfo[[iI]]$indexCode, modelDef$varInfo)
                                           }

                                           if(nrow(BUGSdecl$unrolledIndicesMatrix) > 0)
                                               setupOutputExprs <- setupOutputExprs[ colnames(BUGSdecl$unrolledIndicesMatrix) ]
                                           ## make a unique name
                                           thisNodeGeneratorName <- paste0(nimble:::Rname2CppName(BUGSdecl$targetVarName), '_L', BUGSdecl$sourceLineNumber, '_', nimble:::nimbleUniqueID())
                                           ## create the nimbleFunction generator (i.e. unspecialized nimbleFunction)

                                           nfGenerator <- nimble:::nodeFunctionNew(LHS=LHS, RHS=RHS, name = thisNodeGeneratorName, altParams=altParams, bounds=bounds, parentsSizeAndDims = parentsSizeAndDims, logProbNodeExpr=logProbNodeExpr, type=type, setupOutputExprs=setupOutputExprs, dynamicIndexInfo = dynamicIndexInfo, nodeDim = nodeDim, evaluate=TRUE, where = where)
                                           nodeGenerators[[i]] <<- nfGenerator
                                           names(nodeGenerators)[i] <<- thisNodeGeneratorName
                                           nodeFunctionGeneratorNames[i] <<- thisNodeGeneratorName
                                           nodeFunctions[[i]] <<- nfGenerator(.self, BUGSdecl)
                                           names(nodeFunctions)[i] <<- thisNodeGeneratorName ## not sure what we need here
                                       }
                                   },

                                    buildNodesList = function() {   ## DANGEROUS!!  CAUSES R Studio TO CRASH!!  Unless the option NOT to try to inspect objects is used.
                                        nodes <<- list2env(nodeFunctions)			#trying to speed things up
                                    #    nodes <<- lapply(nodes, function(nf) getFunctionEnvVar(nf, 'nfRefClassObject'))
                                        return(NULL)
                                    },
                                   show = function() {
                                       cat(paste0('Rmodel object with     name: \'', name,    '\'\n'))
                                   }
                               )
)

RMakeCustomModelClass <- function(vars, className, isDataVars, modelDef, where = globalenv()) {
    newvars <- vars
    varnames <- if(is.list(vars)) names(vars) else vars

    inputList <- list(vars = newvars,
                      modelDef = modelDef,
                      isDataVars = isDataVars)

    ## uncomment this line if we want to ensure that every model refClass we generate is uniquely named internally
    className <- paste0(className, '_', nimbleUniqueID())

    eval(substitute(newClass <- setRefClass(
        Class = className,
        contains = 'RmodelBaseClass',
        fields = FIELDS,
        methods = list(
            initialize = function(inputList, ...) {
                nodes <<- new.env()		# list()
                classEnvironment <<- new.env()
                isDataEnv <<- new.env()
                nodeFunctions <<- list()
                nodeGenerators <<- list()
                vars <<- inputList$vars
                isDataVars <<- inputList$isDataVars
                callSuper(modelDef = inputList$modelDef, ...)
                setupDefaultMV()
                init_isDataEnv()
            }
        ), where = where),
                    list(FIELDS = makeBUGSclassFields(varnames, vars)
                         )))
    ans <- function(name = character()) {
        newClass(inputList = inputList, name = name)
    }
    ans
}

MakeCustomModelClass <- function(vars, className, where = globalenv())
    RMakeCustomModelClass(vars, className, where = where)      ## So demos work...

## This builds the list of all fields needed for a reference class definition.
## It is built as an unevaluated list parse tree so that when the class if created the function definitions are evaluated then.
makeBUGSclassFields <- function(vars, varDims) {
    activeBindingDefs <- list()
    envDefs <- as.list(rep('ANY', length(vars)))
    names(envDefs) <- makeEnvName(vars)
    rowDefs <-as.list(rep('ANY', length(vars)))
    nameDefs <- as.list(rep('ANY', length(vars)))
    names(rowDefs) <- makeRowName(vars)
    names(nameDefs) <- makeNameName(vars)
    for(var in vars) {
        activeBindingDefs[[var]] <- makeBUGSactiveBindingDef(makeEnvName(var), makeNameName(var), makeRowName(var), varDims[[var]])
    }
    as.call(c(as.name("list"), activeBindingDefs, envDefs, nameDefs, rowDefs))
}

## This uses the activeBindingTemplate and plugs in the 3 needed names
makeBUGSactiveBindingDef <- function(envVarName, varVarName, rowVarName, dims) {
    if(length(dims) == 0) dims <- 1
    if(prod(dims) == 1) {
        if(length(dims) > 1) {
            template <- activeBindingTemplateLength1NonScalar
        } else {
            template <- activeBindingTemplateLength1Vector
        }
    } else
        template <- activeBindingTemplate

    eval( substitute( substitute(aBT, list(ENVNAME = as.name(envVarName), VARNAME = as.name(varVarName), ROWNAME = as.name(rowVarName), DIMNAME = dims)), list(aBT = template) ) )
}
##e.g.  makeBUGSactiveBindingDef('.env_x','.name_x','.row_x')

## Parse tree template for the active binding functions
activeBindingTemplateLength1NonScalar <- quote( function(value) {
    if(missing(value)) return(if(is.na(ROWNAME)) ENVNAME[[VARNAME]] else ENVNAME[[VARNAME]][[ROWNAME]]) ## commas will get inserted after ROWNAME
    else {
        value <- array(value, dim = DIMNAME)
        if(is.na(ROWNAME)) ENVNAME[[VARNAME]] <- value
        else ENVNAME[[VARNAME]][[ROWNAME]] <- value
        return(invisible(value))
    }
})

activeBindingTemplateLength1Vector <- quote( function(value) {
    if(missing(value)) return(if(is.na(ROWNAME)) ENVNAME[[VARNAME]] else ENVNAME[[VARNAME]][[ROWNAME]]) ## commas will get inserted after ROWNAME
    else {
        value <- value[1]
        if(is.na(ROWNAME)) ENVNAME[[VARNAME]] <- value
        else ENVNAME[[VARNAME]][[ROWNAME]] <- value
        return(invisible(value))
    }
})

activeBindingTemplate <- quote( function(value) {
    if(missing(value)) return(if(is.na(ROWNAME)) ENVNAME[[VARNAME]] else ENVNAME[[VARNAME]][[ROWNAME]]) ## commas will get inserted after ROWNAME
    else {
        if(is.na(ROWNAME)) ENVNAME[[VARNAME]] <- value
        else ENVNAME[[VARNAME]][[ROWNAME]] <- value
        return(invisible(value))
    }
})


createDefault_isDataObj <- function(obj) {
    if(length(obj) == 0) return(FALSE)
    return(array(FALSE, dim = obj))
}

isValid <- function(value) {
    if(is(value, 'try-error')) return(FALSE)
    if(any(is.nan(value))) return (FALSE)
    if(any(is.na(value))) return(FALSE)
    if(any(abs(value)==Inf)) return(FALSE)
    return(TRUE)
}

whyInvalid <- function(value) {
    if(isValid(value)) { warning('Checking why a valid value is invalid.'); return(NULL) }
    if(any(is.nan(value))) return('nan')
    if(any(is.na(value))) return('na')
    if(any(abs(value)==Inf)) return('inf')
    stop('should never happen')
}

# The following roxygen is basically redundant with the method documentation for modelBaseClass::getConditionallyIndependentSets.  Not sure we need both.

#' Get a list of conditionally independent sets of nodes in a nimble model
#'
#' Conditionally independent sets of nodes are typically groups of latent states whose joint probability (density) will not change even if any other non-fixed node is changed.  Default fixed nodes are data nodes and parameter nodes (with no parent nodes), but this can be controlled.
#'
#' @param model A nimble model object (uncompiled or compiled).
#'
#' @param nodes A vector of node names or their graph IDs that are the starting nodes from which conditionally independent sets of nodes should be found.  If omitted, the default will be all latent nodes, defined as stochastic nodes that are not data and have at least one stochastic parent node (possible with determinstic nodes in between).  Note that this will omit latent states that have no hyperparameters.  An example is the first latent state in some state-space (time-series) models, which is sometimes declared with known prior.  See \code{type} because it relates to the interpretation of \code{nodes}.
#'
#' @param givenNodes A vector of node names or their graph IDs that should be considered as fixed and hence can be conditioned on.  If omitted, the default will be all data nodes and all parameter nodes, the latter defined as nodes with no stochastic parent nodes (skipping over deterministic parent nodes).
#'
#' @param omit A vector of node names or their graph IDs that should be omitted and should block further graph exploration. 
#'
#' @param inputType The method of graph exploration depends on what the nodes argument represents.  For \code{latent}, the input \code{nodes} are interpreted as latent states, from which both parent and descendent graph exploration should be done to find nodes in the same set (nodes that are NOT conditionally independent from each other).  For \code{param}, the input \code{nodes} are interpreted as parameters, so graph exploration begins from the  top (input) and explores descendents.  For \code{data}, the input \code{nodes} are interpreted as data nodes, so graph exploration begins from the bottom (input) explores parent nodes.
#' 
#' @param stochOnly Logical for whether only stochastic nodes should be returned (default = TRUE).  If FALSE, both deterministic and stochastic nodes are returned.
#' 
#' @param returnType Either \code{names} for returned nodes to be node names or \code{ids} for returned nodes to be graph IDs.
#' 
#' @param returnScalarComponents If FALSE (default), multivariate nodes are returned as full names (e.g. \code{x[1:3]}).  If TRUE, they are returned as scalar elements (e.g. \code{x[1]},  \code{x[2]},  \code{x[3]}).
#' 
#' @author Perry de Valpine
#'
#' @details This function returns sets of conditionally independent nodes.  Multiple input \code{nodes} might be in the same set or different sets, and other nodes (not in \code{nodes}) will be included.
#'
#' By default, deterministic dependencies of givenNodes are also
#' counted as given nodes.  This is relevant only for parent nodes.
#' This allows the givenNodes to include only stochastic nodes.  Say
#' we have A -> B -> C -> D.  A and D are givenNodes.  C is a latent
#' node.  B is a deterministic node.  By default, B is considered
#' given.  Otherwise, other dependent networks of nodes that that depend on B would be grouped
#' in the same output set as C, but this is usually what is wanted.
#' Any use of the resulting output must ensure that B is calculated when
#' necessary, as usual with nimble's model-generic programming.  To
#' turn off this feature, set
#' \code{nimbleOptions(groupDetermWithGivenInCondIndSets = FALSE)}.
#' 
#' @return List of nodes that are in conditionally independent sets.  With each set, nodes are returned in topologically sorted order.  The sets themselves are returned in topologically sorted order of their first nodes.
#'
#' @seealso There is a non-exported function \code{nimble:::testConditionallyIndependentSets(model, sets, initialize = TRUE)} that tests whether the conditional independence of sets is valid.  It should be the case that \code{nimble:::testConditionallyIndependentSets(model, getConditionallyIndependentSets(model), initialize = TRUE)} returns \code{TRUE}.
#'
getConditionallyIndependentSets <- function(model,
                                            nodes,
                                            givenNodes,
                                            omit = integer(),
                                            inputType = c("latent", "param", "data"),
                                            stochOnly = TRUE,
                                            returnType = 'names',
                                            returnScalarComponents = FALSE) {
  inputType <- match.arg(inputType)
  if(missing(nodes)) { # default to latent nodes
    nodeIDs <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE, returnType = 'ids')
  } else {
    if(is.character(nodes))
      nodeIDs <- model$expandNodeNames(nodes, returnType = 'ids')
    else
      nodeIDs <- nodes
  }
  if(missing(givenNodes)) { # default to top nodes and data nodes. need to be deliberate about end nodes
    givenNodeIDs <- c(model$getNodeNames(topOnly = TRUE, returnType = 'ids'),
                      model$getNodeNames(dataOnly = TRUE, returnType = 'ids'))
  } else {
    if(is.character(givenNodes))
        givenNodeIDs <- model$expandNodeNames(givenNodes, returnType = 'ids')
    else if(is.numeric(givenNodes))
        givenNodeIDs <- givenNodes
  }
  if(isTRUE(nimbleOptions("groupDetermWithGivenInCondIndSets"))) {
    givenNodeIDs <- unique(c(givenNodeIDs, model$getDependencies(givenNodeIDs, determOnly = TRUE, self = FALSE, returnType = 'ids')))
  }
  if(is.character(omit)) {
    # This mimcs getDependencies.  I think it allows omit to include split nodes, whereas getNodeNames would not.
    # It would not make sense for nodes or givenNode to include split nodes.
    elementIDs <- model$modelDef$nodeName2GraphIDs(omit, FALSE)
    omitIDs <- unique(model$modelDef$maps$elementID_2_vertexID[elementIDs],     ## turn into IDs in the graph
                      FALSE,
                      FALSE,
                      NA)
  }
  else if(is.numeric(omit))
    omitIDs <- omit
  
  startUp <- startDown <- TRUE
  if(inputType == "param") startUp <- FALSE
  if(inputType == "data") startDown <- FALSE
  result <- model$modelDef$maps$nimbleGraph$getConditionallyIndependentSets(
    nodeIDs = nodeIDs,
    givenNodeIDs = givenNodeIDs,
    omitIDs = omitIDs,
    startUp,
    startDown)
  if(returnType == 'ids' && returnScalarComponents)
    warning("NIMBLE development warning: calling getConditionallyIndependentSets with returnType = ids and returnScalarComponents may not be meaningful.")
  result <- lapply(result,
                   function(IDs) {
                     if(stochOnly) IDs <- IDs[model$modelDef$maps$types[IDs] == 'stoch']
                     if(returnType == 'ids') IDs
                     if(returnType == 'names') {
                       if(returnScalarComponents)
                         model$modelDef$maps$elementNames[IDs]
                       else
                         model$modelDef$maps$nodeNames[IDs]
                     }
                   })
  result
}

# testConditionallyIndependentSets checks whether a set of nodes are conditionally independent
# model: a nimble model
# sets: a list of node names or IDs
# intialize: should the model be forced into full initialization by full simulation (except data) and calculation?
#
# This works as follows:
#    For each focal set sets[[i]]:
#         Determine the logProb from calculating dependencies of sets[[i]] (which includes sets[[i]], data that depends on it, and deterministic nodes in between)
#         Simulate dependencies of all other sets to change their values.
#         Re-determine the logProb from calculating dependencies of sets[[i]]
#         If sets[[i]] is really conditionally independent of other sets, its logProb should be unchanged by having simulated with all other sets.
#
# Example: testConditionallyIndependentSets(model, getConditionallyIndependentSets(model), TRUE)
testConditionallyIndependentSets <- function(model, sets, initialize = TRUE) {
  if(initialize) { # would be better to use our initializeModel method, but I am doing a quick-and-dirty version here:
    model$simulate()
    model$calculate()
  }
  # sets is a list of stochastic (and optionally deterministic) nodes.
  # This function checks that the nodes in each element are conditionally independent of the others.
  # We check this by simulating all but one set and checking that the logProb of the one set hasn't changed.
  # We do that for each set.
  ok <- TRUE
  # Nodes for calculation/simulation for each set.
  calcNodeSets <- lapply(sets, function(x) model$getDependencies(x))
  for(i in seq_along(sets)) { # i is the set being currently checked
    prevLogProb <- model$calculate(calcNodeSets[[i]]) ## find the logProb for set i
    for(j in seq_along(sets)) {           # Simulate all other sets (with dependencies)
      if(i != j) {
        model$simulate(calcNodeSets[[j]]) # This assumes the bottom nodes of sets are data, which won't be simulated.
      }
    }
    newLogProb <- model$calculate(calcNodeSets[[i]]) # find the logProb for set i again
    if(prevLogProb != newLogProb) {                  # if it has changed, that set is not conditionally independent all the others
      message("Problem: Set ", i, " is not conditionally independent.")
      ok <- FALSE
    }
  }
  ok
}

#' Information on initial values in a NIMBLE model
#'
#'  Having uninitialized nodes in a NIMBLE model can potentially cause some algorithms to fail and can lead to poor performance in others.  Here are some
#'  general guidelines on how non-initialized variables can affect performance:
#'  \itemize{
#'    \item MCMC will auto-initialize but will do so from the prior distribution.  This can cause slow convergence, especially in the case of diffuse priors.
#'    \item Likewise, particle filtering methods will initialize top-level parameters from their prior distributions, which can lead to errors or poor performance in these methods.
#' }
#' Please see this Section (\url{https://r-nimble.org/html_manual/cha-mcmc.html#sec:initMCMC}) of the NIMBLE user manual for further suggestions.
#'
#' @name modelInitialization
#' @rdname modelInitialization
NULL

