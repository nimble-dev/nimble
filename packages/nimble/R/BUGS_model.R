#' Class \code{modelBaseClass}
#' @aliases modelBaseClass getVarNames getNodeNames topologicallySortNodes resetData setData isData getDependencies setInits checkConjugacy newModel [[,modelBaseClass-method [[<-,modelBaseClass-method
#' @export
#' @description
#' This class underlies all NIMBLE model objects: both R model objects created from the return value of nimbleModel(), and compiled model objects.
#' The model object contains a variety of member functions, for providing information about the model structure, setting or querying properties of the model, or accessing various internal components of the model.
#' These member functions of the modelBaseClass are commonly used in the body of the \code{setup} function argument to nimbleFunction(), to aid in preparation of node vectors, nimbleFunctionLists, and other runtime inputs.
#' See documentation for \code{nimbleModel} for details of creating an R model object.
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
modelBaseClass <- setRefClass('modelBaseClass',
                              fields = list(
                                  modelDef = 'ANY',     
                                  nodes = 'ANY',       #list
                                  vars = 'ANY',         
                                  graph = 'ANY',        
                                  defaultModelValues = 'ANY',
                                  name = 'ANY', 		#character  
                                  ##.ModelValuesLookUpName = 'character',
                                  isDataVars = 'ANY', #list           ## list with the dimensions of isData_vars
                                  isDataEnv = 'ANY',	#environment      ## environment holding 'logical' objects, with isData flags
                                  classEnvironment = 'ANY', # environment in which the reference classes will be defined
                                  origData = 'ANY',
                                  origInits = 'ANY',
                                  nimbleProject = 'ANY'
                                  ),
                              methods = list(
                                  calculate = function(nodes) nimble:::calculate(.self, nodes),
                                  calculateDiff = function(nodes) nimble:::calculateDiff(.self, nodes),
                                  getLogProb = function(nodes) nimble:::getLogProb(.self, nodes),
                                  simulate = function(nodes, includeData = FALSE) nimble:::simulate(.self, nodes, includeData),
                                  getParam = function(node, param) nimble:::getParam(.self, node, param),
                                  
                                  getGraph = function() graph,
                                  setGraph = function(value) graph <<- value,
                                  getModelDef = function() modelDef,
                                  setModelDef = function(value) modelDef <<- value,
                                  getMaps = function(mapName, all = FALSE){
                                  	if(all == TRUE)		return(modelDef$maps)
                                  	return(modelDef$maps[[mapName]])
                                   },
                                   isNodeEnd = function(nodeNames){  #Note: it says nodeNames, but graphIDs are fine too. Actually they are better
                                   	if(is.character(nodeNames))
                                   		nodeNames = expandNodeNames(nodeNames, returnType = 'ids')
                                   	return(modelDef$maps$isEndNode_byGID[nodeNames])
                                   },

                                  ## returns the type of one or more node names, e.g., 'stoch' or 'determ'
                                  getNodeType = function(nodes) {
                                      graphIDs <- modelDef$nodeName2GraphIDs(nodes)
                                      types <- getMaps('types')[graphIDs]
                                      return(types)
                                  },

                                  ## returns the declaration ID corresponding to nodes
                                  getDeclID = function(nodes) {
                                      graphIDs <- modelDef$nodeName2GraphIDs(nodes)
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
                                  getNodeDistribution = function(nodes) {
                                      if(length(nodes)==1) getDeclInfo(nodes)[[1]][['distributionName']]
                                      else unlist(lapply(getDeclInfo(nodes), `[[`, 'distributionName'))
                                      ##getDeclInfo(nodes)[[1]]$getDistributionName()
                                  },

                                  ## returns the expr corresponding to 'param' in the distribution of 'node'
                                  getNodeParamExpr = function(node, param) {
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
                                  getNodeValueExpr = function(node) {
                                      expr <- getDeclInfo(node)[[1]]$valueExprReplaced
                                      unrolledIndices <- getUnrolledIndicesList(node)
                                      subExpr <- eval(substitute(substitute(EXPR, unrolledIndices), list(EXPR = expr)))
                                      return(subExpr)
                                  },

                                  isDiscrete = function(node) {
                                      dist <- getNodeDistribution(node)
                                      # explicit reference to namespace needed as class definition objects inheriting from modelBaseClass not in namespace
                                      discrete <- getDistribution(dist)$discrete
                                      return(discrete)
                                  },

                                  isBinary = function(node) {
                                      dist <- getNodeDistribution(node)
                                      if(dist == 'dbern') return(TRUE)
                                      if(dist == 'dbin') {
                                          if(getNodeParamExpr(node, 'size') == 1)
                                              return(TRUE)
                                      }
                                      return(FALSE)
                                  },

                                  isTruncated = function(node) {
                                      di <- getDeclInfo(node)[[1]]
                                      if(di$truncated) return(TRUE) else return(FALSE)
                                  },


                                  getVarNames = function(includeLogProb = FALSE, nodes, includeData = TRUE) {                                  
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
                                          ans <- unique(removeIndexing(nodes))
                                          if(!all(ans %in% modelDef$varNames))
                                              stop(c('invalid node names provided to model$getVarNames') )
                                      }
                                      if(!includeData) {
                                          allData <- unlist(lapply(mget(ans, envir = tm$isDataEnv, inherits = FALSE, ifnotfound = TRUE), all))
                                          ans <- ans[!allData]
                                      }
    	                              return(ans)
                                  },

                                  getNodeFunction = function(nodes) {
                                      gids <- modelDef$nodeName2GraphIDs(nodes)
                                      dclids <- modelDef$graphIDs2indexedNodeInfo(gids)$declIDs
                                      if(length(dclids) == 1)
                                          return(nodeFunctions[[dclids]])
                                      else
                                          return(nodeFunctions[dclids])
                                  },
                                  
                                  getNodeNames = function(determOnly = FALSE, stochOnly = FALSE,
                                                          includeData = TRUE, dataOnly = FALSE, includeRHSonly = FALSE,
                                                          topOnly = FALSE, latentOnly = FALSE, endOnly = FALSE,
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

returnType: Character argument specific type object returned. Options are \'names\' (returns character vector) and \'ids\' (returns numeric graph IDs for model)

returnScalar Componenets: Logical argument specifying whether multivariate nodes should return full node name (i.e. \'x[1:2]\') or should break down into scalar componenets (i.e. \'x[1]\' and \'x[2]\')

Details: Multiple logical input arguments may be used simultaneously.  For example, model$getNodeNames(endOnly = TRUE, dataOnly = TRUE) will return all end-level nodes from the model which are designated as \'data\'.
'
                                      validValues = rep(TRUE, length(modelDef$maps$graphIDs) )
                                      if(!includeRHSonly)		validValues[modelDef$maps$types == 'RHSonly'] <- FALSE
                                      if(determOnly)			validValues[modelDef$maps$types != 'determ']	<- FALSE
                                      if(stochOnly)			validValues[modelDef$maps$types != 'stoch']	<- FALSE
                                      if(!includeData)		validValues[isDataFromGraphID(modelDef$maps$graphIDs)] <- FALSE
                                      if(dataOnly)			validValues[!isDataFromGraphID(modelDef$maps$graphIDs)] <- FALSE
                                      if(topOnly)				validValues[-modelDef$maps$top_IDs] <- FALSE
                                      if(latentOnly)			validValues[-modelDef$maps$latent_IDs] <- FALSE
                                      if(endOnly)				validValues[-modelDef$maps$end_IDs] <- FALSE
                                      
                                      ans <- expandNodeNames(modelDef$maps$graphID_2_nodeName[validValues], 
                                      						 returnScalarComponents = returnScalarComponents,
                                      						 returnType = returnType) 
                                      return(ans)                                      
                                  },
                                  
                                  expandNodeNames = function(nodeNames, env = parent.frame(), returnScalarComponents = FALSE, returnType = 'names', sort = FALSE){
                                      '
Takes a vector of nodeNames and returns the unique and expanded names in the model, i.e. \'x\' expands to \'x[1]\', \'x[2]\', ...

Arguments:

nodeNames: a vector of characters of nodes to be expanded. Alternatively, can be a vector of integer graph IDs, but this use is intended only for advanced users 

returnScalarComponents: should multivariate nodes (i.e. dmnorm or dmulti) be broken up into scalar components?

returnType: return type. Options are \'names\' (character vector) or \'ids\' (graph IDs)

sort: should names be topologically sorted before being returned?
'

                                      if(length(nodeNames) == 0) return(if(returnType=='names') character() else numeric())
                                      graphID <- modelDef$nodeName2GraphIDs(nodeNames, !returnScalarComponents)
                                      if(sort) 
                                          graphID <- sort(graphID)
                                      if(returnType == 'names'){
                                          if(returnScalarComponents) nodeNames <- modelDef$maps$elementNames[graphID] ## these are really elementIDs
                                          else nodeNames <- modelDef$maps$graphID_2_nodeName[graphID]
                                          return(nodeNames)
                                      }
                                      if(returnType == 'ids'){
                                          if(returnScalarComponents) print("NIMBLE development warning: returning IDs of scalar components may not be meaningful.  Checking to see if we ever see this message.") 
                                          return(graphID)
                                      }                  
                                      if(!(returnType %in% c('ids','names')))
                                      	stop('instead expandNodeNames, imporper returnType chosen')
                                  },
                                  
                                  topologicallySortNodes = function(nodeNames, returnType = 'names') {
                                      '
Sorts the input list of node names according to the topological dependence ordering of the model structure. 

Arguments:

nodeNames: A character vector of node names, which is to be topologically sorted. Alternatively can be a numeric vector of graphIDs

returnType: character vector indicating return type. Choices are "names" or "ids"

Details: This function merely reorders its input argument.  This may be inportany prior to calls such as simulate(model, nodes) or calculate(model, nodes), to enforce that the operation is performed in topological order.
'
                                      nodeIDs <- expandNodeNames(nodeNames, returnType = 'ids')			#modelDef$maps$nodeName_2_graphID[nodeNames]
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
                                      list2env(lapply(isDataVars, createDefault_isDataObj), isDataEnv)
                                  },
                                  
                                  resetData = function() {
'
Resets the \'data\' property of ALL model nodes to FALSE.  Subsequent to this call, the model will have no nodes flagged as \'data\'. 
'
                                      ## re-initializes the 'isDataEnv', setting everything to 'FALSE'
                                      init_isDataEnv()
                                      return(invisible(NULL))
                                  },
                                  
                                  setData = function(data, warnAboutMissingNames = TRUE) {
'
Sets the \'data\' flag for specified nodes to TRUE, and also sets the value of these nodes to the value provided.  This is the exclusive method for specifying \'data\' nodes in a model object.  When a \'data\' argument is provided to \'nimbleModel()\', it uses this method to set the data nodes.

Arguments:

data: A named list.  The names of list elements must correspond to model variable names.  The elements of the list must be of class numeric, with size and dimension each matching the corresponding model variable.  These numeric scalars, vectors, arrays, etc, may only contain numeric data, or NAs.

Details: If a list element contains some number of NA values, then the model nodes corresponding to these NAs will not have their value set, and will not be designated as \'data\'.  Only model nodes corresponding to numeric values in the argument list elements will be designated as data.  Designating a deterministic model node as \'data\' will result in an error.  Designating part of a multivariate node as \'data\' and part as non-data (NA) will resilt in an error; multivariate nodes must be entirely data, or entirely non-data.
'
                                      origData <<- data
                                      ## argument is a named list of data values.
                                      ## all nodes specified (except with NA) are set to that value, and have isDataEnv$VAR set to TRUE
                                      if(length(data) == 0) return()
                                      if(is.null(names(data)))   stop('\'data\' argument must by a named list')
                                      for(iData in seq_along(data)) {
                                          varName <- names(data)[iData]
                                          varValue <- data[[iData]]
                                          if(!(varName %in% names(isDataVars))) {
                                              ## when data is from modelDef$constantsList,
                                              ## it is possible that the constants don't exist on LHS of BUGS decls
                                              ## and so are not variables in the model.  In that case we don't want to issue the warning.
                                              if(warnAboutMissingNames) {
                                                      stop('variable name not suitable for setData(): ', varName)
                                                  } else next
                                              }
                                          if(length(dimOrLength(varValue, scalarize = TRUE)) != length(isDataVars[[varName]]))   stop(paste0('incorrect size or dim in data: ', varName))
                                          if(!(all(dimOrLength(varValue, scalarize = TRUE) == isDataVars[[varName]])))   stop(paste0('incorrect size or dim in data: ', varName))
                                          assign(varName, varValue, inherits = TRUE)
                                          isDataVarValue <- !is.na(varValue)
                                          assign(varName, isDataVarValue, envir = isDataEnv)
                                      }
                                   ##   testDataFlags()  ## this is slow for large models.  it could be re-written if we want to use it routinely
                                      return(invisible(NULL))
                                  },
                                  
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
                                  
                                  isData = function(nodeNames) {
'
Returns a vector of logical TRUE / FALSE values, corresponding to the \'data\' flags of the input node names. 

Arguments:

nodeNames: A character vector containing model variable or node names.

Details: The variable or node names specified is expanded into a vector of model node names.  A logical vector is returned, indicating whether each model node has been flagged as containing \'data\'.
'
                                  g_id = modelDef$nodeName2GraphIDs(nodeNames)
                                  		return(isDataFromGraphID(g_id))                                  
                                  },

                                  isDataFromGraphID = function(g_id){
                                   	nodeNames <- modelDef$maps$graphID_2_nodeName[g_id]	
                                  	ret <- unlist(lapply(as.list(nodeNames),
                                                  function(nn)     
                                                    return(as.vector(eval(parse(text=nn, keep.source = FALSE)[[1]],
                                                                                envir=isDataEnv))[[1]])))
                                    if(is.null(ret))   ret <- logical(0)
                                    return(ret)
                                  },

                                  getDependenciesList = function(returnNames = TRUE, sort = TRUE) {
'
Returns a list of all neighbor relationships.  Each list element gives the one-step dependencies of one vertex, and the element name is the vertex label (integer ID or character node name)

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
                                  
                                  getDependencies = function(nodes, omit = character(), self = TRUE,
                                      determOnly = FALSE, stochOnly = FALSE,
                                      includeData = TRUE, dataOnly = FALSE,
                                      includeRHSonly = FALSE, downstream = FALSE,
                                      returnType = 'names', returnScalarComponents = FALSE) {
'
Returns a character vector of the nodes dependent upon the input argument nodes, sorted topologically according to the model graph.  Aditional input arguments provide flexibility in the values returned.

Arguments:

nodes: Character vector of node names, with index blocks allowed, and/or variable names, the dependents of which will be returned.

omit: Character vector of node names, which will be omitted from the nodes returned.  In addition, dependent nodes subsequent to these omitted nodes will not be returned.  The omitted nodes argument serves to stop the downward search within the hierarchical model struture, and excludes the specified node.

self: Logical argument specifying whether to include the input argument nodes in the return vector of dependent nodes.  Default is TRUE.

determOnly: Logical argument specifying whether to return only deterministic nodes.  Default is FALSE.

stochOnly: Logical argument specifying whether to return only stochastic nodes.  Default is FALSE.

includeData: Logical argument specifying whether to include \'data\' nodes (set via the member method setData).  Default is TRUE.

dataOnly: Logical argument specifying whether to return only \'data\' nodes.  Default is FALSE.

includeRHSonly: Logical argument specifying whether to include right-hand-side-only nodes (model nodes which never appear on the left-hand-side of ~ or <- in the model code).  These nodes are neither stochastic nor deterministic, but instead function as variable inputs to the model.  Default is FALSE.

downstream: Logical argument specifying whether the downward search through the model hierarchical structure should continue beyond the first and subsequent stochastic nodes encountered, hence returning all nodes downstream of the input nodes.  Default is FALSE.

returnType: Character argument specific type object returned. Options are \'names\' (returns character vector) and \'ids\' (returns numeric graph IDs for model)

returnScalar Componenets: Logical argument specifying whether multivariate nodes should return full node name (i.e. \'x[1:2]\') or should break down into scalar componenets (i.e. \'x[1]\' and \'x[2]\')

Details: The downward search for dependent nodes propagates through deterministic nodes, but by default will halt at the first level of stochastic nodes encountered.
'

                                      if(inherits(nodes, 'character')) {
                                          elementIDs <- modelDef$nodeName2GraphIDs(nodes, !returnScalarComponents)
                                          if(returnScalarComponents)
                                             # nodeIDs <- .Internal(unique(modelDef$maps$elementID_2_vertexID[elementIDs],     ## turn into IDs in the graph
                                              nodeIDs <- unique(modelDef$maps$elementID_2_vertexID[elementIDs],     ## turn into IDs in the graph
                                                                   FALSE,
                                                                   FALSE,
                                                                   NA)
                                          else
                                              nodeIDs <- elementIDs
                                      }
                                      else if(inherits(nodes, 'numeric'))
                                          nodeIDs <- nodes
                                      ## else if(inherits(nodes, 'nodeVector')){ 
                                      ##     if(!returnScalarComponenets)
                                      ##         nodeIDs <- nodes$getOrigIDs_functions()
                                      ##     else
                                      ##         nodeIDs <- nodes$getOrigIDs_values()
                                      ## }
                                      
                                      if(inherits(omit, 'character')) {
                                          elementIDs <- modelDef$nodeName2GraphIDs(omit, !returnScalarComponents)
                                          if(returnScalarComponents)
#                                              omitIDs <- .Internal(unique(modelDef$maps$elementID_2_vertexID[elementIDs],
                                              omitIDs <- unique(modelDef$maps$elementID_2_vertexID[elementIDs],
                                                                   FALSE,
                                                                   FALSE,
                                                                   NA)
                                          else
                                              omitIDs <- elementIDs
                                      }
                                      else if(inherits(omit, 'numeric'))
                                          omitIDs <- omit
                                      ## else if(inherits(omit, 'nodeVector')){ 
                                      ##     if(!returnScalarComponenets)
                                      ##         omitIDs <- omit$getOrigIDs_functions()
                                      ##     else
                                      ##         omitIDs <- omit$getOrigIDs_values()
                                      ## }
                                      
                                      depIDs <- gd_getDependencies_IDs(graph = getGraph(), maps = getMaps(all = TRUE), nodes = nodeIDs, omit = omitIDs, downstream = downstream)
                                      if(!includeRHSonly) depIDs <- depIDs[modelDef$maps$types[depIDs] != 'RHSonly']
                                      if(determOnly)	depIDs <- depIDs[modelDef$maps$types[depIDs] == 'determ']
                                      if(stochOnly)	depIDs <- depIDs[modelDef$maps$types[depIDs] == 'stoch']
                                      if(!self)	depIDs <- setdiff(depIDs, nodeIDs)
                                      if(!includeData)	depIDs <- depIDs[!isDataFromGraphID(depIDs)]
                                      if(dataOnly)		depIDs <- depIDs[isDataFromGraphID(depIDs)]
                                      
                                      depIDs <- modelDef$nodeName2GraphIDs(modelDef$maps$graphID_2_nodeName[depIDs], !returnScalarComponents)
                                      if(returnScalarComponents)
                                          depIDs = unique(depIDs, FALSE, FALSE, NA)
                                      if(returnType == 'ids'){
                                          if(returnScalarComponents) print("nimble development warning: calling getDependencies with returnType = ids and returnScalarComponents may not be meaningful.")
                                          return(depIDs)
                                      }		                       			 
                                      if(returnType == 'names') {
                                          if(returnScalarComponents)
                                              return(modelDef$maps$elementNames[depIDs])
                                          retVal <- modelDef$maps$nodeNames[depIDs]
                                          ##attr(retVal, 'nodeName') <- TRUE
                                          return(retVal)
                                      }
                                      if(!(returnType %in% c('ids', 'names')))
                                          stop('instead getDependencies, imporper returnType chosen')
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
                                      origInits <<- inits
                                      
                                      for(i in seq_along(inits)) {
                                          dataVals <- .self$isData(names(inits)[[i]])
                                          if(any(dataVals)) {
                                              .self[[names(inits)[i]]][!dataVals] <- inits[[i]][!dataVals]
                                              if(any(!is.na(inits[[i]][dataVals])))
                                                  warning("Ignoring values in inits for data nodes: ", names(inits)[[i]], ".")
                                          } else  .self[[names(inits)[i]]] <- inits[[i]]
                                      }
                                  },
                                  
                                  checkConjugacy = function(nodeVector) {
                                      '
Determines whether or not the input nodes appear in conjugate relationships

Arguments:

nodeVector: A character vector specifying one or more node or variable names.  If omitted, all stochastic non-data nodes are checked for conjugacy.

Details: The return value is a named list, with an element corresponding to each conjugate node.  The list names are the conjugate node names, and list elements are the control list arguments required by the corresponding MCMC conjugate sampler functions.  If no model nodes are conjugate, an empty list is returned.
'
                                      if(missing(nodeVector))
                                          nodeVector <- getNodeNames(stochOnly=TRUE, includeData=FALSE)
                                      nodeVector <- expandNodeNames(nodeVector)
                                      conjugacyRelationshipsObject$checkConjugacy(.self, nodeVector)
                                  },
                                  checkConjugacy2 = function(nodeVector) {
                                      if(missing(nodeVector)) nodeVector <- getNodeNames(stochOnly=TRUE, includeData=FALSE)
                                      nodeIDs <- expandNodeNames(nodeVector, returnType = 'ids')
                                      conjugacyRelationshipsObject$checkConjugacy2(.self, nodeIDs)
                                  },
                                  check = function() {
                                      '
Checks for common errors in model specification, including missing values, inability to calculate/simulate on a node, and dimension/size mismatches
'
                                      # first do size checking; do before check of calculate/simulate so LHS of deterministic expressions are not filled in

                                      for(j in seq_along(.self$modelDef$declInfo)) {
                                              declInfo <- .self$modelDef$declInfo[[j]]
                                              nn <- length(declInfo$nodeFunctionNames)
                                              nfn <- declInfo$nodeFunctionNames[nn]
                                              ## NEWNODEFXNS
                                              nf <- .self$nodeFunctions[[j]]
                                              #context <- as.list(declInfo$unrolledIndicesMatrix[nrow(declInfo$unrolledIndicesMatrix), ])

                                              if(declInfo$type == 'determ') {
                                                  # check LHS and RHS are same size/dim
                                                  # need to eval within nf; constants not present otherwise
                                                  RHSsize <- try(dimOrLength(eval(codeSubstitute(declInfo$valueExprReplaced, as.list(nf)))), silent = TRUE)

                                                  LHSsize <- try(dimOrLength(eval(codeSubstitute(declInfo$targetExprReplaced, as.list(nf)))), silent = TRUE)
                                                  # apparently implicit dropping of size 1 dimensions is ok in determ node calcs
                                                  if(!is(RHSsize, 'try-error') && !is(LHSsize, 'try-error')) {
                                                      if(length(RHSsize) > 1 && any(RHSsize == 1))
                                                          RHSsize <- RHSsize[RHSsize != 1]
                                                      if(length(LHSsize) > 1 && any(LHSsize == 1))
                                                          LHSsize <- LHSsize[LHSsize != 1]
                                                      
                                                      if(!identical(LHSsize, RHSsize))
                                                          stop("Size/dimension mismatch between left-hand side and right-hand size of BUGS expression: ", deparse(declInfo$code))
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
                                                  dist <- deparse(declInfo$valueExprReplaced[[1]])

                                                  distDims <- as.integer(sapply(getDistribution(dist)$types, function(x) x$nDim))
                                                  nms <- names(getDistribution(dist)$types)
                                                  names(distDims) <- nms
                                                  
                                                  sizes <- list(); length(sizes) <- length(nms); names(sizes) <- nms

                                                  for(k in seq_along(nms)) {
                                        # sometimes get_foo not found in env of nf (and doesn't appear in ls(nf) )
                                                      ##fun <- as.call(parse(text = paste0("nf$get_", nms[k])))
                                                      ##e = try(eval(fun))
                                                      ## NEWNODEFXN
                                                      e <- try(.self$getParam(nfn, nms[k]))
                                                      
                                                      if(!is(e, "try-error")) {
                                                          sizes[[nms[k]]] <- dimOrLength(e)
                                                          if(prod(sizes[[nms[[k]]]]) == 1) sizes[[nms[[k]]]] <- numeric()
                                                      } else warning(paste0("Unable to calculate parameter '", nms[k], "'; this may simply reflect that there are missing values in model variables."))
                                                  }
                                        # check dimensions based on varInfo
                                                  if(length(declInfo$targetExprReplaced) > 1) {
                                                      LHSvar <- deparse(declInfo$targetExprReplaced[[2]])
                                                  } else LHSvar <- deparse(declInfo$targetExprReplaced)
                                                  if(.self$modelDef$varInfo[[LHSvar]]$nDim < distDims['value'])
                                                      stop("Dimension of '", LHSvar, "' does not match required dimension for the distribution '", dist, "'. Necessary dimension is ", distDims['value'], ".", ifelse(distDims['value'] > 0, paste0(" You may need to include explicit indexing information, e.g., variable_name", ifelse(distDims['value'] < 2, "[1:2].", "[1:2,1:2].")), ''))
                                                  nms2 <- nms[nms%in%names(declInfo$valueExprReplaced)]
                                                  for(k in seq_along(nms2)) {
                                                      if(!is.numeric(declInfo$valueExprReplaced[[nms2[k]]]) && !(dist == 'dinterval' && nms2[k] == 'c') && ( length(declInfo$valueExprReplaced[[nms2[k]]]) ==1 || deparse(declInfo$valueExprReplaced[[nms2[k]]][[1]]) == '[' )) {  # can only check variables not expressions or constants
                                                          # also dinterval can take 'c' param as scalar or vector, so don't check
                                                          if(length(declInfo$valueExprReplaced[[nms2[k]]]) > 1) {
                                                              var <- deparse(declInfo$valueExprReplaced[[nms2[k]]][[2]])
                                                          } else var <- deparse(declInfo$valueExprReplaced[[nms2[k]]])
                                                          if(var %in% names(.self$modelDef$varInfo) && .self$modelDef$varInfo[[var]]$nDim < distDims[nms2[k]])
                                                              # check less than because variable dim can be bigger than node dim
                                                              stop("Dimension of '", nms2[k], "' does not match required dimension for the distribution '", dist, "'. Necessary dimension is ", distDims[nms2[k]], ".", ifelse(distDims[nms2[k]] > 0, paste0(" You may need to include explicit indexing information, e.g., variable_name", ifelse(distDims[nms2[k]] < 2, "[1:2].", "[1:2,1:2].")), ''))
                                                      }
                                                  }
                                                  
                                        # actual dimensions
                                                  dims <- sapply(sizes, length)
                                                  toCheck <- names(dims[!is.na(sizes) & sapply(sizes, function(x) !is.null(x))])
                                                  if(dist == 'dinterval') toCheck <- toCheck[toCheck != 'c']
                                        # check dimensions based on empirical size of variables
                                                  if(!identical(dims[toCheck], distDims[toCheck])) {
                                                      mismatches <- which(dims[toCheck] != distDims[toCheck])
                                                      stop("Dimension of distribution argument(s) '", paste(names(mismatches), collapse = ","), "' does not match required dimension(s) for the distribution '", dist, "'. Necessary dimension(s) are ", paste(distDims[toCheck][mismatches], collapse = ","), ".", ifelse(any(distDims[toCheck][mismatches] == 1), " You may need to ensure that you have explicit vectors and not one-row or one-column matrices.", ""))                                                          
                                                  }
                                                  
                                        # check sizes
                                                  mats <- dims == 2
                                                  vecs <- dims == 1
                                                  matRows <- unlist(sapply(sizes[mats], `[`, 1))
                                                  matCols <- unlist(sapply(sizes[mats], `[`, 2))
                                                  if(!length(unique(c(matRows, matCols, unlist(sizes[vecs])))) <= 1)
                                                      if(dist %in% names(distributionsInputList)) {
                                                          stop("Size/dimension mismatch amongst vectors and matrices in BUGS expression: ", deparse(declInfo$code))
                                                      } else {
                                                          warning("Possible size/dimension mismatch amongst vectors and matrices in BUGS expression: ", deparse(declInfo$code), ". Ignore this warning if the user-provided distribution has multivariate parameters with distinct sizes.")                                                                                                                                   }
                                                  
                                              }
                                      }

                                      # check for missing values and inability to calculate/simulate
                                      lp <- try(nimble:::calculate(.self))
                                      if(!isValid(lp)) {
                                          varsToCheck <- character()
                                          for(v in .self$getVarNames())
                                              if(!isValid(.self[[v]]) || !isValid(nimble:::getLogProb(.self, setdiff(expandNodeNames(v), modelDef$maps$nodeNamesRHSonly))))
                                                  varsToCheck <- c(varsToCheck, v)
                                          badVars <- list(na=character(), nan=character(), inf=character())
                                      ##nns <- getNodeNames(includeRHSonly = TRUE)
                                          nns <- expandNodeNames(varsToCheck)
                                          nns <- topologicallySortNodes(nns)   ## should be unnecessary; just in case
                                          for(nn in nns) {
                                              val <- .self[[nn]]
                                              type <- getNodeType(nn)
                                              if(length(type) > 1) stop('something wrong with Daniel\'s understanding of nimbleModel')
                                              if(type == 'RHSonly') {
                                                  if(!isValid(val)) badVars[[whyInvalid(val)]] <- c(badVars[[whyInvalid(val)]], nn)
                                              } else if(type == 'determ') {
                                                  test <- try(nimble:::calculate(.self, nn))
                                                  if(class(test) == 'try-error')
                                                      cat("Note: cannot calculate logProb for node ", nn, ".\n")
                                                  val <- .self[[nn]]
                                                  if(!isValid(val)) badVars[[whyInvalid(val)]] <- c(badVars[[whyInvalid(val)]], nn)
                                              } else if(type == 'stoch') {
                                                  if(!isValid(val)) badVars[[whyInvalid(val)]] <- c(badVars[[whyInvalid(val)]], nn)
                                                  test <- try(val <- nimble:::calculate(.self, nn))
                                                  if(class(test) == 'try-error')
                                                      cat("Note: cannot calculate logProb for node ", nn, ".\n")
                                                  
                                                  if(!isValid(val)) badVars[[whyInvalid(val)]] <- c(badVars[[whyInvalid(val)]], paste0('logProb_', nn))
                                              } else stop('unknown node type: ', type)
                                          }
                                          badVars <- lapply(badVars, removeIndexing)
                                          badVars <- lapply(badVars, unique)
                                          badVars <- lapply(badVars, function(nns) if(length(nns>0)) paste0(nns, collapse=', '))
                                          conds <- list(c('na','NAs'), c('nan','NaNs'), c('inf','Infinite values'))
                                          for(i in seq_along(conds)) {
                                              v <- badVars[[conds[[i]][1]]]
                                              m <- conds[[i]][2]
                                              if(!is.null(v)) cat(m, ' were detected in model variable', if(grepl(',',v)) 's' else '', ': ', v, ".\n", sep = '')
                                          }
                                      }

                                  },

                                  newModel = function(data = NULL, inits = NULL, modelName = character(), replicate = FALSE, check = getNimbleOption('checkModel')) {
                                      '
Returns a new R model object, with the same model definiton (as defined from the original model code) as the existing model object.

Arguments:

data: A named list specifying data nodes and values, for use in the newly returned model.  If not provided, the data argument from the creation of the original R model object will be used.

inits: A named list specifying initial values, for use in the newly returned model.  If not provided, the inits argument from the creation of the original R model object will be used.

replicate: Logical specifying whether to repliate all current values and data flags from the current model in the new model.  If TRUE, then the data and inits arguments are not used.  Default value is FALSE.

check: A logical indicating whether to check the model object for missing or invalid values.  Default is given by the NIMBLE option \'checkModel\', see help on \'nimbleOptions\' for details.

modelName: An optional character string, used to set the internal name of the model object.  If provided, this name will propagate throughout the generated C++ code, serving to improve readability.

Details: The newly created model object will be identical to the original model in terms of structure and functionality, but entirely distinct in terms of the internal values.
'
                                      if(replicate) {
                                          newlyCreatedModel <- modelDef$newModel(check = FALSE)
                                          nimCopy(from = .self, to = newlyCreatedModel, logProb = TRUE)
                                          for(var in ls(isDataEnv)) newlyCreatedModel$isDataEnv[[var]] <- isDataEnv[[var]]
                                          if(check) newlyCreatedModel$check()
                                          return(newlyCreatedModel)
                                      }
                                      if(is.null(data)) data <- origData
                                      if(is.null(inits)) inits <- origInits
                                      modelDef$newModel(data = data, inits = inits, modelName = modelName, check = check)
                                  }
                              )
                              )



setMethod('[[', 'modelBaseClass',
          function(x, i) {
              if(!is.indexed(i)) {
                  eval(substitute(x$VAR, list(VAR=i)))
              } else {
                  parsedNode <- parse(text=i)[[1]]
                  parsedNode[[2]] <- substitute(x$VAR, list(VAR=parsedNode[[2]]))
                  eval(parsedNode)
              }
          }
)

setMethod('[[<-', 'modelBaseClass',
          function(x, i, value) {
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
    if(!is.null(code)) message(paste('confused about reaching end of insertSingleBrackets with ', deparse(code)))
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
                                       pointAt(.self, defaultModelValues, index = 1)
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
                                           code <- insertSingleIndexBrackets(code, modelDef$varInfo)
                                           LHS <- code[[2]]
                                           RHS <- code[[3]]
                                           altParams <- BUGSdecl$altParamExprs
                                           altParams <- lapply(altParams, insertSingleIndexBrackets, modelDef$varInfo)
                                           logProbNodeExpr <- BUGSdecl$logProbNodeExpr
                                           logProbNodeExpr <- insertSingleIndexBrackets(logProbNodeExpr, modelDef$logProbVarInfo)
                                           setupOutputExprs <- BUGSdecl$replacementNameExprs
                                           ## ensure they are in the same order as the columns of the unrolledIndicesMatrix, because that is assumed in nodeFunctionNew
                                           ## This can be necessary in a case like for(j in ...) for(i in ...) x[i,j] ~ ...; because x uses inner index first
                                           if(nrow(BUGSdecl$unrolledIndicesMatrix) > 0)
                                               setupOutputExprs <- setupOutputExprs[ colnames(BUGSdecl$unrolledIndicesMatrix) ]
                                           ## make a unique name
                                           thisNodeGeneratorName <- paste0(Rname2CppName(BUGSdecl$targetVarName), '_L', BUGSdecl$sourceLineNumber, '_', nimbleUniqueID())
                                           ## create the nimbleFunction generator (i.e. unspecialized nimbleFunction)
                                           nfGenerator <- nodeFunctionNew(LHS=LHS, RHS=RHS, name = thisNodeGeneratorName, altParams=altParams, logProbNodeExpr=logProbNodeExpr, type=type, setupOutputExprs=setupOutputExprs, evaluate=TRUE, where = where)
                                           nodeGenerators[[i]] <<- nfGenerator
                                           names(nodeGenerators)[i] <<- thisNodeGeneratorName
                                           nodeFunctionGeneratorNames[i] <<- thisNodeGeneratorName
                                           nodeFunctions[[i]] <<- nfGenerator(.self, BUGSdecl)
                                           names(nodeFunctions)[i] <<- thisNodeGeneratorName ## not sure what we need here
                                       }
                                   },
                                   
                                   buildNodeFunctions_old = function(where = globalenv(), debug = FALSE) {
                                       ## This creates the nodeFunctions, which are basically nimbleFunctions, for the model
                                       if(debug) browser()
                                       iNextNodeFunction <- 1
                                       nodeFunctions <<- vector('list', length = modelDef$numNodeFunctions)  ## for the specialized instances
                                       nodeFunctionGeneratorNames <<- character(modelDef$numNodeFunctions) ## possible error: should this be length(modelDef$declInfo)?
                                       nodeGenerators <<- vector('list', length = length(modelDef$declInfo)) ## for the nimbleFunctions
                                       for(i in seq_along(modelDef$declInfo)) {
                                           BUGSdecl <- modelDef$declInfo[[i]]
                                           if(BUGSdecl$numUnrolledNodes == 0) next
                                           ## extract needed pieces
                                           type <- BUGSdecl$type
                                           code <- BUGSdecl$codeReplaced
                                           code <- insertSingleIndexBrackets(code, modelDef$varInfo)
                                           LHS <- code[[2]]
                                           RHS <- code[[3]]
                                           altParams <- BUGSdecl$altParamExprs
                                           altParams <- lapply(altParams, insertSingleIndexBrackets, modelDef$varInfo)
                                           logProbNodeExpr <- BUGSdecl$logProbNodeExpr
                                           logProbNodeExpr <- insertSingleIndexBrackets(logProbNodeExpr, modelDef$logProbVarInfo)
                                           setupOutputExprs <- BUGSdecl$replacementNameExprs

                                           ## make a unique name
                                           thisNodeGeneratorName <- paste0(Rname2CppName(BUGSdecl$targetVarName), '_L', BUGSdecl$sourceLineNumber, '_', nimbleUniqueID())
                                           ## create the nimbleFunction generator (i.e. unspecialized nimbleFunction)
                                           nfGenerator <- nodeFunction(LHS=LHS, RHS=RHS, name = thisNodeGeneratorName, altParams=altParams, logProbNodeExpr=logProbNodeExpr, type=type, setupOutputExprs=setupOutputExprs, evaluate=TRUE, where = where)
                                           nodeGenerators[[i]] <<- nfGenerator

                                           newNodeFunctionNames <- BUGSdecl$nodeFunctionNames
                                           ## We include "_L[source line number]" in the names for the nodeGenerators so we can trace what line of BUGS code they came from
                                           ## This propagates to the C++ class names
                                           names(nodeGenerators)[i] <<- thisNodeGeneratorName

                                           ## If it is a singleton with no replacements, we can build the node simply:
                                           if(length(setupOutputExprs)==0) { 
                                               nodeFunctions[[iNextNodeFunction]] <<- nfGenerator(.self)
                                               nodeFunctionGeneratorNames[iNextNodeFunction] <<- thisNodeGeneratorName
                                               names(nodeFunctions)[iNextNodeFunction] <<- newNodeFunctionNames
                                               iNextNodeFunction <- iNextNodeFunction + 1
                                               next
                                           }

                                           ## It is either within for loops (contexts) and/or has a replacement
                                           ## so we construct code to call the nfGenerator with the needed arguments
                                           

                                           assign('MODEL_UNIQUE_NAME_', .self, envir = BUGSdecl$replacementsEnv)
                                           BUGSdecl$replacementsEnv[['nfGenCall_UNIQUE_NAME_']] <- as.call(c(list(quote(nfGenerator_UNIQUE_NAME_)), list(model = quote(MODEL_UNIQUE_NAME_)), lapply(setupOutputExprs, function(z) substitute(x[[index_UNIQUE_NAME_]], list(x = z)))))
                                           BUGSdecl$replacementsEnv[['nfGenerator_UNIQUE_NAME_']] <- nfGenerator
                                           ## nfGenCall is code for a call to nfGenerator, like nfGenerator(MODEL_UNIQUE_NAME_, j = j[[index_UNIQUE_NAME_]]) 
                                           evalq({
                                               nfGenWrap_UNIQUE_NAME_ <- function(index_UNIQUE_NAME_) x
                                               body(nfGenWrap_UNIQUE_NAME_) <- nfGenCall_UNIQUE_NAME_
                                           }, envir = BUGSdecl$replacementsEnv)
                                           
                                           numNewFunctions <- BUGSdecl$outputSize
                                           nodeFunctions[iNextNodeFunction-1+(1:numNewFunctions)] <<- evalq(lapply(1:outputSize, nfGenWrap_UNIQUE_NAME_), envir = BUGSdecl$replacementsEnv)
                                           nodeFunctionGeneratorNames[iNextNodeFunction-1+(1:numNewFunctions)] <<- thisNodeGeneratorName
                                           rm(list = c('MODEL_UNIQUE_NAME_', 'nfGenCall_UNIQUE_NAME_', 'nfGenerator_UNIQUE_NAME_', 'nfGenWrap_UNIQUE_NAME_'), envir = BUGSdecl$replacementsEnv)
                                           
                                           names(nodeFunctions)[iNextNodeFunction-1+(1:numNewFunctions)] <<- BUGSdecl$nodeFunctionNames
                                           iNextNodeFunction <- iNextNodeFunction + numNewFunctions
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
                # setData(modelDef$constantsList, warnAboutMissingNames = FALSE)
                # removed given new handling of lumped data and constants
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
    if(prod(dims) == 1)
        template <- activeBindingTemplateLength1NonScalar
    else
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
    if(isValid(value)) { warning('checking why a valid value is invalid'); return(NULL) }
    if(any(is.nan(value))) return('nan')
    if(any(is.na(value))) return('na')
    if(any(abs(value)==Inf)) return('inf')
    stop('should never happen')
}





