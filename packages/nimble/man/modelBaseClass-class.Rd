% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BUGS_model.R
\docType{class}
\name{modelBaseClass-class}
\alias{[[,modelBaseClass-method}
\alias{[[<-,modelBaseClass-method}
\alias{checkConjugacy}
\alias{expandNodeNames}
\alias{getCode}
\alias{getDependencies}
\alias{getDependenciesList}
\alias{getDimension}
\alias{getDistribution}
\alias{getDownstream}
\alias{getNodeNames}
\alias{getVarNames}
\alias{initializeInfo}
\alias{isBinary}
\alias{isData}
\alias{isDeterm}
\alias{isDiscrete}
\alias{isEndNode}
\alias{isMultivariate}
\alias{isStoch}
\alias{isTruncated}
\alias{isUnivariate}
\alias{modelBaseClass}
\alias{modelBaseClass-class}
\alias{newModel}
\alias{resetData}
\alias{setData}
\alias{setInits}
\alias{topologicallySortNodes}
\title{Class \code{modelBaseClass}}
\description{
This class underlies all NIMBLE model objects: both R model objects created from the return value of nimbleModel(), and compiled model objects.
The model object contains a variety of member functions, for providing information about the model structure, setting or querying properties of the model, or accessing various internal components of the model.
These member functions of the modelBaseClass are commonly used in the body of the \code{setup} function argument to nimbleFunction(), to aid in preparation of node vectors, nimbleFunctionLists, and other runtime inputs.
See documentation for \code{\link{nimbleModel}} for details of creating an R model object.
}
\section{Methods}{

\describe{
\item{\code{calculate(nodes)}}{See `help(calculate)`}

\item{\code{calculateDiff(nodes)}}{See `help(calculateDiff)`}

\item{\code{check()}}{Checks for errors in model specification and for missing values that prevent use of calculate/simulate on any nodes}

\item{\code{checkBasics()}}{Checks for size/dimension mismatches and for presence of NAs in model variables (the latter is not an error but a note of this is given to the user)}

\item{\code{checkConjugacy(nodeVector, restrictLink = NULL)}}{Determines whether or not the input nodes appear in conjugate relationships

Arguments:

nodeVector: A character vector specifying one or more node or variable names.  If omitted, all stochastic non-data nodes are checked for conjugacy.

Details: The return value is a named list, with an element corresponding to each conjugate node.  The list names are the conjugate node names, and list elements are the control list arguments required by the corresponding MCMC conjugate sampler functions.  If no model nodes are conjugate, an empty list is returned.}

\item{\code{expandNodeNames(nodes, env = parent.frame(), returnScalarComponents = FALSE,
  returnType = "names", sort = FALSE, unique = TRUE)}}{Takes a vector of names of nodes or variables and returns the unique and expanded names in the model, i.e. 'x' expands to 'x[1]', 'x[2]', ...

Arguments:

nodes: a vector of names of nodes (or variables) to be expanded. Alternatively, can be a vector of integer graph IDs, but this use is intended only for advanced users

returnScalarComponents: should multivariate nodes (i.e. dmnorm or dmulti) be broken up into scalar components?

returnType: return type. Options are 'names' (character vector) or 'ids' (graph IDs)

sort: should names be topologically sorted before being returned?

unique: should names be the unique names or should original ordering of nodes (after expansion of any variable names into node names) be preserved}

\item{\code{getBound(node, bound)}}{See `help(getBound)`                                                                                                                                                                                                                               }

\item{\code{getCode()}}{Return the code for a model after processing if-then-else statements, expanding macros, and replacing some keywords (e.g. nimStep for step) to avoid R ambiguity.}

\item{\code{getConditionallyIndependentSets(nodes, givenNodes, omit = integer(),
  inputType = c("latent", "param", "data"), stochOnly = TRUE,
  returnType = "names", returnScalarComponents = FALSE)}}{Get a list of conditionally independent sets of nodes in a nimble model.

Conditionally independent sets of nodes are typically groups of latent states whose joint conditional probability (density) will not change even if any other non-fixed node is changed.  Default fixed nodes are data nodes and parameter nodes (with no parent nodes), but this can be controlled.

model: A nimble model object (uncompiled or compiled).

nodes: A vector of node names or their graph IDs that are the starting nodes from which conditionally independent sets of nodes should be found.  If omitted, the default will be all latent nodes, defined as stochastic nodes that are not data and have at least one stochastic parent node (possible with deterministic nodes in between).  Note that this will omit latent states that have no hyperparameters.  An example is the first latent state in some state-space (time-series) models, which is sometimes declared with known prior.  See type because it relates to the interpretation of nodes.

givenNodes: A vector of node names or their graph IDs that should be considered as fixed and hence can be conditioned on.  If omitted, the default will be all data nodes and all parameter nodes, the latter defined as nodes with no stochastic parent nodes (skipping over deterministic parent nodes).

omit: A vector of node names or their graph IDs that should be omitted and should block further graph exploration. 

intputType: Type of input nodes provided in nodes argument.  For 'latent', the input nodes are interpreted as latent states, from which both downstream and upstream exploration should be done to find nodes in the same set (nodes that are not conditionally independent from each other).  For 'param', the input nodes are interpreted as parameters, so graph exploration begins from the  top (input) and proceeds downstream.  For 'data', the input nodes are interpreted and data nodes, so graph exploration begins from the bottom (input) and proceeds upstream.

stochOnly: Logical for whether only stochastic nodes should be returned (default = TRUE).  If FALSE, both deterministic and stochastic nodes are returned.

returnType: Either 'names' for returned nodes to be node names or 'ids' for returned nodes to be graph IDs.

returnScalarComponents: If FALSE (default), multivariate nodes are returned as full names (e.g. 'x[1:3]').  If TRUE, they are returned as scalar elements (e.g. 'x[1]',  'x[2]',  'x[3]').

Details: This function returns sets of conditionally independent nodes.  Multiple input nodes might be in the same set or different sets, and other nodes (not in codes) will be included.

By default, deterministic dependencies of givenNodes are also counted as given nodes.  This is relevant only for parent nodes. This allows the givenNodes to include only stochastic nodes.  Say we have A -> B -> C -> D.  A and D are givenNodes.  C is a latent node.  B is a deterministic node.  By default, B is considered given.  Otherwise, other branches that depend on B would be grouped in the same output set as C, but this is usually not what is wanted.  Any use of the resulting output must ensure that B is calculated when necessary, as usual with nimble's model-generic programming.  To turn off this feature, set nimbleOptions(groupDetermWithGivenInCondIndSets = FALSE).

There is a non-exported function `nimble:::testConditionallyIndependentSets(model, sets, initialize = TRUE)` that tests whether the conditional independence of sets is valid.  It should be the case that `nimble:::testConditionallyIndependentSets(model, getConditionallyIndependentSets(model), initialize = TRUE)` returns `TRUE`.

Return value: List of nodes that are in conditionally independent sets.  Within each set, nodes are returned in topologically sorted order.  The sets themselves are returned in topologically sorted order of their first nodes.}

\item{\code{getDependencies(nodes, omit = character(), self = TRUE,
  determOnly = FALSE, stochOnly = FALSE, includeData = TRUE,
  dataOnly = FALSE, includeRHSonly = FALSE, downstream = FALSE,
  returnType = "names", returnScalarComponents = FALSE)}}{Returns a character vector of the nodes dependent upon the input argument nodes, sorted topologically according to the model graph. In the genealogical metaphor for a graphical model, this function returns the "children" of the input nodes.  In the river network metaphor, it returns downstream nodes. By default, the returned nodes include the input nodes, include both deterministic and stochastic nodes, and stop at stochastic nodes. Aditional input arguments provide flexibility in the values returned.

Arguments:

nodes: Character vector of node names, with index blocks allowed, and/or variable names, the dependents of which will be returned.

omit: Character vector of node names, which will be omitted from the nodes returned.  In addition, dependent nodes subsequent to these omitted nodes will not be returned.  The omitted nodes argument serves to stop the downward search within the hierarchical model structure, and excludes the specified node.

self: Logical argument specifying whether to include the input argument nodes in the return vector of dependent nodes.  Default is TRUE.

determOnly: Logical argument specifying whether to return only deterministic nodes.  Default is FALSE.

stochOnly: Logical argument specifying whether to return only stochastic nodes.  Default is FALSE.  If both determOnly and stochOnly are TRUE, no nodes will be returned.

includeData: Logical argument specifying whether to include 'data' nodes (set via nimbleModel or the setData method).  Default is TRUE.

dataOnly: Logical argument specifying whether to return only 'data' nodes.  Default is FALSE.

includeRHSonly: Logical argument specifying whether to include right-hand-side-only nodes (model nodes which never appear on the left-hand-side of ~ or <- in the model code).  These nodes are neither stochastic nor deterministic, but instead function as variable inputs to the model.  Default is FALSE.

downstream: Logical argument specifying whether the downward search through the hierarchical model structure should continue beyond the first and subsequent stochastic nodes encountered, hence returning all nodes downstream of the input nodes.  Default is FALSE.

returnType: Character argument specifying type of object returned. Options are 'names' (returns character vector) and 'ids' (returns numeric graph IDs for model).

returnScalarComponenets: Logical argument specifying whether multivariate nodes should be returned as full node names (i.e. 'x[1:2]') or as scalar componenets (i.e. 'x[1]' and 'x[2]').

Details: The downward search for dependent nodes propagates through deterministic nodes, but by default will halt at the first level of stochastic nodes encountered.  Use getDependenciesList for a list of one-step dependent nodes of each node in the model.}

\item{\code{getDependenciesList(returnNames = TRUE, sort = TRUE)}}{Returns a list of all dependent neighbor relationships.  Each list element gives the one-step dependencies of one vertex, and the element name is the vertex label (integer ID or character node name)

Arguments:

returnNames: If TRUE (default), list names and element contents are returns as character node names, e.g. 'x[1]'.  If FALSE, everything is returned using graph IDs, which are unique integer labels for each node.

sort: If TRUE (default), each list element is returned in topologically sorted order.  If FALSE, they are returned in arbitrary order.

Details: This provides a fairly raw representation of the graph (model) structure that may be useful for inspecting what NIMBLE has created from model code.}

\item{\code{getDimension(node, params = NULL, valueOnly = is.null(params) &&
  !includeParams, includeParams = !is.null(params))}}{Determines the dimension of the value and/or parameters of the node

Arguments:

node: A character vector specifying a single node

params: an optional character vector of names of parameters for which dimensions are desired (possibly including 'value' and alternate parameters)

valueOnly: a logical indicating whether to only return the dimension of the value of the node

includeParams: a logical indicating whether to return dimensions of parameters. If TRUE and 'params' is NULL then dimensions of all parameters, including the dimension of the value of the node, are returned

Details: The return value is a numeric vector with an element for each parameter/value requested.}

\item{\code{getDistribution(nodes)}}{Returns the names of the distributions for the requested node or nodes

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.}

\item{\code{getDownstream(...)}}{Identical to getDependencies(..., downstream = TRUE)

Details: See documentation for member method getDependencies.}

\item{\code{getLogProb(nodes)}}{See `help(getLogProb)`}

\item{\code{getNodeNames(determOnly = FALSE, stochOnly = FALSE, includeData = TRUE,
  dataOnly = FALSE, includeRHSonly = FALSE, topOnly = FALSE,
  latentOnly = FALSE, endOnly = FALSE, returnType = "names",
  returnScalarComponents = FALSE)}}{Returns a character vector of all node names in the model, in topologically sorted order.  A variety of logical arguments allow for flexible subsetting of all model nodes.

Arguments:

determOnly: Logical argument specifying whether to return only deterministic nodes.  Default is FALSE.

stochOnly: Logical argument specifying whether to return only stochastic nodes.  Default is FALSE.

includeData: Logical argument specifying whether to include 'data' nodes (set via the member method setData).  Default is TRUE.

dataOnly: Logical argument specifying whether to return only 'data' nodes.  Default is FALSE.

includeRHSonly: Logical argument specifying whether to include right-hand-side-only nodes (model nodes which never appear on the left-hand-side of ~ or <- in the model code).  Default is FALSE.

topOnly: Logical argument specifying whether to return only top-level nodes from the hierarchical model structure.

latentOnly: Logical argument specifying whether to return only latent (mid-level) nodes from the hierarchical model structure.

endOnly: Logical argument specifying whether to return only end nodes from the hierarchical model structure.

returnType: Character argument specific type object returned. Options are 'names' (returns character vector) and 'ids' (returns numeric graph IDs for model)

returnScalar Componenets: Logical argument specifying whether multivariate nodes should return full node name (i.e. 'x[1:2]') or should break down into scalar componenets (i.e. 'x[1]' and 'x[2]')

Details: Multiple logical input arguments may be used simultaneously.  For example, `model$getNodeNames(endOnly = TRUE, dataOnly = TRUE)` will return all end-level nodes from the model which are designated as 'data'.}

\item{\code{getParam(node, param)}}{See `help(getParam)`                                                                                                                                                                                                                               }

\item{\code{getParents(nodes, omit = character(), self = FALSE, determOnly = FALSE,
  stochOnly = FALSE, includeData = TRUE, dataOnly = FALSE,
  includeRHSonly = FALSE, upstream = FALSE, immediateOnly = FALSE,
  returnType = "names", returnScalarComponents = FALSE)}}{ Returns a character vector of the nodes on which the input nodes depend, sorted topologically according to the model graph, by default recursing and stopping at stochastic parent nodes.  In the genealogical metaphor for a graphical model, this function returns the "parents" of the input nodes. In the river network metaphor, it returns upstream nodes.  By default, the returned nodes omit the input nodes, include only stochastic nodes, and stop at stochastic nodes.  Aditional input arguments provide flexibility in the values returned.

Arguments:

nodes: Character vector of node names, with index blocks allowed, and/or variable names, the parents of which will be returned.

omit: Character vector of node names, which will be omitted from the nodes returned.  In addition, parent nodes beyond these omitted nodes will not be returned.  The omitted nodes argument serves to stop the upward search through the hierarchical model structure, and excludes the specified node.

self: Logical argument specifying whether to include the input argument nodes in the return vector of dependent nodes.  Default is FALSE.

determOnly: Logical argument specifying whether to return only deterministic nodes.  Default is FALSE.

stochOnly: Logical argument specifying whether to return only stochastic nodes.  Default is TRUE.  If both determOnly and stochOnly are TRUE, no nodes will be returned.

includeData: Logical argument specifying whether to include 'data' nodes (set via nimbleModel or the setData method).  Default is TRUE.

dataOnly: Logical argument specifying whether to return only 'data' nodes.  Default is FALSE.

includeRHSonly: Logical argument specifying whether to include right-hand-side-only nodes (model nodes which never appear on the left-hand-side of ~ or <- in the model code).  These nodes are neither stochastic nor deterministic, but instead function as variable inputs to the model.  Default is FALSE.

upstream: Logical argument specifying whether the upward search through the hierarchical model structure should continue beyond the first and subsequent stochastic nodes encountered, hence returning all nodes upstream of the input nodes.  Default is FALSE.

immediateOnly: Logical argument specifying whether only the immediate parent nodes should be returned, even if they are deterministic.  If FALSE, getParents recurses and stops at stochastic nodes.  Default is FALSE.

returnType: Character argument specifying type of object returned. Options are 'names' (returns character vector) and 'ids' (returns numeric graph IDs for model).

returnScalarComponenets: Logical argument specifying whether multivariate nodes should be returned as full node names (i.e. 'x[1:2]') or as scalar componenets (i.e. 'x[1]' and 'x[2]').

Details: The upward search for dependent nodes propagates through deterministic nodes, but by default will halt at the first level of stochastic nodes encountered.  Use getParentsList for a list of one-step parent nodes of each node in the model.}

\item{\code{getParentsList(returnNames = TRUE, sort = TRUE)}}{Returns a list of all parent neighbor relationships.  Each list element gives the one-step parents of one vertex, and the element name is the vertex label (integer ID or character node name)

Arguments:

returnNames: If TRUE (default), list names and element contents are returns as character node names, e.g. 'x[1]'.  If FALSE, everything is returned using graph IDs, which are unique integer labels for each node.

sort: If TRUE (default), each list element is returned in topologically sorted order.  If FALSE, they are returned in arbitrary order.

Details: This provides a fairly raw representation of the graph (model) structure that may be useful for inspecting what NIMBLE has created from model code.}

\item{\code{getVarNames(includeLogProb = FALSE, nodes)}}{Returns the names of all variables in a model, optionally including the logProb variables

Arguments:

logProb: Logical argument specifying whether or not to include the logProb variables.  Default is FALSE.

nodes: An optional character vector supplying a subset of nodes for which to extract the variable names and return the unique set of variable names}

\item{\code{initializeInfo()}}{Provides more detailed information on which model nodes are not initialized.}

\item{\code{isBinary(nodes)}}{Determines whether one or more nodes represent binary random variables

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.}

\item{\code{isData(nodes)}}{Returns a vector of logical TRUE / FALSE values, corresponding to the 'data' flags of the input node names.

Arguments:

nodes: A character vector of node or variable names.

Details: The variable or node names specified is expanded into a vector of model node names.  A logical vector is returned, indicating whether each model node has been flagged as containing 'data'.}

\item{\code{isDeterm(nodes)}}{Determines whether one or more nodes are deterministic

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.}

\item{\code{isDiscrete(nodes)}}{Determines whether one or more nodes represent discrete random variables

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.}

\item{\code{isEndNode(nodes)}}{Determines whether one or more nodes are end nodes (nodes with no stochastic dependences)

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is logical vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.}

\item{\code{isMultivariate(nodes)}}{Determines whether one or more nodes represent multivariate nodes

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a logical vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.}

\item{\code{isStoch(nodes)}}{Determines whether one or more nodes are stochastic

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent node names, so the length of the output may be longer than that of the input.}

\item{\code{isTruncated(nodes)}}{Determines whether one or more nodes are truncated

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent nodes names, so the length of the output may be longer than that of the input}

\item{\code{isUnivariate(nodes)}}{Determines whether one or more nodes represent univariate random variables

Arguments:

nodes: A character vector specifying one or more node or variable names.

Details: The return value is a character vector with an element for each node indicated in the input. Note that variable names are expanded to their constituent nodes names, so the length of the output may be longer than that of the input}

\item{\code{newModel(data = NULL, inits = NULL, modelName = character(),
  replicate = FALSE, check = getNimbleOption("checkModel"),
  calculate = TRUE)}}{Returns a new R model object, with the same model definiton (as defined from the original model code) as the existing model object.

Arguments:

data: A named list specifying data nodes and values, for use in the newly returned model.  If not provided, the data argument from the creation of the original R model object will be used.

inits: A named list specifying initial values, for use in the newly returned model.  If not provided, the inits argument from the creation of the original R model object will be used.

modelName: An optional character string, used to set the internal name of the model object.  If provided, this name will propagate throughout the generated C++ code, serving to improve readability.

replicate: Logical specifying whether to replicate all current values and data flags from the current model in the new model.  If TRUE, then the data and inits arguments are not used.  Default value is FALSE.

check: A logical indicating whether to check the model object for missing or invalid values.  Default is given by the NIMBLE option 'checkModel', see help on 'nimbleOptions' for details. 

calculate: A logical indicating whether to run 'calculate' on the model; this will calculate all deterministic nodes and logProbability values given the current state of all nodes. Default is TRUE. For large models, one might want to disable this, but note that deterministic nodes, including nodes introduced into the model by NIMBLE, may be NA. 

Details: The newly created model object will be identical to the original model in terms of structure and functionality, but entirely distinct in terms of the internal values.}

\item{\code{resetData()}}{Resets the 'data' property of ALL model nodes to FALSE.  Subsequent to this call, the model will have no nodes flagged as 'data'.}

\item{\code{setData(..., warnAboutMissingNames = TRUE)}}{Sets the 'data' flag for specified nodes to TRUE, and also sets the value of these nodes to the value provided.  This is the exclusive method for specifying 'data' nodes in a model object.  When a 'data' argument is provided to 'nimbleModel()', it uses this method to set the data nodes.

Arguments:

...:  Arguments may be provided as named elements with numeric values or as character names of model variables.  These may be provided in a single list, a single character vector, or as multiple arguments.  When a named element with a numeric value is provided, the size and dimension must match the corresponding model variable.  This value will be copied to the model variable and any non-NA elements will be marked as data.  When a character name is provided, the value of that variable in the model is not changed but any currently non-NA values are marked as data.  Examples: setData('x', y = 1:10) will mark both x and y as data and will set the value of y to 1:10.  setData(list('x', y = 1:10)) is equivalent.  setData(c('x','y')) or setData('x','y') will mark both x and y as data.

Details: If a provided value (or the current value in the model when only a name is specified) contains some NA values, then the model nodes corresponding to these NAs will not have their value set, and will not be designated as 'data'.  Only model nodes corresponding to numeric values in the argument list elements will be designated as data.  Designating a deterministic model node as 'data' will result in an error.  Designating part of a multivariate node as 'data' and part as non-data (NA) will result in an error; multivariate nodes must be entirely data, or entirely non-data.}

\item{\code{setInits(inits)}}{Sets initial values (or more generally, any named list of value elements) into the model

Arguments:

inits: A named list.  The names of list elements must correspond to model variable names.  The elements of the list must be of class numeric, with size and dimension each matching the corresponding model variable.}

\item{\code{simulate(nodes, includeData = FALSE)}}{See `help(simulate)`                                                                                                                                                                                                                                }

\item{\code{topologicallySortNodes(nodes, returnType = "names")}}{Sorts the input list of node names according to the topological dependence ordering of the model structure.

Arguments:

nodes: A character vector of node or variable names, which is to be topologically sorted. Alternatively can be a numeric vector of graphIDs

returnType: character vector indicating return type. Choices are "names" or "ids"

Details: This function merely reorders its input argument.  This may be important prior to calls such as model$simulate(nodes) or model$calculate(nodes), to enforce that the operation is performed in topological order.}
}}
\examples{
code <- nimbleCode({
    mu ~ dnorm(0, 1)
    x[1] ~ dnorm(mu, 1)
    x[2] ~ dnorm(mu, 1)
})
Rmodel <- nimbleModel(code)
modelVars <- Rmodel$getVarNames()   ## returns 'mu' and 'x'
modelNodes <- Rmodel$getNodeNames()   ## returns 'mu', 'x[1]' and 'x[2]'
Rmodel$resetData()
Rmodel$setData(list(x = c(1.2, NA)))   ## flags only 'x[1]' node as data
Rmodel$isData(c('mu', 'x[1]', 'x[2]'))   ## returns c(FALSE, TRUE, FALSE)
}
\author{
Daniel Turek
}
\seealso{
\code{\link{initializeModel}}
}

