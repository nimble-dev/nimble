nimbleGraphClass <- setRefClass(
    "nimbleGraphClass",
    fields = list(
        graphExtPtr = 'ANY'
    ),
    methods = list(
        setGraph = function(edgesFrom, edgesTo, edgesFrom2ParentExprIDs, nodeFunctionIDs, types, names, numNodes) {
            edgesFrom2ParentExprIDs[ is.na(edgesFrom2ParentExprIDs) ] <- 0
            for(i in list(edgesFrom, edgesTo, edgesFrom2ParentExprIDs, numNodes)) {
                if(length(i) > 0) {
                    if(any(abs(i - round(i)) > 0.1)) {
                      stop('caught something wrong on input to setGraph', call. = FALSE)
                    }
                }
            }
            ## Some nodeFunctionIDs may be zero, if there really is no nodeFunction (e.g. RHSonly)
            ## But on the C++ side we need these IDs to be self, so there is something valid.
            boolZero <- nodeFunctionIDs == 0
            nodeFunctionIDs[boolZero] <- (1:length(nodeFunctionIDs))[boolZero]
            graphExtPtr <<- .Call(C_setGraph, edgesFrom, edgesTo, edgesFrom2ParentExprIDs, nodeFunctionIDs, types, names, numNodes)
        },
        anyStochDependencies = function() {
           .Call(C_anyStochDependencies,graphExtPtr)
        },
        anyStochParents = function() {
            .Call(C_anyStochParents,graphExtPtr)
        },
        getDependencies = function(nodes, omit = integer(), downstream = FALSE) {
            for(i in list(nodes, omit)) {
                if(length(i) > 0) {
                    if(any(abs(i - round(i)) > 0.1)) {
                        stop('caught something wrong on input to getDependencies', call. = FALSE)
                    }
                }
            }
            .Call(C_getDependencies, graphExtPtr, nodes, omit, downstream)
        },
        getParents = function(nodes, omit = integer(), upstream = FALSE, immediateOnly = FALSE) {
          # nodes and omit must be IDs
          for(i in list(nodes, omit)) {
            if(length(i) > 0) {
              if(any(abs(i - round(i)) > 0.1)) {
                stop('caught something wrong on input to getParents', call. = FALSE)
              }
            }
          }
          .Call(C_getParents, graphExtPtr, nodes, omit, upstream, immediateOnly)
        },
        getConditionallyIndependentSets = function(nodeIDs,
                                                   givenNodeIDs,
                                                   omitIDs = integer(),
                                                   startUp = TRUE,
                                                   startDown = TRUE,
                                                   unknownAsGiven = FALSE
                                                   ) {
          .Call(nimble:::C_getConditionallyIndependentSets,
                graphExtPtr,
                nodeIDs,
                givenNodeIDs,
                omitIDs,
                startUp,
                startDown,
                unknownAsGiven)
        },
        getDependencyPathCountOneNode = function(node) {
            if(length(node) > 1)
                stop("getDependencyPathCountOneNode: argument 'node' should provide a single node.")
            .Call(C_getDependencyPathCountOneNode, graphExtPtr, node)
        },
        getDependencyPaths = function(node) {
            if(length(node) > 1)
                stop("getDependencyPaths: argument 'node' should provide a single node.")
            .Call(C_getDependencyPaths, graphExtPtr, node)  
        }
    ))


# The following roxygen is basically redundant with the method documentation for
# modelBaseClass::getConditionallyIndependentSets. Not sure we need both.

#' Get a list of conditionally independent sets of nodes in a nimble model
#'
#' Conditionally independent sets of nodes are typically groups of latent states
#' whose joint probability (density) will not change even if any other non-fixed
#' node is changed. Default fixed nodes are data nodes and parameter nodes (with
#' no parent nodes), but this can be controlled.
#'
#' @param model A nimble model object (uncompiled or compiled).
#'
#' @param nodes A vector of node names or their graph IDs that are the starting
#'   nodes from which conditionally independent sets of nodes should be found.
#'   If omitted, the default will be all latent nodes, defined as stochastic
#'   nodes that are not data and have at least one stochastic parent node
#'   (possible with determinstic nodes in between). Note that this will omit
#'   latent states that have no hyperparameters. An example is the first latent
#'   state in some state-space (time-series) models, which is sometimes declared
#'   with known prior. See \code{type} because it relates to the interpretation
#'   of \code{nodes}.
#'
#' @param givenNodes A vector of node names or their graph IDs that should be
#'   considered as fixed and hence can be conditioned on. If omitted, the
#'   default will be all data nodes and all parameter nodes, the latter defined
#'   as nodes with no stochastic parent nodes (skipping over deterministic
#'   parent nodes).
#'
#' @param omit A vector of node names or their graph IDs that should be omitted
#'   and should block further graph exploration.
#'
#' @param inputType The method of graph exploration depends on what the \code{nodes}
#'   argument represents. For "\code{latent}", the input \code{nodes} are
#'   interpreted as latent states, from which both parent and descendent graph
#'   exploration should be done to find nodes in the same set (nodes that are
#'   NOT conditionally independent from each other). For "\code{param}", the input
#'   \code{nodes} are interpreted as parameters, so graph exploration begins
#'   from the top (input) and explores descendents. For "\code{data}", the input
#'   \code{nodes} are interpreted as data nodes, so graph exploration begins
#'   from the bottom (input) explores parent nodes.
#'
#' @param stochOnly Logical for whether only stochastic nodes should be returned
#'   (default = TRUE). If FALSE, both deterministic and stochastic nodes are
#'   returned.
#'
#' @param unknownAsGiven Logical for whether a model node not in \code{nodes} or
#'   \code{givenNodes} should be treated as given (default = FALSE). Otherwise
#'   (and by default) such a node may be grouped into a conditionally
#'   independent set.
#'
#' @param returnType Either \code{names} for returned nodes to be node names or
#'   \code{ids} for returned nodes to be graph IDs.
#'
#' @param returnScalarComponents If FALSE (default), multivariate nodes are
#'   returned as full names (e.g. \code{x[1:3]}). If TRUE, they are returned as
#'   scalar elements (e.g. \code{x[1]}, \code{x[2]}, \code{x[3]}).
#'
#' @author Perry de Valpine
#'
#' @details This function returns sets of conditionally independent nodes.
#'   Multiple input \code{nodes} might be in the same set or different sets, and
#'   other nodes (not in \code{nodes}) will be included.
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
#' @return List of nodes that are in conditionally independent sets. With each
#'   set, nodes are returned in topologically sorted order. The sets themselves
#'   are returned in topologically sorted order of their first nodes.
#'
#' @seealso There is a non-exported function
#'   \code{nimble:::testConditionallyIndependentSets(model, sets, initialize =
#'   TRUE)} that tests whether the conditional independence of sets is valid. It
#'   should be the case that
#'   \code{nimble:::testConditionallyIndependentSets(model,
#'   getConditionallyIndependentSets(model), initialize = TRUE)} returns
#'   \code{TRUE}.
#'
#' @export
getConditionallyIndependentSets <- function(model,
                                            nodes,
                                            givenNodes,
                                            omit = integer(),
                                            explore = c("both", "down", "up"),
                                            unknownAsGiven = TRUE,
                                            returnType = 'names',
                                            returnScalarComponents = FALSE,
                                            endAsGiven = FALSE) {
  explore <- match.arg(explore)

  # Make default givenNodes if not provided.
  # default to stochastic top nodes and data nodes.
  # NB We do not assume end nodes should be givenNodes.
  if(missing(givenNodes)) {
    givenNodeIDs <- c(model$getNodeNames(topOnly = TRUE, stochOnly = TRUE, returnType = 'ids'),
                      model$getNodeNames(dataOnly = TRUE, stochOnly = TRUE, returnType = 'ids'))
    if(endAsGiven)
      givenNodeIDs <- unique(c(givenNodeIDs, model$getNodeNames(endOnly = TRUE, returnType = 'ids')))
  } else {
    if(is.character(givenNodes))
        givenNodeIDs <- model$expandNodeNames(givenNodes, returnType = 'ids')
    else if(is.numeric(givenNodes))
        givenNodeIDs <- givenNodes
  }

  # Make default nodes (i.e. latents in a typical use case) if not provided.
  # Default to structurally latent nodes that are stochastic, not data, and
  # somewhere have a descendant that is a givenNode
  if(missing(nodes)) {
    nodeIDs <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE, includeData = FALSE, returnType = 'ids')
    allGivenParents <- model$getParents(givenNodeIDs, upstream = TRUE, returnType = 'ids')
    nodeIDs <- intersect(nodeIDs, allGivenParents)
  } else {
    if(is.character(nodes))
      nodeIDs <- model$expandNodeNames(nodes, returnType = 'ids')
    else
      nodeIDs <- nodes
  }

  if(!missing(nodes)) {
    if(missing(givenNodes))
      givenNodesIDs <- setdiff(givenNodeIDs, nodeIDs)
  }
  if(!missing(givenNodes)) {
    nodeIDs <- setdiff(nodeIDs, givenNodeIDs)
  }

  ## if(isTRUE(nimbleOptions("groupDetermWithGivenInCondIndSets"))) {
  ##   givenNodeIDs <- unique(c(givenNodeIDs, model$getDependencies(givenNodeIDs, determOnly = TRUE, self = FALSE, returnType = 'ids')))
  ## }
  if(is.character(omit)) {
    # This mimcs code in getDependencies.  I think it allows omit to include split nodes, whereas getNodeNames would not.
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
  if(explore == "down") startUp <- FALSE
  if(explore == "up")   startDown <- FALSE
  result <- model$modelDef$maps$nimbleGraph$getConditionallyIndependentSets(
    nodeIDs = nodeIDs,
    givenNodeIDs = givenNodeIDs,
    omitIDs = omitIDs,
    startUp = startUp,
    startDown = startDown,
    unknownAsGiven = unknownAsGiven)
  if(returnType == 'ids' && returnScalarComponents)
    warning("NIMBLE development warning: calling getConditionallyIndependentSets with returnType = ids and returnScalarComponents may not be meaningful.")
  result <- lapply(result,
                   function(IDs) {
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
