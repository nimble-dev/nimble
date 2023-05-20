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


# The roxygen for modelBaseClass$getConditionallyIndependentSets refers users to
# here:

#' Get a list of conditionally independent sets of nodes in a nimble model
#'
#' A conditionally independent set of nodes is such that the joint probability
#' (density) of nodes in the set will not change even if any non-given
#' node outside the set is changed. Default given nodes are data nodes and
#' parameter nodes (aka "top-level" nodes, i.e. nodes with no parent nodes), but
#' this can be controlled.
#'
#' @param model A nimble model object (uncompiled or compiled), such as returned
#'   by \code{nimbleModel}.
#'
#' @param nodes A vector of stochastic node names (or their graph IDs) to split
#'   into conditionally independent sets, conditioned on the \code{givenNodes}.
#'   If \code{unknownAsGiven=FALSE}, the \code{nodes} are the starting nodes
#'   from which conditionally independent sets of nodes should be found,
#'   possibly including additional nodes not included in the \code{nodes}
#'   argument. If \code{nodes} is omitted, the default will be all latent nodes
#'   (defined as stochastic nodes that are not data and have at least one
#'   stochastic parent node, possibly with determinstic nodes in-between) that
#'   are a parent of a \code{givenNode} (either provided or default). Note that
#'   this will omit latent states that have no hyperparameters. An example is
#'   the first latent state in some state-space (time-series) models, which is
#'   sometimes declared with known prior.
#'
#' @param givenNodes A vector of node names or their graph IDs that should be
#'   considered as fixed (given) and hence can be conditioned on. If omitted,
#'   the default will be all data nodes and all parameter nodes, the latter
#'   defined as nodes with no stochastic parent nodes (skipping over
#'   deterministic parent nodes). See \code{endAsGiven} for a variant on
#'   defaults.
#'
#' @param omit A vector of node names or their graph IDs that should be omitted
#'   and should block further graph exploration.
#'
#' @param explore The method of graph exploration, which may corresond to what
#'   the \code{nodes} argument represents. For "down", graph exploration starts
#'   only down (towards descendants) from \code{nodes}. For "up", graph
#'   exploration starts only up (towards ancestors) from \code{nodes}. For
#'   "both" (the default and normal setting), both directions are explored.
#'
#' @param unknownAsGiven Logical for whether a model node not in \code{nodes} or
#'   \code{givenNodes} should be treated as given (default = TRUE). Otherwise
#'   (and by default) such a node may be grouped into a conditionally
#'   independent set, resulting in more output nodes than input \code{nodes}.
#'
#' @param returnType Either "names" for returned nodes to be node names or
#'   "ids" for returned nodes to be graph IDs.
#'
#' @param returnScalarComponents If FALSE (default), multivariate nodes are
#'   returned as full names (e.g. \code{x[1:3]}). If TRUE, they are returned as
#'   scalar elements (e.g. \code{x[1]}, \code{x[2]}, \code{x[3]}).
#'
#' @param endAsGiven If TRUE, end nodes (defined as nodes with stochastic
#'   parents but no stochastic children, skipping through deterministic nodes)
#'   are included in the default for \code{givenNodes}.
#'
#' @author Perry de Valpine
#'
#' @details This function returns sets of conditionally independent nodes.
#'   Multiple input \code{nodes} might be in the same set or different sets.
#'
#' The \code{nodes} input and the returned sets include only stochastic nodes
#' because conditional independence is a property of random variables.
#' Deterministic nodes are considered in determining the sets. \code{givenNodes}
#' may contain stochastic or deterministic nodes.
#'
#' @return List of nodes that are in conditionally independent sets. With each
#'   set, nodes are returned in topologically sorted order. The sets themselves
#'   are returned in topologically sorted order of their first nodes.
#'
#' Other nodes (not in \code{nodes}) may be included in the output if
#'   \code{unknownAsGiven=FALSE}.
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
#
# Removed from roxygen draft for getConditionallyIndependentSets:
#' @seealso There is a non-exported function
#'   \code{nimble:::testConditionallyIndependentSets(model, sets, initialize =
#'   TRUE)} that tests whether the conditional independence of sets is valid. It
#'   should be the case that
#'   \code{nimble:::testConditionallyIndependentSets(model,
#'   getConditionallyIndependentSets(model), initialize = TRUE)} returns
#'   \code{TRUE}.
#'
#
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
