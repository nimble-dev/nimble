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
        getParents = function(nodes, omit = integer(), upstream = FALSE) {
          # nodes and omit must be IDs
          for(i in list(nodes, omit)) {
            if(length(i) > 0) {
              if(any(abs(i - round(i)) > 0.1)) {
                stop('caught something wrong on input to getParents', call. = FALSE)
              }
            }
          }
          .Call(C_getParents, graphExtPtr, nodes, omit, upstream)
        },
        getConditionallyIndependentSets = function(nodeIDs,
                                                   givenNodeIDs,
                                                   omitIDs = integer(),
                                                   startUp = TRUE,
                                                   startDown = TRUE) {
          .Call(nimble:::C_getConditionallyIndependentSets,
                graphExtPtr,
                nodeIDs,
                givenNodeIDs,
                omitIDs,
                startUp,
                startDown)
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
