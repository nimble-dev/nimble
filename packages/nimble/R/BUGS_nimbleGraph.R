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
                        cat('caught something wrong on input to setGraph')
                        browser()
                    }
                }
            }
            ## Some nodeFunctionIDs may be zero, if there really is no nodeFunction (e.g. RHSonly)
            ## But on the C++ side we need these IDs to be self, so there is something valid.
            boolZero <- nodeFunctionIDs == 0
            nodeFunctionIDs[boolZero] <- (1:length(nodeFunctionIDs))[boolZero]
            graphExtPtr <<- .Call('setGraph', edgesFrom, edgesTo, edgesFrom2ParentExprIDs, nodeFunctionIDs, types, names, numNodes)
        },
        anyStochDependencies = function() {
           .Call("anyStochDependencies",graphExtPtr)
        },
        anyStochParents = function() {
            .Call("anyStochParents",graphExtPtr)
        },
        getDependencies = function(nodes, omit = integer(), downstream = FALSE) {
            for(i in list(nodes, omit)) {
                if(length(i) > 0) {
                    if(any(abs(i - round(i)) > 0.1)) {
                        cat('caught something wrong on input to setGraph')
                        browser()
                    }
                }
            }

            .Call(C_getDependencies, graphExtPtr, nodes, omit, downstream)
        },
        getDependencyPathCountOneNode = function(node) {
            if(length(node) > 1)
                stop("getDependencyPathCountOneNode: argument 'node' should provide a single node.")
            .Call("getDependencyPathCountOneNode", graphExtPtr, node)
        }
    ))
