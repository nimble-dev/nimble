## CASE 1: nodeFunctionVector
## model$simulate(nodes) becomes:
## model_nodes <- nodeFunctionVector(model, nodes)
## simulate(model_nodes)
## model_nodes$getNodeFunctions() returns a list of the nfRefClassObjets underlying the nodeFunctions

nodeFunctionVector <-
    function(model,
             nodeNames,
             excludeData = FALSE,
             sortUnique = TRUE,
             errorContext = "")
{
    if(length(nodeNames) == 0) {
        gids <- numeric(0)
        indexingInfo <- list(declIDs = integer(), rowIndices = integer())
    } else {
        if(sortUnique) {
            temp_gids <-
                unique(
                    sort(model$modelDef$nodeName2GraphIDs(nodeNames)),
                    FALSE,
                    FALSE,
                    NA) 
            if(excludeData == TRUE)
                temp_gids <-
                    temp_gids[!model$isDataFromGraphID(temp_gids)]
        } else {
            temp_gids <-
                model$modelDef$nodeName2GraphIDs(nodeNames, unique = FALSE)
            if(length(temp_gids) != length(nodeNames))
                stop(paste0("In nodeFunctionVector from a case like model$doSomething(nodes[i]) where nodes may not contain well-defined node names.  Context is ",
                            errorContext))
            if(excludeData) {
                if(sum(model$isDataFromGraphID(temp_gids)) > 0)
                    ## message uses "includeData" instead of "excludeData" b/c that's what users see as argument
                    warning(paste0("In nodeFunctionVector, usage is of the form model$doSomething(nodes[i]) where nodes includes data nodes, but includeData is FALSE.  Set includeData = TRUE if you need to include data nodes in this case.  Context is ",
                                   errorContext))
            }
        }
        gids <- temp_gids
        indexingInfo <-
            model$modelDef$graphIDs2indexedNodeInfo(temp_gids)
    }
    ## Any modification to this list ordering needs to be updated in
    ## populateNodeFxnVectorNew_copyFromRobject in accessorClasses.cpp
    structure(list(gids = gids,
                   indexingInfo = indexingInfo,
                   model = model),
              class = "nodeFunctionVector")
}

