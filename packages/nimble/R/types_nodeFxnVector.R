## CASE 1: nodeFunctionVector
## model$simulate(nodes) becomes:
## model_nodes <- nodeFunctionVector(model, nodes)
## simulate(model_nodes)
## model_nodes$getNodeFunctions() returns a list of the nfRefClassObjets underlying the nodeFunctions

nodeFunctionVector_WithDerivsOutputNodes <- function(model,
                                                 calcNodes,
                                                 excludeData,
                                                 sortUnique) {
  nimDerivsInfo <- nimDerivsInfoClass(calcNodes = calcNodes,
                                      thisModel = model,
                                      case = "outputOnly")
  
  ## Make one dummy node that can be used for set_CppAD_tape_info_for_model
  ## It will never be called from this object
  NFV <- nodeFunctionVector(model = model,
                            nodeNames = calcNodes,
                            excludeData = excludeData,
                            sortUnique = sortUnique)
  class(NFV) <- "nodeFunctionVector_nimDerivs"
  NFV$nimDerivsInfo <- nimDerivsInfo
  NFV
}

nodeFunctionVector_DerivsModelUpdateNodes <- function(model,
                                                      updateNodes = NULL,
                                                      constantNodes = NULL) {
  nimDerivsInfo <- nimDerivsInfoClass(updateNodes = updateNodes,
                                      constantNodes = constantNodes,
                                      thisModel = model,
                                      case = "updateOnly")
  ## Make one dummy node that can be used for set_CppAD_tape_info_for_model
  ## It will never be called from this object
  ## Using the first of updateNodes or constantNodes did not work because
  ## it could be a RHSonly node, which nodes not have a nodeFunction,
  ## which manifests as a crash if not caught.  Instead now we use the first
  ## node in the model (by default, includeRHSnodes = FALSE).
  ## dummyNodeNames <- c(updateNodes, constantNodes)[1]
  dummyNodeNames <- model$getNodeNames()[1]
  if(isTRUE(is.na(dummyNodeNames))) dummyNodeNames <- character()
  if(isTRUE(is.null(dummyNodeNames))) dummyNodeNames <- character()
  NFV <- nodeFunctionVector(model = model,
                            nodeNames = dummyNodeNames,
                            excludeData = FALSE,
                            sortUnique = FALSE)
  class(NFV) <- "nodeFunctionVector_nimDerivs"
  NFV$nimDerivsInfo <- nimDerivsInfo
  NFV
}

nodeFunctionVector <-
    function(model,
             nodeNames,
             wrtNodes = NULL,
             excludeData = FALSE,
             sortUnique = TRUE,
             errorContext = "")
{
    if(!is.null(wrtNodes)){
        nimDerivsInfo <- nimDerivsInfoClass(wrtNodes = wrtNodes, calcNodes = nodeNames, thisModel = model)
    }
    else{
        nimDerivsInfo <- NULL
    }
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
    classLabel <- if(is.null(nimDerivsInfo))
                      "nodeFunctionVector"
                  else
                      "nodeFunctionVector_nimDerivs"
    structure(list(gids = gids,
                   indexingInfo = indexingInfo,
                   model = model,
                   nimDerivsInfo = nimDerivsInfo),
              class = classLabel)
}
