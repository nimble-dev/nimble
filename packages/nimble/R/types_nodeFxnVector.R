## CASE 1: nodeFunctionVector
## simulate(model, nodes) becomes:
## model_nodes <- nodeFunctionVector(model, nodes)
## simulate(model_nodes)
## model_nodes$getNodeFunctions() returns a list of the nfRefClassObjets underlying the nodeFunctions

nodeFunctionVector <- setRefClass(
    Class = 'nodeFunctionVector',
    fields = list(
        model = 'ANY',
        gids = 'ANY',
        indexingInfo = 'ANY'
                  ),
               #   nodes = 'ANY'),
               #   nodeFunctionRefClassObjects = 'ANY'),
    methods = list(
        initialize = function(model, nodeNames, excludeData = FALSE, sortUnique = TRUE, errorContext = "") { ##env = parent.frame()) {
            model <<- model
            if(length(nodeNames) == 0) {
            	gids <<- numeric(0)
                indexingInfo <<- list(declIDs = integer(), rowIndices = integer())
            } else {
                ## This is an old case no longer used
                ## if(is.numeric(nodeNames)) { 		#In case we start with graph ids instead of names
                ##     message("FOUND A CASE WHERE NODEFUNCTIONVECTOR IS CREATED WITH NUMERIC NODENAMES")
                ##     temp_gids <- unique(sort(nodeNames, FALSE),
                ##                                   FALSE,
                ##                                   FALSE,
                ##                                   NA) 
    	        ## } else {
                    if(sortUnique) {
                        temp_gids <- unique(sort(model$modelDef$nodeName2GraphIDs(nodeNames), FALSE),
                                            FALSE,
                                            FALSE,
                                            NA) ##unique( sort(model$modelDef$nodeName2GraphIDs(nodeNames)))
                        if(excludeData == TRUE)
                            temp_gids <- temp_gids[!model$isDataFromGraphID(temp_gids)]
                    } else {
                        temp_gids <- model$modelDef$nodeName2GraphIDs(nodeNames, unique = FALSE)
                        if(length(temp_gids) != length(nodeNames)) stop(paste0("In nodeFunctionVector from a case like model$doSomething(nodes[i]) where nodes may not contain well-defined node names.  Context is ", errorContext))
                        if(excludeData) {
                            if(sum(model$isDataFromGraphID(temp_gids)) > 0) ## message uses includeData instead of excludeData b/c that's what users see as argument
                                warning(paste0("In nodeFunctionVector, usage is of the form model$doSomething(nodes[i]) where nodes includes data nodes, but includeData is FALSE.  Set includeData = TRUE if you need to include data nodes in this case.  Context is ", errorContext))
                        }
                    }
##                }
                gids <<- temp_gids ## old nodeFun system
                indexingInfo <<- model$modelDef$graphIDs2indexedNodeInfo(temp_gids) ## new nodeFun system
            }
        },
        getNodeNames = function(){ ## not used anywhere. provided only for debugging/inspection
        	model$expandNodeNames(gids)	
        },
        show = function() {
            cat('nodeFunctionVector: \n')
            print(cbind(indexingInfo$declIDs, indexingInfo$unrolledIndicesMatrixRows))
        }
    )
)

