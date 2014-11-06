## CASE 1: nodeFunctionVector
## simulate(model, nodes) becomes:
## model_nodes <- nodeFunctionVector(model, nodes)
## simulate(model_nodes)
## model_nodes$getNodeFunctions() returns a list of the nfRefClassObjets underlying the nodeFunctions

nodeFunctionVector <- setRefClass(
    Class = 'nodeFunctionVector',
    fields = list(model = 'ANY',
                  nodes = 'ANY'),
               #   nodeFunctionRefClassObjects = 'ANY'),
    methods = list(
        initialize = function(model, nodeNames, excludeData = FALSE, env = parent.frame()) {
			if(inherits(nodeNames, 'numeric'))
				nodeFunctionNames <- unique(model$modelDef$maps$graphID_2_nodeFunctionName[sort(nodeNames)])
			else if(inherits(nodeNames, 'character'))
        		nodeFunctionNames <- model$expandNodeNames(nodeNames, sort = TRUE)
        	if(excludeData){
        		nodeIsData <- model$isData(nodeFunctionNames)
        		nodeFunctionNames[!nodeIsData]
        	}
            model <<- model
            nodes <<- nodeFunctionNames
        },
        show = function() cat(paste0('nodeFunctionVector: ', paste(nodes, collapse=', '), '\n'))
    )
)

getNodeFxnPtrs <- function(cModel){	
    lapply( cModel$nodes, `[[`, '.basePtr' )
    ## output <- list()
    ## rNodeFxnVecNames = names(cModel$nodes)						
    ## for(vN in rNodeFxnVecNames){
    ##     output[[vN]] <- cModel$nodes[[vN]]$.basePtr		
    ## }
    ## return(output)
}
