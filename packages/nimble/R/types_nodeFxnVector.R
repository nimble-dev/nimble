## CASE 1: nodeFunctionVector
## simulate(model, nodes) becomes:
## model_nodes <- nodeFunctionVector(model, nodes)
## simulate(model_nodes)
## model_nodes$getNodeFunctions() returns a list of the nfRefClassObjets underlying the nodeFunctions

nodeFunctionVector <- setRefClass(
    Class = 'nodeFunctionVector',
    fields = list(model = 'ANY',
                  nodes = 'ANY', 		#'character',
                  nodeFunctionRefClassObjects = 'ANY'),
    methods = list(
        initialize = function(model, nodeNames, env = parent.frame()) {
            nodeNames <- model$expandNodeNames(nodeNames, env)  # expands nodeNames to fully indexed form, including expanding variables using the symbolTable, and only keep ones in the model
            nl_checkNodeNamesInModel(model, nodeNames)           # checks that all nodeNames are present in model
            nodeNames <- model$topologicallySortNodes(nodeNames)    # topological sort of nodeNames
            nodeFunctionNames <- unique(model$getMaps('nodeName_2_nodeFunctionName')[nodeNames])
            model <<- model
            nodes <<- nodeFunctionNames
         ##   nodeFunctionRefClassObjects <<- lapply(model$nodeFunctions[nodeFunctionNames], nf_getRefClassObject)
        },
        getNodeFunctions = function() {  ## We haven't actually needed this.  This will NOT work if the model is a Cmodel 
            if(inherits(nodeFunctionRefClassObjects, 'uninitializedField') ) {
                nodeFunctionRefClassObjects <<- lapply(model$nodeFunctions[nodes], nf_getRefClassObject)
            }
            return(nodeFunctionRefClassObjects)
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
