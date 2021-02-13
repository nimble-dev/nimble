mapsClass <- setRefClass(
    'mapsClass',
    
    fields = list(
        ## set directly from graphNodesList:
        nimbleGraph = 'ANY', ##graphNodeClass
        nodeNames = 'ANY',  ## like vertexID_2_nodeNames
        graphIDs = 'ANY',
        types = 'ANY',
        notStoch = 'ANY',
        
        ## vectors of nodeNames, representing different subsets of 'types'
        nodeNamesLHSall = 'ANY',                 ## names of all nodes with node functions, i.e. all node names *except* RHSonly nodes
        nodeNamesRHSonly = 'ANY',                ## names of all nodes appearing on RHS only
        elementNames = 'ANY',
        
        vertexID_2_nodeID = 'ANY',  ## new
        elementID_2_vertexID = 'ANY',
        
        ## nodeName_2_xxxx maps
        nodeName_2_graphID = 'ANY',                ## named vector of numeric graphIDs
        nodeName_2_nodeFunctionName = 'ANY',     ## named vector of character nodeFunctionNames
        
        ## graphID_2_xxxx maps
        graphID_2_nodeName = 		'ANY',              ## vector of character nodeNames
        graphID_2_logProbName =       'ANY',
        graphID_2_nodeFunctionName ='ANY',      ## vector of character nodeFunctionNames
        graphID_2_declID   =  'ANY',
        graphID_2_unrolledIndicesMatrixRow = 'ANY',
        
        ## varName2GraphID maps
        vars2GraphID_values = 		'ANY',
        vars2GraphID_functions =	'ANY',
        vars2GraphID_functions_and_RHSonly = 'ANY', 
        vars2ID_elements = 'ANY',
        vars2LogProbName =			'ANY',
##        vars2LogProbID = 			'ANY',
        
        ## positions vectors of nodeNames (top, latent, end)
        isEndNode_byGID = 'ANY',

        ## Numeric Vectors containing the graphIDs's for the following node types
        top_IDs = 'ANY',
        latent_IDs = 'ANY',
        end_IDs = 'ANY',

        edgesFrom = 'ANY',
        edgesTo = 'ANY',
        edgesParentExprID = 'ANY',
        edgesFrom2To = 'ANY',
        edgesFrom2ParentExprID = 'ANY'
        
        
    ),
    
    methods = list(
        setPositions3 = function() {}
    )
)

assignLogProbName <- function(nodeInfo, nodeName2LogProbMap){
	allLogProbNames <- as.character(unlist(lapply(nodeInfo, function(ni) ni$logProbNodeReplacedWithValues )))
	allNodeNames <- gsub('logProb_', '', allLogProbNames)
	allLogProbNameswQuotes <- paste0("'", allLogProbNames, "'")
	allNodeCalls <- paste(allNodeNames, " <- " , allLogProbNameswQuotes)
	for(call in allNodeCalls)
		eval(parse(text = call)[[1]], envir= nodeName2LogProbMap)
		
}

mapsClass$methods(setPositions3 = function(graph) { ## graph not actually used any more!
    ## new version to work with XXX3 system from modelDefClass
    ## determine who has any stochastic dependents (descendents)

    boolAnyStochDep <- nimbleGraph$anyStochDependencies()
    boolAnyStochParent <- nimbleGraph$anyStochParents()

    ## end nodes have no stochastic dependents
    ## top nodes have no stochastic ancestor
    ## latent nodes have a stochastic descendent and ancestor

    top_IDs <<- which(!boolAnyStochParent)
    end_IDs <<- which(!boolAnyStochDep)
    latent_IDs <<- which(boolAnyStochParent & boolAnyStochDep)

    isEndNode_byGID <<- rep(FALSE, length(nodeNames))
    isEndNode_byGID[end_IDs] <<- TRUE

    NULL
})



