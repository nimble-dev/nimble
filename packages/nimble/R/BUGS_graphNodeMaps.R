mapsClass <- setRefClass(
    'mapsClass',
    
    fields = list(
        ## set directly from graphNodesList:
        nodeNames = 'ANY',
        graphIDs = 'ANY',
        nodeFunctionNames = 'ANY',
        types = 'ANY',
        
        ## vectors of nodeNames, representing different subsets of 'types'
        nodeNamesStoch = 'ANY',                  ## names of all LHS declared 'stoch' nodes
        nodeNamesDeterm = 'ANY',                 ## names of all LHS declared 'determ' nodes
        nodeNamesLHSinferred = 'ANY',            ## names of all nodes inferred from LHS multivariate stochastic distributions
        nodeNamesLHSall = 'ANY',                 ## names of all nodes with node functions, i.e. all node names *except* RHSonly nodes
        nodeNamesRHSonly = 'ANY',                ## names of all nodes appearing on RHS only
        nodeNamesInModel = 'ANY',                ## names of all nodes in the model; everything *except* LHSinferred (graph-only) nodes
        
        ## nodeName_2_xxxx maps
        nodeName_2_graphID = 'ANY',                ## named vector of numeric graphIDs
        nodeName_2_type = 'ANY',                 ## named vector of character types
        nodeName_2_nodeFunctionName = 'ANY',     ## named vector of character nodeFunctionNames
        nodeName_2_originNodeName = 'ANY',       ## named vector of character nodeNames
        
        ## graphID_2_xxxx maps
        graphID_2_nodeName = 'ANY',              ## vector of character nodeNames
        graphID_2_type = 'ANY',                  ## vector of character types
        graphID_2_nodeFunctionName = 'ANY',      ## vector of character nodeFunctionNames
        graphID_2_originNodeName = 'ANY',        ## vector of character nodeNames
        
        ## varName2GraphID maps
        vars2GID = 'ANY',
        
        ## positions vectors of nodeNames (top, latent, end)
        nodeNamesTop = 'ANY',
        nodeNamesLatent = 'ANY',
        nodeNamesEnd = 'ANY'
    ),
    
    methods = list(
        setup = function() {},
        setPositions = function() {}
    )
)



mapsClass$methods(setup = function(graphNodesList, graph, varInfo) {
    
    nodeNames <<- names(graphNodesList)
    graphIDs <<- unlist(lapply(graphNodesList, function(gn) gn$graphID), use.names = FALSE)
    nodeFunctionNamesRaw <- lapply(graphNodesList, function(gn) gn$nodeFunctionName)
    originNodeNamesRaw <- lapply(graphNodesList, function(gn) gn$originNodeName)
    types <<- unlist(lapply(graphNodesList, function(gn) gn$type), use.names = FALSE)
    nodeFunctionNames <<- unlist(nodeFunctionNamesRaw[types != 'RHSonly'], use.names = FALSE)

    nodeNamesStoch <<- nodeNames[types == 'stoch']
    nodeNamesDeterm <<- nodeNames[types == 'determ']
    nodeNamesLHSinferred <<- nodeNames[types == 'LHSinferred']
    nodeNamesLHSall <<- nodeNames[types != 'RHSonly']
    nodeNamesRHSonly <<- nodeNames[types == 'RHSonly']
    nodeNamesInModel <<- nodeNames[types != 'LHSinferred']
    
    ## nodeName_2_xxxx maps
    nodeName_2_graphID <<- graphIDs
    names(nodeName_2_graphID) <<- nodeNames
    nodeName_2_type <<- types
    names(nodeName_2_type) <<- nodeNames
    nodeName_2_nodeFunctionName <<- unlist(nodeFunctionNamesRaw[types != 'RHSonly']) 
    nodeName_2_originNodeName <<- unlist(originNodeNamesRaw)
    
    ## graphID_2_xxxx maps
    graphID_2_nodeName <<- nodeNames
    graphID_2_type <<- types
    graphID_2_nodeFunctionName <<- unlist(nodeFunctionNamesRaw, use.names = FALSE)
    graphID_2_originNodeName <<- unlist(originNodeNamesRaw, use.names = FALSE)
    
    vars2GID <<- new.env()
    strippedNodeNames <- removeIndexing(nodeNames)
    for(var in varInfo){
    	varName = var[['varName']]
    	if(var$nDim == 0)	vars2GID[[varName]] <<- nodeName_2_graphID[[varName]]
    	else{
	    	vars2GID[[varName]] <<- array(dim = var$maxs)
	    	nodeNames4Var <- nodeNames[strippedNodeNames == varName]
	    	var_GIDs = nodeName_2_graphID[nodeNames4Var]
	    	flatIndices = extractFlatIndices_wVarInfo(nodeNames4Var, var)
	    	vars2GID[[varName]][flatIndices] <<- var_GIDs
	    	}
    }
    
    setPositions(graph)
})


mapsClass$methods(setPositions = function(graph) {
    
    non_top <- numeric(0)
    end <- numeric(0)
    
    for(id in graphIDs) {
        type <- types[id]
        if(type == 'RHSonly')         next
        if(type == 'LHSinferred')     next
        dep_ids <- gd_getDependencies_IDs(graph=graph, maps=.self, nodes=id, omit=numeric(0), downstream=FALSE)
        dep_ids <- setdiff(dep_ids, id)
        dep_types <- types[dep_ids]
        if(type == 'stoch')     non_top <- unique(c(non_top, dep_ids))
        if(!any(dep_types == 'stoch'))     end <- c(end, id)
    }
    
    top <- setdiff(graphIDs, non_top)
    latent <- setdiff(non_top, end)
    
    nodeNamesTop <<- nodeNames[top]
    nodeNamesLatent <<- nodeNames[latent]
    nodeNamesEnd <<- nodeNames[end]
})


