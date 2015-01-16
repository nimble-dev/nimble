mapsClass <- setRefClass(
    'mapsClass',
    
    fields = list(
        ## set directly from graphNodesList:
        nodeNames = 'ANY',
        graphIDs = 'ANY',
        types = 'ANY',
        
        ## vectors of nodeNames, representing different subsets of 'types'

        nodeNamesLHSall = 'ANY',                 ## names of all nodes with node functions, i.e. all node names *except* RHSonly nodes
        nodeNamesRHSonly = 'ANY',                ## names of all nodes appearing on RHS only
        
        ## nodeName_2_xxxx maps
        nodeName_2_graphID = 'ANY',                ## named vector of numeric graphIDs
        nodeName_2_nodeFunctionName = 'ANY',     ## named vector of character nodeFunctionNames
        
        ## graphID_2_xxxx maps
        graphID_2_nodeName = 		'ANY',              ## vector of character nodeNames
        graphID_2_nodeFunctionName ='ANY',      ## vector of character nodeFunctionNames

        ## varName2GraphID maps
        vars2GraphID_values = 		'ANY',
        vars2GraphID_functions =	'ANY',
        vars2LogProbName =			'ANY',
        vars2LogProbID = 			'ANY',
        
        logProbIDs_2_LogProbName =	'ANY',
        ## positions vectors of nodeNames (top, latent, end)
        isEndNode_byGID = 'ANY',

        ## Numeric Vectors containing the graphIDs's for the following node types
        top_IDs = 'ANY',
        latent_IDs = 'ANY',
        end_IDs = 'ANY'
        
        
    ),
    
    methods = list(
        setup = function() {},
        setup2 = function() {},
        setPositions = function() {}
    )
)


mapsClass$methods(setup2 = function(nodeNames, graphIDs, nodeFunctionNamesRaw, originNodeNamesRaw, types, graph, varInfo, nodeInfo) {
    types <<- types
    nodeNames <<- nodeNames
    graphIDs <<- graphIDs


    nodeNamesLHSall <<- nodeNames[types != 'RHSonly']
    nodeNamesRHSonly <<- nodeNames[types == 'RHSonly']

    ## nodeName_2_xxxx maps
    nodeName_2_graphID <<- graphIDs
    names(nodeName_2_graphID) <<- nodeNames
    boolRHSonly <- types != 'RHSonly'
    nodeName_2_nodeFunctionName <<- nodeFunctionNamesRaw[boolRHSonly] # need names
    names(nodeName_2_nodeFunctionName) <<- nodeNames[boolRHSonly]
    
    ## graphID_2_xxxx maps
    graphID_2_nodeName <<- nodeNames
    graphID_2_nodeFunctionName <<- nodeFunctionNamesRaw

    vars2GraphID_values <<- new.env()
    vars2GraphID_functions <<- new.env()
    vars2LogProbName <<- new.env()
    vars2LogProbID <<- new.env()

    isMultiVariateFunction <- grepl(':', nodeNames)
    strippedNodeNames <- removeIndexing(nodeNames)
    for(var in varInfo){
    	varName = var[['varName']]
    	if(var$nDim == 0){
    		vars2GraphID_values[[varName]] <<- nodeName_2_graphID[[varName]]
    		vars2GraphID_functions[[varName]] <<- nodeName_2_graphID[[varName]]
    		vars2LogProbName[[varName]] <<- as.character(NA)
    		vars2LogProbID[[varName]] <<- as.numeric(NA)
    	}
    	else{
	    	vars2GraphID_values[[varName]] <<- array(dim = var$maxs)
	    	vars2LogProbName[[varName]] <<- array(dim = var$maxs)
	    	vars2LogProbID[[varName]] <<- array(dim = var$maxs)
	    	storage.mode(vars2LogProbName[[varName]]) <<- 'character'
	    	nodeNames4Var <- nodeNames[strippedNodeNames == varName & !isMultiVariateFunction]
	    	var_GIDs = as.numeric(nodeName_2_graphID[nodeNames4Var])		#The only reason 'as.numeric' is used is to strip off names
	    	flatIndices = extractFlatIndices_wVarInfo(nodeNames4Var, var)
	    	vars2GraphID_values[[varName]][flatIndices] <<- var_GIDs
	    	vars2GraphID_functions[[varName]] <<- vars2GraphID_values[[varName]]
	    	
	    	nodeNames4Var <- nodeNames[strippedNodeNames == varName & isMultiVariateFunction]
	    	if(length(nodeNames4Var) > 0){
				uniqueNames <- unique(nodeName_2_nodeFunctionName[nodeNames4Var])
				for(funName in uniqueNames){
		    		var_GIDs = nodeName_2_graphID[funName]
					nodeNamesWithCall <- paste0(funName, "<- var_GIDs")
		    		eval(parse(text = nodeNamesWithCall[1])[[1]], envir = vars2GraphID_functions)		    	
		    	}
	    	}
	    }
    }
    assignLogProbName(nodeInfo, vars2LogProbName)
    logProbIDs_2_LogProbName <<- assignLogProbID(vars2LogProbName, vars2LogProbID)
    setPositions(graph)


})



#mapsClass$methods(setup = function(graphNodesList, graph, varInfo, nodeInfo) {
#    
#    nodeNames <<- names(graphNodesList)
#    graphIDs <<- unlist(lapply(graphNodesList, function(gn) gn$graphID), use.names = FALSE)
#    nodeFunctionNamesRaw <- lapply(graphNodesList, function(gn) gn$nodeFunctionName)
#    originNodeNamesRaw <- lapply(graphNodesList, function(gn) gn$originNodeName)
#    types <<- unlist(lapply(graphNodesList, function(gn) gn$type), use.names = FALSE)
#    nodeFunctionNames <<- unlist(nodeFunctionNamesRaw[types != 'RHSonly'], use.names = FALSE)
#
#    nodeNamesStoch <<- nodeNames[types == 'stoch']
#    nodeNamesDeterm <<- nodeNames[types == 'determ']
#    nodeNamesLHSinferred <<- nodeNames[types == 'LHSinferred']
#    nodeNamesLHSall <<- nodeNames[types != 'RHSonly']
#    nodeNamesRHSonly <<- nodeNames[types == 'RHSonly']
#    nodeNamesInModel <<- nodeNames[types != 'LHSinferred']
#    
#    ## nodeName_2_xxxx maps
#    nodeName_2_graphID <<- graphIDs
#    names(nodeName_2_graphID) <<- nodeNames
#    nodeName_2_type <<- types
#    names(nodeName_2_type) <<- nodeNames
#    nodeName_2_nodeFunctionName <<- unlist(nodeFunctionNamesRaw[types != 'RHSonly']) 
#    nodeName_2_originNodeName <<- unlist(originNodeNamesRaw)
#    
#    ## graphID_2_xxxx maps
#    graphID_2_nodeName <<- nodeNames
#    graphID_2_type <<- types
#    graphID_2_nodeFunctionName <<- unlist(nodeFunctionNamesRaw, use.names = FALSE)
#    graphID_2_originNodeName <<- unlist(originNodeNamesRaw, use.names = FALSE)
#    
#    vars2GraphID_values <<- new.env()
#    vars2GraphID_functions <<- new.env()
#    vars2LogProbName <<- new.env()
#    vars2LogProbID <<- new.env()
#    isMultiVariateFunction <- grepl(':', nodeNames)
#    strippedNodeNames <- removeIndexing(nodeNames)
#    for(var in varInfo){
#    	varName = var[['varName']]
#    	if(var$nDim == 0){
#    		vars2GraphID_values[[varName]] <<- nodeName_2_graphID[[varName]]
#    		vars2GraphID_functions[[varName]] <<- nodeName_2_graphID[[varName]]
#    		vars2LogProbName[[varName]] <<- as.character(NA)
#    		vars2LogProbID[[varName]] <<- as.numeric(NA)
#    	}
#    	else{
#	    	vars2GraphID_values[[varName]] <<- array(dim = var$maxs)
#	    	vars2LogProbName[[varName]] <<- array(dim = var$maxs)
#	    	vars2LogProbID[[varName]] <<- array(dim = var$maxs)
#	    	storage.mode(vars2LogProbName[[varName]]) <<- 'character'
#	    	nodeNames4Var <- nodeNames[strippedNodeNames == varName & !isMultiVariateFunction]
#	    	var_GIDs = as.numeric(nodeName_2_graphID[nodeNames4Var])		#The only reason 'as.numeric' is used is to strip off names
#	    	flatIndices = extractFlatIndices_wVarInfo(nodeNames4Var, var)
#	    	vars2GraphID_values[[varName]][flatIndices] <<- var_GIDs
#	    	vars2GraphID_functions[[varName]] <<- vars2GraphID_values[[varName]]
#	    	
#	    	nodeNames4Var <- nodeNames[strippedNodeNames == varName & isMultiVariateFunction]
#	    	if(length(nodeNames4Var) > 0){
#				uniqueNames <- unique(nodeName_2_nodeFunctionName[nodeNames4Var])
#				for(funName in uniqueNames){
#		    		var_GIDs = nodeName_2_graphID[funName]
#					nodeNamesWithCall <- paste0(funName, "<- var_GIDs")
#		    		eval(parse(text = nodeNamesWithCall[1])[[1]], envir = vars2GraphID_functions)		    	
#		    	}
#	    	}
#	    }
#    }
#    assignLogProbName(nodeInfo, vars2LogProbName)
#    logProbIDs_2_LogProbName <<- assignLogProbID(vars2LogProbName, vars2LogProbID)
#    setPositions(graph)
#    is_NodeFunction <<- rep(FALSE, length(graphIDs))
#    for(var in ls(vars2GraphID_functions))
#	    is_NodeFunction[unlist(as.list(vars2GraphID_functions[[var]]))] <<- TRUE
#	is_NodeFunction <<- which(is_NodeFunction)
#})

assignLogProbName <- function(nodeInfo, nodeName2LogProbMap){
	allLogProbNames <- as.character(unlist(lapply(nodeInfo, function(ni) ni$logProbNodeReplacedWithValues )))
	allNodeNames <- gsub('logProb_', '', allLogProbNames)
	allLogProbNameswQuotes <- paste0("'", allLogProbNames, "'")
	allNodeCalls <- paste(allNodeNames, " <- " , allLogProbNameswQuotes)
	for(call in allNodeCalls)
		eval(parse(text = call)[[1]], envir= nodeName2LogProbMap)
		
}

assignLogProbID <- function(nodeName2LogProbMap, nodeName2LogProbIDMap){
	numIDs = 0
	varNames <- ls(nodeName2LogProbMap)
	logProbIDs_2_LogProbName = character(0)
	for(vName in varNames){
		hasLogProb = which(!is.na(nodeName2LogProbMap[[vName]] ) )
		if(length(hasLogProb) > 0){
			theseIDs <- numIDs + 1:length(hasLogProb)
			nodeName2LogProbIDMap[[vName]][hasLogProb] <- theseIDs
			logProbIDs_2_LogProbName[theseIDs] <- nodeName2LogProbMap[[vName]][hasLogProb]
			numIDs = numIDs + length(hasLogProb)
		}
	}
	return(logProbIDs_2_LogProbName)
}

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
    
    top_IDs <<- top
    end_IDs <<- end
    latent_IDs <<- latent
    
	isEndNodebyGID <- rep(FALSE, length(nodeNames))
	isEndNodebyGID[end_IDs] <- TRUE
	isEndNode_byGID <<- isEndNodebyGID
})


