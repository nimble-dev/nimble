
## CASE 3: modelValuesAccessorVector
## copy(modelValues, xxx, from.row=i, nodes) becomes:
## modelValues_nodes_accessors <- modelValuesAccessorVector(modelValues, nodes)
## copy(modelValues_nodes_accessors$row(i), xxx)
## modelValuesAccessorVector$getAccessors() returns a list of the modelValuesAccessor objects

modelValuesAccessor <- setRefClass(
    Class = 'modelValuesAccessor',
    fields = list(modelValues = 'ANY',
    			   id = 'ANY'
   #               var         = 'ANY', 		#'character',
   #               first       = 'ANY', 		#'numeric',
   #               last        = 'ANY', 		#'numeric',
   #               length	  = 'ANY' 		#'numeric'
    ),
    methods = list(toStr = function() paste0(var, '[', first, ':', last, ']'),
                   show  = function() cat(paste0(toStr(), '\n'))
    )
)

makeSingleModelValuesAccessor <- function(ids, modelValues){
	return(modelValuesAccessor(id = id, modelValues = modelValues))
}
modelValuesAccessorVector <- setRefClass(
    Class = 'modelValuesAccessorVector',
    fields = list(modelValues = 'ANY',
                  gids = 'ANY', 		#'character',
                  length = 'ANY'		#'numeric'),
                  ),
    methods = list(
        initialize = function(modelValues, nodeNames, logProb = FALSE, env = parent.frame()) {
        	
        	if(logProb){
        		if(!inherits(modelValues$modelDef,'modelDefClass'))
	        		stop('calling logProb = TRUE on a modelValues object that was not built from a model')
        		nodeNames <- c(nodeNames, modelValues$modelDef$nodeName2LogProbName(nodeNames))
        	}
        	
        	gids <<- modelValues$expandNodeNames(nodeNames, returnType = 'ids')
        	modelValues <<- modelValues
        	length <<- length(gids)
        	        	
        },
       getNodeNames = function() {
       	modelValues$expandNodeNames(gids)
       	},
        getSingleValue_fromGID = function(accessID, row){
        	valueName <- modelValues$expandNodeNames(gids[accessID])
        	parseInfo <- parse(text = valueName)[[1]]
        	if(length(parseInfo) == 1){
        		vName = as.character(parseInfo)
        		index = 1
        	}
        	else{
        		vName = as.character(parseInfo[[2]])
        		index = parseInfo[[3]]
        	}
        	return(modelValues[vName, row][index])
        },
        setSingleValue_fromGID = function(value, accessID, row){
        	valueName <- modelValues$expandNodeNames(gids[accessID])
        	parseInfo <- parse(text = valueName)[[1]]
        	if(length(parseInfo) == 1){
        		vName = as.character(parseInfo)
        		index = 1
        	}
        	else{
        		vName = as.character(parseInfo[[2]])
        		index = parseInfo[[3]]
        	}
        	modelValues[vName, row][index] <<- value
        },
        show = function() cat(paste0('modelValuesAccessorVector: ', paste0(lapply(modelValuesAccessors, function(x) x$toStr()), collapse=', '), '\n'))
    )
)

length.modelValuesAccessorVector <- function(access)
	return(access$length)
