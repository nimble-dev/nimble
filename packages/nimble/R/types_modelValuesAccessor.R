
## CASE 3: modelValuesAccessorVector
## copy(modelValues, xxx, from.row=i, nodes) becomes:
## modelValues_nodes_accessors <- modelValuesAccessorVector(modelValues, nodes)
## copy(modelValues_nodes_accessors$row(i), xxx)
## modelValuesAccessorVector$getAccessors() returns a list of the modelValuesAccessor objects

modelValuesAccessor <- setRefClass(
    Class = 'modelValuesAccessor',
    fields = list(modelValues = 'ANY',
                  var         = 'ANY', 		#'character',
                  first       = 'ANY', 		#'numeric',
                  last        = 'ANY', 		#'numeric',
                  length	  = 'ANY' 		#'numeric'
    ),
    methods = list(toStr = function() paste0(var, '[', first, ':', last, ']'),
                   show  = function() cat(paste0(toStr(), '\n'))
    )
)

modelValuesAccessorVector <- setRefClass(
    Class = 'modelValuesAccessorVector',
    fields = list(modelValues = 'ANY',
                  nodes = 'ANY', 		#'character',
                  modelValuesAccessors = 'ANY',
                  length = 'ANY') ,		#'numeric'),
    methods = list(
        initialize = function(modelValues, nodeNames, logProb = FALSE, env = parent.frame()) {
            nodeNames <- nl_expandNodeNames(nodeNames, modelValues$symTab, env)  
            # expands nodeNames to fully indexed form, including expanding variables using the symbolTable
            
            #if(logProb){
            #    nodeNames <- c(nodeNames, makeLogProbName(nodeNames))
            #    nodeNames <- nl_removeNodeNamesNotInSymbolTable(nodeNames, modelValues$getSymbolTable())
            #}
            
            if(logProb){
	        	if(!inherits(modelValues$modelDef,'modelDefClass'))
	        		stop('calling logProb = TRUE on a modelValues object that was not built from a model')
	        	logProbNames <- modelValues$modelDef$nodeName2LogProbName(nodeNames)
        		nodeNames <- c(nodeNames, logProbNames)
        	}            
            varsAndFlatIndexRanges <- nl_createVarsAndFlatIndexRanges(nodeNames, modelValues$symTab)    # creates a list of variable names, and ranges of the flat index
            modelValues <<- modelValues
            nodes <<- nodeNames
            modelValuesAccessors <<- lapply(varsAndFlatIndexRanges, function(vafir) modelValuesAccessor(modelValues=modelValues, var=vafir$var, first=vafir$ind[1], last=vafir$ind[2], length = vafir$ind[2] - vafir$ind[1] + 1))
            len = 0
            for(mv in modelValuesAccessors)
            	len = len + mv$length
            length <<- len
        },
        getAccessors = function() return(modelValuesAccessors),
        show = function() cat(paste0('modelValuesAccessorVector: ', paste0(lapply(modelValuesAccessors, function(x) x$toStr()), collapse=', '), '\n'))
    )
)

length.modelValuesAccessorVector <- function(access)
	return(access$length)
