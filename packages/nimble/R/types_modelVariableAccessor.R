## CASE 2: modelVariableAccessorVector
## copy(model1, xxx, nodes) becomes:
## model1_nodes_accessors <- modelVariableAccessorVector(model1, nodes)
## copy(model1_nodes_accessors, xxx)
## model1_nodes_accessors$getAccessors() returns a list of the modelVariableAccessor objects

modelVariableAccessor <- setRefClass(
    Class = 'modelVariableAccessor',
    fields = list(model = 'ANY',
                  var   = 'ANY', 		# 'character',
                  first = 'ANY', 		#'numeric',
                  last  = 'ANY', 		#'numeric',
                  length = 'ANY' 		#'numeric'
    ),
    methods = list(toStr = function() paste0(var, '[', first, ':', last, ']'),
                   show  = function() cat(paste0(toStr(), '\n'))
    )
)

modelVariableAccessorVector <- setRefClass(
    Class = 'modelVariableAccessorVector',
    fields = list(model = 'ANY',
                  nodes ='ANY', 		#'character',
                  modelVariableAccessors = 'ANY',
                  length = 'ANY' ),		#'numeric'),
    methods = list(
        initialize = function(model, nodeNames, logProb = FALSE, env = parent.frame()) {

       #     nodeNames <- nl_expandNodeNames(nodeNames, model$getSymbolTable(), env) 
       		nodeNames <- model$expandNodeNames(nodeNames, returnScalarComponents = TRUE)
       # 	expands nodeNames to fully indexed form, including expanding variables using the symbolTable
    
       #    if(logProb){
       #        nodeNames <- c(nodeNames, makeLogProbName(nodeNames))
       #        nodeNames <- nl_removeNodeNamesNotInSymbolTable(nodeNames, model$getSymbolTable())
       #    }
        	
			if(logProb){
	        	logProbNames <- model$modelDef$nodeName2LogProbName(nodeNames)
        		nodeNames <- c(nodeNames, logProbNames)
        	}
            varsAndFlatIndexRanges <- nl_createVarsAndFlatIndexRanges(nodeNames, model$getSymbolTable())  
            	# creates a list of variable names, and ranges of the flat index
            model <<- model
            nodes <<- nodeNames
            modelVariableAccessors <<- lapply(varsAndFlatIndexRanges, function(vafir) modelVariableAccessor(model=model, var=vafir$var, first=vafir$ind[1], last=vafir$ind[2], length = vafir$ind[2] - vafir$ind[1] + 1))
            
            
            len = 0
            for(mv in modelVariableAccessors)
            	len = len + mv$length
            length <<- len
        },
        getAccessors = function() return(modelVariableAccessors),
        show = function() cat(paste0('modelVariableAccessorVector: ', paste0(lapply(modelVariableAccessors, function(x) x$toStr()), collapse=', '), '\n'))
    )
)

# This function allows you to just call "length(access)", rather than access$length
# Motivation for this is to make it easier to generically use accessors
length.modelVariableAccessorVector <- function(access)
	return(access$length)
