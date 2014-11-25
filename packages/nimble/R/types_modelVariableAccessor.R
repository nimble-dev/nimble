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
    fields = list(model = 	'ANY',
    			  gids = 	'ANY',
    			  l_gids = 	'ANY',
                  length = 	'ANY'
                  
                  #nodes ='ANY', 		#'character',
                  #modelVariableAccessors = 'ANY',
                   ),		#'numeric'),
    methods = list(
        initialize = function(model, nodeNames, logProb = FALSE, env = parent.frame()) {
        	gids <<- model$expandNodeNames(nodeNames, returnScalarComponents = TRUE, returnType = 'ids')
        	l_gids <<- numeric(0)
        	if(logProb)
        		l_gids <<- model$modelDef$nodeName2LogProbID(nodeNames)
            length <<- length(gids) + length(l_gids)
   			model <<- model
   
   			if(length < 1)
   				cat('warning: created a modelVariableAccessor of length less than one. Probably a mistake. Original nodeNames = ', nodeNames, '\n')
   
       	#	nodeNames <- model$expandNodeNames(nodeNames, returnScalarComponents = TRUE)
		#	if(logProb){
	    #    	logProbNames <- model$modelDef$nodeName2LogProbName(nodeNames)
        #		nodeNames <- c(nodeNames, logProbNames)
        #	}
#            varsAndFlatIndexRanges <- nl_createVarsAndFlatIndexRanges(nodeNames, model$getSymbolTable())  
            	# creates a list of variable names, and ranges of the flat index
#            model <<- model
#            nodes <<- nodeNames
#            modelVariableAccessors <<- lapply(varsAndFlatIndexRanges, function(vafir) modelVariableAccessor(model=model, var=vafir$var, first=vafir$ind[1], last=vafir$ind[2], length = vafir$ind[2] - vafir$ind[1] + 1))
            
            
 #           len = 0
 #           for(mv in modelVariableAccessors)
 #           	len = len + mv$length
        },
        getAccessors = function() return(modelVariableAccessors),
        show = function() cat(paste0('modelVariableAccessorVector: ', paste0(lapply(modelVariableAccessors, function(x) x$toStr()), collapse=', '), '\n'))
    )
)

# This function allows you to just call "length(access)", rather than access$length
# Motivation for this is to make it easier to generically use accessors
length.modelVariableAccessorVector <- function(access)
	return(access$length)
