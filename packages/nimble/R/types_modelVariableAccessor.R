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
    			  l_gids = 	'ANY',	#graph IDS for log probabilities, which is a different set of graph IDs than for nodes
                  length = 	'ANY'
                   ),		
    methods = list(
        initialize = function(model, nodeNames, logProb = FALSE, env = parent.frame()) {
        	isLogProbName <- grepl('logProb_', nodeNames)
        	nodeOnlyNames <- nodeNames[!isLogProbName]
        	logProbNames <- nodeNames[isLogProbName]
        	gids <<- model$expandNodeNames(nodeOnlyNames, returnScalarComponents = TRUE, returnType = 'ids')
        	l_gids <<- numeric(0)
			if(length(logProbNames) > 0){
				nodeName_fromLogProbName <- gsub('logProb_', '', logProbNames)
				l_gids <<- c(l_gids, model$modelDef$nodeName2LogProbID(nodeName_fromLogProbName))
			}
        	if(logProb)
        		l_gids <<- c(l_gids, model$modelDef$nodeName2LogProbID(nodeOnlyNames))
            length <<- length(gids) + length(l_gids)
   			model <<- model
         },
        getSingleValue_fromGID = function(accessID){
        	len_variables <- length(gids)
        	if(accessID <= len_variables){
	        	thisExpr <- parse(text = model$expandNodeNames(gids[accessID]))[[1]]
	        	return( eval(thisExpr, envir = model) )
	        	}
	        else{
	        	logVarID <- l_gids[accessID - len_variables]
	        	logProbName <- model$modelDef$maps$logProbIDs_2_LogProbName[logVarID]
	        	thisExpr <- parse(text = logProbName)[[1]]
	        	return(eval(thisExpr, envir = model))
	        	}
	        	
        },
        setSingleValue_fromGID = function(value, accessID){
        	len_variables <- length(gids)
        	if(accessID <= len_variables){
	        	thisExpr <- parse(text = paste0(model$expandNodeNames(gids[accessID]), '<-', value))		
	        	eval(thisExpr, envir = model)
	        }
	        else{
	        	logVarID <- l_gids[accessID - len_variables]
	        	logProbName <- model$modelDef$maps$logProbIDs_2_LogProbName[logVarID]
	        	thisExpr <- parse(text = paste0(logProbName, '<-', value))[[1]]
				eval(thisExpr, envir = model)
	        }
        },
        getNodeNames = function(){
        	c(model$expandNodeNames(gids, returnScalarComponents = TRUE), 	model$modelDef$maps$logProbIDs_2_LogProbName[l_gids])
        },
        show = function(){
	        cat(paste0('modelVariableAccessorVector: ', paste0(lapply(modelVariableAccessors, function(x) x$toStr()), collapse=', '), '\n'))
	        }
    )
)

# This function allows you to just call "length(access)", rather than access$length
# Motivation for this is to make it easier to generically use accessors
length.modelVariableAccessorVector <- function(access)
	return(access$length)
