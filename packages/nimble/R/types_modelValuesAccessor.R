
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

valuesAccessorVector <- setRefClass( ## new implementation
    Class = 'valuesAccessorVector',
    fields = list(
        sourceObject = 'ANY',
        nodeNames = 'ANY',
        logProb = 'ANY', 
        length = 'ANY',
        code = 'ANY',
        accessCode = 'ANY',
        setCode = 'ANY',
        mapInfo = 'ANY'
    ),
    methods = list(
        initialize = function(modelOrModelValues, nodeNames, logProb = FALSE) { ## will need to work with IDs too
            sourceObject <<- modelOrModelValues
            inputNodeNames <- nodeNames
             if(!logProb)
                nodeNames <<- inputNodeNames
            else {
                isLogProbName <- grepl('logProb_', inputNodeNames)
                nodeNames <<- c(inputNodeNames, makeLogProbName(inputNodeNames[!isLogProbName]))
            }
            logProb <<- logProb
            code <<- lapply(.self$nodeNames, function(x) parse(text = x, keep.source = FALSE)[[1]])
            length <<- NA
            mapInfo <<- NULL
        },
        makeMapInfo = function() {
            length <<- 0
            mapInfo <<- lapply(code, function(x) {
                varAndIndices <- getVarAndIndices(x)
                varName <- as.character(varAndIndices$varName)
                varSym <- sourceObject$getSymbolTable()$getSymbolObject(varName) ## previously from model$getVarInfo(varName)
                ans <- varAndIndices2mapParts(varAndIndices, varSym$size, varSym$nDim)
                ans$varName <- varName
                anslength <- prod(ans$sizes)
                length <<- length + anslength
                ans$singleton <- anslength == 1
                ans
            }) ## list elements will be offset, sizes, strides, varName, and singleton in that order.  Any changes must be propogated to C++
            invisible(NULL)
        },
        getMapInfo = function() {
            if(is.null(mapInfo)) makeMapInfo()
            mapInfo
        },
        getLength = function() {
            if(is.null(mapInfo)) makeMapInfo()
            writeLines(paste('length = ', length))
            length
        }
    )
    )


makeSingleModelValuesAccessor <- function(ids, modelValues){
	return(modelValuesAccessor(id = id, modelValues = modelValues))
    }

modelValuesAccessorVector <- setRefClass( ## new implementation
   Class = 'modelValuesAccessorVector',
    contains = 'valuesAccessorVector',
    fields = list(row = 'ANY'),
    methods = list(
        initialize = function(...) {
            callSuper(...)
            row <<- as.integer(NA)
            makeAccessAndSetCode()
        },
        makeAccessAndSetCode = function() {
            accessCode <<- lapply(code, function(temp) {
                if(is.name(temp)) return(substitute(sourceObject$B[[localrow]], list(B = temp)))
                temp[[2]] <- substitute(sourceObject$B[[localrow]], list(B = temp[[2]]))
                temp
            })
            setCode <<- lapply(accessCode, function(x) substitute(A <- vals, list(A = x)))
        },
        setRow = function(row) {
            row <<- row
        },
        getValues = function(i, row = NA) {
            ## model version
            localrow <- if(is.na(row)) .self$row else row
            eval(accessCode[[i]])
        },
        setValues = function(i, vals, row = NA) {
            ## model version
            localrow <- if(is.na(row)) .self$row else row
            eval(setCode[[i]])
        }
    ))

modelValuesAccessorVectorOld <- setRefClass( ## old implementation
    Class = 'modelValuesAccessorVectorOld',
    fields = list(modelValues = 'ANY',
                  gids = 'ANY', 		#'character',
                  length = 'ANY'		#'numeric'),
                  ),
    methods = list(
        initialize = function(modelValues, nodeNames, logProb = FALSE, env = parent.frame()) {
	        isLogProbName <- grepl('logProb_', nodeNames)
	        nodeOnlyNames <- nodeNames[!isLogProbName]
	        suppliedLogProbNames <- nodeNames[isLogProbName]
        	inferredLogProbNames <- character(0)
        	if(logProb){
        		if(!inherits(modelValues$modelDef,'modelDefClass'))
	        		stop('calling logProb = TRUE on a modelValues object that was not built from a model')
        		inferredLogProbNames <- modelValues$modelDef$nodeName2LogProbName(nodeOnlyNames)
        	}
        	nodeNames <- c(nodeOnlyNames, suppliedLogProbNames, inferredLogProbNames)
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
