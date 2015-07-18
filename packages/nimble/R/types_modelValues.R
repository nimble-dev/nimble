## Options for requesting a modelValues
## 1. by a model
## 2. by a symbolTable (which can be created with buildSymbolTable)


#' Create a NIMBLE modelValues Object
#' 
#' Builds modelValues object from a model values specification object, which can include a NIMBLE model
#' 
#' @param spec An object which includes information for building modelValues. Can either be a NIMBLE model (see \code{help(modelBaseClass)}) 
#' or the object returned from \code{modelValuesSpec}
#' @param m	The number of rows to create in the modelValues object.  Can later be changed with \code{resize}
#' @author NIMBLE development team
#' @export
#' @details
#' See the User Manual or \code{help(modelValuesBaseClass)} for information about manipulating NIMBLE modelValues object returned by this function
#'
#' @examples
#'	#From model object:
#' code <- nimbleCode({
#'  a ~ dnorm(0,1)
#'  for(i in 1:3){
#'		for(j in 1:3)
#'			b[i,j] ~ dnorm(0,1)
#'		}
#' })
#' Rmodel <- nimbleModel(code)
#' Rmodel_mv <- modelValues(Rmodel, m = 2)
#'	#Custom modelValues object:
#' mvSpec <- modelValuesSpec(vars = c('x', 'y'), types = c('double', 'int'), sizes = list(x = 3, y = c(2,2)))
#' custom_mv <- modelValues(mvSpec, m = 2)
#' custom_mv['y',]
modelValues <- function(spec, m = 1) {
    if(inherits(spec, 'RModelBaseClass')) return(spec$modelDef$modelValuesClass(m))
    if(isModelValuesSpec(spec)) return(spec(m))
    if(inherits(spec, 'symbolTable')) {
        mvClass <- modelValuesSpec(spec) 
        return(mvClass(m))
    }
}	

#' Class \code{modelValuesBaseClass}
#' @export
#' @description
#'	modelValues are NIMBLE containers built to store values from models. They can either be built directly from 
#' a model or be custom built via the \code{modelValuesSpec} function. They consist of rows, where each
#' row can be thought of as a set of values from a model. Like most nimble objects, and unlike most
#' R objects, they are passed by reference instead of by value. 
#'
#'See user manual for more details.
#' @examples
#'mvSpec <- modelValuesSpec(vars = c('a', 'b'), 
#'		types = c('double', 'double'), 
#'		sizes = list(a = 1, b = c(2,2) ) )
#'mv <- modelValues(mvSpec)
#'as.matrix(mv)
#'resize(mv, 2)
#'as.matrix(mv)
#'mv['a',1] <- 1
#'mv['a',2] <- 2
#'mv['b',1] <- matrix(0, nrow = 2, ncol = 2)
#'mv['b',2] <- matrix(1, nrow = 2, ncol = 2)
#'mv['a',]
#'as.matrix(mv)
#'basicModelCode <- nimbleCode({
#'	a ~ dnorm(0,1)
#'	for(i in 1:4)
#'		b[i] ~ dnorm(0,1)
#'})
#'basicModel <- nimbleModel(basicModelCode)
#'basicMV <- modelValues(basicModel, m = 2)	# m sets the number of rows
#'basicMV['b',]
modelValuesBaseClass <- setRefClass('modelValuesBaseClass',
                                    fields = list(
                                        varNames = 'character',
                                        name = 'ANY',
                                        symTab = 'ANY',
                                        sizes = 'ANY',
                                        nrow = 'numeric',  ## nrow is the actually the length of the lists for each variable
                                        CobjectInterface = 'ANY',
                                        mvSpec = 'ANY', 
                                        modelDef = 'ANY', 
                                        GID_map = 'ANY'),
                                    methods = list(
                                        initialize = function(nrow = 1L,...) {
                                            callSuper(...)
                                            nrow <<- nrow
                                    
                                            for(vN in varNames) {
                                                assign(vN, rep(list(array( data = as.numeric(NA), dim = sizes[[vN]])), nrow), inherits = TRUE)
                                            }
                                            GID_map <<- makeMV_GID_Map(.self)
                                        },
                                        getSymbolTable = function() {
                                            return(symTab)
                                        },
                                        getVarNames = function(includeLogProb = FALSE){
                                            if(includeLogProb)
                                                return(varNames)
                                            return(varNames[!grepl('logProb_', varNames)])
                                        },
                                        expandNodeNames = function(nodeNames, returnType = "names", flatIndices = TRUE) 
                                            {
                                                return(GID_map$expandNodeNames(nodeNames = nodeNames, returnType = returnType, flatIndices = flatIndices))
                                            }                          
                                    )
                                    )


setMethod('[', 'modelValuesBaseClass',
	function(x, i, j){
		if(missing(i) )
			i <- x$varNames
		if(is.numeric(i) ) 
			i <- x$varNames[i]
		if(missing(j) )
			j <- 1:x$nrow
			
		if(length(i) == 1 & length(j) == 1)
				return(x[[i]][[j]])
		if(length(i) == 1)
				return(x[[i]][j])
				
		new_nrow = length(j)
		newClassName = is(x)[1]
		newMV = new(newClassName, new_nrow)
		vN = x$varNames
		for(k in vN)
		if(!(k %in% i) ) 
			newMV[[k]] <- list()
		else
			newMV[[k]][1:new_nrow] = x[[k]] [ j[1:new_nrow] ] 
		return(newMV)
	})

setMethod('[<-', 'modelValuesBaseClass',
			function(x, i, j, ..., value) {
				if(missing(i) )
					i <- x$varNames
				if(missing(j) ) 
					j <- 1:x$nrow
				if(length(i) == 1 & length(j) == 1){
					x[[i]][[j]] <- value
					return(x)
					}
				if(length(i) == 1){
					x[[i]][j] <- value[1:length(j) ]
					return(x)
				}
				for(k in i)
					x[[k]][j] <- value[[k]][j]
				return(x)
			})


#' Create the specs for a custom NIMBLE modelValues Object
#' 
#' Builds an R-based modelValues spec object
#' 
#' @param vars	A vector of character strings naming each variable in the modelValues object
#' @param types	A vector of character strings describing the type of data for the modelValues object.
#' Options include `double' (for real-valued variables) and `int'.
#' @param sizes A list in which the named items of the list match the \code{var} arguments and each item is a numeric vector of the dimensions
#' @param symTab For internal use only
#' @param className For internal use only
#' @param where For internal use only
#' @author Clifford Anderson-Bergman
#' @export
#' @details
#' See the User Manual or \code{help(modelValuesBaseClass)} and \code{help(modelValues)} for information
#'
#' @examples
#'	#Custom modelValues object:
#' mvSpec <- modelValuesSpec(vars = c('x', 'y'), 
#' 				types = c('double', 'int'), 
#'				sizes = list(x = 3, y = c(2,2)))
#' custom_mv <- modelValues(mvSpec, m = 2)
#' custom_mv['y',]
#' [[1]]
#'      [,1] [,2]
#' [1,]   NA   NA
#' [2,]   NA   NA
#'
#' [[2]]
#'      [,1] [,2]
#' [1,]   NA   NA
#' [2,]   NA   NA
modelValuesSpec <- function( symTab, className, vars, types, sizes, modelDef = NA, where = globalenv() ) {
    if(missing(className)) className <- 'modelValuesSpec' ## uniqueID will be appended
    makeCustomModelValuesClass(symTab, className, vars, types, sizes, modelDef = modelDef, where)
}

## Generates and evals a call to setRefClass with customized set of variables
## Input is EITHER symtab OR vars (which should be a character vector)
makeCustomModelValuesClass <- function(symTab, className, vars, types, sizes, modelDef, where = globalenv(), addUniqueID = TRUE){

    if(addUniqueID) className <- paste0(className, '_', nimbleUniqueID())
    ## if(exists(className, modelValuesLibrary, inherits = FALSE)) {
    ## 	oldSymTab = modelValuesLibrary[[className]]$symTab
    ## 	if(areMVSymTabsEqual(oldSymTab, SymbolTable) )
    ##             return(modelValuesLibrary[[className]]$modelValuesClass)
	        
    ##        className <- modelValuesClassLabelCreator()   
    ## }
    if(!missing(symTab) ) {    
        vars <- symTab$getSymbolNames()
        sizes <- list()
        for(vN in vars) {
            sizes[[vN]] <- symTab$symbols[[vN]]$size
            if(length(sizes[[vN]]) == 0) sizes[[vN]] <- 1
        }
    } else {
        if(!is.list(sizes) ) {
            sizes <- list(sizes)
            names(sizes) <- vars
        }
        if(length(types) != length(vars))
        	stop('Creation of modelValues aborted: length(types) != length(vars)')
        symTab <- buildSymbolTable(vars, types, sizes)
    }
    if(missing(modelDef))
    	modelDef <- NA
    ans <- eval(substitute(newClass <- setRefClass(
        Class = className,
        contains = 'modelValuesBaseClass',
        fields = FIELDS,
        methods = list(
            show = function() {
                writeLines(paste0("modelValues object with variables: ", paste(varNames, collapse = ", "), "."))
            },
            initialize = function(symTab, vars, sizes, modelDef, ...) {
                symTab <<- symTab
                varNames <<- vars
                sizes <<- sizes
                modelDef <<- modelDef
                callSuper(...)
            },
            resize = function(rows) {
                if(nrow > rows){
                    for(vN in varNames)
                        .self[[vN]] <- .self[[vN]][1:rows]
                }	
                else if(nrow < rows){
                    for(vN in varNames){
                        dims = dimOrLength(.self[[vN]][[1]])	#nimDim(.self[[vN]][[1]])
                        if(length(dims) == 0) dims = 1
                        for(i in (nrow+1):rows){
                            .self[[vN]][[i]] <- array(NA, dim = dims)
                            if( is.integer(.self[[vN]][[1]]) )
                            	.self[[vN]][[i]] <- as.integer(.self[[vN]][[i]])
                            }
                    }	
                }
                nrow <<- rows
            },
            getSize = function(rows){
                return(nrow)
            }            
            
            ), where = where),
                           list(FIELDS = makeModelValuesClassFields(vars),
                                className = className
                                )))

    .isModelValuesSpec <- TRUE ## will be looked in the environment of the returned function
    ## cppClassName <- Rname2CppName(className)
    ## modelValuesLibrary[[className]] <- list(modelValuesClass = ans,
    ##                                         symTab = SymbolTable,
    ##                                         cppClassName = cppClassName,
    ##                                         cppClass = NULL)
    ans <- function(nrow = 1L) {
        res <- newClass(symTab, vars, sizes, nrow = nrow, modelDef = modelDef)
        res$mvSpec <- ans
        res
    }
    ans
}

getModelValuesSpec <- function(mv) {if(!inherits(mv, 'modelValuesBaseClass')) stop('Bad mv in getModelValuesSpec', call. = FALSE); return(mv$mvSpec)}
isModelValuesSpec <- function(x) return(if(is.function(x)) exists('.isModelValuesSpec', envir = environment(x), inherits = FALSE) else FALSE)

## Generates unevaluated code for field definitions for a derived modelValues class
makeModelValuesClassFields <- function(vars) {
    varDefs <- as.list(rep('list', length(vars)  ) )
    names(varDefs) <- c(vars)
    as.call(c(as.name("list"), varDefs))
}
## e.g. makeModelValuesClassFields(c('a','b','c'))



pointAt <- function(model, to, vars = NULL, toVars = NULL, index = NA,  logProb = FALSE) {
    if(is.null(vars)) {
        vars <- model$getVarNames(includeLogProb = TRUE)
    } else {
        if(!all(vars %in% model$getVarNames(includeLogProb = TRUE)))
            stop(paste('Error, trying to do pointAt incorrectly for variable(s):', paste(vars[ which(!(vars %in% model$getVarNames(includeLogProb = TRUE)))], collapse = ', ')))
    }
    
    ## This won't work generally because it needs to determine which vars actually have logProbs
    if(logProb)         { vars <- c(vars, makeLogProbName(vars))  }
    
    if(is.null(toVars)) toVars <- vars
  
  for(i in seq_along(vars)) {
      model[[makeEnvName(vars[i])]] <- to
      model[[makeNameName(vars[i])]] <- toVars[i]
      model[[makeRowName(vars[i])]] <- as.integer(index)
  }
}



makeMV_GID_Map <- function(mv){
	sizeList = mv$sizes
    varNames = sort(mv$varNames)
    nodeNames = NA
    nodeIndex = 0
    nodeNames2GID_maps <- new.env()
    all.names <- character()
    all.flatIndexNames <- character()
    for (i in seq_along(varNames)) {
        baseName = varNames[i]
        dims = sizeList[[baseName]]
        if (length(dims) == 1 & dims[1] == 1) {
            nodeNames[nodeIndex + 1] = baseName
            nodeIndex <- nodeIndex + 1
            nodeNames2GID_maps[[baseName]] <- nodeIndex
            all.names <- c(all.names, baseName)
            all.flatIndexNames <- c(all.flatIndexNames, baseName)
        }
        else {
            mins <- rep(1, length(dims))
            indexStuff <- paste(mins, dims, sep = ":", collapse = ", ")
            compactNodeNames <- paste0(baseName, "[", indexStuff, 
                "]")
            expandedNodeNames <- nl_expandNodeIndex(compactNodeNames)
            nodeNames[nodeIndex + 1:length(expandedNodeNames)] <- expandedNodeNames
            
            nodeNames2GID_maps[[baseName]] <- array(dim = dims)
            nodeNames2GID_maps[[baseName]][1:length(expandedNodeNames)] <- nodeIndex + 1:length(expandedNodeNames)
            nodeIndex = nodeIndex + length(expandedNodeNames)
        	all.names <- c(all.names, expandedNodeNames)
        	all.flatIndexNames <- c(all.flatIndexNames, paste0(baseName, '[', 1:length(expandedNodeNames) ,']' ))
        }
    	
    }
    GID_Map <- new.env()
    GID_Map[['nodeNames2GID_maps']] <- nodeNames2GID_maps
    GID_Map[['graphID_2_nodeName']] <- all.names
    GID_Map[['graphID_2_flatIndexNames']] <- all.flatIndexNames
    GID_Map[['expandNodeNames']] <- function(nodeNames, returnType = 'names', flatIndices = TRUE){
    	if(length(nodeNames) == 0)	return(NULL)
    	gIDs <- unlist(sapply(nodeNames, parseEvalNumeric, env = GID_Map$nodeNames2GID_maps, USE.NAMES = FALSE) )
    	gIDs <- gIDs[!is.na(gIDs)]
    	if(returnType == 'names'){
    		if(flatIndices == TRUE)
    			return(GID_Map$graphID_2_flatIndexNames[gIDs])
    		return(GID_Map$graphID_2_nodeName[gIDs])
    		}
    	return(gIDs)
    }
    return(GID_Map)
}
