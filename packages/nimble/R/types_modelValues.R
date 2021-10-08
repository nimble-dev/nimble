## Options for requesting a modelValues
## 1. by a model
## 2. by a symbolTable (which can be created with buildSymbolTable)


#' Create a NIMBLE modelValues Object
#' 
#' Builds modelValues object from a model values configuration object, which can include a NIMBLE model
#' 
#' @param conf An object which includes information for building modelValues. Can either be a NIMBLE model (see \code{help(modelBaseClass)}) 
#' or the object returned from \code{modelValuesConf}
#' @param m	The number of rows to create in the modelValues object.  Can later be changed with \code{resize}
#' @author NIMBLE development team
#' @export
#' @details
#' See the \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} or \code{help(modelValuesBaseClass)} for information about manipulating NIMBLE modelValues object returned by this function
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
#' mvConf <- modelValuesConf(vars = c('x', 'y'),
#'              types = c('double', 'int'),
#'              sizes = list(x = 3, y = c(2,2)))
#' custom_mv <- modelValues(mvConf, m = 2)
#' custom_mv['y',]
modelValues <- function(conf, m = 1) {
    if(inherits(conf, 'RmodelBaseClass')) return(conf$modelDef$modelValuesClass(m))
    if(isModelValuesConf(conf)) return(conf(m))
    if(inherits(conf, 'symbolTable')) {
        mvClass <- modelValuesConf(conf) 
        return(mvClass(m))
    }
}	

#' Class \code{modelValuesBaseClass}
#' @export
#' @description
#'	modelValues are NIMBLE containers built to store values from models. They can either be built directly from 
#' a model or be custom built via the \code{modelValuesConf} function. They consist of rows, where each
#' row can be thought of as a set of values from a model. Like most nimble objects, and unlike most
#' R objects, they are passed by reference instead of by value. 
#'
#'See the \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} for more details.
#' @aliases [,CmodelValues-method [<-,CmodelValues-method [[,CmodelValues-method [[<-,CmodelValues-method [,CmodelValues-method,character,missing [,modelValuesBaseClass-method [<-,modelValuesBaseClass-method [,CmodelValues-method,character,missing,ANY-method [,CmodelValues-method,ANY,ANY [,CmodelValues-method,ANY,ANY
#' @examples
#'mvConf <- modelValuesConf(vars = c('a', 'b'), 
#'		types = c('double', 'double'), 
#'		sizes = list(a = 1, b = c(2,2) ) )
#'mv <- modelValues(mvConf)
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
                                        mvConf = 'ANY', 
                                        modelDef = 'ANY', 
                                        GID_map = 'ANY'),
                                    methods = list(
                                        initialize = function(nrow = 1L,...) {
                                            callSuper(...)
                                            nrow <<- nrow
                                    
                                            for(vN in varNames) {
                                                if(prod(sizes[[vN]]) == 1)
                                                    if(length(sizes[[vN]]) == 1) {
                                                        assign(vN, rep(list(as.numeric(NA)), nrow), inherits = TRUE)
                                                        next
                                                    }
                                                assign(vN, rep(list(array( data = as.numeric(NA), dim = sizes[[vN]])), nrow), inherits = TRUE)
                                            }
                                           ## GID_map <<- nimble:::makeMV_GID_Map(.self)
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
                                            stop('There was a call to expandNodeNames for a modelValues object.\n  This is deprecated and should not have occurred.\n  Please contact nimble developers at the nimble-users google group or at nimble.stats@gmail.com to let them know this happened.  \n Thank you.')
                                            ##return(GID_map$expandNodeNames(nodeNames = nodeNames, returnType = returnType, flatIndices = flatIndices))
                                        }                          
                                    )
                                    )


setMethod('[', 'modelValuesBaseClass',
	function(x, i, j, ..., drop){
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


#' Create the confs for a custom NIMBLE modelValues object
#' 
#' Builds an R-based modelValues conf object
#' 
#' @param vars	A vector of character strings naming each variable in the modelValues object
#' @param types	A vector of character strings describing the type of data for the modelValues object.
#' Options include `double' (for real-valued variables) and `int'.
#' @param sizes A list in which the named items of the list match the \code{var} arguments and each item is a numeric vector of the dimensions
#' @param symTab For internal use only
#' @param className For internal use only
#' @param where For internal use only
#' @param modelDef For internal use only
#' @author Clifford Anderson-Bergman
#' @export
#' @details
#' See the \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} or \code{help(modelValuesBaseClass)} and \code{help(modelValues)} for information
#'
#' @examples
#'	#Custom modelValues object:
#' mvConf <- modelValuesConf(vars = c('x', 'y'), 
#' 				types = c('double', 'int'), 
#'				sizes = list(x = 3, y = c(2,2)))
#' custom_mv <- modelValues(mvConf, m = 2)
#' custom_mv['y',]
modelValuesConf <- function( symTab, className, vars, types, sizes, modelDef = NA, where = globalenv() ) {
    if(missing(className)) className <- 'modelValuesConf' ## uniqueID will be appended
    makeCustomModelValuesClass(symTab, className, vars, types, sizes, modelDef = modelDef, where)
}

## Generates and evals a call to setRefClass with customized set of variables
## Input is EITHER symtab OR vars (which should be a character vector)
makeCustomModelValuesClass <- function(symTab, className, vars, types, sizes, modelDef, where = globalenv(), addUniqueID = TRUE){

    if(addUniqueID) className <- paste0(className, '_', nimbleUniqueID())
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
                        dims = nimbleInternalFunctions$dimOrLength(.self[[vN]][[1]])
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

    .isModelValuesConf <- TRUE ## will be looked in the environment of the returned function
    ans <- function(nrow = 1L) {
        res <- newClass(symTab, vars, sizes, nrow = nrow, modelDef = modelDef)
        res$mvConf <- ans
        res
    }
    ans
}

getModelValuesConf <- function(mv) {if(!inherits(mv, 'modelValuesBaseClass')) stop('Bad mv in getModelValuesConf', call. = FALSE); return(mv$mvConf)}
isModelValuesConf <- function(x) return(if(is.function(x)) exists('.isModelValuesConf', envir = environment(x), inherits = FALSE) else FALSE)

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
