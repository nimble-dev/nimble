###		These functions are used for calculate/sim/getLP for the nodeFunctionVectors
###		Can either enter model, nodes or model_nodes

#' Turn a numeric vector into a single-row or single-column matrix
#'
#' Turns a numeric vector into a matrix that has 1 row or 1 column.  Part of NIMBLE language executable in R.
#'
#' @aliases asCol
#'
#' @param x       Numeric to be turned into a single row or column matrix
#'
#' @author Perry de Valpine
#' @export
#' @details
#' In the NIMBLE language, some automatic determination of how to turn vectors into single-row or single-column matrices is done.
#' For example, in \code{A \%*\% x}, where A is a matrix and x a vector, x will be turned into a single-column matrix unless
#' it is known at compile time that A is a single column, in which case x will be turned into a single-row matrix.
#' However, if it is desired that x be turned into a single row but A cannot be determined at compile time to be a single column,
#' then one can use \code{A \%*\% asRow(x)} to force this conversion.
asRow <- function(x) {
    matrix(x, nrow = 1)
}

## Aliased in asRow
asCol <- function(x) {
    matrix(x, ncol = 1)
}

rCalcNodes <- function(model, nodes){
	l_Prob = 0
	for(nName in nodes)
		l_Prob = l_Prob + model$nodes[[nName]]$calculate()
	return(l_Prob)
}

#' calculate, simulate, or get the current log probabilities (densities) a set of nodes in a NIMBLE model
#'
#' calculate, simulate, or get the current log probabilities (densities) of one or more nodes of a NIMBLE model and (for calculate and getLogProb) return the sum of their log probabilities (or densities).  Part of R and NIMBLE.
#'
#' @aliases simulate getLogProb
#' 
#' @param model        A NIMBLE model, either the compiled or uncompiled version
#' @param nodes        A character vector of node names, with index blocks allowed, such as 'x', 'y[2]', or 'z[1:3, 2:4]'
#' @author NIMBLE development team
#' @export
#' @details
#' These functions expands the nodes and then process them in the model in the order provided.  Expanding nodes means turning 'y[1:2]' into c('y[1]','y[2]') if y is a vector of scalar nodes.
#' Calculation is defined for a stochastic node as executing the log probability (density) calculation and for a deterministic node as calculating whatever function was provided on the right-hand side of the model declaration.
#'
#' Simulation is defined for a stochastic node as drawing a random value from its distribution, and for deterministic node as equivalent to calculate.
#'
#' getLogProb simply collects the sum of the log probabilities of nodes if they are known to have already been calculated.
#'
#' These functions can be used from R or in NIMBLE run-time functions that will be compiled.  When executed in R (including when an uncompiled nimbleFunction is executed), they can be slow because the nodes are expanded each time.  When compiled in NIMBLE, the nodes are expanded only once during compilation, so execution will be much faster.
#'
#' It is common to want the nodes to be provided in topologically sorted order, so that they will be calculated or simulated following the order of the model graph.  Functions such as model$getDependencies(nodes, ...) return nodes in topologically sorted order.  They can be directly sorted by model$topologicallySortNodes(nodes), but if so it is a good idea to expand names first by model$topologicallySortNodes(model$expandNodeNames(nodes))
#'
#' @return calculate and getLogProb return the sum of the log probabilities (densities) of the calculated nodes, with a contribution of 0 from any deterministic nodes  simulate returns NULL.
#' 
#' @examples
#' calculate(model, c('x', 'y[2:4]', 'z[2:5, 1:10]'))
#' 
calculate <- function(model, nodes, nodeFxnVector)		
{
	if(!missing(nodeFxnVector)){
			model <- nodeFxnVector$model
			nodes <- nodeFxnVector$getNodeNames()
		return(rCalcNodes(model, nodes))
		}
	if(inherits(model, 'modelBaseClass') ){
		if(missing(nodes) ) 
			nodes <- model$getMaps('nodeNamesLHSall')
		nfv <- nodeFunctionVector(model, nodes)
		nodeNames <- nfv$getNodeNames()
		return(rCalcNodes(model, nodeNames))
	}	
}

rGetLogProbsNodes <- function(model, nodes){
	l_Prob = 0
	for(nName in nodes)
		l_Prob = l_Prob + model$nodes[[nName]]$getLogProb()
	return(l_Prob)
}

getLogProb <- function(model, nodes, nodeFxnVector)		
{
	if(!missing(nodeFxnVector)){
		model <- nodeFxnVector$model
		nodes <- nodeFxnVector$getNodeNames()
		return(rGetLogProbsNodes(model, nodes))
	}
	if( inherits(model, "modelBaseClass") ){		
		if(missing(nodes) ) 
                    nodes <- model$getMaps('nodeNamesLHSall')

		nfv <- nodeFunctionVector(model, nodes)
		nodeNames <- nfv$getNodeNames()

    	return(rGetLogProbsNodes(model, nodeNames))
    }        
}


rSimNodes <- function(model, nodes){
	for(nName in nodes)
		model$nodes[[nName]]$simulate()
}

simulate <- function(model, nodes, includeData = FALSE, nodeFxnVector)		
{
	if(!missing(nodeFxnVector)){
		model <- nodeFxnVector$model
		nodes <- nodeFxnVector$getNodeNames()
		rSimNodes(model, nodes)
	}
	else if( inherits(model, "modelBaseClass") ) {
		if(missing(nodes) ) 
			nodes <- model$getMaps('nodeNamesLHSall')
		nfv <- nodeFunctionVector(model, nodes, excludeData = !includeData)
		nodeNames <- model$expandNodeNames(nfv$gids)			
		rSimNodes(nfv$model, nodeNames)
	}
}


# Fill a vector with flattened values from a set of model nodes
#
# Take a vector of node names for a model and fill a vector by concatenating their values.  Works in R and NIMBLE.
#
# @param vals        the variable in the calling function to receive the values
# @param model       a NIMBLE model object, either compiled or uncompiled
# @param nodes       a vector of node names, allowing index blocks that will be expanded
#
# @author NIMBLE development team
# @export
# @details
# Calling \code{getValues(P, model, nodes)} will modify P in the calling function.  This is being deprecated by
# \code{P <- values(model, nodes)}, but the development syntax of \code{getValues(P, model, nodes)} is still supported.
#
# The result P will be the concatenation of values from the nodes requested.  When requested nodes are from matrices or arrays, the values will be flattened into a vector following column-wise order.
#
# The reverse of \code{getValues(P, model, nodes)} is \code{setValues(P, model, nodes)}, which is being deprecated by \code{values(model, nodes) <- P}.
#
# These functions work in R and in NIMBLE run-time code that can be compiled.
#
# @return NULL, but this function works by the side-effect of modifying P in the calling environment.


getValues <- function(vals, model, nodes)
	{
	valsExp = substitute(vals)
	access = modelVariableAccessorVector(model, nodes, logProb = FALSE)
	output = getValuesAccess(vals, access)
	assign(as.character(valsExp), output, envir = parent.frame(n = 1) )
	}

getValuesAccess <- function(vals, access){
#	if(length(vals)!= length(access) ) {
#		writeLines('Length of object to copy into does not match')
#	}
	output = as.numeric(NA)
	totValues <- length(access$gids) + length(access$l_gids)
	for(i in seq_along(access$gids) )
		output[i] <- access$getSingleValue_fromGID(i)
	gid_len = length(access$gids)
	for(i in seq_along(access$l_gids))
		output[i + gid_len] = access$getSingleValue_fromGID(i + gid_len)
	return(output)
}



setValuesAccess <- function(input, access){
	tot_length = length(access$gids) + length(access$l_gids)
	if(length(input)!= tot_length ) 
		writeLines('Length of input does not match accessor')
	else{
		for(i in seq_along(access$gids) )
			access$setSingleValue_fromGID(input[i], i)
		gid_len = length(access$gids)
		for(i in seq_along(access$l_gids) )
			access$setSingleValue_fromGID(input[i+gid_len], i + gid_len)
	}
}	

# Fill a set of nodes in a model from a vector of values
#
# Take a vector of node names for a model and fill them sequentially from the values in a vector.  Works in R and NIMBLE.
#
# @param input        a vector of values to be put in the model's nodes
# @param model       a NIMBLE model object, either compiled or uncompiled
# @param nodes       a vector of node names, allowing index blocks that will be expanded
#
# @author NIMBLE development team
# @export
# @details
# Calling setValues(P, model, nodes) will place values from P, in order, into the nodes provided for the model.  This is being deprecated by
# values(model, nodes) <- P, but the development syntax of setValues(P, model, nodes) is still supported.
#
# When provided nodes are from matrices or arrays, the values will filled following column-wise order.
#
# The reverse of setValues(P, model, nodes) is getValues(P, model, nodes), which is being deprecated by P <- values(model, nodes)
#
# These functions work in R and in NIMBLE run-time code that can be compiled.
#
# @return NULL, but this function works by the side-effect of modifying the model.


setValues <- function(input, model, nodes){
	access = modelVariableAccessorVector(model, nodes, logProb = FALSE)
	setValuesAccess(input, access)
	
}

#' Access values for a set of nodes in a model
#'
#' @aliases values values<-
#' 
#' @description
#' Get or set values for a set of nodes in a model
#'
#' @param model       a NIMBLE model object, either compiled or uncompiled
#' @param nodes       a vector of node names, allowing index blocks that will be expanded
#'
#' @author NIMBLE development team
#' @export
#' @details
#' Calling \code{values(model, nodes)} returns a vector of the concatenation of values from the nodes requested
#' \code{P <- values(model, nodes)} is a newer syntax for \code{getValues(P, model, values)}, which still works and modifies P in the calling environment.
#'
#' Calling \code{values(model, nodes) <- P} sets the value of the nodes in the model, in sequential order from the vector P.
#'
#' In both uses, when requested nodes are from matrices or arrays, the values will be handled following column-wise order.
#'
#' The older function \code{getValues(P, model, nodes)} is equivalent to \code{P <- values(model, nodes)}, and the older function \code{setValues(P, model, nodes)} is equivalent to \code{values(model, nodes) <- P}
#'
#'
#' These functions work in R and in NIMBLE run-time code that can be compiled.
#'
#' @return A vector of values concatenated from the provided nodes in the model
values <- function(model, nodes){
	ans <- NA
	getValues(ans, model, nodes)
	ans
}

`values<-` <- function(model, nodes, value){
	setValues(value, model, nodes)
	return(model)
}

#' Copying function for NIMBLE
#' 
#' Copies values from a NIMBLE model or modelValues object to another NIMBLE model or modelValues. Work in R and NIMBLE.  The NIMBLE keyword \code{copy} is identical to \code{nimCopy}
#' 
#' @param from		Either a NIMBLE model or modelValues object
#' @param to		Either a NIMBLE model or modelValues object
#' @param nodes		The nodes of object \code{from} which will be copied from
#' @param nodesTo	The nodes of object \code{to} which will be copied to. If \code{nodesTo == NA}, will automatically be set to \code{nodes}
#' @param row		If \code{from} is a modelValues, the row which will be copied from
#' @param rowTo		If \code{to} is a modelValues, the row which will be copied to. If \code{rowTo == NA}, will automatically be set to \code{row}
#' @param logProb	A logical value indicating whether the log probabilities of the given nodes should also be copied (i.e. if \code{nodes = 'x'}
#' and \code{logProb = TRUE}, then both \code{'x'} and \code{'logProb_x'} will be copied)
#' 
#' @author Clifford Anderson-Bergman
#' @export
#' @details
#'
#' See the User Manual for more details
#'
#' @examples
#'	# Building model and modelValues object
#' simpleModelCode <- nimbleCode({
#'	for(i in 1:100)
#'		x[i] ~ dnorm(0,1)
#'})
#' rModel <- nimbleModel(simpleModelCode)
#' rModelValues <- modelValues(rModel)
#'
#' #Setting model nodes
#' rModel$x <- rnorm(100)
#' #Using nimCopy in R.
#' nimCopy(from = rModel, to = rModelValues, nodes = 'x', rowTo = 1)
#' 
#' #Use of nimCopy in a simple nimbleFunction
#' cCopyGen <- nimbleFunction(
#' 	setup = function(model, modelValues, nodeNames){},
#' 	run = function(){
#' 		nimCopy(from = model, to = modelValues, nodes = nodeNames, rowTo = 1)
#' 	}
#' )
#' 
#' rCopy <- cCopyGen(rModel, rModelValues, 'x')
#' cModel <- compileNimble(rModel)
#' cCopy <- compileNimble(rCopy, project = rModel)
#' cModel[['x']] <- rnorm(100)
#' 
#' cCopy() ## execute the copy with the compiled function
nimCopy <- function(from, to, nodes, nodesTo = NA, row = NA, rowTo = NA, logProb = FALSE){

    isFromModel = NA
    isToModel = NA
    if(missing(nodes) ) 
        nodes = allNodeNames(from)
    if( inherits(from, "modelBaseClass") ){
        accessFrom = modelVariableAccessorVector(from, nodes, logProb = logProb)
        rowFrom = NA
        isFromModel = TRUE
    }
    else if(inherits(from, "modelValuesBaseClass"))
        {
            accessFrom = modelValuesAccessorVector(from, nodes, logProb = logProb)
            if(is.na(row))
                stop("Error: need to supply 'row' for a modelValues copy")
            rowFrom = row
            isFromModel = FALSE
        }
    if( inherits(to, "modelBaseClass") ){
        if(is.na(nodesTo[[1]]) ) 
            accessTo = modelVariableAccessorVector(to, nodes, logProb = logProb)
        else
            accessTo = modelVariableAccessorVector(to, nodesTo, logProb = logProb)
        rowTo = NA
        isToModel = TRUE
    }
    else if(inherits(to, "modelValuesBaseClass"))
        {
            if(is.na(nodesTo[[1]]) ) 
                accessTo = modelValuesAccessorVector(to, nodes, logProb = logProb)
            else
                accessTo = modelValuesAccessorVector(to, nodesTo, logProb = logProb)
            if(is.na(rowTo))
                rowTo = row
            isToModel = FALSE
        }
    if(is.na(isFromModel))
        stop('argument "from" in nimCopy is neither a model nor modelValues')
    if(is.na(isToModel))
        stop('argument "to" in nimCopy is neither a model nor modelValues')

    
    lengthTo = accessTo$length
    lengthFrom = accessFrom$length
    if(lengthTo != lengthFrom)
        stop('lengths not equal in nimCopy') 
    if(lengthTo > 0){
        for(i in 1:lengthTo){
            if(isFromModel)
                valueFrom = accessFrom$getSingleValue_fromGID(i)
            else
                valueFrom = accessFrom$getSingleValue_fromGID(i, rowFrom)
            
            if(isToModel)
                accessTo$setSingleValue_fromGID(valueFrom, i)
            else
                accessTo$setSingleValue_fromGID(valueFrom, i, rowTo)
        }
    }
}


allNodeNames <- function(object, logProb = FALSE){
	if(inherits(object, 'modelValuesBaseClass') ) {	
		if(!is.null(object$modelDef))
			all.Names <- ls(object$modelDef$nodeInfo)
		else
			all.Names = object$varNames
		if(logProb == TRUE)
				return(all.Names)
		for(i in 1:length(all.Names) ) {
			if(gsub("logProb_", "", all.Names[i]) != all.Names[i])
				all.Names[i] = NA		
		}
		all.Names = all.Names[!is.na(all.Names)]
		return(all.Names)	
		}
	if(inherits(object, 'modelBaseClass') ) {
		all.Names <- object$getNodeNames()	
		if(logProb == TRUE)
				return(all.Names)
		for(i in 1:length(all.Names) ) {
			if(gsub("logProb_", "", all.Names[i]) != all.Names[i])
				all.Names[i] = NA		
		}
		all.Names = all.Names[!is.na(all.Names)]
		return(all.Names)	
	}
}

#' Access or set a member variable of a nimbleFunction created during \code{setup}
#'
#' @aliases nfVar nfVar<-
#'
#' @description
#' Access or set a member variable of a specialized nimbleFunction, i.e. a variable passed to or created during the \code{setup} function that is used in run code or preserved by \code{setupOutputs}.  Works in R for any variable and in NIMBLE for numeric variables.
#'
#' @param nf      A specialized nimbleFunction, i.e. a function returned by executing a function returned from \code{nimbleFunction} with \code{setup} arguments
#' @param varName A character string naming a variable in the \code{setup} function.
#' @author NIMBLE development team
#' @export
#' @details
#' When \code{nimbleFunction} is called and a \code{setup} function is provided, then \code{nimbleFunction} returns a function.  That function is a generator that should be called with arguments to the \code{setup} function and returns another function that will execute the \code{run} function.  The \code{run} function and any other \code{methods} provided can use objects created or passed to \code{setup}.  \code{nfVar} provides direct access to those objects by name.
#'
#' In R, \code{nfVar} can retrieve anything in \code{setup}, including models, other nimbleFunctions, and vectors of node lists.  In NIMBLE, it can only retrieve numeric variables
#'
#' For access to methods of \code{nf}, see \code{nfMethod}.
#' 
#' For more information, see \code{?nimbleFunction} and the NIMBLE User Manual.
#'
#' @return whatever varName is in  the nimbleFunction nf.
#' @examples
#' nfGen1 <- nimbleFunction(
#'     setup = function(A) {
#'     B <- matrix(rnorm(4), nrow = 2)
#'     setupOutputs(B) ## preserves B even though it is not used in run-code
#'    },
#'    run = function() {
#'       print('This is A', A, '\n')
#' })
#'
#' nfGen2 <- nimbleFunction(
#'   setup = function() {
#'     nf1 <- nfGen1(1000)
#'   },
#'   run = function() {
#'       print('accessing A:', nfVar(nf1, 'A'))
#'       nfVar(nf1, 'B')[2,2] <<- -1000
#'       print('accessing B:', nfVar(nf1, 'B'))
#'    })
#'        
#' nf2 <- nfGen2()
#' nf2()
#' Cnf2 <- compileNimble(nf2)
#' Cnf2()
nfVar <- function(nf, varName) {
    refClassObj <- nf_getRefClassObject(nf)
    v <- refClassObj[[varName]]
    return(v)
}

`nfVar<-` <- function(nf, varName, value) {
    refClassObj <- nf_getRefClassObject(nf)
    refClassObj[[varName]] <- value
    return(nf)
}

#' access a member function (or method) of a nimbleFunction
#'
#' access (call) a member function of a nimbleFunction other than \code{run}, which calling the nimbleFunction directly uses.  Works in R and NIMBLE.
#'
#' @param nf          a specialized nimbleFunction, i.e. one that has already had setup parameters processed
#' @param methodName  a character string giving the name of the member function to call
#'
#' @author NIMBLE development team
#' @export
#' @details
#' nimbleFunctions have a default member function called \code{run}, and may have other member functions provided via the \code{methods} argument to \code{nimbleFunction}.
#' \code{nfMethod} is primarily used to call the latter, although it could be used with \code{methodName = 'run'}.  Normally arguments will be provided after \code{nfMethod}, e.g. \code{nfMethod(myNimbleFunction, 'foo')(x)}.
#'
#' @return a function that can be called.
nfMethod <- function(nf, methodName) {
    # refClassObj <- nf_getRefClassObject(nf)
    # methodCall <- substitute(refClassObj$METHODNAME, list(METHODNAME = methodName))
    # eval(methodCall)
    eval(substitute(nf_getRefClassObject(nf)$METHODNAME, list(METHODNAME = methodName)))
}

#' Generates a weighted sample (with replacement) of ranks
#'
#' Takes a set of non-negative \code{weights} (do not need to sum to 1) and 
#' returns a sample with \code{size} elements of the integers \code{1:length(weights)}, where the probability of being sampled is proportional
#' to the value of \code{weights}. An important note is that the output vector
#' will be sorted in ascending order. Also, right now it works slightly odd syntax (see example below). Later releases of NIMBLE will contain more natural syntax. 
#'
#' If invalid weights provided (i.e. negative weights or weights sum to 1), sets output = rep(1, size) and prints warning. 
#' \code{rankSample} can be used inside nimble functions.
#'
#' @param weights		A vector of numeric weights. Does not need to sum to 1, but must be non-negative
#' @param size			size of sample
#' @param output		an R object into which the values will be placed. See example below for proper use
#' @author	Clifford Anderson-Bergman
#' @export
#' @details		

#' \code{rankSample} first samples from the joint distribution \code{size} uniform(0,1) distributions by conditionally sampling from the rank statistics. This leads to 
#' a sorted sample of uniform(0,1)'s. Then, a cdf vector is constructed from weights. Because the sample of uniforms is sorted, \code{rankSample} walks
#' down the cdf in linear time and fills out the sample.
#'  
#' @examples
#' set.seed(1)
#' sampInts = NA	#sampled integers will be placed in sampInts
#' rankSample(weights = c(1, 1, 2), size = 10, sampInts)
#' sampInts
#'# [1] 1 1 2 2 2 3 3 3 3 3
#' rankSample(weights = c(1, 1, 2), size = 10000, sampInts)
#' table(sampInts)
#' #sampInts
#' #   1    2    3 
#' #2429 2498 5073 
#'
#' #Used in a nimbleFunction
#' sampGen <- nimbleFunction(setup = function(){
#' 	x = 1:2
#' },
#' run = function(weights = double(1), k = integer() ){
#' 	rankSample(weights, k, x)
#' 	returnType(integer(1))
#' 	return(x)
#' })
#' rSamp <- sampGen()
#' cSamp <- compileNimble(rSamp)
#' cSamp(1:4, 5)
#' #[1] 1 1 4 4 4
rankSample <- function(weights, size, output){
	if(!is.loaded('rankSample') )
		stop('rankSample does not work because we have not loaded any nimble code yet')
	assign(as.character(substitute(output) ), value = .Call('rankSample', as.numeric(weights), as.integer(size), as.integer(output)), envir = parent.frame(n = 1) ) 
}

#' print function for use in nimbleFunctions, where it is identical to \code{print}, but not R's \code{print}.
#'
#' print function for use in nimbleFunctions, where it is identical to \code{print}.  Works in R or NIMBLE, not quite identically.
#'
#' @param ...  an abitrary set of arguments that will be printed in sequence
#'
#' @details The keyword \code{print} in nimbleFunction run-time code will be automatically turned into \code{nimPrint}.  This is a function that prints its arguments int order using \code{cat} in R, or using \code{std::cout} in C++ code generated from compiling nimbleFunctions.
#' Non-scalar numeric objects can be included, although their output will be formatted slightly different in uncompiled and compiled nimbleFunctions.
#'
#'
#' @examples
#' ans <- matrix(1:4, nrow = 2) ## R code, not NIMBLE code
#' nimPrint('Answer is ', ans, '\n') ## would work in R or NIMBLE
nimPrint <- function(...) {
    items <- list(...)
    for(i in seq_along(items)) {if(is.array(items[[i]])) print(items[[i]]) else cat(items[[i]])}
    cat('\n')
}

#' Explicitly declare a variable in run-time code of a nimbleFunction
#'
#' Explicitly declare a variable in run-time code of a nimbleFunction, for cases when its dimensions cannot be inferred before it is used.  Works in R and NIMBLE.
#'
#' @param name     Name of a variable to declare, without quotes
#' @param def      NIMBLE type declaration, of the form \code{TYPE(nDim, sizes)}, where \code{TYPE} is \code{integer}, \code{double}, or \code{logical}, \code{nDim} is the number of dimensions, and \code{sizes} is an optional vector of sizes concatenated with \code{c}.  If \code{nDim} is omitted, it defaults to 0, indicating a scalar.  If sizes are provided, they should not be changed subsequently in the function, including by assignment.  Omitting \code{nDim} results in a scalar.  For \code{logical}, only scalar is currently supported.
#'
#' @author NIMBLE development team
#' @export
#' @details
#' In a run-time function of a nimbleFunction (either the \code{run} function or a function provided in \code{methods} when calling \code{nimbleFunction}), the dimensionality and numeric type of a variable is inferred when possible from the statement first assigning into it.  E.g. \code{A <- B + C} infers that \code{A} has numeric types, dimensions and sizes taken from \code{B + C}.  However, if the first appearance of \code{A} is e.g. \code{A[i] <- 5}, \code{A} must have been explicitly declared.  In this case, \code{declare(A, double(1))} would make \code{A} a 1-dimensional (i.e. vector) double.
#'
#' When sizes are not set, they can be set by a call to \code{setSize} or by assignment to the whole object.  Sizes are not automatically extended if assignment is made to elements beyond the current sizes.  In compiled nimbleFunctions doing so can cause a segfault and crash the R session.
#'
#'
#' This part of the NIMBLE language is needed for compilation, but it also runs in R.  When run in R, is works by the side effect of creating or modifying \code{name} in the calling environment.
#'
#' @examples
#' declare(A, logical())             ## scalar logical, the only kind allowed
#' declare(B, integer(2, c(10, 10))) ## 10 x 10 integer matrix
#' declare(C, double(3))             ## 3-dimensional double array with no sizes set.
declare <- function(name, def){
    defCode <- substitute(def)
    name <- substitute(name)
    if(exists(as.character(name), parent.frame())) return(invisible(NULL))
    value <- if(defCode[[1]] == 'logical') FALSE else 0
    if(length(defCode) == 1){
        assign(as.character(name), value, envir = parent.frame() )
        return()
    }
    nDim = eval(defCode[[2]], envir = parent.frame() )
    if(nDim == 0 ){
        assign(as.character(name), value, envir = parent.frame() )
        return()
    }
    dims = rep(1, nDim)
    if(length(defCode) == 3)
        dims = eval(defCode[[3]], envir = parent.frame() )
    if(length(dims) != nDim)
        stop('in declare, dimensions are not declared properly')
    assign(as.character(name), array(value, dim = dims), envir = parent.frame() )
}


is.na.vec <- function(x) any(is.na(x))

is.nan.vec <- function(x) any(is.nan(x))
