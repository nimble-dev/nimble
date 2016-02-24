###		These functions are used for calculate/sim/getLP for the nodeFunctionVectors
###		Can either enter model, nodes or model_nodes

#' Halt execution of a nimbleFunction function method.  Part of the NIMBLE language
#'
#' @param msg Character object to be output as an error message
#'
#' @author Perry de Valpine
#' @export
#' @details
#' The NIMBLE stop is similar to the native R stop, but it takes only one argument, the error message to be output.  During uncompiled NIMBLE execution, nimStop simply calls R's stop funtion. During compiled execution it calls the error function from the R headers.  stop is an alias for nimStop in the NIMBLE language  
nimStop <- function(msg) stop(msg, call. = FALSE)
# we use call.=FALSE because otherwise the error msg indicates the
# error itself occurs in nimStop() and not in the calling frame

#' Check for interrupt (e.g. Ctrl-C) during nimbleFunction execution. Part of the NIMBLE language.
#'
#' @author Perry de Valpine
#' @export
#' @details
#' During execution of nimbleFunctions that take a long time, it is nice to occassionally check if the user has entered an interrupt and bail out of execution if so.  This function does that.  During uncompiled nimbleFunction execution, it does nothing.  During compiled execution, it calls R_checkUserInterrupt() of the R headers. 
checkInterrupt <- function() {}

#' Turn a numeric vector into a single-row or single-column matrix
#'
#' Turns a numeric vector into a matrix that has 1 row or 1 column.  Part of NIMBLE language.
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

#' @export
makeParamInfo <- function(model, node, param) {
    distInfo <- getDistribution(model$getNodeDistribution(node))
    ans <- c(list(paramID = distInfo$paramIDs[param]), distInfo$types[[param]])
    class(ans) <- 'getParam_info'
    ans
}

#' Get value of a parameter of a stochastic node in a model
#'
#' Part of the NIMBLE language
#'
#' @param model A NIMBLE model object
#'
#' @param node  The name of a stochastic node in the model
#'
#' @param param The name of a parameter for the node
#' 
#' @export
#' @details For example, suppose node 'x[1:5]' follows a multivariate
#' normal distribution (dmnorm) in a model declared by BUGS code.
#' getParam(model, 'x[1:5]', 'mean') would return the current value of
#' the mean parameter (which may be determined from other nodse).  The
#' parameter requested does not have to be part of the
#' parameterization used to declare the node.  Rather, it can be any
#' parameter known to the distribution.  For example, one can request
#' the scale or rate parameter of a gamma distribution, regardless of
#' which one was used to declare the node.
getParam <- function(model, node, param) {
    if(missing(param)) { ## already converted by keyword conversion
        nodeFunction <- model
        paramInfo <- node
    } else {
        ## not already converted
        nodeFunction <- model$nodes[[node]]
        paramInfo <- makeParamInfo(model, node, param)
    }
    paramID <- paramInfo$paramID
    nDim <- paramInfo$nDim
    type <- paramInfo$type
    funName <- paste0('getParam_',nDim,'D_',type)
    ans <- eval(substitute(nodeFunction$FUNNAME(paramID), list(FUNNAME = as.name(funName))))
    return(ans)
}

#' @export
nimSwitch <- function(paramID, IDoptions, ...) {
    dotsList <- eval(substitute(alist(...)))
    iUse <- which(IDoptions == paramID)
    eval(dotsList[[iUse]], envir = parent.frame())
    invisible(NULL)
}

rCalcNodes <- function(model, nodes){
    l_Prob = 0
    
    if(inherits(model, 'CmodelBaseClass') & getNimbleOption('useMultiInterfaceForNestedNimbleFunctions')) 
        for(nName in nodes)
            l_Prob = l_Prob + model$nodes[[nName]][[1]]$callMemberFunction(model$nodes[[nName]][[2]], 'calculate')
    else
        for(nName in nodes)
            l_Prob = l_Prob + model$nodes[[nName]]$calculate()
    
    return(l_Prob)
}

rCalcDiffNodes <- function(model, nodes){
    l_Prob <- 0
    if(inherits(model, 'CmodelBaseClass') & getNimbleOption('useMultiInterfaceForNestedNimbleFunctions')) 
        for(nName in nodes)
            l_Prob = l_Prob + model$nodes[[nName]][[1]]$callMemberFunction(model$nodes[[nName]][[2]], 'calculateDiff')
    else
        for(nName in nodes)
            l_Prob = l_Prob + model$nodes[[nName]]$calculateDiff()
    return(l_Prob)
}


#' calculate, calculateDiff, simulate, or get the current log probabilities (densities) a set of nodes in a NIMBLE model
#'
#' calculate, calculateDiff, simulate, or get the current log probabilities (densities) of one or more nodes of a NIMBLE model and (for calculate and getLogProb) return the sum of their log probabilities (or densities).  Part of R and NIMBLE.
#' @name nodeFunctions
#' 
#' @param model        A NIMBLE model, either the compiled or uncompiled version
#' @param nodes        A character vector of node names, with index blocks allowed, such as 'x', 'y[2]', or 'z[1:3, 2:4]'
#' @param nodeFxnVector An optional vector of nodeFunctions on which to operate, in lieu of \code{model} and \code{nodes}
#' @param includeData  A logical argument specifying whether \code{data} nodes should be simulated into (only relevant for \link{simulate}
#' @author NIMBLE development team
#' @export
#' @details
#' These functions expands the nodes and then process them in the model in the order provided.  Expanding nodes means turning 'y[1:2]' into c('y[1]','y[2]') if y is a vector of scalar nodes.
#' Calculation is defined for a stochastic node as executing the log probability (density) calculation and for a deterministic node as calculating whatever function was provided on the right-hand side of the model declaration.
#'
#' Difference calculation (calculateDiff) executes the operation(s) on the model as calculate, but it returns the sum of the difference between the new log probabilities and the previous ones.
#' 
#' Simulation is defined for a stochastic node as drawing a random value from its distribution, and for deterministic node as equivalent to calculate.
#'
#' getLogProb simply collects the sum of the log probabilities of nodes if they are known to have already been calculated.
#'
#' These functions can be used from R or in NIMBLE run-time functions that will be compiled.  When executed in R (including when an uncompiled nimbleFunction is executed), they can be slow because the nodes are expanded each time.  When compiled in NIMBLE, the nodes are expanded only once during compilation, so execution will be much faster.
#'
#' It is common to want the nodes to be provided in topologically sorted order, so that they will be calculated or simulated following the order of the model graph.  Functions such as model$getDependencies(nodes, ...) return nodes in topologically sorted order.  They can be directly sorted by model$topologicallySortNodes(nodes), but if so it is a good idea to expand names first by model$topologicallySortNodes(model$expandNodeNames(nodes))
#'
#' @return calculate and getLogProb return the sum of the log probabilities (densities) of the calculated nodes, with a contribution of 0 from any deterministic nodes
#'
#' @return calculateDiff returns the sum of the difference between the new and old log probabilities (densities) of the calculated nodes, with a contribution of 0 from any deterministic nodes.
#'
#' simulate returns NULL.
#' 
NULL

#' @rdname nodeFunctions
#' @export
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

#' @rdname nodeFunctions
#' @export
calculateDiff <- function(model, nodes, nodeFxnVector)		
{
    if(!missing(nodeFxnVector)){
        model <- nodeFxnVector$model
        nodes <- nodeFxnVector$getNodeNames()
        return(rCalcDiffNodes(model, nodes))
    }
    if(inherits(model, 'modelBaseClass') ){
        if(missing(nodes) ) 
            nodes <- model$getMaps('nodeNamesLHSall')
        nfv <- nodeFunctionVector(model, nodes)
        nodeNames <- nfv$getNodeNames()
        return(rCalcDiffNodes(model, nodeNames))
    }	
}

rGetLogProbsNodes <- function(model, nodes){
    l_Prob = 0

    if(inherits(model, 'CmodelBaseClass') & getNimbleOption('useMultiInterfaceForNestedNimbleFunctions')) 
        for(nName in nodes)
            l_Prob = l_Prob + model$nodes[[nName]][[1]]$callMemberFunction(model$nodes[[nName]][[2]], 'getLogProb')
    else
        for(nName in nodes)
            l_Prob = l_Prob + model$nodes[[nName]]$getLogProb()
    return(l_Prob)
}

#' @rdname nodeFunctions
#' @export
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
    if(inherits(model, 'CmodelBaseClass') & getNimbleOption('useMultiInterfaceForNestedNimbleFunctions')) 
        for(nName in nodes)
            model$nodes[[nName]][[1]]$callMemberFunction(model$nodes[[nName]][[2]], 'simulate')
    else 
        for(nName in nodes)
            model$nodes[[nName]]$simulate()
}

#' @rdname nodeFunctions
#' @export
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


getValues <- function(vals, model, nodes, envir = parent.frame()) {
    valsExp = substitute(vals)
    access = modelVariableAccessorVector(model, nodes, logProb = FALSE)
    output = getValuesAccess(access)
    if(is.name(valsExp))
        assign(as.character(valsExp), output, envir = parent.frame(n = 1) ) ## again this assumes valsExp has no indexing -- need to fix
    else {
        varName <- valsExp[[2]]
        orig <- get(varName, envir = envir)
        assignTarget <- valsExp
        assignTarget[[2]] <- quote(orig)
        eval(substitute(AT <- output, list(AT = assignTarget)))
        assign(varName, orig, envir = envir)
    }
}

getValuesAccess <- function(access) {
    fromCode <- makeGetCodeFromAccessorVector(access)
    if(length(fromCode)==0) return(numeric())
    sourceFromObject <- access[[1]]
    unlist(lapply(fromCode, function(i) eval(i)))

##    if(access$numAccessors==0) return(numeric()) ## NEW ACCESSORS
##    unlist(lapply(1:access$numAccessors, function(i) access$getValues(i)))

}


setValuesAccess <- function(input, access) {

    toCode <- makeSetCodeFromAccessorVector(access)
    mapInfo <- makeMapInfoFromAccessorVector(access) ## inefficient to do both of these, but we need the lengths!
    if(length(toCode)==0) return(NULL)
    sourceToObject <- access[[1]]
    nextIndex <- 0
    ## not easy to check lengths any more
    for(i in 1:length(toCode)) {
        nextLength <- mapInfo[[i]]$length
        oneValue <- input[nextIndex + (1:nextLength)]
        eval(toCode[[i]])
        nextIndex <- nextIndex + nextLength
    }
    invisible(NULL)

    ## if(access$numAccessors==0) return(NULL) ## NEW ACCESSORS
    ## nextIndex <- 0
    ## if(access$getLength() != length(input)) {
    ##     writeLines('Length of input does not match accessor')
    ##     if(access$getLength() > length(input)) stop('Bailing out because not enough values were provided for accessor')
    ##     writeLines('Too many input values provided.  Continuing anyway')
    ## }
    ## for(i in 1:length(access$code)) {
    ##     nextLength <- access$getLength(i)
    ##     access$setValues(i, input[nextIndex + (1:nextLength)])
    ##     nextIndex <- nextIndex + nextLength
    ## }
    ## invisible(NULL)
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
	getValues(ans, model, nodes, parent.frame())
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
#' \dontrun{
#' cModel <- compileNimble(rModel)
#' cCopy <- compileNimble(rCopy, project = rModel)
#' cModel[['x']] <- rnorm(100)
#' 
#' cCopy$run() ## execute the copy with the compiled function
#' }
nimCopy <- function(from, to, nodes = NULL, nodesTo = NULL, row = NA, rowTo = NA, logProb = FALSE){
    if(is.null(nodes) )
        nodes = from$getVarNames(includeLogProb = logProb) ## allNodeNames(from)
    if( inherits(from, "modelBaseClass") ){
        accessFrom = modelVariableAccessorVector(from, nodes, logProb = logProb)
    }
    else
        if(inherits(from, "modelValuesBaseClass")) {
            accessFrom = modelValuesAccessorVector(from, nodes, logProb = logProb)
            if(is.na(row))
                stop("Error: need to supply 'row' for a modelValues copy")
            ##accessFrom$setRow(row) ## NEW ACCESSORS
        }
        else stop('argument "from" in nimCopy is neither a model nor modelValues')

    if( inherits(to, "modelBaseClass") ){
        if(is.null(nodesTo) ) 
            accessTo = modelVariableAccessorVector(to, nodes, logProb = logProb)
        else
            accessTo = modelVariableAccessorVector(to, nodesTo, logProb = logProb)
    }
    else
        if(inherits(to, "modelValuesBaseClass")) {
            if(is.null(nodesTo) ) 
                accessTo = modelValuesAccessorVector(to, nodes, logProb = logProb)
            else
                accessTo = modelValuesAccessorVector(to, nodesTo, logProb = logProb)
            if(is.na(rowTo))
                rowTo = row
            ##accessTo$setRow(rowTo) ## NEW ACCESSORS
        }
        else stop('argument "to" in nimCopy is neither a model nor modelValues')

    sourceToObject <- accessTo[[1]]
    sourceFromObject <- accessFrom[[1]]
    setCode <- makeSetCodeFromAccessorVector(accessTo)
    getCode <- makeGetCodeFromAccessorVector(accessFrom)
    lengthTo <- length(setCode)
    if(length(getCode) != lengthTo)
    stop('unequal number of entries in nimCopy') 
    if(lengthTo > 0){
        for(i in 1:lengthTo){
            oneValue <- eval(getCode[[i]]) ## may have row hardwired in
            eval(setCode[[i]]) ## oneValue is hard-wired in. may have rowTo hardwired in
        }
    }
    
    ## lengthTo <- length(accessTo$code) ## NEW ACCESSORS
    ## if(length(accessFrom$code) != lengthTo)
    ##     stop('unequal number of entries in nimCopy') 
    ## if(lengthTo > 0){
    ##     for(i in 1:lengthTo){
    ##         accessTo$setValues( i, accessFrom$getValues(i) )
    ##     }
    ## }
}


#' Internal way to access or set a member variable of a nimbleFunction created during \code{setup}.  Normally in NIMBLE code you would use \code{nf$var} instead of \code{nfVar(nf, var)}. 
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
#' When \code{nimbleFunction} is called and a \code{setup} function is provided, then \code{nimbleFunction} returns a function.  That function is a generator that should be called with arguments to the \code{setup} function and returns another function with \code{run} and possibly other member functions.  The member functions can use objects created or passed to \code{setup}.  During internal processing, the NIMBLE compiler turns some cases of \code{nf$var} into \code{nfVar(nf, var)}. These provide direct access to setup variables (member data).  \code{nfVar} is not typically called by a NIMBLE user or programmer.
#'
#'
#' For internal access to methods of \code{nf}, see \link{nfMethod}.
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
#' nf2$run()
#' Cnf2 <- compileNimble(nf2)
#' Cnf2$run()
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

#' Internal function for accessing a member function (method) of a nimbleFunction.  Normally a user will write \code{nf$method(x)} instead of \code{nfMethod(nf, method)(x)}.
#'
#' access (call) a member function of a nimbleFunction, including \code{run}.
#'
#' @param nf          a specialized nimbleFunction, i.e. one that has already had setup parameters processed
#' @param methodName  a character string giving the name of the member function to call
#'
#' @author NIMBLE development team
#' @export
#' @details
#' nimbleFunctions have a default member function called \code{run}, and may have other member functions provided via the \code{methods} argument to \code{nimbleFunction}.
#' As an internal step, the NIMBLE compiler turns \code{nf$method(x)} into \code{nfMethod(nf, method)(x)}, but a NIMBLE user or programmer would not normally need to use \code{nfMethod} directly.
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
#' @param size			Size of sample
#' @param output		An R object into which the values will be placed. See example below for proper use
#' @param silent Logical indicating whether to suppress logging information
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
#' cSamp$run(1:4, 5)
#' #[1] 1 1 4 4 4
rankSample <- function(weights, size, output, silent = FALSE) {
    ##cat('in R version rankSample\n')
    if(!is.loaded('rankSample'))
        stop('rankSample does not work because we have not loaded any nimble code yet')
    assign(as.character(substitute(output)), .Call('rankSample', as.numeric(weights), as.integer(size), as.integer(output), as.logical(silent)), envir = parent.frame())
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

nimRound <- round

