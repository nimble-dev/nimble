#' Basic nimbleFunctions for calculate, simulate, and getLogProb with a set of nodes
#'
#' simulate, calculate, or get existing log probabilities for the current values in a NIMBLE model
#'
#' @param model       A NIMBLE model
#' @param nodes       A set of nodes. If none are provided, default is all \code{model$getNodeNames()}
#' @author Perry de Valpine
#' @export
#' @aliases calcNodes getLogProbNodes
#' @details
#' These are basic nimbleFunctions that take a model and set of nodes and return a function that will call \code{calculate}, \code{simulate}, or \code{getLogProb} on those nodes.  Each is equivalent to a direct call from R, but in nimbleFunction form they can be be compiled and can be put into a nimbleFunctionList.  For example, \code{myCalc <- calcNodes(model, nodes); ans <- myCalc()} is equivalent to \code{ans <- model$calculate(nodes)}, but one can also do \code{CmyCalc <- compileNimble(myCalc)} to get a faster version.
#'
#' In nimbleFunctions, for only one set of nodes, it is equivalent or slightly better to simply use \code{model$calculate(nodes)} in the run-time code.  The compiler will process the model-nodes combination in the same way as would occur by creating a specialized \code{calcNodes} in the setup code.  However, if there are multiple sets of nodes, one can do the following:
#'
#' Setup code: \code{myCalcs <- nimbleFunctionList(calcNodes); myCalcs[[1]] <- calcNodes(model, nodes[[1]]); myCalcs[[2]] <- calcNodes[[2]]}
#'
#' Run code: \code{for(i in seq_along(myCalcs)) {ans[i] <- myCalcs[[i]]()} }
#' 
simNodes <- nimbleFunction(
    name = 'simNodes',
    setup = function(model, nodes){
        if(missing(nodes) )
            nodes <- model$getNodeNames()
        else {
            nodes <- model$expandNodeNames(nodes)
            nodes <- model$topologicallySortNodes(nodes)
        }
    },
    run = function(){
        model$simulate(nodes)
    })

#' @rdname simNodes
#' @export
calcNodes <- nimbleFunction(
    name = 'calcNodes',
	setup = function(model, nodes){
		if(missing(nodes) )
                    depNodes <- model$getNodeNames()
                else
                    depNodes <- model$getDependencies(nodes)
	},
	run = function(){
            ans <- model$calculate(depNodes)
            return(ans)
            returnType(double())
	})	

#' @rdname simNodes
#' @export
getLogProbNodes <- nimbleFunction(
    name = 'getLogProbNodes',
	setup = function(model, nodes) {
		if(missing(nodes) )
                    depNodes <- model$getNodeNames()
                else
                    depNodes <- model$getDependencies(nodes)
	},
	run = function(){
            ans <- model$getLogProb(depNodes)
            return(ans)
            returnType(double())
	})

#' Basic nimbleFunctions for using a NIMBLE model with sets of stored values
#'
#' simulate, calculate, or get the existing log probabilities for values in a modelValues object using a NIMBLE model
#'
#' @param model		A nimble model. 
#' @param nodes		A set of nodes. If none are provided, default is all \code{model$getNodeNames()}
#' @param mv		A modelValues object in which multiple sets of model variables and their corresponding logProb values are or will be saved. \code{mv} must include the nodes provided
#' @author Clifford Anderson-Bergman
#' @export
#' @details
#' \code{simNodesMV} simulates values in the given nodes and saves them in \code{mv}. \code{calcNodesMV} calculates these nodes for each row of \code{mv} and returns a vector of the total log probabilities (densities) for each row. \code{getLogProbNodesMV} is like \code{calcNodesMV} without actually doing the calculations.
#'
#' Each of these will expand variables or index blocks and topologically sort them so that each node's parent nodes are processed before itself.
#'
#' \code{getLogProbMV} should be used carefully.  It is generally for situations where the logProb values are guaranteed to have already been calculated, and all that is needed is to query them.  The risk is that a program may have changed the values in the nodes, in which case \code{getLogProbMV} would collect logProb values that are out of date with the node values.
#' 
#' @aliases calcNodesMV getLogProbNodesMV
#' 
#' @section Run time arguments:
#' \itemize{
#'	\item \code{m}. (\code{simNodesMV} only). Number of simulations requested.
#'      \item \code{saveLP}. (\code{calcNodesMV}only). Whether to save the logProb values in \code{mv}.  Should be given as \code{TRUE} unless there is a good reason not to.
#' }
#'
#' @return from \code{simNodesMV}: NULL.  from \code{calcNodesMV} and \code{getLogProbMV}: a vector of the sum of log probabilities (densities) from any stochastic nodes in \code{nodes}.  
#' 
#'	
#' @examples
#' code <- nimbleCode({
#'	for(i in 1:5)
#'	x[i] ~ dnorm(0,1)
#' })
#'
#' myModel <- nimbleModel(code)
#' myMV <- modelValues(myModel)
#'
#' Rsim <- simNodesMV(myModel, myMV)
#' Rcalc <- calcNodesMV(myModel, myMV)
#' Rglp <- getLogProbNodesMV(myModel, myMV)
#' \dontrun{
#'   cModel <- compileNimble(myModel)
#'   Csim <- compileNimble(Rsim, project = myModel)
#'   Ccalc <- compileNimble(Rcalc, project = myModel)
#'   Cglp <- compileNimble(Rglp, project = myModel)
#'   Csim$run(10)
#'   Ccalc$run(saveLP = TRUE)
#'   Cglp$run()	#Gives identical answers to Ccalc because logProbs were saved
#'   Csim$run(10)
#'   Ccalc$run(saveLP = FALSE)
#'   Cglp$run()	  #Gives wrong answers because logProbs were not saved
#'   result <- as.matrix(Csim$mv)
#' }
simNodesMV <- nimbleFunction(
    name = 'simNodesMV',
    setup = function(model, mv, nodes) {
        if(missing(nodes) )
            nodes <- model$getNodeNames()
        else {
            nodes <- model$expandNodeNames(nodes)
            nodes <- model$topologicallySortNodes(nodes)
        }
    },
    run = function(m = integer(0)){
        resize(mv, m)
        for(i in 1:m){
            model$simulate(nodes)
            nimCopy(from = model, to = mv, nodes = nodes, row = i)
        } 
    })

#' @rdname simNodesMV
#' @export
calcNodesMV <- nimbleFunction(
    name = 'calcNodesMV',
	setup = function(model, mv, nodes) {
		if(missing(nodes) )
                    nodes <- depNodes <- model$getNodeNames()
                else 
                    depNodes <- model$getDependencies(nodes)
		logPvec <- rep(0,2)
	},
	run = function(saveLP = logical()){
		m <- getsize(mv)
		setSize(logPvec, m)
		for(i in 1:m){
			nimCopy(from = mv, to = model, nodes = nodes, row = i)
			logPvec[i] <<- model$calculate(depNodes)	
			if(saveLP)
				nimCopy(from = model, to = mv, nodes = depNodes, rowTo = i, logProb = TRUE)
		}
	returnType(double(1))
	return(logPvec)
	})

#' @rdname simNodesMV
#' @export
getLogProbNodesMV <- nimbleFunction(
    name = 'getLogProbNodesMV',
	setup = function(model, mv, nodes) {
		if(missing(nodes) )
                    nodes <- depNodes <- model$getNodeNames()
                else
                    depNodes <- model$getDependencies(nodes)
		logPvec <- rep(0,2)
	},
	run = function(){
		m <- getsize(mv)
		setSize(logPvec, m)
		for(i in 1:m){
			nimCopy(from = mv, to = model, nodes = depNodes, row = i, logProb = TRUE)
			logPvec[i] <<- model$getLogProb(depNodes)	
		}
	returnType(double(1))
	return(logPvec)
	})

#' Create an Identity matrix (Deprecated)
#'
#' Returns a d-by-d identity matrix (square matrix of 0's, with 1's on the main diagnol).
#'
#' This function can be used in NIMBLE run code.  It is deprecated because now one can use diag(d) instead.
#'
#' @param d The size of the identity matrix to return, will return a d-by-d matrix
#'
#' @return A d-by-d identity matrix
#'
#' @author Daniel Turek
#'
#' @examples
#' Id <- identityMatrix(d = 3)
#'
#' @export
identityMatrix <- nimbleFunction(
    name = 'identityMatrix',
    run = function(d = double()) {
        declare(arr, double(2, c(d, d)))
        for(i in 1:d)   for(j in 1:d)   arr[i, j] <- 0
        for(i in 1:d)                   arr[i, i] <- 1
        returnType(double(2))
        return(arr)
    })
