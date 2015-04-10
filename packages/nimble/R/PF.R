

pfStepVirtual <- nimbleFunctionVirtual(
    run = function(m = integer())
        returnType(double())
)


pfStep <- nimbleFunction(
    contains = pfStepVirtual,
    setup = function(model, mv, nodes, iNode) {
        notFirst <- iNode != 1
        prevNode <- nodes[if(notFirst) iNode-1 else iNode]
        thisNode <- nodes[iNode]
        prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
        thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
        thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    },
    run = function(m = integer()) {
        browser()
        returnType(double())
        declare(wts, double(1, m))
        declare(ids, integer(1, m))
        for(i in 1:m) {
            if(notFirst) { copy(mv, model, nodes = 'xs', nodesTo = prevNode, row = i)
                           calculate(model, prevDeterm) }
            simulate(model, thisNode)
            copy(model, mv, nodes = thisNode, nodesTo = 'x', row = i)
            calculate(model, thisDeterm)
            wts[i] <- calculate(model, thisData)
        }
        rankSample(wts, m, ids)
        for(i in 1:m)
            copy(mv, mv, nodes = 'x', nodesTo = 'xs', row = ids[i], rowTo = i)
        return(log(mean(exp(wts))))
    }#,  where = getLoadingNamespace()
)


PF <- nimbleFunction(
    setup = function(model, nodes) {
        my_initializeModel <- initializeModel(model)
        nodes <- model$expandNodeNames(nodes, sort = TRUE)
        dims <- lapply(nodes, function(n) nimbleDim(model[[n]]))
        if(length(unique(dims)) > 1) stop('sizes or dimension of latent states varies')
        mv <- modelValues(modelValuesSpec(vars = c('x', 'xs'),
                                          type = c('double', 'double'),
                                          size = list(x = dims[[1]], xs = dims[[1]])))
        pfStepFunctions <- nimbleFunctionList(pfStepVirtual)
        for(iNode in seq_along(nodes))
            pfStepFunctions[[iNode]] <- pfStep(model, mv, nodes, iNode)
    },
    run = function(m = integer(default = 10000)) {
        returnType(double())
        my_initializeModel$run()
        resize(mv, m)
        logL <- 0
        for(iNode in seq_along(pfStepFunctions))
            logL <- logL + pfStepFunctions[[iNode]]$run(m)
        return(logLik)
    }#,  where = getLoadingNamespace()
)








## buildPF <- nimbleFunctionSimple(
    
##     compileArgs = list(model, latentNodes),
##     args = list(m = int()),
##     returnType = double(),
    
##     setupCode = {
##         latentNodes   <- processNodeList(latentNodes, model)
##         numNodes      <- length(latentNodes)
##         endStochNodes <- getPositions(model, type = 'endStoch')
##         calcNodesList <- lapply(as.list(latentNodes), function(node) intersect(getDependencies(model,node), endStochNodes))
        
##         mv <- modelValues(list(original=double(), sample=double(), weight=double()))
##         samplerFunction <- resample(mv=mv, original='original', sample='sample', weight='weight')
##     },
    
##     code = {
##         lik     <- double()
##         likSum  <- double()
##         meanLik <- double(dim=1, size=numNodes)
        
##         resize(mv, m)
        
##         for(k in seq_along(latentNodes)) {
##             likSum <- 0
##             for(i in runtime(1:m)) {
##                 if(k > 1)          {  copy(mv[[i]], model, 'sample', latentNodes[k-1])  }
##                 simulate(model, latentNodes[k])
##                 copy(model, mv[[i]], latentNodes[k], 'original')
##                 lik <- calculate(model, calcNodesList[[k]])
##                 lik <- exp(lik)
##                 mv[[i]]$weight <- lik
##                 likSum <- likSum + lik
##             }
##             if(runtime(likSum == 0))  { return(-Inf) }
##             meanLik[k] <- likSum / m
##             if(k < numNodes)     {  samplerFunction()  }
##         }
        
##         ## calculate logs, add them together, to estimate total marginal log-likelihood
##         ##     for(k in runtime(1:numNodes))     {  meanLik[k] <- log(meanLik[k])  }
##         meanLik <- log(meanLik)
##         lik <- sum(meanLik)
##         return(lik)
##     }
    
## )






## #' nimbleFunction for a basic particle filter
## #'
## #' Builds a particle filter for a scalar state-space model.
## #' 
## #' @param model					A state space model for which the particle filter will be applied
## #' @param orderedNodeVector		A character vector of the hidden state space for which the particles will be simulated for. Must be correctly
## #' ordered, i.e. \code{ c('x[1]', 'x[2]', ...) }
## #' @author Clifford Anderson-Bergman
## #'
## #' @details  A particle filter approximates the log-likelihood of a state-space model using simulations. At each time step, a sample of latent state values is simulated forward, weighted by the probability density of the observation, and resampled according to those weights for the next time step.  The average of the weights is a factor in the likelihood.  This version is for scalar states and observations.
## #' 
## #' 
## #' @section Run time arguments:
## #' \itemize{
## #'	\item{\code{m} } {number of particles to use for each time step}	
## #'	}
## #'
## buildParticleFilter <- nimbleFunction(
##     setup = function(model, orderedNodeVector) { 
##     	sampledHiddenValues <- modelValues(buildSymbolTable(c('XSample', 'XWeightedSample'), c('double', 'double'), list(1,1) ))
##         weights = c(0,0) 
##         numSteps <- length(orderedNodeVector)
##         steppers <- nimbleFunctionList(nfPFStepVir, length = numSteps) ## We need a list of functions for each step
##         steppers[[1]] <- nfPFstep(model = model, nodeNext = orderedNodeVector[1], sampledHiddenValues = sampledHiddenValues, prior = TRUE) ## set up the functions to go from  a "now" to a "next"
##         for(i in 2:numSteps) {
##             steppers[[i]] <- nfPFstep(model= model, nodeNow = orderedNodeVector[i-1], nodeNext = orderedNodeVector[i], sampledHiddenValues)
##         }
##     },
##     run = function(m = integer(0)) {
##         resize(sampledHiddenValues, m)
##         logProb <- 0
##         for(i in 1:numSteps) {
##             logProb <- logProb + steppers[[i]]$run(m)
##         }
##         returnType(double())
##         return(logProb)
##     }, 	where = getLoadingNamespace()	
## )


## nfPFStepVir <- nimbleFunctionVirtual(
##     run = function(m = integer(0)) {
##         returnType(double(0))
##     }
## )

## ## This will be the nimbleFunction for one step of a particle filter
## nfPFstep <- nimbleFunction(
##     contains = nfPFStepVir,
##     setup = function(model, nodeNow, nodeNext, sampledHiddenValues, prior = FALSE) {
##         if(prior){
##             nodeNow <- nodeNext		# I think this is necessary because we are looking to create nodeNow, even if it's not included
##         }
##         dataDepNodes <- model$getDependencies(nodeNext, dataOnly = TRUE, self = FALSE)  		
##         determNodesNext <- model$getDependencies(nodeNext, determOnly = TRUE, self = FALSE)
##         determNodesNow <- model$getDependencies(nodeNow, determOnly = TRUE, self = FALSE)

##         weights <- c(0,0)  
##         log_weights <- c(0,0)      
##         sampledRanks <- 1:2
##     },
##     run = function(m = integer(0)) {
##         setSize(weights, m)
##         setSize(log_weights, m)
##     	tot_prob = 0
##         for(i in 1:m) {
##             if(!prior) {
##             	nimCopy(from = sampledHiddenValues, row = i, nodes = 'XWeightedSample', to = model, nodesTo = nodeNow)
##                 calculate(model, determNodesNow)	# fills in the deterministinc nodes dependent on x[i-1]
##             }
##             simulate(model, nodeNext)			# simulates x[i] conditional on x[i-1] from weighted sample
##             calculate(model, determNodesNext)	# to fill deterministic nodes dependent on x[i]
##             nimCopy(from = model, nodes = nodeNext, to = sampledHiddenValues, rowTo = i, nodesTo = 'XSample')		# saves simulated value
##             log_weights[i] <<- calculate(model, dataDepNodes)		#
##             weights[i] <<- exp(log_weights[i])
##             tot_prob <- tot_prob + weights[i]
##         }
##         rankSample(weights, m, sampledRanks)
##         for(i in 1:m){
##             nimCopy(from = sampledHiddenValues, to = sampledHiddenValues, nodes = 'XSample', nodesTo = 'XWeightedSample', row = sampledRanks[i], rowTo = i)
##         }
##         return(log(tot_prob/m))
##         returnType(double(0))
##     }, 	where = getLoadingNamespace()	
## )


