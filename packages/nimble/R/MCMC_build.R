
#' Create an MCMC function, from an MCMCspec object
#' 
#' Accepts a single required argument, which may be of class MCMCspec, or inherit from class modelBaseClass (a NIMBLE model obejct).  Returns an MCMC function; see details section.
#' 
#' @param obj An object of class MCMCspec, which specifys the model, samplers, monitors, and thinning intervals for the resulting MCMC function.  See \code{configureMCMC} for details of creating MCMCspec objects.  Alternatively, \code{obj} may a NIMBLE model object, in which case an MCMC function corresponding to the defult MCMC specification for this model is returned.
#' 
#' @author Daniel Turek
#' @export
#' @details
#' Calling buildMCMC(obj) will produce an R mcmc function object, say 'Rmcmc'.
#'
#' The Rmcmc function will have arguments:
#'
#' niter: The number of iterations to run the MCMC.
#'
#' reset: Boolean specifying whether to reset the model and stored samples.  This will simulate into any stochastic nodes with value NA,
#' propagate values through any deterministic nodes, and calculate all model probabilities.
#' This will also reset the internal stored MCMC samples.
#' Specifying reset=FALSE allows the MCMC algorithm to continue running from where it left off.
#' Generally, reset=FALSE should only be used when the MCMC has already been run.  See examples.
#'
#' simulateAll: Boolean specifying whether to simulate into all stochastic nodes.  This will overwrite the current values in all stochastic nodes.
#' 
#' Samples corresponding to the 'monitors' and 'monitors2' from the MCMCspec are stored into the interval variables 'mvSamples' and 'mvSamples2', respectively.
#' These may be accessed via:
#' Rmcmc$mvSamples
#' Rmcmc$mvSamples2
#' 
#' The Rmcmc function may be compiled to a C MCMC object, taking care to compile in the same project as the R model object, using:
#' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
#' 
#' The Cmcmc function will function identically to the Rmcmc object, except acting on the C model object.
#' @examples
#' code <- nimbleCode({
#'  mu ~ dnorm(0, 1)
#'  x ~ dnorm(mu, 1)
#' })
#' Rmodel <- nimbleModel(code)
#' spec <- configureMCMC(Rmodel)
#' Rmcmc <- buildMCMC(spec)
#' Rmcmc$run(10)
#' samples <- Rmcmc$mvSamples
#' samples[['x']]
#' Rmcmc$run(100, reset = FALSE)
buildMCMC <- nimbleFunction(
    
    setup = function(mcmcspec, ...) {
    	if(inherits(mcmcspec, 'modelBaseClass'))
    		mcmcspec <- configureMCMC(mcmcspec, ...)
    	
    	else if(!inherits(mcmcspec, 'MCMCspec'))	
            stop('mcmcspec must either be a nimbleModel or a MCMCspec object (created by configureMCMC(...) )')
        
        model <- mcmcspec$model
        
        RHSonlyNodes <- model$getMaps('nodeNamesRHSonly')
        hasRHSonlyNodes <- length(RHSonlyNodes) > 0

        AllDetermNodes <- model$getNodeNames(determOnly = TRUE)
        determNodesNodeFxnVector <- nodeFunctionVector(model = model, nodeNames = AllDetermNodes)
        
        stochNonDataNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        
        initFunctionList <- nimbleFunctionList(mcmcNodeInit_virtual)
        tot_length = hasRHSonlyNodes + length(stochNonDataNodes)

        iter = 1
        if(hasRHSonlyNodes){
            initFunctionList[[iter]] <- mcmcCheckRHS_Init(model = model, node = RHSonlyNodes)
            iter = iter + 1
        }
        
        for(i in seq_along(stochNonDataNodes)){
            initFunctionList[[iter + i - 1]] <- mcmcStochNode_Init(model, stochNonDataNodes[i])
        }
        
        mvSaved <- modelValues(model)
        samplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(mcmcspec$samplerSpecs))     { samplerFunctions[[i]] <- mcmcspec$samplerSpecs[[i]]$buildSampler(model=model, mvSaved=mvSaved) }
        
        monitors  <- mcmcspec$monitors
        monitors2 <- mcmcspec$monitors2
        thin  <- mcmcspec$thin
        thin2 <- mcmcspec$thin2
        mvSamples  <- mcmcspec$newMvSamples()
        mvSamples2 <- mcmcspec$newMvSamples2()
    },
    
    run = function(niter = integer(), reset = logical(default=TRUE), simulateAll = logical(default=FALSE)) {
        if(simulateAll)     simulate(model)    ## default behavior excludes data nodes
        if(reset) {
            for(i in seq_along(initFunctionList))  {   
            	calculate(nodeFxnVector = determNodesNodeFxnVector)
            	initFunctionList[[i]]$run()
   			}
            calculate(model)
            nimCopy(from = model, to = mvSaved, row = 1, logProb = TRUE)
            for(i in seq_along(samplerFunctions))      {   samplerFunctions[[i]]$reset()   }
            mvSamples_offset  <- 0
            mvSamples2_offset <- 0
            resize(mvSamples,  niter/thin)
            resize(mvSamples2, niter/thin2)
        } else {
            mvSamples_offset  <- getsize(mvSamples)
            mvSamples2_offset <- getsize(mvSamples2)
            resize(mvSamples,  mvSamples_offset  + niter/thin)
            resize(mvSamples2, mvSamples2_offset + niter/thin2)
        }
        
        for(iter in 1:niter) {
            for(i in seq_along(samplerFunctions))      {   samplerFunctions[[i]]$run()   }
            if(iter %% thin  == 0) { nimCopy(from = model, to = mvSamples,  row = mvSamples_offset  + iter/thin,  nodes = monitors) }
            if(iter %% thin2 == 0) { nimCopy(from = model, to = mvSamples2, row = mvSamples2_offset + iter/thin2, nodes = monitors2) }
        }
    },  where = getLoadingNamespace()
)






