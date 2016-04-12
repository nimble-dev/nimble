
#' Create an MCMC function, from an MCMCspec object
#'
#' Accepts a single required argument, which may be of class MCMCspec, or inherit from class modelBaseClass (a NIMBLE model obejct).  Returns an MCMC function; see details section.
#'
#' @param mcmcspec An object of class MCMCspec that specifies the model, samplers, monitors, and thinning intervals for the resulting MCMC function.  See \code{configureMCMC} for details of creating MCMCspec objects.  Alternatively, \code{mcmcspec} may a NIMBLE model object, in which case an MCMC function corresponding to the defult MCMC specification for this model is returned.
#' @param ... Additional arguments to be passed to \code{configureMCMC} if \code{mcmcspec} is a NIMBLE model object
#'
#' @author Daniel Turek
#' @export
#' @details
#' Calling buildMCMC(mcmcspec) will produce an R mcmc function object, say 'Rmcmc'.
#'
#' The Rmcmc function will have arguments:
#'
#' \code{niter}: The number of iterations to run the MCMC.
#'
#' \code{reset}: Boolean specifying whether to reset the model and stored samples.  This will simulate into any stochastic nodes with value NA,
#' propagate values through any deterministic nodes, and calculate all model probabilities.
#' This will also reset the internal stored MCMC samples.
#' Specifying \code{reset=FALSE} allows the MCMC algorithm to continue running from where it left off.
#' Generally, \code{reset=FALSE} should only be used when the MCMC has already been run.  See examples.
#'
#' \code{simulateAll}: Boolean specifying whether to simulate into all stochastic nodes.  This will overwrite the current values in all stochastic nodes.
#'
#' Samples corresponding to the \code{monitors} and \code{monitors2} from the MCMCspec are stored into the interval variables \code{mvSamples} and \code{mvSamples2}, respectively.
#' These may be accessed via:
#' \code{Rmcmc$mvSamples}
#' \code{Rmcmc$mvSamples2}
#'
#' The Rmcmc function may be compiled to a C MCMC object, taking care to compile in the same project as the R model object, using:
#' \code{Cmcmc <- compileNimble(Rmcmc, project=Rmodel)}
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

    	else if(!inherits(mcmcspec, 'MCMCspec')) stop('mcmcspec must either be a nimbleModel or a MCMCspec object (created by configureMCMC(...) )')

        model <- mcmcspec$model
        my_initializeModel <- initializeModel(model)
        mvSaved <- modelValues(model)
        samplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(mcmcspec$samplerSpecs))
            samplerFunctions[[i]] <- mcmcspec$samplerSpecs[[i]]$buildSampler(model=model, mvSaved=mvSaved)

        monitors  <- processMonitorNames(model, mcmcspec$monitors)
        monitors2 <- processMonitorNames(model, mcmcspec$monitors2)
        thin  <- mcmcspec$thin
        thin2 <- mcmcspec$thin2
        mvSamplesSpec  <- mcmcspec$getMvSamplesSpec(1)
        mvSamples2Spec <- mcmcspec$getMvSamplesSpec(2)
        mvSamples <- modelValues(mvSamplesSpec)
        mvSamples2 <- modelValues(mvSamples2Spec)
    },

    run = function(niter = integer(), reset = logical(default=TRUE), simulateAll = logical(default=FALSE)) {
        if(simulateAll)     simulate(model)    ## default behavior excludes data nodes
        my_initializeModel$run()
        if(reset) {
            nimCopy(from = model, to = mvSaved, row = 1, logProb = TRUE)
            for(i in seq_along(samplerFunctions))
                samplerFunctions[[i]]$reset()
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
            for(i in seq_along(samplerFunctions))
                samplerFunctions[[i]]$run()
            if(iter %% thin  == 0)
                nimCopy(from = model, to = mvSamples,  row = mvSamples_offset  + iter/thin,  nodes = monitors)
            if(iter %% thin2 == 0)
                nimCopy(from = model, to = mvSamples2, row = mvSamples2_offset + iter/thin2, nodes = monitors2)
        }
    },  where = getLoadingNamespace()
)


# This is a function that will weed out missing indices from the monitors
processMonitorNames <- function(model, nodes){
	isLogProbName <- grepl('logProb', nodes)
	expandedNodeNames <- model$expandNodeNames(nodes[!isLogProbName])
	origLogProbNames <- nodes[isLogProbName]
	expandedLogProbNames <- character()
	if(length(origLogProbNames) > 0){
		nodeName_fromLogProbName <- gsub('logProb_', '', origLogProbNames)
		expandedLogProbNames <- model$modelDef$nodeName2LogProbName(nodeName_fromLogProbName)
	}
	return( c(expandedNodeNames, expandedLogProbNames) )
}


