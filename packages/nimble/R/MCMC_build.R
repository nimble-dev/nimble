
#' Create an MCMC function, from an MCMCconf object
#'
#' Accepts a single required argument, which may be of class MCMCconf, or inherit from class modelBaseClass (a NIMBLE model object).  Returns an MCMC function; see details section.
#'
#' @param conf An object of class MCMCconf that specifies the model, samplers, monitors, and thinning intervals for the resulting MCMC function.  See \code{configureMCMC} for details of creating MCMCconf objects.  Alternatively, \code{MCMCconf} may a NIMBLE model object, in which case an MCMC function corresponding to the default MCMC configuration for this model is returned.
#' @param ... Additional arguments to be passed to \code{configureMCMC} if \code{conf} is a NIMBLE model object
#'
#' @author Daniel Turek
#' @export
#' @details
#' Calling buildMCMC(conf) will produce an uncompiled (R) R mcmc function object, say 'Rmcmc'.
#'
#' The uncompiled MCMC function will have arguments:
#'
#' \code{niter}: The number of iterations to run the MCMC.
#'
#' \code{reset}: Boolean specifying whether to reset the internal MCMC sampling algorithms to their initial state (in terms of self-adapting tuning parameters), and begin recording posterior sample chains anew. Specifying \code{reset=FALSE} allows the MCMC algorithm to continue running from where it left off, appending additional posterior samples to the already existing sample chains. Generally, \code{reset=FALSE} should only be used when the MCMC has already been run (default = TRUE).
#'
#' \code{simulateAll}: Boolean specifying whether to simulate into all stochastic nodes.  This will overwrite the current values in all stochastic nodes (default = FALSE).
#'
#' \code{time}: Boolean specifying whether to record runtimes of the individual internal MCMC samplers.  When \code{time=TRUE}, a vector of runtimes (measured in seconds) can be extracted from the MCMC using the method \code{mcmc$getTimes()} (default = FALSE).
#'
#' \code{progressBar}: Boolean specifying whether to display a progress bar during MCMC execution (default = TRUE).  The progress bar can be permanently disabled by setting the system option \code{nimbleOptions(MCMCprogressBar = FALSE)}.
#'
#' Samples corresponding to the \code{monitors} and \code{monitors2} from the MCMCconf are stored into the interval variables \code{mvSamples} and \code{mvSamples2}, respectively.
#' These may be accessed and converted into R matrix objects via:
#' \code{as.matrix(mcmc$mvSamples)}
#' \code{as.matrix(mcmc$mvSamples2)}
#' 
#' After the MCMC has been run, calling the \code{calculateWAIC()} method of the MCMC object will return the WAIC for the model, calculated using the posterior samples from the MCMC run.  For more information about the
#' WAIC algorithm, see \code{\\help{calculateWAIC()}}.
#'
#' The uncompiled (R) MCMC function may be compiled to a compiled MCMC object, taking care to compile in the same project as the R model object, using:
#' \code{Cmcmc <- compileNimble(Rmcmc, project=Rmodel)}
#'
#' The compiled function will function identically to the uncompiled object, except acting on the compiled model object.
#' 
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' Rmodel <- nimbleModel(code)
#' conf <- configureMCMC(Rmodel)
#' Rmcmc <- buildMCMC(conf)
#' Cmodel <- compileNimble(Rmodel)
#' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
#' Cmcmc$run(10000)
#' samples <- as.matrix(Cmcmc$mvSamples)
#' head(samples)
#' }
buildMCMC <- nimbleFunction(
    name = 'MCMC',
    setup = function(conf, ...) {
    	if(inherits(conf, 'modelBaseClass'))
    		conf <- configureMCMC(conf, ...)

    	else if(!inherits(conf, 'MCMCconf')) stop('conf must either be a nimbleModel or a MCMCconf object (created by configureMCMC(...) )')

        model <- conf$model
        my_initializeModel <- initializeModel(model)
        mvSaved <- modelValues(model)
        samplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(conf$samplerConfs))
            samplerFunctions[[i]] <- conf$samplerConfs[[i]]$buildSampler(model=model, mvSaved=mvSaved)

        monitors  <- processMonitorNames(model, conf$monitors)
        monitors2 <- processMonitorNames(model, conf$monitors2)
        thin  <- conf$thin
        thin2 <- conf$thin2
        mvSamplesConf  <- conf$getMvSamplesConf(1)
        mvSamples2Conf <- conf$getMvSamplesConf(2)
        mvSamples <- modelValues(mvSamplesConf)
        mvSamples2 <- modelValues(mvSamples2Conf)
        samplerTimes <- c(0,0) ## establish as a vector
        progressBarLength <- 52  ## multiples of 4 only
        progressBarDefaultSetting <- getNimbleOption('MCMCprogressBar')
        
        calcWAICFunction <- calculateWAIC(model, mvSamples)
    },

    run = function(niter = integer(), reset = logical(default=TRUE), simulateAll = logical(default=FALSE), time = logical(default=FALSE), progressBar = logical(default=TRUE)) {
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
            samplerTimes <<- numeric(length(samplerFunctions) + 1)       ## default inititialization to zero
        } else {
            mvSamples_offset  <- getsize(mvSamples)
            mvSamples2_offset <- getsize(mvSamples2)
            resize(mvSamples,  mvSamples_offset  + niter/thin)
            resize(mvSamples2, mvSamples2_offset + niter/thin2)
            if(dim(samplerTimes)[1] != length(samplerFunctions) + 1)
                samplerTimes <<- numeric(length(samplerFunctions) + 1)   ## first run: default inititialization to zero
        }
        if(niter < progressBarLength+3 | !progressBarDefaultSetting) progressBar <- progressBar & 0  ## cheap way to avoid compiler warning
        if(progressBar) { for(iPB1 in 1:4) { cat('|'); for(iPB2 in 1:(progressBarLength/4)) cat('-') }; print('|'); cat('|') }
        progressBarIncrement <- niter/(progressBarLength+3)
        progressBarNext <- progressBarIncrement
        progressBarNextFloor <- floor(progressBarNext)
        for(iter in 1:niter) {
            checkInterrupt()
            if(time)
                for(i in seq_along(samplerFunctions))
                    samplerTimes[i] <<- samplerTimes[i] + run.time(samplerFunctions[[i]]$run())
            else 
                for(i in seq_along(samplerFunctions))
                    samplerFunctions[[i]]$run()
            if(iter %% thin  == 0)
                nimCopy(from = model, to = mvSamples,  row = mvSamples_offset  + iter/thin,  nodes = monitors)
            if(iter %% thin2 == 0)
                nimCopy(from = model, to = mvSamples2, row = mvSamples2_offset + iter/thin2, nodes = monitors2)
            if(progressBar & (iter == progressBarNextFloor)) { cat('-')
                                                               progressBarNext <- progressBarNext + progressBarIncrement
                                                               progressBarNextFloor <- floor(progressBarNext) }
        }
        if(progressBar) print('|')
    },
    methods = list(
        getTimes = function() {
            returnType(double(1))
            return(samplerTimes[1:(length(samplerTimes)-1)])
        },
        calculateWAIC = function(){
          burnIn <- .1 * getsize(mvSamples)
          burnIn <- round(burnIn)
          returnType(double(0))
          WAIC <- calcWAICFunction$run(burnIn)
          return(WAIC)
        }),
    where = getLoadingNamespace()
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


