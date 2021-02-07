
#' Create an MCMC function from a NIMBLE model, or an MCMC configuration object
#'
#' First required argument, which may be of class \code{MCMCconf} (an MCMC configuration object), or inherit from class \code{modelBaseClass} (a NIMBLE model object).  Returns an uncompiled executable MCMC function.  See details.
#'
#' @param conf An MCMC configuration object of class \code{MCMCconf} that specifies the model, samplers, monitors, and thinning intervals for the resulting MCMC function.  See \code{configureMCMC} for details of creating MCMC configuration objects.  Alternatively, \code{conf} may a NIMBLE model object, in which case an MCMC function corresponding to the default MCMC configuration for this model is returned.
#' 
#' @param ... Additional arguments to be passed to \code{configureMCMC} if \code{conf} is a NIMBLE model object
#'
#' @details
#' 
#' Calling buildMCMC(conf) will produce an uncompiled MCMC function object.  The uncompiled MCMC function will have arguments:
#'
#' \code{niter}: The number of iterations to run the MCMC.
#'
#' \code{thin}: The thinning interval for the \code{monitors} that were specified in the MCMC configuration.  If this argument is provided at MCMC runtime, it will take precedence over the \code{thin} interval that was specified in the MCMC configuration.  If omitted, the \code{thin} interval from the MCMC configuration will be used.
#'
#' \code{thin2}: The thinning interval for the second set of monitors (\code{monitors2}) that were specified in the MCMC configuration.  If this argument is provided at MCMC runtime, it will take precedence over the \code{thin2} interval that was specified in the MCMC configuration.  If omitted, the \code{thin2} interval from the MCMC configuration will be used.
#'
#' \code{reset}: Boolean specifying whether to reset the internal MCMC sampling algorithms to their initial state (in terms of self-adapting tuning parameters), and begin recording posterior sample chains anew. Specifying \code{reset = FALSE} allows the MCMC algorithm to continue running from where it left off, appending additional posterior samples to the already existing sample chains. Generally, \code{reset = FALSE} should only be used when the MCMC has already been run (default = TRUE).
#'
#' \code{resetMV}: Boolean specifying whether to begin recording posterior sample chains anew. This argument is only considered when using \code{reset = FALSE}.  Specifying \code{reset = FALSE, resetMV = TRUE} allows the MCMC algorithm to continue running from where it left off, but without appending the new posterior samples to the already existing samples, i.e. all previously obtained samples will be erased. This option can help reduce memory usage during burn-in (default = FALSE).
#'
#' \code{nburnin}: Number of initial, pre-thinning, MCMC iterations to discard (default = 0).
#'
#' \code{time}: Boolean specifying whether to record runtimes of the individual internal MCMC samplers.  When \code{time = TRUE}, a vector of runtimes (measured in seconds) can be extracted from the MCMC using the method \code{mcmc$getTimes()} (default = FALSE).
#'
#' \code{progressBar}: Boolean specifying whether to display a progress bar during MCMC execution (default = TRUE).  The progress bar can be permanently disabled by setting the system option \code{nimbleOptions(MCMCprogressBar = FALSE)}.
#'
#' Samples corresponding to the \code{monitors} and \code{monitors2} from the MCMCconf are stored into the interval variables \code{mvSamples} and \code{mvSamples2}, respectively.
#' These may be accessed and converted into R matrix or list objects via:
#' \code{as.matrix(mcmc$mvSamples)}
#' \code{as.list(mcmc$mvSamples)}
#' \code{as.matrix(mcmc$mvSamples2)}
#' \code{as.list(mcmc$mvSamples2)}
#' 
#' The uncompiled MCMC function may be compiled to a compiled MCMC object, taking care to compile in the same project as the R model object, using:
#' \code{Cmcmc <- compileNimble(Rmcmc, project = Rmodel)}
#'
#' The compiled function will function identically to the uncompiled object, except acting on the compiled model object.
#'
#' @section Calculating WAIC:
#' 
#' After the MCMC has been run, calling the \code{calculateWAIC()} method of the MCMC object will return the WAIC for the model, calculated using the posterior samples from the MCMC run.
#' 
#' \code{calculateWAIC()} accepts a single arugment:
#'
#' \code{nburnin}: The number of pre-thinning MCMC samples to remove from the beginning of the posterior samples for WAIC calculation (default = 0). These samples are discarded in addition to any burn-in specified when running the MCMC.
#' 
#' The \code{calculateWAIC} method can only be used if the \code{enableWAIC} 
#' argument to \code{configureMCMC} or to \code{buildMCMC} is set to \code{TRUE}, or if the NIMBLE option
#' \code{enableWAIC} is set to \code{TRUE}.  If a user attempts
#' to call \code{calculateWAIC} without having set \code{enableWAIC = TRUE}
#' (either in the call to \code{configureMCMC}, or \code{buildMCMC}, or as a NIMBLE option),
#' an error will occur.  
#' 
#' The \code{calculateWAIC} method calculates the WAIC of the model that the
#' MCMC was performed on. The WAIC (Watanabe, 2010) is calculated from
#' Equations 5, 12, and 13 in Gelman et al. (2014) (i.e., using \emph{p}WAIC2).
#'
#' Note that there is not a unique value of WAIC for a model. The current version of
#' NIMBLE only provides the conditional WAIC, namely the version of WAIC where all
#' parameters directly involved in the likelihood are treated as \eqn{theta}
#' for the purposes of Equation 5 from Gelman et al. (2014). As a result, the user
#' must set the MCMC monitors (via the \code{monitors} argument) to include all stochastic
#' nodes that are parents of any data nodes; by default the MCMC monitors are only
#' the top-level nodes of the model. For more detail on the use of different predictive
#' distributions, see Section 2.5 from Gelman et al. (2014) or Ariyo et al. (2019).

#' Also note that WAIC relies on a partition of the observations, i.e., 'pointwise'
#' prediction. In NIMBLE the sum over log pointwise predictive density values treats
#' each data node as contributing a single value to the sum. When a data node is
#' multivariate, that data node contributes a single value to the sum based on the
#' joint density of the elements in the node. Note that if one wants the WAIC
#' calculation to be based on the joint predictive density for each group of observations
#' (e.g., grouping the observations from each person or unit in a longitudinal
#' data context), one would need to use a multivariate distribution for the
#' observations in each group (potentially by writing a user-defined distribution).
#' 
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#'     y ~ dnorm(x, 1)
#' })
#' Rmodel <- nimbleModel(code, data = list(y = 0))
#' conf <- configureMCMC(Rmodel, monitors = c('mu', 'x'))
#' Rmcmc <- buildMCMC(conf, enableWAIC = TRUE)
#' Cmodel <- compileNimble(Rmodel)
#' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
#' Cmcmc$run(10000)
#' samples <- as.matrix(Cmcmc$mvSamples)
#' samplesAsList <- as.list(Cmcmc$mvSamples)
#' head(samples)
#' WAIC <- Cmcmc$calculateWAIC(nburnin = 1000)
#' }
#'
#' @seealso \code{\link{configureMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
#' 
#' @author Daniel Turek
#' 
#' @references 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. \emph{Journal of Machine Learning Research} 11: 3571-3594.
#' 
#' Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. \emph{Statistics and Computing} 24(6): 997-1016.
#'
#' Ariyo, O., Quintero, A., Munoz, J., Verbeke, G. and Lesaffre, E. (2019). Bayesian model selection in linear mixed models for longitudinal data. \emph{Journal of Applied Statistics} 47: 890-913.
#' @export
buildMCMC <- nimbleFunction(
    name = 'MCMC',
    setup = function(conf, ...) {
    	if(inherits(conf, 'modelBaseClass'))   conf <- configureMCMC(conf, ...)
    	else if(!inherits(conf, 'MCMCconf')) stop('conf must either be a nimbleModel or a MCMCconf object (created by configureMCMC(...) )')
        dotdotdotArgs <- list(...)
        enableWAICargument <- if(!is.null(dotdotdotArgs$enableWAIC)) dotdotdotArgs$enableWAIC else nimbleOptions('MCMCenableWAIC')    ## accept enableWAIC argument regardless
        model <- conf$model
        my_initializeModel <- initializeModel(model)
        mvSaved <- modelValues(model)
        samplerFunctions <- nimbleFunctionList(sampler_BASE)
        samplerFunctionsHMC <- nimbleFunctionList(sampler_HMC_BASE)
        for(i in seq_along(conf$samplerConfs)) {
            newSF <- conf$samplerConfs[[i]]$buildSampler(model=model, mvSaved=mvSaved)
            samplerFunctions[[i]] <- newSF
            if(conf$samplerConfs[[i]]$name %in% c('HMC', 'HMC2'))
                samplerFunctionsHMC[[length(samplerFunctionsHMC)+1]] <- newSF
        }
        samplerExecutionOrderFromConfPlusTwoZeros <- c(conf$samplerExecutionOrder, 0, 0)  ## establish as a vector
        monitors  <- mcmc_processMonitorNames(model, conf$monitors)
        monitors2 <- mcmc_processMonitorNames(model, conf$monitors2)
        thinFromConfVec <- c(conf$thin, conf$thin2)  ## vector
        thinToUseVec <- c(0, 0)                      ## vector, needs to member data
        mvSamplesConf  <- conf$getMvSamplesConf(1)
        mvSamples2Conf <- conf$getMvSamplesConf(2)
        mvSamples <- modelValues(mvSamplesConf)
        mvSamples2 <- modelValues(mvSamples2Conf)
        samplerTimes <- c(0,0) ## establish as a vector
        progressBarLength <- 52  ## multiples of 4 only
        progressBarDefaultSetting <- getNimbleOption('MCMCprogressBar')
        nimbleVerboseOption <- getNimbleOption('verbose')
        ## WAIC setup:
        dataNodes <- model$getNodeNames(dataOnly = TRUE)
        dataNodeLength <- length(dataNodes)
        sampledNodes <- model$getVarNames(includeLogProb = FALSE, nodes = monitors)
        sampledNodes <- sampledNodes[sampledNodes %in% model$getVarNames(includeLogProb = FALSE)]
        paramDeps <- model$getDependencies(sampledNodes, self = FALSE, downstream = TRUE)
        allVarsIncludingLogProbs <- model$getVarNames(includeLogProb = TRUE)
        enableWAIC <- enableWAICargument || conf$enableWAIC   ## enableWAIC comes from MCMC configuration, or from argument to buildMCMC
        if(enableWAIC) {
            if(dataNodeLength == 0)   stop('WAIC cannot be calculated, as no data nodes were detected in the model.')
            mcmc_checkWAICmonitors_conditional(model = model, monitors = sampledNodes, dataNodes = dataNodes)
        }
    },

    run = function(
        niter                 = integer(),
        reset                 = logical(default = TRUE),
        resetMV               = logical(default = FALSE), ## Allows resetting mvSamples when reset==FALSE
        time                  = logical(default = FALSE),
        progressBar           = logical(default = TRUE),
        ## reinstate samplerExecutionOrder as a runtime argument, once we support non-scalar default values for runtime arguments:
        ##samplerExecutionOrder = integer(1, default = -1)
        nburnin               = integer(default =  0),
        thin                  = integer(default = -1),
        thin2                 = integer(default = -1),
        chain                 = integer(default =  1)) {
        if(niter < 0)       stop('cannot specify niter < 0')
        if(nburnin < 0)     stop('cannot specify nburnin < 0')
        if(nburnin > niter) stop('cannot specify nburnin > niter')
        thinToUseVec <<- thinFromConfVec
        if(thin  != -1)   thinToUseVec[1] <<- thin
        if(thin2 != -1)   thinToUseVec[2] <<- thin2
        my_initializeModel$run()
        nimCopy(from = model, to = mvSaved, row = 1, logProb = TRUE)
        if(reset) {
            samplerTimes <<- numeric(length(samplerFunctions) + 1)       ## default inititialization to zero
            for(i in seq_along(samplerFunctions))   samplerFunctions[[i]]$reset()
            if(length(samplerFunctionsHMC) > 0)   for(i in seq_along(samplerFunctionsHMC))   samplerFunctionsHMC[[i]]$initializeWarmup(niter, chain, nimbleVerboseOption)
            mvSamples_copyRow  <- 0
            mvSamples2_copyRow <- 0
        } else {
            if(nburnin != 0)   stop('cannot specify nburnin when using reset = FALSE.')
            if(dim(samplerTimes)[1] != length(samplerFunctions) + 1)   samplerTimes <<- numeric(length(samplerFunctions) + 1)   ## first run: default inititialization to zero
            if (resetMV) {
                mvSamples_copyRow  <- 0
                mvSamples2_copyRow <- 0                
            } else {
                mvSamples_copyRow  <- getsize(mvSamples)
                mvSamples2_copyRow <- getsize(mvSamples2)
            }
        }
        resize(mvSamples,  mvSamples_copyRow  + floor((niter-nburnin) / thinToUseVec[1]))
        resize(mvSamples2, mvSamples2_copyRow + floor((niter-nburnin) / thinToUseVec[2]))
        ## reinstate samplerExecutionOrder as a runtime argument, once we support non-scalar default values for runtime arguments:
        ##if(dim(samplerExecutionOrder)[1] > 0 & samplerExecutionOrder[1] == -1) {   ## runtime argument samplerExecutionOrder was not provided
        ##    lengthSamplerExecutionOrderFromConf <- dim(samplerExecutionOrderFromConfPlusTwoZeros)[1] - 2
        ##    if(lengthSamplerExecutionOrderFromConf == 0) samplerExecutionOrderToUse <- numeric(0) else samplerExecutionOrderToUse <- samplerExecutionOrderFromConfPlusTwoZeros[1:lengthSamplerExecutionOrderFromConf]
        ##} else {   ## runtime argument samplerExecutionOrder was provided
        ##    samplerExecutionOrderToUse <- samplerExecutionOrder
        ##}
        lengthSamplerExecutionOrderFromConf <- dim(samplerExecutionOrderFromConfPlusTwoZeros)[1] - 2
        if(lengthSamplerExecutionOrderFromConf == 0) samplerExecutionOrderToUse <- numeric(0) else samplerExecutionOrderToUse <- samplerExecutionOrderFromConfPlusTwoZeros[1:lengthSamplerExecutionOrderFromConf]
        if(niter < progressBarLength+3 | !progressBarDefaultSetting) progressBar <- progressBar & 0  ## cheap way to avoid compiler warning
        if(progressBar) { for(iPB1 in 1:4) { cat('|'); for(iPB2 in 1:(progressBarLength/4)) cat('-') }; print('|'); cat('|') }
        progressBarIncrement <- niter/(progressBarLength+3)
        progressBarNext <- progressBarIncrement
        progressBarNextFloor <- floor(progressBarNext)
        if(niter < 1) return()
        for(iter in 1:niter) {
            checkInterrupt()
            if(time) {
                for(i in seq_along(samplerExecutionOrderToUse)) {
                    ind <- samplerExecutionOrderToUse[i]
                    samplerTimes[ind] <<- samplerTimes[ind] + run.time(samplerFunctions[[ind]]$run())
                }
            } else {
                for(i in seq_along(samplerExecutionOrderToUse)) {
                    ind <- samplerExecutionOrderToUse[i]
                    samplerFunctions[[ind]]$run()
                }
            }
            ## adding "accumulators" to MCMC?
            ## https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
            if(iter > nburnin) {
                sampleNumber <- iter - nburnin
                if(sampleNumber %% thinToUseVec[1] == 0) {
                    mvSamples_copyRow  <- mvSamples_copyRow  + 1
                    nimCopy(from = model, to = mvSamples,  row = mvSamples_copyRow,  nodes = monitors)
                }
                if(sampleNumber %% thinToUseVec[2] == 0) {
                    mvSamples2_copyRow <- mvSamples2_copyRow + 1
                    nimCopy(from = model, to = mvSamples2, row = mvSamples2_copyRow, nodes = monitors2)
                }
            }
            if(progressBar & (iter == progressBarNextFloor)) {
                cat('-')
                progressBarNext <- progressBarNext + progressBarIncrement
                progressBarNextFloor <- floor(progressBarNext)
            }
        }
        if(progressBar) print('|')
        if((length(samplerFunctionsHMC) > 0) & nimbleVerboseOption) {
            for(i in seq_along(samplerFunctionsHMC)) {
                maxTreeDepth <- samplerFunctionsHMC[[i]]$getMaxTreeDepth()
                numDivergences <- samplerFunctionsHMC[[i]]$getNumDivergences()
                numTimesMaxTreeDepth <- samplerFunctionsHMC[[i]]$getNumTimesMaxTreeDepth()
                if(numDivergences == 1) print('HMC sampler encountered ', numDivergences, ' divergent path')
                if(numDivergences  > 1) print('HMC sampler encountered ', numDivergences, ' divergent paths')
                if(numTimesMaxTreeDepth == 1) print('HMC sampler reached the maximum search tree depth ', numTimesMaxTreeDepth, ' time')
                if(numTimesMaxTreeDepth  > 1) print('HMC sampler reached the maximum search tree depth ', numTimesMaxTreeDepth, ' times')
            }
        }
        returnType(void())
    },
    methods = list(
        getTimes = function() {
            returnType(double(1))
            return(samplerTimes[1:(length(samplerTimes)-1)])
        },
        calculateWAIC = function(nburnin = integer(default = 0),
            burnIn = integer(default = 0)) {
            if(!enableWAIC) {
                print('Error: must set enableWAIC = TRUE in buildMCMC. See help(buildMCMC) for additional information.')
                return(NaN)
            }
            if(burnIn != 0) {
                print('Warning: \'burnIn\' argument is deprecated and will not be supported in future versions of NIMBLE. Please use the \'nburnin\' argument instead.')
                ## If nburnin has not been changed, replace with burnIn value
                if(nburnin == 0)   nburnin <- burnIn
            }
            nburninPostThinning <- ceiling(nburnin/thinToUseVec[1])
            numMCMCSamples <- getsize(mvSamples) - nburninPostThinning
            if((numMCMCSamples) < 2) {
                print('Error: need more than one post burn-in MCMC samples')
                return(-Inf)
            }
            logPredProbs <- matrix(nrow = numMCMCSamples, ncol = dataNodeLength)
            logAvgProb <- 0
            pWAIC <- 0
            currentVals <- values(model, allVarsIncludingLogProbs)
            
            for(i in 1:numMCMCSamples) {
                copy(mvSamples, model, nodesTo = sampledNodes, row = i + nburninPostThinning)
                model$simulate(paramDeps)
                model$calculate(dataNodes)
                for(j in 1:dataNodeLength)
                    logPredProbs[i,j] <- model$getLogProb(dataNodes[j])
            }
            for(j in 1:dataNodeLength) {
                maxLogPred <- max(logPredProbs[,j])
                thisDataLogAvgProb <- maxLogPred + log(mean(exp(logPredProbs[,j] - maxLogPred)))
                logAvgProb <- logAvgProb + thisDataLogAvgProb
                pointLogPredVar <- var(logPredProbs[,j])
                pWAIC <- pWAIC + pointLogPredVar
            }
            WAIC <- -2*(logAvgProb - pWAIC)
            values(model, allVarsIncludingLogProbs) <<- currentVals
            if(is.nan(WAIC)) print('WAIC was calculated as NaN.  You may need to add monitors to model latent states, in order for a valid WAIC calculation.')
            returnType(double())
            return(WAIC)
        })
)

