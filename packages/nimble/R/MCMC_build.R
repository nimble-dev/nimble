#' Create an MCMC object from a NIMBLE model, or an MCMC configuration object
#'
#' First required argument, which may be of class \code{MCMCconf} (an MCMC configuration object), or inherit from class \code{modelBaseClass} (a NIMBLE model object).  Returns an uncompiled executable MCMC object.  See details.
#'
#' @param conf Either an MCMC configuration object of class \code{MCMCconf} or a NIMBLE model object. An MCMC configuration object would be returned from \code{configureMCMC} and contains information on the model, samplers, monitors, and thinning intervals to be used. Alternatively, \code{conf} may a NIMBLE model object, in which case default configuration from calling \code{configureMCMC(model, ...l)} will be used.
#'
#' @param print A logical argument, specifying whether to print details of the MCMC samplers and monitors.
#' 
#' @param ... Additional arguments to be passed to \code{configureMCMC} if \code{conf} is a NIMBLE model object (see \code{help(configureMCMC)}).
#'
#' @details
#' 
#' Calling buildMCMC(conf) will produce an uncompiled MCMC object. The object contains several methods, including the main \code{run} function for running the MCMC, a \code{getTimes} function for determining the computation time spent in each sampler (see 'getTimes' section below), and functions related to WAIC (\code{getWAIC}, \code{getWAICdetails}, \code{calculateWAIC} (see \code{help(waic)}).
#'
#' The uncompiled \code{run} function will have arguments:
#'
#' \code{niter}: The number of iterations to run the MCMC.
#'
#' \code{nburnin}: Number of initial, pre-thinning, MCMC iterations to discard (default = 0).
#'
#' \code{thin}: The thinning interval for the \code{monitors} that were specified in the MCMC configuration.  If this argument is provided at MCMC runtime, it will take precedence over the \code{thin} interval that was specified in the MCMC configuration.  If omitted, the \code{thin} interval from the MCMC configuration will be used.
#'
#' \code{thin2}: The thinning interval for the second set of monitors (\code{monitors2}) that were specified in the MCMC configuration.  If this argument is provided at MCMC runtime, it will take precedence over the \code{thin2} interval that was specified in the MCMC configuration.  If omitted, the \code{thin2} interval from the MCMC configuration will be used.
#'
#' \code{reset}: Boolean specifying whether to reset the internal MCMC sampling algorithms to their initial state (in terms of self-adapting tuning parameters), and begin recording posterior sample chains anew. Specifying \code{reset = FALSE} allows the MCMC algorithm to continue running from where it left off, appending additional posterior samples to the already existing sample chains. Generally, \code{reset = FALSE} should only be used when the MCMC has already been run (default = TRUE).
#'
#' \code{resetMV}: Boolean specifying whether to begin recording posterior sample chains anew. This argument is only considered when using \code{reset = FALSE}.  Specifying \code{reset = FALSE, resetMV = TRUE} allows the MCMC algorithm to continue running from where it left off, but without appending the new posterior samples to the already existing samples, i.e. all previously obtained samples will be erased. This option can help reduce memory usage during burn-in (default = FALSE).
#'
#' \code{resetWAIC}: Boolean specifying whether to reset the WAIC summary statistics to their initial states and thereby begin the WAIC calculation anew (default = TRUE). Specifying \code{resetWAIC = FALSE} allows the WAIC calculation to continue running from where it left off. 
#' 
#' \code{chain}: Integer specifying the MCMC chain number.  The chain number is passed to each MCMC sampler's before_chain and after_chain methods.  The value for this argument is specified automatically from invocation via runMCMC, and genernally need not be supplied when calling mcmc$run (default = 1).

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
#' The compiled object will function identically to the uncompiled object except acting on the compiled model object.
#'
#' @section Timing the MCMC samplers:
#' 
#' If you want to obtain the computation time spent in each sampler, you can set \code{time=TRUE} as a run-time argument to \code{run()} and then use the method \code{getTimes()} to obtain the times.
#' 
#' @section Calculating WAIC:
#' 
#' Please see \code{help(waic)} for more information.
#' 
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#'     y ~ dnorm(x, 1)
#' })
#' Rmodel <- nimbleModel(code, data = list(y = 0))
#' conf <- configureMCMC(Rmodel, monitors = c('mu', 'x'), enableWAIC = TRUE)
#' Rmcmc <- buildMCMC(conf)
#' Cmodel <- compileNimble(Rmodel)
#' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
#'
#' ## Running the MCMC with `run`
#' Cmcmc$run(10000)
#' samples <- as.matrix(Cmcmc$mvSamples)
#' samplesAsList <- as.list(Cmcmc$mvSamples)
#' head(samples)
#'
#' ## Getting WAIC
#' waicInfo <- Cmcmc$getWAIC()
#' waicInfo$WAIC
#' waicInfo$pWAIC
#'
#' ## Timing the samplers (must set `time = TRUE` when running the MCMC)
#' Cmcmc$run(10000, time = TRUE)
#' Cmcmc$getTimes()
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
    setup = function(conf, print, ...) {
        if(missing(print)) print <- getNimbleOption('verbose')  ## default value defined by getNimbleOption has to be here, inside the setup function code
        dotdotdotArgNames <- names(list(...))
        if(inherits(conf, 'MCMCconf') && ('enableWAIC' %in% dotdotdotArgNames || 'controlWAIC' %in% dotdotdotArgNames))
            stop("buildMCMC: 'enableWAIC' and 'controlWAIC' can only be given as arguments when running 'buildMCMC' directly on a model object, not on an MCMC configuration object. Instead pass these argument(s) directly to 'configureMCMC'.")
        
        if(inherits(conf, 'modelBaseClass')) {
            conf <- configureMCMC(conf, print = FALSE, ...)
            if(print)   conf$show(includeConfGetUnsampledNodes = FALSE)
        } else {
            if(!inherits(conf, 'MCMCconf')) stop('conf must either be a nimbleModel or a MCMCconf object (created by configureMCMC(...) )')
            if(getNimbleOption('MCMCwarnUnsampledStochasticNodes')) {
                conf$setUnsampledNodes()
                conf$warnUnsampledNodes()
            }
        }
        
        enableWAIC <- conf$enableWAIC
        model <- conf$model
        my_initializeModel <- initializeModel(model)
        mvSaved <- modelValues(model)
        
        if(getNimbleOption('MCMCorderPosteriorPredictiveSamplersLast') && length(conf$samplerConfs)) {
            ## put all posterior_predictive  samplers at the end
            samplerNames <- sapply(conf$samplerConfs, `[[`, 'name')
            postPredSamplerBool <- grepl('^posterior_predictive$', samplerNames)
            ppSamplerInd <- which(postPredSamplerBool)
            regularSamplerInd <- which(!postPredSamplerBool)
            if(length(ppSamplerInd) && length(regularSamplerInd) && (max(regularSamplerInd) > min(ppSamplerInd))) {
                messageIfVerbose('  [Note] Reordering posterior_predictive samplers to execute last')
                exOrder <- conf$samplerExecutionOrder
                if((length(exOrder)!=length(conf$samplerConfs)) || !all(exOrder==1:length(conf$samplerConfs))) stop('Halting, rather than reordering samplers in the presence of a modified sampler execution order.  If a modified execution order is needed, then: (1) reorder posterior predictive samplers to be last in the MCMC configuration printSamplers method output, (2) set the desired sampler execution order, and (3) run buildMCMC.')
                conf$samplerConfs <- conf$samplerConfs[c(regularSamplerInd, ppSamplerInd)]
            }
        }

        ## build sampler functions.
        samplerFunctions <- nimbleFunctionList(sampler_BASE)

        predictiveNodeIDs <- conf$model$getPredictiveNodeIDs()
        if(length(predictiveNodeIDs)) {
            ## save current values of 'getDependenciesIncludesPredictiveNodes' system option, then temporarily
            ## change its value to that of 'MCMCusePredictiveDependenciesInCalculations'
            getDependenciesIncludesPredictiveNodes_save <- getNimbleOption('getDependenciesIncludesPredictiveNodes')
            nimbleOptions(getDependenciesIncludesPredictiveNodes = getNimbleOption('MCMCusePredictiveDependenciesInCalculations'))
            e <- try({
                for(i in seq_along(conf$samplerConfs)) {
                    ## if option MCMCusePredictiveDependenciesInCalculations = FALSE, disallowed assignment of joint samplers to *both* PP and non-PP nodes:
                    targetScalarComponentsIsPP <- conf$model$modelDef$nodeName2GraphIDs(conf$samplerConfs[[i]]$targetAsScalar) %in% predictiveNodeIDs
                    samplingPredictiveNode <- isTRUE(any(targetScalarComponentsIsPP))
                    if(samplingPredictiveNode && !all(targetScalarComponentsIsPP) && !getNimbleOption('MCMCusePredictiveDependenciesInCalculations'))
                        stop('Cannot assign a joint sampler to simultaneously update both posterior predictive and non-posterior predictive nodes, when nimble option MCMCusePredictiveDependenciesInCalculations = FALSE')
                    ## caveat: if this sampler is sampling a predictive node, then revert the 'getDependenciesIncludesPredictiveNodes'
                    ## setting back to its original value, for creation of this sampler:
                    if(samplingPredictiveNode)   nimbleOptions(getDependenciesIncludesPredictiveNodes = TRUE)
                    samplerFunctions[[i]] <- conf$samplerConfs[[i]]$buildSampler(model=model, mvSaved=mvSaved)
                    if(samplingPredictiveNode)   nimbleOptions(getDependenciesIncludesPredictiveNodes = getNimbleOption('MCMCusePredictiveDependenciesInCalculations'))
                }},
                silent = TRUE
                )
            ## regardless whether an error occurred during sampler building, restore the original system option value:
            nimbleOptions(getDependenciesIncludesPredictiveNodes = getDependenciesIncludesPredictiveNodes_save)
            ## if an error occurred during sampler building, then quit here:
            if(inherits(e, 'try-error'))   { errorMessage <- sub('^Error.+?: ', '', e[1]); stop(errorMessage) }
        } else {
            for(i in seq_along(conf$samplerConfs))
                samplerFunctions[[i]] <- conf$samplerConfs[[i]]$buildSampler(model=model, mvSaved=mvSaved)
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
        ##nimbleVerboseOption <- getNimbleOption('verbose')   ## not currently used anywhere
        waicFun <- nimbleFunctionList(waicClass_base)
        if(enableWAIC && !('online' %in% names(conf$controlWAIC) && !conf$controlWAIC$online)) {
           waicFun[[1]] <- buildWAIC(model, mvSaved, conf$controlWAIC)
           onlineWAIC <- waicFun[[1]]$online 
           thinWAIC <- waicFun[[1]]$thin
           nburnin_extraWAIC <- waicFun[[1]]$nburnin_extra
        } else {
            if(enableWAIC) {  
                ## Setup for original (offline) WAIC prior to v. 0.12.0, namely cWAIC with no grouping capability.
                ## Retained for backward compatibility.
                ## It could also be useful to allow post hoc calculation of mWAIC given its computational
                ## cost, but that is not implemented.
                waicFun[[1]] <- buildOfflineWAIC(model, mvSamples, conf$controlWAIC, monitors)
            } else waicFun[[1]] <- buildDummyWAIC()
            onlineWAIC <- FALSE
            thinWAIC <- FALSE
            nburnin_extraWAIC <- 0
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
        nburnin               = double(default =  0),
        thin                  = double(default = -1),
        thin2                 = double(default = -1),
        resetWAIC             = logical(default = TRUE),
        chain                 = integer(default =  1)) {
        if(niter < 0)       stop('cannot specify niter < 0')
        if(nburnin < 0)     stop('cannot specify nburnin < 0')
        if(nburnin > niter) stop('cannot specify nburnin > niter')
        thinToUseVec <<- thinFromConfVec
        if(thin  != -1)   thinToUseVec[1] <<- thin
        if(thin2 != -1)   thinToUseVec[2] <<- thin2
        for(iThin in 1:2) {
            if(thinToUseVec[iThin] < 1)   stop('cannot use thin < 1')
            if(thinToUseVec[iThin] != floor(thinToUseVec[iThin]))   stop('cannot use non-integer thin')
        }
        my_initializeModel$run()
        nimCopy(from = model, to = mvSaved, row = 1, logProb = TRUE)
        if(reset) {
            samplerTimes <<- numeric(length(samplerFunctions) + 1)       ## default inititialization to zero
            for(i in seq_along(samplerFunctions))   samplerFunctions[[i]]$reset()
            for(i in seq_along(samplerFunctions))   samplerFunctions[[i]]$before_chain(niter, nburnin, chain)
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
        if(onlineWAIC & resetWAIC)
            waicFun[[1]]$reset()
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
                if(enableWAIC & onlineWAIC & iter > nburnin + nburnin_extraWAIC) {
                    if (!thinWAIC) {
                        waicFun[[1]]$updateStats()
                    } else if (sampleNumber %% thinToUseVec[1] == 0){ 
                        waicFun[[1]]$updateStats()
                    }
                }
            }
            if(progressBar & (iter == progressBarNextFloor)) {
                cat('-')
                progressBarNext <- progressBarNext + progressBarIncrement
                progressBarNextFloor <- floor(progressBarNext)
            }
        }
        if(progressBar) print('|')
        for(i in seq_along(samplerFunctions))   samplerFunctions[[i]]$after_chain()
        returnType(void())
    },
    methods = list(
        getTimes = function() {
            returnType(double(1))
            return(samplerTimes[1:(length(samplerTimes)-1)])
        },
        ## Old-style post-sampling WAIC calculation.
        calculateWAIC = function(nburnin = integer(default = 0)) {
            if(!enableWAIC) {
                print('Error: One must set enableWAIC = TRUE in \'configureMCMC\' or \'buildMCMC\'. See \'help(configureMCMC)\' for additional information.')
                return(NaN)
            }
            result <- waicFun[[1]]$calculateWAIC(nburnin, thinToUseVec[1])
            returnType(double())
            return(result$WAIC)
        },
        getWAIC = function() {
            returnType(waicNimbleList())
            if(enableWAIC) {
                return(waicFun[[1]]$get())
            } else {
                print("WAIC was disabled based on the 'enableWAIC = FALSE'. You may be able to use the \'calculateWAIC\' function.")
                return(waicNimbleList$new(WAIC = NA, lppd = NA, pWAIC = NA))
            }
        },
        getWAICdetails = function(returnElements = logical(default = FALSE)) {
            returnType(waicDetailsNimbleList())
            if(enableWAIC & onlineWAIC) {
                return(waicFun[[1]]$getDetails(returnElements))
            } else {
                print("WAIC details are only available when using online WAIC. Online WAIC was disabled based on the 'onlineWAIC' element of WAIC control list.")
                return(waicDetailsNimbleList$new(marginal = FALSE, niterMarginal = 0, thin = FALSE, online = FALSE, nburnin_extra = 0))
            }
        }
    )
        
)
