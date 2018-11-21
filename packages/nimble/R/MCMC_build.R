
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
#' \code{nburnin}: Number of initial, pre-thinning, MCMC iterations to discard (default = 0).
#'
#' \code{time}: Boolean specifying whether to record runtimes of the individual internal MCMC samplers.  When \code{time = TRUE}, a vector of runtimes (measured in seconds) can be extracted from the MCMC using the method \code{mcmc$getTimes()} (default = FALSE).
#'
#' \code{progressBar}: Boolean specifying whether to display a progress bar during MCMC execution (default = TRUE).  The progress bar can be permanently disabled by setting the system option \code{nimbleOptions(MCMCprogressBar = FALSE)}.
#'
#' Samples corresponding to the \code{monitors} and \code{monitors2} from the MCMCconf are stored into the interval variables \code{mvSamples} and \code{mvSamples2}, respectively.
#' These may be accessed and converted into R matrix objects via:
#' \code{as.matrix(mcmc$mvSamples)}
#' \code{as.matrix(mcmc$mvSamples2)}
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
#' \code{nburnin}: The number of pre-thinning MCMC samples to remove from the beginning of the posterior samples for WAIC calculation (default = 0).
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
#' Equations 5, 12, and 13 in Gelman (2014) (i.e. using \emph{p}WAIC2).  The set
#' of all stochastic nodes monitored by the MCMC object will be treated as
#' \eqn{theta} for the purposes of e.g. Equation 5 from Gelman (2014). 
#' All non-monitored nodes downstream of the monitored nodes that are necessary
#' to calculate \eqn{p(y|theta)} will be simulated from the posterior samples of 
#' \eqn{theta}.  This allows customization of exactly what predictive 
#' distribution \eqn{p(y|theta)} to use for calculations.  For more detail
#' on the use of different predictive distributions, see Section 2.5 from Gelman
#' (2014).  
#' 
#' Note that there exist sets of monitored parameters that do not lead to valid
#' WAIC calculations.  Specifically, for a valid WAIC calculation, every 
#' node that a data node depends on must be either monitored, or be
#' downstream from monitored nodes.  An easy way to ensure this is satisfied
#' is to monitor all top-level parameters in a model (NIMBLE's default).  
#' Another way to guarantee correctness is to monitor all nodes
#' directly upstream from a data node. However, other combinations of monitored
#' nodes are also valid.  If \code{enableWAIC = TRUE}, NIMBLE checks to see if
#' the set of monitored nodes is valid, and returns an error if not.
#' 
#' 
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#'     y ~ dnorm(x, 1)
#' })
#' Rmodel <- nimbleModel(code, data = list(y = 0))
#' conf <- configureMCMC(Rmodel)
#' Rmcmc <- buildMCMC(conf, enableWAIC = TRUE)
#' Cmodel <- compileNimble(Rmodel)
#' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
#' Cmcmc$run(10000)
#' samples <- as.matrix(Cmcmc$mvSamples)
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
#' @export
buildMCMC <- nimbleFunction(
    name = 'MCMC',
    setup = function(conf, ...) {
    	if(inherits(conf, 'modelBaseClass'))   conf <- configureMCMC(conf, ...)
    	else if(!inherits(conf, 'MCMCconf')) stop('conf must either be a nimbleModel or a MCMCconf object (created by configureMCMC(...) )')
        dotdotdotArgs <- list(...)
        enableWAICargument <- if(!is.null(dotdotdotArgs$enableWAIC)) dotdotdotArgs$enableWAIC else nimbleOptions('enableWAIC')    ## accept enableWAIC argument regardless
        model <- conf$model
        my_initializeModel <- initializeModel(model)
        mvSaved <- modelValues(model)
        samplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(conf$samplerConfs))
            samplerFunctions[[i]] <- conf$samplerConfs[[i]]$buildSampler(model=model, mvSaved=mvSaved)
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
        ## WAIC setup:
        dataNodes <- model$getNodeNames(dataOnly = TRUE)
        dataNodeLength <- length(dataNodes)
        sampledNodes <- model$getVarNames(includeLogProb = FALSE, nodes = monitors)
        sampledNodes <- sampledNodes[sampledNodes %in% model$getVarNames(includeLogProb = FALSE)]
        paramDeps <- model$getDependencies(sampledNodes, self = FALSE, downstream = TRUE)
        enableWAIC <- enableWAICargument || conf$enableWAIC   ## enableWAIC comes from MCMC configuration, or from argument to buildMCMC
        if(enableWAIC) {
            if(dataNodeLength == 0)   stop('WAIC cannot be calculated, as no data nodes were detected in the model.')
            mcmc_checkWAICmonitors(model = model, monitors = sampledNodes, dataNodes = dataNodes)
        }
    },

    run = function(
        niter                 = integer(),
        reset                 = logical(default = TRUE),
        time                  = logical(default = FALSE),
        progressBar           = logical(default = TRUE),
        ## reinstate samplerExecutionOrder as a runtime argument, once we support non-scalar default values for runtime arguments:
        ##samplerExecutionOrder = integer(1, default = -1)
        nburnin               = integer(default =  0),
        thin                  = integer(default = -1),
        thin2                 = integer(default = -1)) {
        if(nburnin < 0)   stop('cannot specify nburnin < 0')
        thinToUseVec <<- thinFromConfVec
        if(thin  != -1)   thinToUseVec[1] <<- thin
        if(thin2 != -1)   thinToUseVec[2] <<- thin2
        my_initializeModel$run()
        nimCopy(from = model, to = mvSaved, row = 1, logProb = TRUE)
        if(reset) {
            for(i in seq_along(samplerFunctions))   samplerFunctions[[i]]$reset()
            samplerTimes <<- numeric(length(samplerFunctions) + 1)       ## default inititialization to zero
            mvSamples_copyRow  <- 0
            mvSamples2_copyRow <- 0
        } else {
            if(nburnin != 0)   stop('cannot specify nburnin when using reset = FALSE.')
            if(dim(samplerTimes)[1] != length(samplerFunctions) + 1)   samplerTimes <<- numeric(length(samplerFunctions) + 1)   ## first run: default inititialization to zero
            mvSamples_copyRow  <- getsize(mvSamples)
            mvSamples2_copyRow <- getsize(mvSamples2)
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
            currentVals2 <- values(model, sampledNodes)
            
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
            values(model, sampledNodes) <<- currentVals2
            ##model$calculate(paramDeps)
            model$calculate()
            returnType(double())
            someVals <- model$values()    ## now using model$values() syntax
            return(WAIC)
        }),
    where = getLoadingNamespace()
)

