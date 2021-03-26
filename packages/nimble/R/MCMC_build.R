
#' Create an MCMC function from a NIMBLE model, or an MCMC configuration object
#'
#' First required argument, which may be of class \code{MCMCconf} (an MCMC configuration object), or inherit from class \code{modelBaseClass} (a NIMBLE model object).  Returns an uncompiled executable MCMC function.  See details.
#'
#' @param conf An MCMC configuration object of class \code{MCMCconf} that specifies the model, samplers, monitors, and thinning intervals for the resulting MCMC function.  See \code{configureMCMC} for details of creating MCMC configuration objects.  Alternatively, \code{conf} may a NIMBLE model object, in which case an MCMC function corresponding to the default MCMC configuration for this model is returned.
#' 
#' @param conf A list of arguments for the online WAIC algorithm see the WAIC documentation below for valid inputs
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
#' \code{calculateWAIC()} accepts a single argument:
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
#' The options below are to use the online WAIC 
#' 
#' All of these should be supplied within the list \code{specsWAIC}, and if this 
#' list is non empty then \code{enableWAIC} is set to TRUE by default. 
#' 
#' The option \code{disableOnlineWAIC} will disable the online calculation 
#' of any WAIC, by default this is false, therefore if no arguments are supplied
#' and \code{enableWAIC} is activated then online WAIC calculation is performed. 
#' 
#' \code{getWAIC} takes in no arguments and is wholly tied to the MCMC process 
#' outlined earlier. 
#'
#' \code{groupingWAIC}: A list of which element is a character vector of nodes
#' to specified for that particular group. NIMBLE will then compute joint 
#' values of the PPD based on these groupings
#' 
#' \code{marginalWAIC}: A vector of nodes to integrate out when computing WAIC, 
#' in order to use a marginal likelihood.
#' 
#' \code{mWAICnIts}: Monte Carlo iterations for marginalizing (default is 1000)
#' 
#' \code{thinWAIC}: Boolean for specifying whether to do WAIC calculations
#' only on thinned samples (default is FALSE)
#' 
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
#' @author Daniel Turek & me UwU
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
    setup = function(conf, specsWAIC = NULL,  ...) {
        if(inherits(conf, 'modelBaseClass'))   conf <- configureMCMC(conf, ...)
        else if(!inherits(conf, 'MCMCconf')) stop('conf must either be a nimbleModel or a MCMCconf object (created by configureMCMC(...) )')
        dotdotdotArgs <- list(...)
        enableWAICargument <- if(!is.null(dotdotdotArgs$enableWAIC)) dotdotdotArgs$enableWAIC else nimbleOptions('MCMCenableWAIC')    ## accept enableWAIC argument regardless
        enableWAIC <- enableWAICargument || conf$enableWAIC   ## enableWAIC comes from MCMC configuration, or from argument to buildMCMC
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
        ## enable WAIC if the specsWAIC argument is non empty 
        if (!is.null(specsWAIC)){
            enableWAIC <- TRUE
        }
        ## disable online 
        onlineWAIC <- if(!is.null(specsWAIC$disableOnlineWAIC)) specsWAIC$disableOnlineWAIC else FALSE
        ## thin WAIC 
        thinWAIC <- if(!is.null(specsWAIC$thinWAIC)) specsWAIC$thinWAIC else FALSE
        ## for grouping
        dataNodes <- model$getNodeNames(dataOnly = TRUE)
        dataNodeLength <- length(dataNodes)
        if(!is.null(specsWAIC$groupingWAIC)){
            useGroupsWAIC <- TRUE
            groupNodesWAIC <- unlist(specsWAIC$groupingWAIC)
            groupNodesWAIC <- model$expandNodeNames(groupNodesWAIC)
            if (length(unique(groupNodesWAIC)) != dataNodeLength ) {
                warning("Group nodes supplied do not contain all data nodes")
            }
            groupIndicesWAIC <- sapply(specsWAIC$groupingWAIC,
                                       function(x) length(model$expandNodeNames(x)), USE.NAMES = FALSE)
            groupIndicesWAIC <- cumsum(groupIndicesWAIC)
        } else{
            useGroupsWAIC <- FALSE
            groupNodesWAIC <- dataNodes
            groupIndicesWAIC <- rep(1,dataNodeLength)
        }
        ## for mWAIC
        if(!is.null(specsWAIC$marginalWAIC)){
            mWAIC <- TRUE
            latentNodes <- model$getDependencies(specsWAIC$marginalWAIC, self = TRUE, downstream = TRUE)
            if (!is.null(specsWAIC$mWAICnIts)){
                mMCits <- specsWAIC$mWAICnIts
            }else{
                mMCits <- 1000
            }
        } else {
            mWAIC <- FALSE
            latentNodes <- dataNodes
            mMCits <- 1
        }
        dataNodeLengthWAIC <- length(groupIndicesWAIC)
        ## here to check if the group length is 1 we cannot have a length 1 vector
        if (length(groupIndicesWAIC==1)) {
            groupIndicesWAIC <- c(groupIndicesWAIC,0)
        }
        ## mWAIC variance checks
        ## A vector for checking the mWAIC is converging
        if(!is.null(specsWAIC$mWAICVarCheck)) {
            mWAICsplits <- specsWAIC$mWAICVarCheck
            varCheckItsWAIC <- seq(from = 1, to = mMCits, by = floor(mMCits/mWAICsplits) )
            lengthVarCheck <- length(varCheckItsWAIC)
            mWAICiterations <- c(varCheckItsWAIC, mMCits)
        }else if(mWAIC) {
            mWAICsplits <- 5
            varCheckItsWAIC <- floor(quantile(1:mMCits, names = FALSE))[c(2,3,4)]
            lengthVarCheck <- length(varCheckItsWAIC)
            mWAICiterations <- c(varCheckItsWAIC, mMCits)
        }else {
            mWAICsplits <- 0
            varCheckItsWAIC <- rep(0,3)
            lengthVarCheck <- 0
            mWAICiterations <- 0
        }
        ## A matrix to store the different values of mWAIC values
        mWAIClogprobmat <- matrix(0, nrow = (lengthVarCheck + 1), ncol = dataNodeLengthWAIC )
        ## the logProbVec we work with
        logProbVec <- rep(0, (dataNodeLengthWAIC + 1))
        ## lppd stores
        lppdSumMaxmat <-  matrix(0, nrow = (lengthVarCheck + 1), ncol = dataNodeLengthWAIC )
        lppdCurSummat <-  matrix(0, nrow = (lengthVarCheck + 1), ncol = dataNodeLengthWAIC )
        # pwaic stores
        sspWAICmat <- matrix(0, nrow = (lengthVarCheck + 1), ncol = dataNodeLengthWAIC )
        varCount <- 0
        meanpWAICmat <- matrix(0, nrow = (lengthVarCheck + 1), ncol = dataNodeLengthWAIC )
        delta1pWAICmat <- matrix(0, nrow = (lengthVarCheck + 1), ncol = dataNodeLengthWAIC )
        delta2pWAICmat <- matrix(0, nrow = (lengthVarCheck + 1), ncol = dataNodeLengthWAIC )
        ## mWAIC stores
        ## we put the +1 in here to fix the length 1 vector issue
        margSumMax <- rep(0, (dataNodeLengthWAIC + 1))
        margCurSum <- rep(0, (dataNodeLengthWAIC + 1))
        ## final WAIC outputs
        logAvgProb <- rep(0, (lengthVarCheck + 2))
        pWAIC <- rep(0, (lengthVarCheck + 2))
        WAIC <- rep(0, (lengthVarCheck + 2))
        ## return WAIC list
        WAICList <- nimbleList(mWAICits = double(1), WAIC = double(1), lppd = double(1), pWAIC = double(1))
        ## below is all from old WAIC method
        sampledNodes <- model$getVarNames(includeLogProb = FALSE, nodes = monitors)
        sampledNodes <- sampledNodes[sampledNodes %in% model$getVarNames(includeLogProb = FALSE)]
        paramDeps <- model$getDependencies(sampledNodes, self = FALSE, downstream = TRUE)
        allVarsIncludingLogProbs <- model$getVarNames(includeLogProb = TRUE)
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
        thin2                 = integer(default = -1)) {
        if(niter < 0)       stop('cannot specify niter < 0')
        if(nburnin < 0)     stop('cannot specify nburnin < 0')
        if(nburnin > niter) stop('cannot specify nburnin > niter')
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
                if(enableWAIC & !onlineWAIC) {
                    if (!thinWAIC) {
                        updateWAICStats()
                    } else if (sampleNumber %% thinToUseVec[1] == 0){ 
                        updateWAICStats()
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
            logAvgProbCalc <- 0
            pWAICCalc <- 0
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
                logAvgProbCalc <- logAvgProbCalc + thisDataLogAvgProb
                pointLogPredVar <- var(logPredProbs[,j])
                pWAICCalc <- pWAICCalc + pointLogPredVar
            }
            WAICCalc <- -2*(logAvgProbCalc - pWAICCalc)
            values(model, allVarsIncludingLogProbs) <<- currentVals
            if(is.nan(WAICCalc)) print('WAIC was calculated as NaN.  You may need to add monitors to model latent states, in order for a valid WAIC calculation.')
            returnType(double())
            return(WAICCalc)
        },
        updateWAICStats = function() {
            ticker <- 0
            varCount <<- varCount + 1
            for (w in 1:mMCits) {
                if(mWAIC) {
                    model$simulate(latentNodes)
                }
                model$calculate(groupNodesWAIC)
                ## get the correct logProbs
                for (j in 1:dataNodeLengthWAIC) {
                    if (!useGroupsWAIC) {
                        logProbMarg <- model$getLogProb(groupNodesWAIC[j]) 
                    } else if (j == 1) {
                        logProbMarg <-
                            model$getLogProb(groupNodesWAIC[1:groupIndicesWAIC[1]])
                    } else {
                        logProbMarg <-
                            model$getLogProb(groupNodesWAIC[(groupIndicesWAIC[(j - 1)] + 1):groupIndicesWAIC[j]])
                    }
                    ## online logSumExp for marginalization, averaging over MC samples
                    if (w == 1) {
                        margSumMax[j] <<- logProbMarg
                        margCurSum[j] <<- 1
                    } else if (logProbMarg > margSumMax[j]) {
                        newMax <- logProbMarg - margSumMax[j]
                        margSumMax[j] <<- margSumMax[j] + newMax
                        margCurSum[j] <<- margCurSum[j] * exp(-newMax) + 1
                    } else margCurSum[j] <<- margCurSum[j] + exp(logProbMarg - margSumMax[j])
                    
                }
                ## the variance of mWAIC 
                if(mWAIC & (any(w == varCheckItsWAIC))){
                    ticker <- ticker +1
                    index <- ticker
                    mWAIClogprobmat[index , ] <<- (margSumMax + log(margCurSum) - log(w))[1:dataNodeLengthWAIC]
                }
            }
            logProbVec <<- (margSumMax + log(margCurSum) - log(mMCits))[1:dataNodeLengthWAIC]
            mWAIClogprobmat[(lengthVarCheck + 1), ] <<- logProbVec
            ## update the lppd and pWAIC online 
            for (i in 1:(lengthVarCheck + 1)){
                for (j in 1:dataNodeLengthWAIC) {
                    ## online logSumExp for averaging over MCMC samples
                    #logPredProbs <- logProbVec[j]
                    logPredProbs <- mWAIClogprobmat[i,j]
                    ## lppd
                    if (varCount == 1) {
                        lppdSumMaxmat[i,j] <<- logPredProbs
                        lppdCurSummat[i,j] <<- 1
                    } else if (logPredProbs > lppdSumMaxmat[i,j]) {
                        newMax <- logPredProbs - lppdSumMaxmat[i,j]
                        lppdSumMaxmat[i,j] <<- lppdSumMaxmat[i,j] + newMax
                        lppdCurSummat[i,j] <<- lppdCurSummat[i,j] * exp(-newMax) + 1
                    } else {
                        lppdCurSummat[i,j] <<- lppdCurSummat[i,j] + exp(logPredProbs - lppdSumMaxmat[i,j])
                    }
                    ## Welford's algorithm for pWAIC
                    delta1pWAICmat[i,j] <<- logPredProbs - meanpWAICmat[i,j]
                    meanpWAICmat[i,j] <<- meanpWAICmat[i,j] + delta1pWAICmat[i,j] / varCount
                    delta2pWAICmat[i,j] <<- logPredProbs - meanpWAICmat[i,j]
                    sspWAICmat[i,j] <<- sspWAICmat[i,j] + delta1pWAICmat[i,j] * delta2pWAICmat[i,j]
                }
            }
            ## return the nodes to current proper state
            if(mWAIC) {
                ## reset the nodes if we simulated
                nimCopy(from = mvSaved, to = model, row= 1, nodes = latentNodes, logProb = TRUE)
            }
        },
        finalizeWAIC = function(j = integer(0)) {
            lppd1 <- sum(lppdSumMaxmat[j,] + log(lppdCurSummat[j,]))
            logAvgProb[j] <<- lppd1 - dataNodeLengthWAIC * log(varCount)
            pWAIC[j] <<- sum(sspWAICmat[j,] / (varCount - 1))
            WAIC[j] <<- -2 * (logAvgProb[j] - pWAIC[j])
        },
        getWAIC = function() {
            returnType(WAICList())
            if(!enableWAIC) {
                outs <- rep(NA, (lengthVarCheck + 1))
                WAICListFinal <- WAICList$new(mWAICits = outs,WAIC = outs, lppd = outs, pWAIC = outs)
                print('Error: One must set `enableWAIC = TRUE` in `configureMCMC` or `buildMCMC` in order to calculate WAIC. See `help(buildMCMC)` for additional information.')
                return(WAICListFinal)
            } else {
                for (i in 1:(lengthVarCheck + 1)){
                    finalizeWAIC(i) 
                }
                badpWAIC <- sum(((sspWAICmat[(lengthVarCheck + 1),] / (varCount -1))>0.4))
                if (badpWAIC>0) {
                    print("There are individual pWAIC values that are greater than 0.4, WAIC estimate may be unstable" )
                }
                WAIC <<- WAIC[1:(lengthVarCheck+1)]
                logAvgProb <<- logAvgProb[1:(lengthVarCheck+1)]
                pWAIC <<- pWAIC[1:(lengthVarCheck+1)]
                WAICListFinal <- WAICList$new(mWAICits = mWAICiterations,WAIC = WAIC, lppd = logAvgProb, pWAIC = pWAIC)
                return(WAICListFinal)
            }
        })
)
