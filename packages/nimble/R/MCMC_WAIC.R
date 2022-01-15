waicClass_base <- nimbleFunctionVirtual(
    name = 'waicClass_base',
    methods = list(
        reset = function() {},
        updateStats = function() {},
        get = function() { returnType(waicList()) },
        getDetails = function(returnElements = logical(default = FALSE)) { returnType(waicDetailsList()) },
        calculateWAIC = function(nburnin = integer(default = 0), thin = double(default = 1)) { returnType(waicList()) }
    )
)

buildDummyWAIC <- nimbleFunction(
    name = 'waicClass_dummy',
    contains = waicClass_base,
    setup = function() {},
    run = function() {},
    methods = list(
        reset = function() {},
        updateStats = function() {},
        get = function() { return(waicList$new(WAIC = NA, lppd = NA, pWAIC = NA)); returnType(waicList()) },
        getDetails = function(returnElements = logical(default = FALSE)) {return(waicDetailsList$new(marginal = FALSE, niterMarginal = 0, thin = FALSE, online = FALSE));  returnType(waicDetailsList()) },
        calculateWAIC = function(nburnin = integer(default = 0), thin = double(default = 1)) { return(waicList$new(WAIC = NA, lppd = NA, pWAIC = NA)); returnType(waicList()) }
    )
)

buildOfflineWAIC <- nimbleFunction(
    name = 'waicClass_offline',
    contains = waicClass_base,
    setup = function(model, mvSamples, control, monitors) {
        dataNodes <- model$getNodeNames(dataOnly = TRUE)
        dataNodeLength <- length(dataNodes)
        sampledNodes <- model$getVarNames(includeLogProb = FALSE, nodes = monitors)
        sampledNodes <- sampledNodes[sampledNodes %in% model$getVarNames(includeLogProb = FALSE)]
        paramDeps <- model$getDependencies(sampledNodes, self = FALSE, downstream = TRUE)
        allVarsIncludingLogProbs <- model$getVarNames(includeLogProb = TRUE)
        if(!dataNodeLength)
            stop('Offline WAIC cannot be calculated, as no data nodes were detected in the model.')
        mcmc_checkWAICmonitors_conditional(model = model, monitors = sampledNodes, dataNodes = dataNodes)
        lppd <- -Inf
        pWAIC <- 0
        WAIC <- -Inf
        finalized <- FALSE
    },
    run = function() {},
    methods = list(
        reset = function() {},
        updateStats = function() {},
        getDetails = function(returnElements = logical(default = FALSE)) {return(waicDetailsList$new(marginal = FALSE, niterMarginal = 0, thin = FALSE, online = FALSE));  returnType(waicDetailsList()) },
        calculateWAIC = function(nburnin = integer(default = 0), thin = double(default = 1)) {
            nburninPostThinning <- ceiling(nburnin/thin)
            numMCMCSamples <- getsize(mvSamples) - nburninPostThinning
            if((numMCMCSamples) < 2) {
                print('Error: need more than one post burn-in MCMC samples')
            }
            logPredProbs <- matrix(nrow = numMCMCSamples, ncol = dataNodeLength)
            logAvgProb <- 0
            pWAIC <<- 0
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
                pWAIC <<- pWAIC + pointLogPredVar
            }
            lppd <<- logAvgProb
            WAIC <<- -2*(logAvgProb - pWAIC)
            
            values(model, allVarsIncludingLogProbs) <<- currentVals
            if(is.nan(WAIC)) print('WAIC was calculated as NaN.  You may need to add monitors to model latent states, in order for a valid WAIC calculation.')
            finalized <<- TRUE

            returnType(waicList())
            return(get())
        },
        get = function() {
            ## Extract WAIC summary information.
            returnType(waicList())
            ## Ideally we would call finalize if needed, but to mimic old offline behavior,
            ## we need to pass in nburnin, and finalize method for online WAIC doesn't
            ## take arguments.
            if(!finalized)
                stop("Please execute the 'calculateWAIC' method of the offline WAIC object before running 'get'.")
            output <- waicList$new()

            output$WAIC <- WAIC
            output$lppd <- lppd
            output$pWAIC  <- pWAIC

            return(output)
        }
    )
)

## waicList and waicDetailsList definitions are in nimbleList_core.R.

buildWAIC <- nimbleFunction(
    name = 'waicClass',
    contains = waicClass_base, 
    setup = function(model, mvSaved, control) {
        online              <- extractControlElement(control, 'online',              TRUE)
        thin                <- extractControlElement(control, 'thin',                FALSE)
        dataGroups          <- extractControlElement(control, 'dataGroups',          NULL)
        marginalizeNodes    <- extractControlElement(control, 'marginalizeNodes',    NULL)
        niterMarginal       <- extractControlElement(control, 'niterMarginal',       1000)
        convergenceSet      <- extractControlElement(control, 'convergenceSet',      NULL)

        setupOutputs(online, thin)
        
        ## Partition of data nodes. By default one data node per group.
        dataNodes <- model$getNodeNames(dataOnly = TRUE)
        dataNodeLength <- length(dataNodes)
        if(!dataNodeLength)
            stop('buildWAIC: cannot compute WAIC with no data nodes in the model.')
        if(!is.null(dataGroups)) { 
            useGroups <- TRUE
            dataNodesUser <- lapply(dataGroups, function(x) model$expandNodeNames(x))
            groupIndices <- lapply(dataNodesUser, length)
            groupIndices <- cumsum(groupIndices)
            dataNodesUser <- unlist(dataNodesUser)
            if (length(dataNodesUser) != dataNodeLength || sort(dataNodesUser) != sort(dataNodes)) 
                warning("buildWAIC: Potential problem with data grouping. The nodes included in 'dataGroups' do not match the full set of data nodes in the model.")
            dataNodes <- dataNodesUser
        } else{
            useGroups <- FALSE
            groupIndices <- rep(1, dataNodeLength)
        }
        nGroups <- length(groupIndices)
        ## Cannot have a length 1 vector.
        if(nGroups == 1) 
            groupIndices <- c(groupIndices, 0)
        
        ## Determine marginal versus conditional WAIC.
        if(!is.null(marginalizeNodes)) {
            marginal <- TRUE
            otherLatentNodes <- model$getDependencies(marginalizeNodes, self = FALSE, downstream = TRUE, stochOnly = TRUE, includeData = FALSE)
            marginalizeNodes <- model$getDependencies(marginalizeNodes, self = TRUE, downstream = TRUE, includeData = FALSE)
            if(length(otherLatentNodes))
                warning("buildWAIC: Additional latent nodes found that were not specified in 'marginalizeNodes': ", paste0(otherLatentNodes, collapse = ', '), ". These nodes are also being marginalized over.")
            latentAndDataNodes <- model$getDependencies(marginalizeNodes, self = TRUE, downstream = TRUE)
            if(any(!marginalizeNodes %in% model$getNodeNames(latentOnly = TRUE)))
                warning("buildWAIC: Potential problem with nodes to marginalize over. One or more of the nodes in 'marginalizeNodes' are not latent nodes in the model.")
        } else {
            marginal  <- FALSE 
            niterMarginal <- 1
            ## Next two are unused in conditional case but need to exist for compilation.
            marginalizeNodes <- dataNodes[1]
            latentAndDataNodes <- dataNodes[1]

        }
        
        
        ## For mWAIC, we save relevant quantities for subsets of the full Monte
        ## Carlo sample to allow informal checking that the Monte Carlo estimates
        ## are stable.
        ## By default, we use 25%, 50%, and 75% of the total number of Monte Carlo samples.
        if(marginal) {
            if(!is.null(convergenceSet)) {
                convergenceSet <- convergenceSet[!convergenceSet %in% c(0, 1)]
                if(!length(convergenceSet) || min(convergenceSet) <= 0 || max(convergenceSet) >= 1)
                    stop("buildWAIC: 'convergenceSet' for assessing Monte Carlo error in marginal WAIC should be a set of values between 0 and 1.")
                checkIts <- convergenceSet * niterMarginal
                lengthConvCheck <- length(checkIts) + 1
            } else {  # Default to 25%, 50%, 75%.
                checkIts <- floor(quantile(1:niterMarginal, names = FALSE))[c(2,3,4)]
                lengthConvCheck <- length(checkIts) + 1
            }
            if(length(checkIts) == 1) checkIts <- c(checkIts, 0)
        } else {
            checkIts <- rep(0, 2)
            lengthConvCheck <- 1
        }
        
        logProbMat <- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
        ## lppd stores
        lppdSumMaxMat <-  matrix(0, nrow = lengthConvCheck, ncol = nGroups)
        lppdCurSumMat <-  matrix(0, nrow = lengthConvCheck, ncol = nGroups)
        # pWAIC stores
        sspWAICmat <- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
        meanpWAICmat <- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
        delta1pWAICmat <- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
        delta2pWAICmat <- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
        ## mWAIC stores
        len <- nGroups; if(len == 1) len <- 2  # avoid length one vectors
        margSumMax <- rep(0, len)
        margCurSum <- rep(0, len)
        ## final WAIC outputs
        len <- lengthConvCheck; if(len == 1) len <- 2
        lppd <- rep(0, len)
        pWAIC <- rep(0, len)
        WAIC <- rep(0, len)

        mcmcIter <- 0
        finalized <- FALSE
    },
    run = function() {},
    methods = list(
        updateStats = function() {
            ## Online updating of summary stats, called once per MCMC iteration that is used.
            ticker <- 0  # indexes over MC subsets
            mcmcIter <<- mcmcIter + 1
            for (k in 1:niterMarginal) {  # loop over MC samples; for conditional, there is only one 'iteration'
                if(marginal) {
                    model$simulate(marginalizeNodes)
                    model$calculate(dataNodes) 
                }
                ## Extract logProbs for each group.
                for (j in 1:nGroups) {
                    if (!useGroups) {
                        logProbMarg <- model$getLogProb(dataNodes[j]) 
                    } else if (j == 1) {
                        logProbMarg <-
                            model$getLogProb(dataNodes[1:groupIndices[1]])
                    } else {
                        logProbMarg <-
                            model$getLogProb(dataNodes[(groupIndices[(j - 1)] + 1):groupIndices[j]])
                    }
                    ## online logSumExp for marginalization, averaging over MC samples
                    if (k == 1) {
                        margSumMax[j] <<- logProbMarg
                        margCurSum[j] <<- 1
                    } else if (logProbMarg > margSumMax[j]) {
                        margCurSum[j] <<- margCurSum[j] * exp(-logProbMarg + margSumMax[j]) + 1
                        margSumMax[j] <<- logProbMarg
                    } else margCurSum[j] <<- margCurSum[j] + exp(logProbMarg - margSumMax[j])
                }
                ## Save quantities for each element of convergence check set.
                if(marginal & (any(k == checkIts))){
                    ticker <- ticker + 1
                    logProbMat[ticker, ] <<- (margSumMax + log(margCurSum) - log(k))[1:nGroups]
                }
            }
            logProbMat[lengthConvCheck, ] <<- (margSumMax + log(margCurSum) - log(niterMarginal))[1:nGroups]
            ## update the lppd and pWAIC in an online manner. 
            for (i in 1:lengthConvCheck){
                for (j in 1:nGroups) {
                    ## online logSumExp for averaging over MCMC samples
                    logPredProb <- logProbMat[i, j]
                    ## lppd
                    if (mcmcIter == 1) {
                        lppdSumMaxMat[i, j] <<- logPredProb
                        lppdCurSumMat[i, j] <<- 1
                    } else if (logPredProb > lppdSumMaxMat[i, j]) {
                        lppdCurSumMat[i, j] <<- lppdCurSumMat[i, j] * exp(-logPredProb + lppdSumMaxMat[i, j]) + 1
                        lppdSumMaxMat[i, j] <<- logPredProb
                    } else {
                        lppdCurSumMat[i, j] <<- lppdCurSumMat[i, j] + exp(logPredProb - lppdSumMaxMat[i, j])
                    }
                    ## Welford's algorithm for pWAIC
                    delta1pWAICmat[i, j] <<- logPredProb - meanpWAICmat[i, j]
                    meanpWAICmat[i, j] <<- meanpWAICmat[i, j] + delta1pWAICmat[i, j] / mcmcIter
                    delta2pWAICmat[i, j] <<- logPredProb - meanpWAICmat[i, j]
                    sspWAICmat[i, j] <<- sspWAICmat[i, j] + delta1pWAICmat[i, j] * delta2pWAICmat[i, j]
                }
            }
            ## Return the nodes and logProbs to proper original state, as mWAIC simulates into them.
            if(marginal) {
                nimCopy(from = mvSaved, to = model, row= 1, nodes = latentAndDataNodes, logProb = TRUE)
            }
        },
        finalize = function() {
            ## Calculate WAIC quantities after done updating.
            for(i in 1:lengthConvCheck) {
                tmp <- sum(lppdSumMaxMat[i, ] + log(lppdCurSumMat[i, ]))
                lppd[i] <<- tmp - nGroups * log(mcmcIter)
                pWAIC[i] <<- sum(sspWAICmat[i, ]) / (mcmcIter - 1)
                WAIC[i] <<- -2 * (lppd[i] - pWAIC[i])
            }
            finalized <<- TRUE
        },
        get = function() {
            ## Extract WAIC summary information.
            returnType(waicList())
            if(!finalized)
                finalize()
            badpWAIC <- sum( sspWAICmat[lengthConvCheck, ] / (mcmcIter-1) > 0.4 )
            if(badpWAIC) {  
                print("There are individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes." )
            }
            
            output <- waicList$new()
            
            output$WAIC <- WAIC[lengthConvCheck]
            output$lppd <- lppd[lengthConvCheck]
            output$pWAIC  <- pWAIC[lengthConvCheck]

            return(output)
        },
        getDetails = function(returnElements = logical(default = FALSE)) {
            ## Extract WAIC detailed information.
            returnType(waicDetailsList())
            if(!finalized)
                finalize()

            output <- waicDetailsList$new()
            
            output$marginal <- marginal
            output$thin <- thin
            output$online <- online

            if(marginal) {
                output$niterMarginal <- niterMarginal
                output$WAIC_partialMC <- WAIC[1:(lengthConvCheck - 1)]
                output$lppd_partialMC <- lppd[1:(lengthConvCheck - 1)]
                output$pWAIC_partialMC <- pWAIC[1:(lengthConvCheck - 1)]
                output$niterMarginal_partialMC <- checkIts
            } else output$niterMarginal <- 0
            
            ## Per-observations values useful for (manual) SE calculation when
            ## contrasting WAIC for two models, per Vehtari et al. 2017 (Staistics and Computing)
            ## This would be long if many obs, so have under user control.
            if(returnElements) {  
                output$lppd_elements <- lppdSumMaxMat[lengthConvCheck, ] +
                    log(lppdCurSumMat[lengthConvCheck, ]) - log(mcmcIter)
                output$pWAIC_elements <- sspWAICmat[lengthConvCheck, ] / (mcmcIter - 1)
                output$WAIC_elements <- -2 * (output$lppd_elements - output$pWAIC_elements)
            }
            return(output)
        },
        reset = function() {
            ## Reset values that are accumulated in updateStats().
            mcmcIter <<- 0
            lppdSumMaxMat <<-  matrix(0, nrow = lengthConvCheck, ncol = nGroups)
            lppdCurSumMat <<-  matrix(0, nrow = lengthConvCheck, ncol = nGroups)
            sspWAICmat <<- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
            meanpWAICmat <<- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
            delta1pWAICmat <<- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
            delta2pWAICmat <<- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
            logProbMat <<- matrix(0, nrow = lengthConvCheck, ncol = nGroups)
            finalized <<- FALSE
        },
        calculateWAIC = function(nburnin = integer(default = 0), thin = double(default = 1)) { return(waicList$new(WAIC = NA, lppd = NA, pWAIC = NA)); returnType(waicList()) }
    )
)


#' Calculating WAIC using an offline algorithm
#'
#' In addition to the core online algorithm, NIMBLE implements an offline
#' WAIC algorithm that can be computed on the results of an MCMC. In contrast
#' to NIMBLE's built-in online WAIC, offline WAIC can compute only conditional
#' WAIC and does not allow for grouping data nodes.
#' 
#' @param mcmc An MCMC object (compiled or uncompiled) or matrix or dataframe
#' of MCMC samples as the first argument of \code{calculateWAIC}.
#'
#' @param model A model (compiled or uncompiled) as the second argument of
#' \code{calculateWAIC}. Only required if \code{mcmc} is a matrix/dataframe
#' of samples.
#'
#' @param nburnin The number of pre-thinning MCMC samples to remove from the beginning
#' of the posterior samples for offline WAIC calculation via \code{calculateWAIC}
#' (default = 0). These samples are discarded in addition to any burn-in specified when
#' running the MCMC.
#'
#' @param thin Thinning factor interval to apply to the samples for offline
#' WAIC calculation using \code{calculateWAIC} (default = 1,
#' corresponding to no thinning).
#'
#' @details
#'
#' The ability to calculate WAIC post hoc after all MCMC sampling has been done
#' has certain advantages (e.g., allowing a user to exclude additional burnin
#' samples beyond that specified initially for the MCMC) in addition to
#' providing compatibility with versions of NIMBLE before 0.12.0. This
#' functionality includes the ability to call the \code{calculateWAIC} function
#' on an MCMC object or matrix of samples after running an MCMC and without
#' setting up the MCMC initially to use WAIC, provided the necessary variables
#' are monitored.
#'
#' See \code{help(waic)} for details on using NIMBLE's recommended online
#' algorithm for WAIC.
#' 
#' @section Offline WAIC (WAIC computed after MCMC sampling):
#'
#' As an alternative to online WAIC, NIMBLE also provides a function,
#' \code{calculateWAIC}, that can be called on an MCMC object or a matrix of
#' samples, after running an MCMC. This function does not require that one
#' set \code{enableWAIC = TRUE} nor \code{WAIC = TRUE} when calling
#' \code{runMCMC}. The function checks that the necessary variables were
#' monitored in the MCMC and returns an error if they were not. This function
#' behaves identically to the \code{calculateWAIC} method of an MCMC object.
#' Note that to use this function when using \code{nimbleMCMC} one would
#' need to build the model outside of \code{nimbleMCMC}.
#'
#' The \code{calculateWAIC} function requires either an MCMC object or a matrix
#' (or dataframe) of posterior samples plus a model object. In addition, one
#' can provide optional \code{burnin} and \code{thin} arguments.
#' 
#' In addition, for compatibility with older versions of NIMBLE (prior to
#' v0.12.0), one can also use the \code{calculateWAIC} method of the MCMC
#' object to calculate WAIC after all sampling has been completed.
#'
#' The \code{calculateWAIC()} method accepts a single argument, \code{nburnin},
#' equivalent to the \code{nburnin} argument of the \code{calculateWAIC}
#' function described above.
#' 
#' The \code{calculateWAIC} method can only be used if the \code{enableWAIC} 
#' argument to \code{configureMCMC} or to \code{buildMCMC} is set to \code{TRUE},
#' or if the NIMBLE option \code{enableWAIC} is set to \code{TRUE}.  If a user
#' attempts to call \code{calculateWAIC} without having set
#' \code{enableWAIC = TRUE} (either in the call to \code{configureMCMC}, or
#' \code{buildMCMC}, or as a NIMBLE option), an error will occur.  
#'
#' The \code{calculateWAIC} function and method calculate the WAIC based on
#' Equations 5, 12, and 13 in Gelman et al. (2014) (i.e., using \emph{p}WAIC2).
#'
#' Note that there is not a unique value of WAIC for a model. The 
#' \code{calculateWAIC} function and method only provide the conditional WAIC,
#' namely the version of WAIC where all parameters directly involved in the
#' likelihood are treated as \eqn{theta} for the purposes of Equation 5 from
#' Gelman et al. (2014). As a result, the user must set the MCMC monitors
#' (via the \code{monitors} argument) to include all stochastic nodes that
#' are parents of any data nodes; by default the MCMC monitors are only the
#' top-level nodes of the model. For more detail on the use of different
#' predictive distributions, see Section 2.5 from Gelman et al. (2014) or
#' Ariyo et al. (2019).

#' Also note that WAIC relies on a partition of the observations, i.e.,
#' 'pointwise' prediction. In \code{calculateWAIC} the sum over log pointwise
#' predictive density values treats each data node as contributing a single
#' value to the sum. When a data node is multivariate, that data node contributes
#' a single value to the sum based on the joint density of the elements in the
#' node. Note that if one wants the WAIC calculation via \code{calculateWAIC}
#' to be based on the joint predictive density for each group of observations
#' (e.g., grouping the observations from each person or unit in a longitudinal
#' data context), one would need to use a multivariate distribution for the
#' observations in each group (potentially by writing a user-defined
#' distribution).
#'
#' For more control over and flexibility in how WAIC is calculated, see
#' \code{help(waic)}.
#'
#' @author Joshua Hug and Christopher Paciorek
#' 
#' @seealso \code{\link{waic}} \code{\link{configureMCMC}}
#' \code{\link{buildMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
#'
#' @references 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' \emph{Journal of Machine Learning Research} 11: 3571-3594.
#' 
#' Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models.
#' \emph{Statistics and Computing} 24(6): 997-1016.
#'
#' Ariyo, O., Quintero, A., Munoz, J., Verbeke, G. and Lesaffre, E. (2019).
#' Bayesian model selection in linear mixed models for longitudinal data.
#' \emph{Journal of Applied Statistics} 47: 890-913.
#'
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing} 27: 1413-1432.
#'
#' Hug, J.E.  and Paciorek, C.J. (2021). A numerically stable online
#' implementation and exploration of WAIC through variations of the
#' predictive density, using NIMBLE. \emph{arXiv e-print} <arXiv:2106.13359>.
#'
#' @export
#'
#' @examples
#' code <- nimbleCode({
#'   for(j in 1:J) {
#'     for(i in 1:n) 
#'       y[j, i] ~ dnorm(mu[j], sd = sigma)
#'     mu[j] ~ dnorm(mu0, sd = tau)
#'   }
#'   tau ~ dunif(0, 10)
#'   sigma ~ dunif(0, 10)
#' })
#' J <- 5
#' n <- 10
#' y <- matrix(rnorm(J*n), J, n)
#' Rmodel <- nimbleModel(code, constants = list(J = J, n = n), data = list(y = y),
#'                       inits = list(tau = 1, sigma = 1))
#'
#' ## Make sure the needed variables are monitored.
#' ## Only conditional WAIC without data grouping is available via this approach.
#' conf <- configureMCMC(Rmodel, monitors = c('mu', 'sigma'))
#' \dontrun{
#' Cmodel <- compileNimble(Rmodel)
#' Rmcmc <- buildMCMC(conf)
#' Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#' output <- runMCMC(Cmcmc, niter = 1000)
#' calculateWAIC(Cmcmc)           # Can run on the MCMC object
#' calculateWAIC(output, Rmodel)  # Can run on the samples directly
#'
#' ## Apply additional burnin (additional to any burnin already done in the MCMC.
#' calculateWAIC(Cmcmc, burnin = 500)
#' }
#' @export
calculateWAIC <- function(mcmc, model, nburnin = 0, thin = 1) {
    if((is(mcmc, 'MCMC') || is(mcmc, 'MCMC_refClass')) &&
       identical(nfGetDefVar(mcmc, 'name'), 'MCMC')) {
        ## MCMC is provided
        if(exists('model', mcmc, inherits = FALSE)) compiled <- FALSE else compiled <- TRUE
        if(compiled) {
            if(!exists('Robject', mcmc, inherits = FALSE) || !exists('model', mcmc$Robject, inherits = FALSE))
                stop("calculateWAIC: problem with finding model object in compiled MCMC")
            model <- mcmc$Robject$model
            mvSamples <- mcmc$Robject$mvSamples
        } else {
            model <- mcmc$model
            mvSamples <- mcmc$mvSamples
        }
        samples <- as.matrix(mcmc$mvSamples) # because compilation below is not giving access to existing samples
        usingMCMC <- TRUE
    } else {
        ## Samples matrix is provided
        if(!(is.matrix(mcmc) || is.data.frame(mcmc)))
            stop("calculateWAIC: 'mcmc' must be a matrix/dataframe of samples or a nimble MCMC object.")
        if(is.data.frame(mcmc))
            mcmc <- as.matrix(mcmc)
        if(missing(model)) 
            stop("calculateWAIC: 'model' is required if samples are directly provided rather than an MCMC object.")
        ## Get uncompiled model to avoid warning when passing model into compileNimble
        if(exists('Rmodel', model, inherits = FALSE))
            model <- model$Rmodel

        nm <- model$getVarNames(nodes = colnames(mcmc))
        modelSymbolObjects = model$getSymbolTable()$getSymbolObjects()
        if(!all(nm %in% names(modelSymbolObjects)))
            stop('calculateWAIC: some column names indicate variables that are not in the model.') 
        mvSamplesConf <- modelValuesConf(symbolTable(symbols = modelSymbolObjects[nm]))
        mvSamples <- modelValues(mvSamplesConf)

        usingMCMC <- FALSE
    } 
    waicFun <- buildOfflineWAIC(model, mvSamples, NULL, mvSamples$varNames)
    cwaicFun <- compileNimble(waicFun, project = model, resetFunctions = TRUE)
    if(!usingMCMC) { # copy to compiled mvSamples for speed 
        matrix2mv(mcmc, cwaicFun$mvSamples)
    } else matrix2mv(samples, cwaicFun$mvSamples) # for the moment we can't use the existing mvSamples values directly
    result <- cwaicFun$calculateWAIC(nburnin, thin)
    return(result)
}

#' Using WAIC
#'
#' Details of the WAIC measure for comparing models. NIMBLE implements an online
#' WAIC algorithm, computed during the course of the MCMC iterations.
#' 
#' @name waic
#' 
#' @aliases getWAIC getWAICdetails buildWAIC WAIC enableWAIC
#' 
#' @details
#'
#' To use WAIC, set \code{enableWAIC = TRUE} when configuring or (if not using
#' \code{configureMCMC} building an MCMC) and set \code{WAIC = TRUE} when
#' calling \code{nimbleMCMC} and optionally when calling \code{runMCMC}.
#'
#' By default, NIMBLE calculates WAIC using an online algorithm that updates
#' required summary statistics at each post-burnin iteration of the MCMC.
#'
#' One can also use \code{calculateWAIC} to run an offline version of the
#' WAIC algorithm after all MCMC sampling has been done. This allows calculation
#' of WAIC from a matrix (or dataframe) of posterior samples and also retains
#' compatibility with WAIC in versions of NIMBLE before 0.12.0. However, the
#' offline algorithm is less flexible than the online algorithm and only
#' provides conditional WAIC without the ability to group data points. See
#' \code{help(calculateWAIC)} for details.
#'
#' @section \code{controlWAIC} list:
#'
#' The \code{controlWAIC} argument is a list that controls the behavior of the
#' WAIC algorithm and is passed to either \code{configureMCMC} or (if not using
#' \code{configureMCMC}) \code{buildMCMC}. One can supply any of the following
#' optional components:
#'
#' \code{online}: Logical value indicating whether to calculate WAIC during the
#' course of the MCMC. Default is \code{TRUE} and setting to \code{FALSE} is
#' primarily for backwards compatibility to allow use of the old
#' \code{calculateWAIC} method that calculates WAIC from monitored values after
#' the MCMC finishes.
#'
#' \code{dataGroups}: Optional list specifying grouping of data nodes,
#' one element per group, with each list element containing the node names
#' for the data nodes in that group. If provided, the predictive density values
#' computed will be the joint density values, one joint density per group.
#' Defaults to one data node per 'group'. See details.
#'
#' \code{marginalizeNodes}: Optional set of nodes (presumably latent nodes)
#' over which to marginalize to compute marginal WAIC (i.e., WAIC based on a
#' marginal likelihood), rather than the default conditional WAIC (i.e., WAIC
#' conditioning on all parent nodes of the data nodes). See details.
#'
#' \code{niterMarginal}: Number of Monte Carlo iterations to use when
#' marginalizing (default is 1000).
#'
#' \code{convergenceSet}: Optional vector of numbers between 0 and 1 that
#' specify a set of shorter Monte Carlo simulations for marginal WAIC
#' calculation as fractions of the full (\code{niterMarginal}) Monte Carlo
#' simulation. If not provided, NIMBLE will use 0.25, 0.50, and 0.75.
#' NIMBLE will report the WAIC, lppd, and pWAIC that would have been obtained
#' for these smaller Monte Carlo simulations, allowing assessment of the number
#' of Monte Carlo samples needed for stable calculation of WAIC.
#' 
#' \code{thin}: Logical value for specifying whether to do WAIC calculations
#' only on thinned samples (default is \code{FALSE}). Likely only useful for
#' reducing computation when using marginal WAIC.
#'
#' @section Extracting WAIC:
#'
#' The calculated WAIC and related quantities can be obtained in various ways
#' depending on how the MCMC is run. If using \code{nimbleMCMC} and setting
#' \code{WAIC = TRUE}, see the \code{WAIC} component of the output list. If using
#' \code{runMCMC} and setting \code{WAIC = TRUE}, either see the \code{WAIC}
#' component of the output list or use the \code{getWAIC} method of the MCMC
#' object (in the latter case \code{WAIC = TRUE} is not required). If using
#' the \code{run} method of the MCMC object, use the \code{getWAIC} method of
#' the MCMC object.
#'
#' The output of running WAIC (unless one sets \code{online = FALSE}) is a list
#' containing the following components:
#'
#' \code{WAIC}: The computed WAIC, on the deviance scale. Smaller values are
#' better when comparing WAIC for two models.
#' 
#' \code{lppd}: The log predictive density component of WAIC.
#' 
#' \code{pWAIC}: The pWAIC estimate of the effective number of parameters,
#' computed using the \emph{p}WAIC2 method of Gelman et al. (2014).
#'
#' To get further information, one can use the \code{getWAICdetails} method
#' of the MCMC object.  The result of running \code{getWAICdetails} is a list
#' containing the following components:
#' 
#' \code{marginal}: Logical value indicating whether marginal (\code{TRUE}) or
#' conditional (\code{FALSE}) WAIC was calculated.
#' 
#' \code{niterMarginal}: Number of Monte Carlo iterations used in computing
#' marginal likelihoods if using marginal WAIC.
#' 
#' \code{thin}: Whether WAIC was calculated based only on thinned samples.
#' 
#' \code{online}: Whether WAIC was calculated during MCMC sampling.
#' 
#' \code{WAIC_partialMC}, \code{lppd_partialMC}, \code{pWAIC_partialMC}: The
#' computed marginal WAIC, lppd, and pWAIC based on fewer Monte Carlo
#' simulations, for use in assessing the sensitivity of the WAIC calculation
#' to the number of Monte Carlo iterations.
#'
#' \code{niterMarginal_partialMC}: Number of Monte Carlo iterations used for the
#' values in \code{WAIC_partialMC}, \code{lppd_partialMC}, \code{pWAIC_partialMC}.
#' 
#' \code{WAIC_elements}, \code{lppd_elements}, \code{pWAIC_elements}: Vectors of
#' individual WAIC, lppd, and pWAIC values, one element per data node (or group
#' of nodes in the case of specifying \code{dataGroups}). Of use in computing
#' the standard error of the difference in WAIC between two models, following
#' Vehtari et al. (2017).
#'
#' @section Online WAIC:
#'
#' As of version 0.12.0, NIMBLE provides enhanced WAIC functionality, with user
#' control over whether to use conditional or marginal versions of WAIC and
#' whether to group data nodes. In addition, users are no longer required to
#' carefully choose MCMC monitors. WAIC by default is now calculated in an online
#' manner (updating the required summary statistics at each MCMC iteration),
#' using all post-burnin samples. The WAIC (Watanabe, 2010) is calculated from
#' Equations 5, 12, and 13 in Gelman et al. (2014) (i.e., using 'pWAIC2').
#'
#' Note that there is not a unique value of WAIC for a model. By default, WAIC
#' is calculated conditional on the parent nodes of the data nodes, and the
#' density values used are the individual density values of the data nodes.
#' However, by modifying the \code{marginalizeNodes} and \code{dataGroups}
#' elements of the control list, users can request a marginal WAIC (using a
#' marginal likelihood that integrates over user-specified latent nodes) and/or
#' a WAIC based on grouping observations (e.g., all observations in a cluster)
#' to use joint density values. See the MCMC Chapter of the NIMBLE
#' \href{https://r-nimble.org/html_manual/cha-mcmc.html}{User Manual}
#' for more details.
#'
#' For more detail on the use of different predictive distributions, see Section
#' 2.5 from Gelman et al. (2014) or Ariyo et al. (2019).
#'
#' Note that based on a limited set of simulation experiments in Hug and Paciorek
#' (2021) our tentative recommendation is that users only use marginal WAIC if
#' also using grouping. 
#' 
#' @author Joshua Hug and Christopher Paciorek
#' 
#' @seealso \code{\link{calculateWAIC}} \code{\link{configureMCMC}}
#' \code{\link{buildMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
#'
#' @references 
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' \emph{Journal of Machine Learning Research} 11: 3571-3594.
#' 
#' Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models.
#' \emph{Statistics and Computing} 24(6): 997-1016.
#'
#' Ariyo, O., Quintero, A., Munoz, J., Verbeke, G. and Lesaffre, E. (2019).
#' Bayesian model selection in linear mixed models for longitudinal data.
#' \emph{Journal of Applied Statistics} 47: 890-913.
#'
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing} 27: 1413-1432.
#'
#' Hug, J.E.  and Paciorek, C.J. (2021). A numerically stable online
#' implementation and exploration of WAIC through variations of the
#' predictive density, using NIMBLE. \emph{arXiv e-print} <arXiv:2106.13359>.
#'
#'
#' @examples
#' code <- nimbleCode({
#'   for(j in 1:J) {
#'     for(i in 1:n) 
#'       y[j, i] ~ dnorm(mu[j], sd = sigma)
#'     mu[j] ~ dnorm(mu0, sd = tau)
#'   }
#'   sigma ~ dunif(0, 10)
#'   tau ~ dunif(0, 10)
#' })
#' J <- 5
#' n <- 10
#' groups <- paste0('y[', 1:J, ', 1:', n, ']') 
#' y <- matrix(rnorm(J*n), J, n)
#' Rmodel <- nimbleModel(code, constants = list(J = J, n = n), data = list(y = y),
#'                       inits = list(tau = 1, sigma = 1))
#'
#' ## Various versions of WAIC available via online calculation.
#' ## Conditional WAIC without data grouping:
#' conf <- configureMCMC(Rmodel, enableWAIC = TRUE)
#' ## Conditional WAIC with data grouping
#' conf <- configureMCMC(Rmodel, enableWAIC = TRUE, controlWAIC = list(dataGroups = groups))
#' ## Marginal WAIC with data grouping:
#' conf <- configureMCMC(Rmodel, enableWAIC = TRUE, controlWAIC =
#'             list(dataGroups = groups, marginalizeNodes = 'mu'))
#' \dontrun{
#' Rmcmc <- buildMCMC(conf)
#' Cmodel <- compileNimble(Rmodel)
#' Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#' output <- runMCMC(Cmcmc, niter = 1000, WAIC = TRUE)
#' output$WAIC              # direct access
#' ## Alternatively call via the `getWAIC` method; this doesn't require setting
#' ## `waic=TRUE` in `runMCMC`
#' Cmcmc$getWAIC()          
#' Cmcmc$getWAICdetails()
#' }
NULL
