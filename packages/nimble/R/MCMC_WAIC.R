#' waicList definition
#' 
#' \code{waicList} definition for the \code{nimbleList} type returned by WAIC
#' computation.
#'
#' @details
#'
#' See \code{help(waic)} for details on the elements of the list.
#' 
#' @author NIMBLE development team
#'
#' @export
#' 
waicList <- nimbleList(
    list(
        nimbleType('WAIC', 'double', 0),
        nimbleType('lppd', 'double', 0),
        nimbleType('pWAIC', 'double', 0),
        nimbleType('marginal', 'logical', 0),
        nimbleType('nItsMarginal', 'double', 0),
        nimbleType('thin', 'logical', 0),
        nimbleType('online', 'logical', 0),

        ## values for shorter MC runs to assess convergence for marginal calculation
        nimbleType('WAIC_partialMC', 'double', 1),
        nimbleType('lppd_partialMC', 'double', 1),
        nimbleType('pWAIC_partialMC', 'double', 1),
        nimbleType('nItsMarginal_partialMC', 'double' , 1),  # checkIts

        ## per data group values potentially useful for SE for contrasting WAIC of two models
        nimbleType('WAIC_elements', 'double', 1),
        nimbleType('lppd_elements', 'double', 1),
        nimbleType('pWAIC_elements', 'double', 1)

    ), name = 'waicList'
)
    

buildWAIC <- nimbleFunction(
    name = 'waicClass',
    contains = nimbleFunctionVirtual(), 
    setup = function(model, mvSaved, control) {
        online              <- extractControlElement(control, 'online',              TRUE)
        thin                <- extractControlElement(control, 'thin',                FALSE)
        dataGroups          <- extractControlElement(control, 'dataGroups',          NULL)
        latentNodes         <- extractControlElement(control, 'marginalizeNodes',    NULL)
        nItsMarginal        <- extractControlElement(control, 'nItsMarginal',        1000)
        convergenceSet      <- extractControlElement(control, 'convergenceSet',      NULL)

        setupOutputs(online, thin)
        
        ## Partition of data nodes. By default one data node per group.
        dataNodes <- model$getNodeNames(dataOnly = TRUE)
        dataNodeLength <- length(dataNodes)
        if(!dataNodeLength)
            stop('buildWAIC: cannot compute WAIC with no data nodes in the model.')
        if(!is.null(dataGroups)) { 
            useGroups <- TRUE
            dataNodes <- lapply(dataGroups, function(x) model$expandNodeNames(x))
            groupIndices <- lapply(dataNodes, length)
            groupIndices <- cumsum(groupIndices)
            dataNodes <- unlist(dataNodes)
            if (length(dataNodes) != dataNodeLength ) {
                warning("buildWAIC: Group nodes supplied does not contain all data nodes or contains duplicates.")
            }
        } else{
            useGroups <- FALSE
            groupIndices <- rep(1, dataNodeLength)
        }
        nGroups <- length(groupIndices)
        ## Cannot have a length 1 vector.
        if(nGroups == 1) 
            groupIndices <- c(groupIndices, 0)
        
        ## Determine marginal versus conditional WAIC.
        if(!is.null(latentNodes)) {
            marginal <- TRUE
            latentNodes <- model$getDependencies(latentNodes, self = TRUE, downstream = TRUE, includeData = FALSE)
            latentAndDataNodes <- model$getDependencies(latentNodes, self = TRUE, downstream = TRUE)
        } else {
            marginal  <- FALSE 
            nItsMarginal <- 1
            ## Next two are unused in conditional case but need to exist for compilation.
            latentNodes <- dataNodes[1]
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
                checkIts <- convergenceSet * nItsMarginal
                lengthConvCheck <- length(checkIts) + 1
            } else {  # Default to 25%, 50%, 75%.
                checkIts <- floor(quantile(1:nItsMarginal, names = FALSE))[c(2,3,4)]
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
            ticker <- 0  # iterates over MC subsets
            mcmcIter <<- mcmcIter + 1
            for (k in 1:nItsMarginal) {  # loop over MC samples; for conditional, only one 'iteration'
                if(marginal) {
                    model$simulate(latentNodes)
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
            logProbMat[lengthConvCheck, ] <<- (margSumMax + log(margCurSum) - log(nItsMarginal))[1:nGroups]
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
        },
        get = function(returnElements = logical(default = FALSE)) {
            ## Extract WAIC summary information.
            returnType(waicList())
            if(!finalized)
                finalize()
            badpWAIC <- sum( sspWAICmat[lengthConvCheck, ] / (mcmcIter-1) > 0.4 )
            if(badpWAIC) {   # Check with Josh for how he decided on this cutoff.
                print("There are individual pWAIC values that are greater than 0.4, WAIC estimate may be unstable" )
            }
            
            output <- waicList$new()
            
            output$WAIC <- WAIC[lengthConvCheck]
            output$lppd <- lppd[lengthConvCheck]
            output$pWAIC  <- pWAIC[lengthConvCheck]
            output$marginal <- marginal
            output$thin <- thin
            output$online <- online

            if(marginal) {
                output$nItsMarginal <- nItsMarginal
                output$WAIC_partialMC <- WAIC[1:(lengthConvCheck - 1)]
                output$lppd_partialMC <- lppd[1:(lengthConvCheck - 1)]
                output$pWAIC_partialMC <- pWAIC[1:(lengthConvCheck - 1)]
                output$nItsMarginal_partialMC <- checkIts
            } else output$nItsMarginal <- 0
            
            if(returnElements) {
                ## Per-observations values useful for (manual) SE calculation when
                ## contrasting WAIC for two models, per Vehtari et al. 2017 (Staistics and Computing)
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
        }
    )
)



#' Using WAIC
#'
#' Details of the WAIC measure for comparing models. NIMBLE implements an online
#' WAIC algorithm, computed during the course of the MCMC iterations.
#' 
#' @param enableWAIC A logical argument, specifying whether to enable WAIC
#' calculations for the resulting MCMC algorithm.  Defaults to the value of
#' \code{nimbleOptions('MCMCenableWAIC')}, which in turn defaults to FALSE.
#' 
#' @param waicControl An optional named list of inputs that control the
#' behavior of the WAIC calculation. See `Details`.
#'
#' @name waic
#' 
#' @aliases getWAIC calculateWAIC buildWAIC
#' 
#' @details
#'
#' To use WAIC, set \code{enableWAIC = TRUE} when configuring or building
#' an MCMC and set \code{WAIC = TRUE} when calling \code{nimbleMCMC} and
#' optionally when calling \code{runMCMC}.
#'
#' By default, NIMBLE calculates WAIC using an online algorithm that updates
#' required summary statistics at each post-burnin iteration of the MCMC. The
#' ability to calculate WAIC post hoc after all MCMC sampling has been done has
#' been retained for compatibility with versions of NIMBLE before 0.12.0. See
#' below for details on the two different approaches.
#'
#' @section \code{waicControl} list:
#'
#' The \code{waicControl} argument is a list that controls the behavior of the
#' WAIC algorithm and is passed to either \code{configureMCMC} or
#' \code{buildMCMC}. One can supply any of the following optional components:
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
#' \code{nItsMarginal}: Number of Monte Carlo iterations to use when
#' marginalizing (default is 1000).
#'
#' \code{convergenceSet}: Optional vector of numbers between 0 and 1 that
#' specify a set of shorter Monte Carlo simulations for marginal WAIC
#' calculation as fractions of the full (\code{nItsMarginal}) Monte Carlo
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
#' object. If using the \code{run} method of the MCMC object, use the
#' \code{getWAIC} method of the MCMC object.
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
#' \code{marginal}: Logical value indicating whether marginal (\code{TRUE}) or
#' conditional (\code{FALSE}) WAIC was calculated.
#' 
#' \code{nItsMarginal}: Number of Monte Carlo iterations used in computing
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
#' \code{nItsMarginal_partialMC}: Number of Monte Carlo iterations used for the
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
#' using all post-burnin samples.
#'
#' The \code{calculateWAIC} method calculates the WAIC of the model that the
#' MCMC was performed on using post-burnin MCMC samples.
#' The WAIC (Watanabe, 2010) is calculated from
#' Equations 5, 12, and 13 in Gelman et al. (2014) (i.e., using 'pWAIC2').
#'
#' Note that there is not a unique value of WAIC for a model. By default, WAIC
#' is calculated conditional on the parent nodes of the data nodes, and the
#' density values used are the individual density values of the data nodes.
#' However, by modifying the \code{marginalizeNodes} and \code{dataGroups}
#' elements of the control list, users can request a marginal WAIC (using a
#' marginal likelihood that integrates over user-specified latent nodes) and/or
#' a WAIC based on grouping observations (e.g., all observations in a cluster)
#' to use joint density values.
#'
#' For more detail on the use of different predictive distributions, see Section
#' 2.5 from Gelman et al. (2014) or Ariyo et al. (2019).
#'
#' Note that based on a limited set of simulation experiments in Hug and Paciorek
#' (2021) we recommend users only use marginal WAIC if also using grouping. 
#' 
#' @section WAIC computed after MCMC sampling:
#'
#' For compatibility with older versions of NIMBLE (prior to v0.12.0), one can
#' also calculate WAIC after all sampling has been completed, provided
#' the necessary variables have been monitored in the MCMC.
#'
#' After the MCMC has been run, calling the \code{calculateWAIC()} method of the
#' MCMC object will return the WAIC for the model, calculated using the posterior
#' samples from the MCMC run.
#' 
#' \code{calculateWAIC()} accepts a single argument:
#'
#' \code{nburnin}: The number of pre-thinning MCMC samples to remove from the
#' beginning of the posterior samples for WAIC calculation (default = 0). These
#' samples are discarded in addition to any burn-in specified when running the
#' MCMC.
#' 
#' The \code{calculateWAIC} method can only be used if the \code{enableWAIC} 
#' argument to \code{configureMCMC} or to \code{buildMCMC} is set to \code{TRUE},
#' or if the NIMBLE option \code{enableWAIC} is set to \code{TRUE}.  If a user
#' attempts to call \code{calculateWAIC} without having set
#' \code{enableWAIC = TRUE} (either in the call to \code{configureMCMC}, or
#' \code{buildMCMC}, or as a NIMBLE option), an error will occur.  
#' 
#' The \code{calculateWAIC} method calculates the WAIC of the model that the
#' MCMC was performed on. The WAIC (Watanabe, 2010) is calculated from
#' Equations 5, 12, and 13 in Gelman et al. (2014) (i.e., using \emph{p}WAIC2).
#'
#' Note that there is not a unique value of WAIC for a model.
#' \code{calculateWAIC} only provides the conditional WAIC, namely the version
#' of WAIC where all parameters directly involved in the likelihood are treated
#' as \eqn{theta} for the purposes of Equation 5 from Gelman et al. (2014). As
#' a result, the user must set the MCMC monitors (via the \code{monitors}
#' argument) to include all stochastic nodes that are parents of any data nodes;
#' by default the MCMC monitors are only the top-level nodes of the model. For
#' more detail on the use of different predictive distributions, see Section 2.5
#' from Gelman et al. (2014) or Ariyo et al. (2019).

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
#' `Online WAIC` above.
#'
#' @author Joshua Hug and Christopher Paciorek
#' 
#' @seealso \code{\link{configureMCMC}} \code{\link{buildMCMC}}
#' \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
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
#' Hug, J.E.  and Paciorek, C.J. (2021). A numerically stable online implementation and exploration of WAIC through variations of the predictive density, using NIMBLE. \emph{arXiv} xx.yy.
#'
#'
#' @examples
#' code <- nimbleCode({
#'   for(j in 1:J) {
#'     for(i in 1:n) 
#'       y[j, i] ~ dnorm(mu[j], 1)
#'     mu[j] ~ dnorm(mu0, 1)
#'   }
#' })
#' J <- 5
#' n <- 10
#' groups <- paste0('y[', 1:J, ', 1:', n, ']') 
#' y <- matrix(rnorm(J*n), J, n)
#' Rmodel <- nimbleModel(code, constants = list(J = J, n = n), data = list(y = y)) 
#' ## Conditional WAIC without data grouping:
#' conf <- configureMCMC(Rmodel, enableWAIC = TRUE)
#' ## Conditional WAIC with data grouping
#' conf <- configureMCMC(Rmodel, enableWAIC = TRUE, waicControl = list(dataGroups = groups))
#' ## Marginal WAIC with data grouping:
#' conf <- configureMCMC(Rmodel, enableWAIC = TRUE, waicControl =
#'             list(dataGroups = groups, marginalizeNodes = 'mu'))
NULL
