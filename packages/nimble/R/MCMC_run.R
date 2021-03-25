


#' Run one or more chains of an MCMC algorithm and return samples, summary and/or WAIC
#'
#' Takes as input an MCMC algorithm (ideally a compiled one for speed)
#' and runs the MCMC with one or more chains, any returns any combination
#' of posterior samples, posterior summary statistics, and a WAIC value.
#'
#' @param mcmc A NIMBLE MCMC algorithm.  See details.
#'
#' @param niter Number of iterations to run each MCMC chain.  Default value is 10000.
#'
#' @param nburnin Number of initial, pre-thinning, MCMC iterations to discard.  Default value is 0.
#'
#' @param thin Thinning interval for collecting MCMC samples, corresponding to \code{monitors}.  Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#'
#' @param thin2 Thinning interval for collecting MCMC samples, corresponding to the second, optional set of \code{monitors2}.  Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#'
#' @param nchains Number of MCMC chains to run.  Default value is 1.
#'
#' @param inits Optional argument to specify initial values for each chain.  See details.
#'
#' @param setSeed Logical or numeric argument.  If a single numeric value is provided, R's random number seed will be set to this value at the onset of each MCMC chain.  If a numeric vector of length \code{nchains} is provided, then each element of this vector is provided as R's random number seed at the onset of the corresponding MCMC chain.  Otherwise, in the case of a logical value, if \code{TRUE}, then R's random number seed for the ith chain is set to be \code{i}, at the onset of each MCMC chain.  Note that specifying the argument \code{setSeed = 0} does not prevent setting the RNG seed, but rather sets the random number generation seed to \code{0} at the beginning of each MCMC chain.  Default value is \code{FALSE}.
#'
#' @param progressBar Logical argument.  If \code{TRUE}, an MCMC progress bar is displayed during execution of each MCMC chain.  Default value is defined by the nimble package option MCMCprogressBar.
#'
#' @param samples Logical argument.  If \code{TRUE}, then posterior samples are returned from each MCMC chain.  These samples are optionally returned as \code{coda} \code{mcmc} objects, depending on the \code{samplesAsCodaMCMC} argument.  Default value is \code{TRUE}.  See details.
#'
#' @param samplesAsCodaMCMC Logical argument.  If \code{TRUE}, then a \code{coda} \code{mcmc} object is returned instead of an R matrix of samples, or when \code{nchains > 1} a \code{coda} \code{mcmc.list} object is returned containing \code{nchains} \code{mcmc} objects.  This argument is only used when \code{samples} is \code{TRUE}.  Default value is \code{FALSE}.  See details.
#' 
#' @param samplesAsList Logical argument.  If \code{TRUE}, then samples are returned as a named list object, where each element corresponds to the monitored variable of that name.  List elements are vectors or arrays.  The first dimension of each element corresponds to the number of samples collected, and the subsequent dimensions of each element will match the dimensions of the variable itself.  When multiple chains are run, a list consisting of \code{nchains} named lists is returned, with each nested list containing the samples collected during one chain.  This argument is only used when \code{samples} is \code{TRUE}, and cannot be using concurrently with \code{samplesAsCodaMCMC = TRUE}.  Default value is \code{FALSE}.
#' 
#' @param summary Logical argument.  When \code{TRUE}, summary statistics for the posterior samples of each parameter are also returned, for each MCMC chain.  This may be returned in addition to the posterior samples themselves.  Default value is \code{FALSE}.  See details.
#'
#' @param WAIC Logical argument.  When \code{TRUE}, the WAIC (Watanabe, 2010) of the model is calculated and returned.  Note that in order for the WAIC to be calculated, the \code{mcmc} object must have also been created with the argument `enableWAIC = TRUE`.  If multiple chains are run, then a single WAIC value is calculated using the posterior samples from all chains.  Default value is \code{FALSE}.  See details.
#'
#' @return A list is returned with named elements depending on the arguments passed to \code{nimbleMCMC}, unless this list contains only a single element, in which case only that element is returned.  These elements may include \code{samples}, \code{summary}, and \code{WAIC}, and when the MCMC is monitoring a second set of nodes using \code{monitors2}, also \code{samples2}.  When \code{nchains = 1}, posterior samples are returned as a single matrix, and summary statistics as a single matrix.  When \code{nchains > 1}, posterior samples are returned as a list of matrices, one matrix for each chain, and summary statistics are returned as a list containing \code{nchains+1} matrices: one matrix corresponding to each chain, and the final element providing a summary of all chains, combined.  If \code{samplesAsCodaMCMC} is \code{TRUE}, then posterior samples are provided as \code{coda} \code{mcmc} and \code{mcmc.list} objects.  When \code{WAIC} is \code{TRUE}, a single WAIC value is returned.
#'
#' @details
#'
#' At least one of \code{samples}, \code{summary} or \code{WAIC} must be \code{TRUE}, since otherwise, nothing will be returned.  Any combination of these may be \code{TRUE}, including possibly all three, in which case posterior samples and summary statistics are returned for each MCMC chain, and an overall WAIC value is calculated and returned.
#'
#' When \code{samples = TRUE}, the form of the posterior samples is determined by the \code{samplesAsCodaMCMC} argument, as either matrices of posterior samples, or \code{coda} \code{mcmc} and \code{mcmc.list} objects.
#'
#' Posterior summary statistics are returned individually for each chain, and also as calculated from all chains combined (when \code{nchains > 1}).
#'
#' If provided, the \code{inits} argument can be one of three things:
#' 
#' (1) a function to generate initial values, which will be executed to generate initial values at the beginning of each MCMC chain, or
#' (2) a single named list of initial values which, will be used for each chain, or
#' (3) a list of length \code{nchains}, each element being a named list of initial values which be used for one MCMC chain.
#' 
#' The \code{inits} argument may also be omitted, in which case the current values in the \code{model} object will be used as the initial values of the first chain, and subsequent chains will begin using starting values where the previous chain ended.
#' 
#' Other aspects of the MCMC algorithm, such as the specific sampler assignments, must be specified in advance using the MCMC configuration object (created using \code{configureMCMC}), which is then used to build an MCMC algorithm (using \code{buildMCMC}) argument.
#'
#' The \code{niter} argument specifies the number of pre-thinning MCMC iterations, and the \code{nburnin} argument specifies the number of pre-thinning MCMC samples to discard.  After discarding these burn-in samples, thinning of the remaining samples will take place.  The total number of posterior samples returned will be floor((niter-nburnin)/thin).
#'
#' The MCMC option \code{mcmc$run(..., reset = FALSE)}, used to continue execution of an MCMC chain, is not available through \code{runMCMC()}.
#' 
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, sd = 1000)
#'     sigma ~ dunif(0, 1000)
#'     for(i in 1:10) {
#'         x[i] ~ dnorm(mu, sd = sigma)
#'     }
#' })
#' Rmodel <- nimbleModel(code)
#' Rmodel$setData(list(x = c(2, 5, 3, 4, 1, 0, 1, 3, 5, 3)))
#' Rmcmc <- buildMCMC(Rmodel)
#' Cmodel <- compileNimble(Rmodel)
#' Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#' inits <- function() list(mu = rnorm(1,0,1), sigma = runif(1,0,10))
#' samplesList <- runMCMC(Cmcmc, niter = 10000, nchains = 3, inits = inits)
#' }
#'
#' @seealso \code{\link{configureMCMC}} \code{\link{buildMCMC}} \code{\link{nimbleMCMC}}
#'
#' @author Daniel Turek
#'
#' @export
runMCMC <- function(mcmc,
                    niter = 10000,
                    nburnin = 0,
                    thin,
                    thin2,
                    nchains = 1,
                    inits,
                    ## reinstate samplerExecutionOrder as a runtime argument, once we support non-scalar default values for runtime arguments:
                    ##samplerExecutionOrder,
                    setSeed = FALSE,
                    progressBar = getNimbleOption('MCMCprogressBar'),
                    samples = TRUE,
                    samplesAsCodaMCMC = FALSE,
                    samplesAsList = FALSE,
                    summary = FALSE,
                    WAIC = FALSE) {
    if(missing(mcmc)) stop('must provide a NIMBLE MCMC algorithm')
    if(!identical(nfGetDefVar(mcmc, 'name'), 'MCMC')) stop('mcmc argument must be a NIMBLE MCMC algorithm')
    if(!is.Cnf(mcmc)) message('Warning: running an uncompiled MCMC algorithm, use compileNimble() for faster execution.')
    if(!samples && !summary && !WAIC) stop('no output specified, use samples = TRUE, summary = TRUE, or WAIC = TRUE')
    if(samples && samplesAsCodaMCMC && samplesAsList) stop('cannot specify both samplesAsCodaMCMC = TRUE and samplesAsList = TRUE')
    if(nchains < 1) stop('must have nchains > 0')
    if(!missing(inits)) {
        if(!is.function(inits) && !is.list(inits)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
        if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]]) && (length(inits) != nchains)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
    }
    if(WAIC && !mcmc$enableWAIC) stop('mcmc argument must have been created with "enableWAIC = TRUE" in order to calculate WAIC.')
    model <- if(is.Cnf(mcmc)) mcmc$Robject$model$CobjectInterface else mcmc$model
    if(!is.model(model)) stop('something went wrong')
    hasMonitors2 <- length(if(is.Cnf(mcmc)) mcmc$Robject$monitors2 else mcmc$monitors2) > 0
    samplesList  <- vector('list', nchains); names(samplesList)  <- paste0('chain', 1:nchains)
    samplesList2 <- vector('list', nchains); names(samplesList2) <- paste0('chain', 1:nchains)
    samplesReturnList <- samplesReturnList2 <- vector('list', nchains)
    thinToUseVec <- c(0, 0)
    thinToUseVec[1] <- if(!missing(thin))  thin  else mcmc$thinFromConfVec[1]
    thinToUseVec[2] <- if(!missing(thin2)) thin2 else mcmc$thinFromConfVec[2]
    ## if(thinToUseVec[1] > 1 && nburnin > 0) message("runMCMC's handling of nburnin changed in nimble version 0.6-11. Previously, nburnin samples were discarded *post-thinning*.  Now nburnin samples are discarded *pre-thinning*.  The number of samples returned will be floor((niter-nburnin)/thin).")
    ## reinstate samplerExecutionOrder as a runtime argument, once we support non-scalar default values for runtime arguments:
    ##samplerExecutionOrderToUse <- if(!missing(samplerExecutionOrder)) samplerExecutionOrder else mcmc$samplerExecutionOrderFromConfPlusTwoZeros[mcmc$samplerExecutionOrderFromConfPlusTwoZeros>0]
    for(i in 1:nchains) {
        if(nimbleOptions('verbose')) message('running chain ', i, '...')
        ##if(setSeed) set.seed(i)
        if(is.numeric(setSeed)) {
            if(length(setSeed) == 1) {
                set.seed(setSeed)
            } else { if(length(setSeed) == nchains) set.seed(setSeed[i]) else stop('setSeed argument has different length from nchains') }
        } else if(setSeed) set.seed(i)
        if(!missing(inits)) {
            if(is.function(inits)) {
                theseInits <- inits()
            } else if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]])) {
                theseInits <- inits[[i]]
            } else theseInits <- inits
            model$setInits(theseInits)
        }
        ##model$calculate()   # shouldn't be necessary, since mcmc$run() includes call to my_initializeModel$run()
        mcmc$run(niter, nburnin = nburnin, thin = thinToUseVec[1], thin2 = thinToUseVec[2], progressBar = progressBar) #, samplerExecutionOrder = samplerExecutionOrderToUse)
        samplesList[[i]] <- as.matrix(mcmc$mvSamples)
        if(hasMonitors2)   samplesList2[[i]] <- as.matrix(mcmc$mvSamples2)
        if(samplesAsList) { samplesReturnList[[i]] <- as.list(mcmc$mvSamples)
                            if(hasMonitors2)   samplesReturnList2[[i]] <- as.list(mcmc$mvSamples2) }
    }
    if(WAIC) {
        if(nchains > 1) {
            samplesPerChain <- dim(samplesList[[1]])[1]
            posteriorSamplesMatrix <- matrix(0, nrow = samplesPerChain*nchains, ncol = dim(samplesList[[1]])[2])
            for(i in seq_along(samplesList)) {
                posteriorSamplesMatrix[((i-1)*samplesPerChain + 1):(i*samplesPerChain),] <- samplesList[[i]][,]
            }
            colnames(posteriorSamplesMatrix) <- colnames(samplesList[[1]])
            matrix2mv(posteriorSamplesMatrix, mcmc$mvSamples)  ## transfer all posterior samples into mcmc$mvSamples
        }
        WAICvalue <- mcmc$calculateWAIC()
    }
    if(samplesAsCodaMCMC) {
        samplesList <- as.mcmc.list(lapply(samplesList, as.mcmc))
        if(hasMonitors2)   samplesList2 <- as.mcmc.list(lapply(samplesList2, as.mcmc))
    }
    if(nchains == 1) {
        samplesList <- samplesList[[1]]                       ## returns matrix when nchains = 1
        if(hasMonitors2)   samplesList2 <- samplesList2[[1]]  ## returns matrix when nchains = 1
        if(samplesAsList) { samplesReturnList <- samplesReturnList[[1]]
                            if(hasMonitors2)   samplesReturnList2 <- samplesReturnList2[[1]] }
    }
    if(summary) {
        if(nchains == 1) {
            summaryObject <- samplesSummary(samplesList)
            if(hasMonitors2)   summaryObject <- rbind(summaryObject, samplesSummary(samplesList2))   ## combine summaries
        } else {
            summaryObject <- lapply(samplesList, samplesSummary)
            names(summaryObject) <- paste0('chain', 1:nchains)
            summaryObject$all.chains <- samplesSummary(do.call('rbind', samplesList))
            if(hasMonitors2) {
                summaryObject2 <- lapply(samplesList2, samplesSummary)
                summaryObject2$all.chains <- samplesSummary(do.call('rbind', samplesList2))
                summaryObject <- mapply(rbind, summaryObject, summaryObject2, SIMPLIFY = FALSE)      ## combine summaries
            }
        }
    }
    retList <- list()
    if(samples) { retList$samples <- if(samplesAsList) samplesReturnList else samplesList
                  if(hasMonitors2)   retList$samples2 <- if(samplesAsList) samplesReturnList2 else samplesList2 }
    if(summary)   retList$summary <- summaryObject
    if(WAIC)      retList$WAIC    <- WAICvalue
    if(length(retList) == 1) retList <- retList[[1]]
    return(retList)
}


#' Executes one or more chains of NIMBLE's default MCMC algorithm, for a model specified using BUGS code
#'
#' \code{nimbleMCMC} is designed as the most straight forward entry point to using NIMBLE's default MCMC algorithm.  It provides capability for running multiple MCMC chains, specifying the number of MCMC iterations, thinning, and burn-in, and which model variables should be monitored.  It also provides options to return the posterior samples, to return summary statistics calculated from the posterior samples, and to return a WAIC value.
#'
#' The entry point for this function is providing the \code{code}, \code{constants}, \code{data} and \code{inits} arguments, to create a new NIMBLE model object, or alternatively providing an exisiting NIMBLE model object as the \code{model} argument.
#'
#' @param code The quoted code expression representing the model, such as the return value from a call to \code{nimbleCode}). Not required if \code{model} is provided.  
#' 
#' @param constants Named list of constants in the model.  Constants cannot be subsequently modified. For compatibility with JAGS and BUGS, one can include data values with constants and \code{nimbleModel} will automatically distinguish them based on what appears on the left-hand side of expressions in \code{code}.
#' 
#' @param data Named list of values for the data nodes.  Data values can be subsequently modified.  Providing this argument also flags nodes as having data for purposes of algorithms that inspect model structure. Values that are NA will not be flagged as data.
#'
#' @param inits Argument to specify initial values for the model object, and for each MCMC chain.  See details.
#'
#' @param dimensions Named list of dimensions for variables.  Only needed for variables used with empty indices in model code that are not provided in constants or data.
#'
#' @param model A compiled or uncompiled NIMBLE model object.  When provided, this model will be used to configure the MCMC algorithm to be executed, rather than using the \code{code}, \code{constants}, \code{data} and \code{inits} arguments to create a new model object.  However, if also provided, the \code{inits} argument will still be used to initialize this model prior to running each MCMC chain.
#' 
#' @param monitors A character vector giving the node names or variable names to monitor.  The samples corresponding to these nodes will returned, and/or will have summary statistics calculated. Default value is all top-level stochastic nodes of the model.
#' 
#' @param thin Thinning interval for collecting MCMC samples.  Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#' 
#' @param niter Number of MCMC iterations to run.  Default value is 10000.
#' 
#' @param nburnin Number of initial, pre-thinning, MCMC iterations to discard.  Default value is 0.
#' 
#' @param nchains Number of MCMC chains to run.  Default value is 1.
#' 
#' @param check Logical argument, specifying whether to check the model object for missing or invalid values.  Default value is \code{TRUE}.
#' 
#' @param setSeed Logical or numeric argument.  If a single numeric value is provided, R's random number seed will be set to this value at the onset of each MCMC chain.  If a numeric vector of length \code{nchains} is provided, then each element of this vector is provided as R's random number seed at the onset of the corresponding MCMC chain.  Otherwise, in the case of a logical value, if \code{TRUE}, then R's random number seed for the ith chain is set to be \code{i}, at the onset of each MCMC chain.  Note that specifying the argument \code{setSeed = 0} does not prevent setting the RNG seed, but rather sets the random number generation seed to \code{0} at the beginning of each MCMC chain.  Default value is \code{FALSE}.
#'
#' @param progressBar Logical argument.  If \code{TRUE}, an MCMC progress bar is displayed during execution of each MCMC chain.  Default value is defined by the nimble package option MCMCprogressBar..
#'
#' @param samples Logical argument.  If \code{TRUE}, then posterior samples are returned from each MCMC chain.  These samples are optionally returned as \code{coda} \code{mcmc} objects, depending on the \code{samplesAsCodaMCMC} argument.  Default value is \code{TRUE}.  See details.
#'
#' @param samplesAsCodaMCMC Logical argument.  If \code{TRUE}, then a \code{coda} \code{mcmc} object is returned instead of an R matrix of samples, or when \code{nchains > 1} a \code{coda} \code{mcmc.list} object is returned containing \code{nchains} \code{mcmc} objects.  This argument is only used when \code{samples} is \code{TRUE}.  Default value is \code{FALSE}.  See details.
#' 
#' @param samplesAsList Logical argument.  If \code{TRUE}, then samples are returned as a named list object, where each element corresponds to the monitored variable of that name.  List elements are vectors or arrays.  The first dimension of each element corresponds to the number of samples collected, and the subsequent dimensions of each element will match the dimensions of the variable itself.  When multiple chains are run, a list consisting of \code{nchains} named lists is returned, with each nested list containing the samples collected during one chain.  This argument is only used when \code{samples} is \code{TRUE}, and cannot be using concurrently with \code{samplesAsCodaMCMC = TRUE}.  Default value is \code{FALSE}.
#'
#' @param summary Logical argument.  When \code{TRUE}, summary statistics for the posterior samples of each parameter are also returned, for each MCMC chain.  This may be returned in addition to the posterior samples themselves.  Default value is \code{FALSE}.  See details.
#'z
#' @param WAIC Logical argument.  When \code{TRUE}, the WAIC (Watanabe, 2010) of the model is calculated and returned.  If multiple chains are run, then a single WAIC value is calculated using the posterior samples from all chains.  Default value is \code{FALSE}.  See details.
#' 
#' @return A list is returned with named elements depending on the arguments passed to \code{nimbleMCMC}, unless only one among samples, summary, and WAIC are requested, in which case only that element is returned.  These elements may include \code{samples}, \code{summary}, and \code{WAIC}.  When \code{nchains = 1}, posterior samples are returned as a single matrix, and summary statistics as a single matrix.  When \code{nchains > 1}, posterior samples are returned as a list of matrices, one matrix for each chain, and summary statistics are returned as a list containing \code{nchains+1} matrices: one matrix corresponding to each chain, and the final element providing a summary of all chains, combined.  If \code{samplesAsCodaMCMC} is \code{TRUE}, then posterior samples are provided as \code{coda} \code{mcmc} and \code{mcmc.list} objects.  When \code{WAIC} is \code{TRUE}, a single WAIC value is returned.
#'
#' @details
#'
#' At least one of \code{samples}, \code{summary} or \code{WAIC} must be \code{TRUE}, since otherwise, nothing will be returned.  Any combination of these may be \code{TRUE}, including possibly all three, in which case posterior samples, summary statistics, and WAIC values are returned for each MCMC chain.
#'
#' When \code{samples = TRUE}, the form of the posterior samples is determined by the \code{samplesAsCodaMCMC} argument, as either matrices of posterior samples, or \code{coda} \code{mcmc} and \code{mcmc.list} objects.
#'
#' Posterior summary statistics are returned individually for each chain, and also as calculated from all chains combined (when \code{nchains > 1}).
#'
#' The \code{inits} argument can be one of three things:
#' 
#' (1) a function to generate initial values, which will be executed once to initialize the model object, and once to generate initial values at the beginning of each MCMC chain, or
#' (2) a single named list of initial values which, will be used to initialize the model object and for each MCMC chain, or
#' (3) a list of length \code{nchains}, each element being a named list of initial values.  The first element will be used to initialize the model object, and once element of the list will be used for each MCMC chain.
#' 
#' The \code{inits} argument may also be omitted, in which case the model will not be provided with initial values.  This is not recommended.
#'
#' The \code{niter} argument specifies the number of pre-thinning MCMC iterations, and the \code{nburnin} argument specifies the number of pre-thinning MCMC samples to discard.  After discarding these burn-in samples, thinning of the remaining samples will take place.  The total number of posterior samples returned will be floor((niter-nburnin)/thin).
#' 
#' @examples
#' 
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, sd = 1000)
#'     sigma ~ dunif(0, 1000)
#'     for(i in 1:10) {
#'         x[i] ~ dnorm(mu, sd = sigma)
#'     }
#' })
#' data <- list(x = c(2, 5, 3, 4, 1, 0, 1, 3, 5, 3))
#' inits <- function() list(mu = rnorm(1,0,1), sigma = runif(1,0,10))
#' mcmc.output <- nimbleMCMC(code, data = data, inits = inits,
#'                           monitors = c("mu", "sigma"), thin = 10,
#'                           niter = 20000, nburnin = 1000, nchains = 3,
#'                           summary = TRUE, WAIC = TRUE)
#' }
#'
#' @seealso \code{\link{configureMCMC}} \code{\link{buildMCMC}} \code{\link{runMCMC}}
#' 
#' @author Daniel Turek
#' 
#' @export
nimbleMCMC <- function(code,
                       constants = list(),
                       data = list(),
                       inits,
                       dimensions = list(),
                       model,
                       monitors,
                       thin = 1,
                       niter = 10000,
                       nburnin = 0,
                       nchains = 1,
                       check = TRUE,
                       setSeed = FALSE,
                       progressBar = getNimbleOption('MCMCprogressBar'),
                       samples = TRUE,
                       samplesAsCodaMCMC = FALSE,
                       samplesAsList = FALSE,
                       summary = FALSE,
                       WAIC = FALSE) {
    #### process 'code' argument, to accept a filename, or a function
    ##if(is.character(code) || is.function(code)) {
    ##    if(is.function(code)) modelText <- mergeMultiLineStatements(deparse(body(code)))
    ##    if(is.character(code)) {
    ##        if(!file.exists(code)) stop("'code' argument does not reference an existing file.")
    ##        modelText <- mergeMultiLineStatements(processModelFile(code)$modelLines)
    ##    }
    ##    modelText <- processNonParseableCode(modelText)  ## deal with T() and I() syntax
    ##    code <- parse(text = modelText)[[1]]
    ##}
    if(missing(code) && missing(model)) stop('must provide either code or model argument')
    if(!samples && !summary && !WAIC) stop('no output specified, use samples = TRUE, summary = TRUE, or WAIC = TRUE')
    if(!missing(code) && inherits(code, 'modelBaseClass')) model <- code   ## let's handle it, if model object is provided as un-named first argument to nimbleMCMC
    if(missing(model)) {  ## model object not provided
        if(!missing(inits)) {
            if(!is.function(inits) && !is.list(inits)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
            if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]]) && (length(inits) != nchains)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
            if(is.function(inits)) {
                if(is.numeric(setSeed) || setSeed) { if(is.numeric(setSeed)) set.seed(setSeed[1]) else set.seed(0) }
                theseInits <- inits()
            } else if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]])) {
                theseInits <- inits[[1]]
            } else theseInits <- inits
            Rmodel    <- nimbleModel(code, constants, data, theseInits, dimensions = dimensions, check = check)    ## inits provided
        } else Rmodel <- nimbleModel(code, constants, data,             dimensions = dimensions, check = check)    ## inits not provided
    } else {              ## model object provided
        if(!is.model(model)) stop('model argument must be a NIMBLE model object')
        Rmodel <- if(is.Rmodel(model)) model else model$Rmodel
        if(!is.Rmodel(Rmodel)) stop('something went wrong')
    }
    conf <- configureMCMC(Rmodel, monitors = monitors, thin = thin, enableWAIC = WAIC, print = FALSE)
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(Rmodel, Rmcmc)    ## only one compileNimble() call
    Cmcmc <- compiledList$Rmcmc
    runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits,
            setSeed = setSeed, progressBar = progressBar, samples = samples,
            samplesAsCodaMCMC = samplesAsCodaMCMC, samplesAsList = samplesAsList,
            summary = summary, WAIC = WAIC)
}


