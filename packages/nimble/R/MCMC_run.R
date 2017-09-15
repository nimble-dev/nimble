

#' @export
samplesSummary <- function(samples) {
    cbind(
        `Mean`      = apply(samples, 2, mean),
        `Median`    = apply(samples, 2, median),
        `St.Dev.`   = apply(samples, 2, sd),
        `95%CI_low` = apply(samples, 2, function(x) quantile(x, 0.025)),
        `95%CI_upp` = apply(samples, 2, function(x) quantile(x, 0.975)))
}


#' Run one or more chains of an MCMC algorithm and return samples, summary and/or WAIC
#'
#' Takes as input an MCMC algorithm (ideally a compiled one for speed)
#' and runs the MCMC with one or more chains, any returns any combination
#' of posterior samples, posterior summary statistics, and WAIC values.
#'
#' @param mcmc A NIMBLE MCMC algorithm.  See details.
#'
#' @param niter Number of iterations to run each MCMC chain (default = 10000).
#'
#' @param nburnin Number of initial samples to discard from each MCMC chain (default = 0).
#'
#' @param nchains Number of MCMC chains to run (default = 1).
#'
#' @param inits Optional argument to specify initial values for each chain.  See details.
#'
#' @param setSeed Logical argument.  If \code{TRUE}, then R's random number seed is set to \code{i} (using \code{set.seed(i)}) at the onset of each MCMC chain number \code{i} (default = \code{FALSE}).
#'
#' @param progressBar Logical argument.  If \code{TRUE}, an MCMC progress bar is displayed during execution of each MCMC chain (default = \code{TRUE}).
#'
#' @param returnSamples.  Logical argument.  If \code{TRUE}, then posterior samples are returned from each MCMC chain.  These samples are optionally returned as \code{coda} \code{mcmc} objects, depending on the \code{samplesAsCodaMCMC} argument.  Default value is \code{TRUE}.  See details.
#'
#' @param samplesAsCodaMCMC Logical argument.  If \code{TRUE}, then a \code{coda} \code{mcmc} object is returned instead of an R matrix of samples, or when \code{nchains > 1} a \code{coda} \code{mcmc.list} object is returned containing \code{nchains} \code{mcmc} objects.  This argument is only used when \code{returnSamples} is \code{TRUE}.  Default value is \code{FALSE}.  See details.
#' 
#' @param returnSummary Logical argument.  When \code{TRUE}, summary statistics for the posterior samples of each parameter are also returned, for each MCMC chain.  This may be returned in addition to the posterior samples themselves.  Default value is \code{FALSE}.  See details.
#'
#' @param returnWAIC Logical argument.  When \code{TRUE}, the WAIC (Watanabe, 2010) of the model is calculated and returned.  If multiple chains are run, then WAIC is calculated separately for each chain.  Default value is \code{FALSE}.  See details.
#'
#' @return A list is returned with named elements depending on the arguments passed to \code{nimbleMCMC}, unless only one among samples, summary, and WAIC are requested, in which case only that element is returned.  These elements may include \code{samples}, \code{summary}, and \code{waic}.  When \code{nchains = 1}, posterior samples are returned as a single matrix, and summary statistics as a single matrix.  When \code{nchains > 1}, posterior samples are returned as a list of matrices, one matrix for each chain, and summary statistics are returned as a list containing \code{nchains+1} matrices: one matrix corresponding to each chain, and the final element providing a summary of all chains, combined.  If \code{samplesAsCodaMCMC} is \code{TRUE}, then posterior samples are provided as \code{coda} \code{mcmc} and \code{mcmc.list} objects.  When \code{returnWAIC} is \code{TRUE}, WAIC values are returned as a numeric vector with one element corresponding to each chain.
#'
#' @details
#'
#' At least one of \code{returnSamples}, \code{returnSummary} or \code{returnWAIC} must be \code{TRUE}, since otherwise, nothing will be returned.  Any combination of these may be \cd{TRUE}, including possibly all three, in which case posterior samples, summary statistics, and WAIC values are returned for each MCMC chain.
#'
#' When \code{returnSamples = TRUE}, the form of the posterior samples is determined by the \code{samplesAsCodaMCMC} argument, as either matrices of posterior samples, or \cd{coda} \cd{mcmc} and \cd{mcmc.list} objects.
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
#' Other aspects of the MCMC algorithm, such as sampler assignments and thinning, must be specified in advance using the MCMC configuration object (created using \code{configureMCMC}), which is then used to build the MCMC algorithm (using \code{buildMCMC}) argument.
#'
#' The \code{niter} argument specifies the number of pre-thinning MCMC iterations, and the \code{nburnin} argument will remove post-thinning samples.
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
                    niter = 10000, nburnin = 0, nchains = 1,
                    inits,
                    setSeed = FALSE,
                    progressBar = TRUE,
                    returnSamples = TRUE,
                    samplesAsCodaMCMC = FALSE,
                    returnSummary = FALSE,
                    returnWAIC = FALSE) {
    if(missing(mcmc)) stop('must provide a NIMBLE MCMC algorithm')
    if(!identical(nf_getGeneratorFunction(mcmc), buildMCMC)) stop('mcmc argument must be a NIMBLE MCMC algorithm')
    if(!is.Cnf(mcmc)) message('Warning: running an uncompiled MCMC algorithm, use compileNimble() for faster execution.')
    if(!returnSamples && !returnSummary && !returnWAIC) stop('no output specified, use returnSamples = TRUE, returnSummary = TRUE, or returnWAIC = TRUE')
    if(samplesAsCodaMCMC) require(coda)
    if(nchains < 1) stop('must have nchains > 0')
    if(!missing(inits)) {
        if(!is.function(inits) && !is.list(inits)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
        if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]]) && (length(inits) != nchains)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
    }
    model <- if(is.Cnf(mcmc)) mcmc$Robject$model$CobjectInterface else mcmc$model
    if(!is.model(model)) stop('something went wrong')
    samplesList <- vector('list', nchains)
    names(samplesList) <- paste0('chain', 1:nchains)
    if(returnWAIC) {
        waic <- numeric(nchains)
        if(nchains > 1) names(waic) <- paste0('chain', 1:nchains)
    }
    for(i in 1:nchains) {
        if(nimbleOptions('verbose')) message('running chain ', i, '...')
        if(setSeed) set.seed(i)
        if(!missing(inits)) {
            if(is.function(inits)) {
                theseInits <- inits()
            } else if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]])) {
                theseInits <- inits[[i]]
            } else theseInits <- inits
            model$setInits(theseInits)
        }
        model$calculate()
        mcmc$run(niter, progressBar = progressBar)
        samples <- as.matrix(mcmc$mvSamples)
        if(nburnin > 0) samples <- samples[-(1:nburnin), , drop = FALSE]
        samplesList[[i]] <- samples
        if(returnWAIC) waic[i] <- mcmc$calculateWAIC(nburnin = nburnin)
    }
    if(samplesAsCodaMCMC) samplesList <- coda::as.mcmc.list(lapply(samplesList, coda::as.mcmc))
    if(nchains == 1) samplesList <- samplesList[[1]]  ## returns matrix when nchains = 1
    if(returnSummary) {
        if(nchains == 1) {
            summary <- samplesSummary(samplesList)
        } else {
            summary <- lapply(samplesList, samplesSummary)
            names(summary) <- paste0('chain', 1:nchains)
            summary$all.chains <- samplesSummary(do.call('rbind', samplesList))
        }
    }
    retList <- list()
    if(returnSamples) retList$samples <- samplesList
    if(returnSummary) retList$summary <- summary
    if(returnWAIC)    retList$waic    <- waic
    if(returnSamples + returnSummary + returnWAIC == 1) retList <- retList[[1]]
    return(retList)
}


#' Executes one or more chains of NIMBLE's default MCMC algorithm, for a model specified using BUGS code
#'
#' \code{nimbleMCMC} is designed as the most straight forward entry point to using NIMBLE's default MCMC algorithm.  It provides capability for running multiple MCMC chains, specifying the number of MCMC iterations, thinning, and burn-in, and which model variables should be monitored.  It also provides options to return the posterior samples, to return summary statistics calculated from the posterior samples, and to return the WAIC value calculated from each chain.
#'
#' The entry point for this function is providing the \cd{code}, \cd{constants}, \cd{data} and \cd{inits} arguments, to create a new NIMBLE model object, or alternatively providing an exisiting NIMBLE model object as the \cd{model} argument.
#'
#' @param code The quoted code expression representing the model, such as the return value from a call to \code{nimbleCode}).  No default value, this is a required argument.
#' 
#' @param constants Named list of constants in the model.  Constants cannot be subsequently modified. For compatibility with JAGS and BUGS, one can include data values with constants and \code{nimbleModel} will automatically distinguish them based on what appears on the left-hand side of expressions in \code{code}.
#' 
#' @param data Named list of values for the data nodes.  Data values can be subsequently modified.  Providing this argument also flags nodes as having data for purposes of algorithms that inspect model structure. Values that are NA will not be flagged as data.
#'
#' @param inits Argument to specify initial values for the model object, and for each MCMC chain.  See details.
#'
#' @param model A compiled or uncompiled NIMBLE model object.  When provided, this model will be used to configure the MCMC algorithm to be executed, rather than using the \cd{code}, \cd{constants}, \cd{data} and \cd{inits} arguments to create a new model object.  However, if also provided, the \cd{inits} argument will still be used to initialize this model prior to running each MCMC chain.
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
#' @param setSeed Logical argument.  If \code{TRUE}, then R's random number seed is set to a fixed value at the onset of each MCMC chain, which allows for reproducible results.  Default value is \code{FALSE}.
#'
#' @param progressBar Logical argument.  If \code{TRUE}, an MCMC progress bar is displayed during execution of each MCMC chain.  Default value is \code{TRUE}.
#'
#' @param returnSamples.  Logical argument.  If \code{TRUE}, then posterior samples are returned from each MCMC chain.  These samples are optionally returned as \code{coda} \code{mcmc} objects, depending on the \code{samplesAsCodaMCMC} argument.  Default value is \code{TRUE}.  See details.
#'
#' @param samplesAsCodaMCMC Logical argument.  If \code{TRUE}, then a \code{coda} \code{mcmc} object is returned instead of an R matrix of samples, or when \code{nchains > 1} a \code{coda} \code{mcmc.list} object is returned containing \code{nchains} \code{mcmc} objects.  This argument is only used when \code{returnSamples} is \code{TRUE}.  Default value is \code{FALSE}.  See details.
#' 
#' @param returnSummary Logical argument.  When \code{TRUE}, summary statistics for the posterior samples of each parameter are also returned, for each MCMC chain.  This may be returned in addition to the posterior samples themselves.  Default value is \code{FALSE}.  See details.
#'
#' @param returnWAIC Logical argument.  When \code{TRUE}, the WAIC (Watanabe, 2010) of the model is calculated and returned.  If multiple chains are run, then WAIC is calculated separately for each chain.  Default value is \code{FALSE}.  See details.
#'
#' @return A list is returned with named elements depending on the arguments passed to \code{nimbleMCMC}, unless only one among samples, summary, and WAIC are requested, in which case only that element is returned.  These elements may include \code{samples}, \code{summary}, and \code{waic}.  When \code{nchains = 1}, posterior samples are returned as a single matrix, and summary statistics as a single matrix.  When \code{nchains > 1}, posterior samples are returned as a list of matrices, one matrix for each chain, and summary statistics are returned as a list containing \code{nchains+1} matrices: one matrix corresponding to each chain, and the final element providing a summary of all chains, combined.  If \code{samplesAsCodaMCMC} is \code{TRUE}, then posterior samples are provided as \code{coda} \code{mcmc} and \code{mcmc.list} objects.  When \code{returnWAIC} is \code{TRUE}, WAIC values are returned as a numeric vector with one element corresponding to each chain.
#'
#' @details
#'
#' At least one of \code{returnSamples}, \code{returnSummary} or \code{returnWAIC} must be \code{TRUE}, since otherwise, nothing will be returned.  Any combination of these may be \cd{TRUE}, including possibly all three, in which case posterior samples, summary statistics, and WAIC values are returned for each MCMC chain.
#'
#' When \code{returnSamples = TRUE}, the form of the posterior samples is determined by the \code{samplesAsCodaMCMC} argument, as either matrices of posterior samples, or \cd{coda} \cd{mcmc} and \cd{mcmc.list} objects.
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
#'                           returnSummary = TRUE, returnWAIC = TRUE)
#' }
#'
#' @seealso \code{\link{configureMCMC}} \code{\link{buildMCMC}} \code{\link{runMCMC}}
#' 
#' @author Daniel Turek
#' 
#' @export
nimbleMCMC <- function(code, constants = list(), data = list(), inits, model,
                       monitors, thin = 1, niter = 10000, nburnin = 0, nchains = 1,
                       check = TRUE, setSeed = FALSE, progressBar = TRUE,
                       returnSamples = TRUE, samplesAsCodaMCMC = FALSE, returnSummary = FALSE, returnWAIC = FALSE) {
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
    if(!returnSamples && !returnSummary && !returnWAIC) stop('no output specified, use returnSamples = TRUE, returnSummary = TRUE, or returnWAIC = TRUE')
    if(samplesAsCodaMCMC) require(coda)
    if(missing(model)) {  ## model object not provided
        if(!missing(inits)) {
            if(!is.function(inits) && !is.list(inits)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
            if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]]) && (length(inits) != nchains)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
            if(is.function(inits)) {
                if(setSeed) set.seed(0)
                theseInits <- inits()
            } else if(is.list(inits) && (length(inits) > 0) && is.list(inits[[1]])) {
                theseInits <- inits[[1]]
            } else theseInits <- inits
            Rmodel <- nimbleModel(code, constants, data, theseInits, check = check)    ## inits provided
        } else Rmodel <- nimbleModel(code, constants, data, check = check)             ## inits not provided
    } else {              ## model object provided
        if(!is.model(model)) stop('model argument must be a NIMBLE model object')
        Rmodel <- if(is.Rmodel(model)) model else model$Rmodel
        if(!is.Rmodel(Rmodel)) stop('something went wrong')
    }
    conf <- configureMCMC(Rmodel, monitors = monitors, thin = thin)
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(Rmodel, Rmcmc)    ## only one compileNimble() call
    Cmcmc <- compiledList$Rmcmc
    nburnin <- ceiling(nburnin/thin)    ## accounts for thinning *first*, then dropping burnin
    runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits,
            setSeed = setSeed, progressBar = progressBar, returnSamples = returnSamples,
            samplesAsCodaMCMC = samplesAsCodaMCMC, returnSummary = returnSummary, returnWAIC = returnWAIC)
}


