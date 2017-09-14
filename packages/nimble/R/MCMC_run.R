#' Run one or more chains of an MCMC algorithm and extract samples
#'
#' Takes as input an MCMC algorithm (ideally a compiled one for speed)
#' and runs the MCMC with one or more chains, automatically extracting
#' the samples.
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
#' @param silent Logical argument.  If \code{TRUE}, then all output is suppressed during execution.  This overrides the \code{progressBar} argument (default = \code{FALSE}).
#'
#' @param returnCodaMCMC Logical argument.  If \code{TRUE}, then a \code{coda} \code{mcmc} object is returned instead of an R matrix of samples, or when \code{nchains > 1} a \code{coda} \code{mcmc.list} object is returned containing \code{nchains} \code{mcmc} objects (default = \code{FALSE}).
#'
#' @return When \code{nchains = 1}, a matrix of MCMC samples.  When \code{nchains > 1}, a list of length \code{nchains}, where each list element is a matrix of MCMC samples.  If \code{returnCodaMCMC = TRUE}, then a \code{coda} \code{mcmc} or \code{mcmc.list} object is returned instead.
#'
#' @details
#'
#' The \code{mcmc} argument can be a compiled or uncompiled NIMBLE MCMC algorithm, which is generated using \code{buildMCMC}.  Using a compiled algorithm will give substantially faster execution.
#'
#' If provided, the \code{inits} argument can be one of three things:
#' (1) a function to generate initial values, which will be executed to generate initial values at the beginning of each MCMC chain,
#' (2) a single named list of initial values which, will be used for each chain, or
#' (3) a list of length \code{nchains}, each element being a named list of initial values which be used for one MCMC chain.
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
#' @seealso \code{\link{configureMCMC}} \code{\link{buildMCMC}}
#'
#' @author Daniel Turek
#'
#' @export
runMCMC <- function(mcmc,
                    niter = 10000, nburnin = 0, nchains = 1,
                    inits,
                    setSeed = FALSE,
                    progressBar = TRUE,
                    silent = FALSE,
                    returnCodaMCMC = FALSE) {
    if(missing(mcmc)) stop('must provide a NIMBLE MCMC algorithm')
    if(!identical(nf_getGeneratorFunction(mcmc), buildMCMC)) stop('mcmc argument must be a NIMBLE MCMC algorithm')
    if(!is.Cnf(mcmc)) message('Warning: running an uncompiled MCMC algorithm, use compileNimble() for faster execution.')
    if(nchains < 1) stop('must have nchains > 0')
    if(!missing(inits)) {
        if(!is.function(inits) && !is.list(inits)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
        if(is.list(inits) && is.list(inits[[1]]) && (length(inits) != nchains)) stop('inits must be a function, a list of initial values, or a list (of length nchains) of lists of inital values')
    }
    model <- if(is.Cnf(mcmc)) mcmc$Robject$model$CobjectInterface else mcmc$model
    if(!inherits(model, 'modelBaseClass')) stop('something went wrong')
    samplesList <- vector('list', nchains)
    for(i in 1:nchains) {
        if(!silent) message('running chain ', i, '...')
        if(setSeed) set.seed(i)
        if(!missing(inits)) {
            if(is.function(inits)) {
                theseInits <- inits()
            } else if(is.list(inits[[1]])) {
                theseInits <- inits[[i]]
            } else theseInits <- inits
            model$setInits(theseInits)
        }
        model$calculate()
        mcmc$run(niter, progressBar = progressBar && !silent)
        samples <- as.matrix(mcmc$mvSamples)
        if(nburnin > 0) samples <- samples[-(1:nburnin), , drop = FALSE]
        samplesList[[i]] <- samples
    }
    if(returnCodaMCMC) samplesList <- coda::as.mcmc.list(lapply(samplesList, coda::as.mcmc))
    if(nchains == 1) samplesList <- samplesList[[1]]  ## returns matrix when nchains=1
    return(samplesList)
}



#' Builds R model and compiles C model, configures and compiles MCMC,
#' and can optionally run the MCMC or return the model and MCMC
#' components
#'
#'
#' @details
#'
#' Creates and optionally runs a BUGS model.
#' By default, this will run the provided model, record all samples
#' and generate summary statistics.  This default behavior can be
#' altered via a variety of arguments.  Following execution of the
#' MCMC algorithms, returns a named list containing \code{samples} and
#' \code{summary}. If runMCMC = FALSE, returns the model components
#'     required to run the model outside the function.
#'
#' @param code The quoted code expression representing the model, such
#'     as the return value from a call to \code{nimbleCode}).  No
#'     default value, this is a required argument.
#'
#' @param constants A named list giving values of constants for the
#'     model.  This is the same as the \code{constants} argument which
#'     would be passed to \code{nimbleModel}.  Default value is
#'     \code{list()}.
#'
#' @param data A named list giving the data values for the model.
#'     This is the same as the \code{data} argument which would be
#'     passed to \code{nimbleModel} or \code{model$setData}.  Default
#'     value is \code{list()}.
#'
#' @param inits A named list giving the initial values for the model.
#'     This is the same as the \code{inits} argument which would be
#'     passed to \code{nimbleModel} or \code{model$setInits}.  Default
#'     value is \code{list()}.
#'
#' @param monitors A character vector giving the node names or
#'     variable names to monitor.  The samples corresponding to these
#'     nodes will be stored in the output samples, will have summary
#'     statistics calculated, and density and trace plots generated.
#'     Default value is all top-level stochastic nodes of the model.
#'
#' @param niter Number of MCMC iterations to run. Default value is 10000.
#'
#' @param burnin Number of initial, post-thinning, MCMC iterations to discard.
#' Default value is 100.
#'
#' #' @param nchains Number MCMC chains to run.
#' Default value is 3.
#'
#' @param thin Thinning interval for the MCMC samples.  This applies
#'     to all MCMC algorithms in the suite.  The thinning occurs prior
#'     to the burnin samples being discarded.  Default value is 1.
#'
#' @param check Logical argument, specifying whether to check the
#'     model object for missing or invalid values.  Default value is
#'     FALSE.
#'
#' @param debug Logical argument, specifying whether to enter a
#'     \code{browser()} at the onset of executing each MCMC algrithm.
#'     For use in debugging individual MCMC algorithms, if necessary.
#'     Default value is FALSE.
#'
#' @param runMCMC Logical argument, specifying whether to run the MCMC
#'     or return the setup model objects Default value is FALSE.
#'
#' @return Returns a named list containing elements: samples: A matrix
#'     containing samples from each MCMC algorithm.  summary: A matrix
#'     containing summary statistics for each variable and algorithm.
#'     If runMCMC = FALSE, returns a named list containing elements:
#'     R.model: NIMBLE model created by a call to \code{nimbleModel()}
#'     C.mcmc: compiled NIMBLE model and MCMC created by a call to
#'     \code{compileNimble()} mcmc.spec: defaut MCMC configuration for
#'     a given model created by a call to \code{configureMCMC()}
#'
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' model.samples <- runNimbleMCMC(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     monitors = 'mu')
#'
#' model.components <- runNimbleMCMC(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     runMCMC = FALSE,
#'                     monitors = 'mu')
#'
#'
#' }
#'
#' @author L. Ponisio
#' @export

runNimbleMCMC <- function(code,
                          constants=list(),
                          data=list(),
                          inits=list(),
                          monitors=NULL,
                          niter = 10000,
                          thin=1,
                          burnin=100,
                          nchains=3,
                          check = FALSE,
                          debug = FALSE,
                          runMCMC = TRUE, ...){

    ## Build the R model
    R.model <- nimbleModel(code = code,
                           constants = constants,
                           data = data,
                           inits = inits,
                           check = check,
                           debug = debug,...)
    message('**** R model created ****')
    if(is.null(monitors)){
        monitors <- R.model$getNodeNames(topOnly=TRUE)
    }
    ## Configure and build mcmc
    mcmc.spec <- configureMCMC(R.model,
                               print = FALSE,
                               monitors = monitors,
                               thin = thin)
    mcmc <- buildMCMC(mcmc.spec)
    message('**** MCMC built ****')

    ## Compile model in C++
    C.model <- compileNimble(R.model)
    C.mcmc <- compileNimble(mcmc, project = R.model)
    message('**** NIMBLE model compiled ****')

    if(runMCMC){
        ## Runs the MCMC, returns the samples and summary statistics.
        ## Run the model
        message('**** Running model ****')
        sample.list <- runMCMC(C.mcmc,
                               niter = niter,
                               nburnin = burnin,
                               nchains = nchains)

        calcSummaryStats <- function(samples){
            ## Calculates mean, median, sd, and CI
            means <- mean(samples)
            sds <- sd(samples)
            return(as.matrix(c(means,
                               median(samples), sds,
                               means - sds*1.96, means + sds*1.96)))
        }
        ## If the number of chains is > 1 takes calculates the summary
        ## statistics across all of the chains.
        ##
        ## Extract the MCMC for each parameter from the different
        ## chains, the caclulate the summary statistics.
        getParSamples <- function(samples, par.name){
            samples[, par.name]
        }
        summary.stats.all.chains <- sapply(monitors, function(par.name){
            all.par.samples <- do.call(c,
                                       lapply(sample.list, getParSamples,
                                              par.name))
            summary.stats <- calcSummaryStats(all.par.samples)
            colnames(summary.stats) <- par.name
            return(summary.stats)
        })
        ## Ronames are lost with multi parameter models if not set here
        rownames(summary.stats.all.chains) <-
            c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp')

        return(list(samples=sample.list,
                    summary=summary.stats.all.chains))

    }else{
        ## Return the model components which can be used to run the
        ## model using runMCMC outside this function.
        model.components <- list(R.model = R.model,
                                 mcmc.spec = mcmc.spec,
                                 C.mcmc = C.mcmc)
        return(model.components)
    }
}

