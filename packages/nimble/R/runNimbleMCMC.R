
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

