
#' Builds R model and compiles C model, configures and compiles MCMC,
#' and can optionally run the MCMC or return all the model and MCMC
#' components
#'
#'
#' @details
#'
#' Creates and runs an BUGS model.
#' By default, this will run the provided model, record all samples
#' and generates summary statistics.  This default behavior can be
#' altered via a variety of arguments.  Following execution of the
#' MCMC algorithms, returns a named list containing \code{samples} and
#' \code{summary}
#'
#' @param code The quoted code expression representing the model, such
#'     as the return value from a call to \code{nimbleCode}).  No
#'     default value, this is a required argument.
#'
#' @param constants A named list giving values of constants for the
#'     model.  This is the same as the \code{constants} argument which
#'     would be passed to \code{nimbleModel}.  Default value is
#'     list().
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
#' @param niter Number of MCMC iterations to run.  This applies to all
#'     MCMC algorithms in the suite.  Default value is 1000.
#'
#' @param burnin Number of initial, post-thinning, MCMC iterations to discard.
#' Default value is 10.
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
#' @return Returns a named list containing elements:
#' samples: A matrix containing samples from each MCMC algorithm.
#' summary: A matrix containing summary statistics for each variable and algorithm.
#'
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' output <- runNimbleMCMC(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     niter = 10000,
#'                     monitors = 'mu')
#' }
#'
#' @author L. Ponisio
#' @export

runNimbleMCMC <- function(code, ## model created by nimbleCode()
                          constants,
                          data,
                          inits, ## inital values for parameters to estimate
                          monitors, ## vector of parameters to monitor
                          niter = 1000, ## number of samples
                          thin=1,
                          burin=10,
                          nchains=3,
                          check = FALSE,
                          debug = FALSE,
                          runMCMC = TRUE, ...){ ## thinning rate
    ## build model
    R.model <- nimbleModel(code = code,
                           constants = constants,
                           data = data,
                           inits = inits,
                           check = check,
                           debug = debug,...)
    message('R model created')

    ## configure and build mcmc
    mcmc.spec <- configureMCMC(R.model,
                               print = FALSE,
                               monitors = monitors,
                               thin = thin)
    mcmc <- buildMCMC(mcmc.spec)
    message('MCMC built')

    ## compile model in C++
    C.model <- compileNimble(R.model)
    C.mcmc <- compileNimble(mcmc, project = R.model)
    message('NIMBLE model compiled')

    if(runMCMC){
        ## run model
        message('running model')
        samplesList <- runMCMC(C.mcmc, niter = niter, nburnin = burnin, nchains = nchains)
    }
    return(C.mcmc)
}
