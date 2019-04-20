#' Executes multiple MCMC algorithms and organizes results.
#'
#' Creates, runs, and organizes output from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
#' This can include WinBUGS, OpenBUGS, JAGS and Stan MCMCs, as well as NIMBLE MCMC algorithms.
#' Trace plots and density plots for the MCMC samples may also be generated and saved.
#'
#' @details
#' Creates and runs an MCMC Suite.
#' By default, this will execute the specified MCMCs, record all samples, generate summary statistics, and create and save trace plots and posterior density plots.
#' This default behavior can ben altered via a variety of arguments.
#' Following execution of the MCMC algorithms, returns a named list containing \code{samples}, \code{summary}, and \code{timing} elements.
#' See the NIMBLE User Manual for more information about the organization of the return object.
#' 
#' @param code The quoted code expression representing the model, such as the return value from a call to \code{nimbleCode}).
#' No default value, this is a required argument.
#'
#' @param constants A named list giving values of constants for the model.
#' This is the same as the \code{constants} argument which would be passed to \code{nimbleModel}.
#' Default value is list().
#'
#' @param data A named list giving the data values for the model.
#' This is the same as the \code{data} argument which would be passed to \code{nimbleModel} or \code{model$setData}.
#' Default value is \code{list()}.
#' 
#' @param inits A named list giving the initial values for the model.
#' This is the same as the \code{inits} argument which would be passed to \code{nimbleModel} or \code{model$setInits}.
#' Default value is \code{list()}.
#'
#' @param monitors A character vector giving the node names or variable names to monitor.
#' The samples corresponding to these nodes will be stored in the output samples, will have summary statistics calculated, and density and trace plots generated.
#' Default value is all top-level stochastic nodes of the model.
#' 
#' @param niter Number of MCMC iterations to run.
#' This applies to all MCMC algorithms in the suite.
#' Default value is 10,000.
#'
#' @param burnin Number of initial, post-thinning, MCMC iterations to discard.
#' Default value is 2,000.
#' 
#' @param thin Thinning interval for the MCMC samples.
#' This applies to all MCMC algorithms in the suite.  The thinning occurs prior to the burnin samples being discarded.
#' Default value is 1.
#' 
#' @param summaryStats A character vector, specifying the summary statistics to calculate on the MCMC samples.
#' Each element may be the character name of an exisiting R function (possibly user-defined) which acts on a numeric vector and returns a scalar (e.g., \code{mean} or \code{sd},
#' or a character string which when parsed and evaluted will define such a function (e.g., \code{function(x) mean(sqrt(x))}).
#' Default value is \code{c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp')}, where the final two elements are functions which calculate the limits of a 95 percent Bayesian credible interval.
#' 
#' @param calculateEfficiency A logical, specifying whether to calculate the efficiency for each MCMC algorithm.  Efficiency is defined as the effective sample size (ESS) of each model parameter divided by the algorithm runtime (in seconds).  Default is FALSE.
#'
#' @param MCMCs A character vector specifying the MCMC algorithms to run.
#' \code{'winbugs'} specifies WinBUGS;
#' \code{'openbugs'} specifies OpenBUGS;
#' \code{'jags'} specifies JAGS;
#' \code{'stan'} specifies Stan; in this case, must also provide the \code{'stan_model'} argument;
#' \code{'nimble'} specifies NIMBLE's default MCMC algorithm;
#' \code{'nimble_noConj'} specifies NIMBLE's default MCMC algorithm without the use of any conjugate Gibbs sampling;
#' \code{'nimble_RW'} specifies NIMBLE MCMC algorithm using only random walk Metropolis-Hastings (\code{'RW'}) samplers;
#' \code{'nimble_slice'} specifies NIMBLE MCMC algorithm using only slice (\code{'slice'}) samplers;
#' \code{'autoBlock'} specifies NIMBLE MCMC algorithm with block sampling of dynamically determined parameter groups attempting to maximize sampling efficiency;
#' Anything else will be interpreted as NIMBLE MCMC algorithms, and must have associated entries in the MCMCdefs argument.
#' Default value is \code{'nimble'}, which specifies NIMBLE's default MCMC algorithm.
#' 
#' @param MCMCdefs A named list of MCMC definitions.  The names of list elements should corespond to any custom MCMC algorithms specified in the \code{MCMCs} argument.
#' The list elements should be quoted expressions, enclosed in {} braces.  When executed, the internal code must return an MCMC configuration object, 
#' specifying the corresponding MCMC algorithm; in particular, setting the appropriate samplers.  The code may assume existance of the R model object \code{Rmodel},
#' and must *return* the MCMC configuration object.  Therefore, the final line of such a code block would frequently be a standalone \code{MCMCconf}, to return this object.
#' 
#' @param winbugs_directory A character string giving the directory of the executable WinBUGS program for the WinBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'C:/WinBUGS14'}.
#' 
#' @param winbugs_program A character string giving the name of the WinBUGS program, for the WinBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'WinBUGS'}.
#'
#' @param openbugs_directory A character string giving the directory of the executable OpenBUGS program for the OpenBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'C:/OpenBUGS323'}.
#' 
#' @param openbugs_program A character string giving the name of the OpenBUGS program, for the OpenBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'OpenBUGS'}.
#' 
#' @param stan_model A character string specifying the location and name of the model file (\code{'modelName.stan'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.stan'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' 
#' @param stan_inits A character string specifying the location and name of the inits file (\code{'modelName.init.R'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.init.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' If omitted, it will attempt to locate an inits file in the same directory as the Stan model file.
#'
#' @param stan_data A character string specifying the location and name of the data file (in the form \code{'modelName.data.R'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.data.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' If omitted, it will attempt to locate a data file in the same directory as the Stan model file.
#'
#' @param stanNameMaps A list specifying name mappings between Stan and WinBUGS/OpenBUGS.
#' The syntax for list elements is list(BUGS_PARAM_NAME = list(StanSourceName = 'STAN_PARAM_NAME', transform = function(x) TRANSFORMATION_FUNCTION(x))).
#' The transformation is optional.
#' 
#' @param makePlot Logical argument, specifying whether to generate the trace plots and posterior density plots, for each monitored node.
#' Default value is \code{TRUE}.
#' 
#' @param savePlot Logical argument, specifying whether to save the trace plots and density plots.
#' Plots will be saved into the current working directory.
#' Only used when \code{makePlot == TRUE}.
#' Default value is \code{TRUE}.
#' 
#' @param plotName Character string, giving the file name for saving the trace plots and density plots.
#' Only used when \code{makePlot == TRUE} and \code{savePlot == TRUE}.
#' Default value is \code{'MCMCsuite'}.
#'
#' @param setSeed Logical argument, specifying whether to set.seed(0) prior to MCMC sampling.
#' Default value is \code{TRUE}.
#' 
#' @param check Logical argument, specifying whether to check the model object for missing or invalid values.  Default is given by the NIMBLE option 'checkModel', see help on \code{nimbleOptions} for details.
#' 
#' @param debug Logical argument, specifying whether to enter a \code{browser()} at the onset of executing each MCMC algrithm.
#' For use in debugging individual MCMC algorithms, if necessary.
#' Default value is FALSE.
#'
#' @param ... For internal use only
#'
#' @return Returns a named list containing elements:
#' samples: A 3-dimensional array containing samples from each MCMC algorithm.
#' summary: A 3-dimensional array containing summary statistics for each variable and algorithm.
#' timing: A numeric vector containing timing information.
#' efficiency: Minimum and mean sampling efficiencies for each algorithm (only provided if option calculateEfficiency = TRUE).
#' See the NIMBLE User Manual for more information about the organization of the return object.
#'  
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' output <- MCMCsuite(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     niter = 10000,
#'                     monitors = 'mu',
#'                     MCMCs = c('nimble', 'nimble_RW'),
#'                     summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
#'                     makePlot = FALSE)
#' }
#' 
#' @author Daniel Turek
#' @export
MCMCsuite <- function(
                      code,
            constants           = list(),
            data                = list(),
            inits               = list(),
            monitors            = character(),
            niter               = 10000,
            burnin              = 2000,
            thin                = 1,
            summaryStats        = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp'),
            calculateEfficiency = FALSE,
            MCMCs               = 'nimble',
            MCMCdefs            = list(),
            winbugs_directory   = 'C:/WinBUGS14',
            winbugs_program     = 'WinBUGS',
            openbugs_directory  = 'C:/OpenBUGS323',
            openbugs_program    = 'OpenBUGS',
            stan_model          = '',
            stan_inits          = NULL,
            stan_data           = NULL,
            stanNameMaps        = list(),
            makePlot            = TRUE,
            savePlot            = TRUE,
            plotName            = 'MCMCsuite',
            setSeed             = TRUE,
            check               = getNimbleOption('checkModel'),
            debug               = FALSE) {
    ## aliased in MCMCsuiteClass
    suite <- MCMCsuiteClass(code, constants, data, inits, monitors, niter, burnin, thin, summaryStats, calculateEfficiency,
                            MCMCs, MCMCdefs, winbugs_directory, winbugs_program, openbugs_directory, openbugs_program,
                            stan_model, stan_inits, stan_data, stanNameMaps, makePlot, savePlot, plotName, setSeed,
                            check, debug)
    return(suite$output)
}

#' Class \code{MCMCsuiteClass}
#'
#' @aliases MCMCsuiteClass-class
#'
#' @description
#' Objects of this class create, run, and organize output from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
#' This can include WinBUGS, OpenBUGS, JAGS and Stan MCMCs, as well as NIMBLE MCMC algorithms.
#' Trace plots and density plots for the MCMC samples may also be generated and saved.
#'
#' @seealso \code{\link{MCMCsuite}}
#' 
#' @author Daniel Turek
#' @export
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' output <- MCMCsuite(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     niter = 10000,
#'                     monitors = 'mu',
#'                     MCMCs = c('nimble', 'nimble_RW'),
#'                     summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
#'                     makePlot = FALSE)
#' }
#' 
MCMCsuiteClass <- setRefClass(

    Class = 'MCMCsuiteClass',
    
    fields = list(
        ## set in initialize()
        code = 'ANY',   ## parsed expression for the model code; must be contained in { ... }    --- ORIGINAL ARGUMENT
        constants = 'list',   ## list of the constants (ORIGINAL ARGUMENT)
        data = 'list',   ## list of the data    --- ORIGINAL ARGUMENT
        inits = 'list',  ## named list of initial values used for all MCMC algorithms    --- ORIGINAL ARGUMENT
        constantsAndData = 'list',   ## data list used for WinBUGS, OpenBUGS, JAGS.  is equal to c(constantList, dataList)
        Rmodel = 'ANY',   ## Rmodel object
        
        ## setMonitors()
        monitors = 'character',    ## the original character vector argument to initialize()    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        monitorVars = 'character',    ## character vector of VARIABLE names of parameters to save
        monitorNodesNIMBLE = 'character',  ## character vector of the monitor node names, with spaces as in nimble: 'y[1, 1]'
        monitorNodesBUGS = 'character',    ## same as monitorNodes, except for WinBUGS and OpenBUGS: no spaces in node names: 'y[1,1]'
        nMonitorNodes = 'numeric',   ## number of monitorNodes
        
        ## set in initialize()
        niter = 'numeric',    ## number of MCMC iterations to run    --- ORIGINAL ARGUMENT
        burnin = 'numeric',   ## burn-in period, the number of initial samples to discard, prior to thinning    --- ORIGINAL ARGUMENT
        thin = 'numeric',   ## thinning interval    --- ORIGINAL ARGUMENT
        nkeep = 'numeric',   ## number of samples we'll keep. equal to (niter/thin - burnin)
        burninFraction = 'numeric',  ## fraction of total sampling effort spent on burnin (burnin / (nkeep + burnin))
        
        ## setSummaryStats()
        summaryStats = 'character',    ## character vector of parseable summary statistic functions    --- ORIGINAL ARGUMENT
        calculateEfficiency = 'logical',   ## logical specifying whether to calculate ESS and Efficiency    --- ORIGINAL ARGUMENT
        summaryStatFunctions = 'list',  ## list of the function objects for summary statistics
        summaryStatDimNames = 'character',   ## character vector of the dimension names in output$summary
        nSummaryStats = 'numeric',   ## the number of summary statistics
        
        ## setMCMCs()
        MCMCs = 'character',   ## character vector of the MCMC analyses.  'winbugs', 'openbugs', 'jags', 'stan', or anything else is nimble    --- ORIGINAL ARGUMENT
        winbugsMCMCflag = 'logical',   ## whether 'winbugs' is in MCMCs
        openbugsMCMCflag = 'logical',   ## whether 'openbugs' is in MCMCs
        jagsMCMCflag = 'logical',   ## whether 'jags' is in MCMCs
        stanMCMCflag = 'logical',   ## whether 'stan' is in MCMCs
        nimbleMCMCs = 'character',    ## the names of the remaining (presumably nimble) MCMCs
        nNimbleMCMCs = 'numeric',    ## the number of remaining (nimble) MCMCs
        nimbleMCMCflag = 'logical',   ## a flag indicating whether there are any remaining (nimble) MCMCs
        nMCMCs = 'numeric',   ## the number of MCMC algorithms being ran
        
        ## setMCMCdefs()
        MCMCdefs = 'list',   ## named list of {} expression code blocks, corresponding the setup for nimble MCMCs    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        MCMCdefNames = 'character',   ## names of the MCMCdefs list
        
        ## set in initialize()
        winbugs_directory = 'character',    ## directory for WinBUGS program    --- ORIGINAL ARGUMENT
        winbugs_program = 'character',     ## program for WinBUGS    --- ORIGINAL ARGUMENT
        openbugs_directory = 'character',    ## directory for OpenBUGS program    --- ORIGINAL ARGUMENT
        openbugs_program = 'character',     ## program for OpenBUGS    --- ORIGINAL ARGUMENT
        stan_model = 'character',     ## *.stan model file    --- ORIGINAL ARGUMENT
        makePlot = 'logical',    ## whether to generate plots    --- ORIGINAL ARGUMENT
        savePlot = 'logical',   ## whether or not to save plot PDFs    --- ORIGINAL ARGUMENT
        plotName = 'character',     ## name of the file where we save density and trace plots    --- ORIGINAL ARGUMENT
        setSeed = 'logical',   ## whether to setSeed(0) prior to running each algorithm    --- ORIGINAL ARGUMENT
        debug = 'logical',   ## whether to enter browser() before running each algorithm    --- ORIGINAL ARGUMENT
        modelFileName = 'character',     ## name of the text file where we write the model code, set to a fixed value

        ## Maps with possible transformations from Stan to WinBUGS/OpenBUGS
        ## e.g. for blocker: StanNameMaps <- list(tau = list(StanSourceName = 'sigmasq_delta', transform = function(x) 1/x)) ## transform can be omitted
        StanNameMaps = 'ANY',
        
        ## set in run()
        Cmodel = 'ANY',   ## compiled Cmodel object
        RmcmcFunctionList = 'list',    ## list of the R (nimble) MCMC functions
        CmcmcFunctionList = 'list',    ## list of the C (nimble) MCMC functions
        output = 'list'   ## list of numeric outputs: samples, summary, timing
    ),
    
    methods = list(
        
        initialize = function(
            code,
            constants           = list(),
            data                = list(),
            inits               = list(),
            monitors            = character(),
            niter               = 10000,
            burnin              = 2000,
            thin                = 1,
            summaryStats        = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp'),
            calculateEfficiency = FALSE,
            MCMCs               = 'nimble',
            MCMCdefs            = list(),
            winbugs_directory   = 'C:/WinBUGS14',
            winbugs_program     = 'WinBUGS',
            openbugs_directory  = 'C:/OpenBUGS323',
            openbugs_program    = 'OpenBUGS',
            stan_model          = '',
            stan_inits          = NULL,
            stan_data           = NULL,
            stanNameMaps        = list(),
            makePlot            = TRUE,
            savePlot            = TRUE,
            plotName            = 'MCMCsuite',
            setSeed             = TRUE,
            check               = getNimbleOption('checkModel'),
            debug               = FALSE) {
            
            if(debug) browser()
            code <<- code
            constants <<- constants
            data <<- data
            inits <<- inits
            constantsAndData <<- c(constants, data)
            Rmodel <<- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=check)
            niter <<- niter
            burnin <<- burnin
            thin <<- thin
            nkeep <<- floor(niter/thin) - burnin
            if(nkeep < 0) stop('niter/thin - burnin is negative; this would not retain any samples; try increasing niter, or decreasing burnin')
            burninFraction <<- burnin / (nkeep + burnin)
            setMonitors(monitors)
            setSummaryStats(summaryStats, calculateEfficiency)
            setMCMCs(MCMCs)
            setMCMCdefs(MCMCdefs)
            winbugs_directory <<- winbugs_directory
            winbugs_program <<- winbugs_program
            openbugs_directory <<- openbugs_directory
            openbugs_program <<- openbugs_program
            stan_model <<- stan_model
            if(is.null(stan_inits)) stan_inits <- gsub('stan$', 'init.R', stan_model)
            if(is.null(stan_data))  stan_data  <- gsub('stan$', 'data.R', stan_model)
            StanNameMaps <<- stanNameMaps
            makePlot <<- makePlot
            savePlot <<- savePlot
            plotName <<- plotName
            setSeed <<- setSeed
            debug <<- debug
            modelFileName <<- 'model.txt'

            ## run
            checkMCMCdefNames()
            init_output()
            writeModelFile()
            if(debug)              browser()
            if(winbugsMCMCflag)    run_winbugs()
            if(openbugsMCMCflag)   run_openbugs()
            if(jagsMCMCflag)       run_jags()
            if(stanMCMCflag)       run_stan(stan_data, stan_inits)
            if(nimbleMCMCflag)     run_nimble()
            unlink(modelFileName)
            if(makePlot)           generate_plots()
        },
        
        setMonitors = function(newMonitors) {
            if(length(newMonitors) == 0) newMonitors <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
            newMonitors <- Rmodel$expandNodeNames(newMonitors, returnScalarComponents = TRUE)
            dataFlags <- unlist(lapply(newMonitors, function(mon) eval(parse(text=mon, keep.source=FALSE)[[1]], envir=Rmodel$isDataEnv)))
            newMonitors <- newMonitors[!dataFlags]
            monitors <<- newMonitors
            monitorVars <<- unique(removeIndexing(monitors))
            monitorNodesNIMBLE <<- monitors
            monitorNodesBUGS <<- gsub(' ', '', monitorNodesNIMBLE)
            nMonitorNodes <<- length(monitorNodesNIMBLE)
        },
        
        setSummaryStats = function(summaryStats_arg, calculateEfficiency) {
            calculateEfficiency <<- calculateEfficiency
            if(calculateEfficiency) {
                n <- length
                ess <- effectiveSize
                efficiency <- function(x) return(0)   ## placeholder; calculation done in addToOutput()
                summaryStats_arg <- c(summaryStats_arg, 'n', 'ess', 'efficiency')
            }
            summaryStats <<- summaryStats_arg
            CI95_low <- function(x) quantile(x, probs = 0.025)
            CI95_upp <- function(x) quantile(x, probs = 0.975)
            summaryStatFunctions <<- lapply(summaryStats, function(txt) eval(parse(text=txt)[[1]]))
            summaryStatDimNames <<- gsub('function *\\(.*?\\)', '', summaryStats)
            summaryStatDimNames <<- gsub('^ *', '', summaryStatDimNames)
            nSummaryStats <<- length(summaryStats)
        },
        
        setMCMCs = function(MCMCs) {
            MCMCs <<- unique(MCMCs)
            winbugsMCMCflag <<- 'winbugs' %in% MCMCs
            openbugsMCMCflag <<- 'openbugs' %in% MCMCs
            jagsMCMCflag <<- 'jags' %in% MCMCs
            stanMCMCflag <<- 'stan' %in% MCMCs
            nimbleMCMCs <<- setdiff(MCMCs, c('winbugs', 'openbugs', 'jags', 'stan'))
            nNimbleMCMCs <<- length(nimbleMCMCs)
            nimbleMCMCflag <<- if(nNimbleMCMCs > 0) TRUE else FALSE
            nMCMCs <<- length(MCMCs)
        },
        
        setMCMCdefs = function(newMCMCdefs) {
            MCMCdefs <<- list(nimble        = quote(configureMCMC(Rmodel)),
                              nimble_noConj = quote(configureMCMC(Rmodel, useConjugacy = FALSE)),
                              nimble_RW     = quote(configureMCMC(Rmodel, onlyRW       = TRUE)),
                              nimble_slice  = quote(configureMCMC(Rmodel, onlySlice    = TRUE)),
                              autoBlock     = quote(configureMCMC(Rmodel, autoBlock    = TRUE)))
            MCMCdefs[names(newMCMCdefs)] <<- newMCMCdefs
            MCMCdefNames <<- names(MCMCdefs)
        },
        
        init_output = function() {
            samples <- array(NA, dim = c(nMCMCs, nMonitorNodes, nkeep))
            dimnames(samples) <- list(MCMCs, monitorNodesNIMBLE, NULL)
            summary <- array(NA, dim = c(nMCMCs, nSummaryStats, nMonitorNodes))
            dimnames(summary) <- list(MCMCs, summaryStatDimNames, monitorNodesNIMBLE)
            timing <- rep(NA, nMCMCs+1)
            names(timing) <- c(MCMCs, 'nimble_compile')
            if(stanMCMCflag) timing['stan_compile'] <- NA
            runParams <- c(niter = niter, burnin = burnin, thin = thin, nkeep = nkeep, burninFraction = burninFraction) 
            initialOutput <- list(samples=samples, summary=summary, timing=timing, runParams = runParams)
            if(calculateEfficiency) initialOutput$efficiency <- list(min=NA, mean=NA)
            output <<- initialOutput
        },
        
        run_winbugs = function() {
            if(setSeed) set.seed(0)
            if(requireNamespace('R2WinBUGS', quietly = TRUE)) {
                timeResult <- system.time({
                                              winbugs_out <- R2WinBUGS::bugs(data=constantsAndData, inits=list(inits), parameters.to.save=monitorVars, model.file=modelFileName,
                                                                             n.chains=1, n.iter=niter, n.burnin=0, n.thin=thin, bugs.directory=winbugs_directory, program=winbugs_program)
                                          })
                tempArray <- winbugs_out$sims.array[, 1, ]        ## must use sims.array
                samplesArray <- tempArray[(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
                addToOutput('winbugs', samplesArray, timeResult)
            } else warning("run_winbugs: R2WinBUGS package is required for 'winbugs' option.")
        },
        
        run_openbugs = function() {
            if(requireNamespace('R2WinBUGS', quietly = TRUE)) {
                if(setSeed) set.seed(0)
                timeResult <- system.time({
                                              openbugs_out <- R2WinBUGS::bugs(data=constantsAndData, inits=list(inits), parameters.to.save=monitorVars, model.file=modelFileName,
                                                                              n.chains=1, n.iter=niter, n.burnin=0, n.thin=thin, bugs.directory=openbugs_directory, program=openbugs_program)
                                          })
                tempArray <- openbugs_out$sims.array[, 1, ]        ## must use sims.array
                samplesArray <- tempArray[(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
                addToOutput('openbugs', samplesArray, timeResult)
            } else warning("run_openbugs: R2WinBUGS package is required for 'openbugs' option.")
        },
        
        run_jags = function() {
            if(setSeed) set.seed(0)
            if(requireNamespace('rjags', quietly = TRUE)) {
                jags_mod <- rjags::jags.model(file=modelFileName, data=constantsAndData, inits=inits, n.chains=1, quiet=FALSE)
                timeResult <- system.time({
                                              jags_out <- rjags::coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=thin)
                                          })
                samplesArray <- jags_out[[1]][(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
                addToOutput('jags', samplesArray, timeResult)
            } else warning("run_jags: rjags package is required for 'jags' option.")
        },

        run_stan = function(dataFile, initFile) {
            if(setSeed) set.seed(0)
            if(require('rstan', quietly = TRUE)) {
                 ## warning("MCMCsuite: use of rstan is not yet provided via the CRAN version of NIMBLE because of packaging issues. To use this functionality, please install NIMBLE from <http://r-nimble.org>.")
                 if(stan_model == '') stop('must provide \'stan_model\' argument to run Stan MCMC')
                 ##            dataFile <- gsub('stan$', 'data.R', stan_model)
                 ##            initFile <- gsub('stan$', 'init.R', stan_model)
                 if(!is.list(dataFile)) 
                     constantsAndDataStan <- fileToList(dataFile)
                 else
                     constantsAndDataStan <- dataFile
                
                 if(!is.list(initFile)) {
                     if(file.exists(initFile))
                         initsStan <- fileToList(initFile)
                     else
                         initsStan <- NULL
                 } else
                     initsStan <- initFile

                
                 timeResult <- system.time(stan_mod <- rstan::stan_model(file = stan_model))
                 addTimeResult('stan_compile', timeResult)
                        
                 if(is.null(initsStan)) {
                     ## missing model.init.R file (stan inits file)
                     timeResult <- system.time(stan_out <- rstan::sampling(stan_mod, data=constantsAndDataStan, chains=1, iter=niter, thin=thin))
                 } else {
                       ## we have the model.init.R file
                       ## this one includes inits = ...
                       timeResult <- system.time(stan_out <- rstan::sampling(stan_mod, data=constantsAndDataStan, chains=1, iter=niter, thin=thin, init=list(initsStan)))
                   }
                
                 tempArray <- rstan::extract(stan_out, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
                 for(BUGSname in names(StanNameMaps)) {
                     iCol <- which(StanNameMaps[[BUGSname]]$StanSourceName == colnames(tempArray))
                     if(length(iCol)==1) {
                         if(!is.null(StanNameMaps[[BUGSname]]$transform))
                             tempArray[,iCol] <- StanNameMaps[[BUGSname]]$transform(tempArray[,iCol])
                         colnames(tempArray)[iCol] <- BUGSname
                     }
                 }
                 dimnames(tempArray)[[2]] <- gsub('_', '.', dimnames(tempArray)[[2]])
                 if(!all(monitorNodesBUGS %in% dimnames(tempArray)[[2]])) {
                     missingNames <- setdiff(monitorNodesBUGS, dimnames(tempArray)[[2]])
                     warning(paste0('Stan output is missing values for: ', paste0(missingNames,collapse=', ')))
                 }
                 samplesArray <- array(0, dim = c(nkeep, length(monitorNodesBUGS)))
                 dimnames(samplesArray)[[2]] <- monitorNodesBUGS
                 monitorsWeHave <- intersect(monitorNodesBUGS, dimnames(tempArray)[[2]])
                 samplesArray[, monitorsWeHave] <- tempArray[(burnin+1):floor(niter/thin), monitorsWeHave, drop=FALSE]
                 addToOutput('stan', samplesArray, timeResult)
            } else warning("run_stan: rstan package is required for 'stan' option.")
        },
            
        run_nimble = function() {
            for(iMCMC in seq_along(nimbleMCMCs)) {
                mcmcTag <- nimbleMCMCs[iMCMC]
                mcmcDef <- MCMCdefs[[mcmcTag]]
                mcmcConf <- eval(mcmcDef)
                mcmcConf$addMonitors(monitorVars, print = FALSE)
                mcmcConf$setThin(thin, print = FALSE)
                RmcmcFunctionList[[mcmcTag]] <<- buildMCMC(mcmcConf)
            }
            timeResult <- system.time({
                Cmodel <<- compileNimble(Rmodel)
                CmcmcFunctionList_temp <- compileNimble(RmcmcFunctionList, project = Rmodel)
                if(nNimbleMCMCs == 1) { CmcmcFunctionList[[nimbleMCMCs[1]]] <<- CmcmcFunctionList_temp
                } else                { CmcmcFunctionList                   <<- CmcmcFunctionList_temp }
            })
            addTimeResult('nimble_compile', timeResult)

            ## Record full set of model states
            allInitialModelStates <- list()
            allModelVars <- Cmodel$getVarNames(includeLogProb = TRUE)
            for(var in allModelVars)   allInitialModelStates[[var]] <- Cmodel[[var]]
                
            for(iMCMC in seq_along(nimbleMCMCs)) {
                for(var in allModelVars)   Cmodel[[var]] <<- allInitialModelStates[[var]]
                mcmcTag <- nimbleMCMCs[iMCMC]
                Cmcmc <- CmcmcFunctionList[[mcmcTag]]
                if(setSeed) set.seed(0)
                timeResult <- system.time({ Cmcmc$run(niter) })
                CmvSamples <- Cmcmc$mvSamples
                samplesArray <- as.matrix(CmvSamples, varNames = monitorVars)
                samplesArray <- samplesArray[(burnin+1):floor(niter/thin), monitorNodesNIMBLE, drop=FALSE]
                addToOutput(mcmcTag, samplesArray, timeResult)
            }
        },
        
        addToOutput = function(MCMCtag, samplesArray, timeResult) {
            output$samples[MCMCtag, , ] <<- t(samplesArray) ## makes dim1:monitors, and dim2:iter
            addTimeResult(MCMCtag, timeResult)
            summaryArray <- array(NA, c(nSummaryStats, nMonitorNodes))
            for(iStat in seq_along(summaryStats)) {
                summaryArray[iStat, ] <- apply(samplesArray, 2, summaryStatFunctions[[iStat]])
            }
            if(calculateEfficiency) {
                essDim <- which(summaryStatDimNames == 'ess')
                effDim <- which(summaryStatDimNames == 'efficiency')
                thisTime <- output$timing[[MCMCtag]]
                summaryArray[effDim, ] <- summaryArray[essDim, ] / thisTime
            }
            output$summary[MCMCtag, , ] <<- summaryArray
            if(calculateEfficiency) {
                output$efficiency$min  <<- apply(output$summary[, 'efficiency', , drop=FALSE], 1, min)
                output$efficiency$mean <<- apply(output$summary[, 'efficiency', , drop=FALSE], 1, mean)
            }
        },
        
        addTimeResult = function(MCMCtag, timeResult) {
            output$timing[MCMCtag] <<- timeResult[[3]]
        },
        
        generate_plots = function() {
            cols <- c(2:6, 8:9)
            if(nMCMCs > length(cols))    { message('too many MCMCs to plot'); return() }
            
            ## for each monitorNode, generate traceplot for each MCMC
            for(monitorNode in monitorNodesNIMBLE) {
                dev.new()
                par(mfrow=c(nMCMCs,1), mar=c(3,3,2,1), mgp=c(0,0.6,0), tcl=-0.3)
                for(i in seq_along(MCMCs)) {
                    plot(x=1:nkeep, y=output$samples[MCMCs[i], monitorNode, ],
                         main=paste0(monitorNode, ' traceplot:  ', MCMCs[i]),
                         type='l', col=cols[i], , xlab='', ylab='', xaxt='n', bty='l') }
                filename <- paste0(plotName, '_traceplots_', monitorNode, '.pdf')
                if(savePlot)   { dev.print(device = pdf, file = filename) }
            }
            
            ## density plots
            dev.new()
            par(mfrow = c(nMonitorNodes,1), mar=c(3,3,2,1), mgp=c(0,0.6,0), tcl=-0.3)
            for(monitorNode in monitorNodesNIMBLE) {
                densityList <- apply(output$samples[, monitorNode, , drop=FALSE], 1, density)
                xlim <- range(unlist(lapply(densityList, function(d) d$x)))
                xlim <- mean(xlim) + (xlim-mean(xlim)) * 1.1
                ymax <- max(unlist(lapply(densityList, function(d) d$y))) * 1.1
                plot(-100, -100, xlim=xlim, ylim=c(0,ymax),
                     main=paste0('posterior density:  ', monitorNode),
                     xlab='', ylab='', yaxt='n', bty='n')
                legend(x='topleft', legend=MCMCs, lty=1, lwd=2, col=cols[1:nMCMCs], bty='n')
                for(i in seq_along(MCMCs))     polygon(densityList[[i]], border=cols[i])
                abline(h=0, col='white')
            }
            filename <- paste0(plotName, '_densities.pdf')
            if(savePlot)   { dev.print(device = pdf, file = filename) }
        },
        
        checkMCMCdefNames = function() {
            if(!all(nimbleMCMCs %in% MCMCdefNames)) stop(paste0('missing MCMCdefs for: ', paste0(setdiff(nimbleMCMCs, MCMCdefNames), collapse=', ')))
        },
        
        writeModelFile = function() {
            writeLines(paste0('model\n', paste0(deparse(code), collapse='\n')), con=modelFileName)
        },

        fileToList = function(file) {
            if(!file.exists(file)) {
                warning(paste0('missing Stan input file: \'', file, '\''))
                return(NULL)
            }
            env <- new.env()
            source(file, local = env)
            lst <- list()
            for(name in ls(env))   lst[[name]] <- get(name, env)
            return(lst)
        },
        
        show = function() {
            cat(paste0('MCMCsuite object\n',
                       'algorithms:  ', paste0(MCMCs, collapse=', '), '\n',
                       'monitors:  ', paste0(monitorNodesNIMBLE, collapse=', '), '\n',
                       'model code:\n',
                       paste0(deparse(model), collapse='\n')))
        }
    )
)




