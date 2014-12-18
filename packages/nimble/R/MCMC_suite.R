


#' Class \code{MCMCsuite}
#' @aliases MCMCsuite
#' @export
#' @description
#' Creates, runs, and compiles samples from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
#' This can include BUGS and JAGS MCMCs, as well as NIMBLE MCMC algorithms.
#' Trace plots and density plots for the MCMC samples may also be plotted.
#' See documentation for method \code{initialize()}, for details of creating and running an MCMC Suite.
#' Upon executing the MCMC algorithms, an internal \code{output} variable is created, which may be accessed via \code{suite$output}.
#' \code{output} is a named list containing three elements:
#' \code{output$samples} is a 3-dimensional array containing all MCMC samples.
#' \code{output$summary} is a 3-dimensional array containing all summary statistics.
#' \code{output$timing} is a numeric vector containing timing information.
#' See the NIMBLE User Manual for more information about these \code{output} list elements.
#' @author Daniel Turek
#' @examples
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' suite <- MCMCsuite(code,
#'                    data = list(x=3),
#'                    inits = list(mu=0),
#'                    niter = 10000,
#'                    thin = 10,
#'                    monitors = 'mu',
#'                    MCMCs = c('bugs', 'jags', 'nimble'),
#'                    summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
#'                    plotName = 'example')
MCMCsuite <- setRefClass(

    Class = 'MCMCsuite',
    
    fields = list(
        ## set in initialize()
        model = 'ANY',   ## parsed expression for the model code; must be contained in { ... }    --- ORIGINAL ARGUMENT
        constants = 'list',   ## list of the constants (ORIGINAL ARGUMENT)
        Rmodel = 'ANY',   ## Rmodel object
        
        ## setData()
        data = 'list',   ## list of the data    --- ORIGINAL ARGUMENT
        BUGSdata = 'list',   ## data list used for BUGS and JAGS.  is equal to c(constantList, dataList)
        
        ## setInits()
        inits = 'list',  ## named list of initial values used for all MCMC algorithms    --- ORIGINAL ARGUMENT
        
        ## setMonitors()
        monitors = 'character',    ## the original character vector argument to initialize()    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        monitorVars = 'character',    ## character vector of VARIABLE names of parameters to save
        monitorNodesNIMBLE = 'character',  ## character vector of the monitor node names, with spaces as in nimble: 'y[1, 1]'
        monitorNodesBUGS = 'character',    ## same as monitorNodes, except for BUGS: no spaces in node names: 'y[1,1]'
        nMonitorNodes = 'numeric',   ## number of monitorNodes
        
        ## set in initialize()
        niter = 'numeric',    ## number of MCMC iterations to run    --- ORIGINAL ARGUMENT
        burnin = 'numeric',   ## burn-in period, the number of initial samples to discard, prior to thinning    --- ORIGINAL ARGUMENT
        thin = 'numeric',   ## thinning interval    --- ORIGINAL ARGUMENT

        ## set in setNkeep()
        nkeep = 'numeric',   ## number of samples we'll keep. equal to (niter/thin - burnin)
        
        ## setSummaryStats()
        summaryStats = 'character',    ## character vector of parseable summary statistic functions    --- ORIGINAL ARGUMENT
        summaryStatFunctions = 'list',  ## list of the function objects for summary statistics
        summaryStatDimNames = 'character',   ## character vector of the dimension names in output$summary
        nSummaryStats = 'numeric',   ## the number of summary statistics
        
        ## setMCMCs()
        MCMCs = 'character',   ## character vector of the MCMC analyses.  'bugs', 'jags', or anything else is nimble    --- ORIGINAL ARGUMENT
        bugsMCMCflag = 'logical',   ## whether 'bugs' is in MCMCs
        jagsMCMCflag = 'logical',   ## whether 'jags' is in MCMCs
        nimbleMCMCs = 'character',    ## the names of the remaining (presumably nimble) MCMCs
        nNimbleMCMCs = 'numeric',    ## the number of remaining (nimble) MCMCs
        nimbleMCMCflag = 'logical',   ## a flag indicating whether there are any remaining (nimble) MCMCs
        nMCMCs = 'numeric',   ## the number of MCMC algorithms being ran
        
        ## setMCMCdefs()
        MCMCdefs = 'list',   ## named list of {} expression code blocks, corresponding the setup for nimble MCMCs    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        MCMCdefNames = 'character',   ## names of the MCMCdefs list
        
        ## set in initialize()
        bugs_directory = 'character',    ## directory for BUGS program    --- ORIGINAL ARGUMENT
        bugs_program = 'character',     ## program for BUGS    --- ORIGINAL ARGUMENT
        makePlot = 'logical',    ## whether to generate plots    --- ORIGINAL ARGUMENT
        savePlot = 'logical',   ## whether or not to save plot PDFs    --- ORIGINAL ARGUMENT
        plotName = 'character',     ## name of the file where we save density and trace plots    --- ORIGINAL ARGUMENT
        debug = 'logical',   ## whether to enter browser() before running each algorithm    --- ORIGINAL ARGUMENT
        modelFileName = 'character',     ## name of the text file where we write the model code, set to a fixed value
        cols = 'numeric',   ## vector of the colors used for plotting, set to a fixed value
        
        ## set in run()
        bugs_out = 'ANY',    ## output from bugs()
        jags_mod = 'ANY',   ## jags model object
        jags_out = 'ANY',    ## output from coda.samples() for jags
        Cmodel = 'ANY',   ## compiled Cmodel object
        RmcmcFunctionList = 'list',    ## list of the R (nimble) MCMC functions
        CmcmcFunctionList = 'list',    ## list of the C (nimble) MCMC functions
        output = 'list'   ## list of numeric outputs: samples, summary, timing
    ),
    
    methods = list(
        
        initialize = function(model,
                              constants      = list(),
                              data           = list(),
                              inits          = list(),
                              monitors       = character(),
                              niter          = 10000,
                              burnin         = 2000,
                              thin           = 1,
                              summaryStats   = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp'),
                              MCMCs          = c('bugs', 'jags', 'nimble', 'nimble_RW', 'nimble_slice'),
                              MCMCdefs       = list(),
                              bugs_directory = 'C:/WinBUGS14',
                              bugs_program   = 'WinBUGS',
                              makePlot       = TRUE,
                              savePlot       = TRUE,
                              plotName       = 'MCMCsuite',
                              debug          = FALSE,
                              runNow         = TRUE) {
'
Creates and runs an MCMC Suite.
By default, this will execute the specified MCMCs, record all samples, generate summary statistics, and create and save trace plots and posterior density plots.
This default behavior can ben altered via a variety of input arguments.
Following execution of the MCMC algorithms, the interal \'output\' list is created, which contains named \'samples\', \'summary\', and \'timing\' elements.
See the NIMBLE User Manual for more information about these \'output\' list objects.

Arguments:

model: The quoted code expression representing the model, such as the return value from a call to nimbleCode({...}).
No default value, this is a required argument.

constants: A named list giving values of constants for the model.
This is the same as the \'constants\' argument which would be passed to nimbleModel(...).
Default value is list().

data: A named list giving the data values for the model.
This is the same as the \'data\' argument which would be passed to nimbleModel(...) or model$setData(...).
Default value is list().

inits: A named list giving the initial values for the model.
This is the same as the \'inits\' argument which would be passed to nimbleModel(...) or model$setInits(...).
Default value is list().

monitors: A character vector giving the node names or variable names to monitor.
The samples corresponding to these nodes will be stored in the output samples, will have summary statistics calculated, and density and trace plots generated.
Default value is all top-level stochastic nodes of the model.

niter: Number of MCMC iterations to run.
This applies to all MCMC algorithms in the suite.
Default value is 10,000.

burnin: Number of initial, post-thinning, MCMC iterations to discard.
Default value is 2,000.

thin: Thinning interval for the MCMC samples.
This applies to all MCMC algorithms in the suite.  The thinning occurs prior to the burnin samples being discarded.
Default value is 1.

summaryStats: A character vector, specifying the summary statistics to calculate on the MCMC samples.
Each element may be the character name of an exisiting R function (possibly user-defined) which acts on a numeric vector and returns a scalar (e.g., \'mean\' or \'sd\'),
or a character string which when parsed and evaluted will define such a function (e.g., \'function(x) mean(sqrt(x))\').
Default value is c(\'mean\', \'median\', \'sd\', \'CI95_low\', \'CI95_upp\'), where the final two elements are functions which calculate the limits of a 95% Bayesian credible interval.

MCMCs: A character vector specifying the MCMC algorithms to run.
\'bugs\' specifies a WinBUGS/BUGS MCMC, and \'jags\' specifies a JAGS MCMC.
All other entries will be interpreted as NIMBLE MCMC algorithms.
All NIMBLE algorithms must have an associated entry in the MCMCdefs argument.
The MCMCdefs provides, by default, definitions for \'nimble\', \'nimble_RW\', and \'nimble_slice\'.
Default value is c(\'bugs\', \'jags\', \'nimble\', \'nimble_RW\', \'nimble_slice\').

MCMCdefs: A named list of MCMC definitions.  The names of list elements should corespond to any custom MCMC algorithms specified in the \'MCMCs\' argument.
The list elements should be quoted expressions, enclosed in {} braces.  When executed, the internal code must return an MCMC specification object, 
specifying the corresponding MCMC algorithm; in particular, setting the appropriate samplers.  The code may assume existance of the R model object \'Rmodel\',
and must *return* the MCMC specification object.  Therefore, the final line of such a code block would frequently be a standalong \'mcmcspec\', to return this object.

bugs_directory: A character string giving the directory of the executable BUGS program for the WinBUGS/BUGS MCMC.
This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
Default value is \'C:/WinBUGS14\'.

bugs_program: A character string giving the name of the BUGS program, for the WinBUGS/BUGS MCMC.
This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
Default value is \'WinBUGS\'.

makePlot: Logical argument, specifying whether to generate the trace plots and posterior density plots, for each monitored node.
Default value is TRUE.

savePlot: Logical argument, specifying whether to save the trace plots and density plots.
Plots will be saved into the current working directory.
Only used when makePlot == TRUE.
Default value is TRUE.

plotName: Character string, giving the file name for saving the trace plots and density plots.
Only used when makePlot == TRUE and savePlot == TRUE.
Default value is \'MCMCsuite\'.

debug: Logical argument, specifying whether to enter a broswer() at the onset of executing each MCMC algrithm.
For use in debugging individual MCMC algorithms, if necessary.
Default value is FALSE.

runNow: Logical argument, specifying whether to run the MCMC algorithms straight away.
If FALSE, the MCMC algorithms may be executed by calling suiteObject$run().
Default value is TRUE.
'
            model <<- model
            constants <<- constants
            Rmodel <<- nimbleModel(model, constants = constants)
            setData(data)
            setInits(inits)
            setMonitors(monitors)
            niter <<- niter
            burnin <<- burnin
            thin <<- thin
            setNkeep()
            setSummaryStats(summaryStats)
            setMCMCs(MCMCs)
            setMCMCdefs(MCMCdefs)
            bugs_directory <<- bugs_directory
            bugs_program <<- bugs_program
            makePlot <<- makePlot
            savePlot <<- savePlot
            plotName <<- plotName
            debug <<- debug
            modelFileName <<- 'model.txt'
            cols <<- c(2:6, 8:9)
            if(runNow) run()
        },
        
        setData = function(data) {
            data <<- data
            BUGSdata <<- c(constants, data)
            Rmodel$resetData()
            Rmodel$setData(data = data)
        },
        
        setInits = function(inits) {
            inits <<- inits
            Rmodel$setInits(inits)
        },
        
        setMonitors = function(newMonitors) {
            if(length(newMonitors) == 0) newMonitors <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
            newMonitors <- Rmodel$expandNodeNames(newMonitors)
            dataFlags <- unlist(lapply(newMonitors, function(mon) eval(parse(text=mon, keep.source=FALSE)[[1]], envir=Rmodel$isDataEnv)))
            newMonitors <- newMonitors[!dataFlags]
            monitors <<- newMonitors
            monitorVars <<- unique(removeIndexing(monitors))
            monitorNodesNIMBLE <<- monitors
            monitorNodesBUGS <<- gsub(' ', '', monitorNodesNIMBLE)
            nMonitorNodes <<- length(monitorNodesNIMBLE)
        },
        
        setNiter = function(newNiter) {
            niter <<- newNiter
            setNkeep()
        },

        setBurnin = function(newBurnin) {
            burnin <<- newBurnin
            setNkeep()
        },
        
        setThin = function(newThin) {
            thin <<- newThin
            setNkeep()
        },

        setNkeep = function() {
            nkeep <<- floor(niter/thin) - burnin
        },
        
        setSummaryStats = function(summaryStats) {
            summaryStats <<- summaryStats
            CI95_low <- function(x) quantile(x, probs = 0.025)
            CI95_upp <- function(x) quantile(x, probs = 0.975)
            summaryStatFunctions <<- lapply(summaryStats, function(txt) eval(parse(text=txt)[[1]]))
            summaryStatDimNames <<- gsub('function *\\(.*?\\)', '', summaryStats)
            summaryStatDimNames <<- gsub('^ *', '', summaryStatDimNames)
            nSummaryStats <<- length(summaryStats)
        },
        
        setMCMCs = function(MCMCs) {
            MCMCs <<- unique(MCMCs)
            bugsMCMCflag <<- 'bugs' %in% MCMCs
            jagsMCMCflag <<- 'jags' %in% MCMCs
            nimbleMCMCs <<- setdiff(MCMCs, c('bugs', 'jags'))
            nNimbleMCMCs <<- length(nimbleMCMCs)
            nimbleMCMCflag <<- if(nNimbleMCMCs > 0) TRUE else FALSE
            nMCMCs <<- length(MCMCs)
        },
        
        setMCMCdefs = function(newMCMCdefs) {
            MCMCdefs <<- list(nimble       = quote(configureMCMC(Rmodel)),
                              nimble_RW    = quote(configureMCMC(Rmodel, onlyRW    = TRUE)),
                              nimble_slice = quote(configureMCMC(Rmodel, onlySlice = TRUE)))
            MCMCdefs[names(newMCMCdefs)] <<- newMCMCdefs
            MCMCdefNames <<- names(MCMCdefs)
        },
        
        run = function() {
            checkMCMCdefNames()
            init_output()
            writeModelFile()
            if(bugsMCMCflag)     run_bugs()
            if(jagsMCMCflag)     run_jags()
            if(nimbleMCMCflag)   run_nimble()
            unlink(modelFileName)
            if(makePlot) generate_plots()
            return(output)
        },
        
        init_output = function() {
            samples <- array(NA, dim = c(nMCMCs, nMonitorNodes, nkeep))
            dimnames(samples) <- list(MCMCs, monitorNodesNIMBLE, NULL)
            summary <- array(NA, dim = c(nMCMCs, nSummaryStats, nMonitorNodes))
            dimnames(summary) <- list(MCMCs, summaryStatDimNames, monitorNodesNIMBLE)
            timing <- rep(NA, nMCMCs+1)
            names(timing) <- c(MCMCs, 'nimble_compile')
            output <<- list(samples=samples, summary=summary, timing=timing)
        },
        
        run_bugs = function() {
            library(R2WinBUGS)
            if(debug) browser()
            timeResult <- system.time({
                bugs_out <<- bugs(data=BUGSdata, inits=list(inits), parameters.to.save=monitorVars, model.file=modelFileName,
                                  n.chains=1, n.iter=niter, n.burnin=0, n.thin=thin, bugs.directory=bugs_directory, program=bugs_program)
            })
            tempArray <- bugs_out$sims.array[, 1, ]        ## must use sims.array
            samplesArray <- tempArray[(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
            addToOutput('bugs', samplesArray, timeResult)
        },
        
        run_jags = function() {
            library(rjags)
            if(debug) browser()
            timeResult <- system.time({
                jags_mod <<- jags.model(file=modelFileName, data=BUGSdata, inits=inits, n.chains=1, quiet=FALSE)
                jags_out <<- coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=thin)
            })
            samplesArray <- jags_out[[1]][(burnin+1):floor(niter/thin), monitorNodesBUGS, drop=FALSE]
            addToOutput('jags', samplesArray, timeResult)
        },
        
        run_nimble = function() {
            if(debug) browser()
            for(iMCMC in seq_along(nimbleMCMCs)) {
                mcmcTag <- nimbleMCMCs[iMCMC]
                mcmcDef <- MCMCdefs[[mcmcTag]]
                mcmcspec <- eval(mcmcDef)
                mcmcspec$addMonitors(monitorVars, print = FALSE)
                mcmcspec$setThin(thin, print = FALSE)
                RmcmcFunctionList[[mcmcTag]] <<- buildMCMC(mcmcspec)
            }
            
            timeResult <- system.time({
                Cmodel <<- compileNimble(Rmodel)
                CmcmcFunctionList_temp <- compileNimble(RmcmcFunctionList, project = Rmodel)
                if(nNimbleMCMCs == 1) { CmcmcFunctionList[[nimbleMCMCs[1]]] <<- CmcmcFunctionList_temp
                } else                { CmcmcFunctionList                   <<- CmcmcFunctionList_temp }
            })
            addTimeResult('nimble_compile', timeResult)
            
            for(iMCMC in seq_along(nimbleMCMCs)) {
                Cmodel$setInits(inits);     calculate(Cmodel)
                mcmcTag <- nimbleMCMCs[iMCMC]
                Cmcmc <- CmcmcFunctionList[[mcmcTag]]
                timeResult <- system.time({ Cmcmc(niter) })
                CmvSamples <- Cmcmc$mvSamples
                samplesArray <- as.matrix(CmvSamples, varNames = monitorVars)
                samplesArray <- samplesArray[(burnin+1):floor(niter/thin), monitorNodesNIMBLE, drop=FALSE]
                addToOutput(mcmcTag, samplesArray, timeResult)
            }
        },
        
        addToOutput = function(MCMCtag, samplesArray, timeResult) {
            summaryArray <- array(NA, c(nSummaryStats, nMonitorNodes))
            for(iStat in seq_along(summaryStats)) {
                summaryArray[iStat, ] <- apply(samplesArray, 2, summaryStatFunctions[[iStat]])
            }
            output$samples[MCMCtag, , ] <<- t(samplesArray)      ## makes dim1:monitors, and dim2:iter
            output$summary[MCMCtag, , ] <<- summaryArray
            addTimeResult(MCMCtag, timeResult)
        },
        
        addTimeResult = function(MCMCtag, timeResult) {
            output$timing[MCMCtag] <<- timeResult[[3]] / 60
        },
        
        generate_plots = function() {
            if(nMCMCs > length(cols))    { cat('too many MCMCs to plot'); return() }
            
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
            writeLines(paste0('model\n', paste0(deparse(model), collapse='\n')), con=modelFileName)
        },
        
        show = function() {
            cat(paste0('MCMCsuite object\n',
                       'algorithms:  ', paste0(MCMCs, collapse=', '), '\n',
                       'monitors:  ', paste0(monitorNodesNIMBLE, collapse=', '), '\n',
                       'model:\n',
                       paste0(deparse(model), collapse='\n')))
        }
    )
)






# 
# artificial generation of output$samples
# 
# monitorNodesNIMBLE <- c('x', 'y', 'z')
# MCMCs <- c('bugs', 'jags', 'nimble', 'nimble_slice', 'df')
# niter <- 1000
# plotName <- 'MCMCsuite'
# nMonitorNodes <- length(monitorNodesNIMBLE)
# nMCMCs <- length(MCMCs)
# nkeep <- niter
# samples <- array(NA, dim = c(nMCMCs, nMonitorNodes, nkeep))
# dimnames(samples) <- list(MCMCs, monitorNodesNIMBLE, NULL)
# for(iMonitor in 1:nMonitorNodes) {
#     for(iMCMC in 1:nMCMCs) {
#         samples[iMCMC, iMonitor, ] <- rnorm(nkeep, mean=iMonitor, sd=iMCMC)
#     }
# }
# output <- list(samples=samples)
# 





