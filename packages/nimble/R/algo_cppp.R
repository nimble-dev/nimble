
pppFuncVirtual <- nimbleFunctionVirtual(
    run = function(N = integer(0)) returnType(double(0))
)

discrepancyFunction_BASE <- nimbleFunctionVirtual(
    run = function() returnType(double(0))
)

calcPPP <- nimbleFunction(
    setup = function(model,
                     dataNames,
                     paramNames,
                     MCMCIter,
                     thin,
                     burnIn,
                     discrepancyFunction,
                     nBootstrapSDReps,
                     discrepancyFunctionArgs){
        paramDependencies <- model$getDependencies(paramNames)
        dataDependencies <- model$getDependencies(dataNames)
        discrepancyFunctionList <- nimbleFunctionList(discrepancyFunction_BASE)
        discrepancyFunctionList[[1]] <- discrepancyFunction(model, discrepancyFunctionArgs)
        l <- ceiling(min(1000, (MCMCIter/thin - burnIn)/20))
        ## length of each block, ensures it's not too big total number
        ## of blocks available
        q <- (MCMCIter/thin - burnIn) - l + 1
        ## to sample from number of blocks to use for ppp
        h <- ceiling((MCMCIter/thin - burnIn)/l)
        ##function calculation
        blockMatrix <- matrix(0, nrow = h*l, ncol = length(paramNames))
        ## approximately of size m (number of mc samples)
    },
    run = function(MCMCOutput = double(2), nCalcSamples = integer(0),
                   useBlockMatrix = logical(0)){
        countDeviant <- numeric(nCalcSamples)
        origDataValues <- values(model, dataNames)
        for(i in 1:nCalcSamples){
            randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
            if(useBlockMatrix == TRUE){
                values(model, paramNames) <<- blockMatrix[randNum, ]
            }
            else{
                values(model, paramNames) <<- MCMCOutput[burnIn + randNum, ]
            }
            ## calculate deviance from real data
            calculate(model, paramDependencies)
            obsDeviance <- discrepancyFunctionList[[1]]$run()

            ## simulate data from model, then calculate deviance
            simulate(model, dataNames, includeData = TRUE)
            calculate(model, dataDependencies)
            simDeviance <- discrepancyFunctionList[[1]]$run()

            ## count the number of simulated data sets where the
            ## simulated data has a higher deviance than the observed
            if(simDeviance >= obsDeviance) countDeviant[i] <- 1
            values(model, dataNames) <<- origDataValues
            calculate(model, dataDependencies)
        }
        ppp <- mean(countDeviant)
        returnType(double(0))
        return(ppp)
    },
    methods = list(
        calcBlockBootstrapSD = function(MCMCOutput = double(2),
                                        nCalcSamples=double(0)){
            blockpppValues <- numeric(nBootstrapSDReps)
            ## block estimates of ppp
            for(r in 1:nBootstrapSDReps){
                for(i in 1:h){
                    randNum <- runif(1,0,1)
                    ## random starting index for blocks (post burn-in)
                    randIndex <- ceiling(randNum*q)
                    for(j in 1:l){
                        ## fill in blockMV with chosen blocks
                        blockMatrix[(i - 1)*l + j,] <<-
                            MCMCOutput[burnIn +
                                       randIndex - 1 + j,]
                    }
                }
                ## as per Caffo, calculate both Q functions using the same
                ## samples from the latent variables
                blockpppValues[r]  <- run(MCMCOutput, nCalcSamples, TRUE)
            }
            blockpppValuesSD <- sd(blockpppValues)
            returnType(double())
            return(blockpppValuesSD)
        }
    ))

calcCPPP <- function(MCMCIter,
                     burnIn,
                     nCalcSamples,
                     C.calcPPP,
                     C.mcmc = NULL,
                     samples = NULL,
                     returnSamples){
    ## calculate proportion of deviances that are greater or equal to the
    ## observed value
    if(is.null(samples)){
        C.mcmc$run(MCMCIter, reset = TRUE, simulateAll = TRUE)
        samples <- as.matrix(C.mcmc$mvSamples)
    }
    pre.pp <- C.calcPPP$run(samples, nCalcSamples, FALSE)
    bootSD <- C.calcPPP$calcBlockBootstrapSD(samples, nCalcSamples)
    if(!is.finite(pre.pp))    pre.pp <- NA
    return(list(pre.pp = pre.pp,
                bootSD = bootSD,
                samples = if(returnSamples == T) samples else NULL))
}

#' Calculates the calibrated posterior predictive p-value (CPPP) for a given NIMBLE model.  
#' 
#' EXPERIMENTAL.  Takes a NIMBLE model and calculates the CPPP for it.  The CPPP measures the degree of  "surprise" in the observed data given the specified model.
#' 
#' @param model an uncompiled nimble model 
#' @param dataNames character vector of the names of the data nodes of the model that should be re-simulated to calculate the CPPP.  It is allowed to leave out some names of data elements in the model,
#'  in which case only the data names provided here will be re-simulated for CPPP calculation.
#' @param mcmcCreator (optional) a function that takes the \code{model} as its sole argument, and returns an uncompiled nimble MCMC algorithm created from a call to \code{buildMCMC()}.  
#' This argument can be used to provide a custom MCMC algorithm for use in CPPP calculation (e.g. an algorithm with custom samplers for different model nodes).
#' If this argument not provided, NIMBLE's default MCMC samplers will be used.  
#' @param discrepancyFunction a \code{nimbleFunction} specifying the discrepancy function to be used.  See details. 
#' @param discrepancyFunctionArgs an R \code{list} containing informaton to be used in the \code{discrepancyFunction}.  See details.
#' @param cpppControl	(optional) an R \code{list} with parameters governing the CPPP calculation.  See details for specific parameters.
#' @param mcmcControl (optional) an R \code{list} with parameters governing the MCMC  algorithm,  See details for specific parameters.
#' @param origMCMCOutput (optional) a matrix with samples from the posterior distribution of the model for all model parameters.
#' Each matrix column should contain samples from a different model parameter, and should be named as that parameter.  This argument should be supplied if
#' a user has already run an MCMC and obtained posterior samples from the model, and wishes to use these samples in the CPPP algorithm. 
#' @param returnSamples a logical argument indicating whether to return all posterior samples from all simulated data sets.  This can result in a very large 
#' returned object,  Defaults to \code{FALSE}.
#' @param nCores an integer argument indicating the number of cpu cores to use for CPPP calculation.  In order to use >1 core, the \code{parallel} R package must be 
#' installed.  Only Mac and Linux operating systems support multiple cores at this time.  
#' 
#' @author Nicholas Michaud and Lauren Ponisio
#' @export
#' @details The CPPP algorithm used is based on that introduced by Hjoort et al. (2006), with some modifications.  Specfically, for each simulated data set and posterior predictive p-value
#' calculated from that data set, a bootstrap estimate of the Monte-Carlo error in the p-value is calculated as well.  These errors are used to produce an error-adjusted calibrated posterior predictive p-value. 
#'
#'  @section The \code{discrepancyFunction} Argument:
#' The \code{discrepancyFunction} argument is a \code{nimbleFunction} that will calculate the discrepancy between the model and the observations.  The choice of discrepancy function is highly model-specific, 
#' and as such no default is provided for this argument.  The \code{discrepancyFunction} provided must follow these rules:
#' \enumerate{
#'   \item The \code{discrepancyFunction} must be a \code{nimbleFunction} with \code{setup} code and {run} code.  See the User's Manual for more information on writing \code{nimbleFunctions}.
#'   \item The \code{setup} code to \code{discrepancyFunction} must accept two arguments: \code{model}, 
#'   which is the model the CPPP is being calculated for, and \code{discrepancyFunctionArgs}, 
#'   which is an R list with any needed information for calculating the discrepancy (e.g. node names, tuning parameters, etc.).
#'   \item The \code{run} code of the discrepancy function must take no arguments and return a scalar double, and must \code{contain} the \code{discrepancyFunction_BASE} virtual \code{nimbleFunction}.
#' }
#' The \code{discrepancyFunctionArgs} argument to \code{runCPPP} is passed to the \code{discrepancyFunction}.  See the examples section below for an example of a discrepancy function.
#' 
#' 
#' @section The \code{cpppControl} Argument:
#' The \code{cpppControl} argument is a list with the following elements:
#'	\itemize{
#'	  \item{\code{nSimVals}}	{
#'      an integer argument determining how many simulated posterior predictive p-values should be calculated to estimate the distribution of p-values.  Defaults to 100.
#'	  }
#'	  \item{\code{nCalcSamples}}{
#'	    an integer argument determining how many posterior samples should be used to calculate the posterior predictive p-value for each simulated data set.  Defaults to 1000.
#'	  }
#'	  \item{\code{nBootstrapSDReps}} {
#'	    an integer argument determining how many bootstrap repititions should be used to estimate the standard error of each posterior predictive p-value.  Defaults to 200.
#'	  }
#'  }
#' @section The \code{mcmcControl} Argument:
#' The \code{mcmcControl} argument is a list with the following elements:
#'	\itemize{
#'	  \item{\code{nMCMCIters}}	{
#'      an integer argument determining how many MCMC iterations should be run for each posterior predictive p-value calculation.  Defaults to 10000, but should probably be manually set.
#'	  }
#'	\item{\code{burnInProp}}{
#'	   the proportion of the MCMC chain to discard as burn-in for each posterior predictive p-value calculation.  Must be between 0 and 1.  Defaults to 0.
#'	}
#'  }
#' 
#' @return
#' an R list with four elements:
#' \itemize{
#' \item{\code{CPPP}} {The CPPP value, measuring the amount of "surprise" the data has given the model.  Lower values indicate greater "surprise".  Will be between 0 and 1.}
#' \item{\code{observed.PPP}} {The posterior predictive p-value for the observed data, along with an estimate of the Monte Carlo error for this p-value.}
#' \item{\code{simulated.PPPs}} {Posterior predictive p-values for each of the simulated data sets, along with error estimates for each p-value.  The set of simulated posterior predictive p-values provide an empirical distribution of p-values, which is used to calculate the CPPP.}
#' \item{\code{samples}} {An R list, only returned when \code{returnSamples = TRUE}.  Each element of this list will be a matrix of posterior samples from the model given a set of simulated data.  There will be \code{nSimVals} sets of samples.}
#' }
#'
#' @references 
#' Hjort, N. L., Dahl, F. A., & Steinbakk, G. H. (2006). Post-processing posterior predictive p values. \emph{Journal of the American Statistical Association}, 101(475), 1157-1174.
#' @examples
#' \dontrun{
#'
#'  
#' pumpCode <- nimbleCode({
#'   for (i in 1:N){
#'     theta[i] ~ dgamma(alpha, beta)
#'     lambda[i] <- theta[i]*t[i]
#'     x[i] ~ dpois(lambda[i])
#'   }
#'   alpha ~ dexp(1.0)
#'   beta ~ dgamma(0.1,1.0)
#' })
#' 
#' pumpConsts <- list(N = 10,
#'                    t = c(94.3, 15.7, 62.9, 126, 5.24,
#'                          31.4, 1.05, 1.05, 2.1, 10.5))
#' 
#' pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
#' 
#' pumpModel <- nimbleModel(code=pumpCode,
#'                          constants=pumpConsts,
#'                          data=pumpData)
#' 
#' # Define the discrepancy function.
#' # Note that the discrepancy function must be a nimbleFunction that takes
#' # model and discFunctionArgs as its only two arguments.
#' 
#' pumpDiscFunction <- nimbleFunction(
#'   setup = function(model, discFunctionArgs){ 
#'     dataNames <- discFunctionArgs[['dataNames']]
#'   },
#'   run = function(){
#'     dataVals <- values(model, dataNames)
#'     lambdaVals <- values(model, 'lambda')
#'     disc <- 0
#'     for(i in 1:dim(dataVals)[1]){
#'       disc <- disc + ((dataVals[i] - lambdaVals[i])/sqrt(lambdaVals[i])) # subtract the expected value, divide by the sd.
#'                                                                          # Note here the expected value and sd are the same,
#'                                                                          # as the data comes from a poisson distribution.
#'     }
#'     returnType(double(0))  # must return a double(0)
#'     return(disc)
#'   },
#'   contains = discrepancyFunction_BASE
#' )
#' 
#' pumpCPPP <- runCPPP(pumpModel,
#'                     dataNames = 'x',
#'                     discrepancyFunction = pumpDiscFunction,
#'                     discrepancyFunctionArgs = list(dataNames = 'x'),
#'                     mcmcControl = list(nMCMCIters = 20000))
#' }

runCPPP <-  function(model,
                     dataNames, ## names of the data
                     mcmcCreator = NULL,
                     discrepancyFunction = NULL,
                     discrepancyFunctionArgs = list(),
                     cpppControl = list(),
                     mcmcControl = list(),
                     origMCMCOutput = NULL,
                     returnSamples = FALSE,
                     nCores = 1){
  
    nSimVals <- cpppControl[['nSimVals']] ## number of simulated PPP values
    if(is.null(nSimVals)){
      print("Defaulting to 100 simulated PPP values.")
      nSimVals <- 100
    }
    nCalcSamples <- cpppControl[['nCalcSamples']] ## number of samples from posterior to use for ppp calculation
    if(is.null(nCalcSamples)){
      nCalcSamples <- 1000
    }
    nBootstrapSDReps  <- cpppControl[['nBootstrapSDReps']]  
    if(is.null(nBootstrapSDReps)){
      nBootstrapSDReps <- 200
    }
    nMCMCIters <- mcmcControl[['nMCMCIters']]  
    if(!is.null(origMCMCOutput)){
      origMCMCIter <- nrow(origMCMCOutput)*thin
      if(origMCMCIter <= 1){
        stop("origMCMCOutput has no rows!")
      }
      if(is.null(nMCMCIters)){
        nMCMCIters <- origMCMCIter
      }
      if(nMCMCIters < 2){
        stop("nMCMCIters must be an integer > 1")
      }
    }
    else if(is.null(nMCMCIters)){
      nMCMCIters <- 10000
      warning("Defaulting to 10,000 MCMC iterations for each MCMC run.")
    }
    else{
      origMCMCIter <- nMCMCIters
    }
    burnInProp <-  mcmcControl[['burnInProp']]  
    if(is.null(burnInProp)){
      burnInProp <- 0
    }
    else if(burnInProp >= 1 | burnInProp < 0){
        stop("burnInProp needs to be between 0 and 1")
    }
    if(is.null(discrepancyFunction)){
      stop("No discrepancy function specified.")
    }
    if(!inherits(model, "RmodelBaseClass")){
      stop("model is not an Rmodel")
    }
    if(is(model$CobjectInterface, "uninitializedField")){
      compileNimble(model)
    }

    testDataNames <- try(model[[dataNames]], silent=TRUE)
    if(inherits(testDataNames, "try-error")){
        stop(paste("dataNames", dataNames,
                   "is not the name of the data in model"))
    } else{
        test2DataNames <- all(model$expandNodeNames(dataNames) %in%
                              model$getNodeNames(dataOnly=TRUE))
        if(test2DataNames == FALSE){
            stop(paste("dataNames", dataNames,
                       "is not the name of the data in model"))
        }
    }

    if(returnSamples == TRUE){
      warning("returnSamples = TRUE, so all posterior samples from all simulations will be returned.  This may produce 
              a very large return object.")
    }
   
    ## keep track of the real data
    origData <- values(model, dataNames)

    if(is.null(mcmcCreator)){
      newMCMC.spec <-   configureMCMC(model,
                                      print=FALSE,
                                      nthin=thin)
      newMCMC <- buildMCMC(newMCMC.spec)
    }
    else{
      newMCMC <- mcmcCreator(model)
    }
    thin <- newMCMC$thin
    
    if(is.null(origMCMCOutput)){
      CnewMCMC <- compileNimble(newMCMC, project = model, resetFunctions = TRUE)
      CnewMCMC$run(nMCMCIters)
      origMCMCOutput <- as.matrix(CnewMCMC$mvSamples)
    }
    burnIn <- ceiling(burnInProp*(nMCMCIters/thin))
    
    paramNames <- colnames(origMCMCOutput)
    testParamNames <- lapply(paramNames, function(x){
      test.this.param <- try(model[[x]], silent=TRUE)
      if(inherits(test.this.param, "try-error")){
        stop(paste("paramNames", x,
                   "are not parameters in model"))
      }
    })
    test2ParamNames <- all(model$expandNodeNames(paramNames) %in%
                             model$getNodeNames(includeData=FALSE,
                                                  stochOnly=TRUE))
    if(test2ParamNames == FALSE){
      stop(paste("paramNames", paramNames,
                 "are not parameters in model"))
    }
    

    ## sample posterior, simulate data from sample
    paramDependencies <- model$getDependencies(paramNames)
    modelcalcPPP <- calcPPP(model,
                            dataNames,
                            paramNames,
                            origMCMCIter,
                            thin,
                            burnIn,
                            discrepancyFunction = discrepancyFunction,
                            nBootstrapSDReps=nBootstrapSDReps,
                            discrepancyFunctionArgs)
    C.calcPPP <- compileNimble(modelcalcPPP,
                               project = model)

    ## calculate deviances
    obs.CPPP <- calcCPPP(origMCMCIter,
                         burnIn,
                         nCalcSamples=nCalcSamples,
                         C.calcPPP,
                         samples = origMCMCOutput,
                         returnSamples = returnSamples)


    CPPPTask <- function(coreNum, nBootIter){
        R.newModel <- model$newModel()
        C.newModel <- compileNimble(R.newModel)
        if(is.null(mcmcCreator)){
            newMCMC.spec <-   configureMCMC(R.newModel,
                                            print=FALSE,
                                            nthin=thin)
            newMCMC <- buildMCMC(newMCMC.spec)
        }
        else{
            newMCMC <- mcmcCreator(R.newModel)
        }

        modelcalcPPP <- calcPPP(R.newModel,
                                dataNames,
                                paramNames,
                                nMCMCIters,
                                thin,
                                burnIn,
                                discrepancyFunction = discrepancyFunction,
                                nBootstrapSDReps=nBootstrapSDReps,
                                discrepancyFunctionArgs)

        C.newMcmc <- compileNimble(newMCMC, project = R.newModel)
        C.newPppFunc <- compileNimble(modelcalcPPP,
                                      project = R.newModel)
        out <- list()
        for(iteration in 1:nBootIter){
            message(paste("refitting data iteration", iteration, "for core number", coreNum))
            simulate(C.newModel,  includeData =  TRUE)
            out[[iteration]] <- calcCPPP(nMCMCIters,
                                         burnIn,
                                         nCalcSamples,
                                         C.calcPPP = C.newPppFunc,
                                         C.mcmc = C.newMcmc,
                                         returnSamples = returnSamples)
        }
        return(out)
    }

    libError <-  try(library('parallel'), silent = TRUE)
    if(inherits(libError, 'try-error') && nCores > 1){
      warning("The 'parallel' package must be installed to use multiple cores for CPPP calculation.  Defaulting to one core.")
      nCores <- 1
    }
    
    if(nCores > 1){
    ## simulate the ppp values
      sim.ppp.output <- mclapply(X = 1:nCores, FUN = CPPPTask,
                                 nBootIter =
                                     ceiling(nSimVals/nCores))
    }
    else{
      sim.ppp.output <- lapply(X = 1, FUN = CPPPTask,
                                 nBootIter =
                                   ceiling(nSimVals/nCores))
    }

    ## extract simulated ppp and boot SDs
    sim.ppp <- sapply(sapply(sim.ppp.output,
                             function(x) x), function(x)  x$pre.pp)
    sim.vars <- (sapply(sapply(sim.ppp.output, function(x) x),
                        function(x) x$bootSD))^2

    ## approximate the distribution of observed ppp and simulated ppp
    diff.ppp <- sim.ppp - obs.CPPP$pre.pp 
    diff.vars.ppp <- obs.CPPP$bootSD^2 + sim.vars


    ## calculate the average prob that simulated values are less than
    ## the observed
    CPPP <- mean(pnorm(0, diff.ppp, diff.vars.ppp),
                 na.rm=TRUE)

    ## extract chains and diagnostics
    if(returnSamples){
        sim.samples <- lapply(sapply(sim.ppp.output, function(x) x),
                              function(x) x$samples)
    }
    else sim.samples <- NULL

    ## output real data and model
    values(model, dataNames) <- origData

    out <- list(CPPP=CPPP,
                observed.PPP=c(estimate=obs.CPPP$pre.pp,
                               SE=(obs.CPPP$bootSD)),
                simulated.PPPs=cbind(esimate=sim.ppp,
                                     SE=sqrt(sim.vars)),
                samples=sim.samples)
    return(out)
}
