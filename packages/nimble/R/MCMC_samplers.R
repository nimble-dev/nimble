


####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################

#' @rdname samplers
#' @export
sampler_BASE <- nimbleFunctionVirtual(
    name = 'sampler_BASE',
    methods = list(
        reset        = function() { },
        before_chain = function(MCMCniter = double(), MCMCnburnin = double(), MCMCchain = double()) { },
        after_chain  = function() { }
    ),
    methodControl = list(
        before_chain = list(required = FALSE),
        after_chain  = list(required = FALSE)
    )
)



####################################################################
### prior_samples sampler for using exisiting samples as prior #####
####################################################################

#' @rdname samplers
#' @export
sampler_prior_samples <- nimbleFunction(
    name = 'sampler_prior_samples',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        samples     <- extractControlElement(control, 'samples',     error = '  [Error] prior_samples sampler missing required control argument: samples')
        randomDraws <- extractControlElement(control, 'randomDraws', FALSE)   ## default is taking sequential draws
        ## node list generation
        targetExpanded <- model$expandNodeNames(target)
        targetAsScalar <- model$expandNodeNames(targetExpanded, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, targetExpanded)
        calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
        ## numeric value generation
        k <- length(targetAsScalar)
        if(is.null(dim(samples)))   samples <- matrix(samples, ncol = 1)    ## make vectors into 1-column array
        nSamples <- dim(samples)[1]
        ind <- 0
        ## checks
        if(length(dim(samples)) != 2)   stop('prior_samples sampler \'samples\' control argument must be a 2-dimensional array, but value provided was a ', length(dim(samples)), '-dimensional array')
        if(!(storage.mode(samples) %in% c('integer', 'double')))   stop('prior_samples sampler \'samples\' control argument must be numeric or integer type')
        if(dim(samples)[2] != k)   stop('prior_samples sampler \'samples\' control argument had ', dim(samples)[2], ' columns, but target nodes have ', k, ' scalar elements. These numbers must be equal.')
        if(any(model$getNodeType(target) == 'stoch'))   messageIfVerbose('  [Note] \'prior_samples\' sampler has been assigned to one or more stochastic nodes. The prior distribution for these nodes will be overridden by the prior samples.')
    },
    run = function() {
        if(randomDraws) {
            ind <<- ceiling(runif(1, 0, nSamples))   ## random draws
        } else {
            ind <<- ind + 1                          ## sequential draws (the default)
            if(ind > nSamples)   ind <<- 1           ## recycle sequential draws, if necessary
        }
        values(model, targetExpanded) <<- samples[ind, 1:k]
        model$calculate(calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    },
    methods = list(
        before_chain = function(MCMCniter = double(), MCMCnburnin = double(), MCMCchain = double()) {
            ## issue a note if sequential draws from prior samples will have to be recycled
            if((!randomDraws) & (MCMCniter > nSamples))   cat('  [Note] prior_samples sampler will recycle sequential draws from prior samples, since ', nSamples, ' samples were provided, and ', MCMCniter, ' MCMC iterations will be run.\n')
        },
        reset = function() {
            ind <<- 0
        }
    )
)



####################################################################
### posterior_predictive sampler for trailing non-data networks ####
####################################################################

#' @rdname samplers
#' @export
sampler_posterior_predictive <- nimbleFunction(
    name = 'sampler_posterior_predictive',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        simNodes  <- model$getDependencies(target, downstream = TRUE, includePredictive = TRUE)
        calcNodes <- model$getDependencies(target, downstream = TRUE, includePredictive = TRUE, stochOnly = TRUE)
        ## checks
        if(!all(model$expandNodeNames(target) %in% model$getNodeNames(predictiveOnly = TRUE)))
            stop('posterior_predictive sampler should not be assigned to non-posterior-predictive nodes')
    },
    run = function() {
        model$simulate(simNodes)     ## simulate values of all downstream dependent nodes (stoch and determ both)
        model$calculate(calcNodes)   ## update logProbs of dependent stochastic nodes
        nimCopy(from = model, to = mvSaved, row = 1, nodes = simNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    )
)



####################################################################
### binary Gibbs sampler ###########################################
####################################################################

#' @rdname samplers
#' @export
sampler_binary <- nimbleFunction(
    name = 'sampler_binary',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
        ## checks
        if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
        if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
    },
    run = function() {
        currentLogProb <- model$getLogProb(calcNodes)
        model[[target]] <<- 1 - model[[target]]
        otherLogProbPrior <- model$calculate(target)
        if(otherLogProbPrior == -Inf) {
            otherLogProb <- otherLogProbPrior
        } else {
            otherLogProb <- otherLogProbPrior + model$calculate(calcNodesNoSelf)
        }
        acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
        jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
        if(jump) {
            ##model$calculate(calcNodesPPomitted)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        }
    },
    methods = list(
        reset = function() { }
    )
)



###########################################################################
### categorical sampler for dcat distributions ############################
###########################################################################

#' @rdname samplers
#' @export
sampler_categorical <- nimbleFunction(
    name = 'sampler_categorical',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        length <- extractControlElement(control, 'length', 'prob')
        check <- extractControlElement(control, 'check', TRUE)
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
        ## numeric value generation
        if(is.character(length)) {
            k <- length(model$getParam(target, length))
        } else if(is.numeric(length)) {
            k <- length
        } else stop('Invalid \'length\' control parameter provided for categorical sampler.\nSee help(samplers) for details of the categorical sampler.')
        probs <- numeric(k)
        logProbs <- numeric(k)
        ## checks
        if(length(targetAsScalar) > 1)  stop('cannot use categorical sampler on more than one target node')
        if(check && length == 'prob' && model$getDistribution(target) != 'dcat') stop('Can only use categorical sampler on node with dcat distribution.\nUse control argument \'check = FALSE\' to allow use on other distributions.')
        if(!is.numeric(k) || k <= 0 || k != round(k))   stop('Invalid \'length\' control parameter provided for categorical sampler.\nSee help(samplers) for details of the categorical sampler.')
    },
    run = function() {
        currentValue <- model[[target]]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
        for(i in 1:k) {
            if(i != currentValue) {
                model[[target]] <<- i
                logProbPrior <- model$calculate(target)
                if(logProbPrior == -Inf) {
                    logProbs[i] <<- -Inf
                } else {
                    if(is.nan(logProbPrior)) {
                        logProbs[i] <<- -Inf
                    } else {
                        logProbs[i] <<- logProbPrior + model$calculate(calcNodesNoSelf)
                        if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
                    }
                }
            }
        }
        maxLP <- max(logProbs)
        if(maxLP == Inf | is.nan(maxLP))   cat("Warning: categorical sampler for '", target, "' encountered an invalid model density, and sampling results are likely invalid.\n")
        logProbs <<- logProbs - maxLP
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)   ## rcat normalizes the probabilities internally
        if(!is.na(newValue) & newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)
            ##model$calculate(calcNodesPPomitted)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        }
    },
    methods = list(
        reset = function() { }
    )
)


####################################################################
### noncentered joint sampler on location or scale hyperparameter ##
### and associated set of (random) effects                        ##
####################################################################

## This provides a sampler for the location or scale hyperparameter
## for a set of random effects, where the random effects are then
## shifted or scaled based on the proposed hyperparameter as if
## a noncentered parameterization were used. When added to a model
## parameterized with a centered parameterization, i.e.,
## b[i] ~ dnorm(mu, sd = sigma)`, this is a variation on the ASIS/
## interweaving sampling scheme of Yu and Meng 2011.

## The hyperparameter sampling can either use RW or slice.
## An adjustment to the acceptance probability based on the
## Jacobian of the transformation in the case of sampling
## the scale is included. For normal random effects, this
## adjustment cancels with the density of the random
## effects.

#' @rdname samplers
#' @export
sampler_noncentered <- nimbleFunction(
    name = 'sampler_noncentered',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        type <- extractControlElement(control, 'sampler', 'RW')
        param <- extractControlElement(control, 'param', 'location')
        if(!param %in% c("location", "scale")) 
            stop("noncentered sampler: `sampleParam` control element must be either 'location' or 'scale'")

        if(param == "location") param <- 0
        if(param == "scale") param <- 1

        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)
            stop("sampler noncentered: target for noncentered sampler must be a scalar node")
        effects <- model$getDependencies(target, stochOnly = TRUE, self = FALSE)
        if(!length(effects))
            stop("sampler noncentered: there are no stochastic dependencies of target")
        if(any(model$isMultivariate(effects)))
            stop("sampler noncentered: all dependent nodes must be univariate")
        
        dists <- model$getDistribution(effects)
        uniqueDists <- which(!duplicated(dists))
        for(idx in uniqueDists) {
            result <- try(model$getParam(effects[idx], 'mean'), silent = TRUE)
            if(inherits(result, 'try-error'))
                stop("sampler noncentered: the distribution `", dists[idx], "` does not have a `mean` function, which is required for all dependent nodes")
        }
        
        if(!all(dists %in% c('dnorm','dbeta','dgamma','dinvgamma','dunif')))
            stop("sampler noncentered: dependent nodes must be scalars whose distribution has a function for calculating its mean in nimble")
        
        if(type == "RW") {
            samplerFunction <- sampler_RW_noncentered
        } else {
            if(type == "slice") {
                samplerFunction <- sampler_slice_noncentered
            } else stop("noncentered sampler: `sampler` control element must be either 'RW' or 'slice'")
        }

        ## Pass along control elements to (variant of) RW or slice sampler, with info on the effects to be sampled.
        if(is.list(control)) {
            control[['noncenteredEffects']] <- effects
            control[['noncenteredParam']] <- param
        } else {
            control <- list('noncenteredEffects' = effects, 'noncenteredParam' = param)
        }
        
        samplerList <- nimbleFunctionList(sampler_BASE)
        samplerList[[1]] <- samplerFunction(model = model, mvSaved = mvSaved, target=target, control=control)
    },
    run = function() {
        samplerList[[1]]$run()
    },
    methods = list(
        reset = function() {
            samplerList[[1]]$reset()
        }
    )
)
        

####################################################################
### scalar RW sampler with normal proposal distribution ############
####################################################################

#' @rdname samplers
#' @export
sampler_RW <- nimbleFunction(
    name = 'sampler_RW',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale            <- extractControlElement(control, 'log',                 FALSE)
        reflective          <- extractControlElement(control, 'reflective',          FALSE)
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodes
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory      <- c(0, 0)   ## scaleHistory
        acceptanceHistory <- c(0, 0)   ## scaleHistory
        saveMCMChistory <- getNimbleOption('MCMCsaveHistory')
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
        if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
        if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
                       propValue <- currentValue * exp(propLogScale)
        } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
        if(reflective) {
            lower <- model$getBound(target, 'lower')
            upper <- model$getBound(target, 'upper')
            while(propValue < lower | propValue > upper) {
                if(propValue < lower) propValue <- 2*lower - propValue
                if(propValue > upper) propValue <- 2*upper - propValue
            }
        }
        model[[target]] <<- propValue
        logMHR <- model$calculateDiff(target)
        if(logMHR == -Inf) {
            jump <- FALSE
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
        } else {
            logMHR <- logMHR + model$calculateDiff(calcNodesNoSelf) + propLogScale
            jump <- decide(logMHR)
            if(jump) {
                ##model$calculate(calcNodesPPomitted)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
            }
        }
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)                 ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale                ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                }
                gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                ## If there are upper and lower bounds, enforce a maximum scale of
                ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
                ## Otherwise, for a poorly-informed posterior,
                ## the scale could grow without bound to try to reduce
                ## acceptance probability.  This creates enormous cost of
                ## reflections.
                if(reflective) {
                    lower <- model$getBound(target, 'lower')
                    upper <- model$getBound(target, 'upper')
                    if(scale >= 0.5*(upper-lower)) {
                        scale <<- 0.5*(upper-lower)
                    }
                }
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        setScale = function(newScale = double()) {
            scale         <<- newScale
            scaleOriginal <<- newScale
        },
        getScaleHistory = function() {       ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(scaleHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(acceptanceHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
            }
            gamma1 <<- 0
        }
    )
)

## version for use with noncentered sampler

#' @rdname samplers
#' @export
sampler_RW_noncentered <- nimbleFunction(
    name = 'sampler_RW_noncentered',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale            <- extractControlElement(control, 'log',                 FALSE)
        reflective          <- extractControlElement(control, 'reflective',          FALSE)
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)

        ## sampler_RW_noncentered is used by sampler_noncentered to do joint, noncentered sampling of
        ## target and dependent (random) effects, a version of the ASIS/interweaving sampler of Yu and Meng (2011).
        noncenteredEffects <- extractControlElement(control, 'noncenteredEffects', "")
        noncenteredParam <- extractControlElement(control, 'noncenteredParam', 0)
        
        noncenteredLen <- length(model$expandNodeNames(noncenteredEffects))
        noncenteredMean <- rep(0, noncenteredLen)
        ccList <- mcmc_determineCalcAndCopyNodes(model, c(target, noncenteredEffects))
        
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)

        calcNodesNoSelf <- c(ccList$calcNodesNoSelf, noncenteredEffects)
        copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodes

        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory      <- c(0, 0)   ## scaleHistory
        acceptanceHistory <- c(0, 0)   ## scaleHistory
        saveMCMChistory <- getNimbleOption('MCMCsaveHistory')
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
        if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
        if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
                       propValue <- currentValue * exp(propLogScale)
        } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
        if(reflective) {
            lower <- model$getBound(target, 'lower')
            upper <- model$getBound(target, 'upper')
            while(propValue < lower | propValue > upper) {
                if(propValue < lower) propValue <- 2*lower - propValue
                if(propValue > upper) propValue <- 2*upper - propValue
            }
        }
        model[[target]] <<- propValue

        logMHR <- model$calculateDiff(target)
        if(logMHR == -Inf) {
            jump <- FALSE
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
        } else {
            ## Shift effects and add log-determinant of Jacobian of transformation. 
            logMHR <- logMHR + updateNoncentered(propValue, currentValue)
            logMHR <- logMHR + model$calculateDiff(calcNodesNoSelf) + propLogScale
            jump <- decide(logMHR)
            if(jump) {
                ##model$calculate(calcNodesPPomitted)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = noncenteredEffects, logProb = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = noncenteredEffects, logProb = TRUE)
            }
        }
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        updateNoncentered = function(newValue = double(), oldValue = double()) {
            if(noncenteredParam == 0) {
                values(model, noncenteredEffects) <<- values(model, noncenteredEffects) - oldValue + newValue
                return(0)
            } else {
                ## Mean will often be scalar; should we try to determine if this is the case to avoid work of `getParam` calls?
                for(i in 1:noncenteredLen) 
                    noncenteredMean[i] <<- model$getParam(noncenteredEffects[i], 'mean') 
                values(model, noncenteredEffects) <<- noncenteredMean + (newValue/oldValue) * (values(model, noncenteredEffects) - noncenteredMean)
                return(noncenteredLen * (log(newValue) - log(oldValue)))  # log determinant of Jacobian; accounts for fact we are computing prior for effects, which one doesn't in ASIS re-parameterized formulation. 
            }
            returnType(double())
        },
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)                 ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale                ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                }
                gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                ## If there are upper and lower bounds, enforce a maximum scale of
                ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
                ## Otherwise, for a poorly-informed posterior,
                ## the scale could grow without bound to try to reduce
                ## acceptance probability.  This creates enormous cost of
                ## reflections.
                if(reflective) {
                    lower <- model$getBound(target, 'lower')
                    upper <- model$getBound(target, 'upper')
                    if(scale >= 0.5*(upper-lower)) {
                        scale <<- 0.5*(upper-lower)
                    }
                }
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        setScale = function(newScale = double()) {
            scale         <<- newScale
            scaleOriginal <<- newScale
        },
        getScaleHistory = function() {       ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(scaleHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(acceptanceHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
            }
            gamma1 <<- 0
        }
    )
)

########################################################################
### block RW sampler with multi-variate normal proposal distribution ###
########################################################################

#' @rdname samplers
#' @export
sampler_RW_block <- nimbleFunction(
    name = 'sampler_RW_block',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',      FALSE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        propCov             <- extractControlElement(control, 'propCov',             'identity')
        tries               <- extractControlElement(control, 'tries',               1)
        maxDimCovHistory    <- extractControlElement(control, 'maxDimCovHistory',    10)
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
        finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
        if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in RW_block sampler')
        calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
        calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        d <- length(targetAsScalar)
        scaleHistory <- c(0, 0)                                                                  ## scaleHistory
        acceptanceHistory <- c(0, 0)                                                             ## scaleHistory
        propCovHistory <- if(d <= maxDimCovHistory) array(0, c(2,d,d)) else array(0, c(2,2,2))   ## scaleHistory
        saveMCMChistory <- if(getNimbleOption('MCMCsaveHistory')) TRUE else FALSE
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
        ## nested function and function list definitions
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
        ## checks
        if(any(model$isDiscrete(target)))       stop('cannot use RW_block sampler on discrete-valued target')
        if(!inherits(propCov, 'matrix'))        stop('propCov must be a matrix\n')
        if(!inherits(propCov[1,1], 'numeric'))  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))             stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))               stop('propCov matrix must be symmetric')
    },
    run = function() {
        for(i in 1:tries) {
            propValueVector <- generateProposalVector()
            values(model, targetNodesAsScalar) <<- propValueVector
            lpD <- model$calculateDiff(calcNodesProposalStage)
            if(lpD == -Inf) {
                jump <- FALSE
                nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
            } else {
                lpD <- lpD + model$calculateDiff(calcNodesDepStage)
                jump <- decide(lpD)
                if(jump) {
                    ##model$calculate(calcNodesPPomitted)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                } else {
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
                }
            }
            if(adaptive)     adaptiveProcedure(jump)
        }
    },
    methods = list(
        generateProposalVector = function() {
            propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)                 ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale                ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                    if(d <= maxDimCovHistory) {
                        propCovTemp <- propCovHistory                                           ## scaleHistory
                        setSize(propCovHistory, timesAdapted, d, d)                             ## scaleHistory
                        if(timesAdapted > 1)                                                    ## scaleHistory
                            for(iTA in 1:(timesAdapted-1))                                      ## scaleHistory
                                propCovHistory[iTA, 1:d, 1:d] <<- propCovTemp[iTA, 1:d, 1:d]    ## scaleHistory
                        propCovHistory[timesAdapted, 1:d, 1:d] <<- propCov[1:d, 1:d]            ## scaleHistory
                    }
                }
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$getGamma1()
                    for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                    empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                }
                chol_propCov_scale <<- chol_propCov * scale
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        setScale = function(newScale = double()) {
            scale         <<- newScale
            scaleOriginal <<- newScale
            chol_propCov_scale <<- chol_propCov * scale
        },
        setPropCov = function(newPropCov = double(2)) {
            propCov         <<- newPropCov
            propCovOriginal <<- newPropCov
            chol_propCov <<- chol(propCov)
            chol_propCov_scale <<- chol_propCov * scale
        },
        getScaleHistory = function() {  ## scaleHistory
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
            returnType(double(1))
            return(scaleHistory)
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
            return(acceptanceHistory)
        },                  
        getPropCovHistory = function() { ## scaleHistory
            if(!saveMCMChistory | d > maxDimCovHistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.  Note that to reduce memory use, proposal covariance histories are only saved for parameter vectors of length <= 10; this value can be modified using the 'maxDimCovHistory' control list element.")
            returnType(double(3))
            return(propCovHistory)
        },
        ##getScaleHistoryExpanded = function() {                                                              ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)                         ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]                      ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                                          ## scaleHistory
        ##getPropCovHistoryExpanded = function() {                                                            ## scaleHistory
        ##    propCovHistoryExpanded <- array(dim=c(timesAdapted*adaptInterval,d,d), init=FALSE)              ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
        ##            propCovHistoryExpanded[(iTA-1)*adaptInterval+j,1:d,1:d] <- propCovHistory[iTA,1:d,1:d]  ## scaleHistory
        ##    returnType(double(3)); return(propCovHistoryExpanded) },                                        ## scaleHistory
        reset = function() {
            scale   <<- scaleOriginal
            propCov <<- propCovOriginal
            chol_propCov <<- chol(propCov)
            chol_propCov_scale <<- chol_propCov * scale
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
                if(d <= maxDimCovHistory)
                    propCovHistory <<- nimArray(0, dim = c(2,d,d))
            }
            my_calcAdaptationFactor$reset()
        }
    )
)


#############################################################################
### RW_llFunction, does a RW, but using a generic log-likelihood function ###
#############################################################################

#' @rdname samplers
#' @export
sampler_RW_llFunction <- nimbleFunction(
    name = 'sampler_RW_llFunction',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- extractControlElement(control, 'adaptive',       TRUE)
        adaptInterval  <- extractControlElement(control, 'adaptInterval',  200)
        scale          <- extractControlElement(control, 'scale',          1)
        llFunction     <- extractControlElement(control, 'llFunction',     error = 'RW_llFunction sampler missing required control argument: llFunction')
        includesTarget <- extractControlElement(control, 'includesTarget', error = 'RW_llFunction sampler missing required control argument: includesTarget')
        ## nested function and function list definitions
        mvInternal <- modelValues(model)
        RWControl <- list(adaptive=adaptive, adaptInterval=adaptInterval, scale=scale, log=FALSE, reflective=FALSE)
        targetRWSamplerFunction <- sampler_RW(model, mvInternal, target, RWControl)
        my_setAndCalculateOne <- setAndCalculateOne(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, target = target)
    },
    run = function() {
        modelLP0 <- llFunction$run()
        if(!includesTarget)     modelLP0 <- modelLP0 + model$getLogProb(target)
        propValue <- rnorm(1, mean = model[[target]], sd = scale)
        my_setAndCalculateOne$run(propValue)
        modelLP1 <- llFunction$run()
        if(!includesTarget)     modelLP1 <- modelLP1 + model$getLogProb(target)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
        if(adaptive) {
            targetRWSamplerFunction$adaptiveProcedure(jump)
            scale <<- targetRWSamplerFunction$scale
        }
    },
    methods = list(
        reset = function() {
            targetRWSamplerFunction$reset()
        }
    )
)



####################################################################
### slice sampler (discrete or continuous) #########################
####################################################################

#' @rdname samplers
#' @export
sampler_slice <- nimbleFunction(
    name = 'sampler_slice',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive               <- extractControlElement(control, 'adaptive',               TRUE)
        adaptInterval          <- extractControlElement(control, 'adaptInterval',          200)
        width                  <- extractControlElement(control, 'sliceWidth',             1)
        maxSteps               <- extractControlElement(control, 'sliceMaxSteps',          100)
        maxContractions        <- extractControlElement(control, 'maxContractions',        1000)
        maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning', TRUE)
        eps <- 1e-15
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
        ## numeric value generation
        widthOriginal <- width
        timesRan      <- 0
        timesAdapted  <- 0
        sumJumps      <- 0
        widthHistory  <- c(0, 0)   ## widthHistory
        if(getNimbleOption('MCMCsaveHistory')) {
            saveMCMChistory <- TRUE
        } else saveMCMChistory <- FALSE
        discrete      <- model$isDiscrete(target)
        ## checks
        if(length(targetAsScalar) > 1)     stop('cannot use slice sampler on more than one target node')
    },
    run = function() {
        u <- model$getLogProb(calcNodes) - rexp(1, 1)    # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
        x0 <- model[[target]]    # create random interval (L,R), of width 'width', around current value of target
        L <- x0 - runif(1, 0, 1) * width
        R <- L + width
        maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
        maxStepsR <- maxSteps - 1 - maxStepsL
        lp <- setAndCalculateTarget(L)
        while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
            L <- L - width
            lp <- setAndCalculateTarget(L)
            maxStepsL <- maxStepsL - 1
        }
        lp <- setAndCalculateTarget(R)
        while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
            R <- R + width
            lp <- setAndCalculateTarget(R)
            maxStepsR <- maxStepsR - 1
        }
        x1 <- L + runif(1, 0, 1) * (R - L)
        lp <- setAndCalculateTarget(x1)
        numContractions <- 0
        while((is.nan(lp) | lp < u) & (R-L)/(abs(R)+abs(L)+eps) > eps & numContractions < maxContractions) {   # must be is.nan()
            ## The checks for R-L small and max number of contractions are for cases where model is in
            ## invalid state and lp calculations are NA/NaN or where R and L contract to each other
            ## division by R+L+eps ensures we check relative difference and that contracting to zero is ok
            if(x1 < x0) { L <- x1
                      } else      { R <- x1 }
            x1 <- L + runif(1, 0, 1) * (R - L)           # sample uniformly from (L,R) until sample is inside of slice (with shrinkage)
            lp <- setAndCalculateTarget(x1)
            numContractions <- numContractions + 1
        }
        if((R-L)/(abs(R)+abs(L)+eps) <= eps | numContractions == maxContractions) {
            if(maxContractionsWarning)
                cat("Warning: slice sampler reached maximum number of contractions for '", target, "'. Current parameter value is ", x0, ".\n")
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
            ##model$calculate(calcNodesPPomitted)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
            jumpDist <- abs(x1 - x0)
            if(adaptive)     adaptiveProcedure(jumpDist)
        }
    },
    methods = list(
        setAndCalculateTarget = function(value = double()) {
            if(discrete)     value <- floor(value)
            model[[target]] <<- value
            lp <- model$calculate(target)
            if(lp == -Inf) return(-Inf) 
            lp <- lp + model$calculate(calcNodesNoSelf)
            returnType(double())
            return(lp)
        },
        adaptiveProcedure = function(jumpDist = double()) {
            timesRan <<- timesRan + 1
            sumJumps <<- sumJumps + jumpDist   # cumulative (absolute) distance between consecutive values
            if(timesRan %% adaptInterval == 0) {
                adaptFactor <- (3/4) ^ timesAdapted
                meanJump <- sumJumps / timesRan
                width <<- width + (2*meanJump - width) * adaptFactor   # exponentially decaying adaptation of 'width' -> 2 * (avg. jump distance)
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(widthHistory, timesAdapted)                 ## widthHistory
                    widthHistory[timesAdapted] <<- width                ## widthHistory
                }
                timesRan <<- 0
                sumJumps <<- 0
            }
        },
        getWidthHistory = function() {       ## widthHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(widthHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },
        reset = function() {
            width        <<- widthOriginal
            timesRan     <<- 0
            timesAdapted <<- 0
            sumJumps     <<- 0
            if(saveMCMChistory) {
                widthHistory  <<- c(0, 0)    ## widthHistory
            }
        }
    )
)

## version for use with noncentered sampler

#' @rdname samplers
#' @export
sampler_slice_noncentered <- nimbleFunction(
    name = 'sampler_slice_noncentered',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive               <- extractControlElement(control, 'adaptive',               TRUE)
        adaptInterval          <- extractControlElement(control, 'adaptInterval',          200)
        width                  <- extractControlElement(control, 'sliceWidth',             1)
        maxSteps               <- extractControlElement(control, 'sliceMaxSteps',          100)
        maxContractions        <- extractControlElement(control, 'maxContractions',        1000)
        maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning', TRUE)
        eps <- 1e-15

        ## sampler_slice_noncentered is used by sampler_noncentered to do joint, noncentered sampling of
        ## target and dependent (random) effects, a version of the ASIS/interweaving sampler of Yu and Meng (2011).
        noncenteredEffects <- extractControlElement(control, 'noncenteredEffects', "")
        noncenteredParam <- extractControlElement(control, 'noncenteredParam', 0)

        noncenteredLen <- length(model$expandNodeNames(noncenteredEffects))
        noncenteredMean <- rep(0, noncenteredLen)
        noncenteredEffects0 <- rep(0, noncenteredLen)
        ccList <- mcmc_determineCalcAndCopyNodes(model, c(target, noncenteredEffects))

        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        
        calcNodes <- ccList$calcNodes; calcNodesNoSelf <- c(ccList$calcNodesNoSelf, noncenteredEffects); copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
        ## numeric value generation
        widthOriginal <- width
        timesRan      <- 0
        timesAdapted  <- 0
        sumJumps      <- 0
        widthHistory  <- c(0, 0)   ## widthHistory
        if(getNimbleOption('MCMCsaveHistory')) {
            saveMCMChistory <- TRUE
        } else saveMCMChistory <- FALSE
        discrete      <- model$isDiscrete(target)

        x0 <- 0

        ## checks
        if(length(targetAsScalar) > 1)     stop('cannot use slice sampler on more than one target node')
    },
    run = function() {
        u <- model$getLogProb(calcNodes) - rexp(1, 1)    # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
        x0 <<- model[[target]]    # create random interval (L,R), of width 'width', around current value of target
        noncenteredEffects0 <<- values(model, noncenteredEffects)
        L <- x0 - runif(1, 0, 1) * width
        R <- L + width
        maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
        maxStepsR <- maxSteps - 1 - maxStepsL
        lp <- setAndCalculateTarget(L)
        while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
            L <- L - width
            lp <- setAndCalculateTarget(L)
            maxStepsL <- maxStepsL - 1
        }
        lp <- setAndCalculateTarget(R)
        while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
            R <- R + width
            lp <- setAndCalculateTarget(R)
            maxStepsR <- maxStepsR - 1
        }
        x1 <- L + runif(1, 0, 1) * (R - L)
        lp <- setAndCalculateTarget(x1)
        numContractions <- 0
        while((is.nan(lp) | lp < u) & (R-L)/(abs(R)+abs(L)+eps) > eps & numContractions < maxContractions) {   # must be is.nan()
            ## The checks for R-L small and max number of contractions are for cases where model is in
            ## invalid state and lp calculations are NA/NaN or where R and L contract to each other
            ## division by R+L+eps ensures we check relative difference and that contracting to zero is ok
            if(x1 < x0) { L <- x1
                      } else      { R <- x1 }
            x1 <- L + runif(1, 0, 1) * (R - L)           # sample uniformly from (L,R) until sample is inside of slice (with shrinkage)
            lp <- setAndCalculateTarget(x1)
            numContractions <- numContractions + 1
        }
        if((R-L)/(abs(R)+abs(L)+eps) <= eps | numContractions == maxContractions) {
            if(maxContractionsWarning)
                cat("Warning: slice sampler reached maximum number of contractions for '", target, "'. Current parameter value is ", x0, ".\n")
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = noncenteredEffects, logProb = TRUE)
        } else {
            ##model$calculate(calcNodesPPomitted)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = noncenteredEffects, logProb = TRUE)
            jumpDist <- abs(x1 - x0)
            if(adaptive)     adaptiveProcedure(jumpDist)
        }
    },
    methods = list(
        setAndCalculateTarget = function(value = double()) {
            if(discrete)     value <- floor(value)
            model[[target]] <<- value
            lp <- model$calculate(target)
            if(lp == -Inf) return(-Inf)
            lp <- lp + updateNoncentered(value)
            lp <- lp + model$calculate(calcNodesNoSelf)
            returnType(double())
            return(lp)
        },
        updateNoncentered = function(newValue = double()) {
            if(noncenteredParam == 0) {
                values(model, noncenteredEffects) <<- noncenteredEffects0 - x0 + newValue
                return(0)
            } else {
                ## Mean will often be scalar; should we try to determine if this is the case to avoid work of `getParam` calls?
                for(i in 1:noncenteredLen) 
                    noncenteredMean[i] <<- model$getParam(noncenteredEffects[i], 'mean') 
                values(model, noncenteredEffects) <<- noncenteredMean + (newValue/x0) * (noncenteredEffects0 - noncenteredMean)
                return(noncenteredLen * (log(newValue) - log(x0)))  # log determinant of Jacobian; accounts for fact we are computing prior for effects, which one doesn't do in ASIS re-parameterized formulation. 
            }
            returnType(double())
        },
        adaptiveProcedure = function(jumpDist = double()) {
            timesRan <<- timesRan + 1
            sumJumps <<- sumJumps + jumpDist   # cumulative (absolute) distance between consecutive values
            if(timesRan %% adaptInterval == 0) {
                adaptFactor <- (3/4) ^ timesAdapted
                meanJump <- sumJumps / timesRan
                width <<- width + (2*meanJump - width) * adaptFactor   # exponentially decaying adaptation of 'width' -> 2 * (avg. jump distance)
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(widthHistory, timesAdapted)                 ## widthHistory
                    widthHistory[timesAdapted] <<- width                ## widthHistory
                }
                timesRan <<- 0
                sumJumps <<- 0
            }
        },
        getWidthHistory = function() {       ## widthHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(widthHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
                return(numeric(1, 0))
            }
        },
        reset = function() {
            width        <<- widthOriginal
            timesRan     <<- 0
            timesAdapted <<- 0
            sumJumps     <<- 0
            if(saveMCMChistory) {
                widthHistory  <<- c(0, 0)    ## widthHistory
            }
        }
    )
)


####################################################################
### elliptical slice sampler (dmnorm distributions) ################
####################################################################


essNFList_virtual <- nimbleFunctionVirtual(
    name = 'essNFList_virtual',
    run = function() { returnType(double(1)) }
)

essNF_univariate <- nimbleFunction(
    name = 'essNF_univariate',
    contains = essNFList_virtual,
    setup = function(model, node) {
        if(!(model$getDistribution(node) == 'dnorm'))   stop('something went wrong')
        mean <- numeric(2)
    },
    run = function() {
        setSize(mean, 1)
        mean[1] <<- model$getParam(node, 'mean')
        returnType(double(1))
        return(mean)
    }
)

essNF_multivariate <- nimbleFunction(
    name = 'essNF_multivariate',
    contains = essNFList_virtual,
    setup = function(model, node) {
        if(!(model$getDistribution(node) == 'dmnorm'))   stop('something went wrong')
    },
    run = function() {
        mean <- model$getParam(node, 'mean')
        returnType(double(1))
        return(mean)
    }
)

#' @rdname samplers
#' @export
sampler_ess <- nimbleFunction(
    name = 'sampler_ess',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        maxContractions        <- extractControlElement(control, 'maxContractions',        1000)
        maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning', TRUE)
        eps <- 1e-15
        ## node list generation
        target <- model$expandNodeNames(target)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodes
        ## numeric value generation
        Pi <- pi
        d <- length(model$expandNodeNames(target, returnScalarComponents = TRUE))
        d2 <- max(d, 2)
        target_mean <- f <- nu <- numeric(d2)
        ## nested function and function list definitions
        essNFList <- nimbleFunctionList(essNFList_virtual)
        if(model$getDistribution(target) == 'dnorm')    essNFList[[1]] <- essNF_univariate(model, target)
        if(model$getDistribution(target) == 'dmnorm')   essNFList[[1]] <- essNF_multivariate(model, target)
        ## checks
        if(length(target) > 1)                                           stop('elliptical slice sampler only applies to one target node')
        if(!(model$getDistribution(target) %in% c('dnorm', 'dmnorm')))   stop('elliptical slice sampler only applies to normal distributions')
    },
    run = function() {
        u <- model$getLogProb(calcNodesNoSelf) - rexp(1, 1)
        target_mean[1:d] <<- essNFList[[1]]$run()
        f[1:d] <<- values(model, target) - target_mean[1:d]
        model$simulate(target)
        nu[1:d] <<- values(model, target) - target_mean[1:d]
        theta <- runif(1, 0, 2*Pi)
        theta_min <- theta - 2*Pi
        theta_max <- theta
        values(model, target) <<- f[1:d]*cos(theta) + nu[1:d]*sin(theta) + target_mean[1:d]
        lp <- model$calculate(calcNodesNoSelf)
        numContractions <- 0
        while((is.nan(lp) | lp < u) & theta_max - theta_min > eps & numContractions < maxContractions) {   # must be is.nan()
            ## The checks for theta_max - theta_min small and max number of contractions are
            ## for cases where model is in invalid state and lp calculations are NA/NaN or where
            ## theta interval contracts to zero
            if(theta < 0)   theta_min <- theta   else   theta_max <- theta
            theta <- runif(1, theta_min, theta_max)
            values(model, target) <<- f[1:d]*cos(theta) + nu[1:d]*sin(theta) + target_mean[1:d]
            lp <- model$calculate(calcNodesNoSelf)
            numContractions <- numContractions + 1
        }
        if(theta_max - theta_min <= eps | numContractions == maxContractions) {
            if(maxContractionsWarning)
                cat("Warning: elliptical slice sampler reached maximum number of contractions.\n")
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
            ##model$calculate(calcNodesPPomitted)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        }
    },
    methods = list(
        reset = function() { }
    )
)



#####################################################################################
### Automated factor slice sampler (discrete or continuous) #########################
#####################################################################################

#' @rdname samplers
#' @export
sampler_AF_slice <- nimbleFunction(
    name = 'sampler_AF_slice',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        widthVec               <- extractControlElement(control, 'sliceWidths',              'oneVec')
        maxSteps               <- extractControlElement(control, 'sliceMaxSteps',            100)
        adaptFactorMaxIter     <- extractControlElement(control, 'sliceAdaptFactorMaxIter',  15000)
        adaptFactorInterval    <- extractControlElement(control, 'sliceAdaptFactorInterval', 200)
        adaptWidthMaxIter      <- extractControlElement(control, 'sliceAdaptWidthMaxIter',   512)
        adaptWidthTolerance    <- extractControlElement(control, 'sliceAdaptWidthTolerance', 0.1)
        maxContractions        <- extractControlElement(control, 'maxContractions',          1000)
        maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning',   TRUE)
        eps <- 1e-15
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
        finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
        if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in AF_slice sampler')
        calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
        calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
        ## numeric value generation
        d                  <- length(targetAsScalar)
        discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
        anyDiscrete        <- any(discrete)
        gammaMatrix        <- diag(d)         # matrix of orthogonal bases
        if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
        widthVecOriginal   <- widthVec
        nExpansions        <- rep(0, d)       # number of expansions
        nContracts         <- rep(0, d)       # number of contractions
        adaptFactorMaxIterOriginal <- adaptFactorMaxIter
        factorCounter      <- 0               # number of iterations since last factor adaptation
        factorTimesAdapted <- 0               # number of times factors have adapted
        empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
        empirCov           <- diag(d)
        allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
        widthCounter       <- 0               # number of iterations since last width adaptation
        adaptWidthMaxIterOriginal <- adaptWidthMaxIter
        adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
        widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
        ## checks
        if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
        if(!inherits(widthVec, 'numeric') && !inherits(widthVec, 'integer'))
            stop('sliceWidths must be a numeric vector')
        if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)
    },
    run = function() {
        maxContractionsReached <- FALSE
        for(i in 1:d) {
            eigenVec <- gammaMatrix[, i]
            width <- widthVec[i]
            u <- model$getLogProb(calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
            x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
            Lbound <- -1.0 * runif(1, 0, 1) * width
            Rbound <- Lbound + width
            L <- x0 + Lbound * eigenVec
            R <- x0 + Rbound * eigenVec
            maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
            maxStepsR <- maxSteps - 1 - maxStepsL
            lp <- setAndCalculateTarget(L)
            while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
                Lbound <- Lbound - width
                L <- x0 + Lbound * eigenVec
                lp <- setAndCalculateTarget(L)
                maxStepsL <- maxStepsL - 1
                nExpansions[i] <<- nExpansions[i] + 1
            }
            lp <- setAndCalculateTarget(R)
            while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
                Rbound <- Rbound + width
                R <- x0 + Rbound * eigenVec
                lp <- setAndCalculateTarget(R)
                maxStepsR <- maxStepsR - 1
                nExpansions[i] <<- nExpansions[i] + 1
            }
            prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
            x1 <- x0 + prop * eigenVec
            lp <- setAndCalculateTarget(x1)
            numContractions <- 0
            while((is.nan(lp) | lp < u) & Rbound - Lbound > eps & numContractions < maxContractions) {   # must be is.nan()
                ## The checks for Rbound - Lbound small and max number of contractions are
                ## for cases where model is in invalid state and lp calculations are NA/NaN or where
                ## interval contracts to zero
                if(prop < 0) { Lbound <- prop }
                else         { Rbound <- prop }
                nContracts[i] <<- nContracts[i] + 1
                prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
                x1 <- x0 + prop * eigenVec
                lp <- setAndCalculateTarget(x1)
                numContractions <- numContractions + 1
            }
            if(Rbound - Lbound <= eps | numContractions == maxContractions)
                maxContractionsReached <- TRUE
        }
        if(maxContractionsReached) {
            if(maxContractionsWarning)
                cat("Warning: AF slice sampler reached maximum number of contractions in at least one dimension.\n")
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
            ##model$calculate(calcNodesPPomitted)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        }
        if(allWidthsAdapted == 0)   adaptWidths()
        if(adaptFactorMaxIter > 0)  adaptFactors()
    },
    methods = list(
        setAndCalculateTarget = function(targetValues = double(1)) {
            if(anyDiscrete == 1)
                for(i in 1:d)
                    if(discrete[i] == 1)   targetValues[i] <- floor(targetValues[i])            
            values(model, target) <<- targetValues
            lp <- model$calculate(calcNodesProposalStage)
            if(lp == -Inf) return(lp)
            lp <- lp + model$calculate(calcNodesDepStage)
            returnType(double())
            return(lp)
        },
        adaptFactors = function() {
            adaptFactorMaxIter <<- adaptFactorMaxIter - 1
            factorCounter <<- factorCounter + 1
            empirSamp[factorCounter, 1:d] <<- values(model, target)
            if(factorCounter == adaptFactorInterval) {
                for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
                gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
                factorTimesAdapted <<- factorTimesAdapted + 1
                factorCounter      <<- 0
                nExpansions        <<- rep(0, d)
                nContracts         <<- rep(0, d)
                allWidthsAdapted   <<- 0
                widthCounter       <<- 0
                adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
                adaptWidthInterval <<- 1
                widthIndicatorVec  <<- rep(1, d)
            }
        },
        adaptWidths = function() {
            adaptWidthMaxIter <<- adaptWidthMaxIter - 1
            widthCounter <<- widthCounter + 1
            if(widthCounter == adaptWidthInterval) {
                for(i in 1:d) {
                    if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
                        if(nExpansions[i] == 0)   nExpansions[i] <<- 1
                        widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
                        widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
                        adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
                        nExpansions[i] <<- 0
                        nContracts[i] <<- 0
                        if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
                            widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
                    }
                }
                allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
                widthCounter     <<- 0
            }
            if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
                allWidthsAdapted <<- 1
        },
        reset = function() {
            gammaMatrix        <<- diag(d)
            empirCov           <<- diag(d)
            widthVec           <<- widthVecOriginal
            nExpansions        <<- rep(0, d)
            nContracts         <<- rep(0, d)
            adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
            factorCounter      <<- 0
            factorTimesAdapted <<- 0
            allWidthsAdapted   <<- 0
            widthCounter       <<- 0
            adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
            adaptWidthInterval <<- 1
            widthIndicatorVec  <<- rep(1, d)
        }
    )
)



####################################################################
### crossLevel sampler #############################################
####################################################################

getPosteriorDensityFromConjSampler_virtual <- nimbleFunctionVirtual(
    run = function() { returnType(double()) }
)

getPosteriorDensityFromConjSampler <- nimbleFunction(
    name = 'getPosteriorDensityFromConjSampler',
    contains = getPosteriorDensityFromConjSampler_virtual,
    setup = function(conjugateSamplerFunction) {},
    run = function() {
        posteriorLogDensity <- conjugateSamplerFunction$getPosteriorLogDensity()
        returnType(double())
        return(posteriorLogDensity)
    }
)

#' @rdname samplers
#' @export
sampler_crossLevel <- nimbleFunction(
    name = 'sampler_crossLevel',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive <- extractControlElement(control, 'adaptive', TRUE)
        ## node list generation
        target       <- model$expandNodeNames(target)
        lowNodes     <- model$getDependencies(target, self = FALSE, stochOnly = TRUE, includeData = FALSE)
        lowCalcNodes <- model$getDependencies(lowNodes)
        calcNodes    <- model$getDependencies(c(target, lowNodes))
        ## nested function and function list definitions
        mvInternal <- modelValues(model)
        topRWblockSamplerFunction <- sampler_RW_block(model, mvInternal, target, control)
        lowConjugateSamplerFunctions <- nimbleFunctionList(sampler_BASE)
        lowConjugateGetLogDensityFunctions <- nimbleFunctionList(getPosteriorDensityFromConjSampler_virtual)
        for(iLN in seq_along(lowNodes)) {
            lowNode <- lowNodes[iLN]
            conjugacyResult <- model$checkConjugacy(lowNode)[[lowNode]]
            if(is.null(conjugacyResult))     stop('non-conjugate lowNode \'', lowNode, '\' in crossLevel sampler')
            prior <- conjugacyResult$prior
            dependentCounts <- sapply(conjugacyResult$control, length)
            names(dependentCounts) <- gsub('^dep_', '', names(dependentCounts))
            conjSamplerName <- createDynamicConjugateSamplerName(prior = prior, dependentCounts = dependentCounts)
            if(!dynamicConjugateSamplerExists(conjSamplerName)) {
                conjSamplerDef <- conjugacyRelationshipsObject$generateDynamicConjugateSamplerDefinition(prior = prior, dependentCounts = dependentCounts)
                dynamicConjugateSamplerAdd(conjSamplerName, conjSamplerDef)
            }
            conjSamplerFunction <- dynamicConjugateSamplerGet(conjSamplerName)
            nameToPrint <- gsub('^sampler_', '', conjSamplerName)
            conjugateSamplerConf <- samplerConf(name = nameToPrint, samplerFunction = conjSamplerFunction, target = lowNode, control = conjugacyResult$control, model = model)
            lowConjugateSamplerFunctions[[iLN]] <- conjugateSamplerConf$buildSampler(model, mvInternal)
            lowConjugateGetLogDensityFunctions[[iLN]] <- getPosteriorDensityFromConjSampler(lowConjugateSamplerFunctions[[iLN]])
        }
        my_setAndCalculateTop <- setAndCalculate(model, target)
    },
    run = function() {
        modelLP0 <- model$getLogProb(calcNodes)
        propLP0 <- 0
        for(iSF in seq_along(lowConjugateGetLogDensityFunctions))  { propLP0 <- propLP0 + lowConjugateGetLogDensityFunctions[[iSF]]$run() }
        propValueVector <- topRWblockSamplerFunction$generateProposalVector()
        topLP <- my_setAndCalculateTop$run(propValueVector)
        if(is.na(topLP)) {
            logMHR <- -Inf
            jump <- decide(logMHR)
            if(jump) {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
            }
        }
        else {
            for(iSF in seq_along(lowConjugateSamplerFunctions))
                lowConjugateSamplerFunctions[[iSF]]$run()
            modelLP1 <- model$calculate(calcNodes)
            propLP1 <- 0
            for(iSF in seq_along(lowConjugateGetLogDensityFunctions))
                propLP1 <- propLP1 + lowConjugateGetLogDensityFunctions[[iSF]]$run()
            logMHR <- modelLP1 - modelLP0 - propLP1 + propLP0
            jump <- decide(logMHR)
            if(jump) {
                nimCopy(from = model,   to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodes, logProb = TRUE)
            }
    	}
        if(adaptive)     topRWblockSamplerFunction$adaptiveProcedure(jump)
    },
    methods = list(
        reset = function() {
            topRWblockSamplerFunction$reset()
            for(iSF in seq_along(lowConjugateSamplerFunctions)) {
                lowConjugateSamplerFunctions[[iSF]]$reset()
            }
        }
    )
)



########################################################################################
### RW_llFunctionBlock, does a block RW, but using a generic log-likelihood function ###
########################################################################################

#' @rdname samplers
#' @export
sampler_RW_llFunction_block <- nimbleFunction(
    name = 'sampler_RW_llFunction_block',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',      FALSE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        propCov             <- extractControlElement(control, 'propCov',             'identity')
        llFunction          <- extractControlElement(control, 'llFunction',          error = 'RW_llFunction_block sampler missing required control argument: llFunction')
        includesTarget      <- extractControlElement(control, 'includesTarget',      error = 'RW_llFunction_block sampler missing required control argument: includesTarget')
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        d <- length(targetAsScalar)
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
        ## nested function and function list definitions
        my_setAndCalculate <- setAndCalculate(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, target = target)
        my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
        ## checks
        if(!inherits(propCov, 'matrix'))        stop('propCov must be a matrix\n')
        if(!inherits(propCov[1,1], 'numeric'))  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    },
    run = function() {
        modelLP0 <- llFunction$run()
        if(!includesTarget)     modelLP0 <- modelLP0 + model$getLogProb(target)
        propValueVector <- generateProposalVector()
        my_setAndCalculate$run(propValueVector)
        modelLP1 <- llFunction$run()
        if(!includesTarget)     modelLP1 <- modelLP1 + model$getLogProb(target)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        generateProposalVector = function() {
            propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$getGamma1()
                    for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                    empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                }
                chol_propCov_scale <<- chol_propCov * scale
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        reset = function() {
            scale   <<- scaleOriginal
            propCov <<- propCovOriginal
            chol_propCov <<- chol(propCov)
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            my_calcAdaptationFactor$reset()
        }
    )
)

#########################################################################################
##### RW_multinomial sampler for multinomial distributions ##############################
#########################################################################################
##
## @rdname samplers
## @export
##sampler_RW_multinomial <- nimbleFunction(
##    name = 'sampler_RW_multinomial',
##    contains = sampler_BASE,
##    setup = function(model, mvSaved, target, control) {
##        ## control list extraction
##        adaptive      <- extractControlElement(control, 'adaptive',      TRUE)
##        adaptInterval <- extractControlElement(control, 'adaptInterval', 200)
##        ## node list generation
##        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
##        targetAllNodes <- unique(model$expandNodeNames(target))
##        calcNodes      <- model$getDependencies(target)
##        lTarget        <- length(targetAsScalar)
##        Ntotal         <- sum(values(model,target))
##        NOverL         <- Ntotal / lTarget
##        ## numeric value generation
##        propVector        <- rep(0, lTarget)
##        Zeros             <- matrix(0, lTarget, lTarget)
##        Ones              <- matrix(1, lTarget, lTarget)
##        timesRan          <- Zeros
##        AcceptRates       <- Zeros
##        ScaleShifts       <- Zeros
##        totalAdapted      <- Zeros
##        timesAccepted     <- Zeros
##        ENSwapMatrix      <- Ones
##        ENSwapDeltaMatrix <- Ones
##        RescaleThreshold  <- 0.2 * Ones
##        ## nested function and function list definitions
##        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
##        my_decideAndJump       <- decideAndJump(model, mvSaved, target = target)
##        ## checks
##        if(model$getDistribution(target) != 'dmulti')   stop('can only use RW_multinomial sampler for multinomial distributions')
##        if(length(targetAllNodes) > 1)                  stop('cannot use RW_multinomial sampler on more than one target')
##        if(adaptive & adaptInterval < 100)              stop('adaptInterval < 100 is not recommended for RW_multinomial sampler')
##    },
##    run = function() {
##        for(iFROM in 1:lTarget) {
##            for(iTO in 1:(lTarget-1)) {
##                if(runif(1,0,1) > 0.5) {
##                    iFrom <- iFROM
##                    iTo   <- iTO
##                    if(iFrom == iTo)   iTo <- lTarget
##                } else {
##                    iFrom <- iTO
##                    iTo   <- iFROM
##                    if(iFrom == iTo)   iFrom <- lTarget
##                }
##                ## generate proposal vector
##                propVector <<- values(model,target)
##                pSwap       <- min(1, max(1, ENSwapMatrix[iFrom,iTo]) / propVector[iFrom])
##                nSwap       <- rbinom(n=1,   size=propVector[iFrom], prob=pSwap)
##                lpProp      <- dbinom(nSwap, size=propVector[iFrom], prob=pSwap, log=TRUE)
##                propVector[iFrom] <<- propVector[iFrom] - nSwap
##                propVector[iTo]   <<- propVector[iTo]   + nSwap
##                pRevSwap    <- min(1, max(1, ENSwapMatrix[iTo,iFrom]) / propVector[iTo])
##                lpRev       <- dbinom(nSwap, size=propVector[iTo], prob=pRevSwap, log=TRUE)
##                ## decide and jump
##                lpMHR <- my_setAndCalculateDiff$run(propVector) + lpRev - lpProp
##                jump  <- my_decideAndJump$run(lpMHR, 0, 0, 0)
##                ## adaptation
##                if(adaptive)   adaptiveProcedure(jump=jump, iFrom=iFrom, iTo=iTo)
##            }
##        }
##    },
##    methods = list(
##        adaptiveProcedure = function(jump=logical(), iFrom=integer(), iTo=integer()) {
##            timesRan[iFrom, iTo] <<- timesRan[iFrom, iTo] + 1
##            if(jump)
##                timesAccepted[iFrom, iTo] <<- timesAccepted[iFrom, iTo] + 1
##            if(timesRan[iFrom, iTo] %% adaptInterval == 0) {
##                totalAdapted[iFrom, iTo] <<- totalAdapted[iFrom, iTo] + 1
##                accRate                   <- timesAccepted[iFrom, iTo] / timesRan[iFrom, iTo]
##                AcceptRates[iFrom, iTo]  <<- accRate
##                if(accRate > 0.5) {
##                    ENSwapMatrix[iFrom, iTo] <<- min(Ntotal, ENSwapMatrix[iFrom,iTo] + ENSwapDeltaMatrix[iFrom, iTo] / totalAdapted[iFrom,iTo])
##                } else {
##                    ENSwapMatrix[iFrom, iTo] <<- max(1, ENSwapMatrix[iFrom,iTo] - ENSwapDeltaMatrix[iFrom,iTo] / totalAdapted[iFrom,iTo])
##                }
##                if(accRate<RescaleThreshold[iFrom,iTo] | (accRate > (1-RescaleThreshold[iFrom,iTo]))) {
##                    ## rescale iff ENSwapMatrix[iFrom, iTo] is not set to an upper or lower bound
##                    if(ENSwapMatrix[iFrom, iTo] > 1 & ENSwapMatrix[iFrom, iTo] < Ntotal) {
##                        ScaleShifts[iFrom, iTo]       <<- ScaleShifts[iFrom, iTo] + 1
##                        ENSwapDeltaMatrix[iFrom, iTo] <<- min(NOverL, ENSwapDeltaMatrix[iFrom, iTo] * totalAdapted[iFrom,iTo] / 10)
##                        ENSwapDeltaMatrix[iTo, iFrom] <<- ENSwapDeltaMatrix[iFrom, iTo]
##                        RescaleThreshold[iFrom,iTo]   <<- 0.2 * 0.95^ScaleShifts[iFrom, iTo]
##                    }
##                }
##                ## lower Bound
##                if(ENSwapMatrix[iFrom, iTo] < 1)   ENSwapMatrix[iFrom, iTo] <<- 1
##                ## symmetry in ENSwapMatrix helps maintain good acceptance rates
##                ENSwapMatrix[iTo,iFrom]   <<- ENSwapMatrix[iFrom,iTo]
##                timesRan[iFrom, iTo]      <<- 0
##                timesAccepted[iFrom, iTo] <<- 0
##            }
##        },
##        reset = function() {
##            timesRan          <<- Zeros
##            AcceptRates       <<- Zeros
##            ScaleShifts       <<- Zeros
##            totalAdapted      <<- Zeros
##            timesAccepted     <<- Zeros
##            ENSwapMatrix      <<- Ones
##            ENSwapDeltaMatrix <<- Ones
##            RescaleThreshold  <<- 0.2 * Ones
##        }
##    )
##)



#####################################################################################
### RW_dirichlet sampler for dirichlet distributions ################################
#####################################################################################

#' @rdname samplers
#' @export
sampler_RW_dirichlet <- nimbleFunction(
    name = 'sampler_RW_dirichlet',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scaleOriginal       <- extractControlElement(control, 'scale',               1)
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch  # not used: calcNodes
        ## numeric value generation
        d <- length(targetAsScalar)
        thetaVec         <- rep(0, d)
        scaleVec         <- rep(scaleOriginal, d)
        timesRan         <- 0
        timesAcceptedVec <- rep(0, d)
        timesAdapted     <- 0
        optimalAR        <- 0.44
        gamma1           <- 0
        ## checks
        if(length(model$expandNodeNames(target)) > 1)    stop('RW_dirichlet sampler only applies to one target node')
        if(model$getDistribution(target) != 'ddirch')    stop('can only use RW_dirichlet sampler for dirichlet distributions')
    },
    run = function() {
        if(thetaVec[1] == 0)   thetaVec <<- values(model, target)   ## initialization
        alphaVec <- model$getParam(target, 'alpha')
        for(i in 1:d) {
            currentValue <- thetaVec[i]
            propLogScale <- rnorm(1, mean = 0, sd = scaleVec[i])
            propValue <- currentValue * exp(propLogScale)
            if(propValue != 0) {
                thetaVecProp <- thetaVec
                thetaVecProp[i] <- propValue
                values(model, target) <<- thetaVecProp / sum(thetaVecProp)
                logMHR <- alphaVec[i]*propLogScale + currentValue - propValue + model$calculateDiff(calcNodesNoSelf)
                jump <- decide(logMHR)
            } else jump <- FALSE
            if(adaptive & jump)   timesAcceptedVec[i] <<- timesAcceptedVec[i] + 1
            if(jump) {
                thetaVec <<- thetaVecProp
                ##model$calculate(calcNodesPPomitted)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
            }
            model$calculate(target)                                                             ## update target logProb
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProbOnly = TRUE)    ##
        }
        if(adaptive) {
            timesRan <<- timesRan + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRateVec <- timesAcceptedVec / timesRan
                timesAdapted <<- timesAdapted + 1
                gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
                adaptFactorVec <- exp(10 * gamma1 * (acceptanceRateVec - optimalAR))
                scaleVec <<- scaleVec * adaptFactorVec
                timesRan <<- 0
                timesAcceptedVec <<- numeric(d, 0)
            }
        }
    },
    methods = list(
        reset = function() {
            thetaVec         <<- numeric(d, 0)
            scaleVec         <<- numeric(d, scaleOriginal)
            timesRan         <<- 0
            timesAcceptedVec <<- numeric(d, 0)
            timesAdapted     <<- 0
            gamma1           <<- 0
        }
    )
)


######################################################################################
### RW_wishart sampler for non-conjugate Wishart and inverse-Wishart distributions ###
######################################################################################

#' @rdname samplers
#' @export
sampler_RW_wishart <- nimbleFunction(
    name = 'sampler_RW_wishart',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        propCov             <- extractControlElement(control, 'propCov',             'identity')
        ## node list generation
        target <- model$expandNodeNames(target)
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch  # not used: calcNodesNoSelf
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        d <- sqrt(length(targetAsScalar))
        nTheta <- d*(d+1)/2
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(nTheta)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=nTheta)
        thetaVec      <- numeric(nTheta)
        thetaVec_prop <- numeric(nTheta)
        currentValue      <- array(0, c(d,d))
        currentValue_chol <- array(0, c(d,d))
        propValue         <- array(0, c(d,d))
        propValue_chol    <- array(0, c(d,d))
        ## nested function and function list definitions
        my_calcAdaptationFactor <- calcAdaptationFactor(nTheta, adaptFactorExponent)
        ## checks
        dist <- model$getDistribution(target)
        if(d < 2)                             stop('RW_wishart sampler requires target node dimension to be at least 2x2')
        if(!inherits(propCov, 'matrix'))        stop('propCov must be a matrix')
        if(!inherits(propCov[1,1], 'numeric'))  stop('propCov matrix must be numeric')
        if(!all(dim(propCov) == nTheta))      stop('propCov matrix must have dimension ', d, 'x', d)
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    },
    run = function() {
        currentValue <<- model[[target]]
        currentValue_chol <<- chol(currentValue)
        ## calculate thetaVec
        for(i in 1:d)   thetaVec[i] <<- log(currentValue_chol[i,i])
        nextInd <- d + 1
        for(i in 1:(d-1)) {
            for(j in (i+1):d) {
                thetaVec[nextInd] <<- currentValue_chol[i,j]
                nextInd <- nextInd + 1
            }
        }
        ## generate thetaVec proposal on transformed scale
        thetaVec_prop <<- rmnorm_chol(1, thetaVec, chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
        ## un-transform thetaVec_prop to get propValue_chol
        for(i in 1:d)   propValue_chol[i,i] <<- exp(thetaVec_prop[i])
        nextInd <- d + 1
        for(i in 1:(d-1)) {
            for(j in (i+1):d) {
                propValue_chol[i,j] <<- thetaVec_prop[nextInd]
                nextInd <- nextInd + 1
            }
        }
        ## matrix multiply to get proposal value (matrix)
        model[[target]] <<- t(propValue_chol) %*% propValue_chol
        ## decide and jump
        logMHR <- model$calculateDiff(calcNodes)
        deltaDiag <- thetaVec_prop[1:d]-thetaVec[1:d]
        for(i in 1:d)   logMHR <- logMHR + (d+2-i)*deltaDiag[i]  ## took me quite a while to derive this
        jump <- decide(logMHR)
        if(jump) {
            ##model$calculate(calcNodesPPomitted)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
        }
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!jump)    empirSamp[timesRan, 1:nTheta] <<- thetaVec
            if(jump)     empirSamp[timesRan, 1:nTheta] <<- thetaVec_prop
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                gamma1 <- my_calcAdaptationFactor$getGamma1()
                for(i in 1:nTheta)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
                propCov <<- propCov + gamma1 * (empirCov - propCov)
                chol_propCov <<- chol(propCov)
                chol_propCov_scale <<- chol_propCov * scale
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        reset = function() {
            scale   <<- scaleOriginal
            propCov <<- propCovOriginal
            chol_propCov <<- chol(propCov)
            chol_propCov_scale <<- chol_propCov * scale
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            my_calcAdaptationFactor$reset()
        }
    )
)

######################################################################################
### RW_lkj_corr_cholesky sampler for non-conjugate LKJ correlation Chlesky factor  ###
######################################################################################

#' @rdname samplers
#' @export
sampler_RW_lkj_corr_cholesky <- nimbleFunction(
    name = 'sampler_RW_lkj_corr_cholesky',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scaleOriginal       <- extractControlElement(control, 'scale',               1)
        ## node list generation
        target <- model$expandNodeNames(target)
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)  
        ## numeric value generation
        d <- sqrt(length(targetAsScalar))
        nTheta <- d*(d-1)/2    # this many unconstrained elements to be sampled
        if(length(adaptInterval) > 1 || length(adaptFactorExponent) > 1 || length(scaleOriginal) > 1)
            stop("RW_lkj_corr_cholesky: 'adaptInterval', 'adaptFactorExponent', and 'scaleOriginal' should be single values.")
        if(scaleOriginal < 0)
            stop('Cannot use RW_lkj_corr_cholesky sampler with scale control parameter less than 0.')
        ## adaptation objects
        if(nTheta == 1) nTheta <- 2
        scaleVec            <- rep(scaleOriginal, nTheta)
        timesRan            <- 0
        timesAcceptedVec    <- rep(0, nTheta)
        timesAdapted        <- 0
        optimalAR           <- 0.44
        ##
        z                   <- array(0, c(d, d))  # canonical partial correlations
        diag(z)             <- 1
        partialSums         <- array(0, c(d, d))  # 1-x_{13}^2, 1-x_{14}^2, 1-x_{14}^2-x_{24}^2, ...
        partialSums[1, ]    <- 1
        ## temporary vectors for current column calculations
        partialSumsProp     <- numeric(d)    
        currentValue        <- numeric(d)
        propValue           <- numeric(d)
        ## checks
        dist <- model$getDistribution(target)
        if(dist != 'dlkj_corr_cholesky') stop('RW_lkj_corr_cholesky sampler can only be used with the dlkj_corr_cholesky distribution.')
        if(d < 2)                        stop('RW_lkj_corr_cholesky sampler requires target node dimension to be at least 2x2.')
        if(adaptFactorExponent < 0)      stop('Cannot use RW_lkj_corr_cholesky sampler with adaptFactorExponent control parameter less than 0.')
    },
    run = function() {
        ## calculate transformed values (in unconstrained space) and partial sums in each column
        transform(model[[target]])  # compute z and partialSums
        ## Individual univariate RW on the nTheta elements:
        ## advantage: probably better movement through space than with block update
        ## (plus note only column values of target matrix are recalculated for a given scalar update in a given column)
        ## disadvantage: all downstream dependencies of entire matrix will be recalculated, so much slower than default block sampler
        cnt <- 0
        for(i in 2:d) {   # Iterate over columns
            currentValue <<- model[[target]][1:d, i]
            partialSumsProp <<- partialSums[, i]
            for(j in 1:(i-1)) {
                cnt <- cnt + 1
                ## RW on unconstrained y
                yCurrent <- atanh(z[j, i])
                yProp <- rnorm(1, yCurrent, scaleVec[cnt])
                zProp <- tanh(yProp)   
                propValue[j] <<- zProp * sqrt(partialSumsProp[j])  # proposed value of element of U
                ## Update remainder of column of U so that length of column is 1
                for(jprime in (j+1):i) {
                    partialSumsProp[jprime] <<- partialSumsProp[jprime-1] - propValue[jprime-1]^2
                    propValue[jprime] <<- z[jprime, i] * sqrt(partialSumsProp[jprime])
                }
                model[[target]][j:i, i] <<- propValue[j:i] 
                logMHR <- calculateDiff(model, calcNodesNoSelf) + calculateDiff(model, target)
                ## Adjust MHR to account for non-symmetric proposal by adjusting prior on U to transformed scale (i.e., y).
                ## cosh component is for dz/dy and other component is for du/dz  where 'u' is the corr matrix.
                ## This follows Stan reference manual Section 10.12 (for version 2.27).
                ## There is an important subtlety here - the transformation is from y to U not y to Sigma, so the Jacobian
                ## omits the dSigma/dU, which is in Section 10.9 of the Stan reference manual.
                ## This is because dlkj_corr_cholesky is written directly in terms of a density on U and that density
                ## already accounts for the Jacobian in going from Sigma to U.
                logMHR <- logMHR - 2*(log(cosh(yProp)) - log(cosh(yCurrent))) 
                if(j < i-1) 
                    logMHR <- logMHR + 0.5*sum(log(partialSumsProp[(j+1):(i-1)]) - log(partialSums[(j+1):(i-1), i]))
 
                jump <- decide(logMHR)
                ## Avoid copying entire target matrix as we are modifying one column at a time.
                if(jump) {
                    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelf, logProb = TRUE)
                    timesAcceptedVec[cnt] <<- timesAcceptedVec[cnt] + 1
                    partialSums[(j+1):i, i] <<- partialSumsProp[(j+1):i]
                    currentValue[j:i] <<- propValue[j:i]
                } else {
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelf, logProb = TRUE)
                    model[[target]][j:i, i] <<- currentValue[j:i]
                    model$calculate(target)   ## Update target logProb since not part of calcNodesNoSelf.
                    partialSumsProp[j+1] <<- partialSums[j+1, i]  ## Remaining elements will be modified above anyway.
                }
            }
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
        if(adaptive) {
            timesRan <<- timesRan + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRateVec <- timesAcceptedVec / timesRan
                timesAdapted <<- timesAdapted + 1
                gamma1 <- 1/((timesAdapted + 3)^adaptFactorExponent)
                adaptFactorVec <- exp(10 * gamma1 * (acceptanceRateVec - optimalAR))   
                scaleVec <<- scaleVec * adaptFactorVec
                timesRan <<- 0
                timesAcceptedVec <<- numeric(nTheta, 0)
            }
        }
    },
    methods = list(
        transform = function(x = double(2)) {
            ## Calculate canonical partial correlations (z) and partial sums (remaining lengths of U columns)
            z[1, 2:d] <<- x[1, 2:d]
            partialSums[2, 2] <<- 1 - x[1, 2]^2
            if(d >= 3) {
                for(i in 3:d) {
                    for(j in 2:(i-1)) {
                        partialSums[j, i] <<- partialSums[j-1, i] - x[j-1, i]^2
                        z[j, i] <<- x[j, i] / sqrt(partialSums[j, i])
                    }
                    partialSums[i, i] <<- partialSums[i-1, i] - x[i-1, i]^2
                }
            }
        },
        reset = function() {
            scaleVec         <<- numeric(nTheta, scaleOriginal)
            timesRan         <<- 0
            timesAdapted     <<- 0
            timesAcceptedVec <<- numeric(nTheta, 0)
        }
    )
)

############################################################################################
### RW_block_lkj_corr_cholesky sampler for non-conjugate LKJ correlation Chlesky factor  ###
############################################################################################


#' @rdname samplers
#' @export
sampler_RW_block_lkj_corr_cholesky <- nimbleFunction(
    name = 'sampler_RW_block_lkj_corr_cholesky',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',      FALSE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        propCov             <- extractControlElement(control, 'propCov',             'identity')
        tries               <- extractControlElement(control, 'tries',               1)
        maxDimCovHistory    <- extractControlElement(control, 'maxDimCovHistory',    10)
        ## node list generation
        target <- model$expandNodeNames(target)
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
        if(!is.integer(finalTargetIndex) |
           length(finalTargetIndex) != 1 |
           is.na(finalTargetIndex[1]))
           stop('Problem with target node in sampler_RW_block_lkj_corr_cholesky')
        calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
        calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
##        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        p <- sqrt(length(targetAsScalar))
        d <- p*(p-1)/2    # this many unconstrained elements to be sampled
        scaleHistory <- c(0, 0)                                                                  ## scaleHistory
        acceptanceHistory <- c(0, 0)                                                             ## scaleHistory
        propCovHistory <- if(d <= maxDimCovHistory) array(0, c(2,d,d)) else array(0, c(2,2,2))   ## scaleHistory
        saveMCMChistory <- if(getNimbleOption('MCMCsaveHistory')) TRUE else FALSE
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
        ## nested function and function list definitions
        ##        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)

        z                   <- array(0, c(p, p))  # canonical partial correlations
        diag(z)             <- 1
        y                   <- numeric(d)         # unconstrained parameters; sampling occurs on this scale
        partialSums         <- array(0, c(p, p))  # 1-x_{21}^2, 1-x_{31}^2, 1-x_{31}^2-x_{32}^2, ...
        partialSums[1, ]    <- 1
        logDetJac           <- 0                  # for log of determinant of Jacobian for initial value
        
        ## checks
        if(!inherits(propCov, 'matrix'))       stop('propCov must be a matrix\n')
        if(!inherits(propCov[1,1], 'numeric')) stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))            stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))              stop('propCov matrix must be symmetric')

        dist <- model$getDistribution(target)
        if(dist != 'dlkj_corr_cholesky') stop('RW_block_lkj_corr_cholesky sampler can only be used with the dlkj_corr_cholesky distribution.')
        if(d < 3)                        stop('RW_block_lkj_corr_cholesky sampler requires target node dimension to be at least 3x3.')
        if(adaptFactorExponent < 0)      stop('Cannot use RW_block_lkj_corr_cholesky sampler with adaptFactorExponent control parameter less than 0.')

    },
    run = function() {
        transform(model[[target]])  # compute z and partialSums
        logMHR <- 0
        for(k in 1:tries) {
            yPropValueVector <- generateProposalVector()
            ##        lpMHR <- my_setAndCalculateDiff$run(propValueVector)
            cnt <- 1
            for(i in 2:p) {
                model[[target]][1, i] <<- tanh(yPropValueVector[cnt])
                partialSumsProp <- 1 - model[[target]][1, i]^2
                ## Adjust MHR to account for non-symmetric proposal by adjusting prior on U to transformed scale (i.e., y).
                ## cosh component is for dz/dy and other component is for du/dz  where 'u' is the corr matrix.
                ## This follows Stan reference manual Section 10.12 (for version 2.27).
                ## There is an important subtlety here - the transformation is from y to U not y to Sigma, so the Jacobian
                ## omits the dSigma/dU (which is in Section 10.9 of the Stan reference manual).
                ## This is because dlkj_corr_cholesky is written directly in terms of a density on U and that density
                ## already accounts for the Jacobian in going from Sigma to U.
                logMHR <- logMHR - 2*log(cosh(yPropValueVector[cnt])) 
                cnt <- cnt+1
                if(i > 2) {
                    for(j in 2:(i-1)) {
                        ## Inverse transform from y to U
                        model[[target]][j, i] <<- tanh(yPropValueVector[cnt]) * sqrt(partialSumsProp) 
                        ## Adjust for prior for y using determinant of y to x transformation.
                        logMHR <- logMHR - 2*log(cosh(yPropValueVector[cnt])) +
                            0.5*log(partialSumsProp)
                        cnt <- cnt+1
                        partialSumsProp <- partialSumsProp - model[[target]][j, i]^2
                    }
                }
                model[[target]][i, i] <<- sqrt(partialSumsProp)
            }
            ## Adjust for log determinant term from initial values
            logMHR <- logMHR - logDetJac

            lpD <- calculateDiff(model, calcNodesProposalStage)
            if(lpD == -Inf) {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
                jump <- FALSE
            } else {
                logMHR <- logMHR + lpD + calculateDiff(model, calcNodesDepStage)
                jump <- decide(logMHR)
                if(jump) {
                    nimCopy(from = model,   to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
                    y <<- yPropValueVector
                } else {
                    nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodes, logProb = TRUE)
                }
            }
            if(adaptive)     adaptiveProcedure(jump)
        }
    },
    methods = list(
        transform = function(x = double(2)) {
            ## Compute canonical partial correlations (z), transformed parameters (y), and
            ## remaining lengths of U columns (partialSums).
            z[1, 2:p] <<- x[1, 2:p]
            y[1] <<- atanh(z[1, 2])
            logDetJac <<- -2*log(cosh(y[1]))
            cnt <- 2
            partialSums[2, 2] <<- 1 - x[1, 2]^2
            for(i in 3:p) {
                y[cnt] <<- atanh(z[1, i])
                logDetJac <<- logDetJac - 2*log(cosh(y[cnt]))
                cnt <- cnt+1
                partialSums[2, i] <<- 1 - x[1, i]^2
                for(j in 2:(i-1)) {
                    z[j, i] <<- x[j, i] / sqrt(partialSums[j, i])
                    y[cnt] <<- atanh(z[j, i])
                    logDetJac <<- logDetJac - 2*log(cosh(y[cnt])) + 0.5*log(partialSums[j, i]) 
                    cnt <- cnt+1
                    partialSums[j+1, i] <<- partialSums[j, i] - x[j, i]^2
                }
            }
        },
        generateProposalVector = function() {
            propValueVector <- rmnorm_chol(1, y, chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- y
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)                 ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale                ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                    if(d <= maxDimCovHistory) {
                        propCovTemp <- propCovHistory                                           ## scaleHistory
                        setSize(propCovHistory, timesAdapted, d, d)                             ## scaleHistory
                        if(timesAdapted > 1)                                                    ## scaleHistory
                            for(iTA in 1:(timesAdapted-1))                                      ## scaleHistory
                                propCovHistory[iTA, 1:d, 1:d] <<- propCovTemp[iTA, 1:d, 1:d]    ## scaleHistory
                        propCovHistory[timesAdapted, 1:d, 1:d] <<- propCov[1:d, 1:d]            ## scaleHistory
                    }
                }
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$getGamma1()
                    for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                    empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                }
                chol_propCov_scale <<- chol_propCov * scale
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        getScaleHistory = function() {  ## scaleHistory
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
            returnType(double(1))
            return(scaleHistory)
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.")
            return(acceptanceHistory)
        },                  
        getPropCovHistory = function() { ## scaleHistory
            if(!saveMCMChistory | d > maxDimCovHistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC.  Note that to reduce memory use, proposal covariance histories are only saved for parameter vectors of length <= 10; this value can be modified using the 'maxDimCovHistory' control list element.")
            returnType(double(3))
            return(propCovHistory)
        },
        ##getScaleHistoryExpanded = function() {                                                              ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)                         ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]                      ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                                          ## scaleHistory
        ##getPropCovHistoryExpanded = function() {                                                            ## scaleHistory
        ##    propCovHistoryExpanded <- array(dim=c(timesAdapted*adaptInterval,d,d), init=FALSE)              ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                                      ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                                   ## scaleHistory
        ##            propCovHistoryExpanded[(iTA-1)*adaptInterval+j,1:d,1:d] <- propCovHistory[iTA,1:d,1:d]  ## scaleHistory
        ##    returnType(double(3)); return(propCovHistoryExpanded) },                                        ## scaleHistory
        reset = function() {
            scale   <<- scaleOriginal
            propCov <<- propCovOriginal
            chol_propCov <<- chol(propCov)
            chol_propCov_scale <<- chol_propCov * scale
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
                if(d <= maxDimCovHistory)
                    propCovHistory <<- nimArray(0, dim = c(2,d,d))
            }
            my_calcAdaptationFactor$reset()
        }
    )
)

####################################################################################################
### CAR_scalar_{postPred,conjugate,RW} component samplers for dcar_{normal,proper} distributions ###
####################################################################################################


## posterior predictive sampler for scalar components of dcar_normal() and dcar_proper() nodes,
## i.e. those scalar components which have no stochastic dependencies
CAR_scalar_postPred <- nimbleFunction(
    name = 'CAR_scalar_postPred',
    contains = sampler_BASE,
    setup = function(model, mvSaved, targetScalar, neighborNodes, neighborWeights, Mi, proper) {
        ## node list generation
        calcNodes <- c(targetScalar, model$getDependencies(targetScalar, self = FALSE))
        ## nested function and function list definitions
        dcarList <- nimbleFunctionList(CAR_evaluateDensity_base)
        if(proper) { dcarList[[1]] <- CAR_proper_evaluateDensity(model, targetScalar, neighborNodes, neighborWeights, Mi)
                 } else { dcarList[[1]] <- CAR_normal_evaluateDensity(model, targetScalar, neighborNodes, neighborWeights) }
        ## checks
        targetDCAR <- model$expandNodeNames(targetScalar)
        if(length(targetDCAR) != 1)                                                     stop('something went wrong')
        if(!(model$getDistribution(targetDCAR) %in% c('dcar_normal', 'dcar_proper')))   stop('something went wrong')
    },
    run = function() {
        if(is.na(model[[targetScalar]]) | is.nan(model[[targetScalar]])) {
            model[[targetScalar]] <<- 0
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        priorPrec <- dcarList[[1]]$getPrec()
        if(priorPrec == 0) return()    ## dcar_normal with 0 neighbors, or dcar_proper with 0 neighbors and auto-generated M
        newValue <- rnorm(1, mean = dcarList[[1]]$getMean(), sd = sqrt(1/priorPrec))
        model[[targetScalar]] <<- newValue
        model$calculate(calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    )
)


## dnorm-dnorm conjugate sampler for scalar components of dcar_normal() and dcar_proper() nodes,
## i.e. those scalar components which only have dnorm() dependencies
CAR_scalar_conjugate <- nimbleFunction(
    name = 'CAR_scalar_conjugate',
    contains = sampler_BASE,
    setup = function(model, mvSaved, targetScalar, neighborNodes, neighborWeights, Mi, proper) {
        ## node list generation
        calcNodes <- c(targetScalar, model$getDependencies(targetScalar, self = FALSE))
        calcNodesDeterm <- model$getDependencies(targetScalar, determOnly = TRUE)
        depStochNodes_dnorm <- model$getDependencies(targetScalar, self = FALSE, stochOnly = TRUE)
        nDependents <- length(depStochNodes_dnorm)
        if(nDependents < 1) stop('something went wrong')
        ## numeric value generation
        if(!all(model$getDistribution(depStochNodes_dnorm) == 'dnorm')) stop('something went wrong')
        linearityCheckExprList <- lapply(depStochNodes_dnorm, function(node) model$getParamExpr(node, 'mean'))
        linearityCheckExprList <- lapply(linearityCheckExprList, function(expr) cc_expandDetermNodesInExpr(model, expr, skipExpansionsNode=model$expandNodeNames(targetScalar)))
        if(!all(sapply(linearityCheckExprList, function(expr) cc_nodeInExpr(targetScalar, expr)))) stop('something went wrong')
        linearityCheckResultList <- lapply(linearityCheckExprList, function(expr) cc_checkLinearity(expr, targetScalar))
        if(any(sapply(linearityCheckResultList, function(expr) is.null(expr)))) stop('something went wrong')
        offsetList <- lapply(linearityCheckResultList, '[[', 'offset')
        scaleList <- lapply(linearityCheckResultList, '[[', 'scale')
        allIdentityLinks <- ifelse(all(sapply(offsetList, function(offset) offset == 0)) && sapply(scaleList, function(scale) scale == 1), 1, 0)
        ##if(allIdentityLinks)  cat(paste0(targetScalar, ' conjugate sampler: skipping two-point method\n'))
        ##if(!allIdentityLinks) cat(paste0(targetScalar, ' conjugate sampler: require two-point method\n'))
        ## nested function and function list definitions
        dcarList <- nimbleFunctionList(CAR_evaluateDensity_base)
        if(proper) { dcarList[[1]] <- CAR_proper_evaluateDensity(model, targetScalar, neighborNodes, neighborWeights, Mi)
                 } else { dcarList[[1]] <- CAR_normal_evaluateDensity(model, targetScalar, neighborNodes, neighborWeights) }
        ## checks
        targetDCAR <- model$expandNodeNames(targetScalar)
        if(length(targetDCAR) != 1)                                                     stop('something went wrong')
        if(!(model$getDistribution(targetDCAR) %in% c('dcar_normal', 'dcar_proper')))   stop('something went wrong')
    },
    run = function() {
        prior_mean <- dcarList[[1]]$getMean()
        prior_tau <- dcarList[[1]]$getPrec()
        dependent_values <- values(model, depStochNodes_dnorm)
        dependent_taus <- numeric(nDependents)
        for(i in 1:nDependents)
            dependent_taus[i] <- model$getParam(depStochNodes_dnorm[i], 'tau')
        if(allIdentityLinks) {   ## don't need to calculate coeff and offset
            contribution_mean <- sum(dependent_values * dependent_taus)
            contribution_tau <- sum(dependent_taus)
        } else {                 ## use "two-point" method to calculate coeff and offset
            dependent_offsets <- numeric(nDependents)
            dependent_coeffs <- numeric(nDependents)
            model[[targetScalar]] <<- 0
            model$calculate(calcNodesDeterm)
            for(i in 1:nDependents)
                dependent_offsets[i] <- model$getParam(depStochNodes_dnorm[i], 'mean')
            model[[targetScalar]] <<- 1
            model$calculate(calcNodesDeterm)
            for(i in 1:nDependents)
                dependent_coeffs[i] <- model$getParam(depStochNodes_dnorm[i], 'mean') - dependent_offsets[i]
            contribution_mean <- sum(dependent_coeffs * (dependent_values - dependent_offsets) * dependent_taus)
            contribution_tau <- sum(dependent_coeffs^2 * dependent_taus)
        }
        ## from MCMC_conjugacy.R:
        ## dnorm  = list(param = 'mean', contribution_mean = 'coeff * (value-offset) * tau', contribution_tau = 'coeff^2 * tau'),
        ## posterior = 'dnorm(mean = (prior_mean*prior_tau + contribution_mean) / (prior_tau + contribution_tau),
        ##                    sd   = (prior_tau + contribution_tau)^(-0.5))'
        posterior_mean <- (prior_mean * prior_tau + contribution_mean) / (prior_tau + contribution_tau)
        posterior_sd <- (prior_tau + contribution_tau)^(-0.5)
        newValue <- rnorm(1, mean = posterior_mean, sd = posterior_sd)
        model[[targetScalar]] <<- newValue
        model$calculate(calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    )
)


## RW sampler for non-conjugate scalar components of dcar_normal() and dcar_proper() nodes
CAR_scalar_RW <- nimbleFunction(
    name = 'CAR_scalar_RW',
    contains = sampler_BASE,
    setup = function(model, mvSaved, targetScalar, neighborNodes, neighborWeights, Mi, control, proper) {
        ## control list extraction
        adaptive      <- extractControlElement(control, 'adaptive',      TRUE)
        adaptInterval <- extractControlElement(control, 'adaptInterval', 200)
        scale         <- extractControlElement(control, 'scale',         1)
        ## node list generation
        depNodes <- model$getDependencies(targetScalar, self = FALSE)
        copyNodes <- c(targetScalar, depNodes)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        optimalAR     <- 0.44
        gamma1        <- 0
        ## nested function and function list definitions
        dcarList <- nimbleFunctionList(CAR_evaluateDensity_base)
        if(proper) { dcarList[[1]] <- CAR_proper_evaluateDensity(model, targetScalar, neighborNodes, neighborWeights, Mi)
                 } else { dcarList[[1]] <- CAR_normal_evaluateDensity(model, targetScalar, neighborNodes, neighborWeights) }
        ## checks
        targetDCAR <- model$expandNodeNames(targetScalar)
        if(length(targetDCAR) != 1)                                                     stop('something went wrong')
        if(!(model$getDistribution(targetDCAR) %in% c('dcar_normal', 'dcar_proper')))   stop('something went wrong')
    },
    run = function() {
        if(is.na(model[[targetScalar]]) | is.nan(model[[targetScalar]])) {
            model[[targetScalar]] <<- 0
            model$calculate(copyNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodes, logProb = TRUE)
        }
        lp0 <- dcarList[[1]]$run() + model$getLogProb(depNodes)
        propValue <- rnorm(1, mean = model[[targetScalar]], sd = scale)
        model[[targetScalar]] <<- propValue
        lp1 <- dcarList[[1]]$run() + model$calculate(depNodes)
        logMHR <- lp1 - lp0
        jump <- decide(logMHR)
        if(jump) {
            model$calculate(targetScalar)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodes, logProb = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodes, logProb = TRUE)
        }
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            gamma1 <<- 0
        }
    )
)


#################################################################################################
### CAR_normal sampler for intrinsic conditionally autoregressive (dcar_normal) distributions  ###
#################################################################################################


#' @rdname samplers
#' @export
sampler_CAR_normal <- nimbleFunction(
    name = 'sampler_CAR_normal',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        useConjugacy  <- extractControlElement(control, 'carUseConjugacy', TRUE)
        ## node list generation
        target <- model$expandNodeNames(target)
        targetScalarComponents <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        ## numeric value generation
        zero_mean_caused_a_problem <- 0
        ## checks
        if(length(target) > 1)                               stop('CAR_normal sampler only applies to one target node')
        if(model$getDistribution(target) != 'dcar_normal')   stop('CAR_normal sampler only applies to dcar_normal distributions')
        ## nested function and function list definitions
        adj <- model$getParam(target, 'adj')
        weights <- model$getParam(target, 'weights')
        num <- model$getParam(target, 'num')
        zero_mean <- model$getParam(target, 'zero_mean')
        neighborLists <- CAR_normal_processParams(model, target, adj, weights, num)
        componentSamplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(targetScalarComponents)) {
            targetScalar <- targetScalarComponents[i]
            nDependents <- length(model$getDependencies(targetScalar, self = FALSE, stochOnly = TRUE))
            conjugate <- CAR_checkConjugacy(model, targetScalar, carNode = target)
            neighborNodes <- neighborLists$neighborNodeList[[i]]
            neighborWeights <- neighborLists$neighborWeightList[[i]]
            ##if(length(neighborNodes)==0) cat(paste0('island node detected: ', targetScalar, '\n'))
            if(nDependents == 0) {
                ##cat(paste0('dcar() component node ', targetScalar, ': assigning posterior predictive sampler\n'))
                componentSamplerFunctions[[i]] <- CAR_scalar_postPred(model, mvSaved, targetScalar, neighborNodes, neighborWeights, Mi = 0, proper = FALSE)
            } else if(conjugate && useConjugacy) {
                ##cat(paste0('dcar() component node ', targetScalar, ': assigning conjugate sampler\n'))
                componentSamplerFunctions[[i]] <- CAR_scalar_conjugate(model, mvSaved, targetScalar, neighborNodes, neighborWeights, Mi = 0, proper = FALSE)
            } else {
                ##cat(paste0('dcar() component node ', targetScalar, ': assigning RW sampler\n'))
                componentSamplerFunctions[[i]] <- CAR_scalar_RW(model, mvSaved, targetScalar, neighborNodes, neighborWeights, Mi = 0, control = control, proper = FALSE)
            }
        }
    },
    run = function() {
        for(iSF in seq_along(componentSamplerFunctions))
            componentSamplerFunctions[[iSF]]$run()
        if(zero_mean) {
            targetValues <- values(model, target)
            values(model, target) <<- targetValues - mean(targetValues)
            lp <- model$calculate(calcNodes)
            if(is.na(lp) | is.nan(lp) | abs(lp) == Inf)   zero_mean_caused_a_problem <<- 1
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list(
        after_chain = function() {
            if(zero_mean_caused_a_problem == 1) {
                print('  [Error] Invalid model values were induced during the MCMC, which were caused by centering of the CAR (dcar_normal) process.  Centering takes place because of the argument zero_mean=1 to the dcar_normal distribution.  Results of this MCMC run may be invalid. This can be avoided by absorbing the mean into the CAR process, and omitting the zero_mean argument.')
            }
        },
        reset = function() {
            zero_mean_caused_a_problem <<- 0
            for(iSF in seq_along(componentSamplerFunctions))
                componentSamplerFunctions[[iSF]]$reset()
        }
    )
)


##############################################################################################
### CAR_proper sampler for proper conditionally autoregressive (dcar_proper) distributions  ###
##############################################################################################


#' @rdname samplers
#' @export
sampler_CAR_proper <- nimbleFunction(
    name = 'sampler_CAR_proper',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        useConjugacy  <- extractControlElement(control, 'carUseConjugacy', TRUE)
        ## node list generation
        target <- model$expandNodeNames(target)
        targetScalarComponents <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        ## checks
        if(length(target) > 1)                               stop('CAR_proper sampler only applies to one target node')
        if(model$getDistribution(target) != 'dcar_proper')   stop('CAR_proper sampler only applies to dcar_proper distributions')
        ## nested function and function list definitions
        C <- model$getParam(target, 'C')
        adj <- model$getParam(target, 'adj')
        num <- model$getParam(target, 'num')
        M <- model$getParam(target, 'M')
        neighborLists <- CAR_proper_processParams(model, target, C, adj, num, M)
        componentSamplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(targetScalarComponents)) {
            targetScalar <- targetScalarComponents[i]
            nDependents <- length(model$getDependencies(targetScalar, self = FALSE, stochOnly = TRUE))
            conjugate <- CAR_checkConjugacy(model, targetScalar, carNode = target)
            neighborNodes <- neighborLists$neighborNodeList[[i]]
            neighborCs <- neighborLists$neighborCList[[i]]
            Mi <- M[i]
            ##if(length(neighborNodes)==0) cat(paste0('island node detected: ', targetScalar, '\n'))
            if(nDependents == 0) {
                ##cat(paste0('dcar() component node ', targetScalar, ': assigning posterior predictive sampler\n'))
                componentSamplerFunctions[[i]] <- CAR_scalar_postPred(model, mvSaved, targetScalar, neighborNodes, neighborCs, Mi = Mi, proper = TRUE)
            } else if(conjugate && useConjugacy) {
                ##cat(paste0('dcar() component node ', targetScalar, ': assigning conjugate sampler\n'))
                componentSamplerFunctions[[i]] <- CAR_scalar_conjugate(model, mvSaved, targetScalar, neighborNodes, neighborCs, Mi = Mi, proper = TRUE)
            } else {
                ##cat(paste0('dcar() component node ', targetScalar, ': assigning RW sampler\n'))
                componentSamplerFunctions[[i]] <- CAR_scalar_RW(model, mvSaved, targetScalar, neighborNodes, neighborCs, Mi = Mi, control = control, proper = TRUE)
            }
        }
    },
    run = function() {
        for(iSF in seq_along(componentSamplerFunctions))
            componentSamplerFunctions[[iSF]]$run()
    },
    methods = list(
        reset = function() {
            for(iSF in seq_along(componentSamplerFunctions))
                componentSamplerFunctions[[iSF]]$reset()
        }
    )
)

################################################################################################################
### polyagamma conjugate sampler for parameters appearing additively in logistic regression linear predictor ###
################################################################################################################

## nimbleFunction implementation of Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013.
##   Forecasting High-Dimensional, Time-Varying Variance-Covariance Matrices
##   with High-Frequency Data and Sampling Pólya-Gamma Random Variates for
##   Posterior Distributions Derived from Logistic Likelihoods  
## Adapted from C++ code in R package `pgdraw`.
samplePolyaGamma <- nimbleFunction(
    setup = function(){
        logpi <- log(pi)
        pisq <- pi*pi
        t <- 2/pi ## Best choice of t is near 0.64
        w <- sqrt(pi/2)
        pisq_8 <- pisq/8
        pi_2 <- pi/2
    },
    run = function(){},
    methods = list(
        rpolyagamma = function(size = double(1), y = double(1)) {
            m <- length(size)
            n <- length(y)
            sample <- numeric(value = 0, length = n)
            thisSize <- size[1]

            if(m > 1) {
                for(i in 1:n) {
                    thisSize <- size[i]
                    if(abs(y[i]) < Inf & thisSize > 0) {  ## Default zero if y = Inf or -Inf or if size = 0.
                        ## Note: check on y[i] not strictly needed if we only call this for non-zero-inflated cases, which is now the situation.
                        for (j in 1:thisSize)
                            sample[i] <- sample[i] + rpolyagammaOne(y[i])
                    }
                }
            } else {
                if(thisSize > 0) {  ## Avoid check of size inside loop.
                    for(i in 1:n) {
                        if(abs(y[i]) < Inf) {  ## Default zero if y = Inf or -Inf
                            ## Note: check on y[i] not strictly needed if we only call this for non-zero-inflated cases, which is now the situation.
                            for (j in 1:thisSize)
                                sample[i] <- sample[i] + rpolyagammaOne(y[i])
                        }
                    }
                }
            }
            returnType(double(1))
            return(sample)	
        },
        rpolyagammaOne = function(y = double(0)) {
            z <- 0.5*abs(y)
            K <- z*z*0.5 + pisq_8
            logA <- log(4) - logpi - z
            logK <- log(K)
            Kt <- K * t
            logf1 <- logA + pnorm(w*(t*z - 1), 0, 1, 1, 1) + logK + Kt
            logf2 <- logA + 2*z + pnorm(-w*(t*z+1), 0, 1, 1, 1) + logK + Kt
            p_over_q <- exp(logf1) + exp(logf2)
            ratio <- 1 / (1 + p_over_q)
            
            done1 <- FALSE
            while(!done1){
                u <- runif(1)
                if(u < ratio){
                    sample <- t + rexp(1, 1)/K
                } else {
                    sample <- rtinvgauss(z)
                }
                
                i <- 1
                Sn <- aterm(0, sample)
                U <- runif(1) * Sn
                asgn <- -1
                even <- -1
                
                done2 <- FALSE
                while(!done2) {
                    Sn <- Sn + asgn * aterm(i, sample)
                    ## if odd
                    if(even < 0 & U <= Sn) {
                        sample <- sample * 0.25
                        done2 <- TRUE
                        done1 <- TRUE
                    }
                    ## if even
                    if(even > 0 & U > Sn)
                        done2 <- TRUE
                    even <- -even
                    asgn <- -asgn
                    i <- i + 1
                }
            }
            returnType(double())
            return(sample)
        },
        aterm = function(n = integer(), x = double()) {
            f <- logpi + log(n + 0.5)
            if(x <= t) {
                f <- f + 1.5*(log(2) - logpi - log(x)) - 2*(n + 0.5)*(n + 0.5)/x
            } else {
                f <- f - 0.5 * x * pisq * (n + 0.5)*(n + 0.5)
            }   
            returnType(double())
            return(exp(f))			
        },
        rtinvgauss = function(z = double()) {
            mu <- 1/z
            done <- FALSE
            if(mu > t) {
                while(!done) {
                    u <- runif(1)
                    sample <- 1 / rtruncgamma()
                    if (log(u) < (-z*z*0.5*sample)) done <- TRUE
                }
            } else {
                sample <- t + 1
                while(sample >= t) {
                    sample <- rinvgaussian(mu)
                }
            }
            returnType(double())
            return(sample)
        },
        rtruncgamma = function() {
            done <- FALSE
            while(!done) {
                sample <- rexp(1, 1) * 2 + pi_2
                gX <- w / sqrt(sample)
                if(runif(1) <= gX) done <- TRUE
            }
            returnType(double())
            return(sample)
        },
        rinvgaussian = function(mu = double()) {
            u <- rnorm(1, 0, 1)
            V <- u*u
            sample <- mu + 0.5*mu * ( mu*V - sqrt(4*mu*V + mu*mu * V*V) )
            if(runif(1) > mu /(mu+sample)) {    
                sample <- mu*mu / sample 
            }   
            returnType(double())
            return(sample)
        }
    )
)


## Tooling for dealing with allowing both dnorm and dmnorm nodes.
## Copied from INLA work; eventually consider handling as general tool.

getParam_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
      getMean = function(index = integer()) { returnType(double(1)) },
      getPrecision = function(index = integer()) { returnType(double(2)) }
  )
)

## A place holder to not take up much memory.
emptyParam <- nimbleFunction(
    contains = getParam_BASE,
    setup = function() {},
    run = function() {},
    methods = list(
        getPrecision = function(index = integer()){
            returnType(double(2))
            return(matrix(1, nrow = 1, ncol = 1))
        },
        getMean = function(index = integer()){
            returnType(double(1))
            return(numeric(1, length = 1))
        }
    )
)

## Need at least one dnorm to use this.
## NodeNames relate to node names in the model that are dnorm distributed
## gNodes indicates a 1 if dnorm, 0 o/w. 
## This makes it easy to get the correct indices when pass in the node index in a loop.
gaussParam <- nimbleFunction(
    contains = getParam_BASE,
    setup = function(model, nodeNames, gNodes) {
        indexConvert <- cumsum(gNodes)
        if(length(indexConvert) == 1)
            indexConvert <- c(indexConvert, -1)
    },
    run = function() {},
    methods = list(
        getPrecision = function(index = integer()){
            i <- indexConvert[index]
            return(matrix(model$getParam(nodeNames[i], "tau"), nrow = 1, ncol = 1))
            returnType(double(2))
        },
        getMean = function(index = integer()){
            i <- indexConvert[index]
            return(numeric(model$getParam(nodeNames[i], "mean"), length = 1))
            returnType(double(1))
        }
    )
)

## Need at least one dmnorm to use this.
multiGaussParam <- nimbleFunction(
    contains = getParam_BASE,
    setup = function(model, nodeNames, gNodes) {
        indexConvert <- cumsum(gNodes)
        if(length(indexConvert) == 1)
            indexConvert <- c(indexConvert, -1)
    },
    run = function(){},
    methods = list(
        getPrecision = function(index = integer()){
            i <- indexConvert[index]
            return(model$getParam(nodeNames[i], "prec"))
            returnType(double(2))
        },
        getMean = function(index = integer()){
            i <- indexConvert[index]
            return(model$getParam(nodeNames[i], "mean"))
            returnType(double(1))
        }
    )
)

####################################################################
### Polya-Gamma Data Augmentation ##################################
####################################################################

sampler_polyagamma <- nimbleFunction(
    name = 'sampler_polyagamma',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        
        ## This allows user to override error trapping for unusual model structures.
        check <- extractControlElement(control, 'check', TRUE)   

        ## node list generation
        target <- model$expandNodeNames(target)
        targetDists <- model$getDistribution(target)
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        nonTarget <- model$expandNodeNames(control$nonTargetNodes)
        ccList <- mcmc_determineCalcAndCopyNodes(model, target)
        calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf;
        copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch

        nTarget <- length(target)
        nCoef <- length(targetAsScalar)
     #   if(nCoef == 1)  ## We could/should relax this, though not clear how common it would be.
     #       stop("polyagamma sampler not set up to handle a scalar target node")
        nodeLengths <- sapply(target, function(x) length(model$expandNodeNames(x, returnScalarComponents = TRUE)))
        
        
        if(is.null(control$response)) {
            ## We don't need the user to provide the response and prefer they not, but leaving this flexibility for now.
            yNodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE)
        } else {
            yNodes <- model$expandNodeNames(control$response)
        }
        N <- length(yNodes)

        checkMessage <- "If your model is in a non-standard form and you are sure the Pólya-gamma sampler is appropriate, you can disable this check by setting the control argument `check=FALSE`."

        ## Conjugacy checking, part 1.
        if(check) {
            if(!all(targetDists %in% c("dnorm", "dmnorm")))
                stop("polyagamma sampler: all target nodes must have `dnorm` or `dmnorm` priors. ", checkMessage)
            if(!all(model$getDistribution(yNodes) %in% c("dbern", "dbin")) ) 
                stop("polyagamma sampler: response nodes must be distributed `dbern` or `dbin`. ", checkMessage)
            nodeIDs <- model$expandNodeNames(yNodes, returnType = 'ids')
            if(length(unique(model$modelDef$maps$graphID_2_declID[nodeIDs])) > 1)  # So that we can do conj checking only on one item.
                stop("polyagamma sampler: response nodes should all be part of the same declaration. ", checkMessage)
        }
            
        probAndSizeNodes <- model$getParents(yNodes, immediateOnly = TRUE)
        depNodes <- model$getDependencies(target, self = TRUE)
        probNodes <- intersect(probAndSizeNodes, depNodes)  # These nodes may reflect zero-inflated probabilities.
        sizeNodes <- setdiff(probAndSizeNodes, probNodes)

        zeroInflated <- FALSE
        
        ## Conjugacy checking, part 2.
        ## Make sure any stochastic dependencies between target and y are Bernoulli (i.e. only zero-inflation allowed)
        ## and that zero-inflation variable multiplies the baseline probability.

        ## First we need some processing to make sure that we can simply check inflation based only on `probNodes[1]`,
        ## to avoid costly checking.
        if(check) {
            inflationNodes <- model$getParents(probNodes, omit = c(target, nonTarget), stochOnly = TRUE)
            if(length(inflationNodes)) {
                ## Check that inflation probabilities are directly specified as parents of `probNodes`
                ## to avoid having to check multiple declarations. Seemingly anything otherwise would be
                ## an unusual zero inflation construction.
                inflationNodes <- setdiff(inflationNodes, nonTarget)
                test <- model$getParents(probNodes, omit = c(target, nonTarget), stochOnly = TRUE, immediateOnly = TRUE)
                test <- setdiff(test, nonTarget)
                if(!identical(test, inflationNodes))  # So we need to only consider a single declaration.
                    stop("polyagamma sampler: zero-inflation probabilities should be specified directly as Bernoulli or binomial (with `size=1`) random variables in the response node declaration to enable NIMBLE to efficiently check model validity. ", checkMessage)
                nodeIDs <- model$expandNodeNames(probNodes, returnType = 'ids')
                if(length(unique(model$modelDef$maps$graphID_2_declID[nodeIDs])) > 1)  # If declaration of zero-inflation nodes occurs in one declaration, we can do conj checking only on one item.
                    stop("polyagamma sampler: zero-inflation nodes should all be part of the same declaration to enable NIMBLE to efficiently check model validity. ", checkMessage)
            }
        }
        
        inflationStochNodesOne <- model$getParents(probNodes[1], omit = c(target, nonTarget), stochOnly = TRUE, self = FALSE)
        ## We ask user to provide non-target nodes in the linear predictor as otherwise hard to distinguish from zero-inflation nodes.
        if(length(inflationStochNodesOne)) {
            zeroInflated <- TRUE
            ones <- rep(1, length(model$expandNodeNames(inflationNodes, returnScalarComponents = TRUE)))
            inflationNodesDeps <- model$getDependencies(inflationNodes, determOnly = TRUE, self = FALSE)
            dists <- model$getDistribution(inflationStochNodesOne)

            if(check) {
                if(!all(dists %in% c("dbern", "dbin")))
                    stop("polyagamma sampler: Invalid stochastic nodes found as parents of response. Any such nodes other than the target must specify zero inflation, and any non-target nodes in the linear predictor must be included in `control$nonTargetNodes`. ", checkMessage)
                binomDists <- dists == 'dbin'
                if(any(binomDists)) {
                    if(!all(sapply(inflationStochNodesOne[binomDists], function(x) model$getParamExpr(x, 'size') == 1)))
                        stop("polyagamma sampler: Zero inflation nodes must be `dbern` or `dbin` with `size=1`. ", checkMessage)
                }
            }

            probNodesInflated <- probNodes
            probNodes <- intersect(model$getParents(probNodes, self = FALSE, immediateOnly = TRUE), depNodes)

            ## Check probability is product of inflationNodes and non-inflated probability.
            linearityCheckExprRaw <- model$getValueExpr(probNodesInflated[1])
            for(node in inflationStochNodesOne) {
                linearityCheckExpr <- cc_expandDetermNodesInExpr(model, linearityCheckExprRaw, targetNode = node)
                linearityCheck  <- cc_checkLinearity(linearityCheckExpr, node)
                linkCheck <- cc_linkCheck(linearityCheck, 'multiplicative')
                if(check && (is.null(linkCheck) || linkCheck != 'multiplicative'))
                    stop("polyagamma sampler: with zero inflation, probability must be specified as the product of one or more Bernoulli random variables and the expit-transformed linear predictor. ", checkMessage)
            }
            linearityCheck  <- cc_checkLinearity(linearityCheckExprRaw, probNodes[1])
            linkCheck <- cc_linkCheck(linearityCheck, 'multiplicative')
            if(check && (is.null(linkCheck) || linkCheck != 'multiplicative'))
                stop("polyagamma sampler: with zero inflation, probability must be specified as the product of one or more Bernoulli random variables and the expit-transformed linear predictor. ", checkMessage)
        } else {
            ## Placeholders to allow compilation.
            ones <- rep(1, 2)
            inflationNodes <- probNodes[1]
            inflationNodesDeps <- probNodes[1]
            probNodesInflated <- probNodes[1]
        }
        

        ## At this point, `probNodes` has the nodes for the non-inflated probabilities.

        ## Conjugacy checking, part 3: Check linearity of target nodes in logit link.
        if(check) {
            ## In order to do conjugacy checking only on one item for efficiency, we need all `probNodes` and all
            ## nodes declared in the sequence leading to the target nodes declared in single declarations.
            ## These checks see if anything more complicated is going on.
            nodeIDs <- model$expandNodeNames(probNodes, returnType = 'ids')
            if(length(unique(model$modelDef$maps$graphID_2_declID[nodeIDs])) > 1)  
                stop("polyagamma sampler: probabilities for all response nodes should all be constructed in the same declaration to enable NIMBLE to efficiently check model validity. ", checkMessage)
            nodesToCheck <- model$getParents(probNodes, immediateOnly = TRUE)
            while(!any(target %in% nodesToCheck)) {
                nodeIDs <- model$expandNodeNames(nodesToCheck, returnType = 'ids')
                if(length(unique(model$modelDef$maps$graphID_2_declID[nodeIDs])) > 1)  
                    stop("polyagamma sampler: linear predictors should all be constructed in the same declaration to enable NIMBLE to efficiently check model validity. ", checkMessage)
                nodesToCheck <- model$getParents(nodesToCheck, immediateOnly = TRUE)
            }
             
            if(model$getValueExpr(probNodes[1])[[1]] != 'expit')  
                stop("polyagamma sampler: target must be related to response via logit link. Also note that zero inflation cannot be specified directly in the declaration for the linear predictor to enable NIMBLE to efficiently check model validity. ", checkMessage)   ## `z[i]*expit(b0+b1*x[i])` would be harder to check for validity.
            linearityCheckExprRaw <- model$getValueExpr(probNodes[1])[[2]]
            for(node in targetAsScalar) {
                linearityCheckExpr <- cc_expandDetermNodesInExpr(model, linearityCheckExprRaw, targetNode = node)
                linearityCheck  <- cc_checkLinearity(linearityCheckExpr, node)
                linkCheck <- cc_linkCheck(linearityCheck, "linear")
                if(is.null(linkCheck) || !linkCheck %in% c('identity', 'additive', 'multiplicative', 'linear'))
                    stop("polyagamma sampler: probability must be specified (via logit link) as a linear function of the target nodes. ", checkMessage)
            }
        }
      
        stochSize <- FALSE
        if(length(sizeNodes) && length(model$getParents(sizeNodes, stochOnly = TRUE, self = TRUE)))  
            stochSize <- TRUE  ## This assumes any RHSonly nodes will not change.

        singleSize <- FALSE

        dnormNodes <- targetDists == "dnorm"
        dmnormNodes <- targetDists == "dmnorm"
        n_dnorm <- sum(dnormNodes)
        n_dmnorm <- sum(dmnormNodes)
        
        ## Build nimble function list allowing for both dnorm and dmnorm, as required for compilation to work.
        getParam_nfl <- nimbleFunctionList(getParam_BASE)
        if(n_dnorm > 0) {
            getParam_nfl[[1]] <- gaussParam(model, target[dnormNodes], dnormNodes)
        } else {
            getParam_nfl[[1]] <- emptyParam()
        }
        if(n_dmnorm > 0) {
            getParam_nfl[[2]] <- multiGaussParam(model, target[dmnormNodes], dmnormNodes)
        } else{ 
            getParam_nfl[[2]] <- emptyParam()
        }
        normTypes <- dnormNodes + 2 * dmnormNodes   # vector of values in {1,2} indicating dnorm or dmnorm
        
        ## Build design matrix, which account for all effects (fixed and random) in target.
        if(is.null(control$designMatrix)) {
            X <- matrix(0, nrow = N, ncol = nCoef)
            fixed <- FALSE
            ## Do more on inferring fixed columns?
            ## Generally anything except effects that are stochastically indexed
            ## or where covariates are random (e.g., missing data).
            ## If we can rule out stoch indexing for a model (using model$modelDef$varName$anyDynamicallyIndexed
            ## for the variables contained in `target`, we can probably determine
            ## which columns are fixed by inserting `.Machine$double.xmax` for parameters
            ## (need to think about assuming RHSonly don't change - could ask user?) and figuring out where those are
            ## in the `setDesignMatrix` processing (actually presumably a separate method called once).
            if(is.null(control$fixedDesignColumns)) {  # User has specified which columns are fixed.
                fixedColumns <- rep(FALSE, nCoef)
            } else {
                fixedColumns <- control$fixedDesignColumns
                if(length(fixedColumns) == 1)
                    fixedColumns <- rep(fixedColumns, nCoef)
            }
            if(all(fixedColumns)) 
                fixed <- TRUE
        } else {
            X <- control$designMatrix
            if(ncol(X) != nCoef)
                stop("polyagamma sampler: number of columns of design matrix, ", ncol(X), ", doesn't match number of parameters being sampled, ", nCoef)
            if(nrow(X) != N)
                stop("polyagamma sampler: number of rows of design matrix, ", nrow(X), ", doesn't match number of Bernoulli observations, ", N)
            fixed <- TRUE
            fixedColumns <- rep(TRUE, nCoef)
        }

        initializeSize <- TRUE
        initializeX <- TRUE
        pgSampler <- samplePolyaGamma()

        Q <- matrix(0, nrow = nCoef, ncol = nCoef)
        mu <- numeric(nCoef)				
        b <- rep(0, nCoef)
        bTemp <- rep(0, nCoef)

        if(nCoef == 1) {
            mu <- c(mu, -1)
            b <- c(b, -1)
            bTemp <- c(bTemp, -1)
            fixedColumns <- c(fixedColumns, TRUE)
        }
        if(nTarget == 1) {
            normTypes <- c(normTypes, -1)
            nodeLengths <- c(nodeLengths, -1)
        }
        
        
        probNonZero <- rep(0, N)  ## Track ids where prob == 0 (zero inflated).
        n <- N  ## Number of active (non-zero-inflated) obs.
        
        psi <- numeric(N)
        size <- numeric(N)  
        sizeContig <- numeric(N)  
        
        ## Preallocate storage for sampling. Not clear how much some or all of this helps.
        XW <- matrix(0, nrow = nCoef, ncol = N)
        Xd <- matrix(0, nrow = N, ncol = nCoef)
        kpre <- numeric(N)
        w <- numeric(N)
        one_time_fixes_done <- FALSE
    },
    run = function() {
        if(!one_time_fixes_done) one_time_fixes()
        
        if(initializeSize | stochSize)
            setSizeParam() 
    
        ## Get current values.
        y <- values(model, yNodes)
        if(singleSize) {
            k <- y - size[1]*0.5
        } else k <- y - size*0.5

        start <- 1
        for( i in 1:nTarget ) {
            end <- start + nodeLengths[i] - 1
            Q[start:end,start:end] <<- getParam_nfl[[normTypes[i]]]$getPrecision(i)
            mu[start:end] <<- getParam_nfl[[normTypes[i]]]$getMean(i)
            start <- end + 1
        }
        
        ## Build design matrix.
        if(initializeX | !fixed)
            setDesignMatrix()

        ## Determine logit(probs) and which obs are active (in zero-inflated case).
        setProbParam()

        if(zeroInflated) {
            ## `psi` already has non-zero-prob values in first n elements based on `setProbParam`.
            if(n > 0) {
                if(singleSize) {
                    w[1:n] <<- pgSampler$rpolyagamma(c(size[1]), psi[1:n])
                } else {
                    sizeContig[1:n] <<- size[probNonZero[1:n]] 
                    w[1:n] <<- pgSampler$rpolyagamma(sizeContig, psi[1:n])
                }
            }
        } else {
            if(singleSize) {
                w[1:N] <<- pgSampler$rpolyagamma(c(size[1]), psi)
            } else w[1:N] <<- pgSampler$rpolyagamma(size, psi)  ## w|beta ~ pg(n, x %*% beta)
        }
        
        ## Note that the calculations below involving X don't take advantage of sparsity, including with random intercepts or
        ## spatial processes. For the latter case, one would often have the relevant columns of X be the identity matrix.
        ## We could ask user for this information or detect it and then fill in blocks of XtWX matrix, avoiding
        ## matrix manipulations below for components of X.
        if(zeroInflated){
            if(n > 0) {
                Xd[1:n,] <<- X[probNonZero[1:n],]
                kpre[1:n] <<- k[probNonZero[1:n]]
                for( j in 1:nCoef ) {
                    XW[j,1:n] <<-  Xd[1:n,j]*w[1:n]
                    b[j] <<- sum(Xd[1:n,j] * kpre[1:n]) + sum(Q[j,] * mu)
                }
                Q1 <- XW[,1:n] %*% Xd[1:n,] + Q
            } else {  ## No relevant observations, so drawing from prior.
                for( j in 1:nCoef ) 
                    b[j] <<- sum(Q[j,]*mu)
                Q1 <- Q
            }
        } else {
            for( j in 1:nCoef ){
                XW[j,] <<-  X[,j]*w
                b[j] <<- sum(X[,j] * k) + sum(Q[j,] * mu)
            }
            Q1 <- XW %*% X + Q
        }
        UQ1 <- chol(Q1)
        M <- backsolve( UQ1, forwardsolve(t(UQ1), b) )
	
        values(model, targetAsScalar) <<- rmnorm_chol(n = 1, mean = M, cholesky = UQ1, prec_param = TRUE)
        model$calculate(calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    },
    methods = list(
        setDesignMatrix = function() {
            if(initializeX & zeroInflated) {
                ## Temporarily remove zero inflation in order to fix design matrix.
                ## This avoids repeated rounds of determining matrix values when
                ## zero inflation not initialized all at values of one.
                inflationValuesSaved <- values(model, inflationNodes)
                inflationDepsValuesSaved <- values(model, inflationNodesDeps)
                values(model, inflationNodes) <<- ones
                model$calculate(inflationNodesDeps)
                ## Check for cases like `dbern(z[i]*3*p[i])`
                if(any(values(model, probNodes) != values(model, probNodesInflated)))
                    stop("Zero inflation not specified as multiplying by one or more Bernoulli variables")
            }
            for(j in 1:nCoef) { 
                if(initializeX | !fixedColumns[j]) {
                    bTemp[j] <<- 1
                    values(model, targetAsScalar) <<- bTemp
                    model$calculate(copyNodesDeterm)
                    ## With zero inflation accounted for above, we update every element of a given column.
                    for(i in 1:N) {
                        X[i, j] <<- logit(model$getParam(yNodes[i], 'prob'))
                    }
                    bTemp[j] <<- 0
                }
            }
            if(initializeX & zeroInflated) {
                values(model, inflationNodes) <<- inflationValuesSaved
                values(model, inflationNodesDeps) <<- inflationDepsValuesSaved
            }
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
            initializeX <<- FALSE
        },
        setSizeParam = function() {
            for(i in 1:N) {
                size[i] <<- model$getParam(yNodes[i], 'size')
            }
            ## Not worth checking singleSize if it is stochastic and varies by observation. Perhaps worth considering if
            ## size is constant across observations but stochastic. Not sure if there is a use case of that.
            singleSize <<- FALSE
            if(!stochSize & initializeSize) {
                singleSize <<- TRUE
                i <- 1
                while(i <= N & singleSize) {
                    if(size[i] != size[1]) {
                        singleSize <<- FALSE
                    }
                    i <- i+1
                }
            }
            initializeSize <<- FALSE
        },
        setProbParam = function() {
            ## Note that zero size cases are handled directly in PG sampling;
            ## assumption is that these will be rare so not worth finding them
            ## here and treating as part of zero inflation.
            if(zeroInflated) {
                n <<- 0
                for(i in 1:N) {
                    probi <- model$getParam(yNodes[i], 'prob')
                    if(probi > 0 ) {
                        n <<- n + 1
                        probNonZero[n] <<- i
                        psi[n] <<- logit(probi)
                    }
                }
            } else {  ## Avoid unneeded `if` condition evaluation above for efficiency.
                for(i in 1:N) 
                    psi[i] <<- logit(model$getParam(yNodes[i], 'prob'))
            }
        },
        reset = function() { },
        one_time_fixes = function() {
            ## Run this once after compiling; remove extraneous -1 if necessary.
            if(one_time_fixes_done) return()
            mu <<- fix_one_vec(mu)
            b <<- fix_one_vec(b)
            bTemp <<- fix_one_vec(bTemp)
            one_time_fixes_done <<- TRUE
        },
        fix_one_vec = function(x = double(1)) {
            if(length(x) == 2) {
                if(x[2] == -1) {
                    ans <- numeric(length = 1, value = x[1])
                    return(ans)
                }
            }
            return(x)
            returnType(double(1))
        }
    )
)


    
#' MCMC Sampling Algorithms
#'
#' Details of the MCMC sampling algorithms provided with the NIMBLE MCMC engine; HMC samplers are in the \code{nimbleHMC} package and particle filter samplers are in the \code{nimbleSMC} package.
#'
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target node(s) on which the sampler will be used
#' @param control named list that controls the precise behavior of the sampler, with elements specific to \code{samplertype}.  The default values for control list are specified in the setup code of each sampling algorithm.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the \code{control} argument) appear below.
#'
#' @section Hamiltonian Monte Carlo samplers:
#'
#' Hamiltonian Monte Carlo (HMC) samplers are provided separately in the \code{nimbleHMC} R package.  After loading \code{nimbleHMC}, see \code{help(HMC)} for details.
#'
#' @section Particle filter samplers:
#'
#' As of Version 0.10.0 of NIMBLE, the \code{RW_PF} and \code{RW_PF_block} samplers are provided separately in the `nimbleSMC` package. After loading \code{nimbleSMC}, see \code{help(samplers)} for details.
#' 
#' @section \code{sampler_base}: base class for new samplers
#'
#' When you write a new sampler for use in a NIMBLE MCMC (see \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual}), you must include \code{contains = sampler_BASE}.
#'
#' @section binary sampler:
#'
#' The binary sampler performs Gibbs sampling for binary-valued (discrete 0/1) nodes.  This can only be used for nodes following either a \code{dbern(p)} or \code{dbinom(p, size=1)} distribution.
#'
#' The binary sampler accepts no control list arguments.
#'
#' @section categorical sampler:
#'
#' The categorical sampler performs Gibbs sampling for a single node, which generally would follow a categorical (\code{dcat}) distribution.  The categorical sampler can be assigned to other distributions as well, in which case the number of possible outcomes (1, 2, 3, ..., k) of the distribution must be specified using the 'length' control argument.
#'
#' The categorical sampler accepts the following control list elements:
#' \itemize{
#' \item length. A character string or a numeric argument.  When a character string, this should be the name of a parameter of the distribution of the target node being sampled.  The length of this distribution parameter (considered as a 1-dimensional vector) will be used to determine the number of possible outcomes of the target node's distribution.  When a numeric value, this value will be used as the number of possible outcomes of the target node's distribution.  (default = "prob")
#' \item check. A logical argument.  When FALSE, no check for a 'dcat' prior distribution for the target node takes place. (default = TRUE)
#' }
#' 
#' @section RW sampler:
#'
#' The RW sampler executes adaptive Metropolis-Hastings sampling with a normal proposal distribution (Metropolis, 1953), implementing the adaptation routine given in Shaby and Wells, 2011.  This sampler can be applied to any scalar continuous-valued stochastic node, and can optionally sample on a log scale.
#'
#' The RW sampler accepts the following control list elements:
#' \itemize{
#' \item log. A logical argument, specifying whether the sampler should operate on the log scale. (default = FALSE)
#' \item reflective. A logical argument, specifying whether the normal proposal distribution should reflect to stay within the range of the target distribution. (default = FALSE)
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (proposal standard deviation) throughout the course of MCMC execution to achieve a theoretically desirable acceptance rate. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW sampler will perform its adaptation procedure.  This updates the scale variable, based upon the sampler's achieved acceptance rate over the past adaptInterval iterations. (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the normal proposal standard deviation.  If adaptive = FALSE, scale will never change. (default = 1)
#' }
#'
#' The RW sampler cannot be used with options log=TRUE and reflective=TRUE, i.e. it cannot do reflective sampling on a log scale.
#'
#' After an MCMC algorithm has been configuration and built, the value of the proposal standard deviation of a RW sampler can be modified using the setScale method of the sampler object.  This use the scalar argument to modify the current value of the proposal standard deviation, as well as modifying the initial (pre-adaptation) value to which the proposal standard deviation is reset, at the onset of a new MCMC chain.
#'
#' @section RW_block sampler:
#'
#' The RW_block sampler performs a simultaneous update of one or more model nodes, using an adaptive Metropolis-Hastings algorithm with a multivariate normal proposal distribution (Roberts and Sahu, 1997), implementing the adaptation routine given in Shaby and Wells, 2011.  This sampler may be applied to any set of continuous-valued model nodes, to any single continuous-valued multivariate model node, or to any combination thereof. \cr
#'
#' The RW_block sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (a coefficient for the entire proposal covariance matrix) and propCov (the multivariate normal proposal covariance matrix) throughout the course of MCMC execution.  If only the scale should undergo adaptation, this argument should be specified as TRUE. (default = TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for scale (TRUE) or also for provCov (FALSE).  This argument is only relevant when adaptive = TRUE.  When adaptScaleOnly = FALSE, both scale and propCov undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When adaptScaleOnly = TRUE, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW_block sampler will perform its adaptation procedure, based on the past adaptInterval iterations. (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the scalar multiplier for propCov.  If adaptive = FALSE, scale will never change. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string 'identity', in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default = 'identity')
#' \item tries. The number of times this sampler will repeatedly operate on each MCMC iteration.  Each try consists of a new proposed transition and an accept/reject decision of this proposal.  Specifying tries > 1 can help increase the overall sampler acceptance rate and therefore chain mixing. (default = 1)
#' }
#'
#' After an MCMC algorithm has been configuration and built, the value of the proposal standard deviation of a RW_block sampler can be modified using the setScale method of the sampler object.  This use the scalar argument to will modify the current value of the proposal standard deviation, as well as modifying the initial (pre-adaptation) value which the proposal standard deviation is reset to, at the onset of a new MCMC chain.
#'
#' Operating analogous to the setScale method, the RW_block sampler also has a setPropCov method.  This method accepts a single matrix-valued argument, which will modify both the current and initial (used at the onset of a new MCMC chain) values of the multivariate normal proposal covariance.
#'
#' Note that modifying elements of the control list may greatly affect the performance of this sampler. In particular, the sampler can take a long time to find a good proposal covariance when the elements being sampled are not on the same scale. We recommend providing an informed value for \code{propCov} in this case (possibly simply a diagonal matrix that approximates the relative scales), as well as possibly providing a value of \code{scale} that errs on the side of being too small. You may also consider decreasing \code{adaptFactorExponent} and/or \code{adaptInterval}, as doing so has greatly improved performance in some cases. 
#'
#' @section RW_llFunction sampler:
#'
#' Sometimes it is useful to control the log likelihood calculations used for an MCMC updater instead of simply using the model.  For example, one could use a sampler with a log likelihood that analytically (or numerically) integrates over latent model nodes.  Or one could use a sampler with a log likelihood that comes from a stochastic approximation such as a particle filter, allowing composition of a particle MCMC (PMCMC) algorithm (Andrieu et al., 2010).  The RW_llFunction sampler handles this by using a Metropolis-Hastings algorithm with a normal proposal distribution and a user-provided log-likelihood function.  To allow compiled execution, the log-likelihood function must be provided as a specialized instance of a nimbleFunction.  The log-likelihood function may use the same model as the MCMC as a setup argument, but if so the state of the model should be unchanged during execution of the function (or you must understand the implications otherwise).
#'
#' The RW_llFunction sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (proposal standard deviation) throughout the course of MCMC execution. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the normal proposal standard deviation. (default = 1)
#' \item llFunction. A specialized nimbleFunction that accepts no arguments and returns a scalar double number.  The return value must be the total log-likelihood of all stochastic dependents of the target nodes -- and, if includesTarget = TRUE, of the target node(s) themselves --  or whatever surrogate is being used for the total log-likelihood.  This is a required element with no default.
#' \item includesTarget. Logical variable indicating whether the return value of llFunction includes the log-likelihood associated with target.  This is a required element with no default.
#' }
#'
#' @section slice sampler:
#'
#' The slice sampler performs slice sampling of the scalar node to which it is applied (Neal, 2003).  This sampler can operate on either continuous-valued or discrete-valued scalar nodes.  The slice sampler performs a 'stepping out' procedure, in which the slice is iteratively expanded to the left or right by an amount sliceWidth.  This sampler is optionally adaptive, governed by a control list element, whereby the value of sliceWidth is adapted towards the observed absolute difference between successive samples.
#'
#' The slice sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler will adapt the value of sliceWidth throughout the course of MCMC execution. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item sliceWidth. The initial value of the width of each slice, and also the width of the expansion during the iterative 'stepping out' procedure. (default = 1)
#' \item sliceMaxSteps. The maximum number of expansions which may occur during the 'stepping out' procedure. (default = 100)
#' \item maxContractions. The maximum number of contractions of the interval that may occur during sampling (this prevents infinite looping in unusual situations). (default = 100)
#' \item maxContractionsWarning. A logical argument specifying whether to warn when the maximum number of contractions is reached. (default = TRUE)
#' }
#'
#' @section ess sampler:
#'
#' The ess sampler performs elliptical slice sampling of a single node, which must follow either a univariate or multivariate normal distribution (Murray, 2010).  The algorithm is an extension of slice sampling (Neal, 2003), generalized to context of the Gaussian distribution.  An auxiliary variable is used to identify points on an ellipse (which passes through the current node value) as candidate samples, which are accepted contingent upon a likelihood evaluation at that point.  This algorithm requires no tuning parameters and therefore no period of adaptation, and may result in very efficient sampling from Gaussian distributions.
#'
#' The ess sampler accepts the following control list arguments.
#' \itemize{
#' \item maxContractions. The maximum number of contractions of the interval that may occur during sampling (this prevents infinite looping in unusual situations). (default = 100)
#' \item maxContractionsWarning. A logical argument specifying whether to warn when the maximum number of contractions is reached. (default = TRUE)
#' }
#' 
#' @section AF_slice sampler:
#' 
#' The automated factor slice sampler conducts a slice sampling algorithm on one or more model nodes.  The sampler uses the eigenvectors of the posterior covariance between these nodes as an orthogonal basis on which to perform its 'stepping Out' procedure.  The sampler is adaptive in updating both the width of the slices and the values of the eigenvectors.  The sampler can be applied to any set of continuous or discrete-valued model nodes, to any single continuous or discrete-valued multivariate model node, or to any combination thereof. 
#
#' The automated factor slice sampler accepts the following control list elements:
#' \itemize{
#' \item sliceWidths.  A numeric vector of initial slice widths.  The length of the vector must be equal to the sum of the lengths of all nodes being used by the automated factor slice sampler.  Defaults to a vector of 1's.
#' \item sliceAdaptFactorMaxIter.  The number of iterations for which the factors (eigenvectors) will continue to adapt to the posterior correlation. (default = 15000)
#' \item sliceAdaptFactorInterval.  The interval on which to perform factor adaptation. (default = 200)
#' \item sliceAdaptWidthMaxIter.  The maximum number of iterations for which to adapt the widths for a given set of factors. (default = 512)
#' \item sliceAdaptWidthTolerance. The tolerance for when widths no longer need to adapt, between 0 and 0.5. (default = 0.1)
#' \item sliceMaxSteps.  The maximum number of expansions which may occur during the 'stepping out' procedure. (default = 100)
#' \item maxContractions. The maximum number of contractions of the interval that may occur during sampling (this prevents infinite looping in unusual situations). (default = 100)
#' \item maxContractionsWarning. A logical argument specifying whether to warn when the maximum number of contractions is reached. (default = TRUE)
#' }
#'
#' @section crossLevel sampler:
#'
#' This sampler is constructed to perform simultaneous updates across two levels of stochastic dependence in the model structure.  This is possible when all stochastic descendents of node(s) at one level have conjugate relationships with their own stochastic descendents.  In this situation, a Metropolis-Hastings algorithm may be used, in which a multivariate normal proposal distribution is used for the higher-level nodes, and the corresponding proposals for the lower-level nodes undergo Gibbs (conjugate) sampling.  The joint proposal is either accepted or rejected for all nodes involved based upon the Metropolis-Hastings ratio. This sampler is a conjugate version of Scheme 3 in Knorr-Held and Rue (2002). It can also be seen as a Metropolis-based version of collapsed Gibbs sampling (in particular Sampler 3 of van Dyk and Park (2008)). 
#'
#' The requirement that all stochastic descendents of the target nodes must themselves have only conjugate descendents will be checked when the MCMC algorithm is built.  This sampler is useful when there is strong dependence across the levels of a model that causes problems with convergence or mixing.
#'
#' The crossLevel sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. Logical argument, specifying whether the multivariate normal proposal distribution for the target nodes should be adaptived. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the scalar multiplier for propCov. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string 'identity' or any positive definite matrix of the appropriate dimensions. (default = 'identity')
#' }
#'
#' @section RW_llFunction_block sampler:
#'
#' Sometimes it is useful to control the log likelihood calculations used for an MCMC updater instead of simply using the model.  For example, one could use a sampler with a log likelihood that analytically (or numerically) integrates over latent model nodes.  Or one could use a sampler with a log likelihood that comes from a stochastic approximation such as a particle filter, allowing composition of a particle MCMC (PMCMC) algorithm (Andrieu et al., 2010) (but see samplers listed below for NIMBLE's direct implementation of PMCMC).  The \code{RW_llFunction_block} sampler handles this by using a Metropolis-Hastings algorithm with a multivariate normal proposal distribution and a user-provided log-likelihood function.  To allow compiled execution, the log-likelihood function must be provided as a specialized instance of a nimbleFunction.  The log-likelihood function may use the same model as the MCMC as a setup argument, but if so the state of the model should be unchanged during execution of the function (or you must understand the implications otherwise).
#'
#' The RW_llFunction_block sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the proposal covariance throughout the course of MCMC execution. (default is TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for scale (TRUE) or also for provCov (FALSE).  This argument is only relevant when adaptive = TRUE.  When adaptScaleOnly = FALSE, both scale and propCov undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When adaptScaleOnly = TRUE, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the scalar multiplier for propCov.  If adaptive = FALSE, scale will never change. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string 'identity', in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default = 'identity')
#' \item llFunction. A specialized nimbleFunction that accepts no arguments and returns a scalar double number.  The return value must be the total log-likelihood of all stochastic dependents of the target nodes -- and, if includesTarget = TRUE, of the target node(s) themselves --  or whatever surrogate is being used for the total log-likelihood.  This is a required element with no default.
#' \item includesTarget. Logical variable indicating whether the return value of llFunction includes the log-likelihood associated with target.  This is a required element with no default.
#' }
#'
#' @section RW_dirichlet sampler:
#'
#' This sampler is designed for sampling non-conjugate Dirichlet distributions.  The sampler performs a series of Metropolis-Hastings updates (on the log scale) to each component of a gamma-reparameterization of the target Dirichlet distribution.  The acceptance or rejection of these proposals follows a standard Metropolis-Hastings procedure.
#'
#' The \code{RW_dirichlet} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should independently adapt the scale (proposal standard deviation, on the log scale) for each componentwise Metropolis-Hasting update, to achieve a theoretically desirable acceptance rate for each. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the sampler will perform its adaptation procedure.  (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the proposal standard deviation (on the log scale) for each component of the reparameterized Dirichlet distribution.  If adaptive = FALSE, the proposal standard deviations will never change. (default = 1)
#' }
#'
#' @section RW_wishart sampler:
#'
#' This sampler is designed for sampling non-conjugate Wishart and inverse-Wishart distributions.  More generally, it can update any symmetric positive-definite matrix (for example, scaled covariance or precision matrices).  The sampler performs block Metropolis-Hastings updates following a transformation to an unconstrained scale (Cholesky factorization of the original matrix, then taking the log of the main diagonal elements.
#'
#' The \code{RW_wishart} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale and proposal covariance for the multivariate normal Metropolis-Hasting proposals, to achieve a theoretically desirable acceptance rate for each. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the sampler will perform its adaptation procedure.  (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the scalar multiplier for the multivariate normal Metropolis-Hastings proposal covariance.  If adaptive = FALSE, scale will never change. (default = 1)
#' }
#'
#' @section RW_block_lkj_corr_cholesky sampler:
#'
#' This sampler is designed for sampling non-conjugate LKJ correlation Cholesky factor distributions. The sampler performs a blocked Metropolis-Hastings update following a transformation to an unconstrained scale (using the signed stickbreaking approach documented in Section 10.12 of the Stan Language Reference Manual, version 2.27). 
#'
#' The \code{RW_block_lkj_corr_cholesky} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (a coefficient for the entire proposal covariance matrix) and propCov (the multivariate normal proposal covariance matrix) throughout the course of MCMC execution.  If only the scale should undergo adaptation, this argument should be specified as TRUE. (default = TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for scale (TRUE) or also for provCov (FALSE).  This argument is only relevant when adaptive = TRUE.  When adaptScaleOnly = FALSE, both scale and propCov undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When adaptScaleOnly = TRUE, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW_block sampler will perform its adaptation procedure, based on the past adaptInterval iterations. (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the scalar multiplier for propCov.  If adaptive = FALSE, scale will never change. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string 'identity', in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default = 'identity')
#' }
#'
#' This is the default sampler for the LKJ distribution. However, blocked samplers may perform poorly if the adaptation configuration is poorly chosen. See the comments in the RW_block section of this documentation.
#'
#' @section RW_lkj_corr_cholesky sampler:
#'
#' This sampler is designed for sampling non-conjugate LKJ correlation Cholesky factor distributions. The sampler performs individual Metropolis-Hastings updates following a transformation to an unconstrained scale (using the signed stickbreaking approach documented in Section 10.12 of the Stan Language Reference Manual, version 2.27). 
#'
#' The \code{RW_lkj_corr_cholesky} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scales of the univariate normal Metropolis-Hasting proposals, to achieve a theoretically desirable acceptance rate for each. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the sampler will perform its adaptation procedure.  (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the scalar multiplier for the multivariate normal Metropolis-Hastings proposal covariance.  If adaptive = FALSE, scale will never change. (default = 1)
#' }
#'
#' Note that this sampler is likely run much more slowly than the blocked sampler for the LKJ distribution, as updating each single element will generally incur the full cost of updating all dependencies of the entire matrix. 
#'
#' @section polyagamma sampler:
#'
#' The polyagamma sampler uses Pólya-gamma data augmentation to do conjugate sampling for the parameters in the linear predictor of a logistic regression model (Polson et al., 2013), analogous to the Albert-Chib data augmentation scheme for probit regression. This sampler is not assigned as a default sampler by \code{configureMCMC} and so can only be used if manually added to an MCMC configuration.
#'
#' As an example, consider model code containing:
#' \preformatted{for(i in 1:n) {
#'   logit(prob[i]) <- beta0 + beta1*x1[i] + beta2*x2[i] + u[group[i]] 
#'   y[i] ~ dbinom(prob = prob[i], size = size[i])
#' }
#' for(j in 1:num_groups)
#'   u[j] ~ dnorm(0, sd = sigma_group)
#' }
#' where \code{beta0}, \code{beta1}, and \code{beta2} are fixed effects with normal priors, \code{u[j]} are random effects associated with groups of data, and \code{group[i]} is from the \code{constants} list and gives the group index of the i-th observation, \code{y[i]}. In this model, the parameters \code{beta0}, \code{beta1}, \code{beta2}, and all the \code{u[j]} effects can be jointly sampled by the polyagamma sampler, i.e., they will be its target nodes.
#'
#' After building a model (calling \code{nimbleModel}) containing the above code, one would configure the sampler as follows:
#' \preformatted{MCMCconf <- configureMCMC(model)
#' logistic_nodes <- c("beta0", "beta1", "beta2", "u")
#' MCMCconf$removeSamplers(logistic_nodes) # Optionally, remove default samplers.
#' MCMCconf$addSampler(target = logistic_nodes, type = "polyagamma", control = list(fixedDesignColumns=TRUE))
#' }
#' 
#' As shown here, the stochastic dependencies (\code{y[i]} here) of the target nodes must follow \code{dbin} or \code{dbern} distributions. The logit transformation of their probability parameter must be a linear function (technically an affine function) of the target nodes, which themselves must have `dnorm` or `dmnorm` priors. Zero inflation to account for structural zeroes is also supported, allowed as discussed below. The stochastic dependencies will often but not always be the observations in the logistic regression and will be referred to as 'responses' henceforth. Internally, the sampler draws latent values from the Pólya-gamma distribution, one per response. These latent values are then used to draw from the multivariate normal conditional distribution of the target nodes.
#'
#' Importantly, note that because the Pólya-gamma draws are not retained when an iteration of the sampler finishes, one generally wants to apply the sampler to all parameter nodes involved in the linear predictor of the logistic regression, to avoid duplicative Pólya-gamma draws of the latent values. If there are stochastic indices (e.g., if \code{group[i]} above is stochastic), the Pólya-gamma sampler can still be used, but the stochastic nodes cannot be sampled by it and must have separate sampler(s). It is also possible in some models that regression parameters can be split into (conditionally independent) groups that can be sampled independently, e.g., if one has distinct logistic regression specifications for different sets of responses in the model.
#' 
#' Sampling involves use of the design matrix. The design matrix includes one column corresponding to each regression covariate as well as one columnar block corresponding to each random effect. In the example above, the columns would include a vector of ones (to multiply \code{beta0}), the vectors \code{x1} and \code{x2}, and a vector of indicators for each \code{u[j]} (with a 1 in row \code{i} if \code{group[i]} is \code{j}), resulting in \code{3+num_groups} columns. Note that the polyagamma sampler can determine the design matrix from the model, even when written as above such that the design matrix is not explicitly in the model code. It is also possible to write model code for the linear prediction using matrix multiplication and an explicit design matrix, but that is not necessary.
#'
#' Often the design matrix is fixed in advance, but in some cases elements of the matrix may be stochastic and sampled during the MCMC. That would be the case if there are missing values (e.g., missing covariate values (declared as stochastic nodes)) or stochastic indexing (e.g., unknown assignment of responses to clusters, as mentioned above). Note that changes in the values of any of the \code{beta} or \code{u[j]} target nodes in the example above do not change the design matrix.
#' 
#' Recalculating the elements of the design matrix at every iteration is costly and will likely greatly slow the sampler. To alleviate this, users can specify which columns of the design matrix are fixed (non-stochastic) using the `fixedDesignColumns` control argument. If all columns are fixed, which will often be the case, this argument can be specified simply as `TRUE`. Note that the sampler does not determine if any or all columns are fixed, so users wishing to take advantage of the large speed gains from having fixed columns should provide this control argument. Columns indicated as fixed will be determined (if necessary) when the sampler is first run and retained for subsequent iterations.
#'
#' By default, NIMBLE will determine the design matrix (and as discussed above will do so repeatedly at each iteration for any columns not indicated as being fixed). If the matrix has no stochastic elements, users may choose to provide the matrix directly to the sampler via the `designMatrix` control argument. This will save time in computing the matrix initially but likely will have limited benefit relative to the cost of running many iterations of MCMC and therefore can be omitted in most cases.
#'
#' The sampler allows for binomial responses with stochastic sizes. This would be the case in the above example if the \code{size[i]} values are themselves declared as unobserved stochastic nodes and thus are sampled by MCMC.
#'
#' The sampler allows for zero inflation of the response probability in that the probability determined by the inverse logit transformation of the linear predictor can be multipled by one or more binary scalar nodes to produce the response probability. These binary nodes must be specified via the `dbern` distribution or the `dbin` distribution with size equal to one. This functionality is intended for use in cases where another part of the model introduces structural zeroes, such as in determining occupancy in ecological occupancy models. An example would be if the above were modified by \code{y[i] ~ dbinom(prob = z[i] * p[i], size = size[i])}, where each \code{z[i]} is either 0 or 1.
#'
#' The polyagamma sampler accepts the following control list elements:
#' \itemize{
#' \item fixedDesignColumns. Either a single logical value indicating if the design matrix is fixed (non-stochastic) or a logical vector indicating which columns are fixed. In the latter case, the columns must be ordered exactly as the ordering of target node elements given by \code{model$expandNodeNames(target, returnScalarComponents = TRUE)}, where \code{target} is the same as the \code{target} argument to \code{configureMCMC$addSampler} above. (default = FALSE)
#' \item designMatrix. The full design matrix with rows corresponding to the ordering of the responses and columns ordered exactly as the ordering of target node elements given by \code{model$expandNodeNames(target, returnScalarComponents = TRUE)}, where \code{target} is the same as the \code{target} argument to \code{configureMCMC$addSampler} above. If provided, all columns are assumed to be fixed, ignoring the \code{fixedDesignColumns} control element.
#' \item nonTargetNodes. Additional stochastic nodes involved in the linear predictor that are not to be sampled as part of the sampler. This must include any nodes specifying stochastic indexes (e.g., \code{"group"} if the \code{group[i]} values are stochastic) and any parameters considered known or that for any reason one does not want to sample. Providing \code{nonTargetNodes} is required in order to allow NIMBLE to check for the presence of zero inflation.  
#' \item check. A logical value indicating whether NIMBLE should check various conditions required for validity of the sampler. This is provided for rare cases where the checking may be overly conservative and a user is sure that the sampler is valid and wants to override the checking. (default = TRUE)
#' }
#'
#' @section noncentered sampler:
#'
#' The noncentered sampler is designed to sample the mean or standard deviation of a set of centered random effects while also moving the random effects values to possibly allow better mixing. The noncentered sampler deterministically shifts or scales the dependent node values to be consistent with the proposed value of the target (the mean or the standard deviation) such that the effect is to sample in a noncentered parameterization (Yu and Meng 2011), via an on-the-fly reparameterization. This can improve mixing by updating the target node based on information in the model nodes whose parent nodes are the dependent nodes of the target (i.e,. the "grandchild" nodes of the target; these will often be data nodes). This comes at the extra computational cost of calculating the logProbability of the "grandchild" nodes.
#'
#' It is still necessary to have other samplers on the random effects values.
#' 
#' Mathematically, the noncentered sampler operates in one dimension of a transformed parameter space. When sampling a mean, all random effects will be shifted by the same amount as the mean. When sampling a standard deviation, all random effects (relative to their means) will be scaled by the same factor as the standard deviation.

#' Consider a model that includes the following code:
#' \preformatted{ for(i in 1:n)
#'   y[i] ~ dnorm(beta1*x[i] + u[group[i]], sd = sigma_obs)
#' for(j in 1:num_groups)
#'   u[j] ~ dnorm(beta0, sd = sigma_group)
#' }
#' where \code{u[j]} is a random effect associated with groups of data, and \code{group[i]} gives the group index of the i-th observation. This model has centered random effects, the \code{u[j]}, because these have the intercept \code{beta0} as their mean. In basic univariate sampling, updates to \code{beta0} or to \code{sigma_group} do not change \code{u[j]}, making only small moves possible if \code{num_groups} is large. When the noncentered sampler considers a new value for \code{beta0}, it will shift all the \code{u[j]} so that (in this case) their prior probabilities do not change. If the noncentered sampler considers a new value for \code{sigma_group}, it will rescale all the \code{u[j]} accordingly.
#'
#' The effect of such a sampling strategy is to update \code{beta0} and \code{sigma_group} as if the model had been written in a different (noncentered) way. For updating \code{beta0}, it would be:
#' \preformatted{ for(i in 1:n)
#'   y[i] ~ dnorm(beta0 + beta1*x[i] + u[group[i]], sd = sigma_obs)
#' for(j in 1:num_groups)
#'   u[j] ~ dnorm(0, sd = sigma_group)
#' }
#' For updating \code{sigma_group}, it would be:
#' \preformatted{ for(i in 1:n)
#'   y[i] ~ dnorm(beta1*x[i] + u[group[i]] * sigma_group, sd = sigma_obs)
#' for(j in 1:num_groups)
#'   u[j] ~ dnorm(beta0, sd = 1)
#' }
#'
#' Whether centered or noncentered parameterizations result in better sampling can depend on the model and the data. Therefore Yu and Meng (2011) recommended an "interweaving" strategy of using both kinds of samplers. Adding the noncentered sampler (on either the mean or standard deviation or both) to an existing MCMC configuration for a model specified using the centered parameterization (and with an sampler already assigned to the target node) produces an overall sampling approach that is a variation on the interweaving strategy of Yu and Meng (2011). This provides the benefits of sampling in both the centered and noncentered parameterizations in a single MCMC. 
#'
#' There is a higher computational cost to the noncentered sampler (or to writing the model directly in one of the equivalent ways shown). The cost is that when updating \code{beta0} or \code{sigma_group}, the relevant log probabilities calculations will include (in this case) all the of \code{y[i]}, i.e. the "grandchild" nodes of \code{beta0} or \code{sigma_group}.
#'
#' The noncentered sampler is not assigned by default by \code{configureMCMC} but must be manually added. For example:
#' \preformatted{ MCMCconf <- configureMCMC(model)
#' MCMCconf$addSampler(target = "beta0", type = "noncentered",
#'   control = list(param = "location", sampler = "RW"))
#' MCMCconf$addSampler(target = "sigma_group", type = "noncentered",
#'   control = list(param = "scale", sampler = "RW"))
#' }
#'
#' While the target node will generally be either the mean (location) or standard deviation (scale) of a set of other nodes (e.g., random effects), it could in theory be used in other contexts and one can choose whether the transformation is a shift or a scale operation. In a shift operation (e.g., when the sampling target is a mean), the dependent nodes are set to their previous values plus the difference between the proposed value and previous value for the target. In a scale operation (e.g., when the sampling target is the standard deviation), the dependent nodes minus their means are multiplied by the ratio of the proposed value to the previous value for the target and the previous value for the target. Whether to shift or scale is determined from the \code{param} element of the control list.
#'
#' The sampling algorithm for the target node can either be adaptive Metropolis random walk (which uses NIMBLE's \code{RW} sampler) or slice sampling (which uses NIMBLE's \code{slice} sampler), determined from the \code{sampler} element of the control list. In either case, the underlying sampling accounts for the Jacobian of the deterministic shifting or scaling of the dependent nodes (in the case of shifting, the Jacobian is equal to 1 and has no impact). When the target is the standard deviation of normally-distributed dependent nodes, the Jacobian cancels with the prior distribution for the dependent nodes, and the update is in effect based only on the prior for the target and the distribution of the "grandchild" nodes. 
#'
#' The \code{noncentered} sampler accepts the following control list elements:
#' \itemize{
#' \item sampler. A character string, either \code{"RW"} or \code{"slice"} specifying the type of sampler to be used for the target node. (default = \code{"RW"})
#' \item param. A character string, either \code{"location"} or \code{"scale"} specifying whether sampling is done as shifting or scaling the dependent nodes. (default = \code{"location"})
#' }
#' 
#' @section CAR_normal sampler:
#'
#' The CAR_normal sampler operates uniquely on improper (intrinsic) Gaussian conditional autoregressive (CAR) nodes, those with a \code{dcar_normal} prior distribution.  It internally assigns one of three univariate samplers to each dimension of the target node: a posterior predictive, conjugate, or RW sampler; however these component samplers are specialized to operate on dimensions of a \code{dcar_normal} distribution.
#'
#' The CAR_normal sampler accepts the following control list elements:
#' \itemize{
#' \item \code{carUseConjugacy}. A logical argument, specifying whether to assign conjugate samplers for conjugate components of the target node. If \code{FALSE}, a RW sampler would be assigned instead. (default = TRUE)
#' \item \code{adaptive}. A logical argument, specifying whether any component RW samplers should adapt the scale (proposal standard deviation), to achieve a theoretically desirable acceptance rate. (default = \code{TRUE})
#' \item \code{adaptInterval}. The interval on which to perform adaptation for any component RW samplers.  Every \code{adaptInterval} MCMC iterations (prior to thinning), component RW samplers will perform an adaptation procedure.  This updates the \code{scale} variable, based upon the sampler's achieved acceptance rate over the past \code{adaptInterval} iterations. (default = 200)
#' \item \code{scale}. The initial value of the normal proposal standard deviation for any component RW samplers.  If \code{adaptive = FALSE}, \code{scale} will never change. (default = 1)
#' }
#'
#' @section CAR_proper sampler:
#'
#' The CAR_proper sampler operates uniquely on proper Gaussian conditional autoregressive (CAR) nodes, those with a \code{dcar_proper} prior distribution.  It internally assigns one of three univariate samplers to each dimension of the target node: a posterior predictive, conjugate, or RW sampler, however these component samplers are specialized to operate on dimensions of a \code{dcar_proper} distribution.
#'
#' The CAR_proper sampler accepts the following control list elements:
#' \itemize{
#' \item \code{carUseConjugacy}. A logical argument, specifying whether to assign conjugate samplers for conjugate components of the target node. If \code{FALSE}, a RW sampler would be assigned instead. (default = \code{TRUE})
#' \item \code{adaptive}. A logical argument, specifying whether any component RW samplers should adapt the scale (proposal standard deviation), to achieve a theoretically desirable acceptance rate. (default = \code{TRUE})
#' \item \code{adaptInterval}. The interval on which to perform adaptation for any component RW samplers.  Every adaptInterval MCMC iterations (prior to thinning), component RW samplers will perform an adaptation procedure.  This updates the scale variable, based upon the sampler's achieved acceptance rate over the past adaptInterval iterations. (default = 200)
#' \item \code{scale}. The initial value of the normal proposal standard deviation for any component RW samplers.  If \code{adaptive = FALSE}, \code{scale} will never change. (default = 1)
#' }
#'
#' @section CRP sampler:
#' 
#' The CRP sampler is designed for fitting models involving Dirichlet process mixtures. It is exclusively assigned by NIMBLE's default MCMC configuration to nodes having the Chinese Restaurant Process distribution, \code{dCRP}. It executes sequential sampling of each component of the node (i.e., the cluster membership of each element being clustered). Internally, either of two samplers can be assigned, depending on conjugate or non-conjugate structures within the model. For conjugate and non-conjugate model structures, updates are based on Algorithm 2 and Algorithm 8 in Neal (2000), respectively.
#'
#  The CRP sampler accepts the following control list elements:
#' \itemize{
#' \item \code{checkConjugacy}. A logical argument, specifying whether to assign conjugate samplers if valid. (default = \code{TRUE})
#' \item \code{printTruncation}. A logical argument, specifying whether to print a warning when the MCMC attempts to use more clusters than the maximum number specified in the model. Only relevant where the user has specified the maximum number of clusters to be less than the number of observations. (default = \code{TRUE})
#' }
#' 
#' @section CRP_concentration sampler:
#' 
#' The CRP_concentration sampler is designed for Bayesian nonparametric mixture modeling. It is exclusively assigned to the concentration parameter of the Dirichlet process when the model is specified using the Chinese Restaurant Process distribution, \code{dCRP}. This sampler is assigned by default by NIMBLE's default MCMC configuration and can only be used when the prior for the concentration parameter is a gamma distribution. The assigned sampler is an augmented beta-gamma sampler as discussed in Section 6 in Escobar and West (1995).
#'
#' @section prior_samples sampler:
#'
#' The prior_samples sampler uses a provided set of numeric values (\code{samples}) to define the prior distribution of one or more model nodes.  One every MCMC iteration, the prior_samples sampler takes value(s) from the numeric values provided, and stores these value(s) into the target model node(s).  This allows one to define the prior distribution of model parameters empirically, using a set of numeric \code{samples}, presumably obtained previously using MCMC.  The \code{target} node may be either a single scalar node (scalar case), or a collection of model nodes.
#'
#' The prior_samples sampler provides two options for selection of the value to use on each MCMC iteration.  The default behaviour is to take sequential values from the \code{samples} vector (scalar case), or in the case of multiple dimensions, sequential rows of the \code{samples} matrix are used.  The alternative behaviour, by setting the control argument \code{randomDraws = TRUE}, will instead use random draws from the \code{samples} vector (scalar case), or randomly selected rows of the \code{samples} matrix in the multidimensional case.
#'
#' If the default of sequential selection of values is used, and the number of MCMC iterations exceeds the length of the \code{samples} vector (scalar case) or the number of rows of the \code{samples} matrix, then \code{samples} will be recycled as necessary for the number of MCMC iterations.  A message to this effect is also printed at the beginning of the MCMC chain.
#'
#' Logically, prior_samples samplers might want to operate first, in advance of other samplers, on every MCMC iteration.  By default, at the time of MCMC building, all prior_samples samplers are re-ordered to appear first in the list of MCMC samplers.  This behaviour can be subverted, however, by setting nimbleOptions(MCMCorderPriorSamplesSamplersFirst = FALSE).
#'
#' The prior_samples sampler can be assigned to non-stochastic model nodes (nodes which are not assigned a prior distribution in the model). In fact, it is recommended that nodes being assigned a prior_samples are not provided with a prior distribution in the model, and rather, that these nodes only appear on the right-hand-side of model declaration lines.  In such case that a prior_samples sampler is assigned to a nodes with a prior distribution, the prior distribution will be overridden by the sample values provided to the sampler; however, the node will still be a stochastic node for other purposes, and will contribute to the model joint-density (using the sample values provided relative to the prior distribution), will have an MCMC sampler assigned to it by default, and also may introduce potential for confusion.  In this case, a message is issued at the time of MCMC building.
#'
#' The prior_samples sampler accepts the following control list elements:
#' \itemize{
#' \item \code{samples}. A numeric vector or matrix.  When the \code{target} node is a single scalar-valued node, \code{samples} should be a numeric vector.  When the \code{target} node specifies d > 2 model dimensions, \code{samples} should be a matrix containing d columns.  The \code{samples} control argument is required.
#' \item \code{randomDraws}. A logical argument, specifying whether to use a random draw from \code{samples} on each iteration.  If \code{samples} is a matrix, then a randomly-selected row of the \code{samples} matrix is used.  When \code{FALSE}, sequential values (or sequential matrix rows) are used (default = \code{FALSE}).
#' }
#'
#' @section posterior_predictive sampler:
#'
#' The posterior_predictive sampler operates only on posterior predictive stochastic nodes. A posterior predictive node is a node that is not itself data and has no data nodes in its entire downstream (descendant) dependency network. Note that such nodes play no role in inference for model parameters but have often been included in BUGS models to make predictions, including for posterior predictive checks. As of version 0.13.0, NIMBLE samples model parameters without conditioning on the posterior predictive nodes and samples conditionally from the posterior predictive nodes as the last step of each MCMC iteration.
#'
#' (Also note that NIMBLE allows posterior predictive values to be simulated independently of running MCMC, for example by writing a nimbleFunction to do so.  This means that in many cases where terminal stochastic (posterior predictive) nodes have been included in BUGS models, they are not needed when using NIMBLE.)
#'
#' The posterior_predictive sampler functions by simulating new values for all downstream (dependent) nodes using their conditional distributions, as well as updating the associated model probabilities.  A posterior_predictive sampler will automatically be assigned to all trailing non-data stochastic nodes in a model, or when possible, to any node at a point in the model after which all downstream (dependent) stochastic nodes are non-data.
#'
#' The posterior_predictive sampler accepts no control list arguments.
#'
#' @section RJ_fixed_prior sampler:
#'
#' This sampler proposes addition/removal for variable of interest in the framework of variable selection using reversible jump MCMC, with a specified prior probability of inclusion.  A normal proposal distribution is used to generate proposals for the addition of the variable. This is a specialized sampler used by \code{configureRJ} function, when the model code is written without using indicator variables. See \code{help{configureRJ}} for details. It is not intended for direct assignment.
#'
#' @section RJ_indicator sampler:
#'
#' This sampler proposes transitions of a binary indicator variable, corresponding to a variable of interest, in the framework of variable selection using reversible jump MCMC.  This is a specialized sampler used by \code{configureRJ} function, when the model code is written using indicator variables. See \code{help{configureRJ}} for details. It is not intended for direct assignment.
#'
#' @section RJ_toggled sampler:
#'
#' This sampler operates in the framework of variable selection using reversible jump MCMC.  Specifically, it conditionally performs updates of the target variable of interest using the originally-specified sampling configuration, when variable is "in the model".  This is a specialized sampler used by \code{configureRJ} when adding a reversible jump MCMC . See \code{help{configureRJ}} for details. It is not intended for direct assignment.
#' 
#' @name samplers
#'
#' @aliases sampler binary categorical prior_samples posterior_predictive RW RW_block RW_dirichlet RW_wishart RW_llFunction slice AF_slice crossLevel RW_llFunction_block sampler_prior_samples sampler_posterior_predictive sampler_binary sampler_categorical sampler_RW sampler_RW_block sampler_RW_multinomial sampler_RW_dirichlet sampler_RW_wishart sampler_RW_llFunction sampler_slice sampler_AF_slice sampler_crossLevel sampler_RW_llFunction_block CRP CRP_concentration DPmeasure RJ_fixed_prior RJ_indicator RJ_toggled RW_PF RW_PF_block RW_lkj_corr_cholesky sampler_RW_lkj_corr_cholesky RW_block_lkj_corr_cholesky sampler_RW_block_lkj_corr_cholesky 
#'
#' @examples
#' ## y[1] ~ dbern() or dbinom():
#' # mcmcConf$addSampler(target = 'y[1]', type = 'binary')   
#' 
#' # mcmcConf$addSampler(target = 'a', type = 'RW',
#' #    control = list(log = TRUE, adaptive = FALSE, scale = 3))
#' # mcmcConf$addSampler(target = 'b', type = 'RW',
#' #    control = list(adaptive = TRUE, adaptInterval = 200))
#' # mcmcConf$addSampler(target = 'p', type = 'RW',
#' #    control = list(reflective = TRUE))
#'
#' ## a, b, and c all continuous-valued:
#' # mcmcConf$addSampler(target = c('a', 'b', 'c'), type = 'RW_block')   
#' 
#' # mcmcConf$addSampler(target = 'p', type = 'RW_llFunction',
#' #    control = list(llFunction = RllFun, includesTarget = FALSE))
#' 
#' # mcmcConf$addSampler(target = 'y[1]', type = 'slice',
#' #    control = list(adaptive = FALSE, sliceWidth = 3))
#' # mcmcConf$addSampler(target = 'y[2]', type = 'slice',
#' #    control = list(adaptive = TRUE, sliceMaxSteps = 1))
#' 
#' # mcmcConf$addSampler(target = 'x[1:10]', type = 'ess')   ## x[1:10] ~ dmnorm()
#' 
#' # mcmcConf$addSampler(target = 'p[1:5]', type = 'RW_dirichlet')   ## p[1:5] ~ ddirch()
#'
#' ## y[1] is a posterior predictive node:
#' # mcmcConf$addSampler(target = 'y[1]', type = 'posterior_predictive')   
#' 
#' @seealso \code{\link{configureMCMC}} \code{\link{addSampler}} \code{\link{buildMCMC}} \code{\link{runMCMC}}
#'
#' @author Daniel Turek
#'
#' @references
#'
#' Andrieu, C., Doucet, A., and Holenstein, R. (2010). Particle Markov Chain Monte Carlo Methods. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 72(3), 269-342.
#'
#' Hoffman, Matthew D., and Gelman, Andrew (2014). The No-U-Turn Sampler: Adaptively setting path lengths in Hamiltonian Monte Carlo. \emph{Journal of Machine Learning Research}, 15(1): 1593-1623.
#'
#' Escobar, M. D., and West, M. (1995). Bayesian density estimation and inference using mixtures. \emph{Journal of the American Statistical Association}, 90(430), 577-588.
#'
#' Knorr-Held, L. and Rue, H. (2003). On block updating in Markov random field models for disease mapping. \emph{Scandinavian Journal of Statistics}, 29, 597-614.
#'
#' Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., and Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. \emph{The Journal of Chemical Physics}, 21(6), 1087-1092.
#'
#' Murray, I., Prescott Adams, R., and MacKay, D. J. C. (2010). Elliptical Slice Sampling. \emph{arXiv e-prints}, arXiv:1001.0175.
#'
#' Neal, R. M. (2000). Markov chain sampling methods for Dirichlet process mixture models. \emph{Journal of Computational and Graphical Statistics}, 9(2), 249-265.
#' 
#' Neal, R. M. (2003). Slice Sampling. \emph{The Annals of Statistics}, 31(3), 705-741.
#'
#' Neal, R. M. (2011). MCMC Using Hamiltonian Dynamics. \emph{Handbook of Markov Chain Monte Carlo}, CRC Press, 2011.
#'
#' Pitt, M. K. and Shephard, N. (1999). Filtering via simulation: Auxiliary particle filters. \emph{Journal of the American Statistical Association} 94(446), 590-599.
#'
#' Polson, N.G., Scott, J.G., and J. Windle. (2013). Bayesian inference for logistic models using Pólya-gamma latent variables. Journal of the American Statistical Association, 108(504), 1339–1349. https://doi.org/10.1080/01621459.2013.829001
#'
#' Roberts, G. O. and S. K. Sahu (1997). Updating Schemes, Correlation Structure, Blocking and Parameterization for the Gibbs Sampler. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 59(2), 291-317.
#'
#' Shaby, B. and M. Wells (2011). \emph{Exploring an Adaptive Metropolis Algorithm}. 2011-14. Department of Statistics, Duke University.
#'
#' Stan Development Team (2020). \emph{Stan Language Reference Manual, Version 2.22, Section 10.12}.
#'
#' Tibbits, M. M.,  Groendyke, C.,  Haran, M., and Liechty, J. C. (2014).  Automated Factor Slice Sampling.  \emph{Journal of Computational and Graphical Statistics}, 23(2), 543-563.
#'
#' van Dyk, D.A. and T. Park. (2008). Partially collapsed Gibbs Samplers. \emph{Journal of the American Statistical Association}, 103(482), 790-796.
#'
#' Yu, Y. and Meng, X. L. (2011). To center or not to center: That is not the question - An ancillarity-sufficiency interweaving strategy (ASIS) for boosting MCMC efficiency. Journal of Computational and Graphical Statistics, 20(3), 531–570. https://doi.org/10.1198/jcgs.2011.203main
#' 
NULL




##
## Note: RW_multinomial sampler was removed from package, following version 0.1.0
##
##@section RW_multinomial sampler:
##
##This sampler is designed for sampling multinomial target distributions.  The sampler performs a series of Metropolis-Hastings steps between pairs of groups.  Proposals are generated via a draw from a binomial distribution, whereafter the proposed number density is moved from one group to another group.  The acceptance or rejection of these proposals follows a standard Metropolis-Hastings procedure.  Probabilities for the random binomial proposals are adapted to a target acceptance rate of 0.5.
##
##The \code{RW_multinomial} sampler accepts the following control list elements:
##\itemize{
##\item adaptive.  A logical argument, specifying whether the sampler should adapt the binomial proposal probabilities throughout the course of MCMC execution. (default = TRUE)
##\item adaptInterval.  The interval on which to perform adaptation.  A minimum value of 100 is required. (default = 200)
##}
