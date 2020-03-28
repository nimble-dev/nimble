


####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################

#' @rdname samplers
#' @export
sampler_BASE <- nimbleFunctionVirtual(
    methods = list(
        reset = function() { }
    )
)



####################################################################
### posterior_predictive sampler for trailing stoch. nodes #########
####################################################################

#' @rdname samplers
#' @export
sampler_posterior_predictive <- nimbleFunction(
    name = 'sampler_posterior_predictive',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        calcNodes <- model$getDependencies(target)
    },
    run = function() {
        simulate(model, target)
        calculate(model, calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
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
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## checks
        if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
        if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
    },
    run = function() {
        currentLogProb <- getLogProb(model, calcNodes)
        model[[target]] <<- 1 - model[[target]]
        otherLogProbPrior <- calculate(model, target)
        if(otherLogProbPrior == -Inf) {
            otherLogProb <- otherLogProbPrior
        } else {
            otherLogProb <- otherLogProbPrior + calculate(model, calcNodesNoSelf)
        }
        acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
        if(!is.nan(acceptanceProb) & runif(1,0,1) < acceptanceProb)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
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
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes  <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        k <- length(model$getParam(target, 'prob'))
        probs <- numeric(k)
        logProbs <- numeric(k)
        ## checks
        if(length(targetAsScalar) > 1)  stop('cannot use categorical sampler on more than one target node')
        if(model$getDistribution(target) != 'dcat') stop('can only use categorical sampler on node with dcat distribution')
    },
    run = function() {
        currentValue <- model[[target]]
        logProbs[currentValue] <<- getLogProb(model, calcNodes)
        for(i in 1:k) {
            if(i != currentValue) {
                model[[target]] <<- i
                logProbPrior <- calculate(model, target)
                if(logProbPrior == -Inf) {
                    logProbs[i] <<- -Inf
                } else {
                    if(is.nan(logProbPrior)) {
                        logProbs[i] <<- -Inf
                    } else {
                        logProbs[i] <<- logProbPrior + calculate(model, calcNodesNoSelf)
                        if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
                    }
                }
            }
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs) #rcat normalizes the probs internally
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            calculate(model, calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
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
        logScale            <- if(!is.null(control$log))                 control$log                 else FALSE
        reflective          <- if(!is.null(control$reflective))          control$reflective          else FALSE
        adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
        adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 200
        adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
        scale               <- if(!is.null(control$scale))               control$scale               else 1
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory  <- c(0, 0)   ## scaleHistory
        acceptanceHistory  <- c(0, 0)   ## scaleHistory
        if(nimbleOptions('MCMCsaveHistory')) {
            saveMCMChistory <- TRUE
        } else saveMCMChistory <- FALSE
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
        logMHR <- calculateDiff(model, target)
        if(logMHR == -Inf) {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            ## Drawing a random number is needed during first testing
            ## of this step in order to keep the random numbers identical
            ## to old behavior to see if tests that depend on particular
            ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
            ## we call decide() to ensure same behavior.
            ## jump <- decide(logMHR)
            ## When new behavior is acceptable, we can remove the above line
            ## and uncomment the following:
            jump <- FALSE
        } else {
            logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf) + propLogScale
            jump <- decide(logMHR)
            if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
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
                    setSize(scaleHistory, timesAdapted)         ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale        ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)         ## scaleHistory
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
        getScaleHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(scaleHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
                return(numeric(1, 0))
            }
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(acceptanceHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
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
    ), where = getLoadingNamespace()
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
        adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
        adaptScaleOnly      <- if(!is.null(control$adaptScaleOnly))      control$adaptScaleOnly      else FALSE
        adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 200
        adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
        scale               <- if(!is.null(control$scale))               control$scale               else 1
        propCov             <- if(!is.null(control$propCov))             control$propCov             else 'identity'
        tries               <- if(!is.null(control$tries))               control$tries               else 1
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
        if(!is.integer(finalTargetIndex) |
           length(finalTargetIndex) != 1 |
           is.na(finalTargetIndex[1]))
            stop('Problem with target node in sampler_RW_block')
        calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
        calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
#        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        d <- length(targetAsScalar)
        scaleHistory  <- c(0, 0)                                                 ## scaleHistory
        acceptanceHistory  <- c(0, 0)                                            ## scaleHistory
        propCovHistory <- if(d<=10) array(0, c(2,d,d)) else array(0, c(2,2,2))   ## scaleHistory
        saveMCMChistory <- if(nimbleOptions('MCMCsaveHistory')) TRUE else FALSE
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
        ## nested function and function list definitions
        ##        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
##        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
        ## checks
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    },
    run = function() {
        for(i in 1:tries) {
            propValueVector <- generateProposalVector()
            ##        lpMHR <- my_setAndCalculateDiff$run(propValueVector)
            values(model, targetNodesAsScalar) <<- propValueVector
            lpD <- calculateDiff(model, calcNodesProposalStage)
            if(lpD == -Inf) {
                nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
            ## Drawing a random number is needed during first testing
            ## of this step in order to keep the random numbers identical
            ## to old behavior to see if tests that depend on particular
            ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
            ## we call decide() to ensure same behavior.
            ## jump <- decide(lpD)
            ## When new behavior is acceptable, we can remove the above line
            ## and uncomment the following:
            jump <- FALSE
            } else {
                ##        jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## will use lpMHR - 0
                lpD <- lpD + calculateDiff(model, calcNodesDepStage)
                jump <- decide(lpD)
                if(jump) { nimCopy(from = model,   to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
                } else   { nimCopy(from = mvSaved, to = model,   row = 1, nodes = calcNodes, logProb = TRUE) }
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
                    if(d <= 10) {
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
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
            returnType(double(1))
            return(scaleHistory)
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(!saveMCMChistory)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
            return(acceptanceHistory)
        },                  
        getPropCovHistory = function() { ## scaleHistory
            if(!saveMCMChistory | d > 10)   print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC and note that to reduce memory use we only save the proposal covariance history for parameter vectors of length 10 or less")
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
                if(d <= 10) 
                    propCovHistory <<- nimArray(0, dim = c(2,d,d))
            }
            my_calcAdaptationFactor$reset()
        }
    ), where = getLoadingNamespace()
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
        adaptive       <- if(!is.null(control$adaptive))       control$adaptive       else TRUE
        adaptInterval  <- if(!is.null(control$adaptInterval))  control$adaptInterval  else 200
        scale          <- if(!is.null(control$scale))          control$scale          else 1
        llFunction     <- if(!is.null(control$llFunction))     control$llFunction     else stop('RW_llFunction sampler missing required control argument: llFunction')
        includesTarget <- if(!is.null(control$includesTarget)) control$includesTarget else stop('RW_llFunction sampler missing required control argument: includesTarget')
        ## node list generation
        calcNodes <- model$getDependencies(target)
        ## nested function and function list definitions
        mvInternal <- modelValues(model)
        RWControl <- list(adaptive=adaptive, adaptInterval=adaptInterval, scale=scale, log=FALSE, reflective=FALSE)
        targetRWSamplerFunction <- sampler_RW(model, mvInternal, target, RWControl)
        my_setAndCalculateOne <- setAndCalculateOne(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    },
    run = function() {
        modelLP0 <- llFunction$run()
        if(!includesTarget)     modelLP0 <- modelLP0 + getLogProb(model, target)
        propValue <- rnorm(1, mean = model[[target]], sd = scale)
        my_setAndCalculateOne$run(propValue)
        modelLP1 <- llFunction$run()
        if(!includesTarget)     modelLP1 <- modelLP1 + getLogProb(model, target)
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
    ), where = getLoadingNamespace()
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
        adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
        adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200
        width         <- if(!is.null(control$sliceWidth))    control$sliceWidth    else 1
        maxSteps      <- if(!is.null(control$sliceMaxSteps)) control$sliceMaxSteps else 100
        maxContractions        <- if(!is.null(control$maxContractions))
                                      control$maxContractions else 1000
        maxContractionsWarning <- if(!is.null(control$maxContractionsWarning))
                                      control$maxContractionsWarning else TRUE
        eps <- 1e-15
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        widthOriginal <- width
        timesRan      <- 0
        timesAdapted  <- 0
        sumJumps      <- 0
        discrete      <- model$isDiscrete(target)
        ## checks
        if(length(targetAsScalar) > 1)     stop('cannot use slice sampler on more than one target node')
    },
    run = function() {
        u <- getLogProb(model, calcNodes) - rexp(1, 1)    # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
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
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            jumpDist <- abs(x1 - x0)
            if(adaptive)     adaptiveProcedure(jumpDist)
        }
    },
    methods = list(
        setAndCalculateTarget = function(value = double()) {
            if(discrete)     value <- floor(value)
            model[[target]] <<- value
            lp <- calculate(model, target)
            if(lp == -Inf) return(-Inf) 
            lp <- lp + calculate(model, calcNodesNoSelf)
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
                timesRan <<- 0
                sumJumps <<- 0
            }
        },
        reset = function() {
            width        <<- widthOriginal
            timesRan     <<- 0
            timesAdapted <<- 0
            sumJumps     <<- 0
        }
    ), where = getLoadingNamespace()
)



####################################################################
### elliptical slice sampler (dmnorm distributions) ################
####################################################################


#' @rdname samplers
#' @export
sampler_ess <- nimbleFunction(
    name = 'sampler_ess',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        maxContractions        <- if(!is.null(control$maxContractions))
                                      control$maxContractions else 1000
        maxContractionsWarning <- if(!is.null(control$maxContractionsWarning))
                                      control$maxContractionsWarning else TRUE
        eps <- 1e-15
        ## node list generation
        target <- model$expandNodeNames(target)
        calcNodes <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        Pi <- pi
        ## nested function and function list definitions
        ##target_nodeFunctionList <- nimbleFunctionList(node_stoch_dmnorm)
        ##target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
        ## checks
        if(length(target) > 1)                          stop('elliptical slice sampler only applies to one target node')
        if(model$getDistribution(target) != 'dmnorm')   stop('elliptical slice sampler only applies to multivariate normal distributions')
    },
    run = function() {
        u <- getLogProb(model, calcNodes) - rexp(1, 1)
        target_mean <- model$getParam(target, 'mean') ##target_nodeFunctionList[[1]]$get_mean()
        f <- model[[target]] - target_mean
        simulate(model, target)
        nu <- model[[target]] - target_mean
        theta <- runif(1, 0, 2*Pi)
        theta_min <- theta - 2*Pi
        theta_max <- theta
        model[[target]] <<- f*cos(theta) + nu*sin(theta) + target_mean
        lp <- calculate(model, calcNodes)
        numContractions <- 0
        while((is.nan(lp) | lp < u) & theta_max - theta_min > eps & numContractions < maxContractions) {   # must be is.nan()
            ## The checks for theta_max - theta_min small and max number of contractions are
            ## for cases where model is in invalid state and lp calculations are NA/NaN or where
            ## theta interval contracts to zero
            if(theta < 0)   theta_min <- theta   else   theta_max <- theta
            theta <- runif(1, theta_min, theta_max)
            model[[target]] <<- f*cos(theta) + nu*sin(theta) + target_mean
            lp <- calculate(model, calcNodes)
            numContractions <- numContractions + 1
        }
        if(theta_max - theta_min <= eps | numContractions == maxContractions) {
            if(maxContractionsWarning)
                cat("Warning: elliptical slice sampler reached maximum number of contractions.\n")
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        } else
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
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
        widthVec            <- if(!is.null(control$sliceWidths))              control$sliceWidths              else 'oneVec'
        maxSteps            <- if(!is.null(control$sliceMaxSteps))            control$sliceMaxSteps            else 100
        adaptFactorMaxIter  <- if(!is.null(control$sliceAdaptFactorMaxIter))  control$sliceAdaptFactorMaxIter  else 15000
        adaptFactorInterval <- if(!is.null(control$sliceAdaptFactorInterval)) control$sliceAdaptFactorInterval else 1000
        adaptWidthMaxIter   <- if(!is.null(control$sliceAdaptWidthMaxIter))   control$sliceAdaptWidthMaxIter   else 512
        adaptWidthTolerance <- if(!is.null(control$sliceAdaptWidthTolerance)) control$sliceAdaptWidthTolerance else 0.1
        maxContractions     <- if(!is.null(control$maxContractions))          control$maxContractions else 1000
        maxContractionsWarning <- if(!is.null(control$maxContractionsWarning))
                                      control$maxContractionsWarning else TRUE
        eps <- 1e-15
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target)
        finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
        if(!is.integer(finalTargetIndex) |
           length(finalTargetIndex) != 1 |
           is.na(finalTargetIndex[1]))
            stop('Problem with target node in sampler_AF_slice')
        calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
        calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
##        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
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
        if(class(widthVec) != 'numeric')   stop('sliceWidths must be a numeric vector')
        if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)
    },
    run = function() {
        maxContractionsReached <- FALSE
        for(i in 1:d) {
            eigenVec <- gammaMatrix[, i]
            width <- widthVec[i]
            u <- getLogProb(model, calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
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
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        } else
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
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
    ), where = getLoadingNamespace()
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
    }, where = getLoadingNamespace()
)

#' @rdname samplers
#' @export
sampler_crossLevel <- nimbleFunction(
    name = 'sampler_crossLevel',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive <- if(!is.null(control$adaptive)) control$adaptive else TRUE
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
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    },
    run = function() {
        modelLP0 <- getLogProb(model, calcNodes)
        propLP0 <- 0
        for(iSF in seq_along(lowConjugateGetLogDensityFunctions))  { propLP0 <- propLP0 + lowConjugateGetLogDensityFunctions[[iSF]]$run() }
        propValueVector <- topRWblockSamplerFunction$generateProposalVector()
        topLP <- my_setAndCalculateTop$run(propValueVector)
        if(is.na(topLP))
            jump <- my_decideAndJump$run(-Inf, 0, 0, 0)
        else {
            for(iSF in seq_along(lowConjugateSamplerFunctions))
                lowConjugateSamplerFunctions[[iSF]]$run()
            modelLP1 <- calculate(model, calcNodes)
            propLP1 <- 0
            for(iSF in seq_along(lowConjugateGetLogDensityFunctions))
                propLP1 <- propLP1 + lowConjugateGetLogDensityFunctions[[iSF]]$run()
            jump <- my_decideAndJump$run(modelLP1, modelLP0, propLP1, propLP0)
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
    ), where = getLoadingNamespace()
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
        adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
        adaptScaleOnly      <- if(!is.null(control$adaptScaleOnly))      control$adaptScaleOnly      else FALSE
        adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 200
        adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
        scale               <- if(!is.null(control$scale))               control$scale               else 1
        propCov             <- if(!is.null(control$propCov))             control$propCov             else 'identity'
        llFunction          <- if(!is.null(control$llFunction))          control$llFunction          else stop('RW_llFunction_block sampler missing required control argument: llFunction')
        includesTarget      <- if(!is.null(control$includesTarget))      control$includesTarget      else stop('RW_llFunction_block sampler missing required control argument: includesTarget')
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
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
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
        ## checks
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    },
    run = function() {
        modelLP0 <- llFunction$run()
        if(!includesTarget)     modelLP0 <- modelLP0 + getLogProb(model, target)
        propValueVector <- generateProposalVector()
        my_setAndCalculate$run(propValueVector)
        modelLP1 <- llFunction$run()
        if(!includesTarget)     modelLP1 <- modelLP1 + getLogProb(model, target)
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
    ), where = getLoadingNamespace()
)



#######################################################################################
### RW_PF, does a univariate RW, but using a particle filter likelihood function ######
#######################################################################################

#' @rdname samplers
#' @export
sampler_RW_PF <- nimbleFunction(
    name = 'sampler_RW_PF',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- if(!is.null(control$adaptive))             control$adaptive             else TRUE
        adaptInterval  <- if(!is.null(control$adaptInterval))        control$adaptInterval        else 200
        scale          <- if(!is.null(control$scale))                control$scale                else 1
        m              <- if(!is.null(control$pfNparticles))         control$pfNparticles         else 1000
        existingPF     <- if(!is.null(control$pf))                   control$pf                   else NULL
        filterType     <- if(!is.null(control$pfType))               control$pfType               else 'bootstrap'
        filterControl  <- if(!is.null(control$pfControl))            control$pfControl            else list()
        optimizeM      <- if(!is.null(control$pfOptimizeNparticles)) control$pfOptimizeNparticles else FALSE
        latents        <- if(!is.null(control$latents))              control$latents              else stop('RW_PF sampler missing required control argument: latents')
        
        if(!is.null(control$pfLookahead)) {
          print("Warning, the `pfLookahead` control list argument is deprecated
                and will not be supported in future versions of NIMBLE. Please
                specify the lookahead function via the pfControl argument 
                instead.")
          filterControl$lookahead  <-  control$pfLookahead
        }                    
        else if(is.null(filterControl$lookahead)) {
          filterControl$lookahead  <-  'simulate'
        } 
        
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        latentSamp <- FALSE
        MCMCmonitors <- tryCatch(parent.frame(2)$conf$monitors, error = function(e) e)
        if(identical(MCMCmonitors, TRUE))
            latentSamp <- TRUE
        else if(any(model$expandNodeNames(latents) %in% model$expandNodeNames(MCMCmonitors)))
            latentSamp <- TRUE
        latentDep <- model$getDependencies(latents)
        topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)
        ## numeric value generation
        optimizeM       <- as.integer(optimizeM)
        scaleOriginal   <- scale
        timesRan        <- 0
        timesAccepted   <- 0
        timesAdapted    <- 0
        prevLL          <- 0
        nVarEsts        <- 0
        itCount         <- 0
        optimalAR       <- 0.44
        gamma1          <- 0
        storeParticleLP <- -Inf
        storeLLVar      <- 0
        ## Number of LL estimates to compute to get each LL
        ## variance estimate for m optimization.
        nVarReps        <- 7   
        ## Number of LL variance estimates to compute before deciding optimal m.
        mBurnIn         <- 15   
        d               <- length(targetAsScalar)
        if(optimizeM) m <- 3000
        ## Nested function and function list definitions.
        my_setAndCalculate <- setAndCalculateOne(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        if(!is.null(existingPF)) {
            my_particleFilter <- existingPF
        } else {
            if(latentSamp == TRUE) { 
                filterControl$saveAll <- TRUE
                filterControl$smoothing <- TRUE
            } else {
                filterControl$saveAll <- FALSE
                filterControl$smoothing <- FALSE
            }
            filterControl$initModel <- FALSE
            if(is.character(filterType) && filterType == 'auxiliary') {
                my_particleFilter <- buildAuxiliaryFilter(model, latents, 
                                                          control = filterControl)
            }
            else if(is.character(filterType) && filterType == 'bootstrap') {
                my_particleFilter <- buildBootstrapFilter(model, latents,
                                                          control = filterControl)
            }
            else if(is.nfGenerator(filterType)){
                my_particleFilter <- filterType(model, latents,
                                                control = filterControl)
            }
            else stop('filter type must be either "bootstrap", "auxiliary", or a
                  user defined filtering algorithm created by a call to 
                  nimbleFunction(...).')
        }
        particleMV <- my_particleFilter$mvEWSamples
        ## checks
        if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
        if(length(targetAsScalar) > 1)                      stop('more than one top-level target; cannot use RW_PF sampler, try RW_PF_block sampler')
    },
    run = function() {
        storeParticleLP <<- my_particleFilter$getLastLogLik()
        modelLP0 <- storeParticleLP + getLogProb(model, target)
        propValue <- rnorm(1, mean = model[[target]], sd = scale)
        my_setAndCalculate$run(propValue)
        particleLP <- my_particleFilter$run(m)
        modelLP1 <- particleLP + getLogProb(model, target)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
        if(!jump) {
            my_particleFilter$setLastLogLik(storeParticleLP)
        }
        if(jump & latentSamp){
            ## if we jump, randomly sample latent nodes from pf output and put into model so that they can be monitored
            index <- ceiling(runif(1, 0, m))
            copy(particleMV, model, latents, latents, index)
            calculate(model, latentDep)
            copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
        }
        else if(!jump & latentSamp){
            ## if we don't jump, replace model latent nodes with saved latent nodes
            copy(from = mvSaved, to = model, nodes = latentDep, row = 1, logProb = TRUE)
        }
##        if(jump & !resample)  storeParticleLP <<- particleLP
        if(jump & optimizeM) optimM()
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        optimM = function() {
            tempM <- 15000
            declare(LLEst, double(1, nVarReps))
            if(nVarEsts < mBurnIn) {  # checks whether we have enough var estimates to get good approximation
                for(i in 1:nVarReps)
                    LLEst[i] <- my_particleFilter$run(tempM)
                ## next, store average of var estimates
                if(nVarEsts == 1)
                    storeLLVar <<- var(LLEst)/mBurnIn
                else {
                    LLVar <- storeLLVar
                    LLVar <- LLVar + var(LLEst)/mBurnIn
                    storeLLVar<<- LLVar
                }
                nVarEsts <<- nVarEsts + 1
            }
            else {  # once enough var estimates have been taken, use their average to compute m
                m <<- m*storeLLVar/(0.92^2)
                m <<- ceiling(m)
                storeParticleLP <<- my_particleFilter$run(m)
                optimizeM <<- 0
            }
        },
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
            storeParticleLP <<- -Inf
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)



#######################################################################################
### RW_PF_block, does a block RW, but using a particle filter likelihood function #####
#######################################################################################

#' @rdname samplers
#' @export
sampler_RW_PF_block <- nimbleFunction(
    name = 'sampler_RW_PF_block',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target,  control) {
        ## control list extraction
        adaptive            <- if(!is.null(control$adaptive))             control$adaptive             else TRUE
        adaptScaleOnly      <- if(!is.null(control$adaptScaleOnly))       control$adaptScaleOnly       else FALSE
        adaptInterval       <- if(!is.null(control$adaptInterval))        control$adaptInterval        else 200
        adaptFactorExponent <- if(!is.null(control$adaptFactorExponent))  control$adaptFactorExponent  else 0.8
        scale               <- if(!is.null(control$scale))                control$scale                else 1
        propCov             <- if(!is.null(control$propCov))              control$propCov              else 'identity'
        existingPF          <- if(!is.null(control$pf))                   control$pf                   else NULL
        m                   <- if(!is.null(control$pfNparticles))         control$pfNparticles         else 1000
        filterType          <- if(!is.null(control$pfType))               control$pfType               else 'bootstrap'
        filterControl       <- if(!is.null(control$pfControl))            control$pfControl            else list()
        optimizeM           <- if(!is.null(control$pfOptimizeNparticles)) control$pfOptimizeNparticles else FALSE
        latents             <- if(!is.null(control$latents))              control$latents              else stop('RW_PF sampler missing required control argument: latents')
        
        if(!is.null(control$pfLookahead)) {
          print("Warning, the `pfLookahead` control list argument is deprecated
                and will not be supported in future versions of NIMBLE. Please
                specify the lookahead function via the pfControl argument 
                instead.")
          filterControl$lookahead  <-  control$pfLookahead
        }                    
        else if(is.null(filterControl$lookahead)) {
          filterControl$lookahead  <-  'simulate'
        } 
        
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        latentSamp <- FALSE
        MCMCmonitors <- tryCatch(parent.frame(2)$conf$monitors, error = function(e) e)
        if(identical(MCMCmonitors, TRUE))
            latentSamp <- TRUE
        else if(any(model$expandNodeNames(latents) %in% model$expandNodeNames(MCMCmonitors)))
            latentSamp <- TRUE
        latentDep <- model$getDependencies(latents)
        topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)
        target <- model$expandNodeNames(target)
        ## numeric value generation
        optimizeM     <- as.integer(optimizeM)
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        prevLL        <- 0
        nVarEsts      <- 0
        itCount       <- 0
        d <- length(targetAsScalar)
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
        storeParticleLP <- -Inf
        storeLLVar  <- 0
        nVarReps <- 7    # number of LL estimates to compute to get each LL variance estimate for m optimization
        mBurnIn  <- 15   # number of LL variance estimates to compute before deciding optimal m
        if(optimizeM)   m <- 3000
        ## nested function and function list definitions
        my_setAndCalculate <- setAndCalculate(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)
        if(!is.null(existingPF)) {
            my_particleFilter <- existingPF
        } else {
            if(latentSamp == TRUE) { 
                filterControl$saveAll <- TRUE
                filterControl$smoothing <- TRUE
            } else {
                filterControl$saveAll <- FALSE
                filterControl$smoothing <- FALSE
            }
            filterControl$initModel <- FALSE
            if(is.character(filterType) && filterType == 'auxiliary') {
                my_particleFilter <- buildAuxiliaryFilter(model, latents, 
                                                          control = filterControl)
            }
            else if(is.character(filterType) && filterType == 'bootstrap') {
                my_particleFilter <- buildBootstrapFilter(model, latents,
                                                          control = filterControl)
            }
            else if(is.nfGenerator(filterType)){
                my_particleFilter <- filterType(model, latents,
                                                control = filterControl)
                
            }
            else stop('filter type must be either "bootstrap", "auxiliary", or a
                  user defined filtering algorithm created by a call to 
                  nimbleFunction(...).')
        }
        particleMV <- my_particleFilter$mvEWSamples
        ## checks
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
        if(length(targetAsScalar) < 2)        stop('less than two top-level targets; cannot use RW_PF_block sampler, try RW_PF sampler')
        if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
    },
    run = function() {
        storeParticleLP <<- my_particleFilter$getLastLogLik()
        modelLP0 <- storeParticleLP + getLogProb(model, target)
        propValueVector <- generateProposalVector()
        my_setAndCalculate$run(propValueVector)
        particleLP <- my_particleFilter$run(m)
        modelLP1 <- particleLP + getLogProb(model, target)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
        if(!jump) {
            my_particleFilter$setLastLogLik(storeParticleLP)
        }
        if(jump & latentSamp) {
            ## if we jump, randomly sample latent nodes from pf output and put
            ## into model so that they can be monitored
            index <- ceiling(runif(1, 0, m))
            copy(particleMV, model, latents, latents, index)
            calculate(model, latentDep)
            copy(from = model, to = mvSaved, nodes = latentDep, row = 1, logProb = TRUE)
        }
        else if(!jump & latentSamp) {
            ## if we don't jump, replace model latent nodes with saved latent nodes
            copy(from = mvSaved, to = model, nodes = latentDep, row = 1, logProb = TRUE)
        }
      ##  if(jump & !resample)  storeParticleLP <<- particleLP
        if(jump & optimizeM) optimM()
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        optimM = function() {
            tempM <- 15000
            declare(LLEst, double(1, nVarReps))
            if(nVarEsts < mBurnIn) {  # checks whether we have enough var estimates to get good approximation
                for(i in 1:nVarReps)
                    LLEst[i] <- my_particleFilter$run(tempM)
                ## next, store average of var estimates
                if(nVarEsts == 1)
                    storeLLVar <<- var(LLEst)/mBurnIn
                else {
                    LLVar <- storeLLVar
                    LLVar <- LLVar + var(LLEst)/mBurnIn
                    storeLLVar <<- LLVar
                }
                nVarEsts <<- nVarEsts + 1
            }
            else {  # once enough var estimates have been taken, use their average to compute m
                m <<- m*storeLLVar/(0.92^2)
                m <<- ceiling(m)
                storeParticleLP <<- my_particleFilter$run(m)
                optimizeM <<- 0
            }
        },
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
            chol_propCov_scale <<- chol_propCov * scale
            storeParticleLP <<- -Inf
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            my_calcAdaptationFactor$reset()
        }
    ), where = getLoadingNamespace()
)



#######################################################################################
### RW_multinomial sampler for multinomial distributions ##############################
#######################################################################################

#' @rdname samplers
#' @export
sampler_RW_multinomial <- nimbleFunction(
    name = 'sampler_RW_multinomial',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
        adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        targetAllNodes <- unique(model$expandNodeNames(target))
        calcNodes      <- model$getDependencies(target) 
        lTarget        <- length(targetAsScalar)
        Ntotal         <- sum(values(model,target))
        NOverL         <- Ntotal / lTarget
        ## numeric value generation
        Zeros             <- matrix(0, lTarget, lTarget)
        Ones              <- matrix(1, lTarget, lTarget)
        timesRan          <- Zeros
        AcceptRates       <- Zeros
        ScaleShifts       <- Zeros
        totalAdapted      <- Zeros
        timesAccepted     <- Zeros
        ENSwapMatrix      <- Ones
        ENSwapDeltaMatrix <- Ones
        RescaleThreshold  <- 0.2 * Ones
        lpProp  <- 0
        lpRev   <- 0
        Pi      <- pi 
        PiOver2 <- Pi / 2 ## Irrational number prevents recycling becoming degenerate
        u       <- runif(1, 0, Pi)
        ## nested function and function list definitions
        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        my_decideAndJump       <- decideAndJump(model, mvSaved, calcNodes)
        ## checks
        if(model$getDistribution(target) != 'dmulti')   stop('can only use RW_multinomial sampler for multinomial distributions')
        if(length(targetAllNodes) > 1)                      stop('cannot use RW_multinomial sampler on more than one target')
        if(adaptive & adaptInterval < 100)                  stop('adaptInterval < 100 is not recommended for RW_multinomial sampler')
    },
    run = function() {
        for(iFROM in 1:lTarget) {            
            for(iTO in 1:(lTarget-1)) {
                if(u > PiOver2) {                
                    iFrom <- iFROM
                    iTo   <- iTO
                    if (iFrom == iTo)
                        iTo <- lTarget
                    u <<- 2 * (u - PiOver2)   # recycle u
                } else {
                    iFrom <- iTO
                    iTo   <- iFROM
                    if (iFrom == iTo)
                        iFrom <- lTarget
                    u <<- 2 * (PiOver2 - u)   # recycle u
                }
                propValueVector <- generateProposalVector(iFrom, iTo)
                lpMHR <- my_setAndCalculateDiff$run(propValueVector) + lpRev - lpProp 
                jump  <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## returns lpMHR + 0 - 0 + 0
                if(adaptive)   adaptiveProcedure(jump=jump, iFrom=iFrom, iTo=iTo)
            }
        }
    },
    methods = list(
        generateProposalVector = function(iFrom = integer(), iTo = integer()) { 
            propVector <- values(model,target) 
            pSwap      <- min(1, max(1, ENSwapMatrix[iFrom,iTo]) / propVector[iFrom]) 
            nSwap      <- rbinom(n=1,   size=propVector[iFrom], prob=pSwap) 
            lpProp    <<- dbinom(nSwap, size=propVector[iFrom], prob=pSwap, log=TRUE) 
            propVector[iFrom] <- propVector[iFrom] - nSwap 
            propVector[iTo]   <- propVector[iTo]   + nSwap 
            pRevSwap   <- min(1, max(1, ENSwapMatrix[iTo,iFrom]) / (propVector[iTo] + nSwap)) 
            lpRev     <<- dbinom(nSwap, size=propVector[iTo], prob=pRevSwap, log=TRUE) 
            returnType(double(1)) 
            return(propVector) 
        },
        adaptiveProcedure = function(jump=logical(), iFrom=integer(), iTo=integer()) {
            NVector <- values(model,target) 
            timesRan[iFrom, iTo] <<- timesRan[iFrom, iTo] + 1
            if(jump)
                timesAccepted[iFrom, iTo] <<- timesAccepted[iFrom, iTo] + 1
            if (timesRan[iFrom, iTo] %% adaptInterval == 0) {
                totalAdapted[iFrom, iTo] <<- totalAdapted[iFrom, iTo] + 1
                accRate                   <- timesAccepted[iFrom, iTo] / timesRan[iFrom, iTo]
                AcceptRates[iFrom, iTo]  <<- accRate
                if (accRate > 0.5) {
                    ENSwapMatrix[iFrom, iTo] <<-
                        min(Ntotal,
                            ENSwapMatrix[iFrom,iTo] + ENSwapDeltaMatrix[iFrom, iTo] / totalAdapted[iFrom,iTo])
                } else {
                    ENSwapMatrix[iFrom, iTo] <<-
                        max(1,
                            ENSwapMatrix[iFrom,iTo] - ENSwapDeltaMatrix[iFrom,iTo] / totalAdapted[iFrom,iTo])
                } 
                if(accRate<RescaleThreshold[iFrom,iTo] | accRate>(1-RescaleThreshold[iFrom,iTo])) {
                    ## rescale iff ENSwapMatrix[iFrom, iTo] is not set to an upper or lower bound 
                    if (ENSwapMatrix[iFrom, iTo] > 1 & ENSwapMatrix[iFrom, iTo] < Ntotal) {
                        ScaleShifts[iFrom, iTo]       <<- ScaleShifts[iFrom, iTo] + 1 
                        ENSwapDeltaMatrix[iFrom, iTo] <<- min(NOverL, ENSwapDeltaMatrix[iFrom, iTo] * totalAdapted[iFrom,iTo] / 10)
                        ENSwapDeltaMatrix[iTo, iFrom] <<- ENSwapDeltaMatrix[iFrom, iTo] 
                        RescaleThreshold[iFrom,iTo]   <<- 0.2 * 0.95^ScaleShifts[iFrom, iTo]
                    }
                }
                ## lower Bound 
                if(ENSwapMatrix[iFrom, iTo] < 1)
                    ENSwapMatrix[iFrom, iTo] <<- 1                
                ## symmetry in ENSwapMatrix helps maintain good acceptance rates
                ENSwapMatrix[iTo,iFrom]   <<- ENSwapMatrix[iFrom,iTo]
                timesRan[iFrom, iTo]      <<- 0
                timesAccepted[iFrom, iTo] <<- 0
            }
        },
        reset = function() {
            timesRan          <<- Zeros
            AcceptRates       <<- Zeros
            ScaleShifts       <<- Zeros
            totalAdapted      <<- Zeros
            timesAccepted     <<- Zeros
            ENSwapMatrix      <<- Ones
            ENSwapDeltaMatrix <<- Ones
            RescaleThreshold  <<- 0.2 * Ones
        }
    ), where = getLoadingNamespace()
)



#####################################################################################
### RW_dirichlet sampler for multinomial distributions ##############################
#####################################################################################

#' @rdname samplers
#' @export
sampler_RW_dirichlet <- nimbleFunction(
    name = 'sampler_RW_dirichlet',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
        adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200
        scaleOriginal <- if(!is.null(control$scale))         control$scale         else 1
        ## node list generation
        calcNodes <- model$getDependencies(target)
        depNodes  <- model$getDependencies(target, self = FALSE)
        targetScalarNodes <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        ## numeric value generation
        d <- length(targetScalarNodes)
        thetaVec         <- rep(0, d)
        scaleVec         <- rep(scaleOriginal, d)
        timesRan         <- 0
        timesAcceptedVec <- rep(0, d)
        timesAdapted     <- 0
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
                logMHR <- alphaVec[i]*propLogScale + currentValue - propValue + calculateDiff(model, depNodes)
                jump <- decide(logMHR)
            } else jump <- FALSE
            if(adaptive & jump)   timesAcceptedVec[i] <<- timesAcceptedVec[i] + 1
            if(jump) { thetaVec <<- thetaVecProp
                       nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            } else   { nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE) }
            model$calculate(target)                                                         ## update target logProb
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)    ##
        }
        if(adaptive) {
            timesRan <<- timesRan + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRateVec <- timesAcceptedVec / timesRan
                timesAdapted <<- timesAdapted + 1
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                adaptFactorVec <- exp(10 * gamma1 * (acceptanceRateVec - 0.44))   ## optimalAR = 0.44
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
    ), where = getLoadingNamespace()
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
        adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
        adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 200
        adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
        scale               <- if(!is.null(control$scale))               control$scale               else 1
        propCov             <- if(!is.null(control$propCov))             control$propCov             else 'identity'
        ## node list generation
        target <- model$expandNodeNames(target)
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
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
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric')
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
        logMHR <- calculateDiff(model, calcNodes)
        deltaDiag <- thetaVec_prop[1:d]-thetaVec[1:d]
        for(i in 1:d)   logMHR <- logMHR + (d+2-i)*deltaDiag[i]  ## took me quite a while to derive this
        jump <- decide(logMHR)
        if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
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
    ), where = getLoadingNamespace()
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
        ## node list generation
        target <- model$expandNodeNames(target)
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)  

        d <- sqrt(length(targetAsScalar))
        nTheta <- d*(d-1)/2    # this many unconstrained elements to be sampled
        ## control list extraction
        adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
        adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 200
        adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
        scaleOriginal       <- if(!is.null(control$scale))               control$scale               else 1

        scaleVec            <- rep(scaleOriginal, nTheta)
        timesRan            <- 0
        timesAcceptedVec    <- rep(nTheta, 0)
        timesAdapted        <- 0
        optimalAR           <- 0.44

        z                   <- array(0, c(d, d))
        diag(z)             <- 1
        partialSums         <- array(0, c(d, d))  # 1-x_{21}^2, 1-x_{31}^2, 1-x_{31}^2-x_{32}^2, ...
        partialSums[ , 1]   <- 1

        ## Temporary vectors for current row calculations:
        partialSumsProp     <- numeric(d)    
        currentValue        <- numeric(d)
        propValue           <- numeric(d)
        
        ## checks
        dist <- model$getDistribution(target)
        if(dist != 'dlkj_corr_cholesky') stop('RW_lkj_corr_cholesky sampler can only be used with the dlkj_corr_cholesky distribution.')
        if(d < 2)                        stop('RW_lkj_corr_cholesky sampler requires target node dimension to be at least 2x2.')
        if(adaptFactorExponent < 0)      stop('Cannot use RW_lkj_corr_cholesky sampler with adaptFactorExponent control parameter less than 0.')
        if(any(scaleVec < 0))            stop('Cannot use RW_lkj_corr_cholesky sampler with scale control parameter less than 0.')

    },
    run = function() {
        ## FIXME: work with U not L in all cases without t()
        ## calculate transformed values (in unconstrained space) and partial sums in each column
        transform(model[[target]])  # compute z and partialSums
 
        ## Individual univariate RW on the nTheta elements:
        ## advantage: probably better movement through space than with block update
        ## (plus note only column values of target matrix are recalculated for a given scalar update in a given column)
        ## disadvantage: all downstream dependencies of entire matrix will be recalculated.
        cnt <- 0
        for(i in 2:d) {
            currentValue <<- model[[target]][, i]
            partialSumsProp <<- partialSums[i, ]
            for(j in 1:(i-1)) {
                cnt <- cnt + 1
                ## RW on unconstrained y
                yCurrent <- atanh(z[i, j])
                yProp <- rnorm(1, yCurrent, scaleVec[cnt])
                zProp <- tanh(yProp)
                propValue[j] <<- zProp * sqrt(partialSumsProp[j])
                ## Update remainder of row so that length of row is 1
                for(jprime in (j+1):i) {
                    partialSumsProp[jprime] <<- partialSumsProp[jprime-1] - propValue[i, jprime-1]^2
                    propValue[jprime] <<- z[i, jprime] * sqrt(partialSumsProp[jprime])
                }
                model[[target]][j:i, i] <<- propValue[j:i] 
                logMHR <- calculateDiff(model, calcNodesNoSelf) + calculateDiff(model, target)
                ## Adjust MHR to account for non-symmetric proposal by adjusting prior to transformed scale.
                logMHR <- logMHR + 2*(log(cosh(yCurrent)) - log(cosh(yProp)))
                if(j < i-1) {
                    for(jprime in (j+1):(i-1)) 
                        logMHR <- logMHR + 0.5*(log(partialSumsProp[j]) - log(partialSums[i, j]))
                    ## logMHR <- logMHR + 0.5*sum(log(partialSumsProp[(j+1):(i-1)]) - log(partialSums[i, (j+1):(i-1)]))
                }
                jump <- decide(logMHR)
                ## Avoid copying entire target matrix as we are modifying one column at a time.
                if(jump) {
                    timesAcceptedVec[cnt] <<- timesAcceptedVec[cnt] + 1
                    partialSums[i, (j+1):i] <<- partialSumsProp[(j+1):i]
                    currentValue[j:i] <<- propValue[j:i]
                }
                else {
                    nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelf, logProb = TRUE)
                    model[[target]][j:i, i] <<- currentValue[j:i]
                    model$calculate(target)   ## Update target logProb since not part of calcNodesNoSelf.
                    partialSumsProp[j+1] <<- partialSums[i, j+1]
                }
            }
        }
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
            z[2:d, 1] <<- x[1, 2:d]
            for(i in 3:d) {
                for(j in 2:(i-1)) {
                    partialSums[i, j] <<- partialSums[i, j-1] - x[j-1, i]^2
                    z[i, j] <<- x[j, i] / sqrt(partialSums[i, j])
                }
            }
        },
        reset = function() {
            scaleVec         <<- numeric(nTheta, scaleOriginal)
            timesRan         <<- 0
            timesAdapted     <<- 0
            timesAcceptedVec <<- numeric(nTheta, 0)
        }
    ), where = getLoadingNamespace()
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
    ), where = getLoadingNamespace()
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
    ), where = getLoadingNamespace()
)


## RW sampler for non-conjugate scalar components of dcar_normal() and dcar_proper() nodes
CAR_scalar_RW <- nimbleFunction(
    name = 'CAR_scalar_RW',
    contains = sampler_BASE,
    setup = function(model, mvSaved, targetScalar, neighborNodes, neighborWeights, Mi, control, proper) {
        ## control list extraction
        adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
        adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200
        scale         <- if(!is.null(control$scale))         control$scale         else 1
        ## node list generation
        depNodes <- model$getDependencies(targetScalar, self = FALSE)
        copyNodes <- c(targetScalar, depNodes)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        gamma1        <- 0
        optimalAR     <- 0.44
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
        } else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodes, logProb = TRUE)
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
    ), where = getLoadingNamespace()
)


#################################################################################################
### CAR_normal sampler for intrinsic conitionally autoregressive (dcar_normal) distributions  ###
#################################################################################################


#' @rdname samplers
#' @export
sampler_CAR_normal <- nimbleFunction(
    name = 'sampler_CAR_normal',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        useConjugacy  <- if(!is.null(control$carUseConjugacy)) control$carUseConjugacy else TRUE
        ## node list generation
        target <- model$expandNodeNames(target)
        targetScalarComponents <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
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
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list(
        reset = function() {
            for(iSF in seq_along(componentSamplerFunctions))
                componentSamplerFunctions[[iSF]]$reset()
        }
    ), where = getLoadingNamespace()
)


##############################################################################################
### CAR_proper sampler for proper conitionally autoregressive (dcar_proper) distributions  ###
##############################################################################################


#' @rdname samplers
#' @export
sampler_CAR_proper <- nimbleFunction(
    name = 'sampler_CAR_proper',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        useConjugacy  <- if(!is.null(control$carUseConjugacy)) control$carUseConjugacy else TRUE
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
    ), where = getLoadingNamespace()
)


#' MCMC Sampling Algorithms
#'
#' Details of the MCMC sampling algorithms provided with the NIMBLE MCMC engine
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' 
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target node(s) on which the sampler will be used
#' @param control named list that controls the precise behavior of the sampler, with elements specific to \code{samplertype}.  The default values for control list are specified in the setup code of each sampling algorithm.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the \code{control} argument) appear below.
#'
#' @section \code{sampler_base}: base class for new samplers
#'
#' When you write a new sampler for use in a NIMBLE MCMC (see User Manual), you must include \code{contains = sampler_BASE}.
#'
#' @section binary sampler:
#'
#' The binary sampler performs Gibbs sampling for binary-valued (discrete 0/1) nodes.  This can only be used for nodes following either a \code{dbern(p)} or \code{dbinom(p, size=1)} distribution.
#'
#' The binary sampler accepts no control list arguments.
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
#' }
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
#' The ess sampler performs elliptical slice sampling of a single node, which must follow a multivariate normal distribution (Murray, 2010).  The algorithm is an extension of slice sampling (Neal, 2003), generalized to the multivariate normal context.  An auxilliary variable is used to identify points on an ellipse (which passes through the current node value) as candidate samples, which are accepted contingent upon a likelihood evaluation at that point.  This algorithm requires no tuning parameters and therefore no period of adaptation, and may result in very efficient sampling from multivariate Gaussian distributions.
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
#' \item sliceAdaptFactorInterval.  The interval on which to perform factor adaptation. (default = 1000)
#' \item sliceAdaptWidthMaxIter.  The maximum number of iterations for which to adapt the widths for a given set of factors. (default = 512)
#' \item sliceAdaptWidthTolerance. The tolerance for when widths no longer need to adapt, between 0 and 0.5. (default = 0.1)
#' \item sliceMaxSteps.  The maximum number of expansions which may occur during the 'stepping out' procedure. (default = 100)
#' \item maxContractions. The maximum number of contractions of the interval that may occur during sampling (this prevents infinite looping in unusual situations). (default = 100)
#' \item maxContractionsWarning. A logical argument specifying whether to warn when the maximum number of contractions is reached. (default = TRUE)
#' }
#'
#' @section crossLevel sampler:
#'
#' This sampler is constructed to perform simultaneous updates across two levels of stochastic dependence in the model structure.  This is possible when all stochastic descendents of node(s) at one level have conjugate relationships with their own stochastic descendents.  In this situation, a Metropolis-Hastings algorithm may be used, in which a multivariate normal proposal distribution is used for the higher-level nodes, and the corresponding proposals for the lower-level nodes undergo Gibbs (conjugate) sampling.  The joint proposal is either accepted or rejected for all nodes involved based upon the Metropolis-Hastings ratio.
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
#' Sometimes it is useful to control the log likelihood calculations used for an MCMC updater instead of simply using the model.  For example, one could use a sampler with a log likelihood that analytically (or numerically) integrates over latent model nodes.  Or one could use a sampler with a log likelihood that comes from a stochastic approximation such as a particle filter, allowing composition of a particle MCMC (PMCMC) algorithm (Andrieu et al., 2010) (but see samplers listed below for NIMBLE's direct implementation of PMCMC).  The \code{RW_llFunctionBlock} sampler handles this by using a Metropolis-Hastings algorithm with a multivariate normal proposal distribution and a user-provided log-likelihood function.  To allow compiled execution, the log-likelihood function must be provided as a specialized instance of a nimbleFunction.  The log-likelihood function may use the same model as the MCMC as a setup argument, but if so the state of the model should be unchanged during execution of the function (or you must understand the implications otherwise).
#'
#' The RW_llFunctionBlock sampler accepts the following control list elements:
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
#' @section RW_PF sampler:
#'
#' The particle filter sampler allows the user to perform particle MCMC (PMCMC) (Andrieu et al., 2010), primarily for state-space or hidden Markov models of time-series data. This method uses Metropolis-Hastings samplers for top-level parameters but uses the likelihood approximation of a particle filter (sequential Monte Carlo) to integrate over latent nodes in the time-series.  The \code{RW_PF} sampler uses an adaptive Metropolis-Hastings algorithm with a univariate normal proposal distribution for a scalar parameter.  Note that samples of the latent states can be retained as well, but the top-level parameter being sampled must be a scalar.   A bootstrap, auxiliary, or user defined particle filter can be used to integrate over latent states.
#'
#' For more information about user-defined samplers within a PMCMC sampler, see the NIMBLE User Manual.
#'
#' The \code{RW_PF} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (proposal standard deviation) throughout the course of MCMC execution to achieve a theoretically desirable acceptance rate. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW sampler will perform its adaptation procedure.  This updates the scale variable, based upon the sampler's achieved acceptance rate over the past adaptInterval iterations. (default = 200)
#' \item scale. The initial value of the normal proposal standard deviation.  If \code{adaptive = FALSE}, scale will never change. (default = 1)
#' \item pfNparticles.  The number of particles to use in the approximation to the log likelihood of the data (default = 1000).
#' \item latents.  Character vector specifying the nodes that are latent states over which the particle filter will operate to approximate the log-likelihood function.
#' \item pfType.  Character argument specifying the type of particle filter that should be used for likelihood approximation.  Choose from \code{"bootstrap"} and \code{"auxiliary"}.  Defaults to \code{"bootstrap"}.
#' \item pfControl.  A control list that is passed to the particle filter function.  For details on control lists for bootstrap or auxiliary particle filters, see \code{\link{buildBootstrapFilter}} or \code{\link{buildAuxiliaryFilter}} respectively.  Additionally, this can be used to pass custom arguments into a user-defined particle filter.
#' \item pfOptimizeNparticles.  A logical argument, specifying whether to use an experimental procedure to automatically determine the optimal number of particles to use, based on Pitt and Shephard (2011).  This will override any value of \code{pfNparticles} specified above.
#' \item pf.  A user-defined particle filter object, if a bootstrap or auxiliary particle filter is not adequate.
#' }
#' 
#' @section RW_PF_block sampler:
#'
#' The particle filter block sampler allows the user to perform particle MCMC (PMCMC) (Andrieu et al., 2010) for multiple parameters jointly, primarily for state-space or hidden Markov models of time-series data.  This method uses Metropolis-Hastings block samplers for top-level parameters but uses the likelihood approximation of a particle filter (sequential Monte Carlo) to integrate over latent nodes in the time-series.  The \code{RW_PF} sampler uses an adaptive Metropolis-Hastings algorithm with a multivariate normal proposal distribution.  Note that samples of the latent states can be retained as well, but the top-level parameter being sampled must be a scalar.   A bootstrap, auxiliary, or user defined particle filter can be used to integrate over latent states.
#'
#' For more information about user-defined samplers within a PMCMC sampler, see the NIMBLE User Manual.
#' 
#' The \code{RW_PF_block} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the proposal covariance throughout the course of MCMC execution. (default = TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaptation should be done only for \code{scale} (TRUE) or also for \code{provCov} (FALSE).  This argument is only relevant when \code{adaptive = TRUE}.  When \code{adaptScaleOnly = FALSE}, both \code{scale} and \code{propCov} undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When \code{adaptScaleOnly = TRUE}, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the scalar multiplier for \code{propCov}.  If \code{adaptive = FALSE}, \code{scale} will never change. (default = 1)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the \code{'identity'}, in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default is \code{'identity'})
#' \item pfNparticles.  The number of particles to use in the approximation to the log likelihood of the data (default = 1000).
#' \item latents.  Character vector specifying the nodes that are latent states over which the particle filter will operate to approximate the log-likelihood function.
#' \item pfType.  Character argument specifying the type of particle filter that should be used for likelihood approximation.  Choose from \code{"bootstrap"} and \code{"auxiliary"}.  Defaults to \code{"bootstrap"}.
#' \item pfControl.  A control list that is passed to the particle filter function.  For details on control lists for bootstrap or auxiliary particle filters, see \code{\link{buildBootstrapFilter}} or \code{\link{buildAuxiliaryFilter}} respectively.  Additionally, this can be used to pass custom arguments into a user defined particle filter.
#' \item pfOptimizeNparticles.  A logical argument, specifying whether to automatically determine the optimal number of particles to use, based on Pitt and Shephard (2011).  This will override any value of \code{pfNparticles} specified above.
#' \item pf.  A user-defined particle filter object, if a bootstrap or auxiliary particle filter is not adequate.
#' }
#'
#'
#' @section RW_multinomial sampler:
#'
#' This sampler is designed for sampling multinomial target distributions.  The sampler performs a series of Metropolis-Hastings steps between pairs of groups.  Proposals are generated via a draw from a binomial distribution, whereafter the proposed number density is moved from one group to another group.  The acceptance or rejection of these proposals follows a standard Metropolis-Hastings procedure.  Probabilities for the random binomial proposals are adapted to a target acceptance rate of 0.5.
#'
#' The \code{RW_multinomial} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive.  A logical argument, specifying whether the sampler should adapt the binomial proposal probabilities throughout the course of MCMC execution. (default = TRUE)
#' \item adaptInterval.  The interval on which to perform adaptation.  A minimum value of 100 is required. (default = 200)
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
#' \item scale. The initial value of the proposal standard deviation (on the log scale) for each component of the reparameterized Dirichlet distribution.  If adaptive = FALSE, the proposal standard deviations will never change. (default = 1)
#' }
#'
#' @section RW_wishart sampler:
#'
#' This sampler is designed for sampling non-conjugate Wishart and inverse-Wishart distributions.  More generally, it can update any symmetric positive-definite matrix (for example, scaled covaraiance or precision matrices).  The sampler performs block Metropolis-Hastings updates following a transformation to an unconstrained scale (Cholesky factorization of the original matrix, then taking the log of the main diagonal elements.
#'
#' The \code{RW_wishart} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale and proposal covariance for the multivariate normal Metropolis-Hasting proposals, to achieve a theoretically desirable acceptance rate for each. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the sampler will perform its adaptation procedure.  (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the scalar multiplier for the multivariate normal Metropolis-Hastings proposal covariance.  If adaptive = FALSE, scale will never change. (default = 1)
#' }
#'
#' @section RW_lkj_corr_cholesky sampler:
#'
#' This sampler is designed for sampling non-conjugate LKJ correlation Cholesky factor distributions. The sampler performs individual Metropolis-Hastings updates following a transformation to an unconstrained scale (using the signed stickbreaking approach documented in Section 10.12 of the Stan Language Reference Manual, version 2.22).
#'
#' The \code{RW_lkj_corr_cholesky} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scales of the univariate normal Metropolis-Hasting proposals, to achieve a theoretically desirable acceptance rate for each. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the sampler will perform its adaptation procedure.  (default = 200)
#' \item adaptFactorExponent. Exponent controling the rate of decay of the scale adaptation factor.  See Shaby and Wells, 2011, for details. (default = 0.8)
#' \item scale. The initial value of the scalar multiplier for the multivariate normal Metropolis-Hastings proposal covariance.  If adaptive = FALSE, scale will never change. (default = 1)
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
#' @section CRP_concentration sampler:
#' 
#' The CRP_concentration sampler is designed for Bayesian nonparametric mixture modeling. It is exclusively assigned to the concentration parameter of the Dirichlet process when the model is specified using theChinese Restaurant Process distribution, \code{dCRP}. This sampler is assigned by default by NIMBLE's default MCMC configuration and is and can only be used when the prior for the concentration is a gamma distribution. The assigned sampler is an augmented beta-gamma sampler as discussed in Section 6 in Escobar and West (1995).
#' 
#' @section posterior_predictive sampler:
#'
#' The posterior_predictive sampler is only appropriate for use on terminal stochastic nodes.  Note that such nodes play no role in inference but have often been included in BUGS models to accomplish posterior predictive checks.  NIMBLE allows posterior predictive values to be simulated independently of running MCMC, for example by writing a nimbleFunction to do so.  This means that in many cases where terminal stochastic nodes have been included in BUGS models, they are not needed when using NIMBLE.
#'
#' The posterior_predictive sampler functions by calling the simulate() method of relevant node, then updating model probabilities and deterministic dependent nodes.  The application of a posterior_predictive sampler to any non-terminal node will result in invalid posterior inferences.  The posterior_predictive sampler will automatically be assigned to all terminal, non-data stochastic nodes in a model by the default MCMC configuration, so it is uncommon to manually assign this sampler.
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
#' @aliases sampler posterior_predictive RW RW_block RW_multinomial RW_dirichlet RW_wishart RW_lkj_corr_cholesky RW_llFunction slice AF_slice crossLevel RW_llFunction_block RW_PF RW_PF_block sampler_posterior_predictive sampler_RW sampler_RW_block sampler_RW_multinomial sampler_RW_dirichlet sampler_RW_wishart sampler_RW_lkj_corr_cholesky sampler_RW_llFunction sampler_slice sampler_AF_slice sampler_crossLevel sampler_RW_llFunction_block sampler_RW_PF sampler_RW_PF_block CRP CRP_concentration DPmeasure RJ_fixed_prior RJ_indicator RJ_toggled
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
#' # mcmcConf$addSampler(target = 'x[1:5]', type = 'RW_multinomial')   ## x[1:5] ~ dmulti()
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
#' Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., and Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. \emph{The Journal of Chemical Physics}, 21(6), 1087-1092.
#'
#' Neal, Radford M. (2003). Slice Sampling. \emph{The Annals of Statistics}, 31(3), 705-741.
#'
#' Murray, I., Prescott Adams, R., and MacKay, D. J. C. (2010). Elliptical Slice Sampling. \emph{arXiv e-prints}, arXiv:1001.0175.
#'
#' Pitt, M.K. and Shephard, N. (1999). Filtering via simulation: Auxiliary particle filters. \emph{Journal of the American Statistical Association} 94(446), 590-599.
#'
#' Roberts, G. O. and S. K. Sahu (1997). Updating Schemes, Correlation Structure, Blocking and Parameterization for the Gibbs Sampler. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 59(2), 291-317.
#'
#' Shaby, B. and M. Wells (2011). \emph{Exploring an Adaptive Metropolis Algorithm}. 2011-14. Department of Statistics, Duke University.
#'
#' Stan Development Team (2020). \emph{Stan Language Reference Manual, Version 2.22, Section 10.12}.
#'
#' Tibbits, M. M.,  Groendyke, C.,  Haran, M., and Liechty, J. C. (2014).  Automated Factor Slice Sampling.  \emph{Journal of Computational and Graphical Statistics}, 23(2), 543-563.
#' 
#' Escobar, M. D., and West, M. (1995). Bayesian density estimation and inference using mixtures. \emph{Journal of the American Statistical Association}, 90(430), 577-588.
#' 
#' Neal, R. M. (2000). Markov chain sampling methods for Dirichlet process mixture models. \emph{Journal of Computational and Graphical Statistics}, 9(2), 249-265.
#' 
NULL



