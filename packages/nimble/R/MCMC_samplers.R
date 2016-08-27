


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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        calcNodes  <- model$getDependencies(target)
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes  <- model$getDependencies(target)
        ## checks
        if(length(targetAsScalar) > 1)  stop('cannot use binary sampler on more than one target node')
        if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
    },
    run = function() {
        currentProb <- exp(getLogProb(model, calcNodes))
        model[[target]] <<- 1 - model[[target]]
        otherProb <- exp(calculate(model, calcNodes))
        if(!is.nan(otherProb) & runif(1,0,1) < otherProb/(currentProb+otherProb))
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- control$log
        reflective    <- control$reflective
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes  <- model$getDependencies(target)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        optimalAR     <- 0.44
        gamma1        <- 0
        range <- getDistribution(model$getNodeDistribution(target))$range
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
                       propValue <- currentValue * exp(propLogScale)
                   } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
        if(reflective)   while(propValue < range[1] | propValue > range[2]) {
            if(propValue < range[1]) propValue <- 2*range[1] - propValue
            if(propValue > range[2]) propValue <- 2*range[2] - propValue    }
        model[[target]] <<- propValue
        logMHR <- calculateDiff(model, calcNodes) + propLogScale
        jump <- decide(logMHR)
        if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
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



########################################################################
### block RW sampler with multi-variate normal proposal distribution ###
########################################################################

#' @rdname samplers
#' @export
sampler_RW_block <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
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
        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
        ## checks
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    },
    run = function() {
        propValueVector <- generateProposalVector()
        lpMHR <- my_setAndCalculateDiff$run(propValueVector)
        jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## will use lpMHR - 0
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
                    gamma1 <- my_calcAdaptationFactor$gamma1
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
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        llFunction     <- control$llFunction
        includesTarget <- control$includesTarget
        ## node list generation
        calcNodes <- model$getDependencies(target)
        ## nested function and function list definitions
        mvInternal <- modelValues(model)
        RWControl <- list(adaptive = adaptive, adaptInterval = adaptInterval, scale = scale)
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        width         <- control$sliceWidth
        maxSteps      <- control$sliceMaxSteps
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
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
        while(is.nan(lp) | lp < u) {   # must be is.nan()
            if(x1 < x0) { L <- x1
                      } else      { R <- x1 }
            x1 <- L + runif(1, 0, 1) * (R - L)           # sample uniformly from (L,R) until sample is inside of slice (with shrinkage)
            lp <- setAndCalculateTarget(x1)
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        jumpDist <- abs(x1 - x0)
        if(adaptive)     adaptiveProcedure(jumpDist)
    },
    methods = list(
        setAndCalculateTarget = function(value = double()) {
            if(discrete)     value <- floor(value)
            model[[target]] <<- value
            lp <- calculate(model, calcNodes)
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## node list generation
        target <- model$expandNodeNames(target)
        calcNodes  <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        Pi <- pi
        ## nested function and function list definitions
        ##target_nodeFunctionList <- nimbleFunctionList(node_stoch_dmnorm)
        ##target_nodeFunctionList[[1]] <- model$nodeFunctions[[target]]
        ## checks
        if(length(target) > 1)                              stop('elliptical slice sampler only applies to one target node')
        if(model$getNodeDistribution(target) != 'dmnorm')   stop('elliptical slice sampler only applies to multivariate normal distributions')
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
        while(is.nan(lp) | lp < u) {   # must be is.nan()
            if(theta < 0)   theta_min <- theta   else   theta_max <- theta
            theta <- runif(1, theta_min, theta_max)
            model[[target]] <<- f*cos(theta) + nu*sin(theta) + target_mean
            lp <- calculate(model, calcNodes)
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
)



####################################################################
### crossLevel sampler #############################################
####################################################################

getPosteriorDensityFromConjSampler_virtual <- nimbleFunctionVirtual(
    run = function() { returnType(double()) }
)

getPosteriorDensityFromConjSampler <- nimbleFunction(
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        ## node list generation
        target       <- model$expandNodeNames(target)
        lowNodes     <- model$getDependencies(target, self = FALSE, stochOnly = TRUE, includeData = FALSE)
        lowCalcNodes <- model$getDependencies(lowNodes)
        calcNodes    <- model$getDependencies(c(target, lowNodes))
        ## nested function and function list definitions
        mvInternal <- modelValues(model)
        RWblockControl <- list(adaptive = adaptive, adaptScaleOnly = adaptScaleOnly, adaptInterval = adaptInterval, scale = scale, propCov = propCov)
        topRWblockSamplerFunction <- sampler_RW_block(model, mvInternal, target, RWblockControl)
        lowConjugateSamplerFunctions <- nimbleFunctionList(sampler_BASE)
        lowConjugateGetLogDensityFunctions <- nimbleFunctionList(getPosteriorDensityFromConjSampler_virtual)
        for(iLN in seq_along(lowNodes)) {
            lowNode <- lowNodes[iLN]
            conjugacyResult <- model$checkConjugacy2(lowNode)[[lowNode]]
            if(is.null(conjugacyResult))     stop('non-conjugate lowNode \'', lowNode, '\' in crossLevel updater')
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        llFunction     <- control$llFunction
        includesTarget <- control$includesTarget
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
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
        storeLP0 <- -Inf
        ## checks
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
    },
    run = function() {
        modelLP0 <- storeLP0
        if(!includesTarget)     modelLP0 <- modelLP0 + getLogProb(model, target)
        propValueVector <- generateProposalVector()
        my_setAndCalculate$run(propValueVector)
        modelLP1 <- llFunction$run()+ getLogProb(model, target)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
        if(adaptive)     adaptiveProcedure(jump)
        if(jump) storeLP0 <<- modelLP1
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
                    gamma1 <- my_calcAdaptationFactor$gamma1
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        m              <- control$pfNparticles
        resample       <- control$pfResample
        filterType     <- control$pfType
        lookahead      <- control$pfLookahead
        optimizeM      <- as.integer(control$pfOptimizeNparticles)
        latents        <- control$latents
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
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        prevLL        <- 0
        nVarEsts      <- 0
        itCount       <- 0
        optimalAR     <- 0.44
        gamma1        <- 0
        storeLP0      <- -Inf
        storeLLVar    <- 0
        nVarReps      <- 7    # number of LL estimates to compute to get each LL variance estimate for m optimization
        mBurnIn       <- 15   # number of LL variance estimates to compute before deciding optimal m
        d             <- length(targetAsScalar)
        if(optimizeM)    m <- 3000
        ## nested function and function list definitions
        my_setAndCalculate <- setAndCalculate(model, target)
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        if(filterType == 'auxiliary') {
            my_particleFilter <- buildAuxiliaryFilter(model, latents, control = list(saveAll = TRUE, smoothing = TRUE, lookahead = lookahead))
        }
        else if(filterType == 'bootstrap') {
            my_particleFilter <- buildBootstrapFilter(model, latents, control = list(saveAll = TRUE, smoothing = TRUE))
        }
        else   stop('filter type must be either bootstrap or auxiliary')
        particleMV <- my_particleFilter$mvEWSamples
        ## checks
        if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
        if(length(targetAsScalar) > 1)                      stop('more than one top-level target; cannot use RW_PF sampler, try RW_PF_block sampler')
    },
    run = function() {
        if(resample) {
            modelLP0 <- my_particleFilter$run(m)
            modelLP0 <- modelLP0 + getLogProb(model, target)
        }
        else   modelLP0 <- storeLP0
        propValue <- rnorm(1, mean = model[[target]], sd = scale)
        model[[target]] <<- propValue
        modelLP1 <- my_particleFilter$run(m)
        modelLP1 <- modelLP1 + getLogProb(model, target)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
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
        if(jump & !resample)  storeLP0 <<- modelLP1
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
                if(!resample) {  #reset LL with new m value after burn-in period
                    modelLP0 <- my_particleFilter$run(m)
                    storeLP0 <<- modelLP0 + getLogProb(model, target)
                }
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
            storeLP0 <<- -Inf
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target,  control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        m              <- control$pfNparticles
        resample       <- control$pfResample
        filterType     <- control$pfType
        lookahead      <- control$pfLookahead
        optimizeM      <- as.integer(control$pfOptimizeNparticles)
        latents        <- control$latents
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
        storeLP0    <- -Inf
        storeLLVar  <- 0
        nVarReps <- 7    # number of LL estimates to compute to get each LL variance estimate for m optimization
        mBurnIn  <- 15   # number of LL variance estimates to compute before deciding optimal m
        if(optimizeM)   m <- 3000
        ## nested function and function list definitions
        my_setAndCalculate <- setAndCalculate(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
        if(latentSamp == TRUE) { saveAllVal <- TRUE
                                 smoothingVal <- TRUE
                             } else {
                                 saveAllVal <- FALSE
                                 smoothingVal <- FALSE
                             }
        if(filterType == 'auxiliary') {
            my_particleFilter <- buildAuxiliaryFilter(model, latents, control = list(saveAll = saveAllVal, smoothing = smoothingVal, lookahead = lookahead))
        }
        else if(filterType == 'bootstrap') {
            my_particleFilter <- buildBootstrapFilter(model, latents, control = list(saveAll = saveAllVal, smoothing = smoothingVal))
        }
        else   stop('filter type must be either bootstrap or auxiliary')
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
        if(resample) {
            modelLP0 <- my_particleFilter$run(m)
            modelLP0 <- modelLP0 + getLogProb(model, target)
        }
        else   modelLP0 <- storeLP0
        propValueVector <- generateProposalVector()
        my_setAndCalculate$run(propValueVector)
        modelLP1 <- my_particleFilter$run(m)
        modelLP1 <- modelLP1 + getLogProb(model, target)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
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
        if(jump & !resample)  storeLP0 <<- modelLP1
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
                if(!resample){  #reset LL with new m value after burn-in period
                    modelLP0 <- my_particleFilter$run(m)
                    storeLP0 <<- modelLP0 + getLogProb(model, target)
                }
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
                    gamma1 <- my_calcAdaptationFactor$gamma1
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
            storeLP0 <<- -Inf
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
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
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
        if(model$getNodeDistribution(target) != 'dmulti')   stop('can only use RW_multinomial sampler for multinomial distributions')
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



#' MCMC Sampling Algorithms
#'
#' Details of the MCMC sampling algorithms provided with the NIMBLE MCMC engine
#'
#' @param model (uncompiled) model on which the MCMC is to be run
#' 
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target node(s) on which the sampler will be used
#' @param control named list that controls the precise behavior of the sampler, with elements specific to \code{samplertype}.  The default values for control list elements are determined by the NIMBLE system option \code{'MCMCcontrolDefaultList'}.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the \code{control} argument) appear below.
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
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for scale (TRUE) or also for provCov (FALSE).  This argument is only relevant when adaptive = TRUE.  When adaptScaleOnly = FALSE, both scale and propCov undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When adaptScaleOnly = FALSE, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW_block sampler will perform its adaptation procedure, based on the past adaptInterval iterations. (default = 200)
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
#' }
#'
#' @section ess sampler:
#'
#' The ess sampler performs elliptical slice sampling of a single node, which must follow a multivariate normal distribution (Murray, 2010).  The algorithm is an extension of slice sampling (Neal, 2003), generalized to the multivariate normal context.  An auxilliary variable is used to identify points on an ellipse (which passes through the current node value) as candidate samples, which are accepted contingent upon a likelihood evaluation at that point.  This algorithm requires no tuning parameters and therefore no period of adaptation, and may result in very efficient sampling from multivariate Gaussian distributions.
#'
#' The ess sampler accepts no control list arguments.
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
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for scale (TRUE) or also for provCov (FALSE).  This argument is only relevant when adaptive = TRUE.  When adaptScaleOnly = FALSE, both scale and propCov undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When adaptScaleOnly = FALSE, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the scalar multiplier for propCov.  If adaptive = FALSE, scale will never change. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string 'identity', in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default = 'identity')
#' \item llFunction. A specialized nimbleFunction that accepts no arguments and returns a scalar double number.  The return value must be the total log-likelihood of all stochastic dependents of the target nodes -- and, if includesTarget = TRUE, of the target node(s) themselves --  or whatever surrogate is being used for the total log-likelihood.  This is a required element with no default.
#' \item includesTarget. Logical variable indicating whether the return value of llFunction includes the log-likelihood associated with target.  This is a required element with no default.
#' }
#'
#'
#' @section RW_PF sampler:
#'
#' The particle filter sampler allows the user to perform PMCMC (Andrieu et al., 2010), integrating over latent nodes in the model to sample top-level parameters.  The \code{RW_PF} sampler uses a Metropolis-Hastings algorithm with a univariate normal proposal distribution for a scalar parameter.  Note that latent states can be sampled as well, but the top-level parameter being sampled must be a scalar.   A bootstrap or auxiliary particle filter can be used to integrate over latent states.
#'
#' The \code{RW_PF} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (proposal standard deviation) throughout the course of MCMC execution to achieve a theoretically desirable acceptance rate. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW sampler will perform its adaptation procedure.  This updates the scale variable, based upon the sampler's achieved acceptance rate over the past adaptInterval iterations. (default = 200)
#' \item scale. The initial value of the normal proposal standard deviation.  If \code{adaptive = FALSE}, scale will never change. (default = 1)
#' \item pfNparticles.  The number of particles to use in the approximation to the log likelihood of the data (default = 1000).
#' \item latents.  Character vector specifying the latent model nodes over which the particle filter will stochastically integrate over to estimate the log-likelihood function.
#' \item pfType.  Character argument specifying the type of particle filter that should be used for likelihood approximation.  Choose from \code{"bootstrap"} and \code{"auxiliary"}.  Defaults to \code{"bootstrap"}.
#' \item pfLookahead. Optional character argument specifying the lookahead function for the auxiliary particle filter.  Choose from \code{"simulate"} and \code{"mean"}.  Only applicable if \code{pfType} is set to \code{"auxiliary"}.
#' \item pfResample.  A logical argument, specifying whether to resample log likelihood given current parameters at beginning of each MCMC step, or whether to use log likelihood from previous step.
#' \item pfOptimizeNparticles.  A logical argument, specifying whether to automatically determine the optimal number of particles to use, based on Pitt and Shephard (2011).  This will override any value of \code{pfNparticles} specified above.
#' }
#' 
#' @section RW_PF_block sampler:
#'
#' The particle filter sampler allows the user to perform PMCMC (Andrieu et al., 2010), integrating over latent nodes in the model to sample top-level parameters.  The \code{RW_PF_block} sampler uses a Metropolis-Hastings algorithm with a multivariate normal proposal distribution.  A bootstrap or auxiliary particle filter can be used to integrate over latent states.
#'
#' The \code{RW_PF_block} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the proposal covariance throughout the course of MCMC execution. (default = TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for \code{scale} (TRUE) or also for \code{provCov} (FALSE).  This argument is only relevant when \code{adaptive = TRUE}.  When \code{adaptScaleOnly = FALSE}, both \code{scale} and \code{propCov} undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When \code{adaptScaleOnly = FALSE}, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the scalar multiplier for \code{propCov}.  If \code{adaptive = FALSE}, \code{scale} will never change. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the \code{'identity'}, in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default is \code{'identity'})
#' \item pfNparticles.  The number of particles to use in the approximation to the log likelihood of the data (default = 1000).
#' \item latents.  Character vector specifying the latent model nodes over which the particle filter will stochastically integrate to estimate the log-likelihood function.
#' \item pfResample.  A logical argument, specifying whether to resample log likelihood given current parameters at beginning of each mcmc step, or whether to use log likelihood from previous step.
#' \item pfType.  Character argument specifying the type of particle filter that should be used for likelihood approximation.  Choose from \code{"bootstrap"} and \code{"auxiliary"}.  Defaults to \code{"bootstrap"}.
#' \item pfLookahead. Optional character argument specifying the lookahead function for the auxiliary particle filter.  Choose from \code{"simulate"} and \code{"mean"}.  Only applicable if \code{pfType = "auxiliary"}.
#' \item pfOptimizeNparticles.  A logical argument, specifying whether to automatically determine the optimal number of particles to use, based on Pitt and Shephard (2011).  This will override any value of \code{pfNparticles} specified above.
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
#' @section posterior_predictive sampler:
#'
#' The posterior_predictive sampler is only appropriate for use on terminal stochastic nodes.  Note that such nodes play no role in inference but have often been included in BUGS models to accomplish posterior predictive checks.  NIMBLE allows posterior predictive values to be simulated independently of running MCMC, for example by writing a nimbleFunction to do so.  This means that in many cases where terminal stochastic nodes have been included in BUGS models, they are not needed when using NIMBLE.
#'
#' The posterior_predictive sampler functions by calling the simulate() method of relevant node, then updating model probabilities and deterministic dependent nodes.  The application of a posterior_predictive sampler to any non-terminal node will result in invalid posterior inferences.  The posterior_predictive sampler will automatically be assigned to all terminal, non-data stochastic nodes in a model by the default MCMC configuration, so it is uncommon to manually assign this sampler.
#'
#' The posterior_predictive sampler accepts no control list arguments.
#'
#' @name samplers
#'
#' @aliases sampler posterior_predictive RW RW_block RW_multinomial RW_llFunction slice crossLevel RW_llFunction_block RW_PF RW_PF_block sampler_posterior_predictive sampler_RW sampler_RW_block sampler_RW_multinomial sampler_RW_llFunction sampler_slice sampler_crossLevel sampler_RW_llFunction_block sampler_RW_PF sampler_RW_PF_block
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
NULL



