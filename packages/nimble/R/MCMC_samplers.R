

####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################

sampler_BASE <- nimbleFunctionVirtual(
    methods = list(
        reset = function() { }
    )
)



####################################################################
### end sampler for trailing stochastic (predictive) nodes #########
####################################################################

sampler_end <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  node list generation  ###
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
### scalar RW sampler with normal proposal distribution ############
####################################################################
sampler_RW <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)     stop('more than one target; cannot use RW sampler, try RW_block sampler')
        calcNodes  <- model$getDependencies(target)
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
	# variables previously inside of nested functions:
        optimalAR <- 0.44
        gamma1    <- 0
    ##    nodeID <- which(model$modelDef$maps$graphID_2_nodeName == target)
    },
    
    run = function() {
    ##    print('entering sampler for nodeID ',nodeID,'\n')
        propValue <- rnorm(1, mean = model[[target]], sd = scale)
     	model[[target]] <<- propValue
        logMHR <- calculateDiff(model, calcNodes)

        jump <- decide(logMHR)
        if(jump) {
      ##      print('accepting')
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        else {
       ##     print('rejecting')
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
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
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
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
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)

sampler_RW_log <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)     stop('more than one target; cannot use RW sampler, try RW_block sampler')
        calcNodes  <- model$getDependencies(target)
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
	# variables previously inside of nested functions:
        optimalAR <- 0.44
        gamma1    <- 0
    },
    
    run = function() {
        logCurrentValue <- log(model[[target]])
        logPropValue <- rnorm(1, mean = logCurrentValue, sd = scale)
        propValue <- exp(logPropValue)
        propRatio <- logPropValue - logCurrentValue
     	model[[target]] <<- propValue
        logMHR <- calculateDiff(model, calcNodes) + propRatio

        jump <- decide(logMHR)
        if(jump)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },
    
    methods = list(
    
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
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
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)


sampler_RW_nonDiff <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)     stop('more than one target; cannot use RW sampler, try RW_block sampler')
        calcNodes  <- model$getDependencies(target)
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
	# variables previously inside of nested functions:
        optimalAR <- 0.44
        gamma1    <- 0
    },
    
    run = function() {
        modelLP0 <- getLogProb(model, calcNodes)
        propValue <- rnorm(1, mean = model[[target]], sd = scale)
     	model[[target]] <<- propValue
        modelLP1 <- calculate(model, calcNodes)
        logMHR <- modelLP1 - modelLP0
        jump <- decide(logMHR)
        if(jump)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },
    
    methods = list(
    
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
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
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)



########################################################################
### block RW sampler with multi-variate normal proposal distribution ###
########################################################################
sampler_RW_block <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
        d <- length(targetAsScalar)
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        statSums  <- matrix(0, nrow=1, ncol=d)   # sums of each node, stored as a row-matrix
        statProds <- matrix(0, nrow=d, ncol=d)   # sums of pairwise products of nodes
        ###  nested function and function list definitions  ###
        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
    },
    
    run = function() {
        propValueVector <- generateProposalVector()
        lpMHR <- my_setAndCalculateDiff$run(propValueVector)
        jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## will use lpMHR - 0
        if(adaptive)     adaptiveProcedure(jump)
    },
    
    methods = list(
        
        generateProposalVector = function() {
            propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov * scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly) {
                declare(newValues, double(1, d))

                newValues <- values(model, target)
                statSums  <<- statSums + asRow(newValues)
                statProds <<- statProds + asCol(newValues) %*% asRow(newValues)
            }
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$gamma1
                    empirCov <- (statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                    statSums  <<- statSums  * 0
                    statProds <<- statProds * 0      ##  setAll(statProds, 0)    ## setAll() doesn't work in R, and doesn't work for vectors (only works for dim=2 objects)
                }
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
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            statSums  <<- statSums  * 0
            statProds <<- statProds * 0
            my_calcAdaptationFactor$reset()
        }
    ),  where = getLoadingNamespace()
)

sampler_RW_block_noDiff <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
        d <- length(targetAsScalar)
        if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
        if(class(propCov) != 'matrix')        stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric')  stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov) == d))           stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))             stop('propCov matrix must be symmetric')
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        statSums  <- matrix(0, nrow=1, ncol=d)   # sums of each node, stored as a row-matrix
        statProds <- matrix(0, nrow=d, ncol=d)   # sums of pairwise products of nodes
        ###  nested function and function list definitions  ###
        my_setAndCalculate <- setAndCalculate(model, target)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
    },
    
    run = function() {
        modelLP0 <- getLogProb(model, calcNodes)
        propValueVector <- generateProposalVector()
        modelLP1 <- my_setAndCalculate$run(propValueVector)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
        if(adaptive)     adaptiveProcedure(jump)
    },
    
    methods = list(
        
        generateProposalVector = function() {
            propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov * scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly) {
                declare(newValues, double(1, d))

                newValues <- values(model, target)
                statSums  <<- statSums + asRow(newValues)
                statProds <<- statProds + asCol(newValues) %*% asRow(newValues)
            }
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)
                setSize(acceptanceRateHistory, timesAdapted)
                scaleHistory[timesAdapted] <<- scale
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$gamma1
                    empirCov <- (statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                    statSums  <<- statSums  * 0
                    statProds <<- statProds * 0      ##  setAll(statProds, 0)    ## setAll() doesn't work in R, and doesn't work for vectors (only works for dim=2 objects)
                }
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
            scaleHistory          <<- scaleHistory          * 0
            acceptanceRateHistory <<- acceptanceRateHistory * 0
            statSums  <<- statSums  * 0
            statProds <<- statProds * 0
            my_calcAdaptationFactor$reset()
        }
    ),  where = getLoadingNamespace()
)



#############################################################################
### RW_llFunction, does a RW, but using a generic log-likelihood function ###
#############################################################################

sampler_RW_llFunction <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive       <- control$adaptive
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        llFunction     <- control$llFunction
        includesTarget <- control$includesTarget
        ###  node list generation  ###
        calcNodes <- model$getDependencies(target)
        ###  nested function and function list definitions  ###
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
    ),  where = getLoadingNamespace()
)



####################################################################
### slice sampler (discrete or continuous) #########################
####################################################################

sampler_slice <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        width         <- control$sliceWidth
        maxSteps      <- control$sliceMaxSteps
        ###  node list generation  ###
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        if(length(targetAsScalar) > 1)     stop(paste0('more than one target: ', target, '; cannot use slice sampler'))
        calcNodes <- model$getDependencies(target)
        ###  numeric value generation  ###
        widthOriginal <- width
        timesRan      <- 0
        timesAdapted  <- 0
        sumJumps      <- 0
        discrete      <- model$isDiscrete(targetAsScalar)
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
        while(is.nan(lp) | lp < u) {   # MUST be is.nan()
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
    ),  where = getLoadingNamespace()
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
    },
    where = getLoadingNamespace()
)

sampler_crossLevel <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ###  control list extraction  ###
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        ###  node list generation  ###
        target       <- model$expandNodeNames(target)
        lowNodes     <- model$getDependencies(target, self = FALSE, stochOnly = TRUE, includeData = FALSE)
        lowCalcNodes <- model$getDependencies(lowNodes)
        calcNodes    <- model$getDependencies(c(target, lowNodes))
        ###  nested function and function list definitions  ###
        mvInternal <- modelValues(model)
        RWblockControl <- list(adaptive = adaptive, adaptScaleOnly = adaptScaleOnly, adaptInterval = adaptInterval, scale = scale, propCov = propCov)
        topRWblockSamplerFunction <- sampler_RW_block(model, mvInternal, target, RWblockControl)
        lowConjugateSamplerFunctions <- nimbleFunctionList(sampler_BASE)
        lowConjugateGetLogDensityFunctions <- nimbleFunctionList(getPosteriorDensityFromConjSampler_virtual)
        for(iLN in seq_along(lowNodes)) {
            lowNode <- lowNodes[iLN]
            conjugacyResult <- model$checkConjugacy2(lowNode)[[lowNode]]
            if(is.null(conjugacyResult))     stop('non-conjugate lowNode \'', lowNode, '\' in crossLevel updater')
            samplerType <- conjugacyResult$type
            samplerFunction <- eval(as.name(paste0('sampler_',samplerType)))
            conjugateSamplerSpec <- samplerSpec(name = samplerType, samplerFunction = samplerFunction, target = lowNode, control = conjugacyResult$control, model = model)
            lowConjugateSamplerFunctions[[iLN]] <- conjugateSamplerSpec$buildSampler(model, mvInternal)
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
        else{
	        for(iSF in seq_along(lowConjugateSamplerFunctions))        { lowConjugateSamplerFunctions[[iSF]]$run() }
	        modelLP1 <- calculate(model, calcNodes)
	        propLP1 <- 0
	        for(iSF in seq_along(lowConjugateGetLogDensityFunctions))  { propLP1 <- propLP1 + lowConjugateGetLogDensityFunctions[[iSF]]$run() }
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
    ),  where = getLoadingNamespace()
)



#' MCMC Sampling Algorithms
#'
#' Details of the control argument for customizing sampling algorithms provided with the NIMBLE MCMC engine.
#'
#' The precise behaviour NIMBLE's MCMC sampling algorithms may be customized using the control argument provided to addSampler().  The usage syntax is:
#'
#' mcmcspec$addSampler(target = TARGETNODE, type = SAMPLERTYPE, control = CONTROLLIST)
#'
#' where CONTROLLIST is a named list, with elements specific to SAMPLERTYPE.  The default values for control list elements are determined by the NIMBLE system option \'MCMCcontrolDefaultList\'.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the control argument) appear below.
#'
#' ----- end sampler -----
#'
#' The end sampler is only appropriate for use on terminal stochastic nodes.  Note that such nodes play no role in inference but have often been included in BUGS models to accomplish posterior predictive checks.  NIMBLE allows posterior predictive values to be simulated independently of running MCMC, for example by writing a nimbleFunction to do so.  This means that in many cases where terminal stochastic nodes have been included in BUGS models, they are not needed when using NIMBLE.
#'
#' The end sampler functions by calling the simulate() method of relevant node, then updating model probabilities and deterministic dependent nodes.  The application of an end sampler to any non-terminal node will result in invalid posterior inferences.  The end sampler will automatically be assigned to all terminal, non-data stochastic nodes in a model by the default MCMC configuration, so it is uncommon to manually assign this sampler.
#'
#' The end sampler accepts no control list arguments.
#' 
#' ----- RW sampler -----
#'
#' The RW sampler executes adaptive Metropolis-Hastings sampling with a normal proposal distribution (Metropolis, 1953), implementing the adaptation routine given in Shaby and Wells, 2011.  This sampler can be applied to any scalar continuous-valued stochastic node.  
#'
#' The RW sampler accepts the following control list elements:
#' \begin{itemize}
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (proposal standard deviation) throughout the course of MCMC execution to achieve a theoretically desirable acceptance rate. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW sampler will perform its adaptation procedure.  This updates the scale variable, based upon the sampler's achieved acceptance rate over the past adaptInterval iterations. (default = 200)
#' \item scale. The initial value of the normal proposal standard deviation.  If adaptive = FALSE, scale will never change. (default = 1)
#' \end{itemize}
#' 
#' ----- RW_block sampler -----
#'
#' The RW\_block sampler performs a simultaneous update of one or more model nodes, using an adaptive Metropolis-Hastings algorithm with a multivariate normal proposal distribution (Roberts and Sahu, 1997), implementing the adaptation routine given in Shaby and Wells, 2011.  This sampler may be applied to any set of continuous-valued model nodes, to any single continuous-valued multivariate model node, or to any combination thereof.
#'
#' The RW\_block sampler accepts the following control list elements:
#' \begin{itemize}
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (a coefficient for the entire proposal covariance matrix) and propCov (the multivariate normal proposal covariance matrix) throughout the course of MCMC execution.  If only the scale should undergo adaptation, this argument should be specified as TRUE. (default = TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for scale (TRUE) or also for provCov (FALSE).  This argument is only relevant when adaptive = TRUE.  When adaptScaleOnly = FALSE, both scale and propCov undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When adaptScaleOnly = FALSE, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW\_block sampler will perform its adaptation procedure, based on the past adaptInterval iterations. (default = 200)
#' \item scale. The initial value of the scalar multiplier for propCov.  If adaptive = FALSE, scale will never change. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string \'identity\', in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default = \'identity\')
#' \end{itemize}
#'
#' ----- RW_llFunction sampler -----
#'
#' Sometimes it is useful to control the log likelihood calculations used for an MCMC updater instead of simply using the model.  For example, one could use a sampler with a log likelihood that analytically (or numerically) integrates over latent model nodes.  Or one could use a sampler with a log likelihood that comes from a stochastic approximation such as a particle filter, allowing composition of a particle MCMC (PMCMC) algorithm (Andrieu, 2010).  The RW\_llFunction sampler handles this by using a Metropolis-Hastings algorithm with a normal proposal distribution and a user-provided log-likelihood function.  To allow compiled execution, the log-likelihood function must be provided as a specialized instance of a nimbleFunction.  The log-likelihood function may use the same model as the MCMC as a setup argument, but if so the state of the model should be unchanged during execution of the function (or you must understand the implications otherwise).
#'
#' The RW\_llFunction sampler accepts the following control list elements:
#' \begin{itemize}
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the \scale (proposal standard deviation) throughout the course of MCMC execution. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the normal proposal standard deviation. (default = 1)
#' \item llFunction. A specialized nimbleFunction that accepts no arguments and returns a scalar double number.  The return value must be the total log-likelihood of all stochastic dependents of the target nodes -- and, if includesTarget = TRUE, of the target node(s) themselves --  or whatever surrogate is being used for the total log-likelihood.  This is a required element with no default.
#' \item includesTarget. Logical variable indicating whether the return value of llFunction includes the log-likelihood associated with target.  This is a required element with no default.
#' \end{itemize}
#' 
#' ----- slice sampler -----
#'
#' The slice sampler performs slice sampling of the scalar node to which it is applied (Neal, 2003).  This sampler can operate on either continuous-valued or discrete-valued scalar nodes.  The slice sampler performs a \'stepping out\' procedure, in which the slice is iteratively expanded to the left or right by an amount sliceWidth.  This sampler is optionally adaptive, governed by a control list element, whereby the value of sliceWidth is adapted towards the observed absolute difference between successive samples.
#'
#' The slice sampler accepts the following control list elements:
#' \begin{itemize}
#' \item adaptive. A logical argument, specifying whether the sampler will adapt the value of sliceWidth throughout the course of MCMC execution. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item sliceWidth. The initial value of the width of each slice, and also the width of the expansion during the iterative \'stepping out\' procedure. (default = 1)
#' \item sliceMaxSteps. The maximum number of expansions which may occur during the \'stepping out\' procedure. (default = 100)
#' \end{itemize}
#'
#' ----- crossLevel sampler -----
#'
#' This sampler is constructed to perform simultaneous updates across two levels of stochastic dependence in the model structure.  This is possible when all stochastic descendents of node(s) at one level have conjugate relationships with their own stochastic descendents.  In this situation, a Metropolis-Hastings algorithm may be used, in which a multivariate normal proposal distribution is used for the higher-level nodes, and the corresponding proposals for the lower-level nodes undergo Gibbs (conjugate) sampling.  The joint proposal is either accepted or rejected for all nodes involved based upon the Metropolis-Hastings ratio.
#'
#' The requirement that all stochastic descendents of the target nodes must themselves have only conjugate descendents will be checked when the MCMC algorithm is built.  This sampler is useful when there is strong dependence across the levels of a model that causes problems with convergence or mixing.
#' 
#' The crossLevel sampler accepts the following control list elements:
#' \begin{itemize}
#' \item adaptive. Logical argument, specifying whether the multivariate normal proposal distribution for the target nodes should be adaptived. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item scale. The initial value of the scalar multiplier for propCov. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string \'identity\' or any positive definite matrix of the appropriate dimensions. (default = \'identity\')
#' \end{itemize}
#'
#' @name samplers
#'
#' @aliases sampler end RW RW_block RW_llFunction slice crossLevel sampler_end sampler_RW sampler_RW_block sampler_RW_llFunction sampler_slice sampler_crossLevel
#'
#' @seealso configureMCMC addSampler buildMCMC
#'
#' @references Andrieu, C., Doucet, A., and Holenstein, R. (2010). Particle Markov Chain Monte Carlo Methods. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 72(3), 269-342.
#' @references Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., and Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. \emph{The Journal of Chemical Physics}, 21(6), 1087-1092.
#' @references Neal, Radford M. (2003). Slice Sampling. \emph{The Annals of Statistics}, 31(3), 705-741.
#' @references Roberts, G. O. and S. K. Sahu (1997). Updating Schemes, Correlation Structure, Blocking and Parameterization for the Gibbs Sampler. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 59(2), 291-317.
#' @references Shaby, B. and M. Wells (2011). \emph{Exploring an Adaptive Metropolis Algorithm}. 2011-14. Department of Statistics, Duke University.
NULL



