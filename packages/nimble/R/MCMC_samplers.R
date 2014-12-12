

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
    setup = function(model, mvSaved, control) {
        ###  control list extraction  ###
        targetNode <- control$targetNode
        ###  node list generation  ###
        calcNodes  <- model$getDependencies(targetNode, returnType = 'names')
    },
    run = function() {
        simulate(model, targetNode)
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
    setup = function(model, mvSaved, control) {
        ###  control list extraction  ###
        targetNode    <- control$targetNode
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ###  node list generation  ###
        targetNodeAsScalar <- model$expandNodeNames(targetNode, returnScalarComponents = TRUE)
        if(length(targetNodeAsScalar) > 1)     stop('more than one targetNode; cannot use RW sampler, try RW_block sampler')
        calcNodes  <- model$getDependencies(targetNode)
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
        propValue <- rnorm(1, mean = model[[targetNode]], sd = scale)
     	model[[targetNode]] <<- propValue
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

#sampler_RW <- nimbleFunction(
#    contains = sampler_BASE,
#    setup = function(model, mvSaved, control) {
#        ###  control list extraction  ###
#        targetNode    <- control$targetNode
#        adaptive      <- control$adaptive
#        adaptInterval <- control$adaptInterval
#        scale         <- control$scale
#        ###  node list generation  ###
#        calcNodes  <- model$getDependencies(targetNode, returnType = 'nodeSet')
#        ###  numeric value generation  ###
#        scaleOriginal <- scale
#        timesRan      <- 0
#        timesAccepted <- 0
#        timesAdapted  <- 0
#        scaleHistory          <- c(0, 0)
#        acceptanceRateHistory <- c(0, 0)
#        ###  nested function and function list definitions  ###
#        my_setAndCalculateOne <- setAndCalculateOne(model, targetNode)
#        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
#        my_calcAdaptationFactor <- calcAdaptationFactor(1)
#    },
#    
#    run = function() {
#        modelLP0 <- getLogProb(model, calcNodes)
#        propValue <- generateProposalValue()
#        modelLP1 <- my_setAndCalculateOne(propValue)
#        jump <- my_decideAndJump(modelLP1, modelLP0, 0, 0)
#        if(adaptive)     adaptiveProcedure(jump)
#    },
#    
#    methods = list(
#    
#        generateProposalValue = function() {
#            propValue <- rnorm(1, mean = model[[targetNode]], sd = scale)
#            returnType(double())
#            return(propValue)
#        },
#        
#        adaptiveProcedure = function(jump = logical()) {
#            timesRan <<- timesRan + 1
#            if(jump)     timesAccepted <<- timesAccepted + 1
#            if(timesRan %% adaptInterval == 0) {
#                acceptanceRate <- timesAccepted / timesRan
#                timesAdapted <<- timesAdapted + 1
#                setSize(scaleHistory,          timesAdapted)
#                setSize(acceptanceRateHistory, timesAdapted)
#                scaleHistory[timesAdapted] <<- scale
#                acceptanceRateHistory[timesAdapted] <<- acceptanceRate
#                adaptFactor <- my_calcAdaptationFactor(acceptanceRate)
#                scale <<- scale * adaptFactor
#                timesRan <<- 0
#                timesAccepted <<- 0
#            }
#        },
#        
#        reset = function() {
#            scale <<- scaleOriginal
#            timesRan      <<- 0
#            timesAccepted <<- 0
#            timesAdapted  <<- 0
#            scaleHistory          <<- scaleHistory          * 0
#            acceptanceRateHistory <<- acceptanceRateHistory * 0
#            nfMethod(my_calcAdaptationFactor, 'reset')()
#        }
#    ), where = getLoadingNamespace()
#)


########################################################################
### block RW sampler with multi-variate normal proposal distribution ###
########################################################################

sampler_RW_block <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, control) {
        ###  control list extraction  ###
        targetNodes    <- control$targetNodes
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        ###  node list generation  ###
        targetNodes_asScalars <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(targetNodes)
        ###  numeric value generation  ###
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory          <- c(0, 0)
        acceptanceRateHistory <- c(0, 0)
        d <- length(targetNodes_asScalars)
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
        my_setAndCalculate <- setAndCalculate(model, targetNodes)
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
            declare(origValueVector, double(1, d))
            origValueVector <- values(model, targetNodes)
            ######## propValueVector <- rmnorm(1, mean = origValueVector, chol = chol_propCov * scale)
            ######## propValueVector <- rmnorm(1, mean = origValueVector, sigma = chol_propCov %*% t(chol_propCov) * scale^2, method = 'chol')
            declare(normalVarVector, double(1, d))
            for(i in 1:d)     {   normalVarVector[i] <- rnorm(1, 0, 1)   }
            propValueMatrix <- asRow(origValueVector) + asRow(normalVarVector) %*% chol_propCov * scale
            propValueVector <- propValueMatrix[1, ]
            returnType(double(1))
            return(propValueVector)
        },
        
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(!adaptScaleOnly) {
                declare(newValues, double(1, d))

                newValues <- values(model, targetNodes)
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
                    gamma1 <- nfVar(my_calcAdaptationFactor, 'gamma1')
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
            nfMethod(my_calcAdaptationFactor, 'reset')()
        }
    ),  where = getLoadingNamespace()
)



#############################################################################
### RW_llFunction, does a RW, but using a generic log-likelihood function ###
#############################################################################

sampler_RW_llFunction <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, control) {
        ###  control list extraction  ###
        targetNode     <- control$targetNode
        adaptive       <- control$adaptive
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        llFunction     <- control$llFunction
        includesTarget <- control$includesTarget
        ###  node list generation  ###
        calcNodes <- model$getDependencies(targetNode)
        ###  nested function and function list definitions  ###
        mvInternal <- modelValues(model)
        RWControl <- list(targetNode = targetNode, adaptive = adaptive, adaptInterval = adaptInterval, scale = scale)
        targetRWSamplerFunction <- sampler_RW(model, mvInternal, RWControl)
        my_setAndCalculateOne <- setAndCalculateOne(model, targetNode)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    },
    
    run = function() {
        modelLP0 <- llFunction$run()
        if(!includesTarget)     modelLP0 <- modelLP0 + getLogProb(model, targetNode)
        propValue <- nfMethod(targetRWSamplerFunction, 'generateProposalValue')()
        my_setAndCalculateOne$run(propValue)
        modelLP1 <- llFunction()
        if(!includesTarget)     modelLP1 <- modelLP1 + getLogProb(model, targetNode)
        jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
        if(adaptive)     nfMethod(targetRWSamplerFunction, 'adaptiveProcedure')(jump)
    },
    
    methods = list(
        reset = function() {
            nfMethod(targetRWSamplerFunction, 'reset')()
        }
    ),  where = getLoadingNamespace()
)





####################################################################
### slice sampler (discrete or continuous) #########################
####################################################################

sampler_slice <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, control) {
        ###  control list extraction  ###
        targetNode    <- control$targetNode
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        width         <- control$sliceWidth
        maxSteps      <- control$sliceMaxSteps
        ###  node list generation  ###
        targetNodeAsScalar <- model$expandNodeNames(targetNode, returnScalarComponents = TRUE)
        if(length(targetNodeAsScalar) > 1)     stop(paste0('more than one targetNode: ', targetNode, '; cannot use slice sampler'))
        targetNodeFunction <- model$expandNodeNames(targetNode)
        calcNodes <- model$getDependencies(targetNode)
        ###  numeric value generation  ###
        widthOriginal <- width
        timesRan      <- 0
        timesAdapted  <- 0
        sumJumps      <- 0
        discrete      <- model$getNodeInfo()[[targetNodeFunction]]$isDiscrete()
    },
    
    run = function() {
        u <- getLogProb(model, calcNodes) - rexp(1, 1)    # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
        x0 <- model[[targetNode]]    # create random interval (L,R), of width 'width', around current value of targetNode
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
            model[[targetNode]] <<- value
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
        posteriorLogDensity <- nfMethod(conjugateSamplerFunction, 'getPosteriorLogDensity')()
        returnType(double())
        return(posteriorLogDensity)
    },
    where = getLoadingNamespace()
)

sampler_crossLevel <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, control) {
        ###  control list extraction  ###
        topNodes       <- control$topNodes
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        ###  node list generation  ###
        topNodes     <- model$expandNodeNames(topNodes)
        lowNodes     <- model$getDependencies(topNodes, self = FALSE, stochOnly = TRUE, includeData = FALSE)
        lowCalcNodes <- model$getDependencies(lowNodes)
        calcNodes    <- model$getDependencies(c(topNodes, lowNodes))
        ###  nested function and function list definitions  ###
        mvInternal <- modelValues(model)
        RWblockControl <- list(targetNodes = topNodes, adaptive = adaptive, adaptScaleOnly = adaptScaleOnly, adaptInterval = adaptInterval, scale = scale, propCov = propCov)
        topRWblockSamplerFunction <- sampler_RW_block(model, mvInternal, RWblockControl)
        lowConjugateSamplerFunctions <- nimbleFunctionList(sampler_BASE)
        lowConjugateGetLogDensityFunctions <- nimbleFunctionList(getPosteriorDensityFromConjSampler_virtual)
        for(iLN in seq_along(lowNodes)) {
            conjugacyResult <- model$checkConjugacy(lowNodes[iLN])
            if(is.null(conjugacyResult))     stop('non-conjugate lowNode \'', lowNodes[iLN], '\' in crossLevel updater')
            conjugateSamplerSpec <- samplerSpec(type = conjugacyResult$samplerType, control = conjugacyResult$control)
            lowConjugateSamplerFunctions[[iLN]] <- conjugateSamplerSpec$buildSampler(model, mvInternal)
            lowConjugateGetLogDensityFunctions[[iLN]] <- getPosteriorDensityFromConjSampler(lowConjugateSamplerFunctions[[iLN]])
        }
        my_setAndCalculateTop <- setAndCalculate(model, topNodes)
        my_decideAndJump <- decideAndJump(model, mvSaved, calcNodes)
    },
    
    run = function() {
    	    	
        modelLP0 <- getLogProb(model, calcNodes)
        propLP0 <- 0
        
        for(iSF in seq_along(lowConjugateGetLogDensityFunctions))  { propLP0 <- propLP0 + lowConjugateGetLogDensityFunctions[[iSF]]$run() }
        propValueVector <- nfMethod(topRWblockSamplerFunction, 'generateProposalVector')()
        
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
        if(adaptive)     nfMethod(topRWblockSamplerFunction, 'adaptiveProcedure')(jump)
    },
    
    methods = list(
        reset = function() {
            nfMethod(topRWblockSamplerFunction, 'reset')()
            for(iSF in seq_along(lowConjugateSamplerFunctions)) {
                nfMethod(lowConjugateSamplerFunctions[[iSF]], 'reset')()
            }
        }
    ),  where = getLoadingNamespace()
)
