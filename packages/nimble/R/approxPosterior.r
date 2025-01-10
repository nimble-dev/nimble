## Code to approximate the posterior distribution built on top of inner Laplace.
## Will start to build the buildApproxPosterior functionality from the other branch here 
## as discussed with Chris P.
## Note that we will separate fixed effects that are normally distributed from the hyperparameters.
## * diff from Laplace.
buildApproxPosterior <- nimbleFunction(
  name = 'ApproxPost',
  setup = function(model, hyperParamNodes, fixedEffectsNodes, randomEffectsNodes, calcNodes,  ## *** Figure out latent nodes.
                   calcNodesOther, control = list()) {
    split <- extractControlElement(control, 'split', TRUE)
    check <- extractControlElement(control, 'check', TRUE)
    innerOptimWarning <- extractControlElement(control, 'innerOptimWarning', FALSE)
    
    hyperGrid <- extractControlElement(control, 'hyperGrid', 'CCD')
    nQuadOuter <- extractControlElement(control, 'nQuadOuter', 3)
    nQuadInner <- extractControlElement(control, 'nQuadInner', 1)
    nQuadMarginal <- extractControlElement(control, 'nQuadOuterMarginalHyperparams', 3)
    quadRuleMarginal <- extractControlElement(control, 'outerMarginalQuadRule', 'AGHQ')
    transformMethod <- extractControlElement(control, 'quadTransform', 'spectral')
    
    ## Default starting value for Approx Posterior set here and passed to Laplace.
    ## zero makes sense but should test others on the AGHQ grid and see how they work.
    control$innerOptimStart <- extractControlElement(control, "innerOptimStart", "zero")

    innerMethods <- buildAGHQ(model, nQuadInner, hyperParamNodes, c(fixedEffectsNodes, randomEffectsNodes), 
                                 calcNodes, calcNodesOther, control)
    
    paramNodes <- innerMethods$paramNodes
    randomEffectsNodes <- innerMethods$randomEffectsNodes

    ## Need to check this as it's is now computed in the "buildAGHQ" function:
    scalarRENodes <- model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE)
    nre <- length(scalarRENodes)

    ## Outer optimization settings
    outerOptimControl_   <- nimOptimDefaultControl()
    optimControlArgNames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", 
                              "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    if(!is.null(control$outerOptimControl)){
      validNames <- intersect(names(control$outerOptimControl), optimControlArgNames)
      numValidNames <- length(validNames)
      if(numValidNames > 0){
        for(i in 1:numValidNames){
          outerOptimControl_[[validNames[i]]] <- control$outerOptimControl[[validNames[i]]]
        }   
      }
    }
    outerOptimControl_$fnscale <- -1
    
    paramNodesAsScalars <- model$expandNodeNames(innerMethods$paramNodes, returnScalarComponents = TRUE)
    npar <- length(paramNodesAsScalars)
    paramNodesAsScalars_vec <- paramNodesAsScalars
    if(npar == 1) paramNodesAsScalars_vec <- c(paramNodesAsScalars, "_EXTRA_")
    paramNodesAsScalars_first <- paramNodesAsScalars[1]
    if(npar == 1) p_indices <- c(1, -1)
    else p_indices <- 1:npar
    
    
    ## APPROX SPECIFIC SetUp
    ## Automated transformation for parameters:
    ## ***  Need to connect with Daniel and Chris about parameter transformations for marginals:
    ## *** Key is to do regular transformations on the whole, but individual for the marginals.
    paramsTransform <- parameterTransformI(model, paramNodes, control = list(allowDeterm = FALSE))
    pTransform_length <- paramsTransform$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
    
    ## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    
    ## Default calculation method for AGHQuad
    computeMethod_ <- extractControlElement(control, "computeMethod", 2)
    useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)

    ## Automated transformation for parameters
    paramsTransformIndividual <- parameterTransformI(model, paramNodes, control = list(allowDeterm = allowNonPriors))
    ## Setup specific to approx posteriors:
    pTransform_length <- paramsTransformIndividual$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
   
		## Define the hyperparameter grid:
		## Probably need to say if(pTransform_length > 0) but surely that is implied...
    theta_grid_nfl <- nimbleFunctionList(QUAD_GRID_BASE)
    inner_grid_cache_nfl <- nimbleFunctionList(INNER_CACHE_BASE)

    gridOptions <- c("CCD", "AGHQ", "Custom")

    ## CCD 
		I_CCD <- 1
    theta_grid_nfl[[I_CCD]] <- buildQuadGrid(d = pTransform_length, nQuad = nQuadOuter, quadRule = "CCD")
    inner_grid_cache_nfl_[[I_CCD]] <- inner_cache_methods(nre = 0, nGrid = 1)
    ## AGHQ -
		I_AGHQ <- 2   
    theta_grid_nfl[[I_AGHQ]] <- buildQuadGrid(d = pTransform_length, nQuad = nQuadOuter, quadRule = "AGHQ") 
    inner_grid_cache_nfl_[[I_AGHQ]] <- inner_cache_methods(nre = 0, nGrid = 1)

    ## Custom User Provided -
    # I_USR <- 3
    # theta_grid_nfl[[I_USR]] <- buildQuadGrid(d = pTransform_length, nQuad = nQuadOuter, quadRule = "Custom") 
    # inner_grid_cache_nfl_[[I_USR]] <- inner_cache_methods(nre = 0, nGrid = 1)

    ## Naming convention thought... hyperQuadRule?
    ## gridRule?
    gridMethod <- which(hyperGrid == gridOptions)
    
    ## Default for CCD. Will update when building grid.
		thetaModeIndex <- 1
		nGrid <- 1

		## Build marginal AGHQ grid to compute the hyperparameter marginals (integrate over pT-1 theta values).
    I_GRID <- 1
    I_AGHQ_ONE <- 2
    theta_grid_marg_nfl <- nimbleFunctionList(QUAD_GRID_BASE)
    ## This is set up to add other quad grids in the future. quadRule := "AGHQ" to start.
    theta_grid_marg_nfl[[I_GRID]] <- buildQuadGrid(d = (pTransform_length-1), nQuad = nQuadMarginal, quadRule = quadRuleMarginal)
    ## AGHQ Nodes to approx marginal along (e.g. Stringer)
    theta_grid_marg_nfl[[I_AGHQ_ONE]] <- buildQuadGrid(d = 1, nQuad = nQuadMarginal, quadRule = "AGHQ") ## Only AGHQ currently.
		    
    ## Cached values for convenience:
    ## For marginal distributions in AGHQ over d-1.
    pTransformFix <- 0 
    indexFix <- 0

		## We will want to cache the standard deviation skew terms.
    ## Default will be not to skew (e.g. 1)
		covTheta <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
    cholNegHess <- matrix(0, pTransform_length, pTransform_length)
    eigenVecNegHess <- matrix(0, pTransform_length, pTransform_length)
    eigenValNegHess <- rep(0, pTransform_length)
    if(pTransform_length == 1)
      eigenValNegHess <- c(0,-1)
    calcEigen <- FALSE
    calcChol <- FALSE
   
		skewedStdDev <- matrix(1, nrow = pTransform_length, ncol = 2)
    ## Points for Asymmetric Gaussian Interpolation (integration free...?) Taken from INLA Code:
		extraPoints <- c(-15.0, -10.0, -7.0, 7.0, 10.0, 15.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0)
		aghqPoints <- pracma::gaussHermite(51)$x*sqrt(2)
		zMargGrid <- sort(unique(c(extraPoints, aghqPoints)))
    nzMargGrid <- length(zMargGrid)

    ## Some cached values for summary statistics and reporting:
    marg_P <- matrix(0, nrow=nzMargGrid, ncol = npar)
    marg_theta <- array(0, c(pTransform_length, nzMargGrid, 2) )
    post_sims <- matrix(0, nrow=3, ncol = nre)  ## Just initialize it.

    ## Optim info:
    calcMode <- FALSE
    pTransformMode <- numeric(pTransform_length)
    pTransformNegHess <- matrix(0, nrow = pTransform_length, ncol = pTransform_length)
    logPostProbMode <- 0
    logDetNegHesspTransform <- 0
		
    ## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
  },
  run = function(){},
  methods = list(
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(pTransform_length == 1){
        pTransform_indices <<- numeric(length = 1, value = 1)
        pTransformPostMode <<- numeric(length = 1, value = 0)
        eigenValNegHess <<- numeric(length = 1, value = 0)
      }
      if(nre == 1)
        pTransformPostMode <<- numeric(length = 1, value = 1)
      one_time_fixes_done <<- TRUE
    },
    ## *** Posterior mode for hyperparameters. findMAP
    posteriorMode = function(pStart = double(1, default = Inf),
                       method  = character(0, default = "BFGS"),
                       hessian = logical(0, default = TRUE),
                       parscale = character(0, default = "transformed")){
      optRes <- innerMethods$optimize(pStart = pStart, prior = TRUE, jacobian = TRUE, 
          method  = method, hessian = TRUE, parscale = parscale)
      calcMode <<- TRUE
      pTransformMode <<- optRes$par
      pTransformNegHess <<- optRes$hessian
      logPostProbMode <<- optRes$value
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Build hyper quad grid and cache system.
    buildHyperGrid = function(){
			theta_grid_nfl[[gridMethod]]$buildGrid(nQuadOuter)
      inner_grid_cache_nfl[[gridMethod]]$buildCache(nGridUpdate = theta_grid_nfl[[gridMethod]]$gridSize())
      if(calcMode)
        posteriorMode(rep(Inf, npar), method = "BFGS", hessian = TRUE, parscale = "transformed")
		},
		changeHyperGrid = function(quadRule = character(0, default = "AGHQ"), 
                               nQuadUpdate = integer(0, default = 3)){
      hyperGrid <<- quadRule
      gridMethod <<- which(quadRule == gridOptions)
      nQuadOuter <<- nQuadUpdate
			buildHyperGrid()
    },
    findEigen = function(){
      E <- eigen(pTransformNegHess, symmetric = TRUE) ## Should be symmetric...
      eigenVecNegHess <<- E$vectors
      eigenValNegHess <<- E$values
      logDetNegHesspTransform <<- sum(log(eigenValNegHess))
      calcEigen <<- TRUE
    },
    findCholesky = function(){
      cholNegHess <<- chol(pTransformNegHess)
      logDetNegHesspTransform <<- 2 * sum(log(diag(cholNegHess)))
      calcChol <<- TRUE
    },
    ## Transform from standard (z) to param transform (theta) scale.
    z_to_theta = function(z = double(1), method = character(0, default = "spectral")){
      if(quadTransform_ == "spectral"){
        if(!calcEigen)  ## *** eigenCached
          findEigen()
        theta <- numeric(value = 0, length = pTransform_length)
        for( i in 1:pTransform_length ){
          theta[i] <- pTransformMode[i] + sum(eigenVecNegHess[,i] * z) / sqrt(eigenValNegHess[i]) ## *** cache as sqrt
        }
      }else{
        if(!calcChol) ## *** cholCached
          findCholesky()
        theta <- pTransformMode + backsolve(cholNegHess, z)
      }
      return(theta)
      returnType(double(1))
    },
    ## Transform from param transform (theta) to standard (z) scale.
    theta_to_z = function(theta = double(1), method = character(0, default = "spectral")){
      if(quadTransform_ == "spectral"){
        if(!calcEigen) 
          findEigen()
        z <- numeric(value = 0, length = pTransform_length)
        for( i in 1:pTransform_length ){
          z[i] <- pTransformMode[i] + sum(eigenVecNegHess[,i] * theta) * sqrt(eigenValNegHess[i])
        }
      }else{
        if(!calcChol)
          findCholesky()
        z <- (theta - pTransformMode) %*% t(cholNegHess) 
      }
      return(z)
      returnType(double(1))
    },
    calcSkewedSD = function() {
      ## Require the grid to have been built and the mode found.
      buildHyperGrid() ## *** should be just CCD
      for( i in 1:pTransform_length){
        z <- numeric(value = 0, length = pTransform_length)
        z[i] <- -sqrt(2)
        theta <- z_to_theta(z, transformMethod)
        logDens2Neg <- innerMethods$calcPostLogDens_pTransformed(pTransform = theta)
        skewedStdDev[i, 1] <<- sqrt(2 / (2.0 * (logPostProbMode-logDens2Neg))) 	## numerator (-sqrt(2)) ^2
        z[i] <- sqrt(2)
        theta <- z_to_theta(z, transformMethod)
        logDens2Pos <- innerMethods$calcPostLogDens_pTransformed(pTransform = theta)
        skewedStdDev[i, 2] <<- sqrt(2 / (2.0 * (logPostProbMode-logDens2Pos))) 	## numerator (-sqrt(2)) ^2
      }
    },
    getSkewedStdDev = function(){
      returnType(double(2))
      return(skewedStdDev)
    },
    calcMarginalLogLikApprox = function(){
      ## Line 2748 in r-inla/blob/devel/gmrflib/approx-inference.c
      ## Commit # ef4eb20
      # marg <- logPostProbMode + 0.5*pTransform_length*log(2*pi) -
        # 0.5*(logDetNegHesspTransform) - sum(log(skewedStdDev[,1] * skewedStdDev[,2]))
      ## *** What Paul thinks it should be. ***
      marg <- logPostProbMode + 0.5*pTransform_length*log(2*pi) -
        0.5*(logDetNegHesspTransform) + sum(log((skewedStdDev[,1] + skewedStdDev[,2])/2))
      returnType(double())
			return(marg)
		},
    ### **** PVDB working here next.
		## This is a loop for calculating theta on the grid points.
		## It stores all values we need for simulation inference on the fixed and random-effects.
		## Once this is called, we can then simulate.
		calcHyperGrid = function(){
      buildHyperGrid()
      ## Transform the grid from z to theta.
      theta_grid_nfl[[gridMethod]]$transformGrid(skewedStdDev)
      ## Now fill in the grid values.
      for( i in 1:nGrid ){
				if(theta_grid_nfl[[gridMethod]]$calcCheck(i) == 0 ){
          theta <- theta_grid_nfl[[gridMethod]]$getTheta(i)
					theta_grid_nfl[[gridMethod]]$saveLogDens(i=i, logDensity = calcPostLogProb_pTransformed(theta))

          ## If it didn't converge. Warn and save as -Inf
          # if(checkInnerConvergence() == 0) {
            # print("Warning: Inner marginalization was not successful. Assuming log density = -Inf for that step." )
            # theta_grid_nfl[[gridMethod]]$saveLogDens(i=i, logDensity = -Inf)
          # }
          ## Save values for inference on fixed- random-effects
          theta_grid_nfl[[gridMethod]]$saveInnerCholesky(i, get_inner_cholesky(atOuterMode = 0))
          theta_grid_nfl[[gridMethod]]$saveInnerMode(i, get_inner_mode(atOuterMode = 0))
				}
			}
		},
		## Quadrature based marginal log-likelihood
		## Probably not particularly accurate for CCD.
		calcMarginalLogLikQuad = function(){
      if(gridMethod == I_CCD) 
        print("Warning:  CCD not theoretically supported to compute marginal. Switch to AGHQ if accuracy is required." )
      marg <- theta_grid_nfl[[gridMethod]]$quadSum()
			returnType(double())
			return(marg)
		},
    ##
		getMarginalQuadratureGrid = function(){
			nG <- theta_grid_marg$getGridSize()
			ans <- matrix(0, nrow = nG, ncol = pTransform_length + 1)
			for(i in 1:nG){
				ans[i, 1:(pTransform_length-1)] <- theta_grid_marg$getTheta(i)
				ans[i, pTransform_length] <- theta_grid_marg$getWeights(i)
				ans[i, pTransform_length + 1] <- theta_grid_marg$getLogDensity(i)
			}
			returnType(double(2))
			return(ans)
		},
    ## Marginals AGHQ from Stringer et al.
    ## *** Investigate pruning for AGHQ.
		findMarginalPosteriorDensity = function(pIndex = integer(), nQp = integer(0, default = 3), ## *** update naming of nQp -> nGrid.
                                            nQinner = integer(0, default = 3))  ## *** update naming of nQinner -> nQuad.
		{
      if( gridBuiltMarg == 0 ){
        theta_grid_marg$buildGrid()
        gridBuiltMarg <<- 1
      }
      if(nQinner != theta_grid_marg$getGridSize()){
        theta_grid_marg$resetGrid(nQinner, keepInner = 0) ## Don't need the inner mode and cholesky here.
      }
      ## Make sure we have already calculated the mode:
      if(postModeFound == 0) findPostMode(pStartTransform = rep(0, pTransform_length), method = "BFGS", hessian = TRUE, buildGrid = TRUE)
      if(theta_grid_nfl[[gridMethod]]$calcCheck(0) == 0){
        addModeGridInfo() ## If log dens at the mode hasn't been calculated. Do this to make sure covTheta gets calculated.
      }
      
      ## 1D quadrature to evaluate the theta on.
      zwtheta <- theta_grid_marg$buildAGHQOne(nQp)	## nodes + weights
      stdDev <- sqrt(covTheta[pIndex, pIndex])
      
      # Get indices excluding pIndex.
      inds <- pTransform_indices[pTransform_indices != pIndex]

      # Initial values really matter.
      initPTransform  <- pTransformPostMode[inds]
      subsetNegHess <- negHesspTransformPostMode[inds,inds] ## Can't passByMap
      
      res <- matrix(0, nrow = nQp, ncol = 3)
      ## For each value of theta, we need to do AGHQ which means 
      ## finding the mode of the other parameters.
      for( j in 1:nQp ){
        res[j,1] <- zwtheta[j,1]*stdDev + pTransformPostMode[pIndex]
        indexFix <<- pIndex               ## Need these values fixed for calcPostLogProb_pTransformedj.
        pTransformFix <<- res[j,1] 

        ## If this is the mode then we know some information:
        if(zwtheta[j,1] == 0){
          theta_grid_marg$saveOptimInfo(pTransformMax = initPTransform, maxLogDens = logPostProbMode, negHessian = subsetNegHess)
        }else{
          ## *** Need to make sure this converged. Extreme AGHQ points may be unlikely to.
          findPostModeFixedj(pStartTransform = initPTransform, 
              j=pIndex, pTransformFixed=res[j,1], method = "BFGS", 
              hessian = TRUE, buildGrid = 1)
        }
        theta_grid_marg$transformGrid(skewSD = skewedStdDev) # Won't Skew as this is AGHQ.
        for( i in 1:nQinner ){  ## *** nQinner^d buggy code confuses Chris!!
          if( i != theta_grid_marg$getThetaModeIndex()) {
            theta_grid_marg$saveLogDens(i=i, logDensity = calcPostLogProb_pTransformedj(theta_grid_marg$getTheta(i)))
            ## If it didn't converge. Warn and save as -Inf
            # if(checkInnerConvergence() == 0) {
              # print("Warning: Inner marginalization was not successful. Assuming log density = -Inf for that step." )
              # theta_grid_marg$saveLogDens(i=i, logDensity = -Inf)
            # }
          }
        }
        res[j,2] <- theta_grid_marg$quadSum()
      }
      # Because thetaj values are AGHQ, we can normalize to get the proper posterior prob.
      margj <- sum(exp(res[,2] - logPostProbMode)*zwtheta[,2])
      lognormconst <- log(margj) + logPostProbMode + log(stdDev)
      res[,3] <- res[,2] - lognormconst
      returnType(double(2))
      return(res)
		},
    ## This can't be until I've built the CCD grid. Needs addModeGridInfo() to be run
    ## so that we have covTheta.
    ## Should also ensure that if they plan to skew the grid that is also done.
		findMarginalHyperIntFree = function(pIndex = integer())
		{
			stdDev <- sqrt(covTheta[pIndex, pIndex])
			thetai <- rep(0, pTransform_length)

			for( i in 1:nzMargGrid ){	# Known fixed # of points
				thetai[pIndex] <- pTransformPostMode[pIndex] + zMargGrid[i] * stdDev
        marg_theta[pIndex, i, 1] <<- thetai[pIndex]
				## Find the conditional mean:
				for( j in 1:pTransform_length ){
					if(j != pIndex){
						thetai[j] <- pTransformPostMode[j] + covTheta[pIndex, j]/covTheta[pIndex, pIndex] * 
							(thetai[pIndex] - pTransformPostMode[pIndex])
					}
				}
				## Calculate asymmetric Gaussian:
				zi <- theta_grid_nfl[[gridMethod]]$theta_to_z(thetai)
        logDens <- 0
				for( j in 1:pTransform_length ){
					side <- 2
					if(zi[j] <= 0) side <- 1
					logDens <- logDens - 0.5*(zi[j]/skewedStdDev[j, side])^2
				}
        marg_theta[pIndex, i, 2] <<- logDens
			}
			## Note to self: Explore normalization... Should we use the spline? Can we know it analytically?
			returnType(double(2))
			return(marg_theta[pIndex,,])			
		},
    marginalTransformedSplineDensity = function(pIndex = integer()) {
      returnType(double(2))
      return(marginalSplineR(marg_theta[pIndex, , 1], marg_theta[pIndex, , 2]))
    },
		findApproxPosterior = function(){	
			## Basic approx posterior steps:
			##-------------------------------------
			
			## 1) Find posterior mode + build hyperparameter grid.
			## Values are saved to theta_grid_nfl and locally cached.
			## Should think about initialization.
			findPostMode(pStartTransform = rep(0, pTransform_length), method = "BFGS", hessian = TRUE, buildGrid = TRUE)

			## 2) Calculate skew and if gridMethod == CCD, skew grid.
			calcSkewedSD()	## Cache "skewedStdDev"
			if(gridMethod == I_CCD) theta_grid_nfl[[I_CCD]]$transformGrid(skewSD = skewedStdDev)

			## 3) Calculate the density on the grid points. Saves to theta_grid_nfl
			## This is used for inference on the fixed and random-effects.
			calcHyperGrid()
			
			## 4) Calculate Marginal log-Likelihood
			## Based on asymmetric Gaussian assumption of the marginal of theta.
			marginalAG <- calcMarginalLogLikApprox()
			marginalQuad <- calcMarginalLogLikQuad()	## This one probably make sense only for AGHQ

			## 5) Marginals for theta: 
      ## Automatically do integration free for now. Values are cached. 
      ## User can compare with manually doing aghq.
      for( i in 1:pTransform_length ){
        findMarginalHyperIntFree(i)
      }
      
			## 6) Marginals for Fixed and Random-Effects: 
			## Only simulation based. Will assume 10000? User can add more or do less after testing.
      # simulateLatentEffects(10000) ## I've now added this into the caching system.
		}
    ## Add posterior mode function: 
    ## Add posterior mode function: for holding one fixed. Need to update the marginal likelihood func for that.
    
  )
)