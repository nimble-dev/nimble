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
    theta_length <- paramsTransform$getTransformedLength()
    if(theta_length > 1) theta_indices <- 1:theta_length
    else theta_indices <- c(1, -1)
    
    ## Indicator for removing the redundant index -1 in theta_indices
    one_time_fixes_done <- FALSE
    
    ## Default calculation method for AGHQuad
    computeMethod_ <- extractControlElement(control, "computeMethod", 2)
    useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)
   
		## Define the hyperparameter grid:
		## Probably need to say if(theta_length > 0) but surely that is implied...
    theta_grid_nfl <- nimbleFunctionList(QUAD_GRID_BASE)
    inner_grid_cache_nfl <- nimbleFunctionList(INNER_CACHE_BASE)

    gridOptions <- c("CCD", "AGHQ", "Custom")

    ## CCD 
		I_CCD <- 1
    theta_grid_nfl[[I_CCD]] <- buildQuadGrid(d = theta_length, nQuad = nQuadOuter, quadRule = "CCD")
    inner_grid_cache_nfl[[I_CCD]] <- inner_cache_methods(nre = 0, nGrid = 1)
    ## AGHQ -
		I_AGHQ <- 2   
    theta_grid_nfl[[I_AGHQ]] <- buildQuadGrid(d = theta_length, nQuad = nQuadOuter, quadRule = "AGHQ") 
    inner_grid_cache_nfl[[I_AGHQ]] <- inner_cache_methods(nre = 0, nGrid = 1)

    ## Custom User Provided -
    # I_USR <- 3
    # theta_grid_nfl[[I_USR]] <- buildQuadGrid(d = theta_length, nQuad = nQuadOuter, quadRule = "Custom") 
    # inner_grid_cache_nfl[[I_USR]] <- inner_cache_methods(nre = 0, nGrid = 1)
    
    ## Store the quadrature sums:
    marginalPostDensity <- c(-Inf, -Inf, -Inf)

    ## Naming convention thought... hyperQuadRule?
    ## gridRule?
    gridMethod <- which(hyperGrid == gridOptions)
    
		## Build marginal AGHQ grid to compute the hyperparameter marginals (integrate over pT-1 theta values).
    I_GRID <- 1
    I_AGHQ_ONE <- 2
    theta_grid_marg_nfl <- nimbleFunctionList(QUAD_GRID_BASE)
    ## This is set up to add other quad grids in the future. quadRule := "AGHQ" to start.
    theta_grid_marg_nfl[[I_GRID]] <- buildQuadGrid(d = (theta_length-1), nQuad = nQuadMarginal, quadRule = quadRuleMarginal)
    ## AGHQ Nodes to approx marginal along (e.g. Stringer)
    theta_grid_marg_nfl[[I_AGHQ_ONE]] <- buildQuadGrid(d = 1, nQuad = nQuadMarginal, quadRule = "AGHQ") ## Only AGHQ currently.
		    
    ## Transformation methods for the quadrature grid:
    gridTransform <- quadGridTransformMethods()
        
    ## Cached values for convenience:
    ## For marginal distributions in AGHQ over d-1.
    pTransformFix <- 0 
    indexFix <- 0

		## We will want to cache the standard deviation skew terms.
    ## Default will be not to skew (e.g. 1)
		covTheta <- matrix(0, nrow=theta_length, ncol=theta_length)
    cholNegHess <- matrix(0, theta_length, theta_length)
    eigenVecNegHess <- matrix(0, theta_length, theta_length)
    eigenValNegHess <- rep(0, theta_length)
    if(theta_length == 1)
      eigenValNegHess <- c(0,-1)
    eigenCalculated <- FALSE
    choleskyCalculated <- FALSE
   
		skewedStdDev <- matrix(1, nrow = theta_length, ncol = 2)
    ## Points for Asymmetric Gaussian Interpolation (integration free...?) Taken from INLA Code:
		extraPoints <- c(-15.0, -10.0, -7.0, 7.0, 10.0, 15.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0)
		aghqPoints <- pracma::gaussHermite(51)$x*sqrt(2)
		zMargGrid <- sort(unique(c(extraPoints, aghqPoints)))
    nzMargGrid <- length(zMargGrid)

    ## Some cached values for summary statistics and reporting:
    marg_P <- matrix(0, nrow=nzMargGrid, ncol = npar)
    marg_theta <- array(0, c(theta_length, nzMargGrid, 2) )
    post_sims <- matrix(0, nrow=3, ncol = nre)  ## Just initialize it.

    ## Optim info:
    calcMode <- FALSE
    thetaMode <- numeric(theta_length)
    thetaNegHess <- matrix(0, nrow = theta_length, ncol = theta_length)
    logPostProbMode <- 0
    logDetNegHessTheta <- 0
		
    ## Indicator for removing the redundant index -1 in theta_indices
    one_time_fixes_done <- FALSE
  },
  run = function(){},
  methods = list(
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(theta_length == 1){
        theta_indices <<- numeric(length = 1, value = 1)
        thetaMode <<- numeric(length = 1, value = 0)
        eigenValNegHess <<- numeric(length = 1, value = 0)
      }
      if(nre == 1)
        thetaMode <<- numeric(length = 1, value = 1)
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
      thetaMode <<- optRes$par
      thetaNegHess <<- optRes$hessian
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
    calcEigen = function(){
      E <- eigen(thetaNegHess, symmetric = TRUE) ## Should be symmetric...
      eigenVecNegHess <<- E$vectors
      eigenValNegHess <<- E$values
      logDetNegHessTheta <<- sum(log(eigenValNegHess))
      eigenCalculated <<- TRUE
    },
    calcCholesky = function(){
      cholNegHess <<- chol(thetaNegHess)
      logDetNegHessTheta <<- 2 * sum(log(diag(cholNegHess)))
      choleskyCalculated <<- TRUE
    },
    ## Transform from standard (z) to param transform (theta) scale.
    z_to_theta = function(z = double(1), method = character(0, default = "spectral")){
      if(method == "spectral"){
        if(!eigenCalculated)  ## *** eigenCached
          calcEigen()
        theta <- gridTransform$z_to_theta_spectral(z, thetaMode, eigenVecNegHess, eigenValNegHess)
      }else{
        if(!choleskyCalculated) ## *** cholCached
          calcCholesky()
        theta <- thetaMode + backsolve(cholNegHess, z)
      }
      returnType(double(1))
      return(theta)
    },
    ## Transform from param transform (theta) to standard (z) scale.
    theta_to_z = function(theta = double(1), method = character(0, default = "spectral")){
      if(method == "spectral"){
        if(!eigenCalculated) 
          calcEigen()
        z <- gridTransform$theta_to_z_spectral(theta, thetaMode, eigenVecNegHess, eigenValNegHess)
      }else{
        if(!choleskyCalculated)
          calcCholesky()
        z <- (theta - thetaMode) %*% t(cholNegHess)
      }
      returnType(double(1))
      return(z)
    },
    calcSkewedSD = function() {
      ## Require the grid to have been built and the mode found.
      buildHyperGrid() ## *** should be just CCD
      for( i in 1:theta_length){
        z <- numeric(value = 0, length = theta_length)
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
      # marg <- logPostProbMode + 0.5*theta_length*log(2*pi) -
        # 0.5*(logDetNegHessTheta) - sum(log(skewedStdDev[,1] * skewedStdDev[,2]))
      ## *** What Paul thinks it should be. ***
      marg <- logPostProbMode + 0.5*theta_length*log(2*pi) -
        0.5*(logDetNegHessTheta) + sum(log((skewedStdDev[,1] + skewedStdDev[,2])/2))
      returnType(double())
			return(marg)
		},
		## This is the meat and potatoes for being able to make inference on the latent nodes.
    ## Calculate theta on the quadrature grid points. AGHQ or CCD.
		## Stores all values we need for simulation inference on the latent nodes.
		calcHyperGrid = function(){
      buildHyperGrid()
      ## Transform the grid from z to theta.
      nGrid <- theta_grid_marg_nfl[[gridMethod]]$gridSize()

      ans <- 0
      ## Now fill in the grid values.
      for( i in 1:nGrid ){
        ## Operations at the mode:
        if(i == theta_grid_marg_nfl[[gridMethod]]$modeI()){
          wgt <- theta_grid_marg_nfl[[gridMethod]]$weighti(indx = i)
          inner_grid_cache_nfl[[gridMethod]]$cache_inner_mode(innerMethods$get_inner_mode(atOuterMode = 1), indx = i)
          inner_grid_cache_nfl[[gridMethod]]$cache_inner_negHess(innerMethods$get_inner_cholesky(atOuterMode = 1), indx = i)
          inner_grid_cache_nfl[[gridMethod]]$cache_weights(wgt, indx = i)
          ans <- ans + wgt
        }else{
          node <- theta_grid_marg_nfl[[gridMethod]]$nodei(indx = i)
          wgt <- theta_grid_marg_nfl[[gridMethod]]$weighti(indx = i)
          
          ## Skew the CCD values:
          ## I really think these need to be investigated for quadrature weights... Think more later.
          if(gridMethod == I_CCD){
            for(d in 1:theta_length) {
              node[d] <- node[d]*skewedStdDev[d, step(node[d]) + 1] ## negative skew column 1, positive skew column 2
            }
          }
          ## Transform to theta scale:
          node <- z_to_theta(node, transformMethod)
          thetaLogPostDens <- innerMethods$calcPostLogDens_pTransformed(node)
          wgt_dens <- wgt*exp(thetaLogPostDens - logPostProbMode)
          ## Marginal sum:
          ans <- ans + wgt_dens
          
          ## Cache everything for simulation:
          inner_grid_cache_nfl[[gridMethod]]$cache_inner_mode(innerMethods$get_inner_mode(atOuterMode = 0), indx = i)
          inner_grid_cache_nfl[[gridMethod]]$cache_inner_negHess(innerMethods$get_inner_cholesky(atOuterMode = 0), indx = i)
          inner_grid_cache_nfl[[gridMethod]]$cache_weights(wgt_dens, indx = i)
        }
        ## Add a convergence check.
      }
      ## Marginal log posterior density, a normalizing constant for other methods.
      marginalPostDensity[gridMethod] <<- log(ans) + logPostProbMode - 0.5 * logDetNegHessTheta
		},
		## Quadrature based marginal log-likelihood
		## Probably not particularly accurate for CCD.
		calcMarginalLogLikQuad = function(){
      if(gridMethod == I_CCD) 
        print("Warning:  CCD not theoretically supported to compute marginal. Switch to AGHQ if accuracy is required." )
      ## Need to check if `calcHyperGrid()` has actually been called.
			returnType(double())
			return(marginalPostDensity[gridMethod)
		},
    ## *** PVDB START HERE NEXT!
    ## Marginals AGHQ from Stringer et al.
    ## *** Investigate pruning for AGHQ.
    ## *** Is this a good name?
		findMarginalPosteriorDensity = function(pIndex = integer(), nGrid = integer(0, default = 3), ## *** update naming of nQp -> nGrid.
                                            nQuad = integer(0, default = 3))  ## *** update naming of nQinner -> nQuad.
		{
      ## Build the quadrature grid points:
      theta_grid_marg_nfl[[I_GRID]]$buildGrid()
      theta_grid_marg_nfl[[I_AGHQ_ONE]]$buildGrid()

      if( gridBuiltMarg == 0 ){
        theta_grid_marg$buildGrid()
        gridBuiltMarg <<- 1
      }
      if(nQinner != theta_grid_marg$getGridSize()){
        theta_grid_marg$resetGrid(nQinner, keepInner = 0) ## Don't need the inner mode and cholesky here.
      }
      ## Make sure we have already calculated the mode:
      if(postModeFound == 0) findPostMode(pStartTransform = rep(0, theta_length), method = "BFGS", hessian = TRUE, buildGrid = TRUE)
      if(theta_grid_nfl[[gridMethod]]$calcCheck(0) == 0){
        addModeGridInfo() ## If log dens at the mode hasn't been calculated. Do this to make sure covTheta gets calculated.
      }
      
      ## 1D quadrature to evaluate the theta on.
      zwtheta <- theta_grid_marg$buildAGHQOne(nQp)	## nodes + weights
      stdDev <- sqrt(covTheta[pIndex, pIndex])
      
      # Get indices excluding pIndex.
      inds <- theta_indices[theta_indices != pIndex]

      # Initial values really matter.
      initPTransform  <- thetaMode[inds]
      subsetNegHess <- negHessthetaMode[inds,inds] ## Can't passByMap
      
      res <- matrix(0, nrow = nQp, ncol = 3)
      ## For each value of theta, we need to do AGHQ which means 
      ## finding the mode of the other parameters.
      for( j in 1:nQp ){
        res[j,1] <- zwtheta[j,1]*stdDev + thetaMode[pIndex]
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
			thetai <- rep(0, theta_length)

			for( i in 1:nzMargGrid ){	# Known fixed # of points
				thetai[pIndex] <- thetaMode[pIndex] + zMargGrid[i] * stdDev
        marg_theta[pIndex, i, 1] <<- thetai[pIndex]
				## Find the conditional mean:
				for( j in 1:theta_length ){
					if(j != pIndex){
						thetai[j] <- thetaMode[j] + covTheta[pIndex, j]/covTheta[pIndex, pIndex] * 
							(thetai[pIndex] - thetaMode[pIndex])
					}
				}
				## Calculate asymmetric Gaussian:
				zi <- theta_grid_nfl[[gridMethod]]$theta_to_z(thetai)
        logDens <- 0
				for( j in 1:theta_length ){
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
			findPostMode(pStartTransform = rep(0, theta_length), method = "BFGS", hessian = TRUE, buildGrid = TRUE)

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
      for( i in 1:theta_length ){
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