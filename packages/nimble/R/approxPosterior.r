## Code to approximate the posterior distribution built on top of inner Laplace.
## Will start to build the buildApproxPosterior functionality from the other branch here 
## as discussed with Chris P.
buildApproxPosterior <- nimbleFunction(
  name = 'ApproxPost',
  setup = function(model, hyperParamNodes, randomEffectsNodes, calcNodes,
                   calcNodesOther, control = list()) {
    split <- extractControlElement(control, 'split', TRUE)
    check <- extractControlElement(control, 'check', TRUE)
    innerOptimWarning <- extractControlElement(control, 'innerOptimWarning', FALSE)
    
    hyperGrid <- extractControlElement(control, 'hyperGrid', "ccd")
    nQuadOuter <- extractControlElement(control, 'nQuadOuter', 3)
    nQuadInner <- extractControlElement(control, 'nQuadInner', 1)

    ## Default starting value for Approx Posterior set here and passed to Laplace.
    ## zero makes sense but should test others on the AGHQ grid and see how they work.
    control$innerOptimStart <- extractControlElement(control, "innerOptimStart", "zero")

    innerMethods <- buildAGHQ(model, nQuadInner, hyperParamNodes, randomEffectsNodes, 
                                 calcNodes, calcNodesOther, control)
    
    paramNodes <- innerMethods$paramNodes
    randomEffectsNodes <- innerMethods$randomEffectsNodes

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
    
    ## Automated transformation for parameters
    paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = FALSE))
    pTransform_length <- paramsTransform$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
    
    ## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    
    ## Default calculation method for AGHQuad
    computeMethod_ <- extractControlElement(control, "computeMethod", 2)
    useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)

    # Possible future feature
    #if(is.null(control$allowNonPriors)) allowNonPriors <- FALSE else  allowNonPriors <- control$allowNonPriors
    allowNonPriors <- FALSE

    ## Automated transformation for parameters
    paramsTransformIndividual <- parameterTransformI(model, paramNodes, control = list(allowDeterm = allowNonPriors))
    ## Setup specific to approx posteriors:
    pTransform_length <- paramsTransformIndividual$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
   
		## Define the hyperparameter grid:
		## Probably need to say if(pTransform_length > 0) but surely that is implied...
		I_CCD <- 1
		I_AGHQ <- 2
		## To be able to compare accuracy we will build both grids, or at least an empty version of aghq if not requested.
    theta_grid_nfl <- nimbleFunctionList(GRID_BASE)
		if(approxMethods$hyperGrid == 'ccd'){
			gridMethod <- I_CCD
			theta_grid_nfl[[I_CCD]] <- buildQuadGrid(d = pTransform_length, nQ = 1, method = I_CCD, nre = nre)
			theta_grid_nfl[[I_AGHQ]] <- buildQuadGrid(d = pTransform_length, nQ = 1, method = I_AGHQ, nre = nre)
		}else{
			gridMethod <- I_AGHQ
			theta_grid_nfl[[I_AGHQ]] <- buildQuadGrid(d = pTransform_length, nQ = approxMethods$nQuadAGHQ, method = I_AGHQ, nre = nre)					
		}
		gridBuilt <- c(0,0)
		thetaModeIndex <- 1	## Default for CCD. Will update when building grid.
		nQuadAGHQ <- approxMethods$nQuadAGHQ
		nGrid <- theta_grid_nfl[[gridMethod]]$getGridSize()

		## Build marginal AGHQ grid to integrate over pT-1 theta values.
		theta_grid_marg <- buildQuadGrid(d = (pTransform_length-1), nQ = approxMethods$nQuadAGHQ, method = I_AGHQ, nre = nre)
		gridBuiltMarg <- 0
		
    ## Cache values from optim step.
		pTransformPostMode <- rep(0, pTransform_length)
    logPostProbMode <- 0
		negHesspTransformPostMode <- matrix(0,  nrow = pTransform_length, ncol = pTransform_length)
    logDetNegHesspTransform <- 0
    
		## Cache values for marginal distributions in AGHQ over d-1.
		pTransformFix <- 0
		indexFix <- 0
    
		## Need to have eigen decomp saved.
		# hessEigenVecs <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		# hessEigenVals <- numeric(pTransform_length)
		# ATransform <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		# AInverse <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		covTheta <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		
		## We will want to cache the standard deviation skew terms.
    ## Default will be not to skew (e.g. 1)
		skewedStdDev <- matrix(1, nrow = pTransform_length, ncol = 2)
    ## Points for Asymmetric Gaussian Interpolation (integration free...?)
		extraPoints <- c(-15.0, -10.0, -7.0, 7.0, 10.0, 15.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0)
		aghqPoints <- pracma::gaussHermite(51)$x*sqrt(2)
		zMargGrid <- sort(unique(c(extraPoints, aghqPoints)))
    nzMargGrid <- length(zMargGrid)
    postModeFound <- 0

    ## Some cached values for summary statistics and reporting:
    marg_P <- matrix(0, nrow=nzMargGrid, ncol = npar)
    marg_theta <- array(0, c(pTransform_length, nzMargGrid, 2) )
    post_sims <- matrix(0, nrow=3, ncol = nre)  ## Just initialize it.

		## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    ## Default calculation method for AGHQuad
    methodID <- 2
    ## The nimbleList definitions AGHQuad_params and AGHQuad_summary
    ## have moved to predefined nimbleLists.


  },
