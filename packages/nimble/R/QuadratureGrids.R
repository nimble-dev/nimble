## Build a Grid Base to access all these functions in a nimble function list.
GRID_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
		buildGrid = function(){
    },
    quadSum = function(){
      returnType(double())
    },
		saveLogDens = function(i = integer(0), logDensity = double()){
    },
    transformGrid = function(skewSD = double(2)){
    },
		getWeights = function(i=integer()){
      returnType(double())
    },
    getAllWeights = function(){
      returnType(double(1))
    },
		getNodes = function(i=integer()){
      returnType(double(1))
    },
    getNodesTransformed = function(i=integer()){
      returnType(double(1))
    },
		getLogDensity = function(i=integer()){
      returnType(double())
    },
		setGridSize = function(nQUpdate = integer()){
    },
    skewGridPoints = function(skewSD = double(2)){
    },
		getModeIndex = function(){
      returnType(integer())
    },
		getGridSize = function(){
      returnType(integer())
    },
    buildGrid = function(){
    }
	)
)

## Write a basic quad rule:
quadRule_AGHQ = nimbleFunction(
  setup = function(d = 1, nQuad = 1){
		
    odd <- TRUE
    if(nQuad %% 2 == 0) odd <- FALSE
      
		if(nQuad > 35) {
			print("We don't currently support more than 35 quadrature nodes per dimension. Setting nQuad to 35")
      nQuad <- 35
		}
 
    one_time_fixes_done <- FALSE
    nQ <- nQuad^d
    wgt <- numeric(nQ)
    if(nQ == 1) wgt <- c(0, -1)
    zVals <- matrix(0, nrow = nQ, ncol = d)
    modeIndex <- -1
    ## Need to do a reverse for Eigen Vectors:
    inner_max <- 121
		reverse <- inner_max:1    
  },
  run = function(){},
  methods = list(
    one_time_fixes = function(){
			if(one_time_fixes_done) return()
			if(nQ == 1) {
				wgt <<- numeric(length = 1, value = wgt[1])
			}
			one_time_fixes_done <<- TRUE
    },
    buildAGHQOne = function(nQ1 = integer()){
      res <- matrix(0, nrow = nQ1, ncol = 2)
      if( nQ1 == 1 ){
        ## Laplace Approximation:
        res[,1] <- 0
        res[,2] <- sqrt(2*pi)
      }else{
        i <- 1:(nQ1-1)
        dv <- sqrt(i/2)
        ## Recreate pracma::Diag for this problem.
        y <- matrix(0, nrow = nQ1, ncol = nQ1)
        y[1:(nQ1-1), 1:(nQ1-1) + 1] <- diag(dv)
        y[1:(nQ1-1) + 1, 1:(nQ1-1)] <- diag(dv)
        E <- eigen(y, symmetric = TRUE)
        L <- E$values	# Always biggest to smallest.
        V <- E$vectors
        inds <- reverse[(inner_max-nQ1+1):inner_max]	## Hard coded to maximum 120.
        x <- L[inds]
        ## Make mode hard zero. We know nQ is odd and > 1.
        if(odd) x[ceiling(nQ1 / 2 ) ] <- 0
        V <- t(V[, inds])
        ## Update nodes and weights in terms of z = x/sqrt(2) 
        ## and include Gaussian kernel in weight to integrate an arbitrary function.
        w <- V[, 1]^2  * sqrt(2*pi) * exp(x^2)
        x <- sqrt(2) * x
        res[,1] <- x
        res[,2] <- w
      }
      returnType(double(2))
      return(res)
    },
    buildRule = function(nQuadUpdate = integer(0, default = 0)){ 
      one_time_fixes()
      
      if(nQuadUpdate > 0){
        nQuad <<- nQuadUpdate
        nQ <<- nQuad^d
        zVals <<- matrix(0, nrow = nQ, ncol = d)
      }

      odd <<- TRUE
      if(nQuad %% 2 == 0) odd <<- FALSE
      
      if( nQuad == 1 ){
        ## Laplace Approximation:
        zVals <<- matrix(0, nrow = nQuad, ncol = d)
        wgt <<- numeric(value = exp(0.5 * d * log(2*pi)), length = nQ)
        modeIndex <<- 1
      }else{
        nodes <- buildAGHQOne(nQuad)
        ## If d = 1, then we are done.
        if(d == 1){
          zVals[,1] <<- nodes[,1]
          wgt <<- nodes[,2]
          if(odd) modeIndex <<- which(zVals[,1] == 0)[1]
        }else{
          ## Build the multivariate quadrature rule.
          wgt <<- rep(1, nQ)
          
          ## A counter for when to swap.
          swp <- numeric(value = 0, length = d)
          for( ii in 1:d ) swp[ii] <- nQuad^(ii-1)

          ## Repeat x for each dimension swp times.
          for(j in 1:d ) {
            indx <- 1
            for( ii in 1:nQ ) {
              zVals[ii, j] <<- nodes[indx,1]
              wgt[ii] <<- wgt[ii]*nodes[indx,2]
              k <- ii %% swp[j] 
              if(k == 0) indx <- indx + 1
              if(indx > nQuad) indx <- 1
            }
          }
          ## Assuming mode index is the middle number.
          if(odd) {
            modeIndex <<- ceiling(nQ/2)
            ## Just in case that goes horribly wrong...
            if(sum(abs(zVals[modeIndex,])) != 0) {
              for(ii in 1:nQ) {
                if(sum(abs(zVals[ii,])) == 0) modeIndex <<- ii
              }
            }
          }
        }
      }
    },
    getNodes = function(){
      returnType(double(2))
      return(zVals)
    },
    getWeights = function(){
      returnType(double(1))
      return(wgt)
    }
  )
)

quadRule_CCD <- nimbleFunction(
	setup = function(d = 1){
    if ((d > 120 | d < 1)) stop("Dimension of Theta must be in [1,120]")	
    
		## Walsh Index Assignments for Resolution V Fractional Factorials
		index <- c(1, 2, 4, 8, 15, 16, 32, 51, 64, 85, 106, 128,
			150, 171, 219, 237, 247, 256, 279, 297, 455, 512, 537,
			557, 594, 643, 803, 863, 998, 1024, 1051, 1070, 1112,
			1169, 1333, 1345, 1620, 1866, 2048, 2076, 2085, 2185,
			2372, 2456, 2618, 2800, 2873, 3127, 3284, 3483, 3557,
			3763, 4096, 4125, 4135, 4174, 4435, 4459, 4469, 4497,
			4752, 5255, 5732, 5804, 5915, 6100, 6369, 6907, 7069,
			8192, 8263, 8351, 8422, 8458, 8571, 8750, 8858, 9124,
			9314, 9500, 10026, 10455, 10556, 11778, 11885, 11984,
			13548, 14007, 14514, 14965, 15125, 15554, 16384, 16457,
			16517, 16609, 16771, 16853, 17022, 17453, 17891, 18073,
			18562, 18980, 19030, 19932, 20075, 20745, 21544, 22633,
			23200, 24167, 25700, 26360, 26591, 26776, 28443, 28905,
			29577, 32705)
			
		## Number of grid points for different dimensions of theta.
		nCCD <- index; p <- 1
		for (i in 1:length(index)) {
			if (index[i]>=p) p <- p * 2
			nCCD[i] <- p
		}
		nC <- nCCD[d] ## minimum 2. Note that if d = 1, INLA does a grid approximation instead of CCD. Should do same here.
    nQ <- nC + 2*d + 1
		
    zVals <- matrix(0, nrow = nQ, ncol = d)

		## One time fixes for scalar / vector changes.
		wgt <- numeric(nQ)

		## CCD mode index is 1 always.
		modeIndex <- 1
  },
	run=function(){},
	methods = list(
    ## Taken from Simon Wood's mgcv package.
		## https://github.com/cran/mgcv/blob/master/R/inla.r
		## However, we do scaled design following INLA such that z*zT = 1
		## from https://github.com/hrue/r-inla/blob/devel/gmrflib/design.c
		buildRule = function(){
			## First point is mode.
			design <- matrix(0, nQ, d)
			for (i in 1:d) {
				design[index[i]+2,i] <- 1
				design[2:(nC+1),i] <- fwt(x = design[2:(nC+1),i], n = nC)
			}
			design <- design/sqrt(d)
			
			## Next are the star points on the axes. (scaled)
			design[(nC+2):(nC + d + 1), 1:d] <- diag(d)*1
			design[(nC + d + 2):(nC + 2*d + 1), 1:d] <- diag(d)*-1

			## Weights as defined by Rue 2009. 
			## Note that the paper weights are incorrect: https://groups.google.com/g/r-inla-discussion-group/c/sy2xYin7YJA
			f0 <- 1.1
			wgts <- 1 / ((nQ - 1 ) * ( 1 + exp(- (d * f0^2)/2) * (f0^2 - 1 )) ) 
			wgt0 <- 1 - (nQ-1)*wgts
			
			zVals <<- design
      wgt[1] <<- wgt0
			wgt[2:nQ] <<- rep(wgts, nQ-1)
		},
		## fast Walsh transform taken from Wood MGCV inla.
		fwt = function(x = double(1), n = integer()) {
			lag <- 1
			while (lag < n) {
			offset <-  lag * 2
			ngroups <- length(x)/offset
				for (group in 0:(ngroups-1)) { ## vectorized
					j <- 1:lag + group*offset
					k <- j + lag
					xj <- x[j]; xk <- x[k]
					x[j] <- xj + xk
					x[k] <- xj - xk
				}
			lag <- offset
			} ## while lag
			returnType(double(1))
			return(x)
		},
    getNodes = function(){
      returnType(double(2))
      return(zVals)
    },
    getWeights = function(){
      returnType(double(1))
      return(wgt)
    }
  )
)


## Method for summing likelihoods on real scale with possible small values.
## Returns back on log scale.
#' @export
logSumExp = nimbleFunction(
  run = function(log1 = double(), log2 = double()){
  if(log1 > log2) ans <- log( 1 + exp(log2 - log1)) + log1
  else ans <- log(1 + exp(log1 - log2)) + log2
  returnType(double())
  return(ans)
})

#' @export
buildAGHQGrid <- nimbleFunction(
  contains = GRID_BASE,
	setup = function(d = 1, nQuad = 3){
    
		odd <- TRUE
    if(nQuad %% 2 == 0) odd <- FALSE
      
		if(nQuad > 35) {
			print("We don't currently support more than 35 quadrature nodes per dimension. Setting nQuad to 35")
      nQuad <- 35
		}

		## nQ will be total number of quadrature nodes.
		nQ <- nQuad^d	## Maybe dimension reduced if we prune.
		zVals <- matrix(0, nrow = nQ, ncol = d)
    nodeVals <- matrix(0, nrow = nQ, ncol = d)

		## Need to do a reverse for Eigen Vectors:
    inner_max <- 121
		reverse <- inner_max:1

		## One time fixes if we run into some scalar issues for compilation.
		## This is exclusively if the user requests Laplace (nQ = 1 AGHQ).
		gridFix <- 0
		if(nQ == 1)	{
			gridFix <- 1
		}

		## One time fixes for scalar / vector changes.
		one_time_fixes_done <- FALSE
		wgt <- numeric(nQ + gridFix)
		logDensity <- numeric(nQ + gridFix)
 		logdetNegHessian <- 0
    margDens <- 0
    
    gridBuilt <- FALSE

    ## Transformation methods
    transformation <- "spectral"
    transformationSet <- FALSE
    calcEigen <- FALSE
    ATransform <- matrix(0, d, d)
    AInverseTransform <- matrix(0, d, d)
    cholNegHess_saved <- matrix(0, d, d)

		## AGHQ mode will be in the middle.
		modeIndex <- -1
    maxLogDensity <- 0
  },
	run=function(){},
	methods = list(
		one_time_fixes = function() {
			## Run this once after compiling; remove extraneous -1 if necessary
			if(one_time_fixes_done) return()
			if(nQ == 1) {
				logDensity <<- numeric(length = 1, value = logDensity[1])
				wgt <<- numeric(length = 1, value = wgt[1])
			}
			one_time_fixes_done <<- TRUE
		},	
		buildAGHQOne = function(nQ1 = integer()){
			res <- matrix(0, nrow = nQ1, ncol = 2)
			if( nQ1 == 1 ){
				## Laplace Approximation:
				res[,1] <- 0
				res[,2] <- sqrt(2*pi)
			}else{
				i <- 1:(nQ1-1)
				dv <- sqrt(i/2)
				## Recreate pracma::Diag for this problem.
				y <- matrix(0, nrow = nQ1, ncol = nQ1)
				y[1:(nQ1-1), 1:(nQ1-1) + 1] <- diag(dv)
				y[1:(nQ1-1) + 1, 1:(nQ1-1)] <- diag(dv)
				E <- eigen(y, symmetric = TRUE)
				L <- E$values	# Always biggest to smallest.
				V <- E$vectors
				inds <- reverse[(inner_max-nQ1+1):inner_max]	## Hard coded to maximum 120.
				x <- L[inds]
				## Make mode hard zero. We know nQ is odd and > 1.
				if(odd) x[ceiling(nQ1 / 2 ) ] <- 0
				V <- t(V[, inds])
				## Update nodes and weights in terms of z = x/sqrt(2) 
				## and include Gaussian kernel in weight to integrate an arbitrary function.
				w <- V[, 1]^2  * sqrt(2*pi) * exp(x^2)
				x <- sqrt(2) * x
				res[,1] <- x
				res[,2] <- w
			}
			returnType(double(2))
			return(res)
		},
		buildAGHQ = function(){
			one_time_fixes()
			if( nQuad == 1 ){
				## Laplace Approximation:
				zVals <<- matrix(0, nrow = 1, ncol = d)
				wgt <<- numeric(value = exp(0.5 * d * log(2*pi)), length = nQ)
				modeIndex <<- 1
			}else{
        nodes <- buildAGHQOne(nQuad)
        ## If d = 1, then we are done.
        if(d == 1){
          zVals[,1] <<- nodes[,1]
          wgt <<- nodes[,2]
          if(odd) modeIndex <<- which(zVals[,1] == 0)[1]
        }else{
          ## Build the multivariate quadrature rule.
          wgt <<- rep(1, nQ)
          
          ## A counter for when to swap.
          swp <- numeric(value = 0, length = d)
          for( ii in 1:d ) swp[ii] <- nQuad^(ii-1)

          ## Repeat x for each dimension swp times.
          for(j in 1:d ) {
            indx <- 1
            for( ii in 1:nQ ) {
              zVals[ii, j] <<- nodes[indx,1]
              wgt[ii] <<- wgt[ii]*nodes[indx,2]
              k <- ii %% swp[j] 
              if(k == 0) indx <- indx + 1
              if(indx > nQuad) indx <- 1
            }
          }
          ## Assuming mode index is the middle number.
          if(odd) {
            modeIndex <<- ceiling(nQ/2)
            ## Just in case that goes horribly wrong...
            if(sum(abs(zVals[modeIndex,])) != 0) {
              for(ii in 1:nQ) {
                if(sum(abs(zVals[ii,])) == 0) modeIndex <<- ii
              }
            }
          }
        }
      }
		},
		## Doesn't default to building the grid.
		buildGrid = function(){
			one_time_fixes()
			if(!gridBuilt) buildAGHQ()
      gridBuilt <<- TRUE
		},
    ## Should I update this to the logSumExp code from CCD grid?
    quadSum = function(){
      if(!odd) modeIndex <<- -1 ## Make sure it's negative.
      margDens <<- 0
      for( k in 1:nQ ){
        if(k == modeIndex) margDens <<- margDens + wgt[k]
        else margDens <<- margDens + exp(logDensity[k] - maxLogDensity)*wgt[k]
      }
      ans <- log(margDens) + maxLogDensity - 0.5 * logdetNegHessian
      returnType(double())
      return(ans)
    },
		## Reset the sizes of the storage to change the grid if the user wants more/less AGHQ.
		setGridSize = function(nQUpdate = integer()){
			one_time_fixes()
 
      if(nQuad != nQUpdate){
        nQ <<- nQUpdate^d
        nQuad <<- nQUpdate

        if(nQ %% 2 == 0) {
          odd <<- FALSE
          modeIndex <<- -1
        }else{
          odd <<- TRUE
        }
        ## Update weights and nodes.
        setSize(wgt, nQ)
        setSize(zVals, c(nQ, d))
        setSize(nodeVals, c(nQ, d))
        setSize(logDensity, nQ)

        ## Build the new grid (updates modeIndex).
        buildAGHQ()
      }
    },
		saveLogDens = function(i = integer(0, default = -1), logDens = double()){
      if(i == -1){
        if(odd) logDensity[modeIndex] <<- logDens
        maxLogDensity <<- logDens
			}else{
        logDensity[i] <<- logDens
      }
    },
    transformGrid1D = function(negHess = double(2), inner_mode = double(1)){
      SD <- 1/sqrt(negHess[1,1])
      for( i in 1:nQ) nodeVals[i,] <<- inner_mode + SD*zVals[i,]
      logdetNegHessian <<- log(negHess[1,1])
    }, 
    ## Set the transformation method to be used and compute Eigen Decomp if reqiured.
    setTransformation = function(cholNegHess = double(2), inner_mode = double(1), method = character()){
      cholNegHess_saved <<- cholNegHess
      logdetNegHessian <<- 2*sum(log(diag(cholNegHess)))
      if(method == "spectral"){
        transformation <<- "spectral"      
        negHess <- t(cholNegHess) %*% cholNegHess
        eigenDecomp <- nimEigen(negHess)
        for( i in 1:d ){
            ATransform[,i] <<- eigenDecomp$vectors[,i]/sqrt(eigenDecomp$values[i]) # eigenDecomp$vectors %*% diag(1/sqrt(eigenDecomp$values))
            AInverseTransform[,i] <<- eigenDecomp$vectors[,i]*sqrt(eigenDecomp$values[i]) # eigenDecomp$vectors %*% diag(1/sqrt(eigenDecomp$values))
        }
        calcEigen <<- TRUE
      }else{
        transformation <<- "cholesky"
      }
      transformationSet <<- TRUE
    },
    ## Transform the quadrature grid:
    ## Can change transformation type here.
    transformGrid = function(cholNegHess = double(2), inner_mode = double(1), method = character()){
      if(!transformationSet | (!calcEigen & method == "spectral") ){
        setTransformation(cholNegHess, inner_mode, method)
      }
      if(method == "spectral"){
        transformation <<- "spectral"
        for( i in 1:nQ) nodeVals[i, ] <<- inner_mode + (ATransform %*% zVals[i,])
      }else{
        transformation <<- "cholesky"
        for( i in 1:nQ) nodeVals[i, ] <<- inner_mode + backsolve(cholNegHess, zVals[i,])
      }
    },
    z_to_theta = function(z = double(1)){
      returnType(double(1))
      if(transformation == "spectral"){
        theta <- (nodeVals[modeIndex,] + (ATransform %*% z))[,1]
      }else{
        ## Cholesky transformation.
        theta <- nodeVals[modeIndex,] + backsolve(cholNegHess_saved, z)
      }
      return(theta)
    },
    theta_to_z = function(theta = double(1)){
      returnType(double(1))
      if(transformation == "spectral"){
        z <- ((theta - nodeVals[modeIndex,]) %*% AInverseTransform)[1,]
      }else{
        ## Cholesky transformation.
        z <- (cholNegHess_saved %*% (theta - nodeVals[modeIndex,]))[,1]
      }
      return(z)
    },    
    ## Do nonthing function:
    skewGridPoints = function(skewSD = double(2)){
    },
		getWeights = function(i=integer()){
      returnType(double())
      if(i == -1 & odd)  return(wgt[modeIndex])
      return(wgt[i])
    },
		getAllWeights = function(){
      returnType(double(1))
      return(wgt)
    },    
		getNodesTransformed = function(i=integer()){
      if(i == -1 & odd) return(nodeVals[modeIndex,])
      returnType(double(1)); 
      return(nodeVals[i,])
    },
		getAllNodesTransformed = function(){
      returnType(double(2)); 
      return(nodeVals)
    },
		getNodes = function(i=integer()){
      if(i == -1 & odd) return(zVals[modeIndex,])
      returnType(double(1)); 
      return(zVals[i,])
    },
    getAllNodes = function(){
      returnType(double(2)); 
      return(zVals)
    },    
		getLogDensity = function(i=integer()){
      returnType(double())
      if(i == -1 & odd) return(logDensity[modeIndex])
      return(logDensity[i])
    },
    getModeIndex = function(){
      returnType(integer())
      return(modeIndex)
    },
    getGridSize = function(){
      returnType(double())
      return(nQ)
    }
	)
)## End of buildAGHQGrid

#' Build Adaptive Gauss-Hermite Quadrature Grid
#'
#' Create quadrature grid for use in AGHQuad methods in Nimble.
#'
#' @param d Dimension of quadrature grid being requested.
#'
#' @param nQuad Number of quadrature nodes requested on build.
#'
#' @name buildAGHQGrid
#' 
#' @details
#'
#' This function is used by used by \code{buildOneAGHQuad1D}
#' and \code{buildOneAGHQuad} create the quadrature grid using
#' adaptive Gauss-Hermite quadrature. Handles single or multiple dimension 
#' grids and computes both grid locations and weights. Additionally, acts
#' as a cache system to do transformations, and return marginalized log density.
#'
#' Any of the input node vectors, when provided, will be processed using
#'   \code{nodes <- model$expandNodeNames(nodes)}, where \code{nodes} may be
#'   \code{paramNodes}, \code{randomEffectsNodes}, and so on. This step allows
#'   any of the inputs to include node-name-like syntax that might contain
#'   multiple nodes. For example, \code{paramNodes = 'beta[1:10]'} can be
#'   provided if there are actually 10 scalar parameters, 'beta[1]' through
#'   'beta[10]'. The actual node names in the model will be determined by the
#'   \code{exapndNodeNames} step.
#'
#' Available methods include
#' 
#' \itemize{
#'
#'   \item \code{buildAGHQ}. Builds a adaptive Gauss-Hermite quadrature grid in d dimensions.
#'   Calls \code{buildAGHQOne} to build the one dimensional grid and then expands in each dimension.
#'   Some numerical issues occur in Eigen decomposition making the grid weights only accurate up to 
#'   35 quadrature nodes.
#'
#'   \item Options to get internally cached values are \code{getGridSize},
#'   \code{getModeIndex} for when there are an odd number of quadrature nodes,
#'   \code{getLogDensity} for the cached values, \code{getAllNodes} for the 
#'   quadrature grids, \code{getNodes} for getting a single indexed nodes,
#'   \code{getAllNodesTransformed} for nodes transformed to the parameter scale,
#'   \code{getNodesTransformed} for a single transformed node, \code{getAllWeights} 
#'   to get all quadrature weights, \code{getWeights} single indexed weight.
#'
#'   \item \code{transformGrid(cholNegHess, inner_mode, method)} transforms 
#'   the grid using either cholesky trasnformations,
#'   as default, or spectral that makes use of the Eigen decomposition. For a single
#'   dimension \code{transformGrid1D} is used.
#'
#'   \item As the log density is evaluated externally, it is saved via \code{saveLogDens},
#'   which then is summed via \code{quadSum}.
#'
#'   \item \code{buildGrid} builds the grid the initial time and is only run once in code. After,
#'   the user must choose to \code{setGridSize} to update the grid size.
#'
#'
#'   \item \code{check}. If TRUE (default), a warning is issued if
#'         \code{paramNodes}, \code{randomEffectsNodes} and/or \code{calcNodes}
#'         are provided but seek to have missing elements or unnecessary
#'         elements based on some default inspection of the model. If
#'         unnecessary warnings are emitted, simply set \code{check=FALSE}.
#'
#'   \item \code{innerOptimControl}. A list of control parameters for the inner 
#'         optimization of Laplace approximation using \code{optim}. See 
#'         'Details' of \code{\link{optim}} for further information.
#'
#'   \item \code{innerOptimMethod}. Optimization method to be used in 
#'         \code{optim} for the inner optimization. See 'Details' of 
#'         \code{\link{optim}}. Currently \code{optim} in NIMBLE supports: 
#'         "\code{Nelder-Mead}", "\code{BFGS}", "\code{CG}", and 
#'         "\code{L-BFGS-B}". By default, method "\code{CG}" is used when 
#'         marginalizing over a single (scalar) random effect, and "\code{BFGS}" 
#'         is used for multiple random effects being jointly marginalized over.
#'
#'   \item \code{innerOptimStart}. Choice of starting values for the inner 
#'         optimization. This could be \code{"last"}, \code{"last.best"}, or a 
#'         vector of user provided values. \code{"last"} means the most recent 
#'         random effects values left in the model will be used. When finding 
#'         the MLE, the most recent values will be the result of the most recent 
#'         inner optimization for Laplace. \code{"last.best"} means the random 
#'         effects values corresponding to the largest Laplace likelihood (from 
#'         any call to the \code{calcLaplace} or \code{calcLogLik} method, 
#'         including during an MLE search) will be used (even if it was not the 
#'         most recent Laplace likelihood). By default, the initial random 
#'         effects values will be used for inner optimization.
#'
#'   \item \code{outOptimControl}. A list of control parameters for maximizing
#'         the Laplace log-likelihood using \code{optim}. See 'Details' of
#'         \code{\link{optim}} for further information.
#' }
#'
#' @references
#'
#' Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss Quadrature Rules. 
#' Mathematics of Computation 23 (106): 221-230.
#'
#' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature. Biometrika, 81(3) 624-629.
#'
#' Jackel, P. (2005). A note on multivariate Gauss-Hermite quadrature. London: ABN-Amro. Re.
#'
NULL



#' @export
buildCCDGrid <- nimbleFunction(
  contains = GRID_BASE,
	setup = function(d = 1, nQuad = 3){
    if ((d > 120 | d < 1)) stop("Dimension of Theta must be in [1,120]")	
    
		## Walsh Index Assignments for Resolution V Fractional Factorials
		index <- c(1, 2, 4, 8, 15, 16, 32, 51, 64, 85, 106, 128,
			150, 171, 219, 237, 247, 256, 279, 297, 455, 512, 537,
			557, 594, 643, 803, 863, 998, 1024, 1051, 1070, 1112,
			1169, 1333, 1345, 1620, 1866, 2048, 2076, 2085, 2185,
			2372, 2456, 2618, 2800, 2873, 3127, 3284, 3483, 3557,
			3763, 4096, 4125, 4135, 4174, 4435, 4459, 4469, 4497,
			4752, 5255, 5732, 5804, 5915, 6100, 6369, 6907, 7069,
			8192, 8263, 8351, 8422, 8458, 8571, 8750, 8858, 9124,
			9314, 9500, 10026, 10455, 10556, 11778, 11885, 11984,
			13548, 14007, 14514, 14965, 15125, 15554, 16384, 16457,
			16517, 16609, 16771, 16853, 17022, 17453, 17891, 18073,
			18562, 18980, 19030, 19932, 20075, 20745, 21544, 22633,
			23200, 24167, 25700, 26360, 26591, 26776, 28443, 28905,
			29577, 32705)
			
		## Number of grid points for different dimensions of theta.
		nCCD <- index; p <- 1
		for (i in 1:length(index)) {
			if (index[i]>=p) p <- p * 2
			nCCD[i] <- p
		}
		nC <- nCCD[d] ## minimum 2. Note that if d = 1, INLA does a grid approximation instead of CCD. Should do same here.
    skewed <- 0
    nQ <- nC + 2*d + 1
    
    ## CCD mode is always at 1
    modeIndex <- 1
    optimInfo <- 0
		
    zVals <- matrix(0, nrow = nQ, ncol = d)
    nodeVals <- matrix(0, nrow = nQ, ncol = d)

		## One time fixes for scalar / vector changes.
		wgt <- numeric(nQ)
		logDensity <- numeric(nQ)
 		logdetNegHessian <- 0
    
    gridBuilt <- FALSE
   
		## CCD mode index is 1 always.
		modeIndex <- 1
    maxLogDensity <- 0
    
    ## Transformation methods
    transformation <- "spectral"
    transformationSet <- FALSE
    calcEigen <- FALSE
    ATransform <- matrix(0, d, d)
    AInverseTransform <- matrix(0, d, d)
    cholNegHess_saved <- matrix(0, d, d)
  },
	run=function(){},
	methods = list(
    ## Taken from Simon Wood's mgcv package.
		## https://github.com/cran/mgcv/blob/master/R/inla.r
		## However, we do scaled design following INLA such that z*zT = 1
		## from https://github.com/hrue/r-inla/blob/devel/gmrflib/design.c
		buildCCD = function(){
			## First point is mode.
			design <- matrix(0, nQ, d)
			for (i in 1:d) {
				design[index[i]+2,i] <- 1
				design[2:(nC+1),i] <- fwt(x = design[2:(nC+1),i], n = nC)
			}
			design <- design/sqrt(d)
			
			## Next are the star points on the axes. (scaled)
			design[(nC+2):(nC + d + 1), 1:d] <- diag(d)*1
			design[(nC + d + 2):(nC + 2*d + 1), 1:d] <- diag(d)*-1

			## Weights as defined by Rue 2009. 
			## Note that the paper weights are incorrect: https://groups.google.com/g/r-inla-discussion-group/c/sy2xYin7YJA
			f0 <- 1.1
			wgts <- 1 / ((nQ - 1 ) * ( 1 + exp(- (d * f0^2)/2) * (f0^2 - 1 )) ) 
			wgt0 <- 1 - (nQ-1)*wgts
			
			zVals <<- design
			wgt <<- c(wgt0, rep(wgts, nQ-1))
		},
		## fast Walsh transform taken from Wood MGCV inla.
		fwt = function(x = double(1), n = integer()) {
			lag <- 1
			while (lag < n) {
			offset <-  lag * 2
			ngroups <- length(x)/offset
				for (group in 0:(ngroups-1)) { ## vectorized
					j <- 1:lag + group*offset
					k <- j + lag
					xj <- x[j]; xk <- x[k]
					x[j] <- xj + xk
					x[k] <- xj - xk
				}
			lag <- offset
			} ## while lag
			# returnType(double(1))
			return(x)
		},  
		## Doesn't default to building the grid.
		buildGrid = function(){
			if(!gridBuilt) buildCCD()
      gridBuilt <<- TRUE
		},
    quadSum = function(){
      ans <- 0
      for( k in 1:nQ ){
        ans <- logSumExp(ans, logDensity[k] + log(wgt[k]))
      }
      ans <- ans - 0.5 * logdetNegHessian
      returnType(double())
      return(ans)
    },
		## AGHQ specific do nothing here.
		setGridSize = function(nQUpdate = integer()){
    },
		saveLogDens = function(i = integer(0, default = -1), logDens = double()){
      if(i == -1){
        if(odd) logDensity[modeIndex] <<- logDens
        maxLogDensity <<- logDens
			}else{
        logDensity[i] <<- logDens
      }
    },
    transformGrid1D = function(negHess = double(2), inner_mode = double(1)){
      SD <- 1/sqrt(negHess[1,1])
      for( i in 1:nQ) nodeVals[i,] <<- inner_mode + SD*zVals[i,]
      logdetNegHessian <<- log(negHess[1,1])
    },
    ## Set the transformation method to be used and compute Eigen Decomp if reqiured.
    setTransformation = function(cholNegHess = double(2), inner_mode = double(1), method = character()){
      cholNegHess_saved <<- cholNegHess
      logdetNegHessian <<- 2*sum(log(diag(cholNegHess)))
      if(method == "spectral"){
        transformation <<- "spectral"
        negHess <- t(cholNegHess) %*% cholNegHess
        eigenDecomp <- nimEigen(negHess)
        for( i in 1:d ){
            ATransform[,i] <<- eigenDecomp$vectors[,i]/sqrt(eigenDecomp$values[i]) # eigenDecomp$vectors %*% diag(1/sqrt(eigenDecomp$values))
            AInverseTransform[,i] <<- eigenDecomp$vectors[,i]*sqrt(eigenDecomp$values[i]) # eigenDecomp$vectors %*% diag(1/sqrt(eigenDecomp$values))
        }
        calcEigen <<- TRUE
      }else{
        transformation <<- "cholesky"
      }
      transformationSet <<- TRUE
    },
    ## Transform the quadrature grid:
    ## Can change transformation type here.
    transformGrid = function(cholNegHess = double(2), inner_mode = double(1), method = character()){
      if(!transformationSet | (!calcEigen & method == "spectral") ){
        setTransformation(cholNegHess, inner_mode, method)
      }
      if(method == "spectral"){
        transformation <<- "spectral"
        for( i in 1:nQ) nodeVals[i, ] <<- inner_mode + (ATransform %*% zVals[i,])
      }else{
        transformation <<- "cholesky"
        for( i in 1:nQ) nodeVals[i, ] <<- inner_mode + backsolve(cholNegHess, zVals[i,])
      }
    },
    z_to_theta = function(z = double(1)){
      returnType(double(1))
      if(transformation == "spectral"){
        theta <- (nodeVals[modeIndex,] + (ATransform %*% z))[,1]
      }else{
        ## Cholesky transformation.
        theta <- nodeVals[modeIndex,] + backsolve(cholNegHess_saved, z)
      }
      return(theta)
    },
    theta_to_z = function(theta = double(1)){
      returnType(double(1))
      if(transformation == "spectral"){
        z <- ((theta - nodeVals[modeIndex,]) %*% AInverseTransform)[1,]
      }else{
        ## Cholesky transformation.
        z <- (cholNegHess_saved %*% (theta - nodeVals[modeIndex,]))[,1]
      }
      return(z)
    },
    ## Skew the CCD grid according to the +/- skewed std normal on each side of mode.
    ## Matches with INLA code base. f = 1 as far as I can tell.
    ## z_local[i] = f * design->experiment[k][i]
    ##		    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
    skewGridPoints = function(skewSD = double(2)){
      if(skewed == 0){
        for( i in 1:nQ )
        {
        for( j in 1:d ){
          sdAdj <- skewSD[j, 2]	# Positive Skew
          if(zVals[i,j] <= 0) sdAdj <- skewSD[j, 1]	# Negative Skew
            zVals[i, j] <<- zVals[i, j] * sdAdj
          }
        }
        skewed <<- 1  
      }
    },
		getWeights = function(i=integer()){
      returnType(double())
      if(i == -1 & odd)  return(wgt[modeIndex])
      return(wgt[i])
    },
		getAllWeights = function(){
      returnType(double(1))
      return(wgt)
    },
		getNodesTransformed = function(i=integer()){
      if(i == -1) return(nodeVals[modeIndex,])
      returnType(double(1)); 
      return(nodeVals[i,])
    },
		getAllNodesTransformed = function(){
      returnType(double(2)); 
      return(nodeVals)
    },
		getNodes = function(i=integer()){
      if(i == -1) return(zVals[modeIndex,])
      returnType(double(1)); 
      return(zVals[i,])
    },
    getAllNodes = function(){
      returnType(double(2)); 
      return(zVals)
    },    
		getLogDensity = function(i=integer()){
      returnType(double())
      if(i == -1 & odd) return(logDensity[modeIndex])
      return(logDensity[i])
    },
    getModeIndex = function(){
      returnType(integer())
      return(modeIndex)
    },
    getGridSize = function(){
      returnType(double())
      return(nQ)
    }
	)
)## End of buildCCDGrid


## Create a caching random effects system for simulating the posterior random effect distribution according to Stringer:
## This requires the inner mode, the inner cholesky, the 

## Things I need to add for approx posterior
## save inner mode, save inner cholesky etc. for doing simulation of random-effects
## save outer mode and negHessian. 
## Need wgt*density
## and inner mode
## and inner cholesky for each point:
grid_inner_cache = nimbleFunction(
  setup = function(nre = 0, nQuad = 0){
    innerMode <- matrix(0, nrow = nQuad, ncol = nre)
    innerNegHess <- array(0, c(nQuad, nre, nre))
    if(nQuad < 1) wgts <- c(-1,-1)
    else wgts <- numeric(nQuad)
    
    one_time_fixes_done <- FALSE
  },
  run = function(){},
  methods = list(
    one_time_fixes = function(){
      if(!one_time_fixes_done)
        wgts <<- numeric(value = 0, length = nQuad)
    },
    ## Note to self, this wgt will be density*wgt, strictly for simulating.
    cache_wgts = function(weight = double(), indx = integer()){
      one_time_fixes()
      wgts[indx] <<- weight
    },
    cache_inner_mode = function(mode = double(1), indx = integer()){
      innerMode[indx,] <<- mode
    },
    cache_inner_negHess = function(negHess = double(2), indx = integer()){
      innerNegHess[indx,,] <<- negHess
    },
    update_nQuad = function(nQuadUpdate = integer()){
      nQuad <<- nQuadUpdate
      wgts <<- numeric(value = 0, length = nQuad)
      innerMode <<- matrix(0, nrow = nQuad, ncol = nre)
      innerNegHess <<- array(0, c(nQuad, nre, nre))
    },
    simulate = function(n = integer()){
      val <- matrix(0, nrow = n, ncol = nre)
      simwgt <- wgts/sum(wgts)
      for( i in 1:n ){
        k <- rcat(1, prob = simwgt)
        val[i,] <- rmnorm_chol(n=1, mean = innerMode[k,],  
                                cholesky = innerNegHess[k,,], prec_param = TRUE)
      }
      returnType(double(2))
      return(val)
    }
  )
)
