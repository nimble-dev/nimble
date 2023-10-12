## Build a Grid Base to access all these functions in a nimble function list.
GRID_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
		buildAGHQ = function(){},
		buildCCD = function(){},
		buildGrid = function(){},
		saveValues = function(i = integer(0),theta = double(1), logDensity = double(), 
														innerCholesky = double(2), innerMode = double(1)){},
		calcCheck = function(i=integer()){returnType(integer())},
		getCholesky = function(i=integer()){returnType(double(2))},
		getMode = function(i=integer()){returnType(double(1))},
		getTheta = function(i=integer()){returnType(double(1))},
		getWeights = function(i=integer()){returnType(double())},
		updateWeights = function(i=integer(), weights = double()){},
		getNodes = function(i=integer()){returnType(double(1))},
		updateNodes = function(i=integer(), zNodes = double(1)){},
		getLogDensity = function(i=integer()){returnType(double())},
		resetGrid = function(nQUpdate = integer()){},
		skewGridPoints = function(skewSD = double(2)){},
		getThetaModeIndex = function(){returnType(integer())},
		getGridSize = function(){returnType(integer())},
		buildAGHQOne = function(nQ1 = integer()){returnType(double(2))}
	)
)

## Added in dimensions here so it can be updated to different grid sizes.
## Use an index to pass which is which. Added one time fixes if they want to
## do Laplace or some sort of Empirical Bayes method.
# nQ - Number of quadrature points for AGHQ or custom
# method - Grid method to use, either CCD, AGHQ, or Custom (which will require the user to input the grid in future iterations)
# nre - Number of fixed and random effects.
buildQuadGrid <- nimbleFunction(
	contains = GRID_BASE,
	setup = function(d, nQ, method, nre){
		# if(!method %in% c("CCD", "AGHQ", "Custom")) print("Only CCD, AGHQ, and Custom Grids are supported")
		if (d > 120 | d<1) print("Dimension of Theta must be in [1,120]")		

		## Don't let them use an even quadrature grid. Inefficient.
		if(method == 2 & nQ %% 2 == 0) {
			print("Even Quadrature grids are not recommended: Changing to + 1")
			nQ <- nQ + 1	
		}
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
		nC <- nCCD[d]

		## Need to do a reverse for Eigen Vectors:
		reverse <- 120:1

		## nQ will be total number of quadrature points.
		nQ_aghq <- nQ
		if( method == 1 ) nQ <- nC + 2*d + 1
		if( method == 2 ) {
			nQ <- nQ_aghq^d	## May be dimension reduced if we prune.
		}
		z <- matrix(0, nrow = nQ, ncol = d)

		## One time fixes if we run into some scalar issues for compilation.
		## This is exclusively if the user requests Laplace (nQ = 1 AGHQ).
		gridFix <- 0
		if(nQ == 1)	{
			gridFix <- 1
		}

		thetaVals <- matrix(0, nrow = nQ, ncol = d)
		reMode <- matrix(0, nrow = nQ, ncol = nre)
		cholVals <- array(0, c(nQ, nre, nre))
		
		## One time fixes for scalar / vector changes.
		one_time_fixes_done <- FALSE
		wgt <- numeric(d + gridFix)
		calculated <- numeric(nQ + gridFix)
		logDensTheta <- numeric(nQ + gridFix)
		
		## CCD mode index will be 1, but AGHQ it will be in the middle. Won't allow even nQ.
		modeIndex <- 1
	},
	run=function(){},
	methods = list(
		one_time_fixes = function() {
			## Run this once after compiling; remove extraneous -1 if necessary
			if(one_time_fixes_done) return()
			if(nQ == 1) {
				calculated <<- numeric(length = 1, value = calculated[1])
				logDensTheta <<- numeric(length = 1, value = logDensTheta[1])
				wgt <<- numeric(length = 1, value = wgt[1])
			}
			one_time_fixes_done <<- TRUE
		},	
		## Taken from Simon Wood's mgcv package.
		## https://github.com/cran/mgcv/blob/master/R/inla.r
		## However, we do scaled design following INLA such that z*zT = 1
		## from https://github.com/hrue/r-inla/blob/devel/gmrflib/design.c
		buildCCD = function(){
			one_time_fixes()
			## First point is mode.
			nQ <<- nC + 2*d + 1
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
			## Note that the radius is scaled so we don't need exp(-0.5*d*f0^2)
			f0 <- 1.1 # sqrt(d) + 0.1
			wgts <- 1 / ((nQ - 1) * (f0^2 - 1) * (1 + exp(-0.5*f0^2) ))
			wgt0 <- 1 - (nQ-1)*wgts
			
			z <<- design
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
			returnType(double(1))
			return(x)
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
				inds <- reverse[(120-nQ1+1):120]	## Hard coded to maximum d.
				x <- L[inds]
				## Make mode hard zero. We know nQ is odd and > 1.
				x[ceiling(nQ1 / 2 ) ] <- 0
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
			if( nQ_aghq == 1 ){
				## Laplace Approximation:
				z <<- matrix(0, nrow = 1, ncol = d)
				wgt <<- numeric(value = sqrt(2*pi), length = nQ)
				modeIndex <<- 1
				}else{
					nodes <- buildAGHQOne(nQ_aghq)
					## If d = 1, then we are done.
					if(d == 1) {
						z[,1] <<- nodes[,1]
						wgt <<- nodes[,2]
					}else {
						## Build the multivariate quadrature rule.
						wgt <<- rep(1, nQ)
						
						## A counter for when to swap.
						swp <- numeric(value = 0, length = d)
						for( ii in 1:d ) swp[ii] <- nQ_aghq^(ii-1)

						## Repeat x for each dimension swp times.
						for(j in 1:d ) {
							indx <- 1
							for( ii in 1:nQ )
							{
								z[ii, j] <<- nodes[indx,1]
								wgt[ii] <<- wgt[ii]*nodes[indx,2]
								k <- ii %% swp[j] 
								if(k == 0) indx <- indx + 1
								if(indx > nQ_aghq) indx <- 1
							}
							if(sum(abs(z[ii,])) == 0) modeIndex <<- ii
						}
					}
				}
		},
		## Doesn't default to building the grid.
		buildGrid = function(){
			if(method == 1) buildCCD()
			if(method == 2) buildAGHQ()
		},
		## Reset the sizes of the storage to change the grid if the user wants more/less AGHQ.
		resetGrid = function(nQUpdate = integer()){
			one_time_fixes()
			## If calling this, then method must be AGHQ currently. Can't 'update ccd'.
			method <<- 2
			nQ_aghq <<- nQUpdate
			nQ <<- nQUpdate^d
			
			## Update weights and nodes.
			setSize(wgt, nQ)
			setSize(z, c(nQ, d))			
			buildAGHQ()
			
			## This will resize the storage to be saved later.
			setSize(thetaVals, c(nQ, d))
			setSize(reMode, c(nQ, nre))
			setSize(calculated, nQ)
			setSize(cholVals, c(nQ, nre, nre))
			setSize(logDensTheta, nQ)
			
			## Keep mode information, otherwise reset.
			## I assume this all automatically lines up.
			for( i in 1:nQ ){
				if( abs(sum(z[i,])) == 0 ) modeIndex <<- i
				if( i != modeIndex ) calculated[i] <<- 0
			}
		},
		saveValues = function(i = integer(0, default = 1), theta = double(1), logDensity = double(),
									innerCholesky = double(2), innerMode = double(1)){
			one_time_fixes()
			if(i == 0) i <- modeIndex	## Can quickly update mode with input of zero.
			logDensTheta[i] <<- logDensity
			cholVals[i,,] <<- innerCholesky
			thetaVals[i,] <<- theta
			reMode[i,] <<- innerMode
			calculated[i] <<- 1
		},
		## Skew the CCD grid according to the +/- skewed std normal on each side of mode.
		## Matches with INLA code base. f = 1 as far as I can tell.
		## z_local[i] = f * design->experiment[k][i]
		##		    * (design->experiment[k][i] > 0.0 ? stdev_corr_pos[i] : stdev_corr_neg[i]);
		skewGridPoints = function(skewSD = double(2)){
			for( i in 1:nQ )
			{
				for( j in 1:d ){
					sdAdj <- skewSD[j, 2]	# Positive Skew
					if(z[i,j] <= 0) sdAdj <- skewSD[j, 1]	# Negative Skew
					z[i, j] <<- z[i, j] * sdAdj
				}
			}
		},
		calcCheck = function(i=integer()){returnType(integer()); return(calculated[i])},
		getCholesky = function(i=integer()){returnType(double(2)); return(cholVals[i,,])},
		getMode = function(i=integer()){returnType(double(1)); return(reMode[i,])},
		getTheta = function(i=integer()){returnType(double(1)); return(thetaVals[i,])},
		getWeights = function(i=integer()){returnType(double()); return(wgt[i])},
		updateWeights = function(i=integer(), weight = double()){wgt[i] <<- weight},
		getNodes = function(i=integer()){returnType(double(1)); return(z[i,])},
		updateNodes = function(i=integer(), zNodes = double(1)){z[i,] <<- zNodes},
		getLogDensity = function(i=integer()){returnType(double()); return(logDensTheta[i])},
		getThetaModeIndex = function(){returnType(integer()); return(modeIndex)},
		getGridSize = function(){returnType(integer()); return(nQ)}
	)
)

## Paul's notes for testing:
## 1) Make sure quadrature in 1+2 d integrates common functions correctly.
## 2) Compare weights and nodes to pracma and to mvQuad.

## USE FOR TESTING
# gridPts <- nimCCD(dm = 2)
# zGrid <- as.matrix(gridPts$z_nodes)
# zWeights <- gridPts$weights
# test <- storeGridValues(z=0,1, nre = 5)
# cTest <- compileNimble(test)
# gridPts <- Rget_AGHQ_nodes(2, 5)
# chol <- gridPts
# cTest$resetGrid(gridPts$z, gridPts$weights)

## Testing to make sure these match my existing code.
# innerOpt_nfl <- nimbleFunctionList(GRID_BASE)
# innerOpt_nfl[[1]] <- buildQuadGrid(d = 2, nQ = 3, method = 1, nre = 5)
# tmp <- buildQuadGrid(d = 3, nQ = 5, method = 2, nre = 5)
# tmp$buildGrid()
# tmp$getNodes()
# tmp$getWeights()
# tmp$resetGrid(nQUpdate = 5)
# tmp$buildAGHQOne(51)
p
# np <- pracma::gaussHermite(5)
# mvq <- mvQuad::createNIGrid(dim=3, type = "GHe", level = 5)
# z=mvQuad::getNodes(mvq)
# w=mvQuad::getWeights(mvq)

# tmpc <- compileNimble(tmp)
# innerOpt_nfl[[1]]$buildGrid()
# plot(innerOpt_nfl[[1]]$getNodes())
# innerOpt_nfl[[2]] <- buildQuadGrid(d = 2, nQ = 3, method = 2, nre = 5)
# innerOpt_nfl[[2]]$buildGrid()
# vals <- Rget_AGHQ_nodes(2, 3)
# innerOpt_nfl[[2]]$getNodes()
# innerOpt_nfl[[1]]$saveValues(theta = c(0,0,0), cholesky = matrix(rnorm(25), 5, 5), inner_mode = rnorm(5))
# innerOpt_nfl[[1]]$getMode()