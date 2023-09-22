## R Function to build an AGHQ grid in multiple dimensions.
Rget_AGHQ_nodes <- function(dm = 1, n = 3) {
	## If n=1 or is miswritten, do Laplace.
	if(n <= 1){	
		z <- c(0,-1)	## For one time fixes
		w_star <- c(sqrt(2*pi), -1)
	}else{
		## Gauss-Hermite Nodes
		gh_nodes <- pracma::gaussHermite(n)
		## Convert to z = sqrt(2)*z, and wgt = exp(x^2)*sqrt(2)
		z <- sqrt(2) * gh_nodes$x
		w_star <- exp(gh_nodes$x^2 + log(gh_nodes$w)) * sqrt(2)
		if(n %% 2 != 0) z[ceiling(n/2)] <- 0	## Make it a hard zero.
	}
	if(dm == 1) return(list("z_nodes" = z, "weights" = w_star))
	## Expand those nodes for n dimensions.
	w <- rep(1, n^dm)
	z_nodes <- matrix(0, ncol = dm, nrow = n^dm)
	for (i in 1:dm) {
		z_nodes[ ,i] = rep(z, each = n^(i-1), times = (n^dm/n^i))
		w <- w * rep(w_star, each = n^(i-1), times = (n^dm/n^i))
	}	
	return(list("z_nodes" = as.matrix(z_nodes), "weights" = as.numeric(w)))
}

## Build a CCD grid for hyperparameters in the quadrature object.
CCDGrid <- function(dm)
{
	f0 <- 1.1	## Saw 1.1 on INLA and 1.5 default for mgcv. < sqrt(2) seems to be an issue but they use a larger sphere.

	z <- rsm::ccd(dm, n0 = c(1,rep(0,dm-1)), alpha = "rotatable", 
		inscribed = TRUE, oneblock=TRUE, randomize = FALSE)
	z <-  rsm::decode.data(z)[,-(1:2)]
	
	np <- nrow(z)
	wgts <- 1 / ((np - 1) * (1 + exp(-0.5*dm*f0^2) * (f0^2 - 1)))
	wgt0 <- 1 - (np-1)*wgts
			
	weights <- rep(wgts, np)
	weights[which(rowSums(abs(z)) == 0)] <- wgt0
	
	return(list("z_nodes" = as.matrix(z), "weights" = as.numeric(weights)))
}

GRID_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
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
		resetGrid = function(zNew = double(2), wgtNew = double(1)){},
		skewGridPoints = function(skewSD = double(2)){}
	)
)


## Function to store values from inner optimization to be used for 
## joint inference on the fixed and random effects.
## Regardless of AGHQ or CCD, we need to store all this information to be
## able to simulate for \bm{w}. Additional convenience functions added.

## Added in dimensions here so it can be updated to different grid sizes.
## Use an index to pass which is which. Added one time fixes if they want to
## do Laplace or some sort of Empirical Bayes method.
storeGridValues <- nimbleFunction(
	contains = GRID_BASE,
	setup = function(z, wgt, nre){
		gridFix <- 0
		if(!is.matrix(z)){
			nTheta <- length(z)
			nGrid <- 1
			gridFix <- 1
		}else{
			nTheta <- ncol(z)
			nGrid <- nrow(z)
		}
		if(nGrid == 1) {
			z <- matrix(z, ncol = nTheta, nrow = nGrid)
			wgt <- as.numeric(c(wgt, -1))
		}
		one_time_fixes_done <- FALSE

		thetaVals <- matrix(0, nrow = nGrid, ncol = nTheta)
		reMode <- matrix(0, nrow = nGrid, ncol = nre)
		calculated <- numeric(nGrid + gridFix)
		cholVals <- array(0, c(nGrid + gridFix, nre, nre))
		logDensTheta <- numeric(nGrid + gridFix)
		modeIndex <- 1
		if(nTheta == 1 & nGrid > 1) modeIndex <- which(z == 0)
		if(nTheta > 1 & nGrid > 1) modeIndex <- which(rowSums(abs(z)) == 0)
	},
	run=function(){},
	methods = list(
    one_time_fixes = function() {
      ## Run this once after compiling; remove extraneous -1 if necessary
      if(one_time_fixes_done) return()
      if(nGrid == 1) {
				calculated <<- numeric(length = 1, value = calculated[1])
				logDensTheta <<- numeric(length = 1, value = logDensTheta[1])
				wgt <<- numeric(length = 1, value = wgt[1])
			}
			one_time_fixes_done <<- TRUE
    },	
		## Reset the sizes of the storage to change the grid if the user wants more AGHQ.
		resetGrid = function(zNew = double(2), wgtNew = double(1)){
			one_time_fixes()
			nGrid <<- dim(zNew)[1]
			setSize(thetaVals, c(nGrid, nTheta))
			setSize(reMode, c(nGrid, nre))
			setSize(calculated, nGrid)
			setSize(cholVals, c(nGrid, nre, nre))
			setSize(logDensTheta, nGrid)
			setSize(z, c(nGrid, nTheta))
			setSize(wgt, nGrid)
			
			z <<- zNew
			wgt <<- wgtNew
			## Keep mode information, otherwise reset.
			for( i in 1:nGrid ){
				if( abs(sum(z[i,])) == 0 ) modeIndex <<- i
				calculated[i] <<- 0
			}
		},
		saveValues = function(i= integer(0, default = 1), theta = double(1), logDensity = double(),
									innerCholesky = double(2), innerMode = double(1)){
			one_time_fixes()
			if(i == 0) i <- modeIndex
			logDensTheta[i] <<- logDensity
			cholVals[i,,] <<- innerCholesky
			thetaVals[i,] <<- theta
			reMode[i,] <<- innerMode
			calculated[i] <<- 1
		},
		skewGridPoints = function(skewSD = double(2)){
			for( i in 1:nGrid )
			{
				for( j in 1:nTheta ){
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
		getLogDensity = function(i=integer()){returnType(double()); return(logDensTheta[i])}
	)
)

## USE FOR TESTING
# gridPts <- CCDGrid(dm = 2)
# zGrid <- as.matrix(gridPts$z_nodes)
# zWeights <- gridPts$weights
# test <- storeGridValues(z=0,1, nre = 5)
# cTest <- compileNimble(test)
# gridPts <- Rget_AGHQ_nodes(2, 5)
# chol <- gridPts
# cTest$resetGrid(gridPts$z, gridPts$weights)

# innerOpt_nfl <- nimbleFunctionList(GRID_BASE)
# for( i in 1:10 ) innerOpt_nfl[[i]] <- storeGridValues(zGrid[i,], zWeights[i], nre = 5)
# innerOpt_nfl[[1]]$saveValues(theta = c(0,0,0), cholesky = matrix(rnorm(25), 5, 5), inner_mode = rnorm(5))
# innerOpt_nfl[[1]]$getMode()