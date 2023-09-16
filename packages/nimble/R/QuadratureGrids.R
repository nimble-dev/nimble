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
	return(list("z_nodes" = z_nodes, "weights" = w))
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
	
	return(list("z_nodes" = z, "weights" = weights))
}

GRID_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
		saveValues = function(theta = double(1), 
									cholesky = double(2), inner_mode = double(1)){},
		calcCheck = function(){returnType(integer())},
		getCholesky = function(){returnType(double(2))},
		getMode = function(){returnType(double(1))},
		getTheta = function(){returnType(double(1))},
		updateWeights = function(weights = double()){},
		getZWeights = function(){returnType(double())},
		getZNodes = function(){returnType(double(1))},
		updateZNodes = function(zNodes = double(1)){}	
	)
)


## Function to store values from inner optimization to be used for 
## joint inference on the fixed and random effects.
## Regardless of AGHQ or CCD, we need to store all this information to be
## able to simualte for \bm{w}. Additional convenience functions added.
storeGridValues <- nimbleFunction(
	contains = GRID_BASE,
	setup = function(zVal, zWgt, nre){
		z <- zVal
		wgt <- zWgt
		nTheta <- length(z)
		thetaVals <- numeric(length = nTheta)
		reMode <- numeric(length = nre)
		calculated <- 0
		cholVals <- matrix(0, ncol = nre, nrow = nre)
	},
	run=function(){},
	methods = list(
		saveValues = function(theta = double(1), 
									cholesky = double(2), inner_mode = double(1)){
			cholVals <<- cholesky
			thetaVals <<- theta
			reMode <<- inner_mode
			calculated <<- 1
		},
		calcCheck = function(){returnType(integer()); return(calculated)},
		getCholesky = function(){
			returnType(double(2))
			return(cholVals)
		},
		getMode = function(){
			returnType(double(1))
			return(reMode)
		},
		getTheta = function(){
			returnType(double(1))
			return(thetaVals)	
		},
		updateWeights = function(weights = double()){wgt <<- weights},
		getZWeights = function(){returnType(double()); return(wgt)},
		getZNodes = function(){returnType(double(1)); return(z)},
		updateZNodes = function(zNodes = double(1)){z <<- zNodes}
	)
)	

# zGrid <- as.matrix(gridPts$z_nodes)
# zWeights <- gridPts$weights
# test <- storeGridValues(zGrid[1,], zWeights[1], nre = 50)
# cTest <- compileNimble(test)
# innerOpt_nfl <- nimbleFunctionList(GRID_BASE)
# for( i in 1:10 ) innerOpt_nfl[[i]] <- storeGridValues(zGrid[i,], zWeights[i], nre = 5)
# innerOpt_nfl[[1]]$saveValues(theta = c(0,0,0), cholesky = matrix(rnorm(25), 5, 5), inner_mode = rnorm(5))
# innerOpt_nfl[[1]]$getMode()