## Added in dimensions here so it can be updated to different grid sizes.
## Use an index to pass which is which. Added one time fixes if they want to
## do Laplace or some sort of Empirical Bayes method.
# nQ - Number of quadrature points for AGHQ or custom
# method - Grid method to use, either CCD, AGHQ, or Custom (which will require the user to input the grid in future iterations)
# nre - Number of fixed and random effects.
# require(nimble)

## Do we need references? 
## References
## Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss Quadrature Rules. Mathematics of Computation 23 (106): 221-230.
## Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature. Biometrika, 81(3) 624-629.
library(nimble)
Rget_AGHQ_nodes <- function(n = double()) {
	## If n=1 or is miswritten, do Laplace.
	if(n <= 1){
		z <- c(0,-1)	## For one time fixes
		w_star <- c(sqrt(2*pi), -1)
	}else{
		## Gauss-Hermite Nodes
		gh_nodes <- pracma::gaussHermite(n)
		## Convert to z = sqrt(2)*z, and wgt = exp(x^2)*sqrt(2)
		z <- as.numeric(sqrt(2) * gh_nodes$x)
		w_star <- as.numeric(exp(gh_nodes$x^2 + log(gh_nodes$w)) * sqrt(2))
	}
	return(cbind(z, w_star))
}


quadRule <- nimbleRcall(function(n = double()){},
  Rfun = "Rget_AGHQ_nodes",
  returnType = double(2)
)

buildAGHQGrid <- nimbleFunction(
	# contains = GRID_BASE,
	setup = function(d, nQuad){
    
		odd <- TRUE
    if(nQuad %% 2 == 0) odd <- FALSE
      
		if(nQuad > 35) {
			print("We don't currently support more than 35 quadrature nodes per dimension. Setting nQuad to 35")
      nQuad <- 35
		}

		## nQ will be total number of quadrature points.
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
    transformGrid = function(cholNegHess = double(2), inner_mode = double(1), method = character()){
      if(method == "spectral"){
        ## Spectral transformation.
        negHess <- t(cholNegHess) %*% cholNegHess
        eigenDecomp <- nimEigen(negHess)
        ATransform <- matrix(0, nrow = d, ncol = d)   
        for( i in 1:d ){
            ATransform[,i] <- eigenDecomp$vectors[,i]/sqrt(eigenDecomp$values[i]) # eigenDecomp$vectors %*% diag(1/sqrt(eigenDecomp$values))
        }
        for( i in 1:nQ) nodeVals[i, ] <<- inner_mode + (ATransform %*% zVals[i,])
      }else{
        ## Cholesky transformation.
        for( i in 1:nQ) nodeVals[i, ] <<- inner_mode + backsolve(cholNegHess, zVals[i,])
      }
      logdetNegHessian <<- 2*sum(log(diag(cholNegHess)))
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
)
