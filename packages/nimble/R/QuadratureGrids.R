## Added in dimensions here so it can be updated to different grid sizes.
## Use an index to pass which is which. Added one time fixes if they want to
## do Laplace or some sort of Empirical Bayes method.
# nQ - Number of quadrature points for AGHQ or custom
# method - Grid method to use, either CCD, AGHQ, or Custom (which will require the user to input the grid in future iterations)
# nre - Number of fixed and random effects.
buildAGHQGrid <- nimbleFunction(
	# contains = GRID_BASE,
	setup = function(d, nQuad){
    
		## Don't let them use an even quadrature grid. Inefficient.
		if(nQuad %% 2 == 0) {
			print("Even Quadrature grids are not recommended")
		}

		## nQ will be total number of quadrature points.
		nQ <- nQuad^d	## Maybe dimension reduced if we prune.
		zVals <- matrix(0, nrow = nQ, ncol = d)
    nodeVals <- matrix(0, nrow = nQ, ncol = d)
		## Need to do a reverse for Eigen Vectors:
		reverse <- nQ:1

		## One time fixes if we run into some scalar issues for compilation.
		## This is exclusively if the user requests Laplace (nQ = 1 AGHQ).
		gridFix <- 0
		if(nQ == 1)	{
			gridFix <- 1
		}

		## One time fixes for scalar / vector changes.
		one_time_fixes_done <- FALSE
		wgt <- numeric(d + gridFix)
		logDensity <- numeric(nQ + gridFix)
		logGradient <- numeric(nQ + gridFix)
 		logdetNegHessian <- 0
   
		## AGHQ mode will be in the middle.
		modeIndex <- 1
  },
	run=function(){},
	methods = list(
		one_time_fixes = function() {
			## Run this once after compiling; remove extraneous -1 if necessary
			if(one_time_fixes_done) return()
			if(nQ == 1) {
				logDensity <<- numeric(length = 1, value = logDensity[1])
        logGradient <<- numeric(length = 1, value = logDensity[1])
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
				inds <- reverse[(nQ-nQ1+1):nQ]
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
			if( nQuad == 1 ){
				## Laplace Approximation:
				zVals <<- matrix(0, nrow = 1, ncol = d)
				wgt <<- numeric(value = sqrt(2*pi), length = nQ)
				modeIndex <<- 1
			}else{
        nodes <- buildAGHQOne(nQuad)
        ## If d = 1, then we are done.
        if(d == 1){
          zVals[,1] <<- nodes[,1]
          wgt <<- nodes[,2]
          modeIndex <<- which(zVals[,1] == 0)[1]
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
          modeIndex <<- ceiling(nQ/2)
          ## Just in case that goes horribly wrong...
          if(sum(abs(zVals[modeIndex,])) != 0) {
            for(ii in 1:nQ) {
              if(sum(abs(zVals[ii,])) == 0) modeIndex <<- ii
            }
          }
        }
      }
		},
		## Doesn't default to building the grid.
		buildGrid = function(){
			one_time_fixes()
			buildAGHQ()
		},
    quadSum = function(){
      margProb <- 0
      maxVal <- logDensity[modeIndex]
      for( k in 1:nQ ){
        if(k == modeIndex) margProb <- margProb + wgt[k]
        else margProb <- margProb + exp(logDensity[k] - maxVal)*wgt[k]
      }
      ans <- log(margProb) + maxVal - 0.5 * logdetNegHessian
      returnType(double())
      return(ans)
    },
		## Reset the sizes of the storage to change the grid if the user wants more/less AGHQ.
		resetGrid = function(nQUpdate = integer()){
			one_time_fixes()
      if(nQuad != nQUpdate){
        nQ <<- nQUpdate^d
        
        ## Update weights and nodes.
        setSize(wgt, nQ)
        setSize(zVals, c(nQ, d))
        setSize(nodeVals, c(nQ, d))
        setSize(logDensity, nQ)
        setSize(logGradient, nQ)

        ## Build the new grid (updates modeIndex).
        buildAGHQ()
      }
    },
		saveLogDens = function(i = integer(0, default = 0), logDens = double()){
      if(i == 0) {
        logDensity[modeIndex] <<- logDens
			}else{
        logDensity[i] <<- logDens
      }
    },
		saveLogGrad = function(i = integer(0, default = 0), logGrad = double()){
      if(i == 0) {
         logGradient[modeIndex] <<- logGrad
			}else{
         logGradient[i] <<- logGrad
      }
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
      if(i == 0)  return(wgt[modeIndex])
      return(wgt[i])
    },
		getNodesTransformed = function(i=integer()){
      if(i == 0) return(nodeVals[modeIndex,])
      returnType(double(1)); 
      return(nodeVals[i,])
    },
		getAllNodesTransformed = function(){
      returnType(double(2)); 
      return(nodeVals)
    },
		getNodes = function(i=integer()){
      if(i == 0) return(zVals[modeIndex,])
      returnType(double(1)); 
      return(zVals[i,])
    },
		getLogDensity = function(i=integer()){
      returnType(double())
      if(i == 0) return(logDensity[modeIndex])
      return(logDensity[i])
    }
	)
)

# test <- buildAGHQGrid(d = 2, nQuad = 5)
# testc <- compileNimble(test)
# testc$buildGrid()
# testc$getNodes(i=0)

# testc$transformGrid(cholNegHess, method = "cholesky", inner_mode = c(2,1))
# grd <- testc$getAllNodesTransformed()

# grd2 <- matrix(0, 5^2,2)
# for( i in 1:(5^2)) grd2[i,] <- testc$getNodesTransformed(i)
# plot(grd)
# points(grd2, col = 'red', pch = 4)
# z <- NULL
# for( i in 1:(2^5)) z <- rbind(z, testc$getNodes(i))
# plot(grd)
# plot(z)
# x = seq(0, 5, length = 100)
# y = seq(0, 4, length = 100)
# f <- function(x, y){ apply(cbind(x,y), 1, 
  # FUN = function(z){ dmnorm_chol(z, mean = c(2,1), cholesky = cholNegHess, prec=TRUE)} )}
# fz <- outer(x, y, f)
# contour(x,y,fz)
# points(grd)