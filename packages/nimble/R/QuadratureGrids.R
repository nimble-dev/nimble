## nimbleQuad Quadrature Rules + Grids - Rules are choosing between AGHQ and CCD etc.
## Grids are implementations of the rule to generate the actual nodes and weights.
QUAD_RULE_BASE <- nimbleFunctionVirtual(
  name = 'QUAD_RULE_BASE',
  run = function() {},
  methods = list(
    buildGrid = function(nQuad = integer(0, default = 0),  d = integer(0, default = 1)){
      returnType(quadGridListDef())
    }
  )
)


#' Nimble List Quadrature Data type
#'
#' Creates a quadrature nimble list type to be used internally and for making new custom
#' quadrature rules to marginalize random effects and for posterior approximations.  
#'
#' @details
#' 
#'
#' List is generated with three data types. An integer that is the mode index `modeIndex` that indicates
#' which quadrature node is the mode, values are all zero. A numeric vector, `wgts`, that is a weight for each
#' quadrature node. A matrix, `nodes`, that are the quadrature nodes made by the rule that are of dimension `nQ` rows
#' and `d` dimension columns.
#'
#' @author Paul van Dam-Bates
#' @export
quadGridListDef <- nimbleList(modeIndex = integer(0), 
                              wgts = double(1), 
                              nodes = double(2),
                              name = "quadGridList")
                              
#' @export
# quadGridListDef <- nimbleList(
  # list(
    # nimbleType('modeIndex','integer', 0),
    # nimbleType('wgts', 'double', 1),
    # nimbleType('nodes', 'double', 2)
  # ),
  # name = "quadGridListDef",
  # predefined = TRUE
# )

## Stand alone simple 1D AGHQ function:
## Note this is for convenience due to
## Needing this for sparse grids too.
## AGHQ for multivariate as we do now is the "product rule" version.
AGHQ1D <- nimbleFunction(
  run = function(nQuad = integer(0, default = 1)){
      odd <- TRUE
      if(nQuad %% 2 == 0) 
        odd <- FALSE

      res <- matrix(0, nrow = nQuad, ncol = 2)
      if( nQuad == 1 ){
        ## Laplace Approximation:
        res[,2] <- 0
        res[,1] <- sqrt(2*pi)
      }else{
        i <- 1:(nQuad-1)
        dv <- sqrt(i/2)
        ## Recreate pracma::Diag for this problem.        
        if(nQuad == 2)
          fill_diag <- matrix(dv,1,1)
        else 
          fill_diag <- diag(dv)

        y <- matrix(0, nrow = nQuad, ncol = nQuad)
        y[1:(nQuad-1), 1:(nQuad-1) + 1] <- fill_diag
        y[1:(nQuad-1) + 1, 1:(nQuad-1)] <- fill_diag
        E <- eigen(y, symmetric = TRUE)
        L <- E$values	# Always biggest to smallest.
        V <- E$vectors
        inds <- numeric(value = 0, length = nQuad)
        for( j in seq_along(L) ) inds[j] <- nQuad-j+1 ## Is this an efficient way to do it?
        x <- L[inds]
        ## Make mode hard zero. We know nQ is odd and > 1.
        if(odd) x[ceiling(nQuad / 2 ) ] <- 0
        V <- t(V[, inds])
        ## Update nodes and weights in terms of z = x/sqrt(2) 
        ## and include Gaussian kernel in weight to integrate an arbitrary function.
        w <- V[, 1]^2  * sqrt(2*pi) * exp(x^2)
        x <- sqrt(2) * x
        res[,1] <- w
        res[,2] <- x
      }
      returnType(double(2))
      return(res)    
  }
)

## Write a basic quad rule:
## Maybe not cache here at all?
## This will end up being in a different wrapper.
## This avoids generating too much memory as this function gets called.
quadRule_AGHQ = nimbleFunction(
  contains = QUAD_RULE_BASE,
  name = 'quadRule_AGHQ',  
  setup = function(){},
  run = function(){},
  methods = list(
    buildGrid = function(nQuad = integer(0, default = 0), d = integer(0, default = 1)){
      if(nQuad > 35) {
        print("Warning:  More than 35 quadrature nodes per dimension is not supported. Setting nQuad to 35.")
        nQuad <- 35
      }
      if(nQuad == 0) {
        print("Warning:  No default number of quadrature points given. Assuming nQuad = 3 per dimension.")
        nQuad <- 3
      }
      odd <- TRUE
      if(nQuad %% 2 == 0) 
        odd <- FALSE

      nQ <- nQuad^d
      zVals <- matrix(0, nrow = nQ, ncol = d)
      wgt <- numeric(value = 0, length = nQ)

      if( nQuad == 1 ){
        ## Laplace Approximation:
        wgt <- numeric(value = exp(0.5 * d * log(2*pi)), length = nQ)
        modeIndex <- 1
      }else{
        nodes <- AGHQ1D(nQuad)
        ## If d = 1, then we are done.
        if(d == 1){
          zVals[,1] <- nodes[,2]
          wgt <- nodes[,1]
          if(odd) modeIndex <- which(zVals[,1] == 0)[1]
        }else{
          ## Build the multivariate quadrature rule.
          wgt <- rep(1, nQ)
          
          ## A counter for when to swap.
          swp <- numeric(value = 0, length = d)
          for( ii in 1:d ) swp[ii] <- nQuad^(ii-1)

          ## Repeat x for each dimension swp times.
          for(j in 1:d ) {
            indx <- 1
            for( ii in 1:nQ ) {
              zVals[ii, j] <- nodes[indx,2]
              wgt[ii] <- wgt[ii]*nodes[indx,1]
              k <- ii %% swp[j] 
              if(k == 0) indx <- indx + 1
              if(indx > nQuad) indx <- 1
            }
          }
          ## Assuming mode index is the middle number.
          if(odd) {
            modeIndex <- ceiling(nQ/2)
            ## Just in case that goes horribly wrong...
            if(sum(abs(zVals[modeIndex,])) != 0) {
              for(ii in 1:nQ) {
                if(sum(abs(zVals[ii,])) == 0) modeIndex <- ii
              }
            }
          }
        }
        if(!odd)
          modeIndex <- -1  ## No mode is present.
      }
      returnType(quadGridListDef())
      output <- quadGridListDef$new()
      output$modeIndex <- as.integer(modeIndex)
      output$wgts <- wgt
      output$nodes <- zVals
      return(output)
    }
  )
)

## CCD Grid quadrature from Rue et al 2009, adapted based on some code from MGCV
## for their approximate posterior methods.
quadRule_CCD <- nimbleFunction(
  contains = QUAD_RULE_BASE,
  name = 'quadRule_CCD',  
	setup = function(){    
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
			
  },
	run=function(){},
	methods = list(
    ## Taken from Simon Wood's mgcv package.
    ## https://github.com/cran/mgcv/blob/master/R/inla.r
    ## However, we do scaled design following INLA such that z*zT = 1
    ## from https://github.com/hrue/r-inla/blob/devel/gmrflib/design.c
    ## Can't update nQuad here but makes it general.
    buildGrid = function(nQuad = integer(0, default = 0), d = integer(0, default = 1)){ 
      if ((d > 120 | d < 1)) stop("Dimension of Theta must be in [1,120]")	

      ## Number of grid points for different dimensions of theta.
      nCCD <- index; p <- 1
      for (i in seq_along(index)) {
        if (index[i]>=p) p <- p * 2
        nCCD[i] <- p
      }
      nC <- nCCD[d] ## minimum 2. If 1, choose points c(0,-1,1) but they don't make sense.
      nQ <- nC + 2*d + 1

      ## First point is mode.,
      design <- matrix(0, nQ, d)
      
      if(d > 1){
        for (i in 1:d) {
          design[index[i]+2,i] <- 1
          design[2:(nC+1),i] <- fwt(x = design[2:(nC+1),i], n = nC)
        }
        design <- design/sqrt(d)
        ## Next are the star points on the axes. (scaled)
        design[(nC+2):(nC + d + 1), 1:d] <- diag(d)*1
        design[(nC + d + 2):(nC + 2*d + 1), 1:d] <- diag(d)*-1
      }else{
        design <- matrix(c(0,-1,1), nrow = 3, ncol = 1)
        nQ <- 3
      }

      ## Weights as defined by Rue 2009. 
      ## Note that the paper weights are incorrect: https://groups.google.com/g/r-inla-discussion-group/c/sy2xYin7YJA
      ## See https://github.com/hrue/r-inla/blob/devel/gmrflib/approx-inference.c#L1894
      # w = 1.0 / ((design->nexperiments - 1.0) * (1.0 + exp(-0.5 * SQR(f)) * (SQR(f) / nhyper - 1.0)));
      f0 <- 1.1
      ## From INLA: z_local[i] = f * design->experiment[k][i] where f = f0*sqrt(d)
      design <- design*sqrt(d)*f0
      
      ## Weights that actually make sense: 
      ## Including making the points at distance f0*sqrt(m) on the sphere:
      ## ***This part does not match INLA but the theory***
      wgts <- 1 / ((nQ - 1 ) * f0^2 * (2*pi)^(-d/2)*exp(-d*f0^2/2)) 
      wgt0 <- (2*pi)^(d/2)*(1 - f0^-2)
      ## INLA Weights
      # wgts <- 1 / ((nQ - 1 ) * ( 1 + exp(- (d * f0^2)/2) * (f0^2 - 1 )) ) 
      # wgt0 <- 1 - (nQ-1)*wgts

      ## One time fixes for scalar / vector changes.
      wgt <- numeric(value = 0, length = nQ)
      wgt[1] <- wgt0
      wgt[2:nQ] <- rep(wgts, nQ-1)

      returnType(quadGridListDef())
      output <- quadGridListDef$new()
      output$modeIndex <- 1L
      output$wgts <- wgt
      output$nodes <- design
      return(output)
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
    }
  )
)

permsR <- function(a){pracma::perms(a)}
nimPerms <- nimbleRcall(function(a = double(1)){}, "permsR", returnType = double(2))

quadRule_USER <- nimbleFunction(
  contains = QUAD_RULE_BASE,
  name = 'quadRule_Custom',  
	setup = function(){},
  run = function(){},
  methods = list(
    buildGrid = function(nQuad = integer(0, default = 0), d = integer(0, default = 1)){
      ## This will be a place holder for something others may choose to add.
      ## Can look for quadRule_Custom and check if it's implemented. If it is will try and use it...
      returnType(quadGridListDef())
      output <- quadGridListDef$new()
      output$modeIndex <- 1L
      output$wgts <- numeric(nQuad)
      output$nodes <- matrix(0, nrow=nQuad, d)
      return(output)
    }
  )
)

## Method for summing likelihoods on real scale with possible small values.
## Returns back on log scale.
#' @export
logSumExp = nimbleFunction(
  run = function(log1 = double(), log2 = double()){
  if(log1 > log2) 
    ans <- log( 1 + exp(log2 - log1)) + log1
  else 
    ans <- log(1 + exp(log1 - log2)) + log2
  returnType(double())
  return(ans)
  }, buildDerivs = list(run = list())
)

QUAD_CACHE_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
    cacheQuadGrid = function(nQuad = double(), nodes = double(2), wgts = double(1), modeIndex = integer()){},
    nodes = function(indx = integer(0, default = 0)){returnType(double(2))},
    weights = function(indx = integer(0, default = 0)){returnType(double(1))},
    modeI = function(){returnType(integer())},
    gridSize = function(){returnType(integer())},
    checkGrid = function(nQuad = double(0, default = -1), prune = double(0, default = 0)){returnType(logical())},
    pruneGrid = function(prune = double(0, default = 0)){}
  )
)

quadGridCache <- nimbleFunction(
  contains = QUAD_CACHE_BASE,
  setup = function(){
    quadGridList_internal <- quadGridListDef
    nodes_cached <- matrix(0, nrow = 1, ncol = 1)
    weights_cached <- c(0.0,0.0)
    modeIndex_cached <- -1
    nGrid_cached <- 0L
    gridBuilt <- FALSE
    prune_ <- 0
    nQuad_ <- -1
    numError <- 1e-10 ## ***CJP More precise?
    d <- 1
  },
  run = function(){},
  methods = list(
    cacheQuadGrid = function(nQuad = double(), nodes = double(2), wgts = double(1), modeIndex = integer()){
      nodes_cached <<- nodes
      weights_cached <<- wgts
      modeIndex_cached <<- modeIndex
      nGrid_cached <<- dim(nodes_cached)[1]
      gridBuilt <<- TRUE
      nQuad_ <<- nQuad
      prune_ <<- 0
      d <<- dim(nodes_cached)[2]
    },
    checkGrid = function(nQuad = double(0, default = -1), prune = double(0, default = 0)){
      returnType(logical())
      if(!gridBuilt | ((nQuad > 0) & (nQuad != nQuad_)) | ((prune_ != prune) & (prune_ > 0)))
        return(FALSE)
      else 
        return(TRUE)
    },
    ## Keep the biggest prune proportion weights.
    ## Note that our weights are for arbitrary functions.
    ## As a result, pruning will be weights adjust for a multivariate normal.
    pruneGrid = function(prune = double(0, default = 0)){
      if(prune_ == 0){
        if(!gridBuilt & prune > 0) {
          print("Warning: Cannot prune grid as the quadrature grid isn't built yet.")
        }else{
          ntrim <- 0
          ## Adjust weights:
          weights_adj <- numeric(value = 0, length = nGrid_cached)
          for(i in seq_along(weights_cached)){
            weights_adj[i] <- exp(sum(dnorm(nodes_cached[i,],mean=0,sd=1,log=TRUE)))*weights_cached[i]
          }
          ## Leave 3 points in total.
          while(ntrim/nGrid_cached < prune & ntrim < nGrid_cached - 3) {
            keep <- which(weights_adj > min(weights_adj) + numError)  ## error check as weights might be equal but off by numerical.
            if(dim(keep)[1] > 0){
              weights_adj <- weights_adj[keep]
              weights_cached <<- weights_cached[keep]
              nodes_cached <<- matrix(nodes_cached[keep,], nrow = length(keep), ncol = d)
              ntrim <- nGrid_cached - dim(nodes_cached)[1]
              ## Update mode index:
              if(modeIndex_cached > 0){
                modei <- which(keep == modeIndex_cached)
                if(dim(modei)[1] > 0)
                  modeIndex_cached <<- modei[1]
                else
                  modeIndex_cached <<- -1
              }
            }else{
              ## Exit loop.
              ntrim <- nGrid_cached
            }
          }
          nGrid_cached <<- dim(nodes_cached)[1]
        }
      }
      prune_ <<- prune
    },
    nodes = function(indx = integer(0, default = 0)){
      returnType( double(2) )
      if(indx > 0)
        return(matrix(nodes_cached[indx, ], nrow = 1))
      if(indx == -1 & modeIndex_cached > 0)
        return(matrix(nodes_cached[modeIndex_cached,], nrow = 1))
      return(nodes_cached)
    },
    weights = function(indx = integer(0, default = 0)){
      returnType( double(1) )
      if(indx > 0)
        return(numeric(value = weights_cached[indx], length = 1))
      if(indx == -1 & modeIndex_cached > 0)
        return(numeric(value = weights_cached[modeIndex_cached], length = 1))
      return(weights_cached)
    },
    modeI = function(){
      returnType(integer())
      return(modeIndex_cached)
    },
    gridSize = function(){
      returnType(integer())
      return(nGrid_cached)
    }
  )
)

## Wrapper to make quadrature nodes accesible in a nimble function list.
##***CJP check better naming convention on nQuad_.
#' @export
configureQuadGrid <- nimbleFunction(
  name = "quadGridClass",
  setup = function(d = 1, nQuad_ = 3, quadRule = "AGHQ", control = list()){
    ## Can list all possible quad rules here and set it.
    possibleRules <- c("AGHQ", "CCD", "USER")
    
    quadRules <- extractControlElement(control, "quadRules", quadRule)

    if(!any(quadRule == quadRules))
      quadRules <- c(quadRule, quadRules)

    if(!all(quadRules %in% possibleRules))
      stop("Error:  Only AGHQ or CCD or USER suplied rules are currently implemented.")      

    prune_ <- extractControlElement(control, "prune", 0)
    if(prune_ > 1 | prune_ < 0)
      stop("Can only prune a proportion of quadrature points.")
    
    quadGridList_internal <- quadGridListDef
    quadGridCache_nfl <- nimbleFunctionList(QUAD_CACHE_BASE)
    quadRule_nfl <- nimbleFunctionList(QUAD_RULE_BASE)
    
    I_AGHQ <- I_CCD <- I_USER <- 1
    I_RULE <- which(quadRules == quadRule)[1]

    ## Can I loop through these more efficiently?
    ## I have different names for each function so probably not...
    for( i in seq_along(quadRules) ){
      if(quadRules[i] == "AGHQ") {
        I_AGHQ <- i
        quadRule_nfl[[i]] <- quadRule_AGHQ()
      }
      if(quadRules[i] == "CCD") {
        I_CCD <- i
        quadRule_nfl[[i]] <- quadRule_CCD()
      }
      if(quadRules[i] == "USER"){
        I_USER <- i
        quadRule_nfl[[i]] <- quadRule_USER()
      }
      quadGridCache_nfl[[i]] <- quadGridCache()      
    }
    
    modeIndex <- -1
    nGrid <- 0
    gridBuilt <- FALSE
  },
	run=function(){},
	methods = list(
    ## NOCHNG means keep it as is, and nQuad = -1.
    buildGrid = function(method = character(0, default = "NOCHNG"), nQuad = integer(0, default = -1)){
      if(method != "NOCHNG")
        setRule(method)
      if(nQuad != -1)
        nQuad_ <<- nQuad
      if( !quadGridCache_nfl[[I_RULE]]$checkGrid(nQuad_, prune_) | !gridBuilt) {
        newgrid <- quadGridList_internal$new()
        newgrid <- quadRule_nfl[[I_RULE]]$buildGrid(nQuad = nQuad_, d = d)
        quadGridCache_nfl[[I_RULE]]$cacheQuadGrid(nQuad = nQuad_, nodes = newgrid$nodes, wgts = newgrid$wgts, modeIndex = newgrid$modeIndex)
        gridBuilt <<- TRUE
      }
      
      modeIndex <<- quadGridCache_nfl[[I_RULE]]$modeI()
      nGrid <<- quadGridCache_nfl[[I_RULE]]$gridSize()
    },
    ## Prune grid and then cache it again.
    pruneGrid = function(prune = double(0, default = 0)){
      if(prune > 1 | prune < 0)
        stop("Can only prune a proportion of quadrature points.")

      if(I_RULE == I_CCD)
        print("Warning:  CCD grid cannot be pruned.")

      if(I_RULE != I_CCD) {
        ## Need to rebuild the grid if pruning the grid a second time.
        if(!quadGridCache_nfl[[I_RULE]]$checkGrid(nQuad = nQuad_, prune = prune)){
          gridBuilt <<- FALSE
          buildGrid()
        }
        if(prune > 0)
          quadGridCache_nfl[[I_RULE]]$pruneGrid(prune)
      }
      prune_ <<- prune
    },
    ## Surely there is a better way to do this...
    setRule = function(method = character(0, default = "AGHQ")){
      if(method == "AGHQ")
        I_RULE <<- I_AGHQ
      if(method == "CCD")
        I_RULE <<- I_CCD
      if(method == "USER")
        I_RULE <<- I_USER
    },
    setDim = function(ndim = integer(0, default = 1)){
      if(ndim <= 0)
        stop("Can't input negative dimensions")
      else
        d <<- ndim
      ## Make sure the next grid gets built.
      gridBuilt <<- FALSE  
    },
    weights = function(indx = integer(0, default = 0)){
      if(!gridBuilt) buildGrid()   
      if(indx == -1 & modeIndex > 0)
        indx <- modeIndex      
      returnType(double(1))
      return(quadGridCache_nfl[[I_RULE]]$weights(indx = indx))
    },
    nodes = function(indx = integer(0, default = 0)){
      if(!gridBuilt) buildGrid()    
      if(indx == -1 & modeIndex > 0)
        indx <- modeIndex
      returnType(double(2)); 
      return(quadGridCache_nfl[[I_RULE]]$nodes(indx = indx))
    },
    gridSize = function(){
      if(!gridBuilt) buildGrid()    
      returnType(double())
      return(nGrid)
    },
    modeI = function(){
      if(!gridBuilt) buildGrid()    
      returnType(double())
      return(modeIndex)
    }
  )
)## End of configureQuadGrid


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


## Create a caching random effects system for simulating the posterior random effect distribution according to Stringer:
## This requires the inner mode, the inner cholesky, the 

INNER_CACHE_BASE <- nimbleFunctionVirtual(
  run = function(){},
  methods = list(
    buildCache = function(nGridUpdate = integer(), nLatentNodes = integer()){},
    cache_weights = function(weight = double(), indx = integer()){},
    cache_inner_mode = function(mode = double(1), indx = integer()){},
    cache_inner_negHessChol = function(negHessChol = double(2), indx = integer()){},
    weights = function(){
      returnType(double(1))
    },
    simulate = function(n = integer()){
      returnType(double(2))
    }
  )
)

## Things I need to add for approx posterior
## save inner mode, save inner cholesky etc. for doing simulation of random-effects
## save outer mode and negHessian. 
## Need wgt*density
## and inner mode
## and inner cholesky for each point:
inner_cache_methods = nimbleFunction(
  contains = INNER_CACHE_BASE,
  setup = function(nre = 0, nGrid = 0, condIndptSets = NULL, nCondIndptSets = 1){
    innerMode <- matrix(0, nrow = 1, ncol = 1)
    innerNegHessChol <- array(0, c(1, 1, 1))
    wgtsDens <- c(1,-1)
    cacheBuilt <- FALSE
    if(is.null(condIndptSets)) {
      condInptSets <- nre ## Assuming all one set.
      nCondIndptSets <- 1 ## If NULL then this is not relevant.
    }
    if(length(condIndptSets) == 1){
      condIndptSets <- c(condIndptSets, -1) ##  Make sure it's a vector.
    }
  },
  run = function(){},
  methods = list(
    buildCache = function(nGridUpdate = integer(0, default = -1), nLatentNodes = integer()){
      nre <<- nLatentNodes
      ## If the cond indpt sets don't match up, don't use.
      if(nre != sum(condIndptSets[1:nCondIndptSets])){
        print("  Warning: Not able to simulate latent effects from conditionally independent sets.")
        condIndptSets <<- numeric(value = nre, length = 1)
        nCondIndptSets <<- 1
      }      
    
      if( nGridUpdate > 0 & nGridUpdate != nGrid){
        nGrid <<- nGridUpdate
        cacheBuilt <<- FALSE
      }
      if(cacheBuilt){
        nGrid <<- nGridUpdate
        wgtsDens <<- numeric(value = 0, length = nGrid)
        innerMode <<- matrix(0, nrow = nGrid, ncol = nre)
        innerNegHessChol <<- array(0, c(nGrid, nre, nre))
        cacheBuilt <<- TRUE
      }
    },    
    ## Note to self, this wgt will be density*wgt, strictly for simulating.
    cache_weights = function(weight = double(), indx = integer()){
      wgtsDens[indx] <<- weight
    },
    cache_inner_mode = function(mode = double(1), indx = integer()){
      innerMode[indx,] <<- mode
    },
    ## Note potentially storing a lot of zeros here. Could break it into a list of cond indpt sets.
    cache_inner_negHessChol = function(negHessChol = double(2), indx = integer()){
      innerNegHessChol[indx,,] <<- negHessChol
    },
    weights = function(){
      returnType(double(1))
      return(wgtsDens)
    },
    ## Adding first column to be index for theta.
    simulate = function(n = integer()){
      val <- matrix(0, nrow = n, ncol = nre + 1)
      simwgt <- wgtsDens/sum(wgtsDens) ## Did log sum exp when doing input.

      ## Simulate theta points first. 
      ## Seems efficient to separate to not initiate too many index vectors for cond indpt sets.
      for( i in 1:n ) {
        k <- rcat(1, prob = simwgt)
        val[i, 1] <- k
        jStart <- 1
        for( j in 1:nCondIndptSets ){
          val[i,(jStart+1):(jStart + condIndptSets[j])] <- rmnorm_chol(n=1, 
                                  mean = innerMode[k,jStart:(jStart + condIndptSets[j] - 1)],  
                                  cholesky = innerNegHessChol[k,jStart:(jStart + condIndptSets[j] - 1),jStart:(jStart + condIndptSets[j] - 1)], 
                                  prec_param = TRUE)
        }
        jStart <- jStart + condIndptSets[j]
      }
      returnType(double(2))
      return(val)
    }
  )
)