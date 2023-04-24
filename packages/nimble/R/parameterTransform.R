

#' Automated transformations of model nodes to unconstrained scales
#'
#' `parameterTransform` provides general transformations of constrained continuous-valued model nodes (parameters) to an unconstrained scale.  It handles the cases of interval-bounded parameters (e.g. uniform or beta distributions), semi-interval-bounded parameters (e.g. exponential or gamma distributions), and the multivariate wishart, inverse wishart, dirichlet, and LKJ distributions.  Utilities are provided to transform paramters to an unconstrained scale, back-transform from the unconstrained scale to the original scale of the constrained parameterization, and to calculate the natural logarithm of the determinant of the Jacobian matrix of the inverse transformation, calculated at any location in the transformed (unconstrained) space.
#'
#' @param model A `nimble` model object.  See details.
#' 
#' @param nodes A character vector specifying one or more model node names to undergo transformation.  See details.
#'
#' @details
#' 
#' The `parameterTransform` nimbleFunction is an unspecialized function.  Calling `parameterTransform(model, nodes)` will generate and return a specialized nimbleFunction, which provides transformation functionality for the specificed hierarchical model and set of model nodes.  The `nodes` argument can represent mutliple model nodes arising from distinct prior distributions, which will be simultaneously transformed according to their respective distributions and constraints.
#'
#' This specialized nimbleFunction has the following methods:
#'
#' \code{transform}: Transforms a numeric vector of values from the original constrained model scale to a vector of values on the unconstrained scale.
#'
#' \code{inverseTransform}: Transforms a numeric vector of values from the unconstrained scale to the original constrained parameterization scale.
#'
#' The unconstrained scale may have different dimensionality from the original constrained scale of the model parameters.  For example, a d-dimensional dirichlet distribution is constrained to reside on a simplex in d-dimensional space.  In contrast, the corresponding unconstrained parameterization is unrestrained in (d-1) dimensional space.  The specialized `parameterTransform` nimbleFunction also provides utilities to return the dimensionality of the original (constrained) parameterization, and the transformed (unconstrained) parameterization:
#'
#' \code{getOriginalLength}: Returns the dimensionality (number of scalar elements) of the original constrained parameterization.
#'
#' \code{getTransformedLength}: Returns the dimensionality (number of scalar elements) comprising the transformed unconstrained parameterization.
#'
#' The specialized `parameterTransform` nimbleFunction also provides a method for calculating the natural logarithm of the jacobian of the inverse transformation, calculated at any point in the transformed (unconstrained) space:
#'
#' \code{logDetJacobian}
#'
#' The `parameterTransformation` function has no facility for handling discrete-valued parameters.
#'
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     a ~ dnorm(0, 1)
#'     b ~ dgamma(1, 1)
#'     c ~ dunif(2, 10)
#'     d[1:3] ~ dmnorm(mu[1:3], cov = C[1:3,1:3])
#'     e[1:3,1:3] ~ dwish(R = C[1:3,1:3], df = 5)
#' })
#'  
#' constants <- list(mu=rep(0,3), C=diag(3))
#'  
#' Rmodel <- nimbleModel(code, constants)
#'  
#' ## create a specialized parameterTransform function:
#' nodes <- c('a', 'b', 'c', 'd', 'e')
#' pt <- parameterTransform(Rmodel, nodes)
#'  
#' vals <- c(1, 10, 5,    1,2,3,   as.numeric(diag(3)))
#'  
#' ## transform values to unconstrained scale:
#' transformedVals <- pt$transform(vals)
#'  
#' ## back-transform to original constrained scale of parameterization
#' pt$inverseTransform(transformedVals)  ## return is same as original vals
#'  
#' ## dimensionality of original constrained scale = 1 + 1 + 1 + 3 + 9
#' pt$getOriginalLength()      ## 15
#'  
#' ## dimensionality of transformed (unconstrained) scale = 1 + 1 + 1 + 3 + 6
#' pt$getTransformedLength()   ## 12
#'  
#' ## log of the jacobian of the inverse transformation matrix:
#' pt$logDetJacobian(transformedVals)
#' }
#'
#' @seealso \code{\link{sampler_HMC}} \code{\link{buildLaplace}}
#'
#' @author Daniel Turek
#' 
#' @export
parameterTransform <- nimbleFunction(
    name = 'parameterTransform',
    setup = function(model, nodes, control = list()) {
        nodesExpanded <- model$expandNodeNames(nodes)
        allowDeterm <- if(!is.null(control$allowDeterm)) control$allowDeterm else FALSE
        if(!allowDeterm){
          if(any(model$isDeterm(nodesExpanded)))   stop(paste0('parameterTransform cannot operate on deterministic nodes: ',        paste0(nodesExpanded[model$isDeterm(nodesExpanded)],   collapse = ', ')))
          if(any(model$isDiscrete(nodesExpanded))) stop(paste0('parameterTransform cannot operate on discrete-valued nodes: ',      paste0(nodesExpanded[model$isDiscrete(nodesExpanded)], collapse = ', ')))
        }
        nNodes <- length(nodesExpanded)
        if(nNodes < 1) stop('parameterTransform requires at least one model node')
        ## transformType:
        ## 1: scalar unconstrained
        ## 2: scalar semi-interval (0, Inf)
        ## 3: scalar interval-constrained (0, 1)
        ## 4: scalar semi-interval (-Inf, b) or (a, Inf)
        ## 5: scalar interval-constrained (a, b)
        ## 6: multivariate {normal, t}
        ## 7: multivariate {wishart, inverse-wishart}
        ## 8: multivariate dirichlet
        ## 9: LKJ 
        transformType <- as.integer(rep(NA, (nNodes+1)))  ## must be integer for nimSwitch, and also ensure as a vector
        ## transformData
        transformData <- array(NA, dim = c(nNodes, 6))
        NIND1 <- 1
        NIND2 <- 2
        TIND1 <- 3
        TIND2 <- 4
        DATA1 <- 5
        DATA2 <- 6
        ##
        for(i in 1:nNodes) {
            node <- nodesExpanded[i]
            transformData[i,NIND1] <- if(i==1) 1 else transformData[i-1,NIND2]+1
            transformData[i,TIND1] <- if(i==1) 1 else transformData[i-1,TIND2]+1
            if(allowDeterm) {
              if(model$isDeterm(node)) {
                d <- length(model$expandNodeNames(node, returnScalarComponents = TRUE))
                if(d == 1) {       ## copied from case #1 below
                  transformData[i,NIND2] <- transformData[i,NIND1]
                  transformData[i,TIND2] <- transformData[i,TIND1]
                  transformType[i] <- 1L; next }
                if( d > 1) {       ## copied from case #6 below
                  transformType[i] <- 6L
                  transformData[i,NIND2] <- transformData[i,NIND1] + d - 1
                  transformData[i,TIND2] <- transformData[i,TIND1] + d - 1
                  next }
                stop("parameter transformation system is not able to configure for node: ", node, ".")
              }
            }
            dist <- model$getDistribution(node)
            if(!model$isMultivariate(node)) {   ## univariate
                transformData[i,NIND2] <- transformData[i,NIND1]
                transformData[i,TIND2] <- transformData[i,TIND1]
                bounds <- c(model$getBound(node, 'lower'), model$getBound(node, 'upper'))
                if(bounds[1] == -Inf && bounds[2] == Inf) {       ## 1: scalar unconstrained; also set for scalar determ nodes when allowDeterm is TRUE
                    transformType[i] <- 1L; next }
                if(bounds[1] == 0    && bounds[2] == Inf) {       ## 2: scalar semi-interval (0, Inf)
                    transformType[i] <- 2L; next }
                if(bounds[1] == 0    && bounds[2] == 1  ) {       ## 3: scalar interval-constrained (0, 1)
                    transformType[i] <- 3L
                    if(model$isTruncated(node)) {
                        lParam <- 'lower_'
                        uParam <- 'upper_'
                        lowerBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, lParam))
                        upperBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, uParam))
                        if(length(all.vars(lowerBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant lower bound, which cannot be used in parameterTransform.')
                        if(length(all.vars(upperBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant upper bound, which cannot be used in parameterTransform.')
                    }
                    next }
                if((isValid(bounds[1]) && bounds[2] ==  Inf) ||
                   (isValid(bounds[2]) && bounds[1] == -Inf)) {   ## 4: scalar semi-interval (-Inf, b) or (a, Inf)
                    if(model$isTruncated(node)) {
                        bdName  <- if(isValid(bounds[1])) 'lower'  else 'upper'
                        bdParam <- if(isValid(bounds[1])) 'lower_' else 'upper_'
                        bdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, bdParam))
                        if(length(all.vars(bdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant ', bdName, ' bound, which cannot be used in parameterTransform.')
                    }
                    transformType[i] <- 4L
                    transformData[i,DATA1] <- if(isValid(bounds[1])) bounds[1] else bounds[2]   ## formerly boundValue
                    transformData[i,DATA2] <- if(isValid(bounds[1])) 1 else -1                  ## formerly isLowerBound
                    next }
                if(isValid(bounds[1]) && isValid(bounds[2])) {    ## 5: scalar interval-constrained (a, b)
                    if(dist == 'dunif' || model$isTruncated(node)) {     ## uniform distribution, or a truncated node
                        if(dist == 'dunif')         { lParam <- 'min';    uParam <- 'max'    }
                        if(model$isTruncated(node)) { lParam <- 'lower_'; uParam <- 'upper_' }
                        lowerBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, lParam))
                        upperBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, uParam))
                        if(length(all.vars(lowerBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant lower bound, which cannot be used in parameterTransform.')
                        if(length(all.vars(upperBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant upper bound, which cannot be used in parameterTransform.')
                    } else {   ## some other distribution with finite support
                        message('  [Warning] `parameterTransform` system cannot process the ', dist, ' distribution of node ', node, '.\n         The upper and lower bounds of the ', dist, ' distribution must be constant.\n         If you\'re uncertain about this, please get in touch with the NIMBLE development team.')
                    }
                    transformType[i] <- 5L
                    transformData[i,DATA1] <- bounds[1]               ## formerly lowerBound
                    transformData[i,DATA2] <- bounds[2] - bounds[1]   ## formerly range
                    next }
                stop(paste0('`parameterTransform` system doesn\'t have a transformation for the bounds of node: ', node, ', which are (', bounds[1], ', ', bounds[2], ')'))
            } else {   ## multivariate
                if(dist %in% c('dmnorm', 'dmvt')) {               ## 6: multivariate {normal, t}; also set for non-scalar determ nodes when allowDeterm is TRUE
                    transformType[i] <- 6L
                    d <- length(model$expandNodeNames(node, returnScalarComponents = TRUE))
                    transformData[i,NIND2] <- transformData[i,NIND1] + d - 1
                    transformData[i,TIND2] <- transformData[i,TIND1] + d - 1
                    next }
                if(dist %in% c('dwish', 'dinvwish')) {            ## 7: multivariate {wishart, inverse-wishart}
                    transformType[i] <- 7L
                    dSq <- length(model$expandNodeNames(node, returnScalarComponents = TRUE))
                    d <- sqrt(dSq)
                    transformData[i,NIND2] <- transformData[i,NIND1] + dSq - 1
                    transformData[i,TIND2] <- transformData[i,TIND1] + d*(d+1)/2 - 1
                    transformData[i,DATA1] <- d           ## formerly d
                    transformData[i,DATA2] <- d*(d+1)/2   ## formerly tLengthOne
                    next }
                if(dist == 'ddirch') {                            ## 8: multivariate dirichlet
                    transformType[i] <- 8L
                    d <- length(model$expandNodeNames(node, returnScalarComponents = TRUE))
                    transformData[i,NIND2] <- transformData[i,NIND1] + d - 1
                    transformData[i,TIND2] <- transformData[i,TIND1] + d - 2
                    transformData[i,DATA1] <- d
                    next }
                if(dist == 'dlkj_corr_cholesky') {                ## 9: LKJ
                    transformType[i] <- 9L
                    dSq <- length(model$expandNodeNames(node, returnScalarComponents = TRUE))
                    d <- sqrt(dSq)
                    p <- d * (d-1) / 2  # number of transformed params
                    transformData[i,NIND2] <- transformData[i,NIND1] + dSq - 1
                    transformData[i,TIND2] <- transformData[i,TIND1] + p - 1
                    transformData[i,DATA1] <- d
                    transformData[i,DATA2] <- p
                    next }
                stop(paste0('parameterTransform doesn\'t handle \'', dist, '\' distributions.'), call. = FALSE)
            }
        }
        nLength <- transformData[nNodes,NIND2]
        tLength <- transformData[nNodes,TIND2]
        if(nLength != length(model$expandNodeNames(nodesExpanded, returnScalarComponents = TRUE))) stop('something wrong with nLength')
    },
    run = function() { print('Warning: run method of parameterTransform is not defined') },
    methods = list(
        getOriginalLength    = function() { returnType(double()); return(nLength) },
        getTransformedLength = function() { returnType(double()); return(tLength) },
        transform = function(nodeValuesFromModel = double(1)) {
            ## argument values(model, nodes), return vector on unconstrained scale
            transformed <- nimNumeric(tLength)
            for(iNode in 1:nNodes) {
                theseValues <- nodeValuesFromModel[transformData[iNode,NIND1]:transformData[iNode,NIND2]]
                thisType <- transformType[iNode]
                nimSwitch(thisType, 1:9,
                          theseTransformed <- theseValues,         ## 1: scalar unconstrained
                          theseTransformed <- log(theseValues),    ## 2: scalar semi-interval (0, Inf)
                          theseTransformed <- logit(theseValues),  ## 3: scalar interval-constrained (0, 1)
                          theseTransformed <- log(transformData[iNode,DATA2] * (theseValues - transformData[iNode,DATA1])),    ## 4: scalar semi-interval (-Inf, b) or (a, Inf)
                          theseTransformed <- logit((theseValues - transformData[iNode,DATA1]) / transformData[iNode,DATA2]),  ## 5: scalar interval-constrained (a, b)
                          theseTransformed <- theseValues,         ## 6: multivariate {normal, t}
                          {                                        ## 7: multivariate {wishart, inverse-wishart}
                              ## log-Cholesky transform, values are column-wise
                              dd <- transformData[iNode,DATA1]
                              valueAsMatrix <- nimArray(theseValues, dim = c(dd, dd))
                              U <- chol(valueAsMatrix)
                              ## DT: there has to be a better way to do this procedure, below,
                              ## creating the vector of the log-Cholesky transformed values.
                              theseTransformed <- nimNumeric(transformData[iNode,DATA2])
                              tInd <- 1
                              for(j in 1:dd) {
                                  for(i in 1:dd) {
                                      if(i==j) { theseTransformed[tInd] <- log(U[i,j]); tInd <- tInd+1 }
                                      if(i< j) { theseTransformed[tInd] <-     U[i,j];  tInd <- tInd+1 } } }
                          },
                          {                                        ## 8: multivariate dirichlet
                              dd <- transformData[iNode,DATA1] - 1
                              theseTransformed <- nimNumeric(dd)
                              theseTransformed[1] <- logit( theseValues[1] )
                              if(dd > 1) {
                                  runningSum <- 0
                                  for(j in 2:dd) {
                                      runningSum <- runningSum + theseValues[j-1]
                                      theseTransformed[j] <- logit( theseValues[j] / (1-runningSum) )
                                  }
                              }
                          },
                          {                                        ## 9: LKJ
                              dd <- transformData[iNode,DATA1]  # nrow of matrix
                              pp <- transformData[iNode,DATA2]  # number of transformed params
                              theseTransformed <- nimNumeric(pp)
                              theseValuesMatrix <- nimArray(theseValues, dim = c(dd, dd))  # U in matrix form
                              if(dd > 1) {
                                  cnt <- 1
                                  ## Length of each column of U is 1.
                                  ## We first produce the canonical partial correlations and then apply atanh()
                                  ## to make the unconstrained parameters..
                                  for(j in 2:dd) {
                                      theseTransformed[cnt] <- atanh(theseValuesMatrix[1, j])
                                      cnt <- cnt + 1
                                      if(j > 2) {
                                          partialSum <- 1
                                          for(i in 2:(j-1)) {
                                              partialSum <- partialSum - theseValuesMatrix[i-1, j]^2
                                              ## Transformed value is atanh of the proportion of the
                                              ## remaining correlation (which is in 'partialSum').
                                              theseTransformed[cnt] <- atanh(theseValuesMatrix[i, j] / sqrt(partialSum))
                                              cnt <- cnt + 1
                                          }
                                      }
                                  }
                              }
                          })
                transformed[transformData[iNode,TIND1]:transformData[iNode,TIND2]] <- theseTransformed
            }
            returnType(double(1))
            return(transformed)
        },
        inverseTransform = function(transformedValues = double(1)) {
            ## argument on transformed scale, return vector suitable for values(model,)
            modelValuesVector <- nimNumeric(nLength)
            iNode <- 1L; i <- 1L; j <- 1L; ind1 <- 1L; ind2 <- 1L; dd <- 1L   ## integer types
            for(iNode in 1:nNodes) {
                ind1 <- transformData[iNode,TIND1]
                ind2 <- transformData[iNode,TIND2]
                theseValues <- transformedValues[ind1:ind2]
                thisType <- transformType[iNode]
                nimSwitch(thisType, 1:9,
                          theseInvTransformed <- theseValues,          ## 1: scalar unconstrained
                          theseInvTransformed <- exp(theseValues),     ## 2: scalar semi-interval (0, Inf)
                          theseInvTransformed <- ilogit(theseValues),  ## 3: scalar interval-constrained (0, 1)
                          theseInvTransformed <- transformData[iNode,DATA1] + transformData[iNode,DATA2] * exp(theseValues),    ## 4: scalar semi-interval (-Inf, b) or (a, Inf)
                          theseInvTransformed <- transformData[iNode,DATA1] + transformData[iNode,DATA2] * expit(theseValues),  ## 5: scalar interval-constrained (a, b)
                          theseInvTransformed <- theseValues,          ## 6: multivariate {normal, t}
                          {                                            ## 7: multivariate {wishart, inverse-wishart}
                              dd <- transformData[iNode,DATA1]
                              cholAsMatrix <- nimArray(0, dim = c(dd, dd))
                              ## DT: there has to be a better way to do this procedure, below,
                              ## now creating the vector of the Wishart node values.
                              tInd <- 1L
                              for(j in 1:dd) {
                                  for(i in 1:dd) {
                                      if(i==j) { cholAsMatrix[i,j] <- exp(theseValues[tInd]); tInd <- tInd+1 }
                                      if(i< j) { cholAsMatrix[i,j] <-     theseValues[tInd];  tInd <- tInd+1 } } }
                              valuesAsMatrix <- t(cholAsMatrix) %*% cholAsMatrix
                              theseInvTransformed <- nimNumeric(dd*dd, valuesAsMatrix)
                          },
                          {                                            ## 8: multivariate dirichlet
                              dd <- transformData[iNode,DATA1]
                              ddm1 <- dd - 1L
                              theseInvTransformed <- nimNumeric(dd)
                              theseInvTransformed[1] <- ilogit( theseValues[1] )
                              if(dd > 2) {
                                  runningSum <- 0
                                  for(i in 2:ddm1) {
                                      runningSum <- runningSum + theseInvTransformed[i-1]
                                      theseInvTransformed[i] <- (1-runningSum) * ilogit( theseValues[i] )
                                  }
                              }
                              theseInvTransformed[dd] <- 1 - sum(theseInvTransformed[1:ddm1])
                          },
                          {                                            ## 9: LKJ
                              dd <- transformData[iNode,DATA1]
                              ## Directly fill in the vectorized form of the U matrix
                              theseInvTransformed <- nimNumeric(dd*dd, value = 0, init = TRUE)
                              theseInvTransformed[1] <- 1
                              if(dd > 1) {
                                  cntT <- 1L
                                  ## Fill in j'th column
                                  for(j in 2:dd) {
                                      cntI <- (j-1L)*dd + 1L
                                      ## Get elements of U by going from unconstrained to canonical partial correlations.
                                      theseInvTransformed[cntI] <- tanh(theseValues[cntT])
                                      partialSum <- 1 - theseInvTransformed[cntI]^2
                                      cntI <- cntI + 1L
                                      cntT <- cntT + 1L
                                      if(j > 2) {
                                          for(i in 2:(j-1)) {
                                              ## Value of U based on 'normalizing' partial correlation such that length of column of U is 1.
                                              theseInvTransformed[cntI] <- tanh(theseValues[cntT]) * sqrt(partialSum)
                                              partialSum <- partialSum - theseInvTransformed[cntI]^2
                                              cntI <- cntI + 1L
                                              cntT <- cntT + 1L
                                          }
                                      }
                                      theseInvTransformed[cntI] <- sqrt(partialSum)  # Column must sum to 1.
                                  }
                              }
                          })
                ind1 <- transformData[iNode,NIND1]
                ind2 <- transformData[iNode,NIND2]
                modelValuesVector[ind1:ind2] <- theseInvTransformed
            }
            returnType(double(1))
            return(modelValuesVector)
        },
        logDetJacobian = function(transformedValues = double(1)) {
            ## DT: general intended usage of this method:
            ## values(model, nodes) <- pt$inverseTransform(transformedValues)
            ## lp <- model$calculate(calcNodes) + pt$logDetJacobian(transformedValues)
            lp <- 0
            for(iNode in 1:nNodes) {
                theseValues <- transformedValues[transformData[iNode,TIND1]:transformData[iNode,TIND2]]
                thisType <- transformType[iNode]
                nimSwitch(thisType, 1:9,
                          lpAdd <- 0,               ## 1: scalar unconstrained
                          lpAdd <- theseValues[1],  ## 2: scalar semi-interval (0, Inf)
                          {                         ## 3: scalar interval-constrained (0, 1)
                              x <- theseValues[1]
                              lpAdd <- -log(exp(x)+exp(-x)+2)   ## alternate: -2*log(1+exp(-x))-x)
                          },
                          lpAdd <- theseValues[1],  ## 4: scalar semi-interval (-Inf, b) or (a, Inf)
                          {                         ## 5: scalar interval-constrained (a, b)
                              x <- theseValues[1]
                              lpAdd <- log(transformData[iNode,DATA2]) - log(exp(x)+exp(-x)+2)  ## alternate: -2*log(1+exp(-x))-x)
                          },
                          lpAdd <- 0,               ## 6: multivariate {normal, t}
                          {                         ## 7: multivariate {wishart, inverse-wishart}
                              dd <- transformData[iNode,DATA1]
                              lpAdd <- dd * log(2)
                              for(j in 1:dd)   lpAdd <- lpAdd + (dd+2-j) * theseValues[0.5*j*(j+1)] # "/2" instead of "*0.5" inserts a double cast in C++ that breaks the CppAD system at the time of this writing
                          },
                          {                         ## 8: multivariate dirichlet
                              ## copied from inverseTransform method:
                              dd <- transformData[iNode,DATA1]
                              ddm1 <- dd - 1L
                              theseInvTransformed <- nimNumeric(dd)
                              theseInvTransformed[1] <- ilogit( theseValues[1] )
                              if(dd > 2) {
                                  runningSum <- 0
                                  for(i in 2:ddm1) {
                                      runningSum <- runningSum + theseInvTransformed[i-1]
                                      theseInvTransformed[i] <- (1-runningSum) * ilogit( theseValues[i] )
                                  }
                              }
                              ## copying inverseTransform method ends here
                              x <- theseValues[1]
                              lpAdd <- -log(exp(x)+exp(-x)+2)   ## alternate: -2*log(1+exp(-x))-x)
                              if(dd > 2) {
                                  runningSum <- 0
                                  for(i in 2:ddm1) {
                                      runningSum <- runningSum + theseInvTransformed[i-1]
                                      x <- theseValues[i]
                                      lpAdd <- lpAdd + log(1-runningSum) - log(exp(x)+exp(-x)+2)   ## alternate: -2*log(1+exp(-x))-x)
                                  }
                              }
                          },
                          {                         ## 9: LKJ
                              dd <- transformData[iNode,DATA1]
                              lpAdd <- 0
                              if(dd > 1) {
                                  cntT <- 1L
                                  for(j in 2:dd) {
                                      lpAdd <- lpAdd - 2*log(cosh(theseValues[cntT]))
                                      theseInvTransformedOne <- tanh(theseValues[cntT])
                                      partialSum <- 1
                                      cntT <- cntT + 1L
                                      if(j > 2) {
                                          for(i in 2:(j-1)) {
                                              partialSum <- partialSum - theseInvTransformedOne^2
                                              lpAdd <- lpAdd - 2*log(cosh(theseValues[cntT])) + 0.5*log(partialSum)
                                              theseInvTransformedOne <- tanh(theseValues[cntT]) * sqrt(partialSum)
                                              cntT <- cntT + 1L
                                          }
                                      }
                                  }
                              }
                          })
                lp <- lp + lpAdd
            }
            returnType(double())
            return(lp)
        }
    ),
    buildDerivs = list(inverseTransform = list(),
                       logDetJacobian = list(ignore = c('iNode','j','dd','ddm1','i')))
)
