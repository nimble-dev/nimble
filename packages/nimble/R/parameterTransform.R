

## Is it ok to have this next line in here??
## It seems to be needed for package to build.
## -DT May 2020
nimbleOptions(experimentalEnableDerivs = TRUE)

#' Automated transformations of model nodes to unconstrained scales
#'
#' ADD DETAILS
#' 
#' @author Daniel Turek
#' @export
parameterTransform <- nimbleFunction(
    name = 'parameterTransform',
    setup = function(model, nodes) {
        nodesExpanded <- model$expandNodeNames(nodes)
        if(any(model$isDeterm(nodesExpanded)))   stop(paste0('parameterTransform nodes may not be deterministic: ', paste0(nodesExpanded[model$isDeterm(nodesExpanded)],   collapse = ', ')))
        if(any(model$isDiscrete(nodesExpanded))) stop(paste0('parameterTransform nodes may not be discrete: ',      paste0(nodesExpanded[model$isDiscrete(nodesExpanded)], collapse = ', ')))
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
            dist <- model$getDistribution(node)
            transformData[i,NIND1] <- if(i==1) 1 else transformData[i-1,NIND2]+1
            transformData[i,TIND1] <- if(i==1) 1 else transformData[i-1,TIND2]+1
            if(!model$isMultivariate(node)) {   ## univariate
                transformData[i,NIND2] <- transformData[i,NIND1]
                transformData[i,TIND2] <- transformData[i,TIND1]
                bounds <- c(model$getBound(node, 'lower'), model$getBound(node, 'upper'))
                if(bounds[1] == -Inf && bounds[2] == Inf) {       ## 1: scalar unconstrained
                    transformType[i] <- 1L; next }
                if(bounds[1] == 0    && bounds[2] == Inf) {       ## 2: scalar semi-interval (0, Inf)
                    transformType[i] <- 2L; next }
                if(bounds[1] == 0    && bounds[2] == 1  ) {       ## 3: scalar interval-constrained (0, 1)
                    transformType[i] <- 3L; next }
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
                        message('parameterTransform is not familiar with the ', dist, ' distribution of node ', node, '.')
                        message('It will require the upper and lower bounds of the ', dist, ' distribution to be constant.')
                        message('If you\'re uncertain about this, please get in touch with the NIMBLE development team.')
                    }
                    transformType[i] <- 5L
                    transformData[i,DATA1] <- bounds[1]               ## formerly lowerBound
                    transformData[i,DATA2] <- bounds[2] - bounds[1]   ## formerly range
                    next }
                stop(paste0('parameterTransform doesn\'t have a transformation for the bounds of node: ', node, ', which are (', bounds[1], ', ', bounds[2], ')'))
            } else {   ## multivariate
                if(dist %in% c('dmnorm', 'dmvt')) {               ## 6: multivariate {normal, t}
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
                    transformData[i,TIND2] <- transformData[i,NIND1] + d - 2
                    transformData[i,DATA1] <- d
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
                nimSwitch(thisType, 1:8,
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
                nimSwitch(thisType, 1:8,
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
                              theseInvTransformed <- nimNumeric(dd)
                              theseInvTransformed[1] <- ilogit( theseValues[1] )
                              if(dd > 2) {
                                  runningSum <- 0
                                  for(i in 2:(dd-1)) {
                                      runningSum <- runningSum + theseInvTransformed[i-1]
                                      theseInvTransformed[i] <- (1-runningSum) * ilogit( theseValues[i] )
                                  }
                              }
                              theseInvTransformed[dd] <- 1 - sum(theseInvTransformed[1:(dd-1)])
                          })
                ind1 <- transformData[iNode,NIND1]
                ind2 <- transformData[iNode,NIND2]
                modelValuesVector[ind1:ind2] <- theseInvTransformed
            }
            returnType(double(1))
            return(modelValuesVector)
        },
        calcLogDetJacobian = function(transformedValues = double(1)) {
            ## DT: general intended usage of this method:
            ## values(model, nodes) <- pt$inverseTransform(transformedValues)
            ## lp <- model$calculate(calcNodes) + pt$calcLogDetJacobian(transformedValues)
            lp <- 0
            for(iNode in 1:nNodes) {
                theseValues <- transformedValues[transformData[iNode,TIND1]:transformData[iNode,TIND2]]
                thisType <- transformType[iNode]
                nimSwitch(thisType, 1:8,
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
                              for(j in 1:dd)   lpAdd <- lpAdd + (dd+2-j) * theseValues[j*(j+1)/2]
                          },
                          {                         ## 8: multivariate dirichlet
                              ## copied from inverseTransform method:
                              dd <- transformData[iNode,DATA1]
                              theseInvTransformed <- nimNumeric(dd)
                              theseInvTransformed[1] <- ilogit( theseValues[1] )
                              if(dd > 2) {
                                  runningSum <- 0
                                  for(i in 2:(dd-1)) {
                                      runningSum <- runningSum + theseInvTransformed[i-1]
                                      theseInvTransformed[i] <- (1-runningSum) * ilogit( theseValues[i] )
                                  }
                              }
                              ## copying inverseTransform method ends here
                              x <- theseValues[1]
                              lpAdd <- -log(exp(x)+exp(-x)+2)   ## alternate: -2*log(1+exp(-x))-x)
                              if(dd > 2) {
                                  runningSum <- 0
                                  for(i in 2:(dd-1)) {
                                      runningSum <- runningSum + theseInvTransformed[i-1]
                                      x <- theseValues[i]
                                      lpAdd <- lpAdd + log(1-runningSum) - log(exp(x)+exp(-x)+2)   ## alternate: -2*log(1+exp(-x))-x)
                                  }
                              }
                          })
                lp <- lp + lpAdd
            }
            returnType(double())
            return(lp)
        }
    ),
    enableDerivs = 'inverseTransform',
    where = getLoadingNamespace()
)






###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################



####
#### first version of parameterTransformation:
####
## 
##ptNodeVirtual <- nimbleFunctionVirtual(
##    methods = list(
##        getOriginalLengthOne    = function()                 { returnType(double())  },
##        getTransformedLengthOne = function()                 { returnType(double())  },
##        transformOne            = function(mVal = double(1)) { returnType(double(1)) },
##        inverseTransformOne     = function(tVal = double(1)) { returnType(double(1)) },
##        calcLogDetJacobianOne   = function(tVal = double(1)) { returnType(double())  }
##    )
##)
## 
##ptNodeScalar <- nimbleFunction(
##    name = 'ptNodeScalar',
##    contains = ptNodeVirtual,
##    setup = function(model, node) {
##        if(!identical(node, model$expandNodeNames(node))) stop('more than one node in ptNode')
##        if(length(model$expandNodeNames(node, returnScalarComponents = TRUE)) != 1) stop('non-scalar node in ptNodeScalar')
##    },
##    run = function() { print('Warning: run method of ptNode is not defined') },
##    methods = list(
##        getOriginalLengthOne    = function() { returnType(double()); return(1) },
##        getTransformedLengthOne = function() { returnType(double()); return(1) },
##        transformOne = function(mVal = double(1)) {
##            returnType(double(1))
##            return(mVal)
##        },
##        inverseTransformOne = function(tVal = double(1)) {
##            returnType(double(1))
##            return(tVal)
##        },
##        calcLogDetJacobianOne = function(tVal = double(1)) {
##            returnType(double())
##            return(0)
##        }
##    ), where = getLoadingNamespace()
##)
## 
##ptNodeScalarSemiInterval <- nimbleFunction(
##    name = 'ptNodeScalarSemiInterval',
##    contains = ptNodeVirtual,
##    setup = function(model, node) {
##        if(!identical(node, model$expandNodeNames(node))) stop('more than one node in ptNode')
##        if(length(model$expandNodeNames(node, returnScalarComponents = TRUE)) != 1) stop('non-scalar node in ptNodeScalar')
##        bounds <- c(model$getBound(node, 'lower'), model$getBound(node, 'upper'))
##        if(model$isTruncated(node)) {
##            bdName  <- if(isValid(bounds[1])) 'lower'  else 'upper'
##            bdParam <- if(isValid(bounds[1])) 'lower_' else 'upper_'
##            bdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, bdParam))
##            if(length(all.vars(bdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant ', bdName, ' bound, which cannot be used in parameterTransform.', call. = FALSE)
##        }
##        boundValue <- if(isValid(bounds[1])) bounds[1] else bounds[2]
##        isLowerBound <- if(isValid(bounds[1])) 1 else -1
##    },
##    run = function() { print('Warning: run method of ptNode is not defined') },
##    methods = list(
##        getOriginalLengthOne    = function() { returnType(double()); return(1) },
##        getTransformedLengthOne = function() { returnType(double()); return(1) },
##        transformOne = function(mVal = double(1)) {
##            returnType(double(1))
##            return(log(isLowerBound * (mVal - boundValue)))
##        },
##        inverseTransformOne = function(tVal = double(1)) {
##            returnType(double(1))
##            return(boundValue + isLowerBound * exp(tVal))
##        },
##        calcLogDetJacobianOne = function(tVal = double(1)) {
##            returnType(double())
##            return(tVal[1])
##        }
##    ), where = getLoadingNamespace()
##)
## 
##ptNodeScalarInterval <- nimbleFunction(
##    name = 'ptNodeScalarInterval',
##    contains = ptNodeVirtual,
##    setup = function(model, node) {
##        if(!identical(node, model$expandNodeNames(node))) stop('more than one node in ptNode')
##        if(length(model$expandNodeNames(node, returnScalarComponents = TRUE)) != 1) stop('non-scalar node in ptNodeScalar')
##        bounds <- c(model$getBound(node, 'lower'), model$getBound(node, 'upper'))
##        dist <- model$getDistribution(node)
##        if(dist == 'dunif' || model$isTruncated(node)) {     ## uniform distribution, or a truncated node
##            if(dist == 'dunif')         { lParam <- 'min';    uParam <- 'max'    }
##            if(model$isTruncated(node)) { lParam <- 'lower_'; uParam <- 'upper_' }
##            lowerBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, lParam))
##            upperBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, uParam))
##            if(length(all.vars(lowerBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant lower bound, which cannot be used in parameterTransform.', call. = FALSE)
##            if(length(all.vars(upperBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant upper bound, which cannot be used in parameterTransform.', call. = FALSE)
##        } else {   ## some other distribution with finite support
##            message('parameterTransform is not familiar with the ', dist, ' distribution of node ', node, '.')
##            message('It will require the upper and lower bounds of the ', dist, ' distribution to be constant.')
##            message('If you\'re uncertain about this, please get in touch with the NIMBLE development team.')
##        }
##        lowerBound <- bounds[1]
##        range <- bounds[2] - bounds[1]
##        logRange <- log(range)
##    },
##    run = function() { print('Warning: run method of ptNode is not defined') },
##    methods = list(
##        getOriginalLengthOne    = function() { returnType(double()); return(1) },
##        getTransformedLengthOne = function() { returnType(double()); return(1) },
##        transformOne = function(mVal = double(1)) {
##            returnType(double(1))
##            return(logit((mVal - lowerBound) / range))
##        },
##        inverseTransformOne = function(tVal = double(1)) {
##            returnType(double(1))
##            return(lowerBound + range * expit(tVal))
##        },
##        calcLogDetJacobianOne = function(tVal = double(1)) {
##            x <- tVal[1]
##            returnType(double())
##            return(logRange - log(exp(x)+exp(-x)+2))   ## alternate: -2*log(1+exp(-x))-x)
##        }
##    ), where = getLoadingNamespace()
##)
## 
##ptNodeMultiMNorm <- nimbleFunction(
##    name = 'ptNodeMultiMNorm',
##    contains = ptNodeVirtual,
##    setup = function(model, node) {
##        if(!identical(node, model$expandNodeNames(node))) stop('more than one node in ptNode')
##        if(!(model$getDistribution(node) %in% c('dmnorm', 'dmvt'))) stop('wrong distribution in ptNodeMultiMNorm')
##        nLengthOne <- length(model$expandNodeNames(node, returnScalarComponents = TRUE))
##        tLengthOne <- nLengthOne
##    },
##    run = function() { print('Warning: run method of ptNode is not defined') },
##    methods = list(
##        getOriginalLengthOne    = function() { returnType(double()); return(nLengthOne) },
##        getTransformedLengthOne = function() { returnType(double()); return(tLengthOne) },
##        transformOne = function(mVal = double(1)) {
##            returnType(double(1))
##            return(mVal)
##        },
##        inverseTransformOne = function(tVal = double(1)) {
##            returnType(double(1))
##            return(tVal)
##        },
##        calcLogDetJacobianOne = function(tVal = double(1)) {
##            returnType(double())
##            return(0)
##        }
##    ), where = getLoadingNamespace()
##)
## 
##ptNodeMultiWishart <- nimbleFunction(
##    name = 'ptNodeMultiWishart',
##    contains = ptNodeVirtual,
##    setup = function(model, node) {
##        if(!identical(node, model$expandNodeNames(node))) stop('more than one node in ptNode')
##        if(!(model$getDistribution(node) %in% c('dwish', 'dinvwish'))) stop('wrong distribution in ptNodeMultiWishart')
##        nLengthOne <- length(model$expandNodeNames(node, returnScalarComponents = TRUE))
##        d <- sqrt(nLengthOne)
##        tLengthOne <- d*(d+1)/2
##    },
##    run = function() { print('Warning: run method of ptNode is not defined') },
##    methods = list(
##        getOriginalLengthOne    = function() { returnType(double()); return(nLengthOne) },
##        getTransformedLengthOne = function() { returnType(double()); return(tLengthOne) },
##        transformOne = function(mVal = double(1)) {
##            ## log-Cholesky transform, values are column-wise
##            valueAsMatrix <- nimArray(mVal, dim = c(d, d))
##            U <- chol(valueAsMatrix)
##            ## DT: there has to be a better way to do this procedure, below,
##            ## creating the vector of the log-Cholesky transformed values.
##            ## suggestions or ideas are welcomed.
##            transformedOne <- nimNumeric(tLengthOne)
##            tInd <- 1
##            for(j in 1:d) {
##                for(i in 1:d) {
##                    if(i==j) { transformedOne[tInd] <- log(U[i,j]); tInd <- tInd+1 }
##                    if(i< j) { transformedOne[tInd] <-     U[i,j];  tInd <- tInd+1 }
##                }
##            }
##            if(tInd != tLengthOne+1) stop('something wrong in transformOne method of ptNodeMultiWishart')
##            returnType(double(1))
##            return(transformedOne)
##        },
##        inverseTransformOne = function(tVal = double(1)) {
##            cholAsMatrix <- nimArray(0, dim = c(d, d))
##            ## DT: there has to be a better way to do this procedure, below,
##            ## now creating the vector of the Wishart node values.
##            ## once again, suggestions or ideas are welcomed.
##            tInd <- 1
##            for(j in 1:d) {
##                for(i in 1:d) {
##                    if(i==j) { cholAsMatrix[i,j] <- exp(tVal[tInd]); tInd <- tInd+1 }
##                    if(i< j) { cholAsMatrix[i,j] <-     tVal[tInd];  tInd <- tInd+1 }
##                }
##            }
##            if(tInd != tLengthOne+1) stop('something wrong in inverseTransformOne method of ptNodeMultiWishart')
##            valuesAsMatrix <- t(cholAsMatrix) %*% cholAsMatrix
##            valuesAsVector <- nimNumeric(nLengthOne, valuesAsMatrix)
##            returnType(double(1))
##            return(valuesAsVector)
##        },
##        calcLogDetJacobianOne = function(tVal = double(1)) {
##            lp <- d * log(2)
##            for(i in 1:d)   lp <- lp + (d+2-i) * tVal[i*(i+1)/2]
##            returnType(double())
##            return(lp)
##        }
##    ), where = getLoadingNamespace()
##)
## 
###' Automated transformations of model nodes to unconstrained scales
###'
###' ADD DETAILS
###' 
###' @author Daniel Turek
###' @export
##parameterTransform <- nimbleFunction(
##    name = 'parameterTransform',
##    setup = function(model, nodes) {
##        nodesExpanded <- model$expandNodeNames(nodes)
##        calcNodes <- model$getDependencies(nodesExpanded)   ## DT: calcNodes may be unnecessary, see notes later
##        if(any(model$isDeterm(nodesExpanded)))   stop(paste0('parameterTransform nodes may not be deterministic: ', paste0(nodesExpanded[model$isDeterm(nodesExpanded)],   collapse = ', ')))
##        if(any(model$isDiscrete(nodesExpanded))) stop(paste0('parameterTransform nodes may not be discrete: ',      paste0(nodesExpanded[model$isDiscrete(nodesExpanded)], collapse = ', ')))
##        nNodes <- length(nodesExpanded)
##        if(nNodes < 1) stop('parameterTransform requires at least one model node')
##        ##
##        ptNodeNFL <- nimbleFunctionList(ptNodeVirtual)
##        for(i in 1:nNodes) {
##            node <- nodesExpanded[i]
##            if(!model$isMultivariate(node)) {   ## univariate
##                bounds <- c(model$getBound(node, 'lower'), model$getBound(node, 'upper'))
##                if(bounds[1] == -Inf & bounds[2] == Inf) {
##                    ## unconstrained scalar node (ptNodeScalar)
##                    ptNodeNFL[[i]] <- ptNodeScalar(model, node); next }
##                if((isValid(bounds[1]) && bounds[2] ==  Inf) ||
##                   (isValid(bounds[2]) && bounds[1] == -Inf)) {
##                    ## semi-interval scalar node (ptNodeScalarSemiInterval)
##                    ptNodeNFL[[i]] <- ptNodeScalarSemiInterval(model, node); next }
##                if(isValid(bounds[1]) && isValid(bounds[2])) {
##                    ## interval-constrained scalar node (ptNodeScalarInterval)
##                    ptNodeNFL[[i]] <- ptNodeScalarInterval(model, node); next }
##            stop(paste0('parameterTransform doesn\'t have a transformation for the bounds of node: ', node, ', which are (', bounds[1], ', ', bounds[2], ')'), call. = FALSE)
##            } else {   ## multivariate
##                dist <- model$getDistribution(node)
##                if(dist %in% c('dmnorm', 'dmvt')) {
##                    ptNodeNFL[[i]] <- ptNodeMultiMNorm(model, node); next }
##                if(dist %in% c('dwish', 'dinvwish')) {
##                    ptNodeNFL[[i]] <- ptNodeMultiWishart(model, node); next }
##                stop(paste0('parameterTransform doesn\'t handle \'', dist, '\' distributions.'), call. = FALSE)
##            }
##        }
##        ##
##        nLength <- tLength <- 0
##        for(i in 1:nNodes) { nLength <- nLength + ptNodeNFL[[i]]$getOriginalLengthOne()
##                             tLength <- tLength + ptNodeNFL[[i]]$getTransformedLengthOne() }
##        if(nLength != length(model$expandNodeNames(nodesExpanded, returnScalarComponents = TRUE))) stop('something wrong with nLength', call. = FALSE)
##        ##
##        nInd <- tInd <- array(0, c(nNodes, 2))
##        for(i in 1:nNodes) {
##            nInd[i,1] <- if(i==1) 1 else nInd[i-1,2]+1
##            nInd[i,2] <- nInd[i,1] + ptNodeNFL[[i]]$getOriginalLengthOne() - 1
##            tInd[i,1] <- if(i==1) 1 else tInd[i-1,2]+1
##            tInd[i,2] <- tInd[i,1] + ptNodeNFL[[i]]$getTransformedLengthOne() - 1
##        }
##    },
##    run = function() { print('Warning: run method of parameterTransform is not defined') },
##    methods = list(
##        getOriginalLength    = function() { returnType(double()); return(nLength) },
##        getTransformedLength = function() { returnType(double()); return(tLength) },
##        transform = function(nodeValuesFromModel = double(1)) {
##            ## argument values(model, nodes), return vector on unconstrained scale
##            transformed <- nimNumeric(tLength)
##            for(i in 1:nNodes) {
##                theseValues <- nodeValuesFromModel[nInd[i,1]:nInd[i,2]]
##                transformed[tInd[i,1]:tInd[i,2]] <- ptNodeNFL[[i]]$transformOne(theseValues)
##            }
##            returnType(double(1))
##            return(transformed)
##        },
##        inverseTransform = function(transformedValues = double(1)) {
##            ## argument on transformed scale, return vector suitable for values(model,)
##            modelValuesVector <- nimNumeric(nLength)
##            for(i in 1:nNodes) {
##                theseValues <- transformedValues[tInd[i,1]:tInd[i,2]]
##                modelValuesVector[nInd[i,1]:nInd[i,2]] <- ptNodeNFL[[i]]$inverseTransformOne(theseValues)
##            }
##            returnType(double(1))
##            return(modelValuesVector)
##        },
##        calcLogDetJacobian = function(transformedValues = double(1)) {
##            ## DT: general intended usage of this method:
##            ## values(model, nodes) <- pt$inverseTransform(transformedValues)
##            ## lp <- model$calculate(calcNodes) + pt$calcLogDetJacobian(transformedValues)
##            lp <- 0
##            for(i in 1:nNodes) {
##                theseValues <- transformedValues[tInd[i,1]:tInd[i,2]]
##                lp <- lp + ptNodeNFL[[i]]$calcLogDetJacobianOne(theseValues)
##            }
##            returnType(double())
##            return(lp)
##        }
##        ## 
##        ## DT: I'm not certain if we want the *next two methods*, in whatever form,
##        ## or if it would be left to users to take care of the gradient operations
##        ## (via use of model$calculate(), and use of nimDerivs()).
##        ## Anyway, I'm sketching out how I think it might look.
##        ## The purpose of these *next two methods* is to return a gradient vector:
##        ## - the gradient of model$calculate(calcNodes),
##        ## - where calcNodes is the standard model$getDependencies(nodes), and
##        ## - the derivatives used to form the gradient vector are taken w.r.t.
##        ##   the *elements of the vector of transformed values*
##        ## 
##        ## calculateModelLP = function(transformedValues = double(1)) {   ## method name: open for discussion
##        ##     values(model, nodesExpanded) <<- inverseTransform(transformedValues)
##        ##     lp <- model$calculate(calcNodes)
##        ##     returnType(double(0))
##        ##     return(lp)
##        ## },
##        ## 
##        ## gradientModelCalculate = function(transformedValues = double(1)) {   ## method name: open for discussion
##        ##     currentModelValues <- values(model, nodesExpanded)
##        ##     derivsOutput <- nimDerivs(calculateModelLP(transformedValues), order = 1, wrt = transformedValues)
##        ##     grad <- derivsOutput$gradient[1, 1:tLength]
##        ##     values(model, nodesExpanded) <<- currentModelValues
##        ##     model$calculate(calcNodes)
##        ##     returnType(double(1))
##        ##     return(grad)
##        ## }
##        ## 
##    ), where = getLoadingNamespace()
##)








##
## current HMC sampler (with original parameter transformation system),
## copied from branch hmcAD, March 2020
##
##sampler_HMC <- nimbleFunction(
##    name = 'sampler_HMC',
##    contains = sampler_HMC_BASE,     ## note: different contains for HMC sampler
##    setup = function(model, mvSaved, target, control) {
##        ## control list extraction
##        printTimesRan  <- if(!is.null(control$printTimesRan))  control$printTimesRan  else FALSE
##        printGradient  <- if(!is.null(control$printGradient))  control$printGradient  else FALSE
##        printEpsilon   <- if(!is.null(control$printEpsilon))   control$printEpsilon   else FALSE
##        printJ         <- if(!is.null(control$printJ))         control$printJ         else FALSE
##        messages       <- if(!is.null(control$messages))       control$messages       else TRUE
##        warnings       <- if(!is.null(control$warnings))       control$warnings       else 5
##        initialEpsilon <- if(!is.null(control$initialEpsilon)) control$initialEpsilon else 0
##        gamma          <- if(!is.null(control$gamma))          control$gamma          else 0.05
##        t0             <- if(!is.null(control$t0))             control$t0             else 10
##        kappa          <- if(!is.null(control$kappa))          control$kappa          else 0.75
##        delta          <- if(!is.null(control$delta))          control$delta          else 0.65
##        deltaMax       <- if(!is.null(control$deltaMax))       control$deltaMax       else 1000
##        nwarmup        <- if(!is.null(control$nwarmup))        control$nwarmup        else -1
##        maxTreeDepth   <- if(!is.null(control$maxTreeDepth))   control$maxTreeDepth   else 10
##        ## node list generation
##        targetNodes <- model$expandNodeNames(target)
##        if(length(targetNodes) <= 0) stop('HMC sampler must operate on at least one node', call. = FALSE)
##        calcNodes <- model$getDependencies(targetNodes)
##        originalTargetAsScalars <- model$expandNodeNames(target, returnScalarComponents = TRUE)
##        targetNodesAsScalars <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
##        ## processing of bounds and transformations
##        IND_ID   <- 1   ## transformation ID: 1=identity, 2=log, 3=logit
##        IND_LB   <- 2   ## one-sided-bound: offset; interval-bounded parameters: lower-bound
##        IND_RNG  <- 3   ## one-sided-bound: -1/1;   interval-bounded parameters: range
##        IND_LRNG <- 4   ## one-sided-bound: N/A;    interval-bounded parameters: log(range)
##        maxInd <- max(sapply(grep('^IND_', ls(), value = TRUE), function(x) eval(as.name(x))))
##        d <- length(targetNodesAsScalars)
##        d2 <- max(d, 2) ## for pre-allocating vectors
##        transformNodeNums <- c(0, 0)               ## always a vector
##        transformInfo <- array(0, c(2, maxInd))    ## always an array
##        logTransformNodes   <- character()
##        logitTransformNodes <- character()
##        for(i in 1:d) {
##            node <- targetNodesAsScalars[i]
##            if(model$isDeterm(node))      stop(paste0('HMC sampler doesn\'t operate on deterministic nodes: ', node), call. = FALSE)
##            if(model$isDiscrete(node))    stop(paste0('HMC sampler doesn\'t operate on discrete nodes: ', node), call. = FALSE)
##            dist <- model$getDistribution(node)
##            bounds <- c(model$getBound(node, 'lower'), model$getBound(node, 'upper'))
##            if(!model$isMultivariate(node)) {   ## univariate node
##                if(bounds[1] == -Inf & bounds[2] == Inf) {               ## 1 = identity: support = (-Inf, Inf)
##                    1
##                } else if((isValid(bounds[1]) && bounds[2] ==  Inf) ||   ## 2 = log: support = (a, Inf)
##                          (isValid(bounds[2]) && bounds[1] == -Inf)) {   ##      or: support = (-Inf, b)
##                    if(model$isTruncated(node)) {
##                        bdName  <- if(isValid(bounds[1])) 'lower'  else 'upper'
##                        bdParam <- if(isValid(bounds[1])) 'lower_' else 'upper_'
##                        bdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, bdParam))
##                        if(length(all.vars(bdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant ', bdName, ' bound.  HMC sampler does not yet handle that, please contant the NIMBLE development team.', call. = FALSE)
##                    }
##                    logTransformNodes <- c(logTransformNodes, node)
##                    transformNodeNums <- c(transformNodeNums, i)
##                    newRow <- numeric(maxInd)
##                    newRow[IND_ID]  <- 2
##                    newRow[IND_LB]  <- if(isValid(bounds[1])) bounds[1] else bounds[2]
##                    newRow[IND_RNG] <- if(isValid(bounds[1])) 1 else -1
##                    transformInfo <- rbind(transformInfo, newRow)
##                } else if(isValid(bounds[1]) && isValid(bounds[2])) {    ## 3 = logit: support = (a, b)
##                    if(dist == 'dunif' || model$isTruncated(node)) {     ## uniform distribution, or a truncated node
##                        if(dist == 'dunif')         { lParam <- 'min';    uParam <- 'max'    }
##                        if(model$isTruncated(node)) { lParam <- 'lower_'; uParam <- 'upper_' }
##                        lowerBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, lParam))
##                        upperBdExpr <- cc_expandDetermNodesInExpr(model, model$getParamExpr(node, uParam))
##                        if(length(all.vars(lowerBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant lower bound.  HMC sampler does not yet handle that, please contant the NIMBLE development team.', call. = FALSE)
##                        if(length(all.vars(upperBdExpr)) > 0) stop('Node ', node, ' appears to have a non-constant upper bound.  HMC sampler does not yet handle that, please contant the NIMBLE development team.', call. = FALSE)
##                    } else {   ## some other distribution with finite support
##                        message('HMC sampler is not familiar with the ', dist, ' distribution of node ', node, '.')
##                        message('We\'re going to use a logit-transformation for sampling this node, but')
##                        message('this requires the upper and lower bounds of the ', dist, ' distribution are *constant*.')
##                        message('If you\'re uncertain about this, please get in touch with the NIMBLE development team.')
##                    }
##                    logitTransformNodes <- c(logitTransformNodes, node)
##                    range <- bounds[2] - bounds[1]
##                    if(range <= 0) stop(paste0('HMC sampler doesn\'t have a transformation for the bounds of node: ', node), call. = FALSE)
##                    transformNodeNums <- c(transformNodeNums, i)
##                    newRow <- numeric(maxInd)
##                    newRow[IND_ID]   <- 3
##                    newRow[IND_LB]   <- bounds[1]
##                    newRow[IND_RNG]  <- range
##                    newRow[IND_LRNG] <- log(range)
##                    transformInfo <- rbind(transformInfo, newRow)
##                } else stop(paste0('HMC sampler doesn\'t have a transformation for the bounds of node: ', node, ', which are (', bounds[1], ', ', bounds[2], ')'), call. = FALSE)
##            } else {                            ## multivariate node
##                if(!(node %in% originalTargetAsScalars)) stop(paste0('HMC sampler only operates on complete multivariate nodes. Must specify full node: ', model$expandNodeNames(node), ', or none of it'), call. = FALSE)
##                if(dist %in% c('dmnorm', 'dmvt')) {                      ## dmnorm: identity
##                    message('HMC sampler is waiting for derivatives of dmnorm() to be implemented.')  ## waiting for dmnorm() derivatives
##                    message('otherwise, HMC sampler already works on dmnorm nodes.')                  ## waiting for dmnorm() derivatives
##                    stop()                                                                            ## waiting for dmnorm() derivatives
##                    ## NOTE: no transformation should be necessary for dmnorm and dmvt nodes (??)
##                } else if(dist %in% c('dwish', 'dinvwish')) {            ## wishart: log-cholesky
##                    message('HMC sampler is waiting for derivatives of dwish() to be implemented.')   ## waiting for dwish() derivatives
##                    stop()                                                                            ## waiting for dwish() derivatives
##                    ##transformInfo[xxx, IND_ID] <- 5
##                    ## NOTE: not sure what transformation ID for this (??)  5?
##                    ## NOTE: implementing for dwish() and dinvwish() will require a slightly deeper re-design
##                    ## d <- d + sqrt(len) * (sqrt(len)+1) / 2
##                    ## dmodel <- dmodel + len
##                } else stop(paste0('HMC sampler yet doesn\'t handle \'', dist, '\' distributions.'), call. = FALSE)   ## Dirichlet ?
##            }
##        }
##        if(messages && length(logTransformNodes)   > 0) message('HMC sampler is using a log-transformation for: ',   paste0(logTransformNodes,   collapse = ', '))
##        if(messages && length(logitTransformNodes) > 0) message('HMC sampler is using a logit-transformation for: ', paste0(logitTransformNodes, collapse = ', '))
##        ## numeric value generation
##        timesRan <- 0;   epsilon <- 0;   mu <- 0;   logEpsilonBar <- 0;   Hbar <- 0
##        q <- numeric(d2);   qL <- numeric(d2);   qR <- numeric(d2);   qDiff <- numeric(d2);   qNew <- numeric(d2)
##        p <- numeric(d2);   pL <- numeric(d2);   pR <- numeric(d2);   p2 <- numeric(d2);      p3 <- numeric(d2)
##        grad <- numeric(d2);   gradFirst <- numeric(d2);   gradSaveL <- numeric(d2);   gradSaveR <- numeric(d2)
##        log2 <- log(2)
##        warningsOrig <- warnings
##        nwarmupOrig <- nwarmup
##        numDivergences <- 0
##        numTimesMaxTreeDepth <- 0
##        ## nested function and function list definitions
##        qpNLDef <- nimbleList(q  = double(1), p  = double(1))
##        btNLDef <- nimbleList(q1 = double(1), p1 = double(1), q2 = double(1), p2 = double(1), q3 = double(1), n = double(), s = double(), a = double(), na = double())
##        ## checks
##        if(!nimbleOptions('experimentalEnableDerivs')) stop('must enable NIMBLE derivates, use: nimbleOptions(experimentalEnableDerivs = TRUE)', call. = FALSE)
##        if(initialEpsilon < 0) stop('HMC sampler initialEpsilon must be positive', call. = FALSE)
##    },
##    run = function() {
##        ## No-U-Turm Sampler with Dual Averaging, Algorithm 6 from Hoffman and Gelman (2014)
##        if(timesRan == 0) {
##            if(nwarmup == -1) stop('nwarmup was not set correctly')
##            if(initialEpsilon == 0) { initializeEpsilon()                 ## no initialEpsilon value was provided
##                                  } else { epsilon <<- initialEpsilon }   ## user provided initialEpsilon
##            mu <<- log(10*epsilon)
##        }
##        timesRan <<- timesRan + 1
##        if(printTimesRan) print('============ times ran = ', timesRan)
##        if(printEpsilon)  print('epsilon = ', epsilon)
##        transformValues()              ## sets value of member data 'q'
##        for(i in 1:d)     p[i] <<- rnorm(1, 0, 1)
##        qpLogH <- logH(q, p)
##        logu <- qpLogH - rexp(1, 1)    ## logu <- lp - rexp(1, 1) => exp(logu) ~ uniform(0, exp(lp))
##        qL <<- q;   qR <<- q;   pL <<- p;   pR <<- p;   j  <- 0;   n <- 1;   s <- 1;   qNew <<- q
##        while(s == 1) {
##            v <- 2*rbinom(1, 1, 0.5) - 1    ## -1 or 1
##            if(v == -1) { btNL <- buildtree(qL, pL, logu, v, j, epsilon, qpLogH, 1)        ## first call: first = 1
##                          qL <<- btNL$q1;   pL <<- btNL$p1
##                      } else { btNL <- buildtree(qR, pR, logu, v, j, epsilon, qpLogH, 1)   ## first call: first = 1
##                               qR <<- btNL$q2;   pR <<- btNL$p2 }
##            if(btNL$s == 1)   if(runif(1) < btNL$n / n)   qNew <<- btNL$q3
##            n <- n + btNL$n
##            qDiff <<- qR - qL
##            ##s <- btNL$s * nimStep(inprod(qDiff, pL)) * nimStep(inprod(qDiff, pR))                      ## this line replaced with the next,
##            if(btNL$s == 0) s <- 0 else s <- nimStep(inprod(qDiff, pL)) * nimStep(inprod(qDiff, pR))     ## which acccounts for NaN's in btNL elements
##            if(j >= maxTreeDepth) s <- 0
##            if(printJ) {   if(j == 0) cat('j = ', j) else cat(', ', j)
##                           cat('(');   if(v==1) cat('R') else cat('L');   cat(')')
##                           if(s != 1) print(' ')   }
##            if(j >= maxTreeDepth) { numTimesMaxTreeDepth <<- numTimesMaxTreeDepth + 1 }
##            ##                      if(warnings > 0) { print('HMC sampler encountered maximum search tree depth of ', maxTreeDepth);
##            ##                                         warnings <<- warnings - 1 } }
##            j <- j + 1
##            checkInterrupt()
##        }
##        values(model, targetNodes) <<- inverseTransformValues(qNew)
##        model$calculate(calcNodes)
##        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
##        if(timesRan <= nwarmup) {
##            Hbar <<- (1 - 1/(timesRan+t0)) * Hbar + 1/(timesRan+t0) * (delta - btNL$a/btNL$na)
##            logEpsilon <- mu - sqrt(timesRan)/gamma * Hbar
##            epsilon <<- exp(logEpsilon)
##            timesRanToNegativeKappa <- timesRan^(-kappa)
##            logEpsilonBar <<- timesRanToNegativeKappa * logEpsilon + (1 - timesRanToNegativeKappa) * logEpsilonBar
##            if(timesRan == nwarmup)   epsilon <<- exp(logEpsilonBar)
##        }
##        if(warnings > 0) if(is.nan(epsilon)) { print('HMC sampler value of epsilon is NaN, with timesRan = ', timesRan); warnings <<- warnings - 1 }
##    },
##    methods = list(
##        transformValues = function() {
##            q <<- values(model, targetNodes)
##            if(length(transformNodeNums) > 2) {
##                for(i in 3:length(transformNodeNums)) {
##                    nn <- transformNodeNums[i]
##                    id <- transformInfo[i, IND_ID]                ## 1 = identity, 2 = log, 3 = logit
##                    if(id == 2) q[nn] <<- log(   (q[nn] - transformInfo[i, IND_LB]) / transformInfo[i, IND_RNG] )
##                    if(id == 3) q[nn] <<- logit( (q[nn] - transformInfo[i, IND_LB]) / transformInfo[i, IND_RNG] )
##                }
##            }
##        },
##        inverseTransformValues = function(qArg = double(1)) {
##            transformed <- qArg
##            if(length(transformNodeNums) > 2) {
##                for(i in 3:length(transformNodeNums)) {
##                    nn <- transformNodeNums[i]
##                    id <- transformInfo[i, IND_ID]                ## 1 = identity, 2 = log, 3 = logit
##                    x <- qArg[nn]
##                    if(id == 2) transformed[nn] <- transformInfo[i, IND_LB] + transformInfo[i, IND_RNG]*  exp(x)
##                    if(id == 3) transformed[nn] <- transformInfo[i, IND_LB] + transformInfo[i, IND_RNG]*expit(x)
##                }
##            }
##            returnType(double(1));   return(transformed)
##        },
##        logH = function(qArg = double(1), pArg = double(1)) {
##            values(model, targetNodes) <<- inverseTransformValues(qArg)
##            lp <- model$calculate(calcNodes) - sum(pArg^2)/2
##            if(length(transformNodeNums) > 2) {
##                for(i in 3:length(transformNodeNums)) {
##                    nn <- transformNodeNums[i]
##                    id <- transformInfo[i, IND_ID]                ## 1 = identity, 2 = log, 3 = logit
##                    x <- qArg[nn]
##                    if(id == 2) lp <- lp + x
##                    if(id == 3) lp <- lp + transformInfo[i, IND_LRNG] - log(exp(x)+exp(-x)+2)   ## alternate: -2*log(1+exp(-x))-x
##                }
##            }
##            returnType(double());   return(lp)
##        },
##        gradient = function(qArg = double(1)) {
##            values(model, targetNodes) <<- inverseTransformValues(qArg)
##            if(printGradient) { gradQ <- array(0, c(1,d)); gradQ[1,1:d] <- values(model, targetNodes); print(gradQ) }
##            derivsOutput <- derivs(model$calculate(calcNodes), order = 1, wrt = targetNodes)
##            grad <<- derivsOutput$jacobian[1, 1:d]
##            if(length(transformNodeNums) > 2) {
##                for(i in 3:length(transformNodeNums)) {
##                    nn <- transformNodeNums[i]
##                    id <- transformInfo[i, IND_ID]                ## 1 = identity, 2 = log, 3 = logit
##                    x <- qArg[nn]
##                    if(id == 2) grad[nn] <<- grad[nn]*exp(x) + 1
##                    if(id == 3) grad[nn] <<- grad[nn]*transformInfo[i, IND_RNG]*expit(x)^2*exp(-x) + 2/(1+exp(x)) - 1
##                }
##            }
##        },
##        leapfrog = function(qArg = double(1), pArg = double(1), eps = double(), first = double(), v = double()) {
##            ## Algorithm 1 from Hoffman and Gelman (2014)
##            if(first == 1) { gradient(qArg)     ## member data 'grad' is set in gradient() method
##                         } else { if(v ==  1) grad <<- gradSaveR
##                                  if(v == -1) grad <<- gradSaveL
##                                  if(v ==  2) grad <<- gradSaveL }
##            p2 <<- pArg + eps/2 * grad
##            q2 <-  qArg + eps   * p2
##            gradFirst <<- grad
##            gradient(q2)                        ## member data 'grad' is set in gradient() method
##            p3 <<- p2   + eps/2 * grad
##            if(first == 1) { if(v ==  1) { gradSaveL <<- gradFirst;   gradSaveR <<- grad }
##                             if(v == -1) { gradSaveR <<- gradFirst;   gradSaveL <<- grad }
##                             if(v ==  2) { gradSaveL <<- gradFirst                       }
##                         } else { if(v ==  1) gradSaveR <<- grad
##                                  if(v == -1) gradSaveL <<- grad }
##            if(warnings > 0) if(is.nan.vec(c(q2, p3))) { print('HMC sampler encountered a NaN value in leapfrog routine, with timesRan = ', timesRan); warnings <<- warnings - 1 }
##            returnType(qpNLDef());   return(qpNLDef$new(q = q2, p = p3))
##        },
##        initializeEpsilon = function() {
##            ## Algorithm 4 from Hoffman and Gelman (2014)
##            savedCalcNodeValues <- values(model, calcNodes)
##            transformValues()                   ## sets value of member data 'q'
##            p <<- numeric(d)                    ## keep, sets 'p' to size d on first iteration
##            for(i in 1:d)     p[i] <<- rnorm(1, 0, 1)
##            epsilon <<- 1
##            qpNL <- leapfrog(q, p, epsilon, 1, 2)            ## v = 2 is a special case for initializeEpsilon routine
##            while(is.nan.vec(qpNL$q) | is.nan.vec(qpNL$p)) {              ## my addition
##                if(warnings > 0) { print('HMC sampler encountered NaN while initializing step-size; recommend better initial values')
##                                   print('reducing initial step-size'); warnings <<- warnings - 1 }
##                epsilon <<- epsilon / 1000                                ## my addition
##                qpNL <- leapfrog(q, p, epsilon, 0, 2)                     ## my addition
##            }                                                             ## my addition
##            qpLogH <- logH(q, p)
##            a <- 2*nimStep(exp(logH(qpNL$q, qpNL$p) - qpLogH) - 0.5) - 1
##            if(is.nan(a)) if(warnings > 0) { print('HMC sampler caught acceptance prob = NaN in initializeEpsilon routine'); warnings <<- warnings - 1 }
##            ## while(a * (logH(qpNL$q, qpNL$p) - qpLogH) > -a * log2) {   ## replaced by simplified expression:
##            while(a * (logH(qpNL$q, qpNL$p) - qpLogH + log2) > 0) {
##                epsilon <<- epsilon * 2^a
##                qpNL <- leapfrog(q, p, epsilon, 0, 2)        ## v = 2 is a special case for initializeEpsilon routine
##            }
##            values(model, calcNodes) <<- savedCalcNodeValues
##        },
##        buildtree = function(qArg = double(1), pArg = double(1), logu = double(), v = double(), j = double(), eps = double(), logH0 = double(), first = double()) {
##            ## Algorithm 6 (second half) from Hoffman and Gelman (2014)
##            returnType(btNLDef())
##            if(j == 0) {    ## one leapfrog step in the direction of v
##                qpNL <- leapfrog(qArg, pArg, v*eps, first, v)
##                q <<- qpNL$q;   p <<- qpNL$p;   qpLogH <- logH(q, p)
##                n <- nimStep(qpLogH - logu)          ## step(x) = 1 iff x >= 0, and zero otherwise
##                s <- nimStep(qpLogH - logu + deltaMax)
##                ## lowering the initial step size, and increasing the target acceptance rate may keep the step size small to avoid divergent paths.
##                if(s == 0) { numDivergences <<- numDivergences + 1 }
##                ##           if(warnings > 0) { print('HMC sampler encountered a divergent path on iteration ', timesRan, ', with divergence = ', logu - qpLogH)
##                ##                              warnings <<- warnings - 1 } }
##                a <- min(1, exp(qpLogH - logH0))
##                if(is.nan.vec(q) | is.nan.vec(p)) { n <- 0; s <- 0; a <- 0 }     ## my addition
##                return(btNLDef$new(q1 = q, p1 = p, q2 = q, p2 = p, q3 = q, n = n, s = s, a = a, na = 1))
##            } else {        ## recursively build left and right subtrees
##                btNL1 <- buildtree(qArg, pArg, logu, v, j-1, eps, logH0, 0)
##                if(btNL1$s == 1) {
##                    if(v == -1) { btNL2 <- buildtree(btNL1$q1, btNL1$p1, logu, v, j-1, eps, logH0, 0)   ## recursive calls: first = 0
##                                  btNL1$q1 <- btNL2$q1;   btNL1$p1 <- btNL2$p1
##                              } else {
##                                  btNL2 <- buildtree(btNL1$q2, btNL1$p2, logu, v, j-1, eps, logH0, 0)   ## recursive calls: first = 0
##                                  btNL1$q2 <- btNL2$q2;   btNL1$p2 <- btNL2$p2 }
##                    nSum <- btNL1$n + btNL2$n
##                    if(nSum > 0)   if(runif(1) < btNL2$n / nSum)   btNL1$q3 <- btNL2$q3
##                    qDiff <<- btNL1$q2-btNL1$q1
##                    btNL1$a  <- btNL1$a  + btNL2$a
##                    btNL1$na <- btNL1$na + btNL2$na
##                    btNL1$s  <- btNL2$s * nimStep(inprod(qDiff, btNL1$p1)) * nimStep(inprod(qDiff, btNL1$p2))
##                    btNL1$n  <- nSum
##                }
##                return(btNL1)
##            }
##        },
##        getNwarmup              = function() { returnType(double());   return(nwarmup)              },
##        getMaxTreeDepth         = function() { returnType(double());   return(maxTreeDepth)         },
##        getNumDivergences       = function() { returnType(double());   return(numDivergences)       },
##        getNumTimesMaxTreeDepth = function() { returnType(double());   return(numTimesMaxTreeDepth) },
##        setNwarmup              = function(x = double()) { nwarmup <<- x },
##        reset = function() {
##            timesRan       <<- 0
##            epsilon        <<- 0
##            mu             <<- 0
##            logEpsilonBar  <<- 0
##            Hbar           <<- 0
##            numDivergences <<- 0
##            numTimesMaxTreeDepth <<- 0
##            warnings       <<- warningsOrig
##            nwarmup        <<- nwarmupOrig
##        }
##    ), where = getLoadingNamespace()
##)


