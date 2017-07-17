


#' Convert CAR structural parameters to adjecency, weights, num format
#'
#' This will convert alternate representations of CAR process structure into
#' (adj, weights, num) form required by dcar_normal.  Two alternate
#' representations are handled:
#'
#' A single matrix argument will be interpreted as a matrix of symmetric un-normalized weights;
#'
#' Two lists will be interpreted as (the first) a list of numeric vectors
#' specifying the adjacency (neighboring) indices of each CAR process component,
#' and (the second) a list of numeric vectors giving the un-normalized weights
#' for each of these neighboring relationships.
#'
#' @author Daniel Turek
#' @export
as.carAdjacency <- function(...) {
    args <- list(...)
    if(length(args) == 1) return(CAR_convertWeightMatrix(...))
    if(length(args) == 2) return(CAR_convertNeighborWeightLists(...))
    stop('wrong arguments to as.carAdjacency')
}


CAR_convertNeighborWeightLists <- function(neighborList, weightList) {
    adj <- unlist(neighborList)
    weights <- unlist(weightList)
    num <- sapply(neighborList, length)
    return(list(adj=adj, weights=weights, num=num))
}


CAR_convertWeightMatrix <- function(weightMatrix) {
    neighborList <- apply(weightMatrix, 1, function(row) which(row > 0))
    weightList <- apply(weightMatrix, 1, function(row) row[which(row > 0)])
    return(CAR_convertNeighborWeightLists(neighborList, weightList))
}


## specialized conjugacy-checking of the scalar component nodes of dcar_normal() AND dcar_proper() distributions
CAR_checkConjugacy <- function(model, target) {
    depNodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE)
    for(depNode in depNodes) {
        if(!model$getDistribution(depNode) == 'dnorm')   return(FALSE)
        if(model$isTruncated(depNode))   return(FALSE)
        linearityCheckExpr <- model$getParamExpr(depNode, 'mean')
        linearityCheckExpr <- cc_expandDetermNodesInExpr(model, linearityCheckExpr, skipExpansions=TRUE)
        if(!cc_nodeInExpr(target, linearityCheckExpr))   return(FALSE)
        linearityCheck <- cc_checkLinearity(linearityCheckExpr, target)
        if(!cc_linkCheck(linearityCheck, 'linear'))   return(FALSE)
        if(!cc_otherParamsCheck(model, depNode, target, skipExpansions=TRUE))   return(FALSE)
    }
    return(TRUE)
}


## checks validity of adj, weights, num parameters for dcar_normal distribution
CAR_normal_checkAdjWeightsNum <- function(adj, weights, num) {
    if(any(floor(num) != num)) stop('num argument to dcar_normal() can only contain non-negative integers')
    if(any(num < 0)) stop('num argument to dcar_normal() can only contain non-negative integers')
    if(any(num > length(num))) stop('entries in num argument to dcar_normal() cannot exceed length of dcar_normal() node')
    if(sum(num) == 0) stop('dcar_normal() distribution must specify some neighbors')
    if(sum(num) != length(adj)) stop('length of adj argument to dcar_normal() must be equal to total number of neighbors specified in num argument')
    if(length(adj) != length(weights)) stop('length of adj and weight arguments to dcar_normal() must be the same')
    if(any(weights <= 0)) stop('weights argument to dcar_normal() should only contain positive values')
}


## checks validity of adj, num, C, M parameters for dcar_proper distribution
CAR_proper_checkAdjNumCM <- function(adj, num, C, M) {
    if(any(floor(num) != num)) stop('num argument to dcar_proper() can only contain non-negative integers')
    if(any(num < 0)) stop('num argument to dcar_proper() can only contain non-negative integers')
    if(any(num > length(num))) stop('entries in num argument to dcar_proper() cannot exceed length of dcar_proper() node')
    if(sum(num) == 0) stop('dcar_proper() distribution must specify some neighbors')
    if(length(num) != length(M)) stop('length of num and M arguments to dcar_proper() must be the same')
    if(sum(num) != length(adj)) stop('length of adj argument to dcar_proper() must be equal to total number of neighbors specified in num argument')
    if(length(adj) != length(C)) stop('length of adj and C arguments to dcar_proper() must be the same')
    if(any(C <= 0)) stop('C argument to dcar_proper() should only contain positive values')
    if(any(C > 1)) stop('C argument to dcar_proper() values cannot exceed one (normalised weights)')
    if(any(M <= 0)) stop('M argument to dcar_proper() should only contain positive values (conditional variances)')
    subsetIndList <- CAR_calcSubsetIndList(adj, num)
    for(i in seq_along(subsetIndList)) {
        ind <- subsetIndList[[i]]
        if(length(ind) > 0) {
            if(sum(C[ind]) != 1) stop(paste0('C (normalised weights) for component ', i, 'of dcar_proper distribution must sum to one'))
            for(ind.value in ind) {
                j <- adj[ind.value]
                if(C[ind.value]*M[j] != C[subsetIndList[[j]]][which(subsetIndList[[j]]==i)]*M[i])
                    stop('failing dcar_proper() symmetry constraint Cij*Mjj = Cji*Mii for components i = ', i, ' and j = ', j)
            }
        }
    }
}


## determines indicies of adj (and weights) vectors corresponding to each node
CAR_calcSubsetIndList <- function(adj, num) {
    d <- length(num)
    subsetIndList <- vector('list', d)
    nextInd <- 1
    for(i in 1:d) {
        subsetIndList[[i]] <- numeric()
        if(num[i] != 0) subsetIndList[[i]] <- nextInd:(nextInd+num[i]-1)
        nextInd <- nextInd + num[i]
    }
    if(nextInd != length(adj)+1) stop('something went wrong')
    return(subsetIndList)
}


## checks validity of neighborIndList and neighborWeightList for dcar_normal() distribution
CAR_normal_checkNeighborIndWeightLists <- function(neighborIndList, neighborWeightList, targetAsScalar) {
    for(i in 1:length(neighborIndList)) {
        for(j in seq_along(neighborIndList[[i]])) {
            neighborInd <- neighborIndList[[i]][j]
            if(i == neighborInd)   ## no nodes are their own neighbors
                stop(paste0('node ', targetAsScalar[i], ' in dcar_normal() distribution is neighbors with itself (which is not allowed)'))
            if(!(i %in% neighborIndList[[neighborInd]]))  ## neighbor relationships are symmetric
                stop(paste0('neighbor of node ', targetAsScalar[i], ' in dcar_normal() distribution are not symmetric with node ', targetAsScalar[j]))
            if(neighborWeightList[[neighborInd]][which(neighborIndList[[neighborInd]] == i)] != neighborWeightList[[i]][j])   ## weights symmetric
                stop(paste0('weight of node ', targetAsScalar[i], ' with ', targetAsScalar[j], ' in dcar_normal() distribution is not symmetric'))
        }
    }
}


## checks validity of neighborIndList for dcar_proper() distribution
CAR_proper_checkNeighborIndList <- function(neighborIndList, targetAsScalar) {
    for(i in 1:length(neighborIndList)) {
        for(j in seq_along(neighborIndList[[i]])) {
            neighborInd <- neighborIndList[[i]][j]
            if(i == neighborInd)   ## no nodes are their own neighbors
                stop(paste0('node ', targetAsScalar[i], ' in dcar_proper() distribution is neighbors with itself (which is not allowed)'))
            if(!(i %in% neighborIndList[[neighborInd]]))  ## neighbor relationships are symmetric
                stop(paste0('neighbor of node ', targetAsScalar[i], ' in dcar_proper() distribution are not symmetric with node ', targetAsScalar[j]))
        }
    }
}


##' Calculate number of islands (non-recursive nimbleFunction version).
##'
##' @export
CAR_calcNumIslands <- nimbleFunction(
    name = 'CAR_calcNumIslands',
    run = function(adj = double(1), num = double(1)) {
        N <- dim(num)[1]
        numIslands <- 0
        visited <- rep(0, N)
        for(i in 1:N) {
            if(visited[i] == 0) {
                visited[i] <- 1
                numIslands <- numIslands + 1
                nNeighbors <- num[i]
                if(nNeighbors > 0) {
                    adjStartInd <- 1
                    if(i > 1) adjStartInd <- adjStartInd + sum(num[1:(i-1)])
                    indToVisit <- numeric(nNeighbors)
                    indToVisit[1:nNeighbors] <- adj[adjStartInd:(adjStartInd+nNeighbors-1)]
                    lengthIndToVisit <- nNeighbors
                    l <- 1
                    while(l <= lengthIndToVisit) {
                        nextInd <- indToVisit[l]
                        if(visited[nextInd] == 0) {
                            visited[nextInd] <- 1
                            newNneighbors <- num[nextInd]
                            if(newNneighbors > 0) {
                                newIndToVisit <- numeric(newNneighbors)
                                adjStartInd <- 1
                                if(nextInd > 1) adjStartInd <- adjStartInd + sum(num[1:(nextInd-1)])
                                new_indToVisit <- c(indToVisit, adj[adjStartInd:(adjStartInd+newNneighbors-1)])
                                new_lengthIndToVisit <- lengthIndToVisit + newNneighbors
                                lengthIndToVisit <- new_lengthIndToVisit
                                indToVisit <- new_indToVisit
                            }
                        }
                        l <- l + 1
                    }
                }
            }
        }
        returnType(double())
        return(numIslands)
    }
)


## checking of the validity of the adj, weights, and num parameters to dcar_normal(),
## also processes these arguments into lists, easily usable in nimbleFunction setup code
CAR_normal_processParams <- function(model, target, adj, weights, num) {
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    if(length(targetAsScalar) != length(num)) stop('num argument to dcar_normal() must be same length as dcar_normal() node')
    CAR_normal_checkAdjWeightsNum(adj, weights, num)
    subsetIndList <- CAR_calcSubsetIndList(adj, num)
    neighborIndList <- lapply(subsetIndList, function(ind) adj[ind])
    neighborNodeList <- lapply(neighborIndList, function(ind) targetAsScalar[ind])
    neighborWeightList <- lapply(subsetIndList, function(ind) weights[ind])
    names(neighborIndList) <- names(neighborNodeList) <- names(neighborWeightList) <- targetAsScalar
    CAR_normal_checkNeighborIndWeightLists(neighborIndList, neighborWeightList, targetAsScalar)
    return(list(neighborIndList=neighborIndList, neighborNodeList=neighborNodeList, neighborWeightList=neighborWeightList))
}


## checking of the validity of the C, adj, num parameters to dcar_proper(),
## also processes these arguments into lists, easily usable in nimbleFunction setup code
CAR_proper_processParams <- function(model, target, C, adj, num, M) {
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    if(length(targetAsScalar) != length(num)) stop('num argument to dcar_proper() must be same length as dcar_proper() node')
    CAR_proper_checkAdjNumCM(adj, num, C, M)
    subsetIndList <- CAR_calcSubsetIndList(adj, num)
    neighborIndList <- lapply(subsetIndList, function(ind) adj[ind])
    neighborNodeList <- lapply(neighborIndList, function(ind) targetAsScalar[ind])
    neighborCList <- lapply(subsetIndList, function(ind) C[ind])
    names(neighborIndList) <- names(neighborNodeList) <- names(neighborCList) <- targetAsScalar
    CAR_proper_checkNeighborIndList(neighborIndList, targetAsScalar)
    return(list(neighborIndList=neighborIndList, neighborNodeList=neighborNodeList, neighborCList=neighborCList))
}


CAR_evaluateDensity_base <- nimbleFunctionVirtual(
    run = function() { returnType(double()) },
    methods = list(
        getMean = function() { returnType(double()) },
        getPrec = function() { returnType(double()) }
    )
)


## evaluates the conditional density of scalar components of dcar_normal() nodes:
## p(x_i | x_-i, tau), where x[1:N] ~ dcar_normal()
CAR_normal_evaluateDensity <- nimbleFunction(
    name = 'CAR_normal_evaluateDensity',
    contains = CAR_evaluateDensity_base,
    setup = function(model, targetScalar, neighborNodes, neighborWeights) {
        targetDCAR <- model$expandNodeNames(targetScalar)
        island <- length(neighborNodes)==0
        numNeighbors <- length(neighborWeights)                        ## fix length-1 neighborWeights
        neighborWeights <- array(neighborWeights, c(1, numNeighbors))  ## fix length-1 neighborWeights
        sumWeights <- sum(neighborWeights)
        if(length(targetDCAR) != 1)                              stop('something went wrong')
        if(model$getDistribution(targetDCAR) != 'dcar_normal')   stop('something went wrong')
    },
    run = function() {
        priorMean <- getMean()
        priorSigma <- sqrt(1/getPrec())
        lp <- dnorm(model[[targetScalar]], priorMean, priorSigma, log = TRUE)
        returnType(double())
        if(island) return(0)
        return(lp)
    },
    methods = list(
        getMean = function() {
            if(island) return(0)
            neighborValues <- values(model, neighborNodes)
            mean <- sum(neighborValues*neighborWeights[1,1:numNeighbors]) / sumWeights
            returnType(double())
            return(mean)
        },
        getPrec = function() {
            prec <- model$getParam(targetDCAR, 'tau') * sumWeights
            returnType(double())
            return(prec)
        }
    )
)


## evaluates the conditional density of scalar components of dcar_proper() nodes:
## p(x_i | x_-i, mu, tau, gamma, Mi), where x[1:N] ~ dcar_proper()
CAR_proper_evaluateDensity <- nimbleFunction(
    name = 'CAR_proper_evaluateDensity',
    contains = CAR_evaluateDensity_base,
    setup = function(model, targetScalar, neighborNodes, neighborCs, Mi) {
        targetDCAR <- model$expandNodeNames(targetScalar)
        targetDCARscalarComponents <- model$expandNodeNames(targetScalar, returnScalarComponents = TRUE)
        targetIndex <- which(targetDCARscalarComponents == targetScalar)
        island <- length(neighborNodes)==0
        numNeighbors <- length(neighborCs)                        ## fix length-1 neighborCs
        neighborCs <- array(neighborCs, c(1, numNeighbors))       ## fix length-1 neighborCs
        if(length(targetDCAR) != 1)                              stop('something went wrong')
        if(model$getDistribution(targetDCAR) != 'dcar_proper')   stop('something went wrong')
        if(length(targetIndex) != 1)                             stop('something went wrong')
    },
    run = function() {
        priorMean <- getMean()
        priorSigma <- sqrt(1/getPrec())
        lp <- dnorm(model[[targetScalar]], priorMean, priorSigma, log = TRUE)
        returnType(double())
        if(island) return(0)
        return(lp)
    },
    methods = list(
        getMean = function() {
            if(island) return(0)
            mu <- model$getParam(targetDCAR, 'mu')[targetIndex]
            gamma <- model$getParam(targetDCAR, 'gamma')
            neighborValues <- values(model, neighborNodes)
            mean <- mu + gamma * sum(neighborCs[1,1:numNeighbors] * (neighborValues - mu))
            returnType(double())
            return(mean)
        },
        getPrec = function() {
            prec <- model$getParam(targetDCAR, 'tau') / Mi
            returnType(double())
            return(prec)
        }
    )
)











