

#' Convert CAR lists of neighbors and weights to adj, weights, num format
#' @author Daniel Turek
#' @export
CAR_convertNeighborWeightLists <- function(neighborList, weightList) {
    adj <- unlist(neighborList)
    weights <- unlist(weightList)
    num <- sapply(neighborList, length)
    return(list(adj=adj, weights=weights, num=num))
}


#' Convert CAR weight matrix to adj, weights, num format
#' @author Daniel Turek
#' @export
CAR_convertWeightMatrix <- function(weightMatrix) {
    neighborList <- apply(weightMatrix, 1, function(row) which(row > 0))
    weightList <- apply(weightMatrix, 1, function(row) row[which(row > 0)])
    return(CAR_convertNeighborWeightLists(neighborList, weightList))
}


## specialized conjugacy-checking of the scalar component nodes of dcar_normal() distribution
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


## checks validity of adj, weights, num parameters
CAR_checkAdjWeightsNum <- function(adj, weights, num) {
    if(any(floor(num) != num)) stop('num argument to dcar_normal() can only contain non-negative integers')
    if(any(num < 0)) stop('num argument to dcar_normal() can only contain non-negative integers')
    if(any(num > length(num))) stop('entries in num argument to dcar_normal() cannot exceed length of dcar_normal() node')
    if(sum(num) == 0) stop('dcar_normal() distribution must specify some neighbors')
    if(sum(num) != length(adj)) stop('length of adj argument to dcar_normal() must be equal to total number of neighbors specified in num argument')
    if(length(adj) != length(weights)) stop('length of adj and weight arguments to dcar_normal() must be the same')
    if(any(weights <= 0)) stop('weights argument to dcar_normal() should only contain positive values')
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


## checks validity of neighborIndList and neighborWeightList
CAR_checkNeighborIndWeightLists <- function(neighborIndList, neighborWeightList, targetAsScalar) {
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


## helper to calculate number of islands -- recursive version, won't compile
##CAR_markVisited <- nimbleFunction(
##    run = function(adj = double(1), num = double(1), visited = double(1), index = double()) {
##        if(visited[index] == 0) {
##            visited[index] <- 1
##            adjStartInd <- 1
##            i <- 1
##            while(i < index) {
##                adjStartInd <- adjStartInd + num[i]
##                i <- i + 1
##            }
##            i <- 0
##            while(i < num[index]) {
##                visited <- CAR_markVisited(adj, num, visited, adj[adjStartInd+i])
##                i <- i + 1
##            }
##        }
##        returnType(double(1))
##        return(visited)
##    }
##)
##
#### calculate number of islands -- recursive version, won't compile
##CAR_calcNumIslands <- nimbleFunction(
##    run = function(adj = double(1), num = double(1)) {
##        N <- dim(num)[1]
##        numIslands <- 0
##        visited <- rep(0, N)
##        for(i in 1:N) {
##            if(visited[i] == 0) {
##                numIslands <- numIslands + 1
##                visited <- CAR_markVisited(adj, num, visited, i)
##            }
##        }
##        print('numIslands: ', numIslands)   ## XXXXXX delete!
##        returnType(double())
##        return(numIslands)
##    }
##)


## calculate number of islands
## non-recursive nimbleFunction version
CAR_calcNumIslands <- nimbleFunction(
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


## does (optional) extensive checking of the validity of the adj, weights, and num parameters to dcar_normal(),
## also processes these arguments into lists, easily usable in nimbleFunction setup code
CAR_processParams <- function(model, target, adj, weights, num) {
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    if(length(targetAsScalar) != length(num)) stop('num argument to dcar_normal() must be same length as dcar_normal() node')
    CAR_checkAdjWeightsNum(adj, weights, num)
    subsetIndList <- CAR_calcSubsetIndList(adj, num)
    neighborIndList <- lapply(subsetIndList, function(ind) adj[ind])
    neighborNodeList <- lapply(neighborIndList, function(ind) targetAsScalar[ind])
    neighborWeightList <- lapply(subsetIndList, function(ind) weights[ind])
    names(neighborIndList) <- names(neighborNodeList) <- names(neighborWeightList) <- targetAsScalar
    CAR_checkNeighborIndWeightLists(neighborIndList, neighborWeightList, targetAsScalar)
    return(list(neighborIndList=neighborIndList, neighborNodeList=neighborNodeList, neighborWeightList=neighborWeightList))
}


## evaluates the conditional density of scalar components of dcar_normal() nodes:
## p(x_i | x_-i, tau), where x[1:N] ~ dcar_normal()
CAR_evaluateDensity <- nimbleFunction(
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











