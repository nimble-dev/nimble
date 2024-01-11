
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


#' Convert CAR structural parameters to adjacency, weights, num format
#'
#' This will convert alternate representations of CAR process structure into
#' (adj, weights, num) form required by \code{dcar_normal}.  
#'
#' @details
#' 
#' Two alternate representations are handled:
#'
#' A single matrix argument will be interpreted as a matrix of symmetric unnormalized weights;
#'
#' Two lists will be interpreted as (the first) a list of numeric vectors
#' specifying the adjacency (neighboring) indices of each CAR process component,
#' and (the second) a list of numeric vectors giving the unnormalized weights
#' for each of these neighboring relationships.
#'
#' @param ... Either: a symmetric matrix of unnormalized weights, or two lists specifying adjacency indices and the corresponding unnormalized weights.
#'
#' @seealso \code{\link{CAR-Normal}}
#'
#' @author Daniel Turek
#' @export
as.carAdjacency <- function(...) {
    args <- list(...)
    if(length(args) == 1) return(CAR_convertWeightMatrix(...))
    if(length(args) == 2) return(CAR_convertNeighborWeightLists(...))
    stop('as.carAdjacency: wrong arguments to as.carAdjacency')
}

#' Convert weights vector to parameters of \code{dcar_proper} distributio
#' 
#' Convert weights vector to \code{C} and \code{M} parameters of \code{dcar_proper} distribution
#'
#' @details
#' Given a symmetric matrix of unnormalized weights, this function will calculate corresponding values for the \code{C} and \code{M} arguments suitable for use in the \code{dcar_proper} distribution.  This function can be used to transition between usage of \code{dcar_normal} and \code{dcar_proper}, since \code{dcar_normal} uses the \code{adj}, \code{weights}, and \code{num} arguments, while \code{dcar_proper} requires \code{adj}, \code{num}, and also the \code{C} and \code{M} as returned by this function.
#'
#' Here, \code{C} is a sparse vector representation of the row-normalized adjacency matrix, and \code{M} is a vector containing the conditional variance for each region.  The resulting values of \code{C} and \code{M} are guaranteed to satisfy the symmetry constraint imposed on \eqn{C} and \eqn{M}, that \eqn{M^{-1} C} is symmetric, where \eqn{M} is a diagonal matrix and \eqn{C} is the row-normalized adjacency matrix.
#' 
#' @param adj vector of indices of the adjacent locations (neighbors) of each spatial location.  This is a sparse representation of the full adjacency matrix.
#' @param weights vector of symmetric unnormalized weights associated with each pair of adjacent locations, of the same length as adj.  This is a sparse representation of the full (unnormalized) weight matrix.
#' @param num vector giving the number of neighbors of each spatial location, with length equal to the total number of locations.
#'
#' @return A named list with elements \code{C} and \code{M}.  These may be used as the \code{C} and \code{M} arguments to the \code{dcar_proper} distribution.
#'
#' @seealso \code{\link{CAR-Normal}}, \code{\link{CAR-Proper}}
#'
#' @author Daniel Turek
#' @export
as.carCM <- function(adj, weights, num) {
    CAR_normal_checkAdjWeightsNum(adj, weights, num)
    N <- length(num)
    L <- length(adj)
    M <- rep(-1234, N)
    C <- rep(0, L)
    for(i in 1:N) {
        if(num[i] > 0) {
            startInd <- 1
            if(i > 1) startInd <- startInd + sum(num[1:(i-1)])
            theseInd <- startInd:(startInd+num[i]-1)
            theseWeights <- weights[theseInd]
            wPlus <- sum(theseWeights)
            C[theseInd] <- theseWeights / wPlus
            M[i] <- 1/wPlus
        }
    }
    return(list(C = C, M = M))
}


## specialized conjugacy-checking of the scalar component nodes of dcar_normal() AND dcar_proper() distributions
CAR_checkConjugacy <- function(model, target, carNode) {
    depNodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE)
    for(depNode in depNodes) {
        if(!model$getDistribution(depNode) == 'dnorm')   return(FALSE)
        if(model$isTruncated(depNode))   return(FALSE)
        linearityCheckExprRaw <- model$getParamExpr(depNode, 'mean')
        linearityCheckExpr <- cc_expandDetermNodesInExpr(model, linearityCheckExprRaw, skipExpansionsNode=carNode)
        if(!cc_nodeInExpr(target, linearityCheckExpr))   return(FALSE)
        linearityCheck <- cc_checkLinearity(linearityCheckExpr, target)
        check <- cc_linkCheck(linearityCheck, 'linear')
        if(is.null(check)) return(FALSE)
        if(!cc_otherParamsCheck(model, depNode, target, skipExpansionsNode=carNode, depNodeExprExpanded = linearityCheckExpr, depParamNodeName = 'mean'))   return(FALSE)
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
    ## if(any(weights <= 0)) stop('weights argument to dcar_normal() should only contain positive values')  ## negative OK
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
    ## if(any(C <= 0)) stop('C argument to dcar_proper() should only contain positive values')  ## negative OK
    ## if(any(C > 1)) stop('C argument to dcar_proper() values cannot exceed one')  ## don't enforce
    ## Mi = -1234 indicates M was auto-generated, and region i has 0 neighbors:
    if(any((M < 0) & (M != -1234))) stop('M argument to dcar_proper() should only contain positive values (conditional variances)')
    subsetIndList <- CAR_calcSubsetIndList(adj, num)
    for(i in seq_along(subsetIndList)) {
        ind <- subsetIndList[[i]]
        if(length(ind) > 0) {
            ##if(sum(C[ind]) != 1) stop(paste0('C (weights) for region ', i, 'of dcar_proper() distribution must sum to one'))  ## don't enforce
            for(ind.value in ind) {
                j <- adj[ind.value]
                if(round(C[ind.value]*M[j] - C[subsetIndList[[j]]][which(adj[subsetIndList[[j]]]==i)]*M[i], 10) != 0)
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
    if(nextInd != length(adj)+1) stop('dcar distribution internal error')
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


#' Calculate number of islands based on a CAR adjacency matrix.
#'
#' Calculate number of islands (distinct connected groups) based on a CAR adjacency matrix.
#' 
#' @param adj vector of indices of the adjacent locations (neighbors) of each spatial location.  This is a sparse representation of the full adjacency matrix.
#' @param num vector giving the number of neighbors of each spatial location, with length equal to the total number of locations.
#'
#' @seealso \code{\link{CAR-Normal}}
#' 
#' @author Daniel Turek
#' @export
CAR_calcNumIslands <- nimbleFunction(
    name = 'CAR_calcNumIslands',
    run = function(adj_in = double(1), num_in = double(1)) {
        N <- dim(num_in)[1]
        L <- dim(adj_in)[1]
        adj <- nimInteger(L)
        num <- nimInteger(N)
        for(i in 1:L)
            adj[i] <- ADbreak(adj_in[i])
        for(i in 1:N)
            num[i] <- ADbreak(num_in[i])
        numIslands <- 0L
        visited <- rep(0L, N)
        for(i in 1:N) {
            if(visited[i] == 0) {
                visited[i] <- 1
                numIslands <- numIslands + 1
                nNeighbors <- num[i]
                if(nNeighbors > 0) {
                    adjStartInd <- 1L
                    if(i > 1) adjStartInd <- adjStartInd + sum(num[1:(i-1)])
                    indToVisit <- nimInteger(nNeighbors)
                    indToVisit[1:nNeighbors] <- adj[adjStartInd:(adjStartInd+nNeighbors-1)]
                    lengthIndToVisit <- nNeighbors
                    l <- 1L
                    while(l <= lengthIndToVisit) {
                        nextInd <- indToVisit[l]
                        if(visited[nextInd] == 0) {
                            visited[nextInd] <- 1
                            newNneighbors <- num[nextInd]
                            if(newNneighbors > 0) {
                                newIndToVisit <- numeric(newNneighbors)
                                adjStartInd <- 1L
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
    },
    buildDerivs = list(run  = list(ignore = c("i")))
)


## Generates the \code{M} argument of the \code{dcar_proper} distribution
## 
## Generate the vector of conditional variances, as is required as the \code{M} argument of the \code{dcar_proper} distribution.  This is only used internally, and should never need to be invoked directly.
## 
## @author Daniel Turek
CAR_calcM <- nimbleFunction(
    name = 'CAR_calcM',
    run = function(num = double(1)) {
        N <- dim(num)[1]
        M <- rep(-1234, length = N)   ## Mi = -1234 indicates M was auto-generated, and region i has 0 neighbors
        for(i in 1:N) {
            if(num[i] > 0) {
                M[i] <- ADbreak(1/num[i])
            }
        }
        returnType(double(1))
        return(M)
    },
    buildDerivs = list(run  = list(ignore = c("i")))
)


## Generate the \code{C} argument of the \code{dcar_proper} distribution
##
## Generate sparse vector representation of the normalized adjacency matrix \eqn{C}, as is required as the \code{C} argument of the \code{dcar_proper} distribution.  This is only used internally, and should never need to be invoked directly.
## @author Daniel Turek
CAR_calcC <- nimbleFunction(
    name = 'CAR_calcC',
    run = function(adj = double(1), num_in = double(1)) {
        N <- dim(num_in)[1]
        L <- dim(adj)[1]
        num <- nimInteger(N)
        for(i in 1:N)
            num[i] <- ADbreak(num_in[i])
        C <- rep(0, length = L)
        count <- 1L
        for(i in 1:N) {
            if(num[i] > 0) {
                for(j in 1:num[i]) {
                    C[count] <- 1/num[i]   ##sqrt(M[i]/M[adj[count]])
                    count <- count + 1
                }
            }
        }
        if(count != L+1)   stop('dcar distribution internal error')
        returnType(double(1))
        return(C)
    },
    buildDerivs = list(run  = list(ignore = c("i", "j")))
)


## Generates the normalized adjacency matrix for \code{dcar_proper} distribution
## 
## Using the sparse representation of the normalized adjacency matrix, generate the full matrix.  This uses the \code{C}, \code{adj}, and \code{num} parameters of the \code{dcar_proper} distribution.
## 
## @author Daniel Turek
CAR_calcCmatrix <- nimbleFunction(
    name = 'CAR_calcCmatrix',
    run = function(C = double(1), adj_in = double(1), num_in = double(1)) {
        N <- dim(num_in)[1]
        L <- dim(adj_in)[1]
        adj <- nimInteger(L)
        num <- nimInteger(N)
        for(i in 1:L)
            adj[i] <- ADbreak(adj_in[i])
        for(i in 1:N)
            num[i] <- ADbreak(num_in[i])
        Cmatrix <- array(0, dim = c(N, N))
        count <- 1L
        for(i in 1:N) {
            if(num[i] > 0) {
                for(j in 1:num[i]) {
                    Cmatrix[i, adj[count]] <- ADbreak(C[count])
                    count <- count + 1
                }
            }
        }
        if(count != L+1)   stop('dcar distribution internal error')
        returnType(double(2))
        return(Cmatrix)
    },
    buildDerivs = list(run  = list(ignore = c("i", "j")))
)


#' Calculate bounds for the autocorrelation parameter of the \code{dcar_proper} distribution
#'
#' Calculate the lower and upper bounds for the \code{gamma} parameter of the \code{dcar_proper} distribution
#' 
#' @details Bounds for gamma are the inverse of the minimum and maximum eigenvalues of: \eqn{M^(-0.5) C M^(0.5)}.  The lower and upper bounds are returned in a numeric vector.
#' 
#' @param C vector of the same length as \code{adj}, giving the normalized weights associated with each pair of neighboring locations.
#' @param adj vector of indices of the adjacent locations (neighbors) of each spatial location.  This is a sparse representation of the full adjacency matrix.
#' @param num vector giving the number of neighboring locations of each spatial location, with length equal to the number of locations.
#' @param M vector giving the diagonal elements of the conditional variance matrix, with length equal to the number of locations.
#'
#' @return A numeric vector containing the bounds (minimum and maximum allowable values) for the \code{gamma} parameter of the \code{dcar_proper} distribution.
#' 
#' @seealso \code{\link{CAR-Proper}}, \code{\link{carMinBound}}, \code{\link{carMaxBound}}
#' 
#' @author Daniel Turek
#' @export
carBounds <- nimbleFunction(
    name = 'carBounds',
    run = function(C = double(1), adj = double(1), num = double(1), M = double(1)) {
        N <- dim(num)[1]
        L <- dim(adj)[1]
        Cmatrix <- CAR_calcCmatrix(C[1:L], adj[1:L], num[1:N])
        x <- diag(M^-0.5) %*% Cmatrix %*% diag(M^0.5)
        eigenvalues <- nimEigen(x, only.values = TRUE)$values
        lower <- 1 / max(eigenvalues)
        upper <- 1 / min(eigenvalues)
        bounds <- c(lower, upper)
        orderedBounds <- c(min(bounds), max(bounds))
        returnType(double(1))
        return(orderedBounds)
    }
)


#' Calculate the lower bound for the autocorrelation parameter of the \code{dcar_proper} distribution
#'
#' Calculate the lower bound for the \code{gamma} parameter of the \code{dcar_proper} distribution
#' 
#' @details Bounds for \code{gamma} are the inverse of the minimum and maximum eigenvalues of: \eqn{M^(-0.5) C M^(0.5)}.
#' 
#' @param C vector of the same length as \code{adj}, giving the normalized weights associated with each pair of neighboring locations.
#' @param adj vector of indices of the adjacent locations (neighbors) of each spatial location.  This is a sparse representation of the full adjacency matrix.
#' @param num vector giving the number of neighboring locations of each spatial location, with length equal to the number of locations.
#' @param M vector giving the diagonal elements of the conditional variance matrix, with length equal to the number of locations.
#' 
#' @return The lower bound (minimum allowable value) for the \code{gamma} parameter of the \code{dcar_proper} distribution.
#' 
#' @seealso \code{\link{CAR-Proper}}, \code{\link{carMaxBound}}, \code{\link{carBounds}}
#' 
#' @author Daniel Turek
#' @export
carMinBound <- nimbleFunction(
    name = 'carMinBound',
    run = function(C = double(1), adj = double(1), num = double(1), M = double(1)) {
        bounds <- carBounds(C, adj, num, M)
        returnType(double(0))
        return(bounds[1])
    }
)


#' Calculate the upper bound for the autocorrelation parameter of the \code{dcar_proper} distribution
#'
#' Calculate the upper bound for the \code{gamma} parameter of the \code{dcar_proper} distribution
#' 
#' @details Bounds for \code{gamma} are the inverse of the minimum and maximum eigenvalues of \eqn{M^(-0.5) C M^(0.5)}.
#' 
#' @param C vector of the same length as \code{adj}, giving the normalized weights associated with each pair of neighboring locations.
#' @param adj vector of indices of the adjacent locations (neighbors) of each spatial location.  This is a sparse representation of the full adjacency matrix.
#' @param num vector giving the number of neighboring locations of each spatial location, with length equal to the number of locations.
#' @param M vector giving the diagonal elements of the conditional variance matrix, with length equal to the number of locations.
#' 
#' @return The upper bound (maximum allowable value) for the \code{gamma} parameter of the \code{dcar_proper} distribution.
#' 
#' @seealso \code{\link{CAR-Proper}}, \code{\link{carMinBound}}, \code{\link{carBounds}}
#' 
#' @author Daniel Turek
#' @export
carMaxBound <- nimbleFunction(
    name = 'carMaxBound',
    run = function(C = double(1), adj = double(1), num = double(1), M = double(1)) {
        bounds <- carBounds(C, adj, num, M)
        returnType(double(0))
        return(bounds[2])
    }
)

## Calculate the eigenvalues of the normalized adjacency matrix with all equal weights
## 
## This function calculates the \code{evs} parameter values for the \code{dcar_proper} distribution, when \code{C} is inferred from the \code{adj} parameter as having all equal weights, and should never need to be invoked directly.
## 
## @author Daniel Turek

CAR_calcEVs2 <- nimbleFunction(
    name = 'CAR_calcEVs2',
    run = function(adj = double(1), num = double(1)) {
        N <- dim(num)[1]
        C <- CAR_calcC(adj, num)
        Cmatrix_out <- CAR_calcCmatrix(C, adj, num)
        ## Need to make sure Cmatrix not tracked or have issue with conversion of CppAD double in call to nimEigen.
        Cmatrix <- nimMatrix(nrow = N, ncol = N)
        for(i in 1:N)  
            for(j in 1:N)
                Cmatrix[i,j] <- ADbreak(Cmatrix_out[i,j])
        evs <- nimEigen(Cmatrix, only.values = TRUE)$values
        ## NaN eigenvalues occur when Cmatrix is singular.
        ## Need to make sure arg of `any_nan` and `is.nan` are not CppAD types.
        tmp <- nimNumeric(N)
        for(i in 1:N)
            tmp[i] <- ADbreak(evs[i])
        if(any_nan(tmp)) {
            for(i in 1:N) {
                tmpi <- tmp[i]  # o.w. lifted node is CppAD double.
                if(is.nan(tmpi)) evs[i] <- 0
            }
        }
        returnType(double(1))
        return(evs)
    },
    buildDerivs = list(run  = list(ignore = c("i", "j", "tmp", "tmpi", "Cmatrix")))
)


## Calculate the eigenvalues of the normalized adjacency matrix
## 
## This function calculates the \code{evs} parameter values for the \code{dcar_proper} distribution, and should never need to be invoked directly.
## 
## @author Daniel Turek
CAR_calcEVs3 <- nimbleFunction(
    name = 'CAR_calcEVs3',
    run = function(C = double(1), adj = double(1), num = double(1)) {
        N <- dim(num)[1]
        Cmatrix_out <- CAR_calcCmatrix(C, adj, num)
        Cmatrix <- nimMatrix(nrow = N, ncol = N)
        for(i in 1:N)  
            for(j in 1:N)
                Cmatrix[i,j] <- ADbreak(Cmatrix_out[i,j])
        evs <- nimEigen(Cmatrix, only.values = TRUE)$values
        ## NaN eigenvalues occur when Cmatrix is singular.
        tmp <- numNumeric(N)
        for(i in 1:N)
            tmp[i] <- ADbreak(evs[i])
        if(any_nan(tmp)) {
            for(i in 1:N) {
                tmpi <- tmp[i]  # o.w. lifted node is CppAD double.
                if(is.nan(tmpi)) evs[i] <- 0
            }
        }
        returnType(double(1))
        return(evs)
    },
    buildDerivs = list(run  = list(ignore = c("i", "tmp", "tmpi", "Cmatrix")))
)


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
        if(length(targetDCAR) != 1)                              stop('dcar distribution internal error')
        if(model$getDistribution(targetDCAR) != 'dcar_normal')   stop('dcar distribution internal error')
    },
    run = function() {
        if(island) return(0)
        priorMean <- getMean()
        priorSigma <- sqrt(1/getPrec())
        lp <- dnorm(model[[targetScalar]], priorMean, priorSigma, log = TRUE)
        returnType(double())
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
            if(island) return(0)
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
        targetDCARscalarComponents <- model$expandNodeNames(targetDCAR, returnScalarComponents = TRUE)
        targetIndex <- which(targetDCARscalarComponents == targetScalar)
        numNeighbors <- length(neighborCs)                        ## fix length-1 neighborCs
        island <- numNeighbors == 0
        neighborCs <- array(neighborCs, c(1, numNeighbors))       ## fix length-1 neighborCs
        neighborIndices <- array(0, c(1, numNeighbors))
        neighborIndices[1, ] <- match(neighborNodes, targetDCARscalarComponents)
        if(length(neighborNodes) != numNeighbors)                stop('dcar distribution internal error')
        if(length(targetDCAR) != 1)                              stop('dcar distribution internal error')
        if(model$getDistribution(targetDCAR) != 'dcar_proper')   stop('dcar distribution internal error')
        if(length(targetIndex) != 1)                             stop('dcar distribution internal error')
    },
    run = function() {
        if(Mi == -1234) return(0)   ## Mi = -1234 indicates M was auto-generated, and this region has 0 neighbors => conditional density = 0
        priorMean <- getMean()
        priorSigma <- sqrt(1/getPrec())
        lp <- dnorm(model[[targetScalar]], priorMean, priorSigma, log = TRUE)
        returnType(double())
        return(lp)
    },
    methods = list(
        getMean = function() {
            muValues <- model$getParam(targetDCAR, 'mu')
            if(island) return(muValues[targetIndex])
            gamma <- model$getParam(targetDCAR, 'gamma')
            neighborValues <- values(model, neighborNodes)
            mean <- muValues[targetIndex] + gamma * sum(neighborCs[1,1:numNeighbors] * (neighborValues - muValues[neighborIndices[1,1:numNeighbors]]))
            returnType(double())
            return(mean)
        },
        getPrec = function() {
            if(Mi == -1234) return(0)   ## Mi = -1234 indicates M was auto-generated, and this region has 0 neighbors => 1/Mi = num[i] = 0 => prec = 0
            prec <- model$getParam(targetDCAR, 'tau') / Mi
            returnType(double())
            return(prec)
        }
    )
)







