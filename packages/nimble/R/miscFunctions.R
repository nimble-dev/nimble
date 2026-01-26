#' Placeholder for buildLaplace
#'
#' This function has been moved to the `nimbleQuad` package.
#'
#' @param ... arguments
#'
#' @export
#'
buildLaplace <- function(...)
    cat("NIMBLE's Laplace/AGHQ functionality, including this function, now resides in the 'nimbleQuad' package.\n")

#' Placeholder for buildAGHQ
#'
#' This function has been moved to the `nimbleQuad` package.
#'
#' @param ... arguments
#'
#' @export
#'
buildAGHQ <- function(...)
    cat("NIMBLE's Laplace/AGHQ functionality, including this function, now resides in the 'nimbleQuad' package.\n")

#' Placeholder for runLaplace
#'
#' This function has been moved to the `nimbleQuad` package.
#'
#' @param ... arguments
#'
#' @export
#'
runLaplace <- function(...)
    cat("NIMBLE's Laplace/AGHQ functionality, including this function, now resides in the 'nimbleQuad' package.\n")

#' Placeholder for runAGHQ
#'
#' This function has been moved to the `nimbleQuad` package.
#'
#' @param ... arguments
#'
#' @export
#'
runAGHQ <- function(...)
    cat("NIMBLE's Laplace/AGHQ functionality, including this function, now resides in the 'nimbleQuad' package.\n")

#' Placeholder for summaryLaplace
#'
#' This function has been moved to the `nimbleQuad` package.
#'
#' @param ... arguments
#'
#' @export
#'
summaryLaplace <- function(...)
    cat("NIMBLE's Laplace/AGHQ functionality, including this function, now resides in the 'nimbleQuad' package.\n")

#' Placeholder for summaryAGHQ
#'
#' This function has been moved to the `nimbleQuad` package.
#'
#' @param ... arguments
#'
#' @export
#'
summaryAGHQ <- function(...)
    cat("NIMBLE's Laplace/AGHQ functionality, including this function, now resides in the 'nimbleQuad' package.\n")


#' Placeholder for buildAuxiliaryFilter
#'
#' This function has been moved to the `nimbleSMC` package.
#'
#' @param ... arguments
#'
#' @export
#'
buildAuxiliaryFilter <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildBootstrapFilter
#'
#' This function has been moved to the `nimbleSMC` package.
#'
#' @param ... arguments
#'
#' @export
#'
buildBootstrapFilter <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildEnsembleKF
#'
#' This function has been moved to the `nimbleSMC` package.
#'
#' @param ... arguments
#'
#' @export
#'
buildEnsembleKF <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildIteratedFilter2
#'
#' This function has been moved to the `nimbleSMC` package.
#'
#' @param ... arguments
#'
#' @export
#'
buildIteratedFilter2 <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildLiuWestFilter
#'
#' This function has been moved to the `nimbleSMC` package.
#'
#' @param ... arguments
#'
#' @export
#'
buildLiuWestFilter <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

## used in altParams for dmnorm
## this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dmnormAltParams <- nimbleFunction(
    name = 'calc_dmnormAltParams',
    run = function(cholesky = double(2), prec_param = double(), return_prec = double()) {
        if(prec_param == return_prec) {
            ans <- t(cholesky) %*% cholesky
        } else {
            I <- diag(dim(cholesky)[1])
            ans <- backsolve(cholesky, forwardsolve(t(cholesky), I))
        }
        returnType(double(2))
        return(ans)
    }
)

calc_dmnorm_inv_ld_AltParams <- nimbleFunction(
    name = 'calc_dmnorm_inv_ld_AltParams',
    run = function(mat = double(2), inv_ld = double(1), prec_param = double(), return_prec = double()) {
        ## No need for inversion as desired result is either in `mat` from
        ## original parameter provided or inverse is in `inv_ld`.
        if(prec_param == return_prec)
            return(mat)
        
        n <- sqrt(length(inv_ld)-1)
        nsq <- n * n
        ans <- matrix(inv_ld[1:nsq], nrow = n, ncol = n)
        if(n > 1) {  # Fill lower triangle as PDinverse_logdet only guarantees upper.
            for(i in 2:n)
                for(j in 1:(i-1))
                    ans[i,j] <- ans[j,i]
        }
        return(ans)
        returnType(double(2))
    }
)

## This is used in conjugacy definition for ddirch, to calculate 'contribution'
## terms from dcat dependents.
calc_dcatConjugacyContributions <- nimbleFunction(
    name = 'calc_dcatConjugacyContributions',
    run = function(ncat = double(0), value = double(0)) {
        ans <- numeric(ncat)
        ans[value] <- 1
        return(ans)
        returnType(double(1))
    }
)

## used in conjugacy definition for dmnorm, to calculate 'contribution' terms;
## formerly avoided unnecessary matrix multiplications, when 'coeff' is identity matrix by numerical computation
## now we do via code processing to determine the realized link
calc_dmnormConjugacyContributions <- nimbleFunction(
    name = 'calc_dmnormConjugacyContributions',
    run = function(coeff = double(2), prec = double(2), vec = double(1), order = double(), use_coeff = double()) {
        if(use_coeff == 0) {  ## identity, additive
            if(order == 1) ans <- prec %*% asCol(vec)
            if(order == 2) ans <- prec
        } else {
            if(order == 1) ans <- t(coeff) %*% (prec %*% asCol(vec))
            if(order == 2) ans <- t(coeff) %*% prec %*% coeff
        }
        return(ans)
        returnType(double(2))
    }
)


## used in altParams for dwish and dinvwish
## this needs to be sourced after nimbleFunction() is defined, so can't be done in distributions_inputList.R
calc_dwishAltParams <- nimbleFunction(
    name = 'calc_dwishAltParams',
    run = function(cholesky = double(2), scale_param = double(), return_scale = double()) {
        if(scale_param == return_scale) {
            ans <- t(cholesky) %*% cholesky
        } else {
            I <- diag(dim(cholesky)[1])
            ans <- backsolve(cholesky, forwardsolve(t(cholesky), I))
        }
        returnType(double(2))
        return(ans)
    }
)



## used in conjugacy definition for dgamma, to calculate 'contribution_shape' term:
calc_dcar_normalConjugacyContributionShape <- nimbleFunction(
    name = 'calc_dcar_normalConjugacyContributionShape',
    run = function(num = double(1), c = double()) {
        N <- length(num)
        ans <- (N-c)/2
        return(ans)
        returnType(double())
    }
)



## used in conjugacy definition for dgamma, to calculate 'contribution_rate' term:
calc_dcar_normalConjugacyContributionRate <- nimbleFunction(
    name = 'calc_dcar_normalConjugacyContributionRate',
    run = function(adj = double(1), weights = double(1), num = double(1), value = double(1)) {
        N <- length(num)
        L <- length(weights)
        count <- 1L
        ans <- 0
        for(i in 1:N) {
            if(num[i] > 0) {
                for(j in 1:num[i]) {
                    ## this prevents "double-counting" each pair of neighboring nodes:
                    if(i < adj[count]) {
                        ans <- ans + weights[count] * (value[i] - value[adj[count]])^2
                    }
                    count <- count + 1
                }
            }
        }
        if(count != L+1)   stop('gamma-CAR conjugacy calculation internal error')
        return(ans)
        returnType(double())
    }
)
