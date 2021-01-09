#' Placeholder for compareMCMCs
#'
#' This function has been moved to a separate package
#'
#' @param ... arguments
#'
#' @export
#'
compareMCMCs <- function(...)
    cat("This function now resides in a separate package. Please see https://github.com/nimble-dev/compareMCMCs.\n")

#' Placeholder for MCMCsuite
#'
#' This function has been moved to a separate package
#'
#' @param ... arguments
#'
#' @export
#'
MCMCsuite <- function(...)
    cat("This function now resides in a separate package. Please see https://github.com/nimble-dev/compareMCMCs.\n")


#' Placeholder for buildAuxiliaryFilter
#'
#' This function has been moved to the `nimbleSMC` package
#'
#' @param ... arguments
#'
#' @export
#'
buildAuxiliaryFilter <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildBootstrapFilter
#'
#' This function has been moved to the `nimbleSMC` package
#'
#' @param ... arguments
#'
#' @export
#'
buildBootstrapFilter <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildEnsembleKF
#'
#' This function has been moved to the `nimbleSMC` package
#'
#' @param ... arguments
#'
#' @export
#'
buildEnsembleKF <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildIteratedFilter2
#'
#' This function has been moved to the `nimbleSMC` package
#'
#' @param ... arguments
#'
#' @export
#'
buildIteratedFilter2 <- function(...)
    cat("NIMBLE's sequential Monte Carlo functionality, including this function, now resides in the 'nimbleSMC' package.\n")

#' Placeholder for buildLiuWestFilter
#'
#' This function has been moved to the `nimbleSMC` package
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
            ## Chris suggests:
            ## tmp <- forwardsolve(L, I)
            ## ans <- crossprod(tmp)
        }
        returnType(double(2))
        return(ans)
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
            ## Chris suggests:
            ## tmp <- forwardsolve(L, I)
            ## ans <- crossprod(tmp)
        }
        returnType(double(2))
        return(ans)
    }
)

