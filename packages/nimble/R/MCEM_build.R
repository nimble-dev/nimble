#' @export
build_MCEM_expectation_noAD <- nimbleFunction(
  setup = function(model, mvSamples, burnInDefault, paramNodes,
                   latentNodes, calcNodes, calcNodesOther,
                   useTransform, transformer,
                   lengthParams, lengthParamsTrans,
                   derivsDelta = 0.0001) {
    logLik <- as.numeric(0)
    lastParamsTrans <- rep(-Inf, lengthParamsTrans)
    if(lengthParamsTrans == 1) {
      lastParamsTrans <- c(-Inf, -1) # reset will be called first anyway
    }
    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
    if(length(paramDeps) > 0) {
      keep_paramDeps <- logical(length(paramDeps))
      for(i in seq_along(paramDeps)) {
        if(any(paramDeps[i] == calcNodes)) keep_paramDeps[i] <- FALSE
        else {
          nextDeps <- model$getDependencies(paramDeps[i])
          keep_paramDeps[i] <- any(nextDeps %in% calcNodes)
        }
      }
      paramDeps <- paramDeps[keep_paramDeps]
    }

    allCalcNodes <- model$topologicallySortNodes(c(paramDeps, calcNodes))
    useOther <- length(calcNodesOther) > 0

    burnIn <- burnInDefault
  },
  methods = list(
    reset = function() {
      lastParamsTrans <<- rep(-Inf, lengthParamsTrans)
    },
    set_burnIn = function(new_burnIn = integer()) {
      burnIn <<- new_burnIn
    },
    ##
    inverseTransform = function(tT = double(1)) {
      return(transformer$inverseTransform(tT))
      returnType(double(1))
    },
    llh = function(paramsTrans = double(1)) {
      nSamples <- getsize(mvSamples)
      sum_LL <- 0
      if(useTransform) {
        params <- inverseTransform(paramsTrans)
      } else {
        params <- paramsTrans
      }
      values(model, paramNodes) <<- params
      for(i in (burnIn+1):nSamples){
        nimCopy(from = mvSamples, to = model, nodes = latentNodes, row = i)
        sample_LL <- model$calculate(allCalcNodes)
        if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
          stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')' to see if node 'y[3]' is the cause of the problem). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
        sum_LL <- sum_LL + sample_LL
      }
      logLik <<- sum_LL / (nSamples - burnIn)
      if(is.nan(logLik))
        logLik <<- -Inf

      if(useOther) {
        other_llh_derivs <- model$calculate(calcNodesOther)
        logLik <<- logLik + other_llh_derivs
      }

      lastParamsTrans <<- paramsTrans
      return(logLik)
      returnType(double())
    },
    grad_llh = function(paramsTrans = double(1)) {
      stop("Should never call grad_llh within MCEM if AD is not being used.")
      returnType(double(1))
      return(rep(NA, length(paramsTrans)))
    },
    llh_sample = function(paramsTrans = double(1)) {
      nSamples <- getsize(mvSamples)
      if(useTransform) {
        params <- transformer$inverseTransform(paramsTrans)
      } else {
        params <- paramsTrans
      }
      values(model, paramNodes) <<- params
      #
      logLiks <- numeric(length = nSamples - burnIn, init = FALSE)
      for(i in (burnIn+1):nSamples){
        nimCopy(from = mvSamples, to = model, nodes = latentNodes, row = i)
        logLiks[i - burnIn] <- model$calculate(allCalcNodes)
      }
      if(useOther) {
        otherLogLik <- model$calculate(calcNodesOther)
        logLiks <- logLiks + otherLogLik
      }
      return(logLiks)
      returnType(double(1))
    },
    fisherInfrmtn = function(params = double(1), trans = logical(0, default = 0),
                             returnTrans = integer(0, default = -1),
                             atMLE = logical(0, default = 0)) {
      # trans == TRUE means that params provided is already transformed
      nSamples <- getsize(mvSamples)
      if(returnTrans == -1) returnTrans <- trans
      if(trans) {
        paramsTrans <- params
        if(useTransform) paramsActual <- transformer$inverseTransform(params)
        else paramsActual <- params
      } else {
        paramsActual <- params
        if(useTransform) paramsTrans <- transformer$transform(params)
        else paramsTrans <- params
      }
      values(model, paramNodes) <<- paramsActual
      workParamsActual <- paramsActual
      delta <- derivsDelta
      halfdelta <- 0.5*delta
      twodelta <- 2*delta
      deltasq <- delta*delta
      meanGrad <- numeric(value = 0, length = lengthParams)
      grad <- numeric(lengthParams)
      hessian <- matrix(0, nrow = lengthParams, ncol = lengthParams)
      meanGradGrad <-  matrix(0, nrow = lengthParams, ncol = lengthParams)
      meanHessian <- matrix(0, nrow = lengthParams, ncol = lengthParams)
      fxplus <- numeric(value = 0, length = lengthParams)
      fxminus <- numeric(value = 0, length = lengthParams)
      fxygrid <- matrix(0, nrow = 2, ncol = 2)
      for(iSamp in (burnIn+1):nSamples){
        nimCopy(from = mvSamples, to = model, nodes = latentNodes, row = iSamp)
        values(model, paramNodes) <<- workParamsActual
        fx <- model$calculate(allCalcNodes)
        for(ix in 1:lengthParams) {
          workParamsActual[ix] <- paramsActual[ix] - delta
          values(model, paramNodes) <<- workParamsActual
          fxminus[ix] <- model$calculate(allCalcNodes)
          workParamsActual[ix] <- paramsActual[ix] + delta
          values(model, paramNodes) <<- workParamsActual
          fxplus[ix] <- model$calculate(allCalcNodes)
          # Fill meanGrad[ix]
          grad[ix] <- (fxplus[ix] - fxminus[ix])/twodelta
          hessian[ix, ix] <- (fxplus[ix] -2*fx + fxminus[ix]) / deltasq
          workParamsActual[ix] <- paramsActual[ix]
        }
        if(!atMLE) meanGrad <- meanGrad + grad
        meanGradGrad <- meanGradGrad + (grad%*%t(grad))#gradMatrix %*% t(gradMatrix)

        if(lengthParams > 1) {
          for(ix in 2:lengthParams) {
            for(iy in 1:(ix-1)) {
              workParamsActual[ix] <- paramsActual[ix] - halfdelta
              workParamsActual[iy] <- paramsActual[iy] - halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[1,1] <- model$calculate(allCalcNodes)
              workParamsActual[ix] <- paramsActual[ix] + halfdelta
              workParamsActual[iy] <- paramsActual[iy] - halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[2,1] <- model$calculate(allCalcNodes)
              workParamsActual[ix] <- paramsActual[ix] - halfdelta
              workParamsActual[iy] <- paramsActual[iy] + halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[1,2] <- model$calculate(allCalcNodes)
              workParamsActual[ix] <- paramsActual[ix] + halfdelta
              workParamsActual[iy] <- paramsActual[iy] + halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[2,2] <- model$calculate(allCalcNodes)
              hessian[ix, iy] <- ((fxygrid[2,2] - fxygrid[2,1]) - (fxygrid[1,2] - fxygrid[1,1]))/deltasq
              hessian[iy, ix] <- hessian[ix, iy]
              workParamsActual[ix] <- paramsActual[ix]
              workParamsActual[iy] <- paramsActual[iy]
            }
          }
        }
        meanHessian <- meanHessian + hessian
      }
      if(!atMLE) meanGrad <- meanGrad / (nSamples - burnIn)
      meanHessian <- meanHessian / (nSamples - burnIn)
      meanGradGrad <- meanGradGrad / (nSamples - burnIn)
      returnType(double(2))
      returnMat <- -meanHessian - meanGradGrad
      if(!atMLE) returnMat <- returnMat + (meanGrad%*%t(meanGrad))
      if(useOther) {
        values(model, paramNodes) <<- workParamsActual
        fx <- model$calculate(calcNodesOther)
        for(ix in 1:lengthParams) {
          workParamsActual[ix] <- paramsActual[ix] - delta
          values(model, paramNodes) <<- workParamsActual
          fxminus[ix] <- model$calculate(calcNodesOther)
          workParamsActual[ix] <- paramsActual[ix] + delta
          values(model, paramNodes) <<- workParamsActual
          fxplus[ix] <- model$calculate(calcNodesOther)
          # Fill meanGrad[ix]
          # grad[ix] <- (fxplus[ix] - fxminus[ix])/twodelta
          hessian[ix, ix] <- (fxplus[ix] -2*fx + fxminus[ix]) / deltasq
          workParamsActual[ix] <- paramsActual[ix]
        }
        if(lengthParams > 1) {
          for(ix in 2:lengthParams) {
            for(iy in 1:(ix-1)) {
              workParamsActual[ix] <- paramsActual[ix] - halfdelta
              workParamsActual[iy] <- paramsActual[iy] - halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[1,1] <- model$calculate(calcNodesOther)
              workParamsActual[ix] <- paramsActual[ix] + halfdelta
              workParamsActual[iy] <- paramsActual[iy] - halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[2,1] <- model$calculate(calcNodesOther)
              workParamsActual[ix] <- paramsActual[ix] - halfdelta
              workParamsActual[iy] <- paramsActual[iy] + halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[1,2] <- model$calculate(calcNodesOther)
              workParamsActual[ix] <- paramsActual[ix] + halfdelta
              workParamsActual[iy] <- paramsActual[iy] + halfdelta
              values(model, paramNodes) <<- workParamsActual
              fxygrid[2,2] <- model$calculate(calcNodesOther)
              hessian[ix, iy] <- ((fxygrid[2,2] - fxygrid[2,1]) - (fxygrid[1,2] - fxygrid[1,1]))/deltasq
              hessian[iy, ix] <- hessian[ix, iy]
              workParamsActual[ix] <- paramsActual[ix]
              workParamsActual[iy] <- paramsActual[iy]
            }
          }
        }
        returnMat <- returnMat - hessian

      }
      if(returnTrans & useTransform) {
        paramsJacobian <- matrix(0, nrow = lengthParams, ncol = lengthParamsTrans)
        workParamsTrans <- paramsTrans
        for(ix in 1:lengthParamsTrans) {
          workParamsTrans[ix] <- paramsTrans[ix] + halfdelta
          paramsPlus <- inverseTransform(workParamsTrans)
          workParamsTrans[ix] <- paramsTrans[ix] - halfdelta
          paramsMinus <- inverseTransform(workParamsTrans)
          paramsJacobian[,ix] <- (paramsPlus - paramsMinus)/delta
          workParamsTrans[ix] <- paramsTrans[ix]
        }
        returnMat <- t(paramsJacobian) %*% returnMat %*% paramsJacobian
      }
      return(returnMat)
    }
  ),
  buildDerivs = list(inverseTransform = list())
)

#' @export
build_MCEM_expectation <- nimbleFunction(
  setup = function(model, mvSamples, burnInDefault, paramNodes,
                   latentNodes, calcNodes, calcNodesOther,
                   useTransform, transformer,
                   lengthParams, lengthParamsTrans) {
    logLik <- as.numeric(0)
    gradientLik <- rep(as.numeric(0), lengthParamsTrans)
    lastParamsTrans <- rep(-Inf, lengthParamsTrans)
    if(lengthParamsTrans == 1) {
      lastParamsTrans <- c(-Inf, -1) # reset will be called first anyway
      gradientLik <- as.numeric(c(0, -1))
    }
    transform_wrt <- as.integer(1:lengthParamsTrans)

    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
    if(length(paramDeps) > 0) {
      keep_paramDeps <- logical(length(paramDeps))
      for(i in seq_along(paramDeps)) {
        if(any(paramDeps[i] == calcNodes)) keep_paramDeps[i] <- FALSE
        else {
          nextDeps <- model$getDependencies(paramDeps[i])
          keep_paramDeps[i] <- any(nextDeps %in% calcNodes)
        }
      }
      paramDeps <- paramDeps[keep_paramDeps]
    }

    allCalcNodes <- model$topologicallySortNodes(c(paramDeps, calcNodes))
    useOther <- length(calcNodesOther) > 0

    burnIn <- burnInDefault
  },
  methods = list(
    reset = function() {
      lastParamsTrans <<- rep(-Inf, lengthParamsTrans)
      gradientLik <<- rep(-Inf, lengthParamsTrans) # perhaps unnecessary
    },
    set_burnIn = function(new_burnIn = integer()) {
      burnIn <<- new_burnIn
    },
    ##
    inverseTransform = function(tT = double(1)) {
      return(transformer$inverseTransform(tT))
      returnType(double(1))
    },
    llh = function(paramsTrans = double(1)) {
      nSamples <- getsize(mvSamples)
      sum_LL <- 0
      sum_jacobian <- matrix(value = 0, nrow = 1, ncol = lengthParams)
      if(useTransform) {
        paramsDerivs <- nimDerivs(inverseTransform(paramsTrans),
                                 wrt = transform_wrt, order = c(0,1))
        params <- paramsDerivs$value
      } else {
        params <- paramsTrans
      }
      values(model, paramNodes) <<- params
      for(i in (burnIn+1):nSamples){
        nimCopy(from = mvSamples, to = model, nodes = latentNodes, row = i)
        sample_derivs <- nimDerivs(model$calculate(allCalcNodes), wrt = paramNodes, order = c(0, 1))
        sample_LL <- sample_derivs$value[1]
        if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
          stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')' to see if node 'y[3]' is the cause of the problem). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
        sum_LL <- sum_LL + sample_LL
        sum_jacobian <- sum_jacobian + sample_derivs$jacobian
      }
      logLik <<- sum_LL / (nSamples - burnIn)
      if(is.nan(logLik))
        logLik <<- -Inf
      sum_jacobian <- sum_jacobian / (nSamples - burnIn) # now this is mean of jacobian

      if(useOther) {
        other_llh_derivs <- nimDerivs(model$calculate(calcNodesOther), wrt = paramNodes, order = c(0,1))
        logLik <<- logLik + other_llh_derivs$value[1]
        sum_jacobian <- sum_jacobian + other_llh_derivs$jacobian
      }
      if(useTransform)
        sum_jacobian <- sum_jacobian %*% paramsDerivs$jacobian # should be 1xlength(params) %*% length(params)xlength(paramsTrans)
      gradientLik <<- sum_jacobian[1,]

      lastParamsTrans <<- paramsTrans
      return(logLik)
      returnType(double())
    },
    grad_llh = function(paramsTrans = double(1)) {
      if(all(paramsTrans == lastParamsTrans))
        return(gradientLik)
      llh(paramsTrans)
      return(gradientLik)
      returnType(double(1))
    },
    llh_sample = function(paramsTrans = double(1)) {
      nSamples <- getsize(mvSamples)
      if(useTransform) {
        params <- transformer$inverseTransform(paramsTrans)
      } else {
        params <- paramsTrans
      }
      values(model, paramNodes) <<- params
      #
      logLiks <- numeric(length = nSamples - burnIn, init = FALSE)
      for(i in (burnIn+1):nSamples){
        nimCopy(from = mvSamples, to = model, nodes = latentNodes, row = i)
        logLiks[i - burnIn] <- model$calculate(allCalcNodes)
      }
      if(useOther) {
        otherLogLik <- model$calculate(calcNodesOther)
        logLiks <- logLiks + otherLogLik
      }
      return(logLiks)
      returnType(double(1))
    },
    fisherInfrmtn = function(params = double(1), trans = logical(0, default = 0),
                             returnTrans = integer(0, default = -1),
                             atMLE = logical(0, default = 0)) {
      # trans == TRUE means that params provided is already transformed
      nSamples <- getsize(mvSamples)
      if(returnTrans == -1) returnTrans <- trans
      if(trans) {
        paramsTrans <- params
        if(useTransform) paramsActual <- transformer$inverseTransform(params)
        else paramsActual <- params
      } else {
        paramsActual <- params
        if(useTransform) paramsTrans <- transformer$transform(params)
        else paramsTrans <- params
      }
      values(model, paramNodes) <<- paramsActual

      meanGrad <- numeric(lengthParams)
      meanGradGrad <-  matrix(0, nrow = lengthParams, ncol = lengthParams)
      meanHessian <- matrix(0, nrow = lengthParams, ncol = lengthParams)
      for(i in (burnIn+1):nSamples){
        nimCopy(from = mvSamples, to = model, nodes = latentNodes, row = i)
        sample_derivs <- nimDerivs(model$calculate(allCalcNodes), wrt = paramNodes, order = 1:2)
        grad <- sample_derivs$jacobian
        if(!atMLE) meanGrad <- meanGrad + grad[1,]
        meanHessian <- meanHessian + sample_derivs$hessian[,,1]
        meanGradGrad <- meanGradGrad + (t(grad)%*%grad)
      }
      if(!atMLE) meanGrad <- meanGrad / (nSamples - burnIn)
      meanHessian <- meanHessian / (nSamples - burnIn)
      meanGradGrad <- meanGradGrad / (nSamples - burnIn)
      returnType(double(2))
      returnMat <- -meanHessian - meanGradGrad
      if(!atMLE) returnMat <- returnMat + (meanGrad%*%t(meanGrad))
      if(useOther) {
        calcOtherDerivs <- nimDerivs(model$calculate(calcNodesOther), wrt = paramNodes, order = 2)
        returnMat <- returnMat - calcOtherDerivs$hessian[,,1]
      }
      if(returnTrans & useTransform) {
        paramsDerivs <- nimDerivs(inverseTransform(paramsTrans),
                                 wrt = transform_wrt, order = 1)
        returnMat <- t(paramsDerivs$jacobian) %*% returnMat %*% paramsDerivs$jacobian
      }
      return(returnMat)
    }
  ),
  buildDerivs = list(inverseTransform = list())
)

#' @export
MCEM_mcse <- function(samples, m) {
  # in practice, samples will be a vector of logLik differences
  if(!require(mcmcse)) stop("Must install package mcmcse.")
  ans <- mcse(as.matrix(samples, ncol = 1),
              #size = ceiling(min(1000, (m/20))), ## backward compatible. m here is m-burnin below
              method = "obm",
              # r = 1
              )$se[1]
  ans
}

#' @export
R_MCEM_mcse <- nimbleRcall(function(samples = double(1), m = integer()) {},
                           returnType = double(),
                           Rfun = 'MCEM_mcse')

#' Builds an MCEM algorithm for a given NIMBLE model
#'
#' Takes a NIMBLE model (with some missing data, aka random effects or latent
#' state nodes) and builds a Monte Carlo Expectation Maximization (MCEM)
#' algorithm for maximum likelihood estimation. The user can specify which
#' latent nodes are to be integrated out in the E-Step, or default choices will
#' be made based on model structure. All other stochastic non-data nodes will be
#' maximized over. The E-step is done with a sample from a nimble MCMC
#' algorithm. The M-step is done by a call t \code{optim}.
#'
#' @param paramNodes a character vector of names of parameter nodes in the
#'   model; defaults are provided by \code{\link{setupMargNodes}}.
#'   Alternatively, \code{paramNodes} can be a list in the format returned by
#'   \code{setupMargNodes}, in which case \code{latentNodes}, \code{calcNodes},
#'   and \code{calcNodesOther} are not needed (and will be ignored).
#' @param latentNodes a character vector of names of unobserved
#'   (latent) nodes to marginalize (sum or integrate) over; defaults are provided by
#'   \code{\link{setupMargNodes}}.
#' @param calcNodes a character vector of names of nodes for calculating
#'   components of the full-data likelihood that invole \code{latentNodes};
#'   defaults are provided by \code{\link{setupMargNodes}}. There may be
#'   deterministic nodes between \code{paramNodes} and \code{calcNodes}. These
#'   will be included in calculations automatically and thus do not need to be
#'   included in \code{calcNodes} (but there is no problem if they are).
#' @param calcNodesOther a character vector of names of nodes for calculating
#'   terms in the log-likelihood that do not depend on any \code{latentNodes},
#'   and thus are not part of the marginalization, but should be included for
#'   purposes of finding the MLE. This defaults to stochastic nodes that depend
#'   on \code{paramNodes} but are not part of and do not depend on
#'   \code{latentNodes}. There may be deterministic nodes between
#'   \code{paramNodes} and \code{calcNodesOther}. These will be included in
#'   calculations automatically and thus do not need to be included in
#'   \code{calcNodesOther} (but there is no problem if they are).
#' @param control a named list for providing additional settings used in MCEM.
#'   See \code{control} section below.
#'
#' @details \code{buildMCEM} is a nimbleFunction that creates an MCEM algorithm
#'   for a model and choices (perhaps default) of nodes in different roles in
#'   the model. The MCEM can then be compiled for fast execution with a compiled model.
#'
#' Note that \code{buildMCEM} was re-written for nimble version 1.2.0 and is not
#' backward-compatible with previous versions.
#'
#' Denote data by Y, latent states (or missing data) by X, and parameters by T.
#' MCEM works by the following steps, starting from some T:
#'
#' \enumerate{
#'
#' \item Draw a sample of size M from P(X | Y, T) using MCMC.
#'
#' \item Update T to be the maximizer of E[log P(X, Y | T)] where the
#' expectation is approximated as a Monte Carlo average over the sample from step(1)
#'
#' \item Repeat until converged.
#'
#' }
#'
#' The default version of MCEM is the ascent-based MCEM of Caffo et al. (2015).
#' This attempts to update M when necessary to ensure that step 2 really moves
#' uphill given that it is maximizing a Monte Carlo approximation and could
#' accidently move downhill on the real surface of interest due to Monte Carlo
#' error. The tuning parameters include \code{alpha}, \code{beta}, \code{gamma},
#' \code{C}, \code{tol} (tolerance), and others.
#'
#' If the model supports derivatives via nimble's automatic differentiation (AD)
#' (and \code{buildDerivs=TRUE} in \code{nimbleModel}), the maximization step
#' can use gradients from AD. You must manually set \code{useDerivs=FALSE} in
#' the control list if derivatives aren't supported or if you don't want to use
#' them.
#'
#' After maximization in step 2, we estimate the Monte Carlo standard error of
#' the uphill movement. If the standardized uphill step is bigger than 0 with
#' Type I error rate alpha, the iteration is accepted and the algorithm
#' continues. Otherwise, it is not certain that step 2 really moved uphill due
#' to Monte Carlo error, so the MCMC sample size is incremented by a fixed
#' factor (e.g. 0.33 or 0.5, called \code{Mfactor} in the control list), the
#' additional samples are added by continuing step 1, and step 2 is tried again.
#' If the Monte Carlo noise still overwhelms the magnitude of uphill movement,
#' the sample size is increased again, and so on. alpha should be >0 and <=0.5.
#' A larger value than usually used for inference is recommended so that there
#' is an easy threshold to determine uphill movement, which avoids increasing M
#' prematurely.
#'
#' Convergence is determined in a similar way. After a definite move uphill, we
#' determine if the uphill increment is less than \code{tol}, with Type I error
#' rate gamma. (But if \code{M} hits a maximum value, the convergence criterion
#' changes. See below.)
#'
#' beta is used to help set M to a minimal level based on previous iterations.
#' This is a desired Type II error rate, assuming an uphill move and standard
#' error based on the previous iteration. Set \code{adjustM=FALSE} in the
#' control list if you don't want this behavior.
#'
#' There are some additional controls on convergence for practical purposes. Set
#' \code{C} in the control list to be the number of times the convergence
#' criterion mut be satisfied in order to actually stop. E.g setting \code{C=2}
#' means there will always be a restart after the first convergence.
#'
#' One problem that can occur with ascent-based MCEM is that the final iteration
#' can be very slow if M must become very large to satisfy the convergence
#' criterion. Indeed, if the algorithm starts near the MLE, this can occur. Set
#' \code{maxM} in the control list to set the MCMC sample size that should never
#' be exceeded.
#'
#' If \code{M==maxM}, a softer convergence criterion is used. This second
#' convergence criterion is to stop if we can't be sure we moved uphill using
#' Type I error rate delta. This is a soft criterion because for small delta,
#' Type II errors will be common (e.g. if we really did move uphill but can't be
#' sure from the Monte Carlo sample), allowing the algorithm to terminate. One
#' can continue the algorithm from where it stopped, so it is helpful to not
#' have it get stuck when having a very hard time with the first (stricter)
#' convergence criterion.
#'
#' All of alpha, beta, delta, and gamma are utilized based on asymptotic
#' arguments but in practice must be chosen heuristically. In other words, their
#' theoretical basis does not exactly yield practical advice on good choices for
#' efficiency and accuracy, so some trial and error will be needed.
#'
#' It can also be helpful to set a minimum and maximum of allowed iterations (of
#' steps 1 and 2 above). Setting \code{minIter>1} in the control list can
#' sometimes help avoid a false convergence on the first iteration by forcing at
#' least one more iteration. Setting \code{maxIter} provides a failsafe on a
#' stuck run.
#'
#' If you don't want the ascent-based method at all and simply want to run a set
#' of iterations, set \code{ascent=FALSE} in the control list. This will use the
#' second (softer) convergence criterion.
#'
#' Parameters to be maximized will by default be handled in an unconstrained
#' parameter space, transformed if necessary by an
#' \code{\link{parameterTransform}} object. In that case, the default
#' \code{\link{optim}} method will be "BFGS" and can can be changed by setting
#' \code{optimMehod} in the control list. Set \code{useTransform=FALSE} in the
#' control list if you don't want the parameters transformed. In that case the
#' default \code{optimMethod} will be "L-BFGS-B" if there are any actual
#' constraints, and you can provide a list of \code{boxConstraints} in the
#' control list. (Constraints may be determined by priors written in the model
#' for parameters, even though they priors play no other role in MLE. E.g.
#' \code{sigma ~ halfflat()} indicates \code{sigma > 0}).
#'
#' Most of the control list elements can be over-ridden when calling the
#' \code{findMLE} method. The \code{findMLE} argument \code{continue=TRUE}
#' results in attempting to continue the algorithm where the previous call
#' finished, including whatever settings where in use.
#'
#' See \code{\link{setupMargNodes}} (which is called with the given arguments
#' for \code{paramNodes}, \code{calcNodes}, and \code{calcNodesOther}; and with
#' \code{allowDiscreteLatent=TRUE}, \code{randomEffectsNodes=latentNodes}, and
#' \code{check=check}) for more about how the different groups of nodes are
#' determined. In general, you can provide none, one, or more of the different
#' kinds of nodes and \code{setupMargNodes} will try to determine the others in
#' a sensible way. However, note that this cannot work for all ways of writing a
#' model. One key example is that if random (latent) nodes are written as
#' top-level nodes (e.g. following \code{N(0,1)}), they appear structurally to
#' be parameters and you must tell \code{buildMCEM} that they are
#' \code{latentNodes}. The various "Nodes" arguments will all be passed through
#' \code{model$expandNodeNames}, allowing for example simply "x" to be provided
#' when there are many nodes within "x".
#'
#' Estimating the Monte Carlo standard error of the uphill step is not trivial
#' because the sample was obtained by MCMC and so likely is autocorrelated. This
#' is done by calling whatever function in R's global environment is called
#' "MCEM_mcse", which is required to take two arguments: \code{samples} (which
#' will be a vector of the differences in log(P(Y, X | T)) between the new and
#' old values of T, across the sample of X) and \code{m}, the sample size. It
#' must return an estimate of the standard error of the mean of the sample.
#' NIMBLE provides a default version (exported from the package namespace),
#' which calls \code{mcmcse::mcse} with method "obm". Simply provide a different
#' function with this name in your R session to over-ride NIMBLE's default.
#'
#' @section control list details
#'
#' The control list accepts the following named elements:
#'
#' \itemize{
#'
#' \item \code{initM} initial MCMC sample size, \code{M}. Default=1000.
#'
#' \item \code{Mfactor} Factor by which to increase MCMC sample size when step 2
#' results in noise overwhelming the uphill movement. The new \code{M} will be
#' \code{1+Mfactor)*M} (rounded up). \code{Mfactor} is \code{1/k} of Caffo et
#' al. (2015). Default=1/3.
#'
#' \item \code{maxM} Maximum allowed value of \code{M} (see above). Default=\code{initM*20}.
#'
#' \item \code{burnin} Number of burn-in iterations for the MCMC in step 1. Note
#' that the initial states of one MCMC will be the last states from the previous
#' MCMC, so they will often be good initial values after multiple iterations. Default=500.
#'
#' \item \code{thin} Thinning interval for the MCMC in step 1. Default=1.
#'
#' \item \code{alpha} Type I error rate for determining when step 2 has moved
#' uphill. See above. Default=0.25.
#'
#' \item \code{beta} Used for determining a minimal value of $M$ based on
#' previous iteration, if \code{adjustM} is \code{TRUE}. \code{beta} is a desired Type
#' II error rate for determining uphill moves. Default=0.25.
#'
#' \item \code{delta} Type I error rate for the soft convergence approach
#' (second approach above). Default=0.25.
#'
#' \item \code{gamma} Type I error rate for determining when step 2 has moved
#' less than \code{tol} uphill, in which case ascent-based convergence is
#' achieved (first approach above). Default=0.05.
#'
#' \item \code{buffer} A small amount added to lower box constraints and
#' substracted from upper box constraints for all parameters, relevant only if
#' \code{useTransform=FALSE} and some parameters do have \code{boxConstraints}
#' set or have bounds that can be determined from the model.
#'
#' \item \code{tol} Ascent-based convergence tolerance. Default=0.001.
#'
#' \item \code{ascent} Logical to determine whether to use the ascent-based
#' method of Caffo et al. Default=TRUE.
#'
#' \item \code{C} Number of convergences required to actually stop the
#' algorithm. Default = 1.
#'
#' \item \code{maxIter} Maximum number of MCEM iterations to run.
#'
#' \item \code{minIter} Minimum number of MCEM iterations to run.
#'
#' \item \code{adjustM} Logical indicating whether to see if M needs to be
#' increased based on statistical power argument in each iteration (using
#' \code{beta}). Default=TRUE.
#'
#' \item \code{verbose} Logical indicating whether verbose output is desired.
#' Default=TRUE.
#'
#' \item \code{MCMCprogressBar} Logical indicating whether MCMC progress bars
#' should be shown for every iteration of step 1. This argument is passed to
#' \code{configureMCMC}, or to \code{config} if provided. Default=TRUE.
#'
#' \item \code{derivsDelta} If AD derivatives are not used, then the method
#' \code{vcov} must use finite difference derivatives to implement the method of
#' Louis (1982). The finite differences will be \code{delta} or \code{delta/2}
#' for various steps. This is the same for all dimensions. Default=0.0001.
#'
#' \item \code{mcmcControl} This is passed to \code{configureMCMC}, or
#' \code{config} if provided, as the \code{control} argument. i.e.
#' \code{control=mcmcControl}.
#'
#' \item \code{boxContrainst} List of box constraints for the nodes that will be
#' maximized over, only relevant if \code{useTransform=FALSE} and
#' \code{forceNoConstraints=FALSE} (and ignored otherwise). Each constraint is a
#' list in which the first element is a character vector of node names to which
#' the constraint applies and the second element is a vector giving the lower
#' and upper limits. Limits of \code{-Inf} or \code{Inf} are allowed. Any nodes
#' that are not given constrains will have their constraints automatically
#' determined by NIMBLE. See above. Default=list().
#'
#' \item \code{forceNoConstraints} Logical indicating whether to force ignoring
#' constraints even if they might be necessary. Default=FALSE.
#'
#' \item \code{useTransform} Logical indicating whether to use a parameter
#' transformation (see \code{\link{parameterTransform}}) to create an unbounded
#' parameter space for the paramNodes. This allows unconstrained maximization
#' algorithms to be used. Default=TRUE.
#'
#' \item \code{check} Logical passed as the \code{check} argument to
#' \code{\link{setupMargNodes}}.
#'
#' \item \code{useDerivs} Logical indicating whether to use AD. If TRUE, the
#' model must have been build with `buildDerivs=TRUE`. It is not automatically
#' determined from the model whether derivatives are supported.
#'
#' \item \code{config} Optional function to create the MCMC configuration used
#' for step 1. If missing, the MCMC configuration is created by
#'
#' \code{configureMCMC(model, nodes = latentNodes,
#'                     monitors = latentNodes, thin = thinDefault,
#'                     control = mcmcControl, print = FALSE)}.
#'
#' If provided, the MCMC configuration is created by the same call with
#' \code{config} instead of \code{configureMCMC}.
#'
#' \section Methods in the returned algorithm
#'
#' The object returned by \code{buildMCEM} is a nimbleFunction object with the following methods
#'
#' \itemize{
#'
#' \item \code{findMLE} is the main method of interest, launching the MCEM
#' algorithm. It takes the following arguments:
#'
#'   \itemize{
#'
#'    \item \code{pStart}. Vector of initial parameter values. If omitted, the
#'    values currently in the model object are used.
#'
#'    \item \code{returnTrans}. Logical indicating whether to return parameters
#'    in the transformed space, if a parameter transformation is in use. Default=FALSE.
#'
#'    \item \code{continue}. Logical indicating whether to continue the MCEM
#'    from where the last call stopped. In addition, if TRUE, any other control
#'    setting provided in the last call will be used again. If FALSE, all
#'    control settings are reset to the values provided when \code{buildMCEM}
#'    was called. Any control settings provided in the same call as
#'    \code{continue=TRUE} will over-ride these behaviors and be used in the
#'    continued run.
#'
#'    \item All run-time control settings available in the \code{control} list
#'    for \code{buildMCEM} (except for \code{buffer}, \code{boxConstraints},
#'    \code{forceNoConstraints}, \code{useTransform}, and \code{useDerivs}) are
#'    accepted as individual arguments to over-ride the values provided in the
#'    \code{control} list.
#'
#'   }
#'
#' \code{findMLE} returns on object of class \code{optimResultNimbleList} the
#' with the results of the final optimization of step 2. The \code{par} element
#' of this list is the vector of maximum likelihood (MLE) parameters.
#'
#' \item \code{vcov} computes the approximate variance-covariance matrix of the MLE using
#' the method of Louis (1982). It takes the following arguments:
#'
#'    \itemize{
#'
#'      \item \code{params}. Vector of parameters at which to compute the
#'      Hessian matrix used to obtain the \code{vcov} result. Typically this
#'      will be \code{MLE$par}, if \code{MLE} is the output of \code{findMLE}.
#'
#'      \item \code{trans}. Logical indicating whether \code{params} is on the
#'      transformed prameter scale, if a parameter transformation is in use.
#'      Typically this should be the same as the \code{returnTrans} argument to
#'      \code{findMLE}. Default=FALSE.
#'
#'      \item \code{returnTrans}. Logical indicting whether the \code{vcov}
#'      result should be for the transformed parameter space. Default matches
#'      \code{trans}.
#'
#'      \item \code{M}. Number of MCMC samples to obtain if
#'      \code{resetSamples=TRUE}. Default is the final value of \code{M} from
#'      the last call to \code{findMLE}. It can be helpful to increase \code{M}
#'      to obtain a more accurate \code{vcov} result (i.e. with less Monte Carlo
#'      noise).
#'
#'      \item \code{resetSamples}. Logical indicating whether to generate a new
#'      MCMC sample from P(X | Y, T), where T is \code{params}. If FALSE, the
#'      last sample from \code{findMLE} will be used. If MLE convergence was
#'      reasonable, this sample can be used. However, if the last MCEM step made
#'      a big move in parameter space (e.g. if convergence was not achieved),
#'      the last MCMC sample may not be accurate for obtaining \code{vcov}. Default=FALSE.
#'
#'      \item \code{atMLE}. Logical indicating whether you believe the
#'      \code{params} represents the MLE. If TRUE, one part of the computation
#'      will be skipped because it is expected to be 0 at the MLE. If there are
#'      parts of the model that are not connected to the latent nodes, i.e. of
#'      \code{calcNodesOther} is not empty, then \code{atMLE} will be ignored
#'      and set to FALSE. Default=FALSE. It is not really worth using TRUE
#'      unless you are confident and the time saving is meaningful, which is not
#'      very likely. In other words, this argument is provided for technical
#'      completeness.
#'
#'    }
#'
#' \code{vcov} return a matrix that is the inverse of the negative Hessian of
#' the log likelihood surface, i.e. the usual asymptotic approximation of the
#' parameter variance-covariance matrix.
#'
#' \item \code{doMCMC}. This method runs the MCMC to sample from P(X | Y, T).
#' One does not need to call this, as it is called via the MCEM algorithm in
#' \code{findMLE}. This method is provided for users who want to use the MCMC
#' for latent states directly. Samples should be retrieved by
#' \code{as.matrix(compiled_MCEM$mvSamples)}. This method takes the following arguments:
#'
#'   \itemize{
#'
#'     \item \code{M}. MCMC sample size.
#'
#'     \item \code{thin}. MCMC thinning interval.
#'
#'     \item \code{reset}. Logical indicating whether to reset the MCMC (passed
#'     to the MCMC \code{run} method as \code{reset}).
#'
#'   }
#'
#' \item \code{transform} and \code{inverseTransform}. Convert a parameter
#' vector to an unconstrained parameter space and vice-versa, if
#' \code{useTransform=TRUE} in the call to \code{buildDerivs}.
#'
#' \item \code{resetControls}. Reset all control arguments to the values
#' provided in the call to \code{buildMCEM}. The user does not normally need to
#' call this.
#'
#' }
#'
#' @author Perry de Valpine, Clifford Anderson-Bergman and Nicholas Michaud
#' @export
#'
#' @references
#'
#' Caffo, Brian S., Wolfgang Jank, and Galin L. Jones (2005). Ascent-based Monte Carlo expectation-maximization.  \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 67(2), 235-251.
#'
#'  Louis, Thomas A  (1982). Finding the Observed Information Matrix When Using the EM Algorithm. \emph{Journal of the Royal Statistical Society. Series B (Statistical Methodology)}, 44(2), 226-233.
#'
#' @examples
#' \dontrun{
#' pumpCode <- nimbleCode({
#'  for (i in 1:N){
#'      theta[i] ~ dgamma(alpha,beta);
#'      lambda[i] <- theta[i]*t[i];
#'      x[i] ~ dpois(lambda[i])
#'  }
#'  alpha ~ dexp(1.0);
#'  beta ~ dgamma(0.1,1.0);
#' })
#'
#' pumpConsts <- list(N = 10,
#'               t = c(94.3, 15.7, 62.9, 126, 5.24,
#'                 31.4, 1.05, 1.05, 2.1, 10.5))
#'
#' pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
#'
#' pumpInits <- list(alpha = 1, beta = 1,
#'              theta = rep(0.1, pumpConsts$N))
#' pumpModel <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
#'                          data = pumpData, inits = pumpInits,
#'                          buildDerivs=TRUE)
#'
#' pumpMCEM <- buildMCEM(model = pumpModel)
#'
#' CpumpModel <- compileNimble(pumpModel)
#'
#' CpumpMCEM <- compileNimble(pumpMCEM, project=pumpModel)
#'
#' MLE <- CpumpMCEM$findMLE()
#' vcov <- CpumpMCEM$vcov(MLE$par)
#'
#' }
buildMCEM <- nimbleFunction(
  name = "MCEM",
  setup = function(model, paramNodes, latentNodes, calcNodes, calcNodesOther,
                   control = list(), ...) {
    dotsList <- list(...)
    if(any(c("burnIn", "mcmcControl", "boxConstraints",
             "buffer", "alpha", "beta", "gamma", "C",
             "numReps", "forceNoConstraints", "verbose") %in% names(dotsList)))
      stop("From the arguments provided, it looks like you are trying to use the old version of buildMCEM.",
           " buildMCEM was rewritten for nimble 1.2.0.")
    # Extract control elements or set defaults
    initMuse <- initMdefault <- extractControlElement(control, 'initM', 1000)
    MfactorUse <- MfactorDefault <- extractControlElement(control, 'Mfactor', 1/3)
    maxMuse <- maxMdefault <- extractControlElement(control, 'maxM', initMuse * 20)
    burnInUse <- burnInDefault <- extractControlElement(control, 'burnin', 500)
    thinUse <- thinDefault <- extractControlElement(control, 'thin', 1)
    alphaUse <- alphaDefault <- extractControlElement(control, 'alpha', 0.25)
    betaUse <- betaDefault <- extractControlElement(control, 'beta', 0.25)
    deltaUse <- deltaDefault <- extractControlElement(control, 'delta', 0.25)
    gammaUse <- gammaDefault <- extractControlElement(control, 'gamma', 0.05)
    bufferUse <- bufferDefault <- extractControlElement(control, 'buffer', 1e-6)
    tolUse <- tolDefault <- extractControlElement(control, 'tol', 0.001)
    ascentUse <- ascentDefault <- extractControlElement(control, 'ascent', TRUE)
    Cuse <- Cdefault <- extractControlElement(control, 'C', 1)
    maxIterUse <- maxIterDefault <- extractControlElement(control, 'maxIter', 50)
    minIterUse <- minIterDefault <- extractControlElement(control, 'minIter', 1)
    adjustMuse <- adjustMdefault <- extractControlElement(control, 'adjustM', TRUE)
    verboseUse <- verboseDefault <- extractControlElement(control, 'verbose', getNimbleOption("verbose"))
    MCMCprogressBarUse <- MCMCprogressBarDefault <- extractControlElement(control, 'MCMCprogressBar', getNimbleOption("MCMCprogressBar"))
    derivsDeltaUse <- derivsDeltaDefault <- extractControlElement(control, 'derivsDelta', 0.0001)

    if(verboseUse) message("  [Note] You may need to rebuild and recompile a model for each call to buildMCEM.")

    mcmcControl <- extractControlElement(control, 'mcmcControl', list(adaptInterval = 100))
    # optimControl: see below
    # optimMethod: see below
    boxConstraints <- extractControlElement(control, 'boxConstraints', list())
    forceNoConstraints <- extractControlElement(control, 'forceNoConstraints', FALSE)
    useTransform <- extractControlElement(control, 'transform', TRUE)
    check <- extractControlElement(control, 'check', TRUE)
    useDerivs <- extractControlElement(control, 'useDerivs', TRUE)
    config <- extractControlElement(control, 'config', NULL)

    # check for valid control parameters
    if(useTransform && length(boxConstraints)>0) {
      if(verboseUse)
        message("boxConstraints will be ignored because transform is TRUE.")
      boxConstraints <- list()
    }
    if(forceNoConstraints && length(boxConstraints)>0) {
      if(verboseUse)
        message("boxConstraints will be ignored because forceNoConstraints is TRUE.")
      boxConstraints <- list()
    }
    ## if(bufferDefault == 0)
    ##   warning("'buffer' is zero. This can cause problems if the likelihood",
    ##           " function is degenerate on boundary")
    if(bufferDefault < 0)
      stop("'buffer' must be non-negative.")

    # Set up variables used in updating MCMC sample size
    zAlpha <- qnorm(alphaUse, 0, 1, lower.tail=FALSE)
    zBeta <- qnorm(betaUse, 0, 1, lower.tail=FALSE)
    zDelta <- qnorm(deltaUse, 0, 1, lower.tail=FALSE)
    zGamma <- qnorm(gammaUse, 0, 1, lower.tail=FALSE)

    # Set up information on paramNodes, latentNodes, calcNodes, and calcNodesOther
    MargNodes <- NULL
    if(!missing(paramNodes)) {
      if(is.list(paramNodes)) {
        # The user called setupMargNodes and provided a list of that format to paramNodes.
        MargNodes <- paramNodes
        allStochNonDataNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
        if(length(setdiff(MargNodes$randomEffectsNodes, allStochNonDataNodes)) != 0)
          stop('Some latentNodes (possibly given as paramNodes$randomEffectsNodes) provided',
               ' were not found in model.')
      }
    }
    if(is.null(MargNodes)) {
      MargNodes <- setupMargNodes(model = model, paramNodes = paramNodes,
                                  randomEffectsNodes = latentNodes,
                                  calcNodes = calcNodes,
                                  calcNodesOther = calcNodesOther,
                                  split = FALSE,
                                  check = check,
                                  allowDiscreteLatent = TRUE)
    }
    paramNodes <- MargNodes$paramNodes
    latentNodes <- MargNodes$randomEffectsNodes
    calcNodes <- MargNodes$calcNodes
    calcNodesOther <- MargNodes$calcNodesOther
    num_calcNodesOther <- length(calcNodesOther)
    useOther <- num_calcNodesOther > 0
    if(any(model$isDiscrete(paramNodes)))
      stop(paste0("MCEM cannot optimize over discrete top-level parameters. ",
                  "The following top-level parameters in your model are discrete: ",
                  paste0(paramNodes[model$isDiscrete(paramNodes)], collapse = ', ')))
    latentNodesScalar <- model$expandNodeNames(latentNodes, returnScalarComponents = TRUE)
    lengthLatentNodes <- length(latentNodesScalar)
    if(lengthLatentNodes == 0)
      stop("There are no latent nodes for MCEM to use.")
    paramNodesScalar <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
    lengthParams <- length(paramNodesScalar)
    if(lengthLatentNodes == 0)
      stop("There are no parameter nodes for MCEM to use.")

    # optimMethod
    if(!is.null(control$optimMethod) &&
         (control$optimMethod %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B"))){
      optimMethod <- control$optimMethod
    }
    else optimMethod <- "default" # will changed later to BFGS or L-BFGS-B (if there are bounds)

    # Deal with constraints
    low_limits = rep(-Inf, lengthParams )
    hi_limits  = rep(Inf,  lengthParams )
    needConstraints <- optimMethod == "L-BFGS-B" || optimMethod == "default"
    needConstraints <- needConstraints && (!useTransform) && (!forceNoConstraints)
    if(length(boxConstraints)>0 && !needConstraints)
      if(verboseUse)
        warning("boxConstraints will be ignored because optimMethod must be 'L-BFGS-B'",
                " (or unspecified), and transform must be FALSE, to use boxConstraints.")
    if(needConstraints) {
      for(i in seq_along(paramNodesScalar)) {
        low_limits[i] = model$getBound(paramNodesScalar[i], 'lower') + abs(bufferUse)
        hi_limits[i]  = model$getBound(paramNodesScalar[i], 'upper') - abs(bufferUse)
      }

      constraintNames = list()
      for(i in seq_along(boxConstraints) )
        constraintNames[[i]] = model$expandNodeNames(boxConstraints[[i]][[1]])
      for(i in seq_along(constraintNames) ) {
        limits = boxConstraints[[i]][[2]]
        inds = which(paramNodesScalar %in% constraintNames[[i]])
        if(length(inds) == 0)
          stop(paste("warning: provided a constraint for nodes", constraintNames[[i]],
                     ", but those nodes do not exist in the model!"))
        tooLowNodes <- which(limits[1] + abs(bufferUse) < low_limits[inds])
        tooHighNodes <- which(limits[2] - abs(bufferUse) > hi_limits[inds])
        if(length(tooLowNodes) > 0)
          warning(paste0("User-specified lower bound for ", constraintNames[[i]][tooLowNodes],
                         " is below lower bound detected by NIMBLE.  "))
        if(length(tooHighNodes) > 0)
          warning(paste0("User-specified upper bound for ", constraintNames[[i]][tooHighNodes],
                         " is above upper bound detected by NIMBLE.  "))
        low_limits[inds] = limits[1] + abs(bufferUse)
        hi_limits[inds] = limits[2] - abs(bufferUse)
      }
      if(any(low_limits>=hi_limits))
        stop('lower limits greater than or equal to upper limits!')

      if(optimMethod == "default") {
        optimMethod <- "L-BFGS-B"
        if(all(low_limits == -Inf) && all(hi_limits == Inf))
          optimMethod <- "BFGS"
      }
    } else {
      if(optimMethod == "default")
        optimMethod <- "BFGS"
    }
    useConstraints <- optimMethod == "L-BFGS-B"

    ## All derivative calls are of the form nimDerivs(model$calculate(...), ...),
    ## so these derivative info pieces are not needed.
    ## I am leaving the code here as a reminder in case we need to change.
    ##
    ## # Derivative info
    ## if(length(calcNodesOther)) {
    ##   otherLogLik_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = paramNodes, calcNodes = calcNodesOther)
    ##   otherLogLik_updateNodes   <- otherLogLik_derivsInfo$updateNodes
    ##   otherLogLik_constantNodes <- otherLogLik_derivsInfo$constantNodes
    ## }
    ## else { ## calcNodesOther is empty
    ##   otherLogLik_updateNodes   <- character(0)
    ##   otherLogLik_constantNodes <- character(0)
    ## }

    # optimControl
    optimControl <- nimOptimDefaultControl()
    optimControlArgNames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha",
                              "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    if(!is.null(control$optimControl)){
      validNames <- intersect(names(control$optimControl), optimControlArgNames)
      numValidNames <- length(validNames)
      if(numValidNames > 0){
        for(i in 1:numValidNames){
          optimControl[[validNames[i]]] <- control$optimControl[[validNames[i]]]
        }
      }
    }
    optimControl$fnscale <- -1

    # parameter transformation
    if(useTransform)
      transformer <- parameterTransform(model, paramNodes)
    else
      transformer <- parameterTransform(model, paramNodes[1]) # still need an object for compilation. should allow character() but doesn't
    lengthParamsTrans <- transformer$getTransformedLength()

    nimbleVerbose <- getNimbleOption("verbose")
    nimbleOptions(verbose=FALSE)
    # set up MCMC
    if(!is.null(config)) {
      if(!is.function(config)) {
        stop("If config is provided via control list, it must be a ",
             "function that returns an MCMC configuration (similar to what configureMCMC returns).")
        mcmc_Latent_Conf <- config(model, nodes = latentNodes, monitors = latentNodes,
                                   thin = thinDefault, control = mcmcControl, print = FALSE)
        if(!inherits(mcmc_Latent_Conf, "MCMCconf"))
          stop("If config is provided via control list, it must be a ",
               "function that returns an MCMC configuration (an object of class MCMCconf, ",
               "similar to what configureMCMC returns).")
      }
    } else {
      mcmc_Latent_Conf <- configureMCMC(model, nodes = latentNodes,
                                        monitors = latentNodes,
                                        thin = thinDefault,
                                        control = mcmcControl, print = FALSE)
    }
    mcmc_Latent <- buildMCMC(mcmc_Latent_Conf)
    mvSamples <- mcmc_Latent$mvSamples
    setupOutputs(mvSamples)
    nimbleOptions(verbose = nimbleVerbose)

    if(useDerivs) {
      E <- build_MCEM_expectation(model, mvSamples, burnInDefault, paramNodes,
                                latentNodes, calcNodes, calcNodesOther,
                                useTransform, transformer,
                                lengthParams, lengthParamsTrans)
    } else {
      E <- build_MCEM_expectation_noAD(model, mvSamples, burnInDefault, paramNodes,
                                latentNodes, calcNodes, calcNodesOther,
                                useTransform, transformer,
                                lengthParams, lengthParamsTrans,
                                derivsDeltaUse)
    }
    one_time_fixes_done <- FALSE
    if(length(hi_limits)==1) hi_limits <- c(hi_limits, -1)
    if(length(low_limits)==1) low_limits <- c(low_limits, -1)
    lastParamsTrans <- rep(-Inf, lengthParamsTrans)
    if(length(lastParamsTrans)==1) lastParamsTrans <- c(lastParamsTrans, -1)
    finalM <- 0L
    lastDiffQ <- 0
  },
  methods = list(
    fix_one_vec = function(x = double(1)) {
      if(length(x) == 2) {
        if(x[2] == -1) {
          ans <- numeric(length = 1, value = x[1])
          return(ans)
        }
      }
      return(x)
      returnType(double(1))
    },
    do_one_time_fixes = function() {
      if(!one_time_fixes_done) {
        hi_limits <<- fix_one_vec(hi_limits)
        low_limits <<- fix_one_vec(low_limits)
        lastParamsTrans <<- fix_one_vec(lastParamsTrans)
        one_time_fixes_done <<- TRUE
      }
    },
    vcov = function(params = double(1), trans = logical(0, default = 0),
                    returnTrans = integer(0, default = -1),
                    M = integer(0, default = -1), resetSamples = logical(0, default = FALSE),
                    atMLE = logical(0, default = FALSE)) {
      if(atMLE & useOther) {
        atMLE <- FALSE # atMLE only works if the params are at the MLE in the case of no calcNodesOther (useOther is FALSE)
        # This is tricky enough that we will silently revert it and document that atMLE is ignored if calcNodesOther is FALSE
        # The only point of atMLE is to save a minor bit of computation, and it might be only rarely worth worrying about.
      }
      if(returnTrans == -1) {
        returnTrans <- trans
      }
      if(trans) {
        if(length(params) != lengthParamsTrans) {
          print("  [Warning] For vcov with trans=TRUE, params should be length ", lengthParamsTrans, " but is length ", length(params), ".")
          stop()
        }
        paramsTrans <- params
        paramsActual <- inverseTransform(params)
      } else {
        if(length(params) != lengthParams) {
          print("  [Warning] For vcov, params should be length ", lengthParams, " but is length ", length(params), ".")
          stop()
        }
        paramsActual <- params
        paramsTrans <- transform(params)
      }
      if(resetSamples) { # default to using existing samples | (!all(paramsTrans == lastParamsTrans))) {
        if(M == -1) M <- finalM
        lastParamsTrans <<- paramsTrans
        values(model, paramNodes) <<- paramsActual
        doMCMC(M, TRUE)
      }
      FI <- E$fisherInfrmtn(params, trans, returnTrans, atMLE)
      vc <- inverse(FI)
      return(vc)
      returnType(double(2))
    },
    doMCMC = function(M = integer(), thin = integer(), reset = logical(0, default = TRUE)) {
      mcmc_Latent$run(M, thin = thin, reset = reset, progressBar = MCMCprogressBarUse)
    },
    transform = function(params = double(1)) {
      if(useTransform)
        ans <- transformer$transform(params)
      else ans <- params
      return(ans)
      returnType(double(1))
    },
    inverseTransform = function(paramsTrans = double(1)) {
      if(useTransform)
        ans <- transformer$inverseTransform(paramsTrans)
      else ans <- paramsTrans
      return(ans)
      returnType(double(1))
    },
    resetControls = function() {
      initMuse <<- initMdefault
      finalM <<- initMdefault
      MfactorUse <<- MfactorDefault
      maxMuse <<- maxMdefault
      burnInUse <<- burnInDefault
      thinUse <<- thinDefault
      alphaUse <<- alphaDefault
      betaUse <<- betaDefault
      deltaUse <<- deltaDefault
      gammaUse <<- gammaDefault
#      bufferUse <<- bufferDefault
      tolUse <<- tolDefault
      ascentUse <<- ascentDefault
      Cuse <<- Cdefault
      maxIterUse <<- maxIterDefault
      minIterUse <<- minIterDefault
      adjustMuse <<- adjustMdefault
      verboseUse <<- verboseDefault
      MCMCprogressBarUse <<- MCMCprogressBarDefault
      lastDiffQ <<- 0
    },
    findMLE = function(pStart  = double(1, default = Inf),
                       initM = integer(0, default=-1),
                       Mfactor = double(0, default=-1),
                       maxM = integer(0, default=-1),
                       burnin = integer(0, default=-1),
                       thin = integer(0, default=-1),
                       alpha = double(0, default=-1),
                       beta = double(0, default=-1),
                       delta = double(0, default=-1),
                       gamma = double(0, default=-1),
                       #                 buffer = double(0, default=-1),
                       tol = double(0, default=-1),
                       ascent = integer(0, default = -1),
                       C = integer(0, default=-1),
                       maxIter = integer(0, default=-1),
                       minIter = integer(0, default=-1),
                       adjustM = integer(0, default=-1),
                       verbose = integer(0, default=-1),
                       MCMCprogressBar = integer(0, default=-1),
                       returnTrans = logical(0, default = 0),
                       continue = logical(0, default=FALSE)
                       ) {
      do_one_time_fixes()

      if(any(abs(pStart) == Inf)) pStart <- values(model, paramNodes)
      if(length(pStart) != lengthParams) {
        print("  [Warning] For findMLE, pStart should be length ", lengthParams, " but is length ", length(pStart), ".")
        stop()
      }
      ## In case parameter nodes are not properly initialized
      if(any_na(pStart) | any_nan(pStart) | any(abs(pStart)==Inf))
        stop("  [Warning] For findMLE, pStart has some invalid values.")

      if(!continue) resetControls()
      else
        if(initM==-1) initMuse <<- finalM # different from other tuners here
      if(initM!=-1) initMuse <<- initM
      if(Mfactor!=-1) MfactorUse <<- Mfactor
      if(maxM!=-1) maxMuse <<- maxM
      if(burnin!=-1) burnInUse <<- burnin
      if(thin!=-1) thinUse <<- thin
      if(alpha!=-1) alphaUse <<- alpha
      if(beta!=-1) betaUse <<- beta
      if(delta!=-1) deltaUse <<- delta
      if(gamma!=-1) gammaUse <<- gamma
      #    if(buffer!=-1) bufferUse <<- buffer
      if(tol!=-1) tolUse <<- tol
      if(ascent!=-1) ascentUse <<- (ascent != 0)
      if(C!=-1) Cuse <<- C
      if(maxIter!=-1) maxIterUse <<- maxIter
      if(minIter!=-1) minIterUse <<- minIter
      if(adjustM!=-1) adjustMuse <<- (adjustM != 0)
      if(verbose!=-1) verboseUse <<- (verbose != 0)
      if(MCMCprogressBar!=-1) MCMCprogressBarUse <<- (MCMCprogressBar != 0)

      if(burnInUse >= initMuse)
        stop('mcem quitting: burnIn > initM value')
      E$set_burnIn(burnInUse)
      if(!continue) doMCMC(0, TRUE) # To get valid initial values
      params <- pStart

      if(!continue) {
        if(useConstraints) {
          for(i in 1:lengthParams ) {  # check that initial values satisfy constraints
            if((low_limits[i] == -Inf) & (hi_limits[i] < Inf)){
              if(params[i] > hi_limits[i]){
                params[i] <- hi_limits[i] - 1
              }
            }
            else if((hi_limits[i] == Inf) & (low_limits[i] > -Inf)){
              if(params[i] < low_limits[i]){
                params[i] <- low_limits[i] + 1
              }
            }
            else if((low_limits[i] > -Inf) & (hi_limits[i] < Inf)){
              if(!((params[i] >= low_limits[i]) & (params[i] <= hi_limits[i]))){
                params[i] = (low_limits[i] + hi_limits[i])/2
              }
            }
          }
          values(model, paramNodes) <<- params
        }
      }
      m <- initMuse
      finalM <<- m
      numberOfTimesFlat <- 0 # number of times result does not appear to move uphill
      looksFlat <- FALSE
      #endCrit <- C+1 #ensure that first iteration runs
      varDeltaQ <-0 #use initM as m value for first step
      diffQ <- 1 # any nonzero value can be used here, gets overwritten quickly in algo
      if(continue) diffQ <- lastDiffQ
      itNum <- 0
      zAlpha <<- qnorm(alphaUse, 0, 1, lower.tail=FALSE)
      zBeta <<- qnorm(betaUse, 0, 1, lower.tail=FALSE)
      zDelta <<- qnorm(deltaUse, 0, 1, lower.tail=FALSE)
      zGamma <<- qnorm(gammaUse, 0, 1, lower.tail=FALSE)
      hitMaxM <- FALSE
      paramsTrans <- transform(params)
      while(((numberOfTimesFlat < Cuse) & (itNum < maxIterUse)) |
              (itNum < minIterUse)) {  # || itNum <= 1) { #endCrit > C){ # We forced at least 2 iterations. Why?
                checkInterrupt()
                itNum <- itNum + 1
                acceptThisIter <- FALSE
                #starting sample size calculation for this iteration
                if((itNum > 1) & adjustMuse)
                  m <- burnInUse + ceiling(max(m - burnInUse, varDeltaQ*((zAlpha + zBeta)^2)/((diffQ)^2))) # Arguably the second argument to the max should be multipled by thinUse since it is from iid theory
                if(m >= maxMuse) {
                  m <- maxMuse
                  hitMaxM <- TRUE
                }
                lastParamsTrans <<- paramsTrans # records params used for the latent state MCMC sample
                doMCMC(m, thin = thin, reset = TRUE) # initial mcmc run of size m
                #    paramsPrev <- params  #store previous params value
                paramsTransPrev <- paramsTrans
                ## if(itNum == 1) {
                ##   sample_logLiks_new <- E$llh_sample(paramsTransPrev)
                ## }
                ## sample_logLiks_prev <- sample_logLiks_new
                while(!acceptThisIter){
                  E$reset()
                  if(useDerivs) {
                    if(optimMethod == "L-BFGS-B") {
                      optimOutput = optim(par = paramsTrans, fn = E$llh, gr = E$grad_llh,
                                          control = optimControl, method = 'L-BFGS-B',
                                          lower = low_limits, upper = hi_limits)
                    } else {
                      optimOutput = optim(par = paramsTrans, fn = E$llh, gr = E$grad_llh,
                                          control =  optimControl, method = optimMethod)
                    }
                  } else {
                    if(optimMethod == "L-BFGS-B") {
                      optimOutput = optim(par = paramsTrans, fn = E$llh,
                                          control =  optimControl, method = 'L-BFGS-B',
                                          lower = low_limits, upper = hi_limits)
                    } else {
                      optimOutput = optim(par = paramsTrans, fn = E$llh,
                                          control =  optimControl, method = optimMethod)
                    }
                  }
                  paramsTrans = optimOutput$par
                  # In the future, we could revise this to append llh_sample values when M is being extended.
                  # For now, we simply recalculate them all. This should be minor compared to cost of optim.
                  sample_logLiks_prev <- E$llh_sample(paramsTransPrev)
                  sample_logLiks_new <- E$llh_sample(paramsTrans) # Side effect: this puts the latest params into the model
                  # make nimbleRcall or other solution.
                  sdDeltaQ <- R_MCEM_mcse(sample_logLiks_new - sample_logLiks_prev, m-burnInUse)
                  varDeltaQ <- sdDeltaQ*sdDeltaQ
                  ##        varDeltaQ <- cvarCalc$run(m, params, paramsPrev)
                  ## sdDeltaQ <- sqrt(varDeltaQ) #asymptotic std. error
                  diffQ <- mean(sample_logLiks_new) - mean(sample_logLiks_prev)
                  ## Should have: mean(sample_logLiks_new) == mcse_result$est == optimOutput$value
                  ## diffQ <- cCalc_E_llk$run(params, paramsPrev, 1)
                  lastDiffQ <<- diffQ
                  # criterion for deciding we converged because it really looks flat
                  applyAscentRules <- ascentUse & (!hitMaxM)
                  if(applyAscentRules) {
                    # cat("applying ascent rules\n")
                    # criterion for increasing m
                    movedUphill <- (diffQ - zAlpha*sdDeltaQ) >= 0 # We reject H0:DeltaQ==0 (one-sided to H0:DeltaQ > 0) with (high) Type I error rate alpha. i.e. we believe we true diffQ > 0, so moved uphill, and we have a low threshold for belief. That means if we don't reject, there's too much noise (since we are guaranteed to move uphill, really).
                    # cat("uphill crit: ",(diffQ - zAlpha*sdDeltaQ), " ", movedUphill, "\n")
                    if((!movedUphill) & (m<maxMuse)) { #swamped by mc error
                      if(verboseUse) cat("Monte Carlo error too big: increasing MCMC sample size.\n")
                      mAdd <- ceiling(MfactorUse * (m-burnInUse)) #from section 2.3, additional mcmc samples will be taken if diffQerence is not great enough
                      if(m + mAdd >= maxMuse) {
                        mAdd <- maxMuse - m
                        hitMaxM <- TRUE
                      }
                      m <- m + mAdd
                      doMCMC(mAdd, thin = thin, reset = FALSE)
                    } else {
                      acceptThisIter <- TRUE
                      stronglyFlat <- diffQ + zGamma*sdDeltaQ <= tolUse # We reject H0:DeltaQ==tol (one-sided to H0:DeltaQ < tol) with (low) Type I error rate gamma. i.e. we believe  true diffQ < tol, so we did not move uphill more than tol, so we stayed flat. Using low Type I error rate means we have good evidence that true diffQ < tol.
                      # cat("stongly flat crit ", diffQ + zGamma*sdDeltaQ, " ", stronglyFlat, "\n")
                      if(stronglyFlat) numberOfTimesFlat <- numberOfTimesFlat + 1
                    }
                  } else {
                    acceptThisIter <- TRUE
                    weaklyFlat <- diffQ - zDelta*sdDeltaQ <= 0 # We don't reject H0:DeltaQ==0 (one-sided to H0:DeltaQ > 0) with (high) Type I error rate delta. i.e. we can't see clearly that true diffQ > 0, so we did not clearly move uphill, so we stayed weakly flat. Using high Type I error rate means we have a low threshold for believed we moved uphill, so we lower Type II error rate and are less likely to think we didn't if we really did.
                    # cat("weakly flat crit ", diffQ - zDelta*sdDeltaQ, " ", weaklyFlat, "\n")
                    if(weaklyFlat) numberOfTimesFlat <- numberOfTimesFlat + 1
                  }
                  if(acceptThisIter) {
                    if(verboseUse){
                      cat("  [Note] Iteration Number: ", itNum, ".\n", sep = "")
                      cat("  [Note] Current number of MCMC iterations: ", m, ".\n", sep = "")
                      cat("  [Note] Parameter Estimates: \n", sep = "")
                      params <- inverseTransform(paramsTrans)
                      print(params)
                      cat("  [Note] Convergence Criterion: ", diffQ+zGamma*sdDeltaQ, ".\n", sep = "")
                    }
                  }
                  if(itNum == maxIterUse)
                    cat("  [Note] Stopping MCEM because maximum number of MCEM iterations\n",
                        "(maxIter=", maxIterUse, ") was reached\n.")
                }
              }
      finalM <<- m
      if(!returnTrans) {
        if(!verboseUse) {
          # Note params might have been set above when verbose == TRUE
          params <- inverseTransform(paramsTrans)
        }
        optimOutput$par <- params
      }
      return(optimOutput)
      returnType(optimResultNimbleList())
    }
  ),
  run = function() {cat("  [Note] This run method does nothing. Try using findMLE.\n")}
)
