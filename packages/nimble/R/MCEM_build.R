# asymptotic variance is calculated using the moving-block
# bootstrap method of "Markov Chain Monte Carlo in Statistical Mechanics"
# by Mignani & Rosa, 2001 (p. 350)
calc_asympVar = nimbleFunction(
    name = 'calc_asympVar',
  setup = function(model, fixedNodes, sampledNodes, mvBlock, mvSample, burnIn = 0, numReps){
    calc_E_llk <- calc_E_llk_gen(model, fixedNodes = fixedNodes, sampledNodes = sampledNodes, burnIn = 0, mvSample = mvBlock)
  },
  run = function(nsamps = integer(0), theta = double(1), oldTheta = double(1)){
    ##declare(svals, double(1, numReps))
    svals <- numeric(numReps, init=FALSE)
    l <- ceiling(min(1000, (nsamps - burnIn)/20)) #length of each block, ensures it's not too big
    q <- (nsamps - burnIn) - l + 1 #total number of blocks available to sample from
    h <- ceiling((nsamps - burnIn)/l) #number of blocks to use for q function calculation
    resize(mvBlock, h*l) #size our model value object to be approximately of size m (number of mc samples)
    for(r in 1:numReps){
      for(i in 1:h){
        randNum <- rbeta(1,1,1)
        randIndex <- ceiling(randNum*q) #random starting index for blocks (post burn-in)
        for(j in 1:l){
          copy(mvSample, mvBlock, sampledNodes, sampledNodes, burnIn + randIndex-1+j,  (i-1)*l+j) #fill in mvBlock with chosen blocks
        }
      }
      #as per Caffo, calculate both Q functions using the same samples from the latent variables
      svals[r]  <- calc_E_llk$run(theta, oldTheta, 1) 
    }
    svalsSD <- sd(svals)
    svalsVar <- svalsSD^2
    returnType(double())
    return(svalsVar)
  },where = getLoadingNamespace()
)


# Calculates Q function if diff = 0, calculates difference in Q functions if diff = 1.
calc_E_llk_gen = nimbleFunction(
    name = 'calc_E_llk_gen',
  setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0){
    fixedCalcNodes <- model$getDependencies(fixedNodes)	
    latentCalcNodes <- model$getDependencies(sampledNodes)
    paramDepDetermNodes_fixed <- model$getDependencies(fixedNodes, determOnly = TRUE) 
    paramDepDetermNodes_latent <- model$getDependencies(latentCalcNodes, determOnly = TRUE) 
    areFixedDetermNodes <- length(paramDepDetermNodes_fixed) > 0
    areLatentDetermNodes <- length(paramDepDetermNodes_latent) >0
  },
  run = function(paramValues = double(1), oldParamValues = double(1), diff = integer(0)){
    nSamples = getsize(mvSample)
    mean_LL <- 0
    
    for(i in (burnIn+1):nSamples){
      nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
      values(model, fixedNodes) <<- paramValues  #first use new params, then old ones
      if(areFixedDetermNodes){
        simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
      }
      if(areLatentDetermNodes){
        simulate(model, paramDepDetermNodes_latent)	#	Fills in the deterministic nodes
      }
      sample_LL = calculate(model, latentCalcNodes)
      if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
            stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')' to see if node 'y[3]' is the cause of the problem). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
      mean_LL = mean_LL + sample_LL
      if(diff == 1){
        values(model, fixedNodes) <<- oldParamValues #now old params
        if(areFixedDetermNodes){
          simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
        }
        if(areLatentDetermNodes){
          simulate(model, paramDepDetermNodes_latent)  #	Fills in the deterministic nodes
        }
        sample_LL = calculate(model, latentCalcNodes)
        if(is.na(sample_LL) | is.nan(sample_LL) | sample_LL == -Inf | sample_LL == Inf)
            stop("Non-finite log-likelihood occurred; the MCEM optimization cannot continue. Please check the state of the compiled model (accessible as 'name_of_model$CobjectInterface') and determine which parameter values are causing the invalid log-likelihood by calling 'calculate' with subsets of the model parameters (e.g., 'name_of_model$CobjectInterface$calculate('y[3]')'). Note that if your model is maximizing over parameters whose bounds are not constant (i.e., depend on other parameters), this is one possible cause of such problems; in that case you might try running the MCEM without bounds, by setting 'forceNoConstraints = TRUE'.")
        mean_LL = mean_LL - sample_LL
      }
    }
    values(model, fixedNodes) <<- paramValues  #first use new params, then old ones
    if(areFixedDetermNodes){
      simulate(model, paramDepDetermNodes_fixed)  #	Fills in the deterministic nodes
    }
    mean_LL <- mean_LL / nSamples
    if(is.nan(mean_LL)){
      mean_LL = -Inf	
    }
    returnType(double())
    return(mean_LL)
  },where = getLoadingNamespace())


## helper function to extract ranges of nodes to be maximized
getMCEMRanges <- nimbleFunction(name = 'getMCEMRanges',
 setup = function(model, maxNodes, buffer){
    low_limits = rep(-Inf, length(maxNodes) ) 
    hi_limits  = rep(Inf,  length(maxNodes) )
    nodes <- model$expandNodeNames(maxNodes)
    for(i in seq_along(nodes)) {
      low_limits[i] = getBound(model, nodes[i], 'lower') + abs(buffer)
      hi_limits[i]  = getBound(model, nodes[i], 'upper')  - abs(buffer)
    }
    return(list(low_limits, hi_limits))
  }, where = getLoadingNamespace()
)

#' Builds an MCEM algorithm from a given NIMBLE model
#' 
#' Takes a NIMBLE model and builds an MCEM algorithm for it. The user must specify which latent nodes are to be integrated out in the E-Step.
#' All other stochastic non-data nodes will be maximized over. If the nodes do not have positive density on the entire real line, then box constraints can be used
#' to enforce this. 
#' The M-step is done by a nimble MCMC sampler. The E-step is done by a call to R's \code{optim} with \code{method = 'L-BFGS-B'} if the nodes are constrainted, or \code{method = 'BFGS'} if the nodes are unconstrained.
#' 
#' @param model a nimble model 
#' @param latentNodes character vector of the names of the stochastic nodes to integrated out. Names can be expanded, but don't need to be. For example, if the model contains
#' \code{x[1], x[2] and x[3]} then one could provide either \code{latentNodes = c('x[1]', 'x[2]', 'x[3]')} or \code{latentNodes = 'x'}.
#' @param burnIn burn-in used for MCMC sampler in E step
#' @param mcmcControl	list passed to \code{configureMCMC}, which builds the MCMC sampler. See \code{help(configureMCMC)} for more details
#' @param boxConstraints list of box constraints for the nodes that will be maximized over. Each constraint is a list in which the first element is a character vector of node names to which the constraint applies and the second element is a vector giving the lower and upper limits.  Limits of \code{-Inf} or \code{Inf} are allowed.  Any nodes that are not given constrains will have their constraints automatically determined by NIMBLE
#' @param buffer			A buffer amount for extending the boxConstraints. Many functions with boundary constraints will produce \code{NaN} or -Inf when parameters are on the boundary.  This problem can be prevented by shrinking the boundary a small amount. 
#' @param alpha   probability of a type one error - here, the probability of accepting a parameter estimate that does not increase the likelihood.  Default is 0.25. 
#' @param beta    probability of a type two error - here, the probability of rejecting a parameter estimate that does increase the likelihood.  Default is 0.25.
#' @param gamma   probability of deciding that the algorithm has converged, that is, that the difference between two Q functions is less than C, when in fact it has not.  Default is 0.05.
#' @param C      determines when the algorithm has converged - when C falls above a (1-gamma) confidence interval around the difference in Q functions from time point t-1 to time point t, we say the algorithm has converged. Default is 0.001.
#' @param numReps number of bootstrap samples to use for asymptotic variance calculation.
#' @param forceNoConstraints avoid any constraints even from parameter bounds implicit in the model structure (e.g., from dunif or dgamma distributions); setting this to TRUE might allow MCEM to run when the bounds of a parameter being maximized over depend on another parameter.
#' @param verbose logical indicating whether to print additional logging information.
#' 
#' @author Clifford Anderson-Bergman and Nicholas Michaud
#' @export
#' @details \code{buildMCEM} calls the NIMBLE compiler to create the MCMC and objective function as nimbleFunctions.  If the given model has already been used in compiling other nimbleFunctions, it is possible you will need to create a new copy of the model for buildMCEM to use.
#' Uses an ascent-based MCEM algorithm, which includes rules for automatically increasing the number of MC samples as iterations increase, and for determining when convergence has been reached.  Constraints for parameter values can be provided.  If contstraints are not provided, they will be automatically determined by NIMBLE.
#' @return
#' an R list with two elements:
#' \itemize{
#' \item{\code{run}} {A function that when called runs the MCEM algorithm. This function takes the arguments listed in \code{run} Arguments below.}
#' \item{\code{estimateCov}} {An EXPERIMENTAL function that when called estimates the asymptotic covariance of the parameters.  The covariance is estimated using the method of Louis (1982).
#' This function takes the arguments listed in \code{estimateCov} Arguments below.}
#' }
#'
#'
#' @section \code{run} Arguments:
#'	\itemize{
#'	\item{\code{initM}}	{
#'    starting number of iterations for the algorithm.
#'	}
#'  }
#' @section \code{estimateCov} Arguments:
#'  \itemize{
#'  \item{\code{MLEs}}{
#'    named vector of MLE values.  Must have a named MLE value for each stochastic, non-data, non-latent node.  If the \code{run()} method has alread been called,
#'    MLEs do not need to be provided.
#'  } 
#'  \item{\code{useExistingSamples}}{
#'    logical argument.  If \code{TRUE} and the \code{run()} method has previously been called, the covariance estimation will use MCMC samples from the last step of the MCEM algorithm.
#'    Otherwise, an MCMC algorithm will be run for 10,000 iterations, and those samples will be used.  Defaults to \code{FALSE}.
#'  }
#'  } 
#' @references 
#' 
#' Caffo, Brian S., Wolfgang Jank, and Galin L. Jones (2005). Ascent-based Monte Carlo expectation-maximization.  \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 67(2), 235-251.
#'  
#'  Louis, Thomas A  (1982). Finding the Observed Information Matrix When Using the EM Algorithm. \emph{Journal of the Royal Statistical Society. Series B (Statistical Methodology)}, 44(2), 226-233. 
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
#'                   data = pumpData, inits = pumpInits)
#'
#' # Want to maximize alpha and beta (both which must be positive) and integrate over theta
#' box = list( list(c('alpha','beta'), c(0, Inf)))
#'
#' pumpMCEM <- buildMCEM(model = pumpModel, latentNodes = 'theta[1:10]',
#'                        boxConstraints = box)
#' MLEs <- pumpMCEM$run(initM = 1000)
#' cov <- pumpMCEM$estimateCov()
#' }
#' # Could also use latentNodes = 'theta' and buildMCEM() would figure out this means 'theta[1:10]'
#' 
buildMCEM <- function(model, latentNodes, burnIn = 500 , mcmcControl = list(adaptInterval = 100),
                      boxConstraints = list(), buffer = 10^-6, alpha = 0.25, beta = 0.25, 
                      gamma = 0.05, C = 0.001, numReps = 300, forceNoConstraints = FALSE, verbose = TRUE) {
  latentNodes = model$expandNodeNames(latentNodes)
  latentNodes <- intersect(latentNodes, model$getNodeNames(stochOnly = TRUE))
  dataNodes <- model$getNodeNames(dataOnly = TRUE)
  allStochNonDataNodes = model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
  if(buffer == 0)
    warning("'buffer' is zero. This can cause problems if the likelihood function is degenerate on boundary")
  if(buffer < 0)
    stop("'buffer' must be non-negative.")
  
  if(length(setdiff(latentNodes, allStochNonDataNodes) ) != 0 )
    stop('latentNodes provided not found in model')
  maxNodes = model$expandNodeNames(setdiff(allStochNonDataNodes, latentNodes), returnScalarComponents = TRUE)
  if(any(model$isDiscrete(maxNodes))) stop(paste0("MCEM cannot optimize over discrete top-level parameters. The following top-level parameters in your model are discrete: ",
                                                  paste0(maxNodes[model$isDiscrete(maxNodes)], collapse = ', ')))
  
  limits <- getMCEMRanges(model, maxNodes, buffer)
  low_limits = limits[[1]]
  hi_limits  = limits[[2]]
  
  ## will be assign()'ed from within run() for use in getAsymptoticCov
  paramMLEs <- NA
  mcmcIters <- NA
  
  constraintNames = list()
  for(i in seq_along(boxConstraints) )
    constraintNames[[i]] = model$expandNodeNames(boxConstraints[[i]][[1]])
  for(i in seq_along(constraintNames) ) {
    limits = boxConstraints[[i]][[2]]
    inds = which(maxNodes %in% constraintNames[[i]])
    if(length(inds) == 0)
      stop(paste("warning: provided a constraint for nodes", constraintNames[[i]], ", but those nodes do not exist in the model!"))
    tooLowNodes <- which(limits[1] + abs(buffer) < low_limits[inds])
    tooHighNodes <- which(limits[2] - abs(buffer) > hi_limits[inds])
    if(length(tooLowNodes) > 0) warning(paste0("User-specified lower bound for ", constraintNames[[i]][tooLowNodes],
                                               " is below lower bound detected by NIMBLE.  "))
    if(length(tooHighNodes) > 0) warning(paste0("User-specified upper bound for ", constraintNames[[i]][tooHighNodes],
                                                " is above upper bound detected by NIMBLE.  "))
    low_limits[inds] = limits[1] + abs(buffer)
    hi_limits[inds] = limits[2] - abs(buffer)
  }
  if(any(low_limits>=hi_limits))
    stop('lower limits greater than or equal to upper limits!')
  if(forceNoConstraints ||
     identical(low_limits, rep(-Inf, length(low_limits))) &&
     identical(hi_limits, rep(Inf, length(hi_limits))))
    optimMethod = "BFGS"
  else 
    optimMethod = "L-BFGS-B"

  if(length(latentNodes) == 0)
    stop('no latentNodes')
  
  if(length(maxNodes) == 0)
    stop('no nodes to be maximized over')
  resetFunctions <- FALSE
  if(is(model, "RmodelBaseClass") ){
    Rmodel = model
    if(is(model$CobjectInterface, "uninitializedField")){
      cModel <- compileNimble(model)
    }
    else{
      cModel = model$CobjectInterface
      resetFunctions <- TRUE
    }
  }
  else{
    cModel <- model
    Rmodel <- model$Rmodel
    resetFunctions <- TRUE
  }
  
  zAlpha <- qnorm(alpha, 0, 1, lower.tail=FALSE)
  zBeta <- qnorm(beta, 0, 1, lower.tail=FALSE)
  zGamma <- qnorm(gamma, 0, 1, lower.tail=FALSE)
  
  mcmc_Latent_Conf <- configureMCMC(Rmodel, nodes = latentNodes, monitors = model$getVarNames(), control = mcmcControl, print = FALSE)
  Rmcmc_Latent <- buildMCMC(mcmc_Latent_Conf)
  sampledMV <- Rmcmc_Latent$mvSamples
  mvBlock <- modelValues(Rmodel)
  Rcalc_E_llk <- calc_E_llk_gen(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)
  RvarCalc <- calc_asympVar(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvBlock, mvSample = sampledMV, numReps = numReps)
  RgetCov <- bootstrapGetCov(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)
  
  cmcmc_Latent = compileNimble(Rmcmc_Latent, project = Rmodel, resetFunctions = resetFunctions)
  cGetCov = compileNimble(RgetCov, project = Rmodel)  
  cvarCalc <- compileNimble(RvarCalc, project = Rmodel)
  cCalc_E_llk = compileNimble(Rcalc_E_llk, project = Rmodel)  
  nParams = length(maxNodes)
  run <- function(initM = 1000){
    if(burnIn >= initM)
      stop('mcem quitting: burnIn > initial m value')
    cmcmc_Latent$run(1, reset = TRUE)	# To get valid initial values 
    theta <- values(cModel, maxNodes)
    if(optimMethod == "L-BFGS-B"){
      for(i in seq_along(maxNodes) ) {  # check that initial values satisfy constraints
        if(identical(low_limits[i], -Inf) && (hi_limits[i] < Inf)){
          if(theta[i] > hi_limits[i]){
            theta[i] <- hi_limits[i] - 1
          }
        }
        else if(identical(hi_limits[i], Inf) && (low_limits[i] > -Inf)){
          if(theta[i] < low_limits[i]){
            theta[i] <- low_limits[i] + 1
          }
        }
        else if((low_limits[i] > -Inf) && (hi_limits[i] < Inf)){
          if(!(theta[i] >= low_limits[i] & theta[i] <= hi_limits[i])){
            theta[i] = (low_limits[i] + hi_limits[i])/2
          }
        }
      }
      values(cModel, maxNodes) <<- theta
      simulate(cModel, cModel$getDependencies(maxNodes, self = FALSE))
    }
    m <- initM 
    endCrit <- C+1 #ensure that first iteration runs
    sigSq <-0 #use initM as m value for first step
    diff <- 1 # any nonzero value can be used here, gets overwritten quickly in algo
    itNum <- 0
    while(endCrit > C){ 
      acceptCrit <- 0
      #starting sample size calculation for this iteration
      m <- burnIn + ceiling(max(m - burnIn, sigSq*((zAlpha + zBeta)^2)/((diff)^2)))
      cmcmc_Latent$run(m, reset = TRUE)   #initial mcmc run of size m
      thetaPrev <- theta  #store previous theta value
      itNum <- itNum + 1
      while(acceptCrit == 0){
        if(optimMethod == "L-BFGS-B") 
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, oldParamValues = thetaPrev,
                                diff = 0, control = list(fnscale = -1), method = 'L-BFGS-B', 
                                lower = low_limits, upper = hi_limits)
        if(optimMethod == "BFGS")
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, oldParamValues = thetaPrev,
                                diff = 0, control = list(fnscale = -1), method = 'BFGS')
        
        theta = optimOutput$par    
        sigSq <- cvarCalc$run(m, theta, thetaPrev) 
        ase <- sqrt(sigSq) #asymptotic std. error
        diff <- cCalc_E_llk$run(theta, thetaPrev, 1)
        if((diff - zAlpha*ase)<0){ #swamped by mc error
          cat("Monte Carlo error too big: increasing MCMC sample size.\n")
          mAdd <- ceiling((m-burnIn)/2)  #from section 2.3, additional mcmc samples will be taken if difference is not great enough
          cmcmc_Latent$run(mAdd, reset = FALSE)
          m <- m + mAdd
        }
        else{
          acceptCrit <- 1
          endCrit <- diff + zGamma*ase #evaluate ending criterion
          if(itNum == 1)
            endCrit <- C+1 #ensure that at least two iterations are run
          
          if(verbose == T){
            cat("Iteration Number: ", itNum, ".\n", sep = "")
            cat("Current number of MCMC iterations: ", m, ".\n", sep = "")
            output = optimOutput$par
            names(output) = maxNodes
            cat("Parameter Estimates: \n", sep = "")
            print(output)
            cat("Convergence Criterion: ", endCrit, ".\n", sep = "")
          }
        }
      }
    }
    output <- optimOutput$par
    assign('paramMLEs', output, envir = parent.env(environment()))
    assign('mcmcIters', m, envir = parent.env(environment()))
    names(output) <- maxNodes
    return(output)
  }
  
  estimateCov = function(MLEs = NA, useExistingSamples = FALSE){
    delta <- .0001
    if(!(length(MLEs) == 1 && is.na(MLEs))){
      if(!is.numeric(MLEs)) stop("MLEs argument must be numeric.")
      if(!identical(sort(names(MLEs)), sort(maxNodes))){
        stop(paste('MLEs argument must be a named vector with MLEs for all of the following parameters: ', paste(maxNodes, collapse = ", ")))
      }
      covMLEs <- unname(MLEs[maxNodes])
    }
    else{
      if(length(paramMLEs) == 1 && is.na(paramMLEs)){
        stop(paste('No MLEs argument was provided, and the run() method has not been called yet.  Please call the run() method first or provide a named vector of MLEs.'))
      }
      else{
        covMLEs <- unname(paramMLEs)
      }
    }
    if(dim(as.matrix(cmcmc_Latent$mvSamples))[1]<2){
      if(useExistingSamples == TRUE){
        warning('MCMC over latent states has not been run yet, cannot have useExistingSamples = TRUE')
        useExistingSamples <- FALSE
      }
    }
    values(cModel, maxNodes) <- covMLEs
    calculate(cModel, cModel$getDependencies(maxNodes))
    if(!useExistingSamples){
      if(is.na(mcmcIters)) mcmcIters <- 20000
      cmcmc_Latent$run(mcmcIters)
    }
    FIM <- cGetCov$run(covMLEs, delta)
    cov <- solve(FIM)
    colnames(cov) <- maxNodes
    rownames(cov) <- maxNodes
    return(cov)
  }
  return(list(run = run,
              estimateCov = estimateCov))
}

bootstrapGetCov <- nimbleFunction(
    name = 'bootstrapGetCov',
    setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0){
      fixedCalcNodes <- model$getDependencies(fixedNodes)	
      latentCalcNodes <- model$getDependencies(sampledNodes)
      paramDepDetermNodes_fixed <- model$getDependencies(fixedNodes, determOnly = TRUE) 
      paramDepDetermNodes_latent <- model$getDependencies(latentCalcNodes, determOnly = TRUE) 
      areFixedDetermNodes <- length(paramDepDetermNodes_fixed) > 0
      areLatentDetermNodes <- length(paramDepDetermNodes_latent) > 0
    },
  run = function(theta = double(1), delta = double(0)){
    nSamples <- getsize(mvSample)
    paramLengths <- length(theta)
    fxph <- numeric(paramLengths)
    fxmh <- numeric(paramLengths)
    grad <- numeric(paramLengths)
    derivxy <- matrix(0, nrow = paramLengths, ncol = paramLengths)
    meanGrad <- numeric(paramLengths)
    meanGradGrad <-  matrix(0, nrow = paramLengths, ncol = paramLengths)
    meanDerivxy <- matrix(0, nrow = paramLengths, ncol = paramLengths)
    for(i in (burnIn+1):nSamples){
      values(model, fixedNodes) <<- theta
      if(areFixedDetermNodes){
        simulate(model, paramDepDetermNodes_fixed)  
      }
      if(areLatentDetermNodes){
        simulate(model, paramDepDetermNodes_latent)
      }
      nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
      if(areFixedDetermNodes){
        simulate(model, paramDepDetermNodes_fixed) 
      }
      if(areLatentDetermNodes){
        simulate(model, paramDepDetermNodes_latent)	
      }
      origValue <- calculate(model, latentCalcNodes)
      for(iNode in 1:paramLengths){
        theta[iNode] <- theta[iNode] + delta
        values(model, fixedNodes) <<- theta 
        if(areFixedDetermNodes){
          simulate(model, paramDepDetermNodes_fixed)  
        }
        if(areLatentDetermNodes){
          simulate(model, paramDepDetermNodes_latent)	
        }
        fxph[iNode] <- calculate(model, latentCalcNodes)
        theta[iNode] <- theta[iNode] - 2*delta
        values(model, fixedNodes) <<- theta 
        if(areFixedDetermNodes){
          simulate(model, paramDepDetermNodes_fixed)  
        }
        if(areLatentDetermNodes){
          simulate(model, paramDepDetermNodes_latent)
        }
        fxmh[iNode] <- calculate(model, latentCalcNodes)
        grad[iNode] <- (fxph[iNode] - fxmh[iNode])/(2*delta)
        theta[iNode] <- theta[iNode] + delta
        derivxy[iNode, iNode] <- (fxph[iNode] -2*origValue + fxmh[iNode])/(delta^2)
      }
      for(iNode in 1:paramLengths){
        if(iNode != paramLengths){
          for(jNode in (iNode+1):paramLengths){
            theta[iNode] <- theta[iNode] + delta
            theta[jNode] <- theta[jNode] + delta
            values(model, fixedNodes) <<- theta 
            if(areFixedDetermNodes){
              simulate(model, paramDepDetermNodes_fixed) 
            }
            if(areLatentDetermNodes){
              simulate(model, paramDepDetermNodes_latent)	
            }
            fxyph <- calculate(model, latentCalcNodes)
            theta[iNode] <- theta[iNode] - 2*delta
            theta[jNode] <- theta[jNode] - 2*delta
            values(model, fixedNodes) <<- theta
            if(areFixedDetermNodes){
              simulate(model, paramDepDetermNodes_fixed) 
            }
            if(areLatentDetermNodes){
              simulate(model, paramDepDetermNodes_latent)
            }
            fxymh <- calculate(model, latentCalcNodes)
            derivxy[iNode, jNode] <- (fxyph - fxph[iNode] - fxph[jNode] + 2*origValue - fxmh[iNode] - fxmh[jNode] + fxymh)/(2*delta^2)
            theta[iNode] <- theta[iNode] + delta
            theta[jNode] <- theta[jNode] + delta
          }
        }
      }
      meanGrad <- meanGrad + grad
      meanDerivxy <- meanDerivxy + derivxy
      meanGradGrad <- meanGradGrad + (grad%*%t(grad))
    }
    meanGrad <- meanGrad / (nSamples - burnIn)
    meanDerivxy <- meanDerivxy /  (nSamples - burnIn)
    meanGradGrad <- meanGradGrad /  (nSamples - burnIn)
    for(iNode in 1:paramLengths){
      if(iNode != paramLengths){
        for(jNode in (iNode+1):paramLengths){
          meanDerivxy[jNode, iNode] <- meanDerivxy[iNode, jNode] ## smarter way to do this with matrix mult?
        }
      }
    }
    returnType(double(2))
    returnMat <- -meanDerivxy - meanGradGrad +(meanGrad%*%t(meanGrad))
    return(returnMat)
  },where = getLoadingNamespace())
