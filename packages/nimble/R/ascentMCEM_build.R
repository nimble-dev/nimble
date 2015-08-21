# asymptotic variance is calculated using the moving-block
# bootstrap method of "Markov Chain Monte Carlo in Statistical Mechanics"
# by Mignani & Rosa, 2001 (p. 350)
calc_asympVar = nimbleFunction(
  setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0, numReps){
    mvBlock <- modelValues(Rmodel)
    calc_E_llk <- calc_E_llk_gen(model, fixedNodes = fixedNodes, sampledNodes = sampledNodes, burnIn = burnIn, mvSample = mvBlock)
  },
  run = function(nsamps = integer(0), theta = double(1), oldTheta = double(1)){
    declare(svals, double(1, numReps))
    declare(h, integer(0)) #number of blocks to use for q function calculation
    declare(l, integer(0)) #length of each block
    declare(q, integer(0)) #total number of blocks available to sample from
    l <- ceiling(min(1000, nsamps/10))  #method for choosing l, ensures it's not too big
    q <- nsamps - l + 1
    h <- ceiling(nsamps/l)
    resize(mvBlock, h*l) #size our model value object to be approximately of size m (number of mc samples)
    for(r in 1:numReps){
      for(i in 1:h){
        randNum <- rbeta(1,1,1)
        randIndex <- ceiling(randNum*q) #random starting index for blocks
        for(j in 1:l){
          copy(mvSample, mvBlock, sampledNodes, sampledNodes, randIndex-1+j,  (i-1)*l+j) #fill in mvBlock with chosen blocks
        }
      }
      #as per Caffo, calculate both Q functions using the same samples from the latent variables
      oldQ <- calc_E_llk(oldTheta) 
      newQ <- calc_E_llk(theta)
      svals[r] <- newQ - oldQ
    }
    svalsSD <- sd(svals)
    svalsVar <- svalsSD^2
    returnType(double())
    return(svalsVar)
  },where = getLoadingNamespace()
)
    


calc_E_llk_gen = nimbleFunction(
    setup = function(model, fixedNodes, sampledNodes, mvSample, burnIn = 0){
	fixedCalcNodes <- model$getDependencies(fixedNodes)	
	latentCalcNodes <- model$getDependencies(sampledNodes)
	paramDepDetermNodes_fixed <- model$getDependencies(fixedNodes, determOnly = TRUE) 
	paramDepDetermNodes_latent <- model$getDependencies(latentCalcNodes, determOnly = TRUE) 
	areFixedDetermNodes <- length(paramDepDetermNodes_fixed) > 0
	areLatentDetermNodes <- length(paramDepDetermNodes_latent) >0
    },
    run = function(paramValues = double(1)){
        declare(sample_LL, double() )
        declare(mean_LL, double() )
        declare(nSamples, integer() )
        values(model, fixedNodes) <<- paramValues
        if(areFixedDetermNodes){
            simulate(model, paramDepDetermNodes_fixed)	#	Fills in the deterministic nodes
        }
        nSamples = getsize(mvSample)
        mean_LL <- 0
        for(i in (burnIn+1):nSamples){
            nimCopy(from = mvSample, to = model, nodes = sampledNodes, row = i)
            if(areLatentDetermNodes){
                simulate(model, paramDepDetermNodes_latent)	#	Fills in the deterministic nodes
            }
            sample_LL = calculate(model, latentCalcNodes)
            mean_LL = mean_LL + sample_LL
        }
        mean_LL <- mean_LL / nSamples
        if(is.nan(mean_LL)){
            mean_LL = -Inf	
        }
        returnType(double())
        return(mean_LL)
    },where = getLoadingNamespace())

                                      

#' Builds an MCEM algorithm from a given NIMBLE model
#' 
#' Takes a nimble model and builds an MCEM algorithm for it. The user must specify which latent nodes are to be integrated out in the E-Step.
#' All other stochastic non-data nodes will be maximized over. If the nodes do not have positive density on the entire real line, then box constraints can be used
#' to enforce this. 
#' The M-step is done by a nimble MCMC sampler. The E-step is done by a call to R's \code{optim} with \code{method = 'L-BFGS-B'}.
#' 
#' @param model 			A nimble model 
#' @param latentNodes 		A character vector of the names of the stochastic nodes to integrated out. Names can be expanded, but don't need to be. For example, if the model contains
#' \code{x[1], x[2] and x[3]} then one could provide either \code{latentNodes = c('x[1]', 'x[2]', 'x[3]')} or \code{latentNodes = 'x'}. 
#' @param burnIn			burn-in used for MCMC sampler in E step
#' @param mcmcControl		list passed to \code{MCMCSpec}, a nimble function that builds the MCMC sampler. See \code{help(MCMCSpec)} for more details
#' @param boxConstraints	A list of box constraints for the nodes that will be maximized over. Each constraint is a list in which the first element is a character vector of node names to which the constraint applies and the second element is a vector giving the lower and upper limits.  Limits of \code{-Inf} or \code{Inf} are allowed.
#' @param buffer			A buffer amount for extending the boxConstraints. Many functions with boundary constraints will produce \code{NaN} or -Inf when parameters are on the boundary.  This problem can be prevented by shrinking the boundary a small amount. 
#' @param alpha   probability of a type one error - here, the probability of accepting a parameter estimate that does not increase the likelihood
#' @param beta    probability of a type two error - here, the probability of rejecting a parameter estimate that does increase the likelihood
#' @param gamma   probability of deciding that the algorithm has converged, that is, that the difference between two Q functions is less than C, when in fact it has not.
#' @param C      determines when the algorithm has converged - when C falls above a (1-gamma) confidence interval around the difference in Q functions from time point t-1 to time point t, we say the algorithm has converged.
#' @param numReps number of bootstrap samples to use for asymptotic variance calculation


#'  @author Clifford Anderson-Bergman and Nick Michaud
#' @export
#' @details buildMCEM calls the NIMBLE compiler to create the MCMC and objective function as nimbleFunctions.  If the given model has already been used in compiling other nimbleFunctions, it is possible you will need to create a new copy of the model for buildMCEM to use.
#' Uses the ascent-based mcem algorithm of Caffo '05, which includes rules for automatically increasing the number of MC samples as iterations increase, and for determining when convergence has been reached,
#' @return
#' an R function that when called runs the MCEM algorithm. The function returned takes the arguments listed in Runtime Arguments.
#'
#'
#' See user manual for more 
#' @section Runtime Arguments:
#'	\itemize{
#'	\item{\code{initM}}	{
#'    starting number of iterations for the algorithm.
#'	}
#'
#' @examples
#' 
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
#' pumpMCEM <- buildAscent_MCEM(model = pumpModel, latentNodes = 'theta[1:10]',
#'                        boxConstraints = box)
#' pumpMCEM(initM = 1000)
#'
#' # Could also use latentNodes = 'theta' and buildMCEM would figure out this means 'theta[1:10]'
#' 
buildAscentMCEM <- function(model, latentNodes, burnIn = 100 , mcmcControl = list(adaptInterval = 20),
                      boxConstraints = list(), buffer = 10^-6, alpha = 0.05, beta = 0.05, gamma = 0.05, C = .01, numReps = 200, verbose = T) {
    latentNodes = model$expandNodeNames(latentNodes)
    latentNodes <- intersect(latentNodes, model$getNodeNames(stochOnly = TRUE))
    allStochNonDataNodes = model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    
    if(buffer == 0)
    	cat("warning: buffer 0. Can cause problems if the likelihood function is degenerate on boundary")
    if(buffer < 0)
    	stop('buffer must be non-negative')
    
    if(length(setdiff(latentNodes, allStochNonDataNodes) ) != 0 )
        stop('latentNodes provided not found in model')
    maxNodes = setdiff(allStochNonDataNodes, latentNodes)
    
    low_limits = rep(-Inf, length(maxNodes) ) 
    hi_limits = rep(Inf, length(maxNodes) ) 
    constraintNames = list()
    for(i in seq_along(boxConstraints) )
    	constraintNames[[i]] = model$expandNodeNames(boxConstraints[[i]][[1]])
    for(i in seq_along(constraintNames) ) {
        limits = boxConstraints[[i]][[2]]
        inds = which(maxNodes %in% constraintNames[[i]])
        if(length(inds) == 0)
            stop(  paste("warning: provided a constraint for node '", constraintNames[i], "' but that node does not exist in the model!") )
        low_limits[inds] = limits[1] + abs(buffer)
        hi_limits[inds] = limits[2] - abs(buffer)
    }
    
    if(any(low_limits>=hi_limits))
    	stop('lower limits greater than or equal to upper limits!')
    
    if(length(latentNodes) == 0)
        stop('no latentNodes')
    
    if(length(maxNodes) == 0)
        stop('no nodes to be maximized over')
    

    if(is(model, "RModelBaseClass") ){
    	Rmodel = model
        if(is(model$CobjectInterface, "uninitializedField")){
            cModel <- compileNimble(model)
        }
        else
            cModel = model$CobjectInterface
    }
    else{
        cModel <- model
        Rmodel <- model$Rmodel
    }
    
    zAlpha <- qnorm(alpha, 0, 1, lower.tail=FALSE)
    zBeta <- qnorm(beta, 0, 1, lower.tail=FALSE)
    zGamma <- qnorm(gamma, 0, 1, lower.tail=FALSE)
    


    mcmc_Latent_Spec <- configureMCMC(Rmodel, nodes = latentNodes, monitors = model$getVarNames(), control = mcmcControl) 
    Rmcmc_Latent <- buildMCMC(mcmc_Latent_Spec)
    sampledMV = Rmcmc_Latent$mvSamples
    Rcalc_E_llk <- calc_E_llk_gen(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)
    RvarCalc <- calc_asympVar(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV, numReps = numReps)
    
    cvarCalc <- compileNimble(RvarCalc, project = Rmodel)
    cmcmc_Latent = compileNimble(Rmcmc_Latent, project = Rmodel)
    cCalc_E_llk = compileNimble(Rcalc_E_llk, project = Rmodel)    
    nParams = length(maxNodes)
    run <- function(initM = 1000){
        theta = rep(NA, nParams)
        if(burnIn >= initM)
            stop('mcem quitting: burnIn > inital m value')
        cmcmc_Latent$run(1, reset = TRUE)	# To get valid initial values 
        theta <- values(cModel, maxNodes)
        
        for(i in seq_along(theta) ) {
            if(!(theta[i] >= low_limits[i] & theta[i] <= hi_limits[i]) )
                theta[i] = (low_limits[i] + hi_limits[i])/2				# This is necessary to insure that the initial values respect the constraints
                                        # Would only be a problem if the user supplies bounds that are more strict 
                                        # than necessary to imply a proper node value
        }
        
        m <- initM 
        endCrit <- C+1 #ensure that first iteration runs
        oldQ <- -1000 #ensure we accept first estimate
        sigSq <-0 #use initM as m value for first step
        diff <- 1 # any nonzero value can be used here, gets overwritten quickly in algo
        itNum <- 0
        while(endCrit > C){ 
          acceptCrit <- 0
          #starting sample size calculation for this iteration
          m <- ceiling(max(m, sigSq*((zAlpha + zBeta)^2)/((diff)^2)))
          cmcmc_Latent$run(m, reset = TRUE)   #initial mcmc run of size m
          thetaPrev <- theta  #store previous theta value
          itNum <- itNum + 1
          while(acceptCrit == 0){
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, control = list(fnscale = -1), method = 'L-BFGS-B', lower = low_limits, upper = hi_limits)
            theta = optimOutput$par    
            #varOut has two elements: varOut[1] is the sample variance, varOut[2] is the number of samples used to ccalculate the varianve
            sigSq <- cvarCalc$run(m, theta, thetaPrev) 
            ase <- sqrt(sigSq/numReps) #asymptotic std. error
            newQ <- cCalc_E_llk$run(theta)
            diff <- newQ-oldQ
            if((diff - zAlpha*ase)<0){ #swamped by mc error
              mAdd <- ceiling(m/2)  #from section 2.3, additional mcmc samples will be taken if difference is not great enough
              cmcmc_Latent$run(mAdd, reset = FALSE)
              m <- m + mAdd
            }
            else{
              acceptCrit <- 1
              endCrit <- diff + zGamma*ase #evaluate ending criterion
              cmcmc_Latent$run(m, reset = TRUE) #get Q(theta^(t-1),theta^(t-1)) as in equation 6
              oldQ <- cCalc_E_llk$run(theta) #save old q value 
              if(verbose == T){
                print(paste("Iteration Number:", itNum, sep = " "))
                print(paste("Current number of MCMC iterations:", m, sep = " "))
                output = optimOutput$par
                names(output) = maxNodes
                print("Parameter Estimates:")
                print(output)
                print(paste("Convergence Criterion:", endCrit, sep = " "))
              }
            }
          }
        }
        output = optimOutput$par
        names(output) = maxNodes
        return(output)
    }
    return(run)
}

