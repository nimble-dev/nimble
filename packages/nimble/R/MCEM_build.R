
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
    },
    where = getLoadingNamespace())

                                        #	This is an R function that builds an MCEM algorithm for a model. 
                                        #	The latentNodes will be sampled and the rest of the stochastic non data nodes will be maximized



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
#' @author Clifford Anderson-Bergman
#' @export
#' @details buildMCEM calls the NIMBLE compiler to create the MCMC and objective function as nimbleFunctions.  If the given model has already been used in compiling other nimbleFunctions, it is possible you will need to create a new copy of the model for buildMCEM to use.
#' @return
#' an R function that when called runs the MCEM algorithm. The function returned takes the arguments listed in Runtime Arguments.
#'
#'
#' See user manual for more 
#' @section Runtime Arguments:
#'	\itemize{
#'	\item{\code{maxit}}	{
#'	
#'	 Maximum iterations of the algorithm. Right now, the MCEM runs until \code{maxit}, rather than evaluating a more complex stopping criteria}
#'
#'	\item{\code{m1}} 	{
#'	
#'	 number of MCMC samples taken in each iteration before \code{maxit/2}}
#'
#'	\item{\code{m2}}	{
#'	
#'	number of MCMC samples taken in each iteration after \code{maxit/2} }
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
#' pumpMCEM <- buildMCEM(model = pumpModel, latentNodes = 'theta[1:10]',
#'                        boxConstraints = box)
#' pumpMCEM(maxit = 3, m1 = 100, m2 = 300)
#' # full run (above is for quick illustration only)
#' # pumpMCEM(maxit = 40, m1 = 1000, m2 = 5000)
#'
#' # Could also use latentNodes = 'theta' and buildMCEM would figure out this means 'theta[1:10]'
buildMCEM <- function(model, latentNodes, burnIn = 100 , mcmcControl = list(adaptInterval = 20), boxConstraints = list(), buffer = 10^-6) {
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
    

    if(is(model, "RmodelBaseClass") ){
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


    mcmc_Latent_Spec <- configureMCMC(Rmodel, nodes = latentNodes, monitors = model$getVarNames(), control = mcmcControl) 
    Rmcmc_Latent <- buildMCMC(mcmc_Latent_Spec)
    sampledMV = Rmcmc_Latent$mvSamples
    Rcalc_E_llk <- calc_E_llk_gen(model, fixedNodes = maxNodes, sampledNodes = latentNodes, burnIn = burnIn, mvSample = sampledMV)


    cmcmc_Latent = compileNimble(Rmcmc_Latent, project = Rmodel)
    cCalc_E_llk = compileNimble(Rcalc_E_llk, project = Rmodel)    
    nParams = length(maxNodes)
    run <- function(maxit = 50, m1 = 1000, m2 = 5000){
        theta = rep(NA, nParams)
        if(burnIn >= m1)
            stop('mcem quitting: burnIn > m1')
        if(burnIn >= m2)
            stop('mcem quitting: burnIn > m2')
        cmcmc_Latent$run(1, reset = TRUE)	# To get valid initial values 
        theta <- values(cModel, maxNodes)
        
        for(i in seq_along(theta) ) {
            if(!(theta[i] >= low_limits[i] & theta[i] <= hi_limits[i]) )
                theta[i] = (low_limits[i] + hi_limits[i])/2				# This is necessary to insure that the initial values respect the constraints
                                        # Would only be a problem in the user supplies bounds that are more strict 
                                        # than necessary to imply a proper node value
        }
        
        sampleSize = m1
        for(it in 1:maxit){
            if(it > maxit/2)
                sampleSize = m2
            cmcmc_Latent$run(sampleSize, reset = TRUE)
            optimOutput = optim(par = theta, fn = cCalc_E_llk$run, control = list(fnscale = -1), method = 'L-BFGS-B', lower = low_limits, upper = hi_limits)
            theta = optimOutput$par            
            comp_llk <- cCalc_E_llk$run(theta) 	#actually doing this to set values. wasteful, but optim is going to be several fold longer
        }
        output = optimOutput$par
        names(output) = maxNodes
        return(output)
    }
    return(run)
}

