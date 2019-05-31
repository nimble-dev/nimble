##################################################
## Reversible Jump sampler - covariate selection
##################################################
## functions summary
## sampler_RJ:           does the jump proposal for the variable of interest
## sampler_RJ_indicator: does the jump proposal for an indicator variables
## sampler_toggled:      reassign the default sampler to the variable when is in the model  
## configureRJ:         substitute manual removeSampler/addSampler for each variable for which one wants to perform selections

##################################################
## RJ sampler - no indicatorNodes variable ############
##################################################
## Reversible Jump sampler
##
## This sampler perform a reversible jump MCMC step for the node to which is assigned, using an univariate normal proposal distribution. This is a specialized sampler used by configureRJ function, when the model code is written without using indicator variables. 

#' @rdname samplers
#' @export
#' 
sampler_RJ <- nimbleFunction(
  name = 'sampler_RJ',
  contains = sampler_BASE, 
  setup = function(model, mvSaved, target, control) {
    ## target should be a coefficient to be set to a fixed value (usually zero) or not
    ## control should have
    ## 1 - mean of proposal jump distribution  (default 0)
    ## 2 - sd of proposal jump distribution  (default 1)
    ## 3 - prior prob of taking its fixedValue 
    ## 4 - fixedValue (default 0)
    
    proposalScale <- control$scale
    proposalMean  <- control$mean
    fixedValue    <- control$fixedValue
    
    ## precompute ratio between prior probabilities
    logRatioProbFixedOverProbNotFixed <- log(control$priorProb) - log(1-control$priorProb)
    
    ## get all parameters related to the target
    calcNodes <- model$getDependencies(target)
    ## target not in reduced model so its prior not calculated
    calcNodesReduced <- model$getDependencies(target, self = FALSE)
  },
  
  run = function(){
    ## Reversible-jump move
    ## get current value of the parameter we are interested in
    currentValue <- model[[target]]
    
    ## get current posterior log likelihood
    
    if(currentValue != fixedValue)
    {
      ##----------------------##
      ## remove proposal
      ##----------------------##
      currentLogProb <- model$getLogProb(calcNodes)
      
      logProbReverseProposal <- dnorm(currentValue, mean = proposalMean, sd = proposalScale, log = TRUE)

      model[[target]] <<- fixedValue
      proposalLogProb <- model$calculate(calcNodes)
      logAcceptanceProb <- proposalLogProb - currentLogProb - logRatioProbFixedOverProbNotFixed + logProbReverseProposal
    } else {
      ##----------------------##
      ## add proposal
      ##----------------------##
      currentLogProb <- model$getLogProb(calcNodesReduced)

      proposalValue <- rnorm(1, mean = proposalMean, sd = proposalScale)
      logProbForwardProposal <- dnorm(proposalValue, mean =  proposalMean, sd = proposalScale, log = TRUE)

      model[[target]] <<- proposalValue

      proposalLogProb <- model$calculate(calcNodes)
      logAcceptanceProb <- proposalLogProb - currentLogProb + logRatioProbFixedOverProbNotFixed - logProbForwardProposal  
    }
    
    accept <- decide(logAcceptanceProb)
    
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function() {
  })
)

####################################################################
### RJ sampler - indicatorNodes variable ################################
####################################################################
## Reversible Jump sampler
##
## This sampler perform a reversible jump MCMC step for the node to which is assigned, using an univariate normal proposal distribution. This is a specialized sampler used by configureRJ function, when the model code is written using indicator variables associated with the targetNodes of interest. 

#' @rdname samplers
#' @export
#' 
sampler_RJ_indicator <- nimbleFunction(
  name = 'sampler_RJ_indicator',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control ) {
    ## target should be the name of the indicator variables,
    ## control should have 
    ## 1 - element called coef for the name of the corresponding coefficient 
    ## 2 - mean of proposal jump distribution  (default 0)
    ## 3 - sd of proposal jump distribution  (default 1)
    
    coefNode      <- control$targetNode
    proposalScale <- control$scale
    proposalMean  <- control$mean
    
    ## model with coefficient included
    calcNodes <- model$getDependencies(c(coefNode, target))
    ## coefNode not in reduced model so its prior not calculated
    calcNodesReduced <- model$getDependencies(target)
  },
  run = function( ) {
    currentIndicator <- model[[target]]
    if(currentIndicator == 1) {
      ##----------------------##
      ## remove proposal
      ##----------------------##
      currentLogProb <- model$getLogProb(calcNodes)
      currentCoef <- model[[coefNode]]
      
      ## reverse jumping density
      logProbReverseProposal <- dnorm(currentCoef, mean = proposalMean, sd = proposalScale, log = TRUE)
      
      model[[target]] <<- 0
      model[[coefNode]] <<- 0 
      model$calculate(calcNodes)
      ## avoid including prior for coef not in model
      logAcceptanceProb <- model$getLogProb(calcNodesReduced) - currentLogProb + logProbReverseProposal
    } else {
      ##----------------------##
      ## add proposal
      ##----------------------##
      currentLogProb <- model$getLogProb(calcNodesReduced)
      
      ## proposal value for the coefficient
      proposalCoef <- rnorm(1, mean = proposalMean, sd = proposalScale)
      logProbForwardProposal <- dnorm(proposalCoef, mean =  proposalMean, sd = proposalScale, log = TRUE)

      model[[target]] <<- 1
      model[[coefNode]] <<- proposalCoef
      
      proposalLogProb <- model$calculate(calcNodes)
      logAcceptanceProb <- proposalLogProb - currentLogProb - logProbForwardProposal
    }
    accept <- decide(logAcceptanceProb)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function() {})
)

##-------------------------------##
## helper functions 
##-------------------------------##
# Sample according to default assigned sampler when the target is in the model

#' @rdname samplers
#' @export
sampler_toggled <- nimbleFunction(
  name = 'sampler_toggled',
  ## This nimbleFunction generalizes the role of RW_sampler_nonzero from the web-site example
  ## Sample according to build sampler when the target is in the model (here different from zero)
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {

    fixedValue        <- if(!is.null(control[["fixedValue"]]))   control[["fixedValue"]]  else 0
    # fixedValue <- control[["fixedValue"]]
    original_samplerConf <- control[["samplerConf"]]
    
    if(is.null(original_samplerConf))
      stop("sampler_toggled: 'control$samplerConf' must be provided.")
    
    ## List but only with one element 
    original_sampler <- nimbleFunctionList(sampler_BASE)
    ## This function gets the original sampler configuration 
    original_sampler[[1]] <- original_samplerConf$buildSampler(model=model, mvSaved=mvSaved)
  },
  run = function() {
    if(model[[target]] != fixedValue)
      original_sampler[[1]]$run()
  },
  methods = list(
    reset = function() {
      original_sampler[[1]]$reset()
    }
  )
)


###-------------------------------###
## configureRJ
###-------------------------------###
#' Configure Reversible Jump Sampler
#'
#' Modifies the \code{MCMCconf} object of a specific model to include a reversible jump MCMC sampler for variable selection using an univariate normal proposal distribution. Users can control the mean and scale of the proposal. This function supports two different types of model specification: with and without indicator variables. 
#'
#' @param mcmcConf An \code{MCMCconf} object.
#' @param targetNodes A character vector, specifying the nodes and/or variables for which variable selection is to be performed. Nodes may be specified in their indexed form, \code{'y[1, 3]'}. Alternatively, nodes specified without indexing will be expanded fully, e.g., \code{'x'} will be expanded to \code{'x[1]'}, \code{'x[2]'}, etc.  
#' @param indicatorNodes An optional character vector, specifying the indicator nodes and/or variables paired with \code{targetNodes}. Nodes may be specified in their indexed form, \code{'y[1, 3]'}. Alternatively, nodes specified without indexing will be expanded fully, e.g., \code{'x'} will be expanded to \code{'x[1]'}, \code{'x[2]'}, etc. Nodes must be provided consistently with \code{targetNodes} vector. See details.
#' @param priorProb An optional value or vector of prior probabilities for each node to be in the model. See details.
#' @param control An optional list of control arguments:
#' \itemize{
#' \item mean. The mean of the normal proposal distribution (default = 0).
#' \item scale. The standard deviation of the normal proposal distribution (default  = 1).
#' \item fixedValue. (Optional) Value for the variable when it is out of the model, which can be used only when \code{priorProb} is provided (default = 0). If specified when \code{indicatorNodes} is passed, a warning is given and \code{fixedValue} is ignored. 
#' }
#'
#' @details 
#' 
#' This function modifies the samplers in \code{MCMCconf} for each of the nodes provided in the \code{targetNodes} argument. To these elements two samplers are assigned: a specialized reversible jump sampler and a modified version of the existing sampler that is used only when the target node is already in the model. 
#' \code{configureRJ} can handle two different ways of writing a NIMBLE model, either using indicator variables or not as shown in the examples below and discussed further in the MCMC Chapter in the NIMBLE User Manual. When using indicator variables, the user should provide the \code{indicatorNodes} argument and no \code{fixedValue} argument should be used. When not using indicator variables, the user should provide the \code{priorProb} argument and no \code{indicatorNodes} argument should be provided. In the latter case, the user can provide a non-zero value for \code{fixedValue} if desired. 
#'
#' Note that this functionality is intended for variable selection in regression-style models but may be useful for other situations as well. At the moment, setting a variance component to zero and thereby removing a set of random effects that are explicitly part of a model will not work because MCMC sampling in that case would need to propose values for multiple parameters (the random effects), whereas this functionality only proposes to add a single parameter.
#' 
#' @seealso \code{\link{sampler}} \code{\link{configureRJ}}
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' ## Linear regression with intercept and two covariates, using indicator variables
#' 
#' code <- nimbleCode({
#'   beta0 ~ dnorm(0, sd = 100)
#'   beta1 ~ dnorm(0, sd = 100)
#'   beta2 ~ dnorm(0, sd = 100)
#'   sigma ~ dunif(0, 100) 
#' 
#'   z1 ~ dbern(psi)  ## indicator variable associated with beta1
#'   z2 ~ dbern(psi)  ## indicator variable associated with beta2
#'   psi ~ dunif(0, 1) ## hyperprior on inclusion probability
#'   for(i in 1:N) {
#'     Ypred[i] <- beta0 + beta1 * z1 * x1[i] + beta2 * z2 * x2[i]
#'     Y[i] ~ dnorm(Ypred[i], sd = sigma)
#'   }
#' })
#' 
#' ## simulate some data
#' set.seed(1)
#' N <- 100
#' x1 <- runif(N, -1, 1)
#' x2 <- runif(N, -1, 1) ## this covariate is not included
#' Y <- rnorm(N, 1 + 2.5 * x1, sd = 1)
#' 
#' ## build the model
#' rIndicatorModel <- nimbleModel(code, constants = list(N = N),
#'                                data = list(Y = Y, x1 = x1, x2 = x2), 
#'                                inits=  list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y),
#'                                z1 = 1, z2 = 1, psi = 0.5))
#' 
#' indicatorModelConf <- configureMCMC(rIndicatorModel)
#' 
#' ## Add reversible jump  
#' configureRJ(mcmcConf = indicatorModelConf,     ## model configuration
#'             targetNodes = c("beta1", "beta2"), ## coefficients for selection
#'             indicatorNodes = c("z1", "z2"),    ## indicators paired with coefficients
#'             control = list(mean = 0, scale = 2))
#' 
#' indicatorModelConf$addMonitors("beta1", "beta2", "z1", "z2")
#' 
#' rIndicatorMCMC <- buildMCMC(indicatorModelConf)
#' cIndicatorModel <- compileNimble(rIndicatorModel)
#' cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = rIndicatorModel)
#' 
#' set.seed(1)
#' samples <- runMCMC(cIndicatorMCMC, 10000, nburnin = 6000)
#' 
#' ## posterior probability to be included in the mode
#' mean(samples[ , "z1"])
#' mean(samples[ , "z2"])
#' 
#' ## posterior means when in the model
#' mean(samples[ , "beta1"][samples[ , "z1"] != 0])
#' mean(samples[ , "beta2"][samples[ , "z2"] != 0])
#'
#' 
#' ## Linear regression with intercept and two covariates, without indicator variables
#'
#' code <- nimbleCode({
#'   beta0 ~ dnorm(0, sd = 100)
#'   beta1 ~ dnorm(0, sd = 100)
#'   beta2 ~ dnorm(0, sd = 100)
#'   sigma ~ dunif(0, 100)
#'   for(i in 1:N) {
#'     Ypred[i] <- beta0 + beta1 * x1[i] + beta2 * x2[i]
#'     Y[i] ~ dnorm(Ypred[i], sd = sigma)
#'   }
#' })
#' 
#' rNoIndicatorModel <- nimbleModel(code, constants = list(N = N),
#'                                  data = list(Y = Y, x1 = x1, x2 = x2), 
#'                                  inits=  list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y)))
#' 
#' noIndicatorModelConf <- configureMCMC(rNoIndicatorModel)
#' 
#' ## Add reversible jump  
#' configureRJ(mcmcConf = noIndicatorModelConf,   ## model configuration
#'             targetNodes = c("beta1", "beta2"), ## coefficients for selection   
#'             priorProb = 0.5,                   ## prior probability of inclusion
#'             control = list(mean = 0, scale = 2))
#' 
#' ## add monitors
#' noIndicatorModelConf$addMonitors("beta1", "beta2")
#' rNoIndicatorMCMC <- buildMCMC(noIndicatorModelConf) 
#' 
#' cNoIndicatorModel <- compileNimble(rNoIndicatorModel)
#' cNoIndicatorMCMC <- compileNimble(rNoIndicatorMCMC, project = rNoIndicatorModel)
#' 
#' set.seed(1)
#' samples <- runMCMC(cNoIndicatorMCMC, 10000, nburnin = 6000 )
#' 
#' ## posterior probability to be included in the mode
#' mean(samples[ , "beta1"] != 0)
#' mean(samples[ , "beta2"] != 0)
#' 
#' ## posterior means when in the model
#' mean(samples[ , "beta1"][samples[ , "beta1"] != 0])
#' mean(samples[ , "beta2"][samples[ , "beta2"] != 0])
#' }
#' 
#'  
#' @export
#' 
#' @references
#' 
#' Peter J. Green. (1995). Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. \emph{Biometrika}, 82(4), 711-732.
#' 
configureRJ <- function(mcmcConf, targetNodes, indicatorNodes = NULL, priorProb = NULL, control = list(mean = NULL, scale = NULL, fixedValue = NULL)) {

  nNodes <- length(targetNodes)
  
  ## control list extraction
  fixedValue        <- if(!is.null(control$fixedValue))        control$fixedValue       else 0
  mean              <- if(!is.null(control$mean))              control$mean             else 0
  scale             <- if(!is.null(control$scale))             control$scale            else 1
    
  ## repeat values for multiple nodes if only one value is provided
  if(length(fixedValue) != nNodes) {
      if(length(fixedValue) == 1) fixedValue <- rep(fixedValue, nNodes) else
          stop("configureRJ: inconsistent length of 'fixedValue' argument and specified number of 'targetNodes'.")
  }
  if(length(mean) != nNodes) {
      if(length(mean) == 1) mean <- rep(mean, nNodes) else
          stop("configureRJ: inconsistent length of 'mean' argument and specified number of 'targetNodes'.")
  }
  if(length(scale) != nNodes) {
      if(length(scale) == 1) scale <- rep(scale, nNodes) else
          stop("configureRJ: inconsistent length of 'scale' argument and specified number of 'targetNodes'.")
  }
  
  ## flag for indicators and prior
  indicatorFlag <- (!is.null(indicatorNodes))
  priorFlag     <- (!is.null(priorProb))
  
  ## Check: user must provide indicatorNodes OR priorProb
  if(indicatorFlag == priorFlag) {
    stop("configureRJ: Provide 'indicatorNodes' or 'priorProb' vector")
  }
  
  ## Check: fixedValue can be used only with priorProb 
  if(indicatorFlag & any(fixedValue != 0)) {
   warning("configureRJ: 'fixedValue' can be provided only when using 'priorProb'; it will be ignored.")
  }

  ##---------------------------------------##
  ## No indicator
  ##---------------------------------------##
  if(priorFlag){
    
    ## check that priorProb values are in [0,1]
    if(any(priorProb < 0 | priorProb > 1)){
      stop("configureRJ: Elements in priorProb must be probabilities in [0,1].")
    }

    ## If one value for prior is given, it is used for each variable
    if(length(priorProb) != nNodes) {
      if(length(priorProb) == 1) priorProb <- rep(priorProb, nNodes) else stop("configureRJ: Length of 'priorProb' vector must match 'targetNodes' length.")
    }
    
    for(i in 1:nNodes) {     
      
      ## Create node list
      nodeAsScalar <- mcmcConf$model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
      
      ## if the node is multivariate check that samplers are univariate 
      if(length(nodeAsScalar) > 1) {
        if(length(mcmcConf$getSamplers(targetNodes[i])) == 1)
          stop(paste0("configureRJ: '", targetNodes[i], "' is multivariate and uses a joint sampler; only univariate samplers can be used with reversible jump sampling."))
      }

      ## Create RJ control list for the node 
      nodeControl  = list(
        priorProb  = priorProb[i], 
        mean       = mean[i], 
        scale      = scale[i], 
        fixedValue = fixedValue[i]) 
      
      ## if the node is not a scalar iterate through each element
      
      for(j in 1:length(nodeAsScalar)){
        currentConf <- mcmcConf$getSamplers(nodeAsScalar[j])
        
        ## check on node configuration
        if(length(currentConf) == 0) {
          warning(paste0("configureRJ: There are no samplers for '", nodeAsScalar[j],"'. Skipping it."))
          next
        }
       else if(any(sapply(currentConf, '[[', "name")  ==  "RJ" | 
          sapply(currentConf, '[[', "name")  ==  "RJ_indicator" |
          sapply(currentConf, '[[', "name")  ==  "toggled")) {
           ## check if there is already a RJ or RJ_indicator sampler assigned
          stop(paste0("configureRJ: Node '", nodeAsScalar[j],"' is already configured for reversible jump."))
        }
        else if(length(currentConf) > 1){
          warning(paste0("configureRJ: There is more than one sampler for '", nodeAsScalar[j], "'. Only the first will be used by toggled sampler, and others will be removed."))
        }

        ## substitute node sampler
        mcmcConf$removeSamplers(nodeAsScalar[j])
        mcmcConf$addSampler(type = sampler_RJ,
                            target = nodeAsScalar[j],
                            control = nodeControl)
        
        mcmcConf$addSampler(type = sampler_toggled,
                            target = nodeAsScalar[j],
                            control = list(samplerConf = currentConf[[1]], fixedValue = nodeControl$fixedValue))
      }  
        
    }
  } else {
    
    ##---------------------------##
    ## indicator
    ##---------------------------##

    ## Check that indicatorNodes vector match targetNodes length
    if(length(indicatorNodes) != nNodes){ 
      stop("configureRJ: Length of 'indicatorNodes' vector must match 'targetNodes' length.")
    }
    
    for(i in 1:nNodes) {
      
      ## Create node and indicatorNodes list
      nodeAsScalar <- mcmcConf$model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
      indicatorsAsScalar <- mcmcConf$model$expandNodeNames(indicatorNodes[i], returnScalarComponents = TRUE)
      
      ## if the node is multivariate check that samplers are univariate
      if(length(nodeAsScalar) > 1) {
        if(length(mcmcConf$getSamplers(targetNodes[i])) == 1)
          stop(paste0("configureRJ: '", targetNodes[i], "' is multivariate and uses a joint sampler; only univariate samplers can be used with reversible jump sampling."))
      }

      ## check that length of indicatorNodes matches targetNodes
      if(length(nodeAsScalar) != length(indicatorsAsScalar)) {
        stop(paste0("configureRJ: indicatorNodes node '", indicatorNodes[i] ,"' does not match '", targetNodes[i], "' size."))
      }
      

      nodeControl  = list( 
        mean       = mean[i], 
        scale      = scale[i])

      for(j in 1:length(nodeAsScalar)){
          
        ## add coefficients to control
        nodeControl$targetNode <- nodeAsScalar[j]
        
        currentConf <- mcmcConf$getSamplers(nodeAsScalar[j])
        
        ## check on node configuration
        if(length(currentConf) == 0) {
          warning(paste0("configureRJ: There are no samplers for '", nodeAsScalar[j],"'. Skipping it."))
        } else if(any(sapply(currentConf, '[[', "name")  ==  "RJ" | 
          sapply(currentConf, '[[', "name")  ==  "RJ_indicator" |
          sapply(currentConf, '[[', "name")  ==  "toggled")) {
          ## check if there is already a RJ or RJ_indicator sampler assigned
          stop(paste0("configureRJ: Node '", nodeAsScalar[j],"' is already configured for reversible jump."))
        } else if(length(currentConf) > 1){
          warning(paste0("configureRJ: There is more than one sampler for '", nodeAsScalar[j], "'. Only the first will be used by toggled sampler, and others will be removed."))
        }
        
        ## Add reversible jump sampler for the indicatorNodes variable
        mcmcConf$removeSamplers(indicatorsAsScalar[j])
        mcmcConf$addSampler(type = sampler_RJ_indicator,
                            target = indicatorsAsScalar[j],
                            control = nodeControl)
        
        ## Add sampler for the coefficient variable (when is in the model)
        mcmcConf$removeSamplers(nodeAsScalar[j])
        
        mcmcConf$addSampler(type = sampler_toggled,
                            target = nodeAsScalar[j],
                            control = list(samplerConf = currentConf[[1]])) 
      }  
        
    }
  }
  mcmcConf
}

