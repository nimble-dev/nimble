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
## This sampler perform a Reversible Jump MCMC step for the node to which is assigned, using an univariate normal proposal distribution. This is a specialized sampler used by configureRJ function, when the model code is written without using indicator variables. 

#' @rdname samplers
#' @export
sampler_RJ <- nimbleFunction(
  name = 'sampler_RJ',
  contains = sampler_BASE, 
  setup = function(mvSaved, model, target, control) {
    ## target should be a coefficient to be set to a fixed value (usually zero) or not
    ## control should have
    ## 1 - mean of proposal jump distribution  (default 0)
    ## 2 - sd of proposal jump distribution  (default 1)
    ## 3 - prior prob of taking its fixedValue  (default 0.5)
    ## 4 - fixedValue (default 0)
    
    proposalScale <- control$scale
    proposalMean <- control$mean
    positive <- control$positive
    fixedValue <- control$fixedValue
    
    ## precompute ratio between prior probabilities
    logRatioProbFixedOverProbNotFixed <- log(control$priorProb) - log(1-control$priorProb)
    
    ## get all parameters related to the target
    calcNodes <- model$getDependencies(target)
  },
  
  run = function(){
    ## Reversible-jump move
    
    ## get current value of the parameter we are interested in
    currentValue <- model[[target]]
    
    ## get current posterior log likelihood
    currentLogProb <- model$getLogProb(calcNodes)
    
    if(currentValue != fixedValue)
    {
      ##----------------------##
      ## remove proposal
      ##----------------------##
      ## log probability of the reverse proposal
      if(positive){
        ## Taking the absolute value the proposal is a folded normal
        logProbReverseProposal <- log(dnorm(currentValue, mean = proposalMean, sd = proposalScale) + 
                                        dnorm(-currentValue, mean = proposalMean, sd = proposalScale))
      } else {
        logProbReverseProposal <- dnorm(currentValue, mean = proposalMean, sd = proposalScale, log = TRUE)
      }

      model[[target]] <<- fixedValue
      proposalLogProb <- model$calculate(calcNodes)
      log_accept_prob <- proposalLogProb - currentLogProb - logRatioProbFixedOverProbNotFixed + logProbReverseProposal
    } else {
      ##----------------------##
      ## add proposal
      ##----------------------##
      
      ## Only positive proposal
      if(positive){
        ## Taking the absolute value the proposal is a folded normal
        proposalValue <- abs(rnorm(1, mean = proposalMean, sd = proposalScale))
        logProbForwardProposal <- log(dnorm(proposalValue, mean = proposalMean, sd = proposalScale) + 
                                        dnorm(-proposalValue, mean = proposalMean, sd = proposalScale))
      } else {
        proposalValue <- rnorm(1, mean = proposalMean, sd = proposalScale)
        logProbForwardProposal <- dnorm(proposalValue, mean =  proposalMean, sd = proposalScale, log = TRUE)
      }
      
      model[[target]] <<- proposalValue
      
      proposalLogProb <- model$calculate(calcNodes)
      log_accept_prob <- proposalLogProb - currentLogProb + logRatioProbFixedOverProbNotFixed - logProbForwardProposal  
    }
    
    accept <- decide(log_accept_prob)
    
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
## This sampler perform a Reversible Jump MCMC step for the node to which is assigned, using an univariate normal proposal distribution. This is a specialized sampler used by configureRJ function, when the model code is written using indicator variables associated with the targetNodes of interest. 

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
    
    coefNode <- control$targetNode
    proposalScale <- control$scale
    proposalMean <- control$mean
    positive <- control$positive
    
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
      
      if(positive){
        ## Taking the absolute value the proposal is a folded normal
        logProbReverseProposal <- log(dnorm(currentCoef, mean = proposalMean, sd = proposalScale) + 
                                        dnorm(-currentCoef, mean = proposalMean, sd = proposalScale))
      } else {
        ## reverse jumping density
        logProbReverseProposal <- dnorm(currentCoef, mean = proposalMean, sd = proposalScale, log = TRUE)
      }        
      model[[target]] <<- 0
      model[[coefNode]] <<- 0 ## here there should be fixedValue
      model$calculate(calcNodes)
      ## avoid including prior for coef not in model
      log_accept_prob <- model$getLogProb(calcNodesReduced) - currentLogProb + logProbReverseProposal
    } else {
      ##----------------------##
      ## add proposal
      ##----------------------##
      currentLogProb <- model$getLogProb(calcNodesReduced)
      

      if(positive){
        ## Taking the absolute value the proposal is a folded normal
        proposalCoef <- abs(rnorm(1, mean = proposalMean, sd = proposalScale))
        logProbForwardProposal <- log(dnorm(proposalCoef, mean = proposalMean, sd = proposalScale) + 
                                        dnorm(-proposalCoef, mean = proposalMean, sd = proposalScale))
      } else {
        proposalCoef <- rnorm(1, mean = proposalMean, sd = proposalScale)
        logProbForwardProposal <- dnorm(proposalCoef, mean =  proposalMean, sd = proposalScale, log = TRUE)
      }
      model[[target]] <<- 1
      model[[coefNode]] <<- proposalCoef
      
      proposalLogProb <- model$calculate(calcNodes)
      logAcceptanceProb <- proposalLogProb - currentLogProb - logProbForwardProposal
    }
    accept <- decide(log_accept_prob)
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
    fixedValue <- control[["fixedValue"]]
    original_samplerConf <- control[["samplerConf"]]
    
    if(is.null(original_samplerConf))
      stop("Making wrapper_sampler: Must provide control$samplerConf")
    
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
## helper function to check node configuration (in case of multiple calls)
nodeConfigurationCheck <- function(currentConf, node){
  if(length(currentConf) == 0) {
    warning(paste0("There are no samplers for ", node,". Skipping it."))
  }
  if(length(currentConf) > 1)
    warning(paste0("There is more than one sampler for ", node,". Only the first will be toggled."))
}

###-------------------------------###
#' Configure Reversible Jump sampler
#'
#' Modifies the \code{MCMCconf} object of a specific model, to include a Reversible Jump MCMC sampler for variable selection using an univariate normal proposal distribution. User the mean and scale of the proposal, and whether or not the proposal is strictly positive. This function supports two different ways of writing the model specifications. 
#'
#' @param mcmcConf A \code{MCMCconf} object.
#'
#' @param targetNodes A character vector, specifying the nodes and/or variables interested in the variable selection. Nodes may be specified in their indexed form, \code{y[1, 3]}. Alternatively, nodes specified without indexing will be expanded fully, e.g., \code{x} will be expanded to \code{x[1]]}, \code{x[2]}, etc.  
#' @param indicatorNodes An optional character vector, specifying the indicator nodes and/or variables paired with \code{targetNodes}. Nodes may be specified in their indexed form, \code{y[1, 3]}. Alternatively, nodes specified without indexing will be expanded fully, e.g., \code{x} will be expanded to \code{x[1]]}, \code{x[2]}, etc. Nodes must be provided consistently with \code{targetNodes} vector. (see details?)
#' @param priorProb An optional numeric vector of prior probabilities for each node to be in the model. (see details?)
#' @param control An optional list of control arguments to the Reversible Jump function.
#' \itemize{
#' \item mean. The mean of the normal proposal distribution. (default = 0)
#' \item scale. The standard deviation of the normal proposal distribution. (default  = 1)
#' \item positive. A logical argument specifying whether the proposal is strictly positive. (default = FALSE)
#' \item fixedValue. An optional value taken from the variable when out of the model, used only when \code{priorProb} is provided (default = 0)
#' }
#'
#' @details 
#' 
#' This functions modifies the samplers in \code{MCMCconf} for each of the nodes provided in the \code{targetNodes} argument. To these elements two samplers are assigned: a specialized Reversible Jump sampler and a modified version of the build sampler that is used only when the target node is in the model. 
#' \code{configureRJ} can handle two different ways of writing a nimble model, using or not indicator variables (examples?)
#' \enumerate{
#'   \item \code{targetNodes} and \code{priorProb} probabilities; in this case also a \fixedValue can be provided
#'   \item \code{targetNodes} and \code{indicatorNodes}; no fixedValue, but this choice allow to have random inclusion probabilities
#' }
#' 
#' 
#' 
#' 
#' @seealso \code{\link{sampler_BASE}} \code{\link{sampler_RJ}} \code{\link{configureRJ}}
#' 
#' @examples
#' 
#' @export
#' 
#' @references
#' 
#' Peter J. Green. (1995). Reversible jump markov chain monte carlo computation and bayesian model determination. \emph{Biometrika}, 82(4), 711-732.
configureRJ <- function(mcmcConf, targetNodes, indicatorNodes = NULL, priorProb = NULL, control = list(mean = NULL, scale = NULL, positive = NULL, fixedValue = NULL)) {

  nNodes <- length(targetNodes)
  
  ## control list extraction
  fixedValue        <- if(!is.null(control$fixedValue))        control$fixedValue       else 0
  mean              <- if(!is.null(control$mean))              control$mean             else 0
  scale             <- if(!is.null(control$scale))             control$scale            else 1
  positive          <- if(!is.null(control$positive))          control$positive         else FALSE
  
  
  ## if one value is provided for one element of the list, this value is repeated if there are multiple nodes
  # fixedValue        <- if(length(fixedValue) == 1L & nNodes > 1L)       rep(fixedValue, nNodes)       else fixedValue
  # mean              <- if(length(mean) == 1L & nNodes > 1L)             rep(mean, nNodes)             else mean
  # scale             <- if(length(scale) == 1L & nNodes > 1L)            rep(scale, nNodes)            else scale
  # positive  <- if(length(positive) == 1L & nNodes > 1L) rep(positive, nNodes) else positive
  # 
  ## above set of checks paired with this check
  # ## Check that RJ control arguments match targetNodes length
  # if(any(lengths(list(fixedValue, mean, scale, positive)) != nNodes)){
  #   stop("Arguments in control list must be of length 1 or match targetNodes vector length")
  # }
  
  ## more defensive version
  if(length(fixedValue) != nNodes) {
    if(length(fixedValue) == 1) fixedValue <- rep(fixedValue, nNodes) else stop('inconsistent length of fixedValue argument and specified number of RJ targetNodes')
  }
  if(length(mean) != nNodes) {
    if(length(mean) == 1) mean <- rep(mean, nNodes) else stop('inconsistent length of mean argument and specified number of RJ targetNodes')
  }
  if(length(scale) != nNodes) {
    if(length(scale) == 1) scale <- rep(scale, nNodes) else stop('inconsistent length of scale argument and specified number of RJ targetNodes')
  }
  if(length(positive) != nNodes) {
    if(length(positive) == 1) positive <- rep(positive, nNodes) else stop('inconsistent length of positive argument and specified number of RJ targetNodes')
  }
  
  
  ## flag for indicators and prior
  indicatorFlag <- (!is.null(indicatorNodes))
  priorFlag     <- (!is.null(priorProb))
  
  ## Check: user must provide indicatorNodes OR 
  if(indicatorFlag == priorFlag) {
    stop("Provide indicatorNodes or priorProb vector")
  }
  
  ##---------------------------------------##
  ## No indicator
  ##---------------------------------------##
  if(priorFlag){
    

    ## If one value for prior is given, it is used for each variable
    if(length(priorProb) != nNodes) {
      if(length(priorProb) == 1) priorProb <- rep(priorProb, nNodes) else stop('Length of priorProb vector must match targetNodes length')
    }
    

    for(i in 1:nNodes) {
      
      
      ## Create node list
      nodeAsScalar <- mcmcConf$model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
      
      ## Create RJ control list for the node 
      nodeControl  = list(
        fixedValue = fixedValue[i], 
        priorProb  = priorProb[i], 
        mean       = mean[i], 
        scale      = scale[i], 
        positive   = positive[i]) 
      
      ## if the node is not a scalar iterate through each element
      
      for(j in 1:length(nodeAsScalar)){
        currentConf <- mcmcConf$getSamplers(nodeAsScalar[j])
        
        ## check on node configuration
        nodeConfigurationCheck(currentConf, nodeAsScalar[j])
        
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

    ## Check that indicatorNodes vector  match targetNodes lenght
    if(length(indicatorNodes) != nNodes){ 
      stop("Length of indicatorNodes vector must match targetNodes length")
    }
    
    for(i in 1:nNodes) {
      
      ## Create node and indicatorNodes list
      nodeAsScalar <- mcmcConf$model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
      indicatorsAsScalar <- mcmcConf$model$expandNodeNames(indicatorNodes[i], returnScalarComponents = TRUE)
      
      if(length(nodeAsScalar) != length(indicatorsAsScalar)) {
        stop(paste0("indicatorNodes node ", indicatorNodes[i] ," does not match ", targetNodes[i], " size."))
      }
      
      nodeControl  = list(
        fixedValue = fixedValue[i], 
        mean       = mean[i], 
        scale      = scale[i], 
        positive = positive[i]) 
      
      for(j in 1:length(nodeAsScalar)){
          
        ## add coefficients to control
        nodeControl$targetNode <- nodeAsScalar[j]
        
        
        currentConf <- mcmcConf$getSamplers(nodeAsScalar[j])
        
        ## check on node configuration
        nodeConfigurationCheck(currentConf, nodeAsScalar[j])
        
        ## Add reversible jump sampler for the indicatorNodes variable
        mcmcConf$removeSamplers(indicatorsAsScalar[j])
        mcmcConf$addSampler(type = sampler_RJ_indicator,
                            target = indicatorsAsScalar[j],
                            control = nodeControl)
        
        ## Add sampler for the coefficient variable (when is in the model)
        mcmcConf$removeSamplers(nodeAsScalar[j])
        
        mcmcConf$addSampler(type = sampler_toggled,
                            target = nodeAsScalar[j],
                            control = list(samplerConf = currentConf[[1]], fixedValue = 0)) ## fixedValue is actually not used by sampler_RJ_indicator but passed because of how toggled_sampler is written
      }  
        
    }
  }
  mcmcConf
}

