
## configureRJ:            modifies the MCMC configuration object, to make use of the following RJ samplers:
## sampler_RJ_fixed_prior: proposes addition/removal for variable of interest, using a specified prior probability
## sampler_RJ_indicator:   proposes transitions of a binary indicator variable, corresponding to a variable of interest
## sampler_RJ_toggled:     samples the variable of interest using the original sampling configuration, when variable is "in the model"



#' @rdname RJ_samplers
#' @export
#'
sampler_RJ_fixed_prior <- nimbleFunction(
    name = 'sampler_RJ_fixed_prior',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        proposalMean  <- control$mean
        proposalScale <- control$scale
        priorProb     <- control$priorProb
        fixedValue    <- control$fixedValue
        ## node list generation
        calcNodes <- model$getDependencies(target)
        calcNodesReduced <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        logRatioProbFixedOverProbNotFixed <- log(priorProb) - log(1-priorProb)
    },
    run = function() {
        currentValue <- model[[target]]
        if(currentValue == fixedValue) {  ## propose addition of target
            currentLogProb <- model$getLogProb(calcNodesReduced)
            proposalValue <- rnorm(1, proposalMean, proposalScale)
            logProbForwardProposal <- dnorm(proposalValue, proposalMean, proposalScale, log = TRUE)
            model[[target]] <<- proposalValue
            proposalLogProb <- model$calculate(calcNodes)
            logAcceptanceProb <- proposalLogProb - currentLogProb + logRatioProbFixedOverProbNotFixed - logProbForwardProposal
        } else {                          ## propose removal of target
            currentLogProb <- model$getLogProb(calcNodes)
            logProbReverseProposal <- dnorm(currentValue, proposalMean, proposalScale, log = TRUE)
            model[[target]] <<- fixedValue
            model$calculate(calcNodes)
            proposalLogProb <- model$getLogProb(calcNodesReduced)
            logAcceptanceProb <- proposalLogProb - currentLogProb - logRatioProbFixedOverProbNotFixed + logProbReverseProposal
        }
        accept <- decide(logAcceptanceProb)
        if(accept) { copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else     { copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE) }
    },
    methods = list(
        reset = function() { }
    )
)



#' @rdname samplers
#' @export
#' 
sampler_RJ_indicator <- nimbleFunction(
    name = 'sampler_RJ_indicator',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## note: target is the indicator variable,
        ## control$targetNode is the variable conditionally in the model
        ## control list extraction
        coefNode      <- control$targetNode
        proposalScale <- control$scale
        proposalMean  <- control$mean
        ## node list generation
        calcNodes <- model$getDependencies(c(coefNode, target))
        calcNodesReduced <- model$getDependencies(target)
    },
    run = function() {
        currentIndicator <- model[[target]]
        if(currentIndicator == 0) {   ## propose addition of coefNode
            currentLogProb <- model$getLogProb(calcNodesReduced)
            proposalCoef <- rnorm(1, proposalMean, proposalScale)
            logProbForwardProposal <- dnorm(proposalCoef, proposalMean, proposalScale, log = TRUE)
            model[[target]] <<- 1
            model[[coefNode]] <<- proposalCoef
            proposalLogProb <- model$calculate(calcNodes)
            logAcceptanceProb <- proposalLogProb - currentLogProb - logProbForwardProposal
        } else {                      ## propose removal of coefNode
            currentLogProb <- model$getLogProb(calcNodes)
            currentCoef <- model[[coefNode]]
            logProbReverseProposal <- dnorm(currentCoef, proposalMean, proposalScale, log = TRUE)
            model[[target]] <<- 0
            model[[coefNode]] <<- 0
            model$calculate(calcNodes)
            logAcceptanceProb <- model$getLogProb(calcNodesReduced) - currentLogProb + logProbReverseProposal
        }
        accept <- decide(logAcceptanceProb)
        if(accept) { copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else     { copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE) }
    },
    methods = list(
        reset = function() { }
    )
)



#' @rdname samplers
#' @export
sampler_RJ_toggled <- nimbleFunction(
    name = 'sampler_RJ_toggled',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        fixedValue  <- if(!is.null(control$fixedValue))  control$fixedValue  else 0
        samplerType <- if(!is.null(control$samplerType)) control$samplerType else stop('must provide \'samplerType\' control list element to RJ_toggled sampler')
        ## nested function and function list definitions
        samplerFL <- nimbleFunctionList(sampler_BASE)
        samplerFL[[1]] <- samplerType$buildSampler(model = model, mvSaved = mvSaved)
    },
    run = function() {
        if(model[[target]] != fixedValue)
            samplerFL[[1]]$run()
    },
    methods = list(
        reset = function() {
            samplerFL[[1]]$reset()
        }
    )
)



#' Configure Reversible Jump for Variable Selection
#'
#' Modifies an MCMC configuration object to perform a reversible jump MCMC sampling for variable selection, using a univariate normal proposal distribution.  Users can control the mean and scale of the proposal. This function supports two different types of model specification: with and without indicator variables.
#'
#' @param conf An \code{MCMCconf} object.
#' @param targetNodes A character vector, specifying the nodes and/or variables for which variable selection is to be performed. Nodes may be specified in their indexed form, \code{'y[1, 3]'}. Alternatively, nodes specified without indexing will be expanded, e.g., \code{'x'} will be expanded to \code{'x[1]'}, \code{'x[2]'}, etc.
#' @param indicatorNodes An optional character vector, specifying the indicator nodes and/or variables paired with \code{targetNodes}. Nodes may be specified in their indexed form, \code{'y[1, 3]'}. Alternatively, nodes specified without indexing will be expanded, e.g., \code{'x'} will be expanded to \code{'x[1]'}, \code{'x[2]'}, etc. Nodes must be provided consistently with \code{targetNodes}. See details.
#' @param priorProb An optional value or vector of prior probabilities for each node to be in the model. See details.
#' @param control An optional list of control arguments:
#' \itemize{
#' \item mean. The mean of the normal proposal distribution (default = 0).
#' \item scale. The standard deviation of the normal proposal distribution (default = 1).
#' \item fixedValue. Value for the variable when it is out of the model, which can be used only when \code{priorProb} is provided (default = 0). If specified when \code{indicatorNodes} is passed, a warning is given and \code{fixedValue} is ignored.
#' }
#'
#' @return \code{NULL} \code{configureRJ} modifies the input MCMC configuration object in place.
#'
#' @details
#'
#' This function modifies the samplers in MCMC configuration object for each of the nodes provided in the \code{targetNodes} argument. To these elements two samplers are assigned: a reversible jump sampler to transition the variable in/out of the model, and a modified version of the original sampler, which performs updates only when the target node is already in the model.
#'
#' \code{configureRJ} can handle two different ways of writing a NIMBLE model, either with or without indicator variables. When using indicator variables, the \code{indicatorNodes} argument must be provided. Without indicator variables, the \code{priorProb} argument must be provided. In the latter case, the user can provide a non-zero value for \code{fixedValue} if desired.
#'
#' Note that this functionality is intended for variable selection in regression-style models but may be useful for other situations as well. At the moment, setting a variance component to zero and thereby removing a set of random effects that are explicitly part of a model will not work because MCMC sampling in that case would need to propose values for multiple parameters (the random effects), whereas the current functionality only proposes adding/removing a single model node.
#' 
#' @seealso \code{\link{samplers}} \code{\link{configureMCMC}}
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
#'   z1 ~ dbern(psi)   ## indicator variable associated with beta1
#'   z2 ~ dbern(psi)   ## indicator variable associated with beta2
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
#'                                inits = list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = sd(Y),
#'                                z1 = 1, z2 = 1, psi = 0.5))
#' 
#' indicatorModelConf <- configureMCMC(rIndicatorModel)
#' 
#' ## Add reversible jump  
#' configureRJ(conf = indicatorModelConf,        ## model configuration
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
#' configureRJ(conf = noIndicatorModelConf,      ## model configuration
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
#' samples <- runMCMC(cNoIndicatorMCMC, 10000, nburnin = 6000)
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
#' @author Sally Paganin, Perry de Valpine, Daniel Turek
#'  
#' @export
#' 
#' @references
#' 
#' Peter J. Green. (1995). Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. \emph{Biometrika}, 82(4), 711-732.
#' 
configureRJ <- function(conf, targetNodes, indicatorNodes = NULL, priorProb = NULL, control = list(mean = NULL, scale = NULL, fixedValue = NULL)) {
    model <- conf$model
    nNodes <- length(targetNodes)
    fixedValue <- if(!is.null(control$fixedValue)) control$fixedValue else 0
    mean       <- if(!is.null(control$mean))       control$mean       else 0
    scale      <- if(!is.null(control$scale))      control$scale      else 1
    ## repeat values for multiple nodes if only one value is provided
    if(length(fixedValue) != nNodes)
        if(length(fixedValue) == 1) fixedValue <- rep(fixedValue, nNodes) else stop("configureRJ: inconsistent length of 'fixedValue' argument and specified number of 'targetNodes'.")
    if(length(mean) != nNodes)
        if(length(mean) == 1) mean <- rep(mean, nNodes) else stop("configureRJ: inconsistent length of 'mean' argument and specified number of 'targetNodes'.")
    if(length(scale) != nNodes)
        if(length(scale) == 1) scale <- rep(scale, nNodes) else stop("configureRJ: inconsistent length of 'scale' argument and specified number of 'targetNodes'.")
    ## flag for indicators and prior
    indicatorFlag <- !is.null(indicatorNodes)
    priorFlag     <- !is.null(priorProb)
    ## must provide either indicatorNodes or priorProb
    if(indicatorFlag == priorFlag) stop("configureRJ: Provide 'indicatorNodes' or 'priorProb' vector")
    ## fixedValue can be used only with priorProb
    if(indicatorFlag && any(fixedValue != 0)) warning("configureRJ: 'fixedValue' can be provided only when using 'priorProb'; it will be ignored.")
    ##
    if(priorFlag) {    ## no indicator variables; ues RJ_fixed_prior sampler
        ## check that priorProb values are in [0,1]
        if(any(priorProb < 0 | priorProb > 1)) stop("configureRJ: elements in priorProb must be probabilities in [0,1].")
        ## if one value for prior is given, it is used for each variable
        if(length(priorProb) != nNodes)
            if(length(priorProb) == 1) priorProb <- rep(priorProb, nNodes) else stop("configureRJ: Length of 'priorProb' vector must match 'targetNodes' length.")
        for(i in 1:nNodes) {
            nodeAsScalar <- model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
            ## if the node is multivariate check that samplers are univariate
            if(model$isMultivariate(targetNodes[i]) && length(conf$getSamplers(targetNodes[i])) == 1)
                stop(paste0("configureRJ: '", targetNodes[i], "' is multivariate and uses a joint sampler; only univariate samplers can be used with reversible jump sampling."))
            ## Create RJ control list for the node
            nodeControl <- list(priorProb = priorProb[i], mean = mean[i],
                                scale = scale[i], fixedValue = fixedValue[i])
            for(j in 1:length(nodeAsScalar)) {
                currentConf <- conf$getSamplers(nodeAsScalar[j])
                ## check on node configuration
                if(length(currentConf) == 0) {
                    warning(paste0("configureRJ: There are no samplers for '", nodeAsScalar[j],"'. Skipping it."))
                    next
                }
                else if(any(sapply(currentConf,'[[','name') == 'RJ_fixed_prior' |
                                sapply(currentConf,'[[','name') == 'RJ_indicator' |
                                    sapply(currentConf,'[[','name') == 'RJ_toggled')) {
                    stop(paste0("configureRJ: node '", nodeAsScalar[j],"' is already configured for reversible jump."))
                }
                else if(length(currentConf) > 1)
                    warning(paste0("configureRJ: There is more than one sampler for '", nodeAsScalar[j], "'. Only the first will be used by RJ_toggled sampler, and others will be removed."))
                ## substitute node sampler
                conf$removeSamplers(nodeAsScalar[j])
                conf$addSampler(type = sampler_RJ_fixed_prior,
                                    target = nodeAsScalar[j],
                                    control = nodeControl)
                conf$addSampler(type = sampler_RJ_toggled,
                                    target = nodeAsScalar[j],
                                    control = list(samplerType = currentConf[[1]], fixedValue = nodeControl$fixedValue))
            }
        }
    }
    ##
    if(indicatorFlag) {   ## indicator variables; ues RJ_indicator sampler
        if(length(indicatorNodes) != nNodes)
            stop("configureRJ: Length of 'indicatorNodes' vector must match 'targetNodes' length.")
        for(i in 1:nNodes) {
            nodeAsScalar <- model$expandNodeNames(targetNodes[i], returnScalarComponents = TRUE)
            indicatorsAsScalar <- model$expandNodeNames(indicatorNodes[i], returnScalarComponents = TRUE)
            ## if the node is multivariate check that samplers are univariate
            if(model$isMultivariate(targetNodes[i]) && length(conf$getSamplers(targetNodes[i])) == 1)
                stop(paste0("configureRJ: '", targetNodes[i], "' is multivariate and uses a joint sampler; only univariate samplers can be used with reversible jump sampling."))
            ## check that length of indicatorNodes matches targetNodes
            if(length(nodeAsScalar) != length(indicatorsAsScalar))
                stop(paste0("configureRJ: indicatorNodes node '", indicatorNodes[i] ,"' does not match '", targetNodes[i], "' size."))
            nodeControl  = list(mean = mean[i], scale = scale[i])
            for(j in 1:length(nodeAsScalar)) {
                nodeControl$targetNode <- nodeAsScalar[j]
                currentConf <- conf$getSamplers(nodeAsScalar[j])
                ## check on node configuration
                if(length(currentConf) == 0) {
                    warning(paste0("configureRJ: There are no samplers for '", nodeAsScalar[j],"'. Skipping it."))
                } else if(any(sapply(currentConf,'[[','name') == 'RJ_fixed_prior' |
                                  sapply(currentConf,'[[','name') == 'RJ_indicator' |
                                      sapply(currentConf,'[[','name') == 'RJ_toggled')) {
                    stop(paste0("configureRJ: Node '", nodeAsScalar[j],"' is already configured for reversible jump."))
                } else if(length(currentConf) > 1) {
                    warning(paste0("configureRJ: There is more than one sampler for '", nodeAsScalar[j], "'. Only the first will be used by RJ_toggled sampler, and others will be removed."))
                }
                ## Add reversible jump sampler for the indicatorNodes variable
                conf$removeSamplers(indicatorsAsScalar[j])
                conf$addSampler(type = sampler_RJ_indicator,
                                    target = indicatorsAsScalar[j],
                                    control = nodeControl)
                ## Add sampler for the coefficient variable (when is in the model)
                conf$removeSamplers(nodeAsScalar[j])
                conf$addSampler(type = sampler_RJ_toggled,
                                    target = nodeAsScalar[j],
                                    control = list(samplerType = currentConf[[1]]))
            }
        }
    }
    return(invisible(NULL))
}

