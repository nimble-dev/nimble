## samples from measure G after initial MCMC is run on a CRP-based model
## Used when syntax xi[1:N] ~ dCRP(conc) is used in BUGS.

##-----------------------------------------
##  Wrapper function for sampleDPmeasure
##-----------------------------------------

#' Get posterior samples for a Dirichlet process measure
#'
#' This function obtains posterior samples from a Dirichlet process distributed random measure of a model specified using the \code{dCRP} distribution.
#'
#' @param MCMC an MCMC class object, either compiled or uncompiled.
#' @param epsilon  used for determining the truncation level of the representation of the random measure.
#' 
#' @author Claudia Wehrhahn and Christopher Paciorek
#' 
#' @export
#' @details
#' This function provides samples from a random measure having a Dirichlet process prior. Realizations are almost surely discrete and represented by a (finite) stick-breaking representation (Sethuraman, 1994), whose atoms (or point masses) are independent and identically distributed. This sampler can only be used with models containing a \code{dCRP} distribution . 
#'
#' The \code{MCMC} argument is an object of class MCMC provided by \code{buildMCMC}, or its compiled version. The MCMC should already have been run, as \code{getSamplesDPmeasure} uses the parameter samples to  generates samples for the random measure. Note that the monitors associated with that MCMC must include the cluster membership variable (which has the \code{dCRP} distribution), the cluster parameter variables, all variables directly determining the \code{dCRP} concentration parameter, and any stochastic parent variables of the cluster parameter variables. See \code{help(configureMCMC)} or \code{help(addMonitors)} for information on specifying monitors for an MCMC.
#' 
#' The \code{epsilon} argument is used to determine the truncation level of the random measure. \code{epsilon} is the tail probability of the random measure, which together with posterior samples of the concentration parameter, determines the truncation level (see Section 3 in Gelfand, A.E. and Kottas, A. 2002). The default value is 1e-4.
#'  
#' The returned list contains a matrix with samples from the random measure (one sample per row) and the truncation level. The stick-breaking weights are named \code{weights} and the atoms, or point masses, are named based on the cluster variables in the model.
#' 
#' @seealso \code{\link{buildMCMC}}, \code{\link{configureMCMC}}, 
#' @references
#'
#' Sethuraman, J. (1994). A constructive definition of Dirichlet priors. \emph{Statistica Sinica}, 639-650.
#'
#' Gelfand, A.E. and Kottas, A. (2002). A computational approach for full nonparametric Bayesian inference under Dirichlet process mixture models. \emph{ournal of Computational and Graphical Statistics}, 11(2), 289-305.
#' @examples
#' \dontrun{
#'   conf <- configureMCMC(model)
#'   mcmc <- buildMCMC(conf)
#'   cmodel <- compileNimble(model)
#'   cmcmc <- compileNimble(mcmc, project = model)
#'   runMCMC(cmcmc, niter = 1000)
#'   outputG <- getSamplesDPmeasure(cmcmc)
#'   samples <- outputG$samples
#'   truncation <- output$trunc
#' }
getSamplesDPmeasure <- function(MCMC, epsilon = 1e-4) {
    if(exists('model',MCMC, inherits = FALSE)) compiled <- FALSE else compiled <- TRUE
    if(compiled) {
      if(!exists('Robject', MCMC, inherits = FALSE) || !exists('model', MCMC$Robject, inherits = FALSE))
        stop("getSamplesDPmeasure: problem with finding model object in compiled MCMC")
      model <- MCMC$Robject$model
      mvSamples <- MCMC$Robject$mvSamples
    } else {
      model <- MCMC$model
      mvSamples <- MCMC$mvSamples
    }

    ## Create and run the DPmeasure sampler.
    rsampler <- sampleDPmeasure(model, mvSamples, epsilon)
    if(compiled) {
      csampler <- compileNimble(rsampler, project = model)
      csampler$run()
      samplesMeasure <- csampler$samples
    } else {
      rsampler$run()
      samplesMeasure <- rsampler$samples
    }
    
    namesVars <- rsampler$tildeVars
    p <- length(namesVars)

    lengthData <- rsampler$lengthData
    dimTilde <- rsampler$dimTilde
    dimTildeNim <- rsampler$dimTildeNim
    truncG <- ncol(samplesMeasure) / (sum(dimTilde)+1) 
    namesW <- sapply(seq_len(truncG), function(i) paste0("weight[", i, "]"))
    
    namesAtoms <- NULL
    inames <- 1
    for(j in 1:p){
      if(dimTildeNim[j] == 0) { # scalar cluster parameter
        for(l in 1:truncG) {
          namesAtoms[inames] <- paste0(namesVars[j], "[", l, "]")
          inames <- inames + 1
        }
      }
      if(dimTildeNim[j] == 1) { # vector cluster parameter
        for(l in 1:truncG) {
          for(k in 1:lengthData){
            namesAtoms[inames] <- paste0(namesVars[j], "[", k, ",", l, "]")
            inames <- inames + 1
          }
        }
      }
      if(dimTildeNim[j] == 2) { # matrix cluster parameter
        for(l in 1:truncG){
          for(k in 1:lengthData) {
            for(k1 in 1:lengthData) {
              namesAtoms[inames] <- paste0(namesVars[j], "[", k1, ",", k, ",", l, "]")
              inames <- inames + 1
            }
          }
        }
      }
    }  
    
    colnames(samplesMeasure) <- c(namesW, namesAtoms)
    
    output <- list(samples = samplesMeasure, trunc = truncG)
    return(output)
}


sampleDPmeasure <- nimbleFunction(
  name = 'sampleDPmeasure',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, epsilon){
    ## Determine variables in the mv object and nodes/variables in the model.
    mvSavedVars <- mvSaved$varNames
    
    stochNodes <- model$getNodeNames(stochOnly = TRUE)
    distributions <- model$getDistribution(stochNodes) 
    
    ## Determine if there is a dCRP-distributed node and that it is monitored.
    dcrpIndex <- which(distributions == 'dCRP')
    if(length(dcrpIndex) == 1) {
      dcrpNode <- stochNodes[dcrpIndex] 
      dcrpVar <- model$getVarNames(nodes = dcrpNode)
    } else {
      if(length(dcrpIndex) == 0 ){
        stop('sampleDPmeasure: One node with a dCRP distribution is required.\n')
      }
      stop('sampleDPmeasure: Currently only models with one node with a dCRP distribution are allowed.\n')
    }
    if(sum(dcrpVar == mvSavedVars) == 0)
      stop('sampleDPmeasure: The node having the dCRP distribution has to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    
    ## Find the cluster variables, named tildeVars
    dcrpElements <- model$expandNodeNames(dcrpNode, returnScalarComponents = TRUE)
    clusterVarInfo <- findClusterNodes(model, dcrpVar) 
    tildeVars <- clusterVarInfo$clusterVars
    if( is.null(tildeVars) )  ## probably unnecessary as checked in CRP sampler, but best to be safe
      stop('sampleDPmeasure: The model should have at least one cluster variable.\n')

    ## Check that cluster parameters are IID, as required for random measure G
    isIID <- TRUE
    for(i in seq_along(clusterVarInfo$clusterNodes)) {
        valueExprs <- sapply(clusterVarInfo$clusterNodes[[i]], function(x) model$getValueExpr(x))
        names(valueExprs) <- NULL
        if(length(unique(valueExprs)) != 1) {
            isIID <- FALSE
        }
    }
    
    if(!isIID && length(tildeVars) == 2 && checkNormalInvGammaConjugacy(model, clusterVarInfo))
        isIID <- TRUE
      ## Tricky as MCMC might not be using conjugacy, but presumably ok to proceed regardless of how
      ## MCMC was done, since conjugacy existing would guarantee IID.
    if(!isIID) stop('sampleDPmeasure: cluster parameters have to be independent and identically distributed. \n')

    ## Check that necessary variables are being monitored.
      
    ## Check that cluster variables are monitored.
    counts <- tildeVars %in% mvSavedVars
    if( sum(counts) != length(tildeVars) ) 
      stop('sampleDPmeasure: The node(s) representing the cluster variables must be monitored in the MCMC (and therefore stored in the modelValues object).\n')  
    
    parentNodesTildeVars <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes %in% unlist(clusterVarInfo$clusterNodes)]
    for(i in seq_along(candidateParentNodes)) {
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE)
      for(j in seq_along(tildeVars)) {
        if(sum(aux == clusterVarInfo$clusterNodes[[j]][1]))
          parentNodesTildeVars <- c(parentNodesTildeVars, candidateParentNodes[i])
      }
    }
    if(length(parentNodesTildeVars)) {
      parentNodesTildeVarsDeps <- model$getDependencies(parentNodesTildeVars, self = FALSE)
    } else parentNodesTildeVarsDeps <- NULL
    ## make sure tilde nodes are included (e.g., if a tilde node has no stoch parents) so they get simulated
    parentNodesTildeVarsDeps <- model$topologicallySortNodes(c(parentNodesTildeVarsDeps, unlist(clusterVarInfo$clusterNodes)))
    
    if(!all(model$getVarNames(nodes = parentNodesTildeVars) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the cluster variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    if(is.null(parentNodesTildeVars)) parentNodesTildeVars <- tildeVars  ## to avoid NULL which causes compilation issues

    ## Check that parent nodes of cluster IDs are monitored.   
    parentNodesXi <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes == dcrpNode]
    for(i in seq_along(candidateParentNodes)) {
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE)
      if(sum(aux == dcrpNode)) {
        parentNodesXi <- c(parentNodesXi, candidateParentNodes[i])
      }
    }
    
    if(!all(model$getVarNames(nodes = parentNodesXi) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the membership variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    if(is.null(parentNodesXi)) parentNodesXi <- dcrpNode  ## to avoid NULL which causes compilation issues

    ## End of checks of monitors.
      
    fixedConc <- TRUE # assume that conc parameter is fixed. This will change in the if statement if necessary
    if(length(parentNodesXi)) {
      fixedConc <- FALSE
      parentNodesXiDeps <- model$getDependencies(parentNodesXi, self = FALSE)
      parentNodesXiDeps <- parentNodesXiDeps[!parentNodesXiDeps == dcrpNode]
    } else {
      parentNodesXiDeps <- dcrpNode
    }
    
    dataNodes <- model$getDependencies(dcrpNode, stochOnly = TRUE, self = FALSE)
    N <- length(model$getDependencies(dcrpNode, stochOnly = TRUE, self = FALSE))
    p <- length(tildeVars)
    lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE))
    dimTildeNim <- numeric(p+1) # nimble dimension (0 is scalar, 1 is 2D array, 2 is 3D array)
    dimTilde <- numeric(p+1) # dimension to be used in run code
    for(i in 1:p) {
      dimTildeNim[i] <- model$getDimension(clusterVarInfo$clusterNodes[[i]][1])
      dimTilde[i] <- lengthData^(dimTildeNim[i]) 
    }
    nTilde <- numeric(p+1)
    nTilde[1:p] <- clusterVarInfo$nTilde
    if(any(nTilde[1:p] != nTilde[1])){
      stop('sampleDPmeasure: All cluster parameters must have the same number of parameters.\n')
    }
    
    ## Storage object to be sized in run code based on MCMC output.
    samples <- matrix(0, nrow = 1, ncol = 1)   
    ## Truncation level of the random measure 
    truncG <- 0 
    
    ## control list extraction
    ## The error of approximation G is given by (conc / (conc +1))^{truncG-1}. 
    ## we are going to define an error of approximation and based on the posterior values of the conc parameter define the truncation level of G
    ## the error is between errors that are considered very very small in the folowing papers
    ## Ishwaran, H., & James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. Journal of the American Statistical Association, 96(453), 161-173.
    ## Ishwaran, H., & Zarepour, M. (2000). Markov chain Monte Carlo in approximate Dirichlet and beta two-parameter process hierarchical models. Biometrika, 87(2), 371-390.
    # epsilon <- 1e-4

    setupOutputs(lengthData)
  },
  
  run=function(){
    
    niter <- getsize(mvSaved) # number of iterations in the MCMC
    
    # defining the truncation level of the random measure's representation:
    if( fixedConc ) {
      dcrpAux <- model$getParam(dcrpNode, 'conc') + N
      concSamples <- nimNumeric(length = niter, value = dcrpAux)
    } else {
      concSamples <- numeric(niter)
      for( iiter in 1:niter ) {
        nimCopy(from = mvSaved, to = model, nodes = parentNodesXi, row=iiter) 
        model$calculate(parentNodesXiDeps)
        concSamples[iiter] <- model$getParam(dcrpNode, 'conc')
      }
      dcrpAux <- mean(concSamples) + N
    }
    
    truncG <<- log(epsilon) / log(dcrpAux / (dcrpAux+1)) + 1
    truncG <<- round(truncG)
    
    ## Storage object: matrix with nrow = number of MCMC iterations, and ncol = (1 + p)*truncG, where
    ## truncG the truncation level of the random measure G (an integer given by the values of conc parameter)
    ## (p+1) denoted the number of parameters, each of length truncG, to be stored: p is the number of cluster components (length(tildeVars)) and 1 is for the weights.
    samples <<- matrix(0, nrow = niter, ncol = truncG*(sum(dimTilde)+1)) 
    
    ## computing G(.) = sum_{l=1}^{truncG} w_l delta_{atom_l} (.):
    for(iiter in 1:niter){
      checkInterrupt()
      
      ## sampling weights:
      vaux <- rbeta(1, 1, concSamples[iiter] + N)
      v1prod <- 1
      samples[iiter, 1] <<- vaux  
      for(l1 in 2:truncG) {
        v1prod <- v1prod * (1-vaux)
        vaux <- rbeta(1, 1, concSamples[iiter]+N)
        samples[iiter, l1] <<- vaux * v1prod 
      }
      samples[iiter, 1:truncG] <<- samples[iiter, 1:truncG] / (1 - v1prod * (1-vaux)) # normalizing
      
      ## sampling atoms:
      ## getting the sampled unique values (tilde variables) and their probabilities of being sampled,
      ## need for computing density later.
      probs <- nimNumeric(N)
      uniqueValues <- matrix(0, nrow = N, ncol = sum(dimTilde))  
      xiiter <- mvSaved[dcrpVar, iiter]
      range <- min(xiiter):max(xiiter) 
      index <- 1
      for(i in seq_along(range)){   
        cond <- sum(xiiter == range[i])
        if(cond > 0){
          probs[index] <- cond
          ## slight workaround because can't compile mvSaved[tildeVars[j], iiter]  
          nimCopy(mvSaved, model, tildeVars, row = iiter)
          jcol <- 1
          for(j in 1:p){
            for(l in 1:dimTilde[j]) {
              if(dimTildeNim[j] == 0) { # scalars or matrices
                uniqueValues[index, jcol] <- values(model, tildeVars[j])[dimTilde[j]*(range[i] - 1) + l]     
              }
              if(dimTildeNim[j] == 2) { #  matrices
                uniqueValues[index, jcol] <- values(model, tildeVars[j])[dimTilde[j]*(range[i] - 1) + l]    
              }
              if(dimTildeNim[j] == 1) { # vectors
                uniqueValues[index, jcol] <- values(model, tildeVars[j])[range[i] + (l-1)*nTilde[j]]     
              }
              jcol <- jcol + 1
            }
          }
          index <- index+1
        }
      }
      probs[index] <- concSamples[iiter] 
      newValueIndex <- index 
      
      ## copy tilde parents into model for use in simulation below when simulate atoms of G_0  
      nimCopy(mvSaved, model, parentNodesTildeVars, row = iiter)

      for(l1 in 1:truncG) {
        index <- rcat(prob = probs[1:newValueIndex])
        if(index == newValueIndex){   # sample from G_0
          model$simulate(parentNodesTildeVarsDeps)
          sumCol <- truncG
          for(j in 1:p){
            for(l in 1:dimTilde[j]) {
              if(dimTildeNim[j] == 0) {
                samples[iiter, sumCol + dimTilde[j]*(l1 - 1) +l] <<- values(model, tildeVars[j])[l]     
              }
              if(dimTildeNim[j] == 2) {
                samples[iiter, sumCol + dimTilde[j]*(l1 - 1) +l] <<- values(model, tildeVars[j])[l]     
              }
              if(dimTildeNim[j] == 1) {
                samples[iiter, sumCol + dimTilde[j]*(l1 - 1) +l] <<- values(model, tildeVars[j])[(l-1)*nTilde[j] + 1]  
              }
              
            }
            sumCol <- sumCol + dimTilde[j]*truncG 
          }
        } else {   # sample one of the existing values
          sumCol <- truncG
          jcol <- 1
          for(j in 1:p){
            for(l in 1:dimTilde[j]) {
              samples[iiter, sumCol + dimTilde[j]*(l1 - 1) +l] <<- uniqueValues[index, jcol]    
              jcol <- jcol + 1
            }
            sumCol <- sumCol + dimTilde[j]*truncG 
          }
        }
      }
      
    }
  },
  methods = list( reset = function () {} )
)

## Sampler for concentration parameter, conc, of the dCRP distribution.

#' @rdname samplers
#' @export
sampler_CRP_concentration <- nimbleFunction(
  name = 'sampler_CRP_concentration',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    ## only dependency should be membership vector because configureMCMC checks for only one dependency  
    xiNode <- model$getDependencies(target, self = FALSE) 
    n <- length(model[[xiNode]])
  },
  
  
  run = function() {
    shapeParam <- model$getParam(target, 'shape')
    rateParam <- model$getParam(target, 'rate')
    conc <- model[[target]] 
    xi <- model[[xiNode]]
    
    occupied <- numeric(n)
    for( i in 1:n )
      occupied[xi[i]] <- 1
    k <- sum(occupied)
    
    aux1 <- shapeParam + k
    aux2 <- aux1 - 1
    
    ## generating augmented r.v. and computing the weight.
    x <- rbeta(1, conc+1, n)
    aux3 <- rateParam - log(x)
    w <- aux2/(aux2 + n*aux3)
    
    ## updating the concentration parameter.
    if(runif(1) <= w){
      model[[target]] <<- rgamma(1, aux1, aux3)
    }else{
      model[[target]] <<- rgamma(1, aux2, aux3)
    }
    
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {})
)




############################################################################
## Helper functions for sampling new clusters in sampler_CRP
############################################################################

## We need a base class because it is possible (perhaps unlikely) that
## a user model might have two uses of dCRP samplers that are different versions
## e.g., a nonconjugate and a dnorm_dnorm conjugate.
## To allow for this we need to use a one-element nimbleFunctionList that uses
## this virtual base class.
CRP_helper <- nimbleFunctionVirtual(
  methods = list(
    storeParams = function() { },   ## parameters of base measure, which for conjugate cases should be the same for all clusters
    calculate_prior_predictive = function(i = integer()) {  ## calculate prior predictive for new cluster for conjugate cases
      returnType(double())
    },
    sample = function(i = integer(), j = integer()) {}  ## sample parameter values for new cluster
  )
)

## In these helper functions, 'i' indexes the observation that is associated with xi[k]
## indexing into 'dataNodes.
## 'j' is the potential value of xi[k] under the new cluster and
## indexes into the vector of cluster parameter nodes, 'marginalizedNodes'.

CRP_nonconjugate <- nimbleFunction(
  name = "CRP_nonconjugate",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
      mv <- modelValues(model, m = 1)
  },
  methods = list(
    storeParams = function() {},  ## nothing needed for non-conjugate
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      return(model$getLogProb(dataNodes[i]))
    },
    sample = function(i = integer(), j = integer() ) {
        if(j == 0) {  ## reset to stored values (for case of new cluster not opened)
            nimCopy(mv, model, nodes = marginalizedNodes, row = 1)
        } else {
           ## sample from prior
            nimCopy(model, mv, nodes = marginalizedNodes, row = 1)
            if( p == 1 ) {
                model$simulate(marginalizedNodes[j])
            } else {
                for(l in 1:p) {  ## marginalized nodes should be in correct order based on findClusterNodes.
                    model$simulate(marginalizedNodes[(l-1)*nTilde + j])
                }
            }
        }
    }
  )
)

CRP_conjugate_dnorm_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dnorm_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorMean <- nimNumeric(1)
    priorVar <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorMean <<- model$getParam(marginalizedNodes[1], 'mean')
      priorVar <<- model$getParam(marginalizedNodes[1], 'var')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      dataVar <- model$getParam(dataNodes[i], 'var')
      y <- values(model, dataNodes[i])[1]
      return(dnorm(y, priorMean, sqrt(priorVar + dataVar), log=TRUE))
    },
    sample = function(i = integer(), j = integer()) {
      dataVar <- model$getParam(dataNodes[i], 'var')
      y <- values(model, dataNodes[i])[1]
      postVar <- 1 / (1 / dataVar + 1 / priorVar)
      postMean <- postVar * (y / dataVar + priorMean / priorVar)
      values(model, marginalizedNodes[j]) <<- c(rnorm(1, postMean, sqrt(postVar)))
    }
  )
)


CRP_conjugate_dnorm_invgamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dnorm_invgamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes1, marginalizedNodes2, dataNodes) {
    priorMean <- nimNumeric(1)
    kappa <- nimNumeric(1)
    priorShape <- nimNumeric(1)
    priorScale <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorMean <<- model$getParam(marginalizedNodes1[1], 'mean')
      kappa <<- values(model, marginalizedNodes2[1])[1]/model$getParam(marginalizedNodes1[1], 'var') # construct kappa as sigma2/(sigma2/kappa)
      priorShape <<- model$getParam(marginalizedNodes2[1], 'shape')
      priorScale <<- model$getParam(marginalizedNodes2[1], 'scale')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      c1 <- priorShape * log(priorScale) + lgamma(priorShape + 1/2) + 0.5*log(kappa) -
       lgamma(priorShape) - 0.5*log(2) - 0.5*log(pi) - 0.5*log(1 + kappa)
      c2 <- - (priorShape  + 1/2) * log( (priorScale + kappa * (y - priorMean)^2 / (2*(1+kappa)) ) )
      return(c1 + c2)
    },
    sample = function(i = integer(), j = integer()) {
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes2[j]) <<- c(rinvgamma(1, shape = priorShape + 1/2,
                                                  scale = priorScale + kappa * (y - priorMean)^2 / (2*(1+kappa)) ))
      values(model, marginalizedNodes1[j]) <<- c(rnorm(1, mean = (kappa * priorMean + y)/(1 + kappa), 
                                              sd = sqrt(values(model, marginalizedNodes2[j])[1] / (1+kappa))) )
    }
  )
)


CRP_conjugate_dgamma_dpois <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dpois",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape <- nimNumeric(1)
    priorRate <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape <<- model$getParam(marginalizedNodes[1], 'shape') 
      priorRate <<- model$getParam(marginalizedNodes[1], 'rate') 
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      return(priorShape * log(priorRate) - (priorShape + y) * log(priorRate + 1) +
               lgamma(priorShape + y) - lgamma(priorShape) - lfactorial(y))
    },
    sample = function(i = integer(), j = integer()) {
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rgamma(1, shape = priorShape + y, rate = priorRate + 1))
    }
  )
)


CRP_conjugate_dgamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape <- nimNumeric(1)
    priorRate <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape <<- model$getParam(marginalizedNodes[1], 'shape') 
      priorRate <<- model$getParam(marginalizedNodes[1], 'rate') 
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      dataMean <- model$getParam(dataNodes[i], 'mean')
      y <- values(model, dataNodes[i])[1]
      return(-0.5*log(2*pi) + priorShape * log(priorRate) - lgamma(priorShape) +
               lgamma(priorShape + 0.5) - (priorShape + 0.5)*log(priorRate + 0.5*(y-dataMean)^2))
    },
    sample = function(i = integer(), j = integer()) {
      dataMean <- model$getParam(dataNodes[i], 'mean')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rgamma(1, shape = priorShape + 0.5, rate = priorRate + (y-dataMean)^2/2))
    }
  )
)



CRP_conjugate_dbeta_dbern <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dbern",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape1 <- nimNumeric(1)
    priorShape2 <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape1 <<- model$getParam(marginalizedNodes[1], 'shape1') 
      priorShape2 <<- model$getParam(marginalizedNodes[1], 'shape2')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      return(lgamma(priorShape1+y) + lgamma(priorShape2+1-y) - lgamma(priorShape1) - lgamma(priorShape2) - log(priorShape1+priorShape2))
    },
    sample = function(i = integer(), j = integer()) {
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rbeta(1, shape1=priorShape1+y, shape2=priorShape2+1-y))
    }
  )
)

CRP_conjugate_dbeta_dbin <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dbin",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape1 <- nimNumeric(1)
    priorShape2 <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape1 <<- model$getParam(marginalizedNodes[1], 'shape1') 
      priorShape2 <<- model$getParam(marginalizedNodes[1], 'shape2')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      dataSize <- model$getParam(dataNodes[i], 'size')
      return(lgamma(priorShape1+priorShape2) + lgamma(priorShape1+y) + lgamma(priorShape1+dataSize-y) -
               lgamma(priorShape1) - lgamma(priorShape2) - lgamma(priorShape1+priorShape1+dataSize) +
               lfactorial(dataSize) - lfactorial(y) - lfactorial(dataSize-y))
    },
    sample = function(i = integer(), j = integer()) {
      dataSize <- model$getParam(dataNodes[i], 'size')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rbeta(1, shape1=priorShape1+y, shape2=priorShape2+dataSize-y))
    }
  )
)



CRP_conjugate_dbeta_dnegbin <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dnegbin",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape1 <- nimNumeric(1)
    priorShape2 <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape1 <<- model$getParam(marginalizedNodes[1], 'shape1') 
      priorShape2 <<- model$getParam(marginalizedNodes[1], 'shape2')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      dataSize <- model$getParam(dataNodes[i], 'size')
      return(lgamma(priorShape1+priorShape2) + lgamma(priorShape1+dataSize) + lgamma(priorShape1+y) -
               lgamma(priorShape1) - lgamma(priorShape2) - lgamma(priorShape1+priorShape1+dataSize+y) +
               lfactorial(y+dataSize-1) - lfactorial(y) - lfactorial(dataSize-1))
    },
    sample = function(i = integer(), j = integer()) {
      dataSize <- model$getParam(dataNodes[i], 'size')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rbeta(1, shape1=priorShape1+dataSize, shape2=priorShape2+y))
    }
  )
)

CRP_conjugate_dgamma_dexp <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dexp",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape <- nimNumeric(1)
    priorRate <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape <<- model$getParam(marginalizedNodes[1], 'shape') 
      priorRate <<- model$getParam(marginalizedNodes[1], 'rate') 
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      return(log(priorShape) + priorShape*log(priorRate) - (priorShape+1)*log(priorRate+y))
    },
    sample = function(i = integer(), j = integer()) {
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rgamma(1, shape=priorShape+1, rate=priorRate+y))
    }
  )
)


CRP_conjugate_dgamma_dgamma <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dgamma",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape <- nimNumeric(1)
    priorRate <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape <<- model$getParam(marginalizedNodes[1], 'shape') 
      priorRate <<- model$getParam(marginalizedNodes[1], 'rate')  
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      datashape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      return((datashape-1)*log(y) + priorShape*log(priorRate) + lgamma(datashape+priorShape) -
               lgamma(datashape) - lgamma(priorShape) -(datashape+priorShape)*log(priorRate+y))
    },
    sample = function(i = integer(), j = integer()) {
      datashape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rgamma(1, shape=datashape+priorShape, rate=priorRate+y))
    }
  )
)


CRP_conjugate_dgamma_dweib <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dweib",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape <- nimNumeric(1)
    priorRate <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape <<- model$getParam(marginalizedNodes[1], 'shape') 
      priorRate <<- model$getParam(marginalizedNodes[1], 'rate')  
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      dataShape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      return( log(dataShape) + (dataShape-1)*log(y) + priorShape*log(priorRate) +
                lgamma(priorShape+1) - lgamma(priorShape) - (priorShape + 1)*log(priorRate + y^dataShape))
    },
    sample = function(i = integer(), j = integer()) {
      dataShape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rgamma(1, shape=1+priorShape, rate=priorRate+y^dataShape))
    }
  )
)


CRP_conjugate_dgamma_dinvgamma <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dinvgamma",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    priorShape <- nimNumeric(1)
    priorRate <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorShape <<- model$getParam(marginalizedNodes[1], 'shape') 
      priorRate <<- model$getParam(marginalizedNodes[1], 'rate')  
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      dataShape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      return( -(dataShape+1)*log(y) + priorShape*log(priorRate) + lgamma(priorShape + dataShape) -
                lgamma(dataShape) - lgamma(priorShape) - (dataShape + priorShape)*log(priorRate + 1/y))
    },
    sample = function(i = integer(), j = integer()) {
      dataShape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- c(rgamma(1, shape=dataShape+priorShape, rate=priorRate+1/y))
    }
  )
)


CRP_conjugate_ddirch_dmulti <- nimbleFunction(
  name = "CRP_conjugate_ddirch_dmulti",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, p, nTilde) {
    d <- length(model[[marginalizedNodes[1]]])
    priorAlpha <- nimNumeric(d)
  },
  methods = list(
    storeParams = function() {
      priorAlpha <<- model$getParam(marginalizedNodes[1], 'alpha')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])
      return(lfactorial(d) - sum(lfactorial(y) + lgamma(priorAlpha)) +
               lgamma(sum(priorAlpha)) + sum(lgamma(priorAlpha+y)) - lgamma(sum(priorAlpha+y)))
    },
    sample = function(i = integer(), j = integer()) {
      y <- values(model, dataNodes[i])
      values(model, marginalizedNodes[j]) <<- rdirch(1, alpha = priorAlpha+y)
    }
  )
)


## general dCRP sampler covering nonconjugate and conjugate cases

#' @rdname samplers
#' @export
sampler_CRP <- nimbleFunction(
  name = 'sampler_CRP',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    targetVar <- model$getVarNames(nodes = target)  
    n <- length(targetElements) 

    ## Find nodes indexed by the CRP node.
    clusterVarInfo <- findClusterNodes(model, target)
    tildeVars <- clusterVarInfo$clusterVars

    ##  Various checks that model structure is consistent with our CRP sampler. ##
    
    if(!is.null(clusterVarInfo$targetNonIndex))
        stop("sampler_CRP: Detected that the CRP variable is used in some way not as an index: ", clusterVarInfo$targetNonIndex, ". NIMBLE's CRP sampling not set up to handle this case.")

    if(any(clusterVarInfo$nTilde == 0)) {
        var <- which(clusterVarInfo$nTilde == 0)
        stop("sampler_CRP: Detected unusual indexing in ", deparse(clusterVarInfo$indexExpr[[var[1]]]),
             " . NIMBLE's CRP MCMC sampling is not designed for this case.")
    }
    
    if(is.null(tildeVars))
        stop('sampler_CRP: The model should have at least one cluster variable.\n')

    ## Avoid case like y[i] ~ dnorm(mu[xi[i]], exp(mu[xi[i]])) as it would complicate p>1 case below.
    ## Also catches cases like:  mu[1] <- muTilde[xi[1]]; mu[2] <- exp(muTilde[xi[2]]) that break conjugacy.
   if(length(tildeVars) != length(unique(tildeVars)))
        stop("sampler_CRP: Cluster membership variable used in multiple declarations. NIMBLE's CRP MCMC sampling not designed for this case.")
    
    ## Cases like 'muTilde[xi[n-i+1]]'. sampler_CRP may be ok with this, but when we wrap the cluster node sampling
    ## to avoid sampling empty clusters, this kind of indexing will cause incorrect behavior.
    if(any(clusterVarInfo$targetIndexedByFunction))
        stop("sampler_CRP: Detected that CRP variable is indexed by a function such as 'mu[xi[n-i+1]]' rather than simple indexing such as 'mu[xi[i]]'. NIMBLE's CRP MCMC sampling not designed for this case.")
      
    ## Check for indexing that would cause errors in using sample() with 'j' as the index of the set of cluster nodes, e.g., mu[xi[i]+2] or mu[some_function(xi[i])]
    ## This case should be correctly handled now by use of findClusterNodes to determine ordering of cluster nodes.  
    #if(any(!clusterVarInfo$targetIsIndex))
    #    stop("sampler_CRP: At the moment, NIMBLE's CRP MCMC sampling requires simple indexing and cannot work with indexing such as 'mu[xi[i]+2]'.\n")

    allTildeNodes <- unlist(clusterVarInfo$clusterNodes)
    dataNodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE) # only data
    stochDepsTildeNodes <- model$getDependencies(allTildeNodes, self = FALSE, stochOnly = TRUE)

    ## Make sure tildeNodes as determined from clustering actually are in model.
    if(!all(allTildeNodes %in% model$getNodeNames())) {
        missingNodes <- allTildeNodes[which(!allTildeNodes %in% model$getNodeNames())]
        stop("sampler_CRP: These cluster parameters are not nodes in the model: ", paste(missingNodes, collapse = ','))
    }
    
    ## Check that no other non-data nodes depend on cluster variables. 
    if(!identical(sort(dataNodes), sort(stochDepsTildeNodes)))
      stop("sampler_CRP: Only the variables being clustered can depend on the cluster parameters.")  

    ## Check that membership variable is independent of cluster nodes.
    ## Should be redundant with check that no other non-data nodes depend on cluster variables.
    if(target %in% stochDepsTildeNodes)
      stop("sampler_CRP: Cluster membership variable has to be independent of cluster parameters.")
    
    ## Check that cluster nodes are independent of membership variable
    ## (dataNodes are the dependents of target).
    ## Should be redundant with check that no other non-data nodes depend on cluster variables.
    if(length(intersect(dataNodes, allTildeNodes)))
        stop("sampler_CRP: Cluster parameters have to be independent of cluster membership variable.")

    conjugacyResult <- checkCRPconjugacy(model, target)

    if(is.null(conjugacyResult) || conjugacyResult != "conjugate_dnorm_invgamma_dnorm") {
        ## Check that all cluster parameters are independent as it it tricky to make sure we
        ## are correctly simulating from the prior otherwise. 
        tmp <- sapply(allTildeNodes, function(x) {
            if(any(allTildeNodes %in% model$getDependencies(x, self = FALSE, stochOnly = TRUE)))
                stop("sampler_CRP: Cluster parameters must be conditionally independent except for conjugate settings.")
        })
    }
    
    ## Check that observations are independent of each other.
    sapply(dataNodes, function(x) {
        if(any(dataNodes %in% model$getDependencies(x, self = FALSE, stochOnly = TRUE)))
            stop("sampler_CRP: Variables being clustered must be conditionally independent.")
    })

    ## Check for one "observation" per random index. We will relax this to handle the 'Quinn' model.
    ## Note that cases like mu[xi[i],xi[j]] are being trapped in findClusterNodes().
    if(n != length(dataNodes))
      stop("sampler_CRP: At the moment, NIMBLE can only sample when there is one variable being clustered for each cluster membershipID.")
    
    nTilde <- clusterVarInfo$nTilde
    if(length(unique(nTilde)) != 1)
      stop('sampler_CRP: In a model with multiple cluster parameters, the number of those parameters must all be the same.\n')

    #### End of checks of model structure. ####

    min_nTilde <- min(nTilde) ## we need a scalar for use in run code, but note that given check above, all nTilde values are the same...
    if(min_nTilde < n)
      warning('sampler_CRP: The number of cluster parameters is less than the number of potential clusters. The MCMC is not strictly valid if it ever proposes more components than cluster parameters exist; NIMBLE will warn you if this occurs.\n')
    
    
    ## Determine if concentration parameter is fixed or random (code similar to the one in sampleDPmeasure function).
    ## This is used if truncated case to tell user if model is proper or not.
    fixedConc <- TRUE
    parentNodesTarget <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes == target]
    for(i in seq_along(candidateParentNodes)){
      tmp <- model$getDependencies(candidateParentNodes[i], self = FALSE) 
      if(sum(tmp == target)) 
        parentNodesTarget <- c(parentNodesTarget, candidateParentNodes[i])
    }
    if(length(parentNodesTarget)) {
      fixedConc <- FALSE
    }
    
    ## Here we try to set up some data structures that allow us to do observation-specific
    ## computation, to save us from computing for all observations when a single cluster membership is being proposed.
    ## At the moment, this is the only way we can easily handle dependencies for multiple node elements in a
    ## 'vectorized' way.
    nInterm <- length(model$getDependencies(targetElements[1], determOnly = TRUE))
    dataNodes <- rep(targetElements[1], n) ## this serves as dummy nodes that may be replaced below
    ## needs to be legitimate nodes because run code sets up calculate even if if() would never cause it to be used
    type <- 'indivCalcs'
    intermNodes <- dataNodes
    intermNodes2 <- dataNodes
    intermNodes3 <- dataNodes
    if(nInterm > 3) {
      type <- "allCalcs"  ## give up and do the inefficient approach
    } else {
      for(i in seq_len(n)) {
        ## We need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.
        stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self = FALSE) 
        detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
        if(length(stochDeps) != 1) 
          stop("sampler_CRP: NIMBLE cannot currently assign a sampler to a dCRP node unless each membership element is associated with a single component of the variable to be clustered.\n")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
        if(length(detDeps) != nInterm) {
          type <- 'allCalcs'  # give up again; should only occur in strange situations
        } else {
          dataNodes[i] <- stochDeps[1]         
          if(nInterm >= 1) 
            intermNodes[i] <- detDeps[1]
          if(nInterm >= 2)
            intermNodes2[i] <- detDeps[2]
          if(nInterm >= 3)
            intermNodes3[i] <- detDeps[3]
        }
      }
    }
   
    helperFunctions <- nimbleFunctionList(CRP_helper)
    
    ## use conjugacy to determine which helper functions to use
    if(is.null(conjugacyResult)) {
      sampler <- 'CRP_nonconjugate'
    } else 
      sampler <- switch(conjugacyResult,
                        conjugate_dnorm_dnorm = 'CRP_conjugate_dnorm_dnorm',
                        conjugate_dnorm_invgamma_dnorm = 'CRP_conjugate_dnorm_invgamma_dnorm',
                        conjugate_dbeta_dbern  = 'CRP_conjugate_dbeta_dbern',
                        conjugate_dbeta_dbin = 'CRP_conjugate_dbeta_dbin',
                        conjugate_dbeta_dnegbin = 'CRP_conjugate_dbeta_dnegbin',
                        conjugate_dgamma_dpois = 'CRP_conjugate_dgamma_dpois',
                        conjugate_dgamma_dexp = 'CRP_conjugate_dgamma_dexp',
                        conjugate_dgamma_dgamma = 'CRP_conjugate_dgamma_dgamma',
                        conjugate_dgamma_dnorm = 'CRP_conjugate_dgamma_dnorm',
                        conjugate_dgamma_dweib = 'CRP_conjugate_dgamma_dweib',
                        conjugate_dgamma_dinvgamma = 'CRP_conjugate_dgamma_dinvgamma',
                        conjugate_ddirch_dmulti = 'CRP_conjugate_ddirch_dmulti',
                        'CRP_nonconjugate')  ## default if we don't have sampler set up for a conjugacy

    p <- length(tildeVars)

    if(p == 2 && sampler == "CRP_conjugate_dnorm_invgamma_dnorm") {
      for(i in seq_along(tildeVars)) {
        if(model$getDistribution(clusterVarInfo$clusterNodes[[i]][1]) == 'dnorm') {
          marginalizedNodes1 <- clusterVarInfo$clusterNodes[[i]]
        } 
        if(model$getDistribution(clusterVarInfo$clusterNodes[[i]][1]) == 'dinvgamma') {
          marginalizedNodes2 <- clusterVarInfo$clusterNodes[[i]]
        }
      }
      helperFunctions[[1]] <- eval(as.name(sampler))(model, marginalizedNodes1, marginalizedNodes2, dataNodes)
    } else {
        ## p and nTilde only needed for non-conjugate currently.
        ## Note that the elements of tildeNodes will be in order such that the first element corresponds to the cluster
        ## obtained when xi[i] = 1, the second when xi[i] = 2, etc.
      helperFunctions[[1]] <- eval(as.name(sampler))(model, unlist(clusterVarInfo$clusterNodes), dataNodes, p, min_nTilde)
    }
      
    curLogProb <- numeric(n)
    
    printMessage <- TRUE
  },
  
  
  run = function() {
    
    conc <- model$getParam(target, 'conc')
    helperFunctions[[1]]$storeParams()
    
    xi <- model[[target]]
    
    ## Find unique values in model[[target]].
    ## We don't relabel the unique values, but we do create each new cluster as the lowest unused positive integer.
    ## k denotes the number of unique labels in xi
    
    xiUniques <- numeric(min_nTilde)
    xiCounts <- numeric(n)
    
    aux <- min(xi):max(xi) 
    k <- 1
    for(i in seq_along(aux)) { 
      nMembers <- sum(aux[i] == xi)
      if(nMembers > 0) {
        xiCounts[aux[i]] <- nMembers
        xiUniques[k] <- aux[i]
        k <- k + 1
      }
    }
    k <- k-1 # number of unique labels in xi
    
    kNew <- 1 # kNew is the new label that can be sampled
    while(xiCounts[kNew] > 0 & kNew < n) { 
      kNew <- kNew + 1
    }
    if( kNew == n & xiCounts[kNew] > 0 ) {  # all clusters are filled (with singletons)
      kNew <- 0
    }
    if(kNew > min_nTilde & min_nTilde < n) {
      if( printMessage ) {
        if(fixedConc) {
          nimCat('CRP_sampler: This MCMC is for a parametric model. The MCMC attempted to use more components than the number of cluster parameters. To have a sampler for a nonparametric model increase the number of cluster parameters.\n')
        } else {
          nimCat('CRP_sampler: This MCMC is not for a proper model. The MCMC attempted to use more components than the number of cluster parameters. Please increase the number of cluster parameters.\n')
        }
      }
      kNew <- 0 
      printMessage <<- FALSE
    }
    
    
    for(i in 1:n) { # updates one cluster membership at the time , i=1,...,n
      
      xi <- model[[target]]
      xiCounts[xi[i]] <- xiCounts[xi[i]] - 1
      
      # Computing sampling probabilities and sampling an index.
      if( xiCounts[xi[i]] == 0 ) { # cluster is a singleton.

        ## First, compute probability of sampling an existing label.
        reorderXiUniques <- numeric(min_nTilde) # here we save reordered version of xiUniques when there is a singleton. This is used later for updating xiUniques if a component is deleted.
        iprob <- 1
        for(j in 1:k) {
          if( xiCounts[xiUniques[j]] >= 1 ) { 
            model[[target]][i] <<- xiUniques[j] 
            if(type == 'indivCalcs') {
              if(nInterm >= 1) model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
            curLogProb[iprob] <<- log(xiCounts[xiUniques[j]]) + model$getLogProb(dataNodes[i])
            reorderXiUniques[iprob] <- xiUniques[j]
            iprob <- iprob + 1
          }
        }
          
        ## Second, compute probability of sampling a new cluster, here, new cluster is the current cluster!
        model[[target]][i] <<- xi[i] # label of new component
        if(sampler == 'CRP_nonconjugate'){ # simulate tildeVars[xi[i]] # do this everytime there is a singleton so we ensure this comes always from the prior
          helperFunctions[[1]]$sample(i, model[[target]][i])
          if(type == 'indivCalcs') {
            if(nInterm >= 1) model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
        }
        curLogProb[k] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i) # probability of sampling a new label, only k components because xi_i is a singleton

        ## Sample new cluster.
        index <- rcat( n=1, exp(curLogProb[1:k]-max(curLogProb[1:k])) )
        if(index == k) {
          newLab <- xi[i] 
          newLabCond <- TRUE
        } else {
          newLab <- reorderXiUniques[index]
          newLabCond <- FALSE
        }
        
      } else { # cluster is not a singleton.
        ## First, compute probability of sampling an existing label.
        for(j in 1:k) { 
          model[[target]][i] <<- xiUniques[j] 
          if(type == 'indivCalcs') {
            if(nInterm >= 1) model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
          curLogProb[j] <<- log(xiCounts[xiUniques[j]]) + model$getLogProb(dataNodes[i]) 
        }
        ## Second, compute probability of sampling a new cluster depending on the value of kNew.       
        if(kNew == 0) { # no new cluster can be created 
          curLogProb[k+1] <<- log(0)  # k+1 <= n always because k==n requires all singletons, handled above
        } else { # a new cluster can be created
          model[[target]][i] <<- kNew
          if(sampler == 'CRP_nonconjugate'){
            helperFunctions[[1]]$sample(i, model[[target]][i])
            if(type == 'indivCalcs') {
              if(nInterm >= 1) model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }
          curLogProb[k+1] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i) # probability of sampling a new label
        }
        
        # sample an index from 1 to (k+1)
        index <- rcat( n=1, exp(curLogProb[1:(k+1)]-max(curLogProb[1:(k+1)])) )
        if(index == (k+1)) {
          newLab <- kNew
          newLabCond <- TRUE
        } else {
          newLab <- xiUniques[index]
          newLabCond <- FALSE
        }
      }
      
      ## Update metadata about clustering.
      model[[target]][i] <<- newLab
      
      if( newLabCond ) { # a component is created. It can really create a new component or keep the current label if xi_i is a singleton
        if(sampler != 'CRP_nonconjugate') { # updating the cluster parameters of the new cluster
          helperFunctions[[1]]$sample(i, model[[target]][i])
        }
        if( xiCounts[xi[i]] != 0) { # a component is really created
          k <- k + 1
          xiUniques[k] <- newLab 
          kNew <- kNew + 1
          mySum <- sum(xi == kNew) 
          while(mySum > 0 & kNew < (n+1)) { # need to make sure don't go beyond length of vector
            kNew <- kNew+1
            mySum <- sum(xi == kNew)
          }
          if(kNew > min_nTilde & min_nTilde < n) {
            if( printMessage ) {
              if(fixedConc) {
                nimCat('CRP_sampler: This MCMC is for a parametric model. The MCMC attempted to use more components than the number of cluster parameters. To have a sampler for a nonparametric model increase the number of cluster parameters.\n')
              } else {
                nimCat('CRP_sampler: This MCMC is not for a proper model. The MCMC attempted to use more components than the number of cluster parameters. Please increase the number of cluster parameters.\n')
              }
            }
            kNew <- 0
            printMessage <<- FALSE
          }
        }
        xiCounts[model[[target]][i]] <- 1
      } else { # an existing label is sampled
        if(sampler == 'CRP_nonconjugate')   # reset to previous marginalized node value
            helperFunctions[[1]]$sample(i, 0)
        if( xiCounts[xi[i]] == 0 ) { # xi_i is a singleton, a component was deleted
          k <- k - 1
          xiUniques <- reorderXiUniques
          if( kNew == 0 ) { # the sampler was not nonparametric or xi=1:n
            kNew <- xi[i] 
          } else { # the sampler was and remains nonparametric.
            if( kNew > xi[i] ) {
              kNew <- xi[i]
            }
          }
        }
        xiCounts[model[[target]][i]] <- xiCounts[model[[target]][i]] + 1
      }
    }

    ## We have updated cluster variables but not all logProb values are up-to-date.  
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( 
    reset = function () {
      printMessage <<- TRUE
    }
  ), where = getLoadingNamespace()
)




#' @rdname samplers
#' @export
sampler_CRP_old <- nimbleFunction(
  name = 'sampler_CRP_old',
  contains=sampler_BASE,
    
  setup=function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    targetVar <- model$getVarNames(nodes = target)  
    n <- length(targetElements) 
    
    # first check that the sampler can be used: we need one observation per random index
    nObs <- length(model$getDependencies(targetElements, stochOnly = TRUE, self = FALSE))
    if(n != nObs)
        stop("sampler_CRP: The length of membership variable and observations has to be the same.\n")
    
    ## finding 'tilde' variables (the parameters that are being clustered):
    tildeVars <- NULL
    itildeVar <- 1
    
    dep <- model$getDependencies(targetElements[1], self = FALSE)
    for(i in seq_along(dep)) { 
      expr <- cc_getNodesInExpr(model$getValueExpr(dep[i])) 
      for(j in seq_along(expr)) {
        ## look for cases like thetatilde[xi[i]] to identify 'xi' and extract 'thetaTilde'
        tmpexpr <- parse(text = expr[j])[[1]]
        if(length(tmpexpr) >= 3 && is.call(tmpexpr) && tmpexpr[[1]] == '[') {   
          foundTarget <- all.vars(tmpexpr[[3]]) == targetVar   
          if( length(foundTarget) > 0 && sum(foundTarget) > 0 ) {
            tildeVars[itildeVar] <- deparse(tmpexpr[[2]])
            itildeVar <- itildeVar+1 
          }
        }
      }
    }
    if(is.null(tildeVars))
      stop('sampler_CRP:  The model should have at least one cluster variable.\n')
    
    nTilde <- sapply(tildeVars, function(x) length(model[[x]]))
    
    if(length(unique(nTilde)) != 1)
      stop('sampler_CRP: In a model with multiple cluster parameters, the number of those parameters must all be the same.\n')
    
    min_nTilde <- min(nTilde) ## we need a scalar for use in run code
    if(min_nTilde < n)
      warning('sampler_CRP: The number of cluster parameters is less than the number of potential clusters. The MCMC is not strictly valid if it ever proposes more components than cluster parameters exist; NIMBLE will warn you if this occurs.\n')
    
    ## Here we try to set up some data structures that allow us to do observation-specific
    ## computation, to save us from computing for all observations when a single cluster membership is being proposed.
    ## At the moment, this is the only way we can easily handle dependencies for multiple node elements in a
    ## 'vectorized' way.
    nInterm <- length(model$getDependencies(targetElements[1], determOnly = TRUE))
    dataNodes <- rep(targetElements[1], n) ## this serves as dummy nodes that may be replaced below
    ## needs to be legitimate nodes because run code sets up calculate even if if() would never cause it to be used

    type <- 'indivCalcs'
    intermNodes <- dataNodes
    intermNodes2 <- dataNodes
    intermNodes3 <- dataNodes
    if(nInterm > 3) {
      type <- "allCalcs"  ## give up and do the inefficient approach
    } else {
      for(i in seq_len(n)) {
        stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self = FALSE) 
        detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
        if(length(stochDeps) != 1) 
          stop("sampler_CRP: Nimble cannot currently assign a sampler to a dCRP node unless each membership element is associated with a single observation.\n")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
        if(length(detDeps) != nInterm) {
          type <- 'allCalcs'  # give up again; should only occur in strange situations
        } else {
          dataNodes[i] <- stochDeps[1]
          
          if(nInterm >= 1) 
            intermNodes[i] <- detDeps[1]
          if(nInterm >= 2)
            intermNodes2[i] <- detDeps[2]
          if(nInterm >= 3)
            intermNodes3[i] <- detDeps[3]
        }
      }
    }
    
    helperFunctions <- nimbleFunctionList(CRP_helper)
    
    ## use conjugacy to determine which helper functions to use
    if(control$useConjugacy)
        conjugacyResult <- checkCRPconjugacy(model, target) else conjugacyResult <- NULL
    
    if(is.null(conjugacyResult)) {
      sampler <- 'CRP_nonconjugate'
    } else 
      sampler <- switch(conjugacyResult,
                        conjugate_dnorm_dnorm = 'CRP_conjugate_dnorm_dnorm',
                        conjugate_dbeta_dbern  = 'CRP_conjugate_dbeta_dbern',
                        conjugate_dbeta_dbin = 'CRP_conjugate_dbeta_dbin',
                        conjugate_dbeta_dnegbin = 'CRP_conjugate_dbeta_dnegbin',
                        conjugate_dgamma_dpois = 'CRP_conjugate_dgamma_dpois',
                        conjugate_dgamma_dexp = 'CRP_conjugate_dgamma_dexp',
                        conjugate_dgamma_dgamma = 'CRP_conjugate_dgamma_dgamma',
                        conjugate_dgamma_dnorm = 'CRP_conjugate_dgamma_dnorm',
                        conjugate_dgamma_dweib = 'CRP_conjugate_dgamma_dweib',
                        conjugate_dgamma_dinvgamma = 'CRP_conjugate_dgamma_dinvgamma',
                        conjugate_ddirch_dmulti = 'CRP_conjugate_ddirch_dmulti',
                        'CRP_nonconjugate')  ## default if we don't have sampler set up for a conjugacy
    ## we use [1] here because the 2nd/3rd args only used for conjugate cases and currently that is only setup for
    ## single parameters
    helperFunctions[[1]] <- eval(as.name(sampler))(model, tildeVars[1], model$expandNodeNames(tildeVars[1]), dataNodes)
    
    curLogProb <- numeric(n) 
  },
  
  
  run = function() {
    conc <- model$getParam(target, 'conc')
    helperFunctions[[1]]$storeParams()
    for(i in 1:n) { # updates one cluster membership at the time , i=1,...,n
      xi <- model[[target]]
      cond <- sum(xi[i]==xi) # if cond=1, xi_i is a singleton
      for(j in 1:n) { # calculate probability of sampling indexes 1,...,n   
        if(i==j) { # index i denotes a new indicator xi[i]
          if(cond>1) { # a new parameter has to be created to calculate the prob
            newind <- 1
            mySum <- sum(xi == newind)
            while(mySum>0 & newind < n) { # need to make sure don't go beyond length of vector
              newind <- newind+1
              mySum <- sum(xi == newind)
            }
            if(newind > min_nTilde) {
                nimCat('CRP_sampler: This MCMC is not fully nonparametric. The MCMC attempted to use more components than the number of cluster parameters.\n')
                newind <- xi[i]
            }
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              if(nInterm >= 1) model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          } else { # we keep the old parameter as the "new" one
            newind <- xi[i]
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              if(nInterm >= 1) model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes)  
          }
          curLogProb[j] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i)
        } else { ## consider making the new membership the same as the membership of other element
          model[[target]][i] <<- model[[target]][j]
          if(type == 'indivCalcs') {
            if(nInterm >= 1) model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
          curLogProb[j] <<- model$getLogProb(dataNodes[i])
        }  
      } 
      
      index <- rcat(n=1, exp(curLogProb-max(curLogProb)))#
      if(index==i) { # creates a new component: one that is not used
          if(newind != xi[i]) {
              model[[target]][i] <<- newind
              helperFunctions[[1]]$sample(i, newind)
          }  ## when newind == xi[i], it means we tried to create a cluster beyond min_nTilde, so don't sample new cluster parameters
      } else{
        model[[target]][i] <<- model[[target]][index]
      } 
    }
    model$calculate(calcNodes)
    
    
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {})
)


findClusterNodes <- function(model, target) {
  targetVar <- model$getVarNames(nodes = target)
  targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
  deps <- model$getDependencies(target, self = FALSE)
  declIDs <- sapply(deps, function(x) model$getDeclID(x))
  uniqueIDs <- unique(declIDs)
  depsByDecl <- lapply(uniqueIDs, function(x) deps[which(x == declIDs)])
  ## Find one example dependency per BUGS declaration for more efficient processing
  exampleDeps <- sapply(depsByDecl, `[`, 1)
  

  ## Set up an evaluation environment in which (xi[1],...,xi[n]) = (1,2,...,n)
  ## first try was: e[[targetVar]] <- seq_along(targetElements)
  ## However in first try, that wouldn't handle xi[3:10] ~ dCRP(), but next construction does.
  e <- list()
  idxExpr <- model$getDeclInfo(target)[[1]]$indexExpr[[1]]
  eval(substitute(`<-`(`[`(e$VAR, IDX), seq_along(targetElements)), list(VAR = targetVar, IDX = idxExpr)))
  
  clusterNodes <- indexExpr <- list()
  clusterVars <- indexPosition <- numIndexes <- targetIsIndex <- targetIndexedByFunction <- NULL
  varIdx <- 0

  targetNonIndex <- NULL
  
  for(idx in seq_along(exampleDeps)) {
    ## Pull out expressions, either as RHS of deterministic or parameters of stochastic
    expr <- cc_getNodesInExpr(model$getValueExpr(exampleDeps[idx]))
    for(j in seq_along(expr)) {
      subExpr <- parse(text = expr[j])[[1]]  # individual parameter of stochastic or RHS of deterministic
      len <- length(subExpr)
      ## Look for target variable within expression, but only when used within index
      if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[' &&
         sum(all.vars(subExpr) == targetVar) && subExpr[[2]] != targetVar) {
        varIdx <- varIdx + 1
        clusterVars <- c(clusterVars, deparse(subExpr[[2]]))
        
        ## Determine which index the target variable occurs in.
        k <- whichIndex <- 3
        foundTarget <- FALSE
        while(k <= len) {
            if(sum(all.vars(subExpr[[k]]) == targetVar)) {
                if(foundTarget) {
                    stop("findClusterNodes: CRP variable used multiple times in ", deparse(subExpr),
                         ". NIMBLE's CRP MCMC sampling not designed for this situation.")
                } else {
                    foundTarget <- TRUE
                    whichIndex <- k
                }
            }
            k <- k+1
        }
        if(!foundTarget) stop("findClusterNodes: conflicting information about presence of CRP variable in expression.")

        declInfo <-  model$getDeclInfo(exampleDeps[idx])[[1]]
        
        ## Determine how target variable enters into cluster node definition
        indexPosition[varIdx] <- whichIndex-2
        numIndexes[varIdx] <- len - 2
        indexExpr[[varIdx]] <- subExpr
        ## Is target used directly as index, e.g., "mu[xi[.]]" as opposed to something like "mu[1+xi[.]]".
        targetIsIndex[varIdx] <- length(subExpr[[whichIndex]]) == 3 &&
          subExpr[[whichIndex]][[1]] == '[' &&
          subExpr[[whichIndex]][[2]] == targetVar
        ## Is indexing of target a simple index, e.g. xi[i], as opposed to something like "xi[n-i+1]".
        targetIndexedByFunction[varIdx] <- any(sapply(declInfo$symbolicParentNodes,
                                                    function(x) 
                                                        length(x) >= 3 && x[[1]] == '[' &&
                                                        x[[2]] == targetVar && length(x[[3]]) > 1))
        
        ## Determine all sets of index values so they can be evaluated in context of possible values of target element values.
        unrolledIndices <- declInfo$unrolledIndicesMatrix
        n <- nrow(unrolledIndices)
        if(n > 0) {  # catch cases like use of xi[2] rather than xi[i]
            clusterNodes[[varIdx]] <- rep(NA, n)
            
            ## Determine unevaluated expression, e.g., muTilde[xi[i],j] not muTilde[xi[1],2]
            expr <- declInfo$valueExprReplaced
            expr <- parse(text = cc_getNodesInExpr(expr)[[j]])[[1]]
            templateExpr <- expr   
            
            ## Now evaluate index values for all possible target element values, e.g.,
            ## xi[i] for all 'i' values with xi taking values 1,...,n
            for(i in seq_len(n)) { 
                for(k in 3:len) # this will deal with muTilde[xi[i], j] type cases
                    if(length(all.vars(expr[[k]])))  # prevents eval of things like 1:3, which the as.numeric would change to c(1,3)
                        templateExpr[[k]] <- as.numeric(eval(substitute(EXPR, list(EXPR = expr[[k]])),
                                                             c(as.list(unrolledIndices[i,]), e)))  # as.numeric avoids 1L, 2L, etc.
                clusterNodes[[varIdx]][i] <- deparse(templateExpr)  # convert to node names
            }
        } else clusterNodes[[varIdx]] <- character(0)
      } 
      if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[' && subExpr[[2]] == targetVar)
          targetNonIndex <- deparse(model$getDeclInfo(exampleDeps[idx])[[1]]$codeReplaced)
    }
  }
  
  ## Find the potential cluster nodes that are actually model nodes,
  ## making sure that what we decide are real cluster nodes are the full potential set
  ## or a truncated set that starts with the first cluster node, e.g., muTilde[1], ..., muTilde[3] is ok;
  ## muTilde[2], ..., muTilde[4] is not (unless the model nodes are muTilde[2], ...., muTilde[4]).
  nTilde <- sapply(clusterNodes, length)
  modelNodes <- model$getNodeNames()
  for(varIdx in seq_along(clusterVars)) {
    if(nTilde[[varIdx]]) {
        if(any(is.na(clusterNodes[[varIdx]])))  
            stop("findClusterNodes: fewer cluster IDs in ", target, " than elements being clustered.")
        if(!all(clusterNodes[[varIdx]] %in% modelNodes)) {  # i.e., truncated representation
            cnt <- nTilde[varIdx]
            while(cnt > 0) {
                ## Try to find first nTilde nodes such that are all actual model nodes.
                if(all(clusterNodes[[varIdx]][seq_len(cnt)] %in% modelNodes)) {
                    nTilde[varIdx] <- cnt
                    clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][seq_len(cnt)]
                    break
                }            
                cnt <- cnt - 1
            }
            if(cnt == 0) {
                warning("findClusterNodes: missing cluster parameter ", clusterNodes[[varIdx]][1], ".")
                clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][clusterNodes[[varIdx]] %in% modelNodes]
                if(!length(clusterNodes[[varIdx]]))
                    stop("findClusterNodes: no cluster parameters for variable ", clusterVars[varIdx], ".")
            }
        }
    }
  }
  return(list(clusterNodes = clusterNodes, clusterVars = clusterVars, nTilde = nTilde,
              targetIsIndex = targetIsIndex, indexPosition = indexPosition, indexExpr = indexExpr,
              numIndexes = numIndexes, targetIndexedByFunction = targetIndexedByFunction,
              targetNonIndex = targetNonIndex))
}

checkNormalInvGammaConjugacy <- function(model, clusterVarInfo) {
    if(length(clusterVarInfo$clusterVars) != 2)
        stop("checkNormalInvGammaConjugacy: requires two cluster variables.")
    conjugate <- FALSE
    
    clusterNodes1 <- clusterVarInfo$clusterNodes[[1]]
    clusterNodes2 <- clusterVarInfo$clusterNodes[[2]]
    dists <- c(model$getDistribution(clusterNodes1[1]), model$getDistribution(clusterNodes2[1]))
    if(dists[1] == "dinvgamma" && dists[2] ==  "dnorm") {  ## put in order so dnorm node is first
        dists <- c("dnorm", "dinvgamma")
        tmp <- clusterNodes1; clusterNodes1 <- clusterNodes2; clusterNodes2 <- tmp
    }

    ## Check conjugacy for example nodes.
    exampleNodes <- c(clusterNodes1[1], clusterNodes2[1])
    if(dists[1] == "dnorm" && dists[2] == "dinvgamma") {
        conjugacy_dnorm <- model$checkConjugacy(exampleNodes[1], restrictLink = 'identity')
        conjugacy_dinvgamma <- model$checkConjugacy(exampleNodes[2])
        if(length(conjugacy_dnorm) && length(conjugacy_dinvgamma) &&
           sum(conjugacy_dinvgamma[[1]]$control$dep_dnorm == exampleNodes[1])) 
            conjugate <- TRUE
    }
    if(conjugate) {  ## Now check that conjugacy holds for all cluster nodes.
        ## Check that distribution for cluster mean nodes are same
        if(any(model$getDistribution(clusterNodes1) != "dnorm"))
            conjugate <- FALSE
        ## Check that mean for cluster mean nodes are same
        meanExprs <- sapply(clusterNodes1, function(x) model$getParamExpr(x, 'mean'))
        names(meanExprs) <- NULL
        if(length(unique(meanExprs)) != 1)
            conjugate <- FALSE
        ## Check that variance for cluster mean nodes are same except for dependence on variance
        varExprs <- sapply(clusterNodes1, function(x) model$getParamExpr(x, 'var'))
        names(varExprs) <- NULL
        for(i in seq_along(varExprs)) {
            varText <- deparse(varExprs[[i]])
            if(length(grep(clusterNodes2[i], varText, fixed = TRUE)))  ## remove clusterNodes2[i] from expression so var expressions will be the same
                varExprs[[i]] <- parse(text = gsub(clusterNodes2[i], "a", varText, fixed = TRUE))[[1]]
        }
        if(length(unique(varExprs)) != 1)
            conjugate <- FALSE

        ## Check that cluster variance nodes are IID
        valueExprs <- sapply(clusterNodes2, function(x) model$getValueExpr(x))
        names(valueExprs) <- NULL
        if(length(unique(valueExprs)) != 1)
            conjugate <- FALSE

        ## Check that dependent nodes ('observations') from same declaration.
        ## This should ensure they have same distribution and parameters are being
        ## clustered in same way.
        depNodes <- model$getDependencies(clusterNodes1, stochOnly = TRUE, self = FALSE)
        if(length(unique(model$getDeclID(depNodes))) != 1)
            conjugate <- FALSE

    }
    return(conjugate)
}

sampler_CRP_cluster_wrapper <- nimbleFunction(
    name = "CRP_cluster_wrapper", 
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        regular_sampler <- nimbleFunctionList(sampler_BASE)
        if(!exists('samplerFunction', control))
            control$samplerFunction <- eval(as.name(paste0("sampler_", control$wrapped_type)))
        regular_sampler[[1]] <- eval(control$samplerFunction)(model, mvSaved, target, control$control)
        dcrpNode <- control$dcrpNode
        clusterID <- control$clusterID
    },
    run = function() {
        if(any(model[[dcrpNode]] == clusterID)) regular_sampler[[1]]$run()
    },
    methods = list(
        reset = function() {regular_sampler[[1]]$reset()}
    ))
